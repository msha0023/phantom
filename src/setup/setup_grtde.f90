!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! Setup for general relativistic tidal disruption event
!
! :References: None
!
! :Owner: David Liptai
!
! :Runtime parameters:
!   - beta          : *penetration factor*
!   - dumpsperorbit : *number of dumps per orbit*
!   - ecc           : *eccentricity (1 for parabolic)*
!   - mhole         : *mass of black hole (solar mass)*
!   - norbits       : *number of orbits*
!   - relax         : *relax star into hydrostatic equilibrium*
!   - theta         : *inclination of orbit (degrees)*
!
! :Dependencies: eos, externalforces, gravwaveutils, infile_utils, io,
!   kernel, metric, mpidomain, options, part, physcon, relaxstar,
!   setbinary, setstar, setup_params, systemutils, timestep, units,
!   vectorutils
!
 use setstar, only:star_t
 implicit none
 public :: setpart

 real    :: mhole,beta,ecc,norbits,theta
 integer :: dumpsperorbit
 logical :: relax,provide_coord,calc_orbit
 type(star_t) :: star
 real    :: x_gr, &
         y_gr, &
         z_gr, &
         vx_gr, &
         vy_gr, &
         vz_gr

 private

contains

!----------------------------------------------------------------
!+
!  setup for sink particle binary simulation (no gas)
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use part,      only:nptmass,xyzmh_ptmass,vxyz_ptmass,ihacc,ihsoft,igas,&
                     gravity,eos_vars,rad,gr
 use setbinary, only:set_binary
 use setstar,   only:set_star,shift_star
 use units,     only:set_units,umass,udist
 use physcon,   only:solarm,pi,solarr
 use io,        only:master,fatal,warning
 use timestep,  only:tmax,dtmax
 use metric,    only:mass1,a,metric_type
 use eos,       only:ieos,X_in,Z_in
 use kernel,    only:hfact_default
 use mpidomain, only:i_belong
 use externalforces, only:accradius1,accradius1_hard
 use vectorutils,    only:rotatevec
 use gravwaveutils,  only:theta_gw,calc_gravitwaves
 use setup_params,   only:rhozero,npart_total
 use systemutils,    only:get_command_option
 use options,        only:iexternalforce
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real,              intent(out)   :: vxyzu(:,:)
 character(len=120) :: filename
 integer :: ierr,np_default,i
 logical :: iexist,write_profile,use_var_comp
 real    :: rtidal,rp,semia,period,hacc1,hacc2
 real    :: vxyzstar(3),xyzstar(3)
 real    :: r0,vel,lorentz
 real    :: vhat(3),x0,y0
!
!-- general parameters
!
 hfact = hfact_default
 time  = 0.
 polyk = 1.e-10    ! <== uconst
 gamma = 5./3.
 ieos  = 2
 if (.not.gravity) call fatal('setup','recompile with GRAVITY=yes')
!
!-- space available for injected gas particles
!
 npart          = 0
 npartoftype(:) = 0
 xyzh(:,:)      = 0.
 vxyzu(:,:)     = 0.
 nptmass        = 0
!
!-- Default runtime parameters
!
 mhole           = 1.e6  ! (solar masses)
 call set_units(mass=mhole*solarm,c=1.d0,G=1.d0) !--Set central mass to M=1 in code units
 star%mstar      = 1.*solarm/umass
 star%rstar      = 1.*solarr/udist
 np_default      = 1e6
 star%np         = int(get_command_option('np',default=np_default)) ! can set default value with --np=1e5 flag (mainly for testsuite)
 star%iprofile   = 2
 beta            = 5.
 ecc             = 0.8
 norbits         = 5.
 dumpsperorbit   = 100
 theta           = 0.
 write_profile   = .false.
 use_var_comp    = .false.
 relax           = .true.
 provide_coord   = .false.
!
!-- Read runtime parameters from setup file
!
 if (id==master) print "(/,65('-'),1(/,a),/,65('-'),/)",' Tidal disruption in GR'
 filename = trim(fileprefix)//'.setup'
 inquire(file=filename,exist=iexist)
 if (iexist) call read_setupfile(filename,ieos,polyk,ierr)
 if (.not. iexist .or. ierr /= 0) then
    if (id==master) then
       call write_setupfile(filename)
       print*,' Edit '//trim(filename)//' and rerun phantomsetup'
    endif
    stop
 endif

 !
 !--set up and relax a star
 !
 call set_star(id,master,star,xyzh,vxyzu,eos_vars,rad,npart,npartoftype,&
               massoftype,hfact,xyzmh_ptmass,vxyz_ptmass,nptmass,ieos,polyk,gamma,&
               X_in,Z_in,relax,use_var_comp,write_profile,&
               rhozero,npart_total,i_belong,ierr)
 if (ierr /= 0) call fatal('setup','errors in set_star')
 !
 !--place star into orbit
 !
 rtidal          = star%rstar*(mass1/star%mstar)**(1./3.)
 rp              = rtidal/beta
 accradius1_hard = 5.*mass1
 accradius1      = accradius1_hard
 a               = 0.
 theta           = theta*pi/180.

 print*, 'mstar', star%mstar
 print*, 'rstar', star%rstar
 print*, 'umass', umass
 print*, 'udist', udist
 print*, 'mass1', mass1
 print*, 'tidal radius', rtidal
 print*, 'beta', beta

 xyzstar  = 0.
 vxyzstar = 0.
 period   = 0.

 if (ecc<1.) then
    !
    !-- Set a binary orbit given the desired orbital parameters to get the position and velocity of the star
    !
    semia    = rp/(1.-ecc)
    period   = 2.*pi*sqrt(semia**3/mass1)
    print*, 'period', period
    hacc1    = star%rstar/1.e8    ! Something small so that set_binary doesnt warn about Roche lobe
    hacc2    = hacc1
    ! apocentre = rp*(1.+ecc)/(1.-ecc)
    ! trueanom = acos((rp*(1.+ecc)/r0 - 1.)/ecc)*180./pi
    call set_binary(mass1,star%mstar,semia,ecc,hacc1,hacc2,xyzmh_ptmass,vxyz_ptmass,nptmass,ierr,&
                    posang_ascnode=0.,arg_peri=90.,incl=0.,f=-180.)
    vxyzstar = vxyz_ptmass(1:3,2)
    xyzstar  = xyzmh_ptmass(1:3,2)
    nptmass  = 0

    call rotatevec(xyzstar,(/0.,1.,0./),-theta)
    call rotatevec(vxyzstar,(/0.,1.,0./),-theta)

 elseif (abs(ecc-1.) < tiny(0.)) then
    
    if (.not. provide_coord) then !
       !-- Setup a parabolic orbit
       !
       r0       = 10.*rtidal              ! A default starting distance from the black hole.
       period   = 2.*pi*sqrt(r0**3/mass1) !period not defined for parabolic orbit, so just need some number
       y0       = -2.*rp + r0
       x0       = sqrt(r0**2 - y0**2)
       xyzstar  = (/-x0,y0,0./)
       vel      = sqrt(2.*mass1/r0)
       vhat     = (/2.*rp,-x0,0./)/sqrt(4.*rp**2 + x0**2)
       vxyzstar = vel*vhat

       call rotatevec(xyzstar,(/0.,1.,0./),theta)
       call rotatevec(vxyzstar,(/0.,1.,0./),theta)
   else 
       xyzstar = (/x_gr, y_gr, z_gr/)
       vxyzstar = (/vx_gr, vy_gr, vz_gr/)

       print*, xyzstar,"position vector"
       print*, vxyzstar,"vel vector"
   endif 
 else
    call fatal('setup','please choose a valid eccentricity (0<ecc<=1)',var='ecc',val=ecc)
 endif

 lorentz = 1./sqrt(1.-dot_product(vxyzstar,vxyzstar))
 if (lorentz>1.1) call warning('setup','Lorentz factor of star greater than 1.1, density may not be correct')

 tmax      = norbits*period
 dtmax     = period/dumpsperorbit

 if (id==master) then
    print "(/,a)",       ' STAR SETUP:'
    print "(a,3f10.3)"  ,'         Position = ',xyzstar
    print "(a,3f10.3)"  ,'         Velocity = ',vxyzstar
    print "(a,1f10.3)"  ,' Lorentz factor   = ',lorentz
    print "(a,1f10.3)"  ,' Polytropic gamma = ',gamma
    print "(a,3f10.3,/)",'       Pericentre = ',rp
 endif
 call shift_star(npart,xyzh,vxyzu,x0=xyzstar,v0=vxyzstar)
 print*,"shifted npart"
 read*
 if (nptmass .ne. 0) then
         print*,"shifting nptmass"
     do i = 1,nptmass
        print*, i, "i", nptmass, "nptmass"
        print*, xyzmh_ptmass(1:3,i), " xyzmh_ptmass(1:3,i)"
        print*,xyzstar,"xyzstar"
        xyzmh_ptmass(1:3,i) = xyzstar + xyzmh_ptmass(1:3,i)
        vxyz_ptmass(1:3,i) = vxyzstar + vxyz_ptmass(1:3,i)
     enddo 
 endif
 
 if (id==master) print "(/,a,i10,/)",' Number of particles setup = ',npart

 !
 ! set a few options for the input file
 !
 calc_gravitwaves = .true.
 if (abs(ecc-1.) > epsilon(0.)) then
    theta_gw = theta*180./pi
 else
    theta_gw = -theta*180./pi
 endif

 if (.not.gr) iexternalforce = 1

 if (npart == 0)   call fatal('setup','no particles setup')
 if (ierr /= 0)    call fatal('setup','ERROR during setup')

end subroutine setpart

!
!---Read/write setup file--------------------------------------------------
!
subroutine write_setupfile(filename)
 use infile_utils, only:write_inopt
 use setstar,      only:write_options_star
 use relaxstar,    only:write_options_relax
 character(len=*), intent(in) :: filename
 integer :: iunit

 print "(a)",' writing setup options file '//trim(filename)
 open(newunit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for tidal disruption setup'
 
 call write_inopt(provide_coord,'provide_coord','Do you want to provite initial conditions?', iunit)
 call write_inopt(calc_orbit,'calc_orbit','Do you want to automatically calculate orbit?', iunit)
 
 if (provide_coord) then  
    write(iunit,"(/,a)") '# options for black hole and orbit'
    call write_inopt(mhole,        'mhole',        'mass of black hole (solar mass)',iunit)
    call write_inopt(ecc,          'ecc',          'eccentricity (1 for parabolic)', iunit)
    call write_grfile(filename,iunit)
    call write_inopt(dumpsperorbit,'dumpsperorbit','number of dumps per orbit',      iunit)
    call write_options_star(star,iunit)
    call write_inopt(relax,'relax','relax star into hydrostatic equilibrium',iunit)
    if (relax) call write_options_relax(iunit)
 endif
 if (calc_orbit) then 
    write(iunit,"(/,a)") '# options for black hole and orbit'
    call write_inopt(mhole,        'mhole',        'mass of black hole (solar mass)',iunit)
    call write_inopt(ecc,          'ecc',          'eccentricity (1 for parabolic)', iunit)
    call write_inopt(beta,         'beta',         'penetration factor',             iunit)
    call write_inopt(norbits,      'norbits',      'number of orbits',               iunit)
    call write_inopt(theta,        'theta',        'inclination of orbit (degrees)', iunit)
    call write_inopt(dumpsperorbit,'dumpsperorbit','number of dumps per orbit',      iunit)
    call write_options_star(star,iunit)
    call write_inopt(relax,'relax','relax star into hydrostatic equilibrium',iunit)
    if (relax) call write_options_relax(iunit)
 endif 
 close(iunit)

end subroutine write_setupfile

subroutine read_setupfile(filename,ieos,polyk,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 use io,           only:error
 use setstar,      only:read_options_star
 use relaxstar,    only:read_options_relax
 use physcon,      only:solarm,solarr
 use units,        only:set_units
 character(len=*), intent(in)    :: filename
 integer,          intent(inout) :: ieos
 real,             intent(inout) :: polyk
 integer,          intent(out)   :: ierr
 integer, parameter :: iunit = 21
 integer :: nerr,need_iso
 type(inopts), allocatable :: db(:)

 print "(a)",'reading setup options from '//trim(filename)
 nerr = 0
 ierr = 0
 call open_db_from_file(db,filename,iunit,ierr)
 !
 !--read black hole mass and use it to define code units
 !
 call read_inopt(mhole,'mhole',db,min=0.,errcount=nerr)
 call set_units(mass=mhole*solarm,c=1.d0,G=1.d0) !--Set central mass to M=1 in code units
 !
 !--read star options and convert to code units
 !
 call read_inopt(provide_coord,'provide_coord',db,errcount=nerr)
 call read_inopt(calc_orbit,'calc_orbit',db,errcount=nerr)

 if (provide_coord) then 
    call read_inopt(mhole,'mhole',db,min=0.,errcount=nerr)
    call set_units(mass=mhole*solarm,c=1.d0,G=1.d0) !--Set central mass to M=1 in code units
    call read_grfile(filename,iunit,nerr,db)
    call read_inopt(dumpsperorbit,  'dumpsperorbit',  db,min=0 ,errcount=nerr)

    call read_options_star(star,need_iso,ieos,polyk,db,nerr)
    if (relax) call read_options_relax(db,nerr)
    call read_inopt(relax,'relax',db,errcount=ierr)
    call read_inopt(ecc,            'ecc',            db,min=0.,max=1.,errcount=nerr) 
 
 endif

 if (calc_orbit) then  
    call read_inopt(mhole,'mhole',db,min=0.,errcount=nerr)
    call set_units(mass=mhole*solarm,c=1.d0,G=1.d0) !--Set central mass to M=1 in code units
    call read_inopt(beta,           'beta',           db,min=0.,errcount=nerr)
    call read_inopt(norbits,        'norbits',        db,min=0.,errcount=nerr)
    call read_inopt(theta,          'theta',          db,       errcount=nerr)
    
    call read_inopt(dumpsperorbit,  'dumpsperorbit',  db,min=0 ,errcount=nerr)

    call read_options_star(star,need_iso,ieos,polyk,db,nerr)
    if (relax) call read_options_relax(db,nerr)
    call read_inopt(relax,'relax',db,errcount=ierr)
    call read_inopt(ecc,            'ecc',            db,min=0.,max=1.,errcount=nerr)
 endif 
 call close_db(db)
 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
    ierr = nerr
 endif

end subroutine read_setupfile

!
!---Read/write setup file fot GR E=0 case--------------------------------------------------
!
subroutine write_grfile(filename,iunit)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer, intent(in) :: iunit 

 print "(a)",' writing E=0 orbit params file '//trim(filename)
 
 write(iunit,"(a)") '# parameters file for a E=0 GR orbit TDE'
 call write_inopt(x_gr,  'x',  'x coordinate (code units)', iunit)
 call write_inopt(y_gr,  'y',  'y coordinate  (code units)',iunit)
 call write_inopt(z_gr,  'z',  'z coordinate  (code units)',iunit)
 call write_inopt(vx_gr, 'vx', 'vx (code units)',iunit)
 call write_inopt(vy_gr, 'vy', 'vy (code units)',iunit)
 call write_inopt(vz_gr, 'vz', 'vz (code units', iunit)

end subroutine write_grfile

subroutine read_grfile(filename,iunit,ierr,db)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 use io,           only:error
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, intent(in) :: iunit 
 type(inopts), allocatable,intent(inout) :: db(:)
 integer :: nerr

 print "(a)",'reading setup options from '//trim(filename)
 nerr = 0
 ierr = 0
 call read_inopt(x_gr,   'x',   db,min=-100000.,errcount=nerr)
 call read_inopt(y_gr,    'y',     db,min=-100000.,errcount=nerr)
 call read_inopt(z_gr,     'z',     db,min=-10000.,errcount=nerr)
 call read_inopt(vx_gr,     'vx',     db,min=-100.,errcount=nerr)
 call read_inopt(vy_gr,  'vy',  db,min=-100.,errcount=nerr)
 call read_inopt(vz_gr,     'vz',     db,min=-100.,errcount=nerr)

end subroutine read_grfile

!subroutine write_orbit_params(filename, iunit)
 !use infile_utils, only:write_inopt
 !character(len=*), intent(in) :: filename
 !integer, intent(in) :: iunit
 
 !call write_inopt(beta,         'beta',         'penetration factor',             iunit)
! call write_inopt(norbits,      'norbits',      'number of orbits',               iunit)
! call write_inopt(theta,        'theta',        'inclination of orbit (degrees)', iunit)
  
!end subroutine write_orbit_params

!subroutine read_orbit_params(filename, iunit, ierr, db)
 !use infile_units, only:open_db_from_file,inopts,read_inopt,close_db
 !use io,           only:error
 !character(len=*), intent(in)  :: filename
 !integer,          intent(out) :: ierr
 !integer, intent(in) :: iunit
 !type(inopts), allocatable,intent(inout) :: db(:)
 !integer :: nerr

 !print "(a)",'reading setup options from '//trim(filename)
 !nerr = 0
 !ierr = 0

 !call read_inopt(beta,           'beta',           db,min=0.,errcount=nerr)
 !call read_inopt(norbits,        'norbits',        db,min=0.,errcount=nerr)
 !call read_inopt(theta,          'theta',          db,       errcount=nerr)

!end subroutine read_orbit_params



end module setup
