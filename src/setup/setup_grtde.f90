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
 integer :: nstar
 real    :: mhole,beta,ecc,norbits,theta
 integer :: dumpsperorbit
 type(star_t) :: star(2)
 real    :: a_binary,ecc_binary,inc_binary,O_binary,w_binary,f_binary,m_binary
 real*8  :: x_pos1,y_pos1,z_pos1,vx_pos1,vy_pos1,vz_pos1,x_pos2,y_pos2,z_pos2,vx_pos2,vy_pos2,vz_pos2
 logical :: relax,corotate_binary
 logical :: use_set_binary
 logical :: use_MC_code
 real    :: rp_outer
 private

contains

!----------------------------------------------------------------
!+
!  setup for sink particle binary simulation (no gas)
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use part,      only:nptmass,xyzmh_ptmass,vxyz_ptmass,ihacc,ihsoft,igas,&
                     gravity,eos_vars,rad,nsinkproperties
 use setbinary, only:set_binary
 use setstar,   only:set_star,shift_star
 use units,     only:set_units,umass,udist,unit_density,utime
 use physcon,   only:solarm,pi,solarr
 use io,        only:master,fatal,warning
 use options,   only:iexternalforce,damp
 use timestep,  only:tmax,dtmax
! use metric,    only:mass1,a
 use metric, only:a
 use orbits_data,    only:isco_kerr
 use eos,       only:ieos,X_in,Z_in
 use kernel,    only:hfact_default
 use mpidomain, only:i_belong
 use externalforces, only:accradius1,accradius1_hard,iext_corotate,mass1
 use vectorutils,    only:rotatevec
 use gravwaveutils,  only:theta_gw,calc_gravitwaves
 use setup_params,   only:rhozero,npart_total
 use centreofmass,   only:reset_centreofmass
 use dim,            only:gr
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
 integer :: ierr,nptmass_in,iextern_prev,i,np_default
 logical :: iexist,write_profile,use_var_comp,add_spin
 real    :: rtidal,rp,semia,period,hacc1,hacc2
 real    :: vxyzstar(3),xyzstar(3)
 real    :: r0,vel,lorentz
 real    :: vhat(3),x0,y0
 real :: xyzmh_ptmass_in(nsinkproperties,2),vxyz_ptmass_in(3,2)
!
!-- general parameters
!
 hfact = hfact_default
 time  = 0.
 polyk = 1.e-10    ! <== uconst
 gamma = 5./3.
 ieos  = 2
 nptmass_in = 0
 iextern_prev = iexternalforce
 iexternalforce = 0
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
 if (gr) then
    call set_units(mass=mhole*solarm,c=1.d0,G=1.d0) !--Set central mass to M=1 in code units
 else
    call set_units(mass=solarm,dist=solarr,G=1.d0)
 endif
 star%mstar      = 1.*solarm/umass
 star%rstar      = 1.*solarr/udist
 np_default      = 1e6
 star%np         = int(get_command_option('np',default=np_default)) ! can set default value with --np=1e5 flag (mainly for testsuite)
 star%iprofile   = 2
 print*,udist,"udist",umass,"umass after setting the units"
 a_binary        = 10.*solarr/udist
 ecc_binary      = 0.
 inc_binary      = 0.
 O_binary        = 0.
 w_binary        = 270.
 f_binary        = 180.
 m_binary        = 180.
 beta            = 5.
 ecc             = 0.8
 norbits         = 5.
 dumpsperorbit   = 100
 theta           = 0.
 write_profile   = .false.
 use_var_comp    = .false.
 relax           = .true.
 use_set_binary  = .false.
 use_MC_code     = .false.
!
!-- Read runtime parameters from setup file
!
 if (id==master) print "(/,65('-'),1(/,a),/,65('-'),/)",' Tidal disruption in GR'
 filename = trim(fileprefix)//'.setup'
 inquire(file=filename,exist=iexist)
 if (iexist) call read_setupfile(filename,ieos,polyk,mass1,ierr)
 if (.not. iexist .or. ierr /= 0) then
    if (id==master) then
       call write_setupfile(filename)
       print*,' Edit '//trim(filename)//' and rerun phantomsetup'
    endif
    stop
 endif
 print*,mhole,"MASS OF SMBH",umass,"umass",udist,"udist",solarm,"solarm",unit_density,"unit density",mass1,"mass1"
 a_binary        = a_binary*solarr/udist
 !
 !--set up and relax a star.
 !
 do i=1,nstar
    if (star(i)%iprofile > 0) then
       if (i == 1) then
          print "(/,a,i0,a)",' --- ONLY STAR 1 TO SET ---'
          call set_star(id,master,star(1),xyzh,vxyzu,eos_vars,rad,npart,npartoftype,&
                   massoftype,hfact,xyzmh_ptmass,vxyz_ptmass,nptmass,ieos,polyk,gamma,&
                   X_in,Z_in,relax,use_var_comp,write_profile,&
                   rhozero,npart_total,i_belong,ierr)
          if (ierr /= 0) call fatal('setup','errors in set_star')
       else
          print "(/,a,i0,a)",' --- STAR ',i,' ---'
          call set_star(id,master,star(i),xyzh,vxyzu,eos_vars,rad,npart,npartoftype,&
                     massoftype,hfact,xyzmh_ptmass,vxyz_ptmass,nptmass,ieos,polyk,gamma,&
                     X_in,Z_in,relax,use_var_comp,write_profile,&
                     rhozero,npart_total,i_belong,ierr,itype=i)
          if (ierr /= 0) call fatal('setup','errors in set_star')
       endif
    else
    endif
 enddo

 ! If use_set_binary then we would have had read the orbital paramaters which will be used here
 if (use_set_binary) then
    !
    !--place binary stars around each other
    !
    !--now setup orbit using fake sink particles
    !
    nptmass_in = 0
    star(1)%hacc = 0.5
    star(2)%hacc = 0.5
    print*,a_binary,"a_binary",ecc_binary,"ecc_binary",star(1)%hacc,"star(1)%hacc"
    call set_binary(star(1)%mstar,star(2)%mstar,a_binary,ecc_binary,star(1)%hacc,star(2)%hacc,&
                  xyzmh_ptmass_in,vxyz_ptmass_in,nptmass_in,ierr,&
                  posang_ascnode=O_binary,arg_peri=w_binary,incl=inc_binary,f=f_binary,verbose=(id==master))
    add_spin = corotate_binary
    print*,xyzmh_ptmass_in*udist,"xyzmh_ptmass_in after setting stars around one another"
    print*,"Stars have been set in the binary system",star(1)%mstar,"star(1)%mstar",star(1)%hacc,"star(1)%hacc"
    if (ierr /= 0) call fatal ('setup_binary','error in call to set_binary')
 endif

 ! Add sink particles if we are using iprofile of 0
 do i=1,nstar
    if (star(i)%iprofile == 0) then
       nptmass = nptmass + 1
       print*,"Sinks added",nptmass,"nptmass"
    endif
 enddo
 
 ! Now lets add code for sink
 !
 !--place star into orbit
 !
 ! We only set the orbit this way when we are not using Monte Carlo Code parameters
 if (.not. use_MC_code) then
    rtidal          = star(1)%rstar*(mass1/star(1)%mstar)**(1./3.)
    rp              = rtidal/beta
    if (gr) then 
       accradius1_hard = 5.*mass1
       accradius1      = accradius1_hard
    else 
       accradius1_hard = 100
       accradius1      = accradius1_hard
    endif
    a               = 0.
    theta           = theta*pi/180.

    print*, 'mstar', star(1)%mstar
    print*, 'rstar', star(1)%rstar
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
       print*, 'period in sec', period*utime
       if (gr) then 
          hacc1    = star(1)%rstar/1.e8    ! Something small so that set_binary doesnt warn about Roche lobe
          hacc2    = hacc1
       else
          hacc1 = star(1)%rstar/2.122e8
          hacc2 = hacc1
       endif 
       ! apocentre = rp*(1.+ecc)/(1.-ecc)
       ! trueanom = acos((rp*(1.+ecc)/r0 - 1.)/ecc)*180./pi
       if (nstar == 1) then
          call set_binary(mass1,star(1)%mstar,semia,ecc,hacc1,hacc2,xyzmh_ptmass,vxyz_ptmass,nptmass,ierr,&
                       posang_ascnode=0.,arg_peri=90.,incl=0.,f=-180.)
       else
          call set_binary(mass1,star(1)%mstar,semia,ecc,hacc1,hacc2,xyzmh_ptmass,vxyz_ptmass,nptmass,ierr,&
                       posang_ascnode=0.,arg_peri=90.,incl=0.,f=-180.)
       endif
       vxyzstar = vxyz_ptmass(1:3,2)
       xyzstar  = xyzmh_ptmass(1:3,2)
       nptmass  = 0
       print*,mass1,"mass1",star(1)%mstar,"star(1)%mstar",vxyzstar,"vxyzstar",xyzstar,"xyzstar"
       call rotatevec(xyzstar,(/0.,1.,0./),-theta)
       call rotatevec(vxyzstar,(/0.,1.,0./),-theta)

    elseif (abs(ecc-1.) < tiny(0.)) then
       !
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
       call fatal('setup','please choose a valid eccentricity (0<ecc<=1)',var='ecc',val=ecc)
    endif
    tmax      = norbits*period
    dtmax     = period/dumpsperorbit
 endif

 ! Here we shift the stars around the SMBH
 if (use_MC_code) then
    do i=1,nstar
       if (i == 1) then
          xyzstar = (/x_pos1,y_pos1,z_pos1/)
          vxyzstar = (/vx_pos1,vy_pos1,vz_pos1/)
       else
          xyzstar = (/x_pos2,y_pos2,z_pos2/)
          vxyzstar = (/vx_pos2,vy_pos2,vz_pos2/)
       endif

       ! if no sink are added we do the following
       if (nptmass == 0) then
          print*,"!!!!!!!!!!!!!!!1xyzmh_ptmass(1:3,i)",xyzmh_ptmass(1:3,i)
          print*,xyzmh_ptmass(1:3,i)+xyzstar,"xyzmh_ptmass(1:3,i)+xyzstar Lets shift the star!"
          call shift_star(npart,xyzh,vxyzu,x0=xyzmh_ptmass(1:3,i)+xyzstar,&
                  v0=vxyz_ptmass(1:3,i)+vxyzstar,itype=i)
       else
          print*,star(i)%mstar,"star(i)%mstar",star(i)%hacc,"star(i)%hacc"
          xyzmh_ptmass(4,i) = star(i)%mstar
          xyzmh_ptmass(5,i) = star(i)%hacc
          xyzmh_ptmass(1:3,i) = xyzstar
          vxyz_ptmass(1:3,i) = vxyzstar
          print*,"For use_mc_code we have set the sink ptmass"
       endif
    enddo
    tmax      = 1e4
    dtmax     = 1e2
 else
    do i=1,nstar
       if (nptmass == 0) then
          if (nstar == 1) then
             call shift_star(npart,xyzh,vxyzu,x0=xyzstar,v0=vxyzstar)
          else
             print*,xyzmh_ptmass_in(1:3,i),"xyzmh_ptmass_in(1:3,i) shifitng"
             print*,xyzmh_ptmass_in(1:3,i)*udist,"xyzmh_ptmass_in(1:3,i)",xyzstar*udist,"xyzstar"
             call shift_star(npart,xyzh,vxyzu,x0=xyzmh_ptmass_in(1:3,i)+xyzstar,&
                v0=vxyz_ptmass_in(1:3,i)+vxyzstar,itype=i)
             print*,"Shfting the stars"
          endif
       else
          print*,"USING USE_SET_BINARY and sinks",xyzstar,"xyzstar",xyzmh_ptmass_in(1:3,i),"xyzmh_ptmass_in(1:3,i)"
          xyzmh_ptmass(4,i) = star(i)%mstar
          xyzmh_ptmass(5,i) = star(i)%hacc
          xyzmh_ptmass(:,i) = xyzmh_ptmass_in(:,i)
          vxyz_ptmass(:,i) = vxyz_ptmass_in(:,i)
          xyzmh_ptmass(1:3,i) = xyzstar + xyzmh_ptmass_in(1:3,i)
          vxyz_ptmass(1:3,i) = vxyzstar + vxyz_ptmass_in(1:3,i)
       endif
    enddo
 endif
 
 if (.not. gr) then
    iexternalforce = 1
    damp = 0
    print*,mass1,"mass1",mhole,"mhole"
 endif

 if (gr) then 
    lorentz = 1./sqrt(1.-dot_product(vxyzstar,vxyzstar))
    if (lorentz>1.1) call warning('setup','Lorentz factor of star greater than 1.1, density may not be correct')
 endif
 if (id==master) then
    print "(/,a)",       ' STAR SETUP:'
    print "(a,3f10.3)"  ,'         Position = ',xyzstar
    print "(a,3f10.3)"  ,'         Velocity = ',vxyzstar
    print "(a,1f10.3)"  ,' Lorentz factor   = ',lorentz
    print "(a,1f10.3)"  ,' Polytropic gamma = ',gamma
    print "(a,3f10.3,/)",'       Pericentre = ',rp
 endif

 if (id==master) print "(/,a,i10,/)",' Number of particles setup = ',npart
 if (.not. gr .and. iexternalforce == 0) call fatal('setup','no gravitational potential set')
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

 call write_inopt(nstar, 'nstar', 'number of stars you want?', iunit)
 if (nstar == 1) then
     call write_options_star(star(1),iunit)
     call write_inopt(relax,'relax','relax star into hydrostatic equilibrium',iunit)
     if (relax) call write_options_relax(iunit)
 endif
 call write_inopt(use_MC_code,'use_MC_code','Do you want to provide pos and vel from Alex MC code?',iunit)
 if (nstar == 2) then
     call write_inopt(use_set_binary,'use_set_binary','Do you want to use orbital params for binary?',iunit)

     call write_options_star(star(1),iunit,label='1')
     call write_options_star(star(2),iunit,label='2')

     if (any(star(:)%iprofile > 0)) then
        write(iunit,"(/,a)") '# relaxation options'
        call write_inopt(relax,'relax','relax stars into equilibrium',iunit)
        call write_options_relax(iunit)
     endif

     if (use_set_binary) call write_setupfile_binary(iunit)

 endif
 if (use_MC_code) call write_binary_options(iunit)
 if ( use_set_binary) then

    write(iunit,"(/,a)") '# options for black hole and orbit'
    call write_inopt(mhole,        'mhole',        'mass of black hole (solar mass)',iunit)
    call write_inopt(beta,         'beta',         'penetration factor',             iunit)
    call write_inopt(ecc,          'ecc',          'eccentricity (1 for parabolic)', iunit)
    call write_inopt(norbits,      'norbits',      'number of orbits',               iunit)
    call write_inopt(dumpsperorbit,'dumpsperorbit','number of dumps per orbit',      iunit)
    call write_inopt(theta,        'theta',        'inclination of orbit (degrees)', iunit)
 endif

 close(iunit)

end subroutine write_setupfile

subroutine read_setupfile(filename,ieos,polyk,mass1,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 use io,           only:error
 use setstar,      only:read_options_star
 use relaxstar,    only:read_options_relax
 use physcon,      only:solarm,solarr
 use units,        only:set_units,umass
 use dim,          only:gr
 character(len=*), intent(in)    :: filename
 integer,          intent(inout) :: ieos
 real,             intent(inout) :: polyk
 real,             intent(out)   :: mass1
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
 !call set_units(mass=mhole*solarm,c=1.d0,G=1.d0) !--Set central mass to M=1 in code units
 if (gr) then 
    mass1 = mhole*solarm/umass
 else
    mass1 = mhole*solarm/umass
    print*,solarm,"solarm",umass,"umass"
 endif
 print*,"mass of BH for metric in code units = ", mass1
 !
 !--read star options and convert to code units
 !
 call read_inopt(nstar,'nstar',db,min=0,errcount=nerr)
  call read_inopt(use_MC_code,'use_MC_code',db,errcount=nerr)
  if (nstar==1) then
    call read_options_star(star(1),need_iso,ieos,polyk,db,nerr)
    call read_inopt(relax,'relax',db,errcount=nerr)
    if (relax) call read_options_relax(db,nerr)
  else
    call read_inopt(use_set_binary,'use_set_binary',db,errcount=nerr)
    call read_options_star(star(1),need_iso,ieos,polyk,db,nerr,label='1')
    call read_options_star(star(2),need_iso,ieos,polyk,db,nerr,label='2')


    if (any(star(:)%iprofile > 0)) then
       call read_inopt(relax,'relax',db,errcount=nerr)
       call read_options_relax(db,nerr)
    endif
    if (use_set_binary) call read_setupfile_binary(ierr,db,iunit)
 endif
 if (use_MC_code) call read_binary_options(ierr,db,iunit)
 if (use_set_binary) then
    call read_inopt(beta,           'beta',           db,min=0.,errcount=nerr)
    call read_inopt(ecc,            'ecc',            db,min=0.,max=1.,errcount=nerr)
    call read_inopt(norbits,        'norbits',        db,min=0.,errcount=nerr)
    call read_inopt(dumpsperorbit,  'dumpsperorbit',  db,min=0 ,errcount=nerr)
    call read_inopt(theta,          'theta',          db,       errcount=nerr)
 endif

 call close_db(db)

 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
    ierr = nerr
 endif

end subroutine read_setupfile


!----------------------------------------------------------------
!+
!  write options for the binary orbit
!+
!----------------------------------------------------------------
subroutine write_setupfile_binary(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# binary orbit settings'
 call write_inopt(a_binary,'a','semi-major axis (e.g. 1 au), period (e.g. 10*days) or rp if e=1',iunit)
 call write_inopt(ecc_binary,'ecc','eccentricity',iunit)
 call write_inopt(inc_binary,'inc','inclination (deg)',iunit)
 call write_inopt(O_binary,'O','position angle of ascending node (deg)',iunit)
 call write_inopt(w_binary,'w','argument of periapsis (deg)',iunit)
 call write_inopt(f_binary,'f','initial true anomaly (180=apoastron)',iunit)
 call write_inopt(m_binary,'M','initial mean anomaly (180=apocentre)',iunit)
 call write_inopt(corotate_binary,'corotate','set stars in corotation',iunit)

end subroutine write_setupfile_binary

!----------------------------------------------------------------
!+
!  read options from .setup file
!+
!----------------------------------------------------------------
subroutine read_setupfile_binary(ierr,db,iunit)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 use io,           only:error,fatal
 integer,          intent(inout) :: ierr
 integer, intent(in) :: iunit
 integer :: nerr,need_iso
 type(inopts), allocatable, intent(inout) :: db(:)
 call read_inopt(a_binary,'a',db,errcount=nerr)
 call read_inopt(ecc_binary,'ecc',db,min=0.,errcount=nerr)
 call read_inopt(inc_binary,'inc',db,errcount=nerr)
 call read_inopt(O_binary,'O',db,errcount=nerr)
 call read_inopt(w_binary,'w',db,errcount=nerr)
 call read_inopt(f_binary,'f',db,errcount=nerr)
 call read_inopt(m_binary,'M',db,errcount=nerr)
 call read_inopt(corotate_binary,'corotate',db,errcount=nerr)

end subroutine read_setupfile_binary
!----------------------------------------------------------------
!+
!  write options for the binary orbit of stars by providing pos
!  and velocities based on Alex's code around the black hole
!+
!----------------------------------------------------------------
subroutine write_binary_options(iunit)
   use infile_utils, only:write_inopt
   integer, intent(in) :: iunit

   write(iunit,"(/,a)") '# options for black hole and orbit'
   call write_inopt(mhole,        'mhole',        'mass of black hole (solar units)',iunit)
   call write_inopt(x_pos1,        'x_pos1',        'x-coordinate of star 1',iunit)
   call write_inopt(y_pos1,        'y_pos1',        'y-coordinate of star 1',iunit)
   call write_inopt(z_pos1,        'z_pos1',        'z-coordinate of star 1',iunit)
   call write_inopt(vx_pos1,        'vx_pos1',        'vx-coordinate of star 1',iunit)
   call write_inopt(vy_pos1,        'vy_pos1',        'vy-coordinate of star 1',iunit)
   call write_inopt(vz_pos1,        'vz_pos1',        'vz-coordinate of star 1',iunit)
   call write_inopt(x_pos2,        'x_pos2',        'x-coordinate of star 2',iunit)
   call write_inopt(y_pos2,        'y_pos2',        'y-coordinate of star 2',iunit)
   call write_inopt(z_pos2,        'z_pos2',        'z-coordinate of star 2',iunit)
   call write_inopt(vx_pos2,        'vx_pos2',        'vx-coordinate of star 2',iunit)
   call write_inopt(vy_pos2,        'vy_pos2',        'vy-coordinate of star 2',iunit)
   call write_inopt(vz_pos2,        'vz_pos2',        'vz-coordinate of star 2',iunit)

 end subroutine  write_binary_options

 subroutine read_binary_options(ierr,db,iunit)
   use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
   use io,           only:error,fatal
   integer,          intent(inout) :: ierr
   integer, intent(in) :: iunit
   integer :: nerr,need_iso
   type(inopts), allocatable, intent(inout) :: db(:)
   call read_inopt(x_pos1,            'x_pos1',           db,errcount=nerr)
   call read_inopt(y_pos1,            'y_pos1',            db,errcount=nerr)
   call read_inopt(z_pos1,            'z_pos1',        db,errcount=nerr)
   call read_inopt(vx_pos1,            'vx_pos1',           db,errcount=nerr)
   call read_inopt(vy_pos1,            'vy_pos1',            db,errcount=nerr)
   call read_inopt(vz_pos1,        'vz_pos1',        db,errcount=nerr)
   call read_inopt(x_pos2,            'x_pos2',           db,errcount=nerr)
   call read_inopt(y_pos2,            'y_pos2',            db,errcount=nerr)
   call read_inopt(z_pos2,            'z_pos2',        db,errcount=nerr)
   call read_inopt(vx_pos2,            'vx_pos2',           db,errcount=nerr)
   call read_inopt(vy_pos2,            'vy_pos2',            db,errcount=nerr)
   call read_inopt(vz_pos2,        'vz_pos2',        db,errcount=nerr)

 end subroutine read_binary_options

end module setup
