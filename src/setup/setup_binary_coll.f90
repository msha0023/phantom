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
!   kernel, metric, mpidomain, part, physcon, relaxstar, setbinary,
!   setstar, setup_params, timestep, units, vectorutils
!
 use setstar, only:star_t
 implicit none
 public :: setpart
 integer :: nstar
 real    :: mhole,beta,ecc,norbits,theta
 integer :: dumpsperorbit
 type(star_t) :: star(2)
 real*8 :: x_pos1,y_pos1,z_pos1,vx_pos1,vy_pos1,vz_pos1,x_pos2,y_pos2,z_pos2,vx_pos2,vy_pos2,vz_pos2
 logical :: relax,corotate_binary
 logical :: provide_rp
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
 use units,     only:set_units,umass,udist,unit_density,unit_velocity
 use physcon,   only:solarm,pi,solarr
 use io,        only:master,fatal,warning
 use options,   only:iexternalforce
 use timestep,  only:tmax,dtmax
 use metric,    only:mass1,a
 use eos,       only:ieos,X_in,Z_in
 use kernel,    only:hfact_default
 use mpidomain, only:i_belong
 use externalforces, only:accradius1,accradius1_hard,iext_corotate
 use vectorutils,    only:rotatevec
 use gravwaveutils,  only:theta_gw,calc_gravitwaves
 use setup_params,   only:rhozero,npart_total
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
 integer :: ierr,nptmass_in,iextern_prev,i
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
 call set_units(mass=mhole*solarm,c=1.d0,G=1.d0) !--Set central mass to M=1 in code units
 print*,umass,"umass",1e6*solarm,"1e6*solarm"
 star%mstar      = 1.*solarm/umass
 star%rstar      = 1.*solarr/udist
 star%np         = 1e6
 star%iprofile   = 2
 beta            = 5.
 ecc             = 0.8
 norbits         = 5.
 dumpsperorbit   = 100
 theta           = 0.
 write_profile   = .false.
 use_var_comp    = .false.
 relax           = .false.
 provide_rp      = .false.
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
 print*,mhole,"MASS OF SMBH",umass,"umass",udist,"udist",solarm,"solarm",unit_density,"unit density"
 !
 !--set up and relax a star
 !

 if (nstar == 1) then
         print*,"--------"
         print*,"Setting single star and relaxing it"
 call set_star(id,master,star(1),xyzh,vxyzu,eos_vars,rad,npart,npartoftype,&
               massoftype,hfact,xyzmh_ptmass,vxyz_ptmass,nptmass,ieos,polyk,gamma,&
               X_in,Z_in,relax,use_var_comp,write_profile,&
               rhozero,npart_total,i_belong,ierr)
 endif

 if (nstar == 2) then 
   do i=1,nstar
    if (star(i)%iprofile > 0) then
       print "(/,a,i0,a)",' --- STAR ',i,' ---'
       call set_star(id,master,star(i),xyzh,vxyzu,eos_vars,rad,npart,npartoftype,&
                     massoftype,hfact,xyzmh_ptmass,vxyz_ptmass,nptmass,ieos,polyk,gamma,&
                     X_in,Z_in,relax,use_var_comp,write_profile,&
                     rhozero,npart_total,i_belong,ierr,itype=i)
    endif
 enddo

 endif
 if (ierr /= 0) call fatal('setup','errors in set_star')

  ! use the pos and vel coordinates to shift the stars around the black hole
  do i=1,nstar 
     if (i == 1) then 
         xyzstar = (/x_pos1,y_pos1,z_pos1/)
         vxyzstar = (/vx_pos1,vy_pos1,vz_pos1/)
     else
         xyzstar = (/x_pos2,y_pos2,z_pos2/)
         vxyzstar = (/vx_pos2,vy_pos2,vz_pos2/)
     endif 
     print*,xyzmh_ptmass(1:3,i)+xyzstar,"xyzmh_ptmass(1:3,i)+xyzstar Lets shift the star!"
     call shift_star(npart,xyzh,vxyzu,x0=xyzmh_ptmass(1:3,i)+xyzstar,&
                       v0=vxyz_ptmass(1:3,i)+vxyzstar,itype=i)
  enddo 
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
 
 if (nstar == 2) then
     call write_options_star(star(1),iunit,label='1')
     !call write_inopt(relax,'relax','relax star into hydrostatic equilibrium',iunit)
     !if (relax) call write_options_relax(iunit)

     call write_options_star(star(2),iunit,label='2')
     !call write_inopt(relax,'relax','relax star into hydrostatic equilibrium',iunit)
     !if (relax) call write_options_relax(iunit)
     if (any(star(:)%iprofile > 0)) then
        write(iunit,"(/,a)") '# relaxation options'
        call write_inopt(relax,'relax','relax stars into equilibrium',iunit)
        call write_options_relax(iunit)
     endif

 endif
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
 close(iunit)

end subroutine write_setupfile

subroutine read_setupfile(filename,ieos,polyk,mass1,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 use io,           only:error
 use setstar,      only:read_options_star
 use relaxstar,    only:read_options_relax
 use physcon,      only:solarm,solarr
 use units,        only:set_units,umass,udist,unit_velocity
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
 ! call set_units(mass=mhole*solarm,c=1.d0,G=1.d0) !--Set central mass to M=1 in code units
 mass1 = mhole*solarm/umass
 print*,"mass of BH for metric in code units = ", mass1
 !
 !--read star options and convert to code units
 !
 call read_inopt(nstar,'nstar',db,min=0,errcount=nerr)

 if (nstar == 1) then 
    call read_options_star(star(1),need_iso,ieos,polyk,db,nerr)
    call read_inopt(relax,'relax',db,errcount=nerr)
    if (relax) call read_options_relax(db,nerr)
 endif 

 if (nstar == 2) then
    call read_options_star(star(1),need_iso,ieos,polyk,db,nerr,label='1')
    call read_options_star(star(2),need_iso,ieos,polyk,db,nerr,label='2')
    if (any(star(:)%iprofile > 0)) then
       call read_inopt(relax,'relax',db,errcount=nerr)
       call read_options_relax(db,nerr)
    endif
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
 call close_db(db)
  
 print*,x_pos1, y_pos1,z_pos1,"1st star's position"
 print*,x_pos2,y_pos2,z_pos2,"2nd star's position" 
 endif
 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
    ierr = nerr
 endif

end subroutine read_setupfile




end module setup
