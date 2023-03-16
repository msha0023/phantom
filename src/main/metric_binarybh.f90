!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2023 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module metric
!
! :Runtime parameters:
!   - mass1 : *black hole mass in code units*
!   - mass2 : *2nd black hole mass in code units*
! :Dependencies: infile_utils, io
!
 implicit none
 character(len=*), parameter :: metric_type = 'binarybh'
 integer,          parameter :: imetric     = 4

 real, public  :: mass1 = 1.       ! mass of central object
 real, public  :: mass2 = 1.       ! mass of second black hole
 real, public  :: b     = 20.      ! separation between two BHs
 real, public  :: a     = 0.0      ! spin of central object
contains

!----------------------------------------------------------------
!+
!  The metric tensor in 'CARTESIAN-like form'
!+
!----------------------------------------------------------------
pure subroutine get_metric_cartesian(position,gcov,gcon,sqrtg,tt)
 ! subroutine get_metric_cartesian(position,gcov,gcon,sqrtg,tt)
  use iso_c_binding, only: c_double
  use metric_binarybh_utils
  use inverse4x4,    only: inv4x4

  real, intent(in) :: position(3) !This is the position of the particle
  real, intent(out) :: gcov(0:3,0:3)
  real, intent(in), optional:: tt
  real, intent(out), optional :: gcon(0:3,0:3)
  real, intent(out), optional :: sqrtg
  real :: xx,yy,zz,val
  real :: det
  !print*,tt,"tt in binarybh"
  if (present(sqrtg)) sqrtg = 1.

  gcov = 0.
  xx = position(1)
  yy = position(2)
  zz = position(3)
  
  gcov(0,0) = c_gtt(tt,xx,yy,zz,mass1,mass2,b)
  gcov(0,1) = c_gtx(tt,xx,yy,zz,mass1,mass2,b)
  gcov(0,2) = c_gty(tt,xx,yy,zz,mass1,mass2,b)
  gcov(0,3) = c_gtz(tt,xx,yy,zz,mass1,mass2,b)
  gcov(1,1) = c_gxx(tt,xx,yy,zz,mass1,mass2,b)
  gcov(1,2) = c_gxy(tt,xx,yy,zz,mass1,mass2,b)
  gcov(1,3) = c_gxz(tt,xx,yy,zz,mass1,mass2,b)
  gcov(2,2) = c_gyy(tt,xx,yy,zz,mass1,mass2,b)
  gcov(2,3) = c_gyz(tt,xx,yy,zz,mass1,mass2,b)
  gcov(3,3) = c_gzz(tt,xx,yy,zz,mass1,mass2,b)
  gcov(1,0) = gcov(0,1)
  gcov(2,0) = gcov(0,2)
  gcov(2,1) = gcov(1,2)
  gcov(3,0) = gcov(0,3)
  gcov(3,1) = gcov(1,3)
  gcov(3,2) = gcov(2,3)

  if (present(gcon)) then
      call inv4x4(gcov,gcon,det)
  endif 
end subroutine get_metric_cartesian

!----------------------------------------------------------------
!+
!  The metric tensor in SPHERICAL-like form
!+
!----------------------------------------------------------------
pure subroutine get_metric_spherical(position,gcov,gcon,sqrtg,tt)
! pure subroutine get_metric_spherical(position,gcov,gcon,sqrtg)
 real, intent(in)  :: position(3)
 real, intent(in), optional :: tt
 real, intent(out) :: gcov(0:3,0:3)
 real, intent(out), optional :: gcon(0:3,0:3)
 real, intent(out), optional :: sqrtg

 gcov = 0.
 
 if (present(gcon)) then
   gcon = 0.
 endif

 if (present(sqrtg)) sqrtg = 0.

end subroutine get_metric_spherical

!----------------------------------------------------------------
!+
!  Derivatives of the covariant 'CARTESIAN' metric
!+
!----------------------------------------------------------------
pure subroutine metric_cartesian_derivatives(position,dgcovdx, dgcovdy, dgcovdz)
 real, intent(in) :: position(3)
 real, intent(out), dimension(0:3,0:3) :: dgcovdx,dgcovdy,dgcovdz

 dgcovdx = 0.
 dgcovdy = 0.
 dgcovdz = 0.

end subroutine metric_cartesian_derivatives

!----------------------------------------------------------------
!+
!  Derivatives of the covariant 'SPHERICAL' metric
!+
!----------------------------------------------------------------
pure subroutine metric_spherical_derivatives(position,dgcovdr, dgcovdtheta, dgcovdphi)
 real, intent(in) :: position(3)
 real, intent(out), dimension(0:3,0:3) :: dgcovdr,dgcovdtheta,dgcovdphi

 dgcovdr = 0.
 dgcovdtheta = 0.
 dgcovdphi = 0.

end subroutine metric_spherical_derivatives


!----------------------------------------------------------------
!+
!  (Jacobian tensor) Derivatives of Boyer-Lindquist 'Spherical'
!  with respect to 'Cartesian' coordinates
!+
!----------------------------------------------------------------
pure subroutine get_jacobian(position,dxdx)
 real, intent(in), dimension(3) :: position
 real, intent(out), dimension(0:3,0:3) :: dxdx
 
 dxdx = 0.

end subroutine get_jacobian

!-----------------------------------------------------------------------
!+
!  Boyer-Lindquist coordinate transformations from CARTESIAN to SPHERICAL
!+
!-----------------------------------------------------------------------
pure subroutine cartesian2spherical(xcart,xspher)
 real, intent(in) :: xcart(3)
 real, intent(out) ::xspher(3)

 xspher = 0.

end subroutine cartesian2spherical

!-----------------------------------------------------------------------
!+
!  Boyer-Lindquist coordinate transformations from SPHERICAL to CARTESIAN
!+
!-----------------------------------------------------------------------
pure subroutine spherical2cartesian(xspher,xcart)
 real, intent(in) :: xspher(3)
 real, intent(out) :: xcart(3)

 xcart = 0.

end subroutine spherical2cartesian

!-----------------------------------------------------------------------
!+
!  writes metric options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_metric(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# options relating to the '//trim(metric_type)//' metric'

 call write_inopt(mass1,'mass1','black hole mass in code units',iunit)
 call write_inopt(mass2,'mass2','second black hole mass in code units',iunit)
 call write_inopt(b,'b','separation between black holes',iunit)
end subroutine write_options_metric

!-----------------------------------------------------------------------
!+
!  reads metric options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_metric(name,valstring,imatch,igotall,ierr)
 use io, only:fatal,warning
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 character(len=*), parameter :: tag = 'metric'
 integer, save :: ngot = 0

 imatch  = .true.
 igotall = .false.
 select case(trim(name))
 case('mass1')
    read(valstring,*,iostat=ierr) mass1
    if (mass1 < 0.)  call fatal(tag,'black hole mass: mass1 < 0')
    if (mass1 == 0.) call warning(tag,'black hole mass: mass1 = 0')
    ngot = ngot + 1
 case('mass2')
    read(valstring,*,iostat=ierr) mass2
    if (mass2 < 0.)  call fatal(tag,'black hole mass: mass2 < 0')
    if (mass2 == 0.) call warning(tag,'black hole mass: mass2 = 0')
    ngot = ngot + 1
case('b')
    read(valstring,*,iostat=ierr) b
    if (b < 0.)  call fatal(tag,'black hole separation: b < 0')
    if (b == 0.) call warning(tag,'black holes separation: b = 0')
    ngot = ngot + 1
 case default
    imatch = .false.
 end select

 igotall = (ngot >= 2)

end subroutine read_options_metric

end module metric
