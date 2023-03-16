!This is an interface to read the functions from all_h2PNDm_Fullg.c present in
! ../src/lib/Binary_bh
module metric_binarybh_utils 

 use iso_c_binding, only:c_float,c_double 
 implicit none 

private 

public :: c_gtt,c_gtx,c_gty,c_gtz,c_gxx,c_gxy,c_gxz,c_gyy,c_gyz,c_gzz

interface c_gtt
 !we call the all_h2PNDm_Fullgtt function to get gtt value 

 pure function all_h2PNDm_Fullgtt(tt,xx,yy,zz,m1,m2,b)bind(c,NAME='all_h2PNDm_Fullgtt')

    import 
    implicit none 
    real(kind=c_double),intent(in),value :: tt,xx,yy,zz,m1,m2,b
    real(kind=c_double) :: all_h2PNDm_Fullgtt

    end function all_h2PNDm_Fullgtt

end interface c_gtt

interface c_gtx
pure function all_h2PNDm_Fullgtx(tt,xx,yy,zz,m1,m2,b)bind(c,NAME='all_h2PNDm_Fullgtx')

    import
    implicit none
    real(kind=c_double),intent(in),value :: tt,xx,yy,zz,m1,m2,b
    real(kind=c_double) :: all_h2PNDm_Fullgtx

    end function all_h2PNDm_Fullgtx

end interface c_gtx

interface c_gty
pure function all_h2PNDm_Fullgty(tt,xx,yy,zz,m1,m2,b)bind(c,NAME='all_h2PNDm_FZgty')

    import
    implicit none
    real(kind=c_double),intent(in),value :: tt,xx,yy,zz,m1,m2,b
    real(kind=c_double) :: all_h2PNDm_Fullgty

    end function all_h2PNDm_Fullgty

end interface c_gty

interface c_gtz
pure function all_h2PNDm_Fullgtz(tt,xx,yy,zz,m1,m2,b)bind(c,NAME='all_h2PNDm_Fullgtz')

    import
    implicit none
    real(kind=c_double),intent(in),value :: tt,xx,yy,zz,m1,m2,b
    real(kind=c_double) :: all_h2PNDm_Fullgtz

    end function all_h2PNDm_Fullgtz

end interface c_gtz

interface c_gxx
pure function all_h2PNDm_Fullgxx(tt,xx,yy,zz,m1,m2,b)bind(c,NAME='all_h2PNDm_Fullgxx')

    import
    implicit none
    real(kind=c_double),intent(in),value :: tt,xx,yy,zz,m1,m2,b
    real(kind=c_double) :: all_h2PNDm_Fullgxx

    end function all_h2PNDm_Fullgxx

end interface c_gxx

interface c_gxy
pure function all_h2PNDm_Fullgxy(tt,xx,yy,zz,m1,m2,b)bind(c,NAME='all_h2PNDm_Fullgxy')

    import
    implicit none
    real(kind=c_double),intent(in),value :: tt,xx,yy,zz,m1,m2,b
    real(kind=c_double) :: all_h2PNDm_Fullgxy

    end function all_h2PNDm_Fullgxy

end interface c_gxy

interface c_gxz
pure function all_h2PNDm_Fullgxz(tt,xx,yy,zz,m1,m2,b)bind(c,NAME='all_h2PNDm_Fullgxz')

    import
    implicit none
    real(kind=c_double),intent(in),value :: tt,xx,yy,zz,m1,m2,b
    real(kind=c_double) :: all_h2PNDm_Fullgxz

    end function all_h2PNDm_Fullgxz

end interface c_gxz

interface c_gyy
pure function all_h2PNDm_Fullgyy(tt,xx,yy,zz,m1,m2,b)bind(c,NAME='all_h2PNDm_Fullgyy')

    import
    implicit none
    real(kind=c_double),intent(in),value :: tt,xx,yy,zz,m1,m2,b
    real(kind=c_double) :: all_h2PNDm_Fullgyy

    end function all_h2PNDm_Fullgyy

end interface c_gyy

interface c_gyz
pure function all_h2PNDm_Fullgyz(tt,xx,yy,zz,m1,m2,b)bind(c,NAME='all_h2PNDm_Fullgyz')

    import
    implicit none
    real(kind=c_double),intent(in),value :: tt,xx,yy,zz,m1,m2,b
    real(kind=c_double) :: all_h2PNDm_Fullgyz

    end function all_h2PNDm_Fullgyz

end interface c_gyz

interface c_gzz
pure function all_h2PNDm_Fullgzz(tt,xx,yy,zz,m1,m2,b)bind(c,NAME='all_h2PNDm_Fullgzz')

    import
    implicit none
    real(kind=c_double),intent(in),value :: tt,xx,yy,zz,m1,m2,b
    real(kind=c_double) :: all_h2PNDm_Fullgzz

    end function all_h2PNDm_Fullgzz

end interface c_gzz
end module metric_binarybh_utils   


