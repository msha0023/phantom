double all_h2PNDm_IZgtt(double, double, double, double, double, double, double);
double all_h2PNDm_IZgtx(double, double, double, double, double, double, double);
double all_h2PNDm_IZgty(double, double, double, double, double, double, double);
double all_h2PNDm_IZgtz(double, double, double, double, double, double, double);
double all_h2PNDm_IZgxx(double, double, double, double, double, double, double);
double all_h2PNDm_IZgxy(double, double, double, double, double, double, double);
double all_h2PNDm_IZgxz(double, double, double, double, double, double, double);
double all_h2PNDm_IZgyy(double, double, double, double, double, double, double);
double all_h2PNDm_IZgyz(double, double, double, double, double, double, double);
double all_h2PNDm_IZgzz(double, double, double, double, double, double, double);

double all_h2PNDm_IZ2gtt (double t,double x,double y,double z,double m1,double m2,double b)
{ 
return( 
	   all_h2PNDm_IZgtt(t,-x,-y,z,m2,m1,b) 
	   ); 
}

double all_h2PNDm_IZ2gtx (double t,double x,double y,double z,double m1,double m2,double b)
{ 
	return( 
		   -all_h2PNDm_IZgtx(t,-x,-y,z,m2,m1,b) 
		   ); 
}

double all_h2PNDm_IZ2gty (double t,double x,double y,double z,double m1,double m2,double b)
{ 
	return( 
		   -all_h2PNDm_IZgty(t,-x,-y,z,m2,m1,b) 
		   ); 
}

double all_h2PNDm_IZ2gtz (double t,double x,double y,double z,double m1,double m2,double b)
{ 
	return( 
		   all_h2PNDm_IZgtz(t,-x,-y,z,m2,m1,b) 
		   ); 
}

double all_h2PNDm_IZ2gxx (double t,double x,double y,double z,double m1,double m2,double b)
{ 
	return( 
		   all_h2PNDm_IZgxx(t,-x,-y,z,m2,m1,b) 
		   ); 
}

double all_h2PNDm_IZ2gxy (double t,double x,double y,double z,double m1,double m2,double b)
{ 
	return( 
		   all_h2PNDm_IZgxy(t,-x,-y,z,m2,m1,b) 
		   ); 
}

double all_h2PNDm_IZ2gxz (double t,double x,double y,double z,double m1,double m2,double b)
{ 
	return( 
		   -all_h2PNDm_IZgxz(t,-x,-y,z,m2,m1,b) 
		   ); 
}

double all_h2PNDm_IZ2gyy (double t,double x,double y,double z,double m1,double m2,double b)
{ 
	return( 
		   all_h2PNDm_IZgyy(t,-x,-y,z,m2,m1,b) 
		   ); 
}

double all_h2PNDm_IZ2gyz (double t,double x,double y,double z,double m1,double m2,double b)
{ 
	return( 
		   -all_h2PNDm_IZgyz(t,-x,-y,z,m2,m1,b) 
		   ); 
}

double all_h2PNDm_IZ2gzz (double t,double x,double y,double z,double m1,double m2,double b)
{ 
	return( 
		   all_h2PNDm_IZgzz(t,-x,-y,z,m2,m1,b) 
		   ); 
}

