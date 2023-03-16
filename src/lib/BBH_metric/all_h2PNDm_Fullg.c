#include <math.h>

double ftrans(double,double,double,double,double);

double all_COMTraj1x(double, double, double, double);
double all_COMTraj1y(double, double, double, double);
double all_COMTraj1z(double, double, double, double);
double all_COMTraj2x(double, double, double, double);
double all_COMTraj2y(double, double, double, double);
double all_COMTraj2z(double, double, double, double);

double all_h2PNDm_FZgtt(double,double,double,double,double,double,double);
double all_h2PNDm_FZgtx(double,double,double,double,double,double,double);
double all_h2PNDm_FZgty(double,double,double,double,double,double,double);
double all_h2PNDm_FZgtz(double,double,double,double,double,double,double);
double all_h2PNDm_FZgxx(double,double,double,double,double,double,double);
double all_h2PNDm_FZgxy(double,double,double,double,double,double,double);
double all_h2PNDm_FZgxz(double,double,double,double,double,double,double);
double all_h2PNDm_FZgyy(double,double,double,double,double,double,double);
double all_h2PNDm_FZgyz(double,double,double,double,double,double,double);
double all_h2PNDm_FZgzz(double,double,double,double,double,double,double);

double all_h2PNDm_NZgtt(double,double,double,double,double,double,double);
double all_h2PNDm_NZgtx(double,double,double,double,double,double,double);
double all_h2PNDm_NZgty(double,double,double,double,double,double,double);
double all_h2PNDm_NZgtz(double,double,double,double,double,double,double);
double all_h2PNDm_NZgxx(double,double,double,double,double,double,double);
double all_h2PNDm_NZgxy(double,double,double,double,double,double,double);
double all_h2PNDm_NZgxz(double,double,double,double,double,double,double);
double all_h2PNDm_NZgyy(double,double,double,double,double,double,double);
double all_h2PNDm_NZgyz(double,double,double,double,double,double,double);
double all_h2PNDm_NZgzz(double,double,double,double,double,double,double);

double all_h2PNDm_IZgtt(double,double,double,double,double,double,double);
double all_h2PNDm_IZgtx(double,double,double,double,double,double,double);
double all_h2PNDm_IZgty(double,double,double,double,double,double,double);
double all_h2PNDm_IZgtz(double,double,double,double,double,double,double);
double all_h2PNDm_IZgxx(double,double,double,double,double,double,double);
double all_h2PNDm_IZgxy(double,double,double,double,double,double,double);
double all_h2PNDm_IZgxz(double,double,double,double,double,double,double);
double all_h2PNDm_IZgyy(double,double,double,double,double,double,double);
double all_h2PNDm_IZgyz(double,double,double,double,double,double,double);
double all_h2PNDm_IZgzz(double,double,double,double,double,double,double);

double all_h2PNDm_IZ2gtt(double,double,double,double,double,double,double);
double all_h2PNDm_IZ2gtx(double,double,double,double,double,double,double);
double all_h2PNDm_IZ2gty(double,double,double,double,double,double,double);
double all_h2PNDm_IZ2gtz(double,double,double,double,double,double,double);
double all_h2PNDm_IZ2gxx(double,double,double,double,double,double,double);
double all_h2PNDm_IZ2gxy(double,double,double,double,double,double,double);
double all_h2PNDm_IZ2gxz(double,double,double,double,double,double,double);
double all_h2PNDm_IZ2gyy(double,double,double,double,double,double,double);
double all_h2PNDm_IZ2gyz(double,double,double,double,double,double,double);
double all_h2PNDm_IZ2gzz(double,double,double,double,double,double,double);

double all_h2PNDm_Fullgtt(
						 double tt,
						 double xx,
						 double yy,
						 double zz,
						 double m1,
						 double m2,
						 double b)
{
	double ftransfar, ftransnear, ftransinner1, ftransinner2;
	double lambda, r1T, r2T, PI, r, r1, r2;
	double finalmetric;
	
	PI = 3.1415926535897932384626433832795;
	lambda = PI*sqrt(b*b*b/(m1+m2));
	r1T = pow(m1*m1*m1*b*b*b*b*b/(m1+m2),1.0/7.0);
	r2T = pow(m2*m2*m2*b*b*b*b*b/(m1+m2),1.0/7.0);
	
	r = sqrt(xx*xx + yy*yy + zz*zz);
	
	r1 = sqrt( (xx - all_COMTraj1x(0,m1,m2,b))*(xx - all_COMTraj1x(0,m1,m2,b)) 
			  +(yy - all_COMTraj1y(0,m1,m2,b))*(yy - all_COMTraj1y(0,m1,m2,b)) 
			  +(zz - all_COMTraj1z(0,m1,m2,b))*(zz - all_COMTraj1z(0,m1,m2,b)));
	r2 = sqrt( (xx - all_COMTraj2x(0,m1,m2,b))*(xx - all_COMTraj2x(0,m1,m2,b)) 
			  +(yy - all_COMTraj2y(0,m1,m2,b))*(yy - all_COMTraj2y(0,m1,m2,b)) 
			  +(zz - all_COMTraj2z(0,m1,m2,b))*(zz - all_COMTraj2z(0,m1,m2,b)));	
	
	ftransfar = ftrans(r,lambda/5.0,lambda,1.0,2.50);
	ftransnear = ftrans(xx,2.2*m2 - m1*b/(m1+m2),b-2.2*(m1+m2),1.0,2.50);
	ftransinner1 = ftrans(r1,0.256*r1T,3.17*pow((m1+m2)*(m1+m2)*b*b*b*b*b,1.0/7.0),0.2,b/(m1+m2));
	ftransinner2 = ftrans(r2,0.256*r2T,3.17*pow((m1+m2)*(m1+m2)*b*b*b*b*b,1.0/7.0),0.2,b/(m1+m2));
		
	finalmetric = (1 - ftransfar)*(ftransnear*(ftransinner1*all_h2PNDm_NZgtt(tt,xx,yy,zz, m1,m2,b) 
											   + (1 - ftransinner1)*all_h2PNDm_IZgtt(tt,xx,yy,zz, m1,m2,b))
								   +(1 - ftransnear)*(ftransinner2*all_h2PNDm_NZgtt(tt,xx,yy,zz, m1,m2,b) 
													  + (1 - ftransinner2)*all_h2PNDm_IZ2gtt(tt,xx,yy,zz,m1,m2,b)))
				+ ftransfar*all_h2PNDm_FZgtt(tt,xx,yy,zz, m1,m2,b);
	
	return finalmetric;
}

double all_h2PNDm_Fullgtx(
						 double tt,
						 double xx,
						 double yy,
						 double zz,
						 double m1,
						 double m2,
						 double b)
{
	double ftransfar, ftransnear, ftransinner1, ftransinner2;
	double lambda, r1T, r2T, PI, r, r1, r2;
	double finalmetric;
	
	PI = 3.1415926535897932384626433832795;
	lambda = PI*sqrt(b*b*b/(m1+m2));
	r1T = pow(m1*m1*m1*b*b*b*b*b/(m1+m2),1.0/7.0);
	r2T = pow(m2*m2*m2*b*b*b*b*b/(m1+m2),1.0/7.0);
	
	r = sqrt(xx*xx + yy*yy + zz*zz);
	
	r1 = sqrt( (xx - all_COMTraj1x(0,m1,m2,b))*(xx - all_COMTraj1x(0,m1,m2,b)) 
			  +(yy - all_COMTraj1y(0,m1,m2,b))*(yy - all_COMTraj1y(0,m1,m2,b)) 
			  +(zz - all_COMTraj1z(0,m1,m2,b))*(zz - all_COMTraj1z(0,m1,m2,b)));
	r2 = sqrt( (xx - all_COMTraj2x(0,m1,m2,b))*(xx - all_COMTraj2x(0,m1,m2,b)) 
			  +(yy - all_COMTraj2y(0,m1,m2,b))*(yy - all_COMTraj2y(0,m1,m2,b)) 
			  +(zz - all_COMTraj2z(0,m1,m2,b))*(zz - all_COMTraj2z(0,m1,m2,b)));	
	
	ftransfar = ftrans(r,lambda/5.0,lambda,1.0,2.50);
	ftransnear = ftrans(xx,2.2*m2 - m1*b/(m1+m2),b-2.2*(m1+m2),1.0,2.50);
	ftransinner1 = ftrans(r1,0.256*r1T,3.17*pow((m1+m2)*(m1+m2)*b*b*b*b*b,1.0/7.0),0.2,b/(m1+m2));
	ftransinner2 = ftrans(r2,0.256*r2T,3.17*pow((m1+m2)*(m1+m2)*b*b*b*b*b,1.0/7.0),0.2,b/(m1+m2));
	
	finalmetric = (1 - ftransfar)*(ftransnear*(ftransinner1*all_h2PNDm_NZgtx(tt,xx,yy,zz, m1,m2,b) 
											   + (1 - ftransinner1)*all_h2PNDm_IZgtx(tt,xx,yy,zz, m1,m2,b))
								   +(1 - ftransnear)*(ftransinner2*all_h2PNDm_NZgtx(tt,xx,yy,zz, m1,m2,b) 
													  + (1 - ftransinner2)*all_h2PNDm_IZ2gtx(tt,xx,yy,zz,m1,m2,b)))
	+ ftransfar*all_h2PNDm_FZgtx(tt,xx,yy,zz, m1,m2,b);
	
	return finalmetric;
}

double all_h2PNDm_Fullgty(
						 double tt,
						 double xx,
						 double yy,
						 double zz,
						 double m1,
						 double m2,
						 double b)
{
	double ftransfar, ftransnear, ftransinner1, ftransinner2;
	double lambda, r1T, r2T, PI, r, r1, r2;
	double finalmetric;
	
	PI = 3.1415926535897932384626433832795;
	lambda = PI*sqrt(b*b*b/(m1+m2));
	r1T = pow(m1*m1*m1*b*b*b*b*b/(m1+m2),1.0/7.0);
	r2T = pow(m2*m2*m2*b*b*b*b*b/(m1+m2),1.0/7.0);
	
	r = sqrt(xx*xx + yy*yy + zz*zz);
	
	r1 = sqrt( (xx - all_COMTraj1x(0,m1,m2,b))*(xx - all_COMTraj1x(0,m1,m2,b)) 
			  +(yy - all_COMTraj1y(0,m1,m2,b))*(yy - all_COMTraj1y(0,m1,m2,b)) 
			  +(zz - all_COMTraj1z(0,m1,m2,b))*(zz - all_COMTraj1z(0,m1,m2,b)));
	r2 = sqrt( (xx - all_COMTraj2x(0,m1,m2,b))*(xx - all_COMTraj2x(0,m1,m2,b)) 
			  +(yy - all_COMTraj2y(0,m1,m2,b))*(yy - all_COMTraj2y(0,m1,m2,b)) 
			  +(zz - all_COMTraj2z(0,m1,m2,b))*(zz - all_COMTraj2z(0,m1,m2,b)));	
	
	ftransfar = ftrans(r,lambda/5.0,lambda,1.0,2.50);
	ftransnear = ftrans(xx,2.2*m2 - m1*b/(m1+m2),b-2.2*(m1+m2),1.0,2.50);
	ftransinner1 = ftrans(r1,0.256*r1T,3.17*pow((m1+m2)*(m1+m2)*b*b*b*b*b,1.0/7.0),0.2,b/(m1+m2));
	ftransinner2 = ftrans(r2,0.256*r2T,3.17*pow((m1+m2)*(m1+m2)*b*b*b*b*b,1.0/7.0),0.2,b/(m1+m2));
	
	finalmetric = (1 - ftransfar)*(ftransnear*(ftransinner1*all_h2PNDm_NZgty(tt,xx,yy,zz, m1,m2,b) 
											   + (1 - ftransinner1)*all_h2PNDm_IZgty(tt,xx,yy,zz, m1,m2,b))
								   +(1 - ftransnear)*(ftransinner2*all_h2PNDm_NZgty(tt,xx,yy,zz, m1,m2,b) 
													  + (1 - ftransinner2)*all_h2PNDm_IZ2gty(tt,xx,yy,zz,m1,m2,b)))
	+ ftransfar*all_h2PNDm_FZgty(tt,xx,yy,zz, m1,m2,b);
	
	return finalmetric;
}

double all_h2PNDm_Fullgtz(
						 double tt,
						 double xx,
						 double yy,
						 double zz,
						 double m1,
						 double m2,
						 double b)
{
	double ftransfar, ftransnear, ftransinner1, ftransinner2;
	double lambda, r1T, r2T, PI, r, r1, r2;
	double finalmetric;
	
	PI = 3.1415926535897932384626433832795;
	lambda = PI*sqrt(b*b*b/(m1+m2));
	r1T = pow(m1*m1*m1*b*b*b*b*b/(m1+m2),1.0/7.0);
	r2T = pow(m2*m2*m2*b*b*b*b*b/(m1+m2),1.0/7.0);
	
	r = sqrt(xx*xx + yy*yy + zz*zz);
	
	r1 = sqrt( (xx - all_COMTraj1x(0,m1,m2,b))*(xx - all_COMTraj1x(0,m1,m2,b)) 
			  +(yy - all_COMTraj1y(0,m1,m2,b))*(yy - all_COMTraj1y(0,m1,m2,b)) 
			  +(zz - all_COMTraj1z(0,m1,m2,b))*(zz - all_COMTraj1z(0,m1,m2,b)));
	r2 = sqrt( (xx - all_COMTraj2x(0,m1,m2,b))*(xx - all_COMTraj2x(0,m1,m2,b)) 
			  +(yy - all_COMTraj2y(0,m1,m2,b))*(yy - all_COMTraj2y(0,m1,m2,b)) 
			  +(zz - all_COMTraj2z(0,m1,m2,b))*(zz - all_COMTraj2z(0,m1,m2,b)));	
	
	ftransfar = ftrans(r,lambda/5.0,lambda,1.0,2.50);
	ftransnear = ftrans(xx,2.2*m2 - m1*b/(m1+m2),b-2.2*(m1+m2),1.0,2.50);
	ftransinner1 = ftrans(r1,0.256*r1T,3.17*pow((m1+m2)*(m1+m2)*b*b*b*b*b,1.0/7.0),0.2,b/(m1+m2));
	ftransinner2 = ftrans(r2,0.256*r2T,3.17*pow((m1+m2)*(m1+m2)*b*b*b*b*b,1.0/7.0),0.2,b/(m1+m2));
	
	finalmetric = (1 - ftransfar)*(ftransnear*(ftransinner1*all_h2PNDm_NZgtz(tt,xx,yy,zz, m1,m2,b) 
											   + (1 - ftransinner1)*all_h2PNDm_IZgtz(tt,xx,yy,zz, m1,m2,b))
								   +(1 - ftransnear)*(ftransinner2*all_h2PNDm_NZgtz(tt,xx,yy,zz, m1,m2,b) 
													  + (1 - ftransinner2)*all_h2PNDm_IZ2gtz(tt,xx,yy,zz,m1,m2,b)))
	+ ftransfar*all_h2PNDm_FZgtz(tt,xx,yy,zz, m1,m2,b);
	
	return finalmetric;
}

double all_h2PNDm_Fullgxx(
						 double tt,
						 double xx,
						 double yy,
						 double zz,
						 double m1,
						 double m2,
						 double b)
{
	double ftransfar, ftransnear, ftransinner1, ftransinner2;
	double lambda, r1T, r2T, PI, r, r1, r2;
	double finalmetric;
	
	PI = 3.1415926535897932384626433832795;
	lambda = PI*sqrt(b*b*b/(m1+m2));
	r1T = pow(m1*m1*m1*b*b*b*b*b/(m1+m2),1.0/7.0);
	r2T = pow(m2*m2*m2*b*b*b*b*b/(m1+m2),1.0/7.0);
	
	r = sqrt(xx*xx + yy*yy + zz*zz);
	
	r1 = sqrt( (xx - all_COMTraj1x(0,m1,m2,b))*(xx - all_COMTraj1x(0,m1,m2,b)) 
			  +(yy - all_COMTraj1y(0,m1,m2,b))*(yy - all_COMTraj1y(0,m1,m2,b)) 
			  +(zz - all_COMTraj1z(0,m1,m2,b))*(zz - all_COMTraj1z(0,m1,m2,b)));
	r2 = sqrt( (xx - all_COMTraj2x(0,m1,m2,b))*(xx - all_COMTraj2x(0,m1,m2,b)) 
			  +(yy - all_COMTraj2y(0,m1,m2,b))*(yy - all_COMTraj2y(0,m1,m2,b)) 
			  +(zz - all_COMTraj2z(0,m1,m2,b))*(zz - all_COMTraj2z(0,m1,m2,b)));	
	
	ftransfar = ftrans(r,lambda/5.0,lambda,1.0,2.50);
	ftransnear = ftrans(xx,2.2*m2 - m1*b/(m1+m2),b-2.2*(m1+m2),1.0,2.50);
	ftransinner1 = ftrans(r1,0.256*r1T,3.17*pow((m1+m2)*(m1+m2)*b*b*b*b*b,1.0/7.0),0.2,b/(m1+m2));
	ftransinner2 = ftrans(r2,0.256*r2T,3.17*pow((m1+m2)*(m1+m2)*b*b*b*b*b,1.0/7.0),0.2,b/(m1+m2));
	
	finalmetric = (1 - ftransfar)*(ftransnear*(ftransinner1*all_h2PNDm_NZgxx(tt,xx,yy,zz, m1,m2,b) 
											   + (1 - ftransinner1)*all_h2PNDm_IZgxx(tt,xx,yy,zz, m1,m2,b))
								   +(1 - ftransnear)*(ftransinner2*all_h2PNDm_NZgxx(tt,xx,yy,zz, m1,m2,b) 
													  + (1 - ftransinner2)*all_h2PNDm_IZ2gxx(tt,xx,yy,zz,m1,m2,b)))
	+ ftransfar*all_h2PNDm_FZgxx(tt,xx,yy,zz, m1,m2,b);
	
	return finalmetric;
}

double all_h2PNDm_Fullgxy(
						 double tt,
						 double xx,
						 double yy,
						 double zz,
						 double m1,
						 double m2,
						 double b)
{
	double ftransfar, ftransnear, ftransinner1, ftransinner2;
	double lambda, r1T, r2T, PI, r, r1, r2;
	double finalmetric;
	
	PI = 3.1415926535897932384626433832795;
	lambda = PI*sqrt(b*b*b/(m1+m2));
	r1T = pow(m1*m1*m1*b*b*b*b*b/(m1+m2),1.0/7.0);
	r2T = pow(m2*m2*m2*b*b*b*b*b/(m1+m2),1.0/7.0);
	
	r = sqrt(xx*xx + yy*yy + zz*zz);
	
	r1 = sqrt( (xx - all_COMTraj1x(0,m1,m2,b))*(xx - all_COMTraj1x(0,m1,m2,b)) 
			  +(yy - all_COMTraj1y(0,m1,m2,b))*(yy - all_COMTraj1y(0,m1,m2,b)) 
			  +(zz - all_COMTraj1z(0,m1,m2,b))*(zz - all_COMTraj1z(0,m1,m2,b)));
	r2 = sqrt( (xx - all_COMTraj2x(0,m1,m2,b))*(xx - all_COMTraj2x(0,m1,m2,b)) 
			  +(yy - all_COMTraj2y(0,m1,m2,b))*(yy - all_COMTraj2y(0,m1,m2,b)) 
			  +(zz - all_COMTraj2z(0,m1,m2,b))*(zz - all_COMTraj2z(0,m1,m2,b)));	
	
	ftransfar = ftrans(r,lambda/5.0,lambda,1.0,2.50);
	ftransnear = ftrans(xx,2.2*m2 - m1*b/(m1+m2),b-2.2*(m1+m2),1.0,2.50);
	ftransinner1 = ftrans(r1,0.256*r1T,3.17*pow((m1+m2)*(m1+m2)*b*b*b*b*b,1.0/7.0),0.2,b/(m1+m2));
	ftransinner2 = ftrans(r2,0.256*r2T,3.17*pow((m1+m2)*(m1+m2)*b*b*b*b*b,1.0/7.0),0.2,b/(m1+m2));
	
	finalmetric = (1 - ftransfar)*(ftransnear*(ftransinner1*all_h2PNDm_NZgxy(tt,xx,yy,zz, m1,m2,b) 
											   + (1 - ftransinner1)*all_h2PNDm_IZgxy(tt,xx,yy,zz, m1,m2,b))
								   +(1 - ftransnear)*(ftransinner2*all_h2PNDm_NZgxy(tt,xx,yy,zz, m1,m2,b) 
													  + (1 - ftransinner2)*all_h2PNDm_IZ2gxy(tt,xx,yy,zz,m1,m2,b)))
	+ ftransfar*all_h2PNDm_FZgxy(tt,xx,yy,zz, m1,m2,b);
	
	return finalmetric;
}

double all_h2PNDm_Fullgxz(
						 double tt,
						 double xx,
						 double yy,
						 double zz,
						 double m1,
						 double m2,
						 double b)
{
	double ftransfar, ftransnear, ftransinner1, ftransinner2;
	double lambda, r1T, r2T, PI, r, r1, r2;
	double finalmetric;
	
	PI = 3.1415926535897932384626433832795;
	lambda = PI*sqrt(b*b*b/(m1+m2));
	r1T = pow(m1*m1*m1*b*b*b*b*b/(m1+m2),1.0/7.0);
	r2T = pow(m2*m2*m2*b*b*b*b*b/(m1+m2),1.0/7.0);
	
	r = sqrt(xx*xx + yy*yy + zz*zz);
	
	r1 = sqrt( (xx - all_COMTraj1x(0,m1,m2,b))*(xx - all_COMTraj1x(0,m1,m2,b)) 
			  +(yy - all_COMTraj1y(0,m1,m2,b))*(yy - all_COMTraj1y(0,m1,m2,b)) 
			  +(zz - all_COMTraj1z(0,m1,m2,b))*(zz - all_COMTraj1z(0,m1,m2,b)));
	r2 = sqrt( (xx - all_COMTraj2x(0,m1,m2,b))*(xx - all_COMTraj2x(0,m1,m2,b)) 
			  +(yy - all_COMTraj2y(0,m1,m2,b))*(yy - all_COMTraj2y(0,m1,m2,b)) 
			  +(zz - all_COMTraj2z(0,m1,m2,b))*(zz - all_COMTraj2z(0,m1,m2,b)));	
	
	ftransfar = ftrans(r,lambda/5.0,lambda,1.0,2.50);
	ftransnear = ftrans(xx,2.2*m2 - m1*b/(m1+m2),b-2.2*(m1+m2),1.0,2.50);
	ftransinner1 = ftrans(r1,0.256*r1T,3.17*pow((m1+m2)*(m1+m2)*b*b*b*b*b,1.0/7.0),0.2,b/(m1+m2));
	ftransinner2 = ftrans(r2,0.256*r2T,3.17*pow((m1+m2)*(m1+m2)*b*b*b*b*b,1.0/7.0),0.2,b/(m1+m2));
	
	finalmetric = (1 - ftransfar)*(ftransnear*(ftransinner1*all_h2PNDm_NZgxz(tt,xx,yy,zz, m1,m2,b) 
											   + (1 - ftransinner1)*all_h2PNDm_IZgxz(tt,xx,yy,zz, m1,m2,b))
								   +(1 - ftransnear)*(ftransinner2*all_h2PNDm_NZgxz(tt,xx,yy,zz, m1,m2,b) 
													  + (1 - ftransinner2)*all_h2PNDm_IZ2gxz(tt,xx,yy,zz,m1,m2,b)))
	+ ftransfar*all_h2PNDm_FZgxz(tt,xx,yy,zz, m1,m2,b);
	
	return finalmetric;
}

double all_h2PNDm_Fullgyy(
						 double tt,
						 double xx,
						 double yy,
						 double zz,
						 double m1,
						 double m2,
						 double b)
{
	double ftransfar, ftransnear, ftransinner1, ftransinner2;
	double lambda, r1T, r2T, PI, r, r1, r2;
	double finalmetric;
	
	PI = 3.1415926535897932384626433832795;
	lambda = PI*sqrt(b*b*b/(m1+m2));
	r1T = pow(m1*m1*m1*b*b*b*b*b/(m1+m2),1.0/7.0);
	r2T = pow(m2*m2*m2*b*b*b*b*b/(m1+m2),1.0/7.0);
	
	r = sqrt(xx*xx + yy*yy + zz*zz);
	
	r1 = sqrt( (xx - all_COMTraj1x(0,m1,m2,b))*(xx - all_COMTraj1x(0,m1,m2,b)) 
			  +(yy - all_COMTraj1y(0,m1,m2,b))*(yy - all_COMTraj1y(0,m1,m2,b)) 
			  +(zz - all_COMTraj1z(0,m1,m2,b))*(zz - all_COMTraj1z(0,m1,m2,b)));
	r2 = sqrt( (xx - all_COMTraj2x(0,m1,m2,b))*(xx - all_COMTraj2x(0,m1,m2,b)) 
			  +(yy - all_COMTraj2y(0,m1,m2,b))*(yy - all_COMTraj2y(0,m1,m2,b)) 
			  +(zz - all_COMTraj2z(0,m1,m2,b))*(zz - all_COMTraj2z(0,m1,m2,b)));	
	
	ftransfar = ftrans(r,lambda/5.0,lambda,1.0,2.50);
	ftransnear = ftrans(xx,2.2*m2 - m1*b/(m1+m2),b-2.2*(m1+m2),1.0,2.50);
	ftransinner1 = ftrans(r1,0.256*r1T,3.17*pow((m1+m2)*(m1+m2)*b*b*b*b*b,1.0/7.0),0.2,b/(m1+m2));
	ftransinner2 = ftrans(r2,0.256*r2T,3.17*pow((m1+m2)*(m1+m2)*b*b*b*b*b,1.0/7.0),0.2,b/(m1+m2));
	
	finalmetric = (1 - ftransfar)*(ftransnear*(ftransinner1*all_h2PNDm_NZgyy(tt,xx,yy,zz, m1,m2,b) 
											   + (1 - ftransinner1)*all_h2PNDm_IZgyy(tt,xx,yy,zz, m1,m2,b))
								   +(1 - ftransnear)*(ftransinner2*all_h2PNDm_NZgyy(tt,xx,yy,zz, m1,m2,b) 
													  + (1 - ftransinner2)*all_h2PNDm_IZ2gyy(tt,xx,yy,zz,m1,m2,b)))
	+ ftransfar*all_h2PNDm_FZgyy(tt,xx,yy,zz, m1,m2,b);
	
	return finalmetric;
}

double all_h2PNDm_Fullgyz(
						 double tt,
						 double xx,
						 double yy,
						 double zz,
						 double m1,
						 double m2,
						 double b)
{
	double ftransfar, ftransnear, ftransinner1, ftransinner2;
	double lambda, r1T, r2T, PI, r, r1, r2;
	double finalmetric;
	
	PI = 3.1415926535897932384626433832795;
	lambda = PI*sqrt(b*b*b/(m1+m2));
	r1T = pow(m1*m1*m1*b*b*b*b*b/(m1+m2),1.0/7.0);
	r2T = pow(m2*m2*m2*b*b*b*b*b/(m1+m2),1.0/7.0);
	
	r = sqrt(xx*xx + yy*yy + zz*zz);
	
	r1 = sqrt( (xx - all_COMTraj1x(0,m1,m2,b))*(xx - all_COMTraj1x(0,m1,m2,b)) 
			  +(yy - all_COMTraj1y(0,m1,m2,b))*(yy - all_COMTraj1y(0,m1,m2,b)) 
			  +(zz - all_COMTraj1z(0,m1,m2,b))*(zz - all_COMTraj1z(0,m1,m2,b)));
	r2 = sqrt( (xx - all_COMTraj2x(0,m1,m2,b))*(xx - all_COMTraj2x(0,m1,m2,b)) 
			  +(yy - all_COMTraj2y(0,m1,m2,b))*(yy - all_COMTraj2y(0,m1,m2,b)) 
			  +(zz - all_COMTraj2z(0,m1,m2,b))*(zz - all_COMTraj2z(0,m1,m2,b)));	
	
	ftransfar = ftrans(r,lambda/5.0,lambda,1.0,2.50);
	ftransnear = ftrans(xx,2.2*m2 - m1*b/(m1+m2),b-2.2*(m1+m2),1.0,2.50);
	ftransinner1 = ftrans(r1,0.256*r1T,3.17*pow((m1+m2)*(m1+m2)*b*b*b*b*b,1.0/7.0),0.2,b/(m1+m2));
	ftransinner2 = ftrans(r2,0.256*r2T,3.17*pow((m1+m2)*(m1+m2)*b*b*b*b*b,1.0/7.0),0.2,b/(m1+m2));
	
	finalmetric = (1 - ftransfar)*(ftransnear*(ftransinner1*all_h2PNDm_NZgyz(tt,xx,yy,zz, m1,m2,b) 
											   + (1 - ftransinner1)*all_h2PNDm_IZgyz(tt,xx,yy,zz, m1,m2,b))
								   +(1 - ftransnear)*(ftransinner2*all_h2PNDm_NZgyz(tt,xx,yy,zz, m1,m2,b) 
													  + (1 - ftransinner2)*all_h2PNDm_IZ2gyz(tt,xx,yy,zz,m1,m2,b)))
	+ ftransfar*all_h2PNDm_FZgyz(tt,xx,yy,zz, m1,m2,b);
	
	return finalmetric;
}

double all_h2PNDm_Fullgzz(
						 double tt,
						 double xx,
						 double yy,
						 double zz,
						 double m1,
						 double m2,
						 double b)
{
	double ftransfar, ftransnear, ftransinner1, ftransinner2;
	double lambda, r1T, r2T, PI, r, r1, r2;
	double finalmetric;
	
	PI = 3.1415926535897932384626433832795;
	lambda = PI*sqrt(b*b*b/(m1+m2));
	r1T = pow(m1*m1*m1*b*b*b*b*b/(m1+m2),1.0/7.0);
	r2T = pow(m2*m2*m2*b*b*b*b*b/(m1+m2),1.0/7.0);
	
	r = sqrt(xx*xx + yy*yy + zz*zz);
	
	r1 = sqrt( (xx - all_COMTraj1x(0,m1,m2,b))*(xx - all_COMTraj1x(0,m1,m2,b)) 
			  +(yy - all_COMTraj1y(0,m1,m2,b))*(yy - all_COMTraj1y(0,m1,m2,b)) 
			  +(zz - all_COMTraj1z(0,m1,m2,b))*(zz - all_COMTraj1z(0,m1,m2,b)));
	r2 = sqrt( (xx - all_COMTraj2x(0,m1,m2,b))*(xx - all_COMTraj2x(0,m1,m2,b)) 
			  +(yy - all_COMTraj2y(0,m1,m2,b))*(yy - all_COMTraj2y(0,m1,m2,b)) 
			  +(zz - all_COMTraj2z(0,m1,m2,b))*(zz - all_COMTraj2z(0,m1,m2,b)));	
	
	ftransfar = ftrans(r,lambda/5.0,lambda,1.0,2.50);
	ftransnear = ftrans(xx,2.2*m2 - m1*b/(m1+m2),b-2.2*(m1+m2),1.0,2.50);
	ftransinner1 = ftrans(r1,0.256*r1T,3.17*pow((m1+m2)*(m1+m2)*b*b*b*b*b,1.0/7.0),0.2,b/(m1+m2));
	ftransinner2 = ftrans(r2,0.256*r2T,3.17*pow((m1+m2)*(m1+m2)*b*b*b*b*b,1.0/7.0),0.2,b/(m1+m2));
	
	finalmetric = (1 - ftransfar)*(ftransnear*(ftransinner1*all_h2PNDm_NZgzz(tt,xx,yy,zz, m1,m2,b) 
											   + (1 - ftransinner1)*all_h2PNDm_IZgzz(tt,xx,yy,zz, m1,m2,b))
								   +(1 - ftransnear)*(ftransinner2*all_h2PNDm_NZgzz(tt,xx,yy,zz, m1,m2,b) 
													  + (1 - ftransinner2)*all_h2PNDm_IZ2gzz(tt,xx,yy,zz,m1,m2,b)))
	+ ftransfar*all_h2PNDm_FZgzz(tt,xx,yy,zz, m1,m2,b);
	
	return finalmetric;
}

