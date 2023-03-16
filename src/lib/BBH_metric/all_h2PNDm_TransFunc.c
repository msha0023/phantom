#include <math.h>

double ftrans(
			  double r, 
			  double r0,
			  double w,
			  double q,
			  double s)
{
	double func, chifunc, PI;
	
	PI = 3.1415926535897932384626433832795;
	
	chifunc = tan(PI*(r - r0)/(2.0*w));
	
	if(r <= r0) func = 0.0;
	if(r >= r0 + w) func = 1.0;
	if(r > r0 && r < r0 + w) func = (1.0/2.0)*(1.0 + tanh((s/PI)*(chifunc - q*q/chifunc)));
	
	return func;
}

