double all_h2PNDm_prematchIZgtt_func(double, double, double, double,double,double,double);
double all_h2PNDm_prematchIZgtx_func(double, double, double, double,double,double,double);
double all_h2PNDm_prematchIZgty_func(double, double, double, double,double,double,double);
double all_h2PNDm_prematchIZgtz_func(double, double, double, double,double,double,double);
double all_h2PNDm_prematchIZgxx_func(double, double, double, double,double,double,double);
double all_h2PNDm_prematchIZgxy_func(double, double, double, double,double,double,double);
double all_h2PNDm_prematchIZgxz_func(double, double, double, double,double,double,double);
double all_h2PNDm_prematchIZgyy_func(double, double, double, double,double,double,double);
double all_h2PNDm_prematchIZgyz_func(double, double, double, double,double,double,double);
double all_h2PNDm_prematchIZgzz_func(double, double, double, double,double,double,double);

double all_h2PNDm_dTdt_func(double, double, double, double, double, double, double);
double all_h2PNDm_dTdx_func(double, double, double, double, double, double, double);
double all_h2PNDm_dTdy_func(double, double, double, double, double, double, double);
double all_h2PNDm_dTdz_func(double, double, double, double, double, double, double);
double all_h2PNDm_dXdt_func(double, double, double, double, double, double, double);
double all_h2PNDm_dXdx_func(double, double, double, double, double, double, double);
double all_h2PNDm_dXdy_func(double, double, double, double, double, double, double);
double all_h2PNDm_dXdz_func(double, double, double, double, double, double, double);
double all_h2PNDm_dYdt_func(double, double, double, double, double, double, double);
double all_h2PNDm_dYdx_func(double, double, double, double, double, double, double);
double all_h2PNDm_dYdy_func(double, double, double, double, double, double, double);
double all_h2PNDm_dYdz_func(double, double, double, double, double, double, double);
double all_h2PNDm_dZdt_func(double, double, double, double, double, double, double);
double all_h2PNDm_dZdx_func(double, double, double, double, double, double, double);
double all_h2PNDm_dZdy_func(double, double, double, double, double, double, double);
double all_h2PNDm_dZdz_func(double, double, double, double, double, double, double);

double all_h2PNDm_IZgtt (
  double t,
  double x,
  double y,
  double z,
  double m1,
  double m2,
  double b)
{

return(
   all_h2PNDm_prematchIZgtt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dTdt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dTdt_func(t,x,y,z,m1,m2,b) 
+  all_h2PNDm_prematchIZgxx_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dXdt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dXdt_func(t,x,y,z,m1,m2,b) 
+  all_h2PNDm_prematchIZgyy_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dYdt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dYdt_func(t,x,y,z,m1,m2,b) 
+  all_h2PNDm_prematchIZgzz_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dZdt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dZdt_func(t,x,y,z,m1,m2,b) 
+  all_h2PNDm_prematchIZgtx_func(t,x,y,z,m1,m2,b) * ( all_h2PNDm_dTdt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dXdt_func(t,x,y,z,m1,m2,b) + all_h2PNDm_dXdt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dTdt_func(t,x,y,z,m1,m2,b)) 
+  all_h2PNDm_prematchIZgty_func(t,x,y,z,m1,m2,b) * ( all_h2PNDm_dTdt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dYdt_func(t,x,y,z,m1,m2,b) + all_h2PNDm_dYdt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dTdt_func(t,x,y,z,m1,m2,b)) 
+  all_h2PNDm_prematchIZgtz_func(t,x,y,z,m1,m2,b) * ( all_h2PNDm_dTdt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dZdt_func(t,x,y,z,m1,m2,b) + all_h2PNDm_dZdt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dTdt_func(t,x,y,z,m1,m2,b)) 
+  all_h2PNDm_prematchIZgxy_func(t,x,y,z,m1,m2,b) * ( all_h2PNDm_dXdt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dYdt_func(t,x,y,z,m1,m2,b) + all_h2PNDm_dYdt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dXdt_func(t,x,y,z,m1,m2,b)) 
+  all_h2PNDm_prematchIZgxz_func(t,x,y,z,m1,m2,b) * ( all_h2PNDm_dXdt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dZdt_func(t,x,y,z,m1,m2,b) + all_h2PNDm_dZdt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dXdt_func(t,x,y,z,m1,m2,b)) 
+  all_h2PNDm_prematchIZgyz_func(t,x,y,z,m1,m2,b) * ( all_h2PNDm_dYdt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dZdt_func(t,x,y,z,m1,m2,b) + all_h2PNDm_dZdt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dYdt_func(t,x,y,z,m1,m2,b)) 
);

}

double all_h2PNDm_IZgtx (
  double t,
  double x,
  double y,
  double z,
  double m1,
  double m2,
  double b)
{

return(
   all_h2PNDm_prematchIZgtt_func(t,x,y,z,m1,m2,b)   * all_h2PNDm_dTdt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dTdx_func(t,x,y,z,m1,m2,b) 
+ all_h2PNDm_prematchIZgxx_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dXdt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dXdx_func(t,x,y,z,m1,m2,b) 
+ all_h2PNDm_prematchIZgyy_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dYdt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dYdx_func(t,x,y,z,m1,m2,b) 
+ all_h2PNDm_prematchIZgzz_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dZdt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dZdx_func(t,x,y,z,m1,m2,b) 
+  all_h2PNDm_prematchIZgtx_func(t,x,y,z,m1,m2,b) * ( all_h2PNDm_dTdt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dXdx_func(t,x,y,z,m1,m2,b) + all_h2PNDm_dXdt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dTdx_func(t,x,y,z,m1,m2,b)) 
+  all_h2PNDm_prematchIZgty_func(t,x,y,z,m1,m2,b) * ( all_h2PNDm_dTdt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dYdx_func(t,x,y,z,m1,m2,b) + all_h2PNDm_dYdt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dTdx_func(t,x,y,z,m1,m2,b)) 
+  all_h2PNDm_prematchIZgtz_func(t,x,y,z,m1,m2,b) * ( all_h2PNDm_dTdt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dZdx_func(t,x,y,z,m1,m2,b) + all_h2PNDm_dZdt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dTdx_func(t,x,y,z,m1,m2,b)) 
+  all_h2PNDm_prematchIZgxy_func(t,x,y,z,m1,m2,b) * ( all_h2PNDm_dXdt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dYdx_func(t,x,y,z,m1,m2,b) + all_h2PNDm_dYdt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dXdx_func(t,x,y,z,m1,m2,b)) 
+  all_h2PNDm_prematchIZgxz_func(t,x,y,z,m1,m2,b) * ( all_h2PNDm_dXdt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dZdx_func(t,x,y,z,m1,m2,b) + all_h2PNDm_dZdt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dXdx_func(t,x,y,z,m1,m2,b)) 
+  all_h2PNDm_prematchIZgyz_func(t,x,y,z,m1,m2,b) * ( all_h2PNDm_dYdt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dZdx_func(t,x,y,z,m1,m2,b) + all_h2PNDm_dZdt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dYdx_func(t,x,y,z,m1,m2,b)) 
);
}

double all_h2PNDm_IZgty (
  double t,
  double x,
  double y,
  double z,
  double m1,
  double m2,
  double b)
{
return(
   all_h2PNDm_prematchIZgtt_func(t,x,y,z,m1,m2,b)   * all_h2PNDm_dTdt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dTdy_func(t,x,y,z,m1,m2,b) 
+ all_h2PNDm_prematchIZgxx_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dXdt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dXdy_func(t,x,y,z,m1,m2,b) 
+ all_h2PNDm_prematchIZgyy_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dYdt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dYdy_func(t,x,y,z,m1,m2,b) 
+ all_h2PNDm_prematchIZgzz_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dZdt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dZdy_func(t,x,y,z,m1,m2,b) 
+  all_h2PNDm_prematchIZgtx_func(t,x,y,z,m1,m2,b) * ( all_h2PNDm_dTdt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dXdy_func(t,x,y,z,m1,m2,b) + all_h2PNDm_dXdt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dTdy_func(t,x,y,z,m1,m2,b)) 
+  all_h2PNDm_prematchIZgty_func(t,x,y,z,m1,m2,b) * ( all_h2PNDm_dTdt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dYdy_func(t,x,y,z,m1,m2,b) + all_h2PNDm_dYdt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dTdy_func(t,x,y,z,m1,m2,b)) 
+  all_h2PNDm_prematchIZgtz_func(t,x,y,z,m1,m2,b) * ( all_h2PNDm_dTdt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dZdy_func(t,x,y,z,m1,m2,b) + all_h2PNDm_dZdt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dTdy_func(t,x,y,z,m1,m2,b)) 
+  all_h2PNDm_prematchIZgxy_func(t,x,y,z,m1,m2,b) * ( all_h2PNDm_dXdt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dYdy_func(t,x,y,z,m1,m2,b) + all_h2PNDm_dYdt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dXdy_func(t,x,y,z,m1,m2,b)) 
+  all_h2PNDm_prematchIZgxz_func(t,x,y,z,m1,m2,b) * ( all_h2PNDm_dXdt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dZdy_func(t,x,y,z,m1,m2,b) + all_h2PNDm_dZdt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dXdy_func(t,x,y,z,m1,m2,b)) 
+  all_h2PNDm_prematchIZgyz_func(t,x,y,z,m1,m2,b) * ( all_h2PNDm_dYdt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dZdy_func(t,x,y,z,m1,m2,b) + all_h2PNDm_dZdt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dYdy_func(t,x,y,z,m1,m2,b)) 
);
}

double all_h2PNDm_IZgtz (
  double t,
  double x,
  double y,
  double z,
  double m1,
  double m2,
  double b)
{
return(
   all_h2PNDm_prematchIZgtt_func(t,x,y,z,m1,m2,b)   * all_h2PNDm_dTdt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dTdz_func(t,x,y,z,m1,m2,b) 
+ all_h2PNDm_prematchIZgxx_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dXdt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dXdz_func(t,x,y,z,m1,m2,b) 
+ all_h2PNDm_prematchIZgyy_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dYdt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dYdz_func(t,x,y,z,m1,m2,b) 
+ all_h2PNDm_prematchIZgzz_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dZdt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dZdz_func(t,x,y,z,m1,m2,b) 
+  all_h2PNDm_prematchIZgtx_func(t,x,y,z,m1,m2,b) * ( all_h2PNDm_dTdt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dXdz_func(t,x,y,z,m1,m2,b) + all_h2PNDm_dXdt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dTdz_func(t,x,y,z,m1,m2,b)) 
+  all_h2PNDm_prematchIZgty_func(t,x,y,z,m1,m2,b) * ( all_h2PNDm_dTdt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dYdz_func(t,x,y,z,m1,m2,b) + all_h2PNDm_dYdt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dTdz_func(t,x,y,z,m1,m2,b)) 
+  all_h2PNDm_prematchIZgtz_func(t,x,y,z,m1,m2,b) * ( all_h2PNDm_dTdt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dZdz_func(t,x,y,z,m1,m2,b) + all_h2PNDm_dZdt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dTdz_func(t,x,y,z,m1,m2,b)) 
+  all_h2PNDm_prematchIZgxy_func(t,x,y,z,m1,m2,b) * ( all_h2PNDm_dXdt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dYdz_func(t,x,y,z,m1,m2,b) + all_h2PNDm_dYdt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dXdz_func(t,x,y,z,m1,m2,b)) 
+  all_h2PNDm_prematchIZgxz_func(t,x,y,z,m1,m2,b) * ( all_h2PNDm_dXdt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dZdz_func(t,x,y,z,m1,m2,b) + all_h2PNDm_dZdt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dXdz_func(t,x,y,z,m1,m2,b)) 
+  all_h2PNDm_prematchIZgyz_func(t,x,y,z,m1,m2,b) * ( all_h2PNDm_dYdt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dZdz_func(t,x,y,z,m1,m2,b) + all_h2PNDm_dZdt_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dYdz_func(t,x,y,z,m1,m2,b)) 
);
}

double all_h2PNDm_IZgxx (
  double t,
  double x,
  double y,
  double z,
  double m1,
  double m2,
  double b)
{
return(
   all_h2PNDm_prematchIZgtt_func(t,x,y,z,m1,m2,b)   * all_h2PNDm_dTdx_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dTdx_func(t,x,y,z,m1,m2,b) 
+ all_h2PNDm_prematchIZgxx_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dXdx_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dXdx_func(t,x,y,z,m1,m2,b) 
+ all_h2PNDm_prematchIZgyy_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dYdx_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dYdx_func(t,x,y,z,m1,m2,b) 
+ all_h2PNDm_prematchIZgzz_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dZdx_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dZdx_func(t,x,y,z,m1,m2,b) 
+  all_h2PNDm_prematchIZgtx_func(t,x,y,z,m1,m2,b) * ( all_h2PNDm_dTdx_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dXdx_func(t,x,y,z,m1,m2,b) + all_h2PNDm_dXdx_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dTdx_func(t,x,y,z,m1,m2,b)) 
+  all_h2PNDm_prematchIZgty_func(t,x,y,z,m1,m2,b) * ( all_h2PNDm_dTdx_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dYdx_func(t,x,y,z,m1,m2,b) + all_h2PNDm_dYdx_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dTdx_func(t,x,y,z,m1,m2,b)) 
+  all_h2PNDm_prematchIZgtz_func(t,x,y,z,m1,m2,b) * ( all_h2PNDm_dTdx_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dZdx_func(t,x,y,z,m1,m2,b) + all_h2PNDm_dZdx_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dTdx_func(t,x,y,z,m1,m2,b)) 
+  all_h2PNDm_prematchIZgxy_func(t,x,y,z,m1,m2,b) * ( all_h2PNDm_dXdx_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dYdx_func(t,x,y,z,m1,m2,b) + all_h2PNDm_dYdx_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dXdx_func(t,x,y,z,m1,m2,b)) 
+  all_h2PNDm_prematchIZgxz_func(t,x,y,z,m1,m2,b) * ( all_h2PNDm_dXdx_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dZdx_func(t,x,y,z,m1,m2,b) + all_h2PNDm_dZdx_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dXdx_func(t,x,y,z,m1,m2,b)) 
+  all_h2PNDm_prematchIZgyz_func(t,x,y,z,m1,m2,b) * ( all_h2PNDm_dYdx_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dZdx_func(t,x,y,z,m1,m2,b) + all_h2PNDm_dZdx_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dYdx_func(t,x,y,z,m1,m2,b)) 
);
}

double all_h2PNDm_IZgxy (
  double t,
  double x,
  double y,
  double z,
  double m1,
  double m2,
  double b)
{
return(
   all_h2PNDm_prematchIZgtt_func(t,x,y,z,m1,m2,b)   * all_h2PNDm_dTdx_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dTdy_func(t,x,y,z,m1,m2,b) 
+ all_h2PNDm_prematchIZgxx_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dXdx_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dXdy_func(t,x,y,z,m1,m2,b) 
+ all_h2PNDm_prematchIZgyy_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dYdx_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dYdy_func(t,x,y,z,m1,m2,b) 
+ all_h2PNDm_prematchIZgzz_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dZdx_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dZdy_func(t,x,y,z,m1,m2,b) 
+  all_h2PNDm_prematchIZgtx_func(t,x,y,z,m1,m2,b) * ( all_h2PNDm_dTdx_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dXdy_func(t,x,y,z,m1,m2,b) + all_h2PNDm_dXdx_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dTdy_func(t,x,y,z,m1,m2,b)) 
+  all_h2PNDm_prematchIZgty_func(t,x,y,z,m1,m2,b) * ( all_h2PNDm_dTdx_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dYdy_func(t,x,y,z,m1,m2,b) + all_h2PNDm_dYdx_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dTdy_func(t,x,y,z,m1,m2,b)) 
+  all_h2PNDm_prematchIZgtz_func(t,x,y,z,m1,m2,b) * ( all_h2PNDm_dTdx_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dZdy_func(t,x,y,z,m1,m2,b) + all_h2PNDm_dZdx_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dTdy_func(t,x,y,z,m1,m2,b)) 
+  all_h2PNDm_prematchIZgxy_func(t,x,y,z,m1,m2,b) * ( all_h2PNDm_dXdx_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dYdy_func(t,x,y,z,m1,m2,b) + all_h2PNDm_dYdx_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dXdy_func(t,x,y,z,m1,m2,b)) 
+  all_h2PNDm_prematchIZgxz_func(t,x,y,z,m1,m2,b) * ( all_h2PNDm_dXdx_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dZdy_func(t,x,y,z,m1,m2,b) + all_h2PNDm_dZdx_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dXdy_func(t,x,y,z,m1,m2,b)) 
+  all_h2PNDm_prematchIZgyz_func(t,x,y,z,m1,m2,b) * ( all_h2PNDm_dYdx_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dZdy_func(t,x,y,z,m1,m2,b) + all_h2PNDm_dZdx_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dYdy_func(t,x,y,z,m1,m2,b)) 
);
}

double all_h2PNDm_IZgxz (
  double t,
  double x,
  double y,
  double z,
  double m1,
  double m2,
  double b)
{
return(
   all_h2PNDm_prematchIZgtt_func(t,x,y,z,m1,m2,b)   * all_h2PNDm_dTdx_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dTdz_func(t,x,y,z,m1,m2,b) 
+ all_h2PNDm_prematchIZgxx_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dXdx_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dXdz_func(t,x,y,z,m1,m2,b) 
+ all_h2PNDm_prematchIZgyy_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dYdx_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dYdz_func(t,x,y,z,m1,m2,b) 
+ all_h2PNDm_prematchIZgzz_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dZdx_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dZdz_func(t,x,y,z,m1,m2,b) 
+  all_h2PNDm_prematchIZgtx_func(t,x,y,z,m1,m2,b) * ( all_h2PNDm_dTdx_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dXdz_func(t,x,y,z,m1,m2,b) + all_h2PNDm_dXdx_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dTdz_func(t,x,y,z,m1,m2,b)) 
+  all_h2PNDm_prematchIZgty_func(t,x,y,z,m1,m2,b) * ( all_h2PNDm_dTdx_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dYdz_func(t,x,y,z,m1,m2,b) + all_h2PNDm_dYdx_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dTdz_func(t,x,y,z,m1,m2,b)) 
+  all_h2PNDm_prematchIZgtz_func(t,x,y,z,m1,m2,b) * ( all_h2PNDm_dTdx_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dZdz_func(t,x,y,z,m1,m2,b) + all_h2PNDm_dZdx_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dTdz_func(t,x,y,z,m1,m2,b)) 
+  all_h2PNDm_prematchIZgxy_func(t,x,y,z,m1,m2,b) * ( all_h2PNDm_dXdx_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dYdz_func(t,x,y,z,m1,m2,b) + all_h2PNDm_dYdx_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dXdz_func(t,x,y,z,m1,m2,b)) 
+  all_h2PNDm_prematchIZgxz_func(t,x,y,z,m1,m2,b) * ( all_h2PNDm_dXdx_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dZdz_func(t,x,y,z,m1,m2,b) + all_h2PNDm_dZdx_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dXdz_func(t,x,y,z,m1,m2,b)) 
+  all_h2PNDm_prematchIZgyz_func(t,x,y,z,m1,m2,b) * ( all_h2PNDm_dYdx_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dZdz_func(t,x,y,z,m1,m2,b) + all_h2PNDm_dZdx_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dYdz_func(t,x,y,z,m1,m2,b)) 
);
}

double all_h2PNDm_IZgyy (
  double t,
  double x,
  double y,
  double z,
  double m1,
  double m2,
  double b)
{
return(
   all_h2PNDm_prematchIZgtt_func(t,x,y,z,m1,m2,b)   * all_h2PNDm_dTdy_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dTdy_func(t,x,y,z,m1,m2,b) 
+ all_h2PNDm_prematchIZgxx_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dXdy_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dXdy_func(t,x,y,z,m1,m2,b) 
+ all_h2PNDm_prematchIZgyy_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dYdy_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dYdy_func(t,x,y,z,m1,m2,b) 
+ all_h2PNDm_prematchIZgzz_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dZdy_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dZdy_func(t,x,y,z,m1,m2,b) 
+  all_h2PNDm_prematchIZgtx_func(t,x,y,z,m1,m2,b) * ( all_h2PNDm_dTdy_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dXdy_func(t,x,y,z,m1,m2,b) + all_h2PNDm_dXdy_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dTdy_func(t,x,y,z,m1,m2,b)) 
+  all_h2PNDm_prematchIZgty_func(t,x,y,z,m1,m2,b) * ( all_h2PNDm_dTdy_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dYdy_func(t,x,y,z,m1,m2,b) + all_h2PNDm_dYdy_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dTdy_func(t,x,y,z,m1,m2,b)) 
+  all_h2PNDm_prematchIZgtz_func(t,x,y,z,m1,m2,b) * ( all_h2PNDm_dTdy_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dZdy_func(t,x,y,z,m1,m2,b) + all_h2PNDm_dZdy_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dTdy_func(t,x,y,z,m1,m2,b)) 
+  all_h2PNDm_prematchIZgxy_func(t,x,y,z,m1,m2,b) * ( all_h2PNDm_dXdy_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dYdy_func(t,x,y,z,m1,m2,b) + all_h2PNDm_dYdy_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dXdy_func(t,x,y,z,m1,m2,b)) 
+  all_h2PNDm_prematchIZgxz_func(t,x,y,z,m1,m2,b) * ( all_h2PNDm_dXdy_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dZdy_func(t,x,y,z,m1,m2,b) + all_h2PNDm_dZdy_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dXdy_func(t,x,y,z,m1,m2,b)) 
+  all_h2PNDm_prematchIZgyz_func(t,x,y,z,m1,m2,b) * ( all_h2PNDm_dYdy_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dZdy_func(t,x,y,z,m1,m2,b) + all_h2PNDm_dZdy_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dYdy_func(t,x,y,z,m1,m2,b)) 
);
}

double all_h2PNDm_IZgyz (
  double t,
  double x,
  double y,
  double z,
  double m1,
  double m2,
  double b)
{
return(
   all_h2PNDm_prematchIZgtt_func(t,x,y,z,m1,m2,b)   * all_h2PNDm_dTdy_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dTdz_func(t,x,y,z,m1,m2,b) 
+ all_h2PNDm_prematchIZgxx_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dXdy_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dXdz_func(t,x,y,z,m1,m2,b) 
+ all_h2PNDm_prematchIZgyy_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dYdy_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dYdz_func(t,x,y,z,m1,m2,b) 
+ all_h2PNDm_prematchIZgzz_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dZdy_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dZdz_func(t,x,y,z,m1,m2,b) 
+  all_h2PNDm_prematchIZgtx_func(t,x,y,z,m1,m2,b) * ( all_h2PNDm_dTdy_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dXdz_func(t,x,y,z,m1,m2,b) + all_h2PNDm_dXdy_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dTdz_func(t,x,y,z,m1,m2,b)) 
+  all_h2PNDm_prematchIZgty_func(t,x,y,z,m1,m2,b) * ( all_h2PNDm_dTdy_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dYdz_func(t,x,y,z,m1,m2,b) + all_h2PNDm_dYdy_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dTdz_func(t,x,y,z,m1,m2,b)) 
+  all_h2PNDm_prematchIZgtz_func(t,x,y,z,m1,m2,b) * ( all_h2PNDm_dTdy_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dZdz_func(t,x,y,z,m1,m2,b) + all_h2PNDm_dZdy_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dTdz_func(t,x,y,z,m1,m2,b)) 
+  all_h2PNDm_prematchIZgxy_func(t,x,y,z,m1,m2,b) * ( all_h2PNDm_dXdy_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dYdz_func(t,x,y,z,m1,m2,b) + all_h2PNDm_dYdy_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dXdz_func(t,x,y,z,m1,m2,b)) 
+  all_h2PNDm_prematchIZgxz_func(t,x,y,z,m1,m2,b) * ( all_h2PNDm_dXdy_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dZdz_func(t,x,y,z,m1,m2,b) + all_h2PNDm_dZdy_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dXdz_func(t,x,y,z,m1,m2,b)) 
+  all_h2PNDm_prematchIZgyz_func(t,x,y,z,m1,m2,b) * ( all_h2PNDm_dYdy_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dZdz_func(t,x,y,z,m1,m2,b) + all_h2PNDm_dZdy_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dYdz_func(t,x,y,z,m1,m2,b)) 
);
}

double all_h2PNDm_IZgzz (
  double t,
  double x,
  double y,
  double z,
  double m1,
  double m2,
  double b)
{
return(
   all_h2PNDm_prematchIZgtt_func(t,x,y,z,m1,m2,b)   * all_h2PNDm_dTdz_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dTdz_func(t,x,y,z,m1,m2,b) 
+ all_h2PNDm_prematchIZgxx_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dXdz_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dXdz_func(t,x,y,z,m1,m2,b) 
+ all_h2PNDm_prematchIZgyy_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dYdz_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dYdz_func(t,x,y,z,m1,m2,b) 
+ all_h2PNDm_prematchIZgzz_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dZdz_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dZdz_func(t,x,y,z,m1,m2,b) 
+  all_h2PNDm_prematchIZgtx_func(t,x,y,z,m1,m2,b) * ( all_h2PNDm_dTdz_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dXdz_func(t,x,y,z,m1,m2,b) + all_h2PNDm_dXdz_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dTdz_func(t,x,y,z,m1,m2,b)) 
+  all_h2PNDm_prematchIZgty_func(t,x,y,z,m1,m2,b) * ( all_h2PNDm_dTdz_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dYdz_func(t,x,y,z,m1,m2,b) + all_h2PNDm_dYdz_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dTdz_func(t,x,y,z,m1,m2,b)) 
+  all_h2PNDm_prematchIZgtz_func(t,x,y,z,m1,m2,b) * ( all_h2PNDm_dTdz_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dZdz_func(t,x,y,z,m1,m2,b) + all_h2PNDm_dZdz_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dTdz_func(t,x,y,z,m1,m2,b)) 
+  all_h2PNDm_prematchIZgxy_func(t,x,y,z,m1,m2,b) * ( all_h2PNDm_dXdz_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dYdz_func(t,x,y,z,m1,m2,b) + all_h2PNDm_dYdz_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dXdz_func(t,x,y,z,m1,m2,b)) 
+  all_h2PNDm_prematchIZgxz_func(t,x,y,z,m1,m2,b) * ( all_h2PNDm_dXdz_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dZdz_func(t,x,y,z,m1,m2,b) + all_h2PNDm_dZdz_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dXdz_func(t,x,y,z,m1,m2,b)) 
+  all_h2PNDm_prematchIZgyz_func(t,x,y,z,m1,m2,b) * ( all_h2PNDm_dYdz_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dZdz_func(t,x,y,z,m1,m2,b) + all_h2PNDm_dZdz_func(t,x,y,z,m1,m2,b) * all_h2PNDm_dYdz_func(t,x,y,z,m1,m2,b)) 
);
}
