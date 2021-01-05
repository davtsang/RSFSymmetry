//module for codes from numerical recipies in c
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define TINY 1.0E-30

//function to advance vars[] by 1 step of size dr  (NR in C)
void RK4step(double vars[], int n, double r, double dr, double var_out[], void (*derivatives)(double, double[], double[])){
  double variables[n], d_vars1[n], d_vars2[n], d_vars3[n];
  int i;
  for(i=0;i<n;i++){variables[i]=vars[i];}
  (*derivatives)(r,variables,d_vars1);
  for(i=0;i<n;i++){variables[i]=vars[i]+(dr/2.0)*d_vars1[i];}
  (*derivatives)(r+(dr/2.0),variables,d_vars2);
  for(i=0;i<n;i++){variables[i]=vars[i]+(dr/2.0)*d_vars2[i];}
  (*derivatives)(r+(dr/2.0),variables,d_vars3);
  for(i=0;i<n;i++){
    variables[i]=vars[i]+dr*d_vars3[i];
    d_vars3[i]=d_vars3[i]+d_vars2[i];
  }
  (*derivatives)(r+dr,variables,d_vars2);
  for(i=0;i<n;i++){var_out[i]=vars[i]+(dr/6.0)*(d_vars1[i]+d_vars2[i]+2.0*d_vars3[i]);}
}

//function to solve RK5 and embedded RK4 using the same 6 function evaluations  (NR in C)
void RKCKstep(double vars[], int n, double r, double dr, double var_out[], double err[], void (*derivatives)(double, double[], double[])){
  double variables[n], k1[n], k2[n], k3[n], k4[n], k5[n], k6[n];
  double const a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875;
  double const b21=0.2,b31=0.075,b32=0.225,b41=0.3,b42=-0.9,b43=1.2,b51=-11.0/54.0,b52=2.5,b53=-70.0/27.0,b54=35.0/27.0;
  double const b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,b64=44275.0/110592.0,b65=253.0/4096.0;
  double const c1=37.0/378.0,c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0;
  double dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,dc4=c4-13525.0/55296.0,dc5=-277.0/14336.0,dc6=c6-0.25;
  int i;
  for(i=0;i<n;i++){variables[i]=vars[i];}
  (*derivatives)(r,variables,k1);
  for(i=0;i<n;i++){variables[i]=vars[i]+dr*(b21*k1[i]);}
  (*derivatives)(r+a2*dr,variables,k2);
  for(i=0;i<n;i++){variables[i]=vars[i]+dr*(b31*k1[i]+b32*k2[i]);}
  (*derivatives)(r+a3*dr,variables,k3);
  for(i=0;i<n;i++){variables[i]=vars[i]+dr*(b41*k1[i]+b42*k2[i]+b43*k3[i]);}
  (*derivatives)(r+a4*dr,variables,k4);
  for(i=0;i<n;i++){variables[i]=vars[i]+dr*(b51*k1[i]+b52*k2[i]+b53*k3[i]+b54*k4[i]);}
  (*derivatives)(r+a5*dr,variables,k5);
  for(i=0;i<n;i++){variables[i]=vars[i]+dr*(b61*k1[i]+b62*k2[i]+b63*k3[i]+b64*k4[i]+b65*k5[i]);}
  (*derivatives)(r+a6*dr,variables,k6);
  for(i=0;i<n;i++){
    var_out[i]=vars[i]+dr*(c1*k1[i]+c3*k3[i]+c4*k4[i]+c6*k6[i]);
    err[i]=dr*(dc1*k1[i]+dc3*k3[i]+dc4*k4[i]+dc5*k5[i]+dc6*k6[i]);
  }
}

//function to advance vars[] by 1 step of a size that causes error of less than tol*vars[]  (NR in C)
void RK_var_step(double vars[], int n, double *r, double dr_try, double tol, double *dr_used, double *dr_next, void (*derivatives)(double, double[], double[])){
  int i,j;
  double errmax,errcon=1.89E-4,safety=0.9,dr,drtemp,rnew,err[n],vartemp[n],varscale[n];
  dr=dr_try;
  for (;;){
    RKCKstep(vars,n,*r,dr,vartemp,err,derivatives);
    errmax=0.0;
    for(i=0;i<n;i++){
     // varscale[i]=vars[i]+dr*vartemp[i]+TINY;
      varscale[i]=vartemp[i]+TINY;
      errmax=fmax(errmax,fabs(err[i]/varscale[i]));
    }
    errmax=errmax/tol;
    if (errmax <= 1){break;}
    drtemp=safety*dr*pow(errmax,-0.25);    //safety used, since error may not scale linearly with dr
    if (dr >= 0.0) {
      dr=fmax(drtemp,0.1*dr);
    } else if (dr < 0.0) {
      dr=fmin(drtemp,0.1*dr);
    }
    rnew=(*r)+dr;
    if (rnew == (*r)) { printf("error, dr is too small\nerror at x = %le\n",*r);exit(0);}  //  dr=dr*100000.0;rnew=(*r)+dr;break;} //
  }
  if (errmax > errcon) {
    *dr_next=safety*dr*pow(errmax,-0.2);
  } else {
    *dr_next=5.0*dr;
  }
  *r=(*r)+dr;
  *dr_used=dr;
  for (i=0;i<n;i++) {vars[i]=vartemp[i];}
}


//find d^2_y/d_x^2 using numerical recipies in c code (chapter 3.3)
//assumes natural cubic splines
//solves tridiagonal matrix for  d^2_y/d_x^2
void cubic_spline(int size, double y[size], double x[size], double d_y[size]){
  int i;
  double sig,t;
  double u[size];

  d_y[0]=0.0;
  d_y[size-1]=0.0;
  u[0]=0.0;
  u[size-1]=0.0;
  for (i=1;i<=size-2;i++){
    sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
    t=sig*d_y[i-1]+2.0;
    d_y[i]=(sig-1.0)/t;
    u[i]=(y[i+1]-y[i])/(x[i+1]-x[i])-(y[i]-y[i-1])/(x[i]-x[i-1]);
    u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/t;
  }
  for (i=size-2;i>=0;i--){
    d_y[i]=d_y[i]*d_y[i+1]+u[i];
  }
}

//find y value at given x value using cubic splines
//using numerical recipies in c code (chapter 3.3)
double spline_value(int size, double y[size], double x[size],  double d_y[size], double x_val){
  int klo,khi,k;
  double h,b,a,output_value;

  klo=0;
  khi=size-1;
  while (khi-klo > 1) {
    k=(khi+klo) >> 1;
    if (x[k] > x_val) khi=k;
    else klo=k;
  } 
  h=x[khi]-x[klo];
  if (h==0.0){
    printf("error in spline_value, h==0 %d  %d  %le  %le  %le\n",khi,klo,x[khi],x[klo],x_val);
    if((x[khi]==x_val) && (x[klo]==x_val)){return (y[klo]+y[khi])/2.0;}
    printf("exiting\n");exit(0);
  }

  a=(x[khi]-x_val)/h;
  b=(x_val-x[klo])/h;
  output_value=a*y[klo]+b*y[khi]+((pow(a,3)-a)*d_y[klo]+(pow(b,3)-b)*d_y[khi])*(h*h)/6.0;
  if (y[klo]==y[khi]) {output_value=y[klo];}
  return output_value;
}

//function to linearly interpolate at a given value
double linear_int(int size,double value,double dep_array[],double ind_arr[]){
  double output,ratio;
  int klo,khi,k;
  klo=0;
  khi=size-1;
  while (khi-klo > 1) {
    k=(khi+klo) >> 1;
    if (ind_arr[k] > value) khi=k;
    else klo=k;
  } 

  if (khi == klo) {printf("error,khi=klo\n"); exit(0);}
  ratio=(value-ind_arr[klo])/(ind_arr[khi]-ind_arr[klo]);
  output=dep_array[klo]+ratio*(dep_array[khi]-dep_array[klo]);
  return output;
}
