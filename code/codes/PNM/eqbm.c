//module containing functions for equilibrium properties of the star
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "nr_codes.h"
#define PI 3.14159265359
#define C 2.99792458E+10
#define HBAR 1.0545726663E-27
#define MN 1.674928610E-24
#define G 6.6725985E-8
#define TOL 1.0E-7
#define E_CHARGE 4.8032042510E-10
#define KBOLTZ 1.380649E-16
extern double radius, NB_CC, mu_A, mu_B;
extern double *r_pos, *d_mass_r, *mass_r, *density, *d_density, *pressure, *d_pressure, *bary_dens, *d_bary_dens, *smod_arr, *d_smod;
extern int rsteps, lines;
extern double *Z, *A, *nb_ZA, *ion_ratio, *d_Z, *d_A, *d_ion_ratio;

void smooth_array(int size, double array[size]){
  int i;
  double temp_arr[size];

  for (i=3;i<size-3;i++){
    temp_arr[i]=array[i-3]+array[i+3]+2.0*array[i-2]+2.0*array[i+2]+4.0*array[i-1]+4.0*array[i+1]+8.0*array[i];
    temp_arr[i]=temp_arr[i]/22.0;
  }
  temp_arr[0]=(array[i+3]+2.0*array[i+2]+4.0*array[i+1]+8.0*array[i])/15.0;
  temp_arr[1]=(4.0*array[i-1]+array[i+3]+2.0*array[i+2]+4.0*array[i+1]+8.0*array[i])/19.0;
  temp_arr[2]=(2.0*array[i-2]+4.0*array[i-1]+array[i+3]+2.0*array[i+2]+4.0*array[i+1]+8.0*array[i])/21.0;
  temp_arr[size-3]=(2.0*array[i-2]+4.0*array[i-1]+array[i-3]+2.0*array[i+2]+4.0*array[i+1]+8.0*array[i])/21.0;
  temp_arr[size-2]=(4.0*array[i-1]+array[i-3]+2.0*array[i-2]+4.0*array[i+1]+8.0*array[i])/19.0;
  temp_arr[size-1]=(array[i-3]+2.0*array[i-2]+4.0*array[i-1]+8.0*array[i])/15.0;

  for (i=0;i<size;i++){
    array[i]=temp_arr[i];
  }
}

double Utildevalue(double r){
  double densval,massval,output;
  densval=spline_value(rsteps,density,r_pos,d_density,r);
  massval=spline_value(rsteps,mass_r,r_pos,d_mass_r,r);
  output=(4.0*PI*r*r*densval)*(r/massval);
  return output;
}

double Vtildevalue(double r){
  double densval,massval,presval,output;
  densval=spline_value(rsteps,density,r_pos,d_density,r);
  presval=spline_value(rsteps,pressure,r_pos,d_pressure,r);
  massval=spline_value(rsteps,mass_r,r_pos,d_mass_r,r);

  output=-(G*densval*massval)/(r*r);
  output=output*(1.0+presval/(densval*C*C));
  output=output*(1.0+((4.0*PI*r*r*r*presval)/(massval*C*C)));
  output=output*(1.0/(1.0-(2.0*G*massval)/(r*C*C)));

  output=-output*(r/presval);
  return output;
}

double c1value(double r){
  double massval,output;
  massval=spline_value(rsteps,mass_r,r_pos,d_mass_r,r);
  output=pow(r/radius,3)*(mass_r[rsteps-1]/massval);
  return output;
}

//shear modulus -> Strohmayer 1991
double mu_val(double r){
  double temperature,a,proton,output,nucleon,n_i,gamma,bary_dens_val,ratio;
  bary_dens_val=spline_value(rsteps,bary_dens,r_pos,d_bary_dens,r);
  if (bary_dens_val>NB_CC){return 0.0;}
  output=spline_value(rsteps,smod_arr,r_pos,d_smod,r);
  return output;
}




//method for getting g1 taken from David's code, basically just taking the derivative of logp(logrho) rather than splitting it into rho/p*dp/drho
void setup_g1(double g1[],int lines,double rho_arr[],double pres_arr[],double EoS_bary_dens[]){
  int i;
  double j=0.5;  //0.5 to 1.0 seems okay-ish
  double *logp=malloc(lines*sizeof(double));
  double *logrho=malloc(lines*sizeof(double));
  double *d_logp=malloc(lines*sizeof(double));
  double *tempg1=malloc(lines*sizeof(double));
  double *d_g1=malloc(lines*sizeof(double));

  for (i=0;i<lines;i++){
    logp[i]=log(pres_arr[i]);
    logrho[i]=log(EoS_bary_dens[i]);
  }
  cubic_spline(lines,logp,logrho,d_logp);

  for (i=0;i<lines;i++){
    if (logrho[i]-j < logrho[0]){
      tempg1[i]=(spline_value(lines,logp,logrho,d_logp,logrho[i]+j)-spline_value(lines,logp,logrho,d_logp,logrho[i]))/(j);
    }else if (logrho[i]+j > logrho[lines-1]){
      tempg1[i]=(spline_value(lines,logp,logrho,d_logp,logrho[i])-spline_value(lines,logp,logrho,d_logp,logrho[i]-j))/(j);
    }else{
      tempg1[i]=(spline_value(lines,logp,logrho,d_logp,logrho[i]+j)-spline_value(lines,logp,logrho,d_logp,logrho[i]-j))/(2.0*j);
    }
  }

  for (i=0;i<lines;i++){d_g1[i]=0.0;}  //set to zero to use linear splines
  for (i=0;i<rsteps;i++){
    g1[i]=spline_value(lines,tempg1,EoS_bary_dens,d_g1,bary_dens[i]);//tempg1[i];
  }
  free(tempg1);
}

