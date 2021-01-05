//solves TOV for a 1.4M_sun neutron star. Modified to use newton's EoSs as inputs (9/6/20)

//Variation used in this code: the shear modulus is taken from a different EoS (the mu EoS) than the rest of the EoS properties (the main EoS). 
//This allows us to test the impact of the shear modulus on the mode, since if the resulting freqeuncy is closer to the mu EoS, the shear modulus has a strong impact on the frequency.
//On the other hand, if the frequency is similar to that of the main EoS, the other properties are more important.

//This code is easy to modify to use either the crust-core transition density of the main EoS or the mu EoS. If using the main EoS, the density and pressure of the mu EoS will be rescaled so
//that they have the same values at the transtion as the main EoS

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "nr_codes.h"
#define PI 3.14159265359
#define C 2.99792458E+10
#define G 6.6725985E-8
#define TOL 1E-10
#define HBAR 1.0545726663E-27
#define MN 1.674928610E-24

double dm_dr( double rho, double r);
double dp_dr( double rho, double p, double r, double m);
void derivatives(double rpos,double var_in[],double d_var_out[]);
void smooth_array(int size, double array[size]);

int ARRAY_SIZE,lines_core,lines_crust;
double *rho_arr,*pres_arr,*d_rho,*bary_arr,*d_bary,*smod_arr,*d_smod,*gamma1_arr,*d_gamma1,*pres_arr_2,*smod_arr_2,*d_smod_2,*rho_arr_2,*bary_arr_2,*d_bary_2;
double *core_rho,*crust_rho,*core_bary,*crust_bary,*d_rho_crust,*d_rho_core,*d_bary_crust,*d_bary_core,*pres_core,*pres_crust;

void main (int argc, char *argv[]) {
  double p_c,dp_c,p,dr,dp,r,e_dens,nb,rho,m,k1m,k2m,k1p,k2p,k3m,k4m,k3p,k4p,p_old,bary,mu,nbcc,M_target,gamma1;
  int i,lines,linesmu;
  char c;
  char eosname[80]="",eosmuname[80]="";
  double variables[2],drused;
  int *Z=malloc(0*sizeof(int)),*A=malloc(0*sizeof(int));

  if (argc != 8) {printf("Input: ./star [J (i.0)] [L (j.00)] [K (k.00)] [M* (solar masses)] [J (mu EoS) (i.0)] [L (mu EoS) (j.00)] [K (mu EoS) (k.00)]\n");exit(0);}
/*-----------------------------setup output files and read the input eos file-----------------------------*/
  FILE * mass_output = fopen("NS_mass.txt", "w");
  FILE * pres_output = fopen("NS_pres.txt", "w");
  FILE * dens_output = fopen("NS_dens.txt", "w");
  FILE * grav_output = fopen("NS_grav.txt", "w");
  FILE * dpdr_output = fopen("NS_dpdr.txt", "w");
  FILE * bary_output = fopen("NS_bary.txt", "w");
  FILE * smod_output = fopen("NS_smod.txt", "w");
  FILE * gam1_output = fopen("NS_gam1.txt", "w");

  //get eos file from input J,L values
  strcat(eosname,"EoSs/eos_LQDROP_SkXi450_");
  strcat(eosname,argv[1]);
  strcat(eosname,"_");
  strcat(eosname,argv[2]);
  strcat(eosname,"_");
  strcat(eosname,argv[3]);
  strcat(eosname,"_2.2_Imid_.dat");
  FILE * eos = fopen(eosname, "r");

  M_target=atof(argv[4]);

  lines=0;
  c = getc(eos);
  while (c != EOF){
    if (c == '\n'){
      lines=lines+1;
    }
    c = getc(eos);
  }
  fseek(eos,0,SEEK_SET);
  ARRAY_SIZE=lines;

  rho_arr=realloc(rho_arr,ARRAY_SIZE*sizeof(double));
  pres_arr=realloc(pres_arr,ARRAY_SIZE*sizeof(double));
  bary_arr=realloc(bary_arr,ARRAY_SIZE*sizeof(double));
  smod_arr=realloc(smod_arr,ARRAY_SIZE*sizeof(double));
  gamma1_arr=realloc(gamma1_arr,ARRAY_SIZE*sizeof(double));
  d_rho=realloc(d_rho,ARRAY_SIZE*sizeof(double));
  d_bary=realloc(d_bary,ARRAY_SIZE*sizeof(double));
  d_smod=realloc(d_smod,ARRAY_SIZE*sizeof(double));
  d_gamma1=realloc(d_gamma1,ARRAY_SIZE*sizeof(double));

  //read values from EoS file
  for(i=ARRAY_SIZE-1;i>=0;i--){
    fscanf(eos,"%le  %le  %le  %le  %le",&e_dens,&p,&nb,&gamma1,&mu);
    if (p>pres_arr[i+1] && i!=ARRAY_SIZE-1) {printf("WARNING: -ve dp/dr occurs in %s\n",eosname);}
    rho_arr[i]=e_dens;  //  column 1 of file = energy density/(C*C)   (=rho at NR limit)
    pres_arr[i]=p;       //  column 2 of file = pressure
    bary_arr[i]=nb;       //  column 3 of file = baryon density
    smod_arr[i]=mu;
    gamma1_arr[i]=gamma1;
  }
  fclose(eos);
/*-----------------------------setup output files and read the input eos file-----------------------------*/

/*-----------------------change the shear modulus to be that of the second EoS file-----------------------*/
  double new_mu_rhocc,new_mu_prescc,new_mu_barycc;

  //get EoS to take mu from
  strcat(eosmuname,"EoSs/eos_LQDROP_SkXi450_");
  strcat(eosmuname,argv[5]);
  strcat(eosmuname,"_");
  strcat(eosmuname,argv[6]);
  strcat(eosmuname,"_");
  strcat(eosmuname,argv[7]);
  strcat(eosmuname,"_2.2_Imid_.dat");
  FILE * eosmu = fopen(eosmuname, "r");

  linesmu=0;
  c = getc(eosmu);
  while (c != EOF){
    if (c == '\n'){
      linesmu=linesmu+1;
    }
    c = getc(eosmu);
  }
  fseek(eosmu,0,SEEK_SET);

  pres_arr_2=realloc(pres_arr_2,linesmu*sizeof(double));
  rho_arr_2=realloc(rho_arr_2,linesmu*sizeof(double));
  smod_arr_2=realloc(smod_arr_2,linesmu*sizeof(double));
  d_smod_2=realloc(d_smod_2,linesmu*sizeof(double));
  bary_arr_2=realloc(bary_arr_2,linesmu*sizeof(double));
  d_bary_2=realloc(d_bary_2,linesmu*sizeof(double));

  for(i=linesmu-1;i>=0;i--){
    fscanf(eosmu,"%le  %le  %le  %le  %le",&e_dens,&p,&nb,&gamma1,&mu);
    pres_arr_2[i]=p;
    bary_arr_2[i]=nb;
    rho_arr_2[i]=e_dens;
    smod_arr_2[i]=mu;
  }
  fclose(eosmu);

  //find the crust core transition density in the mu EoS file
  for (i=0;i<linesmu;i++){
    if (smod_arr_2[i]==0){
      new_mu_rhocc=rho_arr_2[i-1];
      new_mu_prescc=rho_arr_2[i-1];
      new_mu_barycc=rho_arr_2[i-1];
      break;
    }
  }

  //use cubic spines for mu, and scale the density in the mu EoS file so that rho_cc is the same as in the main EoS file, while keeping the lowest rho value constant
  //comment out only this bit to get unscaled mu EoS rho values, but still end the crust at the main file's rho_cc
  //(^warning: will result in some of the shear modulus being cut off, or part of the crust having zero shear modulus)
  //always comment this out if you are using the mu EoSs transition point
/*
  for (i=0;i<linesmu;i++){
    d_smod[i]=0;
    rho_arr_2[i]=crust_rho[0]+(crust_rho[lines_crust-1]-crust_rho[0])*((rho_arr_2[i]-rho_arr_2[0])/(new_mu_rhocc-rho_arr_2[0]));
    pres_arr_2[i]=crust_pres[0]+(crust_pres[lines_crust-1]-crust_pres[0])*((pres_arr_2[i]-pres_arr_2[0])/(new_mu_prescc-pres_arr_2[0]));
    bary_arr_2[i]=crust_bary[0]+(crust_bary[lines_crust-1]-crust_bary[0])*((bary_arr_2[i]-bary_arr_2[0])/(new_mu_barycc-bary_arr_2[0]));
  }
*/
/*-----------------------change the shear modulus to be that of the second EoS file-----------------------*/

/*-------------modify eos data to fix repeated lines, which cause errors in the cubic splines-------------*/
  for (i=2;i<ARRAY_SIZE;i++){
    if ((rho_arr[i]==rho_arr[i-1]) && (rho_arr[i]==rho_arr[i-2])){
      rho_arr[i-2]=(rho_arr[i-3]+2.0*rho_arr[i-2])/3.0;
      rho_arr[i]=(rho_arr[i+1]+2.0*rho_arr[i])/3.0;
    }  else if (rho_arr[i]==rho_arr[i-1]  && (rho_arr[i]!=rho_arr[i-2])) {
      rho_arr[i]=(rho_arr[i+1]+2.0*rho_arr[i])/3.0;
      rho_arr[i-1]=(rho_arr[i-2]+2.0*rho_arr[i-1])/3.0;
    }
  }
  for (i=2;i<ARRAY_SIZE;i++){
    if ((pres_arr[i]==pres_arr[i-1]) && (pres_arr[i]==pres_arr[i-2])){
      pres_arr[i-2]=(pres_arr[i-3]+2.0*pres_arr[i-2])/3.0;
      pres_arr[i]=(pres_arr[i+1]+2.0*pres_arr[i])/3.0;
    }  else if (pres_arr[i]==pres_arr[i-1]  && (pres_arr[i]!=pres_arr[i-2])) {
      pres_arr[i]=(pres_arr[i+1]+2.0*pres_arr[i])/3.0;
      pres_arr[i-1]=(pres_arr[i-2]+2.0*pres_arr[i-1])/3.0;
    }
  }
  for (i=2;i<ARRAY_SIZE;i++){
    if ((bary_arr[i]==bary_arr[i-1]) && (bary_arr[i]==bary_arr[i-2])){
      bary_arr[i-2]=(bary_arr[i-3]+2.0*bary_arr[i-2])/3.0;
      bary_arr[i]=(bary_arr[i+1]+2.0*bary_arr[i])/3.0;
    }  else if (bary_arr[i]==bary_arr[i-1]  && (bary_arr[i]!=bary_arr[i-2])) {
      bary_arr[i]=(bary_arr[i+1]+2.0*bary_arr[i])/3.0;
      bary_arr[i-1]=(bary_arr[i-2]+2.0*bary_arr[i-1])/3.0;
    }
  }
  for (i=1;i<ARRAY_SIZE;i++){if ((bary_arr[i]==bary_arr[i-1])){printf("nb   %d\n",i);}}
  for (i=1;i<ARRAY_SIZE;i++){if ((pres_arr[i]==pres_arr[i-1])){printf("p    %d\n",i);}}
  for (i=1;i<ARRAY_SIZE;i++){if ((rho_arr[i]==rho_arr[i-1])){printf("rho  %d\n",i);}}


  for (i=2;i<linesmu;i++){
    if (smod_arr_2[i]==0){  //get the crust-core transition baryon density   //MODIFIED FOR THE mu EOS FILE!
      nbcc=bary_arr_2[i-1];
      break;
    }
  }

  //set up cubic splines
  //cubic_spline(ARRAY_SIZE,smod_arr,pres_arr,d_smod);
  for (i=0;i<ARRAY_SIZE;i++){d_smod[i]=0;}
  cubic_spline(ARRAY_SIZE,rho_arr,pres_arr,d_rho);
  cubic_spline(ARRAY_SIZE,bary_arr,pres_arr,d_bary);
  cubic_spline(ARRAY_SIZE,gamma1_arr,pres_arr,d_gamma1);
  for (i=0;i<ARRAY_SIZE;i++){d_gamma1[i]=0;}
/*-------------modify eos data to fix repeated lines, which cause errors in the cubic splines-------------*/



/*----------------------use one spline in the core, and a different one in the crust----------------------*/

  //Seperate splines in core and crust (at mu=0 and mu!=0)
  lines_core=0; 
  lines_crust=0;

  for (i=0;i<ARRAY_SIZE;i++){
    if (smod_arr[i]==0){lines_core=lines_core+1;}
  }
  lines_crust=ARRAY_SIZE-lines_core;

  core_rho=realloc(core_rho,lines_core*sizeof(double));
  crust_rho=realloc(crust_rho,lines_crust*sizeof(double));
  core_bary=realloc(core_bary,lines_core*sizeof(double));
  crust_bary=realloc(crust_bary,lines_crust*sizeof(double));
  d_rho_core=realloc(d_rho_core,lines_core*sizeof(double));
  d_rho_crust=realloc(d_rho_crust,lines_crust*sizeof(double));
  d_bary_core=realloc(d_bary_core,lines_core*sizeof(double));
  d_bary_crust=realloc(d_bary_crust,lines_crust*sizeof(double));
  pres_core=realloc(pres_core,lines_core*sizeof(double));
  pres_crust=realloc(pres_crust,lines_crust*sizeof(double));

  for (i=0;i<ARRAY_SIZE;i++){
    if (i<lines_crust){
      crust_rho[i]=rho_arr[i];
      crust_bary[i]=bary_arr[i];
      pres_crust[i]=pres_arr[i];
    } else {
      core_rho[i-lines_crust]=rho_arr[i];
      core_bary[i-lines_crust]=bary_arr[i];
      pres_core[i-lines_crust]=pres_arr[i];
    }
  }

  cubic_spline(lines_core,core_rho,pres_core,d_rho_core);
  cubic_spline(lines_crust,crust_rho,pres_crust,d_rho_crust);
  cubic_spline(lines_core,core_bary,pres_core,d_bary_core);
  cubic_spline(lines_crust,crust_bary,pres_crust,d_bary_crust);

  //if (atof(argv[3])<-250){
    for (i=0;i<ARRAY_SIZE;i++){d_rho_core[i]=0;}
    for (i=0;i<ARRAY_SIZE;i++){d_bary_core[i]=0;}   //use linear splines for large Ksym
  //}

/*----------------------use one spline in the core, and a different one in the crust----------------------*/

/*-----------loop over pressure values in the TOV until to get a total mass of 1.4 solar masses-----------*/
  p_c=1.0e+34;
  dp_c=1.0e+34;
  do{
    m=0.0;
    r=0.0;
    p=p_c;
    dr=0.1;
    drused=dr;

    //loop over RK4
    do {
      variables[0]=p;
      variables[1]=m;
      RK_var_step(variables,2,&r,dr,TOL,&drused,&dr,derivatives);
      p=variables[0];
      m=variables[1];
      
      // check for errors
      if (p>=pres_arr[ARRAY_SIZE-lines_core]) {
        rho=spline_value(lines_core,core_rho,pres_core,d_rho_core,p);
        bary=spline_value(lines_core,core_bary,pres_core,d_bary_core,p);
      } else if (p<=pres_arr[ARRAY_SIZE-lines_core-1]) {
        rho=spline_value(lines_crust,crust_rho,pres_crust,d_rho_crust,p);
        bary=spline_value(lines_crust,crust_bary,pres_crust,d_bary_crust,p);
      } else if ((p>pres_arr[ARRAY_SIZE-lines_core-1]) && (p<pres_arr[ARRAY_SIZE-lines_core])) {
        rho=linear_int(ARRAY_SIZE,p,rho_arr,pres_arr);
        bary=linear_int(ARRAY_SIZE,p,bary_arr,pres_arr);
      }
      mu=spline_value(linesmu,smod_arr_2,pres_arr_2,d_smod_2,p);
      gamma1=spline_value(ARRAY_SIZE,gamma1_arr,pres_arr,d_gamma1,p);
      if((m<0) || (p<0) || (rho<0) || ((G*m)/(r*r)<0) || (bary<0) || (mu<0) || (gamma1<0)){printf("uh-oh, negative value at r=%le    %le %le %le %le %le %le\n",r,m,p,rho,bary,mu,gamma1);}
    } while (p > pres_arr[1]);  //leave loop if P is below the boundary (what is the best boundary?)
  
    rho=spline_value(lines_core,core_rho,pres_core,d_rho_core,p_c);  //don't worry about crust vs core spline, this is definately in the core
//    printf("neutron star of mass %lf M_sun and radius %lf Km, calcualted using central pressure %.6e and central density %.6e.\n",m/(1.989E+33),r/(1.0E+05),p_c,rho);
  
    if (m>(M_target*(1.989E+33))){p_c=p_c-dp_c;dp_c=dp_c*0.1;}
    p_c=p_c+dp_c;
  }while (dp_c>1e+29);
/*-----------loop over pressure values in the TOV until to get a total mass of 1.4 solar masses-----------*/

/*----------------solve the TOV equations for 1.4 solar mass star, while printing out data----------------*/
  p_c=p_c-dp_c;
  m=0.0;
  r=0.0;
  p=p_c;
  p_old=p;
  dr=0.1;
  drused=dr;

  //printf out initial values for mass, pressure, etc...
  fprintf(mass_output,"%.15le  %.15le\n",0.0,m);
  fprintf(pres_output,"%.15le  %.15le\n",0.0,p);
  fprintf(dens_output,"%.15le  %.15le\n",0.0,spline_value(ARRAY_SIZE,rho_arr,pres_arr,d_rho,p));
  fprintf(grav_output,"%.15le  %.15le\n",0.0,0.0);
  fprintf(dpdr_output,"%.15le  %.15le\n",0.0,0.0);
  fprintf(bary_output,"%.15le  %.15le\n",0.0,spline_value(ARRAY_SIZE,bary_arr,pres_arr,d_bary,p));
  fprintf(smod_output,"%.15le  %.15le\n",0.0,spline_value(linesmu,smod_arr_2,pres_arr_2,d_smod_2,p));
  fprintf(gam1_output,"%.15le  %.15le\n",0.0,spline_value(ARRAY_SIZE,gamma1_arr,pres_arr,d_gamma1,p));

  //loop over RK4
  do {
    variables[0]=p;
    variables[1]=m;
    RK_var_step(variables,2,&r,dr,TOL,&drused,&dr,derivatives);
    p=variables[0];
    m=variables[1];
      
    //output some data every step
    if (p>=pres_arr[ARRAY_SIZE-lines_core]) {
      rho=spline_value(lines_core,core_rho,pres_core,d_rho_core,p);
      bary=spline_value(lines_core,core_bary,pres_core,d_bary_core,p);
    } else if (p<=pres_arr[ARRAY_SIZE-lines_core-1]) {
      rho=spline_value(lines_crust,crust_rho,pres_crust,d_rho_crust,p);
      bary=spline_value(lines_crust,crust_bary,pres_crust,d_bary_crust,p);
    } else if ((p>pres_arr[ARRAY_SIZE-lines_core-1]) && (p<pres_arr[ARRAY_SIZE-lines_core])) {
      rho=linear_int(ARRAY_SIZE,p,rho_arr,pres_arr);
      bary=linear_int(ARRAY_SIZE,p,bary_arr,pres_arr);
    }
    mu=spline_value(linesmu,smod_arr_2,pres_arr_2,d_smod_2,p);
    gamma1=spline_value(ARRAY_SIZE,gamma1_arr,pres_arr,d_gamma1,p);
  
    fprintf(mass_output,"%.15le  %.15le\n",r,m);
    fprintf(pres_output,"%.15le  %.15le\n",r,p);
    fprintf(dens_output,"%.15le  %.15le\n",r,rho);
    fprintf(grav_output,"%.15le  %.15le\n",r,(G*m)/(r*r));
    fprintf(dpdr_output,"%.15le  %.15le\n",r,(p-p_old)/drused);
    fprintf(bary_output,"%.15le  %.15le\n",r,bary);
    fprintf(smod_output,"%.15le  %.15le\n",r,mu);
    fprintf(gam1_output,"%.15le  %.15le\n",r,gamma1);
    //printf("r=%lf,  TOV GR term 1 = %lf,  TOV GR term 2 =   %lf,  TOV GR term 3 =   %lf\n",r,1.0+p/(rho*C*C),1.0+(4.0*PI*r*r*r*p)/(m*C*C),1.0/(1.0-(2.0*G*m)/(r*C*C)));
    if((m<0) || (p<0) || (rho<0) || ((G*m)/(r*r)<0) || (bary<0) || (mu<0) || (gamma1<0)){printf("uh-oh, negative value at r=%le    %le %le %le %le %le %le\n",r,m,p,rho,bary,mu,gamma1);}
    p_old=p;
  
  } while (p > pres_arr[1]);  //leave loop if P is below the boundary (what is the best boundary?)
  
  rho=spline_value(lines_core,core_rho,pres_core,d_rho_core,p_c);  //don't worry about crust vs core spline, this is definately in the core
  fprintf(dens_output,"%.15le  %.15le\n",r,rho);  //print out extra line just to make sure total radius is written down somewhere
  fprintf(bary_output,"%.15le  %.15le\n",r,nbcc);  //print out extra line to write down the crust-core transition density is written down somewhere

//printf("%lf  %lf  %lf  %lf",atof(argv[1]),(double)((int)atof(argv[2])),0.0,m/(1.989E+33));    //for TOV_repreat
  printf("J=%lf  L=%lf  Ksym=%lf  M*=%lf\n",atof(argv[1]),atof(argv[2]),atof(argv[3]),m/(1.989E+33));
  printf("final neutron star of mass %lf M_sun and radius %lf Km, calcualted using central pressure %.6e and central density %.6e.\n",m/(1.989E+33),r/(1.0E+05),p_c,rho);

  fclose(mass_output);
  fclose(pres_output);
  fclose(dens_output);
  fclose(grav_output);
  fclose(dpdr_output);
  fclose(bary_output);
/*----------------solve the TOV equations for 1.4 solar mass star, while printing out data----------------*/
}





//variables:  0=p,1=m
void derivatives(double rpos,double var_in[],double d_var_out[]){
  double rho_val;
  if (var_in[0]>=pres_arr[ARRAY_SIZE-lines_core]) {
    rho_val=spline_value(lines_core,core_rho,pres_core,d_rho_core,var_in[0]);
  } else if (var_in[0]<=pres_arr[ARRAY_SIZE-lines_core-1]) {
    rho_val=spline_value(lines_crust,crust_rho,pres_crust,d_rho_crust,var_in[0]);
  } else if ((var_in[0]>pres_arr[ARRAY_SIZE-lines_core-1]) && (var_in[0]<pres_arr[ARRAY_SIZE-lines_core])) {
    rho_val=linear_int(ARRAY_SIZE,var_in[0],rho_arr,pres_arr);
  }

  if (var_in[1] == 0.0 && rpos == 0.0) {
    d_var_out[0]=0.0;
  } else if (var_in[1] == 0.0 && rpos != 0.0) {
    d_var_out[0]=-G*rho_val*(1.0/(rpos*rpos))*(1.0+var_in[0]/(rho_val*pow(C,2)))*((4.0*PI*rpos*rpos*rpos*var_in[0])/pow(C,2));
  } else if (var_in[1] != 0.0) {
    d_var_out[0]=-(G*rho_val*var_in[1])/(rpos*rpos);
    d_var_out[0]=d_var_out[0]*(1.0+var_in[0]/(rho_val*C*C));
    d_var_out[0]=d_var_out[0]*(1.0+((4.0*PI*rpos*rpos*rpos*var_in[0])/(var_in[1]*C*C)));
    d_var_out[0]=d_var_out[0]*(1.0/(1.0-(2.0*G*var_in[1])/(rpos*C*C)));
  }
  d_var_out[1]=4.0*PI*rpos*rpos*rho_val;
}

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
