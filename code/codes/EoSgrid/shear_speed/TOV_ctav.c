//solves TOV for a 1.4M_sun neutron star. Modified to use newton's EoSs as inputs (9/6/20). Modified to require input L and K values to match those of the input EoS (2/11/20)
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
double *rho_arr,*pres_arr,*d_rho,*bary_arr,*d_bary,*smod_arr,*d_smod,*gamma1_arr,*d_gamma1;
double *core_rho,*crust_rho,*core_bary,*crust_bary,*d_rho_crust,*d_rho_core,*d_bary_crust,*d_bary_core,*pres_core,*pres_crust;

void main (int argc, char *argv[]) {
  double p_c,dp_c,p,dr,dp,r,e_dens,nb,rho,m,k1m,k2m,k1p,k2p,k3m,k4m,k3p,k4p,p_old,bary,mu,nbcc,M_target,gamma1;
  int i,lines,underscore,j,l,k;
  char c;
  char eosname[80]="", testchar[80]="", testfile[80]="";
  char J[8]="", L[8]="", K[8]="";
  double variables[2],drused;
  int *Z=malloc(0*sizeof(int)),*A=malloc(0*sizeof(int));
  int Kvalue;

 // if (argc != 3) {printf("Input: ./star [file name] [M* (solar masses)]\n");exit(0);}
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
  strcat(eosname,"../EoSs/eos_LQDROP_SkXi450_");
  strcat(eosname,argv[1]);
  strcat(eosname,".0_");
  strcat(eosname,argv[2]);
  strcat(eosname,".00_");
  strcat(eosname,argv[3]);
  strcat(eosname,".00_2.2_Imid_.dat");
//printf("%s\n",eosname);
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
//if (p>pres_arr[i+1] && i!=ARRAY_SIZE-1) {printf("WARNING: -ve dp/drho occurs in %s\n",argv[1]);}
    rho_arr[i]=e_dens;  //  column 1 of file = energy density/(C*C)   (=rho at NR limit)
    pres_arr[i]=p;       //  column 2 of file = pressure
    bary_arr[i]=nb;       //  column 3 of file = baryon density
    smod_arr[i]=mu;
    gamma1_arr[i]=gamma1;
  }
  fclose(eos);
/*-----------------------------setup output files and read the input eos file-----------------------------*/

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
  for (i=2;i<ARRAY_SIZE;i++){
    if ((smod_arr[i]==smod_arr[i-1]) && (smod_arr[i]==smod_arr[i-2]) && (smod_arr[i]!=0)){
      smod_arr[i-2]=(smod_arr[i-3]+2.0*smod_arr[i-2])/3.0;
      smod_arr[i]=(smod_arr[i+1]+2.0*smod_arr[i])/3.0;
    }  else if (smod_arr[i]==smod_arr[i-1]  && (smod_arr[i]!=smod_arr[i-2]) && (smod_arr[i]!=0)) {
      smod_arr[i]=(smod_arr[i+1]+2.0*smod_arr[i])/3.0;
      smod_arr[i-1]=(smod_arr[i-2]+2.0*smod_arr[i-1])/3.0;
    }
    if (smod_arr[i]==0){  //get the crust-core transition baryon density
      nbcc=bary_arr[i-1];
   //   printf("n_bcc=%le\n",bary_arr[i-1]);
      break;
    }
  }
//  for (i=1;i<ARRAY_SIZE;i++){if ((smod_arr[i]==smod_arr[i-1]) && (smod_arr[i]!=0)){printf("mu   %d\n",i);}}
  for (i=1;i<ARRAY_SIZE;i++){if ((bary_arr[i]==bary_arr[i-1])){printf("nb   %d\n",i);}}
  for (i=1;i<ARRAY_SIZE;i++){if ((pres_arr[i]==pres_arr[i-1])){printf("p    %d\n",i);}}
  for (i=1;i<ARRAY_SIZE;i++){if ((rho_arr[i]==rho_arr[i-1])){printf("rho  %d\n",i);}}

//  for (i=1;i<ARRAY_SIZE;i++){gamma1_arr[i]=gamma1_arr[i]+i/(1e+15);}
//  for (i=1;i<ARRAY_SIZE;i++){if ((gamma1_arr[i]==gamma1_arr[i-1])){printf("gamma1  %d\n",i);}}

  //set up cubic splines
  cubic_spline(ARRAY_SIZE,smod_arr,pres_arr,d_smod);
  for (i=0;i<ARRAY_SIZE;i++){d_smod[i]=0;}
  cubic_spline(ARRAY_SIZE,rho_arr,pres_arr,d_rho);
  cubic_spline(ARRAY_SIZE,bary_arr,pres_arr,d_bary);
  cubic_spline(ARRAY_SIZE,gamma1_arr,pres_arr,d_gamma1);
  for (i=0;i<ARRAY_SIZE;i++){d_gamma1[i]=0;}
/*-------------modify eos data to fix repeated lines, which cause errors in the cubic splines-------------*/

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


//if (Kvalue<=-250){
for (i=0;i<ARRAY_SIZE;i++){d_rho_core[i]=0;}  //use linear splines in core
for (i=0;i<ARRAY_SIZE;i++){d_bary_core[i]=0;}
//}
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
      mu=spline_value(ARRAY_SIZE,smod_arr,pres_arr,d_smod,p);
      gamma1=spline_value(ARRAY_SIZE,gamma1_arr,pres_arr,d_gamma1,p);
      if((m<0) || (p<0) || (rho<0) || ((G*m)/(r*r)<0) || (bary<0) || (mu<0) || (gamma1<0)){printf("uh-oh, negative value at r=%le    %le %le %le %le %le %le\n",r,m,p,rho,bary,mu,gamma1);}
    } while (p > pres_arr[1]);  //leave loop if P is below the boundary (what is the best boundary?)
  
    rho=spline_value(lines_core,core_rho,pres_core,d_rho_core,p_c);  //don't worry about crust vs core spline, this is definately in the core
//    printf("neutron star of mass %lf M_sun and radius %lf Km, calcualted using central pressure %.6e and central density %.6e.\n",m/(1.989E+33),r/(1.0E+05),p_c,rho);
  
    if (m>(M_target*(1.989E+33))){p_c=p_c-dp_c;dp_c=dp_c*0.1;}
    p_c=p_c+dp_c;
  }while (dp_c>1e+29);
/*-----------loop over pressure values in the TOV until to get a total mass of 1.4 solar masses-----------*/
/*
strcat(testfile,"EoS_spline_");
strcat(testfile,J);
strcat(testfile,"_");
strcat(testfile,L);
strcat(testfile,"_");
strcat(testfile,K);
strcat(testfile,"_.txt");
  p=pres_arr[0];
  dp=pres_arr[0]/10.0;
  FILE *new=fopen(testfile,"w");
  do{
    if (p>=pres_arr[ARRAY_SIZE-lines_core]) {
      rho=spline_value(lines_core,core_rho,pres_core,d_rho_core,p);
    } else if (p<=pres_arr[ARRAY_SIZE-lines_core-1]) {
      rho=spline_value(lines_crust,crust_rho,pres_crust,d_rho_crust,p);
    } else if ((p>pres_arr[ARRAY_SIZE-lines_core-1]) && (p<pres_arr[ARRAY_SIZE-lines_core])) {
      rho=linear_int(ARRAY_SIZE,p,rho_arr,pres_arr);
    }
    fprintf(new,"%le  %le\n",rho,p);

    if (p>=dp*10000.0){dp=dp*10.0;}
    p=p+dp;
  } while (p<=pres_arr[ARRAY_SIZE-1]);
  exit(0);
*/
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
  fprintf(smod_output,"%.15le  %.15le\n",0.0,spline_value(ARRAY_SIZE,smod_arr,pres_arr,d_smod,p));
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
    mu=spline_value(ARRAY_SIZE,smod_arr,pres_arr,d_smod,p);
    gamma1=spline_value(ARRAY_SIZE,gamma1_arr,pres_arr,d_gamma1,p);
  
    fprintf(mass_output,"%.15le  %.15le\n",r,m);
    fprintf(pres_output,"%.15le  %.15le\n",r,p);
    fprintf(dens_output,"%.15le  %.15le\n",r,rho);
    fprintf(grav_output,"%.15le  %.15le\n",r,(G*m)/(r*r));
    fprintf(dpdr_output,"%.15le  %.15le\n",r,(p-p_old)/drused);
    fprintf(bary_output,"%.15le  %.15le\n",r,bary);
    fprintf(smod_output,"%.15le  %.15le\n",r,mu);
    fprintf(gam1_output,"%.15le  %.15le\n",r,gamma1);
    if((m<0) || (p<0) || (rho<0) || ((G*m)/(r*r)<0) || (bary<0) || (mu<0) || (gamma1<0)){printf("uh-oh, negative value at r=%le    %le %le %le %le %le %le\n",r,m,p,rho,bary,mu,gamma1);}
    p_old=p;
  
  } while (p > pres_arr[1]);  //leave loop if P is below the boundary (what is the best boundary?)
  
  rho=spline_value(lines_core,core_rho,pres_core,d_rho_core,p_c);  //don't worry about crust vs core spline, this is definately in the core
  fprintf(dens_output,"%.15le  %.15le\n",r,rho);  //print out extra line just to make sure total radius is written down somewhere
  fprintf(bary_output,"%.15le  %.15le\n",r,nbcc);  //print out extra line to write down the crust-core transition density is written down somewhere




 // printf("\n");
//  printf("%s  %s  %s  %lf",J,L,K,m/(1.989E+33));



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
