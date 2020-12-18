//code uses fluid core and solid crust. Shoots from r=dr (core) to r=Rcc, and from r=R*-dr (crust) to r=Rcc. Ocean assumed to be negligable. Modified to use newton's EoSs as inputs (9/6/20)
//TODO: add A so that g-modes can be found?
//TODO: do cubic splines of pressure, density, etc..., then interpolate to find fine mesh of points for g1, Utilde,etc...
//made to be used with star.sh, so used several times in a row
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "nr_codes.h"
#include "eqbm.h"
#include "shoot.h"
#define PI 3.14159265359
#define C 2.99792458E+10
#define HBAR 1.0545726663E-27
#define MN 1.674928610E-24
#define G 6.6725985E-8
#define TOL 1.0E-7
#define E_CHARGE 4.8032042510E-10
#define KBOLTZ 1.380649E-16

double *r_pos,*Vtilde,*x,*c1,*Utilde,*Ar,*g1,*d_Vtilde,*d_rpos,*d_c1,*d_Utilde,*d_Ar,*d_g1;
double *pressure,*d_pressure,radius,*gravity,*d_gravity,*density,*d_density,*mass_r,*d_mass_r,*bary_dens,*d_bary_dens,*smod_arr,*d_smod;
double *Z,*A,*nb_ZA,*ion_ratio,*d_Z,*d_A,*d_ion_ratio;
double *EoS_bary_dens,*pres_arr,*rho_arr;
int rsteps,lines,eoslines;
double boundarypos,boundaryposx;
double *logp,*logrho,*d_logp;
double NB_CC;

void main(int argc, char *argv[])
{
  double *dp_dr=malloc(0*sizeof(double));
  double temp,y1,y2,boundary1,boundary2,boundary3,omega,d_omega,sigma,z1,z2,z3,z4,scale;
  char c;
  char eosname[80]="",outfile1[80]="",outfile2[80]="", testchar[80]="";
  int i,l,j,k,underscore,trials,option=0;
  double a,temperature,gamma,proton,nucleon,n_i,mu;
  int trialsz3=0,test=0;
  double z3test,z3testold,z3testoldold,d_z3,z3init,omega_z4,z3old4,z3old3,z3old2,z3old1,omegainit,z3initinit,d_z3init,z3init_z4;
  double omegalow,omegahigh,z3low,z3high;
  char J[8]="", L[8]="", K[8]="";
  int Kvalue;

  if (argc != 7) {printf("Input: ./star [file name] option omega d_omega z3 d_z3\n");exit(0);}
/*--------------get properties as function of radius using TOV--------------*/
  //read star's data from the files output by the TOV
  FILE * NS_dens = fopen("NS_dens.txt", "r");
  FILE * NS_pres = fopen("NS_pres.txt", "r");
  FILE * NS_dpdr = fopen("NS_dpdr.txt", "r");
  FILE * NS_grav = fopen("NS_grav.txt", "r");
  FILE * NS_mass = fopen("NS_mass.txt", "r");
  FILE * NS_bary = fopen("NS_bary.txt", "r");
  FILE * NS_smod = fopen("NS_smod.txt", "r");
  FILE * NS_gam1 = fopen("NS_gam1.txt", "r");

  //get eos file from input J,L values
  strcat(eosname,argv[1]);
  FILE * eos = fopen(eosname, "r");


  //get J,L and K from file name
  strcat(testchar,argv[1]);
  i=0;j=0;l=0;k=0;
  underscore=0;
  while (testchar[i] != '\0') {       // Stop looping when we reach the null-character
    if (testchar[i]=='_') {
      underscore++;
    }
    if (testchar[i]!='_' && underscore==3) {
      J[j]=testchar[i];
      j++;
    }
    if (testchar[i]!='_' && underscore==4) {
      L[l]=testchar[i];
      l++;
    }
    if (testchar[i]!='_' && underscore==5) {
      K[k]=testchar[i];
      k++;
    }
    i++;
  }
  Kvalue=atoi(K);
  if (Kvalue>-250){exit(0);}



  //get length of files
  rsteps=0;
  c = getc(NS_mass);
  while (c != EOF){
    if (c == '\n'){
      rsteps=rsteps+1;
    }
    c = getc(NS_mass);
  }
  fseek(NS_mass,0,SEEK_SET);

  //allocate memory to arrays
  density=realloc(density,rsteps*sizeof(double));
  pressure=realloc(pressure,rsteps*sizeof(double));
  dp_dr=realloc(dp_dr,rsteps*sizeof(double));
  gravity=realloc(gravity,rsteps*sizeof(double));
  mass_r=realloc(mass_r,rsteps*sizeof(double));
  r_pos=realloc(r_pos,rsteps*sizeof(double));
  bary_dens=realloc(bary_dens,rsteps*sizeof(double));
  smod_arr=realloc(smod_arr,rsteps*sizeof(double));
  g1=realloc(g1,rsteps*sizeof(double));

  //read data in
  for(i=0;i<rsteps;i++){
    fscanf(NS_dens,"%le  %le",&r_pos[i],&density[i]);
    fscanf(NS_pres,"%le  %le",&r_pos[i],&pressure[i]);
    fscanf(NS_dpdr,"%le  %le",&r_pos[i],&dp_dr[i]);
    fscanf(NS_grav,"%le  %le",&r_pos[i],&gravity[i]);
    fscanf(NS_mass,"%le  %le",&r_pos[i],&mass_r[i]);
    fscanf(NS_bary,"%le  %le",&r_pos[i],&bary_dens[i]);
    fscanf(NS_smod,"%le  %le",&r_pos[i],&smod_arr[i]);
    fscanf(NS_gam1,"%le  %le",&r_pos[i],&g1[i]);
  }
  //there is an extra line in the density file at r=r_max so that no matter how often the TOV outputs data to file the total radius is given
  fscanf(NS_dens,"%le  %le",&radius,&temp);
  fscanf(NS_bary,"%le  %le",&radius,&NB_CC);
 // printf("baryon density at crust-core boundary = %le\n",NB_CC);

  for (i=0;i<rsteps;i++){
    gravity[i]=fabs(dp_dr[i]/density[i]);
  }

  fclose(NS_dens);
  fclose(NS_pres);
  fclose(NS_dpdr);
  fclose(NS_grav);
  fclose(NS_mass);
  fclose(NS_bary);
  fclose(NS_smod);
  fclose(NS_gam1);
/*--------------get properties as function of radius using TOV--------------*/

/*------------------------------read eos file-------------------------------*/
  FILE *EoSfile=fopen(eosname, "r");

  //get length of files
  eoslines=0;
  c = getc(EoSfile);
  while (c != EOF){
    if (c == '\n'){
      eoslines=eoslines+1;
    }
    c = getc(EoSfile);
  }
  fseek(EoSfile,0,SEEK_SET);

  EoS_bary_dens=realloc(EoS_bary_dens,eoslines*sizeof(double));
  pres_arr=realloc(pres_arr,eoslines*sizeof(double));
  rho_arr=realloc(rho_arr,eoslines*sizeof(double));

  for(i=0;i<eoslines;i++){
    fscanf(EoSfile,"%le  %le  %le  %le  %le",&rho_arr[eoslines-1-i],&pres_arr[eoslines-1-i],&EoS_bary_dens[eoslines-1-i],&gamma,&mu);
  }
  fclose(EoSfile);
/*------------------------------read eos file-------------------------------*/

/*

  for (i=2;i<eoslines;i++){
    if ((rho_arr[i]==rho_arr[i-1]) && (rho_arr[i]==rho_arr[i-2])){
      rho_arr[i-2]=(rho_arr[i-3]+2.0*rho_arr[i-2])/3.0;
      rho_arr[i]=(rho_arr[i+1]+2.0*rho_arr[i])/3.0;
    }  else if (rho_arr[i]==rho_arr[i-1]  && (rho_arr[i]!=rho_arr[i-2])) {
      rho_arr[i]=(rho_arr[i+1]+2.0*rho_arr[i])/3.0;
      rho_arr[i-1]=(rho_arr[i-2]+2.0*rho_arr[i-1])/3.0;
    }
  }
  for (i=2;i<eoslines;i++){
    if ((pres_arr[i]==pres_arr[i-1]) && (pres_arr[i]==pres_arr[i-2])){
      pres_arr[i-2]=(pres_arr[i-3]+2.0*pres_arr[i-2])/3.0;
      pres_arr[i]=(pres_arr[i+1]+2.0*pres_arr[i])/3.0;
    }  else if (pres_arr[i]==pres_arr[i-1]  && (pres_arr[i]!=pres_arr[i-2])) {
      pres_arr[i]=(pres_arr[i+1]+2.0*pres_arr[i])/3.0;
      pres_arr[i-1]=(pres_arr[i-2]+2.0*pres_arr[i-1])/3.0;
    }
  }
  for (i=2;i<eoslines;i++){
    if ((EoS_bary_dens[i]==EoS_bary_dens[i-1]) && (EoS_bary_dens[i]==EoS_bary_dens[i-2])){
      EoS_bary_dens[i-2]=(EoS_bary_dens[i-3]+2.0*EoS_bary_dens[i-2])/3.0;
      EoS_bary_dens[i]=(EoS_bary_dens[i+1]+2.0*EoS_bary_dens[i])/3.0;
    }  else if (EoS_bary_dens[i]==EoS_bary_dens[i-1]  && (EoS_bary_dens[i]!=EoS_bary_dens[i-2])) {
      EoS_bary_dens[i]=(EoS_bary_dens[i+1]+2.0*EoS_bary_dens[i])/3.0;
      EoS_bary_dens[i-1]=(EoS_bary_dens[i-2]+2.0*EoS_bary_dens[i-1])/3.0;
    }
  }
  for (i=1;i<eoslines;i++){if ((EoS_bary_dens[i]==EoS_bary_dens[i-1])){printf("%d\n",i);}}
  for (i=1;i<eoslines;i++){if ((pres_arr[i]==pres_arr[i-1])){printf("%d\n",i);}}
  for (i=1;i<eoslines;i++){if ((rho_arr[i]==rho_arr[i-1])){printf("%d\n",i);}}

*/

/*----------get dimensionless properties of star and setup splines----------*/
  x=realloc(x,rsteps*sizeof(double));
  d_rpos=realloc(d_rpos,rsteps*sizeof(double));
  d_density=realloc(d_density,rsteps*sizeof(double));
  d_gravity=realloc(d_gravity,rsteps*sizeof(double));
  d_pressure=realloc(d_pressure,rsteps*sizeof(double));
  d_mass_r=realloc(d_mass_r,rsteps*sizeof(double));
  d_bary_dens=realloc(d_bary_dens,rsteps*sizeof(double));
  d_smod=realloc(d_smod,rsteps*sizeof(double));
  d_g1=realloc(d_g1,rsteps*sizeof(double));

  //setup_g1(g1,eoslines,rho_arr,pres_arr,EoS_bary_dens);       //remove comment on this line to use numerical derivative to get g1

  for (i=0;i<rsteps;i++){
    x[i]=log(r_pos[i]/pressure[i]);
  }
  x[0]=x[1]+2.0*(log((0.5*(r_pos[1]-r_pos[0]))/((pressure[1]+pressure[0])/2.0))-x[1]);

  cubic_spline(rsteps,r_pos,x,d_rpos);
  cubic_spline(rsteps,density,r_pos,d_density);
  cubic_spline(rsteps,gravity,r_pos,d_gravity);
  cubic_spline(rsteps,pressure,r_pos,d_pressure);
  cubic_spline(rsteps,mass_r,r_pos,d_mass_r);
  cubic_spline(rsteps,bary_dens,r_pos,d_bary_dens);
//  cubic_spline(rsteps,smod_arr,r_pos,d_smod);
  for (i=0;i<rsteps;i++){d_smod[i]=0;}
//  cubic_spline(rsteps,g1,r_pos,d_g1);
  for (i=0;i<rsteps;i++){d_g1[i]=0;}
/*----------get dimensionless properties of star and setup splines----------*/

/*--------------------use shooting method to find sigma---------------------*/
  l=2;
  boundary2=1e+30;
  boundary3=1e+40;
  z3testold=1e+30;
  z3testoldold=1e+40;

  //find Rcc using a pressure value obtained from the EoS file
  boundarypos=r_pos[0];
  do{
    if ((spline_value(rsteps,bary_dens,r_pos,d_bary_dens,boundarypos)>NB_CC) && (spline_value(rsteps,bary_dens,r_pos,d_bary_dens,boundarypos+1)<NB_CC)){
      break;
    }
    boundarypos=boundarypos+1.0;
  } while (boundarypos<radius);
  boundaryposx=log(boundarypos/spline_value(rsteps,pressure,r_pos,d_pressure,boundarypos));
//  printf("n_bcc=%le,   Rcc = %lf\n",NB_CC,boundarypos);


 // printf("trace the minima (3), use this omega and z3 (2), scan over omega and z3 values (1) or find a solution (0)?\n");
 // scanf("%d",&option);
  option=atof(argv[4]);

  //option 3 is to find the minima in the error in z2/z1 and z4
  if (option == 3){

    //----------------------------------------------------------------------------------------------------------------//
    double z4min,z2z1min;int finished=0;
    int temp;
    FILE * z2_z1file = fopen("z2z1_min.txt", "w");
    FILE * z4file = fopen("z4_min.txt", "w");

    strcat(outfile1,"minima/minima_z2z1_");
    strcat(outfile1,argv[1]);
    strcat(outfile1,"_");
    strcat(outfile1,argv[2]);
    strcat(outfile1,"_");
    strcat(outfile1,argv[3]);
    strcat(outfile1,".txt");
    FILE * out1 = fopen(outfile1, "w");

    strcat(outfile2,"minima/minima_z4_");
    strcat(outfile2,argv[1]);
    strcat(outfile2,"_");
    strcat(outfile2,argv[2]);
    strcat(outfile2,"_");
    strcat(outfile2,argv[3]);
    strcat(outfile2,".txt");
    FILE * out2 = fopen(outfile2, "w");

//    printf("Input lowest omega value: ");
//    scanf("%lf",&omegalow);
//    printf("Input highest omega value: ");
//    scanf("%lf",&omegahigh);
//    printf("Input omega step: ");
//    scanf("%lf",&d_omega);
//    printf("Input lowest z3 value: ");
//    scanf("%lf",&z3low);
//    printf("Input highest z3 value: ");
//    scanf("%lf",&z3high);
//    printf("Input z3 step: ");
//    scanf("%lf",&d_z3);
    omegalow=atof(argv[5]);
    omegahigh=atof(argv[6]);
    d_omega=atof(argv[7]);
    z3low=atof(argv[8]);
    z3high=atof(argv[9]);
    d_z3=atof(argv[10]);

//    printf("reset to initial z3 value after each omega value, or just go back a few steps (faster, but only works for small omega ranges) (reset=0/go_back=1): ");
//    scanf("%d",&temp);
    temp=1;

    omega=omegalow;
    if (temp){
      z3init=z3low+2*d_z3;
      z2z1min=z3init;
      z4min=z3init;
    }
    //scan over omega values, getting the z3 values which minimise the BCs
    do{
      if (temp){
        if (z3init <= z3high-0.5*d_z3){
          if (z2z1min<z4min){
            z3init=z2z1min-2*d_z3;
          } else {z3init=z4min-2*d_z3;}
        } else {z3init=z3low;}
      } else {
        z3init=z3low;
      }
      boundary1=1e+30;
      boundary2=1e+30;

      //increment z3 value until both minima have been passed
      do {
      z3=z3init;
        sigma=sqrt(omega*omega*G*mass_r[rsteps-1]/(radius*radius*radius));
        core(&z1,&y1,&y2,omega,l,0,sigma);
        crust(&z1,&z2,&z3,&z4,omega,l,sigma,0,&scale);

        if (boundary1>fabs(z2/z1-Vtildevalue(boundarypos)*(1.0-y2/y1))) {
          boundary1=fabs(z2/z1-Vtildevalue(boundarypos)*(1.0-y2/y1));      //boundary1 is McDermott eqn 20a/20b
          z2z1min=z3init;
          finished=0;
        }
        if (boundary2>fabs(z4)) {
          boundary2=fabs(z4);                                              //boundary2 is McDermott eqn 20c
          z4min=z3init;
          finished=0;
        }  else {finished=finished+1;} if (finished==5) {finished=0;}//{break;}  //stop early if past the minima in both BCs

//        printf("%le  %.12le  %.10le  %.10le\n",z3init,omega,boundary1,boundary2);
        z3init=z3init+d_z3;
      } while(z3init<=z3high);

      fprintf(out1,"%.12le  %.12le\n",omega,z2z1min);
      fprintf(out2,"%.12le  %.12le\n",omega,z4min);
      omega=omega+d_omega;
      printf("\n");
    } while(omega<=omegahigh);
    //----------------------------------------------------------------------------------------------------------------//

  //option 1 is to scan over omega and z3 values
  } else if (option == 1){
    //----------------------------------------------------------------------------------------------------------------//
    FILE * z2_z1file = fopen("z2z1_scan.txt", "w");
    FILE * z4file = fopen("z4_scan.txt", "w");

    strcat(outfile1,"scan/scan_z2z1_");
    strcat(outfile1,argv[1]);
    strcat(outfile1,"_");
    strcat(outfile1,argv[2]);
    strcat(outfile1,"_");
    strcat(outfile1,argv[3]);
    strcat(outfile1,".txt");
    FILE * out1 = fopen(outfile1, "w");

    strcat(outfile2,"scan/scan_z4_");
    strcat(outfile2,argv[1]);
    strcat(outfile2,"_");
    strcat(outfile2,argv[2]);
    strcat(outfile2,"_");
    strcat(outfile2,argv[3]);
    strcat(outfile2,".txt");
    FILE * out2 = fopen(outfile2, "w");

//    printf("Input lowest omega value: ");
//    scanf("%lf",&omegalow);
//    printf("Input highest omega value: ");
//    scanf("%lf",&omegahigh);
//    printf("Input omega step: ");
//    scanf("%lf",&d_omega);
//    printf("Input lowest z3 value: ");
//    scanf("%lf",&z3low);
//    printf("Input highest z3 value: ");
//    scanf("%lf",&z3high);
//    printf("Input z3 step: ");
//    scanf("%lf",&d_z3);
    omegalow=atof(argv[5]);
    omegahigh=atof(argv[6]);
    d_omega=atof(argv[7]);
    z3low=atof(argv[8]);
    z3high=atof(argv[9]);
    d_z3=atof(argv[10]);

    z3init=z3low;
    //loop over z3 values
    do{
      omega=omegalow;
      //loop over the omega range for this z3 value
      do {
        z3=z3init;
        sigma=sqrt(omega*omega*G*mass_r[rsteps-1]/(radius*radius*radius));
        core(&z1,&y1,&y2,omega,l,0,sigma);
        crust(&z1,&z2,&z3,&z4,omega,l,sigma,0,&scale);
        boundary1=fabs(z2/z1-Vtildevalue(boundarypos)*(1.0-y2/y1));      //boundary1 is McDermott eqn 20a/20b
        boundary2=fabs(z4);                                              //boundary2 is McDermott eqn 20c

        fprintf(out1,"%.12le  %.12le  %.10le\n",omega,z3init,boundary1);
        fprintf(out2,"%.12le  %.12le  %.10le\n",omega,z3init,boundary2);
        printf("%le  %.12le  %.10le  %.10le\n",z3init,omega,boundary1,boundary2);
        omega=omega+d_omega;
      } while(omega<=omegahigh);

      z3init=z3init+d_z3;
      printf("\n");
      fprintf(out1,"\n");
      fprintf(out2,"\n");
   } while(z3init<=z3high);
    //----------------------------------------------------------------------------------------------------------------//

  //option 0 is to find the z3(R*) and omega eigenvalues
  } else if (option == 0) {
    //----------------------------------------------------------------------------------------------------------------//
    int loops;
 //   printf("Input initial omega value: ");
 //   scanf("%lf",&omega);
  //  printf("Input omega step: ");
 //   scanf("%lf",&d_omega);
  //  printf("Input initial z3 value: ");
 //   scanf("%lf",&z3initinit);
  //  printf("Input z3 step: ");
 //   scanf("%lf",&d_z3init);
    omega=atof(argv[3]);
    d_omega=atof(argv[4]);
    z3initinit=atof(argv[5]);
    d_z3init=atof(argv[6]);

    do{
      
      //the first part finds the z3(R*) value that minimises z2/z1-z2(y2)/z1(y1) at Rcc, at the current omega value
      //z3=z3(r) used in this crust, z3init=z3(R*) used while searching for the z3(R*) eigenvalue, z3initinit=z3(R*) is the initial z3init=z3(R*) guess
      z3init=z3initinit;
      d_z3=d_z3init;
      boundary2=1e+30;
      loops=0;
      do {  
        if (loops == 100 && fabs(d_z3) > fabs(d_z3init)/20.0) {  //if after 100 loops the step size hasn't been reduced at least twice, try diffrent starting z3
          z3init=-z3init;                            //try negative of previous z3
          d_z3=d_z3init;                             //basically, this is to try to avoid being on the wrong side of the maxima in the omega-z3 plane
          boundary2=1e+30;                           //warning: will also trigger if z3initinit >= 95*d_z3init away from the minima
          //printf("try starting from -z3init\n");
        }
        z3=z3init;
        sigma=sqrt(omega*omega*G*mass_r[rsteps-1]/(radius*radius*radius));
        core(&z1,&y1,&y2,omega,l,0,sigma);
        crust(&z1,&z2,&z3,&z4,omega,l,sigma,0,&scale);
        boundary1=fabs(z2/z1-Vtildevalue(boundarypos)*(1.0-y2/y1));             //test the boundary condition at r=R*, McDermott eqn 20a/20b
        bound_test(boundary1,&boundary2,&boundary3,&z3init,&d_z3,0);            //simple test to determine the next trial z3(R*) value
        loops++;
      } while (boundary1 > 1e-10 && fabs(d_z3)>fabs(z3init)*(1e-15));           //leave loop when the boundary is satisfied or z3 is at a high precision
      z3init=z3init-d_z3;                                                       //undo the last z3 change, since it was never used

      //set initial trial z3(R*) value (z3initinit) to be closer to the solution to speed up future calculations (assuming initial omega value is fairly close to a solution)
      if (test==0){
        z3initinit=round(z3init/d_z3init)*d_z3init;
  //      printf("changing initial z3 to %lf due to finding z3 at %lf\n",z3initinit,z3init);
        test=1;
      }


      //the second part finds the z3(R*) value that minimises |z4| at Rcc, at the current omega value
      z3init_z4=z3initinit;
      d_z3=d_z3init;
      boundary2=1e+30;
      do {
        z3=z3init_z4;
        sigma=sqrt(omega*omega*G*mass_r[rsteps-1]/(radius*radius*radius));
        core(&z1,&y1,&y2,omega,l,0,sigma);
        crust(&z1,&z2,&z3,&z4,omega,l,sigma,0,&scale);
        boundary1=fabs(z4);                                                     //test the boundary condition at r=R*, McDermott eqn 20c
        bound_test(boundary1,&boundary2,&boundary3,&z3init_z4,&d_z3,0);         //simple test to determine the next trial z3(R*) value
      } while (boundary1 > 1e-10 && fabs(d_z3)>fabs(z3init_z4)*(1e-15));        //leave loop when the boundary is staisfied or z3 is at a high precision
      z3init_z4=z3init_z4-d_z3;                                                 //undo the last z3 change, since it was never used


      //the final part is to test the omega value by finding the difference between the z3(R*) values that minimise the boundary conditions at Rcc
      z3test=fabs(z3init-z3init_z4);
      //printf("z3 trial=%d,  omega=%.12le, l=%d, z3 difference=%le, z3(z2/z1)=%le\n",trialsz3+1,omega,l,z3test,z3init);
      trialsz3=trialsz3+1;
      if (trialsz3 > 60) {
  //      printf("exiting loop, z3 not found\n");
        break;
      }
      bound_test(z3test,&z3testold,&z3testoldold,&omega,&d_omega,1);            //simple test to determine the next trial omega value

    } while (z3test > fabs(z3init)/1e+4);  //quit once the z3(R*) eigenvalue has been found to 5 sig. fig.


    //output z values to files using the omega and z3(R*) eigenvalues
    omega=omega-d_omega;
    z3init_z4=z3init_z4-d_z3;
    sigma=sqrt(omega*omega*G*mass_r[rsteps-1]/(radius*radius*radius));

    //printf("total trials=%d, z3 difference=%le\n",trialsz3,z3test);
    printf("  %lf  %lf  %lf  %lf  %lf  %lf\n",sigma/(2.0*PI),r_pos[rsteps-1],boundarypos,omega,z3init,z3test);
//    printf("%lf  %lf  %lf  %lf  %lf  %lf\n",radius/(1e+5),sigma/(2.0*PI),mass_r[rsteps-1]/(1.989e+33),omega,z3init,z3test);
//    printf("%lf  %lf  %lf  %lf  %lf  %lf\n",boundarypos/(1e+5),sigma/(2.0*PI),mass_r[rsteps-1]/(1.989e+33),omega,z3init,z3test);
    //printf("closest omega value found = %.12le, --> sigma=%lf (rad/s), period = %.10lf (s), frequency=%lf (Hz), z3(R*) = %.10lf\n\n",omega,sigma,2.0*PI/sigma,sigma/(2.0*PI),z3init);

    //scale the y1(dr) and y2(dr) values to meet the z values at Rcc 
    if (y1>1e+20){z1=0.0;}
    core(&z1,&y1,&y2,omega,l,1,sigma);
    crust(&z1,&z2,&z3init,&z4,omega,l,sigma,1,&scale);
    //----------------------------------------------------------------------------------------------------------------//

  //option 2 just uses the given z3 and omega values
  } else if (option == 2) {
    //----------------------------------------------------------------------------------------------------------------//
    printf("Input omega value: ");
    scanf("%lf",&omegainit);
    printf("Input z3 value: ");
    scanf("%lf",&z3init);

    omega=omegainit;
    sigma=sqrt(omega*omega*G*mass_r[rsteps-1]/(radius*radius*radius));
    printf("omega value = %.12le, --> sigma=%lf (rad/s), period = %.10lf (s), frequency=%lf (Hz), z3(R*) = %.10lf\n",omega,sigma,2.0*PI/sigma,sigma/(2.0*PI),z3init);

    //initial run to find the values of y1 and z1 at Rcc
    z3=z3init;
    core(&z1,&y1,&y2,omega,l,0,sigma);
    crust(&z1,&z2,&z3,&z4,omega,l,sigma,0,&scale);
    //second run to print out properly scaled values 
//    if (y1>1e+20){z1=0.0;}
    z3=z3init;
    core(&z1,&y1,&y2,omega,l,1,sigma);
    crust(&z1,&z2,&z3,&z4,omega,l,sigma,1,&scale);
    //----------------------------------------------------------------------------------------------------------------//

  }
//  printf("finished.\n");
/*--------------------use shooting method to find sigma---------------------*/
}
