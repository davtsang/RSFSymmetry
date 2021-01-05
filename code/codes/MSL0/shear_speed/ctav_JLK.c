//Outputs mu(rho=10^?) for the input EoS, but only if the input J and L values match the EoSs J and L values (can easily change to J,K or L,K)
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "nr_codes.h"
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
  char c;
  char eosname[80]="";
  char J[8]="", L[8]="", K[8]="";
  int i,lines,underscore,j,l,smod_max;
  double gamma,mu,temp,r_low,r_high,smod_maxval,r,rstep,ctdens_av,csdens_av;
  double *J_fdata=malloc(0*sizeof(double)), *L_fdata=malloc(0*sizeof(double)), *K_fdata=malloc(0*sizeof(double)), *f_fdata=malloc(0*sizeof(double));
  double rho_ratio,bary_ratio,pres_cc,g1_cc;
  double t1,t2,t3,t4,t5,t6;

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
  strcat(eosname,"../EoSs/eos_LQDROP_SkXi450_");
  strcat(eosname,argv[1]);
  strcat(eosname,".0_");
  strcat(eosname,argv[2]);
  strcat(eosname,".00_");
  snprintf(K,8,"%.2lf",atof(argv[3])/100.0);
  strcat(eosname,K);
  strcat(eosname,"_2.2_Imid_.dat");
  FILE * eos = fopen(eosname, "r");

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

/*-----------------------------read J_L_f file------------------------------*/
  FILE * J_L_f = fopen("../starsh_output.txt", "r");

  //get length of input file
  lines=0;
  c = getc(J_L_f);
  while (c != EOF){
    if (c == '\n'){
      lines=lines+1;
    }
    c = getc(J_L_f);
  }
  fseek(J_L_f,0,SEEK_SET);

  J_fdata=realloc(J_fdata,lines*sizeof(double));
  L_fdata=realloc(L_fdata,lines*sizeof(double));
  K_fdata=realloc(K_fdata,lines*sizeof(double));
  f_fdata=realloc(f_fdata,lines*sizeof(double));

  //read input files
  for(i=0;i<lines;i++){
    fscanf(J_L_f,"%le  %le  %le  %le  %le  %le  %le  %le  %le  %le",&J_fdata[i],&L_fdata[i],&K_fdata[i],&t1,&f_fdata[i],&t2,&t3,&t4,&t5,&t6);
  }

  fclose(J_L_f);
/*-----------------------------read J_L_f file------------------------------*/

/*----------get dimensionless properties of star and setup splines----------*/
  d_rpos=realloc(d_rpos,rsteps*sizeof(double));
  d_density=realloc(d_density,rsteps*sizeof(double));
  d_gravity=realloc(d_gravity,rsteps*sizeof(double));
  d_pressure=realloc(d_pressure,rsteps*sizeof(double));
  d_mass_r=realloc(d_mass_r,rsteps*sizeof(double));
  d_bary_dens=realloc(d_bary_dens,rsteps*sizeof(double));
  d_smod=realloc(d_smod,rsteps*sizeof(double));
  d_g1=realloc(d_g1,rsteps*sizeof(double));

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

  //find Rcc using the crust-core transition baryon density output by the TOV
  r_low=r_pos[0];
  r_high=r_pos[rsteps-1];
  boundarypos=(r_low+r_high)/2.0;
  do{
    if (r_high-r_low<1.0){
      break;
    }
    if (spline_value(rsteps,bary_dens,r_pos,d_bary_dens,boundarypos)>1e+37){
      r_low=boundarypos;
    } else if (spline_value(rsteps,bary_dens,r_pos,d_bary_dens,boundarypos)<1e+37){
      r_high=boundarypos;
    }
    boundarypos=(r_low+r_high)/2.0;
  } while (1);




  for (j=0;j<lines;j++){
    if (atof(argv[1]) == J_fdata[j] && atof(argv[2]) == L_fdata[j]){
      //find the max value of smod in this EoS
      smod_maxval=0.0;
      for (i=0;i<rsteps-1;i++){     
        if (smod_arr[i]>=smod_maxval){
          smod_maxval=smod_arr[i];
          smod_max=i;
        }
      }

      //a few useful quantities
      rho_ratio=spline_value(rsteps,density,r_pos,d_density,boundarypos)/density[smod_max];
      bary_ratio=spline_value(rsteps,bary_dens,r_pos,d_bary_dens,boundarypos)/bary_dens[smod_max];
      pres_cc=spline_value(rsteps,pressure,r_pos,d_pressure,boundarypos);
      g1_cc=spline_value(rsteps,g1,r_pos,d_g1,boundarypos);






//      for (i=0;i<rsteps-1;i++){                                                                              //mu at chosen density
//        if (density[i]>=1e+14 && density[i+1]<1e+14){
//          printf("%s  %s  %s  %le  %le  %lf\n",J,L,K_fdata[j],spline_value(rsteps,smod_arr,r_pos,d_smod,r_pos[i]+(density[i+1]-density[i])*(r_pos[i+1]-r_pos[i])/((1e+14)-density[i])),1e+14,f_fdata[j]);
//        }
//      }

//      printf("%s  %s  %s  %le  %le  %lf\n",J,L,K_fdata[j],spline_value(rsteps,smod_arr,r_pos,d_smod,boundarypos+0.05*(radius-boundarypos)),boundarypos+0.05*(radius-boundarypos),f_fdata[j]);                                                                                      //mu at chosen crust fraction 

//      printf("%s  %s  %s  %le  %le  %lf\n",J,L,K_fdata[j],radius-boundarypos,1.0-boundarypos/radius,f_fdata[j]);             //crust thickness and crust as a fraction of the star


 //     printf("%s  %s  %s  %le  %le  %lf\n",J,L,K_fdata[j],smod_maxval,r_pos[smod_max],f_fdata[j]);            //max value of smod in this EoS



      ctdens_av=0.0;
      csdens_av=0.0;
      r=boundarypos-0.5;
      rstep=1.0;
      do{                                                                                        //average density-weighted shear speed in the crust
        if (r+2.5>radius){
          rstep=(radius-r+0.5);
          r=(r+0.5+radius)/2.0;
        } else {
          r=r+1.0;
          rstep=1.0;
        }
        ctdens_av=ctdens_av+sqrt(spline_value(rsteps,smod_arr,r_pos,d_smod,r)*spline_value(rsteps,density,r_pos,d_density,r))*r*r*rstep;
        csdens_av=csdens_av+sqrt(spline_value(rsteps,g1,r_pos,d_g1,r)*spline_value(rsteps,pressure,r_pos,d_pressure,r)*spline_value(rsteps,density,r_pos,d_density,r))*r*r*rstep;
      } while(r+1.5<=radius);
      ctdens_av=ctdens_av/(mass_r[rsteps-1]-spline_value(rsteps,mass_r,r_pos,d_mass_r,boundarypos));
      csdens_av=csdens_av/(mass_r[rsteps-1]-spline_value(rsteps,mass_r,r_pos,d_mass_r,boundarypos));
      printf("%.2lf  %.2lf  %.2lf  %le  %lf\n",J_fdata[j],L_fdata[j],K_fdata[j],sqrt(spline_value(rsteps,smod_arr,r_pos,d_smod,boundarypos)/spline_value(rsteps,density,r_pos,d_density,boundarypos)),f_fdata[j]);
        //((ctavdivr)) -> J, L, K, ctdens_av/Rcc, ctdens_av/Rstar, ctdens_av/(Rstar-Rcc), f
        //((densweight_ct_av)) -> J, L, K, ctdens_av, ?, f
        //((densweight_ctav_rtests)) -> J, L, K, ctdens/r_av, ctdens/Rcc_av, ctdens*Rcc/r^2_av, ctdens/Rstar_av, ctdens*Rstar/r^2_av, ctdens/Rcrust_av, ctdens*Rcrust/r^2_av, f
        //((densweight_ctavdiv_vals)) -> J, L, K, ctdens/csdens,ctdens/cscc,ctdens/pcc, f




//      printf("%s  %s  %s  %le  %le  %lf\n",J,L,K_fdata[j],pressure[0],r_pos[0],f_fdata[j]);               //pressure at the specified location



//      printf("%s  %s  %s  %le  %le  %lf\n",J,L,K_fdata[j],sqrt(smod_maxval/density[smod_max]),r_pos[smod_max],f_fdata[j]);  //shear speed at the specified location
//      printf("%s  %s  %s  %le  %le  %lf\n",J,L,K_fdata[j],sqrt(pres_cc*g1_cc/spline_value(rsteps,density,r_pos,d_density,boundarypos)),boundarypos,f_fdata[j]);  //sound speed at the specified location


//        printf("%s  %s  %s  %le  %le  %lf\n",J,L,K_fdata[j],sqrt(pres_cc*g1_cc/spline_value(rsteps,density,r_pos,d_density,boundarypos))/sqrt(smod_maxval/density[smod_max]),boundarypos,f_fdata[j]);          //  cs(at Rcc) divided by ct(at mumax)


//        printf("%s  %s  %s  %le  %le  %lf\n",J,L,K_fdata[j],(sqrt(pres_cc*g1_cc/spline_value(rsteps,density,r_pos,d_density,boundarypos))/sqrt(smod_maxval/density[smod_max]))/boundarypos,(sqrt(pres_cc*g1_cc/spline_value(rsteps,density,r_pos,d_density,boundarypos))/sqrt(smod_maxval/density[smod_max]))/radius,f_fdata[j]);            //  (cs(at Rcc) divided by ct(at mumax)) divided by Rcc or R*


//      printf("%s  %s  %s  %le  %le  %lf\n",J,L,K_fdata[j],r_pos[smod_max]-boundarypos,(r_pos[smod_max]-boundarypos)/(radius-boundarypos),f_fdata[j]);  //positon of mumax as distance from Rcc and as a fraction of the crust


//        printf("%s  %s  %s  %le  %lf\n",J,L,K_fdata[j],sqrt(pres_cc*g1_cc/spline_value(rsteps,density,r_pos,d_density,boundarypos))/sqrt(smod_maxval/density[smod_max]),f_fdata[j]);          //  cs(at Rcc) divided by ct(at mumax)


//        printf("%s  %s  %s  %le  %le  %le  %le  %le  %lf\n",J,L,K_fdata[j],sqrt(pres_cc*g1_cc/smod_maxval)/boundarypos,spline_value(rsteps,density,r_pos,d_density,boundarypos),pres_cc,pressure[0],spline_value(rsteps,bary_dens,r_pos,d_bary_dens,boundarypos),f_fdata[j]);          //  p times Gamma_1 (at Rcc) divided by mu (at mumax)


//        printf("%s  %s  %s  %le  %le  %le  %le  %le  %lf\n",J,L,K_fdata[j],sqrt(pres_cc*g1_cc/smod_maxval)/boundarypos,spline_value(rsteps,density,r_pos,d_density,boundarypos),pres_cc,pressure[0],spline_value(rsteps,bary_dens,r_pos,d_bary_dens,boundarypos),f_fdata[j]);          // a bunch of different values ((pxgamdivmu_ccvalues))


//        printf("%s  %s  %s  %le  %le  %le  %le  %le  %le  %le  %lf\n",J,L,K_fdata[j],(sqrt(pres_cc*g1_cc/smod_maxval))/spline_value(rsteps,bary_dens,r_pos,d_bary_dens,boundarypos),(sqrt(pres_cc*g1_cc/smod_maxval))/sqrt(spline_value(rsteps,bary_dens,r_pos,d_bary_dens,boundarypos)),(sqrt(pres_cc*g1_cc/smod_maxval))*spline_value(rsteps,bary_dens,r_pos,d_bary_dens,boundarypos),(sqrt(pres_cc*g1_cc/smod_maxval))*sqrt(spline_value(rsteps,bary_dens,r_pos,d_bary_dens,boundarypos)),spline_value(rsteps,bary_dens,r_pos,d_bary_dens,boundarypos)/bary_dens[smod_max],spline_value(rsteps,bary_dens,r_pos,d_bary_dens,boundarypos),bary_dens[smod_max],f_fdata[j]);  //((pxgamdivmu_nbcctests))



//        printf("%s  %s  %s  %le  %le  %le  %le  %le  %le  %lf\n",J,L,K_fdata[j],(sqrt(pres_cc*g1_cc/smod_maxval))/bary_ratio,(sqrt(pres_cc*g1_cc/smod_maxval))/sqrt(bary_ratio),(sqrt(pres_cc*g1_cc/smod_maxval))/rho_ratio,(sqrt(pres_cc*g1_cc/smod_maxval))/sqrt(rho_ratio),bary_ratio,rho_ratio,f_fdata[j]);  //((pxgamdivmu_ratiotests))


//        printf("%s  %s  %s  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %lf\n",J,L,K_fdata[j],(sqrt(pres_cc*g1_cc/smod_maxval))/sqrt(rho_ratio),(sqrt(pres_cc*g1_cc/smod_maxval))/(sqrt(rho_ratio)*sqrt(rho_ratio)),(sqrt(pres_cc*g1_cc/smod_maxval))/(sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio)),(sqrt(pres_cc*g1_cc/smod_maxval))/(sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio)),(sqrt(pres_cc*g1_cc/smod_maxval))/(sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio)),(sqrt(pres_cc*g1_cc/smod_maxval))/(sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio)),(sqrt(pres_cc*g1_cc/smod_maxval))/(sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio)),(sqrt(pres_cc*g1_cc/smod_maxval))/(sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio)),(sqrt(pres_cc*g1_cc/smod_maxval))/(sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio)),(sqrt(pres_cc*g1_cc/smod_maxval))/(sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio)),f_fdata[j]);  //((pxgamdivmu_densitypowers))


//        printf("%s  %s  %s  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %lf\n",J,L,K_fdata[j],1.0/((sqrt(pres_cc*g1_cc/smod_maxval))/sqrt(rho_ratio)),1.0/((sqrt(pres_cc*g1_cc/smod_maxval))/(sqrt(rho_ratio)*sqrt(rho_ratio))),1.0/((sqrt(pres_cc*g1_cc/smod_maxval))/(sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio))),1.0/((sqrt(pres_cc*g1_cc/smod_maxval))/(sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio))),1.0/((sqrt(pres_cc*g1_cc/smod_maxval))/(sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio))),1.0/((sqrt(pres_cc*g1_cc/smod_maxval))/(sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio))),1.0/((sqrt(pres_cc*g1_cc/smod_maxval))/(sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio))),1.0/((sqrt(pres_cc*g1_cc/smod_maxval))/(sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio))),1.0/((sqrt(pres_cc*g1_cc/smod_maxval))/(sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio))),1.0/((sqrt(pres_cc*g1_cc/smod_maxval))/(sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio)*sqrt(rho_ratio))),1.0/(boundarypos),1.0/(radius),f_fdata[j]);  //((all_pxgamdivmu_densitypowers))   (1/pxgamdivmu_densitypowers, with more powers of rho_ratio)

//      printf("%s  %s  %s  %le  %le  %le  %lf\n",J,L,K_fdata[j],boundarypos,r_pos[smod_max],radius,f_fdata[j]);         //radius values (note, 7 columns, not 6)

      break;
    }
  }
}
