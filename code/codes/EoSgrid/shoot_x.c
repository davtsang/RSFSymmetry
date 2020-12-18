//code for the functions for the shooting method from the core to Rcc, and from the crust to Rcc
//uses (1+V)dy/dz
//this file is not currently used, it is just a backup of the version that uses d/dx (see shoot.c for the version currently being used)
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "nr_codes.h"
#include "eqbm.h"
#define TOL 1.0E-7
extern int rsteps;
extern double boundarypos, radius, boundaryposx;
extern double *g1, *d_g1, *x, *d_rpos, *r_pos, *d_mass_r, *mass_r, *density, *d_density, *pressure, *d_pressure, *gravity, *d_gravity;


//variables:  0=z1,1=z2,2=z3,3=z4,4=omega2,5=l
void derivatives_crust(double xpos,double var_in[],double d_var_out[]){
  double rval,c1val,Arval,Vtildeval,Utildeval,g1val,shear_modval,presval,densval,massval,gravval;
  double a1,a2,a3;
  rval=spline_value(rsteps,r_pos,x,d_rpos,xpos);
  g1val=spline_value(rsteps,g1,r_pos,d_g1,rval);
  shear_modval=mu_val(rval);

  densval=spline_value(rsteps,density,r_pos,d_density,rval);
  massval=spline_value(rsteps,mass_r,r_pos,d_mass_r,rval);
  gravval=spline_value(rsteps,gravity,r_pos,d_gravity,rval);
  presval=spline_value(rsteps,pressure,r_pos,d_pressure,rval);

  Utildeval=Utildevalue(rval);
  Vtildeval=Vtildevalue(rval);
  Arval=0.0;
  c1val=c1value(rval);

  a1=shear_modval/presval;
  a2=g1val-(2.0/3.0)*(shear_modval/presval);
  a3=g1val+(4.0/3.0)*(shear_modval/presval);

  d_var_out[0]=(1.0/(1.0+Vtildeval))*(-(1.0+2.0*(a2/a3))*var_in[0]+(1.0/a3)*var_in[1]+var_in[5]*(var_in[5]+1.0)*(a2/a3)*var_in[2]);
  d_var_out[1]=(1.0/(1.0+Vtildeval))*((-c1val*Vtildeval*var_in[4]-4.0*Vtildeval+Utildeval*Vtildeval+12.0*g1val*a1/a3)*var_in[0]+(Vtildeval-4.0*a1/a3)*var_in[1]+var_in[5]*(var_in[5]+1.0)*(Vtildeval-6.0*g1val*a1/a3)*var_in[2]+var_in[5]*(var_in[5]+1.0)*var_in[3]);
  if (a1 != 0.0){
    d_var_out[2]=(1.0/(1.0+Vtildeval))*(-var_in[0]+(1.0/a1)*var_in[3]);
  } else {
    d_var_out[2]=(1.0/(1.0+Vtildeval))*(-var_in[0]);  //if a1=0.0, shear_mod=0.0 --> if omega is correct then z4=0.0
  }
  d_var_out[3]=(1.0/(1.0+Vtildeval))*((Vtildeval-6.0*g1val*a1/a3)*var_in[0]-(a2/a3)*var_in[1]+(-c1val*Vtildeval*var_in[4]+(2.0/a3)*((2.0*var_in[5]*(var_in[5]+1.0)-1.0)*a1*a2+2.0*(var_in[5]*(var_in[5]+1.0)-1.0)*a1*a1))*var_in[2]+(Vtildeval-3.0)*var_in[3]);
  d_var_out[4]=0.0;
  d_var_out[5]=0.0;
}

//variables:  0=y1,1=y2,2=omega2,3=l
void derivatives(double xpos,double var_in[],double d_var_out[]){
  double rval,c1val,Arval,Vtildeval,Utildeval,g1val,presval,densval,massval,gravval;
  rval=spline_value(rsteps,r_pos,x,d_rpos,xpos);
  g1val=spline_value(rsteps,g1,r_pos,d_g1,rval);
  densval=spline_value(rsteps,density,r_pos,d_density,rval);
  massval=spline_value(rsteps,mass_r,r_pos,d_mass_r,rval);
  gravval=spline_value(rsteps,gravity,r_pos,d_gravity,rval);
  presval=spline_value(rsteps,pressure,r_pos,d_pressure,rval);
  Utildeval=Utildevalue(rval);
  Vtildeval=Vtildevalue(rval);
  Arval=0.0;
  c1val=c1value(rval);

  d_var_out[0]=(1.0/(1.0+Vtildeval))*((Vtildeval/g1val-3.0)*var_in[0]+((var_in[3]*(var_in[3]+1)/(c1val*var_in[2]))-Vtildeval/g1val)*var_in[1]);
  d_var_out[1]=(1.0/(1.0+Vtildeval))*((c1val*var_in[2]+Arval)*var_in[0]+(1.0-Utildeval-Arval)*var_in[1]);
  d_var_out[2]=0.0;
  d_var_out[3]=0.0;
}

void core(double *z1,double *y1,double *y2,double omega,int l,int print,double sigma){
  int i,small_step;
  double dx,c1omega,omega2,dxused,dxnext,xpos,r,scale,Vtilde_ss,Arval=0.0,Utilde_ss;
  double variables[4];
  for (i=0;i<rsteps;i++){
    if (r_pos[i]>=1.0){small_step=i;break;}                                   //take a small step away from the core to avoid numerical errors
  }

  if (print){
    scale=1.0/(*y1);                   //if writing data out, find scale for y1 so that U(Rcc)=1
 //   printf("scale= %le\n",scale);
  }

  omega2=omega*omega;
  c1omega=c1value(r_pos[small_step])*omega2;
  Vtilde_ss=Vtildevalue(r_pos[small_step]);
  Utilde_ss=Utildevalue(r_pos[small_step]);

  *y1=1;                                                                      // at r=radius, define y1=1
  *y2=(4.0-(Vtilde_ss/g1[small_step])+l+Arval*l/c1omega)*c1omega/l;           //get y2 value from the core boundary condition
  *y2=(*y2)/(1.0+l+Utilde_ss+Arval-(Vtilde_ss/g1[small_step])*(c1omega/l));
  *y2=(*y2)*(*y1);

  FILE * y1file = fopen("y1.txt", "w");
  FILE * y2file = fopen("y2.txt", "w");
  FILE * Ucorefile = fopen("Ucore.txt", "w");
  FILE * Vcorefile = fopen("Vcore.txt", "w");
  if (print){                                                                 //write out initial values
    fprintf(y1file,"%.10lf  %.10le\n",r_pos[small_step],*y1*scale);
    fprintf(y2file,"%.10lf  %.10le\n",r_pos[small_step],*y2*scale);
    fprintf(Ucorefile,"%.10lf  %.10le\n",r_pos[small_step]/100000,*y1*r_pos[small_step]*scale/boundarypos);
    fprintf(Vcorefile,"%.10lf  %.10le\n",r_pos[small_step]/100000,*y2/c1omega*r_pos[small_step]*scale/10.0/boundarypos);
  }

  variables[0]=*y1;
  variables[1]=*y2;
  variables[2]=omega2;
  variables[3]=l;

  xpos=x[small_step];
  dxnext=(x[small_step+1]-x[small_step])/10.0;
  dxused=dxnext;
  do {
    if (xpos+dxnext>boundaryposx){dxnext=boundaryposx-xpos;}
    RK_var_step(variables,4,&xpos,dxnext,TOL,&dxused,&dxnext,derivatives);    //take RK4 step

    if (print) {
      r=spline_value(rsteps-1,r_pos,x,d_rpos,xpos);
      fprintf(y1file,"%.10lf  %.10le\n",r,variables[0]*scale);
      fprintf(y2file,"%.10lf  %.10le\n",r,variables[1]*scale);
      fprintf(Ucorefile,"%.10lf  %.10le\n",r/100000,variables[0]*r*scale/boundarypos);
      fprintf(Vcorefile,"%.10lf  %.10le\n",r/100000,variables[1]*r/(c1value(r)*omega2)*scale/10.0/boundarypos);
    }
  } while (xpos<boundaryposx);

  *y1=variables[0];
  *y2=variables[1];
//  if (print)printf("y1_cc= %le\n",*y1*scale);
  if (print)*y1=*y1*scale;
//  if (print)printf("y2_cc= %le\n",*y2*scale);
  if (print)*y2=*y2*scale;
  fclose(y1file);
  fclose(y2file);
  fclose(Ucorefile);
  fclose(Vcorefile);
}

void crust(double *z1,double *z2,double *z3,double *z4,double omega,int l,double sigma,int print,double *scale){
  int i,small_step;
  double dx,dxused,dxnext,xpos,r,c1omega,omega2,Vtilde_ss,Utilde_ss;
  double variables[6];
  for (i=rsteps-1;i>1;i--){
    if (r_pos[i]<=r_pos[rsteps-1]-1.0){small_step=i;break;}                   //take a small step away from the surface to avoid numerical errors
  }

  if (print){
    *scale=1.0/(*z1);                   //if writing data out, find scale for z1 so that U(Rcc)=1
 //   printf("scale= %le\n",scale);
  }

  omega2=omega*omega;
  c1omega=c1value(r_pos[small_step])*omega2;
  Vtilde_ss=Vtildevalue(r_pos[small_step]);
  Utilde_ss=Utildevalue(r_pos[small_step]);

  *z1=1.0;                                                                    // at r=radius, define y1=1
  *z2=Vtilde_ss-c1omega-4.0+Utilde_ss;                                        //get y2 value from the surface boundary condition
  *z2=(*z2)/(l*(l+1.0)/c1omega-Vtilde_ss)+1.0;
  *z2=(*z2)*Vtilde_ss*(*z1);
  *z4=0.0;

  FILE * z1file = fopen("z1.txt", "w");
  FILE * z2file = fopen("z2.txt", "w");
  FILE * z3file = fopen("z3.txt", "w");
  FILE * z4file = fopen("z4.txt", "w");
  FILE * Ucrustfile = fopen("Ucrust.txt", "w");
  FILE * Vcrustfile = fopen("Vcrust.txt", "w");
  if (print) {                                                                //write out initial values
    fprintf(z1file,"%.10lf  %.10le\n",r_pos[small_step],*z1*(*scale));
    fprintf(z2file,"%.10lf  %.10le\n",r_pos[small_step],*z2*(*scale));
    fprintf(z3file,"%.10lf  %.10le\n",r_pos[small_step],*z3*(*scale));
    fprintf(z4file,"%.10lf  %.10le\n",r_pos[small_step],*z4*(*scale));
    fprintf(Ucrustfile,"%.10lf  %.10le\n",r_pos[small_step]/100000,*z1*r_pos[small_step]*(*scale)/boundarypos);
    fprintf(Vcrustfile,"%.10lf  %.10le\n",r_pos[small_step]/100000,*z3*r_pos[small_step]*(*scale)/10.0/boundarypos);
  }

  variables[0]=*z1;
  variables[1]=*z2;
  variables[2]=*z3;
  variables[3]=*z4;
  variables[4]=omega*omega;
  variables[5]=l;

  xpos=x[small_step];
  dxnext=(x[small_step-1]-x[small_step])/10.0;
  dxused=dxnext;
  do {
    if (xpos+dxnext<boundaryposx){dxnext=boundaryposx-xpos;}
    RK_var_step(variables,6,&xpos,dxnext,TOL,&dxused,&dxnext,derivatives_crust);    //take RK4 step

    if (print){
      r=spline_value(rsteps-1,r_pos,x,d_rpos,xpos);
      fprintf(z1file,"%.10lf  %.10le\n",r,variables[0]*(*scale));
      fprintf(z2file,"%.10lf  %.10le\n",r,variables[1]*(*scale));
      fprintf(z3file,"%.10lf  %.10le\n",r,variables[2]*(*scale));
      fprintf(z4file,"%.10lf  %.10le\n",r,variables[3]*(*scale));
      fprintf(Ucrustfile,"%.10lf  %.10le\n",r/100000,variables[0]*r*(*scale)/boundarypos);
      fprintf(Vcrustfile,"%.10lf  %.10le\n",r/100000,variables[2]*r*(*scale)/10.0/boundarypos);
    }
  } while (xpos>boundaryposx);

  *z1=variables[0];
  *z2=variables[1];
  *z3=variables[2];
  *z4=variables[3];
//  if (print)printf("z1_cc= %le\n",*z1*(*scale));
  if (print)*z1=*z1*(*scale);
//  if (print)printf("z2_cc= %le\n",*z2*(*scale));
  if (print)*z2=*z2*(*scale);
//  if (print)printf("z3_cc= %le\n",*z3*(*scale));
  if (print)*z3=*z3*(*scale);
//  if (print)printf("z4_cc= %le\n",*z4*(*scale));
  if (print)*z4=*z4*(*scale);
  fclose(z1file);
  fclose(z2file);
  fclose(z3file);
  fclose(z4file);
  fclose(Ucrustfile);
  fclose(Vcrustfile);
//  if (print)printf("crust scale=  %le\n",(*scale));
  if (!print){(*scale)=1.0/(variables[0]*spline_value(rsteps-1,r_pos,x,d_rpos,boundaryposx)/boundarypos);}
}

void bound_test(double boundary1,double *boundary2,double *boundary3,double *value,double *d_value,int option){
  if (boundary1>(*boundary2)) (*d_value)=-(*d_value);                          //change search direction if it has gone past the minima
  if (((*boundary3)==boundary1) && ((*d_value)>0)) (*d_value)=(*d_value)*0.1;  //change step size if close to the minima
  if ((fabs(*value)<2.0*fabs(*d_value)) && option) (*d_value)=(*d_value)*0.1;  //option to avoid negative value (used to avoid negative omega, but not used for z3)
  (*value)=(*value)+(*d_value);
  (*boundary3)=(*boundary2);   //boundary3 is the boundary 2 steps ago
  (*boundary2)=boundary1;      //boundary2 is the boundary 1 step ago
}
