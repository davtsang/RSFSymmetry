#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "nr_codes.h"

double *J,*L,*K,*f;
int ffound,ffound2;
double *J_ffound, *L_ffound, *K_ffound, *J_ffound2, *L_ffound2, *K_ffound2;

void main(int argc, char *argv[])
{
  char c;
  int lines,i,j,k,m;
  double M,t1,t2,t3,t4,t5;
  double f_target,K_input,Jval,Lval,Kval;

  double *L_temp=malloc(0*sizeof(double)),*J_temp=malloc(0*sizeof(double));
  int i_min,i_sort=0;
  double L_min;

  double *d_J_ffound2=malloc(0*sizeof(double));
  double testL=1000.0,L_max=0.0;

  f_target=atof(argv[1]);
  K_input=atof(argv[2]);

  FILE *in_file=fopen("starsh_output.txt","r"); 
  //get length of input file
  lines=0;
  c = getc(in_file);
  while (c != EOF){
    if (c == '\n'){
      lines=lines+1;
    }
    c = getc(in_file);
  }
  fseek(in_file,0,SEEK_SET);

  J=realloc(J,lines*sizeof(double));
  L=realloc(L,lines*sizeof(double));
  K=realloc(K,lines*sizeof(double));
  f=realloc(f,lines*sizeof(double));

  //read input files
  for(i=0;i<lines;i++){                                      //change to L,K,J (J,K,L) to get L vs K (J vs K) contours
    fscanf(in_file,"%le  %le  %le  %le  %le  %le  %le  %le  %le  %le",&J[i],&L[i],&K[i],&M,&f[i],&t1,&t2,&t3,&t4,&t5);
  }

  //linear interpolation between neighbouring EoSs to find the target f value
  ffound=0;
  for(i=0;i<lines-1;i++){  
    for(j=i+1;j<lines;j++){  
      if (((J[i]==J[j]) && (L[i]+2 > L[j]) && (L[i]-2 < L[j]) && (K[i]+15 > K[j]) && (K[i]-15 < K[j])) || ((J[i]+1.4==J[j]) && (L[i]+9.6 > L[j]) && (L[i]+9.2 < L[j]) && (K[i]+27.3 > K[j]) && (K[i]+24.3 < K[j])) || ((J[i]-1.4==J[j]) && (L[i]-9.2 > L[j]) && (L[i]-9.6 < L[j]) && (K[i]-24.3 > K[j]) && (K[i]-27.3 < K[j]))){
        if (((f[i] > f_target) && (f[j] < f_target)) || ((f[i] < f_target) && (f[j] > f_target))){
          Jval=J[j]+(f_target-f[j])/(f[i]-f[j])*(J[i]-J[j]);
          Lval=L[j]+(f_target-f[j])/(f[i]-f[j])*(L[i]-L[j]);
          Kval=K[j]+(f_target-f[j])/(f[i]-f[j])*(K[i]-K[j]);
          ffound=ffound+1;
          J_ffound=realloc(J_ffound,ffound*sizeof(double));
          L_ffound=realloc(L_ffound,ffound*sizeof(double));
          K_ffound=realloc(K_ffound,ffound*sizeof(double));
          J_ffound[ffound-1]=Jval;
          L_ffound[ffound-1]=Lval;
          K_ffound[ffound-1]=Kval;
        }
      }      
    }
  }


/*
  //PRINT OUT HERE TO GET PLANES IN J,L,Ksym
  for(i=0;i<ffound;i++){ 
    printf("%le  %le  %le\n",J_ffound[i],L_ffound[i],K_ffound[i]);
  } exit(0);
*/



  if (K_input<1000){
/*
    //put interpolated points that are close to the chosen Ksym value in an array (could interpolate to find the chosen Ksym instead?)
    ffound2=0;
    for(i=0;i<ffound;i++){ 
      if (fabs(K_ffound[i]-K_input) < 10.0){   //change to allow a wider range around the chosen Ksym
        ffound2=ffound2+1;
        J_ffound2=realloc(J_ffound2,ffound2*sizeof(double));
        L_ffound2=realloc(L_ffound2,ffound2*sizeof(double));
        J_ffound2[ffound2-1]=J_ffound[i];
        L_ffound2[ffound2-1]=L_ffound[i];
      }
    }
*/

  //linear interpolation between neighbouring interpolated points to find the target K value
  ffound2=0;
  for(i=0;i<ffound-1;i++){ 
    for(j=i+1;j<ffound;j++){  
      if ((L_ffound[i]+2.0 > L_ffound[j]) && (L_ffound[i]-2.0 < L_ffound[j]) && (J_ffound[i]+0.5 > J_ffound[j]) && (J_ffound[i]-0.5 < J_ffound[j])){
        if (((K_ffound[i] > K_input) && (K_ffound[j] < K_input)) || ((K_ffound[i] < K_input) && (K_ffound[j] > K_input))){
          Jval=J_ffound[j]+(K_input-K_ffound[j])/(K_ffound[i]-K_ffound[j])*(J_ffound[i]-J_ffound[j]);
          Lval=L_ffound[j]+(K_input-K_ffound[j])/(K_ffound[i]-K_ffound[j])*(L_ffound[i]-L_ffound[j]);
          m=1;
          for (k=0;k<ffound2;k++){if (fabs(Jval-J_ffound2[k]) < 0.03) m=0;}  //ignore values that are the same as previously found ones
          if (m){
            ffound2=ffound2+1;
            J_ffound2=realloc(J_ffound2,ffound2*sizeof(double));
            L_ffound2=realloc(L_ffound2,ffound2*sizeof(double));
            J_ffound2[ffound2-1]=Jval;
            L_ffound2[ffound2-1]=Lval;
          }
        }
      }      
    }
  }

















    //sort by L value (lowest L value first)
    L_temp=realloc(L_temp,ffound2*sizeof(double));
    J_temp=realloc(J_temp,ffound2*sizeof(double));
    do {
      L_min=1e+10;
      for (i=0;i<ffound2;i++) {
        if (L_ffound2[i]<L_min){L_min=L_ffound2[i];i_min=i;}
      }
      if (L_min<1e+10) {
        L_temp[i_sort]=L_ffound2[i_min];
        J_temp[i_sort]=J_ffound2[i_min];
        i_sort=i_sort+1;
        L_ffound2[i_min]=1e+10;
      } else {break;}    
    } while (1);
  
    for (i=0;i<ffound2;i++) {
      L_ffound2[i]=L_temp[i];
      J_ffound2[i]=J_temp[i];
    }

    //if cubic splines are worth using, use them to interpolate the J,L values that give f=target
    if (ffound2>2000){  //set to (large) to use linear interpolation, or (2) to do cubic splines
      for (i=0;i<ffound2;i++){
        J_ffound2[i]=J_ffound2[i]+i/1e+4;  //make small adjustments to each J value to make sure splines don't become linear due to J1=J2
      }
    
      d_J_ffound2=malloc(ffound2*sizeof(double));
      cubic_spline(ffound2,J_ffound2,L_ffound2,d_J_ffound2);
      // find maximum and minimum L values that have the chosen frequency (that I have data for) (for any J value)
      for(i=0;i<ffound2;i++){
        if(testL>L_ffound2[i]) testL=L_ffound2[i];
        if(L_max<L_ffound2[i]) L_max=L_ffound2[i];
      }
  
      //print out interpolated J,L values that give the chosen frequency
      do {
        printf("%le  %le\n",spline_value(ffound2,J_ffound2,L_ffound2,d_J_ffound2,testL),testL);
        testL=testL+0.1;
      } while (testL<L_max);
    }  
  
    //else print out the J,L,K=Kval values without any interpolation
    else {        
      for (i=0;i<ffound2;i++){
        printf("%le  %le  %le\n",J_ffound2[i],L_ffound2[i],K_input);
      }
      if (ffound2>0){printf("\n");}
    }





  } else {
    for(i=0;i<ffound;i++){ 
      printf("%lf  %lf  %lf\n",J_ffound[i],L_ffound[i],K_ffound[i]);
    }
  }

}
