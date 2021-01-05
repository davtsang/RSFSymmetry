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
  double f_target,J_input,Jval,Lval,Kval;

  double *L_temp=malloc(0*sizeof(double)),*K_temp=malloc(0*sizeof(double));
  int i_min,i_sort=0;
  double L_min;

  double *d_J_ffound2=malloc(0*sizeof(double));
  double testL=1000.0,L_max=0.0;

  f_target=atof(argv[1]);
  J_input=atof(argv[2]);

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
      if ((J[i]==J_input) && (J_input==J[j]) && (L[i]+2 > L[j]) && (L[i]-2 < L[j]) && (K[i]+15 > K[j]) && (K[i]-15 < K[j])){
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


    //sort by K value (lowest K value first)
    K_temp=realloc(K_temp,ffound*sizeof(double));
    L_temp=realloc(L_temp,ffound*sizeof(double));
    do {
      L_min=1e+10;
      for (i=0;i<ffound;i++) {
        if (K_ffound[i]<L_min){L_min=K_ffound[i];i_min=i;}
      }
      if (L_min<1e+10) {
        K_temp[i_sort]=K_ffound[i_min];
        L_temp[i_sort]=L_ffound[i_min];
        i_sort=i_sort+1;
        K_ffound[i_min]=1e+10;
      } else {break;}    
    } while (1);
  
    for (i=0;i<ffound;i++) {
      K_ffound[i]=K_temp[i];
      L_ffound[i]=L_temp[i];
    }



  //PRINT OUT HERE TO GET PLANES IN J,L,Ksym
  for(i=0;i<ffound;i++){ 
    printf("%le  %le  %le\n",J_ffound[i],L_ffound[i],K_ffound[i]);
  }
}
