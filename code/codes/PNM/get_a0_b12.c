#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "nr_codes.h"

double *J, *L, *K, *f, *a0, *b12;

void main()
{
double Cl, Ck,M,t1,t2,t3,t4,t5;
int lines,i,j;
char c;

  FILE *in_file=fopen("J_L_K_f.txt","r"); 
  FILE *out_file=fopen("test1.txt","w"); 
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
  a0=realloc(a0,lines*sizeof(double));
  b12=realloc(b12,lines*sizeof(double));

  //read input files
  for(i=0;i<lines;i++){                                      //change to L,K,J (J,K,L) to get L vs K (J vs K) contours
    fscanf(in_file,"%le  %le  %le  %le  %le  %le  %le  %le  %le  %le",&J[i],&L[i],&K[i],&M,&f[i],&t1,&t2,&t3,&t4,&t5);
  }

  
  for(i=0;i<lines;i++){                                      //change to L,K,J (J,K,L) to get L vs K (J vs K) contours
    Cl=L[i]-6.7*J[i];
    Ck=K[i]-18.4*J[i];
    b12[i]=(Ck-5.0*Cl+7.2)/7.79;
    a0[i]=-(Cl-1.56*b12[i]+59.22)/19.47;
    for(j=0;j<i;j++){
      if ((a0[i]+0.05>a0[j]) && (a0[i]-0.05<a0[j])) {printf("%le  %le  %d  %d\n",a0[i],a0[j],i,j); a0[i]=a0[j]; break;}
    }
    for(j=0;j<i;j++){
      if ((b12[i]+0.5>b12[j]) && (b12[i]-0.5<b12[j])) {printf("%le  %le  %d  %d\n",b12[i],b12[j],i,j); b12[i]=b12[j]; break;}
    }
    fprintf(out_file,"%le  %.3le  %.3le  %le\n",J[i],a0[i],b12[i],f[i]);
  }




}
