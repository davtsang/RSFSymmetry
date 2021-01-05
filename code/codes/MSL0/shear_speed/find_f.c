#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "nr_codes.h"

void find_JL(int lines, double f_J_K[], double L_J_K[], double f_target, double Jval, double Kval);
void find_LJ(int lines, double f_L_K[], double J_L_K[], double f_target, double Lval, double Kval);
void sort_JL(int file_lines, int *lines,double **f_J_K, double **L_J_K, double Jval, double Kval);
void sort_LJ(int file_lines, int *lines,double **f_L_K, double **J_L_K, double Lval, double Kval);
void contour(double f_target, double Kval, int lines);

int ffound;
double *J_ffound, *L_ffound;
double *J, *L, *K, *f;

void main(int argc, char *argv[])
{
  double avct, Kval;
  int i, lines;
  char c;
  double M,f_temp;

  if (argc != 3) {printf("Input: ./find_f [average shear speed (cm/s)] [Ksym (MeV)]\n");exit(0);}

  avct=atof(argv[1]);
  Kval=atof(argv[2]);

  FILE *in_file=fopen("all_ct_nb37.txt","r"); 
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
    fscanf(in_file,"%le  %le  %le  %le  %le",&J[i],&L[i],&K[i],&f[i],&f_temp);
  }

  contour(avct,Kval,lines);
}


//find the J,L values which give the target frequency by varying L along cubic splines
void contour(double f_target, double Kval, int lines){
  double **L_J_Kval=malloc(0*sizeof(double)), **f_J_Kval=malloc(0*sizeof(double)), **d_f_J_Kval=malloc(0*sizeof(double));
  int *lines_J=malloc(0*sizeof(int));

  double **J_L_Kval=malloc(0*sizeof(double)), **f_L_Kval=malloc(0*sizeof(double)), **d_f_L_Kval=malloc(0*sizeof(double));
  int *lines_L=malloc(0*sizeof(int));

  double *Jvalues=malloc(1*sizeof(double)),*Lvalues=malloc(1*sizeof(double));
  int Jvals=1,Lvals=1;

  double *L_temp=malloc(0*sizeof(double)),*J_temp=malloc(0*sizeof(double)),*K_temp=malloc(0*sizeof(double));
  int i,i_min,i_sort=0,j;
  double L_min;

  double *d_J_ffound=malloc(0*sizeof(double));
  double testL=1000.0,L_max=0.0;

  ffound=0;

  //find the J and L grid in the input data
  Jvalues[0]=J[0];
  Lvalues[0]=L[0];

  for (j=1;j<lines;j++){
    for (i=0;i<j;i++) {
      if (J[i]==J[j]) i=10000;
    }
    if (i<10000) {
      Jvals=Jvals+1;
      Jvalues=realloc(Jvalues,Jvals*sizeof(double));
      Jvalues[Jvals-1]=J[j];
    }

    for (i=0;i<j;i++) {
      if (L[i]==L[j]) i=10000;
    }
    if (i<10000) {
      Lvals=Lvals+1;
      Lvalues=realloc(Lvalues,Lvals*sizeof(double));
      Lvalues[Lvals-1]=L[j];
    }
  }

//  for (i=0;i<Jvals;i++){printf("%d  %lf\n",Jvals,Jvalues[i]);}for (i=0;i<Lvals;i++){printf("%d  %lf\n",Lvals,Lvalues[i]);}  //Prints out the J,L grid

  lines_J=realloc(lines_J,Jvals*sizeof(int));
  L_J_Kval=realloc(L_J_Kval,Jvals*sizeof(double));
  f_J_Kval=realloc(f_J_Kval,Jvals*sizeof(double));
  d_f_J_Kval=realloc(d_f_J_Kval,Jvals*sizeof(double));
  for (i=0;i<Jvals;i++){
    sort_JL(lines,&lines_J[i],&f_J_Kval[i],&L_J_Kval[i],Jvalues[i],Kval);    //sort EoSs into arrays for constant J,K=Kval
    find_JL(lines_J[i],f_J_Kval[i],L_J_Kval[i],f_target,Jvalues[i],Kval);    //interpolate along L to find the target frequency for each J,Kval
  }

  lines_L=realloc(lines_L,Lvals*sizeof(int));
  J_L_Kval=realloc(J_L_Kval,Lvals*sizeof(double));
  f_L_Kval=realloc(f_L_Kval,Lvals*sizeof(double));
  d_f_L_Kval=realloc(d_f_L_Kval,Lvals*sizeof(double));
  for (i=0;i<Lvals;i++){
    sort_LJ(lines,&lines_L[i],&f_L_Kval[i],&J_L_Kval[i],Lvalues[i],Kval);    //sort EoSs into arrays for constant L,K=Kval
    find_LJ(lines_L[i],f_L_Kval[i],J_L_Kval[i],f_target,Lvalues[i],Kval);    //interpolate along J to find the target frequency for each L,Kval
  }

  //sort by L value (lowest L value first)
  L_temp=realloc(L_temp,ffound*sizeof(double));
  J_temp=realloc(J_temp,ffound*sizeof(double));
  K_temp=realloc(K_temp,ffound*sizeof(double));
  do {
    L_min=1e+10;
    for (i=0;i<ffound;i++) {
      if (J_ffound[i]<L_min){L_min=J_ffound[i];i_min=i;}
    }
    if (L_min<1e+10) {
      L_temp[i_sort]=L_ffound[i_min];
      J_temp[i_sort]=J_ffound[i_min];
      i_sort=i_sort+1;
      J_ffound[i_min]=1e+10;
    } else {break;}    
  } while (1);

  for (i=0;i<ffound;i++) {
    L_ffound[i]=L_temp[i];
    J_ffound[i]=J_temp[i];
  }

//  for (i=0;i<ffound;i++) {
//    printf("%le  %le  %le\n",J_ffound[i],L_ffound[i],f_target);
//  }printf("\n");


  //if cubic splines are worth using, use them to interpolate the J,L values that give f=target
  if (ffound>20000){  
    for (i=0;i<ffound;i++){
      J_ffound[i]=J_ffound[i]+i/1e+4;  //make small adjustments to each J value to make sure splines don't become linear due to J1=J2
    }
  
    d_J_ffound=malloc(ffound*sizeof(double));
    cubic_spline(ffound,J_ffound,L_ffound,d_J_ffound);
    // find maximum and minimum L values that have the chosen frequency (that I have data for) (for any J value)
    for(i=0;i<ffound;i++){
      if(testL>L_ffound[i]) testL=L_ffound[i];
      if(L_max<L_ffound[i]) L_max=L_ffound[i];
    }

    //print out interpolated J,L values that give the chosen frequency
    do {
      printf("%le  %le\n",spline_value(ffound,J_ffound,L_ffound,d_J_ffound,testL),testL);
      testL=testL+0.1;
    } while (testL<L_max);
  }  

  //else print out the J,L,K=Kval values without any interpolation
  else {        
    for (i=0;i<ffound;i++){
      printf("%le  %le\n",J_ffound[i],L_ffound[i]);
    }
  }
//  printf("\n");


  FILE *out_file1=fopen("J_f_mergetime.txt","w");
  for(j=0;j<Lvals;j++){
    for(i=0;i<lines;i++){
      if(L[i]==Lvalues[j] && f[i]!=0){fprintf(out_file1,"%le  %le\n",J[i],f[i]);}
    }
    fprintf(out_file1,"\n");
  }
  fclose(out_file1);

  free(L_J_Kval);free(f_J_Kval);free(d_f_J_Kval);free(lines_J);free(J_L_Kval);free(f_L_Kval);free(d_f_L_Kval);
  free(lines_L);free(Jvalues);free(Lvalues);free(L_temp);free(J_temp);free(K_temp);free(d_J_ffound);
}



// function for finding the J,L,K=x values that give the target frequency by interpolating between L values for each J value
void find_JL(int lines, double f_J_K[], double L_J_K[], double f_target, double Jval, double Kval){
  double Lval,f_val,f_val_next,d_Lval=0.01;
  double *d_f_J_K=malloc(lines*sizeof(double));
  int i,j;

  if (lines>1){
    //setup cubic splines for f(L) with constant J,K
    cubic_spline(lines,f_J_K,L_J_K,d_f_J_K);

    //loop over L values, interpolating to see what frequencies they give
    Lval=L_J_K[0];
    f_val=spline_value(lines,f_J_K,L_J_K,d_f_J_K,Lval);
    do{
      f_val_next=spline_value(lines,f_J_K,L_J_K,d_f_J_K,Lval+d_Lval);
      //printf("%le  %le  %le\n",Jval,Lval,f_val);
      if ((f_val>f_target && f_val_next<f_target) || (f_val<f_target && f_val_next>f_target)) {  //if the target is between this L and L+L_step...
        //printf("%le  %le  %le  %le  %le\n",Jval,Lval,Kval,f_val,f_val_next);
        j=1;
        for (i=0;i<ffound;i++){
          if (fabs(J_ffound[i]-Jval)<0.5 && fabs(L_ffound[i]-Lval)<1.0) j=0;  //Ignore interpolated values that are very close to other points, 
        }                                                                     //because it causes problems for the splines
        if (j) {
          ffound=ffound+1;
          J_ffound=realloc(J_ffound,ffound*sizeof(double));
          L_ffound=realloc(L_ffound,ffound*sizeof(double));
          J_ffound[ffound-1]=Jval;
          L_ffound[ffound-1]=Lval;
        }
      }
  
      f_val=f_val_next;
      Lval=Lval+d_Lval;
    } while (Lval<=L_J_K[lines-1]-d_Lval);
  }
  //printf("\n");
  free(d_f_J_K);
}

// function for finding the J,L,K=x values that give the target frequency by interpolating between L values for each J value
void find_LJ(int lines, double f_L_K[], double J_L_K[], double f_target, double Lval, double Kval){
  double Jval,f_val,f_val_next,d_Jval=0.01;
  double *d_f_L_K=malloc(lines*sizeof(double));
  int i,j;

  if (lines>1){
    //setup cubic splines for f(J) with constant L,K
    d_f_L_K=realloc(d_f_L_K,lines*sizeof(double));
    cubic_spline(lines,f_L_K,J_L_K,d_f_L_K);

    //loop over J values, interpolating to see what frequencies they give
    Jval=J_L_K[0];
    f_val=spline_value(lines,f_L_K,J_L_K,d_f_L_K,Jval);
    do{
      f_val_next=spline_value(lines,f_L_K,J_L_K,d_f_L_K,Jval+d_Jval);
      //printf("%le  %le  %le\n",Jval,Lval,f_val);
      if ((f_val>f_target && f_val_next<f_target) || (f_val<f_target && f_val_next>f_target)) {  //if the target is between this L and L+L_step...
        //printf("%le  %le  %le  %le  %le\n",Jval,Lval,Kval,f_val,f_val_next);
        j=1;
        for (i=0;i<ffound;i++){
          if (fabs(J_ffound[i]-Jval)<0.5 && fabs(L_ffound[i]-Lval)<1.0) j=0;  //Ignore interpolated values that are very close to other points, 
        }                                                                     //because it causes problems for the splines
        if (j) {
          ffound=ffound+1;
          J_ffound=realloc(J_ffound,ffound*sizeof(double));
          L_ffound=realloc(L_ffound,ffound*sizeof(double));
          J_ffound[ffound-1]=Jval;
          L_ffound[ffound-1]=Lval;
        }
      }
  
      f_val=f_val_next;
      Jval=Jval+d_Jval;
    } while (Jval<=J_L_K[lines-1]-d_Jval);
  }
  //printf("\n");
  free(d_f_L_K);
}


//function for sorting EoSs into arrays with the same J and K values
void sort_JL(int file_lines, int *lines,double **f_J_K, double **L_J_K, double Jval, double Kval){
  int i;

  (*lines)=0;
  for(i=0;i<file_lines;i++){
    if (J[i]==Jval && f[i]>0.0) {
      (*lines)=(*lines)+1;
      *L_J_K=realloc(*L_J_K,(*lines)*sizeof(double));
      *f_J_K=realloc(*f_J_K,(*lines)*sizeof(double));
      (*L_J_K)[(*lines)-1]=L[i];
      (*f_J_K)[(*lines)-1]=f[i];
    }
  }
}
//COMBINE THESE 2 BY INPUTTING J/L ARRAYS RATHER THAN USING GLOBAL VARIABLES
//function for sorting EoSs into arrays with the same L and K values
void sort_LJ(int file_lines, int *lines,double **f_L_K, double **J_L_K, double Lval, double Kval){
  int i;

  (*lines)=0;
  for(i=0;i<file_lines;i++){
    if (L[i]==Lval && f[i]>0.0) {
      (*lines)=(*lines)+1;
      *J_L_K=realloc(*J_L_K,(*lines)*sizeof(double));
      *f_L_K=realloc(*f_L_K,(*lines)*sizeof(double));
      (*J_L_K)[(*lines)-1]=J[i];
      (*f_L_K)[(*lines)-1]=f[i];
    }
  }
}
