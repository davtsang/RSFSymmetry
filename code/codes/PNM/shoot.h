void derivatives_crust(double xpos,double var_in[],double d_var_out[]);
void derivatives(double xpos,double var_in[],double d_var_out[]);
void core(double *z1,double *y1,double *y2,double omega,int l,int print,double sigma);
void crust(double *z1,double *z2,double *z3,double *z4,double omega,int l,double sigma,int print, double *scale);
void bound_test(double boundary1,double *boundary2,double *boundary3,double *value,double *d_value, int option);
