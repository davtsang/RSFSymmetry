void RK4step(double vars[], int n, double r, double dr, double var_out[], void (*derivatives)(double, double[], double[]));
void RKCKstep(double vars[], int n, double r, double dr, double var_out[], double err[], void (*derivatives)(double, double[], double[]));
void RK_var_step(double vars[], int n, double *r, double dr_try, double tol, double *dr_used, double *dr_next, void (*derivatives)(double, double[], double[]));
void cubic_spline(int size, double y[size], double x[size], double d_y[size]);
double spline_value(int size, double y[size], double x[size],  double d_y[size], double x_val);
double linear_int(int size,double value,double dep_array[],double ind_arr[]);
