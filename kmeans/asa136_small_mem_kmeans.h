void kmns ( double *a[], int m, int n, double *c[], int k, int ic1[], int nc[], 
  int iter, double wss[], int *ifault );
void optra ( double *a[], int m, int n, double *c[], int k, int ic1[], 
  int ic2[], int nc[], double an1[], double an2[], int ncp[], double d[], 
  int itran[], int live[], int *indx );
void qtran ( double *a[], int m, int n, double *c[], int k, int ic1[], 
  int ic2[], int nc[], double an1[], double an2[], int ncp[], double d[], 
  int itran[], int *indx );
double r8_huge ( void );
void timestamp ( void );
double distance(double a[], double c[], int len_a, int len_c);
