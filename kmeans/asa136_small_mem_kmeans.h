void kmns ( int m, int n, double *c[], int k, /*int ic1[], int nc[], */
  int iter, int *ifault );
void optra ( int m, int n, double *c[], int k, int ic1[], 
  int ic2[], int nc[], double an1[], double an2[], int ncp[], double d[], 
  int itran[], int live[], int *indx );
void qtran ( int m, int n, double *c[], int k, int ic1[], 
  int ic2[], int nc[], double an1[], double an2[], int ncp[], double d[], 
  int itran[], int *indx );
double r8_huge ( void );
void timestamp ( void );
double distance(double a[], double c[], int len_a, int len_c);
//void distance_test(void * args);
void *find_closest(void * args);
void * optran_p(void * args);
void * qtran_p(void * args);
void qtran_o ( int m, int n, double *c[], int k, int ic1[], 
  int ic2[], int nc[], double an1[], double an2[], int ncp[], double d[], 
  int itran[], int *indx );
struct thread_args{
	int begin;
	int end;
	double ** center;
	int * ic1;
	int * ic2;
//	double * dt;
	int n;
	int k;	
};
struct thread_optran{
	int begin;
	int end;
	double ** center;
	int *ic1, *ic2, *nc, *ncp, *itran, *live, *indx;
	double *an1, *an2, *d;
	int m,n,k;	
};
struct thread_qtran{
	int begin;
	int end;
	double ** center;
	int *ic1, *ic2, *nc, *ncp, *itran, *indx;
	double *an1, *an2, *d;
	int m,n,k;	
};
