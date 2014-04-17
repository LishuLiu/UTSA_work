// Kmeans_SmallMem.cpp : Defines the entry point for the console application.
//

//# include "stdafx.h"
# include "asa136_small_mem_kmeans.h"
# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <fstream>
# include <ctime>
# include <stdio.h>
# include <sstream>
# include <string>
# include <malloc.h>
# include <vector>
# include <string>
# include <string.h>
# include <algorithm>
# include <functional>

using namespace std;

void * test01 (void*);
void * test02 (void*);
void test03 (int m, int n, int k, int iter, char* filename);

int *length;
double **matrix;
int print_flag = 0;
clock_t start,end;
int finish = 0;

double distance(double *,double *, int, int);
//****************************************************************************80

int main(int argc, char* argv[])
//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for ASA136_PRB.
//
{
	timestamp ( );
	if(argc<6)
	{
	cout << "Usage: Kmeans m n k filename iter.max"<<endl;
	return -1;
	}
	int k,m,n,iter;
	char filename[128];
	m = atoi(argv[1]);
	n = atoi(argv[2]);
	k = atoi(argv[3]);
	memset(filename,0,sizeof(filename));
	strcpy(filename,argv[4]);
	iter = atoi(argv[5]);
	length = new int[m];//length of each row
	matrix = new double*[m];

	test03 (m,n,k,iter,filename); 
	/*
 	// Begin another thread to run test, the main will wait for signal 
  	pthread_t test_thread;
  	
  	struct test_argv my_thread;
  	my_thread.m = m;
  	my_thread.n = n;
  	my_thread.k = k;
  	my_thread.iter = iter;
  	my_thread.filename = filename;
  	
  	pthread_create(&test_thread, NULL,test01, &my_thread);
  	
  	char c;		// to get key stroke || 'p' will print, 'k' will kill the loop
  	int loop = 1;		// to control loop
  	
  	
  	while(loop)
  	{
  		if(finish==1)
  		{	
  			cout << finish << endl;
  			loop = 0;
  		}
  		c = getchar();
  		if ( c == 'p')
  		{
  			print_flag = 1;
  			//print out
  		}
  		else if ( c == 'k')
  		{
  			loop = 0;
  		}
  	}
  	*/
  
//	test01 (m,n,k,filename,iter);  

  
	timestamp ( );
	return 0;
}
//****************************************************************************80
void * test02 (void* argv)
{
	
	struct test_argv * my_args;
	my_args = (struct test_argv *)argv;
	int m = my_args -> m;
	int n = my_args -> n;
	int k = my_args -> k;
	int iter = my_args -> iter;
	char * filename = my_args -> filename;
	cout << "m: " << m << endl;
	cout << "n: " << n << endl;
	cout << "k: " << k << endl;
	cout << "iter: " << m << endl;
	cout << "filename: " << filename << endl;
	
	int j;
	while(1){
		if(print_flag){
			for (j=0;j<10;j++)
			{
				cout << "print" << endl;
			}
			print_flag = 0;
		}
	}
	return NULL;
	
}

void * test01 (void* argv)

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tries out the ASA136 routine.
//
//
{
	struct test_argv * my_args;
	my_args = (struct test_argv *)argv;
	int m = my_args -> m;
	int n = my_args -> n;
	int k = my_args -> k;
	int iter = my_args -> iter;
	char * filename = my_args -> filename;

	int i,j,ifault,nw,nd,len,word;
	int *ic1,*nc;

	string line;

	double *wss;
	double value;
	ic1 = new int[m];
	nc = new int[k];
	wss = new double[k];

	cout << "********************************************\n";
	cout << "TEST01\n";
	cout << "  Test the KMNS algorithm,\n";
	cout << "  Applied Statistics Algorithm #136.\n";
	cout << "********************************************\n";

	
//
//  Read the data.
//
	
	FILE *fp;
	nw=0;nd=0;

  	printf("reading data from %s\n", filename);
	fp = fopen(filename,"r");
	if(!fp){
		printf("Unable to open file.\n");
		return NULL;
	}
	matrix = (double**)malloc(sizeof(double*)*m);
	length = (int*)malloc(sizeof(int)*m);
	while((fscanf(fp,"%d",&len)!=EOF))
	{

	    matrix[nd] = (double*)malloc(sizeof(double)*len*2);
		length[nd] = len*2;

//		cout << "***test*** "<< nd << " " << len << endl;

		
		for (j=0;j<len*2;j=j+2)
		{
//			cout << j << " " << word << " " << value <<endl;
		    fscanf(fp,"%d:%lf",&word,&value);
			matrix[nd][j] = (double)word;
			matrix[nd][j+1] = (double)value;
		}
		
		nd++;
	}
	fclose(fp);
	cout << "Data read." << endl;
	
/*	
	cout << "\nHere's a few data:" << endl;
	for (i=0;i<5;i++){
		for (j=0;j<10;j++)
			cout << matrix[i][j] << " ";
		cout << endl;
		}
*/	
	
//
//  Get random k points. Shuffle entire m points, and use first k points.
//
	srand(time(NULL));		//initialize random seed
	vector<int> rand;
	for ( i = 0; i < m; i++) rand.push_back(i); // 0 1 2 3 ... m-1
	random_shuffle (rand.begin(),rand.end());

	
//
//  Initialize the cluster centers.
//  Here, we arbitrarily make the first K data points cluster centers.
//
	double** center = new double*[k];
	for ( i = 1; i <= k; i++ )
	{
//		cout << i << endl;
		center[i-1] = new double[n];
		int row = rand[i-1];		//0~m-1
		len = length[row];
		for ( j = 0; j < len; j = j+2 )
		{
			int dim = int(matrix[row][j])-1;
			double value = matrix[row][j+1];
			center[i-1][dim] = value;
		}

		for (int tt = 0; tt < n; tt ++){
			if (center[i-1][tt] < 0){
			center[i-1][tt] = 0.0;
			}
		}
	}

	cout << "Center chosen!" << endl;

//
//  Compute the clusters.
//
	start = clock();
	finish = kmns ( m, n, center, k, ic1, nc, iter, wss, &ifault );
	cout << "finish = " << finish << endl;

	if ( ifault != 0 )
	{
		cout << "\n";
		cout << "TEST01 - Fatal error!\n";
		cout << "  KMNS returned IFAULT = " << ifault << "\n";
		return NULL;
	}

	/*
	cout << "\n";
	cout << "  Cluster  Population  Energy\n";
	cout << "\n";

	nc_sum = 0;
	wss_sum = 0.0;

	for ( i = 1; i <= k; i++ )
	{
		cout << "  " << setw(8) << i
				<< "  " << setw(8) << nc[i-1]
				<< "  " << setw(14) << wss[i-1] << "\n";
		nc_sum = nc_sum + nc[i-1];
		wss_sum = wss_sum + wss[i-1];
	}

	cout << "\n";
	cout << "     Total"
		<< "  " << setw(8) << nc_sum
		<< "  " << setw(14) << wss_sum << "\n";
	*/

	delete [] ic1;
	delete [] nc;
	delete [] wss;
    
    /*
	//print center
	for( i= 1; i <= k; i++)
	{
		cout << "==============center " << i << "==============" << endl;
		for (j = 1; j <+ n; j++)
		{
  			cout << center[i-1][j-1] << " ";
  			if (j%100 == 0)
  			cout << "\n";
		}
		cout<<endl;
	}
    */
  
	for (j = 0; j < k; j++)
		delete [] center[j];
	delete [] center;

	return NULL;
}

void test03 (int m, int n, int k, int iter, char* filename)
{


	int i,j,ifault,nw,nd,len,word;
	int *ic1,*nc;

	string line;

	double *wss;
	double value;
	ic1 = new int[m];
	nc = new int[k];
	wss = new double[k];

	cout << "********************************************\n";
	cout << "TEST01\n";
	cout << "  Test the KMNS algorithm,\n";
	cout << "  Applied Statistics Algorithm #136.\n";
	cout << "********************************************\n";

	
//
//  Read the data.
//
	
	FILE *fp;
	nw=0;nd=0;

  	printf("reading data from %s\n", filename);
	fp = fopen(filename,"r");
	if(!fp){
		printf("Unable to open file.\n");
		return;
	}
	matrix = (double**)malloc(sizeof(double*)*m);
	length = (int*)malloc(sizeof(int)*m);
	while((fscanf(fp,"%d",&len)!=EOF))
	{

	    matrix[nd] = (double*)malloc(sizeof(double)*len*2);
		length[nd] = len*2;


		
		for (j=0;j<len*2;j=j+2)
		{
//			cout << j << " " << word << " " << value <<endl;
		    fscanf(fp,"%d:%lf",&word,&value);
			matrix[nd][j] = (double)word;
			matrix[nd][j+1] = (double)value;
		}
		
		nd++;
	}
	fclose(fp);
	cout << "Data read." << endl;
	
	
//
//  Get random k points. Shuffle entire m points, and use first k points.
//
	srand(time(NULL));		//initialize random seed
	vector<int> rand;
	for ( i = 0; i < m; i++) rand.push_back(i); // 0 1 2 3 ... m-1
	random_shuffle (rand.begin(),rand.end());

	
//
//  Initialize the cluster centers.
//  Here, we arbitrarily make the first K data points cluster centers.
//
	double** center = new double*[k];
	for ( i = 1; i <= k; i++ )
	{
//		cout << i << endl;
		center[i-1] = new double[n];
		int row = rand[i-1];		//0~m-1
		len = length[row];
		for ( j = 0; j < len; j = j+2 )
		{
			int dim = int(matrix[row][j])-1;
			double value = matrix[row][j+1];
			center[i-1][dim] = value;
		}

		for (int tt = 0; tt < n; tt ++){
			if (center[i-1][tt] < 0){
			center[i-1][tt] = 0.0;
			}
		}
	}

	cout << "Center chosen!" << endl;

//
//  Compute the clusters.
//
	start = clock();
	finish = kmns ( m, n, center, k, ic1, nc, iter, wss, &ifault );
	cout << "finish = " << finish << endl;

	if ( ifault != 0 )
	{
		cout << "\n";
		cout << "TEST01 - Fatal error!\n";
		cout << "  KMNS returned IFAULT = " << ifault << "\n";
		return;
	}


	delete [] ic1;
	delete [] nc;
	delete [] wss;

	for (j = 0; j < k; j++)
		delete [] center[j];
	delete [] center;

	return;
}
