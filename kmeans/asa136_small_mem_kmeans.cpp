//# include "stdafx.h"
# include "asa136_small_mem_kmeans.h"
# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <stdio.h>
# include <fstream>
# include <iostream>
# include <sstream>
# include <string>
# include <vector>
# include <string.h>
# include <cstdlib>

using namespace std;

extern int *length;
extern double **matrix;
extern int print_flag;
extern clock_t start,end;

int kmns ( int m, int n, double *center[], int k, int ic1[], int nc[], 
	int iter, double wss[], int *ifault )

//
//****************************************************************************80
//
//  Purpose:
//
//    KMNS carries out the K-means algorithm.
//
//    Output, int *IFAULT, error indicator.
//    0, no error was detected.
//    1, at least one cluster is empty after the initial assignment.  A better
//       set of initial cluster centers is needed.
//    2, the allowed maximum number off iterations was exceeded.
//    3, K is less than or equal to 1, or greater than or equal to M.
//
{
	double aa;
	double *an1;
	double *an2;
	double *d;
	double da;
	double db;
	double dc;
	double dt[2];
	int i;
	int *ic2;
	int ii;
	int ij;
	int il;
	int indx;
	int *itran;
	int j;
	int l;
	int *live;
	int *ncp;
	double temp;
	int len_a;
	int u,t;

	cout << "begin kmeans" << endl;

	*ifault = 0;

	if ( k <= 1 || m <= k ){
		*ifault = 3;
		return 1;
	}
	ic2 = new int[m];
	an1 = new double[k];
	an2 = new double[k];
	ncp = new int[k];
	d = new double[m];
	itran = new int[k];
	live = new int[k];
  
//
//  For each point I, find its two closest centers, IC1(I) and
//  IC2(I).  Assign the point to IC1(I).
//
	for ( i = 1; i <= m; i++ ){

		ic1[i-1] = 1;
		ic2[i-1] = 2;

		for  (il = 1; il <= 2; il++ ){
			dt[il-1] = 0.0;
			len_a = length[i-1];
			dt[il-1] = distance(matrix[i-1],center[il-1],len_a,n);
		}


		if ( dt[1] < dt[0] ){
			ic1[i-1] = 2;
			ic2[i-1] = 1;
			temp = dt[0];
			dt[0] = dt[1];
			dt[1] = temp;
		}

		for ( l = 3; l <= k; l++ ){
			db = 0.0;
			len_a = length[i-1];
			db = distance(matrix[i-1],center[l-1],len_a,n);

			if ( db < dt[1] ){
				if ( dt[0] <= db ){
					dt[1] = db;
					ic2[i-1] = l;
				}
				else{
					dt[1] = dt[0];
					ic2[i-1] = ic1[i-1];
					dt[0] = db;
					ic1[i-1] = l;
				}
			}
		}
	}


//
//  Update cluster centers to be the average of points contained within them.
//
	memset(nc, 0, sizeof(int)*k);
	for (i = 1; i <= k; i++){
		memset(center[i-1],0, sizeof(double)*n);
	}


	for ( i = 1; i <= m; i++ ){
		l = ic1[i-1];
		nc[l-1] = nc[l-1] + 1;  // l=1~k: each cluster center
		len_a = length[i-1];
		for ( j = 0; j < len_a; j = j+2 ){
			int dim = int(matrix[i-1][j])-1;
			double value = matrix[i-1][j+1];
			center[l-1][dim] = center[l-1][dim] + value;
		}
	}


//
//  Check to see if there is any empty cluster at this stage.
//
	/*
	*ifault = 1;
	for ( l = 1; l <= k; l++ ){
		if ( nc[l-1] == 0 ){
			*ifault = 1;
			return;
		}
	}
	*/
	*ifault = 0;

	cout << "begin transfer" << endl;

	for ( l = 1; l <= k; l++ ){
		aa = ( double ) ( nc[l-1] );
		for ( j = 1; j <= n; j++ ){
			center[l-1][j-1] = center[l-1][j-1] / aa;
		}

//  Initialize AN1, AN2, ITRAN and NCP.
//
//  AN1(L) = NC(L) / (NC(L) - 1)
//  AN2(L) = NC(L) / (NC(L) + 1)
//  ITRAN(L) = 1 if cluster L is updated in the quick-transfer stage,
//           = 0 otherwise
//
//  In the optimal-transfer stage, NCP(L) stores the step at which
//  cluster L is last updated.
//
//  In the quick-transfer stage, NCP(L) stores the step at which
//  cluster L is last updated plus M.
//
		an2[l-1] = aa / ( aa + 1.0 );

		if ( 1.0 < aa ){
			an1[l-1] = aa / ( aa - 1.0 );
		}
		else{
			an1[l-1] = r8_huge ( );
		}
		itran[l-1] = 1;
		ncp[l-1] = -1;
	}
  
//
// transfer
//
	indx = 0;
	*ifault = 2;

	for ( ij = 1; ij <= iter; ij++ ){
//
//  In this stage, there is only one pass through the data.   Each
//  point is re-allocated, if necessary, to the cluster that will
//  induce the maximum reduction in within-cluster sum of squares.
//
		optra ( m, n, center, k, ic1, ic2, nc, an1, an2, ncp, d, itran, live, &indx );	
		cout << "optra finished" << endl;
//
//  Stop if no transfer took place in the last M optimal transfer steps.
//
		if ( indx == m ){
			*ifault = 0;
			break;
		}
//
//  Each point is tested in turn to see if it should be re-allocated
//  to the cluster to which it is most likely to be transferred,
//  IC2(I), from its present cluster, IC1(I).   Loop through the
//  data until no further change is to take place.
//
		qtran ( m, n, center, k, ic1, ic2, nc, an1, an2, ncp, d, itran, &indx );
		cout << "qtran finished" << endl;

//
//  If there are only two clusters, there is no need to re-enter the
//  optimal transfer stage.
//
		if ( k == 2 ){
			*ifault = 0;
			break;
		}
//
//  NCP has to be set to 0 before entering OPTRA.
//
		for ( l = 1; l <= k; l++ ){
			ncp[l-1] = 0;
		}
		
		timestamp ( );
		
		if(print_flag){
			end = clock();
			double elapsed = ((double)(end-start));
			elapsed = (int)(elapsed);
			
			
			cout << "printing..." << endl;
			double * dist = new double[m];	//dist will store distance of each point to its center.
			double *tmp = new double[n]; 
			for (i = 1; i <= m; i++ ){
				ii = ic1[i-1];
				// tmp is to store and calculate difference in center[ii-1] and matrix[i-1]

				for (int t = 0; t < n; t++){
					tmp[t] = center[ii-1][t];
				}

				// for each dim in matrix[i-1], deduct it from center[ii-1]
				len_a = length[i-1];
				for ( j = 0; j < len_a; j = j+2 ){
					int dim = int(matrix[i-1][j])-1;
					double value = matrix[i-1][j+1];
					tmp[dim] = value - tmp[dim];
				}
				for(int t = 0; t < n; t++){
					wss[ii-1] = wss[ii-1] + tmp[t]*tmp[t];
					dist[i-1] = dist[i-1] + tmp[t]*tmp[t];
				}   

			}
			delete [] tmp;
		
			// Output cluster distribution to file.
			// Output distance of each point to its center.
//			string cmd = strcat("mkdir ",atio(elasped))
//			mkdir(cmd.c_str(),0750);
			ofstream output1;
//			string cluster = "cluster." + Str(elapsed) + ".txt";
			output1.open("cluster");	
			ofstream output2;
//			string distance = strncat("distance.", elapsed);
			output2.open("distance");//.c_str());
			ofstream output3;			
//			string distri = strncat("cluster_distribution.", elapsed);
			output3.open("distri");//.c_str());
			for(i = 1 ; i <= m ; i++){
				output1 <<  ic1[i-1] << " " << endl;
				output2 << dist[i-1] << " " << endl;
			}
			cout << "\nCluster saved into 'cluster.txt'." << endl;
			cout << "\nDistance saved into 'distance.txt'." << endl;	
			for(u = 0 ; u < k ; u++){
				output3 << i << ": ";
				for(t = 1 ; t <= m ; t++){
					if(ic1[t-1] == u){
						output3 << ic1[i-1] << " ";
					}
				} 
				output3 << endl;
			}	
			cout << "\nCluster distribution saved into 'cluster_distribution.txt'." << endl;
			output1.close();
			output2.close();
			output3.close();
			
			print_flag = 0;
		}
		
	}
//
//  If the maximum number of iterations was taken without convergence,
//  IFAULT is 2 now.  This may indicate unforeseen looping.
//
	if ( *ifault == 2 )
	{
		cout << "\n";
		cout << "KMNS - Warning!\n";
		cout << "  Maximum number of iterations reached\n";
		cout << "  without convergence.\n";
	}
//
//  Compute the within-cluster sum of squares for each cluster.
//
	memset(wss, 0, sizeof(double)*k);
	for (i = 1; i <= k; i++){
		memset(center[i-1],0, sizeof(double)*n);
	}

	for ( i = 1; i <= m; i++ ){
		ii = ic1[i-1];
		len_a = length[i-1];
		for ( j = 0; j < len_a; j = j+2 ){
			int dim = int(matrix[i-1][j])-1;
			double value = matrix[i-1][j+1];
			center[ii-1][dim] = center[ii-1][dim] + value;
		}
	}

	for ( j = 1; j <= n; j++ ){
		for ( l = 1; l <= k; l++ ){
			center[l-1][j-1] = center[l-1][j-1] / ( double ) ( nc[l-1] );
		}
	}


	double * dist = new double[m];	//dist will store distance of each point to its center.
	double *tmp = new double[n]; 
	for (i = 1; i <= m; i++ ){
		ii = ic1[i-1];
		// tmp is to store and calculate difference in center[ii-1] and matrix[i-1]

		for (int t = 0; t < n; t++){
			tmp[t] = center[ii-1][t];
		}

		// for each dim in matrix[i-1], deduct it from center[ii-1]
		len_a = length[i-1];
		for ( j = 0; j < len_a; j = j+2 ){
			int dim = int(matrix[i-1][j])-1;
			double value = matrix[i-1][j+1];
			tmp[dim] = value - tmp[dim];
		}
		for(int t = 0; t < n; t++){
			wss[ii-1] = wss[ii-1] + tmp[t]*tmp[t];
			dist[i-1] = dist[i-1] + tmp[t]*tmp[t];
		}   

	}
	delete [] tmp; 


// Output cluster distribution to file.
// Output distance of each point to its center.
	ofstream output1;
	output1.open("cluster.txt");	
	ofstream output2;
	output2.open("distance.txt");
	ofstream output3;
	output3.open("cluster_distribution.txt");
	for(i = 1 ; i <= m ; i++){
		output1 <<  ic1[i-1] << " " << endl;
		output2 << dist[i-1] << " " << endl;
	}
	cout << "\nCluster saved into 'cluster.txt'." << endl;
	cout << "\nDistance saved into 'distance.txt'." << endl;
	for(j = 1 ; j <= k ; j++){
		output3 << j << ": ";
		for(i = 1 ; i <= m ; i++){
			if(ic1[i-1] == j){
				output3 << i << " ";
			}
		} 
		output3 << endl;
	}	
	cout << "\nCluster distribution saved into 'cluster_distribution.txt'." << endl;
	output1.close();
	output2.close();
	output3.close();

	

	delete [] ic2;
	delete [] an1;
	delete [] an2;
	delete [] ncp;
	delete [] d;
	delete [] itran;
	delete [] live;

	return 1;
}
//****************************************************************************80

void optra ( int m, int n, double *center[], int k, int ic1[], 
  int ic2[], int nc[], double an1[], double an2[], int ncp[], double d[], 
  int itran[], int live[], int *indx )

//****************************************************************************80
//
//  Purpose:
//
//    OPTRA carries out the optimal transfer stage.
//
{
	double al1;
	double al2;
	double alt;
	double alw;
	double da;
	double db;
	double dc;
	double dd;
	double de;
	double df;
	int i;
	int j;
	int l;
	int l1;
	int l2;
	int ll;
	double r2;
	double rr;
	int len_a;


//
//  If cluster L is updated in the last quick-transfer stage, it
//  belongs to the live set throughout this stage.   Otherwise, at
//  each step, it is not in the live set if it has not been updated
//  in the last M optimal transfer steps.
//
	for ( l = 1; l <= k; l++ )
	{
		if ( itran[l-1] == 1)
		{
			live[l-1] = m + 1;
		}
	}

	for ( i = 1; i <= m; i++ )
	{
		*indx = *indx + 1;
		l1 = ic1[i-1];
		l2 = ic2[i-1];
		ll = l2;
//
//  If point I is the only member of cluster L1, no transfer.
//

		if ( 1 < nc[l1-1]  )
		{
//
//  If L1 has not yet been updated in this stage, no need to
//  re-compute D(I).
//
			if ( ncp[l1-1] != 0 )
			{
				de = 0.0;
				len_a = length[i-1];
				de = distance(matrix[i-1],center[l1-1],len_a,n);
				d[i-1] = de * an1[l1-1];
			}
//
//  Find the cluster with minimum R2.
//
			da = 0.0;
			len_a = length[i-1];
			da = distance(matrix[i-1],center[l2-1],len_a,n);
			r2 = da * an2[l2-1];

			for ( l = 1; l <= k; l++ )
			{

//
//  If LIVE(L1) <= I, then L1 is not in the live set.   If this is
//  true, we only need to consider clusters that are in the live set
//  for possible transfer of point I.   Otherwise, we need to consider
//  all possible clusters.
//
				if ( ( i < live[l1-1] || i < live[l2-1] ) && l != l1 && l != ll )
				{
					rr = r2 / an2[l-1];
					dc = 0.0;
					len_a = length[i-1];
					dc = distance(matrix[i-1],center[l-1],len_a,n);

					if ( dc < rr )
					{
						r2 = dc * an2[l-1];
						l2 = l;
					} 
				}
			}
//
//  If no transfer is necessary, L2 is the new IC2(I).
//
			if ( d[i-1] <= r2 )
			{
				ic2[i-1] = l2;
			}
//
//  Update cluster centers, LIVE, NCP, AN1 and AN2 for clusters L1 and
//  L2, and update IC1(I) and IC2(I).
//
			else
			{
				*indx = 0;
				live[l1-1] = m + i;
				live[l2-1] = m + i;
				ncp[l1-1] = i;
				ncp[l2-1] = i;
				al1 = ( double ) ( nc[l1-1] );
				alw = al1 - 1.0;
				al2 = ( double ) ( nc[l2-1] );
				alt = al2 + 1.0;

				len_a = length[i-1];

				bool *flag = new bool[n];
				memset(flag, 0, n);
				for ( j = 0; j < len_a; j = j+2 )
				{
					int dim = int(matrix[i-1][j])-1;
					flag[dim] = true;
					double value = matrix[i-1][j+1];
					center[l1-1][dim] = (center[l1-1][dim] * al1 - value ) / alw;
					center[l2-1][dim] = (center[l2-1][dim] * al2 + value ) / alt;
				}
// deal with the rest dim in center which is not in matrix[i-1]
				for (j = 0; j < n; j++ ){
					if(!flag[j]){
						center[l1-1][j] = center[l1-1][j] * al1 / alw;
						center[l2-1][j] = center[l2-1][j] * al2 / alt;
					}
				}
				delete [] flag;

				nc[l1-1] = nc[l1-1] - 1;
				nc[l2-1] = nc[l2-1] + 1;
				an2[l1-1] = alw / al1;
				if ( 1.0 < alw )
				{
					an1[l1-1] = alw / ( alw - 1.0 );
				}
				else
				{
					an1[l1-1] = r8_huge ( );
				}
				an1[l2-1] = alt / al2;
				an2[l2-1] = alt / ( alt + 1.0 );
				ic1[i-1] = l2;
				ic2[i-1] = l1;
			}
		}

		if ( *indx == m )
		{
			return;
		}
	}
//
//  ITRAN(L) = 0 before entering QTRAN.   Also, LIVE(L) has to be
//  decreased by M before re-entering OPTRA.
//
	for ( l = 1; l <= k; l++ )
	{
		itran[l-1] = 0;
		live[l-1] = live[l-1] - m;
	}
	return;
}
//****************************************************************************80

void qtran ( int m, int n, double *center[], int k, int ic1[], 
int ic2[], int nc[], double an1[], double an2[], int ncp[], double d[], 
int itran[], int *indx )

//****************************************************************************80
//
//  Purpose:
//
//    QTRAN carries out the quick transfer stage.
//
{
	double al1;
	double al2;
	double alt;
	double alw;
	double da;
	double db;
	double dd;
	double de;
	int i;
	int icoun;
	int istep;
	int j;
	int l1;
	int l2;
	double r2;
	int len_a;
  
//
//  In the optimal transfer stage, NCP(L) indicates the step at which
//  cluster L is last updated.   In the quick transfer stage, NCP(L)
//  is equal to the step at which cluster L is last updated plus M.
//
	icoun = 0;
	istep = 0;

	for ( ; ; )
	{
		for ( i = 1; i <= m; i++ )
		{
			icoun = icoun + 1;
			istep = istep + 1;
			l1 = ic1[i-1];
			l2 = ic2[i-1];
//
//  If point I is the only member of cluster L1, no transfer.
//
			if ( 1 < nc[l1-1] )
			{
//
//  If NCP(L1) < ISTEP, no need to re-compute distance from point I to
//  cluster L1.   Note that if cluster L1 is last updated exactly M
//  steps ago, we still need to compute the distance from point I to
//  cluster L1.
//
				if ( istep <= ncp[l1-1] )
				{
					da = 0.0;
					len_a = length[i-1];
					da = distance(matrix[i-1], center[l1-1],len_a,n);
					d[i-1] = da * an1[l1-1];
				}
//
//  If NCP(L1) <= ISTEP and NCP(L2) <= ISTEP, there will be no transfer of
//  point I at this step.
//
				if ( istep < ncp[l1-1] || istep < ncp[l2-1] )
				{
					r2 = d[i-1] / an2[l2-1];
					dd = 0.0;
					len_a = length[i-1];
					dd = distance(matrix[i-1], center[l2-1],len_a,n);
//
//  Update cluster centers, NCP, NC, ITRAN, AN1 and AN2 for clusters
//  L1 and L2.   Also update IC1(I) and IC2(I).   Note that if any
//  updating occurs in this stage, INDX is set back to 0.
//
					if ( dd < r2 )
					{
						icoun = 0;
						*indx = 0;
						itran[l1-1] = 1;
						itran[l2-1] = 1;
						ncp[l1-1] = istep + m;
						ncp[l2-1] = istep + m;
						al1 = ( double ) ( nc[l1-1] );
						alw = al1 - 1.0;
						al2 = ( double ) ( nc[l2-1] );
						alt = al2 + 1.0;		
						len_a = length[i-1];						
						bool *flag = new bool[n];//temperary store dim in a
						memset(flag, 0, n);
						
						for ( j = 0; j < len_a; j = j+2 )
						{
							int dim = int(matrix[i-1][j])-1;
							flag[dim] = true;
							double value = matrix[i-1][j+1];
							center[l1-1][dim] = (center[l1-1][dim] * al1 - value ) / alw;
							center[l2-1][dim] = (center[l2-1][dim] * al2 + value ) / alt;
						}
						
// deal with the rest dim in center which is not in matrix[i-1]
						for (j = 0; j < n; j++ ){
							if(!flag[j]){    
								center[l1-1][j] = center[l1-1][j] * al1 / alw;
								center[l2-1][j] = center[l2-1][j] * al2 / alt;
							}
						}
						delete [] flag;	
						nc[l1-1] = nc[l1-1] - 1;
						nc[l2-1] = nc[l2-1] + 1;
						an2[l1-1] = alw / al1;
						if ( 1.0 < alw ){
							an1[l1-1] = alw / ( alw - 1.0 );
						}
						else{
							an1[l1-1] = r8_huge ( );
						}
						an1[l2-1] = alt / al2;
						an2[l2-1] = alt / ( alt + 1.0 );
						ic1[i-1] = l2;
						ic2[i-1] = l1;
					}
				}
			}
//
//  If no re-allocation took place in the last M steps, return.
//
			if ( icoun == m )
			{ 
				return;
			}
		}
	}
}
//****************************************************************************80

double r8_huge ( void )

//****************************************************************************80
//
//  Purpose:
//
//    R8_HUGE returns a "huge" R8.
//
{
	double value;

	value = 1.0E+30;

	return value;
}
//****************************************************************************80

void timestamp ( )

//****************************************************************************80
//
//  Purpose:
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//    May 31 2001 09:45:54 AM
//
{
	# define TIME_SIZE 40

	static char time_buffer[TIME_SIZE];
	const struct tm *tm;
	size_t len;
	time_t now;

	now = time ( NULL );
	tm = localtime ( &now );

	len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

	cout << time_buffer << "\n";

	return;
	# undef TIME_SIZE
}

//****************************************************************************80

double distance(double a[],double c[], int len_a, int len_c)
{

//****************************************************************************80
// Calculate distance of a dynamic array to a static array 
    int i;
    int j;


    double *tmp = new double[len_c]; 
    // initiate tmp
    for( i = 0; i < len_c; i++){
        tmp[i] = c[i];
    }
    for ( j = 0; j < len_a; j = j+2){ 
        int dim = int(a[j])-1;
        double value = a[j+1];
        tmp[dim] = tmp[dim] - value;
    }
    
    double dt = 0;
    for (i = 0; i < len_c; i ++){
        dt = dt + tmp[i] * tmp[i];
    }
    
    delete [] tmp;
    return dt;
}
