//# include "stdafx.h"
# include "asa136_kmeans_v4.h"
# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <fstream>
# include <ctime>
# include <stdio.h>
# include <fstream>
# include <sstream>
# include <string>
# include <vector>

using namespace std;
extern int *length;
extern double **matrix;

void kmns ( int m, int n, double *center[], int k, int ic1[], int nc[], 
	int iter, double wss[], int *ifault )

//****************************************************************************80
//
//  Purpose:
//
//    KMNS carries out the K-means algorithm.
//
//  Discussion:
//
//    This routine attempts to divide M points in N-dimensional space into 
//    K clusters so that the within cluster sum of squares is minimized.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 February 2008
//
//  Author:
//
//    Original FORTRAN77 version by John Hartigan, Manchek Wong.
//    C++ version by John Burkardt
//
//  Reference:
//
//    John Hartigan, Manchek Wong,
//    Algorithm AS 136:
//    A K-Means Clustering Algorithm,
//    Applied Statistics,
//    Volume 28, Number 1, 1979, pages 100-108.
//
//  Parameters:
//
//    Input, double A(M,N), the points.
//
//    Input, int M, the number of points.
//
//    Input, int N, the number of spatial dimensions.
//
//    Input/output, double C(K,N), the cluster centers.
//
//    Input, int K, the number of clusters.
//
//    Output, int IC1(M), the cluster to which each point 
//    is assigned.
//
//    Output, int NC(K), the number of points in each cluster.
//
//    Input, int ITER, the maximum number of iterations allowed.
//
//    Output, double WSS(K), the within-cluster sum of squares
//    of each cluster.
//
//    Output, int *IFAULT, error indicator.
//    0, no error was detected.
//    1, at least one cluster is empty after the initial assignment.  A better
//       set of initial cluster centers is needed.
//    2, the allowed maximum number off iterations was exceeded.
//    3, K is less than or equal to 1, or greater than or equal to M.
//
{
//  double *a;
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

  *ifault = 0;

  if ( k <= 1 || m <= k )
  {
    *ifault = 3;
    return;
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
  for ( i = 1; i <= m; i++ )
  {

    ic1[i-1] = 1;
    ic2[i-1] = 2;

    for  (il = 1; il <= 2; il++ )
    {

        dt[il-1] = 0.0;
        len_a = length[i-1];
        dt[il-1] = distance(matrix[i-1],center[il-1],len_a,n);
    }


    if ( dt[1] < dt[0] )
    {
      ic1[i-1] = 2;
      ic2[i-1] = 1;
      temp = dt[0];
      dt[0] = dt[1];
      dt[1] = temp;
    }

    for ( l = 3; l <= k; l++ )
    {
      db = 0.0;
      len_a = length[i-1];
      db = distance(matrix[i-1],center[l-1],len_a,n);

      if ( db < dt[1] )
      {
        if ( dt[0] <= db )
        {
          dt[1] = db;
          ic2[i-1] = l;
        }
        else
        {
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
  double **n_center = new double *[k];
  for (i = 1; i <= k; i++){
    nc[i-1] = 0;
    n_center[i-1] = new double[n];
    for ( j = 1; j <= n; j++ ){
      n_center[i-1][j-1] = 0.0;
    }  
  }


  for ( i = 1; i <= m; i++ )
  {
    l = ic1[i-1];
    nc[l-1] = nc[l-1] + 1;  // l=1~k: each cluster center
//    a = matrix[i-1];       // matrix[i-1]: each point
    len_a = length[i-1];
    for ( j = 0; j < len_a; j = j+2 )
    {
      int dim = int(matrix[i-1][j])-1;
      double value = matrix[i-1][j+1];
      n_center[l-1][dim] = n_center[l-1][dim] + value;
    }
  }



//
//  Check to see if there is any empty cluster at this stage.
//
  *ifault = 1;

  for ( l = 1; l <= k; l++ )
  {
    if ( nc[l-1] == 0 )
    {
      *ifault = 1;
      return;
    }

  }

  *ifault = 0;

  for ( l = 1; l <= k; l++ )
  {
    aa = ( double ) ( nc[l-1] );
    for ( j = 1; j <= n; j++ )
    {
      n_center[l-1][j-1] = n_center[l-1][j-1] / aa;
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

    if ( 1.0 < aa )
    {
      an1[l-1] = aa / ( aa - 1.0 );
    }
    else
    {
      an1[l-1] = r8_huge ( );
    }
    itran[l-1] = 1;
    ncp[l-1] = -1;
  }
  
// n_center -> center
    for (i = 0; i < k; i++){
        for (j = 0; j < n; j++){
            center[i][j] = n_center[i][j];
        }    
    }
  
// delete n_center
  for (i = 0; i < k; ++i)
  	delete [] n_center[i];
  delete [] n_center;

  cout << "Initial cluster distribution: ";
  for(i = 0 ; i < m ; i++)
  {
	cout <<  ic1[i] << " ";
  }
  cout << endl;
  
//
// transfer
//
  indx = 0;
  *ifault = 2;

  for ( ij = 1; ij <= iter; ij++ )
  {
  	cout << "iter " << ij << endl;
//
//  In this stage, there is only one pass through the data.   Each
//  point is re-allocated, if necessary, to the cluster that will
//  induce the maximum reduction in within-cluster sum of squares.
//
    optra ( m, n, center, k, ic1, ic2, nc, an1, an2, ncp, d, itran, live, &indx );
		cout << "optra" << endl;
		/*
    cout << "Cluster distribution after optra: ";
		for(i = 0 ; i < m ; i++)
		{
			cout <<  ic1[i] << " ";
		}
		cout << endl;*/
		/*
		//print center
		for( i= 1; i <= k; i++)
		{
			cout << "\n~~~~~~~~~~center " << i << "~~~~~~~~~~\n" << endl;
			for (j = 1; j <+ n; j++)
			{
			  cout << center[i-1][j-1] << " ";
			  if (j%100 == 0)
			  	cout << "\n";
			}
		}*/
//
//  Stop if no transfer took place in the last M optimal transfer steps.
//
    if ( indx == m )
    {
      *ifault = 0;
      cout << "\nProgram Break!\n" << endl;
      break;
    }
//
//  Each point is tested in turn to see if it should be re-allocated
//  to the cluster to which it is most likely to be transferred,
//  IC2(I), from its present cluster, IC1(I).   Loop through the
//  data until no further change is to take place.
//
    qtran ( m, n, center, k, ic1, ic2, nc, an1, an2, ncp, d, itran, &indx );
    cout << "qtran" << endl;
    /*
    cout << "Cluster distribution after qtran: ";
		for(i = 0 ; i < m ; i++)
		{
			cout <<  ic1[i] << " ";
		}
		cout << endl;*/
		/*
		//print center
		for( i= 1; i <= k; i++)
		{
			cout << "\n~~~~~~~~~~center " << i << "~~~~~~~~~~\n" << endl;
			for (j = 1; j <+ n; j++)
			{
			  cout << center[i-1][j-1] << " ";
			  if (j%100 == 0)
			  	cout << "\n";
			}
		}*/
//
//  If there are only two clusters, there is no need to re-enter the
//  optimal transfer stage.
//
    if ( k == 2 )
    {
      *ifault = 0;
      break;
    }
//
//  NCP has to be set to 0 before entering OPTRA.
//
    for ( l = 1; l <= k; l++ )
    {
      ncp[l-1] = 0;
    }

  cout << "Final cluster distribution: ";
  for(i = 0 ; i < m ; i++)
  {
	cout <<  ic1[i] << " ";
  }
  cout << endl;

  timestamp ( );
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
  cout << "compute within-cluster sum of squares" << endl;
  for ( l = 1; l <= k; l++ )
  {
    wss[l-1] = 0.0;
    for ( j = 1; j <= n; j++ )
    {
      center[l-1][j-1] = 0.0;
    }
  }

  cout << "Initialized." << endl;
  for ( i = 1; i <= m; i++ )
  {
    ii = ic1[i-1];
//    a = matrix[i-1];
    len_a = length[i-1];
    for ( j = 0; j < len_a; j = j+2 )
    {
      int dim = int(matrix[i-1][j])-1;
      double value = matrix[i-1][j+1];
      center[ii-1][dim] = center[ii-1][dim] + value;
    }
  }

  cout << "test" << endl;
  for ( j = 1; j <= n; j++ )
  {
    for ( l = 1; l <= k; l++ )
    {
      center[l-1][j-1] = center[l-1][j-1] / ( double ) ( nc[l-1] );
    }
  }
  
	for (i = 1; i <= m; i++ )
	{
	  ii = ic1[i-1];
	  // tmp is to store and calculate difference in center[ii-1] and matrix[i-1]
	  double *tmp = new double[n]; 
	  for (int t = 0; t < n; t++){
			tmp[t] = center[ii-1][t];
	  }

	  // for each dim in matrix[i-1], deduct it from center[ii-1]
	  len_a = length[i-1];
	  for ( j = 0; j < len_a; j = j+2 )
	  {
			int dim = int(matrix[i-1][j])-1;
			double value = matrix[i-1][j+1];
			tmp[dim] = value - center[ii-1][j-1];
	  }
	  for(int t = 0; t < n; t++){
			wss[ii-1] = wss[ii-1] + tmp[t]*tmp[t];
	  }   
	  delete [] tmp; 
	}


  cout << "Final cluster distribution: ";
  for(i = 0 ; i < m ; i++)
  {
    cout <<  ic1[i] << " ";
  }
  cout << endl;

  delete [] ic2;
  delete [] an1;
  delete [] an2;
  delete [] ncp;
  delete [] d;
  delete [] itran;
  delete [] live;

  return;
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
//  Discussion:
//
//    This is the optimal transfer stage.
//
//    Each point is re-allocated, if necessary, to the cluster that
//    will induce a maximum reduction in the within-cluster sum of
//    squares.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 February 2008
//
//  Author:
//
//    Original FORTRAN77 version by John Hartigan, Manchek Wong.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    John Hartigan, Manchek Wong,
//    Algorithm AS 136:
//    A K-Means Clustering Algorithm,
//    Applied Statistics,
//    Volume 28, Number 1, 1979, pages 100-108.
//
//  Parameters:
//
//    Input, double A(M,N), the points.
//
//    Input, int M, the number of points.
//
//    Input, int N, the number of spatial dimensions.
//
//    Input/output, double C(K,N), the cluster centers.
//
//    Input, int K, the number of clusters.
//
//    Input/output, int IC1(M), the cluster to which each 
//    point is assigned.
//
//    Input/output, int IC2(M), used to store the cluster 
//    which each point is most likely to be transferred to at each step.
//
//    Input/output, int NC(K), the number of points in 
//    each cluster.
//
//    Input/output, double AN1(K).
//
//    Input/output, double AN2(K).
//
//    Input/output, int NCP(K).
//
//    Input/output, double D(M).
//
//    Input/output, int ITRAN(K).
//
//    Input/output, int LIVE(K).
//
//    Input/output, int *INDX, the number of steps since a 
//    transfer took place.
//
{
//  double *a;
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
          //cout << "rr" << rr << "\n";

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

//        a = matrix[i-1];
        len_a = length[i-1];
        int *tmp = new int[len_a/2];//temperary store dim in a
        int tt = 0;
        for ( j = 0; j < len_a; j = j+2 )
        {
            int dim = int(matrix[i-1][j])-1;
            tmp[tt] = dim;
            tt++;
            double value = matrix[i-1][j+1];
            center[l1-1][dim] = (center[l1-1][dim] * al1 - value ) / alw;
            center[l2-1][dim] = (center[l2-1][dim] * al2 + value ) / alt;
        }
        // deal with the rest dim in center which is not in matrix[i-1]
        for (j = 0; j < n; j++ ){
            bool flag = true; 
            for (tt = 0; tt < len_a/2; tt++){
                if (j == tmp[tt]){
                    flag = false;
                    break;
                }
            }
            if(flag){    
                center[l1-1][j] = center[l1-1][j] * al1 / alw;
                center[l2-1][j] = center[l2-1][j] * al2 / alw;
            }
        }
        delete [] tmp;

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
//  Discussion:
//
//    This is the quick transfer stage.
//
//    IC1(I) is the cluster which point I belongs to.
//    IC2(I) is the cluster which point I is most likely to be
//    transferred to.
//
//    For each point I, IC1(I) and IC2(I) are switched, if necessary, to
//    reduce within-cluster sum of squares.  The cluster centers are
//    updated after each step.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 February 2008
//
//  Author:
//
//    Original FORTRAN77 version by John Hartigan, Manchek Wong.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    John Hartigan, Manchek Wong,
//    Algorithm AS 136:
//    A K-Means Clustering Algorithm,
//    Applied Statistics,
//    Volume 28, Number 1, 1979, pages 100-108.
//
//  Parameters:
//
//    Input, double A(M,N), the points.
//
//    Input, int M, the number of points.
//
//    Input, int N, the number of spatial dimensions.
//
//    Input/output, double C(K,N), the cluster centers.
//
//    Input, int K, the number of clusters.
//
//    Input/output, int IC1(M), the cluster to which each 
//    point is assigned.
//
//    Input/output, int IC2(M), used to store the cluster 
//    which each point is most likely to be transferred to at each step.
//
//    Input/output, int NC(K), the number of points in 
//    each cluster.
//
//    Input/output, double AN1(K).
//
//    Input/output, double AN2(K).
//
//    Input/output, int NCP(K).
//
//    Input/output, double D(M).
//
//    Input/output, int ITRAN(K).
//
//    Input/output, int INDX, counts the number of steps 
//    since the last transfer.
//
{
//  double *a;
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
			cout << "~~~~~~~~~~" << i << "~~~~~~~~~~" << endl;
      icoun = icoun + 1;
      istep = istep + 1;
      l1 = ic1[i-1];
      l2 = ic2[i-1];
//
//  If point I is the only member of cluster L1, no transfer.
//
      if ( 1 < nc[l1-1] )
      {
      	cout << "in the if" << endl;
//
//  If NCP(L1) < ISTEP, no need to re-compute distance from point I to
//  cluster L1.   Note that if cluster L1 is last updated exactly M
//  steps ago, we still need to compute the distance from point I to
//  cluster L1.
//
        if ( istep <= ncp[l1-1] )
        {
         	cout << "test1" << endl;
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
         	cout << "test2" << endl;
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
          	cout << "dd<r2" << endl;
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

         		cout << "test3-1" << endl;					
						len_a = length[i-1];
						int *tmp = new int[len_a/2];//temperary store dim in a
						int tt = 0;
						for ( j = 0; j < len_a; j = j+2 )
						{
						    int dim = int(matrix[i-1][j])-1;
						    tmp[tt] = dim;
						    tt++;
						    double value = matrix[i-1][j+1];
						    center[l1-1][dim] = (center[l1-1][dim] * al1 - value ) / alw;
						    center[l2-1][dim] = (center[l2-1][dim] * al2 + value ) / alt;
						}
						
						// deal with the rest dim in center which is not in matrix[i-1]
						for (j = 0; j < n; j++ ){
						    bool flag = true; 
						    for (tt = 0; tt < len_a/2; tt++){
						        if (j == tmp[tt]){
						            flag = false;
						            break;
						        }
						    }
						    if(flag){    
						        center[l1-1][j] = center[l1-1][j] * al1 / alw;
						        center[l2-1][j] = center[l2-1][j] * al2 / alw;
						    }
						}
						delete [] tmp;
         		cout << "test3-2" << endl;	
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
            cout << "test3-3" << endl;	
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
			cout << "icoun = " << icoun << endl;
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
//  Discussion:
//
//    The value returned by this function is NOT required to be the
//    maximum representable R8.  This value varies from machine to machine,
//    from compiler to compiler, and may cause problems when being printed.
//    We simply want a "very large" but non-infinite number.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 October 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_HUGE, a "huge" R8 value.
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
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    May 31 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
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
