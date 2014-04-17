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

void test01 (int,int,int,char*);
int *length;
double **matrix;

double distance(double *,double *, int, int);
//****************************************************************************80

int main(int argc, char* argv[])
//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for ASA136_PRB.
//
//  Discussion:
//
//    ASA136_PRB calls the ASA136 routines.
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
//    John Burkardt
//
{
	timestamp ( );
	if(argc<5)
	{
	cout << "Usage: Kmeans m n k filename iter.max"<<endl;
	return -1;
	}
	int k;
	int m;
	int n;
//	int iter;
	char filename[128];
	m = atoi(argv[1]);
	n = atoi(argv[2]);
	k = atoi(argv[3]);
	memset(filename,0,sizeof(filename));
	strcpy(filename,argv[4]);
//	iter = atoi(argv[5]);
	length = new int[m];//length of each row
	matrix = new double*[m];
  
	test01 (m,n,k,filename);
  
	timestamp ( );
	return 0;
}
//****************************************************************************80

void test01 (int m, int n, int k, char* filename)

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tries out the ASA136 routine.
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
//    John Burkardt
//
{
	int i;
	int *ic1;
	int ifault;
	string line;
	int iter = 50;
	int j;
	int *nc;
	int nc_sum;
	double *wss;
	double wss_sum;
	int len_a;

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
	ifstream input (filename);
	if(!input.is_open()){
		cout<<"Unable to open file."<<endl;
		return;
	}
  

	for(i=0 ; i<m ; ++i){
		int len = 0;
		int iSpace = 0;
		getline(input, line);
		len = line.length();
		for(int l = 0; l< len;l++){
			if(isspace(int(line[l]))){
				iSpace++;
			}
		} 
		len = iSpace;       // len of row. Note: This applies for cases where each line end up with a space.
		length[i] = len;
		matrix[i] = new double[len];
		// split line and read into matrix
		istringstream iss(line);
		for(j=0 ; j<n; ++j){           
			iss >> matrix[i][j];
		}
	}
	
	input.close ( );
	cout<<"Data read into memory."<<endl;
//
//  Get random k points. Shuffle entire m points, and use first k points.
//
	srand(time(NULL));		//initialize random seed
	vector<int> rand;
	for ( i = 0; i < m; i++) rand.push_back(i); // 0 1 2 3 ... m-1
//	random_shuffle (rand.begin(),rand.end());

	
//
//  Initialize the cluster centers.
//  Here, we arbitrarily make the first K data points cluster centers.
//
	double** center = new double*[k];
	for ( i = 1; i <= k; i++ )
	{
		center[i-1] = new double[n];
		int row = rand[i-1];		//0~m-1
		len_a = length[row];
		for ( j = 0; j < len_a; j = j+2 )
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
	cout<<"Random center chosen."<<endl;

//
//  Compute the clusters.
//
	kmns ( m, n, center, k, ic1, nc, iter, wss, &ifault );

	if ( ifault != 0 )
	{
		cout << "\n";
		cout << "TEST01 - Fatal error!\n";
		cout << "  KMNS returned IFAULT = " << ifault << "\n";
		return;
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

	return;
}

