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

void test01 (int,int,int,char*,int);
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
//
{
	timestamp ( );
	if(argc<5)
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
  
	test01 (m,n,k,filename,iter);
  
	timestamp ( );
	return 0;
}
//****************************************************************************80

void test01 (int m, int n, int k, char* filename, int iter)

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tries out the ASA136 routine.
//
//
{
	int i,j,ifault,nw,nd,len,word;
//	int *ic1,*nc;
	double value;

//	ic1 = new int[m];
//	nc = new int[m];

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
	while((fscanf(fp,"%10d",&len)!=EOF))
	{
	    matrix[nd] = (double*)malloc(sizeof(double)*len*2);
		length[nd] = len*2;
		for (j=0;j<len*2;j=j+2)
		{
		    fscanf(fp,"%10d:%lf",&word,&value);
			matrix[nd][j] = (double)word;
			matrix[nd][j+1] = (double)value;
		}
		nd++;
	}
	fclose(fp);
	
	//test
	printf("Data read.\n");
	/*
	printf("\nHere's a few data:\n");
	for (i=0;i<5;i++){
		for (j=0;j<20;j++)
			printf("%.10f ",matrix[i][j]);
		printf("\n");
		}
	*/	

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
	
	for (i = 0; i < m ; i ++){
		cout << endl;
		for (j = 0 ; j < k ; j ++){
			cout << "[" << i << "][" << j << "]: " << distance(matrix[i],center[j], length[i],n) << "\t";
			}
		}

//
//  Compute the clusters.
//
	kmns ( m, n, center, k, /* ic1, nc, */iter, &ifault );

	if ( ifault != 0 )
	{
		cout << "\n";
		cout << "TEST01 - Fatal error!\n";
		cout << "  KMNS returned IFAULT = " << ifault << "\n";
		return;
	}


//	delete [] ic1;

    
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

