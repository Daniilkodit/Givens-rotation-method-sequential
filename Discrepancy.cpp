#include "header.h"
double Discrepancy(double*A, double*X, int n,int m,double* block,double*blockE,double* Res,double*rmax, MPI_Comm com, double* buf /*???*/)
{
        int i, j, k =n/m, l = n-m*k,t;
        double r0 = 0;
        if(n > 11000) return 0.0;
	for(t = 0;t<k;t++)
	{
		for(i=0;i<m;i++)
                	rmax[i] = -1.0;

		for(i=0;i<k;i++)
		{

			j = 0;
                        get_block(A,block,n,m,i,j);
                        get_block(X,blockE,n,m,j,t);
                        Block_multiplication(Res,block,blockE,m,m,m,m);
			for(j=1;j<k;j++)
			{
			
				get_block(A,block,n,m,i,j);
				get_block(X,blockE,n,m,j,t);
				Block_mult_add(Res,block,blockE,m,m,m,m);
			}
			if(l!=0)
			{
				j = k;
				get_block(A,block,n,m,i,j);
                        	get_block(X,blockE,n,m,j,t);
                        	Block_mult_add(Res,block,blockE,m,m,l,m);
			}
			Disr_block(Res,rmax,m,m,m);
		}
		if(l!=0)
		{
			i = k;
			j = 0;
			get_block(A,block,n,m,i,j);
                        get_block(X,blockE,n,m,j,t);
                        Block_multiplication(Res,block,blockE,l,m,m,m);
			for(j=1;j<k;j++)
                        {

                                get_block(A,block,n,m,i,j);
                                get_block(X,blockE,n,m,j,t);
                                Block_mult_add(Res,block,blockE,l,m,m,m);
                        }
                        if(l!=0)
                        {
                                j = k;
                                get_block(A,block,n,m,i,j);
                                get_block(X,blockE,n,m,j,t);
                                Block_mult_add(Res,block,blockE,l,m,l,m);
                        }
                        Disr_block(Res,rmax,m,m,l);
		}
		for(int w  =0;w<m;w++)
			if(r0<rmax[w]) r0 = rmax[w]; 
	}
	if(l!=0)
	{
		t=k;
		for(i=0;i<m;i++)
                	rmax[i] = -1.0;

		for(i=0;i<k;i++)
		{
			j = 0;
                        get_block(A,block,n,m,i,j);
                        get_block(X,blockE,n,m,j,t);
                        Block_multiplication(Res,block,blockE,m,l,m,m);
			for(j=1;j<k;j++)
			{
			
				get_block(A,block,n,m,i,j);
				get_block(X,blockE,n,m,j,t);
				Block_mult_add(Res,block,blockE,m,l,m,m);
			}
			if(l!=0)
			{
				j = k;
				get_block(A,block,n,m,i,j);
                        	get_block(X,blockE,n,m,j,t);
                        	Block_mult_add(Res,block,blockE,m,l,l,m);
			}
			Disr_block(Res,rmax,m,l,m);
		}
		if(l!=0)
		{
			i = k;
			j = 0;
                        get_block(A,block,n,m,i,j);
                        get_block(X,blockE,n,m,j,t);
                        Block_multiplication(Res,block,blockE,l,l,m,m);
			for(j=1;j<k;j++)
                        {

                                get_block(A,block,n,m,i,j);
                                get_block(X,blockE,n,m,j,t);
                                Block_mult_add(Res,block,blockE,l,l,m,m);
                        }
                        if(l!=0)
                        {
                                j = k;
                                get_block(A,block,n,m,i,j);
                                get_block(X,blockE,n,m,j,t);
                                Block_mult_add(Res,block,blockE,l,l,l,m);
                        }
                        Disr_block(Res,rmax,m,l,l);
		}
		for(int w  =0;w<m;w++)
			if(r0<rmax[w]) r0 = rmax[w];
	}
        return r0;
}
double Norm_matrix(double*A,int n)
{
	int col,row;
	double r0=0,r1=0;
	for( col = 0; col < n; col++)
        {
                for( row = 0; row < n; row++)
                {
			r1 += fabs(A[row*n+col]);
		}
		r0 = r0 > r1 ? r0 : r1;
                r1 = 0;
	}
	return r0;
}
void Disr_block(double*Res,double*rmax,int m,int m1,int m2)
{
	int row,col;
	double r1=0;
	for( col = 0;col<m1;col++)
	{
		for( row = 0;row<m2;row++)
		{
			
			r1+=fabs(Res[row*m+col]);
		}
		rmax[col]+=r1;
		r1 = 0;
	}
}
