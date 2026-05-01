#ifndef HEADER
#define HEADER
#include<cstdio>
#include<fstream>
#include<string>
#include<chrono>
#include<cmath>
#include <fenv.h>
int feenableexcept(int excepts);
int fedisableexcept(int excepts);
int fegetexcept(void);
int Ainit(double*, int, int, const std::string&);
void Block_mult_add(double* Result, double* Block_A, double* Block_B, int row_A, int col_B, int col_row, const int m);
int Rotation(int, int,double, double*, double*, double*, double*, double*, double*, double*, double*);
void mul_rotate2(double*block, double*blockE, double*blockC, double*blockCE, int m1, int m2, int j, double*rcos, double*rsin, int l, int k);
void Xmul_rotate2(double*blockE,double*blockCE, int m1, int m2, int j, double*rcos, double*rsin, int l, int k);
int mul_fill_rotate(double*block, double*blockE, int m,int m2, double*rcos, double*rsin, double eps, double Anorm);
void mul_rotate(double*block, double*blockE, int m, int j, double*rcos, double*rsin, int l, int k);
void Xmul_rotate(double*blockE, int m,int m1, int j, double*rcos, double*rsin, int l, int k);
int mul_fill_rotate2(double*block, double*blockE, double*blockC, double*blockCE, int m1, int m2, double*rcos, double*rsin, double eps, double Anorm);
void Print_matrix(int, int, int, double*);
void get_block(double *, double *, int, int, int, int);
void put_block(double *, double *, int, int, int, int);
void Disr_block(double*Res,double*rmax,int m,int m1,int m2);
double Discrepancy(double*A, double*X, int n,int m,double* block,double*blockE,double* Res,double*rmax);
double Norm_matrix(double*,int);
void Block_multiplication(double* Result, double* Block_A, double* Block_B,int row_A,int col_B,int col_row,const int m);
int Inverse(double*X,double*A,int n,int l,double Anorm,double eps);
void Block_mult_diff(double* Result, double* Block_A, double* Block_B, int row_A, int col_B,int col_row,const int m);
#endif
