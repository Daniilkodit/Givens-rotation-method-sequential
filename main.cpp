#include "header.h"
int main(int argc, char*argv[])
{
        int n, m, r, s, i, out = 0, task = 22;
        std::string name;
	double eps = 1.e-15;
        double*A = nullptr;
        double*X = nullptr;
        double*block = nullptr;
        double*rcos = nullptr;
        double*rsin = nullptr;
        double*blockE = nullptr;
        double*blockC = nullptr;
        double*blockCE = nullptr;
	double *rmax = nullptr;
        double t1 = 0, t2=0, r1 = -1, r2 =-1;
	feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);
        if(!((argc == 5 || argc == 6) &&
                        sscanf(argv[1], "%d", &n) == 1 &&
                        sscanf(argv[2], "%d", &m) == 1 &&
                        sscanf(argv[3], "%d", &r) == 1 &&
                        sscanf(argv[4], "%d", &s) == 1))
        {
                printf("Usage: %s n m r s\n", argv[0]);
                return 0;
        }
        if(n <= 0 || m <= 0 || r <= 0 || m > n)
        {
                printf("Usage: n>0 m>0 r>0 s and m<=n\n");
                return 0;
        }
        if(s == 0)
        {
                if(argc != 6)
                {
                        printf("Usage: %s n m r 0 namefile.txt\n", argv[0]);
                        return 0;
                }
                else name = argv[5];
        }
        A = new double[n * n];
        if(Ainit(A, n, s, name) < 0)
        {
                delete[] A;
                return -1;
        }
        Print_matrix(n, n, r, A);
        X = new double[n * n] {};
        block = new double[m * m];
        blockE = new double[m * m];
        blockCE = new double[m * m];
        blockC = new double[m * m];
        rcos = new double[m * m]; //содержит наши косинусы
        rsin = new double[m * m]; //содержит наши косинусы
	rmax = new double[m];
        for(i = 0; i < n; i++)
        {
                X[i * n + i] = 1;
        }
	if(s==4) eps = 1.e-20;
#if 1
        auto start1 = std::chrono::high_resolution_clock::now();
        out = Rotation(n, m,eps, A, X, rcos, rsin, block, blockE, blockC, blockCE); // Метод вращений
        auto end1 = std::chrono::high_resolution_clock::now();
        t1 = std::chrono::duration<double>(end1 - start1).count();
#endif
        if(out == 0)
        {
                Print_matrix(n, n, r, X);
                Ainit(A, n, s, name);
#if 1
                auto start2 = std::chrono::high_resolution_clock::now();
                r1 = Discrepancy(A, X, n,m,blockCE,blockE,block,rmax);
                auto end2 = std::chrono::high_resolution_clock::now();
                t2 = std::chrono::duration<double>(end2 - start2).count();
                r2 = Discrepancy(A, X, n,m,blockCE,blockE,block,rmax);
#endif
        }
        else
        {
                r1 = -1;
                r2 = -1;
        }
        printf ("%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d\n", argv[0], task, r1, r2, t1, t2, s, n, m);
        delete[] A;
        delete[] X;
        delete[] block;
        delete[] blockC;
        delete[] blockCE;
        delete[] blockE;
        delete[] rcos;
        delete[] rsin;
	delete[] rmax;
        return 0;

}
