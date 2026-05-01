#include"header.h"
int Inverse(double* X, double* A, int n, int l, double Anorm, double eps)
{
    int i, j, t;
    for (i = 0; i < l; i++) {
        for (j = 0; j < l; j++) {
            X[i * n + j] = (i == j) ? 1.0 : 0.0;
        }
    }

    for (t = 0; t < l; t++) {
        for (i = l - 1; i >= 0; i--) {
            double sum = 0.0;
	    
            if (fabs(A[i * n + i]) < eps * Anorm) return -1;

            j = i + 1;
            int remainder = (l - j) % 4;

            for (int k = 0; k < remainder; k++) {
                sum += A[i * n + j] * X[j * n + t];
                j++;
            }

            for (; j < l; j += 4) {
                sum += A[i * n + j] * X[j * n + t];
                sum += A[i * n + j + 1] * X[(j + 1) * n + t];
                sum += A[i * n + j + 2] * X[(j + 2) * n + t];
                sum += A[i * n + j + 3] * X[(j + 3) * n + t];
            }

            X[i * n + t] = (X[i * n + t] - sum) / A[i * n + i];
        }
    }
    return 0;
}
void Xmul_rotate (double *blockE, int m,int m1, int j, double *rcos,
            double *rsin, int l, int k)
{
        double cosk, sink, x, y;
        int p, q, c;
        int z = (j == k) ? l : m;
        if(j==k &&l==0)return;
        for (q = 0; q < m1 - 1; q++)
        {
                for (p = q + 1; p < m1; p++)
                {
                        cosk = rcos[m * q + p];
                        sink = rsin[m * q + p];
                        for (c = 0; c < z; c++)
                        {
                                x = cosk * blockE[q * m + c] - sink * blockE[p * m + c];
                                y = sink * blockE[q * m + c] + cosk * blockE[p * m + c];
                                blockE[q * m + c] = x;
                                blockE[p * m + c] = y;
                        }

                }
        }
}
void Xmul_rotate2 (double *blockE, double *blockCE,
             int m1, int m2, int j, double *rcos, double *rsin, int l, int k)
{
        double cosk, sink, x, y;
        int p, q, c;
        int z = (j == k) ? l : m1;
        if(j==k && l==0) return;
        for (q = 0; q < m1; q++)
        {
                for (p = 0; p < m2; p++)
                {
                        cosk = rcos[m1 * q + p];
                        sink = rsin[m1 * q + p];
                        for (c = 0; c < z; c++)
                        {
                                x = cosk * blockE[q * m1 + c] - sink * blockCE[p * m1 + c];
                                y = sink * blockE[q * m1 + c] + cosk * blockCE[p * m1 + c];
                                blockE[q * m1 + c] = x;
                                blockCE[p * m1 + c] = y;
                        }

                }
        }
}
int Rotation (int n, int m,double eps, double *A, double *X, double *rcos, double *rsin,
              double *block, double *blockE, double *blockC, double *blockCE)
{
        int k = n / m, l = n-k*m;
        int r, j,i,t,c;
        double Anorm = Norm_matrix(A,n);
	if(Anorm<=1.e-64) return -1;
        for (r = 0; r < k; r++)
        {
                j = r;			// первые двe формулы алгоритма, заполняем rotator потом перейдем к j=r+1
                get_block (A, block, n, m, r, j);
                get_block (X, blockE, n, m, r, j);
                if (mul_fill_rotate (block, blockE, m,m, rcos, rsin, eps, Anorm) != 0)
                        return -1;
                put_block (A, block, n, m, r, j);
                put_block (X, blockE, n, m, r, j);
		for(c = 0;c<r;c++)
		{
			get_block(X,blockE,n,m,r,c);
			Xmul_rotate(blockE,m,m,c,rcos,rsin,l,k);
			put_block(X,blockE,n,m,r,c);
		}
                for (j = r + 1; j < k; j++)	// насчитали матрицу rotaror теперь первая формула для j=r+1
                {
                        get_block (A, block, n, m, r, j);
                        get_block (X, blockE, n, m, r, j);
                        mul_rotate (block, blockE, m, j, rcos, rsin, l, k);
                        put_block (A, block, n, m, r, j);
                        put_block (X, blockE, n, m, r, j);
                }
		if(l!=0)
		{
			j=k;
			get_block (A, block, n, m, r, j);
                        get_block (X, blockE, n, m, r, j);
                        mul_rotate (block, blockE, m, j, rcos, rsin, l, k);
                        put_block (A, block, n, m, r, j);
                        put_block (X, blockE, n, m, r, j);
		}
                //формулы 3,4,5,6
                for (i = r + 1; i < k; i++)
                {
                        j = r;		//cначала насчитаем матрицу rotation
                        get_block (A, block, n, m, r, j);
                        get_block (A, blockC, n, m, i, j);
                        get_block (X, blockE, n, m, r, j);
                        get_block (X, blockCE, n, m, i, j);
                        if (mul_fill_rotate2
                                        (block, blockE, blockC, blockCE, m, m, rcos, rsin, eps,
                                         Anorm) != 0)
                                return -1;
		
                        put_block (A, block, n, m, r, j);
                        put_block (A, blockC, n, m, i, j);
                        put_block (X, blockE, n, m, r, j);
                        put_block (X, blockCE, n, m, i, j);
			for(c=0;c<r;c++)
			{
				get_block (X, blockE, n, m, r, c);
                        	get_block (X, blockCE, n, m, i, c);
				Xmul_rotate2(blockE,blockCE,m,m,c,rcos,rsin,l,k);
				put_block (X, blockE, n, m, r, c);
                        	put_block (X, blockCE, n, m, i, c);
			}
                        for (j = r + 1; j < k; j++)	//получили матрицу rotation теперь формула 3
                        {

                                get_block (A, block, n, m, r, j);
                                get_block (A, blockC, n, m, i, j);
                                get_block (X, blockE, n, m, r, j);
                                get_block (X, blockCE, n, m, i, j);
                                mul_rotate2 (block, blockE, blockC, blockCE, m, m, j, rcos,
                                             rsin, l, k);
                                put_block (A, block, n, m, r, j);
                                put_block (A, blockC, n, m, i, j);
                                put_block (X, blockE, n, m, r, j);
                                put_block (X, blockCE, n, m, i, j);
                        }
			if(l!=0)
			{
				j = k;
				get_block (A, block, n, m, r, j);
                                get_block (A, blockC, n, m, i, j);
                                get_block (X, blockE, n, m, r, j);
                                get_block (X, blockCE, n, m, i, j);
                                mul_rotate2 (block, blockE, blockC, blockCE, m, m, j, rcos,
                                             rsin, l, k);
                                put_block (A, block, n, m, r, j);
                                put_block (A, blockC, n, m, i, j);
                                put_block (X, blockE, n, m, r, j);
                                put_block (X, blockCE, n, m, i, j);
			}
                }
		if(l==0) continue;
                i = k;			// реализация формул 5,6
                j = r;			//cначала насчитаем матрицу rotation
                get_block (A, block, n, m, r, j);
                get_block (A, blockC, n, m, i, j);
                get_block (X, blockE, n, m, r, j);
                get_block (X, blockCE, n, m, i, j);
                if (mul_fill_rotate2
                                (block, blockE, blockC, blockCE, m, l, rcos, rsin, eps, Anorm) != 0)
                        return -1;
	
                put_block (A, block, n, m, r, j);
                put_block (A, blockC, n, m, i, j);
                put_block (X, blockE, n, m, r, j);
                put_block (X, blockCE, n, m, i, j);
		for(c=0;c<r;c++)
                {
                       get_block (X, blockE, n, m, r, c);
                       get_block (X, blockCE, n, m, i, c);
                       Xmul_rotate2(blockE,blockCE,m,l,c,rcos,rsin,l,k);
                       put_block (X, blockE, n, m, r, c);
                       put_block (X, blockCE, n, m, i, c);
                }	
                for (j = r + 1; j < k + 1; j++)	//получили матрицу rotation
                {

                        get_block (A, block, n, m, r, j);
                        get_block (A, blockC, n, m, i, j);
                        get_block (X, blockE, n, m, r, j);
                        get_block (X, blockCE, n, m, i, j);
                        mul_rotate2 (block, blockE, blockC, blockCE, m, l, j, rcos, rsin, l,
                                     k);
                        put_block (A, block, n, m, r, j);
                        put_block (A, blockC, n, m, i, j);
                        put_block (X, blockE, n, m, r, j);
                        put_block (X, blockCE, n, m, i, j);
                }

        }
	if(l!=0)
	{
		r = k;			// последний шаг алгоритма с блоком l*l
        	j = k;
        	get_block (A, block, n, m, r, j);
        	get_block (X, blockE, n, m, r, j);
        	if (mul_fill_rotate (block, blockE, m,l, rcos, rsin, eps, Anorm) != 0)
                	return -1;
        	put_block (A, block, n, m, r, j);
        	put_block (X, blockE, n, m, r, j);
		for(c = 0;c<r;c++)
                {
                        get_block(X,blockE,n,m,r,c);
                        Xmul_rotate(blockE,m,l,c,rcos,rsin,l,k);//??? не факт
                        put_block(X,blockE,n,m,r,c);
                }

	}
        //Обратный ход метода Гаусса
#if 0
	printf("Обратный\n");
	for ( t = 0; t < n; t++) {
        for (i = n - 1; i >= 0; i--) {
            double sum = 0.0;

            for (j = i + 1; j < n; j++) {
                sum += A[i * n + j] * X[j * n + t];
            }

            X[i * n + t] = (X[i * n  +t] - sum) / A[i * n + i];
        }
    }
	return 0;
#endif
#if 1
	for(i = 0;i<k;i++)
	{
		get_block(A,block,n,m,i,i);
		if(Inverse(blockCE,block,m,m,Anorm,eps)!=0) return -1;
		for(j = i+1;j<k;j++)
		{
			get_block(A,blockE,n,m,i,j);
			Block_multiplication(block,blockCE,blockE,m,m,m,m);
			put_block(A,block,n,m,i,j);
		}
		if(l!=0)
		{
			j = k;
			get_block(A,blockE,n,m,i,j);
                        Block_multiplication(block,blockCE,blockE,m,l,m,m);
                        put_block(A,block,n,m,i,j);
		}
		for(t = 0;t<k;t++)
		{
			get_block(X,blockE,n,m,i,t);
			Block_multiplication(block,blockCE,blockE,m,m,m,m);
                        put_block(X,block,n,m,i,t);
		}
		if(l!=0)
		{
			t = k;
			get_block(X,blockE,n,m,i,t);
                        Block_multiplication(block,blockCE,blockE,m,l,m,m);
                        put_block(X,block,n,m,i,t);
		}
	}
	if(l!=0)
	{
		i = k;
		get_block(A,block,n,m,i,i);
        	if(Inverse(blockCE,block,m,l,Anorm,eps)!=0) return -1;
		for(t = 0;t<k;t++)
        	{
                  get_block(X,blockE,n,m,i,t);
                  Block_multiplication(block,blockCE,blockE,l,m,l,m);
                  put_block(X,block,n,m,i,t);
        	}
		if(l!=0)
                {
                        t = k;
                        get_block(X,blockE,n,m,i,t);
                        Block_multiplication(block,blockCE,blockE,l,l,l,m);
                        put_block(X,block,n,m,i,t);
                }
	}
	for (i = k - 1; i >= 0; i--) {
        	for (j = 0; j < k; j++) {
            		get_block(X, blockE, n, m, i, j);
            		for (t = i + 1; t < k; t++) {
                		get_block(A, block, n, m, i, t);
                		get_block(X, blockCE, n, m, t, j);
                		Block_mult_diff(blockE, block, blockCE, m, m, m,m);
            		}
			if (l != 0) {
                		get_block(A, block, n, m, i, k);
                		get_block(X, blockCE, n, m, k, j);
                		Block_mult_diff(blockE, block, blockCE, m, m, l,m);
            		}
			put_block(X,blockE,n,m,i,j);
        	}
        	if (l!= 0) {
            		get_block(X, blockE, n, m, i, k);
            		for (t = i + 1; t < k; t++) {
                		get_block(A, block, n, m, i, t);
                		get_block(X, blockCE, n, m, t, k);
                		Block_mult_diff(blockE, block, blockCE, m, l, m,m);
            		}
            		get_block(A, block, n, m, i, k);
            		get_block(X, blockCE, n,m, k,k);
            		Block_mult_diff(blockE, block, blockCE, m, l, l,m);
            		put_block(X, blockE, n, m, i, k);
        	}
    	}
#endif
        return 0;
}
void mul_rotate (double *block, double *blockE, int m, int j, double *rcos,
            double *rsin, int l, int k)
{
        double cosk, sink, x, y;
        int p, q, c;
	int z = (j == k) ? l : m;
	if(j==k &&l==0)return;
        for (q = 0; q < m - 1; q++)
        {
                for (p = q + 1; p < m; p++)
                {
                        cosk = rcos[m * q + p];
                        sink = rsin[m * q + p];
                        for (c = 0; c < z; c++)	//произведение матриц поворота на Ar,j
                        {
                                x = cosk * block[q * m + c] - sink * block[p * m + c];
                                y = sink * block[q * m + c] + cosk * block[p * m + c];
                                block[q * m + c] = x;
                                block[p * m + c] = y;
                                x = cosk * blockE[q * m + c] - sink * blockE[p * m + c];
                                y = sink * blockE[q * m + c] + cosk * blockE[p * m + c];
                                blockE[q * m + c] = x;
                                blockE[p * m + c] = y;
                        }

                }
        }
}

int mul_fill_rotate2 (double *block, double *blockE, double *blockC,
                  double *blockCE, int m1, int m2, double *rcos, double *rsin,
                  double eps, double Anorm)
{
        double cosk, sq, sink, x, y;
        int p, q, c;
        for (q = 0; q < m1; q++)
        {
                for (p = 0; p < m2; p++)
                {
                        x = block[q * m1 + q];	//насчитываем cos и sin в матрицах поворота
                        y = blockC[p * m1 + q];
                        sq = (sqrt (x * x + y * y));
                        if (sq< eps * Anorm)
                                return -1;		// x=0 и y=0
                        rcos[m1 * q + p] = x / sq;
                        rsin[m1 * q + p] = -y / sq;
                        cosk = rcos[m1 * q + p];
                        sink = rsin[m1 * q + p];
                        for (c = 0; c < m1; c++)
                        {
                                x = cosk * block[q * m1 + c] - sink * blockC[p * m1 + c];
                                y = sink * block[q * m1 + c] + cosk * blockC[p * m1 + c];

                                block[q * m1 + c] = x;
                                blockC[p * m1 + c] = y;
                                x = cosk * blockE[q * m1 + c] - sink * blockCE[p * m1 + c];
                                y = sink * blockE[q * m1 + c] + cosk * blockCE[p * m1 + c];
                                blockE[q * m1 + c] = x;
                                blockCE[p * m1 + c] = y;
                        }

                }
        }
        return 0;

}

void mul_rotate2 (double *block, double *blockE, double *blockC, double *blockCE,
             int m1, int m2, int j, double *rcos, double *rsin, int l, int k)
{
        double cosk, sink, x, y;
        int p, q, c;
	int z = (j == k) ? l : m1;
        if(j==k && l==0) return;
        for (q = 0; q < m1; q++)
        {
                for (p = 0; p < m2; p++)
                {
                        cosk = rcos[m1 * q + p];
                        sink = rsin[m1 * q + p];
                        for (c = 0; c < z; c++)
                        {
                                x = cosk * block[q * m1 + c] - sink * blockC[p * m1 + c];
                                y = sink * block[q * m1 + c] + cosk * blockC[p * m1 + c];
                                block[q * m1 + c] = x;
                                blockC[p * m1 + c] = y;
                                x = cosk * blockE[q * m1 + c] - sink * blockCE[p * m1 + c];
                                y = sink * blockE[q * m1 + c] + cosk * blockCE[p * m1 + c];
                                blockE[q * m1 + c] = x;
                                blockCE[p * m1 + c] = y;
                        }

                }
        }
}

int mul_fill_rotate (double *block, double *blockE, int m,int m2, double *rcos,
                 double *rsin, double eps, double Anorm)
{
        int q, p, c;
        double cosk, sink, x, y, sq;
        for (q = 0; q < m2 - 1; q++)
        {
                for (p = q + 1; p < m2; p++)
                {
                        x = block[q * m + q];	//насчитываем cos и sin в матрицах поворота
                        y = block[p * m + q];
                        sq = (sqrt (x * x + y * y));
                        if (sq < eps * Anorm)
                                return -1;		// x=0 и y=0
                        rcos[m * q + p] = x / sq;
                        rsin[m * q + p] = -y / sq;
                        cosk = rcos[m * q + p];
                        sink = rsin[m * q + p];
                        for (c = 0; c < m2; c++)	//произведение матриц поворота на Ar,j
                        {
                                x = cosk * block[q * m + c] - sink * block[p * m + c];
                                y = sink * block[q * m + c] + cosk * block[p * m + c];
                                block[q * m + c] = x;
                                block[p * m + c] = y;
                                x = cosk * blockE[q * m + c] - sink * blockE[p * m + c];
                                y = sink * blockE[q * m + c] + cosk * blockE[p * m + c];
                                blockE[q * m + c] = x;
                                blockE[p * m + c] = y;

                        }

                }
        }
        return 0;

}

void get_block (double *matr, double *block, int real_n, int block_n, int i, int j)
{
        int k = real_n / block_n;
        int l = real_n - k * block_n;
        int width = (j < k ? block_n : l);
        int height = (i < k ? block_n : l);
        int row, col;

        double *source_block = matr + i * real_n * block_n + j * block_n;

        for (row = 0; row < height; row++)
        {
                for (col = 0; col < width; col++)
                {
                        block[row * block_n + col] = source_block[row * real_n + col];
                }
        }
}

void put_block (double *matr, double *block, int real_n, int block_n, int i, int j)
{
        int k = real_n / block_n;
        int l = real_n - k * block_n;
        int width = (j < k ? block_n : l);
        int height = (i < k ? block_n : l);
        int row, col;

        double *target_block = matr + i * real_n * block_n + j * block_n;

        for (row = 0; row < height; row++)
        {
                for (col = 0; col < width; col++)
                {
                        target_block[row * real_n + col] = block[row * block_n + col];
                }
        }
}
