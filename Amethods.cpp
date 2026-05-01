#include "header.h"
int Ainit(double* a, int n, int s, const std::string& name)
{
        int i, j;
        if (s == 1)
        {
                for (i = 0; i < n; i++)
                {
                        for (j = 0; j < i; j++)
                                a[i * n + j] = n - i;
                        for (j = i; j < n; j++)
                                a[i * n + j] = n - j;
                }
                return 0;
        }
        else if (s == 2)
        {
                for (i = 0; i < n; i++)
                {
                        for (j = 0; j < i; j++)
                                a[i * n + j] = i + 1;
                        for (j = i; j < n; j++)
                                a[i * n + j] = j + 1;
                }
                return 0;
        }
        else if (s == 3)
        {
                for (i = 0; i < n; i++)
                {
                        for (j = 0; j < i; j++)
                                a[i * n + j] = i - j;
                        for (j = i; j < n; j++)
                                a[i * n + j] = j - i;
                }
                return 0;
        }
        else if (s == 4)
        {
                for (i = 0; i < n; i++)
                {
                        for (j = 0; j < n; j++)
                                a[i * n + j] = 1.0 / (i + j + 1);
                }
                return 0;
        }
        else if( s == 0)
        {
                std::ifstream in;
                in.open(name);
                if (!in.is_open())
                {
                        printf("Ошибка: невозможно открыть файл %s\n", name.c_str());
                        return -1;
                }
                for ( i = 0; i < n; i++)
                {
                        for (j = 0; j < n; j++)
                        {
                                if (!(in >> a[i * n + j]))
                                {
                                        printf("Ошибка: недостаточно данных или неверный формат в файле %s\n", name.c_str());

                                        return -2;
                                }
                        }
                }
                in.close();
                return 0;
        }
        printf("s == %d ?!\n", s);
        return -3;
}
void Print_matrix(int l, int n, int r, double*A)
{
        int i, j;
        int row, col;
        row = l > r ? r : l;
        col = n > r ? r : n;
        for(i = 0; i < row; i++)
        {
                for(j = 0; j < col; j++) printf(" %10.3e", A[i * n + j]);
                printf("\n");
        }
        printf("\n");
}
