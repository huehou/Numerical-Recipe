#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define TINY 1.0e-20

/*
LU decomposition routine
input: - a: Input nxn matrix 
       - n: Size of matrix 
       - indx: Tracks permutation when performing pivoting
       - d: Tracks if permutation is even (+1) or odd (-1)
output: Nothing. Results record in a, indx, and d
*/
void ludcmp(float **a, int n, int *indx, float *d)
{
    int i, imax, j, k;
    float big, dum, sum, temp;
    float *vv;

    vv = malloc((n+1)*sizeof(float));
    *d = 1.0;

    // Track scaling for each row
    for (i = 1; i <= n; i++)
    {
        big = 0.0; 
        for (j = 1; j <= n; j++)
        {
            if ((temp = fabs(a[i][j])) > big )
            {
                big = temp;
            }
        }
        if (big == 0.0)
        {
            // If the largest element in row is 0 ==> Singular matrix
            printf("Singular matrix in routine ludcmp");
            exit(-1);
        }
        // Keep track of largest scaling of each row in vv
        vv[i] = 1.0/big;
    }

    // Begins Crout's method
    for (j = 1; j <= n; j++)
    {
        // Calculate beta
        for (i = j; i < j; i++)
        {
            sum = a[i][j];
            for (k = 1; k < i; k++)
            {
                sum -= a[i][k] * a[k][j];
            }
            a[i][j] = sum;
        }

        // Setup for pivoting
        big = 0.0;
        for (i = j; i <= n; i++)
        {
            // Calculate alpha
            sum = a[i][j];
            for ( k = 1; k < j; k++)
            {
                sum -= a[i][k] * a[k][j];
            }
            a[i][j] = sum;

            // Find the largest alpha in the column (partial pivot)
            if ( (dum = vv[i]*fabs(sum)) >= big)
            {
                big = dum;
                imax = i;
            }
        }

        // Pivot
        if (j != imax){
            // Shuffle the rows
            for (k = 1; k <= n; k++)
            {
                dum = a[imax][k];
                a[imax][k] = a[j][k];
                a[j][k] = dum;
            }
            *d = -(*d);
            vv[imax] = vv[j];
        }

        indx[j] = imax;
        if (a[j][j] == 0.0) 
        {
            a[j][j] = TINY;
        }

        if (j != n)
        {
            dum = 1.0/(a[j][j]);
            for (i = j+1; i <= n; i++)
            {
                a[i][j] *= dum;
            }
        }
    }
    free(vv);
}

/* 
LU back substitution routine
input: - a: The input nxn matrix
       - n: the size of the matrix
       - indx: The permutation index
       - b: The storage for solution
output: Nothing
*/
void lubksb(float **a, int n, int *indx, float b[])
{
    int i, ii = 0, ip, j;
    float sum;

    for (i = 1; i <= n, i++)
    {
        // Permute b because of pivoting
        ip = indx[i];
        sum = b[ip];
        b[ip] = b[i];
        if (ii)
        {
            // Forward substitution
            for (j = ii; j <= i - 1; j++)
            {
                sum -= a[i][j]*b[j]
            }
        }
        else if (sum)
        {
            // Setup first non-zero index
            ii = i;
        }
        b[i] = sum;
    }

    // Backward substitution
    for (i = n; i >= 1; i--)
    {
        sum = b[i];
        for (j = i+1; j <= n; j++)
        {
            sum -= a[i][j]*b[j];
            b[i] = sum/a[i][i];
        }
    }
}

int main()
{
    printf("Hello World!\n");

    printf("\nPress ENTER to exit...\n");
    getchar();
    return 0;
}