#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#define TINY 1.0e-20

/*
LU decomposition routine
input: - a: Input nxn matrix 
       - n: Size of matrix 
       - indx: Tracks permutation when performing pivoting
       - d: Tracks if permutation is even (+1) or odd (-1)
output: Nothing. Results record in a, indx, and d
*/
void ludcmp(double **a, int n, int *indx, double *d)
{
    int i, imax, j, k;
    double big, dum, sum, temp;
    double *vv;
    
    vv = malloc((n+1)*sizeof(double));
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
        for (i = 1; i < j; i++)
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
void lubksb(double **a, int n, int *indx, double b[])
{
    int i, ii = 0, ip, j;
    double sum;

    for (i = 1; i <= n; i++)
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
                sum -= a[i][j]*b[j];
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
        }
        b[i] = sum/a[i][i];
    }
}


int main()
{
    /* Question 1 */
    // double **a;
    // double mat[5][5] = {
    //     {0., 0., 0., 0., 0.},
    //     {0., 1., 3., 3., -5.},
    //     {0., 2., -4., 7., -1.},
    //     {0., 7., 1./2., 3., -6.},
    //     {0., 9., -2., 3., 8.}
    // };
    // double b[5] = {0., 0., 2., 3., -10.};
    // int n = 4, *indx;
    // double d;
    // double cons[5]; // Container for consistency check

    // // Setup matrix
    // a = (double **) malloc((n+1)*sizeof(double *));
    // for (int i = 0; i <= n; ++i)
    // {
    //     a[i] = (double *) malloc((n+1)*sizeof(double));
    // }

    // for (int i = 1; i <= 4; i++)
    // {
    //     for (int j = 1; j <= 4; j++)
    //     {
    //         a[i][j] = mat[i][j];
    //     }
    // }

    // // Setup index array
    // indx = (int *) malloc((n+1)*sizeof(int));

    // // LU decomposition
    // ludcmp(a, n, indx, &d);

    // // Backward Forward substitution
    // lubksb(a, n, indx, b);

    // // Readout solution
    // printf("Solution is: ");
    // for(int i = 1; i <= 4; i++)
    // {
    //     printf("%f\t", b[i]);
    // }
    // printf("\n");

    // // Consistency check
    // printf("Consistency check: b = ");
    // for(int i = 1; i <= 4; i++)
    // {
    //     cons[i] = 0.;
    //     for(int j = 1; j <= 4; j++)
    //     {
    //         cons[i] += mat[i][j]*b[j];
    //     }
    //     printf("%f\t",cons[i]);
    // }
    // printf("\n");



    // free(a);
    // free(indx);

    /* Question 2 
       We setup the solution such that effectively we are solving the equation
       M*V = b
       The solution of V is eventually saved in b so we make no distinction in the declaration
    */
    
    clock_t begin = clock(); // Start of clock
    int L=2; // Size of grid
    int N = (L+1)*(L+1); // Size of matrix
    int *indx; // Track permutation order
    indx = (int *) malloc((N+1)*sizeof(int));
    double d; // Track permutation is odd or even
    double r = 1/sqrt(L); // resistance
    double *V; //
    V = (double *) malloc((N+1)*sizeof(double));

    double **M; // Potential at each grid point
    M = (double **) malloc((N+1)*sizeof(double *));
    for (int i = 1; i <= N; i++)
    {
        M[i] = (double *) malloc((N+1)*sizeof(double));
    }
    

    // Setup the values for V
    for (int i = 1; i <= N-1; i++)
    {
        V[i] = 0.;
    }
    V[N] = 1.;

    // Setup the values for M
    for(int i = 1; i <= N; i++)
    {
        // Initialise everything as zero and change accordingly
        for(int j = 1; j <= N; j++)
        {
            M[i][j] = 0.;
        }

        if (i == 1) // VA
        {
            M[i][1] = 1.;
        }
        else if (i == N) // VB
        {
            M[i][N] = 1.;
        }
        else if (i == L+1) // Lower right corner
        {
            M[i][i-1] = 1.;
            M[i][i] = -2.;
            M[i][i+L+1] = 1.;
        }
        else if (i == N-L) // Upper left corner
        {
            M[i][i] = -2.;
            M[i][i+1] = 1.;
            M[i][i-L-1] = 1.;
        }
        else if (i < L+1) // Lower edge
        {
            M[i][i] = -3.;
            M[i][i-1] = 1.;
            M[i][i+1] = 1.;
            M[i][i+L+1] = 1.;
        }
        else if (i%(L+1) == 1) // Left edge 
        {
            M[i][i] = -3.;
            M[i][i+1] = 1.;
            M[i][i+L+1] = 1.;
            M[i][i-L-1] = 1.;
        }
        else if (i%(L+1) == 0) // Right edge
        {
            M[i][i] = -3.;
            M[i][i-1] = 1.;
            M[i][i+L+1] = 1.;
            M[i][i-L-1] = 1.;
        }
        else if (i > N-L) // Upper edge
        {
            M[i][i] = -3.;
            M[i][i-1] = 1.;
            M[i][i+1] = 1.;
            M[i][i-L-1] = 1.;
        }
        else // Rest of the middle grid points
        {
            M[i][i] = -4.;
            M[i][i-1] = 1.;
            M[i][i+1] = 1.;
            M[i][i-L-1] = 1.;
            M[i][i+L+1] = 1.;
        }
    }

    // LU decomposition
    ludcmp(M, N, indx, &d);

    // LU back substitution
    lubksb(M, N, indx, V);

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    // To calculate total current
    printf("Total current is: %f\n", (V[2]+V[L+2])/r);
    printf("Effective resistance is: %f\n", r/(V[2]+V[L+2]));
    printf("This take time %f seconds", time_spent);



    free(V);
    free(M);
    free(indx);
    printf("\nPress ENTER to exit...\n");
    getchar();
    return 0;
}