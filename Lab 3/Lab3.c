#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// Random number generator
static unsigned long long int x = 1;
double drand64(void) 
{
   x = 6364136223846793005ll * x + (long long int) 1;
   return (double) x * 5.4210108624275218e-20; 
}

/* 
Monte Carlo Procedure for Honeycomb lattice
input: - size: the number of rows and columns of A and B parallelogram, which we assume to be size*size
*/
void MonteCarlo(int size)
{
    double A[size*size + 1], B[size*size + 1];
    // Setup the initial grid
    for(int i = 1; i <= size*size; i++)
    {
        A[i] = drand64() > 0.5 ? 1. : -1;
        B[i] = drand64() > 0.5 ? 1. : -1;
    }

    for(int i = 1; i <= size*size; i++)
    {
        printf("%f\t%f\n", A[i], B[i]);
    }
}

int main()
{
    MonteCarlo(3);
    printf("Hello World!\n");

    printf("\nPress ENTER to exit...\n");
    getchar();
    return 0;
}