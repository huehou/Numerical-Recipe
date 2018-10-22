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
    double J = 1.;
    double energy = 0.;

    // Setup the initial grid
    for(int i = 1; i <= size*size; i++)
    {
        A[i] = drand64() > 0.5 ? 1. : -1;
        B[i] = drand64() > 0.5 ? 1. : -1;
    }

    // Calculate initial energy
    for(int i = 1; i <= size*size; i++)
    {
        if(i == size)
        {
            energy += -J*A[i]*B[i];
            energy += -J*A[1]*B[i];
            energy += -J*A[size*size - size + 1]*B[i];
        }
        else if(i < size)
        {
            energy += -J*A[i]*B[i];
            energy += -J*A[i+1]*B[i];
            energy += -J*A[size*size-size+1+i]*B[i];
        }
        else if(i == size*size)
        {
            energy += -J*A[i]*B[i];
            energy += -J*A[size*size-size+1]*B[i];
            energy += -J*A[size*size-2*size+1]*B[i];
        }
        else if (i%size == 0)
        {
            energy += -J*A[i]*B[i];
            energy += -J*A[i-size+1];
            energy += -J*A[i-2*size+1];
        }
        else
        {
            energy += -J*A[i]*B[i];
            energy += -J*A[i+1]*B[i];
            energy += -J*A[i-size+1]*B[i];
        }
    }

    for(int i = 1; i <= size*size; i++)
    {
        printf("%f\t%f\n", A[i], B[i]);
    }
    printf("Energy is: %f\n", energy);
}

int main()
{
    MonteCarlo(3);
    printf("Hello World!\n");

    printf("\nPress ENTER to exit...\n");
    getchar();
    return 0;
}