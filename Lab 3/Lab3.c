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
    double J = 1., sign;
    double energy = 0., Ediff;
    int index, flag;
    double T = 1;
    double rate;

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
            // Lower right corner for B
            energy += -J*A[i]*B[i];
            energy += -J*A[1]*B[i];
            energy += -J*A[size*size - size + 1]*B[i];
        }
        else if(i < size)
        {
            // Lower edge for B
            energy += -J*A[i]*B[i];
            energy += -J*A[i+1]*B[i];
            energy += -J*A[size*size-size+1+i]*B[i];
        }
        else if (i%size == 0)
        {
            // Right edge for B
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
    printf("Energy: %f\n",energy);

    for(int i = 0; i < 1000; i++)
    {
        // Pick a random site to flip
        index = (int) (drand64()*size*size + 1);

        // Then find out energy difference
        if(drand64() > 0.5)
        {
            // For B
            flag = 1;
            sign = -B[index];
            Ediff = 0.;
            if(index == size)
            {
                // Lower right corner for B
                Ediff += 2*J*A[index]*sign;
                Ediff += 2*J*A[1]*sign;
                Ediff += 2*J*A[size*size - size + 1]*sign;
            }
            else if(index < size)
            {
                // Lower edge for B
                Ediff += 2*J*A[index]*sign;
                Ediff += 2*J*A[index+1]*sign;
                Ediff += 2*J*A[size*size-size+1+index]*sign;
            }
            else if (index%size == 0)
            {
                // Right edge for B
                Ediff += 2*J*A[index]*sign;
                Ediff += 2*J*A[index-size+1];
                Ediff += 2*J*A[index-2*size+1];
            }
            else
            {
                Ediff += 2*J*A[index]*sign;
                Ediff += 2*J*A[index+1]*sign;
                Ediff += 2*J*A[index-size+1]*sign;
            }
            // printf("B: %f\t%d\n",Ediff, index);
        }
        else
        {
            // For A
            flag = 0;
            sign = -A[index];
            Ediff = 0.;
            if(index == size*size - size + 1)
            {
                // Upper left corner for A
                Ediff += 2*J*sign*B[index];
                Ediff += 2*J*sign*B[size*size];
                Ediff += 2*J*sign*B[size];
            }
            else if(index%size == 1)
            {
                // Right edge for A
                Ediff += 2*J*sign*B[index];
                Ediff += 2*J*sign*B[index + size - 1];
                Ediff += 2*J*sign*B[index + 2*size - 1];
            }
            else if(index > size*size-size+1)
            {
                // Upper edge for A
                Ediff += 2*J*sign*B[index];
                Ediff += 2*J*sign*B[index-1];
                Ediff += 2*J*sign*B[index%size - 1];
            }
            else
            {
                Ediff += 2*J*sign*B[index];
                Ediff += 2*J*sign*B[index-1];
                Ediff += 2*J*sign*B[index-1+size];
            }
            // printf("A: %f\n",Ediff);
        }

        // Calculate acceptance rate
        rate = fmin(1, exp(-Ediff/T));
        if(drand64()<rate)
        {
            
            if(flag == 0)
            {
                printf("%d\t%f\t%f\n",index, A[index], sign);
                A[index] = sign;
            }
            else
            {
                printf("%d\t%f\t%f\n",index, B[index], sign);
                B[index] = sign;
            }
            energy += Ediff;
        }
        printf("Energy: %f\t Ediff: %f\n",energy, Ediff);
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