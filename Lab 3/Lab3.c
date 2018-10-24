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
       - steps: Number of steps in Monte Carlo
       - T: Temperature of the system
       - mean: storage for <E>
       - var: storage for <E^2>
       - wait: Number of flips before including the energy for mean and var
output: Nothing. Results are recorded in mean and var
*/
void MonteCarlo(int size, int steps, double T, double* mean, double* var, int wait)
{
    double A[size*size + 1], B[size*size + 1];
    double J = 1.;
    double energy = 0., Ediff = 0.;
    int index, flag;
    double count = 0.;
    double rate;

    *mean = 0.;
    *var = 0.;

    // Setup the initial grid
    for(int i = 1; i <= size*size; i++)
    {
        A[i] = drand64() > 0.5 ? 1. : -1.;
        B[i] = drand64() > 0.5 ? 1. : -1.;
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
            energy += -J*A[i-size+1]*B[i];
            energy += -J*A[i-2*size+1]*B[i];
        }
        else
        {
            energy += -J*A[i]*B[i];
            energy += -J*A[i+1]*B[i];
            energy += -J*A[i-size+1]*B[i];
        }
    }
    

    for(int i = 0; i < steps; i++)
    {
        // Pick a random site to flip
        index = (int) (drand64()*size*size + 1);

        // Then find out energy difference
        if(drand64() > 0.5)
        {
            // For B
            flag = 1;
            Ediff = 0.;
            if(index == size)
            {
                // Lower right corner for B
                Ediff += 2*J*A[index]*B[index];
                Ediff += 2*J*A[1]*B[index];
                Ediff += 2*J*A[size*size - size + 1]*B[index];
            }
            else if(index < size)
            {
                // Lower edge for B
                Ediff += 2*J*A[index]*B[index];
                Ediff += 2*J*A[index+1]*B[index];
                Ediff += 2*J*A[size*size-size+1+index]*B[index];
            }
            else if (index%size == 0)
            {
                // Right edge for B
                Ediff += 2*J*A[index]*B[index];
                Ediff += 2*J*A[index-size+1]*B[index];
                Ediff += 2*J*A[index-2*size+1]*B[index];
            }
            else
            {
                Ediff += 2*J*A[index]*B[index];
                Ediff += 2*J*A[index+1]*B[index];
                Ediff += 2*J*A[index-size+1]*B[index];
            }
            // printf("B: %f\t%d\n",Ediff, index);
        }
        else
        {
            // For A
            flag = 0;
            Ediff = 0.;
            if(index == size*size - size + 1)
            {
                // Upper left corner for A
                Ediff += 2*J*A[index]*B[index];
                Ediff += 2*J*A[index]*B[size*size];
                Ediff += 2*J*A[index]*B[size];
            }
            else if(index%size == 1)
            {
                // Left edge for A
                Ediff += 2*J*A[index]*B[index];
                Ediff += 2*J*A[index]*B[index + size - 1];
                Ediff += 2*J*A[index]*B[index + 2*size - 1];
            }
            else if(index > size*size-size+1)
            {
                // Upper edge for A
                Ediff += 2*J*A[index]*B[index];
                Ediff += 2*J*A[index]*B[index-1];
                if (index % size == 0)
                {
                    Ediff += 2*J*A[index]*B[size - 1];
                }
                else
                {
                    Ediff += 2*J*A[index]*B[index % size - 1];
                }
                
            }
            else
            {
                Ediff += 2*J*A[index]*B[index];
                Ediff += 2*J*A[index]*B[index-1];
                Ediff += 2*J*A[index]*B[index-1+size];
            }
            // printf("A: %f\n",Ediff);
        }

        // Calculate acceptance rate
        rate = fmin(1, exp(-Ediff/T));
        if(drand64()<rate)
        {
            if(flag == 0)
            {
                A[index] *= -1;
            }
            else
            {
                B[index] *= -1;
            }
            energy += Ediff;
        }
        // Calculate mean and variance
        // Only consider the energy after sufficient amount of flips
        if(i % wait == 0)
        {
            // Discard the first 500 points because the system is not yet in equilibrium
            if(count > 500)
            {
                *mean += energy;
                *var += energy*energy;
            }
            count += 1;
        }
    }
    *mean = *mean/(count-500);
    *var = *var/(count-500);
    printf("count = %f\t energy = %f\n", count, energy);
}

int main()
{
    double mean, var;
    int size = 16;
    int N = 2*size*size; //number of particles
    double T = 1.; //temperature
    double C, Cmean, Cvar;

    // Save data to file
    FILE *ofp;
    ofp = fopen("Lab3size16.dat", "w");

    // For various temperature
    for(double i = 1.; i <= 2.; i = i+0.025)
    {
        mean = 0.; var = 0.;
        C = 0.;
        Cvar = 0.;
        Cmean = 0.;
        // Run Monte Carlo 5 times for error bars
        for(int j = 0; j < 5; j++)
        {
            MonteCarlo(size, 20500*2000, i, &mean, &var, N);
            printf("Running T = %f, run = %d\n", i, j);
            C = 1/(i*i*N)*(var-mean*mean);
            Cmean += C;
            Cvar += C*C;
        }
        Cmean /= 5;
        Cvar /= 5;
        Cvar = Cvar - Cmean*Cmean;
        Cvar = sqrt(Cvar);

        fprintf(ofp, "%f\t%f\t%f\n", i, C, Cvar);
    }
    
    fclose(ofp);

    printf("\nPress ENTER to exit...\n");
    getchar();
    return 0;
}