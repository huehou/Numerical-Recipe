#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define TYPE float // TYPE is defined for easy switching (between float and double)

/* 
Polynomial interpolation routine
input: - xa: input array of x values
       - ya: input array of y values
       - n: Length of input array
       - x: The point of interest
       - y: The interpolated value
       - dy: The error of fitting
output: Nothing. Values and errors are stored in y and dy respectively
*/
void polint(TYPE xa[], TYPE ya[], int n, TYPE x, TYPE *y, TYPE *dy)
{
    int i, m, ns = 1; 
    TYPE den, dif, dift, ho, hp, w;

    TYPE *c, *d; // The difference vector

    dif = fabs(x-xa[1]); // Compute abs difference between x and the first point in xa

    // Initialise c and d
    c = malloc((n+1)*sizeof(TYPE));
    d = malloc((n+1)*sizeof(TYPE));

    // Find the index in xa to which x is closest
    for (i = 1; i <= n; i++)
    {
        if ( (dift = fabs(x-xa[i])) < dif )
        {
            ns = i;
            dif = dift;
        }
        // Initialise c and d
        c[i] = ya[i];
        d[i] = ya[i];
    }

    // Initialise y to value closest to x
    *y = ya[ns--];

    for (m = 1; m < n; m++)
    {
        for (i = 1; i <= n - m; i++)
        {
            ho = xa[i] - x;
            hp = xa[i+m] - x;
            w = c[i+1] - d[i];
            if ( (den = ho-hp) == 0.0)
            {
                printf("Error in routine polint");
            }
            den = w/den;
            d[i] = hp*den;
            c[i] = ho*den;
        }

        *y += (*dy = (2*ns < (n-m) ? c[ns+1] : d[ns--]) );
    }
    
    free(c);
    free(d);
}

int main()
{
    // Question 1

    // Define the array
    TYPE xa[5] = {0., -1., 1., 2., 4.};
    TYPE ya[5] = {0., 1.25, 2., 3., 0.};

    // Build the x-axis and y-axis for interpolation
    int pointsNum = 100;
    TYPE x[pointsNum + 1];
    TYPE y[pointsNum + 1];
    TYPE dy[pointsNum + 1];
    TYPE dx = 5./(pointsNum-1);

    for(int i = 1; i <= pointsNum; i++)
    {
        x[i] = -1 + (i-1)*dx;
    }

    // Save data to file
    FILE *ofp;
    ofp = fopen("Lab2Q1.dat", "w");

    for(int i = 1; i <= pointsNum; i++)
    {
        polint(xa, ya, 4, x[i], &y[i], &dy[i]);
        fprintf(ofp, "%f\t%f\n", x[i], y[i]);
    }


    // for(int i = 1; i <= pointsNum; i++)
    // {
    //     printf("%f\n", x[i]);
    // }
    // printf("\n");

    printf("\nPress Enter to exit...");
    getchar();
    return 0;
}