#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define TYPE float // TYPE is defined for easy switching (between float and double)
#define FUNC(x) ((*func)(x))
#define EPS 1.0e-6
#define JMAX 20
#define JMAXP (JMAX+1)
#define K 5

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

/*
Trapezoidal Rule Integration
input: - func: Integrand as a function
       - a: Lower limit of integration
       - b: Upper limit of integration
       - n: Number of calls
output: Result of integration
*/
TYPE trapzd(TYPE (*func)(TYPE), TYPE a, TYPE b, int n)
{
    TYPE x, tnm, sum, del;
    static TYPE s;
    int it, j;

    if (n == 1)
    {
        return (s = 0.5*(b-a)*(FUNC(a)+FUNC(b)));
    }
    else
    {
        for (it = 1, j = 1; j < n-1; j++) 
        {
            it <<= 1;
        }
        tnm = it; 
        del = (b-a)/tnm;
        x = a + 0.5*del;
        for (sum = 0.0, j = 1; j <= it; j++, x += del)
        {
            sum += FUNC(x);
        }
        s = 0.5*(s+(b-a)*sum/tnm);
        return s;
    }
}

/* 
Romberg integration routine
input: - func: Integrand as a function
       - a: Lower integration limit
       - b: Upper integration limit
output: Result of integration
*/
TYPE qromb(TYPE (*func)(TYPE), TYPE a, TYPE b)
{
    TYPE ss, dss;
    TYPE s[JMAXP], h[JMAXP+1];
    int j;

    h[1] = 1.0;
    for (j = 1; j <= JMAX; j++)
    {
        s[j] = trapzd(func, a, b, j);
        if (j >= K)
        {
            polint(&h[j-K], &s[j-K], K, 0.0, &ss, &dss);
            if (fabs(dss) <= EPS*fabs(ss)) return ss;
        }
        h[j+1] = 0.25*h[j];
    }
    printf("Too many steps in routine qromb");
    return 0.0;
}

// Integrand of question 2
TYPE integrand(TYPE x)
{
    return x*x*x*x * log(x + sqrt(x*x + 1));
}

int Question1()
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

    return 0;
}

int Question2()
{
    TYPE s;
    for(int i = 1; i <= 20; i++)
    {
        s = trapzd(integrand, 0, 2, i);
        printf("%.16f\n", s);   
    }

    s = qromb(integrand, 0, 2);
    printf("%.16f\n",s);
    
    return 0;
}

int main()
{
    Question2();

    printf("\nPress Enter to exit...");
    getchar();
    return 0;
}