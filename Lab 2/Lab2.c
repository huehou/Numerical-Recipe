#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define TYPE float // TYPE is defined for easy switching (between float and double)
#define FUNC(x) ((*func)(x))
#define EPS 1.0e-6
#define JMAX 20
#define JMAXP (JMAX+1)
#define K 9
#define PI 3.1415926535897932385

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
            printf("%.24f\t%.24f\n", ss, dss);
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

int Question2()
{
    TYPE s;
    s = qromb(integrand, 0, 2);
    printf("%.24f\n",s);
    
    return 0;
}

static unsigned long long int x = 1;

// Random number generator
double drand64(void) 
{
   x = 6364136223846793005ll * x + (long long int) 1;
   return (double) x * 5.4210108624275218e-20; 
}

void MonteCarlo(double d, double steps, double *mean)
{
    double r1, r2, theta1, theta2, phi1, phi2;
    double x1, x2, y1, y2, z1, z2;
    double dr1, dr2, dtheta1, dtheta2, dphi1, dphi2;
    double dx1, dx2, dy1, dy2, dz1, dz2;
    double x1new, x2new, y1new, y2new, z1new, z2new;
    double rA1, rA2, rB1, rB2, r12, psi, E;
    double rA1new, rA2new, rB1new, rB2new, r12new, psinew, Enew;
    *mean = 0;
    int counter=0;

    // Starting point of r1 and r2
    r1 = 2.*drand64() - 1.;
    theta1 = drand64()*PI;
    phi1 = drand64()*2*PI;
    r2 = 2.*drand64() - 1.;
    theta2 = drand64()*PI;
    phi2 = drand64()*2*PI;
    x1 = r1*sin(theta1)*cos(phi1);
    y1 = r1*sin(theta1)*sin(phi1);
    z1 = r1*cos(theta1);
    x2 = r2*sin(theta2)*cos(phi2);
    y2 = r2*sin(theta2)*sin(phi2);
    z2 = r2*cos(theta2);

    // Other quantities with rA = (0,0,0) and rB = (0,0,d)
    rA1 = sqrt(x1*x1+y1*y1+z1*z1);
    rA2 = sqrt(x2*x2+y2*y2+z2*z2);
    rB1 = sqrt(x1*x1+y1*y1+(z1-d)*(z1-d));
    rB2 = sqrt(x2*x2+y2*y2+(z2-d)*(z2-d));
    r12 = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));
    psi = exp(-rA1)*exp(-rB2) + exp(-rA2)*exp(-rB1);
    E = -((1.-1./rA1-1./rB2)*exp(-rA1-rB2) + (1. - 1./rA2 - 1./rB1)*exp(-rA2-rB1))/psi;
    E = E - (1./rA1 + 1./rB1 + 1./rA2 + 1./rB2) + (1./d + 1./r12);

    for (int i = 0; i < steps; i++)
    {
        // Generate step
        dr1 = (2.*drand64() - 1.)*3;
        dtheta1 = drand64()*PI;
        dphi1 = drand64()*2*PI;
        dr2 = (2.*drand64() - 1.)*3;
        dtheta2 = drand64()*PI;
        dphi2 = drand64()*2*PI;
        dx1 = dr1*sin(dtheta1)*cos(dphi1);
        dy1 = dr1*sin(dtheta1)*sin(dphi1);
        dz1 = dr1*cos(dtheta1);
        dx2 = dr2*sin(dtheta2)*cos(dphi2);
        dy2 = dr2*sin(dtheta2)*sin(dphi2);
        dz2 = dr2*cos(dtheta2);

        // New place
        x1new = x1 + dx1;
        x2new = x2 + dx2;
        y1new = y1 + dy1;
        y2new = y2 + dy2;
        z1new = z1 + dz1;
        z2new = z2 + dz2;

        // Other quantities with rA = (0,0,0) and rB = (0,0,d)
        rA1new = sqrt(x1new*x1new+y1new*y1new+z1new*z1new);
        rA2new = sqrt(x2new*x2new+y2new*y2new+z2new*z2new);
        rB1new = sqrt(x1new*x1new+y1new*y1new+(z1new-d)*(z1new-d));
        rB2new = sqrt(x2new*x2new+y2new*y2new+(z2new-d)*(z2new-d));
        psinew = exp(-rA1new)*exp(-rB2new) + exp(-rA2new)*exp(-rB1new);
        
        double rate = fmin(1, (psinew*psinew)/(psi*psi));
        if (drand64() <= rate)
        {
            counter += 1;
            x1 = x1new;
            y1 = y1new;
            z1 = z1new;
            x2 = x2new;
            y2 = y2new;
            z2 = z2new;
            rA1 = rA1new;
            rA2 = rA2new;
            rB1 = rB1new;
            rB2 = rB2new;
            psi = psinew;
            r12 = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));
            E = -((1.-1./rA1-1./rB2)*exp(-rA1-rB2) + (1. - 1./rA2 - 1./rB1)*exp(-rA2-rB1))/psi;
            E = E - (1./rA1 + 1./rB1 + 1./rA2 + 1./rB2) + (1./d + 1./r12);
        }
        *mean += E;
    }
    
    *mean /= steps;
    // printf("Acceptance rate: %f\n", counter/steps);
}

int Question3()
{
    // Unit conversion:
    double hartree = 27.211385;
    double bohr = 0.52917721092;
    double mean;
    double res, error;
    double times = 100.;
    double d;

    // Save data to file
    FILE *ofp;
    ofp = fopen("Lab2Q3.dat", "w");

    for (int j = 0; j < 108; j++)
    {
        d = 0.8 + j*0.05;
        res = 0;
        error = 0;

        for(int i = 0; i < times; i++)
        {
            MonteCarlo(d, 100000, &mean);
            // printf("%f\n", mean);
            res += mean;
            error += mean*mean;
        }

        res /= times;
        error /= times;
        error -= res*res;
        error = sqrt(error/(times-1));
        printf("rAB = %f A\n", d*bohr);
        printf("Result is %f with error %f \n", res*hartree, error*hartree);

        fprintf(ofp, "%f\t%f\t%f\n", d*bohr, res*hartree, error*hartree);
    }
    
}

int main()
{
    Question3();

    printf("\nPress Enter to exit...");
    getchar();
    return 0;
}