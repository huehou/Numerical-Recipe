#include <stdio.h>
#include <complex.h>
#include <fftw3.h>
#include <math.h>
#define PI 3.1415926535897932385

/* 
Perform FFT on wave function in position space
input: - N: length of the position space
       - x: array of the position space
       - psi: The position wave function
       - k: the array of the momentum space
       - psik: the momentum wave function
output: None. The momentum space and wave function is recorded in k and psik respectively
*/
void fft(int N, double* x, double complex* psi, double* k, double complex* psik)
{
    fftw_complex *in, *out;
    fftw_plan p;

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);

    // FFTSHIFT: Shift the domain for FFT
    for(int i = 0; i < N; i++)
    {
        in[i] = psi[(i+N/2)%N];
    }

    // Perform FFT
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);
    p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);
    
    // IFFTSHIFT: Shift the domain back to normal
    double dx = x[1]-x[0];
    for(int i = 0; i < N; i++)
    {
        psik[i] = out[(i+N/2)%N]*dx/sqrt(2*PI);
    }

    // Momentum domain
    for(int i = 0; i < N; i++)
    {
        k[i] = 2 * PI * (i-N/2) / (N*dx);
    }

    fftw_destroy_plan(p);
    fftw_free(in); fftw_free(out);
}

/* 
Perform IFFT on wave function in momentum space
input: - N: length of the position space
       - x: array of the position space
       - psi: The position wave function
       - k: the array of the momentum space
       - psik: the momentum wave function
output: None. The position space and wave function is recorded in x and psi respectively
*/
void ifft(int N, double* x, double complex* psi, double* k, double complex* psik)
{
    fftw_complex *in, *out;
    fftw_plan p;

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);

    // FFTSHIFT: Shift the domain for IFFT
    for(int i = 0; i < N; i++)
    {
        in[i] = psik[(i+N/2)%N];
    }

    // Perform IFFT
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);
    p = fftw_plan_dft_1d(N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(p);
    
    // IFFTSHIFT: Shift the domain back to normal
    double dk = k[1]-k[0];
    for(int i = 0; i < N; i++)
    {
        psi[i] = out[(i+N/2)%N]*dk/sqrt(2*PI);
    }

    // Momentum domain
    for(int i = 0; i < N; i++)
    {
        x[i] = 2 * PI * (i-N/2) / (N*dk);
    }

    fftw_destroy_plan(p);
    fftw_free(in); fftw_free(out);
}

/* 
Integration by Trapezoidal Rule
input: - psi: The integrand array
       - start: The starting index
       - N: The length of the array
       - dx: The spacing of the position
output: The value of the integration
*/
double integrate(double *psi, int start, int N, double dx)
{
    double sum = 0.;
    for(int i = start+1; i < N-1; i++)
    {
        sum += psi[i]*dx;
    }
    sum += (psi[start] + psi[N-1])*dx/2.;
    return sum;
}

void Problem1b()
{
    int N = 10000;
    double x[N]; // position
    double k[N]; // momentum
    double complex psi[N]; // position wave function
    double complex psik[N]; // momentum wave function
    double V[N]; // potential energy
    double x0 = -500., dx = 1600./N;
    double a = sqrt(50.);

    // Initial position space
    for(int i = 0; i < N; i ++)
    {
        x[i] = -800. + dx*i;
    }

    // Potential energy
    int index;
    for(int i = 0; i < N; i ++)
    {
        if (0 <= x[i] && x[i] < a)
        {
            V[i] = 1;
            index = i;
        }
        else
        {
            V[i] = 0;
        }
    }

    // Initial wave function
    double complex m = 1., hbar = 1.;
    double sigma = 90.;
    double dE = 0.01;
    double complex p0;
    double prob0[N];

    FILE *ofp;
    ofp = fopen("Problem1b.dat","w");


    for(double E = 2.; E > 0.6; E -= dE)
    {
        printf("Doing E = %f \n", E);
        p0 = sqrt(2*m*E);
        for(int i = 0; i < N; i++)
        {
            psi[i] = cexp(-1./(2.*sigma*sigma)*(x[i]-x0)*(x[i]-x0) + I/hbar*p0*(x[i]-x0));
            prob0[i] = cabs(psi[i])*cabs(psi[i]);
        }

        double norm = integrate(prob0, 0, N, dx);

        // Suzuki-Trotter procedure
        double T = 800./p0;
        double dt = dx/p0;
        for(double t = 0; t <= T; t += dt)
        {
            for(int i = 0; i < N; i++)
            {
                psi[i] = cexp(-I/hbar*dt*V[i])*psi[i];
            }

            fft(N, x, psi, k, psik);

            for(int i = 0; i < N; i++)
            {
                psik[i] = cexp(-I/hbar*dt*k[i]*k[i]/2/m)*psik[i];
            }

            ifft(N, x, psi, k, psik);
        }
        
        double prob[N];
        for(int i = 0; i < N; i++)
        {
            prob[i] = cabs(psi[i])*cabs(psi[i]);
        }
        double trans = integrate(prob, index, N, dx);
        trans = trans/norm;
        printf("Transmission probability for E = %f is %f.\n", E, trans);
        fprintf(ofp, "%f\t%f\n", E, trans);
    }
    

    fclose(ofp);
}

void Problem1c()
{
    int N = 10000;
    double x[N]; // position
    double k[N]; // momentum
    double complex psi[N]; // position wave function
    double complex psik[N]; // momentum wave function
    double V[N]; // potential energy
    double x0 = -500., dx = 1600./N;
    double a = sqrt(50.);

    // Initial position space
    for(int i = 0; i < N; i ++)
    {
        x[i] = -800. + dx*i;
    }

    // Potential energy
    int index;
    for(int i = 0; i < N; i ++)
    {
        if (0 <= x[i] && x[i] < a)
        {
            V[i] = 4/(a*a)*(a-x[i])*x[i];
            index = i;
        }
        else
        {
            V[i] = 0;
        }
    }

    // Initial wave function
    double complex m = 1., hbar = 1.;
    double sigma = 90.;
    double dE = 0.01;
    double complex p0;
    double prob0[N];

    FILE *ofp;
    ofp = fopen("Problem1c.dat","w");


    for(double E = 2.; E > 0.6; E -= dE)
    {
        printf("Doing E = %f \n", E);
        p0 = sqrt(2*m*E);
        for(int i = 0; i < N; i++)
        {
            psi[i] = cexp(-1./(2.*sigma*sigma)*(x[i]-x0)*(x[i]-x0) + I/hbar*p0*(x[i]-x0));
            prob0[i] = cabs(psi[i])*cabs(psi[i]);
        }

        double norm = integrate(prob0, 0, N, dx);

        // Suzuki-Trotter procedure
        double T = 800./p0;
        double dt = dx/p0;
        for(double t = 0; t <= T; t += dt)
        {
            for(int i = 0; i < N; i++)
            {
                psi[i] = cexp(-I/hbar*dt*V[i])*psi[i];
            }

            fft(N, x, psi, k, psik);

            for(int i = 0; i < N; i++)
            {
                psik[i] = cexp(-I/hbar*dt*k[i]*k[i]/2/m)*psik[i];
            }

            ifft(N, x, psi, k, psik);
        }
        
        double prob[N];
        for(int i = 0; i < N; i++)
        {
            prob[i] = cabs(psi[i])*cabs(psi[i]);
        }
        double trans = integrate(prob, index, N, dx);
        trans = trans/norm;
        printf("Transmission probability for E = %f is %f.\n", E, trans);
        fprintf(ofp, "%f\t%f\n", E, trans);
    }
    

    fclose(ofp);
}

int main()
{
    Problem1b();
    
    return 0;
}