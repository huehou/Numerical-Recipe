#include <stdio.h>
#include <complex.h>
#include <fftw3.h>
#include <math.h>
#define PI 3.1415926535897932385

/* 
Perform FFT on wave function
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
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);

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

void ifft(int N, double* x, double complex* psi, double* k, double complex* psik)
{
    fftw_complex *in, *out;
    fftw_plan p;

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);

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

void test_fft()
{
    int N = 10000;
    double x[N];
    double k[N];
    double complex psi[N];
    double complex psik[N];
    double x0 = -100., dx = 200./N;
    FILE *ofp;
    ofp = fopen("trial.dat", "w");
    for(int i = 0; i < N; i++)
    {
        x[i] = x0 + dx*i;
        psi[i] = exp(-(x[i]-5)*(x[i]-5));
        // printf("%f\t%f\n", x[i], cabs(psi[i]));
        fprintf(ofp, "%f\t%f\n", x[i], cabs(psi[i]));
    }

    fclose(ofp);

    fft(N, x, psi, k, psik);
    ifft(N, x, psi, k, psik);

    // FILE *ofp;
    ofp = fopen("fft.dat", "w");
    for(int i = 0; i < N; i++)
    {
        printf("%f\t%f\n", k[i], cabs(psi[i]));
        fprintf(ofp, "%f\t%f\t%f\n", x[i], creal(psi[i]),cimag(psi[i]));
    }

    fclose(ofp);

}

void test2()
{
    double x[10];
    double complex psi[10];

    for(int i = -5 ; i < 5; i ++)
    {
        x[i+5] = i;
        psi[i+5] = i;
    }
    
    fft(10, x, psi, x, psi);
    fft(10, x, psi, x, psi);
}

int main()
{
    test_fft();
    
    return 0;
}