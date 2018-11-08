#include <stdio.h>
#include <complex.h>
#include <fftw3.h>
#include <math.h>

void fft(int N, double* x, double complex* psi, double* k, double complex* psik)
{
    fftw_complex *in, *out;
    fftw_plan p;

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);

    for(int i = 0; i < N; i++)
    {
        in[i] = psi[(i+N/2)%N];
    }
    for(int i = 0; i < N; i++)
    {
        printf("%f\t",cabs(in[i]));
    }
    printf("\n");
    for(int i = 0; i < N; i++)
    {
        printf("%f\t",cabs(psi[i]));
    }
}

void test_fft()
{
    double x[10000];
    double complex psi[10000];
    double x0 = -100., dx = 200./10000.;
    FILE *ofp;
    ofp = fopen("trial.dat", "w");
    for(int i = 0; i < 10000; i++)
    {
        x[i] = x0 + dx*i;
        psi[i] = exp(-x[i]*x[i]);
        // printf("%f\t%f\n", x[i], cabs(psi[i]));
        fprintf(ofp, "%f\t%f\n", x[i], cabs(psi[i]));
    }

    fclose(ofp);

    fftw_complex *in, *out;
    fftw_plan p;
    int N = 10000;

    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);
    p = fftw_plan_dft_1d(N, psi, out, FFTW_FORWARD, FFTW_ESTIMATE);
    
    fftw_execute(p);
    
    

    // FILE *ofp;
    ofp = fopen("fft.dat", "w");
    for(int i = 0; i < 10000; i++)
    {
        printf("%f\n", cabs(out[i]));
        fprintf(ofp, "%f\n", cabs(out[i]));
    }

    fclose(ofp);

    fftw_destroy_plan;
    fftw_free(in); fftw_free(out);

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
    test2();
    
    return 0;
}