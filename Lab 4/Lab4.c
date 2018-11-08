#include <stdio.h>
#include <complex.h>
#include <fftw3.h>
#include <math.h>

void fft(double* x[], double complex* psi[], double* k[], double complex* psik[])
{
    fftw_complex *in, *out;
    fftw_plan p;
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

int main()
{
    test_fft();
    
    return 0;
}