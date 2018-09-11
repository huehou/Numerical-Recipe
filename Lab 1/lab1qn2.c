#include <math.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <time.h>
#define NR_END 1
#define FREE_ARG char*
#define TINY 1.0e-20

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

float *vector(long nl, long nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
	float *v;

	v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl+NR_END;
}

void ludcmp(float **a, int n, int *indx, float *d)
{
   int i,imax,j,k;
   float big,dum,sum,temp;
   float *vv;
   vv=malloc((n+1)*sizeof(float));
   *d=1.0;
   
   for (i=1;i<=n;i++) {
      big=0.0;
      for (j=1;j<=n;j++)
         if ((temp=fabs(a[i][j])) > big) big=temp;
      if (big == 0.0) nrerror("Singular matrix in routine ludcmp");
      vv[i]=1.0/big;
   }

   for (j=1;j<=n;j++) {
      for (i=1;i<j;i++) {
         sum=a[i][j];
         for (k=1;k<i;k++) 
            sum -= a[i][k]*a[k][j];
            a[i][j]=sum;
      }

      big=0.0;
      for (i=j;i<=n;i++) {
         sum=a[i][j];
         for (k=1;k<j;k++)
            sum -= a[i][k]*a[k][j];
         a[i][j]=sum;
         if ( (dum=vv[i]*fabs(sum)) >= big) {
            big=dum;
            imax=i;
         }
      }
      if (j != imax) {
         for (k=1;k<=n;k++) {
            dum=a[imax][k];
            a[imax][k]=a[j][k];
            a[j][k]=dum;
         }
         *d = -(*d);
         vv[imax]=vv[j];
      }
      indx[j]=imax;
      if (a[j][j] == 0.0) a[j][j]=TINY;
      if (j != n) {
         dum=1.0/(a[j][j]);
         for (i=j+1;i<=n;i++) a[i][j] *= dum;
      }
   }

   free(vv);
}


//declare lubksb
void lubksb(float **a, int n, int *indx, float b[]){
   int i,ii=0,ip,j;
   float sum;

   for (i=1;i<=n;i++){
      ip=indx[i];
      sum=b[ip];
      b[ip]=b[i];
      if (ii)
         for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
      else if (sum) ii=i;
    b[i]=sum;
   }

   for (i=n;i>=1;i--) {
      sum=b[i];
      for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
      b[i]=sum/a[i][i];
   }
}


//call ludcmp, lubksb
int main(){

   float **a;
   int i, j, l, n;
   int *indx;
   float *b;
   float d;
   b = (float *)malloc((n+1)*sizeof(float));

   clock_t start = clock(); /*start timing*/

   //input matrix size
   printf("Enter L\n");
   scanf("%d", &l);
   n = (l+1)*(l+1);
   indx = (int *)malloc((n+1)*sizeof(int));
   a = (float **)malloc((n+1)*sizeof(float *));
   for(i=1; i<=n; i++)
   {
      a[i] = (float *)malloc((n+1)*sizeof(float));
   }
   
   //input coefficient matrix elements and corresponding col vector
  
  for(i=1; i<=n; i++)
  {
     b[i] = 0;
  }
  b[n] = 1;
  
   /*set up V_A and V_B*/
   a[l+1][l+1] = 1;
   a[n-l][n-l] = 1;
 
   /*upper left corner*/
   a[1][1] = -2;
   a[1][2] = 1;
   a[1][l+2] = 1;
 
   /*lower right corner*/
   a[n][n-(l+1)] = 1;
   a[n][n-1] = 1;
   a[n][n] = -2;

   /*top side: T shape*/
   //i initialised at 2 because exclude top left corner.
   //i stops at l because at l+1 is Vb
  for(i=2; i<=l; i++)
  {
     a[i][i-1] = 1;
     a[i][i] = -3;
     a[i][i+1] = 1;
     a[i][i+l+1] = 1;
  }

  /*bottom side: inverted T shape*/
  //i initialised at n-l+1 because exclude V_a.
  //i stops at n-1 because n corresponds to bottom right corner.
  for(i=n-l+1; i<=n-1; i++)
  {
     a[i][i-(l+1)] = 1;
     a[i][i-1] = 1;
     a[i][i] = -3;
     a[i][i+1] = 1;
  }

  /*left side: left T shape*/
  //i initialised at 2 because exclude top left hand corner.
  //i stops at n-l-1 because n-l corresponds to V_a.
  for(i=2; i<=n-l-1; i++)
  {
    printf("i:%d\n",i);
    printf("1\n");
     a[i-1][i] = 1; 
     printf("2\n");
     a[i+l][i] = -3;
     printf("3\n");
     a[i+l+1][i] = 1;
     printf("4\n");
     a[i+2*l+1][i] = 1;
     printf("5\n");
  }

  /*right side: right T shape*/
  //i initialised at l+2 because l+1 is V_B
  //i stops at n-1 because n corresponds to bottom right corner.
  for(i=l+2; i<=n-1; i++)
  {
     a[i-1][i] = 1;
     a[i+(l-1)][i] = 1;
     a[i+l][i] = -3;
     a[i+(2*l+1)][i] = 1;
  }

  /*all the middle ones shaped like a vertical cross*/
  for(i=2; i<=l; i++)
  {
     for(j=2; j<=l; j++)
     {
        a[i][j] = 1;
        a[i][j+l] = 1;
        a[i][j+(l+1)] = -4;
        a[i][j+(l+2)] = 1;
        a[i][j+2*(l+1)] = 1;
     }
  }

  //using the functions
  ludcmp(a, n, indx, &d);
  lubksb(a, n, indx, b);

  clock_t stop = clock(); /*end timing*/
  double time = (double)(stop - start)/CLOCKS_PER_SEC;
  //output

  //using kirchoff's law
  printf("Total current:%f", (b[2]+b[l+2])*sqrt(l));//r = 1/sqrt(l)
  
  printf("Time taken: %f s", time);

  free(b);
  free(a);
  free(indx);

}
