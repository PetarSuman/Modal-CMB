#include <nag.h>
#include <math.h>
#include <stdio.h>
#include <nag_stdlib.h>
#include <nagf03.h>
#define NMAX 8
#define TDA NMAX
main()
{
double detf, determ, a[NMAX][TDA], p[NMAX];
Integer i, dete, j, n;
static NagError fail;
Vprintf("f03aec Example Program Results\n");
/* Skip heading in data file */
Vscanf("%*[^\n]");
Vscanf("%ld\n",&n);
if (n<1 || n>NMAX)
{
Vfprintf(stderr,"n is out of range: n = %5ld\n",n);
exit(EXIT_FAILURE);
}
for (i=0; i<n; i++)
for (j=0; j<n; j++)
Vscanf("%lf",&a[i][j]);
fail.print = TRUE;
f03aec(n, (double *)a, (Integer)TDA, p, &detf, &dete, &fail);
if (fail.code != NE_NOERROR)
exit(EXIT_FAILURE);
Vprintf("Array A after factorization\n");
for (i=0; i<n; i++)
for (j=0; j<n; j++)
Vprintf("%9.4f%s", a[i][j], (j%8==7 || j==n-1) ? "\n" : " ");
Vprintf("\nArray p\n");
for (i=0; i<n; i++)
Vprintf("%9.4f%s", p[i], (i%8==7 || i==n-1) ? "\n" : " ");
Vprintf("\ndetf = %9.4f dete = %2ld\n\n", detf, dete);
determ = detf*pow(2.0,(double)dete);
Vprintf("Value of determinant = %9.4f\n", determ);
exit(EXIT_SUCCESS);
}