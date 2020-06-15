#include "omp.h"
#include <cstdio>
#include <cassert>
#include <ctime>

double startTime, endTime;
double t1, t2, t3;

void HinesAlgo(double *u, double *l, double *d, double *rhs, int *p, int N)
{
    int i;
    double factor;
    t1 = omp_get_wtime();
    // Backward sweep
    for (i = N - 1; i > 0; i--)
    {
        factor = u[i] / d[i];
        d[p[i]] -= factor * l[i];
        rhs[p[i]] -= factor * rhs[i];
    }
    t2 = omp_get_wtime();
    // Forward sweep
    rhs[0] /= d[0];
    for (i = 1; i < N; i++)
    {
        rhs[i] -= l[i] * rhs[p[i]];
        rhs[i] /= d[i];
    }
    t3 = omp_get_wtime();
}

int main(int argc, char const *argv[])
{
    int N, i, idx;
    double *upper, *lower, *diag, *rhs;
    int *parent;
    FILE *fin, *fout;

    // read data from file
    assert(argc == 3);
    fin = fopen(argv[1], "r");
    fscanf(fin, "%d", &N);
    upper = new double[N];
    lower = new double[N];
    diag = new double[N];
    rhs = new double[N];
    parent = new int[N];

    for (i = 0; i < N; i++)
    {
        fscanf(fin, "%d", &idx);
        fscanf(fin, "%lf %lf %lf %lf %d",
               &upper[idx], &lower[idx], &rhs[idx], &diag[idx], &parent[idx]);
    }
    fclose(fin);

    // Hines algorithm
    startTime = omp_get_wtime();
    HinesAlgo(upper, lower, diag, rhs, parent, N);
    endTime = omp_get_wtime();
    printf("Execution time of serial Hines: %lf s\n", endTime - startTime);
    printf("backward sweep time: %lf s\n", t2 - t1);
    printf("forward sweep time: %lf s\n", t3 - t2);

    // write result to file
    fout = fopen(argv[2], "w");
    for (i = 0; i < N; i++)
        fprintf(fout, "%d %.8e %.8e %.8e %.8e\n",
                i, upper[i], lower[i], rhs[i], diag[i]);
    fclose(fout);

    delete[] upper;
    delete[] lower;
    delete[] diag;
    delete[] rhs;
    delete[] parent;
    return 0;
}
