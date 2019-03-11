/*
 MD_main.c
 
 Created by Anders Lindman on 2013-10-31.
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "initfcc.h"
#include "alpotential.h"

#define kB 8.6173303e-05 // [eV/K]
#define MASS 0.002796439 // eV ps^2/A^2

void sys(int N, double pos[N][3], double *e_pot, double *e_kin, double dt, int T, double L, double m);
void write_array(char *file_name, double **arr, int nrows, int ncols, double dt, double t0);
double **allocate2d(int nrows, int ncols);

/* Main program */
int main()
{

    srand(time(NULL));
    double rando;

    int T = 100000;
    double dt = 0.005; 

    int Nc = 4;
    double perturbation = 0.13;
    double v0 = 63;
    double m = MASS;

    int N = 4 * Nc * Nc * Nc;
    double a0 = pow(v0, 1.0 / 3.0);
    double L = Nc * a0;

    double pos[N][3];

    double e_pot[T];
    double e_kin[T];

    double **e = allocate2d(3, T);
    e[0] = e_kin;
    e[1] = e_pot;

    init_fcc(pos, Nc, a0);

    // perturb system
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            rando = (double)rand() / (double)RAND_MAX;
            pos[i][j] += perturbation * a0 * (rando - 0.5);
        }
    }

    // velovity-verlet
    sys(N, pos, e_pot, e_kin, dt, T, L, m);

    
    for(int i = 0; i < T; i++)
    {
        e[2][i] = e[0][i] + e[1][i];
    }

    write_array("plote.dat", e, 3, T, dt, 0);

    // temperature
    double temp = 0;
    int start = (int) round(T*0.5);
    for(int i = start; i < T; i++)
    {
        temp += e_kin[i]*2.0/(3.0*N*kB);
    }
    temp /= (double) T - start;
    printf("Average temperature: %.4fK\n", temp);
}

void sys(int N, double pos[N][3], double *e_pot, double *e_kin, double dt, int T, double L, double m)
{
    double vs[N][3];
    double f[N][3];

    for (int i = 0; i < N; i++)
        for (int j = 0; j < 3; j++)
            vs[i][j] = 0;

    get_forces_AL(f, pos, L, N);

    for (int t = 0; t < T; t++)
    {
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                vs[i][j] += 0.5 * f[i][j] * dt / m;
                pos[i][j] += vs[i][j] * dt;
            }
        }
        get_forces_AL(f, pos, L, N);
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                vs[i][j] += 0.5 * f[i][j] * dt / m;
            }
        }
        e_pot[t] = get_energy_AL(pos, L, N);
        e_kin[t] = get_kin_energy_AL(vs, N, m);
    }
}


double **allocate2d(int nrows, int ncols)
{
    double **res;
    const size_t pointers = nrows * sizeof *res;
    const size_t elements = nrows * ncols * sizeof **res;
    res = malloc(pointers + elements);

    size_t i;
    double *const data = (double *)&res[0] + nrows;
    for (i = 0; i < nrows; i++)
        res[i] = data + i * ncols;

    return res;
}

void write_array(char *file_name, double **arr, int nrows, int ncols, double dt, double t0)
{

    FILE *fp;
    fp = fopen(file_name, "w");
    for (int i = 0; i < ncols; ++i)
    {
        fprintf(fp, "%.6f ", i * dt + t0);
        for (int j = 0; j < nrows; j++)
        {
            fprintf(fp, "%.6f ", arr[j][i]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}