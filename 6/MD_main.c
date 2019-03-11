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

#define kB 8.6173303e-05       // [eV/K]
#define kappa 2.219            // A^3/eV
#define MASS 0.002796439       // eV ps^2/A^2
#define eV_A3_to_BAR 1602176.6 // conversion factor
#define K_for_0_C 273.15

typedef struct
{
    double a;
    double b;
} tuple;

tuple prog(double temp_eq, double delta_temp);
void sys(int N, double pos[N][3], double vs[N][3], double *e_pot, double *e_kin, double *temps, double *pressures, double dt, int T, double L, double m);
double equilibriate(int N, double pos[N][3], double vs[N][3], double *e_pot, double *e_kin, double *temps, double *pressures, double dt, int T, double L, double m, double tau_T, double temp_eq, double tau_P, double press_eq, int press_flag);
void write_array(char *file_name, double **arr, int nrows, int ncols, double dt, double t0);
double **allocate2d(int nrows, int ncols);

/* Main program */
int main()
{
    printf("Specific heat of Al: 0.0904 å^2/(ps^2*K)\n");
    double delta_temp = 5.0;
    printf("T = %.2f C\n", 500 - delta_temp);
    tuple r1 = prog(500 + K_for_0_C, -delta_temp);
    double e1 = r1.a;
    double t1 = r1.b;
    printf("T = %.2f C\n", 500 + delta_temp);
    tuple r2 = prog(500 + K_for_0_C, + delta_temp);
    double e2 = r2.a;
    double t2 = r2.b;
    double C_V;
    C_V = (e2 - e1) / (t2 - t1);
    printf("C_V = %.4f \n", C_V);

    printf("Specific heat of Al liquid: 0.118 å^2/(ps^2*K)\n");
    printf("T = %.2f C\n", 700 - delta_temp);
    r1 = prog(700 + K_for_0_C, -delta_temp);
    e1 = r1.a;
    t1 = r1.b;
    printf("T = %.2f C\n", 700 + delta_temp);
    r2 = prog(700 + K_for_0_C, delta_temp);
    e2 = r2.a;
    t2 = r2.b;
    C_V = (e2 - e1) / (t2 - t1);
    printf("C_V = %.4f \n", C_V);
}

//returns mean energy
tuple prog(double temp_eq, double delta_temp)
{
    srand(time(NULL));
    double rando;

    double total_time = 30; // ps
    double dt = 0.005;
    int T = (int)total_time / dt;

    int Nc = 4;
    double perturbation = 0.13;
    double v0 = 67.18; // volume of single cell: angström^3
    double m = MASS;

    double tau_T = 4;
    double tau_P = 4;
    double press_eq = 1.0 / eV_A3_to_BAR;

    int N = 4 * Nc * Nc * Nc;
    double a0 = pow(v0, 1.0 / 3.0);
    double L = a0 * Nc;
    double pos[N][3];
    double vs[N][3];

    double e_pot_eq[T];
    double e_kin_eq[T];
    double **e_eq = allocate2d(3, T);
    e_eq[0] = e_kin_eq;
    e_eq[1] = e_pot_eq;

    double temps_eq[T];
    double pressures_eq[T];
    double **TP_eq = allocate2d(2, T);
    TP_eq[0] = temps_eq;
    TP_eq[1] = pressures_eq;

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

    // #### equilibriate ####
    // melt (only ex. 4)
    if (temp_eq > 600.0 + K_for_0_C)
        L = equilibriate(N, pos, vs, e_pot_eq, e_kin_eq, temps_eq, pressures_eq, dt, T, L, m, tau_T, temp_eq + 1000, tau_P, press_eq, 1);
    // equil
    L = equilibriate(N, pos, vs, e_pot_eq, e_kin_eq, temps_eq, pressures_eq, dt, T, L, m, tau_T, temp_eq, tau_P, press_eq, 1);
    L = equilibriate(N, pos, vs, e_pot_eq, e_kin_eq, temps_eq, pressures_eq, dt, T, L, m, tau_T, temp_eq + delta_temp, tau_P, press_eq, 0);

    for (int i = 0; i < T; i++)
    {
        e_eq[2][i] = e_eq[0][i] + e_eq[1][i];
    }

    write_array("plote_eq.dat", e_eq, 3, T, dt, 0);
    write_array("plottp_eq.dat", TP_eq, 2, T, dt, 0);

    printf("Equilibration done.\nV: %.4f\n", L * L * L);

    free(e_eq);
    free(TP_eq);

    // #### velocity verlet ####
    total_time = 120; // ps
    T = (int)total_time / dt;

    double e_pot[T];
    double e_kin[T];
    double **e = allocate2d(3, T);
    e[0] = e_kin;
    e[1] = e_pot;

    double temps[T];
    double pressures[T];
    double **TP = allocate2d(2, T);
    TP[0] = temps;
    TP[1] = pressures;

    // run system
    sys(N, pos, vs, e_pot, e_kin, temps, pressures, dt, T, L, m);

    double mean_e = 0;
    for (int i = 0; i < T; i++)
    {
        mean_e += e_kin[i] + e_pot[i];
    }
    mean_e /= T;

    double temp = 0;
    for (int i = 0; i < T; i++)
    {
        temp += temps[i];
    }
    temp /= T;
    printf("Mean T = %.4f K\n", temp);

    tuple t;
    t.a = mean_e;
    t.b = temp;
    return t;
}

void sys(int N, double pos[N][3], double vs[N][3], double *e_pot, double *e_kin, double *temps, double *pressures, double dt, int T, double L, double m)
{
    double f[N][3];
    double V = L * L * L;

    double temp;
    double press;
    double virial;

    double pos0[N][3];
    for (int i = 0; i < N; i++)
        for (int j = 0; j < 3; j++)
            pos0[i][j] = pos[i][j];

    get_forces_AL(f, pos, L, N);

    e_pot[0] = get_energy_AL(pos, L, N);
    e_kin[0] = get_kin_energy_AL(vs, N, m);

    for (int t = 1; t < T; t++)
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

        // temperature
        temp = e_kin[t] * 2.0 / (3.0 * N * kB);
        temps[t] = temp;
    }
}

double equilibriate(int N, double pos[N][3], double vs[N][3], double *e_pot, double *e_kin, double *temps, double *pressures, double dt, int T, double L, double m, double tau_T, double temp_eq, double tau_P, double press_eq, int press_flag)
{
    double f[N][3];
    double V = L * L * L;

    double alpha_T;
    double temp;
    double alpha_P;
    double press;
    double virial;

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

        // scale temperature
        temp = e_kin[t] * 2.0 / (3.0 * N * kB);
        alpha_T = 1 + 2 * dt / tau_T * (temp_eq - temp) / temp;
        double alpha_T_sqrt = sqrt(alpha_T);
        for (int i = 0; i < N; i++)
            for (int j = 0; j < 3; j++)
                vs[i][j] *= alpha_T_sqrt;
        temps[t] = temp;

        // scale pressure
        virial = get_virial_AL(pos, L, N);
        press = (N * kB * temp + virial) / V;
        if (press_flag == 1)
        {
            alpha_P = 1 - kappa * dt / tau_P * (press_eq - press);
            double alpha_P_1_3 = pow(alpha_P, 1.0 / 3.0);
            L *= alpha_P_1_3;
            for (int i = 0; i < N; i++)
                for (int j = 0; j < 3; j++)
                    pos[i][j] *= alpha_P_1_3;
        }
        V = L * L * L;
        pressures[t] = press * eV_A3_to_BAR;
        //printf("%.8f\n",V);
    }
    return L;
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