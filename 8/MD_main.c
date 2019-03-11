/*
 MD_main.c
 
 Created by Anders Lindman on 2013-10-31.
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <complex.h>
#include "initfcc.h"
#include "alpotential.h"

#define kB 8.6173303e-05       // [eV/K]
#define kappa 2.219            // A^3/eV
#define MASS 0.002796439       // eV ps^2/A^2
#define eV_A3_to_BAR 160217.66 // conversion factor bar/(eV/A^3)
#define K_for_0_C 273.15

void prog(double temp_eq);
void sys(int N, double pos[N][3], double vs[N][3], double *e_pot, double *e_kin, double *temps,
         double *pressures, double dt, int T, double L, double m);
double equilibriate(int N, double pos[N][3], double vs[N][3], double *e_pot, double *e_kin,
                    double *temps, double *pressures, double dt, int T, double L, double m,
                    double tau_T, double temp_eq, double tau_P, double press_eq);
void add_pair_distr(int N, double pos[N][3], double *g, int len_g, double delta_r, double L);
void write_array(char *file_name, double **arr, int nrows, int ncols, double dt, double t0);
void write_vector(char *file_name, double *vec, int n, double dt, double t0);
double **allocate2d(int nrows, int ncols);

/* Main program */
int main()
{
    printf("Structural\n");
    printf("T = 700 C\n");
    prog(700 + K_for_0_C);
}

void prog(double temp_eq)
{
    srand(time(NULL));
    double rando;

    double total_time = 20; // ps
    double dt = 0.01;
    int T = (int)(total_time / dt);

    int Nc = 4;
    double perturbation = 0.13;
    double v0 = 66; // volume of single cell: angström^3
    double m = MASS;

    double tau_T = 2;
    double tau_P = 2;
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
        L = equilibriate(N, pos, vs, e_pot_eq, e_kin_eq, temps_eq, pressures_eq,
                         dt, T, L, m, tau_T, temp_eq + 500, tau_P, press_eq * 20);
    // equil
    L = equilibriate(N, pos, vs, e_pot_eq, e_kin_eq, temps_eq, pressures_eq, dt, T,
                     L, m, tau_T, temp_eq, tau_P, press_eq);

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
    total_time = 1; // ps
    T = (int)(total_time / dt);

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

    // temperature
    double temp = 0;
    int start = (int)round(T * 0);
    for (int i = start; i < T; i++)
    {
        temp += temps[i];
    }
    temp /= (double)T - start;
    printf("Average temperature: %.4f K\n", temp);

    // pressure
    double press = 0;
    start = (int)round(T * 0);
    for (int i = start; i < T; i++)
    {
        press += pressures[i];
    }
    press /= (double)T - start;
    printf("Average pressure: %.4f bar\n", press);

    write_array("plote.dat", e, 3, T, dt, 0);
    write_array("plottp.dat", TP, 2, T, dt, 0);
}

void sys(int N, double pos[N][3], double vs[N][3], double *e_pot, double *e_kin,
         double *temps, double *pressures, double dt, int T, double L, double m)
{
    double f[N][3];
    double V = L * L * L;
    double Nc = N / 4;
    Nc = pow(Nc, 1.0 / 3.0);
    double a0 = L / Nc;

    double temp;
    double press;
    double virial;
    double pos0[N][3];
    for (int i = 0; i < N; i++)
        for (int j = 0; j < 3; j++)
            pos0[i][j] = pos[i][j];

    double q_grid_length = 2 * M_PI / L;
    int q_grid_size_half = 32;
    int q_grid_size = 2 * q_grid_size_half + 1;
    double S_q_grid[q_grid_size][q_grid_size][q_grid_size];

    get_forces_AL(f, pos, L, N);

    e_pot[0] = get_energy_AL(pos, L, N);
    e_kin[0] = get_kin_energy_AL(vs, N, m);

    for (int t = 1; t < T; t++)
    {
        printf("%.4f\n", (double)t / T);
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

        // pressure
        virial = get_virial_AL(pos, L, N);
        press = (N * kB * temp + virial) / V;
        pressures[t] = press * eV_A3_to_BAR;

        // Finding S(q) grid
        double qx, qy, qz, q_dot_r;
        double _Complex sum;
        for (int i = 0; i < q_grid_size; i++)
        {
            for (int j = 0; j < q_grid_size; j++)
            {
                for (int k = 0; k < q_grid_size; k++)
                {
                    qx = q_grid_length * (i - q_grid_size_half);
                    qy = q_grid_length * (j - q_grid_size_half);
                    qz = q_grid_length * (k - q_grid_size_half);

                    sum = 0;
                    for (int ind = 0; ind < N; ind++)
                    {
                        q_dot_r = qx * pos[ind][0] + qy * pos[ind][1] + qz * pos[ind][2];
                        sum += cexp(I * q_dot_r);
                    }

                    S_q_grid[i][j][k] += (pow(creal(sum), 2) + pow(cimag(sum), 2));
                }
            }
        }
    }

    // spherical averaging
    double qx, qy, qz, q_norm;
    int q_index;
    double max_q = 20;
    double delta_q = 0.1;
    int S_q_length = (int)(max_q / delta_q);
    double S_q[S_q_length];
    double S_q_counter[S_q_length]; // to normalize in the end

    for (int i = 0; i < S_q_length; i++)
    {
        S_q[i] = 0;
        S_q_counter[i] = 0;
    }

    for (int i = 0; i < q_grid_size; i++)
    {
        for (int j = 0; j < q_grid_size; j++)
        {
            for (int k = 0; k < q_grid_size; k++)
            {
                qx = q_grid_length * (i - q_grid_size_half);
                qy = q_grid_length * (j - q_grid_size_half);
                qz = q_grid_length * (k - q_grid_size_half);
                if (qx != 0.0 && qy != 0.0 && qz != 0.0)
                {
                    q_norm = sqrt(qx * qx + qy * qy + qz * qz);
                    q_index = (int)floor(q_norm / delta_q);
                    if (q_index < S_q_length)
                    {
                        S_q[q_index] += S_q_grid[i][j][k];
                        S_q_counter[q_index] += 1;
                    }
                }
            }
        }
    }

    // normalize
    for (int i = 0; i < S_q_length; i++)
    {
        if (S_q_counter[i] > 0.0)
            S_q[i] /= S_q_counter[i] * N * T;
    }

    write_vector("sq.dat", S_q, S_q_length, delta_q, 0);
}

double equilibriate(int N, double pos[N][3], double vs[N][3], double *e_pot, double *e_kin,
                    double *temps, double *pressures, double dt, int T, double L, double m,
                    double tau_T, double temp_eq, double tau_P, double press_eq)
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
        alpha_P = 1 - kappa * dt / tau_P * (press_eq - press);
        double alpha_P_1_3 = pow(alpha_P, 1.0 / 3.0);
        L *= alpha_P_1_3;
        for (int i = 0; i < N; i++)
            for (int j = 0; j < 3; j++)
                pos[i][j] *= alpha_P_1_3;
        V = L * L * L;
        pressures[t] = press * eV_A3_to_BAR;
        //printf("%.8f\n",V);
    }
    return L;
}

void add_pair_distr(int N, double pos[N][3], double *g, int len_g, double delta_r, double L)
{
    double delta; //temp variable
    double dist;
    int r_index;

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            if (i != j)
            {
                dist = 0;
                for (int d = 0; d < 3; d++)
                {
                    delta = pos[i][d] - pos[j][d];
                    delta = delta - L * round(delta / L);
                    dist += pow(delta, 2);
                }
                dist = sqrt(dist);

                r_index = (int)floor(dist / delta_r);
                if (r_index < len_g)
                    g[r_index] += 1;
            }
        }
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

void write_vector(char *file_name, double *vec, int n, double dt, double t0)
{

    FILE *fp;
    fp = fopen(file_name, "w");
    for (int i = 0; i < n; ++i)
    {
        fprintf(fp, "%.6f ", i * dt + t0);
        fprintf(fp, "%.6f ", vec[i]);
        fprintf(fp, "\n");
    }
    fclose(fp);
}