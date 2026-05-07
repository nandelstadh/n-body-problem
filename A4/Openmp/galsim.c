#include "../../displaygal.h"
#include "../../graphics.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "omp.h"

/* Accept arguments ./galsim N filename nsteps delta_t graphics */
/* where the input arguments have the following meaning: */
/* N: the number of stars/particles to simulate */
/* filename: the filename of the file to read the initial configuration from */
/* nsteps: the number of timesteps */
/* delta_t: the timestep ∆t */
/* graphics: 1 or 0 meaning graphics on/off */

int main(int argc, char* argv[]) {
    int N = atoi(argv[1]);
    char* filename = argv[2];
    int nsteps = atoi(argv[3]);
    double dt = atof(argv[4]);
    const char graphics = atoi(argv[5]);
    int threadcnt = atoi(argv[6]);

    omp_set_num_threads(threadcnt);

    const int windowWidth = 800;

    if (graphics) {
        InitializeGraphics(argv[5], windowWidth, windowWidth);
        SetCAxes(0, 1);
    }

    FILE* file = fopen(filename, "r");
    double* buffer = malloc(sizeof(double) * 6 * N);
    double* px = malloc(sizeof(double) * N);
    double* py = malloc(sizeof(double) * N);
    double* mass = malloc(sizeof(double) * N);
    double* vx = malloc(sizeof(double) * N);
    double* vy = malloc(sizeof(double) * N);
    double* brightness = malloc(sizeof(double) * N);
    double* axup = calloc((size_t)threadcnt * N, sizeof(double));
    double* ayup = calloc((size_t)threadcnt * N, sizeof(double));

    if (file == NULL || buffer == NULL || px == NULL || py == NULL || mass == NULL || vx == NULL || vy == NULL || brightness == NULL || axup == NULL || ayup == NULL) {
        perror("alloc/open");
        return 1;
    }

    /*Read in the initial state*/
    fread(buffer, sizeof(double), 6 * N, file);
    fclose(file);

#pragma omp parallel for
    for (int i = 0; i < N; i++) {
        px[i] = buffer[6 * i];
        py[i] = buffer[6 * i + 1];
        mass[i] = buffer[6 * i + 2];
        vx[i] = buffer[6 * i + 3];
        vy[i] = buffer[6 * i + 4];
        brightness[i] = buffer[6 * i + 5];
    }

    double e0 = 0.001;
    double G = 100 / (double)N;
    double* ax = calloc(N, sizeof(double));
    double* ay = calloc(N, sizeof(double));

#pragma omp parallel
    {
        int tid = omp_get_thread_num();
        double* ax_local = axup + ((size_t)tid * N);
        double* ay_local = ayup + ((size_t)tid * N);

        for (int time = 0; time < nsteps; time++) {
#pragma omp for
            for (int n = 0; n < N; n++) {
                double px_n = px[n];
                double py_n = py[n];
                double mass_n = mass[n];
                double axn = 0, ayn = 0;
                for (int m = n + 1; m < N; m++) {
                    double xdiff = px_n - px[m];
                    double ydiff = py_n - py[m];
                    double r = sqrt(xdiff * xdiff + ydiff * ydiff) + e0;
                    double r3 = r * r * r;
                    double force_x = -G * xdiff / r3;
                    double force_y = -G * ydiff / r3;
                    axn += mass[m] * force_x;
                    ayn += mass[m] * force_y;
                    ax_local[m] -= mass_n * force_x;
                    ay_local[m] -= mass_n * force_y;
                }
                ax_local[n] += axn;
                ay_local[n] += ayn;
            }

#pragma omp single
            {
                for (int n = 0; n < N; n++) {
                    double ax_total = 0;
                    double ay_total = 0;
                    for (int k = 0; k < threadcnt; k++) {
                        ax_total += axup[(size_t)k * N + n];
                        ay_total += ayup[(size_t)k * N + n];
                    }
                    ax[n] = ax_total;
                    ay[n] = ay_total;
                    vx[n] += dt * ax[n];
                    vy[n] += dt * ay[n];
                    px[n] += dt * vx[n];
                    py[n] += dt * vy[n];
                }

                if (graphics) {
                    for (int i = 0; i < N; i++) {
                        buffer[6 * i] = px[i];
                        buffer[6 * i + 1] = py[i];
                        buffer[6 * i + 2] = mass[i];
                        buffer[6 * i + 3] = vx[i];
                        buffer[6 * i + 4] = vy[i];
                        buffer[6 * i + 5] = brightness[i];
                    }
                    display(buffer, 6 * N);
                }

                memset(axup, 0, sizeof(double) * (size_t)threadcnt * N);
                memset(ayup, 0, sizeof(double) * (size_t)threadcnt * N);
            }
        }
    }

#pragma omp parallel for
    for (int i = 0; i < N; i++) {
        buffer[6 * i] = px[i];
        buffer[6 * i + 1] = py[i];
        buffer[6 * i + 2] = mass[i];
        buffer[6 * i + 3] = vx[i];
        buffer[6 * i + 4] = vy[i];
        buffer[6 * i + 5] = brightness[i];
    }
    FILE* output = fopen("result.gal", "wb");
    fwrite(buffer, sizeof(double), 6 * N, output);
    fclose(output);

    if (graphics) {
        FlushDisplay();
        CloseDisplay();
    }

    free(buffer);
    free(px);
    free(py);
    free(mass);
    free(vx);
    free(vy);
    free(brightness);
    free(ax);
    free(ay);
    free(axup);
    free(ayup);
}
