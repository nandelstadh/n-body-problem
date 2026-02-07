#include "displaygal.h"
#include "graphics.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"

/* Accept arguments ./galsim N filename nsteps delta_t graphics */
/* where the input arguments have the following meaning: */
/* N: the number of stars/particles to simulate */
/* filename: the filename of the file to read the initial configuration from */
/* nsteps: the number of timesteps */
/* delta_t: the timestep ∆t */
/* graphics: 1 or 0 meaning graphics on/off */

int main(int argc, char *argv[]) {
    int N = atoi(argv[1]);
    char *filename = argv[2];
    int nsteps = atoi(argv[3]);
    double dt = atof(argv[4]);
    const char graphics = atoi(argv[5]);

    const int windowWidth = 800;

    if (graphics) {
        InitializeGraphics(argv[5], windowWidth, windowWidth);
        SetCAxes(0, 1);
    }

    FILE *file;
    file = fopen(filename, "r");
    double buffer[6 * N];

    /*Read in the initial state*/
    fread(buffer, sizeof(double), 6 * N, file);
    fclose(file);

    double e0 = 0.001;
    double G = 100 / (double)N;

    for (int time = 0; time < nsteps; time++) {
        // printf("Timestep: %d \n", time);
        for (int n = 0; n < 6 * N; n += 6) {
            /* Compute ax and ay values*/
            double ax = 0, ay = 0;
            for (int m = 0; m < 6 * N; m += 6) {
                if (m == n) continue;
                double xdiff = buffer[n] - buffer[m];
                double ydiff = buffer[n + 1] - buffer[m + 1];
                double diff = sqrt(xdiff * xdiff + ydiff * ydiff);
                ax += -G * (buffer[m + 2] / pow((diff + e0), 3)) * xdiff;
                ay += -G * (buffer[m + 2] / pow((diff + e0), 3)) * ydiff;
            }

            /*Update buffer velocities*/
            buffer[n + 3] += dt * ax;
            buffer[n + 4] += dt * ay;
        }
        /*Update buffer positions for next time step*/
        for (int n = 0; n < 6 * N; n += 6) {
            buffer[n] += dt * buffer[n + 3];
            buffer[n + 1] += dt * buffer[n + 4];
        }
        if (graphics) {
            display(buffer, 6 * N);
        }
    }
    FILE *output = fopen("result.gal", "wb");
    fwrite(buffer, sizeof(double), 6 * N, output);
    fclose(output);

    if (graphics) {
        FlushDisplay();
        CloseDisplay();
    }
}
