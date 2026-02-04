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
    int graphics = atoi(argv[5]);

    FILE *file;
    file = fopen(filename, "r");
    double buffer[6 * N];

    /*Read in the initial state*/
    fread(buffer, sizeof(double), 6 * N, file);
    fclose(file);

    double e0 = 0.001;
    double G = 100 / (double)N;

    for (int n = 0; n < 6 * N; n += 6) {
        /* printf("Element: %i", n); */
        /* printf("xpos: %f", buffer[n]); */
        /* printf("ypos: %f", buffer[n + 1]); */
        /* printf("mass: %f", buffer[n + 2]); */
        /* printf("xvel: %f", buffer[n + 3]); */
        /* printf("yvel: %f", buffer[n + 5]); */
        /* printf("brightness: %f", buffer[n + 6]); */
        /* Compute Fx and Fy values*/
        double Fx = 0, Fy = 0;
        for (int m = 0; m < 6 * N; m += 6) {
            if (m == n) break;
            double xdiff = buffer[n] - buffer[m];
            double ydiff = buffer[n + 1] - buffer[m + 1];
            Fx += -G * buffer[m + 2] / pow((xdiff + e0), 3) * xdiff;
            Fy += -G * buffer[m + 2] / pow((ydiff + e0), 3) * ydiff;
        }

        /*Update buffer values*/
        double ax = Fx / buffer[2];
        double ay = Fy / buffer[2];
        buffer[n + 3] += dt * ax;
        buffer[n + 4] += dt * ay;
        buffer[n] += dt * buffer[n + 3];
        buffer[n + 1] += dt * buffer[n + 4];

        /* printf("\n"); */
    }
    FILE *output = fopen("result.gal", "wb");
    fwrite(buffer, sizeof(double), 6 * N, output);
    fclose(output);
}
