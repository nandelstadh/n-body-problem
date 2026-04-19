#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

int main(int argc, char *argv[]) {
    if (argc != 7) {
        printf("usage ./galsim N filename nsteps delta_t graphics n_threads");
        return 1;
    }
    int N = atoi(argv[1]);
    char *filename = argv[2];
    int nsteps = atoi(argv[3]);
    double dt = atof(argv[4]);
    const int nthreads = atoi(argv[6]);
    FILE *file;
    file = fopen(filename, "rb");
    double buffer[6 * N];
    fread(buffer, sizeof(double), 6 * N, file);
    fclose(file);

    double px[N], py[N], mass[N], vx[N], vy[N], brightness[N], ax[N], ay[N];
    double *axnup[nthreads], *aynup[nthreads];
    for (int k=0; k< nthreads;k++) {  
        axnup[k] = calloc(N,sizeof(double)); 
        aynup[k] = calloc(N, sizeof(double));
    }  
    for (int i = 0; i < N; i++) {
        px[i] = buffer[6 * i];
        py[i] = buffer[6 * i + 1];
        mass[i] = buffer[6 * i + 2];
        vx[i] = buffer[6 * i + 3];
        vy[i] = buffer[6 * i + 4];
        brightness[i] = buffer[6 * i + 5];
        ax[i] = ay[i] = 0;
    }
    #pragma omp parallel num_threads(nthreads)
    {
    int id = omp_get_thread_num();
    for (int time = 0; time < nsteps; time++) { 
        #pragma omp for schedule(dynamic)
        for (int n = 0; n < N; n++) {        
            double px_n = px[n];
            double mass_n = mass[n];
            double py_n = py[n];
            double axn = 0, ayn = 0;
            for (int m = n + 1; m < N; m++) {
                double xdiff = px_n-px[m];
                double ydiff = py_n-py[m];
                double e0 = 0.001;
                double G = 100 / (double)N;
                double r = sqrt(xdiff * xdiff + ydiff * ydiff) + e0;
                double r3 = r * r * r;
                double force_x = -G * xdiff / r3;
                double force_y = -G * ydiff / r3;
                axn += mass[m] * force_x;
                ayn += mass[m] * force_y;
                axnup[id][m] -= mass_n * force_x;
                aynup[id][m] -= mass_n * force_y;
            }
            axnup[id][n] += axn; 
            aynup[id][n] += ayn;
        }
        #pragma omp single
        {
        for (int n=0;n<N;n++){ 
            for(int k=0;  k<nthreads; k++) { 
                ax[n] += axnup[k][n];
                ay[n] += aynup[k][n];
            }
        }
        for (int n  =0; n<N; n++) { 
            vx[n] += dt*ax[n];
            vy[n] += dt*ay[n]; 
            px[n] += dt*vx[n];
            py[n] += dt*vy[n];
        }
        for (int k =0;  k< nthreads; k++) { 
            memset(axnup[k], 0, N*sizeof(double));
            memset(aynup[k], 0, N*sizeof(double));
        }
        memset(ax, 0, N*sizeof(double));
        memset(ay, 0, N*sizeof(double));
        }
        }  
    }
    for (int i = 0; i < N; i++) {
        buffer[6 * i]      = px[i];
        buffer[6 * i + 1] = py[i];
        buffer[6 * i + 2] = mass[i];
        buffer[6 * i + 3] = vx[i];
        buffer[6 * i + 4] = vy[i];
        buffer[6 * i + 5] = brightness[i];
    }
    for (int k  =0;  k <nthreads; k ++) { 
        free(axnup[k]);
        free(aynup[k]);
    }
    FILE *output = fopen("result.gal", "wb");
    if (output == NULL) {
        perror("fopen");
        return 1;
    }
    fwrite(buffer, sizeof(double), 6 * N, output);
    fclose(output);
    return 0;
}
