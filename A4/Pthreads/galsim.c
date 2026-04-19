#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
typedef struct { 
    int id; 
    int nthreads;
    int nsteps;
    double dt; 
    int N;
    double* px; 
    double* py; 
    double* mass; 
    double* vx; 
    double* vy;
    double* axn; 
    double* ayn; 
    double** axnup; 
    double** aynup;
    pthread_barrier_t *pbarr;
} GalParams;



void* pairloop(void* args) { 
    GalParams params = *(GalParams*)(args); 
    int N = params.N;
    double* axnup = params.axnup[params.id]; 
    double* aynup = params.aynup[params.id]; 
    for (int time = 0; time < params.nsteps; time++) { 
        for (int n =  params.id; n < N; n+=params.nthreads) {        
            double px_n = params.px[n];
            double mass_n = params.mass[n];
            double py_n = params.py[n];
            double axn = 0, ayn = 0;
            for (int m = n + 1; m < N; m++) {
                double xdiff = px_n-params.px[m];
                double ydiff = py_n-params.py[m];
                double e0 = 0.001;
                double G = 100 / (double)N;
                double r = sqrt(xdiff * xdiff + ydiff * ydiff) + e0;
                double r3 = r * r * r;
                double force_x = -G * xdiff / r3;
                double force_y = -G * ydiff / r3;
                axn += params.mass[m] * force_x;
                ayn += params.mass[m] * force_y;
                axnup[m] -= mass_n * force_x;
                aynup[m] -= mass_n * force_y;
            }
            axnup[n] += axn; 
            aynup[n] += ayn;
        }
        pthread_barrier_wait(params.pbarr);
        if (params.id == 0) { 
            for (int n=0;n<N;n++){ 
                for(int k=0;  k<params.nthreads; k++) { 
                    params.axn[n] += params.axnup[k][n];
                    params.ayn[n] += params.aynup[k][n];
                }
            }
            for (int n  =0; n<N; n++) { 
                params.vx[n] += params.dt*params.axn[n];
                params.vy[n] += params.dt*params.ayn[n]; 
                params.px[n] += params.dt*params.vx[n];
                params.py[n] += params.dt*params.vy[n];
            }
            for (int k =0;  k< params.nthreads; k++) { 
                memset(params.axnup[k], 0, params.N*sizeof(double));
                memset(params.aynup[k], 0, params.N*sizeof(double));
            }
            memset(params.axn, 0, params.N*sizeof(double));
            memset(params.ayn, 0, params.N*sizeof(double));
        }
        pthread_barrier_wait(params.pbarr);
    }  
    return NULL; 
}

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
    GalParams args[nthreads];
    for (int i = 0; i < N; i++) {
        px[i] = buffer[6 * i];
        py[i] = buffer[6 * i + 1];
        mass[i] = buffer[6 * i + 2];
        vx[i] = buffer[6 * i + 3];
        vy[i] = buffer[6 * i + 4];
        brightness[i] = buffer[6 * i + 5];
        ax[i] = ay[i] = 0;
    }
    pthread_barrier_t pbarr;  
    pthread_barrier_init(&pbarr, NULL, nthreads);  
    for (int k=0; k<nthreads; k++) { 
        args[k].id = k; 
        args[k].nthreads = nthreads;
        args[k].nsteps = nsteps;
        args[k].dt = dt;  
        args[k].px = px; 
        args[k].py = py; 
        args[k].pbarr =  &pbarr;
        args[k].N = N;
        args[k].mass = mass; 
        args[k].axn = ax; 
        args[k].ayn = ay; 
        args[k].vx = vx;  
        args[k].vy = vy;
        args[k].axnup = axnup;
        args[k].aynup = aynup;
    }
    pthread_t threads[nthreads];
    for (int k  =0 ; k <nthreads; k ++) {  
        pthread_create(&threads[k], NULL, pairloop, &args[k]);
    }
    for (int  k=0;  k<nthreads;k++) { 
        pthread_join(threads[k], NULL);  
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
        pthread_barrier_destroy(&pbarr);
        return 1;
    }
    fwrite(buffer, sizeof(double), 6 * N, output);
    fclose(output);
    pthread_barrier_destroy(&pbarr); 
    return 0;
}
