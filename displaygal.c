#include "displaygal.h"

const float radius = 0.02, color = 0;

void display(double* buffer, int N){
    int L = 1, W = 1; // Skärmdimensioner
    ClearScreen();
    // Ritar en cirkel för varje objekt
    for (int i = 0; i < N; i += 6){
        DrawCircle(buffer[i], buffer[i + 1], L, W, radius, color);
    }
    Refresh();
    usleep(3000);
}