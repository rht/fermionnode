#include <stdio.h>
#include <math.h>

float phi2D(int nx, int ny, float length, float x, float y)
{
    return sqrt(2. / length) * cos(nx * 2 * pi() * x / length) *
           cos(ny * 2 * pi() * y / length);
}

