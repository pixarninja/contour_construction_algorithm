#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define SIZE 100

struct point_t {
    double x;
    double y;
};

struct point_t cardioid(double r, double theta);

int main(void)
{
    FILE *shape = fopen("../datapoints/tao_distance/rose.dat", "w");
    struct point_t *point = malloc(sizeof(struct point_t) * SIZE);
    int i = 0;
    for(; i < 32; i++) {
        /* stores the point corresponding to i/32ths of a circle*/
        point[i] = cardioid(1, 2*i*M_PI/SIZE);
        fprintf(shape, "%d: (%lf, %lf)\n", i, point[i].x, point[i].y);
    }
    fclose(shape);
    return 0;
}

struct point_t cardioid(double r, double theta)
{
    struct point_t point;
    point.x = cos(5 * theta) * cos(theta);
    point.y = sin(5 * theta) * sin(theta);
    return point;
}
