#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define SIZE 100

struct point_t {
    double x;
    double y;
};

struct point_t circle(double r, double theta);

int main(void)
{
    FILE *shape = fopen("../datapoints/tao_distance/circle.dat", "w");
    struct point_t *point = malloc(sizeof(struct point_t) * SIZE);
    int i = 0;
    for(; i < SIZE; i++) {
        /* stores the point corresponding to i/32ths of a circle*/
        point[i] = circle(4, 2*i*M_PI/SIZE);
        fprintf(shape, "%d: (%lf, %lf)\n", i, point[i].x, point[i].y);
    }
    fclose(shape);
    return 0;
}

struct point_t circle(double r, double theta)
{
    struct point_t point;
    point.x = r*cos(theta);
    point.y = r*sin(theta);
    return point;
}
