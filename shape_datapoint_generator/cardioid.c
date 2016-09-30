#include <stdio.h>
#include <stdlib.h>
#include <math.h>

struct point_t {
    double x;
    double y;
};

struct point_t cardioid(double r, double theta);

int main(void)
{
    struct point_t *point = malloc(sizeof(struct point_t) * 32);
    int i = 0;
    for(; i < 32; i++) {
        /* stores the point corresponding to i/32ths of a circle*/
        point[i] = cardioid(4, 2*i*M_PI/32);
        printf("%d: (%lf, %lf)\n", i, point[i].x, point[i].y);
    }
    return 0;
}

struct point_t cardioid(double r, double theta)
{
    struct point_t point;
    point.x = r * cos(theta) * (1 - cos(theta));
    point.y = r * sin(theta) * (1 - cos(theta));
    return point;
}
