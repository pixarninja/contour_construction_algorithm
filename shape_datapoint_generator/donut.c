#include <stdio.h>
#include <stdlib.h>
#include <math.h>

struct point_t {
    double x;
    double y;
};

struct point_t circle(double r, double theta);

int main(int argc, char *argv[])
{
    int size;
    int size_a;
    FILE *shape = fopen("./datapoints/tao_distance/donut.dat", "w");
    if(shape == NULL) {
        printf("File not found. Exiting Program. Good Day.\n");
        exit(EXIT_FAILURE);
    }
    struct point_t *point;
    int i = 0;
    if(argc == 1) {
        size = 100;
    }
    else {
        size = atoi(argv[argc - 1]);
    }
    if((size % 2) == 0) {
        size_a = size / 2;
    }
    else {
        size_a = (size / 2) + 1;
    }
    point = malloc(sizeof(struct point_t) * size);
    for(; i < size_a; i++) {
        point[i] = circle(8, 2*(i)*M_PI/size_a);
        fprintf(shape, "%d: (%lf, %lf)\n", i, point[i].x, point[i].y);
    }
    for(; i < size; i++) {
        point[i] = circle(4, 2*(i - size_a + 1)*M_PI/size_a);
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
