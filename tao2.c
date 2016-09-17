#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#define SIZE 7

struct point_t {
    double x;
    double y;
    int index;
};

struct vector_t {
    struct point_t point[2];
    double length;
    double i;
    double j;
};

struct function_t {
    struct point_t point;
    double x;
    double y;
    double m;
    double b;
};

struct triangle_t {
    char name[3];
    char *points;
    struct function_t f1;
    struct function_t f2;
    struct function_t f3;
    double x;
    double y;
};

int global_count = 0;

double shortest_path(struct point_t *point, int size);
int **construct_curve(struct point_t *point, int size);
double calculate_theta(struct vector_t V1, struct vector_t V2);
struct triangle_t *construct_triangles(struct point_t **curve, int size,int divisions);
double distance_p(struct point_t start, struct point_t end);
double distance_v(struct vector_t start, struct vector_t end);
double length_v(struct vector_t v);
double dot_product(struct vector_t start, struct vector_t end);

int main(void)
{
    struct point_t *point = malloc(sizeof(struct point_t) * SIZE);
    int size = SIZE - 1;
    double distance = 0;
    /* O(0,0) */
    point[6].x = 0;
    point[6].y = 0;
    point[6].index = 6;
    /* A(0,2) */
    point[3].x = 0;
    point[3].y = 2;
    point[3].index = 3;
    /* B(2,2) */
    point[5].x = 2;
    point[5].y = 2;
    point[5].index = 5;
    /* C(2,4) */
    point[4].x = 2;
    point[4].y = 4;
    point[4].index = 4;
    /* D(-1,2) */
    point[2].x = -1;
    point[2].y = 2;
    point[2].index = 2;
    /* E(-2,-2) */
    point[1].x = -2;
    point[1].y = -2;
    point[1].index = 1;
    /* F(0,-3) */
    point[0].x = 0;
    point[0].y = -3;
    point[0].index = 0;
    /* runs tao-distance algorithm on data */
    distance = shortest_path(point, size);
    printf("\n");
    printf("Total Permutations: %d\n\n", global_count);
    printf("Shortest Path: %s\n", "derp");
    printf("Distance: %lf\n", distance);
    printf("\n");
    return 0;
}

/* calculates the shortest path */
double shortest_path(struct point_t *point, int size)
{
    double total = size;
    return total;
}

/* returns a 2D array of (0, 1, or 2) neighbors for each point */
int **construct_curve(struct point_t *point, int size)
{
    int *curve = malloc(sizeof(int) * size);
    double range = 0;
    struct point_t center;
    double shortest = DBL_MAX;
    double longest = 0;
    double tmp = 0;
    double tao = 0;
    double division_constant = 0;
    int i = 0;
    int sum_x = 0;
    int sum_y = 0;
    /* calculate average point */
    for(i = 0; i <= size; i++) {
        sum_x += point[i].x;
        sum_y += point[i].y;
    }
    center.x = sum_x / size;
    center.y = sum_y / size;
    /* find shortest and longest distances */
    for(i = 0; i <= size; i++) {
        tmp = distance_p(center, point[i]);
        if(tmp < shortest) {
            shortest = tmp;
        }
        if(tmp > longest) {
            longest = tmp;
        }
    }
    /* return the difference in length of the position vectors */
    tao = longest - shortest;
    division_constant = (tao/((int)sqrt(tao) + 1));
    /* find the points in the current range */
    range = division_constant;
    return &curve;
}

/* constructs triangles given the curve they are on*/
struct triangle_t *construct_triangles(struct point_t **curve, int size,int divisions)
{
    struct triangle_t *triangle = malloc(sizeof(struct triangle_t) * 100);
    return triangle;
}

/* calculates theta given two vectors */
double calculate_theta(struct vector_t V1, struct vector_t V2)
{
    return acos(dot_product(V1, V2)/(V1.length * V2.length));
}

/* calculates distance given two points */
double distance_p(struct point_t start, struct point_t end)
{
    double x_1 = start.x;
    double y_1 = start.y;
    double x_2 = end.x;
    double y_2 = end.y;
    return sqrt((x_2 - x_1) * (x_2 - x_1) + (y_2 - y_1) * (y_2 - y_1));
}

/* calculates distance given two vectors */
double distance_v(struct vector_t start, struct vector_t end)
{
    double i_1 = start.i;
    double j_1 = start.j;
    double i_2 = end.i;
    double j_2 = end.j;
    return sqrt((i_2 - i_1) * (i_2 - i_1) + (j_2 - j_1) * (j_2 - j_1));
}

/* calculates length of a single vectors */
double length_v(struct vector_t v)
{
    double i = v.i;
    double j = v.j;
    return sqrt((i * i) + (j * j));
}

/* calculates dot product of two vectors */
double dot_product(struct vector_t start, struct vector_t end)
{
    double i_1 = start.i;
    double j_1 = start.j;
    double i_2 = end.i;
    double j_2 = end.j;
    return ((i_2 * i_1) + (j_2 * j_1));
}
