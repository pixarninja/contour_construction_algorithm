#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#define SIZE 14

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
double tao_distance(double theta, double length, double curvature);
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
    /* A(0,1) */
    point[0].x = 0;
    point[0].y = 1;
    point[0].index = 0;
    /* B(-1,-1) */
    point[1].x = -1;
    point[1].y = -1;
    point[1].index = 1;
    /* C(1,-1) */
    point[2].x = 1;
    point[2].y = -1;
    point[2].index = 2;
    /* D(0,-4) */
    point[3].x = 0;
    point[3].y = -4;
    point[3].index = 3;
    /* E(4,0) */
    point[4].x = 4;
    point[4].y = 0;
    point[4].index = 4;
    /* F(3,3) */
    point[5].x = 3;
    point[5].y = 3;
    point[5].index = 5;
    /* G(-3,3) */
    point[6].x = -3;
    point[6].y = 3;
    point[6].index = 6;
    /* H(-4,0) */
    point[7].x = -4;
    point[7].y = 0;
    point[7].index = 7;
    /* I(-7,0) */
    point[8].x = -7;
    point[8].y = 0;
    point[8].index = 8;
    /* J(-6,6) */
    point[9].x = -6;
    point[9].y = 6;
    point[9].index = 9;
    /* K(0,7) */
    point[10].x = 0;
    point[10].y = 7;
    point[10].index = 10;
    /* L(7,0) */
    point[11].x = 7;
    point[11].y = 0;
    point[11].index = 11;
    /* M(5,-5) */
    point[12].x = 5;
    point[12].y = -5;
    point[12].index = 12;
    /* N(0,-7) */
    point[13].x = 0;
    point[13].y = -7;
    point[13].index = 13;
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
    int **curve = construct_curve(point, size);
    return total;
}

/* creates a 2D array of curves from the point data */
int **construct_curve(struct point_t *point, int size)
{
    struct point_t center;
    double *distance = malloc(sizeof(double) * size);
    double range[2] = {0};
    double shortest = DBL_MAX;
    double longest = 0;
    double tao = 0;
    double division_constant = 0;
    int **curve = malloc(sizeof(int *) * size);
    int **neighbors = malloc(sizeof(int *) * size);
    int shortest_path[2];
    int i = 0;
    int j = 0;
    int k = 0;
    int tmp = 0;
    int sum_x = 0;
    int sum_y = 0;
    int division_number = 0;
    /* calculate average point */
    for(i = 0; i <= size; i++) {
        sum_x += point[i].x;
        sum_y += point[i].y;
    }
    center.x = sum_x / size;
    center.y = sum_y / size;
    /* find shortest and longest distances */
    for(i = 0; i <= size; i++) {
        distance[i] = distance_p(center, point[i]);
        if(distance[i] < shortest) {
            shortest = distance[i];
        }
        if(distance[i] > longest) {
            longest = distance[i];
        }
    }
    /* return the difference in length of the position vectors */
    tao = longest - shortest;
    division_number = (int)sqrt(tao) + 1;
    division_constant = (tao/division_number);
    /* find the points in the current range */
    range[0] = shortest;
    for(i = 0; i < division_number; i++) {
        range[1] = (division_constant * (i + 1)) + shortest;
        /* initialize curve */
        curve[i] = malloc(sizeof(int) * size);
        for(j = 0; j <= size; j++) {
            curve[i][j] = SIZE;
        }
        for(j = 0; j <= size; j++) {
            /* range[0] < distance[j] < range[1] */
            if((distance[j] >= range[0]) && (distance[j] <= range[1])) {
                curve[i][j] = j;
                //k++;
            }
        }
        range[0] = range[1];
        k = 0;
    }
    /* printing for debug */
    for(i = 0; i < division_number; i++) {
        printf("curve[%d]: ", i);
        for(j = 0; j <= size; j++) {
            if(curve[i][j] != SIZE) {
                printf("%d ", curve[i][j]);
            }
        }
        printf("\n");
    }//*/
    /* fills neighbors in a 2D array of length (size x 2)
       with the indices of the two neighbors (index of SIZE means
       a neighbor wasn't found)
       -- initializes neighbors[] array */
    for(i = 0; i <= size; i++) {
        /* each point can have 0, 1, or 2 neighbors */
        neighbors[i] = malloc(sizeof(int) * 2);
        /* initialize to a non-existant point */
        neighbors[i][0] = SIZE;
        neighbors[i][1] = SIZE;
    }
    for(i = 0; i < division_number; i++) {
        /* calculate distances between all paths within 
           curve[i] */
        for(j = 0; j <= size; j++) {
            /* initializes shortest_path[] array to a non-existant
               node */
            shortest_path[0] = SIZE;
            shortest_path[1] = SIZE;
            /* first find a starting point (index j) */
            if(curve[i][j] != SIZE) {
                /* then compare distances between all the other
                   points (index k) */
                for(k = 0; k <= size; k++) {
                    if((curve[i][k] != SIZE) && (k != j)) {
                        /* stores the two shortest paths
                           -- checks to see if it is the first path
                           -- or then the second path
                           -- stores paths if the current distance is
                              less than the distance to the previous
                              shortest path */
                        if(shortest_path[0] == 14) {
                            shortest_path[0] = k;
                        }
                        else if(shortest_path[1] == 14) {
                            shortest_path[1] = k;
                            /* keeps the smallest value in index 0 */
                            if (shortest_path[1] < shortest_path[0]) {
                                tmp = shortest_path[0];
                                shortest_path[0] = shortest_path[1];
                                shortest_path[1] = tmp;
                            }
                        }
                        else if(distance_p(point[j], point[k]) <= distance_p(point[j], point[shortest_path[1]])) {
                            /* stores new shortest path */
                            shortest_path[1] = k;
                            /* keeps the smallest value in index 0 */
                            if (shortest_path[1] < shortest_path[0]) {
                                tmp = shortest_path[0];
                                shortest_path[0] = shortest_path[1];
                                shortest_path[1] = tmp;
                            }
                        }
                        else if(distance_p(point[j], point[k]) <= distance_p(point[j], point[shortest_path[1]])) {
                            /* stores new shortest path */
                            shortest_path[1] = k;
                        }
                    }
                }
                neighbors[j][0] = shortest_path[0];
                neighbors[j][1] = shortest_path[1];
            }
        }
    }
    /* printing for debug */
    for(i = 0; i <= size; i++) {
        printf("node %d: %d, %d\n", i, neighbors[i][0], neighbors[i][1]);
    }
    return curve;
}

/* calculates theta given two vectors */
double calculate_theta(struct vector_t V1, struct vector_t V2)
{
    return acos(dot_product(V1, V2)/(V1.length * V2.length));
}

/* calculates tao_distance given theta, length and curvature */
double tao_distance(double theta, double length, double curvature)
{
    return abs((tan(theta) + 1) * length / curvature);
}

/* constructs triangles given the curve they are on*/
struct triangle_t *construct_triangles(struct point_t **curve, int size,int divisions)
{
    struct triangle_t *triangle = malloc(sizeof(struct triangle_t) * 100);
    return triangle;
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
