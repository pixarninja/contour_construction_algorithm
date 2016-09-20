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
    struct point_t point1;
    struct point_t point2;
    double m;
};

struct triangle_t {
    char name[3];
    char points[SIZE - 1];
    struct function_t f1;
    struct function_t f2;
    struct function_t f3;
    double x;
    double y;
};

int global_count = 0;

double shortest_path(struct point_t *point, int size);
int **construct_curve(struct point_t *point, int size);
struct triangle_t *construct_triangles(struct point_t **curve, int size,int divisions);
double calculate_theta(struct vector_t V1, struct vector_t V2);
double calculate_curvature(struct point_t prev, struct point_t curr, struct point_t next);
double tao_distance(double theta, double length, double curvature);
double distance_p(struct point_t start, struct point_t end);
double distance_v(struct vector_t start, struct vector_t end);
double length_v(struct vector_t v);
struct vector_t subtract_v(struct vector_t V2, struct vector_t V1);
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
    struct point_t prev;
    struct vector_t *position = malloc(sizeof(struct vector_t) * size);
    double *distance = malloc(sizeof(double) * size);
    double range[2] = {0};
    double shortest = DBL_MAX;
    double longest = 0;
    double tao = 0;
    double division_constant = 0;
    double sum_x = 0;
    double sum_y = 0;
    int **curve = malloc(sizeof(int *) * size);
    int **neighbors = malloc(sizeof(int *) * size);
    int shortest_curvature[2];
    int i = 0;
    int j = 0;
    int k = 0;
    int tmp = 0;
    int division_number = 0;
    /* calculate average point */
    for(i = 0; i <= size; i++) {
        sum_x += point[i].x;
        sum_y += point[i].y;
    }
    center.x = sum_x / size;
    center.y = sum_y / size;
    /* calculate position vectors */
    for(i = 0; i <= size; i++) {
        position[i].point[0].x = center.x;
        position[i].point[0].y = center.y;
        position[i].point[0].index = SIZE;
        position[i].point[1].x = point[i].x;
        position[i].point[1].y = point[i].y;
        position[i].point[1].index = SIZE;
        position[i].i = point[i].x - center.x;
        position[i].j = point[i].y - center.y;
        position[i].length = length_v(position[i]);
    }
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
            if((distance[j] >= (range[0] - 0.01)) && (distance[j] <= (range[1] + 0.01))) {
                curve[i][j] = j;
            }
        }
        range[0] = range[1];
        k = 0;
    }
    /* printing for debug
    for(i = 0; i < division_number; i++) {
        printf("curve[%d]: ", i);
        for(j = 0; j <= size; j++) {
            if(curve[i][j] != SIZE) {
                printf("%d ", curve[i][j]);
            }
        }
        printf("\n");
    }*/
    /* fills neighbors in a 2D array of length (size x 2)
       with the indices of the two neighbors (index of SIZE means
       a neighbor wasn't found)
       -- initializes neighbors[] array to a non-existant node */
    for(i = 0; i <= size; i++) {
        /* each point can have 0, 1, or 2 neighbors */
        neighbors[i] = malloc(sizeof(int) * 2);
        neighbors[i][0] = SIZE;
        neighbors[i][1] = SIZE;
    }
    for(i = 0; i < division_number; i++) {
        /* calculate distances between all paths within curve[i] */
        for(j = 0; j <= size; j++) {
            /* initializes shortest_curvature[] array to a non-existant
               node */
            shortest_curvature[0] = SIZE;
            shortest_curvature[1] = SIZE;
            /* first find a starting point (index j) */
            if(curve[i][j] != SIZE) {
                /* then compare distances between all the other
                   points (index k)
                */
                /* initializing previous point for j */
                prev.x = point[j].x + (position[j].i / position[j].length);
                prev.y = point[j].y + (position[j].j / position[j].length);
                prev.index = SIZE;
                for(k = 0; k <= size; k++) {
                    if((curve[i][k] != SIZE) && (k != j)) {
                        /* stores the two shortest paths
                           -- stores paths if the current curvature is
                              less than the curvature of the last
                              point */
                        if(shortest_curvature[0] == SIZE) {
                            shortest_curvature[0] = k;
                        }
                        else if(shortest_curvature[1] == SIZE) {
                            shortest_curvature[1] = k;
                            /* keeps the smallest value in index 0 */
                            if (calculate_curvature(prev, point[j], point[shortest_curvature[1]]) <= calculate_curvature(prev, point[j], point[shortest_curvature[0]])) {
                                tmp = shortest_curvature[0];
                                shortest_curvature[0] = shortest_curvature[1];
                                shortest_curvature[1] = tmp;
                            }
                        }
                        else if(calculate_curvature(prev, point[j], point[k]) <= calculate_curvature(prev, point[j], point[shortest_curvature[1]])) {
                            /* stores new shortest path */
                            shortest_curvature[1] = k;
                            /* keeps the smallest value in index 0 */
                            if (calculate_curvature(prev, point[j], point[shortest_curvature[1]]) <= calculate_curvature(prev, point[j], point[shortest_curvature[0]])) {
                                tmp = shortest_curvature[0];
                                shortest_curvature[0] = shortest_curvature[1];
                                shortest_curvature[1] = tmp;
                            }
                        }
                    }
                }
                /* printing for debug
                printf("curv[%d] = %0.2lf, : curv[%d] = %0.2lf\n", shortest_curvature[0], calculate_curvature(prev, point[j], point[shortest_curvature[0]]), shortest_curvature[1], calculate_curvature(prev, point[j], point[shortest_curvature[1]])); */
                /* sets neighbors of node i */
                neighbors[j][0] = shortest_curvature[0];
                neighbors[j][1] = shortest_curvature[1];
            }
        }
    }
    /* printing for debug */
    printf("NEIGHBORS ARRAY:\n");
    for(i = 0; i <= size; i++) {
        printf("node %d: %d, %d\n", i, neighbors[i][0], neighbors[i][1]);
    }
    return curve;
}

/* constructs triangles given the curve they are on*/
struct triangle_t *construct_triangles(struct point_t **curve, int size,int divisions)
{
    struct triangle_t *triangle = malloc(sizeof(struct triangle_t) * 100);
    return triangle;
}

/* calculates curvature given the previous, current, and next points */
double calculate_curvature(struct point_t prev, struct point_t curr, struct point_t next)
{
    /* initializing point data
       -- initializing vector V <prev, curr>
    struct vector_t V1;
    V1.point[0].x = prev.x;
    V1.point[0].y = prev.y;
    V1.point[0].index = prev.index;
    V1.point[1].x = curr.x;
    V1.point[1].y = curr.y;
    V1.point[1].index = curr.index;
    V1.i = V1.point[1].x - V1.point[0].x;
    V1.j = V1.point[1].y - V1.point[0].y;
    V1.length = length_v(V1);*/
    /* -- initializing vector T1 */
    struct vector_t T1;
    T1.point[0].x = curr.x;
    T1.point[0].y = curr.y;
    T1.point[0].index = curr.index;
    T1.point[1].x = prev.x;
    T1.point[1].y = prev.y;
    T1.point[1].index = prev.index;
    T1.i = T1.point[1].x - T1.point[0].x;
    T1.j = T1.point[1].y - T1.point[0].y;
    T1.length = length_v(T1);
    /* initializing vector V2, <curr, next> */
    struct vector_t V2;
    V2.point[0].x = curr.x;
    V2.point[0].y = curr.y;
    V2.point[0].index = curr.index;
    V2.point[1].x = next.x;
    V2.point[1].y = next.y;
    V2.point[1].index = next.index;
    V2.i = V2.point[1].x - V2.point[0].x;
    V2.j = V2.point[1].y - V2.point[0].y;
    V2.length = length_v(V2);
    /* -- initializing vector T2 */
    struct vector_t T2;
    T2.point[0].x = curr.x;
    T2.point[0].y = curr.y;
    T2.point[0].index = curr.index;
    T2.point[1].x = (V2.i / V2.length) + T2.point[0].x;
    T2.point[1].y = (V2.j / V2.length) + T2.point[0].y;
    T2.point[1].index = SIZE;
    T2.i = T2.point[1].x - T2.point[0].x;
    T2.j = T2.point[1].y - T2.point[0].y;
    T2.length = length_v(T2);
    //return (length_v(subtract_v(T1, T2)) / calculate_theta(T1, T2));
    return (calculate_theta(T1, T2));
}

/* calculates theta given two vectors */
double calculate_theta(struct vector_t V1, struct vector_t V2)
{
    //printf("dot_product<%0.2lf,%0.2lf>,<%0.2lf,%0.2lf> = %0.2lf\n", V1.i, V1.j, V2.i, V2.j, dot_product(V1, V2));
    //printf("length<%0.2lf,%0.2lf> = %0.2lf,length<%0.2lf,%0.2lf> = %0.2lf\n", V1.i, V1.j, V1.length, V2.i, V2.j, V2.length);
    return acos(dot_product(V1, V2)/(V1.length * V2.length));
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

/* calculates tao_distance given theta, length and curvature */
double tao_distance(double theta, double length, double curvature)
{
    return abs((tan(theta) + 1) * length / curvature);
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

/* calculates length of a single vector */
double length_v(struct vector_t V)
{
    double i = V.i;
    double j = V.j;
    return sqrt((i * i) + (j * j));
}

/* calculates difference of a two vectors */
struct vector_t subtract_v(struct vector_t V2, struct vector_t V1)
{
    struct vector_t S;
    S.point[0].x = V1.point[1].x;
    S.point[0].y = V1.point[1].y;
    S.point[1].x = V2.point[1].x;
    S.point[1].y = V2.point[1].y;
    S.i = S.point[1].x - S.point[0].x;
    S.j = S.point[1].y - S.point[0].y;
    S.length = length_v(S);
    return S;
}
