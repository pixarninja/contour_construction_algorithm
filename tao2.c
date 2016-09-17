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
}

struct triangle_t {
    char name[3];
    char *points;
    struct function_t f1;
    struct function_t f2;
    struct function_t f3;
    double x;
    double y;
}

int global_count = 0;

double shortest_path(struct point_t start, int n, struct point_t *search, int size);
double calculate_theta(struct curvature_t k);
double calculate_curvature(struct curvature_t k);
double tao_distance(struct curvature_t k);
double distance_p(struct point_t start, struct point_t end);
double distance_v(struct vector_t start, struct vector_t end);
double length_v(struct vector_t v);
double dot_product(struct vector_t start, struct vector_t end);
void print_k(struct curvature_t k);

int main(void)
{
    struct point_t *point = malloc(sizeof(struct point_t) * SIZE);
    struct point_t *shortest = malloc(sizeof(struct point_t) * SIZE);
    shortest = point;
    int start = 0;
    int size = SIZE - 1;
    double total = 0;
    /* O(0,0) */
    point[6].x = 0;
    point[6].y = 0;
    point[6].pt = 'O';
    point[6].tao_d = DBL_MAX;
    /* A(0,2) */
    point[3].x = 0;
    point[3].y = 2;
    point[3].pt = 'A';
    point[3].tao_d = DBL_MAX;
    /* B(2,2) */
    point[5].x = 2;
    point[5].y = 2;
    point[5].pt = 'B';
    point[5].tao_d = DBL_MAX;
    /* C(2,4) */
    point[4].x = 2;
    point[4].y = 4;
    point[4].pt = 'C';
    point[4].tao_d = DBL_MAX;
    /* D(-1,2) */
    point[2].x = -1;
    point[2].y = 2;
    point[2].pt = 'D';
    point[2].tao_d = DBL_MAX;
    /* E(-2,-2) */
    point[1].x = -2;
    point[1].y = -2;
    point[1].pt = 'E';
    point[1].tao_d = DBL_MAX;
    /* F(0,-3) */
    point[0].x = 0;
    point[0].y = -3;
    point[0].pt = 'F';
    point[0].tao_d = DBL_MAX;
    /* runs tao-distance algorithm on data */
    total = shortest_path(point[start], start, point, size);
    printf("\n");
    printf("Total Permutations: %d\n\n", global_count);
    printf("Shortest Path: %s\n", "derp");
    printf("Distance: %lf\n", total);
    printf("\n");
    return 0;
}

/* calculates the shortest path */
double shortest_path(struct point_t start, int n, struct point_t *search, int size)
{
    int i;
    int j;
    int count = 0;
    int total_size = size;
    double theta = 0;
    double curvature = 0;
    double total = 0;
    struct point_t best;
    best.x = DBL_MAX;
    best.y = DBL_MAX;
    best.tao_d = DBL_MAX;
    best.pt = '\0';
    /* remove start pt from search */
    for(i = 0; i <= size; i++) {
        if(search[i].pt == start.pt) {
            for(j = i; j < size; j++) {
                search[j] = search[j + 1];
            }
            break;
        }
    }
    size--;
    /* initializing structure curr */
    struct point_t *curr = malloc(sizeof(struct point_t) * size);
    for(i = 0; i <= size; i++) {
        curr[i].x = DBL_MAX;
        curr[i].y = DBL_MAX;
        curr[i].tao_d = DBL_MAX;
        curr[i].pt = '\0';
    }
    struct vector_t T1;
    T1.point[0].x = start.x;
    T1.point[0].y = start.y;
    T1.point[0].pt = start.pt;
    T1.point[1].x = start.x + 0.0;
    T1.point[1].y = start.y + 1.0;
    T1.point[1].pt = 'T';
    T1.i = 0.0;
    T1.j = 1.0;
    T1.length = 1;
    struct curvature_t k;
    /* initializing point 1 data for structure k
    -- initializing vector V */
    k.V.point[0].x = start.x;
    k.V.point[0].y = start.y;
    k.V.point[0].pt = start.pt;
    k.V.i = 0;
    k.V.j = 0;
    k.V.length = 0;
    /* -- initializing vector T1 */
    k.T1 = T1;
    /* -- initializing vector T2 */
    k.T2.point[0].x = start.x;
    k.T2.point[0].y = start.y;
    k.T2.point[0].pt = start.pt;
    /* normalizes n, if n is not already an index within the new size */
    n %= total_size;
    /* outer loop, calculates total distance */
    while(global_count <= total_size) {
        i = 0;
        /* refreshing best pt */
        best.tao_d = DBL_MAX;
        /* loops through all possible pts from start */
        while(count <= size) {
            /* initializing point 2 data for structure k
               -- initializing vector V */
            k.V.point[1].x = search[i].x;
            k.V.point[1].y = search[i].y;
            k.V.point[1].pt = search[i].pt;
            k.V.i = k.V.point[1].x - k.V.point[0].x;
            k.V.j = k.V.point[1].y - k.V.point[0].y;
            k.V.length = length_v(k.V);
            /* -- initializing vector T2 */
            k.T2.point[1].x = k.V.point[1].x;
            k.T2.point[1].y = k.V.point[1].y;
            k.T2.point[1].pt = 'T';
            k.T2.i = (k.T2.point[1].x - k.T2.point[0].x) / k.V.length;
            k.T2.j = (k.T2.point[1].y - k.T2.point[0].y) / k.V.length;
            k.T2.point[1].x = k.V.point[0].x + k.T2.i;
            k.T2.point[1].y = k.V.point[0].y + k.T2.j;
            k.T2.length = length_v(k.T2);
            /* -- initializing tao, theta, and curvature */
            k.tao = (dot_product(k.T1, k.T2)); //length of T1 and T2 is always 1
            if((k.tao >= 0.0) && (k.tao <= 0.0)) {
                k.tao = 0.1;
            }
            k.theta = calculate_theta(k);
            k.curvature = calculate_curvature(k);
            /* initializing structure curr */
            curr[i].x = k.V.point[1].x;
            curr[i].y = k.V.point[1].y;
            curr[i].pt = k.V.point[1].pt;
            /* calculating tao-distance */
            curr[i].tao_d = tao_distance(k);
            k.V.point[1].tao_d = curr[i].tao_d;
            print_k(k);
            i++;
            count++;
        }
        /* find point with the lowest tao-distance */
        for(i = 0; i <= size; i++) {
            if(best.tao_d > curr[i].tao_d) {
                best.x = curr[i].x;
                best.y = curr[i].y;
                best.pt = curr[i].pt;
                best.tao_d = curr[i].tao_d;
                n = i;
            }
        }
        printf("\nChosen path: %c --> %c \n\n", k.V.point[0].pt, best.pt);
        /* remove chosen point from search array */
        for(i = n; i < size; i++) {
            search[i].x = search[i + 1].x;
            search[i].y = search[i + 1].y;
            search[i].pt = search[i + 1].pt;
            search[i].tao_d = search[i + 1].tao_d;
        }
        size--;
        /* reinitializing structure k
           -- reinitializing vector V */
        k.V.point[1].x = best.x;
        k.V.point[1].y = best.y;
        k.V.point[1].pt = best.pt;
        k.V.i = k.V.point[1].x - k.V.point[0].x;
        k.V.j = k.V.point[1].y - k.V.point[0].y;
        k.V.length = length_v(k.V);
        /* -- reinitializing vector T1 */
        k.T2.point[1].x = best.x;
        k.T2.point[1].y = best.y;
        k.T2.point[1].pt = 'T';
        k.T2.i = (k.T2.point[1].x - k.T2.point[0].x) / k.V.length;
        k.T2.j = (k.T2.point[1].y - k.T2.point[0].y) / k.V.length;
        k.T2.length = length_v(k.T2);
        k.T1.point[0].x = best.x;
        k.T1.point[0].y = best.y;
        k.T1.point[0].pt = best.pt;
        k.T1.point[1].x = best.x + k.T2.i;
        k.T1.point[1].y = best.y + k.T2.j;
        k.T1.point[1].pt = 'T';
        k.T1.i = (k.T1.point[1].x - k.T1.point[0].x);
        k.T1.j = (k.T1.point[1].y - k.T1.point[0].y);
        k.T1.length = length_v(k.T1);
        /* -- shift starting point to best point */
        start = best;
        /* -- initializing vector T2 */
        k.T2.point[0].x = start.x;
        k.T2.point[0].y = start.y;
        k.T2.point[0].pt = start.pt;
        k.T2.i = 0;
        k.T2.j = 0;
        k.T2.length = 0;
        /* -- initializing vector V */
        k.V.point[0].x = start.x;
        k.V.point[0].y = start.y;
        k.V.point[0].pt = start.pt;
        k.V.i = 0;
        k.V.j = 0;
        k.V.length = 0;
        /* normalizes n, if n is not already an index within the new size */
        n %= total_size;
        count = 0;
        global_count++;
    }
    return total;
}

/* calculates theta given structure k */
double calculate_theta(struct curvature_t k)
{
    //return (2 * acos(k.tao) - (M_PI / 180));
    return (2 * acos(k.tao));
}

/* calculates curvature given structure k */
double calculate_curvature(struct curvature_t k)
{
    //return (sqrt(((k.T2.i - k.T1.i) * (k.T2.i - k.T1.i)) + ((k.T2.j - k.T1.j) * (k.T2.j - k.T1.j))) / (sqrt((k.T1.i * k.T1.i) + (k.T1.j * k.T1.j)) * sqrt((k.T2.i * k.T2.i) + (k.T2.j * k.T2.j))));
    //return (sqrt(((k.T2.i - k.T1.i) * (k.T2.i - k.T1.i)) + ((k.T2.j - k.T1.j) * (k.T2.j - k.T1.j))) / (sqrt((k.T1.i * k.T1.i) + (k.T1.j * k.T1.j)) * sqrt((k.T2.i * k.T2.i) + (k.T2.j * k.T2.j))) + 0.000001);
    return (sqrt(((k.T2.i - k.T1.i) * (k.T2.i - k.T1.i)) + ((k.T2.j - k.T1.j) * (k.T2.j - k.T1.j))) / calculate_theta(k));
}

/* calculates distance given index and structure */
double tao_distance(struct curvature_t k)
{
    return abs((tan(k.theta) + 1) * k.V.length / k.curvature);
    //return abs((tan(k.theta) * k.theta * k.tao) / (k.curvature + 29)) + k.V.length;
    //return abs((abs(tan(k.theta) * k.curvature) - k.curvature - k.V.length) * abs(k.tao));
    //return abs(k.curvature * k.V.length);
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

/* prints structure k for debugging */
void print_k(struct curvature_t k)
{
    printf("%c:(%lf, %lf), %c:(%lf, %lf), V:<%lf, %lf>, |V| = %lf\n", k.V.point[0].pt, k.V.point[0].x, k.V.point[0].y, k.V.point[1].pt, k.V.point[1].x, k.V.point[1].y, k.V.i, k.V.j, k.V.length);
    printf("%c1[0]:(%lf, %lf), %c1[1](%lf, %lf), T1:<%lf, %lf>, |T1| = %lf\n", k.T1.point[0].pt, k.T1.point[0].x, k.T1.point[0].y, k.T1.point[1].pt, k.T1.point[1].x, k.T1.point[1].y, k.T1.i, k.T1.j, k.T1.length);
    printf("%c2[0]:(%lf, %lf), %c2[1](%lf, %lf), T2:<%lf, %lf>, |T2| = %lf\n", k.T2.point[0].pt, k.T2.point[0].x, k.T2.point[0].y, k.T2.point[1].pt, k.T2.point[1].x, k.T2.point[1].y, k.T2.i, k.T2.j, k.T2.length);
    printf("curvature: %lf; ", k.curvature);
    printf("angle = %lf; ", k.theta * 90 / M_PI);
    printf("tao = %lf; ", k.tao);
    printf("tao-distance = %lf\n\n", k.V.point[1].tao_d);
}
