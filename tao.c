#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#define SIZE 7

struct point_t {
    double x;
    double y;
    double tao_d;
    char pt;
};

struct vector_t {
    struct point_t point[2];
    double length;
    double i;
    double j;
};

struct curvature_t {
    struct vector_t T1;
    struct vector_t T2;
    struct vector_t V;
    double tao;
};

double shortest_path(struct point_t start, int n, struct vector_t T1, struct point_t *search, int size);
double calculate_tao_distance(double tao, double vector_length);
double distance_p(struct point_t start, struct point_t end);
double distance_v(struct vector_t start, struct vector_t end);

char global_shortest[9] = {'\0'};
int global_count = 0;

int main(void)
{
    struct point_t *point = malloc(sizeof(struct point_t) * SIZE);
    struct point_t *shortest = malloc(sizeof(struct point_t) * SIZE);
    struct vector_t T0;
    T0.point[0].x = 0.0;
    T0.point[0].y = 0.0;
    T0.point[0].pt = 'T';
    T0.point[1].x = 0.0;
    T0.point[1].y = 1.0;
    T0.point[1].pt = 'T';
    T0.i = 0.0;
    T0.j = 1.0;
    T0.length = 1;
    shortest = point;
    int start = 0;
    int size = SIZE - 1;
    double total = 0;
    /* O(0,0) */
    point[0].x = 0;
    point[0].y = 0;
    point[0].pt = 'O';
    point[0].tao_d = DBL_MAX;
    /* A(0,2) */
    point[1].x = 0;
    point[1].y = 2;
    point[1].pt = 'A';
    point[1].tao_d = DBL_MAX;
    /* B(2,2) */
    point[2].x = 2;
    point[2].y = 2;
    point[2].pt = 'B';
    point[2].tao_d = DBL_MAX;
    /* C(2,4) */
    point[3].x = 2;
    point[3].y = 4;
    point[3].pt = 'C';
    point[3].tao_d = DBL_MAX;
    /* D(-1,2) */
    point[4].x = -1;
    point[4].y = 2;
    point[4].pt = 'D';
    point[4].tao_d = DBL_MAX;
    /* E(-2,-2) */
    point[5].x = -2;
    point[5].y = -2;
    point[5].pt = 'E';
    point[5].tao_d = DBL_MAX;
    /* F(0,-3) */
    point[6].x = 0;
    point[6].y = -3;
    point[6].pt = 'F';
    point[6].tao_d = DBL_MAX;
    /* runs tao-distance algorithm on data */
    total = shortest_path(point[start], start, T0, point + 1, size);
    printf("\n");
    printf("Total Permutations: %d\n\n", global_count);
    printf("Shortest Path: %s\n", "derp");
    printf("Distance: %lf\n", total);
    printf("\n");
    return 0;
}

/* calculates the shortest path */
double shortest_path(struct point_t start, int n, struct vector_t T1, struct point_t *search, int size)
{
    int i;
    int j;
    double total = 0;
    struct point_t best;
    best.x = DBL_MAX;
    best.y = DBL_MAX;
    best.tao_d = DBL_MAX;
    best.pt = '\0';
    struct point_t *curr = malloc(sizeof(struct point_t) * size);
    for(i = 0; i < size; i++) {
        curr[i].x = DBL_MAX;
        curr[i].y = DBL_MAX;
        curr[i].tao_d = DBL_MAX;
        curr[i].pt = '\0';
    }
    struct curvature_t k;
    /* base case */
    if(n == size) {
        i = n;
        /* initializing structure k
        -- initializing vector V */
        k.V.point[0].x = start.x;
        k.V.point[0].y = start.y;
        k.V.point[0].pt = start.pt;
        k.V.point[1].x = search[i].x;
        k.V.point[1].y = search[i].y;
        k.V.point[1].pt = search[i].pt;
        k.V.i = k.V.point[1].x - k.V.point[0].x;
        k.V.j = k.V.point[1].y - k.V.point[0].y;
        k.V.length = distance_p(start, search[i]);
        /* initializing vector T1 */
        k.T1 = T1;
        /* initializing vector T2 */
        k.T2.point[0].x = k.V.point[0].x/k.V.length;
        k.T2.point[0].y = k.V.point[0].y/k.V.length;
        k.T2.point[0].pt = 'T';
        k.T2.point[1].x = k.V.point[1].x/k.V.length;
        k.T2.point[1].y = k.V.point[1].y/k.V.length;
        k.T2.point[1].pt = 'T';
        k.T2.i = k.V.i/k.V.length;
        k.T2.j = k.V.j/k.V.length;
        k.T2.length = 1;
        /* initializing tao */
        k.tao = ((k.T2.i - k.T1.i) * (k.T2.i - k.T1.i)) + ((k.T2.j - k.T1.j) * (k.T2.j - k.T1.j));
        /* calculating tao-distance */
        printf("%c --> %c: %lf, ", k.V.point[0].pt, k.V.point[1].pt, k.tao);
        total += calculate_tao_distance(k.tao, k.V.length);
        global_count++;
        printf("%lf\n", total);
        return total;
    }
    else {
        i = n;
        for(; i < size; i++) {
            /* initializing structure k
               -- initializing vector V */
            k.V.point[0].x = start.x;
            k.V.point[0].y = start.y;
            k.V.point[0].pt = start.pt;
            k.V.point[1].x = search[i].x;
            k.V.point[1].y = search[i].y;
            k.V.point[1].pt = search[i].pt;
            k.V.i = k.V.point[1].x - k.V.point[0].x;
            k.V.j = k.V.point[1].y - k.V.point[0].y;
            k.V.length = distance_p(start, search[i]);
            /* -- initializing vector T1 */
            k.T1 = T1;
            /* -- initializing vector T2 */
            k.T2.point[0].x = k.V.point[0].x/k.V.length;
            k.T2.point[0].y = k.V.point[0].y/k.V.length;
            k.T2.point[0].pt = 'T';
            k.T2.point[1].x = k.V.point[1].x/k.V.length;
            k.T2.point[1].y = k.V.point[1].y/k.V.length;
            k.T2.point[1].pt = 'T';
            k.T2.i = k.V.i/k.V.length;
            k.T2.j = k.V.j/k.V.length;
            k.T2.length = 1;
            /* -- initializing tao */
            k.tao = ((k.T2.i - k.T1.i) * (k.T2.i - k.T1.i)) + ((k.T2.j - k.T1.j) * (k.T2.j - k.T1.j));
            printf("%c --> %c: %lf, ", k.V.point[0].pt, k.V.point[1].pt, k.tao);
            /* calculating tao-distance */
            curr[i].tao_d = calculate_tao_distance(k.tao, k.V.length);
        }
        for(i = n + 1; i < size; i++) {
            if(best.tao_d > curr[i].tao_d) {
                best = curr[i];
            }
        }
        total = best.tao_d;
        /* remove pt for best_tao_d */
        for(i = 0; i < size; i++) {
            if(search[i].pt == best.pt) {
                for(j = i; j < size; j++) {
                    search[j] = search[j + 1];
                }
                break;
            }
        }
        /* shifting T2 vector, to be T1 of next point */
        k.T2.point[0].x = k.V.point[1].x;
        k.T2.point[0].y = k.V.point[1].y;
        /* -- calculating second point of T2 vector */
        k.T2.point[1].x = k.V.point[1].x + k.T2.i;
        k.T2.point[1].y = k.V.point[1].y + k.T2.j;
        k.T2.length = 1;
        total += shortest_path(best, n + 1, k.T2, search + 1, size - 1);
        global_count++;
        printf("%lf\n", total);
        return total;
    }
}

/* calculates distance given index and structure */
double calculate_tao_distance(double tao, double vector_length)
{
    return ((M_PI * (atan(sqrt(tao)/sqrt(tao - 0.25))) * vector_length)/(90 * sqrt(tao - 0.25)));
    //return ((M_PI * vector_length)/(90 * sqrt(tao - 0.25)));
}

/* calculates distance given two points */
double distance_p(struct point_t start, struct point_t end)
{
    double x_1 = start.x;
    double y_1 = start.y;
    double x_2 = end.x;
    double y_2 = end.y;
    return sqrt((x_2 - x_1)*(x_2 - x_1)+(y_2 - y_1)*(y_2 - y_1));
}

/* calculates distance given two vectors */
double distance_v(struct vector_t start, struct vector_t end)
{
    double i_1 = start.i;
    double j_1 = start.j;
    double i_2 = end.i;
    double j_2 = end.j;
    return sqrt((i_2 - i_1)*(i_2 - i_1)+(j_2 - j_1)*(j_2 - j_1));
}
