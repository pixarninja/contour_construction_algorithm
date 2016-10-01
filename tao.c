/**
 * Calculates the roundest path between the points using only curvature;
 * equates to the shortest path for a dataset of a "round" shape. I have
 * named this method of computation as "tao-distance", because it deals
 * with calculating a constant outlined below:
 *
 *    tao = ...
 *
 * Thus
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <ctype.h>
#include <unistd.h>

#define NUM_FILES 3

struct point_t {
    double x;
    double y;
    double tao_d;
    int index;
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
    double theta;
    double curvature;
};

int global_count = 0;

double shortest_path(struct point_t start, int n, struct point_t *search, int size, FILE *gnu_files[NUM_FILES]);
double calculate_curvature(struct curvature_t k);
double calculate_theta(struct curvature_t k);
double tao_distance(struct curvature_t k);
double angle_v(struct vector_t start, struct vector_t end);
double distance_p(struct point_t start, struct point_t end);
double distance_v(struct vector_t start, struct vector_t end);
double length_v(struct vector_t v);
double dot_product(struct vector_t start, struct vector_t end);
void print_k(struct curvature_t k);

int main(int argc, char *argv[])
{
    FILE *data;
    FILE *gnu_files[NUM_FILES];
    struct point_t *point;
    struct point_t *shortest;
    char buf[1024];
    double distance = 0.0;
    double range = 0.0;
    int size = 0;
    int start = 0;
    int i = 0;
    int c;
    int flag;
    if(argc == 1) {
        printf("Shape option [c s p] or help screen [h] not chosen.\nExiting program. Good day.\n");
        exit(EXIT_FAILURE);
    }
    while ((c = getopt(argc, argv, "csph")) != -1) {
        switch (c) {
        case 'c':
            data = fopen("./datapoints/tao_distance/circle.dat", "r");
            flag = 'c';
            break;
        case 's':
            data = fopen("./datapoints/tao_distance/square.dat", "r");
            flag = 's';
            break;
        case 'p':
            data = fopen("./datapoints/tao_distance/cardioid.dat", "r");
            flag = 'p';
            break;
        case 'h':
            printf("Enter -c for a circle, -s for a square, or -p for a cardioid\n");
            return 0;
        }
    }
    gnu_files[0] = fopen ("./gnu_files/commands.tmp", "w+");
    gnu_files[1] = fopen("./gnu_files/points.tmp", "w+");
    gnu_files[2] = fopen("./gnu_files/lines.tmp", "w+");
    while(fgets(buf, 1024, data)) {
        size++;
    }
    fclose(data);
    point = malloc(sizeof(struct point_t) * size);
    shortest = malloc(sizeof(struct point_t) * size);
    point = malloc(sizeof(struct point_t) * size);
    switch(flag) {
        case 'c':
            data = fopen("./datapoints/tao_distance/circle.dat", "r");
            break;
        case 's':
            data = fopen("./datapoints/tao_distance/square.dat", "r");
            break;
        case 'p':
            data = fopen("./datapoints/tao_distance/cardioid.dat", "r");
            break;
    }
    while(fscanf(data, "%d: (%lf, %lf)", &point[i].index, &point[i].x, &point[i].y) > 0) {
        if(abs(point[i].x) > range) {
            range = abs(point[i].x);
        }
        if(abs(point[i].y) > range) {
            range = abs(point[i].y);
        }
        i++;
    }
    fclose(data);
    /* stores data for gnu_points */
    for(i = 0; i < size; i++) {
        fprintf(gnu_files[1], "%lf %lf %d\n", point[i].x, point[i].y, point[i].index);
    }
    /* plot setup */
    fprintf(gnu_files[0], "set xrange [%lf:%lf]\n", -(range + 1), range + 1);
    fprintf(gnu_files[0], "set yrange [%lf:%lf]\n", -(range + 1), range + 1);
    fprintf(gnu_files[0], "set size ratio 1\n");
    fprintf(gnu_files[0], "set grid\n");
    fprintf(gnu_files[0], "set title \"Shortest Path of a 'Round' Dataset\"\n");
    fprintf(gnu_files[0], "set style line 1 lc rgb \"black\" lw 1\n");
    shortest = point;
    /* runs tao-distance algorithm on data */
    distance = shortest_path(point[start], start, point, size - 1, gnu_files);
    printf("\n");
    printf("Total Permutations: %d\n\n", global_count);
    printf("Shortest Path: %s\n", "TODO");
    printf("Distance: %lf\n", distance);
    printf("\n");
    /* plot */
    fclose(gnu_files[0]);
    fclose(gnu_files[1]);
    fclose(gnu_files[2]);
    system("gnuplot -persistent ./gnu_files/commands.tmp");
    return 0;
}

/* calculates the shortest path */
double shortest_path(struct point_t start, int n, struct point_t *search, int size, FILE *gnu_files[NUM_FILES])
{
    struct vector_t initial;
    struct vector_t check;
    struct point_t best;
    struct point_t end = start;
    struct point_t center;
    double total = 0.0;
    double sum_x = 0.0;
    double sum_y = 0.0;
    double theta = DBL_MAX;
    double tmp = 0.0;
    int i;
    int j;
    int index;
    int count = 0;
    int total_size = size - 1;
    best.x = DBL_MAX;
    best.y = DBL_MAX;
    best.tao_d = DBL_MAX;
    best.index = INT_MAX;
    /* calculate average point */
    for(i = 0; i <= size; i++) {
        sum_x += search[i].x;
        sum_y += search[i].y;
    }
    center.x = sum_x / size;
    center.y = sum_y / size;
    /* initialize initial vector */
    initial.point[0].x = start.x;
    initial.point[0].y = start.y;
    initial.point[0].index = start.index;
    initial.i = (start.x - center.x) / distance_p(start, center);
    initial.j = (start.y - center.y) / distance_p(start, center);
    initial.point[1].x = start.x + initial.i;
    initial.point[1].y = start.y + initial.j;
    initial.point[1].index = INT_MAX;
    initial.length = length_v(initial);
    /* remove start index from search */
    for(i = 0; i <= size; i++) {
        if(search[i].index == start.index) {
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
        curr[i].index = INT_MAX;
    }
    struct curvature_t k;
    /* initializing point 1 data for structure k
    -- initializing vector V */
    k.V.point[0].x = start.x;
    k.V.point[0].y = start.y;
    k.V.point[0].index = start.index;
    k.V.i = 0;
    k.V.j = 0;
    k.V.length = 0;
    /* -- initializing vector T1 */
    k.T1.point[0].x = start.x;
    k.T1.point[0].y = start.y;
    k.T1.point[0].index = start.index;
    /* checks for the point with the smallest deviation in angle from
       the position vector of the starting point (selects the "second"
       point) */
    for(i = 0; i <= size; i ++) {
        check.point[0].x = start.x;
        check.point[0].y = start.y;
        check.point[0].index = start.index;
        check.point[1].x = search[i].x;
        check.point[1].y = search[i].y;
        check.point[1].index = search[i].index;
        check.i = check.point[1].x - check.point[0].x;
        check.j = check.point[1].y - check.point[0].y;
        check.length = length_v(check);
        tmp = angle_v(initial, check);
        if(tmp < theta) {
            theta = tmp;
            index = i;
        }
    }
    /* points T1 in the direction of the "second point" */
    k.T1.i = (search[index].x - start.x) / distance_p(search[index], start);
    k.T1.j = (search[index].y - start.y) / distance_p(search[index], start);
    k.T1.point[1].x = start.x + k.T1.i;
    k.T1.point[1].y = start.y + k.T1.j;
    k.T1.point[1].index = INT_MAX;
    k.T1.length = 1;
    /* -- initializing vector T2 */
    k.T2.point[0].x = start.x;
    k.T2.point[0].y = start.y;
    k.T2.point[0].index = start.index;
    /* normalizes n, if n is not already an index within the new size */
    n %= total_size;
    /* plot */
    fprintf(gnu_files[2], "%lf %lf\n", start.x, start.y);
    /* outer loop, calculates total distance */
    while(global_count <= total_size) {
        i = 0;
        /* refreshing best index */
        best.tao_d = DBL_MAX;
        /* loops through all possible indices from start */
        while(count <= size) {
            /* initializing point 2 data for structure k
               -- initializing vector V */
            k.V.point[1].x = search[i].x;
            k.V.point[1].y = search[i].y;
            k.V.point[1].index = search[i].index;
            k.V.i = k.V.point[1].x - k.V.point[0].x;
            k.V.j = k.V.point[1].y - k.V.point[0].y;
            k.V.length = length_v(k.V);
            /* -- initializing vector T2 */
            k.T2.point[1].x = k.V.point[1].x;
            k.T2.point[1].y = k.V.point[1].y;
            k.T2.point[1].index = INT_MAX;
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
            curr[i].index = k.V.point[1].index;
            /* calculating tao-distance */
            curr[i].tao_d = tao_distance(k);
            k.V.point[1].tao_d = curr[i].tao_d;
            /* for debugging tao-distance function */
            //print_k(k);
            i++;
            count++;
        }
        /* find point with the lowest tao-distance */
        for(i = 0; i <= size; i++) {
            if(best.tao_d > curr[i].tao_d) {
                best.x = curr[i].x;
                best.y = curr[i].y;
                best.index = curr[i].index;
                best.tao_d = curr[i].tao_d;
                n = i;
            }
        }
        /* plot */
        fprintf(gnu_files[2], "%lf %lf\n", best.x, best.y);
        /* remove chosen point from search array */
        for(i = n; i < size; i++) {
            search[i].x = search[i + 1].x;
            search[i].y = search[i + 1].y;
            search[i].index = search[i + 1].index;
            search[i].tao_d = search[i + 1].tao_d;
        }
        size--;
        /* reinitializing structure k
           -- reinitializing vector V */
        k.V.point[1].x = best.x;
        k.V.point[1].y = best.y;
        k.V.point[1].index = best.index;
        k.V.i = k.V.point[1].x - k.V.point[0].x;
        k.V.j = k.V.point[1].y - k.V.point[0].y;
        k.V.length = length_v(k.V);
        /* -- reinitializing vector T1 */
        k.T2.point[1].x = best.x;
        k.T2.point[1].y = best.y;
        k.T2.point[1].index = 'T';
        k.T2.i = (k.T2.point[1].x - k.T2.point[0].x) / k.V.length;
        k.T2.j = (k.T2.point[1].y - k.T2.point[0].y) / k.V.length;
        k.T2.length = length_v(k.T2);
        k.T1.point[0].x = best.x;
        k.T1.point[0].y = best.y;
        k.T1.point[0].index = best.index;
        k.T1.point[1].x = best.x + k.T2.i;
        k.T1.point[1].y = best.y + k.T2.j;
        k.T1.point[1].index = 'T';
        k.T1.i = (k.T1.point[1].x - k.T1.point[0].x);
        k.T1.j = (k.T1.point[1].y - k.T1.point[0].y);
        k.T1.length = length_v(k.T1);
        /* -- shift starting point to best point */
        start = best;
        /* -- initializing vector T2 */
        k.T2.point[0].x = start.x;
        k.T2.point[0].y = start.y;
        k.T2.point[0].index = start.index;
        k.T2.i = 0;
        k.T2.j = 0;
        k.T2.length = 0;
        /* -- initializing vector V */
        k.V.point[0].x = start.x;
        k.V.point[0].y = start.y;
        k.V.point[0].index = start.index;
        k.V.i = 0;
        k.V.j = 0;
        k.V.length = 0;
        /* normalizes n, if n is not already an index within the new size */
        n %= total_size;
        count = 0;
        global_count++;
    }
    /* final point */
    fprintf(gnu_files[2], "%lf %lf\n", end.x, end.y);
    /* plot */
    fprintf(gnu_files[0], "plot './gnu_files/lines.tmp' with lines ls 1 title \"shortest path\",");
    fprintf(gnu_files[0], "'./gnu_files/points.tmp' using 1:2 with points pt 7 notitle,");
    fprintf(gnu_files[0], "'' using 1:2:3 with labels point pt 7 offset char -1,-1 notitle\n");
    return total;
}

/* calculates curvature given structure k */
double calculate_curvature(struct curvature_t k)
{
    return (sqrt(pow((k.T2.i - k.T1.i), 2) + pow((k.T2.j - k.T1.j), 2)) / calculate_theta(k));
}

/* calculates theta given structure k */
double calculate_theta(struct curvature_t k)
{
    return (acos(k.tao) + (M_PI / 180));
}

/* calculates distance given index and structure */
double tao_distance(struct curvature_t k)
{
    return (k.V.length * (k.curvature + 0.000001));
}

/* calculates angle between two vectors */
double angle_v(struct vector_t start, struct vector_t end)
{
    return (acos(dot_product(start, end) / (start.length * end.length)));
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
    printf("%d:(%lf, %lf), %d:(%lf, %lf), V:<%lf, %lf>, |V| = %lf\n", k.V.point[0].index, k.V.point[0].x, k.V.point[0].y, k.V.point[1].index, k.V.point[1].x, k.V.point[1].y, k.V.i, k.V.j, k.V.length);
    printf("%d[0]:(%lf, %lf), %d[1](%lf, %lf), T1:<%lf, %lf>, |T1| = %lf\n", k.T1.point[0].index, k.T1.point[0].x, k.T1.point[0].y, k.T1.point[1].index, k.T1.point[1].x, k.T1.point[1].y, k.T1.i, k.T1.j, k.T1.length);
    printf("%d[0]:(%lf, %lf), %d[1](%lf, %lf), T2:<%lf, %lf>, |T2| = %lf\n", k.T2.point[0].index, k.T2.point[0].x, k.T2.point[0].y, k.T2.point[1].index, k.T2.point[1].x, k.T2.point[1].y, k.T2.i, k.T2.j, k.T2.length);
    printf("curvature: %lf; ", k.curvature);
    printf("angle = %lf; ", k.theta * 90 / M_PI);
    printf("tao = %lf; ", k.tao);
    printf("tao-distance = %lf\n\n", k.V.point[1].tao_d);
}
