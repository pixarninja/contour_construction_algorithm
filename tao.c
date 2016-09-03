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

int global_count = 0;

double shortest_path(struct point_t start, int n, struct vector_t T1, struct point_t *search, int size);
//double calculate_tao_distance(double tao, double vector_length);
double calculate_tao_distance(struct curvature_t k);
double distance_p(struct point_t start, struct point_t end);
double distance_v(struct vector_t start, struct vector_t end);
double length_v(struct vector_t v);
double dot_product(struct vector_t start, struct vector_t end);
void print_k(struct curvature_t k);

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
    total = shortest_path(point[start], start, T0, point, size);
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
    int count = 0;
    int total_size = size;
    char shortest[9] = {'\0'};
    double total = 0;
    struct point_t best;
    best.x = DBL_MAX;
    best.y = DBL_MAX;
    best.tao_d = DBL_MAX;
    best.pt = '\0';
    struct point_t *curr = malloc(sizeof(struct point_t) * size);
    /* remove start pt from search */
    for(i = 0; i <= size; i++) {
        if(search[i].pt == start.pt) {
            for(j = i; j < size; j++) {
                search[j] = search[j + 1];
            }
            break;
        }
    }
    size -= 1;
    for(i = 0; i <= size; i++) {
        curr[i].x = DBL_MAX;
        curr[i].y = DBL_MAX;
        curr[i].tao_d = DBL_MAX;
        curr[i].pt = '\0';
    }
    struct curvature_t k;
    /* initializing point 1 data for structure k
    -- initializing vector V */
    k.V.point[0].x = start.x;
    k.V.point[0].y = start.y;
    k.V.point[0].pt = start.pt;
    k.V.point[1].x = start.x;
    k.V.point[1].y = start.y;
    k.V.point[1].pt = start.pt;
    k.V.i = 0;
    k.V.j = 0;
    k.V.length = 0;
    /* -- initializing vector T1 */
    k.T1 = T1;
    /* normalizes n, if n is not already an index of start */
    n %= total_size;
    /* outer loop, calculates total distance */
    while(global_count <= total_size) {
        //i = n;
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
            k.T2.point[0].x = k.V.point[0].x/k.V.length;
            k.T2.point[0].y = k.V.point[0].y/k.V.length;
            k.T2.point[0].pt = 'T';
            k.T2.point[1].x = k.V.point[1].x/k.V.length;
            k.T2.point[1].y = k.V.point[1].y/k.V.length;
            k.T2.point[1].pt = 'T';
            k.T2.i = k.V.i/k.V.length;
            k.T2.j = k.V.j/k.V.length;
            k.T2.length = length_v(k.T2);
            /* -- initializing tao */
            //k.tao = ((k.T2.i - k.T1.i) * (k.T2.i - k.T1.i)) + ((k.T2.j - k.T1.j) * (k.T2.j - k.T1.j));
            k.tao = (dot_product(k.T1, k.T2) / (length_v(k.T1) * length_v(k.T2)));
            /* debug for tao less than 1 */
            //k.tao += 1;
            printf("%c --> %c; tao: %lf, ", k.V.point[0].pt, k.V.point[1].pt, k.tao);
            /* initializing current point */
            curr[i] = k.V.point[1];
            /* -- calculating tao-distance */
            curr[i].tao_d = calculate_tao_distance(k);
            printf("tao_d: %lf\n", curr[i].tao_d);
            /* increment and modulous */
            count++;
            i++;
            i %= total_size;
        }
        /* finds the path with the lowest tao-distance */
        for(i = 0; i <= size; i++) {
            if(best.tao_d > curr[i].tao_d) {
                best = curr[i];
                n = i;
            }
        }
        printf("Selected path: %c --> %c\n\n", k.V.point[0].pt, best.pt);
        /* remove selected point from search */
        for(i = n; i < size; i++) {
            search[i] = search[i + 1];
        }
        size -= 1;
        /* initializing structure k
           -- setting chosen endpoint */
        k.V.point[1].x = best.x;
        k.V.point[1].y = best.y;
        k.V.point[1].pt = best.pt;
        k.V.i = k.V.point[1].x - k.V.point[0].x;
        k.V.j = k.V.point[1].y - k.V.point[0].y;
        k.V.length = length_v(k.V);
        /* -- shifting T1 vector */
        k.T1.point[0].x = k.V.point[1].x;
        k.T1.point[0].y = k.V.point[1].y;
        /* -- calculating second point of T1 vector */
        k.T1.point[1].x = k.V.point[1].x + k.T2.i;
        k.T1.point[1].y = k.V.point[1].y + k.T2.j;
        k.T1.length = 1;
        /* -- shifting starting point */
        k.V.point[0].x = k.V.point[1].x;
        k.V.point[0].y = k.V.point[1].y;
        k.V.point[0].pt = k.V.point[1].pt;
        print_k(k); //for debugging purposes
        printf("New search: ");
        for(i = 0; i <= size; i++) {
            printf("%c ", search[i].pt);
        }
        printf("\n");
        /* initializing new values
           -- setting starting pt */
        start = best;
        /* -- setting counting variables */
        global_count++;
        count = 0;
    }
    return total;
}

/* calculates distance given index and structure */
double calculate_tao_distance(struct curvature_t k)
{
    double theta = 2 * acos(k.tao) + (M_PI / 180);
    double curvature = sqrt(((k.T2.i - k.T1.i) * (k.T2.i - k.T1.i)) + ((k.T2.j - k.T1.j) * (k.T2.j - k.T1.j))) / (sqrt((k.T1.i * k.T1.i) + (k.T1.j * k.T1.j)) * sqrt((k.T2.i * k.T2.i) + (k.T2.j * k.T2.j))) + 0.000001;
    printf("curvature: %lf; ", curvature);
    printf("angle = %lf; ", theta * 90 / M_PI);
    //return ((M_PI * vector_length)/(90 * sqrt(tao - 0.25))); //first implementation
    //return ((M_PI * (atan(sqrt(tao)/sqrt(tao - 0.25))) * vector_length)/(90 * sqrt(tao - 0.25)) + 1); //second implementation
    return abs(tan(theta) * k.V.length / curvature); //third implementation
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
        printf("%c:(%lf, %lf), %c(%lf, %lf), V:<%lf, %lf>, |V| = %lf\n", k.V.point[0].pt, k.V.point[0].x, k.V.point[0].y, k.V.point[1].pt, k.V.point[1].x, k.V.point[1].y, k.V.i, k.V.j, k.V.length);
        printf("%c1[0]:(%lf, %lf), %c1[1](%lf, %lf), T1:<%lf, %lf>, |T1| = %lf\n", k.T1.point[0].pt, k.T1.point[0].x, k.T1.point[0].y, k.T1.point[1].pt, k.T1.point[1].x, k.T1.point[1].y, k.T1.i, k.T1.j, k.T1.length);
        printf("%c2[0]:(%lf, %lf), %c2[1](%lf, %lf), T2:<%lf, %lf>, |T2| = %lf\n", k.T2.point[0].pt, k.T2.point[0].x, k.T2.point[0].y, k.T2.point[1].pt, k.T2.point[1].x, k.T2.point[1].y, k.T2.i, k.T2.j, k.T2.length);
	printf("tao = %lf\n\n", k.tao);
}

/*
        if(i == size) {
            / initializing structure k
            k.V.point[0].x = start.x;
            k.V.point[0].y = start.y;
            k.V.point[0].pt = start.pt;
            k.V.point[1].x = search[i].x;
            k.V.point[1].y = search[i].y;
            k.V.point[1].pt = search[i].pt;
            k.V.i = k.V.point[1].x - k.V.point[0].x;
            k.V.j = k.V.point[1].y - k.V.point[0].y;
            k.V.length = distance_p(start, search[i]);
            k.T1 = T1;
            k.T2.point[0].x = k.V.point[0].x/k.V.length;
            k.T2.point[0].y = k.V.point[0].y/k.V.length;
            k.T2.point[0].pt = 'T';
            k.T2.point[1].x = k.V.point[1].x/k.V.length;
            k.T2.point[1].y = k.V.point[1].y/k.V.length;
            k.T2.point[1].pt = 'T';
            k.T2.i = k.V.i/k.V.length;
            k.T2.j = k.V.j/k.V.length;
            k.T2.length = 1;
            k.tao = ((k.T2.i - k.T1.i) * (k.T2.i - k.T1.i)) + ((k.T2.j - k.T1.j) * (k.T2.j - k.T1.j));
            printf("%c --> %c: %lf, ", k.V.point[0].pt, k.V.point[1].pt, k.tao);
            //total += calculate_tao_distance(k.tao, k.V.length);
            global_count++;
            printf("%lf\n", total);
            return total;
        }
        else {
            //i = n;
            for(; i < size; i++) {
*/
