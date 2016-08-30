#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#define SIZE 7

struct point_t {
    double x;
    double y;
    char pt;
    double curr;
    char shortest[9];
} point;

void swap(struct point_t *x, struct point_t *y);
void permute(struct point_t *point, struct point_t *shortest, int index, int n);
double distance(struct point_t *point, int first, int last);
int factorial(int n);

//char global[9] = {'\0'};
int global_count = 0;

int main(void)
{
    struct point_t *point = malloc(sizeof(struct point_t) * SIZE);
    struct point_t *shortest = malloc(sizeof(struct point_t) * SIZE);
    int i = 0;
    int j = 0;
    for(; i <= (SIZE - 1); i++) {
        for(; j <= 9; j++)
        point[i].shortest[j] = '\0';
    }
    shortest = point;
    int start = 0;
    int n = SIZE - 1;
    /* A(0,2) */
    point[0].x = 0;
    point[0].y = 2;
    point[0].pt = 'A';
    point[0].curr = DBL_MAX;
    /* B(2,2) */
    point[1].x = 2;
    point[1].y = 2;
    point[1].pt = 'B';
    point[1].curr = DBL_MAX;
    /* C(2,4) */
    point[2].x = 2;
    point[2].y = 4;
    point[2].pt = 'C';
    point[2].curr = DBL_MAX;
    /* D(-1,2) */
    point[3].x = -1;
    point[3].y = 2;
    point[3].pt = 'D';
    point[3].curr = DBL_MAX;
    /* E(-2,-2) */
    point[4].x = -2;
    point[4].y = -2;
    point[4].pt = 'E';
    point[4].curr = DBL_MAX;
    /* F(0,-3) */
    point[5].x = 0;
    point[5].y = -3;
    point[5].pt = 'F';
    point[5].curr = DBL_MAX;
    /* O(0,0) */
    point[6].x = 0;
    point[6].y = 0;
    point[6].pt = 'O';
    point[6].curr = DBL_MAX;
    /* fills array with permutations */
    permute(point, shortest, start, n);
    printf("\n");
    printf("n! = %d\n", factorial(SIZE));
    printf("Total Permutations: %d\n\n", global_count);
    printf("Shortest Path: %s\n", point[0].shortest);
    printf("Distance: %lf\n", shortest[0].curr);
    printf("\n");
    return 0;
}

/* Function to swap values at two pointers */
void swap(struct point_t *x, struct point_t *y)
{
    struct point_t tmp;// = malloc(sizeof(struct point_t));
    tmp.x = x->x;
    tmp.y = x->y;
    tmp.pt = x->pt;
    x->x = y->x;
    x->y = y->y;
    x->pt = y->pt;
    y->x = tmp.x;
    y->y = tmp.y;
    y->pt = tmp.pt;
}

/* calculates all permutations */
void permute(struct point_t *point, struct point_t *shortest, int index, int n)
{
    int i = 0;
    double total = 0;
    double curr = point[0].curr;
    double segment = 0;
    /* base case */
    if(index == n) {
        global_count++;
        /* calculating distance of segments from start to end */
        printf("%c", point[n].pt);
        for(i = 0; i < n; i++) {
            segment = distance(point, i, i + 1);
            printf("%c", point[i].pt);
            total += segment;
        }
        /* calculating the final segment (start and end nodes) */
        segment = distance(point, n, 0);
        total += segment;
        printf("%c: %lf\n", point[n].pt, total);
        /* checking to see if the most recent total is less than
           the current shortest length */
        if(total < curr) {
            point[0].curr = total;
            /* storing new shortest path
               -- workaround for using the global variable
               global[0] = shortest[n].pt;
               global[i + 1] = shortest[i].pt; */
            point[0].shortest[0] = point[n].pt;
            for(i = 0; i <= n; i++) {
                point[0].shortest[i + 1] = point[i].pt;
            }
            shortest = point;
            //printf("Calculated Path: %s\n", global);
            //printf("Distance: %lf\n", shortest[0].curr);
            curr = total;
        }
    }
    else {
        for(i = index; i <= n; i++) {
            swap(point+index, point+i);
            permute(point, shortest, index+1, n);
            swap(point+index, point+i); //backtrack
        }
    }
}

/* calculates distance given index and structure */
double distance(struct point_t *point, int first, int last)
{
    double x_1 = point[first].x;
    double y_1 = point[first].y;
    double x_2 = point[last].x;
    double y_2 = point[last].y;
    double dist = sqrt((x_2 - x_1)*(x_2 - x_1)+(y_2 - y_1)*(y_2 - y_1));
    printf(" %lf ", dist);
    return dist;
}

/* calculates factorial of an integer */
int factorial(int n)
{
    int result = 1;
    int i = 2;
    for(; i <= n; i++) {
        result *= i;
    }
    return result;
}
