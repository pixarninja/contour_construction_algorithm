CC=gcc
CFLAGS=-g -Wall

all: circles cardioid ellipse square triangle two_circles

circles:
	$(CC) $(CFLAGS) circle.c -o circle -lm
	$(CC) $(CFLAGS) two_circles.c -o two_circles -lm
	$(CC) $(CFLAGS) three_circles.c -o three_circles -lm
	$(CC) $(CFLAGS) donut.c -o donut -lm

cardioid:
	$(CC) $(CFLAGS) cardioid.c -o cardioid -lm

ellipse:
	$(CC) $(CFLAGS) ellipse.c -o ellipse -lm

square:
	$(CC) $(CFLAGS) square.c -o square -lm

triangle:
	$(CC) $(CFLAGS) triangle.c -o triangle -lm

clean: 
	rm circle cardioid donut ellipse square triangle two_circles
