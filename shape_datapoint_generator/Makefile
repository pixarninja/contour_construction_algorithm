CC=gcc
CFLAGS=-g -Wall

all: circle cardioid ellipse square triangle

circle:
	$(CC) $(CFLAGS) circle.c -o circle -lm

cardioid:
	$(CC) $(CFLAGS) cardioid.c -o cardioid -lm

ellipse:
	$(CC) $(CFLAGS) ellipse.c -o ellipse -lm

square:
	$(CC) $(CFLAGS) square.c -o square -lm

triangle:
	$(CC) $(CFLAGS) triangle.c -o triangle -lm

clean: 
	rm circle cardioid ellipse square triangle
