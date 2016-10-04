CC=gcc
CFLAGS=-g -Wall

all: circle cardioid ellipse

circle:
	$(CC) $(CFLAGS) circle.c -o circle -lm

ellipse:
	$(CC) $(CFLAGS) ellipse.c -o ellipse -lm

cardioid:
	$(CC) $(CFLAGS) cardioid.c -o cardioid -lm

clean: 
	rm circle cardioid ellipse
