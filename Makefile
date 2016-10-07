CC=gcc
CFLAGS=-g -Wall

all: contours

contours:
	$(CC) $(CFLAGS) construct_contours.c -o construct_contours -lm

clean: 
	rm construct_contours
