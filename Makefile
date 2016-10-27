CC=gcc
CPP=g++
CFLAGS=-g -Wall

all: contours_cpp

contours_c:
	$(CC) $(CFLAGS) construct_contours.c -o construct_contours -lm

contours_cpp:
	$(CPP) $(CFLAGS) construct_contours.cpp -o construct_contours -lm

clean: 
	rm construct_contours
