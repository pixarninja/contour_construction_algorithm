CC=gcc
CFLAGS=-g -Wall

all: contours tao2 shortest

shortest:
	$(CC) $(CFLAGS) shortest.c -o shortest -lm

contours:
	$(CC) $(CFLAGS) construct_contours.c -o construct_contours -lm

tao2:
	$(CC) $(CFLAGS) tao2.c -o tao2 -lm

clean: 
	rm shortest construct_contours tao2
