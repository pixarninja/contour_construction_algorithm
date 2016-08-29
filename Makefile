CC=gcc
CFLAGS=-g -Wall

all: clean shortest

shortest:
	$(CC) $(CFLAGS) shortest.c -o shortest -lm

clean: 
	rm shortest
