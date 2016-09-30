CC=gcc
CFLAGS=-g -Wall

all: tao tao2

shortest:
	$(CC) $(CFLAGS) shortest.c -o shortest -lm

tao:
	$(CC) $(CFLAGS) tao.c -o tao -lm

tao2:
	$(CC) $(CFLAGS) tao2.c -o tao2 -lm

clean: 
	rm shortest tao tao2
