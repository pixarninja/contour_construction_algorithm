CC=gcc
CFLAGS=-g -Wall

all: clean shortest tao

shortest:
	$(CC) $(CFLAGS) shortest.c -o shortest -lm

tao:
	$(CC) $(CFLAGS) tao.c -o tao -lm

clean: 
	rm shortest tao
