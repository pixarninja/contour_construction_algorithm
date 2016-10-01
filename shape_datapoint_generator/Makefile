CC=gcc
CFLAGS=-g -Wall

all: circle cardioid rose

circle:
	$(CC) $(CFLAGS) circle.c -o circle -lm

cardioid:
	$(CC) $(CFLAGS) cardioid.c -o cardioid -lm

rose:
	$(CC) $(CFLAGS) rose.c -o rose -lm

clean: 
	rm circle cardioid rose
