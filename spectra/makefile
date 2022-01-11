# Makefile for interface example
  
CC=mpiicpc
CFLAGS= -mkl -O3

organics.out : ed_ssh.cpp
	${CC} $(CFLAGS) $< -o $@

clean:
	rm -f *.o *.x
