CC= gcc
OBJECTS= pratt.o sequence.o block.o scan.o pattern.o hit.o search.o menu.o swiss.o al.o tree.o mst.o pam_dist.o help.o
CFLAGS= -O2
 
all: pratt 
 
pratt: $(OBJECTS)
	$(CC) -O4  -o pratt $(OBJECTS) -lm 

clean:
	touch *.c; make pratt

