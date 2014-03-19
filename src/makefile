CC=g++ 
CCOPT=g++ -O3
CFLAGS=-W -Wall -pedantic 
LDFLAGS= -lX11 -lpthread
EXEC=inpaint

all: inpainting.o 
	$(CCOPT) -o $(EXEC) inpainting.o $(CFLAGS) $(LDFLAGS)

inpainting.o: inpainting.cpp inpainting.h progressbar.h
	$(CCOPT) -o inpainting.o -c inpainting.cpp $(CFLAGS)
	
clean:
	rm -f *.c~ *.h~ *.cpp~ *.ini~ *.o makefile~ $(EXEC)
