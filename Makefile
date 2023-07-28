all: test

CC = g++-10
OPT= -g -flto --std=c++2a
CFLAGS = $(OPT) -Wall
LIBS = -lssl -lcrypto 

test: test.cc zipf.c hashutil.c
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

clean:
	rm -f test test.o
