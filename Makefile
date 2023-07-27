all: test

CC = g++
OPT= -g -flto -Ofast --std=c++2a
CFLAGS = $(OPT) -Wall
LIBS = -lssl -lcrypto 

test: test.cc sketch.cc zipf.c hashutil.c
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

clean:
	rm -f test test.o
