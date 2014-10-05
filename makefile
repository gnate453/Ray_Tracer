CC=gcc
CFLAGS=-w -Wall -Wextra -Wpedantic -g
MAIN=main.o
INCLUDE=pipe.o io.o
TARGET=HW3 

all: $(INCLUDE) $(MAIN) $(TARGET)

%.o: %.cc
	$(CC) -c $? -o $@ $(CFLAGS)

$(TARGET): $(MAIN) $(OBJS)
	$(CC) -o $@ $? $(FLAGS)

clean:
	rm -v *.o HW3
