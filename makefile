CC=g++
CFLAGS=-w -Wall -Wextra -Wpedantic -g
MAIN=main.o
INCLUDE=pipe.o io.o objects.o
TARGET=HW4 

all: $(INCLUDE) $(MAIN) $(TARGET)

%.o: %.cc
	$(CC) -c $? -o $@ $(CFLAGS)

$(TARGET): $(MAIN) $(INCLUDE)
	$(CC) -o $@ $? $(FLAGS)

clean:
	rm -v *.ppm *.o HW4
