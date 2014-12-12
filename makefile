CC=g++
CFLAGS=-w -g
MAIN=main.o
INCLUDE=pipe.o io.o objects.o
TARGET=HW5 

all: $(INCLUDE) $(MAIN) $(TARGET)

%.o: %.cc
	$(CC) -c $? -o $@ $(CFLAGS)

$(TARGET): $(MAIN) $(INCLUDE)
	$(CC) -o $@ $? $(FLAGS)

clean:
	rm -v *.ppm *.o HW5
