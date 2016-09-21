SHELL = /bin/sh
COMPILER= clang++ 
CFLAGS= -lstdc++ -O0 -I$(shell root-config --incdir) -std=c++11 -ggdb -Wall -Wextra

SOURCES = $(shell echo src/*.cpp)
COMMON = 
HEADERS = $(shell echo include/*.h)
OBJECTS = $(SOURCES:.c=.o)

ALLOBJECTS := $(patsubst %.cpp,%.o,$(wildcard src/*.cpp))
LINKOPTION = $(shell root-config --libs)

all: driver sumRun sumAll

driver: $(OBJECTS) $(COMMON)
	$(COMPILER) $(CFLAGS) -o driver $(OBJECTS) $(LINKOPTION)

sumRun: sumRun.o
	$(COMPILER) $(CFLAGS) -o sumRun sumRun.o $(LINKOPTION)

sumAll: sumAll.o
	$(COMPILER) $(CFLAGS) -o sumAll sumAll.o $(LINKOPTION)

%.o : %.cpp $(HEADERS) $(COMMON)
	$(COMPILER) $(CFLAGS) $< -o $@

#text: text.o $(OBJECTS)
#	$(COMPILER) -o text text.o $(OBJECTS) $(LINKOPTION)

clean:
	rm -f *.o
