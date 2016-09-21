OBJECTS = 
ALLOBJECTS := $(patsubst %.cpp,%.o,$(wildcard *.cpp))
CFLAGS= -c -O2 -I$(shell root-config --incdir) -std=c++11
COMPILER= clang++
LINKOPTION = $(shell root-config --libs)

all: plotLitData

plotLitData: plotLitData.o $(OBJECTS)
	$(COMPILER) -o plotLitData plotLitData.o $(OBJECTS) $(LINKOPTION)

$(ALLOBJECTS): %.o : %.cpp
	$(COMPILER) $(CFLAGS) $< -o $@

clean:
	rm -f *.o
