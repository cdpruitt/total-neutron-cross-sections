SHELL = /bin/sh
COMPILER= clang++ 
CFLAGS= -lstdc++ -O0 -I$(shell root-config --incdir) -std=c++11 -ggdb -Wall -Wextra

SOURCES = src/branches.cpp src/crossSection.cpp src/dataPoint.cpp src/driver.cpp src/histos.cpp src/plots.cpp src/raw.cpp src/resort.cpp src/target.cpp src/waveform.cpp src/waveformFitting.cpp
COMMON = 
HEADERS = $(shell echo include/*.h)
OBJECTS = $(SOURCES:.c=.o)

ALLOBJECTS := $(patsubst %.cpp,%.o,$(wildcard src/*.cpp))
LINKOPTION = $(shell root-config --libs)

all: driver sumRun sumAll readLitData text

driver: $(OBJECTS) $(COMMON)
	$(COMPILER) $(CFLAGS) -o driver $(OBJECTS) $(LINKOPTION)

sumRun: src/sumRun.cpp
	$(COMPILER) $(CFLAGS) -o sumRun src/sumRun.cpp $(LINKOPTION)

sumAll: src/sumAll.cpp
	$(COMPILER) $(CFLAGS) -o sumAll src/sumAll.cpp $(LINKOPTION)

readLitData: src/readLitData.cpp src/dataSet.cpp
	$(COMPILER) $(CFLAGS) -o readLitData src/readLitData.cpp src/dataSet.cpp $(LINKOPTION)

text: src/text.cpp
	$(COMPILER) $(CFLAGS) -o text src/text.cpp $(LINKOPTION)

%.o : %.cpp $(HEADERS) $(COMMON)
	$(COMPILER) $(CFLAGS) $< -o $@

clean:
	rm -f *.o
