SHELL = /bin/sh
COMPILER= clang++ 
CFLAGS= -lstdc++ -O0 -I$(shell root-config --incdir) -std=c++11 -ggdb -Wall -Wextra

SOURCES = src/branches.cpp src/crossSection.cpp src/CSPrereqs.cpp src/dataPoint.cpp src/dataSet.cpp src/driver.cpp src/basicHistos.cpp src/plots.cpp src/raw.cpp src/altSeparate.cpp src/target.cpp src/vetoEvents.cpp src/waveform.cpp src/waveformFitting.cpp
COMMON = 
HEADERS = $(shell echo include/*.h)
OBJECTS = $(SOURCES:.c=.o)

ALLOBJECTS := $(patsubst %.cpp,%.o,$(wildcard src/*.cpp))
LINKOPTION = $(shell root-config --libs)

all: driver sumAll text finetimeCheck readLitData

driver: $(OBJECTS) $(COMMON)
	$(COMPILER) $(CFLAGS) -o driver $(OBJECTS) $(LINKOPTION)

sumAll: src/sumAll.cpp src/dataSet.cpp src/dataPoint.cpp src/CSPrereqs.cpp src/target.cpp src/crossSection.cpp src/plots.cpp
	$(COMPILER) $(CFLAGS) -o sumAll src/sumAll.cpp src/dataSet.cpp src/dataPoint.cpp src/CSPrereqs.cpp src/target.cpp src/crossSection.cpp src/plots.cpp $(LINKOPTION)

readLitData: src/readLitData.cpp src/dataSet.cpp src/dataPoint.cpp
	$(COMPILER) $(CFLAGS) -o readLitData src/readLitData.cpp src/dataSet.cpp src/dataPoint.cpp $(LINKOPTION)

text: src/text.cpp src/raw.cpp src/branches.cpp
	$(COMPILER) $(CFLAGS) -o text src/text.cpp src/raw.cpp src/branches.cpp $(LINKOPTION)

finetimeCheck: src/finetimeCheck.cpp src/branches.cpp
	$(COMPILER) $(CFLAGS) -o finetimeCheck src/finetimeCheck.cpp src/branches.cpp $(LINKOPTION)

makeCSText: src/makeCSText.cpp src/dataPoint.cpp
	$(COMPILER) $(CFLAGS) -o makeCSText src/makeCSText.cpp src/dataPoint.cpp $(LINKOPTION)

%.o : %.cpp $(HEADERS) $(COMMON)
	$(COMPILER) $(CFLAGS) $< -o $@

#.PHONY: all

clean:
	rm -f *.o
