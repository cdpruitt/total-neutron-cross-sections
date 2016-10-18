SHELL = /bin/sh
COMPILER= clang++ 
CFLAGS= -lstdc++ -O0 -I$(shell root-config --incdir) -std=c++11 -ggdb -Wall -Wextra

SOURCES = src/branches.cpp src/crossSection.cpp src/dataPoint.cpp src/dataSet.cpp src/driver.cpp src/basicHistos.cpp src/plots.cpp src/raw.cpp src/separate.cpp src/assignMacro.cpp src/target.cpp src/waveform.cpp src/waveformFitting.cpp
COMMON = 
HEADERS = $(shell echo include/*.h)
OBJECTS = $(SOURCES:.c=.o)

ALLOBJECTS := $(patsubst %.cpp,%.o,$(wildcard src/*.cpp))
LINKOPTION = $(shell root-config --libs)

all: driver sumRun sumAll text finetimeCheck readLitData
	#install -m 755 driver bin
	#install -m 755 sumRun bin
	#install -m 755 sumAll bin
	#install -m 755 text bin
	#install -m 755 finetimeCheck bin
	#install -m 755 readLitData bin

driver: $(OBJECTS) $(COMMON)
	$(COMPILER) $(CFLAGS) -o driver $(OBJECTS) $(LINKOPTION)

sumRun: src/sumRun.cpp
	$(COMPILER) $(CFLAGS) -o sumRun src/sumRun.cpp $(LINKOPTION)

sumAll: src/sumAll.cpp src/dataPoint.cpp src/dataSet.cpp
	$(COMPILER) $(CFLAGS) -o sumAll src/sumAll.cpp src/dataPoint.cpp src/dataSet.cpp $(LINKOPTION)

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
