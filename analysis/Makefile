SHELL = /bin/sh
COMPILER= clang++ 
CFLAGS= -lstdc++ -O2 -I$(shell root-config --incdir) -std=c++11 -ggdb
SOURCES = src/branches.cpp src/crossSection.cpp src/CSPrereqs.cpp src/dataPoint.cpp src/dataSet.cpp src/driver.cpp src/experiment.cpp src/histos.cpp src/plots.cpp src/raw.cpp src/separate.cpp src/target.cpp src/veto.cpp src/waveform.cpp src/waveformFitting.cpp
HEADERS = $(shell echo include/*.h)
OBJECTS = $(SOURCES:.c=.o)

LINKOPTION = $(shell root-config --libs)

all: driver sumSingle sumAll text finetimeCheck readLitData subtractCS relativeCS relative

driver: $(OBJECTS)
	$(COMPILER) $(CFLAGS) -o driver $(OBJECTS) $(LINKOPTION)

sumSingle: src/sumSingle.cpp src/dataSet.cpp src/dataPoint.cpp src/CSPrereqs.cpp src/experiment.cpp src/target.cpp src/crossSection.cpp src/plots.cpp
	$(COMPILER) $(CFLAGS) -o sumSingle src/sumSingle.cpp src/dataSet.cpp src/dataPoint.cpp src/CSPrereqs.cpp src/experiment.cpp src/target.cpp src/crossSection.cpp src/plots.cpp $(LINKOPTION)

sumAll: src/sumAll.cpp src/dataSet.cpp src/dataPoint.cpp src/CSPrereqs.cpp src/experiment.cpp src/target.cpp src/crossSection.cpp src/plots.cpp
	$(COMPILER) $(CFLAGS) -o sumAll src/sumAll.cpp src/dataSet.cpp src/dataPoint.cpp src/CSPrereqs.cpp src/experiment.cpp src/target.cpp src/crossSection.cpp src/plots.cpp $(LINKOPTION)

readLitData: src/readLitData.cpp src/dataSet.cpp src/dataPoint.cpp
	$(COMPILER) $(CFLAGS) -o readLitData src/readLitData.cpp src/dataSet.cpp src/dataPoint.cpp $(LINKOPTION)

text: src/text.cpp src/raw.cpp src/branches.cpp
	$(COMPILER) $(CFLAGS) -o text src/text.cpp src/raw.cpp src/branches.cpp $(LINKOPTION)

finetimeCheck: src/finetimeCheck.cpp src/branches.cpp
	$(COMPILER) $(CFLAGS) -o finetimeCheck src/finetimeCheck.cpp src/branches.cpp $(LINKOPTION)

makeCSText: src/makeCSText.cpp src/dataPoint.cpp
	$(COMPILER) $(CFLAGS) -o makeCSText src/makeCSText.cpp src/dataPoint.cpp $(LINKOPTION)

subtractCS: src/subtractCS.cpp src/dataSet.cpp src/dataPoint.cpp src/crossSection.cpp src/target.cpp src/plots.cpp
	$(COMPILER) $(CFLAGS) -o subtractCS src/subtractCS.cpp src/dataSet.cpp src/dataPoint.cpp src/crossSection.cpp src/target.cpp src/plots.cpp $(LINKOPTION)

relativeCS: src/relativeCS.cpp src/dataSet.cpp src/dataPoint.cpp src/crossSection.cpp src/target.cpp src/plots.cpp
	$(COMPILER) $(CFLAGS) -o relativeCS src/relativeCS.cpp src/dataSet.cpp src/dataPoint.cpp src/crossSection.cpp src/target.cpp src/plots.cpp $(LINKOPTION)

relative: src/dataPoint.cpp src/dataSet.cpp src/experiment.cpp src/crossSection.cpp src/relative.cpp
	$(COMPILER) $(CFLAGS) -o relative src/dataPoint.cpp src/dataSet.cpp src/experiment.cpp src/crossSection.cpp src/relative.cpp $(LINKOPTION)

%.o : %.cpp $(HEADERS)
	$(COMPILER) $(CFLAGS) $< -o $@

#.PHONY: all

clean:
	rm -f *.o
