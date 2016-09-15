SHELL = /bin/sh
COMPILER= clang++ 
CFLAGS= -lstdc++ -O0 -I$(shell root-config --incdir) -std=c++11 -ggdb -Wall -Wextra

TARGET = analyze
SOURCES = $(shell echo src/*.cpp)
COMMON = 
HEADERS = $(shell echo include/*.h)
OBJECTS = $(SOURCES:.c=.o)

ALLOBJECTS := $(patsubst %.cpp,%.o,$(wildcard src/*.cpp))
LINKOPTION = $(shell root-config --libs)

all: $(TARGET)

$(TARGET): $(OBJECTS) $(COMMON)
	$(COMPILER) $(CFLAGS) -o $(TARGET) $(OBJECTS) $(LINKOPTION)

%.o : %.cpp $(HEADERS) $(COMMON)
	$(COMPILER) $(CFLAGS) $< -o $@

#driver DataPoint crossSection target DataSet raw text
#resort histos waveform DPPwaveform sumRun sumAll plotMonitor text
#
#
#
#target: target.o $(OBJECTS)
#	$(COMPILER) -o target target.o $(OBJECTS) $(LINKOPTION)
#
#driver: driver.o $(OBJECTS)
#	$(COMPILER) -o driver driver.o $(OBJECTS) $(LINKOPTION)
#
#raw: raw.o $(OBJECTS)
#	$(COMPILER) -o raw raw.o $(OBJECTS) $(LINKOPTION)
#
#crossSection: crossSection.o $(OBJECTS)
#	$(COMPILER) -o crossSection crossSection.o $(OBJECTS) $(LINKOPTION)
#
#DataPoint: DataPoint.o $(OBJECTS)
#	$(COMPILER) -o DataPoint DataPoint.o $(OBJECTS) $(LINKOPTION)
#
#DataSet: DataSet.o $(OBJECTS)
#	$(COMPILER) -o DataSet DataSet.o $(OBJECTS) $(LINKOPTION)
#
#resort: resort.o $(OBJECTS)
	#$(COMPILER) -o resort resort.o $(OBJECTS) $(LINKOPTION)

#histos: histos.o $(OBJECTS)
#	$(COMPILER) -o histos histos.o $(OBJECTS) $(LINKOPTION)

#finetime: finetime.o $(OBJECTS)
#	$(COMPILER) -o finetime finetime.o $(OBJECTS) $(LINKOPTION)

#plotMonitor: plotMonitor.o $(OBJECTS)
	#$(COMPILER) -o plotMonitor plotMonitor.o $(OBJECTS) $(LINKOPTION)

#sumRun: sumRun.o $(OBJECTS)
#	$(COMPILER) -o sumRun sumRun.o $(OBJECTS) $(LINKOPTION)

#sumAll: sumAll.o $(OBJECTS)
#	$(COMPILER) -o sumAll sumAll.o $(OBJECTS) $(LINKOPTION)

#waveform: waveform.o $(OBJECTS)
#	$(COMPILER) -o waveform waveform.o $(OBJECTS) $(LINKOPTION)

#DPPwaveform: DPPwaveform.o $(OBJECTS)
#	$(COMPILER) -o DPPwaveform DPPwaveform.o $(OBJECTS) $(LINKOPTION)

#text: text.o $(OBJECTS)
#	$(COMPILER) -o text text.o $(OBJECTS) $(LINKOPTION)

clean:
	rm -f *.o
