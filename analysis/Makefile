OBJECTS = 
ALLOBJECTS := $(patsubst %.cpp,%.o,$(wildcard *.cpp))
CFLAGS= -c -O0 -I$(shell root-config --incdir) -std=c++11 -Wall -Wextra -ggdb
COMPILER= clang++ 
LINKOPTION = $(shell root-config --libs)

all: raw resort histos waveform sumRun sumAll plotMonitor

raw: raw.o $(OBJECTS)
	$(COMPILER) -o raw raw.o $(OBJECTS) $(LINKOPTION)

resort: resort.o $(OBJECTS)
	$(COMPILER) -o resort resort.o $(OBJECTS) $(LINKOPTION)

histos: histos.o $(OBJECTS)
	$(COMPILER) -o histos histos.o $(OBJECTS) $(LINKOPTION)

#finetime: finetime.o $(OBJECTS)
#	$(COMPILER) -o finetime finetime.o $(OBJECTS) $(LINKOPTION)

plotMonitor: plotMonitor.o $(OBJECTS)
	$(COMPILER) -o plotMonitor plotMonitor.o $(OBJECTS) $(LINKOPTION)

sumRun: sumRun.o $(OBJECTS)
	$(COMPILER) -o sumRun sumRun.o $(OBJECTS) $(LINKOPTION)

sumAll: sumAll.o $(OBJECTS)
	$(COMPILER) -o sumAll sumAll.o $(OBJECTS) $(LINKOPTION)

waveform: waveform.o $(OBJECTS)
	$(COMPILER) -o waveform waveform.o $(OBJECTS) $(LINKOPTION)

$(ALLOBJECTS): %.o : %.cpp
	$(COMPILER) $(CFLAGS) $< -o $@

clean:
	rm -f *.o
