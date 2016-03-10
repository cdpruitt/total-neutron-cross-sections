OBJECTS = 
ALLOBJECTS := $(patsubst %.cpp,%.o,$(wildcard *.cpp))
CFLAGS= -c -O2 -I$(shell root-config --incdir) --std=c++0x
COMPILER= g++
LINKOPTION = $(shell root-config --libs)

all: raw resort histos

raw: raw.o $(OBJECTS)
	$(COMPILER) -o raw raw.o $(OBJECTS) $(LINKOPTION)

resort: resort.o $(OBJECTS)
	$(COMPILER) -o resort resort.o $(OBJECTS) $(LINKOPTION)

histos: histos.o $(OBJECTS)
	$(COMPILER) -o histos histos.o $(OBJECTS) $(LINKOPTION)

#finetime: finetime.o $(OBJECTS)
#	$(COMPILER) -o finetime finetime.o $(OBJECTS) $(LINKOPTION)

SumHistos: SumHistos.o $(OBJECTS)
	$(COMPILER) -o SumHistos SumHistos.o $(OBJECTS) $(LINKOPTION)

$(ALLOBJECTS): %.o : %.cpp
	$(COMPILER) $(CFLAGS) $< -o $@

clean:
	rm -f *.o
