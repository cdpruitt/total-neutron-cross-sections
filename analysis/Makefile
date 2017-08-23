# Define important locations
SHELL = /bin/sh
BIN=./bin/
SOURCE=./src/
INCLUDE=./include/

# Define compilation flags
COMPILER= clang++ 
CFLAGS= -lstdc++ -O2 -I$(shell root-config --incdir) -std=c++11 -ggdb

# For integration with ROOT data analysis framework
LINKOPTION = $(shell root-config --libs)

# List all targets
TARGETS = driver sumSingle sumAll text finetimeCheck readLitData subtractCS multiplyCS relativeCS relative generateCSCorrection scaledownCS produceRunningRMS createRelativeCS
all: $(addprefix $(BIN), $(TARGETS))

# Build driver (main data analysis engine)
DRIVER_SOURCES = branches.cpp dataPoint.cpp dataSet.cpp driver.cpp experiment.cpp histos.cpp plots.cpp raw.cpp separate.cpp target.cpp veto.cpp waveform.cpp waveformFitting.cpp

$(BIN)driver: $(addprefix $(SOURCE), $(driver_SOURCES))
	$(COMPILER) $(CFLAGS) -o $(BIN)driver $(addprefix $(SOURCE), $(driver_SOURCES)) $(LINKOPTION)

# Build sumSingle (for generating cross sections using data from a single run) 
SUMSINGLE_SOURCES = sumSingle.cpp dataSet.cpp dataPoint.cpp CSPrereqs.cpp experiment.cpp target.cpp crossSection.cpp plots.cpp

$(BIN)sumSingle: $(addprefix $(SOURCE), $(SUMSINGLE_SOURCES))
	$(COMPILER) $(CFLAGS) -o $(BIN)sumSingle $(addprefix $(SOURCE), $(SUMSINGLE_SOURCES)) $(LINKOPTION)

# Build sumAll (for generating cross sections using data from all available runs)
SUMALL_SOURCES = sumAll.cpp dataSet.cpp dataPoint.cpp CSPrereqs.cpp experiment.cpp target.cpp crossSection.cpp plots.cpp

$(BIN)sumAll: $(addprefix $(SOURCE), $(SUMALL_SOURCES))
	$(COMPILER) $(CFLAGS) -o $(BIN)sumAll $(addprefix $(SOURCE), $(SUMALL_SOURCES)) $(LINKOPTION)

# Build readLitData (for reading pre-formatted literature data into ROOT Graphs)
$(BIN)readLitData: src/readLitData.cpp src/dataSet.cpp src/dataPoint.cpp
	$(COMPILER) $(CFLAGS) -o $(BIN)readLitData src/readLitData.cpp src/dataSet.cpp src/dataPoint.cpp $(LINKOPTION)

# Build text (for producing human-readable dump of raw event file data)
$(BIN)text: src/text.cpp src/raw.cpp src/branches.cpp
	$(COMPILER) $(CFLAGS) -o $(BIN)text src/text.cpp src/raw.cpp src/branches.cpp $(LINKOPTION)

# Build finetimeCheck (for comparing the timestamps of the same event, but recorded by different digitizer channels)
$(BIN)finetimeCheck: src/finetimeCheck.cpp src/branches.cpp
	$(COMPILER) $(CFLAGS) -o $(BIN)finetimeCheck src/finetimeCheck.cpp src/branches.cpp $(LINKOPTION)

# Build makeCSText (for writing graphed cross sections to formatted text files)
$(BIN)makeCSText: src/makeCSText.cpp src/dataPoint.cpp
	$(COMPILER) $(CFLAGS) -o $(BIN)makeCSText src/makeCSText.cpp src/dataPoint.cpp $(LINKOPTION)

# Build subtractCS (for taking the difference of two cross section graphs)
$(BIN)subtractCS: src/subtractCS.cpp src/dataSet.cpp src/dataPoint.cpp src/crossSection.cpp src/target.cpp src/plots.cpp
	$(COMPILER) $(CFLAGS) -o $(BIN)subtractCS src/subtractCS.cpp src/dataSet.cpp src/dataPoint.cpp src/crossSection.cpp src/target.cpp src/plots.cpp $(LINKOPTION)

# Build multiplyCS (for multiplying two cross section graphs)
$(BIN)multiplyCS: src/multiplyCS.cpp 
	$(COMPILER) $(CFLAGS) -o $(BIN)multiplyCS src/multiplyCS.cpp src/dataSet.cpp src/dataPoint.cpp src/crossSection.cpp src/target.cpp src/plots.cpp $(LINKOPTION)

# Build relativeCS (for calculating the relative difference of two cross section graphs)
$(BIN)relativeCS: src/relativeCS.cpp
	$(COMPILER) $(CFLAGS) -o $(BIN)relativeCS src/relativeCS.cpp src/dataSet.cpp src/dataPoint.cpp src/crossSection.cpp src/target.cpp src/plots.cpp $(LINKOPTION)

$(BIN)relative: src/dataPoint.cpp src/dataSet.cpp src/experiment.cpp src/crossSection.cpp src/relative.cpp
	$(COMPILER) $(CFLAGS) -o $(BIN)relative src/dataPoint.cpp src/dataSet.cpp src/experiment.cpp src/crossSection.cpp src/relative.cpp $(LINKOPTION)

$(BIN)createRelativeCS: src/createRelativeCS.cpp
	$(COMPILER) $(CFLAGS) -o $(BIN)createRelativeCS src/createRelativeCS.cpp src/dataSet.cpp src/dataPoint.cpp src/crossSection.cpp src/target.cpp src/plots.cpp $(LINKOPTION)

# Build generateCSCorrection (???)
$(BIN)generateCSCorrection: src/generateCSCorrection.cpp
	$(COMPILER) $(CFLAGS) -o $(BIN)generateCSCorrection src/generateCSCorrection.cpp src/dataSet.cpp src/dataPoint.cpp src/crossSection.cpp src/target.cpp src/plots.cpp $(LINKOPTION)

# Build scaledownCS (for rebinning a cross section with a coarser bin size)
$(BIN)scaledownCS: src/scaledownCS.cpp src/dataSet.cpp
	$(COMPILER) $(CFLAGS) -o $(BIN)scaledownCS src/scaledownCS.cpp src/dataSet.cpp src/dataPoint.cpp src/crossSection.cpp src/target.cpp src/plots.cpp $(LINKOPTION)

# Build produceRunningRMS (for plotting the running root-mean-squared difference between two cross
# section graphs)
$(BIN)produceRunningRMS: src/produceRunningRMS.cpp
	$(COMPILER) $(CFLAGS) -o $(BIN)produceRunningRMS src/produceRunningRMS.cpp src/dataSet.cpp src/dataPoint.cpp src/crossSection.cpp src/target.cpp src/plots.cpp $(LINKOPTION)

# Rule for building all targets
%.o : %.cpp
	$(COMPILER) $(CFLAGS) $< -o $@ $(LINKOPTION)

clean:
	rm -f *.o
