################################### SETUP ######################################

# Define important locations
SHELL = /bin/sh
BIN = bin/
SOURCE = src/

# Define compilation flags
COMPILER = clang++ 
CFLAGS = -lstdc++ -I$(shell root-config --incdir) -std=c++11 -ggdb

# For integration with ROOT data analysis framework
LINKOPTION = $(shell root-config --libs)

############################### DEFINE TARGETS #################################

# List all targets
TARGETS = driver sumSingle sumAll text fineTimeCheck readLitData subtractCS multiplyCS relativeCS applyCSCorrectionFactor scaledownCS produceRunningRMS
all: $(addprefix $(BIN), $(TARGETS))

# Build driver (main data analysis engine)
DRIVER_SOURCES = branches.cpp dataPoint.cpp dataSet.cpp driver.cpp experiment.cpp histos.cpp plots.cpp raw.cpp assignMacropulses.cpp target.cpp veto.cpp waveform.cpp waveformFitting.cpp
$(BIN)driver: $(addprefix $(SOURCE), $(DRIVER_SOURCES))
	$(COMPILER) $(CFLAGS) -o $(BIN)driver $(addprefix $(SOURCE), $(DRIVER_SOURCES)) $(LINKOPTION)

# Build sumSingle (for generating cross sections using data from a single run) 
SUMSINGLE_SOURCES = sumSingle.cpp dataSet.cpp dataPoint.cpp CSPrereqs.cpp experiment.cpp target.cpp crossSection.cpp plots.cpp
$(BIN)sumSingle: $(addprefix $(SOURCE), $(SUMSINGLE_SOURCES))
	$(COMPILER) $(CFLAGS) -o $(BIN)sumSingle $(addprefix $(SOURCE), $(SUMSINGLE_SOURCES)) $(LINKOPTION)

# Build sumAll (for generating cross sections using data from all available runs)
SUMALL_SOURCES = sumAll.cpp dataSet.cpp dataPoint.cpp CSPrereqs.cpp experiment.cpp target.cpp crossSection.cpp plots.cpp
$(BIN)sumAll: $(addprefix $(SOURCE), $(SUMALL_SOURCES))
	$(COMPILER) $(CFLAGS) -o $(BIN)sumAll $(addprefix $(SOURCE), $(SUMALL_SOURCES)) $(LINKOPTION)

# Build readLitData (for reading pre-formatted literature data into ROOT Graphs)
READLITDATA_SOURCES = readLitData.cpp dataSet.cpp dataPoint.cpp 
$(BIN)readLitData: $(addprefix $(SOURCE), $(READLITDATA_SOURCES))
	$(COMPILER) $(CFLAGS) -o $(BIN)readLitData $(addprefix $(SOURCE), $(READLITDATA_SOURCES)) $(LINKOPTION)

# Build text (for producing human-readable dump of raw event file data)
TEXT_SOURCES = text.cpp raw.cpp branches.cpp
$(BIN)text: $(addprefix $(SOURCE), $(TEXT_SOURCES))
	$(COMPILER) $(CFLAGS) -o $(BIN)text $(addprefix $(SOURCE), $(TEXT_SOURCES)) $(LINKOPTION)

# Build finetimeCheck (for comparing the timestamps of the same event, but recorded by different digitizer channels)
FINETIMECHECK_SOURCES = fineTimeCheck.cpp raw.cpp assignMacropulses.cpp branches.cpp experiment.cpp
$(BIN)fineTimeCheck: $(addprefix $(SOURCE), $(FINETIMECHECK_SOURCES))
	$(COMPILER) $(CFLAGS) -o $(BIN)fineTimeCheck $(addprefix $(SOURCE), $(FINETIMECHECK_SOURCES)) $(LINKOPTION)

# Build makeCSText (for writing graphed cross sections to formatted text files)
MAKECSTEXT_SOURCES = makeCSText.cpp dataPoint.cpp
$(BIN)makeCSText: $(addprefix $(SOURCE), $(MAKECSTEXT_SOURCES))
	$(COMPILER) $(CFLAGS) -o $(BIN)makeCSText $(addprefix $(SOURCE), $(MAKECSTEXT_SOURCES)) $(LINKOPTION)

# Build subtractCS (for taking the difference of two cross section graphs)
SUBTRACTCS_SOURCES = subtractCS.cpp dataSet.cpp dataPoint.cpp crossSection.cpp target.cpp plots.cpp
$(BIN)subtractCS: $(addprefix $(SOURCE), $(SUBTRACTCS_SOURCES)) 
	$(COMPILER) $(CFLAGS) -o $(BIN)subtractCS $(addprefix $(SOURCE), $(SUBTRACTCS_SOURCES)) $(LINKOPTION)

# Build multiplyCS (for multiplying two cross section graphs)
MULTIPLYCS_SOURCES = multiplyCS.cpp dataSet.cpp dataPoint.cpp crossSection.cpp target.cpp plots.cpp
$(BIN)multiplyCS: $(addprefix $(SOURCE), $(MULTIPLYCS_SOURCES)) 
	$(COMPILER) $(CFLAGS) -o $(BIN)multiplyCS $(addprefix $(SOURCE), $(MULTIPLYCS_SOURCES)) $(LINKOPTION)

# Build relativeCS (for calculating the relative difference of two cross section graphs)
RELATIVECS_SOURCES = relativeCS.cpp dataSet.cpp dataPoint.cpp crossSection.cpp target.cpp plots.cpp
$(BIN)relativeCS: $(addprefix $(SOURCE), $(RELATIVECS_SOURCES)) 
	$(COMPILER) $(CFLAGS) -o $(BIN)relativeCS $(addprefix $(SOURCE), $(RELATIVECS_SOURCES)) $(LINKOPTION)

# Build applyCSCorrectionFactor (for scaling each point in a cross section by a factor)
APPLYCSCORRECTIONFACTOR_SOURCES = applyCSCorrectionFactor.cpp dataSet.cpp dataPoint.cpp crossSection.cpp target.cpp plots.cpp
$(BIN)applyCSCorrectionFactor: $(addprefix $(SOURCE), $(APPLYCSCORRECTIONFACTOR_SOURCES)) 
	$(COMPILER) $(CFLAGS) -o $(BIN)applyCSCorrectionFactor $(addprefix $(SOURCE), $(APPLYCSCORRECTIONFACTOR_SOURCES)) $(LINKOPTION)

# Build scaledownCS (for rebinning a cross section with a coarser bin size)
SCALEDOWNCS_SOURCES = scaledownCS.cpp dataSet.cpp dataPoint.cpp crossSection.cpp target.cpp plots.cpp
$(BIN)scaledownCS: $(addprefix $(SOURCE), $(SCALEDOWNCS_SOURCES)) 
	$(COMPILER) $(CFLAGS) -o $(BIN)scaledownCS $(addprefix $(SOURCE), $(SCALEDOWNCS_SOURCES)) $(LINKOPTION)

# Build produceRunningRMS (for plotting the running root-mean-squared difference between two cross
# section graphs)
PRODUCERUNNINGRMS_SOURCES = produceRunningRMS.cpp dataSet.cpp dataPoint.cpp crossSection.cpp target.cpp plots.cpp
$(BIN)produceRunningRMS: $(addprefix $(SOURCE), $(PRODUCERUNNINGRMS_SOURCES)) 
	$(COMPILER) $(CFLAGS) -o $(BIN)produceRunningRMS $(addprefix $(SOURCE), $(PRODUCERUNNINGRMS_SOURCES)) $(LINKOPTION)

############################### DEFINE OBJECTS #################################

# Rule for building all targets
%.o : %.cpp
	$(COMPILER) $(CFLAGS) $< -o $@ $(LINKOPTION)

clean:
	rm -f *.o
