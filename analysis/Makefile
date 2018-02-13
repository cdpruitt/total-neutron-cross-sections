################################### SETUP ######################################

# Define important locations
SHELL = /bin/sh
BIN = bin/
SOURCE = src/

# Define compilation flags
COMPILER = clang++ 
CFLAGS = -lstdc++ -O3 -I$(shell root-config --incdir) -std=c++11

# For integration with ROOT data analysis framework
LINKOPTION = $(shell root-config --libs)

############################### DEFINE TARGETS #################################

# List all targets
TARGETS = driver text sumAll eachSubrun sumChunk readLitData subtractCS mergeCS shiftCS multiplyCS relativeDiffCS relativeCS applyCSCorrectionFactor scaledownCS produceRunningRMS detTimeCheck rateHisto #plotCSPrereqs
all: $(addprefix $(BIN), $(TARGETS))

# Build driver (main data analysis engine)
DRIVER_SOURCES = dataPoint.cpp dataSet.cpp driver.cpp config.cpp experiment.cpp fillBasicHistos.cpp fillCSHistos.cpp plots.cpp raw.cpp identifyMacropulses.cpp assignEventsToMacropulses.cpp calculateGammaCorrection.cpp correctForDeadtime.cpp produceEnergyHistos.cpp target.cpp veto.cpp softwareCFD.cpp identifyGoodMacros.cpp correctForBackground.cpp

$(BIN)driver: $(addprefix $(SOURCE), $(DRIVER_SOURCES))
	$(COMPILER) $(CFLAGS) -o $(BIN)driver $(addprefix $(SOURCE), $(DRIVER_SOURCES)) $(LINKOPTION)

# Build sumAll (for generating cross sections using data from all available runs)
SUMALL_SOURCES = sumAll.cpp dataSet.cpp dataPoint.cpp CSPrereqs.cpp config.cpp experiment.cpp target.cpp crossSection.cpp plots.cpp CSUtilities.cpp
$(BIN)sumAll: $(addprefix $(SOURCE), $(SUMALL_SOURCES))
	$(COMPILER) $(CFLAGS) -o $(BIN)sumAll $(addprefix $(SOURCE), $(SUMALL_SOURCES)) $(LINKOPTION)

# Build eachSubrun (for generating cross sections using data from all available runs)
EACHSUBRUN_SOURCES = eachSubrun.cpp dataSet.cpp dataPoint.cpp CSPrereqs.cpp config.cpp experiment.cpp target.cpp crossSection.cpp plots.cpp CSUtilities.cpp
$(BIN)eachSubrun: $(addprefix $(SOURCE), $(EACHSUBRUN_SOURCES))
	$(COMPILER) $(CFLAGS) -o $(BIN)eachSubrun $(addprefix $(SOURCE), $(EACHSUBRUN_SOURCES)) $(LINKOPTION)

# Build sumChunk (for generating cross sections using data from a select set of subruns)
SUMCHUNK_SOURCES = sumChunk.cpp dataSet.cpp dataPoint.cpp CSPrereqs.cpp config.cpp experiment.cpp target.cpp crossSection.cpp plots.cpp CSUtilities.cpp
$(BIN)sumChunk: $(addprefix $(SOURCE), $(SUMCHUNK_SOURCES))
	$(COMPILER) $(CFLAGS) -o $(BIN)sumChunk $(addprefix $(SOURCE), $(SUMCHUNK_SOURCES)) $(LINKOPTION)

# Build plotCSPrereqs (for generating cross sections using data from a select set of subruns)
PLOTCSPREREQS_SOURCES = plotCSPrereqs.cpp dataSet.cpp dataPoint.cpp CSPrereqs.cpp config.cpp experiment.cpp target.cpp crossSection.cpp plots.cpp CSUtilities.cpp
$(BIN)plotCSPrereqs: $(addprefix $(SOURCE), $(PLOTCSPREREQS_SOURCES))
	$(COMPILER) $(CFLAGS) -o $(BIN)plotCSPrereqs $(addprefix $(SOURCE), $(PLOTCSPREREQS_SOURCES)) $(LINKOPTION)

# Build readLitData (for reading pre-formatted literature data into ROOT Graphs)
READLITDATA_SOURCES = readLitData.cpp dataSet.cpp dataPoint.cpp 
$(BIN)readLitData: $(addprefix $(SOURCE), $(READLITDATA_SOURCES))
	$(COMPILER) $(CFLAGS) -o $(BIN)readLitData $(addprefix $(SOURCE), $(READLITDATA_SOURCES)) $(LINKOPTION)

# Build text (for producing human-readable dump of raw event file data)
TEXT_SOURCES = text.cpp raw.cpp softwareCFD.cpp
$(BIN)text: $(addprefix $(SOURCE), $(TEXT_SOURCES))
	$(COMPILER) $(CFLAGS) -o $(BIN)text $(addprefix $(SOURCE), $(TEXT_SOURCES)) $(LINKOPTION)

# Build detTimeCheck (for comparing the timestamps of the same event, but recorded by different digitizer channels)
DETTIMECHECK_SOURCES = detTimeCheck.cpp raw.cpp config.cpp experiment.cpp softwareCFD.cpp
$(BIN)detTimeCheck: $(addprefix $(SOURCE), $(DETTIMECHECK_SOURCES))
	$(COMPILER) $(CFLAGS) -o $(BIN)detTimeCheck $(addprefix $(SOURCE), $(DETTIMECHECK_SOURCES)) $(LINKOPTION)

# Build makeCSText (for writing graphed cross sections to formatted text files)
MAKECSTEXT_SOURCES = makeCSText.cpp dataPoint.cpp
$(BIN)makeCSText: $(addprefix $(SOURCE), $(MAKECSTEXT_SOURCES))
	$(COMPILER) $(CFLAGS) -o $(BIN)makeCSText $(addprefix $(SOURCE), $(MAKECSTEXT_SOURCES)) $(LINKOPTION)

# Build subtractCS (for taking the difference of two cross section graphs)
SUBTRACTCS_SOURCES = subtractCS.cpp dataSet.cpp dataPoint.cpp crossSection.cpp target.cpp plots.cpp CSUtilities.cpp
$(BIN)subtractCS: $(addprefix $(SOURCE), $(SUBTRACTCS_SOURCES)) 
	$(COMPILER) $(CFLAGS) -o $(BIN)subtractCS $(addprefix $(SOURCE), $(SUBTRACTCS_SOURCES)) $(LINKOPTION)

# Build mergeCS (for taking the difference of two cross section graphs)
MERGECS_SOURCES = mergeCS.cpp dataSet.cpp dataPoint.cpp crossSection.cpp target.cpp plots.cpp CSUtilities.cpp
$(BIN)mergeCS: $(addprefix $(SOURCE), $(MERGECS_SOURCES)) 
	$(COMPILER) $(CFLAGS) -o $(BIN)mergeCS $(addprefix $(SOURCE), $(MERGECS_SOURCES)) $(LINKOPTION)

# Build shiftCS (for taking the difference of two cross section graphs)
SHIFTCS_SOURCES = shiftCS.cpp dataSet.cpp dataPoint.cpp crossSection.cpp target.cpp plots.cpp CSUtilities.cpp
$(BIN)shiftCS: $(addprefix $(SOURCE), $(SHIFTCS_SOURCES)) 
	$(COMPILER) $(CFLAGS) -o $(BIN)shiftCS $(addprefix $(SOURCE), $(SHIFTCS_SOURCES)) $(LINKOPTION)

# Build multiplyCS (for multiplying two cross section graphs)
MULTIPLYCS_SOURCES = multiplyCS.cpp dataSet.cpp dataPoint.cpp crossSection.cpp target.cpp plots.cpp CSUtilities.cpp
$(BIN)multiplyCS: $(addprefix $(SOURCE), $(MULTIPLYCS_SOURCES)) 
	$(COMPILER) $(CFLAGS) -o $(BIN)multiplyCS $(addprefix $(SOURCE), $(MULTIPLYCS_SOURCES)) $(LINKOPTION)

# Build relativeCS (for calculating the absolute difference of two cross section graphs)
RELATIVECS_SOURCES = relativeCS.cpp dataSet.cpp dataPoint.cpp crossSection.cpp target.cpp plots.cpp CSUtilities.cpp
$(BIN)relativeCS: $(addprefix $(SOURCE), $(RELATIVECS_SOURCES)) 
	$(COMPILER) $(CFLAGS) -o $(BIN)relativeCS $(addprefix $(SOURCE), $(RELATIVECS_SOURCES)) $(LINKOPTION)

# Build relativeDiffCS (for calculating the relative difference of two cross section graphs)
RELATIVEDIFFCS_SOURCES = relativeDiffCS.cpp dataSet.cpp dataPoint.cpp crossSection.cpp target.cpp plots.cpp CSUtilities.cpp
$(BIN)relativeDiffCS: $(addprefix $(SOURCE), $(RELATIVEDIFFCS_SOURCES)) 
	$(COMPILER) $(CFLAGS) -o $(BIN)relativeDiffCS $(addprefix $(SOURCE), $(RELATIVEDIFFCS_SOURCES)) $(LINKOPTION)

# Build applyCSCorrectionFactor (for scaling each point in a cross section by a factor)
APPLYCSCORRECTIONFACTOR_SOURCES = applyCSCorrectionFactor.cpp dataSet.cpp dataPoint.cpp crossSection.cpp target.cpp plots.cpp CSUtilities.cpp
$(BIN)applyCSCorrectionFactor: $(addprefix $(SOURCE), $(APPLYCSCORRECTIONFACTOR_SOURCES)) 
	$(COMPILER) $(CFLAGS) -o $(BIN)applyCSCorrectionFactor $(addprefix $(SOURCE), $(APPLYCSCORRECTIONFACTOR_SOURCES)) $(LINKOPTION)

# Build scaledownCS (for rebinning a cross section with a coarser bin size)
SCALEDOWNCS_SOURCES = scaledownCS.cpp dataSet.cpp dataPoint.cpp crossSection.cpp target.cpp plots.cpp CSUtilities.cpp
$(BIN)scaledownCS: $(addprefix $(SOURCE), $(SCALEDOWNCS_SOURCES)) 
	$(COMPILER) $(CFLAGS) -o $(BIN)scaledownCS $(addprefix $(SOURCE), $(SCALEDOWNCS_SOURCES)) $(LINKOPTION)

# Build produceRunningRMS (for plotting the running root-mean-squared difference between two cross
# section graphs)
PRODUCERUNNINGRMS_SOURCES = produceRunningRMS.cpp dataSet.cpp dataPoint.cpp crossSection.cpp target.cpp plots.cpp CSUtilities.cpp
$(BIN)produceRunningRMS: $(addprefix $(SOURCE), $(PRODUCERUNNINGRMS_SOURCES)) 
	$(COMPILER) $(CFLAGS) -o $(BIN)produceRunningRMS $(addprefix $(SOURCE), $(PRODUCERUNNINGRMS_SOURCES)) $(LINKOPTION)

# Build rateHisto
RATEHISTO_SOURCES = rateHisto.cpp
$(BIN)rateHisto: $(addprefix $(SOURCE), $(RATEHISTO_SOURCES)) 
	$(COMPILER) $(CFLAGS) -o $(BIN)rateHisto $(addprefix $(SOURCE), $(RATEHISTO_SOURCES)) $(LINKOPTION)

# Build TOFSimulation
TOFSIMULATION_SOURCES = TOFSimulation.cpp correctForDeadtime.cpp
$(BIN)TOFSimulation: $(addprefix $(SOURCE), $(TOFSIMULATION_SOURCES)) 
	$(COMPILER) $(CFLAGS) -o $(BIN)TOFSimulation $(addprefix $(SOURCE), $(TOFSIMULATION_SOURCES)) $(LINKOPTION)

############################### DEFINE OBJECTS #################################

# Rule for building all targets
%.o : %.cpp
	$(COMPILER) $(CFLAGS) $< -o $@ $(LINKOPTION)

clean:
	rm -f *.o
