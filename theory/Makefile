################################### SETUP ######################################

# Define important locations
SHELL = /bin/sh
BIN = bin/
SOURCE = src/
ANALYSIS_SOURCE = ../analysis/src/

# Define compilation flags
COMPILER = clang++ 
CFLAGS = -lstdc++ -O3 -I$(shell root-config --incdir) -std=c++11

# For integration with ROOT data analysis framework
LINKOPTION = $(shell root-config --libs)

############################### DEFINE TARGETS #################################

# List all targets
TARGETS = opticalModel
all: $(addprefix $(BIN), $(TARGETS))

ANALYSIS_SOURCES = crossSection.cpp dataSet.cpp dataPoint.cpp CSUtilities.cpp

# Build optical model
OPTICAL_MODEL_SOURCES = opticalModel.cpp calculateCS.cpp scat.cpp reaction.cpp imaginaryForm.cpp expInt.cpp compound.cpp waves.cpp volume.cpp merson.cpp spinOrbit.cpp hartreeFock.cpp potPara.cpp surVolume.cpp OpticalPotential.cpp twoFermi.cpp asy.cpp gaussInteg.cpp disperse.cpp legendre.cpp surfaceTF.cpp surfaceSTD.cpp sphericalB.cpp standard.cpp coul.cpp whit.cpp channel.cpp level.cpp levelD.cpp compoundException.cpp gauss16.cpp Nucleus.cpp

$(BIN)opticalModel: $(addprefix $(SOURCE), $(OPTICAL_MODEL_SOURCES)) $(addprefix $(ANALYSIS_SOURCE), $(ANALYSIS_SOURCES))
	$(COMPILER) $(CFLAGS) -o $(BIN)opticalModel $(addprefix $(SOURCE), $(OPTICAL_MODEL_SOURCES)) $(addprefix $(ANALYSIS_SOURCE), $(ANALYSIS_SOURCES)) $(LINKOPTION)

# Build simple strongly-absorbing model
RAMSAUER_SOURCES = ramsauer.cpp
$(BIN)ramsauer: $(addprefix $(SOURCE), $(RAMSAUER_SOURCES)) $(addprefix $(ANALYSIS_SOURCE), $(ANALYSIS_SOURCES))
	$(COMPILER) $(CFLAGS) -o $(BIN)ramsauer $(addprefix $(SOURCE), $(RAMSAUER_SOURCES)) $(addprefix $(ANALYSIS_SOURCE), $(ANALYSIS_SOURCES)) $(LINKOPTION)


############################### DEFINE OBJECTS #################################

# Rule for building all targets
%.o : %.cpp
	$(COMPILER) $(CFLAGS) $< -o $@ $(LINKOPTION)

clean:
	rm -f *.o

