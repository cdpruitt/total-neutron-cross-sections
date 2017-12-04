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
TARGETS = OM_tl
all: $(addprefix $(BIN), $(TARGETS))

# Build driver (main data analysis engine)
OM_TL_SOURCES = OM_tl.cpp scat.cpp reaction.cpp imaginaryForm.cpp expInt.cpp compound.cpp waves.cpp volume.cpp merson.cpp spinOrbit.cpp hartreeFock.cpp potPara.cpp surVolume.cpp OpticalPotential.cpp twoFermi.cpp asy.cpp gaussInteg.cpp disperse.cpp legendre.cpp surfaceTF.cpp surfaceSTD.cpp sphericalB.cpp standard.cpp coul.cpp whit.cpp channel.cpp level.cpp levelD.cpp compoundException.cpp gauss16.cpp

$(BIN)OM_tl: $(addprefix $(SOURCE), $(OM_TL_SOURCES))
	$(COMPILER) $(CFLAGS) -o $(BIN)OM_tl $(addprefix $(SOURCE), $(OM_TL_SOURCES)) $(LINKOPTION)

############################### DEFINE OBJECTS #################################

# Rule for building all targets
%.o : %.cpp
	$(COMPILER) $(CFLAGS) $< -o $@ $(LINKOPTION)

clean:
	rm -f *.o

