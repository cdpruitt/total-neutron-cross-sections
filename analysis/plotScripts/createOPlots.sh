#!/bin/bash

./createLitWaterCS.sh # using the binning listed in the experimental
                      # configuration, create literature H2O/D2O CSs

./shiftO18Lit.sh # move up literature O18 data by 1 barn for visibility

./produceOTotalCS.sh  # using literature H and D data, create exp NatO/O18 CSs

./relativeDiffLit_O.sh # calculate relative difference between literature O/exp
                       # O cross sections

./relativeDiff_D2OH2O.sh # calculate relative difference between raw experimental
                         # D2O/H2O cross sections

./relativeDiff_O18O16.sh # calculate relative difference between
                         # experimental NatO/O18 cross sections

#./absoluteCS_O_correctedByLitCPb.sh # apply Sn dataset C/Pb correction to O CSs

#./produceOTotalCS_Corrected.sh # using literature H and D data, create exp
                               # NatO/O18 CSs

#./relativeDiffLit_CorrectedO.sh # calculate relative difference between
                                # literature O/exp O cross sections, with C/Pb
                                # correction applied to experimental data

#./relativeDiff_CorrectedO18O16.sh # calculate relative difference between
                                     # NatO/O18 cross sections, with C/Pb
                                     # correction applied to experimental data
