#!/bin/bash

expFile="/data2/analysis/total.root"
expH2OGraphName="H2O"
expD2OGraphName="D2O"
expH2O18GraphName="H2O18"

relativeFile="/data2/analysis/relative.root"
correctionAverage="CPb_correctionAverage"

correctedFile="/data2/analysis/corrected.root"
correctedH2OGraphName="H2O"
correctedD2OGraphName="D2O"
correctedH2O18GraphName="H2O18"

../bin/applyCSCorrectionFactor "$relativeFile" "$correctionAverage" "$expFile" "$expH2OGraphName" "$correctedFile" "$correctedH2OGraphName"
../bin/applyCSCorrectionFactor "$relativeFile" "$correctionAverage" "$expFile" "$expD2OGraphName" "$correctedFile" "$correctedD2OGraphName"
../bin/applyCSCorrectionFactor "$relativeFile" "$correctionAverage" "$expFile" "$expH2O18GraphName" "$correctedFile" "$correctedH2O18GraphName"
