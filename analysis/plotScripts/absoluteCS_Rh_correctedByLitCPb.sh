#!/bin/bash

expFile="/data1/analysis/total.root"
expRhNatGraphName="RhNat"
relativeFile="/data1/analysis/relative.root"
correctionAverage="CPb_correctionAverage"

correctedFile="/data1/analysis/corrected.root"
correctedRhNatGraphName="RhNat"

../bin/applyCSCorrectionFactor "$relativeFile" "$correctionAverage" "$expFile" "$expRhNatGraphName" "$correctedFile" "$correctedRhNatGraphName"
