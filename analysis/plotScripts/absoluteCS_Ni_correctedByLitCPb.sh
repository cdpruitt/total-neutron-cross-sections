#!/bin/bash

expFile="/data1/analysis/total.root"
expCGraphName="CNat"
expNiNatGraphName="NiNat"
expPbGraphName="PbNat"
expNi58GraphName="Ni58"
expNi64GraphName="Ni64"

litFile="/data1/analysis/literatureData.root"
litPbGraphName="Natural Pb (n,tot)"
litCGraphName="Natural C (n,tot)"

relativeFile="/data1/analysis/relative.root"
relPbGraphName="PbNat, expLit, absolute"
relCGraphName="CNat, expLit, absolute"

correctionAverage="CPb_correctionAverage"

../bin/subtractCS "$relativeFile" "$relCGraphName" "$relativeFile" "$relPbGraphName" -1 2 "$correctionAverage"

correctedFile="/data1/analysis/corrected.root"
correctedCNatGraphName="CNat"
correctedNiNatGraphName="NiNat"
correctedPbNatGraphName="PbNat"
correctedNi58GraphName="Ni58"
correctedNi64GraphName="Ni64"

../bin/applyCSCorrectionFactor "$relativeFile" "$correctionAverage" "$expFile" "$expCGraphName" "$correctedFile" "$correctedCNatGraphName"
../bin/applyCSCorrectionFactor "$relativeFile" "$correctionAverage" "$expFile" "$expNiNatGraphName" "$correctedFile" "$correctedNiNatGraphName"
../bin/applyCSCorrectionFactor "$relativeFile" "$correctionAverage" "$expFile" "$expPbGraphName" "$correctedFile" "$correctedPbNatGraphName"
../bin/applyCSCorrectionFactor "$relativeFile" "$correctionAverage" "$expFile" "$expNi58GraphName" "$correctedFile" "$correctedNi58GraphName"
../bin/applyCSCorrectionFactor "$relativeFile" "$correctionAverage" "$expFile" "$expNi64GraphName" "$correctedFile" "$correctedNi64GraphName"
