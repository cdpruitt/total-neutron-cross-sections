#!/bin/bash

expFile="/data1/analysis/total.root"
expCGraphName="CNat"
expNiNatGraphName="NiNat"
expPbGraphName="PbNat"

litFile="/data1/analysis/literatureData.root"
litPbGraphName="Natural Pb (n,tot)"
litCGraphName="Natural C (n,tot)"

relativeFile="/data1/analysis/relative.root"
relPbGraphName="NatPb, expLit absolute"
relCGraphName="NatC, expLit absolute"

correctionAverage="CPb_correctionAverage"

../bin/subtractCS "$relativeFile" "$relCGraphName" "$relativeFile" "$relPbGraphName" -1 2 "$correctionAverage"

correctedFile="/data1/analysis/corrected.root"
correctedCNatGraphName="CNat"
correctedNiNatGraphName="NiNat"
correctedPbNatGraphName="PbNat"

../bin/applyCSCorrectionFactor "$relativeFile" "$correctionAverage" "$expFile" "$expCGraphName" "$correctedFile" "$correctedCNatGraphName"
../bin/applyCSCorrectionFactor "$relativeFile" "$correctionAverage" "$expFile" "$expNiNatGraphName" "$correctedFile" "$correctedNiNatGraphName"
../bin/applyCSCorrectionFactor "$relativeFile" "$correctionAverage" "$expFile" "$expPbGraphName" "$correctedFile" "$correctedPbNatGraphName"
