#!/bin/bash

expFile="/data2/analysis/total.root"
expCGraphName="CNat"
expSnNatGraphName="SnNat"
expPbGraphName="PbNat"
expSn112GraphName="Sn112"
expSn124GraphName="Sn124"

litFile="/data2/analysis/literatureData.root"
litPbGraphName="Natural Pb (n,tot)"
litCGraphName="Natural C (n,tot)"

relativeFile="/data2/analysis/relative.root"
relPbGraphName="NatPb, expLit absolute"
relCGraphName="NatC, expLit absolute"

correctionAverage="CPb_correctionAverage"

../bin/subtractCS "$relativeFile" "$relCGraphName" "$relativeFile" "$relPbGraphName" -1 2 "$correctionAverage"

correctedFile="/data2/analysis/corrected.root"
correctedCNatGraphName="CNat"
correctedSnNatGraphName="SnNat"
correctedPbNatGraphName="PbNat"
correctedSn112GraphName="Sn112"
correctedSn124GraphName="Sn124"

../bin/applyCSCorrectionFactor "$relativeFile" "$correctionAverage" "$expFile" "$expCGraphName" "$correctedFile" "$correctedCNatGraphName"
../bin/applyCSCorrectionFactor "$relativeFile" "$correctionAverage" "$expFile" "$expSnNatGraphName" "$correctedFile" "$correctedSnNatGraphName"
../bin/applyCSCorrectionFactor "$relativeFile" "$correctionAverage" "$expFile" "$expPbGraphName" "$correctedFile" "$correctedPbNatGraphName"
../bin/applyCSCorrectionFactor "$relativeFile" "$correctionAverage" "$expFile" "$expSn112GraphName" "$correctedFile" "$correctedSn112GraphName"
../bin/applyCSCorrectionFactor "$relativeFile" "$correctionAverage" "$expFile" "$expSn124GraphName" "$correctedFile" "$correctedSn124GraphName"



