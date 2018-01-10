#!/bin/bash

expFile="/data2/analysis/total.root"
expCGraphName="CNat"
expPbGraphName="PbNat"

expSn112GraphName="Sn112"
expSnNatGraphName="SnNat"
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
correctedSn112GraphName="Sn112_corrected"
correctedSnNatGraphName="SnNat_corrected"
correctedSn124GraphName="Sn124_corrected"

../bin/applyCSCorrectionFactor "$relativeFile" "$correctionAverage" "$expFile" "$expSn112GraphName" "$correctedFile" "$correctedSn112GraphName"
../bin/applyCSCorrectionFactor "$relativeFile" "$correctionAverage" "$expFile" "$expSnNatGraphName" "$correctedFile" "$correctedSnNatGraphName"
../bin/applyCSCorrectionFactor "$relativeFile" "$correctionAverage" "$expFile" "$expSn124GraphName" "$correctedFile" "$correctedSn124GraphName"
