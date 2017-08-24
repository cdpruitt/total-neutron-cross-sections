#!/bin/bash

scaledown=$1

expFile="/data2/analysis/total.root"
expPbGraphName="NatPb"
expCGraphName="longCarbon"

litFile="/data2/analysis/literatureData.root"
litPbGraphName="Natural Pb (n,tot)"
litCGraphName="Natural carbon (n,tot)"

relPbGraphName="relDiff_natPb_expLit"
relCGraphName="relDiff_natC_expLit"

../relativeCS "$expFile" "$expCGraphName" "$litFile" "$litCGraphName" "$relCGraphName"
../relativeCS "$expFile" "$expPbGraphName" "$litFile" "$litPbGraphName" "$relPbGraphName"

correctionAverage="CPb_correctionAverage"

../subtractCS "$expFile" "$relCGraphName" "$expFile" "$relPbGraphName" -1 2 "$correctionAverage"

expSn112GraphName="Sn112"
expSnNatGraphName="NatSn"
expSn124GraphName="Sn124"

correctedFile="/data2/analysis/corrected.root"
correctedSn112GraphName="Sn112_corrected"
correctedSnNatGraphName="SnNat_corrected"
correctedSn124GraphName="Sn124_corrected"

../applyCSCorrectionFactor "$expFile" "$correctionAverage" "$expFile" "$expSn112GraphName" "$correctedFile" "$correctedSn112GraphName"
../applyCSCorrectionFactor "$expFile" "$correctionAverage" "$expFile" "$expSnNatGraphName" "$correctedFile" "$correctedSnNatGraphName"
../applyCSCorrectionFactor "$expFile" "$correctionAverage" "$expFile" "$expSn124GraphName" "$correctedFile" "$correctedSn124GraphName"

litFile="/data2/analysis/literatureData.root"
litSnGraphName="Natural Sn (n,tot)"

scaledownSn112GraphName="Sn112_corrected"
scaledownSnNatGraphName="SnNat_corrected"
scaledownSn124GraphName="Sn124_corrected"

if [ "$scaledown" -gt 1 ]
then
    printf "Scaling down relative CS plots by factor of $scaledown.\n"

    scaledownFile="/data2/analysis/scaledown.root"
    ../scaledownCS "$correctedFile" "$correctedSn112GraphName" "$scaledown" "$scaledownFile" "$scaledownSn112GraphName"
    ../scaledownCS "$correctedFile" "$correctedSnNatGraphName" "$scaledown" "$scaledownFile" "$scaledownSnNatGraphName"
    ../scaledownCS "$correctedFile" "$correctedSn124GraphName" "$scaledown" "$scaledownFile" "$scaledownSn124GraphName"
fi
