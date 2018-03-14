#!/bin/bash

expFile="/data1/analysis/total.root"

NiNatGraphName="NiNat"
Ni64GraphName="Ni64"

NiNatGraphSysErrorsName="${NiNatGraphName}SysErrors"
Ni64GraphSysErrorsName="${Ni64GraphName}SysErrors"

purity=0.922
impurity=$(bc <<< "1-$purity")

../bin/subtractCS "$expFile" "$Ni64GraphName" "$expFile" "$NiNatGraphName" "$impurity" "$purity" "Ni64PurityCorrected"

../bin/subtractCS "$expFile" "$Ni64GraphSysErrorsName" "$expFile" "$NiNatGraphSysErrorsName" "$impurity" "$purity" "Ni64PurityCorrectedSysErrors"
