#!/bin/bash

lowEFile="/data1/analysis/total_lowE.root"
highEFile="/data1/analysis/total_highE.root"
outputFile="/data1/analysis/merged.root"

CGraphName="CNat"
RhGraphName="RhNat"
PbGraphName="PbNat"
Ni58GraphName="Ni58"
NiNatGraphName="NiNat"
Ni64GraphName="Ni64"

juncture="6"

./bin/mergeCS "$lowEFile" "$CGraphName" "$highEFile" "$CGraphName" "$juncture" "$outputFile" "$CGraphName"
./bin/mergeCS "$lowEFile" "$RhGraphName" "$highEFile" "$RhGraphName" "$juncture" "$outputFile" "$RhGraphName"
./bin/mergeCS "$lowEFile" "$PbGraphName" "$highEFile" "$PbGraphName" "$juncture" "$outputFile" "$PbGraphName"
./bin/mergeCS "$lowEFile" "$Ni58GraphName" "$highEFile" "$Ni58GraphName" "$juncture" "$outputFile" "$Ni58GraphName"
./bin/mergeCS "$lowEFile" "$NiNatGraphName" "$highEFile" "$NiNatGraphName" "$juncture" "$outputFile" "$NiNatGraphName"
./bin/mergeCS "$lowEFile" "$Ni64GraphName" "$highEFile" "$Ni64GraphName" "$juncture" "$outputFile" "$Ni64GraphName"
