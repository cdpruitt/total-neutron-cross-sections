#include "../include/config.h"
#include "../include/CSUtilities.h"
#include "../include/crossSection.h"
#include "../include/dataSet.h"
#include "../include/dataPoint.h"
#include "../include/experiment.h"
#include "../include/physicalConstants.h"
#include "../include/plots.h"

#include <iostream>
#include <fstream>

using namespace std;

CrossSection mergeCrossSections(CrossSection firstCS, double juncture,
        CrossSection secondCS)
{
    CrossSection outputCS;
    DataSet outputDS;

    DataSet firstDS = firstCS.getDataSet();
    DataSet secondDS = secondCS.getDataSet();

    for(int i=0; i<firstDS.getNumberOfPoints(); i++)
    {
        if(firstDS.getPoint(i).getXValue()>juncture)
        {
            break;
        }

        outputDS.addPoint(firstDS.getPoint(i));
    }

    for(int i=0; i<secondDS.getNumberOfPoints(); i++)
    {
        if(secondDS.getPoint(i).getXValue()<juncture)
        {
            continue;
        }

        outputDS.addPoint(secondDS.getPoint(i));
    }

    outputCS.addDataSet(outputDS);

    return outputCS;
}

int readLitData(string litDirectory, string litOutputName, const Config& config)
{
    string inFileName = litDirectory + "/filesToRead.txt";
    ifstream inFile(inFileName);
    if(!inFile.is_open())
    {
        cout << "Failed to open " << inFileName << endl;
        return(1);
    }

    string dummy;
    vector<string> fileNames;

    while(inFile >> dummy)
    {
        dummy = litDirectory + "/" + dummy;
        fileNames.push_back(dummy);
    }

    // read energy bins from config
    TH1D* TOFHisto = new TH1D("","",config.plot.TOF_BINS, config.plot.TOF_LOWER_BOUND, config.plot.TOF_UPPER_BOUND);
    TH1D* energyHisto = timeBinsToRKEBins(TOFHisto, "");

    vector<double> energyBins;
    
    int numberOfBins = energyHisto->GetNbinsX();

    for(int i=1; i<=numberOfBins; i++)
    {
        energyBins.push_back(energyHisto->GetBinLowEdge(i));
    }

    energyBins.push_back(
            energyHisto->GetBinLowEdge(numberOfBins)
           +energyHisto->GetBinWidth(numberOfBins));

    // recreate output file
    TFile* outFile = new TFile(litOutputName.c_str(),"RECREATE");

    vector<DataSet> allData;

    for(string s : fileNames)
    {
        cout << "Creating plot for " << s << endl;
        allData.push_back(DataSet(s, energyBins));
        outFile->cd();
        allData.back().getPlot()->Write();
    }

    outFile->Close();

    return 0;
}

int correctForBlank(CrossSection& rawCS, CSPrereqs& targetData, string expName)
{
    string blankDataLocation = "../" + expName + "/targetData/blankCorrection.txt";
    ifstream blankDataFile(blankDataLocation.c_str());
    if(!blankDataFile.is_open())
    {
        std::cout << "Attempted to read blank data, but failed to find data in " << blankDataLocation << std::endl;
        exit(1);
    }

    vector<Target> blankComposition; // for holding multiple elements that comprise the blank

    string str;

    while(getline(blankDataFile,str))
    {
        // parse into tokens

        vector<string> tokens;
        istringstream iss(str);
        copy(istream_iterator<string>(iss),
                istream_iterator<string>(),
                back_inserter(tokens));

        if(tokens[0]=="Name:")
        {
            blankComposition.push_back(Target());
            blankComposition.back().setName(tokens[1]);
        }

        else if(tokens[0]=="Length:")
        {
            blankComposition.back().setLength(atof(tokens[1].c_str()));
        }

        else if(tokens[0]=="Diameter:")
        {
            blankComposition.back().setDiameter(atof(tokens[1].c_str()));
        }

        else if(tokens[0]=="Mass:")
        {
            blankComposition.back().setMass(atof(tokens[1].c_str()));
        }

        else if(tokens[0]=="Molar")
        {
            blankComposition.back().setMolarMass(atof(tokens[2].c_str()));
        }

        else
        {
            cerr << "Error - couldn't parse a line in a targetData text file" << endl;
            exit(1);
        }
    }

    for(Target t : blankComposition)
    {
        double blankNumberDensity = (t.getMass()/t.getMolarMass())*AVOGADROS_NUMBER/(t.getLength()*M_PI*pow((t.getDiameter()/2),2));

        // calculate number of atoms in this target
        double numberTargetAtoms =
            (targetData.target.getMass()/targetData.target.getMolarMass())*AVOGADROS_NUMBER;

        // calculate number density (atoms/cm^3) in target
    double targetNumberDensity =
        numberTargetAtoms/(targetData.target.getLength()*pow(targetData.target.getDiameter()/2,2)*M_PI); // area of cylinder end

        double ratioNumberDensities = blankNumberDensity/targetNumberDensity; // ratio of atoms of this element in blank, compared to target
        ratioNumberDensities *= -1; // the correction should be additive, not subtractive

        string graphFileName = "../" + expName + "/literatureData/literatureData.root";
        string blankCSGraphName = t.getName() + "(n,tot)";

        string name = rawCS.name + "blankCorrected";
        rawCS = subtractCS(rawCS, graphFileName, blankCSGraphName, ratioNumberDensities, 1, "blankCorrected");
    }

    return 0;
}

/*int produceTotalCSPlots(string dataLocation, vector<CrossSection>& crossSections)
{
    //  Create total cross section plots

    string outFileName = dataLocation + "/total.root";
    TFile* outFile = new TFile(outFileName.c_str(), "UPDATE");

    outFile->Close();

    return 0;
}*/

/*int produceRelativeCSPlots(string dataLocation, string string expName,
 * vector<CrossSection>& crossSections)
  {

    //  Create relative cross section plots

    string relativeFileName = dataLocation + "/relative.root";
    TFile* relativeFile = new TFile(relativeFileName.c_str(), "UPDATE");

    // read which relative cross section plots to make from the experimental
    // directory
    vector<pair<string,string>> relativeTargets = getRelativePlotNames(expName,"relativePlots.txt");

    for(pair<string,string> p : relativeTargets)
    {
        CrossSection larger;
        CrossSection smaller;

        for(auto& cs : crossSections)
        {
            if(cs.name==p.first)
            {
                larger = cs;
            }
            else if(cs.name==p.second)
            {
                smaller = cs;
            }
        }

        if(larger.name!="" && smaller.name!="")
        {
            // found cross section plots for both individual targets
            cout << "Producing relative cross section plot of " << larger.name << " and " << smaller.name << endl;

            CrossSection relative = calculateRelativeDiff(larger,smaller);
            string relativeTitle = "#frac{#sigma_{" + larger.name +
                "}-#sigma_{" + smaller.name +
                "}}{#sigma_{" + larger.name +
                "}+#sigma_{" + smaller.name + "}}";
            string relativeName = "relative" + larger.name + smaller.name;
            relative.createGraph(relativeName.c_str(), relativeTitle.c_str());
            cout << "Relative plot " << relative.getDataSet().getReference() <<
                " RMS error: " << relative.calculateRMSError() << endl;
        }

        else
        {
            cout << "Failed to find cross section plot for either " << larger.name << " or " << smaller.name << endl;
        }
    }

    relativeFile->Close();

    return 0;
}*/

CrossSection calculateRelativeDiff(CrossSection a, CrossSection b)
{
    CrossSection relative; 

    DataSet aData = a.getDataSet();
    DataSet bData = b.getDataSet();

    if(aData.getNumberOfPoints()!=bData.getNumberOfPoints())
    {
        cerr << "Error: can't calculate relative cross section from "
             << "data sets of different sizes. Returning empty cross section."
             << endl;
        return relative;
    }

    DataSet relativeDataSet;

    // for each point, calculate the relative cross section, including error
    for(int i=0; i<aData.getNumberOfPoints(); i++)
    {
        DataPoint aPoint = aData.getPoint(i);
        DataPoint bPoint = bData.getPoint(i);

        double aXValue = aPoint.getXValue();
        double bXValue = bPoint.getXValue();
        
        if(aXValue != bXValue)
        {
            cerr << "Error: can't calculate relative cross section from "
                 << "data points with different x values. Returning cross section."
                 << endl;
            return relative;
        }

        double aYValue = aPoint.getYValue();
        double bYValue = bPoint.getYValue();

        double yValue = (aYValue-bYValue)/(aYValue+bYValue);

        // calculate cross section error
        double aError = getPartialError(aPoint, bPoint, a.getArealDensity());
        double bError = getPartialError(bPoint, aPoint, b.getArealDensity());
        double totalError = pow(pow(aError,2)+pow(bError,2),0.5);

        relativeDataSet.addPoint(
                DataPoint(aXValue, aPoint.getXError(), 100*yValue, 100*totalError, /* convert to % */
                          aPoint.getBlankMonitorCounts()+bPoint.getBlankMonitorCounts(),
                          aPoint.getTargetMonitorCounts()+bPoint.getTargetMonitorCounts(),
                          aPoint.getBlankDetCounts()+bPoint.getBlankDetCounts(),
                          aPoint.getTargetDetCounts()+bPoint.getTargetDetCounts()));
    }

    relative.addDataSet(relativeDataSet);
    return relative;
}

CrossSection subtractCS(string rawCSFileName, string rawCSGraphName,
        string subtrahendFileName, string subtrahendGraphName,
        double factor, double divisor, string name)
{
    // get rawCS graph
    TFile* rawCSFile = new TFile(rawCSFileName.c_str(),"UPDATE");
    TGraphErrors* rawCSGraph = (TGraphErrors*)rawCSFile->Get(rawCSGraphName.c_str());
    if(!rawCSGraph)
    {
        cerr << "Error: failed to find " << rawCSGraphName << " in " << rawCSFileName << endl;
        exit(1);
    }

    // get subtrahend graph
    TFile* subtrahendFile = new TFile(subtrahendFileName.c_str(),"READ");
    TGraphErrors* subtrahendGraph = (TGraphErrors*)subtrahendFile->Get(subtrahendGraphName.c_str());
    if(!subtrahendGraph)
    {
        cerr << "Error: failed to find " << subtrahendGraphName << " in " << subtrahendFileName << endl;
        exit(1);
    }

    DataSet subtrahendData = DataSet();
    DataSet rawCSData = DataSet(rawCSGraph, rawCSGraphName);

    // for each y-value of the raw CS graph, read the y-value of the subtrahend
    // and the y-error
    for(int i=0; i<rawCSData.getNumberOfPoints(); i++)
    {
        subtrahendData.addPoint(
                DataPoint(rawCSData.getPoint(i).getXValue(),
                    rawCSData.getPoint(i).getXError(),
                    subtrahendGraph->Eval(rawCSData.getPoint(i).getXValue()),
                    subtrahendGraph->GetErrorY(rawCSData.getPoint(i).getXValue()))); 
    }

    // perform the subtraction
    CrossSection differenceCS = CrossSection();
    differenceCS.addDataSet((rawCSData-subtrahendData*factor)/divisor);

    // create graph of difference
    rawCSFile->cd();
    differenceCS.createGraph(name, name);

    return differenceCS;
}

CrossSection subtractCS(CrossSection rawCS, string subtrahendFileName,
        string subtrahendGraphName, double factor, double divisor, string name)
{
    // get subtrahend graph
    TFile* subtrahendFile = new TFile(subtrahendFileName.c_str(),"READ");
    TGraphErrors* subtrahendGraph = (TGraphErrors*)subtrahendFile->Get(subtrahendGraphName.c_str());
    if(!subtrahendGraph)
    {
        cerr << "Error: failed to find " << subtrahendGraphName << " in " << subtrahendFileName << endl;
        subtrahendFile->Close();
        exit(1);
    }

    DataSet subtrahendData = DataSet();
    DataSet rawCSData = rawCS.getDataSet();

    // for each y-value of the raw CS graph, read the y-value of the subtrahend
    // and the y-error
    for(int i=0; i<rawCSData.getNumberOfPoints(); i++)
    {
        subtrahendData.addPoint(
                DataPoint(rawCSData.getPoint(i).getXValue(),
                    rawCSData.getPoint(i).getXError(),
                    subtrahendGraph->Eval(rawCSData.getPoint(i).getXValue()),
                    subtrahendGraph->GetErrorY(rawCSData.getPoint(i).getXValue()))); 
    }

    // perform the subtraction
    CrossSection differenceCS = CrossSection();
    differenceCS.addDataSet((rawCSData-subtrahendData*factor)/divisor);

    differenceCS.name = rawCS.name;

    subtrahendFile->Close();
    return differenceCS;
}

CrossSection shiftCS(string rawCSFileName, string rawCSGraphName,
        double shift, string outputFileName, string outputGraphName)
{
    // get rawCS graph
    TFile* rawCSFile = new TFile(rawCSFileName.c_str(),"READ");
    TGraphErrors* rawCSGraph = (TGraphErrors*)rawCSFile->Get(rawCSGraphName.c_str());
    if(!rawCSGraph)
    {
        cerr << "Error: failed to find " << rawCSGraphName << " in " << rawCSFileName << endl;
        exit(1);
    }

    DataSet rawCSData = DataSet(rawCSGraph, rawCSGraphName);

    // perform the shift
    CrossSection shiftedCS = CrossSection();
    shiftedCS.addDataSet(rawCSData+shift);

    // write output graph
    TFile* outputFile = new TFile(outputFileName.c_str(),"UPDATE");
    shiftedCS.createGraph(outputGraphName, outputGraphName);

    rawCSFile->Close();
    outputFile->Close();

    return shiftedCS;
}

CrossSection multiplyCS(string rawCSFileName, string rawCSGraphName,
        double factor, string name)
{
    // get rawCS graph
    TFile* rawCSFile = new TFile(rawCSFileName.c_str(),"UPDATE");
    TGraphErrors* rawCSGraph = (TGraphErrors*)rawCSFile->Get(rawCSGraphName.c_str());
    if(!rawCSGraph)
    {
        cerr << "Error: failed to find " << rawCSGraphName << " in " << rawCSFileName << endl;
        exit(1);
    }

    DataSet rawDS = DataSet(rawCSGraph,rawCSGraphName);

    // perform the multiplication
    DataSet productDS = rawDS*factor;

    CrossSection productCS = CrossSection();
    productCS.addDataSet(productDS);

    // create graph of difference
    rawCSFile->cd();
    productCS.createGraph(name, name);

    return productCS;
}

double getPartialError(DataPoint aPoint, DataPoint bPoint, double aArealDensity)
{
    double aYValue = aPoint.getYValue();
    double bYValue = bPoint.getYValue();

    double aBMon = aPoint.getBlankMonitorCounts();
    double aTMon = aPoint.getTargetMonitorCounts();
    double aBDet = aPoint.getBlankDetCounts();
    double aTDet = aPoint.getTargetDetCounts();

    if(aBMon<=0 || aTMon<=0 || aBDet<=0 || aTDet<=0)
    {
        //cerr << "Error: could not calculate partial error with non-finite count ratio. " << endl;
        return 0;
    }

    // first calculate ratio of monitor/detector counts for each data set
    double aCountRatio = (aBDet/aTDet)/(aBMon/aTMon);

    // calculate error of aCountRatio
    double errorACountRatio = aCountRatio * pow(1/aBDet+1/aTDet+1/aBMon+1/aTMon,0.5);

    // calculate error of data point from data set a partial derivative
    double leftExpression = errorACountRatio/(aCountRatio*aArealDensity*pow(10,-24)*(aYValue+bYValue));
    double rightExpression = 1-(aYValue/(aYValue+bYValue));

    return leftExpression*rightExpression;
}

void produceRunningRMS(DataSet firstDS, DataSet secondDS, string name)
{
    DataSet rms;

    for(int i=0; i<firstDS.getNumberOfPoints(); i++)
    {
        double relDiff = 
            (firstDS.getPoint(i).getYValue() -
            secondDS.getPoint(i).getYValue())/
            (firstDS.getPoint(i).getYValue() +
             secondDS.getPoint(i).getYValue());

        if(i==0)
        {
            rms.addPoint(
                    DataPoint(firstDS.getPoint(i).getXValue(),
                        0,
                        pow(pow(relDiff,2),0.5),
                        0)
                    );
        }

        else
        {
            rms.addPoint(DataPoint(firstDS.getPoint(i).getXValue(),
                        0,
                        pow((pow(relDiff,2)+pow(rms.getPoint(i-1).getYValue(),2)*i)/(i+1),0.5),
                        0));
        }
    }

    cout << "Total RMS at " << rms.getPoint(rms.getNumberOfPoints()-1).getXValue()
        << " MeV = " << rms.getPoint(rms.getNumberOfPoints()-1).getYValue() << endl;

    CrossSection rmsPlot = CrossSection();
    rmsPlot.addDataSet(rms);

    string n = name + "rms";
    rmsPlot.createGraph(n.c_str(), n.c_str());
}

CrossSection relativeCS(string firstCSFileName, string firstCSGraphName,
        string secondCSFileName, string secondCSGraphName,
        string outputFileName, string name)
{
    // get firstCS graph
    TFile* firstCSFile = new TFile(firstCSFileName.c_str(),"READ");
    TGraphErrors* firstCSGraph = (TGraphErrors*)firstCSFile->Get(firstCSGraphName.c_str());
    if(!firstCSGraph)
    {
        cerr << "Error: failed to find " << firstCSGraphName << " in " << firstCSFileName << endl;
        exit(1);
    }

    // get second graph
    TFile* secondCSFile = new TFile(secondCSFileName.c_str(),"READ");
    TGraphErrors* secondCSGraph = (TGraphErrors*)secondCSFile->Get(secondCSGraphName.c_str());
    if(!secondCSGraph)
    {
        cerr << "Error: failed to find " << secondCSGraphName << " in " << secondCSFileName << endl;
        exit(1);
    }

    DataSet firstCSData = DataSet();
    DataSet secondCSData = DataSet();

    // find maximum value of second CS dataset
    DataSet firstCSDataRaw = DataSet(firstCSGraph, firstCSGraphName);
    DataSet secondCSDataRaw = DataSet(secondCSGraph, secondCSGraphName);
    double maxEnergyValue = secondCSDataRaw.getPoint(secondCSDataRaw.getNumberOfPoints()-1).getXValue();

    // for each y-value of the first CS graph, read the y-value of the second
    // and the y-error
    for(int i=0; i<firstCSDataRaw.getNumberOfPoints(); i++)
    {
        if(firstCSDataRaw.getPoint(i).getXValue()>maxEnergyValue)
        {
            break;
        }

        firstCSData.addPoint(firstCSDataRaw.getPoint(i));

        secondCSData.addPoint(
                DataPoint(firstCSDataRaw.getPoint(i).getXValue(),
                    firstCSDataRaw.getPoint(i).getXError(),
                    secondCSGraph->Eval(firstCSDataRaw.getPoint(i).getXValue()),
                    0
                    /*secondCSGraph->GetErrorY(firstCSDataRaw.getPoint(i).getXValue())*/)); 
    }

    // perform the division
    CrossSection relCS = CrossSection();
    relCS.addDataSet(firstCSData/secondCSData);

    // create graph of relative difference
    TFile* outputFile = new TFile(outputFileName.c_str(), "UPDATE");
    relCS.createGraph(name, name);

    // create running RMS plot
    CrossSection firstCS = CrossSection();
    firstCS.addDataSet(firstCSData);

    CrossSection secondCS = CrossSection();
    secondCS.addDataSet(secondCSData);

    produceRunningRMS(firstCS.getDataSet(), secondCS.getDataSet(), name);

    outputFile->Close();

    firstCSFile->Close();
    secondCSFile->Close();

    return relCS;
}

CrossSection relativeDiffCS(string firstCSFileName, string firstCSGraphName,
        string secondCSFileName, string secondCSGraphName,
        string outputFileName, string name)
{
    // get firstCS graph
    TFile* firstCSFile = new TFile(firstCSFileName.c_str(),"READ");
    TGraphErrors* firstCSGraph = (TGraphErrors*)firstCSFile->Get(firstCSGraphName.c_str());
    if(!firstCSGraph)
    {
        cerr << "Error: failed to find " << firstCSGraphName << " in " << firstCSFileName << endl;
        exit(1);
    }

    // get secondCS graph
    TFile* secondCSFile = new TFile(secondCSFileName.c_str(),"READ");
    TGraphErrors* secondCSGraph = (TGraphErrors*)secondCSFile->Get(secondCSGraphName.c_str());
    if(!secondCSGraph)
    {
        cerr << "Error: failed to find " << secondCSGraphName << " in " << secondCSFileName << endl;
        exit(1);
    }

    DataSet firstCSData = DataSet();
    DataSet secondCSData = DataSet();

    DataSet firstCSDataRaw = DataSet(firstCSGraph, firstCSGraphName);
    DataSet secondCSDataRaw = DataSet(secondCSGraph, secondCSGraphName);

    // find maximum value of datasets
    double maxEnergyValue = min(
            firstCSDataRaw.getPoint(firstCSDataRaw.getNumberOfPoints()-1).getXValue(),
            secondCSDataRaw.getPoint(secondCSDataRaw.getNumberOfPoints()-1).getXValue());

    // find minimum value of datasets
    double minEnergyValue = max(
            firstCSDataRaw.getPoint(0).getXValue(),
            secondCSDataRaw.getPoint(0).getXValue());

    // for each y-value of the first CS graph, read the y-value of the second
    // and the y-error
    for(int i=0; i<firstCSDataRaw.getNumberOfPoints(); i++)
    {
        if(firstCSDataRaw.getPoint(i).getXValue()<minEnergyValue)
        {
            continue;
        }

        if(firstCSDataRaw.getPoint(i).getXValue()>maxEnergyValue)
        {
            break;
        }

        firstCSData.addPoint(firstCSDataRaw.getPoint(i));

        secondCSData.addPoint(
                DataPoint(firstCSDataRaw.getPoint(i).getXValue(),
                    firstCSDataRaw.getPoint(i).getXError(),
                    secondCSGraph->Eval(firstCSDataRaw.getPoint(i).getXValue()),
                    0
                    /*secondCSGraph->GetErrorY(firstCSDataRaw.getPoint(i).getXValue())*/)); 
    }

    // perform the difference
    CrossSection differenceCS = CrossSection();
    differenceCS.addDataSet(firstCSData-secondCSData);

    // perform the sum
    CrossSection sumCS = CrossSection();
    sumCS.addDataSet(firstCSData+secondCSData);

    // perform the division
    CrossSection relDiffCS = CrossSection();
    relDiffCS.addDataSet(differenceCS.getDataSet()/sumCS.getDataSet());

    // create graph of relative difference
    TFile* outputFile = new TFile(outputFileName.c_str(), "UPDATE");
    relDiffCS.createGraph(name, name);

    // create running RMS plot
    CrossSection firstCS = CrossSection();
    firstCS.addDataSet(firstCSData);

    CrossSection secondCS = CrossSection();
    secondCS.addDataSet(secondCSData);

    produceRunningRMS(firstCS.getDataSet(), secondCS.getDataSet(), name);

    outputFile->Close();

    firstCSFile->Close();
    secondCSFile->Close();

    return relDiffCS;
}

void applyCSCorrectionFactor(string CSCorrectionFileName, string CSCorrectionGraphName, string CSToBeCorrectedFileName, string CSToBeCorrectedGraphName, string outputFileName, string outputGraphName)
{
    // get CS correction graph
    TFile* CSCorrectionFile = new TFile(CSCorrectionFileName.c_str(),"READ");
    TGraphErrors* CSCorrectionGraph = (TGraphErrors*)CSCorrectionFile->Get(CSCorrectionGraphName.c_str());
    if(!CSCorrectionGraph)
    {
        cerr << "Error: failed to find " << CSCorrectionGraphName << " in " << CSCorrectionFileName << endl;
        exit(1);
    }

    // get CSToBeCorrected graph
    TFile* CSToBeCorrectedFile = new TFile(CSToBeCorrectedFileName.c_str(),"READ");
    TGraphErrors* CSToBeCorrectedGraph = (TGraphErrors*)CSToBeCorrectedFile->Get(CSToBeCorrectedGraphName.c_str());
    if(!CSToBeCorrectedGraph)
    {
        cerr << "Error: failed to find " << CSToBeCorrectedGraphName << " in " << CSToBeCorrectedFileName << endl;
        exit(1);
    }

    TFile* outputFile = new TFile(outputFileName.c_str(),"UPDATE");

    DataSet CSToBeCorrectedData = DataSet(CSToBeCorrectedGraph, CSToBeCorrectedGraphName);
    DataSet CSCorrectionData = DataSet();

    // for each y-value of the CSToBeCorrected graph, read the y-value of the CSCorrection graph
    for(int i=0; i<CSToBeCorrectedData.getNumberOfPoints(); i++)
    {
        CSCorrectionData.addPoint(
                DataPoint(CSToBeCorrectedData.getPoint(i).getXValue(),
                    CSToBeCorrectedData.getPoint(i).getXError(),
                    CSCorrectionGraph->Eval(CSToBeCorrectedData.getPoint(i).getXValue()),
                    CSCorrectionGraph->GetErrorY(CSToBeCorrectedData.getPoint(i).getXValue())));
    }

    // perform the correction
    CrossSection correctedCS = CrossSection();
    correctedCS.addDataSet(correctCSUsingControl(CSToBeCorrectedData,CSCorrectionData));

    // create graph of correctedCS
    outputFile->cd();
    correctedCS.createGraph(outputGraphName, outputGraphName);

    outputFile->Close();
}

void scaledownCS(string CSToBeCorrectedFileName, string CSToBeCorrectedGraphName, int scaledown, string outputFileName, string outputGraphName)
{
    // get CSToBeCorrected graph
    TFile* CSToBeCorrectedFile = new TFile(CSToBeCorrectedFileName.c_str(),"READ");
    TGraphErrors* CSToBeCorrectedGraph = (TGraphErrors*)CSToBeCorrectedFile->Get(CSToBeCorrectedGraphName.c_str());
    if(!CSToBeCorrectedGraph)
    {
        cerr << "Error: failed to find " << CSToBeCorrectedGraphName << " in " << CSToBeCorrectedFileName << endl;
        exit(1);
    }

    TFile* outputFile = new TFile(outputFileName.c_str(),"UPDATE");

    DataSet CSToBeCorrectedData = DataSet(CSToBeCorrectedGraph, CSToBeCorrectedGraphName);
    DataSet CorrectedCSData = DataSet();

    double rebinnedXValue = 0;
    double rebinnedYValue = 0;
    double rebinnedYError = 0;
    int numberOfPoints = 0;

    // create downscaled CSToBeCorrected graph
    for(int i=0; i<CSToBeCorrectedData.getNumberOfPoints(); i++)
    {
        rebinnedXValue += CSToBeCorrectedData.getPoint(i).getXValue();
        rebinnedYValue += CSToBeCorrectedData.getPoint(i).getYValue();
        rebinnedYError += CSToBeCorrectedData.getPoint(i).getYError();
        numberOfPoints++;

        if(i%scaledown==scaledown-1)
        {
            rebinnedYError /= numberOfPoints;
            CorrectedCSData.addPoint(
                DataPoint(rebinnedXValue/numberOfPoints,
                    CSToBeCorrectedData.getPoint(i).getXError(),
                    rebinnedYValue/numberOfPoints,
                    rebinnedYError));

            rebinnedXValue = 0;
            rebinnedYValue = 0;
            numberOfPoints = 0;
        }
    }

    // perform the correction
    CrossSection correctedCS = CrossSection();
    correctedCS.addDataSet(CorrectedCSData);

    // create graph of correctedCS
    outputFile->cd();
    correctedCS.createGraph(outputGraphName, outputGraphName);

    outputFile->Close();
}


