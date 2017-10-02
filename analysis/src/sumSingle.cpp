#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TLatex.h"

#include "../include/physicalConstants.h"
#include "../include/target.h"
#include "../include/dataSet.h"
#include "../include/dataPoint.h"
#include "../include/CSPrereqs.h"
#include "../include/crossSection.h"
#include "../include/experiment.h"

using namespace std;

struct Plots
{
    vector<TGraphErrors*> CSGraphs;
    vector<TGraphErrors*> CSGraphsWaveform;

    vector<TGraphErrors*> relativeCSGraphs;
    vector<TGraphErrors*> relativeCSGraphsWaveform;
} plots;

// Extract each point from a graph and store their positions in two vectors,
// xValues and yValues
void extractGraphData(
        TGraphErrors* graph,
        vector<double>* xValues,
        vector<double>* xError,
        vector<double>* yValues,
        vector<double>* yError)
{
    int numPoints = graph->GetN();
    
    xValues->resize(numPoints);
    yValues->resize(numPoints);
    xError->resize(numPoints);
    yError->resize(numPoints);

    for(int k=0; k<numPoints; k++)
    {
        graph->GetPoint(k,xValues->at(k),yValues->at(k));
        xError->at(k) = graph->GetErrorX(k);
        yError->at(k) = graph->GetErrorY(k);
    }
}

// Extract each point from a graph and store their positions in a DataSet
void extractGraphData(
        TGraphErrors* graph,
        DataSet& dataSet)
{
    int numPoints = graph->GetN();

    double xValue;
    double xError;
    double yValue;
    double yError;

    for(int k=0; k<numPoints; k++)
    {
        graph->GetPoint(k,xValue,yValue);
        xError = graph->GetErrorX(k);
        yError = graph->GetErrorY(k);

        DataPoint dataPoint(xValue, xError, yValue, yError);
        dataSet.addPoint(dataPoint);
    }
}


// Calculate the root-mean-squared difference between two vectors
double calculateRMS(vector<double> graph1Data, vector<double> graph2Data)
{
    double rms = 0;
    for(int i=0; (size_t)i<graph1Data.size(); i++)
    {
        rms += pow(graph1Data.at(i)-graph2Data.at(i),2);
    }
    rms /= graph1Data.size();
    rms = pow(rms,0.5);
    return rms;
}

// Scale a dataset's y-values by the ratio of two other datasets
DataSet scale(DataSet setToScale, DataSet expReference, DataSet litReference)
{
    DataSet scaleFactor;

    int n = setToScale.getNumberOfPoints();
    for(int i=0; i<n; i++)
    {
        DataPoint p = litReference.getPoint(i)/expReference.getPoint(i);
        scaleFactor.addPoint(p);
    }

    DataSet scaledSet = scaleFactor*setToScale;

    return scaledSet;
}

CrossSection correctForBlank(CrossSection rawCS, double targetNumberDensity, string expName, string graphFileName)
{
    string blankDataLocation = "../" + expName + "/targetData/" + "blank.txt";
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

    CrossSection correctedCS = rawCS;

    for(Target t : blankComposition)
    {
        double factor = (t.getMass()/t.getMolarMass())*AVOGADROS_NUMBER/(t.getLength()*M_PI*pow((t.getDiameter()/2),2));

        if(targetNumberDensity==0)
        {
            continue;
        }

        factor /= targetNumberDensity; // ratio of atoms of this element in blank, compared to target
        factor *= -1; // the correction should be additive, not subtractive

        string graphFileLocation = "../" + expName + "/literatureData/" + graphFileName;
        string graphFileName = t.getName() + "(n,tot)";
        correctedCS.subtractCS(graphFileLocation, graphFileName, factor);
    }

    return correctedCS;
}

CrossSection calculateCS(string CSFileName, CSPrereqs& targetData, CSPrereqs& blankData, string expName)
{
    // make sure the monitor recorded counts for both the blank and the target
    // of interest so we can normalize the flux between them
    if(targetData.monitorCounts == 0 || blankData.monitorCounts == 0)
    {
        cerr << "Error - didn't find any monitor counts for target while trying to calculate cross sections. Returning empty cross section..." << endl;
        return CrossSection();
    }

    // define variables to hold cross section information
    CrossSection crossSection;
    double energyValue;
    double energyError;
    double crossSectionValue;
    double crossSectionError;

    int numberOfBins = targetData.energyHisto->GetNbinsX();

    // calculate the ratio of target/blank monitor counts (normalize flux)
    double tMon = targetData.monitorCounts;
    double bMon = blankData.monitorCounts;
    double monitorRatio = tMon/bMon;

    // calculate number of atoms in this target
    long double numberOfAtoms =
        (targetData.target.getMass()/targetData.target.getMolarMass())*AVOGADROS_NUMBER;

    // calculate areal density (atoms/cm^2) in target
    double arealDensity =
        numberOfAtoms/(pow(targetData.target.getDiameter()/2,2)*M_PI); // area of cylinder end

    // calculate volume density (atoms/cm^3) in target
    double volumeDensity =
        arealDensity/targetData.target.getLength();

    // save this areal density for later error propagation
    crossSection.setArealDensity(arealDensity);

    // loop through each bin in the energy histo, calculating a cross section
    // for each bin
    for(int i=1; i<=numberOfBins-1; i++) // skip the overflow and underflow bins
    {
        // read data from detector histograms for target and blank
        TH1D* bCounts = blankData.energyHisto;
        TH1D* tCounts = targetData.energyHisto;

        energyValue = tCounts->GetBinCenter(i);
        energyError = tCounts->GetBinWidth(i)/2;

        long bDet = bCounts->GetBinContent(i);
        long tDet = tCounts->GetBinContent(i);

        // if any essential values are 0, return an empty DataPoint
        if(bMon <=0 || tMon <=0 || bDet <=0 || tDet <=0)
        {
            crossSectionValue = 0;
            crossSectionError = 0;
            crossSection.addDataPoint(
                DataPoint(energyValue, energyError, crossSectionValue, crossSectionError,
                          bMon, tMon, bDet, tDet));
            continue;
        }

        // calculate the ratio of target/blank counts in the detector
        double detectorRatio = (double)tDet/bDet;

        crossSectionValue =
            -log(detectorRatio/monitorRatio)/arealDensity; // in cm^2

        crossSectionValue *= pow(10,24); // in barns 
            
        // calculate the statistical error
        crossSectionError =
            pow(1/tDet+1/bDet+1/bMon+1/tMon,0.5)/arealDensity; // in cm^2

        crossSectionError *= pow(10,24); // in barns

        crossSection.addDataPoint(
                DataPoint(energyValue, energyError, crossSectionValue, crossSectionError,
                    bMon, tMon, bDet, tDet));
    }

    CrossSection corrected = crossSection;
    //CrossSection corrected = correctForBlank(crossSection, volumeDensity, expName, "literatureData.root");

    TFile* CSFile = new TFile(CSFileName.c_str(),"UPDATE");
    string name = targetData.target.getName();
    crossSection.createCSGraph(name, name);
    name += "target0CorrectedEnergy";
    corrected.createCSGraph(name, name);

    CSFile->Write();
    CSFile->Close();

    return corrected;
}

int main(int, char* argv[])
{
    string dataLocation = argv[1]; // name of directory where analysis is stored
                                   // (omit trailing slash)
    string expName = argv[2];      // experiment directory where runs to-be-sorted
                                   // are listed
    string ROOTFileName = argv[3]; // name of ROOT files that contain data
                                   // used to calculate cross sections

    string runNumber = argv[4];

    int subRun = stoi(argv[5]);

    // Create a CSPrereqs for each target to hold data from all the runs
    vector<CSPrereqs> allData;

    vector<string> targetNames = readExperimentConfig(expName,"targetNames");
    for(auto targetName : targetNames)
    {
        string targetDataLocation = "../" + expName + "/targetData/" + targetName + ".txt";
        allData.push_back(CSPrereqs(targetDataLocation));
    }

    cout << endl;

    stringstream subRunFormatted;
    subRunFormatted << setfill('0') << setw(4) << subRun;

    // open subrun
    stringstream inFileName;
    inFileName << dataLocation << "/" << runNumber << "/"
        << subRunFormatted.str() << "/" << ROOTFileName << ".root";
    ifstream f(inFileName.str());
    if(!f.good())
    {
        cerr << "Failed to open subrun " << inFileName.str() << endl;
        exit(1);
    }

    f.close();

    TFile* inFile = new TFile(inFileName.str().c_str(),"READ");

    // get target order for this run
    vector<string> targetOrder = getTargetOrder(expName, stoi(runNumber));

    cout << "Adding " << runNumber << ", subrun " << subRun << "\r";
    fflush(stdout);

    // Loop through all target positions in this subrun
    for(int j=0; (size_t)j<targetOrder.size(); j++)
    {
        // pull data needed for CS calculation from subrun 
        string targetDataLocation = "../" + expName + "/targetData/" + targetOrder[j] + ".txt";
        CSPrereqs subRunData(targetDataLocation);

        subRunData.readData(inFile, "summedDet", j);

        // find the correct CSPrereqs to add this target's data to
        for(int k=0; (size_t)k<allData.size(); k++)
        {
            if(allData[k].target.getName() == subRunData.target.getName())
            {
                // add subrun data to total
                if(!allData[k].energyHisto)
                {
                    // this is the first subrun to be added
                    allData[k].monitorCounts = subRunData.monitorCounts;
                    allData[k].energyHisto = (TH1D*)subRunData.energyHisto->Clone();
                    // prevent the cloned histogram from being closed
                    // when the subrun is closed
                    allData[k].energyHisto->SetDirectory(0);
                }

                else
                {
                    allData[k] = allData[k] + subRunData;
                }

                break;
            }

            if((size_t)k+1==allData.size())
            {
                cerr << "Failed to find a CSPrereqs to add this subrun to." << endl;
                continue;
            }
        }
    }

    // Close the sub-run input files
    inFile->Close();

    string outFileName = dataLocation + "/total.root";

    vector<CrossSection> crossSections;

    cout << endl << "Total statistics over all runs: " << endl << endl;

    for(CSPrereqs p : allData)
    {
        long totalCounts = 0;
        for(int i=0; i<p.energyHisto->GetNbinsX(); i++)
        {
            totalCounts += p.energyHisto->GetBinContent(i);
        }

        cout << p.target.getName() << ": total events in energy histo = "
             << totalCounts << ", total monitor events = "
             << p.monitorCounts << endl;
        crossSections.push_back(calculateCS(outFileName, p, allData[0], expName));
        cout << "Target " << crossSections.back().getDataSet().getReference() <<
                " RMS error: " << crossSections.back().calculateRMSError() << endl << endl;
    }

    /**************************************************************************
    Create relative cross section plots
    **************************************************************************/

    string relativeFileName = dataLocation + "/relative.root";
    TFile* relativeFile = new TFile(relativeFileName.c_str(), "RECREATE");

    // read which relative cross section plots to make from the experimental
    // directory
    vector<pair<string,string>> relativeTargets = getRelativePlotNames(expName,"relativePlots.txt");

    for(pair<string,string> p : relativeTargets)
    {
        int largerTarget = -1;
        int smallerTarget = -1;
        for(int i=0; (size_t)i<targetNames.size(); i++)
        {
            if(targetNames[i]==p.first)
            {
                largerTarget = i;
            }
            else if(targetNames[i]==p.second)
            {
                smallerTarget = i;
            }
        }

        if(largerTarget>=0 && smallerTarget>=0)
        {
            // found cross section plots for both individual targets,
            // so use them to make a relative cross section plot
            cout << "Producing relative cross section plot of " << targetNames[largerTarget] << " and " << targetNames[smallerTarget] << endl;

            //CrossSection sum = crossSections[largerTarget]+crossSections[smallerTarget];
            //cout << "sum plot " << sum.getDataSet().getReference() <<
            //        " RMS error: " << sum.calculateRMSError() << endl << endl;

            //CrossSection difference = crossSections[largerTarget]-crossSections[smallerTarget];
            //cout << "difference plot " << difference.getDataSet().getReference() <<
            //        " RMS error: " << difference.calculateRMSError() << endl << endl;

            CrossSection relative = calculateRelative(crossSections[largerTarget],crossSections[smallerTarget]);
            string relativeName = "#frac{#sigma_{" + allData[largerTarget].target.getName() +
                                      "}-#sigma_{" + allData[smallerTarget].target.getName() + 
                                     "}}{#sigma_{" + allData[largerTarget].target.getName() +
                                      "}+#sigma_{" + allData[smallerTarget].target.getName() + "}}";
            relative.createCSGraph(relativeName.c_str(), relativeName.c_str());
            cout << "Relative plot " << relative.getDataSet().getReference() <<
                    " RMS error: " << relative.calculateRMSError() << endl;
        }

        else
        {
            cout << "Failed to find cross section plot for either " << targetNames[largerTarget] << " or " << targetNames[smallerTarget] << endl;
        }
    }

    relativeFile->Close();

    /**************************************************************************
    Create subtracted cross section plots
    **************************************************************************/
    /*
    string subtractedFileName = dataLocation + "/subtracted.root";
    TFile* subtractedFile = new TFile(subtractedFileName.c_str(), "RECREATE");

    // read which subtracted cross section plots to make from the experimental
    // directory
    vector<pair<string,string>> subtractedTargets = getSubtractedPlotNames(expName,"subtractedPlots.txt");

    for(pair<string,string> p : subtractedTargets)
    {
        int largerTarget = -1;
        int smallerTarget = -1;
        for(int i=0; (size_t)i<targetNames.size(); i++)
        {
            if(targetNames[i]==p.first)
            {
                largerTarget = i;
            }
            else if(targetNames[i]==p.second)
            {
                smallerTarget = i;
            }
        }

        if(largerTarget>=0)
        {
            // found experimental cross section plot for this target
            cout << "Producing subtracted cross section plot of " << targetNames[largerTarget] << " and " << targetNames[smallerTarget] << endl;

            CrossSection litData = 

            CrossSection subtracted = calculateSubtracted(crossSections[largerTarget],litData);
            string subtractedName = "#sigma_{" + allData[largerTarget].target.getName() + "}-" +
                                    "#sigma_{" + allData[smallerTarget].target.getName() + "}";
                                      
            subtracted.createCSGraph(relativeName.c_str());
            cout << "subtracted plot " << relative.getDataSet().getReference() <<
                    " RMS error: " << subtracted.calculateRMSError() << endl;
        }

        else
        {
            cout << "Failed to find cross section plot for either " << targetNames[largerTarget] << " or " << targetNames[smallerTarget] << endl;
        }
    }

    subtractedFile->Close();
*/
}
