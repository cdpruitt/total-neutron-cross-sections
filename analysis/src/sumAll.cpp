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

Config config;

const int MAX_SUBRUN_NUMBER = 200;

const int FIRST_CS_ENERGY = 2; // in MeV
const int LAST_CS_ENERGY = 600; // in MeV

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

CrossSection calculateCS(CSPrereqs& targetData, CSPrereqs& blankData, string expName)
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
    //monitorRatio *= 1.004;

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

        // ignore bins outside the energies of interest
        if(FIRST_CS_ENERGY>energyValue ||
                LAST_CS_ENERGY<energyValue)
        {
            continue;
        }

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

    string name = targetData.target.getName();
    crossSection.createCSGraph(name, name);

    //CrossSection corrected = correctForBlank(crossSection, volumeDensity, expName, "literatureData.root");
    //name += "blankCorrected";

    //corrected.createCSGraph(name, name);

    //return corrected;
    return crossSection;
}

int main(int, char* argv[])
{
    string dataLocation = argv[1];

    string expName = argv[2]; // experiment directory where runs to-be-sorted
                              // are listed

    vector<CSPrereqs> allCSPrereqs;

    // Open run list
    string runListName = "../" + expName + "/runsToSort.txt";
    ifstream runList(runListName);
    if(!runList.is_open())
    {
        cerr << "Error: couldn't find runlist at " << runListName << endl;
        exit(1);
    }
    cout << endl;

    // Ingest data from every run in the run list
    string runNumber;
    while (runList >> runNumber)
    {
        // read in run config file
        config = Config(expName, atoi(runNumber.c_str()));

        // check to see if data with this target has already been read in
        // (from previous runs). If not, create a new CSPrereq to hold the new
        // target's data
        for(string targetName : config.targetConfig.TARGET_ORDER)
        {
            bool CSPAlreadyExists = false;

            for(CSPrereqs csp : allCSPrereqs)
            {
                if(csp.target.getName() == targetName)
                {
                    // CSPrereq for this target already exists
                    CSPAlreadyExists = true;
                    break;
                }
            }

            if(CSPAlreadyExists)
            {
                continue;
            }

            string targetDataLocation = "../" + expName + "/targetData/" + targetName + ".txt";
            allCSPrereqs.push_back(CSPrereqs(targetDataLocation));
        }

        // Loop through all subruns of this run
        for(int subRun=0; subRun<=MAX_SUBRUN_NUMBER; subRun++)
        {
            stringstream subRunFormatted;
            subRunFormatted << setfill('0') << setw(4) << subRun;

            // open subrun
            stringstream monitorFileName;
            monitorFileName << dataLocation << "/" << runNumber << "/"
                       << subRunFormatted.str() << "/histos.root";
            ifstream f(monitorFileName.str());
            if(!f.good())
            {
                // failed to open this sub-run - skip to the next one
                cerr << "Couldn't open " << monitorFileName.str() << "; continuing.\r";
                fflush(stdout);
                continue;
            }

            TFile* monitorFile = new TFile(monitorFileName.str().c_str(),"READ");

            stringstream energyFileName;
            energyFileName << dataLocation << "/" << runNumber << "/"
                       << subRunFormatted.str() << "/energy.root";
            ifstream g(energyFileName.str());
            if(!g.good())
            {
                // failed to open this sub-run - skip to the next one
                cerr << "Couldn't open " << energyFileName.str() << "; continuing.\r";
                fflush(stdout);
                continue;
            }

            g.close();

            TFile* energyFile = new TFile(energyFileName.str().c_str(),"READ");

            // get target order for this run
            vector<string> targetOrder = getTargetOrder(expName, stoi(runNumber));

            cout << "Adding " << runNumber << ", subrun " << subRun << endl;

            // Loop through all target positions in this subrun
            for(int j=0; (size_t)j<targetOrder.size(); j++)
            {
                // pull data needed for CS calculation from subrun 
                string targetDataLocation = "../" + expName + "/targetData/" + targetOrder[j] + ".txt";
                CSPrereqs subRunData(targetDataLocation);

                // for histos
                subRunData.readEnergyData(energyFile, "summedDet", j);
                subRunData.readMonitorData(monitorFile, "monitor", j);

                // find the correct CSPrereqs to add this target's data to
                for(CSPrereqs& csp : allCSPrereqs)
                {
                    if(csp.target.getName() == subRunData.target.getName())
                    {
                        // add subrun data to total
                        csp = csp + subRunData;
                    }
                }
            }

            // Close the sub-run input files
            energyFile->Close();
            monitorFile->Close();
        }
    }

    string outFileName = dataLocation + "/total.root";
    TFile* outFile = new TFile(outFileName.c_str(), "UPDATE");

    vector<CrossSection> crossSections;

    cout << endl << "Total statistics over all runs: " << endl << endl;

    for(CSPrereqs& p : allCSPrereqs)
    {
        long totalCounts = 0;
        for(int i=0; i<p.energyHisto->GetNbinsX(); i++)
        {
            long tempCounts = p.energyHisto->GetBinContent(i);
            if(tempCounts < 0)
            {
                continue;
            }

            totalCounts += tempCounts;
        }

        cout << p.target.getName() << ": total events in energy histo = "
             << totalCounts << ", total monitor events = "
             << p.monitorCounts << endl;
        crossSections.push_back(calculateCS(p, allCSPrereqs[0], expName));
        cout << "Target " << crossSections.back().getDataSet().getReference() <<
                " RMS error: " << crossSections.back().calculateRMSError() << endl << endl;

        //p.energyHisto->SetDirectory(outFile);
        string name = p.target.getName() + "TOF";
        p.TOFHisto->SetNameTitle(name.c_str(),name.c_str());
        p.TOFHisto->SetDirectory(outFile);
    }

    outFile->Write();
    outFile->Close();

    /**************************************************************************
    Create relative cross section plots
    **************************************************************************/

    string relativeFileName = dataLocation + "/relative.root";
    TFile* relativeFile = new TFile(relativeFileName.c_str(), "UPDATE");

    // read which relative cross section plots to make from the experimental
    // directory
    vector<pair<string,string>> relativeTargets = getRelativePlotNames(expName,"relativePlots.txt");

    for(pair<string,string> p : relativeTargets)
    {
        int largerTarget = -1;
        int smallerTarget = -1;
        for(int i=0; (size_t)i<config.targetConfig.TARGET_ORDER.size(); i++)
        {
            if(config.targetConfig.TARGET_ORDER[i]==p.first)
            {
                largerTarget = i;
            }
            else if(config.targetConfig.TARGET_ORDER[i]==p.second)
            {
                smallerTarget = i;
            }
        }

        if(largerTarget>=0 && smallerTarget>=0)
        {
            // found cross section plots for both individual targets,
            // so use them to make a relative cross section plot
            cout << "Producing relative cross section plot of " << config.targetConfig.TARGET_ORDER[largerTarget] << " and " << config.targetConfig.TARGET_ORDER[smallerTarget] << endl;

            //CrossSection sum = crossSections[largerTarget]+crossSections[smallerTarget];
            //cout << "sum plot " << sum.getDataSet().getReference() <<
            //        " RMS error: " << sum.calculateRMSError() << endl << endl;

            //CrossSection difference = crossSections[largerTarget]-crossSections[smallerTarget];
            //cout << "difference plot " << difference.getDataSet().getReference() <<
            //        " RMS error: " << difference.calculateRMSError() << endl << endl;

            CrossSection relative = calculateRelative(crossSections[largerTarget],crossSections[smallerTarget]);
            string relativeTitle = "#frac{#sigma_{" + allCSPrereqs[largerTarget].target.getName() +
                                      "}-#sigma_{" + allCSPrereqs[smallerTarget].target.getName() + 
                                     "}}{#sigma_{" + allCSPrereqs[largerTarget].target.getName() +
                                      "}+#sigma_{" + allCSPrereqs[smallerTarget].target.getName() + "}}";
            string relativeName = "relative" + allCSPrereqs[largerTarget].target.getName() + allCSPrereqs[smallerTarget].target.getName();
            relative.createCSGraph(relativeName.c_str(), relativeTitle.c_str());
            cout << "Relative plot " << relative.getDataSet().getReference() <<
                    " RMS error: " << relative.calculateRMSError() << endl;
        }

        else
        {
            cout << "Failed to find cross section plot for either " << config.targetConfig.TARGET_ORDER[largerTarget] << " or " << config.targetConfig.TARGET_ORDER[smallerTarget] << endl;
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
        for(int i=0; (size_t)i<config.targetConfig.TARGET_ORDER.size(); i++)
        {
            if(config.targetConfig.TARGET_ORDER[i]==p.first)
            {
                largerTarget = i;
            }
            else if(config.targetConfig.TARGET_ORDER[i]==p.second)
            {
                smallerTarget = i;
            }
        }

        if(largerTarget>=0)
        {
            // found experimental cross section plot for this target
            cout << "Producing subtracted cross section plot of " << config.targetConfig.TARGET_ORDER[largerTarget] << " and " << config.targetConfig.TARGET_ORDER[smallerTarget] << endl;

            CrossSection litData = 

            CrossSection subtracted = calculateSubtracted(crossSections[largerTarget],litData);
            string subtractedName = "#sigma_{" + allCSPrereqs[largerTarget].target.getName() + "}-" +
                                    "#sigma_{" + allCSPrereqs[smallerTarget].target.getName() + "}";
                                      
            subtracted.createCSGraph(relativeName.c_str());
            cout << "subtracted plot " << relative.getDataSet().getReference() <<
                    " RMS error: " << subtracted.calculateRMSError() << endl;
        }

        else
        {
            cout << "Failed to find cross section plot for either " << config.targetConfig.TARGET_ORDER[largerTarget] << " or " << config.targetConfig.TARGET_ORDER[smallerTarget] << endl;
        }
    }

    subtractedFile->Close();
*/
}
