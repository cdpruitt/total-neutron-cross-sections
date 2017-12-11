#include "../include/CSUtilities.h"
#include "../include/crossSection.h"
#include "../include/dataSet.h"
#include "../include/dataPoint.h"
#include "../include/experiment.h"

#include <iostream>

using namespace std;

/*CrossSection correctForBlank(CrossSection rawCS, double targetNumberDensity, string expName, string graphFileName)
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
}*/

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
