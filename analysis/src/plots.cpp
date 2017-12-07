#include "TRandom3.h"

#include "../include/target.h"
#include "../include/plots.h"
#include "../include/physicalConstants.h"
#include "../include/config.h"

#include <iostream>

using namespace std;

TH1D* convertTOFtoEnergy(TH1D* tof, string name)
{
    TH1D* energy = timeBinsToRKEBins(tof, name); 

    if(!tof)
    {
        cerr << "Error: cannot convert empty TOF histogram to energy units in convertTOFtoEnergy()" << endl;
        return energy;
    }

    unsigned int tofBins = tof->GetNbinsX();

    TRandom3 *randomizeBin = new TRandom3();

    for(unsigned int j=1; j<=tofBins; j++)
    {
        // convert time into neutron velocity based on flight path distance
        double velocity = pow(10.,7.)*(config.facility.FLIGHT_DISTANCE)
            /(tof->GetBinCenter(j)
                    /*+randomizeBin->Uniform(
                        -(1/(double)(2*config.plot.TOF_BINS_PER_NS)),
                         (1/(double)(2*config.plot.TOF_BINS_PER_NS)))
                         */
             ); // in meters/sec 

        // convert velocity to relativistic kinetic energy
        double rKE = (pow((1.-pow((velocity/C),2.)),-0.5)-1.)*NEUTRON_MASS; // in MeV

        energy->Fill(rKE,tof->GetBinContent(j));
    }

    unsigned int energyBins = energy->GetNbinsX();
    for(unsigned int j=1; j<=energyBins; j++)
    {
        energy->SetBinError(j,pow(energy->GetBinContent(j),0.5));
    }

    return energy;
}

vector<double> scaleBins(vector<double> inputBins, int scaledown)
{
    vector<double> outputBins;

    for(size_t i=0; i<inputBins.size(); i+=scaledown)
    {
        outputBins.push_back(inputBins[i]);
    }

    // add the max edge of the input bins
    outputBins.push_back(inputBins[inputBins.size()-1]);

    return outputBins;
}

double tofToRKE(double TOF)
{
    double velocity = pow(10.,7.)*config.facility.FLIGHT_DISTANCE/TOF; // in meters/sec 

    if (velocity>C)
    {
        return -1;
    }

    // convert velocity to relativistic kinetic energy
    double RKE = (pow((1.-pow((velocity/C),2.)),-0.5)-1.)*NEUTRON_MASS; // in MeV
    if(RKE<0)
    {
        return -1;
    }

    return RKE;
}

double RKEToTOF(double RKE)
{
    // convert relativistic kinetic energy to velocity
    double velocity = pow(1-pow((1/((RKE/NEUTRON_MASS)+1)),2),0.5)*C;

    if(velocity<0 || velocity>C)
    {
        return -1;
    }

    double TOF = pow(10.,7.)*config.facility.FLIGHT_DISTANCE/velocity; // in meters/sec 

    return TOF; // in ns
}

TH1D* timeBinsToRKEBins(TH1D* inputHisto, string name)
{
    if(!inputHisto)
    {
        cerr << "Error: tried to convert time bins to energy bins for histo " << name << ", but histo pointer was null." << endl;
    }

    int nOldBins = inputHisto->GetNbinsX();
    TAxis* oldAxis = inputHisto->GetXaxis();

    double minimumTime;
    int minimumBin;

    for(int i=1; i<=nOldBins; i++)
    {
        minimumTime = inputHisto->GetBinLowEdge(i);
        minimumBin = i;

        if(tofToRKE(minimumTime)>0 && tofToRKE(minimumTime)<config.plot.ENERGY_UPPER_BOUND)
        {
            break;
        }
    }

    if(tofToRKE(minimumTime)==-1)
    {
        cerr << "Error: energy of old min time " << minimumTime << " was not finite." << endl;
        exit(1);
    }

    double maximumTime;
    int maximumBin;

    for(int i=nOldBins; i>=1; i--)
    {
        maximumTime = inputHisto->GetBinLowEdge(i) + inputHisto->GetBinWidth(i);
        maximumBin = i;

        if(tofToRKE(maximumTime)>config.plot.ENERGY_LOWER_BOUND)
        {
            break;
        }
    }

    if(tofToRKE(maximumTime)==-1)
    {
        cerr << "Error: energy of old maximum time " << maximumTime << " was not finite." << endl;
        exit(1);
    }

    // Remap bins from old histo to new histo
    int nUnscaledEnergyBins = maximumBin-minimumBin+1;
    vector<double> unscaledEnergyBins;

    // Reorder bins to go from lowest energy (shortest time) to highest energy (longest time)
    // n bins are defined for n+1 points (like fence sections and fence posts)
    for(int i=0; i<nUnscaledEnergyBins; i++)
    {
        double newBin = tofToRKE(oldAxis->GetBinLowEdge(maximumBin-i)+oldAxis->GetBinWidth(maximumBin-i));

        if(newBin<=0)
        {
            cerr << "Error: tried to make negative energy bin." << endl;
            exit(1);
        }

        unscaledEnergyBins.push_back(newBin);
    }

    unscaledEnergyBins.push_back(tofToRKE(oldAxis->GetBinLowEdge(minimumBin)));
    
    // Downscale bins to desired granularity
    //vector<double> scaledEnergyBins = unscaledEnergyBins;
    vector<double> scaledEnergyBins = scaleBins(unscaledEnergyBins, unscaledEnergyBins.size()/config.plot.NUMBER_ENERGY_BINS);

    TH1D* outputHisto = new TH1D(name.c_str(),
            name.c_str(),
            scaledEnergyBins.size()-1,
            &scaledEnergyBins[0]);
            //newXMin,
            //newXMax);

    return outputHisto;
}

TH1D* RKEBinsToTimeBins(TH1D *inputHisto, string name)
{
    // extract the total number of bins in the input Histo (minus the
    // overflow and underflow bins)
    int nOldBins = inputHisto->GetSize()-2;

    double minimumEnergy = (((TAxis*)inputHisto->GetXaxis())->GetXmin());
    int minimumBin = 0;

    for(int i=0; i<nOldBins; i++)
    {
        if(RKEToTOF(minimumEnergy)>0 && RKEToTOF(minimumEnergy)<config.plot.TOF_UPPER_BOUND)
        {
            break;
        }

        minimumEnergy = inputHisto->GetBinLowEdge(i);
        minimumBin = i;
    }

    double tentativeTime = RKEToTOF(minimumEnergy);
    if(tentativeTime==-1)
    {
        cerr << "Error: time of old min energy " << minimumEnergy << " was not finite: " << tentativeTime << " (ns)" << endl;
        exit(1);
    }

    //double newXMax = tentativeEnergy;

    double maximumEnergy = (((TAxis*)inputHisto->GetXaxis())->GetXmax());
    int maximumBin = nOldBins;

    for(int i=nOldBins; i>0; i--)
    {
        if(RKEToTOF(maximumEnergy)>config.plot.TOF_LOWER_BOUND)
        {
            break;
        }
        maximumEnergy = inputHisto->GetBinLowEdge(i);
        maximumBin = i;
    }

    tentativeTime = RKEToTOF(maximumEnergy);
    if(tentativeTime==-1)
    {
        cerr << "Error: time of old maximum energy " << maximumEnergy << " was not finite: " << tentativeTime << " (ns)" << endl;
        exit(1);
    }

    //double newXMin = tentativeEnergy;

    TAxis* oldAxis = inputHisto->GetXaxis();

    // Remap bins from old histo to new histo
    int nUnscaledTimeBins = maximumBin-minimumBin;
    vector<double> unscaledTimeBins;

    // Reorder bins to go from lowest time (highest energy) to highest time (lowest energy)
    // n bins are defined n+1 points (like fence sections and fence posts)
    for(int i=0; i<nUnscaledTimeBins+1; i++)
    {
        double newBin = RKEToTOF(oldAxis->GetBinLowEdge(maximumBin-i));
        if(newBin<=0)
        {
            continue;
        }
        unscaledTimeBins.push_back(newBin);
    }

    // Downscale bins to desired granularity
    vector<double> scaledTimeBins = unscaledTimeBins;
    //scaleBins(unscaledTimeBins, scaledTimeBins, nUnscaledTimeBins/NUMBER_TOF_BINS);

    TH1D* outputHisto = new TH1D(name.c_str(),
            name.c_str(),
            scaledTimeBins.size()-1,
            &scaledTimeBins[0]);
            //newXMin,
            //newXMax);

    // Assign the remapped bins to the new histo
    //TH1* outputHistoNonZero = outputHisto->Rebin(scaledEnergyBins.size()-2,"outputHistoNonZero",&scaledEnergyBins[0]);

    //double test = outputHistoNonZero->GetXaxis()->GetBinLowEdge(scaledEnergyBins.size()-2);
    //double test2 = outputHistoNonZero->GetXaxis()->GetBinLowEdge(0);

    //return outputHistoNonZero;
    return outputHisto;
}

// create new plots
Plots::Plots(string name)
{
    string tofName = name + "TOF";
    string rawTOFName = name + "rawTOF";
    string energyName = name + "Energy";
    string deadtimeName = name + "Deadtime";

    TOFHisto = new TH1D(tofName.c_str(),tofName.c_str(),config.plot.TOF_BINS,config.plot.TOF_LOWER_BOUND,config.plot.TOF_UPPER_BOUND);
    rawTOFHisto = new TH1D(rawTOFName.c_str(),rawTOFName.c_str(),config.plot.TOF_BINS,config.plot.TOF_LOWER_BOUND,config.plot.TOF_UPPER_BOUND);
    deadtimeHisto = new TH1D(deadtimeName.c_str(),deadtimeName.c_str(),config.plot.TOF_BINS,config.plot.TOF_LOWER_BOUND,config.plot.TOF_UPPER_BOUND);

    energyHisto = timeBinsToRKEBins(TOFHisto,energyName);

    //energyHisto = new TH1D(energyName.c_str(),energyName.c_str(),ENERGY_BINS,0,ENERGY_RANGE);
    //TOFHisto = RKEBinsToTimeBins(energyHisto,tofName);
    //deadtimeHisto = RKEBinsToTimeBins(energyHisto,deadtimeName);
}

TH1D* Plots::getTOFHisto()
{
    return TOFHisto;
}

TH1D* Plots::getRawTOFHisto()
{
    return rawTOFHisto;
}

TH1D* Plots::getEnergyHisto()
{
    return energyHisto;
}

TH1D* Plots::getDeadtimeHisto()
{
    return deadtimeHisto;
}

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


/*CrossSection calculateCS(const CSPrereqs& targetData, const CSPrereqs& blankData)
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

    // calculate the ratio of target/blank monitor counts (normalize
    // flux/macropulse)
    double tMon = targetData.monitorCounts;
    double bMon = blankData.monitorCounts;
    double monitorRatio = tMon/bMon;

    // calculate the ratio of target/blank good macropulse ratio (normalize
    // macropulse number)
    double tGMN = targetData.goodMacroNumber;
    double bTMN = blankData.totalMacroNumber;
    double macroNumberRatio = (targetData.goodMacroNumber/targetData.totalMacroNumber)
        /(blankData.goodMacroNumber/blankData.totalMacroNumber);
    //macroNumberRatio = 1;

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
    for(int i=1; i<=numberOfBins; i++) // skip the overflow and underflow bins
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
            -log(detectorRatio/(macroNumberRatio*monitorRatio))/arealDensity; // in cm^2

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
}*/

/*int producePlots(string dataLocation, const vector<CSPrereqs>& allCSPrereqs)
{
    string outFileName = dataLocation + "/total.root";
    TFile* outFile = new TFile(outFileName.c_str(), "UPDATE");

    vector<CrossSection> crossSections;

    cout << endl << "Total statistics over all runs: " << endl << endl;

    CSPrereqs blank;
    for(auto& p : allCSPrereqs)
    {
        if(p.target.getName()=="blank")
        {
            blank = p;
            break;
        }
    }

    for(auto& p : allCSPrereqs)
    {
        cout << "all CS Prereqs name = " << p.target.getName() << endl;

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
            << p.monitorCounts << ", good macro number = "
            << p.goodMacroNumber << ", total macro number = "
            << p.totalMacroNumber << endl;
        crossSections.push_back(calculateCS(p, blank));
        cout << "Target " << crossSections.back().getDataSet().getReference() <<
            " RMS error: " << crossSections.back().calculateRMSError() << endl << endl;

        p.energyHisto->SetDirectory(outFile);
        p.energyHisto->Write();

        string name = p.target.getName() + "TOF";
        p.TOFHisto->SetNameTitle(name.c_str(),name.c_str());
        p.TOFHisto->SetDirectory(outFile);

        p.TOFHisto->Write();
    }

    outFile->Close();

    //  Create relative cross section plots

    string relativeFileName = dataLocation + "/relative.root";
    TFile* relativeFile = new TFile(relativeFileName.c_str(), "UPDATE");

    // read which relative cross section plots to make from the experimental
    // directory
    vector<pair<string,string>> relativeTargets = getRelativePlotNames(expName,"relativePlots.txt");

    for(pair<string,string> p : relativeTargets)
    {
        int largerTarget = -1;
        int smallerTarget = -1;
        for(int i=0; (size_t)i<config.target.TARGET_ORDER.size(); i++)
        {
            if(config.target.TARGET_ORDER[i]==p.first)
            {
                largerTarget = i;
            }
            else if(config.target.TARGET_ORDER[i]==p.second)
            {
                smallerTarget = i;
            }
        }

        if(largerTarget>=0 && smallerTarget>=0)
        {
            // found cross section plots for both individual targets,
            // so use them to make a relative cross section plot
            cout << "Producing relative cross section plot of " << config.target.TARGET_ORDER[largerTarget] << " and " << config.target.TARGET_ORDER[smallerTarget] << endl;

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
            cout << "Failed to find cross section plot for either " << config.target.TARGET_ORDER[largerTarget] << " or " << config.target.TARGET_ORDER[smallerTarget] << endl;
        }
    }

    relativeFile->Close();

       string subtractedFileName = dataLocation + "/subtracted.root";
       TFile* subtractedFile = new TFile(subtractedFileName.c_str(), "RECREATE");

    // read which subtracted cross section plots to make from the experimental
    // directory
    vector<pair<string,string>> subtractedTargets = getSubtractedPlotNames(expName,"subtractedPlots.txt");

    for(pair<string,string> p : subtractedTargets)
    {
    int largerTarget = -1;
    int smallerTarget = -1;
    for(int i=0; (size_t)i<config.target.TARGET_ORDER.size(); i++)
    {
    if(config.target.TARGET_ORDER[i]==p.first)
    {
    largerTarget = i;
    }
    else if(config.target.TARGET_ORDER[i]==p.second)
    {
    smallerTarget = i;
    }
    }

    if(largerTarget>=0)
    {
    // found experimental cross section plot for this target
    cout << "Producing subtracted cross section plot of " << config.target.TARGET_ORDER[largerTarget] << " and " << config.target.TARGET_ORDER[smallerTarget] << endl;

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
    cout << "Failed to find cross section plot for either " << config.target.TARGET_ORDER[largerTarget] << " or " << config.target.TARGET_ORDER[smallerTarget] << endl;
    }
    }

    subtractedFile->Close();

    return 0;
}*/
