// Map an input histogram with bins in the time domain to equivalent bins in the relativistic
// kinetic energy domain

#include "../include/target.h"
#include "../include/helperFunctions.h"
#include "../include/analysisConstants.h"
#include "../include/physicalConstants.h"
#include "../include/plottingConstants.h"
#include "../include/plots.h"
#include <iostream>
#include <vector>
#include "TH1.h"
#include "TFile.h"
#include "TAxis.h"

using namespace std;

CrossSection calculateCS(vector<Target*>& targets, vector<Plots*> plots, TFile* histoFile)
{
    // Find number of events in the monitor for each target to use in scaling
    // cross-sections

    // switch to the monitor directory
    histoFile->cd();
    histoFile->cd("/");
    histoFile->cd(dirs[1].c_str());

    // to normalize flux between all channels, we'll need to keep track of the
    // number of counts that recorded by the monitor paddle for each target

    for(int i=0; i<6; i++)
    {
        plots[i]->setMonitorCounts(((TH1I*)gDirectory->Get("targetPosH"))->GetBinContent(i+2));
        cout << "target position " << i << " monitor counts = " << plots[i]->getMonitorCounts() << endl;
    }

    // switch to the detector directory
    histoFile->cd("/");
    histoFile->cd(dirs[2].c_str());

    // Loop through the relativistic kinetic energy histograms and use them
    // to populate cross-section for each target

    // calculate cross sections for each bin of each target
    double crossSectionValue;
    double crossSectionError;
    double energyValue;
    double energyError;

    Target* blank = targets[0];
    Plots* blankPlots = plots[0];
    //TH1I* blankEnergy = blank->getEnergyHisto();
    TH1I* blankEnergyC = plots[0]->getEnergyHistoCorrected();

    for(int i=0; (size_t)i<targets.size(); i++)
    {
        Target* t = targets[i];
        Plots* p = plots[i];
        int numberOfBins = p->getEnergyHisto()->GetNbinsX();

        //TH1I* tof = p->getTOFHisto();
        //TH1I* tofC = p->getTOFHistoCorrected();
        //TH1I* en = p->getEnergyHisto();
        TH1I* enC = p->getEnergyHistoCorrected();

        long targetMonCounts = p->getMonitorCounts();
        long blankMonCounts = blankPlots->getMonitorCounts();
        if(targetMonCounts == 0 || blankMonCounts == 0)
        {
            cout << "Error - didn't find any monitor counts for target while trying to calculate cross sections. Exiting..." << endl;
            exit(1);
        }

        for(int j=1; j<numberOfBins; j++)
        {
            energyValue = enC->GetBinCenter(j);
            energyError = 0;

            // avoid "divide by 0" and "log of 0" errors
            if(blankEnergyC->GetBinContent(j) <= 0 || enC->GetBinContent(j) <= 0)
            {
                crossSectionValue = 0;
                crossSectionError = 0;
            }

            else
            {
                // calculate the cross section
                crossSectionValue =
                    -log(
                            ((double)enC->GetBinContent(j) // counts in target
                             /blankEnergyC->GetBinContent(j))// counts in blank
                            *(blankMonCounts/(double)targetMonCounts) // scale by monitor counts
                        )
                    /
                        (
                            t->getMass()
                            *AVOGADROS_NUMBER
                            *pow(10.,-24) // convert cm^2 to barns 
                            /
                                (pow(t->getDiameter()/2,2)*M_PI // area of cylinder end
                                *t->getMolMass())
                        );

                // calculate the statistical error
                crossSectionError =
                    pow((1/(double)enC->GetBinContent(j) 
                        +1/(double)blankEnergyC->GetBinContent(j)
                        +1/(double)blankMonCounts
                        +1/(double)targetMonCounts
                        ),0.5)
                        /((t->getMass())
                         *AVOGADROS_NUMBER
                         *pow(10.,-24) // convert cm^2 to barns 
                         *(p->getCrossSection().getDataPoint(j).getXValue()) // error of log(x) ~ (errorOfX)/x
                         /
                         ((pow(t->getDiameter()/2,2)*M_PI // area of cylinder end
                           *t->getMolMass())));
            }

            p->getCrossSection().addDataPoint(
                    DataPoint(crossSectionValue,crossSectionError,energyValue,energyError));

        }

        p->createCSGraph();
    }
}
