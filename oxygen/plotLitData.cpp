#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
/*#include "TMath.h"
#include "TAxis.h"
#include "TFile.h"
#include "TColor.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TStyle.h"
*/

/*****************************************************************************/
// Style section: edit to change style of plotted datasets
const std::string colors[4] = {"kRed","kRed","kRed","kRed"};
const std::string markers[4] = {"kPlus","kPlus","kPlus","kPlus"};
/*****************************************************************************/

class DataSet
{
    public:
        // constructors
        DataSet();
        DataSet(std::string dataSetLocation);

        TGraphErrors* getPlot();
        char* getReference();

    private:
        TGraphErrors* dataPlot;
        TColor* color;
        TAttMarker* marker;

        char reference[200];

        std::vector<float> energy;
        std::vector<float> xsection;
        std::vector<float> error;
};

DataSet::DataSet()
{
    std::cout << "Error: attempted to create DataSet without providing an input data file." << std::endl;
    exit(1);
}

DataSet::DataSet(std::string dataSetLocation)
{
    std::ifstream dataFile(dataSetLocation.c_str());
    if(!dataFile.is_open())
    {
        std::cout << "Attempted to create DataSet, but failed to find " << dataSetLocation << std::endl;
        exit(1);
    }

    char dummy[200];
    dataFile.getline(dummy,200);
    dataFile.getline(reference,200);
    dataFile.getline(dummy,200);

    float dum,dum2,dum3;

    while(dataFile >> dum >> dum2 >> dum3)
    {
        energy.push_back(dum);
        xsection.push_back(dum2);
        error.push_back(dum3);
    }

    dataPlot = new TGraphErrors(energy.size(),&energy[0],&xsection[0],0,&error[0]);

    dataPlot->GetXaxis()->SetTitle("Energy [MeV]");
    dataPlot->GetXaxis()->CenterTitle();
    dataPlot->GetXaxis()->SetRangeUser(0,700.);

    dataPlot->GetYaxis()->SetTitle("sigma [b]");
    dataPlot->GetYaxis()->CenterTitle();
}

TGraphErrors* DataSet::getPlot()
{
    return dataPlot;
}

char* DataSet::getReference()
{
    return reference;
}

int plotLitData(string input1, string input2/*, string input3, string input4*/)
{
    gStyle->SetOptStat(0);
    TStyle * Sty = (TStyle*)gROOT->FindObject("MyStyle");
    if(!Sty)      
    {
        Sty = new TStyle("MyStyle","MyStyle");
    }

    Sty->SetOptTitle(0);    
    Sty->SetOptStat(0);
    Sty->SetPalette(1,0);
    Sty->SetCanvasColor(10);      
    Sty->SetCanvasBorderMode(0);    
    Sty->SetFrameLineWidth(3);
    Sty->SetFrameFillColor(10);
    Sty->SetPadColor(10);
    Sty->SetPadTickX(1);
    Sty->SetPadTickY(1);
    Sty->SetPadBottomMargin(.15);
    Sty->SetPadTopMargin(.03);
    Sty->SetPadLeftMargin(.14);
    Sty->SetPadRightMargin(.06);
    Sty->SetHistLineWidth(3);
    Sty->SetHistLineColor(kBlue);
    Sty->SetFuncWidth(3);
    Sty->SetMarkerColor(kBlue);
    Sty->SetLineWidth(1);
    Sty->SetLabelSize(0.06,"xyz");
    Sty->SetLabelOffset(0.02,"y");
    Sty->SetLabelOffset(0.02,"x");
    Sty->SetLabelColor(kBlack,"xyz");
    Sty->SetMarkerSize(1);
    Sty->SetMarkerStyle(21);
    Sty->SetTitleSize(0.06,"xyz");
    Sty->SetTitleOffset(1.15,"y");
    Sty->SetTitleOffset(1.1,"x");
    Sty->SetTitleFillColor(10);
    Sty->SetTitleTextColor(kBlack);
    Sty->SetTickLength(.03,"xz");
    Sty->SetTickLength(.02,"y");
    Sty->SetNdivisions(5,"x");
    Sty->SetNdivisions(10,"yz");
    Sty->SetEndErrorSize(0);
    gROOT->SetStyle("MyStyle");
    gROOT->ForceStyle();

    // create output file
    std::stringstream fileOutName;
    fileOutName << "plotLitData.root";
    TFile* fileOut = new TFile(fileOutName.str().c_str(),"RECREATE");

    // read in data from each inputted dataset location
    std::vector<DataSet*> allData;
    allData.push_back(new DataSet(input1));
    ((TGraphErrors*)allData.back()->getPlot())->SetMarkerColor(kRed);

    /*allData.push_back(new DataSet(input2));
    ((TGraphErrors*)allData.back()->getPlot())->SetMarkerColor(kRed-7);

    allData.push_back(new DataSet(input3));
    ((TGraphErrors*)allData.back()->getPlot())->SetMarkerColor(kBlue);

    allData.push_back(new DataSet(input4));
    ((TGraphErrors*)allData.back()->getPlot())->SetMarkerColor(kBlue-7);
    */

    // plot all data together
    TCanvas* c1 = new TCanvas("xsPlots","Plots of literature data for total neutron cross-sections",1800,900);
    c1->SetLogx();

    TMultiGraph *allPlots = new TMultiGraph();
    for(int i=0; i<allData.size(); i++)
    {
        allPlots->Add(((TGraphErrors*)allData[i]->getPlot()),"p");
    }
    allPlots->Draw("a");

    allPlots->GetXaxis()->SetTitle("MeV");
    allPlots->GetXaxis()->CenterTitle();

    allPlots->GetYaxis()->SetTitle("Barns");
    allPlots->GetYaxis()->CenterTitle();

    legend = new TLegend(0.7,0.7,0.9,0.9);
    legend->SetHeader("Literature data");
    for(int i=0; i<allData.size(); i++)
    {
        legend->AddEntry(((TGraphErrors*)allData[i]->getPlot()),allData[i]->getReference(),"lep");
    }
    legend->Draw();

    allPlots->Write();

    //c1->Update();

    // clean up
    fileOut->Write();
    fileOut->Close();
    return 0;
}
