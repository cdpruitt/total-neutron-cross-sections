
/*****************************************************************************/
// Style section: edit to change style of plotted datasets
const std::string colors[4] = {"kRed","kRed","kRed","kRed"};
const std::string markers[4] = {"kPlus","kPlus","kPlus","kPlus"};
/*****************************************************************************/

// ROOT macro for formatting and combining plots

struct PlotParameters {
    string name = "default name";
    string title = "default title";
    bool plotLine = true;
    int lineColor = 0;
    double lineWidth = 1;
    bool plotPoint = false;;
    int pointColor = 0;
    double pointSize = 1;
    string XAxisTitle = "default x-axis title";
    string YAxisTitle = "default y-axis title";
};

void prettyPrintHisto(string histoFileName, string directoryName, string histoName, string formatFileName)
{
    // open file with formatting information
    ifstream formatFile(formatFileName);
    if(!formatFile.is_open())
    {
        cout << "Failed to open " << formatFileName << endl;
        exit(1);
    }

    string dummy;
    string dummy2;

    PlotParameters p;

    while(formatFile >> dummy >> dummy2)
    {
        if(dummy=="name")
        {
            p.name = dummy2;
        }

        else if(dummy=="title")
        {
            p.title = dummy2;
        }

        else if(dummy=="plotline")
        {
            if(dummy2=="true")
            {
                p.plotLine = true;
            }

            else if (dummy2=="false")
            {
                p.plotLine = false;
            }

            else
            {
                cerr << "Error: expected 'true' or 'false' for plotLine value, but received different string." << endl;
                p.plotLine = false;
            }
        }

        else if(dummy=="lineColor")
        {
            p.lineColor = stoi(dummy2);
        }

        else if(dummy=="lineWidth")
        {
            p.lineWidth = stod(dummy2);
        }

        else if(dummy=="plotPoint")
        {
            if(dummy2=="true")
            {
                p.plotPoint = true;
            }

            else if (dummy2=="false")
            {
                p.plotPoint = false;
            }

            else
            {
                cerr << "Error: expected 'true' or 'false' for plotPoint value, but received different string." << endl;
                p.plotPoint = false;
            }       
        }

        else if(dummy=="pointColor")
        {
            p.pointColor = stoi(dummy2);
        }

        else if(dummy=="pointSize")
        {
            p.pointSize = stod(dummy2);
        }

        else if(dummy=="XAxisTitle")
        {
            p.XAxisTitle = dummy2;
        }

        else if(dummy=="YAxisTitle")
        {
            p.YAxisTitle = dummy2;
        }
    }

    TFile* histoFile = new TFile(histoFileName.c_str(),"UPDATE");

    if(directoryName!="")
    {
        histoFile->cd(directoryName.c_str());
    }

    TH1I* histoToFormat;

    if(gDirectory->Get(histoName.c_str()))
    {
        histoToFormat = (TH1I*)gDirectory->Get(histoName.c_str());
    }

    else
    {
        cerr << "Error: failed to find " << histoName << " in " << directoryName << ". Exiting..." << endl;
        exit(1);
    }

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
    /*Sty->SetMarkerColor(kBlue);
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
    */
    gROOT->SetStyle("MyStyle");
    gROOT->ForceStyle();

    if(p.title!="")
    {
        histoToFormat->SetTitle(p.title.c_str());
    }

    if(p.plotLine)
    {
        histoToFormat->SetLineColor(p.lineColor);
        histoToFormat->SetLineWidth(p.lineWidth);
    }

    if(p.plotPoint)
    {
        histoToFormat->SetMarkerColor(p.pointColor);
        histoToFormat->SetMarkerSize(p.pointSize);
    }

    if(p.XAxisTitle!="")
    {
        histoToFormat->GetXaxis()->SetTitle(p.XAxisTitle.c_str());
        histoToFormat->GetXaxis()->CenterTitle();
        histoToFormat->GetXaxis()->SetTitleSize(0.05);
        histoToFormat->GetXaxis()->SetTitleOffset(1.5);
    }

    if(p.YAxisTitle!="")
    {
        histoToFormat->GetYaxis()->SetTitle(p.YAxisTitle.c_str());
        histoToFormat->GetYaxis()->CenterTitle();
        histoToFormat->GetYaxis()->SetTitleSize(0.04);
        histoToFormat->GetYaxis()->SetTitleOffset(3);
    }

    histoToFormat->Write();

    //TLine *line = new TLine(3,0.0345,300,0.0345);
    //line->SetLineStyle(7);
    //line->SetLineWidth(3);
    //line->SetLineColor(kBlack);
    //line->Draw();

    //TLegend *legend = new TLegend(0.2,0.7,0.4,0.9);
    //legend->SetHeader("data");
    //legend->AddEntry(relativeCSLogGraph,"experimental (preliminary)","p");
    //legend->AddEntry(plot1,"DOM calculation","l");
    //legend->AddEntry(line,"expected from size scaling","l");
    //legend->Draw();

    histoFile->Close();
}
