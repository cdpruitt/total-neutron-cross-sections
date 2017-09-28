// ROOT macro for pretty-printing histograms

// how to use:
//
// root
// .x prettyPrintHisto.cpp("histoFileName","directoryName","histoName","formatFileName")

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

    bool hasLegend = false;
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

    PlotParameters p;

    for(string dummy; getline(formatFile, dummy); )
    {
        // parse dummy string into label and content
        int endOfFirstWord = dummy.find_first_of(' ');
        string firstWord = dummy.substr(0,endOfFirstWord);

        string dummyRemainder = dummy.substr(endOfFirstWord+1);
        int startOfSecondWord = dummyRemainder.find_first_not_of(' ');
        string secondWord = dummyRemainder.substr(startOfSecondWord);

        cout << firstWord << " " << secondWord << endl;

        if(firstWord=="name")
        {
            p.name = secondWord;
        }

        else if(firstWord=="title")
        {
            p.title = secondWord;
        }

        else if(firstWord=="plotline")
        {
            if(secondWord=="true")
            {
                p.plotLine = true;
            }

            else if (secondWord=="false")
            {
                p.plotLine = false;
            }

            else
            {
                cerr << "Error: expected 'true' or 'false' for plotLine value, but received different string." << endl;
                p.plotLine = false;
            }
        }

        else if(firstWord=="lineColor")
        {
            p.lineColor = stoi(secondWord);
        }

        else if(firstWord=="lineWidth")
        {
            p.lineWidth = stod(secondWord);
        }

        else if(firstWord=="plotPoint")
        {
            if(secondWord=="true")
            {
                p.plotPoint = true;
            }

            else if (secondWord=="false")
            {
                p.plotPoint = false;
            }

            else
            {
                cerr << "Error: expected 'true' or 'false' for plotPoint value, but received different string." << endl;
                p.plotPoint = false;
            }       
        }

        else if(firstWord=="pointColor")
        {
            p.pointColor = stoi(secondWord);
        }

        else if(firstWord=="pointSize")
        {
            p.pointSize = stod(secondWord);
        }

        else if(firstWord=="XAxisTitle")
        {
            p.XAxisTitle = secondWord;
        }

        else if(firstWord=="YAxisTitle")
        {
            p.YAxisTitle = secondWord;
        }
    }

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
    Sty->SetFuncColor(kRed);
    Sty->SetLineWidth(1);
    Sty->SetLabelSize(0.06,"xyz");
    Sty->SetLabelOffset(0.02,"y");
    Sty->SetLabelOffset(0.02,"x");
    Sty->SetLabelColor(kBlack,"xyz");
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
        histoToFormat->GetXaxis()->SetTitleOffset(1.1);
    }

    if(p.YAxisTitle!="")
    {
        histoToFormat->GetYaxis()->SetTitle(p.YAxisTitle.c_str());
        histoToFormat->GetYaxis()->CenterTitle();
        histoToFormat->GetYaxis()->SetTitleSize(0.05);
        histoToFormat->GetYaxis()->SetTitleOffset(1.15);
    }

    // Define legend format and contents
    /*TLegend *legend = new TLegend(0.2,0.7,0.4,0.9);
    legend->SetHeader("data");
    legend->AddEntry(relativeCSLogGraph,"experimental (preliminary)","p");
    legend->AddEntry(plot1,"DOM calculation","l");
    legend->AddEntry(line,"expected from size scaling","l");
    legend->Draw();
    */

    histoToFormat->Write();

    histoFile->Close();
}
