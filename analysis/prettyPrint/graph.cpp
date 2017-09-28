// ROOT macro for pretty-printing graphs

// how to use:
//
// root
// .x prettyPrintGraph.cpp("graphFileName","directoryName","graphName","formatFileName")

struct PlotParameters {
    string name = "default name";

    bool hasTitle = false;
    string title = "default title";

    bool plotLine = true;
    int lineColor = 0;
    double lineWidth = 1;
    int lineStyle = 0;

    bool plotPoint = true;
    int pointColor = 0;
    double pointSize = 1;
    int pointStyle = 5;

    double padLeftEdge = 0.005;
    double padRightEdge = 0.995;
    double padTopEdge = 0.995;
    double padBottomEdge = 0.005;

    double padLeftMargin = 0.1;
    double padRightMargin = 0.01;
    double padTopMargin = 0.01;
    double padBottomMargin = 0.1;

    string XAxisTitle = "default x-axis title";
    string YAxisTitle = "default y-axis title";

    double XAxisTitleOffset = 1;
    double YAxisTitleOffset = 1;

    double XAxisTitleSize = 0.05;
    double YAxisTitleSize = 0.05;

    double XAxisLabelOffset = 0.02;
    double YAxisLabelOffset = 0.02;

    double XAxisLabelSize = 0.04;
    double YAxisLabelSize = 0.03;

    int XAxisLabelFont = 2;
    int YAxisLabelFont = 2;

    int XAxisNDivisions = 10;
    int YAxisNDivisions = 10;

    double XAxisTickLength = 0.03;
    double YAxisTickLength = 0.02;

    bool hasLegend = false;
};

void prettyPrintGraph(string graphFileName, string directoryName, string graphName, string formatFileName, string outputFileName)
{
    // open file with formatting information
    ifstream formatFile(formatFileName);
    if(!formatFile.is_open())
    {
        cout << "Failed to open " << formatFileName << endl;
        exit(1);
    }

    // read in formatting information
    PlotParameters p;

    for(string dummy; getline(formatFile, dummy); )
    {
        // parse each line into "label" and "value"
        int endOfFirstWord = dummy.find_first_of(' ');
        string firstWord = dummy.substr(0,endOfFirstWord);

        string dummyRemainder = dummy.substr(endOfFirstWord+1);
        int startOfSecondWord = dummyRemainder.find_first_not_of(' ');
        string secondWord = dummyRemainder.substr(startOfSecondWord);

        if(firstWord=="name")
        {
            p.name = secondWord;
        }

        else if(firstWord=="hasTitle")
        {
            if(secondWord=="true")
            {
                p.hasTitle = true;
            }

            else if (secondWord=="false")
            {
                p.hasTitle = false;
            }

            else
            {
                cerr << "Error: expected 'true' or 'false' for hasTitle value, but received different string." << endl;
                p.hasTitle = false;
            }
        }

        else if (firstWord=="title")
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

        else if(firstWord=="lineStyle")
        {
            p.lineStyle = stoi(secondWord);
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

        else if(firstWord=="pointStyle")
        {
            p.pointStyle = stoi(secondWord);
        }

        else if(firstWord=="padLeftEdge")
        {
            p.padLeftEdge = stod(secondWord);
        }

        else if(firstWord=="padRightEdge")
        {
            p.padRightEdge = stod(secondWord);
        }

        else if(firstWord=="padBottomEdge")
        {
            p.padBottomEdge = stod(secondWord);
        }

        else if(firstWord=="padTopEdge")
        {
            p.padTopEdge = stod(secondWord);
        }

        else if(firstWord=="padLeftMargin")
        {
            p.padLeftMargin = stod(secondWord);
        }

        else if(firstWord=="padRightMargin")
        {
            p.padRightMargin = stod(secondWord);
        }

        else if(firstWord=="padBottomMargin")
        {
            p.padBottomMargin = stod(secondWord);
        }

        else if(firstWord=="padTopMargin")
        {
            p.padTopMargin = stod(secondWord);
        }

        else if(firstWord=="XAxisTitle")
        {
            p.XAxisTitle = secondWord;
        }

        else if(firstWord=="YAxisTitle")
        {
            p.YAxisTitle = secondWord;
        }

        else if(firstWord=="XAxisTitleOffset")
        {
            p.XAxisTitleOffset = stod(secondWord);
        }

        else if(firstWord=="YAxisTitleOffset")
        {
            p.YAxisTitleOffset = stod(secondWord);
        }

        else if(firstWord=="XAxisTitleSize")
        {
            p.XAxisTitleSize = stod(secondWord);
        }

        else if(firstWord=="YAxisTitleSize")
        {
            p.YAxisTitleSize = stod(secondWord);
        }

        else if(firstWord=="XAxisLabelOffset")
        {
            p.XAxisLabelOffset = stod(secondWord);
        }

        else if(firstWord=="YAxisLabelOffset")
        {
            p.YAxisLabelOffset = stod(secondWord);
        }

        else if(firstWord=="XAxisLabelSize")
        {
            p.XAxisLabelSize = stod(secondWord);
        }

        else if(firstWord=="YAxisLabelSize")
        {
            p.YAxisLabelSize = stod(secondWord);
        }

        else if(firstWord=="XAxisLabelFont")
        {
            p.XAxisLabelFont = stoi(secondWord);
        }

        else if(firstWord=="YAxisLabelFont")
        {
            p.YAxisLabelFont = stoi(secondWord);
        }

        else if(firstWord=="XAxisNDivisions")
        {
            p.XAxisNDivisions = stod(secondWord);
        }

        else if(firstWord=="YAxisNDivisions")
        {
            p.YAxisNDivisions = stod(secondWord);
        }

        else if(firstWord=="XAxisTickLength")
        {
            p.XAxisTickLength = stod(secondWord);
        }

        else if(firstWord=="YAxisTickLength")
        {
            p.YAxisTickLength = stod(secondWord);
        }
    }

    TStyle * Sty = (TStyle*)gROOT->FindObject("MyStyle");

    if(!Sty)      
    {
        Sty = new TStyle("MyStyle","MyStyle");
    }

    TCanvas* c = new TCanvas("c1","",1000,800);

    TFile* outputFile = new TFile(outputFileName.c_str(),"UPDATE");
    TFile* graphFile = new TFile(graphFileName.c_str(),"READ");

    if(directoryName!="")
    {
        graphFile->cd(directoryName.c_str());
    }

    TGraphErrors* graphToFormat;

    if(gDirectory->Get(graphName.c_str()))
    {
        graphToFormat = (TGraphErrors*)gDirectory->Get(graphName.c_str());
    }

    else
    {
        cerr << "Error: failed to find " << graphName << " in " << directoryName << ". Exiting..." << endl;
        exit(1);
    }

    if(p.hasTitle)
    {
        if(p.title!="")
        {
            graphToFormat->SetTitle(p.title.c_str());
        }
    }

    else
    {
        Sty->SetOptTitle(0);    
    }

    if(p.plotLine)
    {
        graphToFormat->SetLineColor(p.lineColor);
        graphToFormat->SetLineWidth(p.lineWidth);
        graphToFormat->SetLineStyle(p.lineStyle);
    }

    if(p.plotPoint)
    {
        graphToFormat->SetMarkerColor(p.pointColor);
        graphToFormat->SetMarkerSize(p.pointSize);
        graphToFormat->SetMarkerStyle(p.pointStyle);
    }

    if(p.padLeftEdge && p.padRightEdge && p.padTopEdge && p.padBottomEdge)
    {
        gPad->SetPad(p.padLeftEdge, p.padRightEdge, p.padTopEdge, p.padBottomEdge);
    }

    if(p.padLeftMargin && p.padRightMargin && p.padTopMargin && p.padBottomMargin)
    {
        gPad->SetLeftMargin(p.padLeftMargin);
        gPad->SetRightMargin(p.padRightMargin);
        gPad->SetTopMargin(p.padTopMargin);
        gPad->SetBottomMargin(p.padBottomMargin);
    }

    if(p.XAxisTitle!="")
    {
        graphToFormat->GetXaxis()->SetTitle(p.XAxisTitle.c_str());
        graphToFormat->GetXaxis()->CenterTitle();

        if(p.XAxisTitleSize)
        {
            graphToFormat->GetXaxis()->SetTitleSize(p.XAxisTitleSize);
        }

        if(p.XAxisTitleOffset)
        {
            graphToFormat->GetXaxis()->SetTitleOffset(p.XAxisTitleOffset);
        }
    }

    if(p.YAxisTitle!="")
    {
        graphToFormat->GetYaxis()->SetTitle(p.YAxisTitle.c_str());
        graphToFormat->GetYaxis()->CenterTitle();

        if(p.YAxisTitleSize)
        {
            graphToFormat->GetYaxis()->SetTitleSize(p.YAxisTitleSize);
        }

        if(p.YAxisTitleOffset)
        {
            graphToFormat->GetYaxis()->SetTitleOffset(p.YAxisTitleOffset);
        }
    }

    if(p.XAxisLabelOffset)
    {
        graphToFormat->GetXaxis()->SetLabelOffset(p.XAxisLabelOffset);
    }

    if(p.YAxisLabelOffset)
    {
        graphToFormat->GetYaxis()->SetLabelOffset(p.XAxisLabelOffset);
    }

    if(p.XAxisLabelSize)
    {
        graphToFormat->GetXaxis()->SetLabelSize(p.XAxisLabelSize);
    }

    if(p.YAxisLabelSize)
    {
        graphToFormat->GetYaxis()->SetLabelSize(p.XAxisLabelSize);
    }

    if(p.XAxisLabelFont)
    {
        graphToFormat->GetXaxis()->SetLabelFont(p.XAxisLabelFont);
    }

    if(p.YAxisLabelFont)
    {
        graphToFormat->GetYaxis()->SetLabelFont(p.XAxisLabelFont);
    }

    if(p.XAxisNDivisions)
    {
        graphToFormat->GetXaxis()->SetNdivisions(p.XAxisNDivisions);
    }

    if(p.YAxisNDivisions)
    {
        graphToFormat->GetYaxis()->SetNdivisions(p.YAxisNDivisions);
    }

    if(p.XAxisTickLength)
    {
        graphToFormat->GetXaxis()->SetTickLength(p.XAxisTickLength);
    }

    if(p.YAxisTickLength)
    {
        graphToFormat->GetYaxis()->SetTickLength(p.YAxisTickLength);
    }

    Sty->SetOptStat(0);
    Sty->SetPalette(1,0);
    Sty->SetCanvasColor(10);      
    Sty->SetCanvasBorderMode(0);    
    Sty->SetFrameLineWidth(3);
    Sty->SetFrameFillColor(10);
    Sty->SetPadColor(10);
    Sty->SetPadTickX(1);
    Sty->SetPadTickY(1);
    Sty->SetHistLineWidth(3);
    Sty->SetHistLineColor(kBlue);
    Sty->SetFuncWidth(3);
    Sty->SetFuncColor(kRed);
    Sty->SetLabelColor(kBlack,"xyz");
    //Sty->SetTitleSize(0.06,"xyz");
    Sty->SetTitleFillColor(10);
    Sty->SetTitleTextColor(kBlack);
    Sty->SetEndErrorSize(0);
    gROOT->SetStyle("MyStyle");
    gROOT->ForceStyle();

    graphToFormat->Draw();
    gPad->SetLogx(1);

    // Define legend format and contents
    TLegend *legend = new TLegend(0.2,0.7,0.4,0.9);
    legend->SetHeader("data");
    legend->AddEntry(graphToFormat,"{nat}^Pb","l");
    legend->Draw();

    //outputFile->cd();
    //graphToFormat->Write();

    graphFile->Close();
    outputFile->Close();
}
