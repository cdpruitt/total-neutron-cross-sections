// Plotting program for viewing multiple sub-runs in a single canvas (for identifying bad runs)
//
// To run, start ROOT and run the following command (including quotes):
//
// .x examineRuns.cpp("filepath","run","fileType","subDir","histoType")
//
// filepath = drive source (i.e., /data3)
// run = run to sort (i.e., 170)
// fileType = type of root file (i.e., histos)
// subDir = directory within ROOT file to look for plots (i.e., detS)
// histoType = name of histogram to plot (i.e., target1CS)

void examineRuns(string filePath, string run, string fileType, string subDir, string histoType)
{
    const unsigned int CANVAS_WIDTH = 6;
    const unsigned int CANVAS_HEIGHT = 4;

    //gROOT->SetStyle("Plain");
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

    // Create canvas
    TCanvas *canvas = new TCanvas("mycan","mycan",1800,900);
    canvas->Divide(CANVAS_WIDTH,CANVAS_HEIGHT,0.0001,0.0001);

    // Access existing histograms
    TFile* fileIn;
    stringstream fileName;

    TH1I* histoToAdd;
    TLatex latex;
    stringstream label;

    for(int i=0; i<CANVAS_WIDTH*CANVAS_HEIGHT; i++)
    {
        fileName.str("");
        fileName << filePath.c_str() << "/analysis/" << run.c_str() << "/"
            << std::setfill('0') << std::setw(4) << i << "/"
            << fileType << ".root";
        fileIn = new TFile(fileName.str().c_str(),"READ");
        if(!fileIn->IsOpen())
        {
            cout << "Can't open input file " << i << "." << endl;
            break;
        }

        canvas->cd(i+1);
        gPad->SetLogx();

        if(subDir.c_str())
        {
            gDirectory->cd("/");
            gDirectory->GetDirectory(subDir.c_str())->cd();
        }

        histoToAdd = (TH1I*)gDirectory->Get(histoType.c_str());

        histoToAdd->SetTitleSize(0.5,"t");
        histoToAdd->Draw();

        label.str("");
        label << std::setfill('0') << std::setw(4) << i;
        latex.SetTextSize(0.07);
        latex.DrawLatexNDC(0.2,0.05,label.str().c_str());
        latex.Draw();
    }
}
