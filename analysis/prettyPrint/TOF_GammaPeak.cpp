void TOF_GammaPeak()
{
    TStyle * style = (TStyle*)gROOT->FindObject("histoStyle");

    if(!style)      
    {
        style = new TStyle("graphStyle","graphStyle");
    }

    TCanvas* c = new TCanvas("c1","",1200,1200);

    style->SetOptStat(0);
    style->SetOptTitle(0);    
    style->SetPalette(1,0);
    style->SetCanvasColor(10);      
    style->SetCanvasBorderMode(0);    
    style->SetFrameLineWidth(3);
    style->SetFrameFillColor(10);
    style->SetPadColor(10);
    style->SetHistLineWidth(3);
    style->SetHistLineColor(kBlue);
    style->SetMarkerSize(0.9);
    style->SetMarkerStyle(8);
    style->SetFuncWidth(3);
    style->SetFuncColor(kRed);
    style->SetLabelColor(kBlack,"xyz");
    style->SetTitleSize(0.06,"xyz");
    style->SetTitleFillColor(10);
    style->SetTitleTextColor(kBlack);
    style->SetEndErrorSize(0);

    gROOT->SetStyle("graphStyle");
    gROOT->ForceStyle();

    string TOFFileName = "/data2/analysis/66/0001/histos.root";

    TFile* TOFFile = new TFile(TOFFileName.c_str(),"READ");

    if(!TOFFile)
    {
        cerr << "Error: failed to open " << TOFFileName << endl;
        exit(1);
    }

    string blankTOFName = "blankTOF";

    TOFFile->cd("lowThresholdDet");
     
    TH1D* blankTOF = (TH1D*)gDirectory->Get(blankTOFName.c_str());
    
    if(!blankTOF)
    {
        cout << "Error: failed to find " << blankTOFName << " in file " << TOFFileName << endl;
        exit(1);
    }

    // Pad dimensions and margins
    gPad->SetPad(0.005, 0.995, 0.995, 0.005);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.03);
    gPad->SetBottomMargin(0.15);
    gPad->SetTicky(2);

    // Set histo point and line characteristics
    blankTOF->SetLineColor(kRed);
    blankTOF->SetLineWidth(5);
    blankTOF->SetLineStyle(0);
    blankTOF->SetMarkerColor(kRed);

    // X-axis parameters
    blankTOF->GetXaxis()->SetTitle("Time of flight (ns)");
    blankTOF->GetXaxis()->SetTitleSize(0.05);
    blankTOF->GetXaxis()->SetTitleFont(2);
    blankTOF->GetXaxis()->SetTitleOffset(1.5);
    blankTOF->GetXaxis()->CenterTitle();

    blankTOF->GetXaxis()->SetLabelOffset(0.01);
    blankTOF->GetXaxis()->SetLabelSize(0.05);
    blankTOF->GetXaxis()->SetLabelFont(2);

    blankTOF->GetXaxis()->SetNdivisions(10);
    blankTOF->GetXaxis()->SetTickLength(0.03);

    // Y-axis parameters
    blankTOF->GetYaxis()->SetTitle("Counts");
    blankTOF->GetYaxis()->SetTitleSize(0.06);
    blankTOF->GetYaxis()->SetTitleFont(2);
    blankTOF->GetYaxis()->SetTitleOffset(1.3);
    blankTOF->GetYaxis()->CenterTitle();

    blankTOF->GetYaxis()->SetLabelOffset(0.01);
    blankTOF->GetYaxis()->SetLabelSize(0.05);

    blankTOF->GetYaxis()->SetLabelFont(2);
    blankTOF->GetYaxis()->SetNdivisions(10);
    blankTOF->GetYaxis()->SetTickLength(0.02);

    blankTOF->Draw();
        
    blankTOF->GetYaxis()->SetRangeUser(0,7000);
    blankTOF->GetXaxis()->SetRange(1640, 1780);

    //TLatex latex;
    //latex.SetNDC();
    //latex.SetTextSize(0.05);
    //latex.SetTextAlign(13); // align at top
    //latex.DrawLatex(0.65,0.65,"Pb");
    //latex.DrawLatex(0.47,0.52,"Sn");
    //latex.DrawLatex(0.32,0.4,"C");

    //gPad->SetLogy(1);

    // Define legend format and contents
    //TLegend *legend = new TLegend(0.75,0.7,0.9,0.9);
    //legend->SetHeader("data","C");
    //legend->AddEntry(blankTOF,"blank","l");
    //legend->Draw();
}
