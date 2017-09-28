void absoluteCS_SnCorrected_expLit()
{
    TStyle * style = (TStyle*)gROOT->FindObject("graphStyle");

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

    // Pad dimensions and margins
    gPad->SetPad(0.005, 0.995, 0.995, 0.005);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.03);
    gPad->SetBottomMargin(0.15);
    gPad->SetTicky(2);

    string expFileName = "/data2/analysis/scaledown.root";
    string litFileName = "/data2/analysis/literatureData.root";

    TFile* expFile = new TFile(expFileName.c_str(),"READ");
    TFile* litFile = new TFile(litFileName.c_str(),"READ");

    if(!expFile)
    {
        cerr << "Error: failed to open " << expFileName << endl;
        exit(1);
    }

    if(!litFile)
    {
        cerr << "Error: failed to open " << litFileName << endl;
        exit(1);
    }
    
    string expSn112GraphName = "Sn112_corrected";
    string expSnNatGraphName = "SnNat_corrected";
    string expSn124GraphName = "Sn124_corrected";
    string litSnNatGraphName = "NatSn_sd10";
 
    TGraphErrors* expSn112Graph = (TGraphErrors*)expFile->Get(expSn112GraphName.c_str());
    TGraphErrors* expSnNatGraph = (TGraphErrors*)expFile->Get(expSnNatGraphName.c_str());
    TGraphErrors* expSn124Graph = (TGraphErrors*)expFile->Get(expSn124GraphName.c_str());
    TGraphErrors* litSnNatGraph = (TGraphErrors*)litFile->Get(litSnNatGraphName.c_str());

    if(!expSn112Graph)
    {
        cerr << "Error: failed to open " << expSn112GraphName << " in " << expFile << endl;
        exit(1);
    }

    if(!expSn124Graph)
    {
        cerr << "Error: failed to open " << expSn124GraphName << " in " << expFile << endl;
        exit(1);
    }

    // Set graph point and line characteristics
    expSn112Graph->SetLineColor(kRed-9);
    expSn112Graph->SetLineWidth(5);
    expSn112Graph->SetLineStyle(0);
    expSn112Graph->SetMarkerColor(kRed-9);

    expSnNatGraph->SetLineColor(kRed);
    expSnNatGraph->SetLineWidth(5);
    expSnNatGraph->SetLineStyle(0);
    expSnNatGraph->SetMarkerColor(kRed);

    expSn124Graph->SetLineColor(kRed+3);
    expSn124Graph->SetLineWidth(5);
    expSn124Graph->SetLineStyle(0);
    expSn124Graph->SetMarkerColor(kRed+3);

    litSnNatGraph->SetLineColor(kBlack);
    litSnNatGraph->SetLineWidth(5);
    litSnNatGraph->SetLineStyle(2);
    litSnNatGraph->SetMarkerColor(kBlack);

    // X-axis parameters
    expSn112Graph->GetXaxis()->SetTitle("Energy (MeV)");
    expSn112Graph->GetXaxis()->SetTitleSize(0.05);
    expSn112Graph->GetXaxis()->SetTitleFont(2);
    expSn112Graph->GetXaxis()->SetTitleOffset(1.4);
    expSn112Graph->GetXaxis()->CenterTitle();

    expSn112Graph->GetXaxis()->SetLabelOffset(0.01);
    expSn112Graph->GetXaxis()->SetLabelSize(0.05);
    expSn112Graph->GetXaxis()->SetLabelFont(2);

    expSn112Graph->GetXaxis()->SetNdivisions(10);
    expSn112Graph->GetXaxis()->SetTickLength(0.03);

    // Y-axis parameters
    expSn112Graph->GetYaxis()->SetTitle("#sigma_{tot} (barns)");
    expSn112Graph->GetYaxis()->SetTitleSize(0.06);
    expSn112Graph->GetYaxis()->SetTitleFont(2);
    expSn112Graph->GetYaxis()->SetTitleOffset(0.8);
    expSn112Graph->GetYaxis()->CenterTitle();

    expSn112Graph->GetYaxis()->SetLabelOffset(0.01);
    expSn112Graph->GetYaxis()->SetLabelSize(0.05);

    expSn112Graph->GetYaxis()->SetLabelFont(2);
    expSn112Graph->GetYaxis()->SetNdivisions(10);
    expSn112Graph->GetYaxis()->SetTickLength(0.02);

    expSn112Graph->Draw();
    litSnNatGraph->Draw("same");
    expSn112Graph->Draw("same");
    expSn124Graph->Draw("same");
    expSnNatGraph->Draw("same");

    gPad->SetLogx(1);
    
    expSn112Graph->GetYaxis()->SetRangeUser(1.25,5.75);

    //TLatex latex;
    //latex.SetNDC();
    //latex.SetTextSize(0.05);
    //latex.SetTextAlign(13); // align at top
    //latex.DrawLatex(0.65,0.65,"Pb");
    //latex.DrawLatex(0.47,0.52,"Sn");
    //latex.DrawLatex(0.32,0.4,"C");

    // Define legend format and contents
    TLegend *legend = new TLegend(0.2,0.2,0.5,0.45);
    //legend->SetHeader("data","C");
    legend->AddEntry(expSn112Graph,"{}^{112}Sn (C+Pb corrected) ","l");
    legend->AddEntry(expSnNatGraph,"{}^{nat}Sn (C+Pb corrected) ","l");
    legend->AddEntry(expSn124Graph,"{}^{124}Sn (C+Pb corrected) ","l");
    legend->AddEntry(litSnNatGraph,"{}^{nat}Sn (lit.) ","l");
    legend->Draw();

    expFile->Close();
    litFile->Close();
}
