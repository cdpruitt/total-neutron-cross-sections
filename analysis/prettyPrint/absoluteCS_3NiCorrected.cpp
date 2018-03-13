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

    // read graphs
    string expFileName = "/data1/analysis/corrected.root";
    string litFileName = "/data1/analysis/literatureData.root";

    TFile* expFile = new TFile(expFileName.c_str(),"READ");
    TFile* litFile = new TFile(litFileName.c_str(),"READ");
    
    string exp58GraphName = "Ni58_corrected";
    string expNatGraphName = "NiNat_corrected";
    string exp64GraphName = "Ni64_corrected";

    string litNiGraphName = "Natural Ni (n,tot)";
    
    TGraphAsymmErrors* exp58Graph = (TGraphAsymmErrors*)expFile->Get(exp58GraphName.c_str());
    TGraphAsymmErrors* expNatGraph = (TGraphAsymmErrors*)expFile->Get(expNatGraphName.c_str());
    TGraphAsymmErrors* exp64Graph = (TGraphAsymmErrors*)expFile->Get(exp64GraphName.c_str());

    TGraphAsymmErrors* litNiGraph = (TGraphAsymmErrors*)litFile->Get(litNiGraphName.c_str());

    // Set graph point and line characteristics
    exp58Graph->SetLineColor(kRed);
    exp58Graph->SetLineWidth(5);
    exp58Graph->SetLineStyle(0);
    exp58Graph->SetMarkerColor(kRed);

    expNatGraph->SetLineColor(kBlack);
    expNatGraph->SetLineWidth(5);
    expNatGraph->SetLineStyle(0);
    expNatGraph->SetMarkerColor(kBlack);

    exp64Graph->SetLineColor(kBlue);
    exp64Graph->SetLineWidth(5);
    exp64Graph->SetLineStyle(0);
    exp64Graph->SetMarkerColor(kBlue);

    litNiGraph->SetLineColor(kGray);
    litNiGraph->SetLineWidth(5);
    litNiGraph->SetLineStyle(2);
    litNiGraph->SetMarkerColor(kGray);

    // Pad dimensions and margins
    gPad->SetPad(0.005, 0.995, 0.995, 0.005);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.03);
    gPad->SetBottomMargin(0.15);
    gPad->SetTicky(2);

    // X-axis parameters
    exp58Graph->GetXaxis()->SetTitle("Energy (MeV)");
    exp58Graph->GetXaxis()->SetTitleSize(0.05);
    exp58Graph->GetXaxis()->SetTitleFont(2);
    exp58Graph->GetXaxis()->SetTitleOffset(1.4);
    exp58Graph->GetXaxis()->CenterTitle();

    exp58Graph->GetXaxis()->SetLabelOffset(0.01);
    exp58Graph->GetXaxis()->SetLabelSize(0.05);
    exp58Graph->GetXaxis()->SetLabelFont(2);

    exp58Graph->GetXaxis()->SetNdivisions(10);
    exp58Graph->GetXaxis()->SetTickLength(0.03);

    // Y-axis parameters
    exp58Graph->GetYaxis()->SetTitle("#sigma_{tot} (barns)");
    exp58Graph->GetYaxis()->SetTitleSize(0.06);
    exp58Graph->GetYaxis()->SetTitleFont(2);
    exp58Graph->GetYaxis()->SetTitleOffset(0.8);
    exp58Graph->GetYaxis()->CenterTitle();

    exp58Graph->GetYaxis()->SetLabelOffset(0.01);
    exp58Graph->GetYaxis()->SetLabelSize(0.05);

    exp58Graph->GetYaxis()->SetLabelFont(2);
    exp58Graph->GetYaxis()->SetNdivisions(10);
    exp58Graph->GetYaxis()->SetTickLength(0.02);

    exp58Graph->Draw("");
    expNatGraph->Draw("same");
    exp64Graph->Draw("same");
    litNiGraph->Draw("same");
    exp58Graph->Draw("same");
    expNatGraph->Draw("same");
    exp64Graph->Draw("same");

    gPad->SetLogx(1);
    
    exp58Graph->GetYaxis()->SetRangeUser(0.5,4.5);

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.035);
    latex.SetTextAlign(13); // align at top
    //latex.DrawLatex(0.65,0.65,"Pb (elem.)");
    //latex.DrawLatex(0.35,0.52,"Ni (elem.)");
    //latex.DrawLatex(0.32,0.4,"C (elem.)");

    // Define legend format and contents
    TLegend *legend = new TLegend(0.63, 0.68,0.9,0.9);
    legend->AddEntry(litNiGraph,"{}^{nat}Ni, literature","l");
    legend->AddEntry(exp58Graph,"{}^{58}Ni, lit corrected","l");
    legend->AddEntry(expNatGraph,"{}^{nat}Ni, lit corrected","l");
    legend->AddEntry(exp64Graph,"{}^{64}Ni, lit corrected","l");
    legend->Draw();

    expFile->Close();
    litFile->Close();
}
