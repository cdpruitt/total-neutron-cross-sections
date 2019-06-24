{
    string fileName = "/data2/analysis/relative.root";
    TFile* file = new TFile(fileName.c_str(),"READ");
    
    string expGraphName = "DtoH_exp, percent";
    string litGraphName = "DtoH_lit, percent";
        
    TGraphAsymmErrors* expGraph = (TGraphAsymmErrors*)file->Get(expGraphName.c_str());
    TGraphAsymmErrors* litGraph = (TGraphAsymmErrors*)file->Get(litGraphName.c_str());

    TStyle* style = (TStyle*)gROOT->FindObject("graphStyle");

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

    // Set graph point and line characteristics
    expGraph->SetLineColor(kRed);
    expGraph->SetLineWidth(4);
    expGraph->SetLineStyle(0);
    expGraph->SetMarkerColor(kRed);
    expGraph->SetMarkerStyle(21);
    expGraph->SetMarkerSize(2);

    litGraph->SetLineColor(kBlack);
    litGraph->SetLineWidth(4);
    litGraph->SetLineStyle(0);
    litGraph->SetMarkerColor(kBlack);
    litGraph->SetMarkerStyle(20);
    litGraph->SetMarkerSize(2);

    // Pad dimensions and margins
    gPad->SetPad(0.005, 0.995, 0.995, 0.005);
    gPad->SetLeftMargin(0.2);
    gPad->SetRightMargin(0.005);
    gPad->SetTopMargin(0.002);
    gPad->SetBottomMargin(0.14);
    gPad->SetTickx(2);
    gPad->SetTicky(2);

    TMultiGraph* mg = new TMultiGraph();

    mg->Add(expGraph, "lp");
    mg->Add(litGraph, "lp");

    mg->Draw("alp");

    // X-axis parameters
    mg->GetXaxis()->SetTitle("Energy (MeV)");
    mg->GetXaxis()->SetTitleSize(0.05);
    mg->GetXaxis()->SetTitleFont(2);
    mg->GetXaxis()->SetTitleOffset(1.4);
    mg->GetXaxis()->CenterTitle();

    mg->GetXaxis()->SetLabelOffset(0.01);
    mg->GetXaxis()->SetLabelSize(0.05);
    mg->GetXaxis()->SetLabelFont(2);

    mg->GetXaxis()->SetNdivisions(10);
    mg->GetXaxis()->SetTickLength(0.03);

    // Y-axis parameters
    mg->GetYaxis()->SetTitle("(#frac{#sigma_{D} - #sigma_{H}}{#sigma_{D} + #sigma_{H}}) [%]");
    mg->GetYaxis()->SetTitleSize(0.06);
    mg->GetYaxis()->SetTitleFont(2);
    mg->GetYaxis()->SetTitleOffset(1.3);
    mg->GetYaxis()->CenterTitle();

    mg->GetYaxis()->SetLabelOffset(0.01);
    mg->GetYaxis()->SetLabelSize(0.05);

    mg->GetYaxis()->SetLabelFont(2);
    mg->GetYaxis()->SetNdivisions(10);
    mg->GetYaxis()->SetTickLength(0.02);

    gPad->SetLogx(1);
    
    //expGraph->GetYaxis()->SetRangeUser(-5,20);
    mg->GetXaxis()->SetLimits(2.5,600);

    //TLatex latex;
    //latex.SetNDC();
    //latex.SetTextSize(0.05);
    //latex.SetTextAlign(13); // align at top
    //latex.DrawLatex(0.47,0.52,"Ni");
    //latex.DrawLatex(0.32,0.4,"C");

    // Define legend format and contents
    TLegend *legend = new TLegend(0.55,0.25,0.95,0.4);
    legend->SetTextSize(0.03);
    legend->SetTextAlign(12);
    legend->AddEntry(expGraph,"Present work (D_{2}O, H_{2}O)","lp");
    legend->AddEntry(litGraph,"Abfalterer, 2001 (CH_{2}, C_{8}H_{18}, D_{2}O)","lp");
    legend->Draw();

    file->Close();
}
