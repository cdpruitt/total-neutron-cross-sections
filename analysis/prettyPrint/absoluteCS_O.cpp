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
    string expFileName = "/data2/analysis/total.root";
    string shiftedFileName = "/data2/analysis/shifted.root";

    string litFileName = "/data2/analysis/literatureData.root";

    TFile* expFile = new TFile(expFileName.c_str(),"READ");
    TFile* shiftedFile = new TFile(shiftedFileName.c_str(),"READ");
    TFile* litFile = new TFile(litFileName.c_str(),"READ");
    
    string expNatGraphName = "ONat";
    string exp18GraphName = "O18+1Barn";

    string litNatGraphName = "NatO(n,tot)";
    string lit18GraphName = "18O(n,tot)+1Barn";
    
    TGraphErrors* exp18Graph = (TGraphErrors*)shiftedFile->Get(exp18GraphName.c_str());
    TGraphErrors* expNatGraph = (TGraphErrors*)expFile->Get(expNatGraphName.c_str());

    TGraphErrors* litNatGraph = (TGraphErrors*)litFile->Get(litNatGraphName.c_str());
    TGraphErrors* lit18Graph = (TGraphErrors*)litFile->Get(lit18GraphName.c_str());

    // Set graph point and line characteristics
    exp18Graph->SetLineColor(kBlue);
    exp18Graph->SetLineWidth(5);
    exp18Graph->SetLineStyle(0);
    exp18Graph->SetMarkerColor(kBlue);

    expNatGraph->SetLineColor(kRed);
    expNatGraph->SetLineWidth(5);
    expNatGraph->SetLineStyle(0);
    expNatGraph->SetMarkerColor(kRed);

    litNatGraph->SetLineColor(kRed);
    litNatGraph->SetLineWidth(5);
    litNatGraph->SetLineStyle(2);
    litNatGraph->SetMarkerColor(kRed);

    lit18Graph->SetLineColor(kBlue-5);
    lit18Graph->SetLineWidth(3);
    lit18Graph->SetLineStyle(2);
    lit18Graph->SetMarkerColor(kBlue-5);

    // Pad dimensions and margins
    gPad->SetPad(0.005, 0.995, 0.995, 0.005);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.03);
    gPad->SetBottomMargin(0.15);
    gPad->SetTicky(2);

    // X-axis parameters
    exp18Graph->GetXaxis()->SetTitle("Energy (MeV)");
    exp18Graph->GetXaxis()->SetTitleSize(0.05);
    exp18Graph->GetXaxis()->SetTitleFont(2);
    exp18Graph->GetXaxis()->SetTitleOffset(1.4);
    exp18Graph->GetXaxis()->CenterTitle();

    exp18Graph->GetXaxis()->SetLabelOffset(0.01);
    exp18Graph->GetXaxis()->SetLabelSize(0.05);
    exp18Graph->GetXaxis()->SetLabelFont(2);

    exp18Graph->GetXaxis()->SetNdivisions(10);
    exp18Graph->GetXaxis()->SetTickLength(0.03);

    // Y-axis parameters
    exp18Graph->GetYaxis()->SetTitle("#sigma_{tot} (barns)");
    exp18Graph->GetYaxis()->SetTitleSize(0.06);
    exp18Graph->GetYaxis()->SetTitleFont(2);
    exp18Graph->GetYaxis()->SetTitleOffset(0.8);
    exp18Graph->GetYaxis()->CenterTitle();

    exp18Graph->GetYaxis()->SetLabelOffset(0.01);
    exp18Graph->GetYaxis()->SetLabelSize(0.05);

    exp18Graph->GetYaxis()->SetLabelFont(2);
    exp18Graph->GetYaxis()->SetNdivisions(10);
    exp18Graph->GetYaxis()->SetTickLength(0.02);

    exp18Graph->Draw("");
    expNatGraph->Draw("same");
    litNatGraph->Draw("same");
    lit18Graph->Draw("same");
    exp18Graph->Draw("same");
    expNatGraph->Draw("same");

    gPad->SetLogx(1);
    
    exp18Graph->GetYaxis()->SetRangeUser(0,4);

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.05);
    latex.SetTextAlign(13); // align at top
    latex.DrawLatex(0.37,0.80,"#color[4]{+1 barn}");
    //latex.DrawLatex(0.35,0.52,"Sn (elem.)");
    //latex.DrawLatex(0.32,0.4,"C (elem.)");

    // Define legend format and contents
    TLegend *legend = new TLegend(0.7,0.7,0.9,0.9);
    legend->AddEntry(litNatGraph,"{}^{nat}O (lit. data)","l");
    legend->AddEntry(lit18Graph,"{}^{18}O (lit. data)","l");
    legend->AddEntry(expNatGraph,"{}^{nat}O","l");
    legend->AddEntry(exp18Graph,"{}^{18}O","l");
    legend->Draw();

    expFile->Close();
    litFile->Close();
}
