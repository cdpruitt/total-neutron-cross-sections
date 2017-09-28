void absoluteCS_112Sn124Sn_localDOM()
{
    string Sn112FileName = "~/DOM/DOM48/sn/nsn112.root";
    string Sn124FileName = "~/DOM/DOM48/sn/nsn124.root";

    TFile* Sn112File = new TFile(Sn112FileName.c_str(),"READ");
    TFile* Sn124File = new TFile(Sn124FileName.c_str(),"READ");

    if(!Sn112File)
    {
        cerr << "Error: failed to open " << Sn112FileName << endl;
        exit(1);
    }

    if(!Sn124File)
    {
        cerr << "Error: failed to open " << Sn124FileName << endl;
        exit(1);
    }
    
    string expSn112GraphName = "react_nsn112";
    string expSn124GraphName = "react_nsn124";
    string theSn112GraphName = "react_nsn112";
    string theSn124GraphName = "react_nsn124";
 
    TGraphErrors* expSn112Graph = (TGraphErrors*)Sn112File->Get(expSn112GraphName.c_str());
    TGraphErrors* expSn124Graph = (TGraphErrors*)Sn124File->Get(expSn124GraphName.c_str());
    TGraphErrors* theSn112Graph = (TGraphErrors*)Sn112File->Get(theSn112GraphName.c_str());
    TGraphErrors* theSn124Graph = (TGraphErrors*)Sn124File->Get(theSn124GraphName.c_str());

    if(!expSn112Graph)
    {
        cerr << "Error: failed to open " << expSn112GraphName << " in " << Sn112File << endl;
        exit(1);
    }

    if(!expSn124Graph)
    {
        cerr << "Error: failed to open " << expSn124GraphName << " in " << Sn124File << endl;
        exit(1);
    }

    TStyle * style = (TStyle*)gROOT->FindObject("graphStyle");

    if(!style)      
    {
        style = new TStyle("graphStyle","graphStyle");
    }

    exit(0);

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

    exit(0);

    // Set graph point and line characteristics
    expSn112Graph->SetLineColor(kRed);
    expSn112Graph->SetLineWidth(5);
    expSn112Graph->SetLineStyle(0);
    expSn112Graph->SetMarkerColor(kRed);

    expSn124Graph->SetLineColor(kBlack);
    expSn124Graph->SetLineWidth(3);
    expSn124Graph->SetLineStyle(2);
    expSn124Graph->SetMarkerColor(kBlack);

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
    //expSn124Graph->Draw("same");

    gPad->SetLogx(1);
    
    expSn112Graph->GetYaxis()->SetRangeUser(0,9);

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.05);
    latex.SetTextAlign(13); // align at top
    latex.DrawLatex(0.65,0.65,"Pb");
    latex.DrawLatex(0.47,0.52,"Sn");
    latex.DrawLatex(0.32,0.4,"C");

    // Define legend format and contents
    //TLegend *legend = new TLegend(0.75,0.7,0.9,0.9);
    //legend->SetHeader("data","C");
    //legend->AddEntry(graph,"{}^{nat}Pb","l");
    //legend->Draw();

    Sn112File->Close();
    Sn124File->Close();
}
