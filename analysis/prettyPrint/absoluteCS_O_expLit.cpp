void absoluteCS_O_expLit()
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

    string expFileName = "/data3/analysis/total.root";
    string litFileName = "/data3/analysis/literatureData.root";

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
    
    string expNatOGraphName = "NatO";
    string litNatOGraphName = "NatO(n,tot)";
    //string exp18OGraphName = "O18";
    //string lit18OGraphName = "18O(n,tot)";
 
    TGraphErrors* expNatOGraph = (TGraphErrors*)expFile->Get(expNatOGraphName.c_str());
    TGraphErrors* litNatOGraph = (TGraphErrors*)litFile->Get(litNatOGraphName.c_str());
    //TGraphErrors* exp18OGraph = (TGraphErrors*)expFile->Get(exp18OGraphName.c_str());
    //TGraphErrors* lit18OGraph = (TGraphErrors*)litFile->Get(lit18OGraphName.c_str());

    if(!expNatOGraph)
    {
        cerr << "Error: failed to open " << expNatOGraphName << " in " << expFile << endl;
        exit(1);
    }

    /*if(!exp18OGraph)
    {
        cerr << "Error: failed to open " << exp18OGraphName << " in " << expFile << endl;
        exit(1);
    }*/

    // Set graph point and line characteristics
    expNatOGraph->SetLineColor(kRed);
    expNatOGraph->SetLineWidth(5);
    expNatOGraph->SetLineStyle(0);
    expNatOGraph->SetMarkerColor(kRed);

    litNatOGraph->SetLineColor(kBlack);
    litNatOGraph->SetLineWidth(5);
    litNatOGraph->SetLineStyle(2);
    litNatOGraph->SetMarkerColor(kBlack);

    //exp18OGraph->SetLineColor(kRed);
    //exp18OGraph->SetLineWidth(5);
    //exp18OGraph->SetLineStyle(0);
    //exp18OGraph->SetMarkerColor(kRed);

    //lit18OGraph->SetLineColor(kBlack);
    //lit18OGraph->SetLineWidth(5);
    //lit18OGraph->SetLineStyle(2);
    //lit18OGraph->SetMarkerColor(kBlack);

    // X-axis parameters
    expNatOGraph->GetXaxis()->SetTitle("Energy (MeV)");
    expNatOGraph->GetXaxis()->SetTitleSize(0.05);
    expNatOGraph->GetXaxis()->SetTitleFont(2);
    expNatOGraph->GetXaxis()->SetTitleOffset(1.4);
    expNatOGraph->GetXaxis()->CenterTitle();

    expNatOGraph->GetXaxis()->SetLabelOffset(0.01);
    expNatOGraph->GetXaxis()->SetLabelSize(0.05);
    expNatOGraph->GetXaxis()->SetLabelFont(2);

    expNatOGraph->GetXaxis()->SetNdivisions(10);
    expNatOGraph->GetXaxis()->SetTickLength(0.03);

    // Y-axis parameters
    expNatOGraph->GetYaxis()->SetTitle("#sigma_{tot} (barns)");
    expNatOGraph->GetYaxis()->SetTitleSize(0.06);
    expNatOGraph->GetYaxis()->SetTitleFont(2);
    expNatOGraph->GetYaxis()->SetTitleOffset(0.8);
    expNatOGraph->GetYaxis()->CenterTitle();

    expNatOGraph->GetYaxis()->SetLabelOffset(0.01);
    expNatOGraph->GetYaxis()->SetLabelSize(0.05);

    expNatOGraph->GetYaxis()->SetLabelFont(2);
    expNatOGraph->GetYaxis()->SetNdivisions(10);
    expNatOGraph->GetYaxis()->SetTickLength(0.02);

    expNatOGraph->Draw();
    //lit18OGraph->Draw("same");
    litNatOGraph->Draw("same");
    //exp18OGraph->Draw("same");
    expNatOGraph->Draw("same");

    gPad->SetLogx(1);
    
    expNatOGraph->GetYaxis()->SetRangeUser(0,4.5);

    //TLatex latex;
    //latex.SetNDC();
    //latex.SetTextSize(0.05);
    //latex.SetTextAlign(13); // align at top
    //latex.DrawLatex(0.65,0.65,"Pb");
    //latex.DrawLatex(0.47,0.52,"Sn");
    //latex.DrawLatex(0.32,0.4,"C");

    // Define legend format and contents
    TLegend *legend = new TLegend(0.6,0.7,0.9,0.9);
    //legend->SetHeader("data","C");
    legend->AddEntry(expNatOGraph,"H_{2}O^{exp.} - H_{2}^{lit.} ","l");
    legend->AddEntry(litNatOGraph,"O^{lit.} (derived) ","l");
    //legend->AddEntry(exp18OGraph,"{}^{18}O (exp.) ","l");
    //legend->AddEntry(lit18OGraph,"{}^{18}O (lit.) ","l");
    legend->Draw();
}
