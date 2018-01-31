{
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
    style->SetHistLineColor(kBlack);
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
    gPad->SetLeftMargin(0.1);
    gPad->SetRightMargin(0.1);
    gPad->SetTopMargin(0.03);
    gPad->SetBottomMargin(0.2);
    gPad->SetTicky(2);

    TMultiGraph* mg = new TMultiGraph();

    for(int i=0; i<=40; i++)
    {
        string fileName = "/data1/analysis/subRunCS/15_" + to_string(i) + ".root";

        TFile* file = new TFile(fileName.c_str(),"READ");
        if(!file->IsOpen())
        {
            continue;
        }

        string CGraphName = "CNatLitDifference";

        TGraphErrors* CGraph = (TGraphErrors*)file->Get(CGraphName.c_str());
        if(!CGraph)
        {
            continue;
        }

        mg->Add((TGraphErrors*)(CGraph->Clone()), "AP");

        file->Close();
    }

    mg->Draw("AP");

    gPad->SetLogx(1);
}
