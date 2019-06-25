void CanvasPartition(TCanvas *C, const Int_t Nx=2, const Int_t Ny=2,
        Float_t lMargin=0.15, Float_t rMargin=0.05,
        Float_t bMargin=0.15, Float_t tMargin=0.05);

void ThreePanelRelDiff() {
    TStyle * style = (TStyle*)gROOT->FindObject("graphStyle");

    if(!style)      
    {
        style = new TStyle("graphStyle","graphStyle");
    }

    TCanvas* canvas = new TCanvas("canvas","canvas", 2550, 450);
    //canvas->Divide(3, 2, 0.2, 0);

    int xPads = 3;
    int yPads = 1;

    double lMargin = 0.07;
    double rMargin = 0.01;
    double bMargin = 0.20;
    double tMargin = 0.01;

    CanvasPartition(canvas, xPads, yPads, lMargin, rMargin, bMargin, tMargin);

    style->SetOptStat(0);
    style->SetOptTitle(0);    
    //style->SetPalette(1,0);
    style->SetCanvasColor(10);      
    style->SetCanvasBorderMode(0);    
    style->SetFrameLineWidth(3);
    style->SetFrameFillColor(10);
    style->SetPadColor(10);
    style->SetHistLineWidth(4);
    style->SetHistLineColor(kRed);
    style->SetMarkerSize(0.6);
    style->SetMarkerStyle(8);
    //style->SetFuncWidth(3);
    //style->SetFuncColor(kRed);
    style->SetLabelColor(kBlack,"xyz");
    style->SetTitleSize(0.08,"xyz");
    style->SetTitleFillColor(10);
    style->SetTitleTextColor(kBlack);
    style->SetEndErrorSize(0);

    gROOT->SetStyle("graphStyle");
    gROOT->ForceStyle();

    // o16/o18
    {
        string fileName = "/data2/analysis/relative.root";
        string ramsauerFileName = "../../theory/ramsauer.root";

        TFile* file = new TFile(fileName.c_str(),"READ");
        TFile* ramsauerFile = new TFile(ramsauerFileName.c_str(), "READ");

        string relGraphName = "O18O16, percent";
        string relGraphSEName = "O18O16SysErrors, percent";

        string SARelDiffGraphThirdName = "RelDiff124_112Third";
        string SARelDiffGraphSixthName = "RelDiff124_112Sixth";

        string RamsauerRelDiffGraphName = "RelDiffRamsauer18_16";

        TGraphAsymmErrors* relGraph = (TGraphAsymmErrors*)file->Get(relGraphName.c_str());
        TGraphAsymmErrors* relGraphSE = (TGraphAsymmErrors*)file->Get(relGraphSEName.c_str());

        TGraph* SARelDiffGraphThird = (TGraph*)ramsauerFile->Get(SARelDiffGraphThirdName.c_str());
        TGraph* SARelDiffGraphSixth = (TGraph*)ramsauerFile->Get(SARelDiffGraphSixthName.c_str());
        TGraph* RamsauerRelDiffGraph = (TGraph*)ramsauerFile->Get(RamsauerRelDiffGraphName.c_str());

        // Set graph point and line characteristics
        relGraph->SetLineColor(kRed);
        relGraph->SetLineWidth(5);
        relGraph->SetLineStyle(0);
        relGraph->SetMarkerColor(kRed);
        relGraph->SetFillColor(kRed);
        relGraph->SetFillStyle(3002);

        SARelDiffGraphThird->SetLineStyle(9);
        SARelDiffGraphThird->SetLineWidth(3);
        SARelDiffGraphThird->SetLineColor(kBlack);

        SARelDiffGraphSixth->SetLineStyle(7);
        SARelDiffGraphSixth->SetLineWidth(3);
        SARelDiffGraphSixth->SetLineColor(kGray+2);

        RamsauerRelDiffGraph->SetLineStyle(7);
        RamsauerRelDiffGraph->SetLineWidth(5);
        RamsauerRelDiffGraph->SetLineColor(kGray+2);

        relGraphSE->SetFillColor(kBlue);
        relGraphSE->SetFillStyle(3001);

        // first panel
        {
            canvas->cd(0);

            TPad* pad = (TPad*)gROOT->FindObject("pad_0_0");
            pad->Draw();
            pad->cd();

            pad->SetFrameLineWidth(3);
            pad->SetTicky(1);
            pad->SetTickx(1);

            TMultiGraph* mg = new TMultiGraph();

            mg->Add(SARelDiffGraphThird, "l");
            mg->Add(SARelDiffGraphSixth, "l");
            mg->Add(relGraph,"3l");
            mg->Add(relGraphSE, "3");

            mg->Draw("al");

            // X-axis parameters
            mg->GetXaxis()->SetTitle("Energy [MeV]");
            mg->GetXaxis()->SetTitleSize(0.07);
            mg->GetXaxis()->SetTitleFont(2);
            mg->GetXaxis()->SetTitleOffset(1.4);
            mg->GetXaxis()->CenterTitle();

            //mg->GetXaxis()->SetLabelOffset(0.01);
            mg->GetXaxis()->SetLabelSize(0.07);
            mg->GetXaxis()->SetLabelFont(2);

            mg->GetXaxis()->SetNdivisions(10);
            mg->GetXaxis()->SetTickLength(0.03);

            // Y-axis parameters
            mg->GetYaxis()->SetTitle("#frac{#sigma_{H} - #sigma_{L}}{#sigma_{H} + #sigma_{L}} [%]");
            mg->GetYaxis()->SetTitleSize(0.07);
            mg->GetYaxis()->SetTitleFont(2);
            mg->GetYaxis()->SetTitleOffset(1.1);
            mg->GetYaxis()->CenterTitle();

            mg->GetYaxis()->SetLabelOffset(0.01);
            mg->GetYaxis()->SetLabelSize(0.07);

            mg->GetYaxis()->SetLabelFont(2);
            mg->GetYaxis()->SetNdivisions(8);
            mg->GetYaxis()->SetTickLength(0.02);

            mg->GetXaxis()->SetLimits(8,600);
            mg->GetYaxis()->SetRangeUser(0,4.9);

            gPad->SetLogx(1);

            TLatex latex;
            latex.SetNDC();
            latex.SetTextSize(0.10);
            latex.SetTextAlign(13); // align at top

            latex.DrawLatex(0.40, 0.35, "^{16,18}O");

            latex.SetTextSize(0.10);
            latex.DrawLatex(0.90, 0.33, "(a)");

            // Define legend format and contents
            TLegend *legend = new TLegend(0.32, 0.15, 0.57, 0.37);
            legend->SetNColumns(1);
            legend->AddEntry(relGraph,"Exp data, sys + stat   ","f");
            //legend->AddEntry(SARelDiffGraph,"SAS, r #alpha A^{1/3} ","l");
            legend->AddEntry(relGraphSE,"Exp data, sys only   ","f");
            //legend->AddEntry(RamsauerRelDiffGraph,"Ramsauer","l");
            //legend->Draw();
        }
    }

    // ni58/ni64
    {
        string fileName = "/data1/analysis/relative.root";
        string ramsauerFileName = "../../theory/ramsauer.root";

        TFile* file = new TFile(fileName.c_str(),"READ");
        TFile* ramsauerFile = new TFile(ramsauerFileName.c_str(), "READ");

        string relGraphName = "Ni64Ni58, percent";
        string relGraphSEName = "Ni64Ni58SysErrors, percent";

        string SARelDiffGraphThirdName = "RelDiff124_112Third";
        string SARelDiffGraphSixthName = "RelDiff124_112Sixth";

        string RamsauerRelDiffGraphName = "RelDiffRamsauer64_58";

        TGraphAsymmErrors* relGraph = (TGraphAsymmErrors*)file->Get(relGraphName.c_str());
        TGraphAsymmErrors* relGraphSE = (TGraphAsymmErrors*)file->Get(relGraphSEName.c_str());
        TGraph* SARelDiffGraphThird = (TGraph*)ramsauerFile->Get(SARelDiffGraphThirdName.c_str());
        TGraph* SARelDiffGraphSixth = (TGraph*)ramsauerFile->Get(SARelDiffGraphSixthName.c_str());
        TGraph* RamsauerRelDiffGraph = (TGraph*)ramsauerFile->Get(RamsauerRelDiffGraphName.c_str());

        // Set graph point and line characteristics
        relGraph->SetLineColor(kRed);
        relGraph->SetLineWidth(5);
        relGraph->SetLineStyle(0);
        relGraph->SetMarkerColor(kRed);
        relGraph->SetFillColor(kRed);
        relGraph->SetFillStyle(3002);

        SARelDiffGraphThird->SetLineStyle(9);
        SARelDiffGraphThird->SetLineWidth(3);
        SARelDiffGraphThird->SetLineColor(kBlack);

        SARelDiffGraphSixth->SetLineStyle(7);
        SARelDiffGraphSixth->SetLineWidth(3);
        SARelDiffGraphSixth->SetLineColor(kGray+2);

        RamsauerRelDiffGraph->SetLineStyle(7);
        RamsauerRelDiffGraph->SetLineWidth(5);
        RamsauerRelDiffGraph->SetLineColor(kGray+2);

        relGraphSE->SetFillColor(kBlue);
        relGraphSE->SetFillStyle(3001);

        // second panel
        {
            canvas->cd(0);

            TPad* pad = (TPad*)gROOT->FindObject("pad_1_0");
            pad->Draw();
            pad->cd();

            pad->SetFrameLineWidth(3);
            pad->SetTicky(1);
            pad->SetTickx(1);

            TMultiGraph* mg = new TMultiGraph();

            mg->Add(SARelDiffGraphThird, "l");
            mg->Add(SARelDiffGraphSixth, "l");
            mg->Add(relGraph,"3l");
            mg->Add(relGraphSE, "3");

            mg->Draw("al");

            // X-axis parameters
            mg->GetXaxis()->SetTitle("Energy [MeV]");
            mg->GetXaxis()->SetTitleSize(0.07);
            mg->GetXaxis()->SetTitleFont(2);
            mg->GetXaxis()->SetTitleOffset(1.4);
            mg->GetXaxis()->CenterTitle();

            mg->GetXaxis()->SetLabelOffset(0.01);
            mg->GetXaxis()->SetLabelSize(0.07);
            mg->GetXaxis()->SetLabelFont(2);

            mg->GetXaxis()->SetNdivisions(10);
            mg->GetXaxis()->SetTickLength(0.03);

            // Y-axis parameters
            mg->GetYaxis()->SetTitle("#frac{#sigma_{H} - #sigma_{L}}{#sigma_{H} + #sigma_{L}} [%]");
            mg->GetYaxis()->SetTitleSize(0.07);
            mg->GetYaxis()->SetTitleFont(2);
            mg->GetYaxis()->SetTitleOffset(1.1);
            mg->GetYaxis()->CenterTitle();

            mg->GetYaxis()->SetLabelOffset(0.01);
            mg->GetYaxis()->SetLabelSize(0.07);

            mg->GetYaxis()->SetLabelFont(2);
            mg->GetYaxis()->SetNdivisions(8);
            mg->GetYaxis()->SetTickLength(0.02);

            mg->GetXaxis()->SetLimits(8,600);
            mg->GetYaxis()->SetRangeUser(0,4.9);

            gPad->SetLogx(1);

            TLatex latex;
            latex.SetNDC();
            latex.SetTextSize(0.10);
            latex.SetTextAlign(13); // align at top

            latex.DrawLatex(0.20, 0.88, "^{58,64}Ni");

            latex.SetTextSize(0.10);
            latex.DrawLatex(0.88, 0.33, "(b)");

            // Define legend format and contents
            TLegend *legend = new TLegend(0.22, 0.15, 0.47, 0.33);
            legend->SetNColumns(1);
            legend->AddEntry(relGraph,"Exp data, sys + stat   ","f");
            //legend->AddEntry(SARelDiffGraph,"SAS, r #alpha A^{1/3} ","l");
            legend->AddEntry(relGraphSE,"Exp data, sys only   ","f");
            //legend->AddEntry(RamsauerRelDiffGraph,"Ramsauer","l");
            //legend->Draw();
        }
    }

    // sn112/sn124
    {
        string fileName = "/data2/analysis/relative.root";
        string ramsauerFileName = "../../theory/ramsauer.root";

        TFile* file = new TFile(fileName.c_str(),"READ");
        TFile* ramsauerFile = new TFile(ramsauerFileName.c_str(), "READ");

        string relGraphName = "Sn124Sn112, percent";
        string relGraphSEName = "Sn124Sn112SysErrors, percent";

        string SARelDiffGraphThirdName = "RelDiff124_112Third";
        string SARelDiffGraphSixthName = "RelDiff124_112Sixth";

        string RamsauerRelDiffGraphName = "RelDiffRamsauer124_112";
        string RamsauerRelDiffGraphSixthName = "RelDiffRamsauerSixth124_112";

        TGraphAsymmErrors* relGraph = (TGraphAsymmErrors*)file->Get(relGraphName.c_str());
        TGraphAsymmErrors* relGraphSE = (TGraphAsymmErrors*)file->Get(relGraphSEName.c_str());

        TGraph* SARelDiffGraphThird = (TGraph*)ramsauerFile->Get(SARelDiffGraphThirdName.c_str());
        TGraph* SARelDiffGraphSixth = (TGraph*)ramsauerFile->Get(SARelDiffGraphSixthName.c_str());
        TGraph* RamsauerRelDiffGraph = (TGraph*)ramsauerFile->Get(RamsauerRelDiffGraphName.c_str());
        TGraph* RamsauerRelDiffGraphSixth = (TGraph*)ramsauerFile->Get(RamsauerRelDiffGraphSixthName.c_str());

        // Set graph point and line characteristics
        relGraph->SetLineColor(kRed);
        relGraph->SetLineWidth(3);
        relGraph->SetLineStyle(0);
        relGraph->SetMarkerColor(kRed);
        relGraph->SetFillColor(kRed);
        relGraph->SetFillStyle(3002);

        SARelDiffGraphThird->SetLineStyle(9);
        SARelDiffGraphThird->SetLineWidth(3);
        SARelDiffGraphThird->SetLineColor(kBlack);

        SARelDiffGraphSixth->SetLineStyle(7);
        SARelDiffGraphSixth->SetLineWidth(3);
        SARelDiffGraphSixth->SetLineColor(kGray+2);

        RamsauerRelDiffGraph->SetLineStyle(7);
        RamsauerRelDiffGraph->SetLineWidth(3);
        RamsauerRelDiffGraph->SetLineColor(kGray+2);

        RamsauerRelDiffGraphSixth->SetLineStyle(7);
        RamsauerRelDiffGraphSixth->SetLineWidth(3);
        RamsauerRelDiffGraphSixth->SetLineColor(kGray+2);

        relGraphSE->SetFillColor(kBlue);
        relGraphSE->SetFillStyle(3001);

        // second panel
        {
            canvas->cd(0);

            TPad* pad = (TPad*)gROOT->FindObject("pad_2_0");
            pad->Draw();
            pad->cd();

            pad->SetFrameLineWidth(3);
            pad->SetTicky(1);
            pad->SetTickx(1);

            TMultiGraph* mg = new TMultiGraph();

            mg->Add(SARelDiffGraphThird, "l");
            mg->Add(SARelDiffGraphSixth, "l");
            mg->Add(relGraph,"3l");
            mg->Add(relGraphSE, "3");

            mg->Draw("al");

            // X-axis parameters
            mg->GetXaxis()->SetTitle("Energy [MeV]");
            mg->GetXaxis()->SetTitleSize(0.07);
            mg->GetXaxis()->SetTitleFont(2);
            mg->GetXaxis()->SetTitleOffset(1.4);
            mg->GetXaxis()->CenterTitle();

            mg->GetXaxis()->SetLabelOffset(0.01);
            mg->GetXaxis()->SetLabelSize(0.07);
            mg->GetXaxis()->SetLabelFont(2);

            mg->GetXaxis()->SetNdivisions(10);
            mg->GetXaxis()->SetTickLength(0.03);

            // Y-axis parameters
            mg->GetYaxis()->SetTitle("#frac{#sigma_{H} - #sigma_{L}}{#sigma_{H} + #sigma_{L}} [%]");
            mg->GetYaxis()->SetTitleSize(0.07);
            mg->GetYaxis()->SetTitleFont(2);
            mg->GetYaxis()->SetTitleOffset(1.1);
            mg->GetYaxis()->CenterTitle();

            mg->GetYaxis()->SetLabelOffset(0.01);
            mg->GetYaxis()->SetLabelSize(0.07);

            mg->GetYaxis()->SetLabelFont(2);
            mg->GetYaxis()->SetNdivisions(8);
            mg->GetYaxis()->SetTickLength(0.02);

            mg->GetXaxis()->SetLimits(8,600);
            mg->GetYaxis()->SetRangeUser(0,4.9);

            gPad->SetLogx(1);

            TLatex latex;
            latex.SetNDC();
            latex.SetTextSize(0.10);
            latex.SetTextAlign(13); // align at top

            latex.DrawLatex(0.13, 0.88, "^{112,124}Sn");

            latex.SetTextSize(0.10);
            latex.DrawLatex(0.87, 0.33, "(c)");

            // Define legend format and contents
            TLegend *legend = new TLegend(0.22, 0.15, 0.47, 0.33);
            legend->SetNColumns(1);
            legend->AddEntry(relGraph,"Exp data, sys + stat   ","f");
            //legend->AddEntry(SARelDiffGraphThird,"SAS, r #alpha A^{1/3} ","l");
            legend->AddEntry(relGraphSE,"Exp data, sys only   ","f");
            //legend->AddEntry(SARelDiffGraphSixth,"SAS, r #alpha A^{1/6} ","l");
            //legend->AddEntry(RamsauerRelDiffGraph,"Ramsauer","l");
            //legend->Draw();
        }
    }
}

void CanvasPartition(TCanvas *C,const Int_t Nx,const Int_t Ny,
                     Float_t lMargin, Float_t rMargin,
                     Float_t bMargin, Float_t tMargin)
{
   if (!C) return;

   // Setup Pad layout:
   Float_t vSpacing = 0.0;
   Float_t vStep  = (1.- bMargin - tMargin - (Ny-1) * vSpacing) / Ny;

   Float_t hSpacing = 0.0;
   Float_t hStep  = (1.- lMargin - rMargin - (Nx-1) * hSpacing) / Nx;

   Float_t vposd,vposu,vmard,vmaru,vfactor;
   Float_t hposl,hposr,hmarl,hmarr,hfactor;

   for (Int_t i=0;i<Nx;i++) {

      if (i==0) {
         hposl = 0.0;
         hposr = lMargin + hStep;
         hfactor = hposr-hposl;
         hmarl = lMargin / hfactor;
         hmarr = 0.0;
      } else if (i == Nx-1) {
         hposl = hposr + hSpacing;
         hposr = hposl + hStep + rMargin;
         hfactor = hposr-hposl;
         hmarl = 0.0;
         hmarr = rMargin / (hposr-hposl);
      } else {
         hposl = hposr + hSpacing;
         hposr = hposl + hStep;
         hfactor = hposr-hposl;
         hmarl = 0.0;
         hmarr = 0.0;
      }

      for (Int_t j=0;j<Ny;j++) {

         if (j==0) {
            vposd = 0.0;
            vposu = bMargin + vStep;
            vfactor = vposu-vposd;
            vmard = bMargin / vfactor;
            vmaru = 0.0;
         } else if (j == Ny-1) {
            vposd = vposu + vSpacing;
            vposu = vposd + vStep + tMargin;
            vfactor = vposu-vposd;
            vmard = 0.0;
            vmaru = tMargin / (vposu-vposd);
         } else {
            vposd = vposu + vSpacing;
            vposu = vposd + vStep;
            vfactor = vposu-vposd;
            vmard = 0.0;
            vmaru = 0.0;
         }

         C->cd(0);

         char name[16];
         sprintf(name,"pad_%i_%i",i,j);
         TPad *pad = (TPad*) gROOT->FindObject(name);
         if (pad) delete pad;
         pad = new TPad(name,"",hposl,vposd,hposr,vposu);
         pad->SetLeftMargin(hmarl);
         pad->SetRightMargin(hmarr);
         pad->SetBottomMargin(vmard);
         pad->SetTopMargin(vmaru);

         pad->SetFrameBorderMode(0);
         pad->SetBorderMode(0);
         pad->SetBorderSize(0);

         pad->Draw();
      }
   }
}
