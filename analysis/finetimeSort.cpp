#include "TFile.h"
#include "TH2S.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

void main() {

    // a 2D plot to show fine time (FT) of events
    TH2S* outFTLR = new TH2S("outFTLR","outFTLR",1023,0,1023,1023,0,1023);
    outFTLR->SetMarkerStyle(20);

    TFile *events =
        TFile::Open("events.root");

    /*TTreeReader reader("tree", events);
    TTreeReader reader2("tree", events);

    TTreeReaderValue<Int_t> channel(reader,"event.chNo");
    TTreeReaderValue<Int_t> timetag(reader,"event.timetag");
    TTreeReaderValue<Int_t> finetime(reader, "event.finetime");

    TTreeReaderValue<Int_t> channel2(reader2,"event.chNo");
    TTreeReaderValue<Int_t> timetag2(reader2,"event.timetag");
    TTreeReaderValue<Int_t> finetime2(reader2, "event.finetime");
    */

    while (reader.Next()) {
        while (reader2.Next()) {
            if(timetag == timetag2) {
                if(channel == 6 && channel2 == 7) {
                    outFTLR->Fill(finetime,finetime2)
                }
                if(channel ==7 && channel2 == 6) {
                    outFTLR->Fill(finetime2,finetime)
                }
            }
        }
    }
    outFTLR->Draw();
}
