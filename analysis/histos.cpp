using namespace std;

// ROOT file directory structure 
string dirs[4] = {"targetChanger","monitor","detS","scavenger"};

vector<TTree*> orchard; // holds all channel-specific trees
// so that fillHistos can loop through all the trees and make the same histos
// for each tree

vector<vector<TH1I*>> histos; // holds all histograms for the run
// split into sub-vectors on a per-channel basis

void setBranches(TTree* tree)
{
    tree->SetBranchAddress("macroNo",&macroNo);
    tree->SetBranchAddress("evtNo",&evtNo);
    tree->SetBranchAddress("macroTime",&macroTime);
    tree->SetBranchAddress("completeTime",&completeTime);
    tree->SetBranchAddress("targetPos",&targetPos);
    tree->SetBranchAddress("sgQ",&sgQ);
    tree->SetBranchAddress("lgQ",&lgQ);
    //tree->SetBranchAddress("waveform",&waveform);
}

void fillHistos(TTree* tree)
{

    setBranches(tree);

    int totalEntries = tree->GetEntries();

    gDirectory->cd("/");

    /*TH1I *blankRaw = new TH1I("blank","blank",20000,0,700);
    TH1I *carbonSRaw = new TH1I("carbonS","carbonS",20000,0,700);
    TH1I *carbonLRaw = new TH1I("carbonL","carbonL",20000,0,700);
    TH1I *Sn112Raw = new TH1I("Sn112","Sn112",20000,0,700);
    TH1I *NatSnRaw = new TH1I("NatSn","NatSn",20000,0,700);
    TH1I *Sn124Raw = new TH1I("Sn124","Sn124",20000,0,700);
    TH1I *totalRaw = new TH1I("total","total",20000,0,700);
    

    TH1I *blankRawLog = new TH1I("blankLog","blank",20000,0,TMath::Log10(700));
    TH1I *carbonSRawLog = new TH1I("carbonSLog","carbonS",20000,0,TMath::Log10(700));
    TH1I *carbonLRawLog = new TH1I("carbonLLog","carbonL",20000,0,TMath::Log10(700));
    TH1I *Sn112RawLog = new TH1I("Sn112Log","Sn112",20000,0,TMath::Log10(700));
    TH1I *NatSnRawLog = new TH1I("NatSnLog","NatSn",20000,0,TMath::Log10(700));
    TH1I *Sn124RawLog = new TH1I("Sn124Log","Sn124",20000,0,TMath::Log10(700));
    TH1I *totalRawLog = new TH1I("totalLog","total",20000,0,TMath::Log10(700));

    */
    TH1I *TOF = new TH1I("TTOF","All events time of flight",10000,-MICRO_PERIOD*1.1,MICRO_PERIOD*1.1);
    
    /*TH1I *firstSTOF = new TH1I("firstSTOF","first in micro time of flight",10000,-MICRO_PERIOD*1.1,MICRO_PERIOD*1.1);
    TH1I *secondSTOF = new TH1I("secondSTOF","second in micro time of flight",10000,-MICRO_PERIOD*1.1,MICRO_PERIOD*1.1);
    TH1I *thirdSTOF = new TH1I("thirdSTOF","third in micro time of flight",10000,-MICRO_PERIOD*1.1,MICRO_PERIOD*1.1);

    TH1I *firstSTOFblank = new TH1I("firstSTOFblank","first in micro time of flight, blank",10000,-MICRO_PERIOD*1.1,MICRO_PERIOD*1.1);
    TH1I *firstSTOFcs = new TH1I("firstSTOFcs","first in micro time of flight, carbon s",10000,-MICRO_PERIOD*1.1,MICRO_PERIOD*1.1);
    TH1I *firstSTOFcl = new TH1I("firstSTOFcl","first in micro time of flight, carbon l",10000,-MICRO_PERIOD*1.1,MICRO_PERIOD*1.1);
    TH1I *firstSTOFsn112 = new TH1I("firstSTOFsn112","first in micro time of flight, sn112",10000,-MICRO_PERIOD*1.1,MICRO_PERIOD*1.1);
    TH1I *firstSTOFsnnat = new TH1I("firstSTOFsnnat","first in micro time of flight, snnat",10000,-MICRO_PERIOD*1.1,MICRO_PERIOD*1.1);
    TH1I *firstSTOFsn124 = new TH1I("firstSTOFsn124","first in micro time of flight, sn124",10000,-MICRO_PERIOD*1.1,MICRO_PERIOD*1.1);

    TH1I *secondSTOFblank = new TH1I("secondSTOFblank","second in micro time of flight, blank",10000,-MICRO_PERIOD*1.1,MICRO_PERIOD*1.1);
    TH1I *secondSTOFcs = new TH1I("secondSTOFcs","second in micro time of flight, carbon s",10000,-MICRO_PERIOD*1.1,MICRO_PERIOD*1.1);
    TH1I *secondSTOFcl = new TH1I("secondSTOFcl","second in micro time of flight, carbon l",10000,-MICRO_PERIOD*1.1,MICRO_PERIOD*1.1);
    TH1I *secondSTOFsn112 = new TH1I("secondSTOFsn112","second in micro time of flight, sn112",10000,-MICRO_PERIOD*1.1,MICRO_PERIOD*1.1);
    TH1I *secondSTOFsnnat = new TH1I("secondSTOFsnnat","second in micro time of flight, snnat",10000,-MICRO_PERIOD*1.1,MICRO_PERIOD*1.1);
    TH1I *secondSTOFsn124 = new TH1I("secondSTOFsn124","second in micro time of flight, sn124",10000,-MICRO_PERIOD*1.1,MICRO_PERIOD*1.1);

    TH1I *thirdSTOFblank = new TH1I("thirdSTOFblank","third in micro time of flight, blank",10000,-MICRO_PERIOD*1.1,MICRO_PERIOD*1.1);
    TH1I *thirdSTOFcs = new TH1I("thirdSTOFcs","third in micro time of flight, carbon s",10000,-MICRO_PERIOD*1.1,MICRO_PERIOD*1.1);
    TH1I *thirdSTOFcl = new TH1I("thirdSTOFcl","third in micro time of flight, carbon l",10000,-MICRO_PERIOD*1.1,MICRO_PERIOD*1.1);
    TH1I *thirdSTOFsn112 = new TH1I("thirdSTOFsn112","third in micro time of flight, sn112",10000,-MICRO_PERIOD*1.1,MICRO_PERIOD*1.1);
    TH1I *thirdSTOFsnnat = new TH1I("thirdSTOFsnnat","third in micro time of flight, snnat",10000,-MICRO_PERIOD*1.1,MICRO_PERIOD*1.1);
    TH1I *thirdSTOFsn124 = new TH1I("thirdSTOFsn124","third in micro time of flight, sn124",10000,-MICRO_PERIOD*1.1,MICRO_PERIOD*1.1);
    */

    //TH2I *triangle = new TH2I("triangle","Pulse integral vs. TOF",500,0,MICRO_PERIOD+1,2000,0,65536);

    bool firstInMicro = false;
    bool secondInMicro = false;
    bool thirdInMicro = false;
    bool gammaInMicro = false;

    double trueTime = -1;
    int microNo = -1;
    int microNoPrev = -1;
    int microNo2Prev = -1;
    int microNo3Prev = -1;

    for(int j = 0; j<channelList.size(); j++)
    {
        // loop through tree once per channel number 
        cout << "Populating " << dirs[j] << " histograms..." << endl;
        nE = 0;

        gDirectory->cd("/");
        gDirectory->GetDirectory(dirs[j].c_str())->cd();

        DPPWaveformsDir = (TDirectory*)gDirectory->Get("DPPWaveformsDir");
        WaveWaveformsDir = (TDirectory*)gDirectory->Get("WaveWaveformsDir");

        int targetCounter = 0;
        targetTime = get<1>(targetTimeList[targetCounter]);

        tree->SetEntryList(channelList[j]);
        totalEntries = tree->GetEntries();

        rng = new TRandom3();

        double fullTimeP = 0;

        for (int i=0; i<totalEntries; i++)

        {
            tree->GetEntry(i);

            if (chNo == chNoList[j])
            {
                // prepare for filling basic histos
                fullTime = timetag+pow(2,32)*extTime;

                if (chNo == 4 || chNo == 6 || chNo == 7)
                {
                    fullTime += fineTime*(2./1024.);
                }

                if (evtType == 1)
                {

                    /*if (fullTime < targetTime-1000000)
                      {
                      cout << fullTime << " " << targetTime << endl;
                      }*/

                    while (fullTime-get<1>(targetTimeList[targetCounter+1])+TIME_OFFSET > 0)
                    {
                        // if it's been too long since the last target changer event,
                        // step to the next target changer event - provided
                        // we haven't reset the time because of a recent switch
                        // to waveform mode

                        if ((get<1>(targetTimeList[targetCounter]) < get<1>(targetTimeList[targetCounter+1])) || fullTimeP > fullTime)
                        {
                            targetCounter++;
                            targetTime = get<1>(targetTimeList[targetCounter]);
                            targetPos = get<2>(targetTimeList[targetCounter]);
                            targetType = get<3>(targetTimeList[targetCounter]);
                            fill(evtNo.begin(),evtNo.end(),0);

                            fill(waveformStart.begin(),waveformStart.end(),0); // prepare for next waveform mode

                            fullTimeP = fullTime; // update the time of the last event

                        }

                        else
                        {
                            break;
                        }
                    }

                    /*if (chNo==4 && targetCounter > 23388)
                      {
                      cout << "fullTime " << fullTime << " targetTime " << targetTime << " fineTime " << fineTime << "target counter " << targetCounter << endl;
                      abort();
                      }*/

                    // if event has associated target changer event, fill DPP histo
                    if (fullTime-targetTime+TIME_OFFSET < 650000 && fullTime-targetTime+TIME_OFFSET > 0) 
                    {
                        // within macropulse window; fill histos
                        outMacro = (TH1I*)(gDirectory->Get("outMacro"));
                        outMacro->Fill(get<0>(targetTimeList[targetCounter]));

                        outEvt = (TH1I*)(gDirectory->Get("outEvt"));
                        outEvt->Fill(evtNo[chNo]);

                        outExtTime = (TH1I*)(gDirectory->Get("outExtTime"));
                        outExtTime->Fill(extTime);

                        outTime = (TH1I*)(gDirectory->Get("outTime"));
                        outTime->Fill(timetag);

                        outSGQ = (TH1I*)(gDirectory->Get("outSGQ"));
                        outSGQ->Fill(sgQ);

                        outLGQ = (TH1I*)(gDirectory->Get("outLGQ"));
                        outLGQ->Fill(lgQ);

                        outFT = (TH1I*)(gDirectory->Get("outFT"));
                        outFT->Fill(fineTime + 16*rng->Rndm());

                        if (dummyWaveform->size() > 0)
                        {
                            DPPWaveformsDir->cd();

                            stringstream temp;
                            temp << "macroNo " << get<0>(targetTimeList[targetCounter]) << "evtNo " << evtNo[chNo];
                            DPPWaveform = new TH1I(temp.str().c_str(),temp.str().c_str(),dummyWaveform->size(),0,dummyWaveform->size()*2);

                            for(int i=0;i<dummyWaveform->size();i++)
                            {
                                DPPWaveform->SetBinContent(i,dummyWaveform->at(i));
                            }

                            gDirectory->cd("/");
                            gDirectory->GetDirectory(dirs[j].c_str())->cd();

                        }

                        if (chNo==2 || chNo==4 || chNo==6 || chNo==7)
                        {
                            microNo3Prev = microNo2Prev;
                            microNo2Prev = microNoPrev;
                            microNoPrev = microNo;
                            trueTime = fmod(((double)extTime*pow(2,32)+(double)timetag+((double)fineTime*2./1024.)-targetTime+TIME_OFFSET),MICRO_PERIOD);
                            microNo = floor(((double)extTime*pow(2,32)+(double)timetag+((double)fineTime*2./1024.)-targetTime+TIME_OFFSET)/MICRO_PERIOD);

                            firstInMicro = false;
                            secondInMicro = false;
                            thirdInMicro = false;
                            gammaInMicro = false;

                            if (microNo != microNoPrev)
                            {
                                firstInMicro = true;
                            }

                            if (microNo != microNo2Prev && microNo==microNoPrev)
                            {
                                secondInMicro = true;
                            }

                            if (microNo !=microNo3Prev && microNo==microNoPrev && microNoPrev==microNo2Prev)
                            {
                                thirdInMicro = true;
                            }

                            if (trueTime < 100)
                            {
                                gammaInMicro = true;
                            }

                            timeDiff << trueTime << " " << timetag << " " << targetTime << endl;

                            // convert trueTime into neutron velocity based on flight path distance
                            double velocity = pow(10.,7.)*FLIGHT_DISTANCE/trueTime; // in meters/sec 

                            // convert velocity to relativistic kinetic energy
                            double rKE = (pow((1.-pow((velocity/C),2.)),-0.5)-1.)*NEUTRON_MASS; // in MeV

                            if (trueTime > 100 && chNo != 2/* && !gammaInMicro*/) // gate disallowing gammas and monitors
                            {
                                switch (targetPos)
                                {
                                    case 1:
                                        // BLANK
                                        blankRaw->Fill(rKE);
                                        blankRawLog->Fill(TMath::Log10(rKE));
                                        break;

                                    case 2:
                                        // SHORT CARBON
                                        carbonSRaw->Fill(rKE);
                                        carbonSRawLog->Fill(TMath::Log10(rKE));
                                        break;

                                    case 3:
                                        // LONG CARBON
                                        carbonLRaw->Fill(rKE);
                                        carbonLRawLog->Fill(TMath::Log10(rKE));
                                        break;

                                    case 4:
                                        // Sn112
                                        Sn112Raw->Fill(rKE);
                                        Sn112RawLog->Fill(TMath::Log10(rKE));
                                        break;

                                    case 5:
                                        // Natural Sn
                                        NatSnRaw->Fill(rKE);
                                        NatSnRawLog->Fill(TMath::Log10(rKE));
                                        break;

                                    case 6:
                                        // Sn124
                                        Sn124Raw->Fill(rKE);
                                        Sn124RawLog->Fill(TMath::Log10(rKE));
                                        break;

                                    default:
                                        break;
                                }

                                totalRaw->Fill(rKE);
                                totalRawLog->Fill(TMath::Log10(rKE));
                            }

                            if (trueTime > 50 && chNo == 2) // gate disallowing gammas and non-monitors
                            {
                                switch (targetPos)
                                {
                                    case 1:
                                        // BLANK
                                        monBlankRaw->Fill(rKE);
                                        break;
                                    case 2:
                                        // SHORT CARBON
                                        monCarbonSRaw->Fill(rKE);
                                        break;
                                    case 3:
                                        // LONG CARBON
                                        monCarbonLRaw->Fill(rKE);
                                        break;
                                    case 4:
                                        // Sn112
                                        monSn112Raw->Fill(rKE);
                                        break;
                                    case 5:
                                        // Natural Sn
                                        monNatSnRaw->Fill(rKE);
                                        break;
                                    case 6:
                                        // Sn124
                                        monSn124Raw->Fill(rKE);
                                        break;
                                    default:
                                        break;
                                }
                                monTotalRaw->Fill(rKE);
                            }

                            if (targetPos > 0/* && gammaInMicro*/)
                            {
                                switch (chNo)
                                {
                                    case 2:
                                        MTOF->Fill(trueTime);
                                        break;
                                    case 4:
                                        STOF->Fill(trueTime);

                                        if (firstInMicro)
                                        {

                                            firstSTOF->Fill(trueTime);
                                            switch (targetPos)
                                            {
                                                case 1:
                                                    firstSTOFblank->Fill(trueTime);
                                                    break;
                                                case 2:
                                                    firstSTOFcs->Fill(trueTime);
                                                    break;
                                                case 3:
                                                    firstSTOFcl->Fill(trueTime);
                                                    break;
                                                case 4:
                                                    firstSTOFsn112->Fill(trueTime);
                                                    break;
                                                case 5:
                                                    firstSTOFsnnat->Fill(trueTime);
                                                    break;
                                                case 6:
                                                    firstSTOFsn124->Fill(trueTime);
                                                    break;
                                            }
                                        }

                                        else if (secondInMicro)
                                        {
                                            secondSTOF->Fill(trueTime);
                                            switch (targetPos)
                                            {
                                                case 1:
                                                    secondSTOFblank->Fill(trueTime);
                                                    break;
                                                case 2:
                                                    secondSTOFcs->Fill(trueTime);
                                                    break;
                                                case 3:
                                                    secondSTOFcl->Fill(trueTime);
                                                    break;
                                                case 4:
                                                    secondSTOFsn112->Fill(trueTime);
                                                    break;
                                                case 5:
                                                    secondSTOFsnnat->Fill(trueTime);
                                                    break;
                                                case 6:
                                                    secondSTOFsn124->Fill(trueTime);
                                                    break;
                                            }
                                        }

                                        else if (thirdInMicro)
                                        {
                                            thirdSTOF->Fill(trueTime);
                                            switch (targetPos)
                                            {
                                                case 1:
                                                    thirdSTOFblank->Fill(trueTime);
                                                    break;
                                                case 2:
                                                    thirdSTOFcs->Fill(trueTime);
                                                    break;
                                                case 3:
                                                    thirdSTOFcl->Fill(trueTime);
                                                    break;
                                                case 4:
                                                    thirdSTOFsn112->Fill(trueTime);
                                                    break;
                                                case 5:
                                                    thirdSTOFsnnat->Fill(trueTime);
                                                    break;
                                                case 6:
                                                    thirdSTOFsn124->Fill(trueTime);
                                                    break;
                                            }
                                        }

                                        triangle->Fill(trueTime,lgQ);
                                        break;

                                    case 6:
                                        LTOF->Fill(trueTime);
                                        break;
                                    case 7:
                                        RTOF->Fill(trueTime);
                                }
                            } 
                        }
                    }

                    prevTarget=1;
                }

                else if (evtType == 2)
                {
                    WaveWaveformsDir->cd();

                    TH1I* waveformHolder;

                    if (fullTime >= waveformStart[chNo]+650000 || prevTarget==1)
                    {
                        // new macropulse in waveform mode - create new plot

                        stringstream temp;
                        //cout << "waveform mode targetCounter" << targetCounter << " and supposed macropulse number " << get<0>(targetTimeList[targetCounter]) << endl;
                        temp << "macropulse " << get<0>(targetTimeList[targetCounter]) << " event no " << evtNo[chNo];
                        waveformStart[chNo] = fullTime;

                        switch (chNo)
                        {
                            case 0:
                                WaveWaveform0 = new TH1I(temp.str().c_str(),temp.str().c_str(),400000,0,800000);
                                waveformHolder = WaveWaveform0;
                                break;

                            case 2:
                                WaveWaveform2 = new TH1I(temp.str().c_str(),temp.str().c_str(),400000,0,800000);
                                waveformHolder = WaveWaveform2;
                                break;

                            case 4:
                                WaveWaveform4 = new TH1I(temp.str().c_str(),temp.str().c_str(),400000,0,800000);
                                waveformHolder = WaveWaveform4;
                                break;

                            case 6:
                                WaveWaveform6 = new TH1I(temp.str().c_str(),temp.str().c_str(),400000,0,800000);
                                waveformHolder = WaveWaveform6;
                                break;

                            case 7:
                                WaveWaveform7 = new TH1I(temp.str().c_str(),temp.str().c_str(),400000,0,800000);
                                waveformHolder = WaveWaveform7;
                                break;
                        }
                    }

                    for(int i=0;i<dummyWaveform->size();i++)
                    {
                        waveformHolder->SetBinContent(i+(fullTime-waveformStart[chNo])/2,dummyWaveform->at(i));
                    }

                    gDirectory->cd("/");
                    gDirectory->GetDirectory(dirs[j].c_str())->cd();

                    prevTarget=2;
                }

                nE++;
                evtNo[chNo]++;

                if(nE%100==0)
                {
                    cout << nE << " events\r";
                    fflush(stdout);
                }
            }
        }
        cout << endl;
    }
}

int main(int argc, char* argv[])
{
    // read in the raw tree name
    string runDir = argv[1];
    string runNo = argv[2];

    // Open the raw tree from the initial sort. If it doesn't exist, exit.
    TFile *file;
    TTree *ch0tree;
    TTree *ch2tree;
    TTree *ch4tree;
    TTree *ch6tree;
    TTree *ch0treeW;
    TTree *ch2treeW;
    TTree *ch4treeW;
    TTree *ch6treeW;

    stringstream treeName;
    stringstream fileName;

    treeName << runDir << "-" << runNo; 
    fileName << analysispath <<"analysis/" << runDir << "/" << treeName.str() << ".root";

    file = new TFile(fileName.str().c_str(),"UPDATE");

    ch0tree = (TTree*)file->Get("ch0tree");
    ch2tree = (TTree*)file->Get("ch2tree");
    ch4tree = (TTree*)file->Get("ch4tree");
    ch6tree = (TTree*)file->Get("ch6tree");
    ch0treeW = (TTree*)file->Get("ch0treeW");
    ch2treeW = (TTree*)file->Get("ch2treeW");
    ch4treeW = (TTree*)file->Get("ch4treeW");
    ch6treeW = (TTree*)file->Get("ch6treeW");

    orchard.push_back(ch0tree);
    orchard.push_back(ch2tree);
    orchard.push_back(ch4tree);
    orchard.push_back(ch6tree);
    orchard.push_back(ch0treeW);
    orchard.push_back(ch2treeW);
    orchard.push_back(ch4treeW);
    orchard.push_back(ch6treeW);

    // prepare the root file with 4 directories, one for each channel
    // these directories will hold basic variable histograms showing the
    // raw data in each tree, plus TOF, x-sections, etc histograms

    for(int i=0; i<4; i++)
        {
	  vector<TH1I*> tempVec;
	  histos.push_back(tempVec); // create sub-vector for this channel

	  gDirectory->mkdir(dirs[i].c_str(),dirs[i].c_str());
	  gDirectory->GetDirectory(dirs[i].c_str())->cd();

	  // instantiate histograms

	  histos.back().push_back(new TH1I("outMacro","outMacro",10000,0,1000000));
	  histos.back().back()->GetXaxis()->SetTitle("macropulse number");

	  histos.back().push_back(new TH1I("outEvt","outEvt",1000,0,1000));
	  histos.back().push_back(new TH1I("outExtTime","outExtTime",1000,0,1000));
	  histos.back().push_back(new TH1I("outTime","outTime",250000,0,2500000000));
	  histos.back().push_back(new TH1I("outSGQ","outSGQ",35000,0,35000));
	  histos.back().push_back(new TH1I("outLGQ","outLGQ",70000,0,70000));
	  histos.back().push_back(new TH1I("outFT","outFT",1023,0,1023));

	  gDirectory->mkdir("DPPWaveformsDir","raw DPP waveforms");
	  gDirectory->mkdir("WaveWaveformsDir","concatenated waveform waveforms");

	  gDirectory->cd("/");
        }

    for(int i = 0; i<orchard.size(); i++)
    {
        fillHistos(orchard[i]);
    }

    file->Write();

    file->Close();
}
