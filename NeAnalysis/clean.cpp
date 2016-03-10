int main(int argc, char* argv[])
{
    // read in data run location
    string runDir = argv[1];
    string runNo = argv[2];

    // Create a tree for this run
    TFile *file;

    stringstream treeName;
    stringstream fileName;
    treeName << runDir << "-" << runNo; 
    fileName << analysispath <<"analysis/" << runDir << "/" << treeName.str() << ".root";

    file = new TFile(fileName.str().c_str(),"UPDATE");

    if(file->Get("tree"))
    {
        tree = (TTree*)file->Get("tree");
        cout << "Found previously existing tree - skipping sorting " << treeName << endl;
    }

    else
    {
        tree = new TTree("tree","");
        cout << "Created ROOT tree " << treeName.str() << endl;

        tree->Branch("evtType",&ev.evtType,"evtType/i");
        tree->Branch("chNo",&ev.chNo,"chNo/i");
        tree->Branch("extTime",&ev.extTime,"extTime/i");
        tree->Branch("timetag",&ev.timetag,"timetag/d");
        tree->Branch("fineTime",&ev.fineTime,"fineTime/i");
        tree->Branch("sgQ",&ev.sgQ,"sgQ/i");
        tree->Branch("lgQ",&ev.lgQ,"lgQ/i");
        tree->Branch("waveform",&ev.waveform);
    }

    stringstream runName;
    runName << analysispath <<"output/" << runDir << "/data-" << runNo << ".evt";
    processRun(runName.str());

    file->Write();

    file->Close();

}
