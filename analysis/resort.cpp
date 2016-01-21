// Takes a ROOT tree with raw event data and processes it into a new tree useful
// for analysis

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include "TFile.h"
#include "TTree.h"
#include "TEntryList.h"

using namespace std;

// Experimental constants


// Time delay of target changer coarse time after real macropulse start time
const double TIME_OFFSET = 836; // in ns

// Period of micropulses
const double MICRO_PERIOD = 1788.820; // in ns

// Period of macropulses
const double MACRO_PERIOD = 8330000; // in ns

const string analysispath =  "/media/Drive3/";


// variables for holding raw tree event data
unsigned int chNo, evtType, extTime, fineTime, sgQ, lgQ;
double timetag;

// additional variables for holding new tree event data
unsigned int evtNo, macroNo, microNo, targetPos;
double completeTime, macroTime, microTime;

vector<int> waveform; // for holding one event's waveform data
vector<int> *dummyWaveform; // for transferring waveform data from clean tree to subtrees

// event structure for resorted trees (different from raw tree event structure)
struct event
{
  unsigned int macroNo; // label each event by macropulse
  unsigned int microNo; // label each event by micropulse
  unsigned int evtNo; // uniquely label each event in a macropulse

  double macroTime; // =extTime+timetag for target changer
  double completeTime; // =extTime+timetag+fineTime for detector events
  double microTime; // =completeTime-macroTime+TIME_OFFSET for detector events

  unsigned int targetPos; // target changer position;

  unsigned int sgQ, lgQ; // event charge gates

  vector<int> waveform; // transfer waveforms from raw to cleaned tree
} ev;

TTree* cTree; // holds raw data with beam-off periods removed

// channel-specific trees for DPP and waveform mode
// data from cTree will be reprocessed and divided into these trees for later
// analysis/histogramming
TTree* ch0Tree;
TTree* ch2Tree;
TTree* ch4Tree;
TTree* ch6Tree;
TTree* ch0TreeW;
TTree* ch2TreeW;
TTree* ch4TreeW;
TTree* ch6TreeW;


// Target Changer lgQ gates for setting integral target positions 1-6
const int tarGate[12] = {5000,7500,11000,13500,17000,20000,23000,27000,29000,33000,35000,39000};

void fillTree(TTree* tree)
{

    ev.macroNo = macroNo;
    ev.microNo = microNo;
    ev.evtNo = evtNo;
    ev.macroTime = macroTime;
    ev.completeTime = completeTime;
    ev.microTime = microTime;
    ev.targetPos = targetPos;
    ev.sgQ = sgQ;
    ev.lgQ = lgQ;

    if ((int)completeTime%10000 != 0 && evtType==1)
    {
        waveform.clear(); 
    }

    ev.waveform = waveform; 

    tree->Fill();

    // reset waveform in preparation for next event
    waveform.clear();
}

void cleanTree(TTree* tree)
{
    // This function excises periods of beam off from the raw tree
    // The approach is:
    //  - to loop through all entries in the tree in the order that
    //    they were read out of the digitizer
    //  - when a target changer event is found that doesn't line up
    //    with the macropulse structure (i.e., timestamps 8.3 ms or
    //    16.6 ms after each other), this indicates that beam went
    //    off in the previous macropulse.
    //  - Move backward in the tree until the previous chunk of
    //    target changer events is found, and save the event indices
    //    of the beginning and end of this dead time
    //  - Using this info, fill a vector with the periods when the
    //    beam is on, and recreate the raw tree using these periods
    //    as gates for when data should be allowed. The resulting
    //    tree, "cTree", is used for all analysis.

    // link raw tree branches to variables-to-be-read
    tree->SetBranchAddress("chNo",&chNo);
    tree->SetBranchAddress("evtType",&evtType);
    tree->SetBranchAddress("timetag",&timetag);
    tree->SetBranchAddress("extTime",&extTime);

    cTree = tree->CloneTree(); // COMMENT to start using clean tree method again

    /*
    // prepare variables for evaluating whether a macropulse has been skipped
    // (thus indicating that beam went off)
    double timeDiff = 0; // difference between times of adjacent macropulses
    double prevMacroTime = 0; // placeholder for previous macropulse time
    double currentMacroTime = 0; // holds complete timestamp of current macro
    int init = 0; // for holding the first event index of a 'beam on' period
    int prevEvtType = 2; // keeps track of when acquisition switches modes;
                         // initialized to 2 in order to accept the first
                         // macropulse of each new acquisition pariod
    vector<pair<int,int>> beamOn; // holds event indices for when beam is on
                         // first index is 'start of beam on', second index
                         // is 'end of beam on'

    // loop through raw tree looking for target changer events
    int totalEntries = tree->GetEntries();

    for(int i = 0; i<totalEntries; i++)
    {
        tree->GetEntry(i);

        if (chNo == 0)
        {
            // found a target changer event - now test for beam on/off
            currentMacroTime = ((double)extTime*pow(2,32)+timetag);
            timeDiff = currentMacroTime-prevMacroTime;
            cout << "macropulse time = " << currentMacroTime << "timeDiff = " << timeDiff << "\r";
            fflush(stdout);

            if (fmod(timeDiff,MACRO_PERIOD) > 10000 || evtType == 2 || prevEvtType == 2)
            {
                // the target changer event came within the expected window of
                // 8.3 or 16.6 ms, or beam at 80 Hz (24.9 ms gaps),
                // or was the start of a new acquisition period
                // ...beam is ON so accept this event and continue
            }

            else
            {
                // target changer time was OUTSIDE acceptable bounds
                // thus beam was off for some period, so figure out what that period is

                // we'll need to go backwards to find the previous set of
                // target changer events when beam was still good
                // shift the event index backwards by 10 so we can escape the current
                // chunk of target changer events where beam off was detected
                for (int j = i-10; j>init; j--)
                {
                    //cout << "j = " << j << ", chNo = " << chNo << "\r";
                    //fflush(stdout);

                    tree->GetEntry(j);
                    if(chNo==0)
                    {
                        // found the previous chunk of target changer events;
                        // record their index in beamOn
                        beamOn.push_back(make_pair(init,j));

                        // shift init to the next event in preparation for next
                        // beam on period
                        init = i+1;

                        //cout << "beam off period" << endl;
                        
                        // ... and keep looping forward through raw
                        // tree where we left off (that is, at index i)
                        break;
                    }
                }
            }

            // reset placeholder variables for next event
            prevEvtType = evtType;
            prevMacroTime = currentMacroTime;
        }

        if(i%100==0)
        {
            cout << "Sorting through raw tree, event #" << i << "\r";
            fflush(stdout);
        }
    }

    cout << endl << endl;

    // after loop finishes, add the last leg of the run to the beamOn list
    beamOn.push_back(make_pair(init,totalEntries));

    // create a new empty tree to hold the cleaned data from the raw tree
    
    //UNCOMMENT to start using clean tree method again
    cTree = tree->CloneTree(0);

    // Add only events when beam was on to the cleaned tree
    for(int i = 0; i<beamOn.size(); i++)
    {
        int resume = beamOn[i].first;
        int omit = beamOn[i].second;

        cout << "Filling clean tree with events " << beamOn[i].first << " to " << beamOn[i].second << endl;

        for(int j = resume; j<omit+1; j++)
        {
            tree->GetEntry(j);
            cTree->Fill();

            if(j%100==0)
            {
                cout << "Populating cleaned tree, event #" << j << "\r";
                fflush(stdout);
            }
        }

        //cout << endl << "Beam off starting at event " << omit << endl;

    }

    cout << endl << beamOn.size()-1 << " periods of beam anomaly removed." << endl;
    */
}

void populateTrees()
{
    // This method creates two trees for each channel, one for DPP mode and one
    // for waveform mode (except for ch6). Each event to be added to these trees
    // will be assigned to a macropulse for easier analysis later.

    cTree->SetBranchAddress("chNo",&chNo);
    cTree->SetBranchAddress("evtType",&evtType);
    cTree->SetBranchAddress("timetag",&timetag);
    cTree->SetBranchAddress("extTime",&extTime);
    cTree->SetBranchAddress("sgQ",&sgQ);
    cTree->SetBranchAddress("lgQ",&lgQ);
    cTree->SetBranchAddress("fineTime",&fineTime);
    cTree->SetBranchAddress("waveform",&dummyWaveform);

    for(int j = 0; j<8; j=j+2)
    {
        // looping through tree once per channel number
        cout << "Populating trees, ch = " << j << endl;

        //cTree->SetEntryList(channelList[j]);
        int totalEntries = cTree->GetEntries();

        int prevEvtType = 0; // keeps track of when acquisition switches modes
        // so that macropulse counting is done correctly

        int prevTime = 0;         // keep track of previous macro time 

        switch (j)
        {
            case 0:
                // Target changer (ch0)
                // Populate a tree of only target changer events in preparation
                // for assigning times and target changer positions to the
                // detector channels (ch2, ch4, ch6)

                macroNo = 0;

                for (int i=0; i<totalEntries; i++)
                {
                    cTree->GetEntry(i);

                    if (chNo == j)
                    {
                        macroTime = (double)extTime*pow(2,32)+timetag;
                        completeTime = macroTime;

                        if (evtType==1)
                        {
                            // new macropulse in DPP mode

                            // assign an integral target position based on lgQ value
                            if (lgQ>tarGate[0] && lgQ<tarGate[1])
                            {
                                targetPos = 1;
                            }

                            else if (lgQ>tarGate[2] && lgQ<tarGate[3])
                            {
                                targetPos = 2;
                            }

                            else if (lgQ>tarGate[4] && lgQ<tarGate[5])
                            {
                                targetPos = 3;
                            }

                            else if (lgQ>tarGate[6] && lgQ<tarGate[7])
                            {
                                targetPos = 4;
                            }

                            else if (lgQ>tarGate[8] && lgQ<tarGate[9])
                            {
                                targetPos = 5;
                            }

                            else if (lgQ>tarGate[10] && lgQ<tarGate[11])
                            {
                                targetPos = 6;
                            }

                            else
                            {
                                targetPos = 0;
                            }

                            for(int k = 0; k<dummyWaveform->size(); k++)
                            {
                                waveform.push_back(dummyWaveform->at(k));
                            }

                            if(prevEvtType==2)
                            {
                                // start of a new DPP period - indicate by
                                // making evtNo 1 so we'll know when we sort
                                // through this tree to assign macropulse times
                                // to each event
                                evtNo = 1;
                            }

                            else
                            {
                                evtNo = 0;
                            }

                            fillTree(ch0Tree);

                            //cout << macroTime << endl;

                            macroNo++; // increment macropulse counter
                        }

                        else if (evtType==2)
                        {
                            // waveform mode
                            // don't increment the macropulse number, but fill
                            // a waveform-only tree with timestamp, evtNo, and waveform data
                            targetPos = 0;

                            if (prevEvtType==1)
                            {
                                // the previous event was a DPP event, so this is
                                // the first waveform mode event in the wavelet

                                prevTime = macroTime;
                            }

                            else if (macroTime>prevTime+7000000)
                            {
                                evtNo = 0; // clear the evtNo counter
                            }

                            // transfer the waveform over from the old tree to
                            // where the new tree points
                            for(int k = 0; k<dummyWaveform->size(); k++)
                            {
                                waveform.push_back(dummyWaveform->at(k));
                            }

                            fillTree(ch0TreeW);
                        }

                        prevEvtType = evtType;

                        if (i%100==0)
                        {
                            cout << "Evt number " << i << ", macroTime " << macroTime << "\r";
                            fflush(stdout);
                        }
                    }
                }

                cout << endl << endl;

                // link tree branches to macropulse variables for filling
                // ch2, ch4, ch6 trees
                ch0Tree->SetBranchAddress("macroNo",&macroNo);
                ch0Tree->SetBranchAddress("evtNo",&evtNo);
                ch0Tree->SetBranchAddress("macroTime",&macroTime);
                ch0Tree->SetBranchAddress("targetPos",&targetPos);

                break;

            case 2:
                {
                // Monitor (ch2)
                // Populate a tree of only monitor events

                // prepare the currentMacroTime and currentMacroNo variables
                // for keeping track of the macropulse we're using to assign
                // macroTime and macroNo to events
                //ch0Tree->GetEntry(0);

                double prevCompleteTime = 0;

                evtNo = 0; // first macropulse; reset event counter

                // move the ch0Tree forward one step to prepare for moving to
                // the next macropulse
                ch0Tree->GetEntry(0);

                // the macropulse tree (ch0Tree) is now pointing at a future
                // event; once monitor events 'catch up' to the future
                // macropulse, we'll shift to that macropulses' data for
                // assigning to monitor events
                // start looping through monitor events
                for (int i=0; i<totalEntries; i++)
                {
                    cTree->GetEntry(i);

                    if (chNo == j)
                    {
                        completeTime = (double)extTime*pow(2,32)+timetag;

                        if (evtType == 1)
                        {
                            // DPP mode

                            // first, check to see if we've wrapped around
                            // to a new DPP mode period
                            if (completeTime<prevCompleteTime)
                            {
                                // new DPP mode period
                                // update the ch0Tree to point at the new
                                // macropulse
                                double prevMacroTime = macroTime;

                                cout << "new DPP mode period (from channel wraparound)" << endl;

                                while (prevMacroTime < macroTime)
                                {
                                    prevMacroTime = macroTime;
                                    ch0Tree->GetEntry(macroNo+1);
                                }
                                // ch0Tree is now pointing at the first
                                // macropulse of the next DPP mode period

                                // reset the event counter for this channel
                                evtNo = 0;
                            }

                            // check to see if enough time has elapsed such that
                            // we should update the current macropulse to a new
                            // one
                            else if (completeTime-macroTime-TIME_OFFSET > 650000)
                            {
                                // macropulse expired
                                // time to examine the next macropulse

                                // hold the old macro time so we'll be able to
                                // compare it with the new macro time and make
                                // sure there's still normal beam structure
                                double prevMacroTime = macroTime;

                                // pull the next macropulse 
                                ch0Tree->GetEntry(macroNo+1);

                                // first, we should check to make sure we're not
                                // about to run off the end of ch0Tree and get
                                // to the end of the run
                                if(macroNo+1>=ch0Tree->GetEntries())
                                {
                                    // reached the end of the ch0Tree
                                    // throw away all subsequent events because
                                    // it's the end of the run
                                    break;
                                }

                                // we're not at the end of the run, so we
                                // should examine the next macropulse to see
                                // what to do with it

                                // evtNo was updated when we pulled the next
                                // macropulse - let's check to see if we've
                                // reached the end of a DPP period
                                if (evtNo==1)
                                {
                                    // new DPP period found because the
                                    // macropulse has evtNo==1, so move the
                                    // cleaned tree indices forward to the new
                                    // DPP period

                                    cout << "new DPP mode period (from ch0Tree wraparound)" << endl;

                                    do
                                    {
                                        cout << "prevCompleteTime = " << prevCompleteTime << " completeTime = " << completeTime << endl;
                                        prevCompleteTime = completeTime;
                                        do
                                        {
                                            i++;
                                            cTree->GetEntry(i);
                                            completeTime = (double)extTime*pow(2,32)+timetag;
                                        }
                                        while (chNo != j);
                                    }
                                    while (prevCompleteTime < completeTime);

                                    cout << "completeTime = " << completeTime << " macroTime = " << macroTime << endl;
                                    //exit(0);
                                }

                                // not a new DPP period, so lastly let's check
                                // to see if there are any beam anomalies (that
                                // is, the new macropulse is out-of-sync with
                                // the old one)
                                else if ((macroTime-prevMacroTime > 8300000 && macroTime-prevMacroTime < 8360000) || (macroTime-prevMacroTime > 16600000 && macroTime-prevMacroTime < 17200000) || (macroTime-prevMacroTime > 24900000 && macroTime-prevMacroTime < 25080000))
                                {
                                    // no beam anomaly - keep filling the
                                    // channel-specific tree and continue as
                                    // normal
                                }

                                else
                                {
                                    // there IS a beam anomaly - cycle through
                                    // the cleaned and ch0Trees until we find
                                    // another period of clean beam
                                    while(!(macroTime-prevMacroTime > 8300000 && macroTime-prevMacroTime < 8360000) && !(macroTime-prevMacroTime > 16600000 && macroTime-prevMacroTime < 17200000) && !(macroTime-prevMacroTime > 24900000 && macroTime-prevMacroTime < 25080000))
                                    {
                                        // update the ch0Tree to point at the
                                        // next macropulse so we can check to
                                        // see if it's out-of-sync with the
                                        // previous macropulse 
                                        prevMacroTime = macroTime;
                                        ch0Tree->GetEntry(macroNo+1);

                                        cout << "macroTime is " << macroTime << " prevMacroTime is " << prevMacroTime << endl;
                                    };

                                    // OK - ch0Tree is now in a period of
                                    // good beam. We need to move the cleaned
                                    // tree index forward to point at this new
                                    // area of good beam

                                    while (completeTime < macroTime-TIME_OFFSET)
                                    {
                                        do
                                        {
                                            i++;
                                            cTree->GetEntry(i);
                                            completeTime = (double)extTime*pow(2,32)+timetag;
                                        }
                                        while (chNo != j);

                                        cout << "completeTime is " << completeTime << " macroTime is " << macroTime << " i is " << i << "\r";
                                        fflush(stdout);
                                    }

                                    cout << endl;

                                    // now the cleaned tree index is pointing
                                    // to the area of good beam, so continue
                                    // filling the channel-specific tree
                                }

                                // lastly, because we shifted to a new
                                // macropulse at the start of this conditional,
                                // we should reset the event counter for this
                                // channel
                                evtNo = 0;
                            }

                            // pointing at the correct macropulse - fill the
                            // channel-specific tree with data from this channel
                        
                            // calculate the time elapsed since the start of the
                            // current macropulse and the time elapsed since the start
                            // of the micropulse
                            double timeDiff = completeTime-macroTime-TIME_OFFSET;

                            microTime = fmod(timeDiff,MICRO_PERIOD);
                            microNo = floor(timeDiff/MICRO_PERIOD);

                            for(int k = 0; k<dummyWaveform->size(); k++)
                            {
                                waveform.push_back(dummyWaveform->at(k));
                            }

                            fillTree(ch2Tree);
                        }

                        else
                        {
                            // waveform mode data: fill ch2TreeW
                            
                            for(int k = 0; k<dummyWaveform->size(); k++)
                            {
                                waveform.push_back(dummyWaveform->at(k));
                            }

                            // populate waveform mode tree for monitor
                            fillTree(ch2TreeW);

                            if (prevEvtType==1)
                            {
                                ch0Tree->GetEntry(macroNo+1); // shift to the next
                                // DPP-mode macropulse in preparation for the next DPP mode
                                evtNo=0;
                            }
                        }

                        evtNo++; // increment event counter for channel 2
                        prevEvtType = evtType;

                        if (i%100==0)
                        {
                            cout << "Evt number " << i << ", completeTime " << completeTime << "\r";
                            fflush(stdout);
                        }
                    }
                    // prepare for next iteration of loop
                    prevCompleteTime = completeTime;
                }

                cout << endl << endl;
                break;
                }

            case 4: // summed detector (ch4)
            case 6: // scavenger (ch6)
                {
                // Populate a tree of only monitor events

                // prepare the currentMacroTime and currentMacroNo variables
                // for keeping track of the macropulse we're using to assign
                // macroTime and macroNo to events
                //ch0Tree->GetEntry(0);

                double prevCompleteTime = 0;

                evtNo = 0; // first macropulse; reset event counter

                // move the ch0Tree forward one step to prepare for moving to
                // the next macropulse
                ch0Tree->GetEntry(0);

                // the macropulse tree (ch0Tree) is now pointing at a future
                // event; once monitor events 'catch up' to the future
                // macropulse, we'll shift to that macropulses' data for
                // assigning to monitor events
                // start looping through monitor events
                for (int i=0; i<totalEntries; i++)
                {
                    cTree->GetEntry(i);

                    if (chNo == j)
                    {
                        completeTime = (double)extTime*pow(2,32)+timetag;

                        if (evtType == 1)
                        {
                            // DPP mode

                            // first, check to see if we've wrapped around
                            // to a new DPP mode period
                            if (completeTime<prevCompleteTime)
                            {
                                // new DPP mode period
                                // update the ch0Tree to point at the new
                                // macropulse
                                double prevMacroTime = macroTime;

                                cout << "new DPP mode period (from channel wraparound)" << endl;

                                while (prevMacroTime < macroTime)
                                {
                                    prevMacroTime = macroTime;
                                    ch0Tree->GetEntry(macroNo+1);
                                }
                                // ch0Tree is now pointing at the first
                                // macropulse of the next DPP mode period

                                // reset the event counter for this channel
                                evtNo = 0;
                            }

                            // check to see if enough time has elapsed such that
                            // we should update the current macropulse to a new
                            // one
                            else if (completeTime-macroTime-TIME_OFFSET > 650000)
                            {
                                // macropulse expired
                                // time to examine the next macropulse

                                // hold the old macro time so we'll be able to
                                // compare it with the new macro time and make
                                // sure there's still normal beam structure
                                double prevMacroTime = macroTime;

                                // pull the next macropulse 
                                ch0Tree->GetEntry(macroNo+1);

                                // first, we should check to make sure we're not
                                // about to run off the end of ch0Tree and get
                                // to the end of the run
                                if(macroNo+1>=ch0Tree->GetEntries())
                                {
                                    // reached the end of the ch0Tree
                                    // throw away all subsequent events because
                                    // it's the end of the run
                                    break;
                                }

                                // we're not at the end of the run, so we
                                // should examine the next macropulse to see
                                // what to do with it

                                // evtNo was updated when we pulled the next
                                // macropulse - let's check to see if we've
                                // reached the end of a DPP period
                                if (evtNo==1)
                                {
                                    // new DPP period found because the
                                    // macropulse has evtNo==1, so move the
                                    // cleaned tree indices forward to the new
                                    // DPP period

                                    cout << "new DPP mode period (from ch0Tree wraparound)" << endl;

                                    do
                                    {
                                        cout << "prevCompleteTime = " << prevCompleteTime << " completeTime = " << completeTime << endl;
                                        prevCompleteTime = completeTime;
                                        do
                                        {
                                            i++;
                                            cTree->GetEntry(i);
                                            completeTime = (double)extTime*pow(2,32)+timetag;
                                        }
                                        while (chNo != j);
                                    }
                                    while (prevCompleteTime < completeTime);

                                    cout << "completeTime = " << completeTime << " macroTime = " << macroTime << endl;
                                    //exit(0);
                                }

                                // not a new DPP period, so lastly let's check
                                // to see if there are any beam anomalies (that
                                // is, the new macropulse is out-of-sync with
                                // the old one)
                                else if ((macroTime-prevMacroTime > 8300000 && macroTime-prevMacroTime < 8360000) || (macroTime-prevMacroTime > 16600000 && macroTime-prevMacroTime < 17200000) || (macroTime-prevMacroTime > 24900000 && macroTime-prevMacroTime < 25080000))
                                {
                                    // no beam anomaly - keep filling the
                                    // channel-specific tree and continue as
                                    // normal
                                }

                                else
                                {
                                    // there IS a beam anomaly - cycle through
                                    // the cleaned and ch0Trees until we find
                                    // another period of clean beam
                                    while(!(macroTime-prevMacroTime > 8300000 && macroTime-prevMacroTime < 8360000) && !(macroTime-prevMacroTime > 16600000 && macroTime-prevMacroTime < 17200000) && !(macroTime-prevMacroTime > 24900000 && macroTime-prevMacroTime < 25080000))
                                    {
                                        // update the ch0Tree to point at the
                                        // next macropulse so we can check to
                                        // see if it's out-of-sync with the
                                        // previous macropulse 
                                        prevMacroTime = macroTime;
                                        ch0Tree->GetEntry(macroNo+1);

                                        cout << "macroTime is " << macroTime << " prevMacroTime is " << prevMacroTime << endl;
                                    };

                                    // OK - ch0Tree is now in a period of
                                    // good beam. We need to move the cleaned
                                    // tree index forward to point at this new
                                    // area of good beam

                                    while (completeTime < macroTime-TIME_OFFSET)
                                    {
                                        do
                                        {
                                            i++;
                                            cTree->GetEntry(i);
                                            completeTime = (double)extTime*pow(2,32)+timetag;
                                        }
                                        while (chNo != j);

                                        cout << "completeTime is " << completeTime << " macroTime is " << macroTime << " i is " << i << "\r";
                                        fflush(stdout);
                                    }

                                    cout << endl;

                                    // now the cleaned tree index is pointing
                                    // to the area of good beam, so continue
                                    // filling the channel-specific tree
                                }

                                // lastly, because we shifted to a new
                                // macropulse at the start of this conditional,
                                // we should reset the event counter for this
                                // channel
                                evtNo = 0;
                            }

                            // pointing at the correct macropulse - fill the
                            // channel-specific tree with data from this channel
                        
                            // calculate the time elapsed since the start of the
                            // current macropulse and the time elapsed since the start
                            // of the micropulse
                            double timeDiff = completeTime-macroTime-TIME_OFFSET;

                            microTime = fmod(timeDiff,MICRO_PERIOD);
                            microNo = floor(timeDiff/MICRO_PERIOD);

                            for(int k = 0; k<dummyWaveform->size(); k++)
                            {
                                waveform.push_back(dummyWaveform->at(k));
                            }

                            if (chNo == 4)
                            {
                                fillTree(ch4Tree);
                            }

                            else
                            {
                                fillTree(ch6Tree);
                            }
                        }

                        else
                        {
                            // waveform mode data: fill ch2TreeW
                            
                            for(int k = 0; k<dummyWaveform->size(); k++)
                            {
                                waveform.push_back(dummyWaveform->at(k));
                            }

                            // populate waveform mode tree for monitor
                            fillTree(ch4TreeW);

                            if (prevEvtType==1)
                            {
                                ch0Tree->GetEntry(macroNo+1); // shift to the next
                                // DPP-mode macropulse in preparation for the next DPP mode
                                evtNo=0;
                            }
                        }

                        evtNo++; // increment event counter for channel 2
                        prevEvtType = evtType;

                        if (i%100==0)
                        {
                            cout << "Evt number " << i << ", completeTime " << completeTime << "\r";
                            fflush(stdout);
                        }
                    }
                    // prepare for next iteration of loop
                    prevCompleteTime = completeTime;
                }

                cout << endl << endl;
                break;
                }

            default:
                break;
        }
    }
}

void branch(TTree* tree)
{
    tree->Branch("macroNo",&ev.macroNo,"macroNo/i");
    tree->Branch("microNo",&ev.microNo,"microNo/i");
    tree->Branch("evtNo",&ev.evtNo,"evtNo/i");
    tree->Branch("macroTime",&ev.macroTime,"macroTime/d");
    tree->Branch("completeTime",&ev.completeTime,"completeTime/d");
    tree->Branch("microTime",&ev.microTime,"microTime/d");
    tree->Branch("targetPos",&ev.targetPos,"targetPos/i");
    tree->Branch("sgQ",&ev.sgQ,"sgQ/i");
    tree->Branch("lgQ",&ev.lgQ,"lgQ/i");
    tree->Branch("waveform",&ev.waveform);
}

void branchW(TTree* tree)
{
    tree->Branch("macroNo",&ev.macroNo,"macroNo/i");
    tree->Branch("evtNo",&ev.evtNo,"evtNo/i");
    tree->Branch("completeTime",&ev.completeTime,"completeTime/d");
    tree->Branch("waveform",&ev.waveform);
}

int main(int argc, char* argv[])
{

    // read in the raw tree name
    string runDir = argv[1];
    string runNo = argv[2];

    // Open the raw tree from the initial sort. If it doesn't exist, exit.
    TFile *file;
    TTree *tree;

    stringstream treeName;
    stringstream fileName;

    treeName << runDir << "-" << runNo; 
    fileName << analysispath <<"analysis/" << runDir << "/" << treeName.str() << ".root";

    file = new TFile(fileName.str().c_str(),"UPDATE");

    if(file->Get("tree"))
    {
        // Found the raw tree; start sorting the raw trees into new trees for
        // each channel.
        cout << "Found previous raw tree; creating cleaned subtrees for " << treeName << endl;

        tree = (TTree*)file->Get("tree");

        ch0Tree = new TTree("ch0Tree","");
        branch(ch0Tree);

        ch2Tree = new TTree("ch2Tree","");
        branch(ch2Tree);

        ch4Tree = new TTree("ch4Tree","");
        branch(ch4Tree);

        ch6Tree = new TTree("ch6Tree","");
        branch(ch6Tree);

        ch0TreeW = new TTree("ch0TreeW","");
        branchW(ch0TreeW);

        ch2TreeW = new TTree("ch2TreeW","");
        branchW(ch2TreeW);

        ch4TreeW = new TTree("ch4TreeW","");
        branchW(ch4TreeW);

        ch6TreeW = new TTree("ch6TreeW","");
        branchW(ch6TreeW);
    }

    else
    {
        cout << "Failed to find raw tree - exiting " << fileName << endl;
        exit(1);
    }

    cleanTree(tree);

    populateTrees();

    file->Write();

    file->Close();
}
