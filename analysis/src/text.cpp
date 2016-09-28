/******************************************************************************
  text.cpp
 ******************************************************************************/

// This is called by ./sort.sh and is the first step in analyzing the total
// cross section data.
// It reads raw hexidecimal data from an .evt file and creates a ROOT file,
// runX-Y_raw.root (where X is the run number, and Y is the sub-run number),
// with one ROOT tree containing all events from the input .evt file.
// 
// This file expects .evt file data with the following structure:
//
//      ---------------------------------------------------------------------------
//      EVENT HEADER:
//
//      Size            |   uint32; number of bytes in the event (self-inclusive)
//      Event Type      |   uint32; possible values are 1 (DPP) or 2 (waveform)
//      Channel         |   uint32; channel number of event (range: 0-7)
//      Time Tag        |   uint32; coarse time (2 ns units) of event trigger
//
//      (in every event, there is an EVENT HEADER)
//      ---------------------------------------------------------------------------
//
//      ---------------------------------------------------------------------------
//      DPP EVENT BODY:
//
//      Extra select    |   uint16; describes the contents of Extras
//      Extras          |   uint32; additional PSD data (see DPP Event body section)
//      Short Gate Q    |   uint16; charge integrated in the short gate
//      Long Gate Q     |   uint16; charge integrated in the long gate
//      Baseline        |   uint16; baseline level of ADC
//      Pile up rej.    |   uint16; pile-up rejection flag (not yet implemented)
//      Probe info      |   uint16; turns on an analog probe for examining waveform
//      Num. of samples |   uint32; number of waveform samples collected in Samples
//      Samples         |   series of uint16 giving raw waveform samples
//
//      (only events of type 1 have DPP EVENT BODY)
//      ---------------------------------------------------------------------------
//
//      ---------------------------------------------------------------------------
//
//      WAVEFORM EVENT BODY
//
//      Num. of samples |   uint32; number of waveform samples to follow (mixed mode)
//      Samples         |   series of uint16 giving raw waveform samples
//
//      (only events of type 2 have WAVEFORM EVENT BODY)
//      ---------------------------------------------------------------------------

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

#include "include/physicalConstants.h"

using namespace std;

// To populate events from the raw file into a ROOT tree, we provide
// an event structure listing all the variables that comprise an event.
// When a new event is added to the tree, a subset of these variables is
// included in the new event being saved.
struct Event
{
    // Event header variables:
    unsigned int size;    // size of event, in bytes
    unsigned int evtType; // "event type", either 1 (DPP) or 2 (waveform)
    unsigned int chNo;    // "channel number" indicates event origin (i.e.,
    // detector, monitor, target changer)
    double timetag;       // "coarse timestamp of event" includes 32 bits of time
    // information. Units are the same as the sample period of
    // the digitizer (e.g. 2 ns or 5 ns)

    // Event body variables:
    unsigned int extraSelect; // indicates the meaning of the "extras" words
    unsigned int sgQ;     // "short gate integrated charge" provides the charge
    // integral over an adjustable range of the event's peak 
    unsigned int lgQ;     // "long gate integrated charge" provides the charge
    // integral over an adjustable range of the event's peak 
    unsigned int nSamp; // number of samples in the event's waveform
    vector<int> waveform; // "digital waveform of event" is a series of waveform
    // samples for each event

    // Variables extracted from "extras"

    unsigned int baseline;
    unsigned int flags;
    unsigned int extTime; // "extended timestamp of event" extends timetag with
    // 16 additional bits for times greater than 2^32
    // sample periods.
    unsigned int fineTime;// "fine timestamp of event" sub-divides timetag with
    // 10 additional bits of time granularity, with units of
    // (sample period)/2^10 units.
    unsigned int PZC; // positive zero-crossing (for manual CFD calculation)
    unsigned int NZC; // negative zero-crossing (for manual CFD calculation)
} ev;

// If enabled, organize formatted text output by channel number:
struct TextOutput {
    ofstream targetChanger;
    ofstream monitor;
    ofstream detectorT;

    // print left and right detector channels independently (for diagnostics
    // only)
    ofstream detectorL;
    ofstream detectorR;

    // conglomerate of all channels
    ofstream total;

} textOutput;

const string targetChangerName = "textOutput/targetChanger.txt";
const string monitorName = "textOutput/monitor.txt";
const string detectorTName = "textOutput/detectorT.txt";
const string detectorLName = "textOutput/detectorL.txt";
const string detectorRName = "textOutput/detectorR.txt";

// Raw data is stored as hexadecimal words (16 bits long) in the .evt files
int const BufferWords = 1; // number of chars per buffer word
int const BufferBytes = BufferWords*2; // number of bytes per buffer word

// track number of processed events of each type
long numberOfEvents = 0;
long numberOfDPPs = 0;
long numberOfWaveforms = 0;
long numberOfCh0Waveforms = 0;
long numberOfCh2Waveforms = 0;
long numberOfCh4Waveforms = 0;

// extract a single event from the raw event file and fill its data into the
// tree
void readEvent(ifstream& evtfile)
{
    // clear all DPP-specific event variables to prepare for reading event
    ev.extTime = 0;
    ev.fineTime = 0;
    ev.sgQ = 0;
    ev.lgQ = 0;

    unsigned short buffer[BufferWords];

    // start reading event header (common to all events)

    evtfile.read((char*)buffer,BufferBytes);
    unsigned short size1 = buffer[0];
    evtfile.read((char*)buffer,BufferBytes);
    unsigned short size2 = buffer[0];
    ev.size = (size2 << 16) | size1;

    evtfile.read((char*)buffer,BufferBytes);
    unsigned short evtType1 = buffer[0];
    evtfile.read((char*)buffer,BufferBytes);
    unsigned short evtType2 = buffer[0];
    ev.evtType = (evtType2 << 16) | evtType1;

    evtfile.read((char*)buffer,BufferBytes);
    unsigned short chNo1 = buffer[0];
    evtfile.read((char*)buffer,BufferBytes);
    unsigned short chNo2 = buffer[0];
    ev.chNo = (chNo2 << 16) | chNo1;

    evtfile.read((char*)buffer,BufferBytes);
    unsigned short timetag1 = buffer[0];
    evtfile.read((char*)buffer,BufferBytes);
    unsigned short timetag2 = buffer[0];
    ev.timetag = (timetag2<< 16) | timetag1;
    ev.timetag *= SAMPLE_PERIOD; // timetag converted from samples to ns
    // finished reading event header

    if(ev.evtType==1)
    {
        // start reading event data specific to DPP events

        // EXTRAS is a 32-bit word whose content varies with the value of EXTRA_SELECT.
        // [DEFAULT]    0: extended timestamp (bits 16-31) and baseline*4 (bits 0-15)
        //              1: extended timestamp (bits 16-31) and flags (bits 0-15?)
        //              2: extended timestamp (bits 16-31), flags (bits 10-15), and fine timestamp
        //                  (bits 0-9)
        //              3: pulse peak value (bits 0-15) [UNTESTED FEATURE]
        //              5: CFD positive ZC (bits 16-31) and negative ZC (bits 0-15)
        //              7: fixed value of 0x12345678

        evtfile.read((char*)buffer,BufferBytes);
        ev.extraSelect = buffer[0];

        evtfile.read((char*)buffer,BufferBytes);
        unsigned int extras1 = buffer[0];

        evtfile.read((char*)buffer,BufferBytes);
        unsigned int extras2 = buffer[0];

        switch(ev.extraSelect)
        {
            case 2:
                // retrieve extended time from bits 16-31 (extras2)
                ev.extTime = extras2;
                // retrieve configuration file flags from bits 10:15 (0xfc00)
                ev.flags = (extras1 & 0xfc00);
                // retrieve fine time from bits 0:9 (0x03ff)
                ev.fineTime = (extras1 & 0x03ff);
                break;
            case 5:
                // retrieve positive and negative zero-crossings
                ev.PZC = extras2;
                ev.NZC = extras1;

                // other cases not currently implemented
            default:
                cout << "Error: encountered unimplemented value of extraSelect" << endl;
                exit(1);
                break;
        }

        evtfile.read((char*)buffer,BufferBytes);
        ev.sgQ = buffer[0];

        evtfile.read((char*)buffer,BufferBytes);
        ev.lgQ = buffer[0];

        // "pile-up rejection" flag (not implemented in current acquisition software,
        // so discard this data)
        //puRej = buffer[0];
        evtfile.read((char*)buffer,BufferBytes);

        // "probe" indicates whether an additional diagnostic waveform will be captured along with the
        // normal waveform. The top bit turns on this probe waveform; the bottom two bits describe the
        // probe type.
        // Note that this is a diagnostic tool for acquisition turned off during production runs
        //probe = buffer[0];
        evtfile.read((char*)buffer,BufferBytes);

        numberOfDPPs++;
    }

    else if (ev.evtType==2)
    {
        // waveform mode event
        numberOfWaveforms++;
    }

    else
    {
        cout << "Error: evtType must be either 1 (DPP mode) or 2 (Waveform mode)." << endl;
        exit(1);
    }

    // the DPP-mode only data has been read out, so now let's read out the data
    // that's always present (the number of waveform samples and the waveforms)

    evtfile.read((char*)buffer,BufferBytes);
    unsigned short nSamp1 = buffer[0];
    evtfile.read((char*)buffer,BufferBytes);
    unsigned short nSamp2 = buffer[0];
    ev.nSamp = (nSamp2 << 16) | nSamp1;

    ev.waveform.clear();
    for(int i=0;(size_t)i<ev.nSamp;i++)
    {
        evtfile.read((char*)buffer,BufferBytes);
        ev.waveform.push_back(buffer[0]);
    }

    // Update statistics on events
    numberOfEvents++;
    if(ev.evtType==2)
    {
        if(ev.chNo==0)
        {
            numberOfCh0Waveforms++;
        }

        if(ev.chNo==2)
        {
            numberOfCh2Waveforms++;
        }

        if(ev.chNo==4)
        {
            numberOfCh4Waveforms++;
        }
    }
}

// Pretty-print event data into a text file
void printEvent(Event& event, TextOutput& text)
{
    ofstream* out;

    // segregate text output by channel
    switch(event.chNo)
    {
        case 0:
            out = &text.targetChanger;
            break;
        case 2:
            out = &text.monitor;
            break;
        case 4:
            out = &text.detectorT;
            break;
        case 6:
            out = &text.detectorL;
            break;
        case 7:
            out = &text.detectorR;
            break;
    }

    // Print event's header data (i.e., data that exist in both DPP and waveform mode)
    *out << setfill('*') << setw(63) << "*" << endl;
    *out << "| EVENT " << left << setfill(' ') << setw(54) << numberOfEvents << "|" << endl;
    *out << "|" << right << setfill('-') << setw(62) << "|" << endl;
    //*out << "| run " << runInfo.runNumber << "-" << runInfo.subrunNumber << ", macro " << left << setfill(' ') << setw(44) << "|" << endl;
    *out << "| channel " << event.chNo;

    if(event.evtType==1)
    {
        *out << left << setfill(' ') << setw(51) << ", DPP mode" << "|" << endl;
    }

    else if (event.evtType==2)
    {
        *out << left << setfill(' ') << setw(51) << ", waveform mode " << "|" << endl;
    }

    else
    {
        *out << left << setfill(' ') << setw(50) << ", ERROR DETERMINING MODE" << " |" << endl;
        cout << "Error: event type value out-of-range (DPP=1, waveform=2)" << endl;
        //cout << "Event number = " << event.evtNo[chNo] << endl;
    }

    stringstream temp; // used for formatting strings with units
    temp << event.size << " bytes";
    *out << "| size = " << left << setfill(' ') << setw(53) << temp.str() << "|" << endl;

    temp.str("");
    temp << event.timetag << " ns";
    *out << "| timetag = " << left << setw(50) << temp.str() << "|" << endl;

    *out << "|" << right << setfill('-') << setw(62) << "|" << endl;

    // Print event's body data (DPP-mode/waveform-mode dependent)
    if(event.evtType==1)
    {
        // DPP mode
        // Determine meaning of 'extras' word based on extraSelect:
        switch(event.extraSelect)
        {
            /*case 0: 
                // Extended timestamp + baseline
                temp.str("");
                temp << event.extTime << " *8.59 s";
                *out << "| extended time stamp = " << left << setfill(' ') << setw(38) << temp.str() << "|" << endl;

                // readout gives baseline*4, and we want just baseline
                *out << "| baseline = " << left << setfill(' ') << setw(49) << event.baseline << "|" << endl;
                break;

            case 1:
                // Extended timestamp + configuration file flags
                temp << event.extras2 << " *8.59 s";
                *out << "| extended time stamp = " << left << setfill(' ') << setw(34) << temp.str() << "|" << endl;
                // extract flags from bits 10:15 (0xfc00)
                *out << "| flags = " << left << setfill(' ') << setw(49) << (event.extras1 & 0xfc00) << "|" << endl;
                break; 
                */

            case 2:
                // Extended timestamp + configuration file flags + fine timestamp
                temp.str("");
                temp << event.extTime << " *8.59 s";
                *out << "| extended time stamp = " << left << setfill(' ') << setw(38) << temp.str() << "|" << endl;
                // extract flags from bits 10:15 (0xfc00)
                *out << "| flags = " << left << setfill(' ') << setw(52) << (event.flags & 0xfc00) << "|" << endl;
                *out << "| fine time stamp = " << left << setfill(' ') << setw(42) << event.fineTime << "|" << endl;
                break;

            /*case 3:
                // Maximum amplitude of pulse
                *out << "| maximum amplitude of pulse = " << left << setfill(' ') << setw(49) << event.extras1 << "|" << endl;
                break;

            case 5:
                // Positive zero-crossing and negative zero-crossing (for manual
                // constant fraction discimination calculation)
                *out << "| PZC = " << left << setfill(' ') << setw(54) << event.extras2 << "|" << endl;
                *out << "| NZC = " << left << setfill(' ') << setw(54) << event.extras1 << "|" << endl;
                break;

            case 7:
                // fixed value of 0x12345678 (for diagnostics)
                const unsigned int extras = (extras2 << 16) | extras1;
                *out << "| Fixed value of 305419896 (0x12345678) outputted: " << left << setfill(' ') << setw(17) << event.extras << "|" << endl;
                break;
                */
        }

        *out << "| short gate charge = " << left << setw(40) << event.sgQ << "|" << endl;
        *out << "| long gate charge = " << left << setw(41) << event.lgQ << "|" << endl;

        // Pile-up detection not operational in current washudaq - ignore it.
        //out << "| pile-up detected = " << left << setw(40) << puRej << " |" << endl;

        temp.str("");
        temp << event.nSamp << " samples";
        *out << "| waveform length = " << left << setw(41) << temp.str() << " |" << endl;
        *out << "|" << right << setfill('-') << setw(62) << "|" << endl;

        if(event.nSamp > 0)
            // print wavelet data
        {
            *out << left << setfill(' ') << setw(62) << "| Waveform samples" << "|" << endl;
            *out << "|" << right << setfill(' ') << setw(62) << "|" << endl;

            for(std::vector<int>::size_type i = 0; i != event.waveform.size(); i++)
            {
                if(i%100 == 0 && i>0)
                {
                    *out << "|" << right << setfill(' ') << setw(62) << "|" << endl;
                }

                if(i%10 == 0)
                {
                    temp.str("");
                    temp << "|";
                }

                temp << right << setfill(' ') << setw(6) << event.waveform[i];

                if(i%10==9 || i==event.nSamp-1)
                {
                    *out << left << setw(62) << temp.str() << "|" << endl;
                } 
            }

            *out << "|" << right << setfill('-') << setw(62) << "|" << endl;
        }

        /*if((probe & 0x8000)==0x8000)
        {
            // print analog probe data
            out << "| Analog probe enabled" << right << setfill(' ') << setw(41) << "|" << endl;

            temp.str("");
            temp << nSamp << " samples";
            out << "| waveform length = " << left << setw(41) << temp.str() << " |" << endl;
            out << "|" << right << setfill('-') << setw(62) << "|" << endl;

            if((probe & 0x0003)==0x0002)
            {
                // analog probe is CFD
                if(anSamp > 0)
                {
                    out << left << setfill(' ') << setw(62) << "| CFD samples" << "|" << endl;
                    out << "|" << right << setfill(' ') << setw(62) << "|" << endl;

                    for(std::vector<int>::size_type i = 0; i != anProbe.size(); i++)
                    {
                        if(i%100 == 0 && i>0)
                        {
                            out << "|" << right << setfill(' ') << setw(62) << "|" << endl;
                        }

                        if(i%10 == 0)
                        {
                            temp.str("");
                            temp << "|";
                        }

                        temp << right << setfill(' ') << setw(6) << anProbe[i];

                        if(i%10==9 || i==anSamp-1)
                        {
                            out << left << setw(62) << temp.str() << "|" << endl;
                        } 

                    }

                    // done with this event; increment the wavelet counter to get ready for the next
                    // wavelet.
                    //nCFDs++;
                }
            }

            if((probe & 0x0003)==0x0001)
            {
                // analog probe is baseline
                if(anSamp > 0)
                {

                    out << left << setfill(' ') << setw(62) << "| Baseline samples" << "|" << endl;
                    out << "|" << right << setfill(' ') << setw(62) << "|" << endl;

                    for(int i=0;i<anSamp;i++)
                    {
                        //listBaselines[nBaselines]->SetBinContent(i,buffer[0]);

                        if(i%100 == 0 && i>0)
                        {
                            out << "|" << right << setfill(' ') << setw(62) << "|" << endl;
                        }

                        if(i%10 == 0)
                        {
                            temp.str("");
                            temp << "|";
                        }

                        temp << right << setfill(' ') << setw(6) << buffer[0];

                        if(i%10==9 || i==anSamp-1)
                        {
                            out << left << setw(62) << temp.str() << "|" << endl;
                        } 
                    }
                }
            }
        }

        else
        {
            out << "| Analog probe disabled" << right << setfill(' ') << setw(40) << "|" << endl;
            out << "|" << right << setfill(' ') << setw(62) << "|" << endl;
            out << setfill('*') << setw(63) << "*" << endl;
        }
        */
    }

    else if(event.evtType==2)
    {
        stringstream temp;
        temp << event.nSamp << " samples";
        *out << "| waveform length = " << left << setfill(' ') << setw(41) << temp.str() << " |" << endl;
        *out << "|" << right << setfill('-') << setw(62) << "|" << endl;

        for(int i=0;(size_t)i<event.nSamp;i++)
        {

            if(event.nSamp>1000 && (i%1000 == 0))
            {
                *out << "|" << right << setfill(' ') << setw(62) << "|" << endl;
                stringstream temp;
                temp << "Samples " << i << "-" << i+1000;
                *out << "| " << left << setw(60) << temp.str() << "|" << endl;
            }

            if(i%100 == 0)
            {
                *out << "|" << right << setfill(' ') << setw(62) << "|" << endl;
            }

            if(i%10 == 0)
            {
                *out << "|";
            } 

            *out << right << setfill(' ') << setw(6) << event.waveform[i];

            if(i%10 == 9)
            {
                *out << " |" << endl;
            } 
        }

        *out << "|" << right << setfill(' ') << setw(62) << "|" << endl;
    }

    else
    {
        cout << "ERROR: unknown value for event type" << endl;
        return;
    }
    *out << endl;
}

int main(int argc, char* argv[])
{
    cout << endl << "Entering ./text..." << endl;

    string inFileName = argv[1];

    // attempt to open input file
    ifstream inFile;
    inFile.open(inFileName,ios::binary);
    if (!inFile)
    {
        cout << "Failed to open " << inFileName << ". Please check that the file exists" << endl;
        exit(1);
    }

    cout << inFileName << " opened successfully. Start reading events..." << endl;

    textOutput.targetChanger.open(targetChangerName);
    textOutput.monitor.open(monitorName);
    textOutput.detectorT.open(detectorTName);
    textOutput.detectorL.open(detectorLName);
    textOutput.detectorR.open(detectorRName);

    // we're now pointing at the first 16-bit word in the data stream
    // start looping through the evtfile to extract events
    while(!inFile.eof() /* use to truncate sort && numberOfEvents<1000000*/)
    {
        readEvent(inFile);
        printEvent(ev, textOutput);

        if (numberOfEvents%10000 == 0)
        {
            cout << "Processed " << numberOfEvents << " events\r";
            fflush(stdout);
        }
    }

    // reached end of input file
    cout << "Finished processing event file" << endl;
    cout << "Total events: " << numberOfEvents << endl;
    cout << "Total number of DPP-mode events processed = " << numberOfDPPs << endl;
    cout << "Total number of waveform-mode events processed = " << numberOfWaveforms << endl;

    inFile.close();
    textOutput.targetChanger.close();
    textOutput.monitor.close();
    textOutput.detectorT.close();
    textOutput.detectorL.close();
    textOutput.detectorR.close();
}
