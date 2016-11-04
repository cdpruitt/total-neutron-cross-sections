#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

#include "../include/raw.h"
#include "../include/physicalConstants.h"
#include "../include/dataStructures.h"

extern RawEvent rawEvent;

// keep track of event statistics
long numberOfEvents = 0;

using namespace std;

// Organize formatted text output by channel number:
struct TextOutput {
    ofstream ch0;
    ofstream ch2;
    ofstream ch4;

    // print left and right detector channels independently (for diagnostics
    // only)
    ofstream channel6;
    ofstream channel7;

    // conglomerate of all channels
    ofstream total;

} textOutput;

const string ch0Name = "textOutput/ch0.txt";
const string ch2Name = "textOutput/ch2.txt";
const string ch4Name = "textOutput/ch4.txt";
const string channel6Name = "textOutput/channel6.txt";
const string channel7Name = "textOutput/channel7.txt";

// Pretty-print event data into a text file
void printEvent(RawEvent& rawEvent, TextOutput& text)
{
    ofstream* out;

    // segregate text output by channel
    switch(rawEvent.chNo)
    {
        case 0:
            out = &text.ch0;
            break;
        case 2:
            out = &text.ch2;
            break;
        case 4:
            out = &text.ch4;
            break;
        case 6:
            out = &text.channel6;
            break;
        case 7:
            out = &text.channel7;
            break;
    }

    // Print event's header data (i.e., data that exist in both DPP and waveform mode)
    *out << setfill('*') << setw(63) << "*" << endl;
    *out << "| EVENT " << left << setfill(' ') << setw(54) << numberOfEvents << "|" << endl;
    *out << "|" << right << setfill('-') << setw(62) << "|" << endl;
    //*out << "| run " << runInfo.runNumber << "-" << runInfo.subrunNumber << ", macro " << left << setfill(' ') << setw(44) << "|" << endl;
    *out << "| channel " << rawEvent.chNo;

    if(rawEvent.evtType==1)
    {
        *out << left << setfill(' ') << setw(51) << ", DPP mode" << "|" << endl;
    }

    else if (rawEvent.evtType==2)
    {
        *out << left << setfill(' ') << setw(51) << ", waveform mode " << "|" << endl;
    }

    else
    {
        *out << left << setfill(' ') << setw(50) << ", ERROR DETERMINING MODE" << " |" << endl;
        cerr << "Error: event type value out-of-range (DPP=1, waveform=2)" << endl;
        //cout << "Event number = " << rawEvent.evtNo[chNo] << endl;
    }

    stringstream temp; // used for formatting strings with units
    temp << rawEvent.size << " bytes";
    *out << "| size = " << left << setfill(' ') << setw(53) << temp.str() << "|" << endl;

    temp.str("");
    temp << rawEvent.timetag << " ns";
    *out << "| timetag = " << left << setw(50) << temp.str() << "|" << endl;

    *out << "|" << right << setfill('-') << setw(62) << "|" << endl;

    // Print event's body data (DPP-mode/waveform-mode dependent)
    if(rawEvent.evtType==1)
    {
        // DPP mode
        // Determine meaning of 'extras' word based on extraSelect:
        switch(rawEvent.extraSelect)
        {
            /*case 0: 
                // Extended timestamp + baseline
                temp.str("");
                temp << rawEvent.extTime << " *8.59 s";
                *out << "| extended time stamp = " << left << setfill(' ') << setw(38) << temp.str() << "|" << endl;

                // readout gives baseline*4, and we want just baseline
                *out << "| baseline = " << left << setfill(' ') << setw(49) << rawEvent.baseline << "|" << endl;
                break;

            case 1:
                // Extended timestamp + configuration file flags
                temp << rawEvent.extras2 << " *8.59 s";
                *out << "| extended time stamp = " << left << setfill(' ') << setw(34) << temp.str() << "|" << endl;
                // extract flags from bits 10:15 (0xfc00)
                *out << "| flags = " << left << setfill(' ') << setw(49) << (rawEvent.extras1 & 0xfc00) << "|" << endl;
                break; 
                */

            case 2:
                // Extended timestamp + configuration file flags + fine timestamp
                temp.str("");
                temp << rawEvent.extTime << " *8.59 s";
                *out << "| extended time stamp = " << left << setfill(' ') << setw(38) << temp.str() << "|" << endl;
                // extract flags from bits 10:15 (0xfc00)
                *out << "| flags = " << left << setfill(' ') << setw(52) << (rawEvent.flags & 0xfc00) << "|" << endl;
                *out << "| fine time stamp = " << left << setfill(' ') << setw(42) << rawEvent.fineTime << "|" << endl;
                break;

            case 3:
                // Maximum amplitude of pulse
                *out << "| maximum amplitude of pulse = " << left << setfill(' ') << setw(49) << rawEvent.extTime << "|" << endl;
                break;

            case 5:
                // Positive zero-crossing and negative zero-crossing (for manual
                // constant fraction discimination calculation)
                *out << "| PZC = " << left << setfill(' ') << setw(54) << rawEvent.PZC << "|" << endl;
                *out << "| NZC = " << left << setfill(' ') << setw(54) << rawEvent.NZC << "|" << endl;
                break;
            /*case 7:
                // fixed value of 0x12345678 (for diagnostics)
                const unsigned int extras = (extras2 << 16) | extras1;
                *out << "| Fixed value of 305419896 (0x12345678) outputted: " << left << setfill(' ') << setw(17) << rawEvent.extras << "|" << endl;
                break;
             */
            default:
                cerr << "Error: found unsupported value of EXTRAS. Exiting..." << endl;
                exit(1);
        }

        *out << "| short gate charge = " << left << setw(40) << rawEvent.sgQ << "|" << endl;
        *out << "| long gate charge = " << left << setw(41) << rawEvent.lgQ << "|" << endl;

        // Pile-up detection not operational in current washudaq - ignore it.
        //out << "| pile-up detected = " << left << setw(40) << puRej << " |" << endl;

        temp.str("");
        temp << rawEvent.nSamp << " samples";
        *out << "| waveform length = " << left << setw(41) << temp.str() << " |" << endl;
        *out << "|" << right << setfill('-') << setw(62) << "|" << endl;

        if(rawEvent.nSamp > 0)
            // print wavelet data
        {
            *out << left << setfill(' ') << setw(62) << "| Waveform samples" << "|" << endl;
            *out << "|" << right << setfill(' ') << setw(62) << "|" << endl;

            for(std::vector<int>::size_type i = 0; i != rawEvent.waveform->size(); i++)
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

                temp << right << setfill(' ') << setw(6) << rawEvent.waveform->at(i);

                if(i%10==9 || i==rawEvent.nSamp-1)
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

    else if(rawEvent.evtType==2)
    {
        stringstream temp;
        temp << rawEvent.nSamp << " samples";
        *out << "| waveform length = " << left << setfill(' ') << setw(41) << temp.str() << " |" << endl;
        *out << "|" << right << setfill('-') << setw(62) << "|" << endl;

        for(int i=0;(size_t)i<rawEvent.nSamp;i++)
        {

            if(rawEvent.nSamp>1000 && (i%1000 == 0))
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

            *out << right << setfill(' ') << setw(6) << rawEvent.waveform->at(i);

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

int main(int, char* argv[])
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

    textOutput.ch0.open(ch0Name);
    textOutput.ch2.open(ch2Name);
    textOutput.ch4.open(ch4Name);
    textOutput.channel6.open(channel6Name);
    textOutput.channel7.open(channel7Name);

    rawEvent.waveform = new vector<int>;

    // we're now pointing at the first 16-bit word in the data stream
    // start looping through the evtfile to extract events
    while(!inFile.eof() /* use to truncate sort && numberOfEvents<1000000*/)
    {
        if(readEvent(inFile))
        {
            printEvent(rawEvent, textOutput);
        }

        numberOfEvents++;

        if (numberOfEvents%10000 == 0)
        {
            cout << "Processed " << numberOfEvents << " events\r";
            fflush(stdout);
        }
    }

    // reached end of input file
    cout << "Finished processing event file" << endl;
    cout << "Total events: " << numberOfEvents << endl;

    inFile.close();
    textOutput.ch0.close();
    textOutput.ch2.close();
    textOutput.ch4.close();
    textOutput.channel6.close();
    textOutput.channel7.close();
}
