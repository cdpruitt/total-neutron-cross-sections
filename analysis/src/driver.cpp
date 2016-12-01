// project-specific classes
#include "../include/analysisConstants.h"
#include "../include/plottingConstants.h"

#include "../include/raw.h"
#include "../include/separate.h"
#include "../include/histos.h"
#include "../include/plots.h"
#include "../include/waveform.h"
#include "../include/crossSection.h"
#include "../include/veto.h"

// ROOT library classes
#include "TFile.h"
#include "TTree.h"

// STL classes
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

using namespace std;

vector<string> getChannelMap(string expName, int runNumber)
{
    string channelMapLocation = "../" + expName + "/channelMap.txt";
    ifstream dataFile(channelMapLocation.c_str());
    if(!dataFile.is_open())
    {
        std::cout << "Failed to find channel mapping in " << channelMapLocation << std::endl;
        exit(1);
    }

    string str;
    vector<string> channelMap;

    while(getline(dataFile,str))
    {
        // ignore comments in data file
        string delimiter = "-";
        string token = str.substr(0,str.find(delimiter));
        if(!atoi(token.c_str()))
        {
            // This line starts with a non-integer and is thus a comment; ignore
            continue;
        }

        // parse data lines into space-delimited tokens
        vector<string> tokens;
        istringstream iss(str);
        copy(istream_iterator<string>(iss),
                istream_iterator<string>(),
                back_inserter(tokens));

        // extract run numbers from first token
        string lowRun = tokens[0].substr(0,tokens[0].find(delimiter));
        tokens[0] = tokens[0].erase(0,tokens[0].find(delimiter) + delimiter.length());

        delimiter = "\n";
        string highRun = tokens[0].substr(0,tokens[0].find(delimiter));
        
        if(atoi(lowRun.c_str()) <= runNumber && runNumber <= atoi(highRun.c_str()))
        {
            for(int i=1; (size_t)i<tokens.size(); i++)
            {
                channelMap.push_back(tokens[i]);
            }
            break;
        }
    }

    return channelMap;
}



int main(int, char* argv[])
{
    /*************************************************************************/
    /* Set up input/output filenames */
    /*************************************************************************/

    // location of digitizer-produced data file
    string rawDataFileName = argv[1];

    // location of directory where all analysis output will be stored
    string analysisDirectory = argv[2];
 
    // location of experimental configuration files
    string experimentName = argv[3];

    int runNumber = atoi(argv[4]);

    // create names of all output files
    string rawTreeFileName = analysisDirectory + "raw.root";
    string sortedFileName = analysisDirectory + "sorted.root";
    string vetoedFileName = analysisDirectory + "vetoed.root";
    string histoFileName = analysisDirectory + "histos.root";
    string waveformFileName = analysisDirectory + "waveform.root";
    string errorFileName = analysisDirectory + "error.txt";

    vector<string> channelMap = getChannelMap(experimentName, runNumber);

    /*************************************************************************/
    /* Start analysis:
     * Populate raw event data into a tree */
    /*************************************************************************/
    ifstream f(rawTreeFileName);
    if(!f.good())
    {
        // create a raw data tree for this subrun
        readRawData(rawDataFileName,rawTreeFileName);
    }

    else
    {
        cout << "Raw data tree already exists." << endl;
        f.close();
    }

    /*************************************************************************/
    /* "Process" raw events by assigning time and target data */
    /*************************************************************************/
    ifstream p(sortedFileName);
    if(!p.good())
    {
        // separate all data by channel and event type
        separateByChannel(rawTreeFileName, sortedFileName, channelMap);
    }

    else
    {
        cout << "Sorted data tree already exists." << endl;
        p.close();
    }

    /*************************************************************************/
    /* Eliminate detector events using the veto paddle */
    /*************************************************************************/
    ifstream v(vetoedFileName);
    if(!v.good())
    {
        vetoEvents(sortedFileName, vetoedFileName, detectorNames, "veto");
    }

    else
    {
        cout << "Vetoed data tree already exists." << endl;
        v.close();
    }

    /*************************************************************************/
    /* Process events into histograms in preparation for cross section
     * calculation */
    /*************************************************************************/
    ifstream h(histoFileName);
    if(!h.good())
    {
        histos(sortedFileName, vetoedFileName, histoFileName, channelMap);
    }

    else
    {
        cout << "Histos already exist." << endl;
        h.close();
    }

    /*************************************************************************/
    /* Process events into histograms in preparation for cross section
     * calculation */
    /*************************************************************************/
    correctForDeadtime(histoFileName, histoFileName, detectorNames);

    return 0;
}
