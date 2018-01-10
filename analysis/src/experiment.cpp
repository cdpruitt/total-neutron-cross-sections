#include "../include/experiment.h"
#include "../include/config.h"

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <utility> // for pair

using namespace std;

vector<pair<string,string>> getRelativePlotNames(string expName, string fileName)
{
    string relativePlotNameLocation = "../" + expName + "/" + fileName;
    ifstream dataFile(relativePlotNameLocation.c_str());
    if(!dataFile.is_open())
    {
        std::cout << "Failed to find relative plots configuration in " << relativePlotNameLocation << std::endl;
        exit(1);
    }

    vector<pair<string,string>> relativePlotNames;
    string str;
    string largerTargetName;
    string smallerTargetName;

    while(getline(dataFile,str))
    {
        // parse the two target names from the line in the file
        istringstream(str) >> largerTargetName >> smallerTargetName;
        relativePlotNames.push_back(pair<string,string> (largerTargetName,smallerTargetName));
    }

    return relativePlotNames;
}

// read the channel mapping of the current run
vector<pair<unsigned int, string>> getChannelMap(string expName, unsigned int runNumber)
{
    string channelMapLocation = "../" + expName + "/channelMap.txt";
    ifstream dataFile(channelMapLocation.c_str());
    if(!dataFile.is_open())
    {
        std::cout << "Failed to find channel mapping in " << channelMapLocation << std::endl;
        exit(1);
    }

    string str;
    unsigned int run;
    vector<pair<unsigned int, string>> channelMap;

    while(getline(dataFile,str))
    {
        // ignore comments in data file
        string delimiter = "-";
        stringstream tokenStream(str.substr(0,str.find(delimiter)));

        if(!(tokenStream>>run))
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
                channelMap.push_back(make_pair(i-1,tokens[i]));
            }

            break;
        }
    }

    if(channelMap.size()==0)
    {
        cerr << "Error: failed to recover a channel mapping from " << channelMapLocation << "." << endl;
        exit(1);
    }

    return channelMap;
}

// Determine the target order for a given run
vector<string> getTargetOrder(string expName, int runNumber)
{
    string targetOrderLocation = "../" + expName + "/TargetOrder.txt";
    ifstream dataFile(targetOrderLocation.c_str());
    if(!dataFile.is_open())
    {
        std::cout << "Failed to find target order data in " << targetOrderLocation << std::endl;
        exit(1);
    }

    string str;
    unsigned int run;
    vector<string> targetOrder;

    while(getline(dataFile,str))
    {
        // ignore comments in data file
        string delimiter = "-";
        stringstream tokenStream(str.substr(0,str.find(delimiter)));

        if(!(tokenStream>>run))
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
                targetOrder.push_back(tokens[i]);
            }
            break;
        }
    }

    return targetOrder;
}

// extract configuration data from an experiment directory
vector<string> readExperimentConfig(string expName, string fileName)
{
    string filePath = "../" + expName + "/" + fileName + ".txt";
    ifstream dataFile(filePath.c_str());
    if(!dataFile.is_open())
    {
        std::cout << "Failed to find experimental config data in " << filePath << std::endl;
        exit(1);
    }

    string str;
    vector<string> configContents;

    while(getline(dataFile,str))
    {
        configContents.push_back(str);
    }

    return configContents;
}

// Read analysis parameters
AnalysisConfig readAnalysisConfig(string expName)
{
    string analysisConfigLocation = "../" + expName + "/AnalysisConfig.txt";
    ifstream dataFile(analysisConfigLocation.c_str());
    if(!dataFile.is_open())
    {
        std::cout << "Failed to find analysis configuration in " << analysisConfigLocation << std::endl;
        exit(1);
    }

    string str;
    unsigned int run;

    AnalysisConfig analysisConfig;

    while(getline(dataFile,str))
    {
        // parse into tokens
        vector<string> tokens;
        istringstream iss(str);
        copy(istream_iterator<string>(iss),
                istream_iterator<string>(),
                back_inserter(tokens));

        if(!tokens.size())
        {
            continue;
        }

        if(tokens[0]=="Raw")
        {
            analysisConfig.RAW_TREE_FILE_NAME = tokens.back();
        }

        else if(tokens[0]=="Macropulse-assigned")
        {
            analysisConfig.MACROPULSE_ASSIGNED_FILE_NAME = tokens.back();
        }

        else if(tokens[0]=="Survived-veto")
        {
            analysisConfig.PASSED_VETO_FILE_NAME = tokens.back();
        }
        
        else if(tokens[0]=="Histogram")
        {
            analysisConfig.HISTOGRAM_FILE_NAME = tokens.back();
        }
        
        else if(tokens[0]=="Energy")
        {
            analysisConfig.ENERGY_PLOTS_FILE_NAME = tokens.back();
        }
        
        else if(tokens[0]=="DPP")
        {
            analysisConfig.DPP_TREE_NAME = tokens.back();
        }
        
        else if(tokens[0]=="Waveform")
        {
            analysisConfig.WAVEFORM_TREE_NAME = tokens.back();
        }
        
        else if(tokens[0]=="Macropulse")
        {
            analysisConfig.MACROPULSE_TREE_NAME = tokens.back();
        }

        else if(tokens[0]=="Monitor")
        {
            analysisConfig.MONITOR_TREE_NAME = tokens.back();
        }
        
        else if(tokens[0]=="Gamma")
        {
            analysisConfig.GAMMA_CORRECTION_TREE_NAME = tokens.back();
        }

        else if(tokens[0]=="Low")
        {
            analysisConfig.Q_RATIO_LOW_THRESHOLD = stod(tokens.back());
        }

        else if(tokens[0]=="High")
        {
            analysisConfig.Q_RATIO_HIGH_THRESHOLD = stod(tokens.back());
        }

        else if(tokens[0]=="lgQ_low")
        {
            analysisConfig.CHARGE_GATE_LOW_THRESHOLD = stod(tokens.back());
        }

        else if(tokens[0]=="lgQ_high")
        {
            analysisConfig.CHARGE_GATE_HIGH_THRESHOLD = stod(tokens.back());
        }
    }

    return analysisConfig;
}

// Read deadtime correction parameters
DeadtimeConfig readDeadtimeConfig(string expName)
{
    string deadtimeConfigLocation = "../" + expName + "/DeadtimeConfig.txt";
    ifstream dataFile(deadtimeConfigLocation.c_str());
    if(!dataFile.is_open())
    {
        std::cout << "Failed to find deadtime configuration in " << deadtimeConfigLocation << std::endl;
        exit(1);
    }

    string str;
    unsigned int run;

    DeadtimeConfig deadtimeConfig;

    while(getline(dataFile,str))
    {
        // parse into tokens
        vector<string> tokens;
        istringstream iss(str);
        copy(istream_iterator<string>(iss),
                istream_iterator<string>(),
                back_inserter(tokens));

        if(!tokens.size())
        {
            continue;
        }

        if(tokens[0]=="Logistic_k")
        {
            deadtimeConfig.LOGISTIC_K = stod(tokens.back());
        }

        else if(tokens[0]=="Logistic_mu")
        {
            deadtimeConfig.LOGISTIC_MU = stod(tokens.back());
        }
    }

    return deadtimeConfig;
}

// Read facility  parameters
FacilityConfig readFacilityConfig(string expName, int runNumber)
{
    string facilityConfigLocation = "../" + expName + "/FacilityConfig.txt";
    ifstream dataFile(facilityConfigLocation.c_str());
    if(!dataFile.is_open())
    {
        std::cout << "Failed to find facility configuration in " << facilityConfigLocation << std::endl;
        exit(1);
    }

    string str;
    unsigned int run;
    vector<string> facilityConfig;

    while(getline(dataFile,str))
    {
        // ignore comments in data file
        string delimiter = "-";
        stringstream tokenStream(str.substr(0,str.find(delimiter)));

        if(!(tokenStream>>run))
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
                facilityConfig.push_back(tokens[i]);
            }

            break;
        }
    }

    return FacilityConfig(facilityConfig);
}

// Read software CFD parameters
SoftwareCFDConfig readSoftwareCFDConfig(string expName, int runNumber)
{
    string softwareCFDConfigLocation = "../" + expName + "/softwareCFDConfig.txt";
    ifstream dataFile(softwareCFDConfigLocation.c_str());
    if(!dataFile.is_open())
    {
        std::cout << "Failed to find software CFD configuration in " << softwareCFDConfigLocation << std::endl;
        exit(1);
    }

    string str;
    unsigned int run;
    vector<string> softwareCFDConfig;

    while(getline(dataFile,str))
    {
        // ignore comments in data file
        string delimiter = "-";
        stringstream tokenStream(str.substr(0,str.find(delimiter)));

        if(!(tokenStream>>run))
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
                softwareCFDConfig.push_back(tokens[i]);
            }

            break;
        }
    }

    return SoftwareCFDConfig(softwareCFDConfig);
}


// Read time offset parameters
TimeConfig readTimeConfig(string expName, int runNumber)
{
    string timeConfigLocation = "../" + expName + "/TimeConfig.txt";
    ifstream dataFile(timeConfigLocation.c_str());
    if(!dataFile.is_open())
    {
        std::cout << "Failed to find time configuration in " << timeConfigLocation << std::endl;
        exit(1);
    }

    string str;
    unsigned int run;
    vector<string> timeConfig;

    while(getline(dataFile,str))
    {
        // ignore comments in data file
        string delimiter = "-";
        stringstream tokenStream(str.substr(0,str.find(delimiter)));

        if(!(tokenStream>>run))
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
                timeConfig.push_back(tokens[i]);
            }

            break;
        }
    }

    return TimeConfig(timeConfig);
}

// Read time offset parameters
DigitizerConfig readDigitizerConfig(string expName, int runNumber)
{
    string digitizerConfigLocation = "../" + expName + "/DigitizerConfig.txt";
    ifstream dataFile(digitizerConfigLocation.c_str());
    if(!dataFile.is_open())
    {
        std::cout << "Failed to find digitizer configuration in " << digitizerConfigLocation << std::endl;
        exit(1);
    }

    string str;
    unsigned int run;
    vector<string> digitizerConfig;

    while(getline(dataFile,str))
    {
        // ignore comments in data file
        string delimiter = "-";
        stringstream tokenStream(str.substr(0,str.find(delimiter)));

        if(!(tokenStream>>run))
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
                digitizerConfig.push_back(tokens[i]);
            }

            break;
        }
    }

    vector<pair<unsigned int, string>> channelMap = getChannelMap(expName, runNumber);

    return DigitizerConfig(channelMap, digitizerConfig);
}


// Read CS parameters 
CSConfig readCSConfig(string expName, int runNumber)
{
    string csConfigLocation = "../" + expName + "/CSConfig.txt";
    ifstream dataFile(csConfigLocation.c_str());
    if(!dataFile.is_open())
    {
        std::cout << "Failed to find CS configuration in " << csConfigLocation << std::endl;
        exit(1);
    }

    string str;
    unsigned int run;
    vector<string> csConfig;

    while(getline(dataFile,str))
    {
        // ignore comments in data file
        string delimiter = "-";
        stringstream tokenStream(str.substr(0,str.find(delimiter)));

        if(!(tokenStream>>run))
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
                csConfig.push_back(tokens[i]);
            }

            break;
        }
    }

    return CSConfig(csConfig);
}

// Read plot parameters 
PlotConfig readPlotConfig(string expName, int runNumber)
{
    string plotConfigLocation = "../" + expName + "/PlotConfig.txt";

    ifstream dataFile(plotConfigLocation.c_str());
    if(!dataFile.is_open())
    {
        std::cerr << "Failed to find plot configuration in " << plotConfigLocation << std::endl;
        exit(1);
    }

    string str;
    unsigned int run;
    vector<string> plotConfig;

    while(getline(dataFile,str))
    {
        // ignore comments in data file
        string delimiter = "-";
        stringstream tokenStream(str.substr(0,str.find(delimiter)));

        if(!(tokenStream>>run))
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
                plotConfig.push_back(tokens[i]);
            }

            break;
        }
    }

    return PlotConfig(plotConfig);
}

// Read target configuration parameters 
TargetConfig readTargetConfig(string expName, int runNumber)
{
    string targetConfigLocation = "../" + expName + "/TargetConfig.txt";
    ifstream dataFile(targetConfigLocation.c_str());
    if(!dataFile.is_open())
    {
        std::cout << "Failed to find target configuration in " << targetConfigLocation << std::endl;
        exit(1);
    }

    string str;
    unsigned int run;
    vector<pair<int,int>> targetConfig;

    while(getline(dataFile,str))
    {
        // ignore comments in data file
        string delimiter = "-";
        stringstream tokenStream(str.substr(0,str.find(delimiter)));

        if(!(tokenStream>>run))
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
                delimiter = "-";
                string lowGate = tokens[i].substr(0,tokens[i].find(delimiter));
                tokens[i] = tokens[i].erase(0,tokens[i].find(delimiter) + delimiter.length());

                delimiter = "\n";
                string highGate = tokens[i].substr(0,tokens[i].find(delimiter));

                targetConfig.push_back(pair<int,int>(atoi(lowGate.c_str()),atoi(highGate.c_str())));
            }

            break;
        }
    }

    vector<string> targetPositions = getTargetOrder(expName, runNumber);

    return TargetConfig(targetPositions, targetConfig);
}
