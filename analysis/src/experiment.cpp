#include "../include/experiment.h"

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

// extract configuration data from an experiment directory
vector<string> readExperimentConfig(string expName, string fileName)
{
    string filePath = "../" + expName + "/" + fileName + ".txt";
    ifstream dataFile(filePath.c_str());
    if(!dataFile.is_open())
    {
        std::cout << "Failed to find target names in " << filePath << std::endl;
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

// Determine the target order for a given run
vector<string> getTargetOrder(string expName, int runNumber)
{
    string targetOrderLocation = "../" + expName + "/targetOrder.txt";
    ifstream dataFile(targetOrderLocation.c_str());
    if(!dataFile.is_open())
    {
        std::cout << "Failed to find target order data in " << targetOrderLocation << std::endl;
        exit(1);
    }

    string str;
    vector<string> targetOrder;

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
                targetOrder.push_back(tokens[i]);
            }
            break;
        }
    }

    return targetOrder;
}
