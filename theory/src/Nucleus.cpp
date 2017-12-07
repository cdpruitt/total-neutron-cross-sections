#include "../include/Nucleus.h"

#include <vector>
#include <string>

#include <fstream>
#include <iostream>
#include <sstream>
#include <iterator>

using namespace std;

Nucleus::Nucleus() {}

Nucleus::Nucleus(std::string configFileName)
{
    ifstream configFile(configFileName);

    if(!configFile.is_open())
    {
        cerr << "Error: failed to open nucleus config file "
            << configFileName << endl;
        
        exit(1);
    }

    string buffer;

    while(getline(configFile,buffer))
    {
        vector<string> tokens;
        istringstream iss(buffer);
        copy(istream_iterator<string>(iss),
                istream_iterator<string>(),
                back_inserter(tokens));

        if(tokens[0]=="Name")
        {
            Name = tokens[2];
            continue;
        }

        if(tokens[0]=="A")
        {
            A = stod(tokens[2]);
            continue;
        }

        if(tokens[0]=="Z")
        {
            Z = stod(tokens[2]);
            continue;
        }

    }
}
