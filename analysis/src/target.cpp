#include <string>
#include <iostream>
#include <iterator>
#include <fstream>
#include <sstream>
//#include <stdlib.h>

#include "../include/target.h"

using namespace std;

Target::Target()
{
}

Target::Target(string targetDataLocation)
{
    ifstream dataFile(targetDataLocation.c_str());
    if(!dataFile.is_open())
    {
        std::cout << "Attempted to create Target, but failed to find target data in " << targetDataLocation << std::endl;
        exit(1);
    }

    string str;

    while(getline(dataFile,str))
    {
        // parse into tokens

        vector<string> tokens;
        istringstream iss(str);
        copy(istream_iterator<string>(iss),
                istream_iterator<string>(),
                back_inserter(tokens));

        if(tokens[0]=="Name:")
        {
            name = tokens[1];
        }

        else if(tokens[0]=="Length:")
        {
            length = atof(tokens[1].c_str());
        }

        else if(tokens[0]=="Diameter:")
        {
            diameter = atof(tokens[1].c_str());
        }

        else if(tokens[0]=="Mass:")
        {
            mass = atof(tokens[1].c_str());
        }

        else if(tokens[0]=="Molar")
        {
            molMass = atof(tokens[2].c_str());
        }

        else
        {
            cerr << "Error - couldn't parse a line in a targetData text file" << endl;
            exit(1);
        }
    }
}

string Target::getName()
{
    return name;
}

double Target::getLength()
{
    return length;
}

double Target::getDiameter()
{
    return diameter;
}

double Target::getMass()
{
    return mass;
}

double Target::getMolMass()
{
    return molMass;
}
