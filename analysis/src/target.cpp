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

        if((tokens[0]=="Diameter") && (tokens[1] == "uncertainty:"))
        {
            diameterUncertainty = atof(tokens[2].c_str());
        }

        else if((tokens[0]=="Mass") && (tokens[1] == "uncertainty:"))
        {
            massUncertainty = atof(tokens[2].c_str());
        }

        else if((tokens[0]=="Molar") && (tokens[2] == "uncertainty:"))
        {
            molMassUncertainty = atof(tokens[3].c_str());
        }

        else if(tokens[0]=="Name:")
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

string Target::getName() const
{
    return name;
}

double Target::getLength() const
{
    return length;
}

double Target::getDiameter() const
{
    return diameter;
}

double Target::getMass() const
{
    return mass;
}

double Target::getMolarMass() const
{
    return molMass;
}

double Target::getDiameterUncertainty() const
{
    return diameterUncertainty;
}

double Target::getMassUncertainty() const
{
    return massUncertainty;
}

double Target::getMolarMassUncertainty() const
{
    return molMassUncertainty;
}

void Target::setName(string n)
{
    name = n;
}

void Target::setLength(double l)
{
    length = l;
}

void Target::setDiameter(double d)
{
    diameter = d;
}

void Target::setMass(double m)
{
    mass = m;
}

void Target::setMolarMass(double mm)
{
    molMass = mm;
}

void Target::setDiameterUncertainty(double d)
{
    diameterUncertainty = d;
}

void Target::setMassUncertainty(double m)
{
    massUncertainty = m;
}

void Target::setMolarMassUncertainty(double mm)
{
    molMassUncertainty = mm;
}
