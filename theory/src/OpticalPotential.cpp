#include "../include/OpticalPotential.h"

#include <vector>
#include <string>

#include <fstream>
#include <iostream>
#include <sstream>
#include <iterator>

using namespace std;

OpticalPotential::OpticalPotential() {}

OpticalPotential::OpticalPotential(string configFileName)
{
    ifstream configFile(configFileName);

    if(!configFile.is_open())
    {
        cerr << "Error: failed to open optical potential config file "
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

        if(tokens[0]=="VHF")
        {
            VHF = stod(tokens[2]);
            continue;
        }

        if(tokens[0]=="VHF_Edep")
        {
            VHF_Edep = stod(tokens[2]);
            continue;
        }

        if(tokens[0]=="RHF")
        {
            RHF = stod(tokens[2]);
            continue;
        }

        if(tokens[0]=="RHF_Edep")
        {
            RHF_Edep = stod(tokens[2]);
            continue;
        }

        if(tokens[0]=="RHF_Eshift")
        {
            RHF_Eshift = stod(tokens[2]);
            continue;
        }

        if(tokens[0]=="aHF")
        {
            aHF = stod(tokens[2]);
            continue;
        }

        if(tokens[0]=="Avolume")
        {
            Avolume = stod(tokens[2]);
            continue;
        }

        if(tokens[0]=="Avolume_Edep")
        {
            Avolume_Edep = stod(tokens[2]);
            continue;
        }

        if(tokens[0]=="Avolume_Eshift")
        {
            Avolume_Eshift = stod(tokens[2]);
            continue;
        }

        if(tokens[0]=="Rvolume")
        {
            Rvolume = stod(tokens[2]);
            continue;
        }

        if(tokens[0]=="avolume")
        {
            avolume = stod(tokens[2]);
            continue;
        }

        if(tokens[0]=="Asurface")
        {
            Asurface = stod(tokens[2]);
            continue;
        }

        if(tokens[0]=="Asurface_Edep")
        {
            Asurface_Edep = stod(tokens[2]);
            continue;
        }

        if(tokens[0]=="Asurface_Eshift")
        {
            Asurface_Eshift = stod(tokens[2]);
            continue;
        }

        if(tokens[0]=="CoulombRadiusConstant")
        {
            CoulombRadiusConstant = stod(tokens[2]);
            continue;
        }
    }

    /*cout << "OP.VHF = " << OP.VHF << endl;
          cout << "OP.VHF_Edep = " << OP.VHF_Edep << endl;
          cout << "OP.RHF = " << OP.RHF << endl;
          cout << "OP.RHF_Edep = " << OP.RHF_Edep << endl;
          cout << "OP.RHF_Eshift = " << OP.RHF_Eshift << endl;
          cout << "OP.aHF = " << OP.aHF << endl;
          cout << "OP.Avolume = " << OP.Avolume << endl;
          cout << "OP.Avolume_Edep = " << OP.Avolume_Edep << endl;
          cout << "OP.Avolume_Eshift = " << OP.Avolume_Eshift << endl;
          cout << "OP.Rvolume = " << OP.Rvolume << endl;
          cout << "OP.avolume = " << OP.avolume << endl;
          cout << "OP.Asurface = " << OP.Asurface << endl;
          cout << "OP.Asurface_Edep = " << OP.Asurface_Edep << endl;
          cout << "OP.Asurface_Eshift = " << OP.Asurface_Eshift << endl;
          cout << "OP.CoulombRadiusConstant = " << OP.CoulombRadiusConstant << endl;
          */
}
