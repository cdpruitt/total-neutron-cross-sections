#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "TMath.h"

#include "../include/dataSet.h"

using namespace std;

DataSet::DataSet()
{
    cerr << "Error: attempted to create DataSet without providing an input file or input data." << std::endl;
    exit(1);
}

DataSet::DataSet(std::vector<double> var1, std::vector<double> var2, std::vector<double> var3, std::string ref)
{
    energy = var1;
    xsection = var2;
    error = var3;

    reference = ref;

    createPlot();
}

DataSet::DataSet(string dataSetLocation)
{
    std::ifstream dataFile(dataSetLocation.c_str());
    if(!dataFile.is_open())
    {
        std::cout << "Attempted to create DataSet, but failed to find " << dataSetLocation << std::endl;
        exit(1);
    }

    char dummy[200];
    dataFile.getline(dummy,200);
    getline(dataFile,reference);
    dataFile.getline(dummy,200);

    double dum,dum2,dum3;

    while(dataFile >> dum >> dum2 >> dum3)
    {
        energy.push_back(dum);
        xsection.push_back(dum2);
        error.push_back(dum3);
    }

    createPlot();
}

TGraphErrors* DataSet::getPlot() const
{
    return dataPlot;
}

std::string DataSet::getReference()
{
    return reference;
}

const DataSet operator+(const DataSet& set1, const DataSet& set2)
{
    std::vector<double> sumEnergy;
    std::vector<double> sumXsection;
    std::vector<double> sumError;

    sumEnergy = set1.energy;
    sumXsection = set1.xsection;
    sumError = set1.error;

    for (int i=0; (size_t)i<set2.energy.size(); i++)
    {
        if (set1.energy[i] == set2.energy[i])
        {
            sumXsection[i] += set2.xsection[i];
            sumError[i] += set2.error[i];
        }

        else
        {
            cerr << "Error: set 1 energies != set 2 energies" << endl;
            exit(1);
        }
    }

    return DataSet(sumEnergy, sumXsection, sumError, set1.reference + set2.reference);
}

const DataSet DataSet::plus(const DataSet& set2,const string& name)
{
    std::vector<double> sumEnergy;
    std::vector<double> sumXsection;
    std::vector<double> sumError;

    sumEnergy = energy;
    sumXsection = xsection;
    sumError = error;

    for (int i=0; (size_t)i<energy.size(); i++)
    {
        if (energy[i] == set2.energy[i])
        {
            sumXsection[i] += set2.xsection[i];
            sumError[i] += set2.error[i];
        }

        else
        {
            cerr << "Error: set 1 energies != set 2 energies" << endl;
            exit(1);
        }
    }

    return DataSet(sumEnergy, sumXsection, sumError, name);
}

const DataSet operator-(const DataSet& set1, const DataSet& set2)
{
    std::vector<double> sumEnergy;
    std::vector<double> sumXsection;
    std::vector<double> sumError;

    sumEnergy = set1.energy;
    sumXsection = set1.xsection;
    sumError = set1.error;

    for (int i=0; (size_t)i<set2.energy.size(); i++)
    {
        if (set1.energy[i] == set2.energy[i])
        {
            sumXsection[i] -= set2.xsection[i];
            sumError[i] -= set2.error[i];
        }

        else
        {
            cerr << "Error: set 1 energies != set 2 energies" << endl;
            exit(1);
        }
    }

    return DataSet(sumEnergy, sumXsection, sumError, set1.reference + set2.reference);
}

const DataSet DataSet::minus(const DataSet& set2,const string& name)
{
    std::vector<double> sumEnergy;
    std::vector<double> sumXsection;
    std::vector<double> sumError;

    sumEnergy = energy;
    sumXsection = xsection;
    sumError = error;

    for (int i=0; (size_t)i<set2.energy.size(); i++)
    {
        if (energy[i] == set2.energy[i])
        {
            sumXsection[i] -= set2.xsection[i];
            sumError[i] -= set2.error[i];
        }

        else
        {
            cerr << "Error: set 1 energies != set 2 energies" << endl;
            exit(1);
        }
    }

    return DataSet(sumEnergy, sumXsection, sumError, name);
}

const DataSet operator/(const DataSet& set1, const DataSet& set2)
{
    std::vector<double> sumEnergy;
    std::vector<double> sumXsection;
    std::vector<double> sumError;

    sumEnergy = set1.energy;
    sumXsection = set1.xsection;
    sumError = set1.error;

    for (int i=0; (size_t)i<set2.energy.size(); i++)
    {
        if (set1.energy[i] == set2.energy[i])
        {
            if(set2.xsection[i] != 0)
            {
                sumXsection[i] /= set2.xsection[i];
            }

            if(set2.error[i] != 0)
            {
                sumError[i] /= set2.error[i];
            }
        }

        else
        {
            cerr << "Error: set 1 energies != set 2 energies" << endl;
            exit(1);
        }
    }
    return DataSet(sumEnergy, sumXsection, sumError, set1.reference + "_" + set2.reference);
}

const DataSet DataSet::divideBy(const DataSet& set2, const string& name)
{
    std::vector<double> sumEnergy;
    std::vector<double> sumXsection;
    std::vector<double> sumError;

    sumEnergy = energy;
    sumXsection = xsection;
    sumError = error;

    for (int i=0; (size_t)i<set2.energy.size(); i++)
    {
        if (energy[i] == set2.energy[i])
        {
            if(set2.xsection[i] != 0)
            {
                sumXsection[i] /= set2.xsection[i];
            }

            if(set2.error[i] != 0)
            {
                sumError[i] /= set2.error[i];
            }
        }

        else
        {
            cerr << "Error: set 1 energies != set 2 energies" << endl;
            exit(1);
        }
    }

    return DataSet(sumEnergy, sumXsection, sumError, name);
}

void DataSet::createPlot()
{
    dataPlot = new TGraphErrors(energy.size(),&energy[0],&xsection[0],0,&error[0]);
    dataPlot->SetName(reference.c_str());
    dataPlot->Write();
}
