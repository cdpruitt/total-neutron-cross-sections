#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "TMath.h"

#include "../include/dataPoint.h"
#include "../include/dataSet.h"

using namespace std;

DataSet::DataSet()
{
}

DataSet::DataSet(std::vector<double> var1, std::vector<double> var2, std::vector<double> var3, std::string ref)
{
    energy = var1;
    xsection = var2;
    error = var3;

    reference = ref;

    createPlot(reference);
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
}

void DataSet::addPoint(DataPoint dataPoint)
{
    data.push_back(dataPoint);
}

DataPoint DataSet::getPoint(int i) const
{
    return data[i];
}

TGraphErrors* DataSet::getPlot() const
{
    return dataPlot;
}

std::string DataSet::getReference() const
{
    return reference;
}

void DataSet::setReference(string reference)
{
    this->reference = reference;
}

int DataSet::getNumberOfPoints() const
{
    return data.size();
}

const DataSet operator+(const DataSet& set1, const DataSet& set2)
{
    std::vector<double> sumEnergy;
    std::vector<double> sumXsection;
    std::vector<double> sumError;

    std::vector<double> energies1 = set1.getXValues();
    std::vector<double> xsecs1 = set1.getYValues();
    std::vector<double> errors1 = set1.getYErrors();

    std::vector<double> energies2 = set2.getXValues();
    std::vector<double> xsecs2 = set2.getYValues();
    std::vector<double> errors2 = set2.getYErrors();

    DataSet summedDataSet;
    for (int i=0; i<set2.getNumberOfPoints(); i++)
    {
        if (energies1[i] == energies2[i])
        {
            sumEnergy.push_back(energies1[i]);
            sumXsection.push_back(xsecs1[i] + xsecs2[i]);
            sumError.push_back(pow(pow(errors1[i],2) + pow(errors2[i],2),2));
        }

        else
        {
            cerr << "Error: set 1 energies != set 2 energies" << endl;
            exit(1);
        }
        DataPoint d = DataPoint(sumEnergy[i],0,sumXsection[i],sumError[i]);
        summedDataSet.addPoint(d);
    }

    summedDataSet.setReference(set1.getReference() + set2.getReference() + "sum");
    return summedDataSet;
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
    std::vector<double> diffEnergy;
    std::vector<double> diffXsection;
    std::vector<double> diffError;

    std::vector<double> energies1 = set1.getXValues();
    std::vector<double> xsecs1 = set1.getYValues();
    std::vector<double> errors1 = set1.getYErrors();

    std::vector<double> energies2 = set2.getXValues();
    std::vector<double> xsecs2 = set2.getYValues();
    std::vector<double> errors2 = set2.getYErrors();

    DataSet diffDataSet;
    for (int i=0; i<set2.getNumberOfPoints(); i++)
    {
        if (energies1[i] == energies2[i])
        {
            diffEnergy.push_back(energies1[i]);
            diffXsection.push_back(xsecs1[i] - xsecs2[i]);
            diffError.push_back(pow(pow(errors1[i],2) + pow(errors2[i],2),2));
        }

        else
        {
            cerr << "Error: set 1 energies != set 2 energies" << endl;
            exit(1);
        }
        DataPoint d = DataPoint(diffEnergy[i],0,diffXsection[i],diffError[i]);
        diffDataSet.addPoint(d);
    }

    diffDataSet.setReference(set1.getReference() + set2.getReference() + "diff");
    return diffDataSet;
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

const DataSet operator*(const DataSet& set1, const DataSet& set2)
{
    std::vector<double> multEnergy;
    std::vector<double> multXsection;
    std::vector<double> multError;

    std::vector<double> energies1 = set1.getXValues();
    std::vector<double> xsecs1 = set1.getYValues();
    std::vector<double> errors1 = set1.getYErrors();

    std::vector<double> energies2 = set2.getXValues();
    std::vector<double> xsecs2 = set2.getYValues();
    std::vector<double> errors2 = set2.getYErrors();

    DataSet multDataSet;
    for (int i=0; i<set2.getNumberOfPoints(); i++)
    {
        if (energies1[i] == energies2[i])
        {
            multEnergy.push_back(energies1[i]);
            if(xsecs1[i]!=0 && xsecs2[i]!=0)
            {
                multXsection.push_back(xsecs1[i]*xsecs2[i]);
                multError.push_back(pow(pow(errors1[i]/xsecs1[i],2) + pow(errors2[i]/xsecs2[i],2),2));
            }

            else
            {
                multXsection.push_back(0);
                multError.push_back(0);
            }
        }

        else
        {
            cerr << "Error: set 1 energies != set 2 energies" << endl;
            exit(1);
        }
        DataPoint d = DataPoint(multEnergy[i],0,multXsection[i],multError[i]);
        multDataSet.addPoint(d);
    }

    multDataSet.setReference(set1.getReference() + set2.getReference() + "mult");
    return multDataSet;
}

const DataSet operator/(const DataSet& set1, const DataSet& set2)
{
    std::vector<double> divEnergy;
    std::vector<double> divXsection;
    std::vector<double> divError;

    std::vector<double> energies1 = set1.getXValues();
    std::vector<double> xsecs1 = set1.getYValues();
    std::vector<double> errors1 = set1.getYErrors();

    std::vector<double> energies2 = set2.getXValues();
    std::vector<double> xsecs2 = set2.getYValues();
    std::vector<double> errors2 = set2.getYErrors();

    DataSet divDataSet;
    for (int i=0; i<set2.getNumberOfPoints(); i++)
    {
        if (energies1[i] == energies2[i])
        {
            divEnergy.push_back(energies1[i]);
            if(xsecs1[i]!=0 && xsecs2[i]!=0)
            {
                divXsection.push_back(xsecs1[i]/xsecs2[i]);
                divError.push_back(pow(pow(errors1[i]/xsecs1[i],2) + pow(errors2[i]/xsecs2[i],2),2));
            }

            else
            {
                divXsection.push_back(0);
                divError.push_back(0);
            }
        }

        else
        {
            cerr << "Error: set 1 energies != set 2 energies" << endl;
            exit(1);
        }
        DataPoint d = DataPoint(divEnergy[i],0,divXsection[i],divError[i]);
        divDataSet.addPoint(d);
    }

    divDataSet.setReference(set1.getReference() + set2.getReference() + "div");
    return divDataSet;
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

DataSet DataSet::merge(DataSet set2)
{
    DataSet mergedSet;
    for(int i=0; (size_t)i<data.size(); i++)
    {
        DataPoint point1 = this->getPoint(i);
        DataPoint point2 = set2.getPoint(i);
        mergedSet.addPoint(point1.mergePoints(point2));
    }
    return mergedSet;
}

vector<double> DataSet::getXValues() const
{
    vector<double> xValues;
    for(const DataPoint point : data)
    {
        xValues.push_back(point.getXValue());
    }
    return xValues;
}

vector<double> DataSet::getYValues() const
{
    vector<double> yValues;
    for(const DataPoint point : data)
    {
        yValues.push_back(point.getYValue());
    }
    return yValues;
}

vector<double> DataSet::getXErrors() const
{
    vector<double> xErrors;
    for(const DataPoint point : data)
    {
        xErrors.push_back(point.getXError());
    }
    return xErrors;
}

vector<double> DataSet::getYErrors() const
{
    vector<double> yErrors;
    for(const DataPoint point : data)
    {
        yErrors.push_back(point.getYError());
    }
    return yErrors;
}

TGraphErrors* DataSet::createPlot(string name)
{
    dataPlot = new TGraphErrors(this->getXValues().size(),&this->getXValues()[0],&this->getYValues()[0],&this->getXErrors()[0],&this->getYErrors()[0]);
    dataPlot->SetName(name.c_str());
    dataPlot->Write();
    return dataPlot;
}
