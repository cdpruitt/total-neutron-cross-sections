#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include "TMath.h"
#include "TFile.h"

#include "../include/dataPoint.h"
#include "../include/dataSet.h"

using namespace std;

DataSet::DataSet()
{
}

DataSet::DataSet(std::vector<double> var1, std::vector<double> var2, std::vector<double> var3, std::string ref)
{
    for(int i=0; (size_t)i<var1.size(); i++)
    {
        data.push_back(DataPoint(var1[i],0,var2[i],var3[i]));
    }

    reference = ref;

    //createPlot(reference);
}

DataSet::DataSet(string dataSetLocation)
{
    std::ifstream dataFile(dataSetLocation.c_str());
    if(!dataFile.is_open())
    {
        std::cout << "Attempted to create DataSet, but failed to find " << dataSetLocation << std::endl;
        exit(1);
    }

    string dummy;
    std::getline(dataFile,dummy);
    std::getline(dataFile,reference);
    std::getline(dataFile,dummy);

    double dum,dum2,dum3;

    //cout << reference;

    while(dataFile >> dum >> dum2 >> dum3)
    {
        data.push_back(DataPoint(dum, 0, dum2, dum3));
    }

    createPlot(reference);
    dataFile.close();
}

DataSet::DataSet(TGraphErrors* graph, string ref)
{

    int numberOfPoints = graph->GetN();

    vector<double> x (numberOfPoints);
    vector<double> y (numberOfPoints);
    vector<double> xError (numberOfPoints);
    vector<double> yError (numberOfPoints);

    for(int i=0; i<numberOfPoints; i++)
    {
        graph->GetPoint(i,x[i],y[i]);
        xError[i] = graph->GetErrorX(i);
        yError[i] = graph->GetErrorY(i);

        this->addPoint(DataPoint(x[i],xError[i],y[i],yError[i]));
    }

    reference = ref; 
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
    DataSet summedDataSet;
    if(set1.getNumberOfPoints() != set2.getNumberOfPoints())
    {
        cerr << "Error: tried to add data sets with different number of points. Returning empty data set..." << endl;
        return summedDataSet;
    }

    for (int i=0; i<set1.getNumberOfPoints(); i++)
    {
        if (set1.getPoint(i).getXValue()!=set2.getPoint(i).getXValue())
        {
            cerr << "Error: there was an energy mismatch between a data point in set 1\
                and a data point in set 2. Returning summed dataSet." << endl;
            return summedDataSet;
        }

        double xValue = set1.getPoint(i).getXValue();
        double xError = set1.getPoint(i).getXError();
        double yValue = set1.getPoint(i).getYValue() +
                        set2.getPoint(i).getYValue();
        double yError = pow(pow(set1.getPoint(i).getYError(),2) +
                            pow(set2.getPoint(i).getYError(),2),0.5);

        summedDataSet.addPoint(DataPoint(xValue,xError,yValue,yError));
    }

    summedDataSet.setReference(set1.getReference() + set2.getReference() + "sum");
    return summedDataSet;
}

const DataSet operator-(const DataSet& set1, const DataSet& set2)
{
    DataSet differenceDataSet;
    if(set1.getNumberOfPoints() != set2.getNumberOfPoints())
    {
        cerr << "Error: tried to subtract data sets with different number of points. Returning empty data set..." << endl;
        return differenceDataSet;
    }

    for (int i=0; i<set1.getNumberOfPoints(); i++)
    {
        if (set1.getPoint(i).getXValue()!=set2.getPoint(i).getXValue())
        {
            cerr << "Error: there was an energy mismatch between a data point in set 1\
                and a data point in set 2. Returning difference dataSet." << endl;
            return differenceDataSet;
        }

        differenceDataSet.addPoint(set1.getPoint(i)-set2.getPoint(i));
    }

    differenceDataSet.setReference(set1.getReference() + set2.getReference() + "difference");
    return differenceDataSet;
}

const DataSet operator*(const DataSet& set1, const DataSet& set2)
{
    DataSet productDataSet;
    if(set1.getNumberOfPoints() != set2.getNumberOfPoints())
    {
        cerr << "Error: tried to multiply data sets with different number of points. Returning empty data set..." << endl;
        return productDataSet;
    }

    for (int i=0; i<set1.getNumberOfPoints(); i++)
    {
        if (set1.getPoint(i).getXValue()!=set2.getPoint(i).getXValue())
        {
            cerr << "Error: there was an energy mismatch between a data point in set 1\
                and a data point in set 2. Returning product dataSet." << endl;
            return productDataSet;
        }

        double xValue = set1.getPoint(i).getXValue();
        double xError = set1.getPoint(i).getXError();
        double yValue = set1.getPoint(i).getYValue() *
                        set2.getPoint(i).getYValue();
        double yError = abs(yValue)*
                        pow(pow(set1.getPoint(i).getYError()/set1.getPoint(i).getYValue(),2) +
                            pow(set2.getPoint(i).getYError()/set2.getPoint(i).getYValue(),2),0.5);

        productDataSet.addPoint(DataPoint(xValue,xError,yValue,yError));
    }

    productDataSet.setReference(set1.getReference() + set2.getReference() + "product");
    return productDataSet;
}

const DataSet operator*(const DataSet& multiplicand, const double multiplier)
{
    DataSet product;

    for(int i=0; i<multiplicand.getNumberOfPoints(); i++)
    {
        product.addPoint(multiplicand.getPoint(i)*multiplier);
    }

    return product;
}

const DataSet operator/(const DataSet& set1, const DataSet& set2)
{
    DataSet quotientDataSet;
    if(set1.getNumberOfPoints() != set2.getNumberOfPoints())
    {
        cerr << "Error: tried to divide data sets with different number of points. Returning empty data set..." << endl;
        return quotientDataSet;
    }

    for (int i=0; i<set1.getNumberOfPoints(); i++)
    {
        if (set1.getPoint(i).getXValue()!=set2.getPoint(i).getXValue())
        {
            cerr << "Error: there was an energy mismatch between a data point in set 1\
                and a data point in set 2. Returning quotient dataSet." << endl;
            return quotientDataSet;
        }

        double xValue = set1.getPoint(i).getXValue();
        double xError = set1.getPoint(i).getXError();
        double yValue = set1.getPoint(i).getYValue() /
                        set2.getPoint(i).getYValue();
        double yError = abs(yValue)*
                        pow(pow(set1.getPoint(i).getYError()/set1.getPoint(i).getYValue(),2) +
                            pow(set2.getPoint(i).getYError()/set2.getPoint(i).getYValue(),2),0.5);

        quotientDataSet.addPoint(DataPoint(xValue,xError,yValue,yError));
    }

    quotientDataSet.setReference(set1.getReference() + set2.getReference() + "quotient");
    return quotientDataSet;
}

const DataSet correctCSUsingControl(const DataSet& dataSetToCorrect, const DataSet& correction)
{
    DataSet correctedDataSet;
    if(dataSetToCorrect.getNumberOfPoints() != correction.getNumberOfPoints())
    {
        cerr << "Error: tried to divide data sets with different number of points. Returning empty data set..." << endl;
        return correctedDataSet;
    }

    for (int i=0; i<dataSetToCorrect.getNumberOfPoints(); i++)
    {
        if (dataSetToCorrect.getPoint(i).getXValue()!=correction.getPoint(i).getXValue())
        {
            cerr << "Error: there was an energy mismatch between a data point in set 1\
                and a data point in set 2. Returning quotient dataSet." << endl;
            return correctedDataSet;
        }

        double xValue = dataSetToCorrect.getPoint(i).getXValue();
        double xError = dataSetToCorrect.getPoint(i).getXError();
        double yValue = dataSetToCorrect.getPoint(i).getYValue()/
                        (1+correction.getPoint(i).getYValue());
        double yError = dataSetToCorrect.getPoint(i).getYError()/
                        (1+correction.getPoint(i).getYValue());

        correctedDataSet.addPoint(DataPoint(xValue,xError,yValue,yError));
    }

    correctedDataSet.setReference(dataSetToCorrect.getReference() + correction.getReference() + "correctedByControl");
    return correctedDataSet;
}

const DataSet operator/(const DataSet& dividend, const double divisor)
{
    DataSet quotient;

    for(int i=0; i<dividend.getNumberOfPoints(); i++)
    {
        quotient.addPoint(dividend.getPoint(i)/divisor);
    }

    return quotient;
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
    dataPlot->SetNameTitle(name.c_str(),name.c_str());
    return dataPlot;
}
