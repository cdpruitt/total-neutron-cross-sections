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

DataSet::DataSet(string dataSetLocation, const vector<double>& energyBins)
{
    std::ifstream dataFile(dataSetLocation.c_str());
    if(!dataFile.is_open())
    {
        std::cout << "Attempted to create DataSet, but failed to find " << dataSetLocation << std::endl;
        exit(1);
    }

    vector<DataPoint> rawDataPoints;

    string dummy;
    std::getline(dataFile,dummy);
    std::getline(dataFile,reference);
    std::getline(dataFile,dummy);

    double dum,dum2,dum3;

    while(dataFile >> dum >> dum2 >> dum3)
    {
        rawDataPoints.push_back(DataPoint(dum, 0, dum2, dum3));
    }

    vector<vector<DataPoint>> binnedDataPoints;

    for(int i=0; i<energyBins.size()-1; i++)
    {
        binnedDataPoints.push_back(vector<DataPoint>());
    }

    int energyBinCounter = 0;
    double lowEnergy = energyBins[energyBinCounter];
    double highEnergy = energyBins[energyBinCounter+1];

    bool endSort = false;

    // average data in rawDataPoints according to energy bins
    for(int i=0; i<rawDataPoints.size(); i++)
    {
        DataPoint dp = rawDataPoints[i];

        if(dp.getXValue() < lowEnergy)
        {
            continue;
        }

        while(dp.getXValue() > highEnergy)
        {
            energyBinCounter++;
            if(energyBinCounter>=energyBins.size()-1)
            {
                endSort = true;
                break;
            }

            lowEnergy = energyBins[energyBinCounter];
            highEnergy = energyBins[energyBinCounter+1];
        }

        if(endSort)
        {
            break;
        }
        
        if(dp.getXValue()>=lowEnergy
                && dp.getXValue()<=highEnergy)
        {
            // found correct bin for this data point
            binnedDataPoints[energyBinCounter].push_back(dp);
        }

        else
        {
            cerr << "Error: assignment error while reading literature dataset. Exiting..." << endl;
            exit(1);
        }

        if(i%10==0)
        {
            cout << "Read " << i << " data points in from literature file...\r";
            fflush(stdout);
        }
    }

    for(int i=0; i<binnedDataPoints.size(); i++)
    {
        double averageCrossSection = 0;
        double averageError = 0;

        for(int j=0; j<binnedDataPoints[i].size(); j++)
        {
            averageCrossSection += binnedDataPoints[i][j].getYValue();
            averageError += binnedDataPoints[i][j].getYError();
        }

        if(binnedDataPoints[i].size()==0)
        {
            continue;
        }

        averageCrossSection /= binnedDataPoints[i].size();
        averageError /= binnedDataPoints[i].size();

        double energy = (energyBins[i] + energyBins[i+1])/2;

        data.push_back(DataPoint(energy, 0, averageCrossSection, averageError));
    }

    createPlot(reference);
    dataFile.close();
}

DataSet::DataSet(TGraphAsymmErrors* graph, string ref)
{
    if(!graph)
    {
        cerr << "Error: cannot create a dataset from a non-existent graph" << endl;
        exit(1);
    }

    int numberOfPoints = graph->GetN();

    vector<double> x (numberOfPoints, 0);
    vector<double> y (numberOfPoints, 0);
    vector<double> xError (numberOfPoints, 0);
    vector<double> yError (numberOfPoints, 0);

    for(int i=0; i<numberOfPoints; i++)
    {
        graph->GetPoint(i,x[i],y[i]);
        xError[i] = graph->GetErrorX(i);
        yError[i] = graph->GetErrorY(i);

        addPoint(DataPoint(x[i],xError[i],y[i],yError[i]));
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

TGraphAsymmErrors* DataSet::getPlot() const
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

const DataSet operator+(const DataSet& rawDataSet, const double addend)
{
    DataSet sum;

    for(int i=0; i<rawDataSet.getNumberOfPoints(); i++)
    {
        sum.addPoint(rawDataSet.getPoint(i)+addend);
    }

    return sum;
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
                        correction.getPoint(i).getYValue();
        double yError = dataSetToCorrect.getPoint(i).getYError()/
                        correction.getPoint(i).getYValue();

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

vector<double> DataSet::getXErrorsL() const
{
    vector<double> xErrorsL;
    for(const DataPoint point : data)
    {
        xErrorsL.push_back(point.getXErrorL());
    }
    return xErrorsL;
}

vector<double> DataSet::getXErrorsR() const
{
    vector<double> xErrorsR;
    for(const DataPoint point : data)
    {
        xErrorsR.push_back(point.getXErrorR());
    }
    return xErrorsR;
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

TGraphAsymmErrors* DataSet::createPlot(string name)
{
    dataPlot = new TGraphAsymmErrors(this->getXValues().size(),&this->getXValues()[0],&this->getYValues()[0],&this->getXErrors()[0],&this->getYErrors()[0]);
    dataPlot->SetNameTitle(name.c_str(),name.c_str());
    return dataPlot;
}
