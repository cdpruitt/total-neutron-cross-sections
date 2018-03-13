#ifndef __DATASET_H__
#define __DATASET_H__

#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include <sstream>
#include "../include/dataPoint.h"

class DataSet
{
    public:
        // constructors
        DataSet();
        DataSet(std::string dataSetLocation);
        DataSet(std::string dataSetLocation, const std::vector<double>& energyBins);
        DataSet(std::vector<double> energy, std::vector<double> xsection, std::vector<double> error, std::string reference);
        DataSet(TGraphAsymmErrors* graph, std::string reference);

        TGraphAsymmErrors* getPlot() const;
        int getNumberOfPoints() const;

        std::string getReference() const;
        void setReference(std::string reference);

        void addPoint(DataPoint dataPoint);
        DataPoint getPoint(int i) const;

        friend const DataSet operator+(const DataSet& set1, const DataSet& set2);
        friend const DataSet operator-(const DataSet& set1, const DataSet& set2);
        friend const DataSet operator*(const DataSet& set1, const DataSet& set2);
        friend const DataSet operator/(const DataSet& set1, const DataSet& set2);
        friend const DataSet correctCSUsingControl(const DataSet& dataSetToCorrect, const DataSet& correction);

        friend const DataSet operator+(const DataSet& augend, const double addend);
        friend const DataSet operator*(const DataSet& multiplicand, const double multiplier);
        friend const DataSet operator/(const DataSet& dividend, const double divisor);

        const DataSet plus(const DataSet& set2, const std::string& name);
        const DataSet minus(const DataSet& set2, const std::string& name);
        const DataSet divideBy(const DataSet& set2, const std::string& name);
        DataSet merge(DataSet set2);

        std::vector<double> getXValues() const;
        std::vector<double> getXErrors() const;
        std::vector<double> getXErrorsL() const;
        std::vector<double> getXErrorsR() const;

        std::vector<double> getYValues() const;
        std::vector<double> getYErrors() const;

        TGraphAsymmErrors* createPlot(std::string name);
        TGraphAsymmErrors* createPlot(std::string name, TFile* file);
    private:

        TGraphAsymmErrors* dataPlot;

        std::string reference;

        std::vector<DataPoint> data;
};

#endif
