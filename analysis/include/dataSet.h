#ifndef __DATASET_H__
#define __DATASET_H__

#include "TFile.h"
#include "TGraphErrors.h"
#include <sstream>
#include "../include/dataPoint.h"

class DataSet
{
    public:
        // constructors
        DataSet();
        DataSet(std::string dataSetLocation);
        DataSet(std::vector<double> energy, std::vector<double> xsection, std::vector<double> error, std::string reference);
        DataSet(TGraphErrors* graph, std::string reference);

        TGraphErrors* getPlot() const;
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

        friend const DataSet operator*(const DataSet& multiplicand, const double multiplier);
        friend const DataSet operator/(const DataSet& dividend, const double divisor);

        const DataSet plus(const DataSet& set2, const std::string& name);
        const DataSet minus(const DataSet& set2, const std::string& name);
        const DataSet divideBy(const DataSet& set2, const std::string& name);
        DataSet merge(DataSet set2);

        std::vector<double> getXValues() const;
        std::vector<double> getXErrors() const;
        std::vector<double> getYValues() const;
        std::vector<double> getYErrors() const;

        TGraphErrors* createPlot(std::string name);
        TGraphErrors* createPlot(std::string name, TFile* file);
    private:

        TGraphErrors* dataPlot;

        std::string reference;

        std::vector<DataPoint> data;
};

#endif
