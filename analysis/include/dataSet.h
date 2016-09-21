#ifndef __DATASET_H__
#define __DATASET_H__

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

        TGraphErrors* getPlot() const;

        std::string getReference();
        void setReference(std::string reference);

        friend const DataSet operator+(const DataSet& set1, const DataSet& set2);
        friend const DataSet operator-(const DataSet& set1, const DataSet& set2);
        friend const DataSet operator/(const DataSet& set1, const DataSet& set2);

        const DataSet plus(const DataSet& set2, const std::string& name);
        const DataSet minus(const DataSet& set2, const std::string& name);
        const DataSet divideBy(const DataSet& set2, const std::string& name);

    private:
        void createPlot();

        TGraphErrors* dataPlot;

        std::string reference;

        std::vector<DataPoint> data;
        std::vector<double> energy;
        std::vector<double> xsection;
        std::vector<double> error;
};

#endif
