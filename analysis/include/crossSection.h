#ifndef CROSS_SECTION_H
#define CROSS_SECTION_H

#include <vector>

#include "../include/dataPoint.h"
#include "../include/dataSet.h"
#include "../include/CSPrereqs.h"

class CrossSection
{
    public:
        CrossSection();
        void addDataPoint(DataPoint dataPoint);
        void addDataSet(DataSet dataSet);
        DataPoint getDataPoint(int i) const;
        DataSet getDataSet();
        void createCSGraph(std::string name);

        int getNumberOfPoints() const;
        std::vector<double> getEnergyValues() const;
        std::vector<double> getEnergyErrors() const;
        std::vector<double> getCrossSectionValues() const;
        std::vector<double> getCrossSectionErrors() const;

    private:
        DataSet dataSet;
};

int calculateCS(std::string CSFileName, CSPrereqs& targetData, CSPrereqs& blankData);
void correctForDeadtime(std::string histoFileName, std::string deadtimeFileName, std::string directory);

#endif
