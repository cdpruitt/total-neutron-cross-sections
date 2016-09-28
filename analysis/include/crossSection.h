#ifndef CROSS_SECTION_H
#define CROSS_SECTION_H

#include <vector>

#include "../include/dataPoint.h"

class CrossSection
{
    public:
        CrossSection();
        void addDataPoint(DataPoint dataPoint);
        DataPoint getDataPoint(int i) const;
        void createCSGraph(std::string name);

        int getNumberOfPoints() const;
        std::vector<double> getEnergyValues() const;
        std::vector<double> getEnergyErrors() const;
        std::vector<double> getCrossSectionValues() const;
        std::vector<double> getCrossSectionErrors() const;

    private:
        std::vector<DataPoint> data;
};

void calculateCS(std::string histoFileName, std::string CSFileName, int runNumber);


#endif
