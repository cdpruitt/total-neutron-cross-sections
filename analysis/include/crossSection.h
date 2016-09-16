#ifndef CROSS_SECTION_H
#define CROSS_SECTION_H

#include <vector>

#include "DataPoint.h"

class CrossSection
{
    public:
        CrossSection();
        void addDataPoint(DataPoint dataPoint);
        DataPoint getDataPoint(int i);
        void createCSGraph();

        int getNumberOfPoints();
        std::vector<double> getEnergyValues();
        std::vector<double> getEnergyErrors();
        std::vector<double> getCrossSectionValues();
        std::vector<double> getCrossSectionErrors();

    private:
        std::vector<DataPoint> data;
};

void calculateCS(const std::vector<std::string>& targetOrder, std::string histoFileName, std::string CSFileName);


#endif
