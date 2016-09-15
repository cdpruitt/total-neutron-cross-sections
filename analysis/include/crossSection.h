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

        int getNumberOfPoints();
        std::vector<double> getEnergyValues();
        std::vector<double> getEnergyErrors();
        std::vector<double> getCrossSectionValues();
        std::vector<double> getCrossSectionErrors();

    private:
        std::vector<DataPoint> data;
};

#endif
