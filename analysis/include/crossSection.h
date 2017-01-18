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

        double getArealDensity() const;
        void setArealDensity(double arealDensity);

        friend CrossSection operator+(const CrossSection& augend, const CrossSection& addend);
        friend CrossSection operator-(const CrossSection& minuend, const CrossSection& subtrahend);
        friend CrossSection operator/(const CrossSection& dividend, const CrossSection& divisor);

        double calculateRMSError();

        void subtractCS(std::string subtrahendFileName,
                        std::string subtrahendGraphName, double factor);

    private:
        DataSet dataSet;
        double arealDensity;
};

CrossSection calculateCS(std::string CSFileName, CSPrereqs& targetData, CSPrereqs& blankData);
void correctForDeadtime(std::string histoFileName, std::string deadtimeFileName, std::vector<std::string> detectorChannels);

double getPartialError(DataPoint aPoint, DataPoint bPoint, double aArealDensity);
CrossSection calculateRelative(CrossSection a, CrossSection b);

 CrossSection subtractCS(std::string rawCSFileName, std::string rawCSGraphName,
                        std::string subtrahendFileName, std::string subtrahendGraphName,
                        double factor,  // multiplies the subtrahend
                        double divisor, // divides the final difference
                        std::string name // name given to output graph
                       );

#endif
