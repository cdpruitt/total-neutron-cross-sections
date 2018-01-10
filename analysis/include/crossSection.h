#ifndef CROSS_SECTION_H
#define CROSS_SECTION_H

#include <vector>
#include <string>

#include "../include/dataPoint.h"
#include "../include/dataSet.h"
#include "../include/CSPrereqs.h"

class CrossSection
{
    public:
        CrossSection();
        CrossSection(std::string name);

        void addDataPoint(DataPoint dataPoint);
        void addDataSet(DataSet dataSet);
        DataPoint getDataPoint(int i) const;
        DataSet getDataSet();
        void createGraph(std::string name, std::string title);

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
        friend CrossSection operator*(const CrossSection& cs, double factor);

        double calculateRMSError();

        void calculateCS(const CSPrereqs& targetData, const CSPrereqs& blankData);

        std::string name;

    private:
        DataSet dataSet;
        double arealDensity;
};

double getPartialError(DataPoint aPoint, DataPoint bPoint, double aArealDensity);
#endif
