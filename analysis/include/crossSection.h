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
        void createCSGraph(std::string name, std::string title);

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

        void subtractCS(std::string subtrahendFileName,
                        std::string subtrahendGraphName, double factor);

    private:
        DataSet dataSet;
        double arealDensity;
};

void correctForDeadtime(std::string histoFileName, std::string deadtimeFileName, std::vector<std::string> detectorChannels);

double getPartialError(DataPoint aPoint, DataPoint bPoint, double aArealDensity);
CrossSection calculateRelative(CrossSection a, CrossSection b);

CrossSection subtractCS(std::string rawCSFileName, std::string rawCSGraphName,
                        std::string subtrahendFileName, std::string subtrahendGraphName,
                        double factor,  // multiplies the subtrahend
                        double divisor, // divides the final difference
                        std::string name // name given to output graph
                       );

CrossSection multiplyCS(std::string rawCSFileName, std::string rawCSGraphName,
                        double factor,  // multiplies the subtrahend
                        std::string name // name given to output graph
                       );

CrossSection relativeCS(std::string rawCSFileName, std::string rawCSGraphName,
                        std::string subtrahendFileName, std::string subtrahendGraphName,
                        std::string name // name given to output graph
                       );

CrossSection correctForBlank(CrossSection rawCS, double targetNumberDensity, std::string expName, std::string graphFileName);

CrossSection calculateCS(const CSPrereqs& targetData, const CSPrereqs& blankData);

void applyCSCorrectionFactor(std::string CSCorrectionFilename, std::string CSCorrectionGraphName, std::string CSToBeCorrectedFilename, std::string CSToBeCorrectedGraphName, std::string outputFileName, std::string outputGraphName);

void scaledownCS(std::string CSToBeCorrectedFileName, std::string CSToBeCorrectedGraphName, int scaledown, std::string outputFileName, std::string outputGraphName);

void produceRunningRMS(DataSet firstDS, DataSet secondDS, std::string name);

#endif
