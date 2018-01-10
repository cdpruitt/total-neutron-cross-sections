#ifndef CS_UTILITIES_H
#define CS_UTILITIES_H

#include <string>
#include <vector>

#include "../include/crossSection.h"
#include "../include/CSPrereqs.h"

CrossSection calculateRelativeDiff(CrossSection a, CrossSection b);

int producePlots(
        std::string dataLocation,
        std::string expName,
        std::vector<CrossSection>& crossSections);

CrossSection correctForBlank(CrossSection rawCS, double targetNumberDensity, std::string expName, std::string graphFileName);

CrossSection calculateRelative(CrossSection a, CrossSection b);

CrossSection subtractCS(std::string rawCSFileName, std::string rawCSGraphName,
                        std::string subtrahendFileName, std::string subtrahendGraphName,
                        double factor,  // multiplies the subtrahend
                        double divisor, // divides the final difference
                        std::string name // name given to output graph
                       );

CrossSection shiftCS(std::string rawCSFileName, std::string rawCSGraphName,
                        double shift,  // shift to be added, in barns
                        std::string outputFileName, std::string outputGraphName // name given to output graph
                       );


CrossSection multiplyCS(std::string rawCSFileName, std::string rawCSGraphName,
                        double factor,  // multiplies the subtrahend
                        std::string name // name given to output graph
                       );

CrossSection relativeCS(std::string rawCSFileName, std::string rawCSGraphName,
                        std::string subtrahendFileName, std::string subtrahendGraphName,
                        std::string outputFileName, std::string name // name given to output graph
                       );

CrossSection relativeDiffCS(std::string rawCSFileName, std::string rawCSGraphName,
                        std::string subtrahendFileName, std::string subtrahendGraphName,
                        std::string outputFileName, std::string name // name given to output graph
                       );

CrossSection calculateCS(const CSPrereqs& targetData, const CSPrereqs& blankData);

void applyCSCorrectionFactor(std::string CSCorrectionFilename, std::string CSCorrectionGraphName, std::string CSToBeCorrectedFilename, std::string CSToBeCorrectedGraphName, std::string outputFileName, std::string outputGraphName);

void scaledownCS(std::string CSToBeCorrectedFileName, std::string CSToBeCorrectedGraphName, int scaledown, std::string outputFileName, std::string outputGraphName);

void produceRunningRMS(DataSet firstDS, DataSet secondDS, std::string name);

int produceTotalCSPlots(std::string dataLocation, std::vector<CrossSection>& crossSections);

#endif /* CS_UTILITIES_H */
