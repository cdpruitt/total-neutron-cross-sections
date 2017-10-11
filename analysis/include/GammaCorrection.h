#ifndef GAMMA_CORRECTION__H
#define GAMMA_CORRECTION__H

#include <vector>
#include <string>

struct GammaCorrection
{
    unsigned int numberOfGammas = 0;
    double averageGammaTime = 0;
};

int calculateGammaCorrection(std::string inputFileName, std::string treeName,
        std::vector<GammaCorrection>& gammaCorrectionList);

#endif /* GAMMA_CORRECTION_H */
