#ifndef GAMMA_CORRECTION__H
#define GAMMA_CORRECTION__H

#include <vector>
#include <string>

struct GammaEvent
{
    GammaEvent() {}
    GammaEvent(double t, double e, double w) : time(t), energy(e), weight(w) {}
    double time = 0;
    double energy = 0;
    double weight = 0;
};

struct GammaCorrection
{
    std::vector<GammaEvent> gammaList;

    double averageGammaTime = 0;
    unsigned int numberOfGammas = 0;
    double correction = 0;
};

int calculateGammaCorrection(std::string inputFileName, std::string treeName,
        std::vector<GammaCorrection>& gammaCorrectionList, std::string outputFileName);

#endif /* GAMMA_CORRECTION_H */
