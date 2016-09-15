#ifndef HISTOS_H
#define HISTOS_H

#include "plots.h"

int histos(std::string sortedFileName, TFile*& histoFile, std::string crossSectionFileName, std::vector<Plots*>& plots);

#endif
