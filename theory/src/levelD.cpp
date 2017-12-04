#include "../include/levelD.h"

/**
 * Constructor
 \param file is the name of the file (without *.ld extension) containg the input data
 */
levelD::levelD(string filename)
{
    string name = filename+".ld";
    ifstream file(name.c_str());
    if (!file)
    {
        cout << " could not open file " << name << " in levelD.cpp" << endl;
        hasLevelD = false;
        return;
    }

    hasLevelD = true;

    getline(file,name);
    getline(file,name);
    getline(file,name);
    float stuff1,stuff2,stuff3,stuff4;
    for (int i=0;i<nE;i++)
    {
        file >> Ex0[i]>> stuff1 >> stuff2 >> stuff3 >> stuff4;
        for (int j=0;j<nJ;j++) file >> ld[i][j];
        getline(file,name);
    }
    file.close();
    file.clear();

}
/**
 * returns the spin-dependent level density
 \param Ex is the excitation energy in MeV
 \param J is the spin (odd mass drop the 1/2 and give just integer value)
 */
double levelD::getLD(double Ex, int J)
{
    if(!hasLevelD)
    {
        return 0;
    }

    if (J > 29) return 0.;
    if (Ex > Ex0[nE-1]) return 0.;
    if (Ex < Ex0[0]) return 0.;

    int i = 1;
    for (;;)
    {
        if (Ex < Ex0[i]) break;
        i++;
    }

    if (ld[i-1][J] == 0.) return  (Ex-Ex0[i-1])/(Ex0[i]-Ex0[i-1])*ld[i][J];
    //otherwise extrapolate in a log plot.
    double out =  log(ld[i-1][J]) + (Ex-Ex0[i-1])/(Ex0[i]-Ex0[i-1])
        *(log(ld[i][J])-log(ld[i-1][J]));
    return exp(out);  
}
