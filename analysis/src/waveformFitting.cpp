#include <iostream>
#include <math.h>

using namespace std;

// For functions that use a Fermi function, define the x-axis offset
const double FERMI_OFFSET = 11.5; // in ns

/*****************************************************************************/
/* Define the peak-fitting functional form here:
 *
 *  onePeakForm = A * (t-t0)^n * exp[-((t-t0)^1)/w] + [C + m*(t-t0)]
 *       + B * (t-t1)^n * exp[-((t-t1)^1)/w]
 *
 *  par[0] = A = amplitude of first peak
 *  par[1] = B = amplitude of second peak
 *  par[2] = t1 = starting time of first peak
 *  par[3] = t2 = starting time of second peak
 *  par[4] = n = polynomial order of peak (rise rate)
 *  par[5] = w = exponential decay constant of peak (fall rate)
 *  par[6] = C = baseline offset
 *  par[7] = m = slope of baseline
 *
 *  During the first attempt at fitting, set B = 0 and allow A, t1, C, and m
 *  to vary. If this produces a satisfactory fit (chi squared < ERROR_LIMIT),
 *  then we end fitting and accept the fit's trigger time, amplitude, etc for
 *  this peak.
 *  If this failes to produce a satisfactory fit, then we reset the parameters
 *  to their initial states and allow B and t2 to vary as well (i.e., allow
 *  two peaks in the fit). If this fitting also fails, we skip the peak region
 *  and indicate an error.

double onePeakForm(double *x, double *par)
{
    double fitval;    // fitted value of function
    double arg1 = 0;   // argument of exponential 1
    double arg2 = 0;   // argument of exponential 2

    if (par[5]!=0)
    {
        arg1 = pow((x[0]-par[2]),2)/par[5];
        arg2 = pow((x[0]-par[3]),2)/par[5];
    }

    fitval = par[6] + par[7]*(x[0]-par[2]);
    fitval += par[0]*pow((x[0]-par[2]),par[4])*exp(-arg1);
    fitval += par[1]*pow((x[0]-par[3]),par[4])*exp(-arg2);

    // If before first peak start, set to background + peak's gaussian
    if (x[0]<par[2])
    {
        fitval = par[6] + par[7]*(x[0]-par[2]);
    }

    // If before second peak start, set to background + first peak
    if (x[0]<par[3] && x[0]>=par[2])
    {
        fitval = par[6] + par[7]*(x[0]-par[2]);
        fitval += par[0]*pow((x[0]-par[2]),par[4])*exp(-arg1);
    }

    return fitval;
}
*/

/*****************************************************************************/
/* The peak-fitting functional form is defined as a linear background plus
 * (up to) two peaks, with each peak a convolution of a Maxwell-Boltzmann and a
 * Gaussian.
 *
 * Thus, convolutedPeakFunc =
 *   Integral (A * (i-t1)^n * exp[-((i-t1)^2)/d + (t-i-t1)^2/(2w^2)]) di
 * + Integral (B * (i-t2)^n * exp[-((i-t2)^2)/d + (t-i-t2)^2/(2w^2)]) di
 * + C + m*(t-t0)
 *
 * This form is realized as a TF1Convolution of a user-defined function,
 * fittingFunc, and a generic Gaussian provided by ROOT.
 *
 * fittingFunc = Maxwell-Boltzmann-like distribution + background
 *          = A * (t-t1)^n * exp[-((t-t1)^2)/d] + C + m*(t-t1)
 *
 * Two fittingFuncs are present in the final convolutedPeakFunc expression,
 * and they each have independent amplitudes and time zeroes.
 *
 *  par[0] = A = amplitude of first peak
 *  par[1] = t1 = starting time of first peak
 *  par[2] = n = polynomial order of peak (rise rate)
 *  par[3] = d = exponential decay constant of peak (fall rate)
 *  par[4] = C = baseline offset
 *  par[5] = m = slope of baseline
 *
 *  During the first attempt at fitting, set B = 0 and allow A, t1, C, and m
 *  to vary. If this produces a satisfactory fit (chi squared < ERROR_LIMIT),
 *  then we end fitting and accept the fit's trigger time, amplitude, etc for
 *  this peak.
 *  If this fails to produce a satisfactory fit, allow B and t2 to vary as well
 *  (i.e., allow two peaks in the fit). If this fitting also fails, we ignore
 *  this trigger.
 */



double onePeakForm(double *x, double *par)
{
    double fitval;    // fitted value of function

    if (par[2]==0 || par[3]==0)
    {
        cout << "Error: divide by 0 in onePeakForm" << endl;
        exit(1);
    }

    // Define exponential*fermi function as peak shape
    double arg1 = pow((x[0]-par[1]),1)/par[2];
    double arg2 = pow((x[0]-(par[1]+FERMI_OFFSET)),1)/par[3];

    // Start with linear background
    fitval = par[4] + par[5]*(x[0]-par[1]);

    // add peak
    fitval += par[0]*exp(-arg1)/(1+exp(-arg2));

    return fitval;
}

double onePeakExpBackForm(double *x, double *par)
{
    double fitval;    // fitted value of function

    if (par[2]==0 || par[3]==0)
    {
        cout << "Error: divide by 0 in onePeakOnExpForm" << endl;
        exit(1);
    }

    // Define exponential*fermi function as peak shape
    double arg1 = pow((x[0]-par[1]),1)/par[2];
    double arg2 = pow((x[0]-(par[1]+FERMI_OFFSET)),1)/par[3];

    // Define exponential decay as previous peak background
    double arg3 = pow((x[0]-(par[7])),1)/par[3];

    // Start with linear background
    fitval = par[4] + par[5]*(x[0]-par[1]);

    // Add exponential background of previous peak
    fitval += par[6]*exp(-arg3);

    // add peak
    fitval += par[0]*exp(-arg1)/(1+exp(-arg2));

    return fitval;
}

double twoPeakForm(double *x, double *par)
{
    double fitval;    // fitted value of function

    if (par[4]==0 || par[5]==0)
    {
        cout << "Error: divide by 0 in twoPeakForm" << endl;
        exit(1);
    }

    // Define exponential*fermi function as peak shape
    double arg1 = pow((x[0]-par[2]),1)/par[4];
    double arg2 = pow((x[0]-(par[2]+FERMI_OFFSET)),1)/par[5];

    double arg3 = pow((x[0]-par[3]),1)/par[4];
    double arg4 = pow((x[0]-(par[3]+FERMI_OFFSET)),1)/par[5];


    // Start with linear background
    fitval = par[6] + par[7]*(x[0]-par[2]);

    // add peaks
    fitval += par[0]*exp(-arg1)/(1+exp(-arg2));
    fitval += par[1]*exp(-arg3)/(1+exp(-arg4));

    return fitval;
}

/*double detForm(double *x, double *par)
{
    double fitval;    // fitted value of function
    double arg = 0;   // argument of exponential

    // Uncomment for 'real' detForm
    if (par[2]==0)
    {
        cout << "Error: divide by 0 in detForm" << endl;
        exit(1);
    }

    fitval = par[0]*exp(-pow((x[0]-par[1]),2)/(2*pow(par[2],2)));

    return fitval;
}*/

/*double cfdForm(double *x, double *par)
{
    double fitval;    // fitted value of function

    // uncomment for gaussian peakform definition
    //fitval = par[0]*exp(-pow((x[0]-par[1]),2)/(2*pow(par[2],2)));

    // uncomment for 'correct' peakform definition
    double arg1 = 0;   // argument of peak 1 exponential
    double arg2 = 0;   // argument of peak 2 exponential

    if (par[5]==0)
    {
        cout << "Error: divide by 0 in onePeakForm" << endl;
        exit(1);
    }

    arg1 = pow((x[0]-par[2]),1)/par[4];
    arg2 = pow((x[0]-(par[2]+11.5)),1)/par[5];

    // define background
    fitval = par[6] + par[7]*(x[0]-par[2]);

    // add scaled-down part of CFD form
    fitval += (par[0]/CFD_SCALEDOWN)*exp(-arg1)/(1+par[1]*exp(-arg2));

    // add delayed part of CFD form
    arg1 = pow((x[0]-(par[2]+CFD_DELAY)),1)/par[4];
    arg2 = pow((x[0]-(par[2]+11.5+CFD_DELAY)),1)/par[5];
    fitval += -par[0]*exp(-arg1)/(1+par[1]*exp(-arg2));

    return fitval;
}*/
