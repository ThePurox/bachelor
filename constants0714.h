#ifndef CONSTANTS_H
#define CONSTANTS_H

extern ofstream& psiout;
extern ofstream& out2;
extern ofstream& densityout;
extern ofstream& density2file;
extern ofstream& anout;
extern ofstream& barycenterout;
extern ofstream& relativeout;
extern const char* relativefile;
extern const char* barycenterfile;
extern ofstream& currentout;
extern ofstream& vextout;
extern ofstream& superAdout;
extern ofstream& currentAdout;
extern ofstream& timeout;
extern ofstream& powerAdout;
extern ofstream& powerout;
extern ofstream& dissipationout;
extern ofstream& extPowerout;
extern ofstream& daniout;
extern ofstream& superCurrout;
extern ofstream& squaresout;
extern ofstream& velout;
extern ofstream& velAdout;
extern ofstream& superPowerout;

extern const double sigWWPot;
extern const double eps;
extern const double L;
extern const int N;
extern const double extForceFactor;
extern const double extPotFactor0;
extern long count;
extern const double T;
extern const double tFinal;
extern const int printNumb;
extern const int N_part;
extern const int overscan;
extern const int N_space;
extern const int N_PS;
extern long N_psi;
extern const double dt;
extern const int steps;
extern int NN;
extern double norm;
extern int forceInt;
extern int version;
extern int step;
extern int printStep;
extern double speedRPS;
extern double powderPS,powderNN;

extern double* psi;
extern double* initial;
extern double* anal;
extern double*** extForce;
extern double** M;// last collum is for alpha and beta and gamma in Numerical Recepies in C page 75
extern double* rhs;
extern double* Vext;
//extern double* WWPot;
extern double*** current;
extern double** density;
extern double* psiAd;
extern double*** currentAd;
extern double** freeEnergy;
extern double** power;
extern double** powerAd;
extern double** dissipation;
extern double** extPower;
extern double** daniA;
extern double** squares;
extern double*** vel;
extern double*** velAd;

extern double dx;
extern double alpha;
extern double sigVext0;
extern double initialNorm;
extern bool VextTimeDependent;

#endif
