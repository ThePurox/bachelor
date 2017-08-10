#ifndef FUNCTIONS_H
#define FUNCTION_H

int LinInd(int a[]);
int IndNN(int x, int y);
void combi(int a[]);
void calcF(double t);
double calcWWF(int a[],int select,bool pot);
double getInitial(int a[]);
void initArrays();
void init();
double PMSDiff(int a[],int select);
void fillMDiff();
void fillrhsDiff(int a[],int select);
void fillMDrift(int a[],int select);
void fillrhsDrift(int a[],int select);
void tridag( double r[],double x[]);
void tridag(double r[],int a[],int select);
void cicle(int a[],int select);
double intPsi(double func[]);
void printPsi(double func[],int n,ofstream &out);
void printDensity2(int dim1,int dim2);
void calcRelative(double* psi, double relaA[],double baryA[]);
void analyticfactor(double t);
void density1(double* psi,double rho[],int select);
void density2(double* psi,double rho[],int dim1, int dim2);
void densitySum(double* psi,double rho[]);
double sumCurrent(int posx, int posy, int dir);
double residuum(double a[],double b[],int n,int dim);
void printWWF();
void printCurrent(double*** current,ofstream &currentout);
double diffPsi(double* psi,int a[],int select);
void calcCurrent(double* psi,double*** current);
void doStep(double t);
void doSuperAd();
void calcFreeEnergy();
double diffVext(int a[],int part,int space);
int factorial(int n);
double residuumFreeEnergy();
void calcPower(double* psi,double**vSquared);
void printTime(double** arr, ofstream &file);
void calcExtPower();
double calcdtVext(int x, int y, double t);
double diffDens(int x,int y, int dir);
void dani();
void printSuperCurr(ofstream &superCurrout);
void calcVel(double*** vel);
double powder(int exp);
void printSuperPower(ofstream &superPowerout);
void factorPsi();
void printCorrelation(ofstream &file,ofstream &file1);

#endif

