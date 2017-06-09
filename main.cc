#include <iostream>
#include <math.h>
#include <fstream>
//#include <cmath>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <string>

using namespace std;

#include "constants.h"

const int N_part=2;
const int N_space=1;
const int N=50;
const double tFinal=15;//*sigWWPot^2*gamma/eps
const double dt=5e-3;
const double T=0.40;
const double L=10;//*sigWWPot
const int printNumb=1000;


const bool VextTimeDependent=false;
const bool superAdBool=true;
const bool analyticBool=false;

const bool printPsiBool=true;
const bool printDensityBool=true;
const bool printCurrentBool=true;
const bool printRelativeBool=true;
const bool printWWFBool=false;
const bool printPowerBool=true;

vector<string> filenames={
	"psi",
	"psiAd",
	"density",
	"densityAd",
	"2dout",
	"analytic",
	"barycenter",
	"barycenterAd",
	"relative",
	"relativeAd",
	"current",
	"currentAd",
	"adiabPotential",
	"superAd",
	"dcurrent",
	"ddensity",
	"timeDependent",
	"power",
	"powerAd",
	"psi2aus1",
	"dissipation"
};

vector<ofstream> ofstreams(filenames.size());

ofstream& psiout=(ofstreams[0]);
ofstream& psiAdout=ofstreams[1];

ofstream& densityout=ofstreams[2];
ofstream& densityAdout=ofstreams[3];

ofstream& density2file=ofstreams[4];
ofstream& anout=ofstreams[5];

ofstream& barycenterout=ofstreams[6];
ofstream& barycenterAdout=ofstreams[7];

ofstream& relativeout=ofstreams[8];
ofstream& relativeAdout=ofstreams[9];

ofstream& currentout=ofstreams[10];
ofstream& currentAdout=ofstreams[11];

ofstream& vextout=ofstreams[12];
ofstream& superAdout=ofstreams[13];

// test environment
ofstream& dcurrent=ofstreams[14];
ofstream& ddensity=ofstreams[15];
ofstream& timeout=ofstreams[16];
ofstream& powerout=ofstreams[17];
ofstream& powerAdout=ofstreams[18];
ofstream& out2=ofstreams[19];
ofstream& dissipationout=ofstreams[20];


const double sigWWPot=1;
const double eps=1;
const double extForceFactor=1;
const double extPotFactor0=1;
const int overscan=0;
long count=0;
const int N_PS=N_space*N_part;
long N_psi;
const int steps=int(tFinal/dt);
double dx=L/N;
double alpha=T*dt/dx/dx;
double sigVext0=L/6.*L/6.;
double initialNorm=1;
double norm;
int NN=1;

double* psi;
double* initial;
double* anal;
double*** extForce;
double** M;// last collum is for alpha and beta and gamma in Numerical Recepies in C page 75
double* rhs;
double* Vext;
//double* WWPot;
double*** current;
double** density;
double* psiAd;
double*** currentAd;
double** freeEnergy;
double** power;
double** powerAd;
double** dissipation;


#include "functions.h"

int LinInd(int a[N_PS]){
	int r=0;
	int n=1;
	for(int i=0;i<N_PS;i++){ r+=n*((a[i]+N*(overscan+1))%N); n*=N;}
	return r;
}

int IndNN(int x, int y){
	if(N_space<2) y=0;
	return (x+N*(overscan+1))%N + N*((y+N*(overscan+1))%N);
}

//substitutes most combinatoric sections
void combi(int a[]){
	a[0]++;
	for(int j=0; j<N_PS-1;j++){
		if(a[j]>=N){
			a[j]=0;
			a[j+1]++;
		}
	}
}

void calcF(double t){
	for(int i=0;i<N;i++){
		double x=(i)*dx;//-L/2;
		if(N_space==2){
			for(int j=0;j<N;j++){
				double y=(j)*dx;
				extForce[IndNN(i,j)][0][count%2]=0+extForceFactor*(x)*(x-L)*(x-L/2.)*64./L/L/L/L*(L*L*y*y-2*y*y*y*L+y*y*y*y)*64./L/L/L/L;//*x/L*exp(-(x*x+y*y*(1%N_space))*20/L/L);
				extForce[IndNN(i,j)][1][count%2]=0+extForceFactor*(y)*(y-L)*(y-L/2.)*64./L/L/L/L*(L*L*x*x-2*x*x*x*L+x*x*x*x)*64./L/L/L/L;//*x/L*exp(-(x*x+y*y*(1%N_space))*20/L/L);
			}
		}else{
			//extForce[IndNN(i,0)][0][count%2]=1*M_PI/L*sin(2*M_PI*(x/L+1./1.*t/tFinal));
			extForce[IndNN(i,0)][0][count%2]=0/*-2*x*exp(-x*x);/*-2*(0.5-t/tFinal)*/+extForceFactor*(x)*(x-L)*(x-L/2.)*64./L/L/L/L;//*x/L*exp(-(x*x+y*y*(1%N_space))*20/L/L);
		}
	}
	count++;
	//return T/dx/100;
//	return -2*dx*(a-N/2);
}

double calcWWF(int a[N_PS],int select,bool pot){
//	return 0;	
	double sum=0;
	double distx=0,disty=0,dist,sqdist,temp;
	int shift = select%N_space;
	if(N_space==2){
		for(int i=0;i<(N_PS-1)/2;i++){
			if(2*i==select||2*i+1==select) continue;
			distx=(a[select-shift]-a[2*i])*dx;
			distx-= L*(int)(2*distx/L);
//			distx/=L;
			disty=(a[select-shift+1]-a[2*i+1])*dx;
			disty-= L*(int)(2*disty/L);
//			disty/=L;
			sqdist=distx*distx+disty*disty;  //squared distance
//			dist=sqrt(sqdist);
//			sum+=8*dist/(1+8*8*dist*dist/L/L)/(1+8*8*dist*dist/L/L)/L*T/dx;
			temp=exp(-sqdist);
			if(!pot){
				if(select%2==0) sum+=2*distx*temp;
				else sum+=2*disty*temp;
			}else{
				sum+=temp;
			}
		}	
	}
	else{
		for(int i=0;i<N_PS;i++){
			if(i==select) continue;
			dist=(a[select]-a[i])*dx;
			dist-= L*(int)(2*dist/L);
			temp=exp(-dist*dist);
			if(!pot) sum+=2*dist*temp;
			else sum += temp;
		}
	}
	return sum;
}

void initPrint(){
	ifstream versionfile("/home/nico/version.log");
	int version;
	versionfile>>version;
	version++;
	cout<<"Version: "<<version<<"Bool: "<<versionfile.is_open()<<endl;
	versionfile.close();
	ofstream versionfile1("/home/nico/version.log");
	versionfile1<<version<<endl;
	versionfile1.close();
	string vString;
	string sep=" | ";
	string directory="/home/nico/data/";
	if(version<10) vString="00"+to_string(version);
	else if(version<100) vString="0"+to_string(version);
	else if(version>=100) vString=to_string(version);

//	sprintf(vString,"%03d",version);
	for(int i=0;i<filenames.size();i++){
		filenames[i]=directory + filenames[i] + vString + ".dat";
		ofstreams[i].open(filenames[i]);
		ofstreams[i]<<"#N_part: "<<N_part<<sep<<"N_space: "<<N_space<<sep<<"N: "<<N<<sep<<"t_final: "<<tFinal<<sep<<"dt: "<<dt<<sep<<"Boxsize: "<<L<<sep<<"printNumb: "<<printNumb<<endl;
	}



}

void initArrays(){
	N_psi=1;
	for(int i=0;i<N_PS;i++) N_psi*=N;
	psi = (double*) malloc(N_psi*sizeof(double));
	anal = (double*) malloc(N_psi*sizeof(double));
	psiAd = (double*) malloc(N_psi*sizeof(double));

//	WWPot = (double*) malloc(N*sizeof(double));
	rhs = (double*) malloc(N*sizeof(double));
//	WWPot = (double*) malloc(N*sizeof(double));

	M = (double**) malloc(3*sizeof(double*));
	for(int i=0;i<3;i++)
		M[i] = (double*) malloc((N+1)*sizeof(double));

	Vext = (double*) malloc(NN*sizeof(double));

	freeEnergy = (double**) malloc(4*sizeof(double*));
	for(int i=0;i<4;i++)
		freeEnergy[i]=(double*) malloc(printNumb*sizeof(double));

	// arrays with N*N: current, extForce

	initial = (double*) malloc(NN*sizeof(double));

	density = (double**) malloc(3*sizeof(double*));
	for(int i=0;i<3;i++) 
		density[i]= (double*) malloc(NN*sizeof(double));
	
	
	power = (double**) malloc(5*sizeof(double*));
	for(int i=0;i<5;i++) 
		power[i]= (double*) malloc(printNumb*sizeof(double));

	dissipation = (double**) malloc(3*sizeof(double*));
	for(int i=0;i<3;i++) 
		dissipation[i]= (double*) malloc(printNumb*sizeof(double));

	powerAd = (double**) malloc(5*sizeof(double*));
	for(int i=0;i<5;i++) 
		powerAd[i]= (double*) malloc(printNumb*sizeof(double));
	
	current = (double***) malloc(NN*sizeof(double**));
	for(int j=0;j<NN;j++){
		current[j] = (double**) malloc(3*sizeof(double*));
		for(int k=0;k<3;k++)
			current[j][k] = (double*) malloc(N_space*sizeof(double));
	}

	currentAd = (double***) malloc(NN*sizeof(double**));
	for(int j=0;j<NN;j++){
		currentAd[j] = (double**) malloc(3*sizeof(double*));
		for(int k=0;k<3;k++)
			currentAd[j][k] = (double*) malloc(N_space*sizeof(double));
	}

	extForce = (double***) malloc(NN*sizeof(double**));
	for(int j=0;j<NN;j++){
		extForce[j] = (double**) malloc(N_space*sizeof(double*));
		for(int k=0;k<N_space;k++)
			extForce[j][k] = (double*) malloc(2*sizeof(double));
	}
}

double getInitial(int a[N_PS]){
	double sum=0;
	int index;
	int shift;
	for(int select=0;select<N_PS;select++){
		shift = select%N_space;
		index = IndNN(a[select-shift],a[(select+1-shift)]);
		sum+=(Vext[index]+0.5*calcWWF(a,select,true));
	}
	return exp(-1./T*sum);
}

void init(){
	NN=1;
	for(int i=0; i<N_space; i++) NN*=N;
	initArrays();
	initPrint();
//	double sigVext=L*L/6./6.;
/*	double pos=0;
	for(int j=0;j<N_PS;j++){		////Anfangszustandgenerierendes externes Potential initialisieren
		for(int i=0;i<N;i++){
			pos=(i-N/4)*dx;
			Vext[j][i]=-extPotFactor0*(exp(-pos*pos/sigVext0)+exp(-(pos+L)*(pos+L)/sigVext0)+exp(-(pos-L)*(pos-L)/sigVext0));
		}
	}*/
	for(int i=0;i<NN;i++) Vext[i]=0;
		
	for(int k=-1;k<=1;k++){
		for(int i=0;i<N;i++){
			double x=(i-N/4)*dx+L*k;
			if(N_space==2){
				for(int j=0;j<N;j++){
					double y=(j-N/4)*dx+L*k;
					Vext[IndNN(i,j)]-=extPotFactor0*exp(-(x*x+y*y)/sigVext0);
				}
			}else{
				Vext[IndNN(i,0)]-=extPotFactor0*exp(-(x)*(x)/sigVext0);
			}
		}
	}

//	double sigWWPot=L*L/6./6.;		////Anfangszustandgenerierendes WW-Potential initialisieren, muss mit calcWWF zusammenpassen!!!
/*	for(int i=0;i<N;i++){ 
		double x=i*dx;
		if(i<=N/2) WWPot[i]=0* (1./4.*x*x*x*x-1./3.*x*x*x*(L+L/2.)+L*L/2./2.*x*x)/L/L/Lexp(-i*i*dx*dx);
		else       WWPot[i]=0*exp(-(i-N)*(i-N)*dx*dx);
	}*/
	int a[N_PS];
	for(int i=0;i<N_PS;i++) a[i]=0;
	for(int i=0;i<N_psi;i++){
//		double pro=1;
//		for(int j=0;j<N_PS;j++){pro*=exp(-(a[j]-N/2)*(a[j]-N/2)*dx*dx/sigVext0);}	
//		psi[i]=pro;
		psi[i]=getInitial(a);
		combi(a);		
//		psi[i]=0.5*(1+sin(2*M_PI/L*x));	
	}
	double initialNorm=1./intPsi(psi);
	for(int i=0;i<N_psi;i++) psi[i]*=initialNorm;
	density2(psi,initial,0,1);
//	if(N_PS==1) out2.open("psi2aus1.dat");	
	if(N_space==2){ 
		relativeout.open(relativefile); 
	//	relativeout<<N_psi<<endl<<((printNumb>=steps)?steps:printNumb)<<endl;;
		barycenterout.open(barycenterfile);
	//	barycenterout<<N_psi<<endl;
	}
	currentout<<"#position WWF Vext Entropy"<<endl;
	currentAdout<<"#position WWF Vext Entropy"<<endl;
}



double intPsi(double func[]){
	double sum=0;
	for(int i=0;i<N_psi;i++) sum+=func[i];
	for(int i=0;i<N_PS;i++) sum*=dx;
	return sum;
}


/*
void analytic(double t){
	int a[N_PS];
	int powder = round(pow(N*(2*overscan+1),N_PS));
	for(int i=0; i<N_PS; i++)a[i]=0;
	for(int l=0; l<N_psi; l++){
	
		int b[N_PS];
		for(int j=0; j<N_PS; j++)b[j]=-N*overscan;
		double sum=0;
		for(int i=0; i<powder; i++){
			double pro=1;
			for(int k=0; k<N_PS; k++) pro*=dx*exp(-1.*(b[k]-a[k])*(b[k]-a[k])*dx*dx/4./t/T);
			sum+=initial[LinInd(b)]*pro;
			b[0]++;
			for(int j=0; j<N_PS-1;j++){
				if(b[j]>=N*(overscan+1)){
					b[j]=0;
					b[j+1]++;
				}
			}
		}
//		for(int k=N_PS-1; k>=0; k--) anout << a[k]*dx<< "\t";
//		anout << sum*T/pow(4*M_PI*t*T,N_PS/2.) << endl;
		anal[LinInd(a)]=sum/pow(4*M_PI*t*T,N_PS/2.);
		a[0]++;
		for(int j=0; j<N_PS-1;j++){
			if(a[j]>=N){
				a[j]=0;
				a[j+1]++;
//				if(j==N_PS-2) anout << endl;
			}
		}
	}
//	anout << endl;
}
*/



void analyticfactor(double t){
	double norm=dx/pow(4*M_PI*t*T,0.5);	
	int a[N_PS];
	for(int j=0; j<N_PS;j++) a[j]=0;
	double* analtemp = (double*) malloc(NN*sizeof(double));

	if(N_space>1){
		for(int i=0; i<N; i++){
			for(int j=0; j<N; j++){
				double sum=0;
				for(int b=-N*overscan; b<N*(overscan+1); b++){
					for(int c=-N*overscan; c<N*(overscan+1); c++)
						sum+=initial[IndNN(b,c)]*exp(-((b-i)*(b-i)+(c-j)*(c-j))*dx*dx/4./t/T);
				}
				analtemp[IndNN(i,j)]=sum*norm;
			}
		}
	}else{
		for(int i=0; i<N; i++){
			double sum=0;
			for(int b=-N*overscan; b<N*(overscan+1); b++){
				sum+=initial[IndNN(b,0)]*exp(-(b-i)*(b-i)*dx*dx/4./t/T);
			}
			analtemp[IndNN(i,0)]=sum*norm;
		}	

	}
	for(int l=0; l<N_psi; l++){
		double pro = 1;
		for(int j=0; j<N_part;j++){
			pro*=analtemp[IndNN(a[2*j],a[2*j+1])];
		}
		anal[LinInd(a)] = pro;
		combi(a);
	}
	free(analtemp);
}


/*
void densitiy(double roh[N]){
	int a[N_PS];
	for(int i=0;i<N;i++) a[i]=0;
	for(int p=0;p<N_part;p++){
		for(int i=0;i<N;i++){
			a[N_space*p]=i;
			double sum=0;

	
	
	
	
}*/

void densitySum(double* psi,double rho[]){
	double temp[NN];
	for(int i=0;i<NN;i++)
		rho[i]=0;
	for(int i=0;i<N_part;i++){
		if(N_space==1)density1(psi,temp,i);
		else density2(psi,temp,2*i,2*i+1);
		for(int j=0;j<NN;j++)
			rho[j]+=temp[j];
	}
}

void density1(double* psi,double rho[],int select){
	int a[N_PS];
	for(int j=0; j<N_PS;j++) a[j]=0;
	for(int i=0; i<N; i++){
		 a[select]=i;
		 double sum=0;
		 for(int l=0; l<N_psi/N; l++){
			 sum+=psi[LinInd(a)];
			 if(select==0) a[1]++;
			 else a[0]++;
			 for(int j=0; j<N_PS-1; j++){
				 if(a[j]>=N){
					 a[j]=0;
					 if(j+1==select) a[j+2]++;
					 else a[j+1]++;
				 }
			 }
		 }
		 for(int l=0; l<N_PS-1; l++) sum*=dx;
		rho[i]=sum;
	}
}

void density2(double* psi,double rho[],int dim1,int dim2){
	int a[N_PS];
	for(int l=0; l<N_PS; l++) a[l]=0;
	if(N_space==2){
		if(dim2-dim1<0)cout<<"Fehler in printDensity2"<<endl;
		for(int i=0; i<N; i++){
			a[dim1]=i;
			for(int k=0; k<N; k++){
				a[dim2]=k;
				double sum=0;
				for(int l=0; l<N_psi/N/N; l++){
					sum+=psi[LinInd(a)];
					if(dim2-dim1==1){
						if(dim1==0) a[2]++;
						else a[0]++;
						for(int j=0;j<N_PS-1;j++){
							if(a[j]>=N){
								a[j]=0;
								if(j+1==dim1) a[j+3]++;
								else a[j+1]++;
							}
						}
					}else{
						if(dim1==0) a[1]++;
						else a[0]++;
						for(int j=0;j<N_PS-1;j++){
							if(a[j]>=N){
								a[j]=0;
								if(j+1 != dim1 && j+1 != dim2) a[j+1]++;
								else a[j+2]++;
							}
						}
					}
				}
				for(int l=0; l<N_PS-2; l++) sum*=dx;
				rho[IndNN(i,k)]=sum;
			}
		}
	}else{
		for(int i=0; i<N; i++){
			a[dim1]=i;
			double sum=0;
			for(int l=0; l<N_psi/N; l++){
				sum+=psi[LinInd(a)];
				if(dim1==0) a[1]++;
				else a[0]++;
				for(int j=0; j<N_PS-1; j++){
					if(a[j]>=N){
						a[j]=0;
						if(j+1==dim1) a[j+2]++;
						else a[j+1]++;
					}
				}
			}
			for(int l=0; l<N_PS-1; l++) sum*=dx;
			rho[IndNN(i,0)]=sum;
		}
	}
}

double residuum(double a[],double b[],int n,int dim){
	double sum=0;
	for(int i=0;i<n;i++) sum+=fabs(b[i]-a[i]);
	return sum*pow(dx,dim);
}

double residuumCurrent(int step){
	double sum=0;
	double tempDens[NN];
	densitySum(psi,density[(step+1)%3]);
	if(N_space==1){
		if(step>1){
			for(int i=0;i<N;i++){
				double test;
				//test environment
				dcurrent << i << "\t" << (sumCurrent(i+1,0,0)-sumCurrent(i-1,0,0))/(2*dx) << endl;
				ddensity << i << "\t" << -(density[(step+1)%3][i]-density[(step-1+3)%3][i])/(2*dt) << endl;
				test= fabs((sumCurrent(i+1,0,0)-sumCurrent(i-1,0,0))/(2*dx) + (density[(step+1)%3][i]-density[(step-1+3)%3][i])/(2*dt));
				//cout<<i<<"\t"<<test<<endl;
				sum+=test;
			}
			dcurrent << endl << endl;
			ddensity << endl << endl;

			sum*=dx;
		}
	}
	else{
		if(step>1){
			for(int posx=0;posx<N;posx++){
				for(int posy=0;posy<N;posy++){
					double test;
					//test environment
					dcurrent << posx << "\t" << posy << "\t"  << (sumCurrent(posx+1,posy,0)-sumCurrent(posx-1,posy,0)+sumCurrent(posx,posy+1,1)-sumCurrent(posx,posy-1,1))/(2*dx) << endl;
					ddensity << posx << "\t" << posy << "\t"  << -(density[(step+1)%3][IndNN(posx,posy)]-density[(step-1+3)%3][IndNN(posx,posy)])/(2*dt) << endl;
					test= fabs((sumCurrent(posx+1,posy,0)-sumCurrent(posx-1,posy,0)+sumCurrent(posx,posy+1,1)-sumCurrent(posx,posy-1,1))/(2*dx) + (density[(step+1)%3][IndNN(posx,posy)]-density[(step-1+3)%3][IndNN(posx,posy)])/(2*dt));
					//cout<<i<<"\t"<<test<<endl;
					sum+=test;
				}
				dcurrent << endl;
				ddensity << endl;
			}
			dcurrent << endl;
			ddensity << endl;
			sum*=dx*dx;
		}
	}
	calcCurrent(psi,current);
	return sum;
}


double sumCurrent(int posx, int posy, int dir){
	double sum=0;
	if(N_space==1) dir=0;
	for(int i=0; i<3; i++){
		sum+=current[IndNN(posx,posy)][i][dir];
	}
	return sum;
}

double diffPsi(double* psi,int a[],int select){
	double sum=0;
	a[select]++;
	sum+=psi[LinInd(a)];
	a[select]-=2;
	sum-=psi[LinInd(a)];
	a[select]++;
	return sum/(2*dx);

}


void calcCurrent(double* psi,double*** current){
	double powder = pow(dx,N_part-1);
	if(N_space==2){
		for(int i=0;i<N;i++){
			for(int j=0;j<N;j++){
				int a[N_PS];
				double sum[3][2]={0,0,0,0,0,0};
				for(int select=0;select<N_part;select++){
					for(int k=0;k<N_PS;k++)a[k]=0;
					a[2*select]=i;
					a[2*select+1]=j;
					for(int l=0; l<N_psi/N/N; l++){
						for(int space=0;space<N_space;space++){
							sum[0][space]+=psi[LinInd(a)]*calcWWF(a,2*select+space,false);
							if(current==::current) sum[1][space]+=psi[LinInd(a)]*extForce[IndNN(a[2*select],a[2*select+1])][space][(count+1)%2];
							if(current==::currentAd) sum[1][space]+=psi[LinInd(a)]*diffVext(a,2*select,space);
							sum[2][space]-=T*diffPsi(psi,a,2*select+space);
						}
						if(2*select==0) a[2]++;
						else a[0]++;
						for(int k=0;k<N_PS-1;k++){
							if(a[k]>=N){
								a[k]=0;
								if(k+1==2*select) a[k+3]++;
								else a[k+1]++;
							}
						}
					}
				}
				for(int k=0;k<3;k++){
					current[IndNN(i,j)][k][0]=sum[k][0]*powder*powder;
					current[IndNN(i,j)][k][1]=sum[k][1]*powder*powder;
				}
			}
		}
	}else{
			for(int j=0;j<N;j++){
				int a[N_PS];
				double sum[3]={0,0,0};
				for(int select=0;select<N_part;select++){
					for(int k=0;k<N_PS;k++)a[k]=0;
					a[select]=j;
					for(int l=0; l<N_psi/N; l++){
							sum[0]+=psi[LinInd(a)]*calcWWF(a,select,false);
							if(current==::current) sum[1]+=psi[LinInd(a)]*extForce[IndNN(a[select],0)][0][(count+1)%2];
							if(current==::currentAd) sum[1]+=psi[LinInd(a)]*diffVext(a,select,0);
							sum[2]-=T*diffPsi(psi,a,select);
						if(select==0) a[1]++;
						else a[0]++;
						for(int k=0;k<N_PS-1;k++){
							if(a[k]>=N){
								a[k]=0;
								if(k+1==select) a[k+2]++;
								else a[k+1]++;
							}
						}
					}
				}
				for(int k=0;k<3;k++){
					current[IndNN(j,0)][k][0]=sum[k]*powder;
				}
			}
		}

}

double diffVext(int a[],int part,int space){
	return -(Vext[IndNN(a[part]+1-space,a[part+1]+space)]-Vext[IndNN(a[part]-1+space,a[part+1]-space)])/(2*dx);
}

void calcRelative(double* psi, double relaA[],double baryA[]){
	if(N_part!=2){
		cout<<"Fehler in calcRelative"<<endl;
		return;
	}
	int a[N_PS];
	int dist[N_space],bary[N_space];
	for(int i=0;i<NN;i++){
		relaA[i]=0;
		baryA[i]=0;
	}
	for(int i=0;i<N_PS;i++) a[i]=0;
	for(int i=0;i<N_psi;i++){
		for(int j=0;j<N_space;j++){
			dist[j]=a[j]-a[j+N_space];
			dist[j]-=N*(2*dist[j]/N);//-N/2;
			bary[j]=round(a[j]-dist[j]/2);//(a[j]+a[j+N_space])/2;
	//		cout<<bary[j]<<endl;
			dist[j]+=N/2;
		}
		relaA[IndNN(dist[0],dist[1])]+=psi[LinInd(a)];
		baryA[IndNN(bary[0],bary[1])]+=psi[LinInd(a)];
		combi(a);
	}
	for(int i=0;i<NN;i++){
		for(int j=0;j<N_space;j++){
			relaA[i]*=dx;
			baryA[i]*=dx;
		}
	}
}

void calcPower(double* psi,double**vSquared,int step,int printStep){
	double powder = pow(dx,N_PS);
	for(int k=0; k<5;k++)
		vSquared[k][printStep]=0;
	int a[N_PS];
	for(int select=0;select<N_part;select++){
		for(int k=0;k<N_PS;k++)a[k]=0;
		for(int l=0; l<N_psi; l++){
			double sum[3][N_space];
			double v[N_space];
			for(int space=0;space<N_space;space++){
				v[space]=0;
				sum[0][space]=calcWWF(a,N_space*select+space,false);
				if(psi==::psi) sum[1][space]=extForce[IndNN(a[N_space*select],a[N_space*select+1])][space][(count+1)%2];
				if(psi==::psiAd) sum[1][space]=diffVext(a,N_space*select,space);
				sum[2][space]=-T*diffPsi(psi,a,N_space*select+space)/psi[LinInd(a)];
				for(int k=0; k<3;k++)
					v[space]+=sum[k][space];	
			}
			for(int k=0; k<3;k++){
				for(int space=0;space<N_space;space++)
					vSquared[k+1][printStep]+=0.5*v[space]*psi[LinInd(a)]*sum[k][space];
				
			}
			combi(a);	
		}
	}
	for(int k=0;k<3;k++)
		vSquared[4][printStep]+=vSquared[k+1][printStep];
	for(int k=0; k<4;k++)
		vSquared[k+1][printStep]*=powder;
	vSquared[0][printStep]=step*dt;

}

void calcDiss(int step,int printStep){
	if(printStep!=0 && printStep!=printNumb){
		dissipation[0][printStep]=step*dt;
		dissipation[1][printStep-1]=power[4][printStep-1]-(freeEnergy[3][printStep]-freeEnergy[3][printStep-2])/(dt*steps/printNumb*2.);
		dissipation[2][printStep]=0;
		double curr;
		for(int i=0;i<NN;i++){
			for(int space=0;space<N_space;space++){
				curr=sumCurrent(i%N,i/N,space);
				dissipation[2][printStep]+=curr*curr/density[(step+1)%3][i];
			}
		}
		dissipation[2][printStep]*=pow(dx,N_space)/2.;
	}
		
}

int main(){
	init();
	double psi0=intPsi(psi);
	double baryA[NN],relaA[NN];
	cout<<psi0<<endl;
	time_t t_0;
	time(&t_0);
	if(printWWFBool)	printWWF();

	{calcF(0);calcF(0);}
	cout<<"N Space: "<<N_space<<endl<<"N Part: "<<N_part<<endl;
	for(int i=0;i<steps;i++){
		calcCurrent(psi,current);
		//		printPsi(density[i%3],N,densityout);
		if(i%int(ceil(steps*1./printNumb))==0) {
			cout<<(100.*i)/steps<<"%"<<endl;
			if(superAdBool){
				doSuperAd(i,i*printNumb/steps);
				cout<<"Resdiuum SuperAdiabatic: \t" << residuum(psi,psiAd,N_psi,N_PS)<<endl;
				if(N_PS<=2||printPsiBool){ 
					printPsi(psiAd,N_psi,psiAdout); 
				}
				if(printDensityBool){
					densitySum(psiAd,baryA);
					printPsi(baryA,NN,densityAdout);
				}
				if(printRelativeBool){
					calcRelative(psiAd,relaA,baryA);
					printPsi(relaA,NN,relativeAdout);
					printPsi(baryA,NN,barycenterAdout);
				}
				if(printCurrentBool){
					calcCurrent(psiAd,currentAd);
					printCurrent(currentAd,currentAdout);
				}
				if(printPowerBool){
					calcPower(psiAd,powerAd,i,i*printNumb/steps);
				}
			}
			if(N_PS<=2||printPsiBool){ 
				printPsi(psi,N_psi,psiout); 
			}
			if(printDensityBool){
				densitySum(psi,baryA);
				printPsi(baryA,NN,densityout);
			}
			if(printRelativeBool){
				calcRelative(psi,relaA,baryA);
				printPsi(relaA,NN,relativeout);
				printPsi(baryA,NN,barycenterout);
			}
			if(printCurrentBool){
				calcCurrent(psi,current);
				printCurrent(current,currentout);
				residuumCurrent(i);
			}
			if(printPowerBool)
				calcPower(psi,power,i,i*printNumb/steps);
			if(analyticBool){
				analyticfactor((i+1)*dt);  ////////////////////////////////////////////hier noch +1? glaube nicht weil doStep() jetzt weiter unten ist
				double div=1./intPsi(anal);
				cout<<"Analytic: \t" << 1.-intPsi(anal)<<endl;
				for(int j=0;j<N_psi;j++) anal[j]*=div;
				printPsi(anal,N_psi,anout);
				cout<<"Resdiuum Analytic: \t" << residuum(psi,anal,N_psi,N_PS)<<endl;
			}
			double psi1=intPsi(psi); 
			cout<<"Deviation: \t"<<(psi1-psi0)/psi0<<endl;
			calcDiss(i,i*printNumb/steps);
		}else{
			if(printCurrentBool){
				residuumCurrent(i);
			}
		}
		doStep(i*dt);
	}
	if(printPowerBool){
		if(superAdBool)
			printTime(powerAd,powerAdout);
		printTime(power,powerout);
	}
	if(superAdBool)
		printTime(freeEnergy, timeout);

	printTime(dissipation,dissipationout);

	cout<<"Free Energy deviation: "<< residuumFreeEnergy()<<endl;
	time_t t_1;
	time(&t_1);
	cout<<"Minutes: " <<difftime(t_1,t_0)/60<<endl;

}




