#include <iostream>
#include <math.h>
#include <fstream>
//#include <cmath>
#include <cstdlib>
#include <ctime>

using namespace std;

#include "constants.h"
#include "functions.h"


double PMSDiff(int a[],int select){
	double sum=0;
	a[select]++;
	sum+=psi[LinInd(a)];
	a[select]-=2;
	sum+=psi[LinInd(a)];
	a[select]++;
	return sum;

}

void fillMDiff(){
//double alpha=D*dt/dx/dx;
	for(int i=0;i<N;i++){
		M[0][i]=M[2][i]=0.5*alpha;
		M[1][i]=-alpha-1;
	}	
	M[0][N]=M[1][N]=0.5*alpha;
}

void fillrhsDiff(int a[],int select){
	int aTemp=a[select];
	for(int i=0;i<N;i++){
		a[select]=i;	
		rhs[i]=(alpha-1)*psi[LinInd(a)];
		rhs[i]+=-alpha/2.*PMSDiff(a,select);
	}
	a[select]=aTemp;
}

void fillMDrift(int a[],int select){
	double alpha=4*dx/dt;
	int aTemp=a[select];
	int index;
	int shift=select%N_space;
	for(int i=0;i<N;i++) M[1][i]=alpha;
	for(int i=1;i<N;i++){
		a[select]=i-1;
		index = IndNN(a[select-shift],a[(select+1-shift)]);
		M[0][i]=-(extForce[index][shift][(count+1)%2]+calcWWF(a,select,false)); 
		a[select]++; 
		index = IndNN(a[select-shift],a[(select+1-shift)]);
		M[2][i-1]=+(extForce[index][shift][(count+1)%2]+calcWWF(a,select,false));
	}
	a[select]=N-1;
	index = IndNN(a[select-shift],a[(select+1-shift)]);
	M[1][N]=-(extForce[index][shift][(count+1)%2]+calcWWF(a,select,false));
	a[select]=0;
	index = IndNN(a[select-shift],a[(select+1-shift)]);
	M[0][N]=+(extForce[index][shift][(count+1)%2]+calcWWF(a,select,false));
	a[select]=aTemp;
}

void fillrhsDrift(int a[],int select){
	double alpha=4*dx/dt;
	int shift = select%N_space;
	int index;
	int aTemp=a[select];
	for(int i=0;i<N;i++){
		a[select]=i;
		rhs[i]=alpha*psi[LinInd(a)];
		a[select]++;
		index = IndNN(a[select-shift],a[(select+1-shift)]);
		rhs[i]-=(extForce[index][shift][(count)%2]+calcWWF(a,select,false))*psi[LinInd(a)];
		a[select]-=2;
		index = IndNN(a[select-shift],a[(select+1-shift)]);
		rhs[i]+=(extForce[index][shift][(count)%2]+calcWWF(a,select,false))*psi[LinInd(a)];
	}
	a[select]=aTemp;
}


void solveLin( double rhs[],double sol[]){
	if(M[1][0]==0) {cout<<"linke obere Ecke = 0"<<endl; return;}
	double gamma[N];
	double beta=M[1][0];
	sol[0]=rhs[0]/beta;
	for(int i=1;i<N;i++){
		gamma[i]=M[2][i-1]/beta;
		beta=M[1][i]-M[0][i]*gamma[i];
		if(beta==0){cout<<"error in tridag"<<endl; return;}
		sol[i]=(rhs[i]-M[0][i]*sol[i-1])/beta;
	}
	for(int i=N-2;i>=0;i--){
		sol[i]-=gamma[i+1]*sol[i+1];
	}
}

void solveOff(int a[],int select){
	double u[N],z[N],sol[N];
	int aTemp=a[select];
	M[2][N]=-M[1][0];
	M[1][0]-=M[2][N];
	M[1][N-1]-=M[0][N]*M[1][N]/M[2][N];
//	tridag(rhs,a,select);
	solveLin(rhs,sol);
	u[0]=M[2][N];
	u[N-1]=M[0][N];
	for(int i=1;i<N-1;i++) u[i]=0;
	solveLin(u,z);
//	a[select]=0;
//	double x=psi[(shift+1)%2][LinInd(a)];
//	a[select]=N-1;
	double pro=(sol[0]+M[1][N]*sol[N-1]/M[2][N])/(1.+z[0]+M[1][N]*z[N-1]/M[2][N]);
	for(int i=0;i<N;i++){/*a[select]=i;  psi[(shift+1)%2][LinInd(a)]*/
		sol[i]-=pro*z[i];
	}
	for(int i=0;i<N;i++){
		a[select]=i;
		psi[LinInd(a)]=sol[i];
	}
	M[1][0]+=M[2][N];
	M[1][N-1]+=M[0][N]*M[1][N]/M[2][N];
	a[select]=aTemp;
}

void doStep(double t){
	for(int k=0;k<2;k++){
		int a[N_PS+1];
		if(k==0) fillMDiff();	
		if(k==1&&VextTimeDependent)calcF(t);
		//	cout<< "MatrixGefüllt"<<endl;
		for(int select=0;select<N_PS;select++){
			for(int i=0;i<N_PS;i++) a[i]=0;
//			if(k==1&&VextTimeDependent)calcF(t+dt);  /////////Stimmt nicht mehr!!!!!!!!!!!!
			for(int j=0;j<N_psi/N;j++){
				if(k==0){
					fillrhsDiff(a,select);
				}else{
					fillMDrift(a,select);
					fillrhsDrift(a,select);
				}
				solveOff(a,select);
				//			cout<<"rechte seite gefüllt"<<endl;
				//			cout<<select<<"   "<< j<<endl;

				if(select==0){a[1]++;}
				else {a[0]++;}
				for(int l=0;l<N_PS-1;l++){

					if(a[l]>=N){

						a[l]=0;
						if(l+1==select)l++;
						a[l+1]++;
						if(l+1>N_PS) cout<<"error!"<<endl;
					}
				}

			}
		}
	}
}
