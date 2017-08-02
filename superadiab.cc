#include <iostream>
#include <math.h>
#include <fstream>
//#include <cmath>
#include <cstdlib>
#include <ctime>

using namespace std;

#include "constants.h"
#include "functions.h"

/*
void calcSuperAdFF(){
	int a[N_PS];
	for(int i=0;i<N_psi;i++){
		superAdFF[LinInd(a)]=log(psi[LinInd(a)]);
		superAdFF[LinInd(a)]*=T;

		superAdFF[LinInd(a)]-=calcWWF(a,select,true);
		combi(a);

	
	}
}*/
void shiftVext(){
	double offset=1e6;
	for(int i=0;i<NN;i++){
		if(Vext[i]<offset)
			offset=Vext[i];
	}
	for(int i=0;i<NN;i++){
		Vext[i]-=offset;
	}
	int a[N_PS+1];
		for(int i=0;i<N_PS;i++) a[i]=0;
	for(int i=0;i<N_psi;i++){
			psiAd[LinInd(a)]=getInitial(a);
			combi(a);
		}

		norm=1./intPsi(psiAd);	
for(int i=0;i<N_psi;i++)
			psiAd[i]*=norm;

}

void doSuperAd(){
	double densAd[NN];
	double res=1;
	int counter=0;
	int a[N_PS+1];
	do{
		for(int i=0;i<N_PS;i++) a[i]=0;
		for(int i=0;i<N_psi;i++){
			psiAd[LinInd(a)]=getInitial(a);
			combi(a);
		}
		norm=1./intPsi(psiAd);	
		for(int i=0;i<N_psi;i++)
			psiAd[i]*=norm;
		densitySum(psiAd,densAd);
		res=residuum(density[(step)%3],densAd,NN,N_part-1);
		for(int i=0;i<NN;i++){
			Vext[i]+=0.351*log(densAd[i]/density[(step)%3][i]);
			if(densAd[i]<=0)cout<<"Fehler " <<densAd[i]<<endl;
		}
		counter++;
//		cout<<res<<endl;
		//			if(counter>=1000) cout<<"super ad did not converge"<<endl;
	}while(fabs(res)>=1e-6 && counter<10000);
	shiftVext();
	cout<<"After "<<counter<<" iterations Super Ad did converge"<<endl;
	cout<<"With Residuum "<<res<<endl;
	calcFreeEnergy();
	printPsi(Vext,NN,vextout);
}



void calcFreeEnergy(){
	double wavelength=3;
	double fE=T*log(pow(wavelength,N_PS)*norm*factorial(N_part));
	freeEnergy[3][printStep]=fE;
	double fEtemp=0;
	freeEnergy[5][printStep]=0;
	for(int i=0;i<NN;i++){
		fEtemp+=Vext[i]*density[(step)%3][i];
		freeEnergy[5][printStep]+=density[(step)%3][i]*(log(wavelength*density[(step)%3][i])-1);
	}
	fE-=fEtemp*powderNN;
	double dtF=0;
	for(int x=0;x<N;x++){
		for(int y=0;y<=(N-1)*(N_space-1);y++){
			for(int space=0;space<N_space;space++)
				dtF-=(Vext[IndNN(x+1-space,y+space)]-Vext[IndNN(x-1+space,y-space)])/(2)*sumCurrent(x,y,space);
		}
	}
	dtF*=powder(N_space-1);
	freeEnergy[0][printStep]=step*dt;
	freeEnergy[1][printStep]=fE;
	freeEnergy[2][printStep]=dtF;
	freeEnergy[5][printStep]*=T*powderNN;
}

double residuumFreeEnergy(){
	double sum=0,sum1=0;
	for(int i=1;i<printNumb-1;i++){
		freeEnergy[4][i]=(freeEnergy[1][i+1]-freeEnergy[1][i-1])/(freeEnergy[0][i+1]-freeEnergy[0][i-1]);
		sum+=fabs(freeEnergy[2][i]-freeEnergy[4][i])*(freeEnergy[0][i+1]-freeEnergy[0][i-1]);
		sum1+=fabs(freeEnergy[2][i])*(freeEnergy[0][i+1]-freeEnergy[0][i-1]);
	}
	
	freeEnergy[4][0]=(freeEnergy[1][1]-freeEnergy[1][0])/(freeEnergy[0][1]-freeEnergy[0][0]);
	freeEnergy[4][printNumb-1]=(freeEnergy[1][printNumb-1]-freeEnergy[1][printNumb-2])/(freeEnergy[0][printNumb-1]-freeEnergy[0][printNumb-2]);
	return sum/sum1;
}

int factorial(int n){
	int x=1;
	for(int i=2;i<=n;i++)
		x*=i;
	return x;
}


