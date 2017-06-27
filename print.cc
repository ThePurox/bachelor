
#include <iostream>
#include <math.h>
#include <fstream>
//#include <cmath>
#include <cstdlib>
#include <ctime>

using namespace std;

#include "constants.h"
#include "functions.h"

void printPsi(double func[],int n,ofstream &out){
	if(n==N){
		for(int i=0;i<N;i++) out<<i*dx<<"\t"<<func[i]<<endl;	
		out<<endl<<endl;
		return;
	}
	if(n==NN){
		for(int i=0;i<N;i++){
			for(int j=0;j<N;j++){
				out<<i*dx<<"\t"<<j*dx<<"\t"<<func[IndNN(i,j)]<<endl;
			}
			out<<endl;
		}
		out<<endl;
		return;
	}
	if(n==N_psi){
		int a[N_PS];
		for(int i=0;i<N_PS;i++) a[i]=0;
		for(int i=0;i<N_psi;i++){
			for(int j=N_PS-1;j>=0;j--) out<<a[j]*dx<<"\t";
			out<<func[LinInd(a)]<<endl;			
			a[0]++;
			for(int j=0;j<N_PS-1;j++){
				if(a[j]>=N){a[j+1]++;a[j]=0;if(j==N_PS-2)out<<endl;}
			}	
		}	
		out<<endl;

		if(N_PS==1){
			out<<endl;
			for(int i=0;i<N;i++){
				for(int j=0;j<N;j++){
					out2<<i*dx<<"\t"<<j*dx<<"\t"<<func[i]*func[j]<<endl;
				}
				out2<<endl;
			}
			out2<<endl;

		}
	}
}

void printPsi(double func[],int n,ofstream &out,double t){
	if(n==N){
		for(int i=0;i<N;i++) out<<i*dx<<"\t"<<func[(i-int(n*4*t/tFinal)+5*n)%n]<<endl;	
		out<<endl<<endl;
		return;
	}
	if(n==NN){
		for(int i=0;i<N;i++){
			for(int j=0;j<N;j++){
				out<<i*dx<<"\t"<<j*dx<<"\t"<<func[IndNN(i,j)]<<endl;
			}
			out<<endl;
		}
		out<<endl;
		return;
	}
	if(n==N_psi){
		int a[N_PS];
		for(int i=0;i<N_PS;i++) a[i]=0;
		for(int i=0;i<N_psi;i++){
			for(int j=N_PS-1;j>=0;j--) out<<a[j]*dx<<"\t";
			out<<func[LinInd(a)]<<endl;			
			a[0]++;
			for(int j=0;j<N_PS-1;j++){
				if(a[j]>=N){a[j+1]++;a[j]=0;if(j==N_PS-2)out<<endl;}
			}	
		}	
		out<<endl;

		if(N_PS==1){
			out<<endl;
			for(int i=0;i<N;i++){
				for(int j=0;j<N;j++){
					out2<<i*dx<<"\t"<<j*dx<<"\t"<<func[i]*func[j]<<endl;
				}
				out2<<endl;
			}
			out2<<endl;

		}
	}
}
/*
   void printDensity2(int dim1,int dim2){
   if(N_PS<=1) return;
   if(dim2-dim1<0)cout<<"Fehler in printDensity2"<<endl;
   int a[N_PS];
   for(int l=0; l<N_PS; l++) a[l]=0;
   for(int i=0; i<N; i++){
   a[dim1]=i;
   for(int k=0; k<N; k++){
   a[N_PS-2]=k;
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
   density2file<<i*dx<<"\t"<<k*dx<<"\t"<<sum<<endl;
   }
   density2file << endl;
   }
   density2file << endl;
   }
 
void printRelative(){
	if(N_part!=2)return;
	int a[N_PS];
	double relaA[N][N],baryA[2*N][2*N];
	for(int i=0;i<N;i++)for(int j=0;j<N;j++)relaA[i][j]=0;
	for(int i=0;i<2*N;i++)for(int j=0;j<2*N;j++)baryA[i][j]=0;
	int dist[N_space],bary[N_space];
	for(int i=0;i<N_PS;i++) a[i]=0;
	for(int i=0;i<N_psi;i++){
		for(int j=0;j<N_space;j++){
			dist[j]=a[j]-a[j+N_space];
			dist[j]-=N*(2*dist[j]/N)-N/2;
			bary[j]=a[j]+a[j+N_space];
		}
		relaA[dist[0]][dist[1]]+=psi[LinInd(a)];
		baryA[bary[0]][bary[1]]+=psi[LinInd(a)];
		combi(a);
	}
	for(int i=0;i<N;i++){
		for(int j=0;j<N;j++){
			relativeout<<i*dx<<"\t"<<j*dx<<"\t"<<relaA[i][j]*dx*dx<<endl;
		}
		relativeout<<endl;
	}
	for(int i=0;i<2*N;i++){
		for(int j=0;j<2*N;j++){
			barycenterout<<i*dx/2<<"\t"<<j*dx/2<<"\t"<<baryA[i][j]*dx*dx<<endl;
		}
		barycenterout<<endl;
	}
	barycenterout<<endl;
	relativeout<<endl;
}
*/
void printWWF(){
	ofstream wwfout1("WWF1.dat");
	ofstream wwfout2("WWF2.dat");
	ofstream wwfout3("WWF3.dat");
	int a[N_PS];
	for(int i=0;i<N_PS;i++) a[i]=0;
	for(int i=0;i<N_psi;i++){
		for(int j=0;j<N_PS;j++) wwfout1<<a[j]<<"\t";
		wwfout1<<calcWWF(a,0,false)<<endl;
		combi(a);
	}
	for(int i=0;i<N_PS;i++) a[i]=0;
	for(int i=0;i<N_psi;i++){
		for(int j=0;j<N_PS;j++) wwfout2<<a[j]<<"\t";
		wwfout2<<calcWWF(a,1,false)<<endl;
		combi(a);
	}
	for(int i=0;i<N_PS;i++) a[i]=0;
	for(int i=0;i<N;i++){
		a[0]=i;
		wwfout3<<i*dx<<"\t"<<calcWWF(a,1,false)<<endl;
	}
}

void printCurrent(double*** current,ofstream &currentout){
	if(N_space==1){
		for(int i=0;i<N;i++){
			currentout<<i*dx;
			for(int k=0;k<3;k++) currentout<<"\t"<<current[IndNN(i,0)][k][0];
			currentout<<endl;
		}
		currentout<<endl;
	}else{
		for(int i=0;i<N;i++){
			for(int j=0;j<N;j++){
				currentout<<i*dx<<"\t"<<j*dx;
				for(int k=0;k<3;k++) currentout<<"\t"<<current[IndNN(i,j)][k][0]<<"\t"<<current[IndNN(i,j)][k][1];
				currentout<<endl;
			}
			currentout<<endl;
		}
	}
	currentout<<endl;
}

void printSuperCurr(ofstream &superCurrout){
	for(int x=0;x<N;x++){
		int Ny=N;
		if(N_space==1) Ny=1;
		for(int y=0;y<Ny;y++){
			superCurrout << x*dx << "\t" << y*dx << "\t";
			for(int space=0;space<N_space;space++)
				superCurrout<<current[IndNN(x,y)][0][space] - currentAd[IndNN(x,y)][0][space]<<"\t";
			superCurrout<<endl;
		}
		if(N_space==2)superCurrout<<endl;
	}
	if(N_space==1)superCurrout<<endl;
	superCurrout<<endl;
}

void printTime(double**arr,ofstream &file){
	int col=0;
	if(arr==::power||arr==::powerAd)
		col=6;
	if(arr==::freeEnergy)
		col=5;
	if(arr==::dissipation)
		col=3;
	if(arr==::extPower)
		col=3;	
	if(arr==::daniA)
		col=2;	
	if(arr==::squares)
		col=10;	
	for(int i=0;i<printNumb;i++){
		for(int k=0; k<col;k++)
			file<<arr[k][i]<<"\t";
		file<<endl;
	}
}
