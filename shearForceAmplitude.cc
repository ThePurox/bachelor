#include <iostream>
#include <math.h>
#include <fstream>
//#include <cmath>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <string>

using namespace std;

string directory = "/home/nico/data/";
string filename1 = directory + "superCurr580.dat";
string filename2 = directory + "superCurr581.dat";
const int N=26;
int numbers[N]={580,581,599,688,600,689,690,691,601,692,693,694,695,653,696,654,655,710,711,712,713,714,715,716,717,718};
int F[N];
string filenames[N];
ofstream out("sheardataC.dat");
const int dim = 2;
const int coll1 = 4;
const int coll2 = 7;
const int collumns1=8;
const int lines=50;
int blocks[N];
const int skiplines=2;
const double dnx=1e-1;

int main(){
	for(int i=0;i<N;i++)
		blocks[i]=200;
	blocks[3]=100;
	blocks[5]=100;
	blocks[6]=100;
	F[0]=1;
	for(int i=1;i<=10;i++)
		F[i]=(i)*2;
	F[11]=24;
	F[12]=28;
	F[13]=32;
	F[14]=48;
	F[15]=64;
	F[16]=128;
	F[17]=0.4;
	F[18]=0.6;
	F[19]=0.8;
	F[20]=1.5;
	F[21]=3;
	F[22]=5;
	F[23]=7;
	F[24]=40;
	F[25]=90;
	for(int i=0;i<N;i++){
		filenames[i]=directory + "superCurr" + to_string(numbers[i]) + ".dat";
		ifstream file1(filenames[i]);
		cout<<file1.is_open()<<endl;
		double min1=5e5,min2=5e5,max1=-5e5,max2=-5e5;
		double a1[collumns1];
		string str;
		for(int s=0;s<skiplines;s++){
			getline(file1,str);
		}

		for(int b=0;b<blocks[i];b++){
			double dif=0,abs=0;
			for(int d=0;d<pow(lines,dim-1);d++){
				for(int l=0;l<lines;l++){
					for(int c=0;c<collumns1;c++){
						file1>>a1[c];
					}
					if(b==blocks[i]-1){
						if(a1[coll1]>max1)
							max1=a1[coll1];
						if(a1[coll1]<min1)
							min1=a1[coll1];
						if(a1[coll2]>max2)
							max2=a1[coll2];
						if(a1[coll2]<min2)
							min2=a1[coll1];
					}
					//cout<<"diff: "<<dif<<endl;
					//cout<<"Absolute: "<<abs<<endl;
				}
				getline(file1,str);
			}
			//cout<<"Block: "<<b<<endl;
			//cout<<"Relative: "<<dif/abs<<endl;
			//cout<<"Absolute: "<<dif*dnx<<endl;
			getline(file1,str);

		}
		file1.close();
		out<<F[i]<<"\t"<<max1<<"\t"<<max2<<endl;
	}
}
