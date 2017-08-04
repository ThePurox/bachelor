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
string filename1 = directory + "density500.dat";
string filename2 = directory + "density500.dat";
const int dim = 1;
const int coll1 = 2;
const int coll2 = 2;
const int collumns1=2;
const int collumns2=2;
const int lines=100;
const int blocks=100;
const int skiplines=2;
const double dnx=1e-1;

int main(){
	ifstream file1(filename1);
	ifstream file2(filename2);

	double a1[collumns1];
	double a2[collumns2];
	string str;
	for(int s=0;s<skiplines;s++){
		getline(file1,str);
		getline(file2,str);
	}
	
	for(int b=0;b<blocks;b++){
		double dif=0,abs=0;
		for(int d=0;d<pow(lines,dim-1);d++){
			for(int l=0;l<lines;l++){
				for(int c=0;c<collumns1;c++){
					file1>>a1[c];
				}for(int c=0;c<collumns2;c++){
					file2>>a2[c];
				}
				dif+=fabs(a2[coll2-1]-a1[coll1-1]);
				abs+=fabs(a1[coll1-1]);
				//cout<<"diff: "<<dif<<endl;
				//cout<<"Absolute: "<<abs<<endl;
			}
			getline(file1,str);
			getline(file2,str);
		}
		cout<<"Block: "<<b<<endl;
		cout<<"Relative: "<<dif/abs<<endl;
		cout<<"Absolute: "<<dif*dnx<<endl;
		getline(file1,str);
		getline(file2,str);
	
		}
	cout << file1.is_open();
}
