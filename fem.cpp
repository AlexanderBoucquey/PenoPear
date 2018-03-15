#include <iostream>
#include <string>
#include <iomanip>
#include <math.h>       /* pow */
#include <fstream>	// Open file
using namespace std;

int mostFrequentElement(int arr[], int x);
int numberOccurences(int arr[], int n, int x);

int main()
{	// VARIABELEN
	// Variabelen voor mesh
	
	int number_of_linesp;
	int nump{0};
	int numt{0};
	string line,temp;
	ifstream fp,ft;
	int i;

	// TEST CONSTANTEN
	float T = 272.15;
	float nuu = 2/100;
	float nuv = 0.7/100;

	// CONSTANTEN
	float Dur = 2.8e-10;
	float Duz = 1.1e-9;
	float Dvr = 2.32e-9;
	float Dvz = 6.97e-9;
	float R_g = 8.314;	
	float K_mu = 0.4103;
	float K_mv = 27.2438;
	float K_mfu = 0.1149;
	float r_q = 0.97;
	float hu = 7e-7;
	float hv = 7.5e-7;
	float p_atm = 101300;

	// UITREKENBARE CONSTANTEN
	float C_uamb = p_atm*nuu/(R_g*T);
	float C_vamb = p_atm*nuv/(R_g*T);
	float V_mu = 2.39e-4*exp(80200/R_g*(1/293.15-1/T));
	float V_mfv = 1.61e-4*exp(56700/R_g*(1/293.15-1/T));

	// MESH
	// p zijn de x en y waarden van alle mogelijke hoeken van de triangels.
	// t zijn alle geconnecteerde triangles.	
	fp.open("meshp.txt");
	ft.open("mesht.txt");
	

	// Tel aantal lijnen code.
	int lines_p = 0;
	int lines_t = 0;

	while(!getline(fp, line).eof())
		lines_p++;
	
	while(!getline(ft, line).eof())
		lines_t++;

	fp.close();
	ft.close();

	// Initialiseer de p en t arrays.
	float **p = new float*[lines_p];
	for(int i = 0; i < lines_p; ++i) {
	    p[i] = new float[2];
	}
	int **t = new int*[lines_t];
	for(int i = 0; i < lines_t; ++i) {
	    t[i] = new int[3];
	}
	
	// Vul de p en de t arrays in.
	fp.open("meshp.txt");
	ft.open("mesht.txt");
	
	for(int i = 0;i<lines_p; ++i){
	   getline(fp, temp,',');		
	   p[i][0] = stof(temp);	
	   getline(fp, temp);
	   p[i][1] = stof(temp); 
	}
	
	for(int i = 0;i<lines_t; i++){
	  getline(ft,temp,',');
	  t[i][0] = stoi(temp);
	  getline(ft,temp,',');
	  t[i][1] = stoi(temp);
	  getline(ft,temp);
	  t[i][2] = stoi(temp);	  
	}

	cout << "Mesh is correct aangemaakt!" << endl;
	//(!getline(ft,line).eof()) && 
	// testen of mesh correct werkt
	/*
	ofstream out("test_mesh_p.txt");
	for(int i=0;i<100;i++){
		for (int j=0;j<2;j++){
			out<<p[i][j];
		}
		out<< "\n";	
	}
	out.close();

	ofstream out2("test_mesh_t.txt");
	for(int i=0;i<100;i++){
		for (int j=0;j<3;j++){
			out2<<t[i][j];
		}
		out2<< "\n";	
	}
	out2.close();*/


	// OPSTELLEN VAN DE MATRICES
	// TODO: matlabcode omzetten naar cpp.
	float r[3] = {};
	float z[3] = {};
	float Ku[lines_p][lines_p] = {};
	float Kv[lines_p][lines_p] = {};
	float C[lines_p][lines_p] = {};
	cout <<"Aantal lijnen p:"<< lines_p<<endl;
	for (int i = 0; i< lines_t; i++){
	  for (int j = 0; j <3; ++j){
	    r[j]= p[t[i][j]-1][0];
	    z[j]= p[t[i][j]-1][1];
	  }
          float a[3] = {r[2]*z[1]-r[1]*z[2], r[0]*z[2]-r[2]*z[0], r[1]*z[0]-r[0]*z[1]}; 
       	  float b[3] = {z[2]-z[1],z[0]-z[2],z[1]-z[0]};
          float c[3] = {r[1]-r[2],r[2]-r[0],r[0]-r[1]};
          float area = r[1]*z[2]+r[0]*z[1] + r[2]*z[0] - r[1]*z[0] - r[0]*z[2] - r[2]*z[1]; 
          float K[3][3] = {};
          for (int i = 0; i < 3; i++){
            for (int j = 0; j < 3; ++j){
              Ku[i][j] = (r[0] + r[1] + r[2])/(12*area)*(Dur*b[i]*b[j] + Duz*c[i]*c[j]);
            }
          }
          for (int i = 0; i < 3; i++){
            for (int j = 0; j < 3; ++j){
              Kv[i][j] = (r[0] + r[1] + r[2])/(12*area)*(Dvr*b[i]*b[j] + Dvz*c[i]*c[j]);
            }
          }
	float C[3][3]= {{area/60*6*r[1]+2*r[2]+2*r[3], area/60*2*r[1]+2*r[2]+r[3], area/60*2*r[1]+r[2]+2*r[3]},{area/60*2*r[1]+2*r[2]+r[3],area/60* 2*r[1]+6*r[2]+2*r[3], area/60*r[1]+2*r[2]+2*r[3]},{area/60*2*r[1]+r[2]+2*r[3], area/60*r[1]+2*r[2]+2*r[3], area/60*2*r[1]+2*r[2]+6*r[3]}};
	  for (int m = 0; m <3; ++m){
	    for (int j = 0; j <3; ++j){
	      Ku[t[m][j]][t[m][j]] = Ku[t[m][j]][t[m][j]] + Ku[m][j];
    	      Kv[t[m][j]][t[m][j]] = Kv[t[m][j]][t[m][j]] + Kv[m][j];
    	      C[t[m][j]][t[m][j]] = C[t[m][j]][t[m][j]] + C[m][j];
	    }
	  }	  
	}

	// BEPALEN BOUNDARY
	int x = mostFrequentElement(t[0], lines_t);
	int n =	0;
	for (int i = 0; i < 3; ++i)
	  n = n + numberOccurences(t[i], lines_t, x);
	cout << "Een meest voorkomend element: " << x<<endl;
	cout << "Aantal voorkomens: "<< n << endl;
	/*l = sqrt((rp(1:end-1,1)-rp(2:end,1)).^2+(rp(1:end-1,2)-rp(2:end,2)).^2);
l(length(rn)) = sqrt((rp(1,1)-rp(end,1)).^2+(rp(1,2)-rp(end,2)).^2);
r = rp(:,1);
K_h = zeros(length(p),length(p));
R_q = zeros(length(p),1);

for o = 1:length(rn)-1
    if((r(o)>1E-13) || (r(o+1)>1E-13)) 
        K_h([rn(o) rn(o+1)],[rn(o) rn(o+1)])= K_h([rn(o) rn(o+1)],[rn(o) rn(o+1)])+Kh(r(o:o+1),l(o));
       % disp(o);
        R_q([rn(o) rn(o+1)]) = R_q([rn(o) rn(o+1)]) + Rq(r(o:o+1),l(o));
    end
end


if((r(end)>1E-13) || (r(1)>1E-13)) 
  K_h([rn(end) rn(1)],[rn(length(rn)) rn(1)])= K_h([rn(end) rn(1)],[rn(end) rn(1)])+Kh(r([length(rn) 1]),l(length(rn)));
  R_q([rn(end) rn(1)]) = R_q([rn(end) rn(1)]) + Rq(r([length(rn) 1]),l(length(rn)));
end */


	// LINEAIRE OPLOSSING
	// TODO: lineaire oplossing berekenen.

	// NIET LINEAIRE OPLOSSING
	// TODO: Juiste functies ingeven en Newton implementeren.

	// PLOT OPLOSSINGEN
	// TODO: De oplossingen plotten.
	/*Matplotlib in python heeft kei veel coole plotters en redelijk deftige documentatie
Ik denk dat wij nu trisurf gebruiken als functie
	*/
	return 0;
}

int mostFrequentElement(int arr[], int x){
	int current_candidate = arr[0];
	int counter = 0;
	int i = 0;
	for (i = 0; i < x; ++i){
	  if(current_candidate == arr[i]){
	    ++counter;
	} else if (counter ==0){
	    current_candidate = arr[i];
	    ++counter;
	} else {
	    --counter;
	} cout << current_candidate << endl;
	}
	return current_candidate;
}

int numberOccurences(int arr[], int n, int x){
	int res = 0;
	for (int i=0; i<n; i++)
	  if (x == arr[i])
	    res++;
	return res;
}
