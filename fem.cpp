#include <iostream>
#include <string>
#include <iomanip>
#include <math.h>       /* pow */
#include <fstream>	// Open file
using namespace std;

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
	
	i = 0;
	while(i<lines_p){
	   getline(fp, temp,',');		
	   p[i][0] = stof(temp);	
	   getline(fp, temp);
	   p[i][1] = stof(temp); 
	   i = i+1;
	}
	
	i = 0;
	while(i<lines_t){
	  getline(ft,temp,',');
	  t[i][0] = stoi(temp);
	  getline(ft,temp,',');
	  t[i][1] = stoi(temp);
	  getline(ft,temp);
	  t[i][2] = stoi(temp);
	  i = i+1;
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

	// LINEAIRE OPLOSSING
	// TODO: lineaire oplossing berekenen.

	// NIET LINEAIRE OPLOSSING
	// TODO: Juiste functies ingeven en Newton implementeren.

	// PLOT OPLOSSINGEN
	// TODO: De oplossingen plotten.
	return 0;
}
