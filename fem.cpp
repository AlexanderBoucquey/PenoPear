#include <iostream>
#include <math.h>       /* pow */
#include <fstream>	// Open file

int main()
{	// VARIABELEN
	// Variabelen voor mesh
	vector<float> p, t;
	int nump{0};
	int numt{0};
	ifstream fp,ft;

	// TEST CONSTANTEN
	float T = 272.15;
	float nuu = 2/100;
	float nuv = 0.7/100;

	// CONSTANTEN
	float Dur = 2.8*pow(10.0,-10.0);
	float Duz = 1.1*pow(10.0,-9.0);
	float Dvr = 2.32*pow(10.0,-9.0);
	float Dvz = 6.97*pow(10.0,-9.0);
	float R_g = 8.314;	
	float K_mu = 0.4103;
	float K_mv = 27.2438;
	float K_mfu = 0.1149;
	float r_q = 0.97;
	float hu = 7*pow(10.0,-7.0);
	float hv = 7.5*pow(10.0,-7.0);
	float p_atm = 101300;

	// UITREKENBARE CONSTANTEN
	float C_uamb = p_atm*nuu/(R_g*T);
	float C_vamb = p_atm*nuv/(R_g*T);
	float V_mu = 2.39*pow(10.0,-4.0)*exp(80200/R_g*(1/293.15-1/T));
	float V_mfv = 1.61*pow(10.0,-4.0)*exp(56700/R_g*(1/293.15-1/T));

	// MESH
	// p zijn de x en y waarden van alle mogelijke hoeken van de triangels.
	// t zijn alle geconnecteerde triangles.	
	fp.open("meshp.txt");
	ft.open("mesht.txt");
	while (fin >> nump)
	  p.emplace_back(nump);
	while (ft >> numt)
	  t.emplace_back(numt);

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
