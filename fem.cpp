#include <iostream>
#include <math.h>       /* pow */

int main()
{	// Test Constanten
	T = 272.15;
	nuu = 2/100;
	nuv = 0.7/100;

	// Constanten
	Dur = 2.8*pow(10.0,-10.0);
	Duz = 1.1*pow(10.0,-9.0);
	Dvr = 2.32*pow(10.0,-9.0);
	Dvz = 6.97*pow(10.0,-9.0);
	R_g = 8.314;	
	K_mu = 0.4103;
	K_mv = 27.2438;
	K_mfu = 0.1149;
	r_q = 0.97;
	hu = 7*pow(10.0,-7.0);
	hv = 7.5*pow(10.0,-7.0);
	p_atm = 101300;

	// Uitrekenbare constanten
	C_uamb = p_atm*nuu/(R_g*T);
	C_vamb = p_atm*nuv/(R_g*T);
	V_mu = 2.39*pow(10.0,-4.0)*exp(80200/R_g*(1/293.15-1/T));
	V_mfv = 1.61*pow(10.0,-4.0)*exp(56700/R_g*(1/293.15-1/T));

	// Mesh + randen
	// TODO: mesh inladen en randen bepalen.

	// Opstellen van de Matrices
	// TODO: matlabcode omzetten naar cpp.

	// Lineaire startoplossing
	// TODO: lineaire oplossing berekenen.

	// Niet lineaire oplossing
	// TODO: Juiste functies ingeven en Newton implementeren.

	// Plot
	// TODO: De oplossingen plotten.
	return 0;
}
