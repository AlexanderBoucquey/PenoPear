#include <iostream>
#include <string>
#include <iomanip>
#include <math.h>       /* pow */
#include <fstream>	// Open file
#include <Eigen/Dense>  // Eigen Library
#include <Eigen/IterativeLinearSolvers> // Conjugate gradient

using namespace std;
using namespace Eigen; // Eigen library

int mostFrequentElement(int arr[], int x);
int numberOccurences(int arr[], int n, int x);

int main()
{	// VARIABELEN
	// Variabelen voor mesh
	typedef double precision;
	int number_of_linesp;
	int nump{0};
	int numt{0};
	string line,temp;
	ifstream fp,ft, frn, frp;
	precision Kh, Rq;	

	// TEST CONSTANTEN
	precision T = 272.15;
	precision nuu = 0.02;
	precision nuv = 0.007;

	// CONSTANTEN
	precision Dur = 2.8e-10;
	precision Duz = 1.1e-9;
	precision Dvr = 2.32e-9;
	precision Dvz = 6.97e-9;
	precision R_g = 8.314;	
	precision K_mu = 0.4103;
	precision K_mv = 27.2438;
	precision K_mfu = 0.1149;
	precision r_q = 0.97;
	precision hu = 7e-7;
	precision hv = 7.5e-7;
	precision p_atm = 101300;

	// UITREKENBARE CONSTANTEN
	precision C_uamb = ((precision)((precision)p_atm*(precision)nuu)/(precision)((precision)R_g*(precision)T));
	precision C_vamb = p_atm*nuv/(R_g*T);
	precision V_mu = 2.39e-4*exp((precision)80200.0/R_g*((precision)1.0/293.15-(precision)1.0/T));
	precision V_mfv = 1.61e-4*exp((precision)56700.0/R_g*((precision)1.0/293.15-(precision)1.0/T));

	// MESH
	// p zijn de x en y waarden van alle mogelijke hoeken van de triangels.
	// t zijn alle geconnecteerde triangles.	
	// rn zijn randwaarden.
	// rp zijn de randpunten.
	fp.open("meshp.txt");
	ft.open("mesht.txt");
	frn.open("rn.txt");
	frp.open("rp.txt");

	// Tel aantal lijnen code.
	int lines_p= 0;
	int lines_t= 0;
	int lines_rn= 0;
	int lines_rp= 0;

	while(!getline(fp, line).eof())
		lines_p++;
	
	while(!getline(ft, line).eof())
		lines_t++;
	
	while(!getline(frn, line).eof())
		lines_rn++;
	
	while(!getline(frp, line).eof())
		lines_rp++;

	fp.close();
	ft.close();
	frn.close();
	frp.close();

	// Initialiseer de p,t, rn en rp arrays.
	precision **p = new precision*[lines_p];
	for(int i = 0; i < lines_p; ++i) {
	    p[i] = new precision[2];
	}

	int **t = new int*[lines_t];
	for(int i = 0; i < lines_t; ++i) {
	    t[i] = new int[3];
	}
	
	int *rn = new int[lines_rn];

	precision **rp = new precision*[lines_rp];
	for(int i = 0; i < lines_rp; ++i) {
	    rp[i] = new precision[2];
	}

	// Vul de p,t, rn en rp arrays in.
	fp.open("meshp.txt");
	ft.open("mesht.txt");
	frn.open("rn.txt");
	frp.open("rp.txt");
	
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

	for(int i = 0;i<lines_rn; ++i){
	   getline(frn, temp);
	   rn[i] = stof(temp); 
	}

	for(int i = 0;i<lines_rp; ++i){
	   getline(frp, temp,',');		
	   rp[i][0] = stof(temp);	
	   getline(frp, temp);
	   rp[i][1] = stof(temp); 
	}

	cout << "Mesh is aangemaakt!" << endl;

	cout <<"Aantal lijnen p:"<< lines_p<<endl;

	// OPSTELLEN VAN DE MATRICES
	// Initialiseren
	precision r[3] = {};
	precision z[3] = {};
	Eigen::Matrix<precision, Eigen::Dynamic, Eigen::Dynamic> Ku, Kv, C,K_h, Ku_temp, Kv_temp, l, r2, R_q,u,A1,A2,b1,b2,v;
	Ku.resize(lines_p,lines_p);
	Kv.resize(lines_p,lines_p);
	C.resize(lines_p,lines_p);
	K_h.resize(lines_p,lines_p);
	Ku_temp.resize(3,3);
	Kv_temp.resize(3,3);
	l.resize(lines_rn,1);
	r2.resize(lines_rp,1);
	R_q.resize(lines_p,1);
	u.resize(lines_p,1);
	A1.resize(lines_p,lines_p);
	b1.resize(lines_p,1);
	v.resize(lines_p,1);
	A2.resize(lines_p,lines_p);
	b2.resize(lines_p,1);
	// Invullen Ku, Kv en C.
	for (int i = 0; i< lines_t; i++){
	  for (int j = 0; j <3; ++j){
	    r[j]= p[t[i][j]-1][0];
	    z[j]= p[t[i][j]-1][1];
	  }

          precision a[3] = {r[1]*z[2]-r[2]*z[1], r[2]*z[0]-r[0]*z[2], r[0]*z[1]-r[1]*z[0]}; 
       	  precision b[3] = {z[1]-z[2],z[2]-z[0],z[0]-z[1]};
          precision c[3] = {r[2]-r[1],r[0]-r[2],r[1]-r[0]};
          precision area = r[1]*z[2]+r[0]*z[1] + r[2]*z[0] - r[1]*z[0] - r[0]*z[2] - r[2]*z[1]; 

          for (int l = 0; l < 3; l++){
            for (int j = 0; j < 3; ++j){
              Ku_temp(l,j) = (r[0] + r[1] + r[2])/(12*area)*(Dur*b[l]*b[j] + Duz*c[l]*c[j]);
            }
          }

          for (int l = 0; l < 3; l++){
            for (int j = 0; j < 3; ++j){
              Kv_temp(l,j) = (r[0] + r[1] + r[2])/(12*area)*(Dvr*b[l]*b[j] + Dvz*c[l]*c[j]);
            }
          }

	precision C_temp[3][3] = {{area/60*6*r[1]+2*r[2]+2*r[3], area/60*2*r[1]+2*r[2]+r[3], area/60*2*r[1]+r[2]+2*r[3]},{area/60*2*r[1]+2*r[2]+r[3],area/60* 2*r[1]+6*r[2]+2*r[3], area/60*r[1]+2*r[2]+2*r[3]},{area/60*2*r[1]+r[2]+2*r[3], area/60*r[1]+2*r[2]+2*r[3], area/60*2*r[1]+2*r[2]+6*r[3]}};
	  for (int m = 0; m <3; ++m){
	    for (int j = 0; j <3; ++j){
	      Ku(t[i][m]-1,t[i][j]-1) = Ku(t[i][m]-1,t[i][j]-1) + Ku_temp(m,j);
    	      Kv(t[i][m]-1,t[i][j]-1) = Kv(t[i][m]-1,t[i][j]-1) + Kv_temp(m,j);
    	      C(t[i][m]-1,t[i][j]-1) = C(t[i][m]-1,t[i][j]-1) + C_temp[m][j];

	    }
	  }

	}
	
	// Invullen l en r.	
	for (int i = 0; i<lines_rp-1; i++){
	  l(i) = sqrt(pow(rp[i][0]-rp[i+1][0],2)+pow(rp[i][1]-rp[i+1][1],2));
	  r2(i) = rp[i][0];
	}
	//r2(lines_rp) = rp[lines_rp][0];

	// Invullen K_h en R_q.
	for (int i = 0; i < lines_rn-1; i++){
	  if ((r2(i) >1e-13) || (r2(i+1)>1e-13)){
	    Kh = l(i)/12*(3*r2(i)+r2(i+1));
	    K_h(rn[i]-1,rn[i]-1) = K_h(rn[i]-1,rn[i]-1) + Kh;

	    Kh = l(i)/12*(r2(i)+r2(i+1));
	    K_h(rn[i]-1,rn[i+1]-1) = K_h(rn[i]-1,rn[i+1]-1) + Kh;

	    Kh = l(i)/12*(r2(i)+r2(i+1));
	    K_h(rn[i+1]-1,rn[i]-1) = K_h(rn[i+1]-1,rn[i]-1) + Kh;

	    Kh = l(i)/12*(r2(i)+3*r2(i+1));
	    K_h(rn[i+1]-1,rn[i+1]-1) = K_h(rn[i+1]-1,rn[i+1]-1) + Kh;

	    Rq = l(i)/6*(2*r2(i)+r2(i+1));
	    R_q(rn[i]-1) = R_q(rn[i]-1) + Rq;

	    Rq = l(i)/6*(r2(i)+2*r2(i+1));
	    R_q(rn[i+1]-1) = R_q(rn[i+1]-1) + Rq;

	  }
	}

	
	// LINEAIRE OPLOSSING
	// TODO: lineaire oplossing berekenen.
//
	A1 = Ku + (V_mu/K_mu)*C + hu*K_h;
	b1 = hu*C_uamb*R_q;
	ConjugateGradient<Matrix<precision, Eigen::Dynamic, Eigen::Dynamic>, Lower|Upper> cg;
	cg.compute(A1);
	u = cg.solve(b1);
	std::cout << "#iterations:     " << cg.iterations() << std::endl;
	std::cout << "estimated error: " << cg.error()      << std::endl;
	A2 = Kv + hv*K_h;
	b2 = r_q*(V_mu/K_mu)*C*u + hv*C_vamb*R_q;
	cg.compute(A2);
	v = cg.solve(b2);
	std::cout << "#iterations:     " << cg.iterations() << std::endl;
	std::cout << "estimated error: " << cg.error()      << std::endl;
	cout<<v<<endl;
	cout<<u<<endl;
//	SparseMatrix<float> A(lines_p,lines_p);
//	u = Map<VectorXf>(R_q);//*
//	u = hu*C_uamb*R_q;	

	// NIET LINEAIRE OPLOSSING
	// TODO: Juiste functies ingeven en Newton implementeren.

	// PLOT OPLOSSINGEN
	// TODO: De oplossingen plotten. 
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
