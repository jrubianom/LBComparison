#include <algorithm>
#include <cstdlib>
#include<iostream>
#include<cmath>
#include <fstream>
#include <vector>
#include <functional>
#include "LBEM.h"

using namespace std;

//------------------------CONSTANTS-------------------------------
const int Lx = 100;   //
const int Ly = 100;   //
const int Lz = 100; //
const int Qr = 2, Qi = 6;
//-------------------
//-------------------
const double Epsilon0=3, Mu0=3;
const double Sigma0=0.0;
const double C=1.0/sqrt(Epsilon0*Mu0);

const double E00=0.001,B00=E00/C;
const double J0=0.0001;

const double T = 17.68/100.0*Lz/C, omega = 2*M_PI/T;
const double lambda = C*T, k = omega/C, alpha = 0.5;
const double Z0 = sqrt(Mu0/Epsilon0);

int main(){

  Parameter Params(Lx,Ly,Lz,Epsilon0,Mu0,Sigma0,E00,B00,J0,alpha,T);
  LatticeBoltzmann Dipole(Params);
  double R = lambda;
  int t, tmax=T*(70/25.0);//(R+lambda)/C;
  int N=100;

  int Ntheta = 100, Nphi = 100;
  vector<vector<double>> Ss,SsCurrent,SsPlus;
  Set(Ss,Ntheta,Nphi);
  Set(SsCurrent,Ntheta,Nphi);
  Set(SsPlus,Ntheta,Nphi);

  Dipole.Start();
  cout << "maximum number of iteration:" << tmax << endl;
  for(t=0;t<tmax;t++){
    Dipole.UpdateTime(t);
    cout<< "Iteracion:\t"<< t+1 <<endl;
    Dipole.Collision();
    Dipole.Advection();

    if(R <= C*t && C*t <= R+lambda)
      Dipole.PowerByPlanes(R,Ntheta,Nphi,Ss,SsCurrent,SsPlus);
  }

  Dipole.Print();
  //Print Power data
  Dipole.PrintPowerPlanes(Ss,Ntheta,Nphi,R);

  return 0;
}
