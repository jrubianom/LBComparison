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
const int Qr = 2, Qp = 3, Qi = 4, Qj = 2;
//-------------------
const double Tau = 0.5;
const double UTau = 1/Tau;
const double UmUTau=1-1/Tau;
//-------------------
const double Epsilon0=1, Mu0=2;
const double Sigma0=0.0;
const double C=1.0/sqrt(2.0);

const double E00=0.001,B00=E00/C;
const double J0=0.0001;

const double Z0=sqrt(Mu0/Epsilon0);

int main(){

  int Lx, Ly, Lz;

  ofstream fileTime("fileTime.dat");

  // Loop to refine the mesh
  for(int i=1; i<4; i++){
    // Domain size
    int size = 100+10*(i-1);
    Lx = Ly = Lz = size;
    // Constants
    double T = 17.68/100.0*Lz/C, omega = 2*M_PI/T;
    double lamd = C*T, k = omega/C, alpha = 0.5;

    for(int iteration=0; iteration<1; iteration++){
      Parameter Params(Lx,Ly,Lz,Tau,Epsilon0,Mu0,Sigma0,E00,B00,J0,alpha,T);
      LatticeBoltzmann Dipole(Params);

      int t, tmax=T*(70/25.0);

      // Time measurement
      struct timespec begin, end;

      Dipole.Start();
      cout << "maximum number of iteration:" << tmax << endl;

      clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &begin);
      for(t=0;t<tmax;t++){
        Dipole.UpdateTime(t);
        cout<< "Iteracion:\t"<< t+1 <<endl;
        Dipole.Collision();
        Dipole.Advection();
      }
      clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);

      long seconds = end.tv_sec - begin.tv_sec;
      long nanoseconds = end.tv_nsec - begin.tv_nsec;
      double ellapsed = seconds + nanoseconds*1e-9;
      fileTime << i << " " << ellapsed << endl;

      // Print fields
      string stringB = "fileB_"+to_string(i)+".dat";
      ofstream fileB(stringB);
      string stringE = "fileE_"+to_string(i)+".dat";
      ofstream fileE(stringE);
      string stringContour = "contour_"+to_string(i)+".dat";
      ofstream contour(stringContour);

      Dipole.Print(fileB, fileE, contour);

    }
  }
  
  fileTime.close();
  
  return 0;
}