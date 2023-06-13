#include <iostream>
#include <cmath>
#include <string>
#include "fstream"
#include "Vector.h"

using namespace std;

const double Tau = 0.5;
const double UTau = 1/Tau;
const double UmUTau=1-1/Tau;
//--------
const double Epsilon0=1, Mu0=2;
const double Sigma0=0.0;
const double C=1.0/sqrt(Epsilon0*Mu0);

const double E00=1,B00=E00/C;

int main(){
    int V[12][3], V0[3]; /*V[xyz][p][i]*/  vector3D v[12],v0; //v[p][i]
    vector3D e[12], e0; //e[p][i][j]
    vector3D h[12], h0; //b[p][i][j]
    int i;
    double alpha0=5.0;


    v[0].cargue(1,0,0);
    v[1].cargue(0,-1,0);
    v[2].cargue(0,0,-1);
    v[3].cargue(-1,0,0);
    v[4].cargue(0,1,0);
    v[5].cargue(0,0,1);
    //Electric vectors
    e[0].cargue(0,-1,0);
    e[1].cargue(0,0,1);
    e[2].cargue(-1,0,0);
    e[3].cargue(0,0,-1);
    e[4].cargue(1,0,0);
    e[5].cargue(0,1,0);
    //Magnetic vectors
    h[0].cargue(0,0,-1);
    h[1].cargue(-1,0,0);
    h[2].cargue(0,1,0);
    h[3].cargue(0,-1,0);
    h[4].cargue(0,0,-1);
    h[5].cargue(-1,0,0);
    for(i=6; i < 12;i++){
        v[i] = -1*v[i-6];
        h[i] = -1*h[i-6];
        e[i] = e[i-6];
    }
    V0[0] = 0; V0[1] = 0; V0[2]=0;
    for(i=0;i < 12; i++){
        V[i][0] = v[i].x();
        V[i][1] = v[i].y();
        V[i][2] = v[i].z();
    }

    vector3D f[10];

    vector3D E0,B,D;
    B.cargue(0,B00*exp(-0.25*pow(0,2)/(alpha0*alpha0)),0);
    E0.cargue(E00*exp(-0.25*pow(0,2)/(alpha0*alpha0)),0,0);
    D = E0*Epsilon0;
    vector3D Aux1,Aux2;
    Aux1 = D - ((v[2]^B)*1/(Mu0));
    Aux2 = B + ((v[2]^D)*1/(Epsilon0));
    Aux1.show();
    Aux2.show();



    return 0;
}
