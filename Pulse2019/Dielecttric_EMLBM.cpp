#include <iostream>
#include <cmath>
#include <string>
#include "fstream"
#include "Vector.h"

using namespace std;


const int Lx=1;
const int Ly=1;
const int Lz=200;
///-----------
const double Tau = 0.5;
const double UTau = 1/Tau;
const double UmUTau=1-1/Tau;
//--------
const double Epsilon0=1, Mu0=1;
const double Sigma0=0.0;
const double C=1.0/sqrt(Epsilon0*Mu0);

const double E00=1,B00=E00/(C*Mu0);

//------------------Electromagnetic Constants for the Media------------------------------
double mur(int ix,int iy,int iz){
    return 1.0;
}
double epsilonr(int ix,int iy, int iz){
    return 1.5+0.5*tanh((double)(iz-(Lz/2.0)));//3-tanh((double)(iz-(Lz/2.0)));
}
double sigma(int ix,int iy,int iz){
    return 0.0;
}

//--------------------- class LatticeBoltzmann ------------
class LatticeBoltzmann{
  private:
    int n = 12;
    vector3D v[12];
    vector3D e[12];
    vector3D h[12];
    double f[Lx][Ly][Lz][12],fnew[Lx][Ly][Lz][12];
    vector3D P[Lx][Ly][Lz];
public:
    LatticeBoltzmann(void);
    //Fields from direct sums
    double  rhoc(int ix,int iy,int iz,bool UseNew);
    vector3D E(int ix,int iy,int iz,bool UseNew);
    vector3D H(int ix,int iy,int iz,bool UseNew);
    //Equilibrium Functions
    double feq(vector3D & E,vector3D & H,int i);
    //Simulation Functions
    void Start(void);
    void Collision(void);
    void ImposeFields(int t);
    void Advection(void);
    void Print(void);
};

LatticeBoltzmann::LatticeBoltzmann(void){
  int ix,iy,iz,alpha,r,p,i,j;
  //Velocity vectors V[p][i]=V^p_i (in components)
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
  for(int i=6; i < 12;i++){
    v[i] = -1*v[i-6];
    h[i] = -1*h[i-6];
    e[i] = e[i-6];
  }
}

//-----------------MACROSCOPIC FIELDS------------------
//Fields from direct sums
double LatticeBoltzmann::rhoc(int ix,int iy,int iz,bool UseNew){
  double sum = 0;
  //Start for the distribution for the central (zero) vector
  for(int i=0; i < n; i++)
    if(UseNew)
      sum+=fnew[ix][iy][iz][i];
    else
      sum+=f[ix][iy][iz][i];
  return sum;
}

vector3D LatticeBoltzmann::E(int ix,int iy,int iz,bool UseNew){
  vector3D sum; sum.cargue(0,0,0);
  for(int i=0;i < n;i++)
    if(UseNew)
	    sum += e[i]*fnew[ix][iy][iz][i];
   else
	    sum += e[i]*f[ix][iy][iz][i];
  return sum;
}
vector3D LatticeBoltzmann::H(int ix,int iy,int iz,bool UseNew){
  vector3D sum; sum.cargue(0,0,0);
  for(int i=0;i < n;i++)
    if(UseNew)
	    sum += h[i]*fnew[ix][iy][iz][i];
   else
	    sum += h[i]*f[ix][iy][iz][i];
  return sum;
}
//---------------EQUILIBRIUM FUNCTIONS-------------
double LatticeBoltzmann::feq(vector3D & E,vector3D & H,int i){
  double aux = E*e[i] + H*h[i];
  return aux/4.0;
}


//-------------------SIMULATION FUNCTIONS ----------------------------
void LatticeBoltzmann::Start(void){
  int ix,iy,iz,i;
  vector3D E0,H0,Jprima0,P0;
  double sigma0,mur0,epsilonr0;
  double alpha0=5.0,iz0=40;
  for(ix=0;ix<Lx;ix++) //para cada celda
    for(iy=0;iy<Ly;iy++)
      for(iz=0;iz<Lz;iz++){
        //Compute the constants
        sigma0=sigma(ix,iy,iz); mur0=mur(ix,iy,iz); epsilonr0=epsilonr(ix,iy,iz);
        //Impose the fields
        //rhoc0=0; Jprima0.cargue(0,0,0);
        H0.cargue(0,B00*exp(-0.25*pow(iz-iz0,2)/(alpha0*alpha0)),0);
        E0.cargue(E00*exp(-0.25*pow(iz-iz0,2)/(alpha0*alpha0)),0,0);
        //H0.cargue(0,sin(iz*2*M_PI/Lz),0);
        //E0.cargue(sin(iz*2*M_PI/Lz),0,0);
        P0 = (epsilonr0 - 1)*E0;
        P[ix][iy][iz] = P0;
        E0 = (E0 + P0)/epsilonr0;
        //Impose f=fnew=feq with the desired fields
        for(i=0;i < n; i++)
          fnew[ix][iy][iz][i]=f[ix][iy][iz][i]=
            feq(E0,H0,i);
      }
}


void LatticeBoltzmann::Collision(void){
    int ix,iy,iz,i; double sigma0,mur0,epsilonr0,prefactor0;
    double rhoc0; vector3D D0,B0,E0,H0,Jprima0,Eprima0,P0;
    for(ix=0;ix<Lx;ix++) //para cada celda
        for(iy=0;iy<Ly;iy++)
            for(iz=0;iz<Lz;iz++){
                //Compute the constants
                sigma0=sigma(ix,iy,iz); mur0=mur(ix,iy,iz); epsilonr0=epsilonr(ix,iy,iz);
                //Compute the fields
                E0=E(ix,iy,iz,false); H0=H(ix,iy,iz,false);
                P0 = 2*(epsilonr0-1)*E0 - P[ix][iy][iz];
                E0 = (E0 + P0)/epsilonr0;
                P[ix][iy][iz] = P0;
                //BGK evolution rule
                for(i=0; i < n; i++)
                  fnew[ix][iy][iz][i]=UmUTau*f[ix][iy][iz][i]
                    +UTau*feq(E0,H0,i);
            }
}

void LatticeBoltzmann::ImposeFields(int t){

}

void LatticeBoltzmann::Advection(void){
  int ix,iy,iz,i,ixnew,iynew,iznew;
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++)
      for(iz=0;iz<Lz;iz++){//for each cell
        for(i=0;i<n;i++){
          ixnew=(ix+int(v[i].x())+Lx)%Lx; iynew=(iy+int(v[i].y())+Ly)%Ly; iznew=(iz+int(v[i].z())+Lz)%Lz;
          f[ixnew][iynew][iznew][i]=fnew[ix][iy][iz][i];
        }
      }
}



void LatticeBoltzmann::Print(void){
  int ix=0,iy=0,iz,r,p,i,j; double sigma0,mur0,epsilonr0,prefactor0;
  double rhoc0; vector3D D0,B0,E0,H0,Jprima0,Eprima0,P0; double E2,B2;
  for(iz=0;iz<Lz;iz++){
    //Compute the electromagnetic constants
    sigma0=sigma(ix,iy,iz); mur0=mur(ix,iy,iz); epsilonr0=epsilonr(ix,iy,iz);
    //Compute the Fields
    E0=E(ix,iy,iz,true); H0=H(ix,iy,iz,true);
    P0 = 2*(epsilonr0-1)*E0 - P[ix][iy][iz];
    //P[ix][iy][iz] = P0;
    //E0 = (E0 + P0)/epsilonr0;
    //Print
    E2=norma2(E0); B2=norma2(H0);
    cout<<iz<<" "<<0.5*(epsilonr0*E2+B2/mur0)<<endl;
  }
}

int main(){
  LatticeBoltzmann DielectricPulse;
  int t, tmax= 200;

  DielectricPulse.Start();
  DielectricPulse.ImposeFields(0);

  for(t=0;t<tmax;t++){
    DielectricPulse.Collision();
    DielectricPulse.ImposeFields(t);
    DielectricPulse.Advection();
  }

  DielectricPulse.Print();

  return 0;
}
