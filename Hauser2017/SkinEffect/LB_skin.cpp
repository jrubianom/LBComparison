#include <fstream>
#include <iostream>
#include <cmath>
#include <string>
#include "fstream"
#include "Vector.h"

using namespace std;

//------------------------CONSTANTS-------------------------------
const int scale = 1;

const int Lx=1;
const int Ly=1;
const int Lz=1000*scale;

const int Qr = 2, Qi = 6;
//-------------------
const double Epsilon0 = 3, Mu0= 3;
const double Sigma0=0.0125;
const double C=1.0/sqrt(Epsilon0*Mu0);
const double E00=1,B00=E00/C;

const double epsr1 = 2, epsr2 = 2.5;

//------------------Electromagnetic Constants for the Media------------------------------
double mur(int ix,int iy,int iz){
  return 1.0;
}
double epsilonr(int ix,int iy, int iz){
  return 1;//1.75 + 0.75*tanh(iz-Lz/2);
}
double sigma(int ix,int iy,int iz){
  return Sigma0*(1+tanh(iz-Lz/4.0));
}


//--------------------- class LatticeBoltzmann ------------
class LatticeBoltzmann{
  private:
    int V[Qi][3], V0[3];   vector3D v[Qi];
    vector3D *f=nullptr,*fnew=nullptr;//f[ix][iy][iz][r][i]
  public:
    LatticeBoltzmann();
    ~LatticeBoltzmann();
    double Ccell(double epsilonr0,double mur0){return C/sqrt(epsilonr0*mur0);};
    int index(int ix,int iy,int iz,int r,int i);
    int index0(int ix,int iy,int iz);
    //Fields from direct sums
    vector3D D(int ix,int iy,int iz,bool UseNew);
    vector3D B(int ix,int iy,int iz,bool UseNew);
    //Fields deduced from the first ones through electromagnetic constants
    vector3D E(vector3D & D0,double Epsilonr);
    vector3D H(vector3D & B0,double Mur);
    //Equilibrium Functions
    vector3D feq(vector3D & D,vector3D & B,int r,int i,
                             double epsr,double mur0);
    //Simulation Functions
    void Start(void);
    void Collision(void);
    void ImposeFields(int t);
    void Advection(void);
    void Print(void);
};

LatticeBoltzmann::LatticeBoltzmann(){
  int i;
  //Velocity vectors V[p][i]=V^p_i (in components)

  //Velocity vectors V[p][i]=V^p_i (in components)
  v[0].cargue(1,0,0);
  v[2].cargue(-1,0,0);
  v[1].cargue(0,1,0);
  v[3].cargue(0,-1,0);
  v[5].cargue(0,0,1);
  v[4].cargue(0,0,-1);

  V0[0] = 0; V0[1] = 0; V0[2]=0;
  for(i=0;i < Qi; i++){
    V[i][0] = v[i].x();
    V[i][1] = v[i].y();
    V[i][2] = v[i].z();
  }

  f = new vector3D[Lx*Ly*Lz*Qr*Qi]; fnew = new vector3D[Lx*Ly*Lz*Qr*Qi];
}

LatticeBoltzmann::~LatticeBoltzmann(void){
  delete[] f;  delete[] fnew;
}

int LatticeBoltzmann::index(int ix,int iy,int iz,int r,int i){
  return (iz*Lx*Ly+iy*Lx+ix)*Qr*Qi + (r*Qi + i);
}

int LatticeBoltzmann::index0(int ix,int iy,int iz){
  return (iz*Lx*Ly+iy*Lx+ix);
}

//-----------------MACROSCOPIC FIELDS------------------
//Fields from direct sums
vector3D LatticeBoltzmann::D(int ix,int iy,int iz,bool UseNew){
  int i,id; vector3D sum; sum.cargue(0,0,0);

  for(i=0;i<Qi;i++){
    id = index(ix,iy,iz,0,i);
    if(UseNew)
      sum+=fnew[id];
    else
      sum+=f[id];
  }
  return sum;
}
vector3D LatticeBoltzmann::B(int ix,int iy,int iz,bool UseNew){
  int i,id; vector3D sum; sum.cargue(0,0,0);
  for(i=0;i<Qi;i++){
    id = index(ix,iy,iz,1,i);
    if(UseNew)
      sum+=fnew[id];
    else
      sum+=f[id];
  }
  return sum;
}
//Fields deduced from the first ones through electromagnetic constants
vector3D LatticeBoltzmann::E(vector3D & D0,double Epsilonr){
  return D0*(1.0/(Epsilonr*Epsilon0));
}
vector3D LatticeBoltzmann::H(vector3D & B0,double Mur){
  return B0*(1.0/(Mur*Mu0));
}
//---------------EQUILIBRIUM FUNCTIONS-------------
vector3D LatticeBoltzmann::feq(vector3D & D,vector3D & B,int r,int i,
                               double epsr,double mur0){
  vector3D Aux;
  if(r == 0){
    Aux = D - 3*((v[i]^B)*1/(mur0*Mu0));
  }
  if(r == 1){
    Aux = B + 3*((v[i]^D)*1/(epsr*Epsilon0));
  }
  return Aux*(1.0/6.0);
}


//-------------------SIMULATION FUNCTIONS ----------------------------
void LatticeBoltzmann::Start(void){
  int ix,iy,iz,i,r,id;
  vector3D E0,B0,D0;
  double sigma0,mur0,epsilonr0,speedCell;
  double alpha0=Lz*0.5/(10*sqrt(2)),iz0=Lz/2.0 - Lz/6.0;
  for(ix=0;ix<Lx;ix++) //para cada celda
    for(iy=0;iy<Ly;iy++)
      for(iz=0;iz<Lz;iz++){
        //Compute the constants
        sigma0=sigma(ix,iy,iz); mur0=mur(ix,iy,iz); epsilonr0=epsilonr(ix,iy,iz);
        speedCell = Ccell(epsilonr0,mur0);

        B0.cargue(0,0,0);
        E0.cargue(0,0,0);
        D0 = epsilonr0*Epsilon0*E0*exp(-sigma0/(epsilonr0*Epsilon0));

        for(r=0;r<Qr;r++)
          for(i=0;i<Qi;i++){
            id = index(ix,iy,iz,r,i);
            fnew[id]=f[id]=feq(D0,B0,r,i,epsilonr0,mur0);
          }
      }
}

void LatticeBoltzmann::Collision(void){
  int ix,iy,iz,i,r,id;
  double mur0,epsilonr0,sigma0;
  vector3D D0,B0;
  for(ix=0;ix<Lx;ix++) //para cada celda
    for(iy=0;iy<Ly;iy++)
      for(iz=0;iz<Lz;iz++){
        //Compute the constants
        sigma0=sigma(ix,iy,iz); mur0=mur(ix,iy,iz); epsilonr0=epsilonr(ix,iy,iz);
        //Compute the fields
        D0=D(ix,iy,iz,false)*exp(-sigma0/(epsilonr0*Epsilon0));
        B0=B(ix,iy,iz,false);
        //E0 = E(D0,epsilonr0); H0 = H(B0,mur0);
        //BGK evolution rule
        for(r = 0; r < Qr; r++)
          for(i=0; i < Qi; i++){
            id = index(ix,iy,iz,r,i);
            fnew[id]=2*feq(D0,B0,r,i,epsilonr0,mur0)-f[id];
          }

      }
}

void LatticeBoltzmann::ImposeFields(int t){
  int ix,iy,iz,i,r,id;
  iz = 0;
  vector3D E0,B0,D0,H0,Jprima0,P0;
  double sigma0,mur0,epsilonr0,speedCell;
  double T,omega;
  T = 17.68/C; omega = 2*M_PI/T;
  for(ix=0;ix<Lx;ix++) //para cada celda
    for(iy=0;iy<Ly;iy++){
        //Compute the constants
        sigma0=sigma(ix,iy,iz); mur0=mur(ix,iy,iz); epsilonr0=epsilonr(ix,iy,iz);
        speedCell = Ccell(epsilonr0,mur0);
        //Impose the fields
        //rhoc0=0; Jprima0.cargue(0,0,0);
        B0.cargue(0,(E00/speedCell)*sin(omega*t),0);
        E0.cargue(E00*sin(omega*t),0,0);
        D0 = epsilonr0*Epsilon0*E0*exp(-sigma0/(epsilonr0*Epsilon0));
        //H0.cargue(0,sin(iz*2*M_PI/Lz),0);
        //E0.cargue(sin(iz*2*M_PI/Lz),0,0);
        //Impose f=fnew=feq with the desired fields
        for(r=0; r < Qr; r++)
          for(i=0;i < Qi; i++){
            id = index(ix,iy,iz,r,i);
            fnew[id]=f[id]=feq(D0,B0,r,i,epsilonr0,mur0);
          }
    }
}

void LatticeBoltzmann::Advection(void){
  int ix,iy,iz,r,i,ixnew,iynew,iznew;
  int id,idnew;
  for(ix=0;ix<Lx;ix++) 
    for(iy=0;iy<Ly;iy++)
      for(iz=0;iz<Lz;iz++){//for each cell
        for(r=0;r<Qr;r++)
          for(i=0;i<Qi;i++){
            ixnew=(ix+V[i][0]+Lx)%Lx; iynew=(iy+V[i][1]+Ly)%Ly; iznew=(iz+V[i][2]+Lz)%Lz;
            id = index(ix,iy,iz,r,i);idnew=index(ixnew,iynew,iznew,r,i);
            f[idnew]=fnew[id];
          }
      }
}

void LatticeBoltzmann::Print(void){
  int ix=0,iy=0,iz; double sigma0,mur0,epsilonr0;
  vector3D D0,B0,E0; double E2,B2,eps,mus;
  ofstream Myfile("data.dat");
  for(iz=0;iz<Lz/2.0;iz++){
    //Compute the electromagnetic constants
    sigma0=sigma(ix,iy,iz); mur0=mur(ix,iy,iz); epsilonr0=epsilonr(ix,iy,iz);
    //Compute the Fields
    D0=D(ix,iy,iz,true); B0=B(ix,iy,iz,true);
    E0=E(D0,epsilonr0)*exp(sigma0/(epsilonr0*Epsilon0));
    //Print
    E2=norma2(E0); B2=norma2(B0);
    eps = Epsilon0*epsilonr0;
    mus = Mu0*mur0;
    //cout<<iz<<" "<<0.5*(eps*E2+B2/mus)<<endl;
    Myfile<<iz<<" "<<E0.x()<<endl;
  }
  Myfile.close();
}


int main(){
  LatticeBoltzmann DielectricPulse;
  int t, tmax=Lz/(3*C);

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
