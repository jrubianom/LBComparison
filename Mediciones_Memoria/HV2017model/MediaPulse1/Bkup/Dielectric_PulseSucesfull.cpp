#include <iostream>
#include <cmath>
#include <string>
#include "fstream"
#include "Vector.h"

using namespace std;

const int scale = 1;

const int Lx=1;
const int Ly=1;
const int Lz=200*scale;
//--------
const double Epsilon0 = 3, Mu0= 3;
const double Sigma0=0.0;
const double C=1.0/sqrt(Epsilon0*Mu0);

const double E00=0.001,B00=E00/C;

const double epsr1 = 2, epsr2 = 2.5;

//------------------Electromagnetic Constants for the Media------------------------------
double mur(int ix,int iy,int iz){
  return 1.0;
}
double epsilonr(int ix,int iy, int iz){
  return (epsr1 + epsr2)/2.0 + (epsr2 - epsr1)/2*tanh(iz-Lz/2.0);
}
double sigma(int ix,int iy,int iz){
  return 0.0;
}

//--------------------- class LatticeBoltzmann ------------
class LatticeBoltzmann{
  private:
    int n=6,r0 = 2;
    int V[6][3], V0[3]; /*V[xyz][p][i]*/  vector3D v[6],v0; //v[p][i]
    vector3D e[6], e0; //e[p][i][j]
    vector3D h[6], h0; //b[p][i][j]
    // 12 = non zero directions. 2 = electric or magnetic
    vector3D f[Lx][Ly][Lz][2][6],fnew[Lx][Ly][Lz][2][6];//f[ix][iy][iz][r][p][i][j]
    double f0[Lx][Ly][Lz][2],f0new[Lx][Ly][Lz][2];//f0[ix][iy][iz] (r=0)
  public:
    LatticeBoltzmann(void);
    //Speed of light in cell
    double Ccell(double epsilonr0,double mur0){return C/sqrt(epsilonr0*mur0);};
    //Fields from direct sums
    double  rhoc(int ix,int iy,int iz,bool UseNew);
    vector3D D(int ix,int iy,int iz,bool UseNew);
    vector3D B(int ix,int iy,int iz,bool UseNew);
    vector3D E(vector3D & D0,double Epsilonr){return D0*(1.0/(Epsilonr*Epsilon0));};
    vector3D H(vector3D & B0,double Mur){return B0*(1.0/(Mur*Mu0));};
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

LatticeBoltzmann::LatticeBoltzmann(void){
  int i;
  //Velocity vectors V[p][i]=V^p_i (in components)
  v[0].cargue(1,0,0);
  v[2].cargue(-1,0,0);
  v[1].cargue(0,1,0);
  v[3].cargue(0,-1,0);
  v[5].cargue(0,0,1);
  v[4].cargue(0,0,-1);

  V0[0] = 0; V0[1] = 0; V0[2]=0;
  for(i=0;i < n; i++){
    V[i][0] = v[i].x();
    V[i][1] = v[i].y();
    V[i][2] = v[i].z();
  }
}

//-----------------MACROSCOPIC FIELDS------------------
//Fields from direct sums
vector3D LatticeBoltzmann::D(int ix,int iy,int iz,bool UseNew){
  vector3D sum; sum.cargue(0,0,0);
  for(int i=0;i < n;i++)
    if(UseNew)
	    sum += fnew[ix][iy][iz][0][i];
   else
	    sum += f[ix][iy][iz][0][i];
  return sum;
}
vector3D LatticeBoltzmann::B(int ix,int iy,int iz,bool UseNew){
  vector3D sum; sum.cargue(0,0,0);
  for(int i=0;i < n;i++)
    if(UseNew)
	    sum += fnew[ix][iy][iz][1][i];
   else
	    sum += f[ix][iy][iz][1][i];
  return sum;
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
  int ix,iy,iz,i,r;
  vector3D E0,B0,D0,H0,Jprima0,P0;
  double sigma0,mur0,epsilonr0,speedCell;
  double alpha0=Lz*0.5/(10*sqrt(2)),iz0=Lz/2.0 - Lz/6.0;
  for(ix=0;ix<Lx;ix++) //para cada celda
    for(iy=0;iy<Ly;iy++)
      for(iz=0;iz<Lz;iz++){
        //Compute the constants
        sigma0=sigma(ix,iy,iz); mur0=mur(ix,iy,iz); epsilonr0=epsilonr(ix,iy,iz);
        speedCell = Ccell(epsilonr0,mur0);
        //Impose the fields
        //rhoc0=0; Jprima0.cargue(0,0,0);
        B0.cargue(0,(E00/speedCell)*exp(-pow(iz-iz0,2)/(2*pow(alpha0,2))),0);
        E0.cargue(E00*exp(-pow(iz-iz0,2)/(2*pow(alpha0,2))),0,0);
        D0 = epsilonr0*Epsilon0*E0;
        //H0.cargue(0,sin(iz*2*M_PI/Lz),0);
        //E0.cargue(sin(iz*2*M_PI/Lz),0,0);
        //Impose f=fnew=feq with the desired fields
        for(r=0; r < 2; r++)
          for(i=0;i < n; i++)
            fnew[ix][iy][iz][r][i]=f[ix][iy][iz][r][i]=feq(D0,B0,r,i,epsilonr0,mur0);
      }
}


void LatticeBoltzmann::Collision(void){
    int ix,iy,iz,i,r; double sigma0,mur0,epsilonr0,prefactor0;
    double rhoc0; vector3D D0,B0,E0,H0,Jprima0,Eprima0,P0;
    vector3D aux;
    for(ix=0;ix<Lx;ix++) //para cada celda
        for(iy=0;iy<Ly;iy++)
            for(iz=0;iz<Lz;iz++){
                //Compute the constants
                sigma0=sigma(ix,iy,iz); mur0=mur(ix,iy,iz); epsilonr0=epsilonr(ix,iy,iz);
                //Compute the fields
                D0=D(ix,iy,iz,false); B0=B(ix,iy,iz,false);
                //E0 = E(D0,epsilonr0); H0 = H(B0,mur0);
                //BGK evolution rule
                for(r = 0; r < 2; r++)
                  for(i=0; i < n; i++)
                    fnew[ix][iy][iz][r][i]=2*feq(D0,B0,r,i,epsilonr0,mur0)-f[ix][iy][iz][r][i];
                    //aux=2*feq(D0,B0,r,i,epsilonr0,mur0)-f[ix][iy][iz][r][i];
                //aux.show();
                    //fnew[ix][iy][iz][r][i]=f[ix][iy][iz][r][i];
            }
}

void LatticeBoltzmann::ImposeFields(int t){

}

void LatticeBoltzmann::Advection(void){
  int ix,iy,iz,i,r,ixnew,iynew,iznew;
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++)
      for(iz=0;iz<Lz;iz++){//for each cell
        for(r=0;r < 2;r++)
          for(i=0;i<n;i++){
            ixnew=(ix+V[i][0]+Lx)%Lx; iynew=(iy+V[i][1]+Ly)%Ly; iznew=(iz+V[i][2]+Lz)%Lz;
            f[ixnew][iynew][iznew][r][i]=fnew[ix][iy][iz][r][i];
          }
      }
}


void LatticeBoltzmann::Print(void){
  int ix=0,iy=0,iz,r,p,i,j; double sigma0,mur0,epsilonr0,prefactor0;
  double rhoc0; vector3D D0,B0,E0,H0,Jprima0,Eprima0,P0; double E2,B2,eps,mus;
  ofstream Myfile("data.dat");
  double Et, Er,EAmp;
  for(iz=0;iz<Lz;iz++){
    //Compute the electromagnetic constants
    sigma0=sigma(ix,iy,iz); mur0=mur(ix,iy,iz); epsilonr0=epsilonr(ix,iy,iz);
    //Compute the Fields
    D0=D(ix,iy,iz,true); B0=B(ix,iy,iz,true);
    E0=E(D0,epsilonr0);
    //Print
    E2=norma2(E0); B2=norma2(B0);
    eps = Epsilon0*epsilonr0;
    mus = Mu0*mur0;
    //cout<<iz<<" "<<0.5*(eps*E2+B2/mus)<<endl;
    Myfile<<iz<<" "<<E0.x()/E00<<endl;

    EAmp = abs(E0.x());

    if(iz==0){
      Et = Er = EAmp;
    }
    else{
      if(iz < Lz/2 && EAmp > Er)
        Er = EAmp;
      if(iz > Lz/2 && EAmp > Et)
        Et = EAmp;
    }
  }
  Myfile.close();
  Er = Er/E00; Et = Et/E00;
  double Ertheo,Ettheo,ratio = sqrt(epsr2/epsr1);
  Ertheo = abs((ratio -1 )/(ratio + 1));
  Ettheo = abs(2/(ratio + 1));
  cout << "Simulated:\n";
  cout << Er << "\t"  << Et << endl;
  cout << "Theoretical:\n";
  cout << Ertheo << "\t"  << Ettheo << endl;
  cout << "Realtive Error Er and Et %" << endl;
  cout << abs(100*(Ertheo - Er)/Ertheo) << "\t" << abs(100*(Ettheo - Et)/Ettheo) << endl;
}

int main(){
  LatticeBoltzmann DielectricPulse;
  int t, tmax=Lz/(2*C);

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
