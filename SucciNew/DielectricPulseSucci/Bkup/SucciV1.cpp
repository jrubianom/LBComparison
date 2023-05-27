#include <iostream>
#include <cmath>
#include <string>
#include "fstream"
#include "Vector.h"

using namespace std;

//------------------------CONSTANTS-------------------------------
const int scale = 3;

const int Lx=1;
const int Ly=1;
const int Lz=200*scale;

const int Qr = 2, Qi = 7;
//-------------------
const double Epsilon0 = 1, Mu0= 2;
const double Sigma0=0.0;
const double C=1.0/sqrt(Epsilon0*Mu0);
const double E00=0.001,B00=E00/C;

const double epsr1 = 1, epsr2 = 2.5;

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

double c2(int ix, int iy, int iz){
  double epsr0 = epsilonr(ix,iy,iz);
  double mur0 = mur(ix,iy,iz);
  return 1.0/(Mu0*Epsilon0*mur0*epsr0);
}

//---------Auxiliar funcion
double levi(int i,int j,int k){
  return ((i-j)*(j-k)*(k-i))/2.0;
}

//--------------------- class LatticeBoltzmann ------------
class LatticeBoltzmann{
  private:
    int D=3;
    int V[Qi][3];   vector3D v[Qi];
    double *f=nullptr,*fnew=nullptr;//f[ix][iy][iz][r][i]
  public:
    LatticeBoltzmann();
    ~LatticeBoltzmann();
    double Ccell(double epsilonr0,double mur0){return C/sqrt(epsilonr0*mur0);};
    int index(int ix,int iy,int iz,int i,int alpha,int beta);
    double weight_0(int ix, int iy, int iz){ return 1.0 - 3*c2(ix, iy, iz);};
    double weight_i(int ix, int iy, int iz){return (1.0 - weight_0(ix, iy, iz))/(2*D);};
    double Lambda(vector3D &E, int alpha,int beta);

    //Fields from direct sums
    vector3D E(int ix,int iy,int iz,bool UseNew);
    vector3D B(int ix,int iy,int iz,bool UseNew);
    //Equilibrium Functions
    double feq(vector3D & E,vector3D & B,double c_sqrd,
                               double w_i,int i,int alpha,int beta);
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
  v[0].cargue(0,0,0);
  v[1].cargue(1,0,0);
  v[2].cargue(-1,0,0);
  v[3].cargue(0,1,0);
  v[4].cargue(0,-1,0);
  v[5].cargue(0,0,1);
  v[6].cargue(0,0,-1);


  for(i=0;i < Qi; i++){
    V[i][0] = v[i].x();
    V[i][1] = v[i].y();
    V[i][2] = v[i].z();
  }

  f = new double[Lx*Ly*Lz*Qi*D*D]; fnew = new double[Lx*Ly*Lz*Qi*D*D];
}

LatticeBoltzmann::~LatticeBoltzmann(void){
  delete[] f;  delete[] fnew;
}

int LatticeBoltzmann::index(int ix,int iy,int iz,int i,int alpha,int beta){
  return (iz*Lx*Ly+iy*Lx+ix)*Qi*D*D + (i*D*D + alpha*D + beta);
}

double LatticeBoltzmann::Lambda(vector3D &E, int alpha,int beta){
  double sumaux = 0;
  for (int gamma = 0; gamma < 3; gamma++){
    sumaux += levi(gamma,alpha,beta)*E.entry(gamma);
  }
  return -sumaux;
}

//-----------------MACROSCOPIC FIELDS------------------
//Fields from direct sums
vector3D LatticeBoltzmann::E(int ix,int iy,int iz,bool UseNew){
  int i,id;
  double sum;
  vector3D Eaux;
  double Evec[3];
  double c_sqrd = c2(ix,iy,iz);
  for(int gamma=0; gamma < 3; gamma++){
    sum = 0;
    for(int alpha=0; alpha < 3; alpha++)
      for(int beta=0;beta < 3; beta++)
        for(i=0;i<Qi;i++){
          id = index(ix,iy,iz,i,alpha,beta);
          if(UseNew)
            sum+=levi(gamma,alpha,beta)*fnew[id];
          else
            sum+=levi(gamma,alpha,beta)*f[id];
        }
    Evec[gamma] = -c_sqrd*0.5*sum;
  }
  Eaux.cargue(Evec[0],Evec[1],Evec[2]);
  return Eaux;
}
vector3D LatticeBoltzmann::B(int ix,int iy,int iz,bool UseNew){
  int i,id;
  double sum;
  vector3D Baux;
  double Bvec[3];
  for(int tau = 0; tau < 3; tau++){
    sum = 0;
    for(int lamb = 0; lamb < 3; lamb++)
      for(int gamma=0; gamma < 3; gamma++)
        for(int alpha=0; alpha < 3; alpha++)
          for(int beta=0;beta < 3; beta++)
            for(i=0;i<Qi;i++){
              id = index(ix,iy,iz,i,alpha,beta);
              if(UseNew)
                sum+=levi(tau,lamb,gamma)*levi(gamma,alpha,beta)*V[i][lamb]*fnew[id];
              else
                sum+=levi(tau,lamb,gamma)*levi(gamma,alpha,beta)*V[i][lamb]*f[id];
            }
    Bvec[tau] = -sum/4.0;
  }
  Baux.cargue(Bvec[0],Bvec[1],Bvec[2]);
  return Baux;
}

//---------------EQUILIBRIUM FUNCTIONS-------------
double LatticeBoltzmann::feq(vector3D & E,vector3D & B,double c_sqrd,
                               double w_i,int i,int alpha,int beta){
  return w_i/c_sqrd*(Lambda(E,alpha,beta) + V[i][alpha]*B.entry(beta)-V[i][beta]*B.entry(alpha));
}


//-------------------SIMULATION FUNCTIONS ----------------------------
void LatticeBoltzmann::Start(void){
  int ix,iy,iz,i,r,id;
  vector3D E0,B0,D0;
  double sigma0,mur0,epsilonr0,speedCell,c_sqrd,w_i;
  double alpha0=Lz*0.5/(10*sqrt(2)),iz0= Lz/2.0 - Lz/6.0;
  for(ix=0;ix<Lx;ix++) //para cada celda
    for(iy=0;iy<Ly;iy++)
      for(iz=0;iz<Lz;iz++){
        //Compute the constants
        sigma0=sigma(ix,iy,iz); mur0=mur(ix,iy,iz); epsilonr0=epsilonr(ix,iy,iz);
        //speedCell = Ccell(epsilonr0,mur0);
        c_sqrd = c2(ix,iy,iz);
        speedCell = sqrt(c_sqrd*2);
        B0.cargue(0,(E00/speedCell)*exp(-pow(iz-iz0,2)/(2*pow(alpha0,2))),0);
        E0.cargue(E00*exp(-pow(iz-iz0,2)/(2*pow(alpha0,2))),0,0);

        for(int alpha=0; alpha < 3; alpha++)
          for(int beta=0; beta < 3; beta++)
            for(i=0;i<Qi;i++){
              if(i==0)
                w_i = weight_0(ix,iy,iz);
              else
                w_i = weight_i(ix,iy,iz);

              id = index(ix,iy,iz,i,alpha,beta);
              fnew[id]=f[id]=feq(E0,B0,c_sqrd,w_i,i,alpha,beta);
            }
      }
}

void LatticeBoltzmann::Collision(void){
  int ix,iy,iz,i,r,id;
  double mur0,epsilonr0,sigma0,c_sqrd,w_i;;
  vector3D E0,B0;
  for(ix=0;ix<Lx;ix++) //para cada celda
    for(iy=0;iy<Ly;iy++)
      for(iz=0;iz<Lz;iz++){
        //Compute the constants
        sigma0=sigma(ix,iy,iz); mur0=mur(ix,iy,iz); epsilonr0=epsilonr(ix,iy,iz);
        c_sqrd = c2(ix,iy,iz);
        //Compute the fields
        E0=E(ix,iy,iz,false); B0=B(ix,iy,iz,false);
        //E0 = E(D0,epsilonr0); H0 = H(B0,mur0);
        //BGK evolution rule
        for(int alpha=0; alpha < 3; alpha++)
          for(int beta=0; beta < 3; beta++)
            for(i=0; i < Qi; i++){
              if(i==0)
                w_i = weight_0(ix,iy,iz);
              else
                w_i = weight_i(ix,iy,iz);

              id = index(ix,iy,iz,i,alpha,beta);
              fnew[id]=2*feq(E0,B0,c_sqrd,w_i,i,alpha,beta)-f[id];
            }

      }
}

void LatticeBoltzmann::ImposeFields(int t){

}

void LatticeBoltzmann::Advection(void){
  int ix,iy,iz,r,i,ixnew,iynew,iznew;
  int id,idnew;
  for(ix=0;ix<Lx;ix++) 
    for(iy=0;iy<Ly;iy++)
      for(iz=0;iz<Lz;iz++){//for each cell
        for(int alpha=0; alpha < 3; alpha++)
          for(int beta=0; beta < 3; beta++)
            for(i=0;i<Qi;i++){
              ixnew=(ix+V[i][0]+Lx)%Lx; iynew=(iy+V[i][1]+Ly)%Ly; iznew=(iz+V[i][2]+Lz)%Lz;
              id = index(ix,iy,iz,i,alpha,beta);
              idnew=index(ixnew,iynew,iznew,i,alpha,beta);
              f[idnew]=fnew[id];
            }
      }
}
/*
void LatticeBoltzmann::Print(void){
  int ix=0,iy=0,iz;
  vector3D B0,E0;;
  ofstream Myfile("data.dat");
  for(iz=0;iz<Lz;iz++){
    E0=E(ix,iy,iz,true); B0=B(ix,iy,iz,true);
    Myfile<<iz<<" "<<E0.entry(0)/E00<<endl;
  }
  Myfile.close();
}
*/

void LatticeBoltzmann::Print(void){
  int ix=0,iy=0,iz;
  vector3D B0,E0;
  ofstream Myfile("data.dat");
  double Et, Er,EAmp;
  for(iz=0;iz<Lz;iz++){
    //Compute the Fields
    E0=E(ix,iy,iz,true); B0=B(ix,iy,iz,true);
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
  int t, tmax=Lz/(2*C)*1.4;

  DielectricPulse.Start();

  for(t=0;t<tmax;t++){
    DielectricPulse.Collision();
    DielectricPulse.Advection();
  }

  DielectricPulse.Print();

  return 0;
}
