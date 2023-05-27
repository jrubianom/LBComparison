#include <iostream>
#include <cmath>
#include <string>
#include "fstream"
#include "Vector.h"
#include <algorithm>
#include <vector>

using namespace std;

//------------------------CONSTANTS-------------------------------
const int scale = 1;

const int Lx=1;
const int Ly=1;
const int Lz=500*scale;

const int Qr = 2, Qi = 7;
//-------------------
const double Epsilon0 = 1, Mu0= 2;
const double Sigma0=0.0125/3;
const double C=1.0/sqrt(Epsilon0*Mu0);
const double E00=1,B00=E00/C;

const double epsr1 = 1, epsr2 = 2.5;

const double Xfactor=sqrt(2);

const double T = 17.68/C*Xfactor, omega = 2*M_PI/T;

//------------------Electromagnetic Constants for the Media------------------------------
double mur(int ix,int iy,int iz){
  return 1.0;
}
double epsilonr(int ix,int iy, int iz){
  return 1;//(epsr1 + epsr2)/2.0 + (epsr2 - epsr1)/2*tanh(iz-Lz/2.0);
}
double sigma(int ix,int iy,int iz){
  return Sigma0*(1+tanh(iz-Lz/4.0));
}
double c2(int ix,int iy,int iz){
  double epsr0 = epsilonr(ix,iy,iz);
  double mur0 = mur(ix,iy,iz);
  return 1.0/(Mu0*Epsilon0*mur0*epsr0);
}
double c2Prime(int ix, int iy, int iz){
  return 3*c2(ix,iy,iz);
}

//---------Auxiliar funcion
double levi(int i,int j,int k){
  return ((i-j)*(j-k)*(k-i))/2.0;
}

void Amplitud(vector<double>& Amplitud, vector<double>& oldAmplitud);

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
    double weight_0(int ix, int iy, int iz){ return 1.0 - c2Prime(ix, iy, iz);};
    double weight_i(int ix, int iy, int iz){return (1.0 - weight_0(ix, iy, iz))/(2*D);};
    double Lambda(vector3D &E, int alpha,int beta);

    //Fields from direct sums
    vector3D E(int ix,int iy,int iz,bool UseNew);
    vector3D B(int ix,int iy,int iz,bool UseNew);
    //Force terms
    double S(vector3D & J,vector3D & B,double c_sqrdPrime,int alpha);
    double T(vector3D & J,vector3D & B,double c_sqrdPrime,
             double w_i,int i,int alpha,int beta);
    //Equilibrium Functions
    double feq(vector3D & E,vector3D & B,double c_sqrdPrime,
               double w_i,int i,int alpha,int beta);
    //Simulation Functions
    void Start(void);
    void Collision(void);
    void ImposeFields(int t);
    void Advection(void);

    //Print
    void Print(vector<double>& EAmplit);
    void MeasureAmplitude(vector<double> &EAmplit);

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
  double c_sqrdPrime = c2Prime(ix,iy,iz);
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
    Evec[gamma] = -c_sqrdPrime*0.5*sum;
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
    Bvec[tau] = -3*sum/4.0;
  }
  Baux.cargue(Bvec[0],Bvec[1],Bvec[2]);
  return Baux;
}
//---------------FORCE TERMS-----------------------
double LatticeBoltzmann::S(vector3D & J,vector3D & B,double c_sqrdPrime,
                           int alpha){
  return -J.entry(alpha)/Epsilon0;
}

double LatticeBoltzmann::T(vector3D & J,vector3D & B,double c_sqrdPrime,
                           double w_i,int i,int alpha,int beta){
  double sum=0;
  for(int gamma=0; gamma < 3; gamma++){
    sum += levi(gamma,alpha,beta)*S(J,B,c_sqrdPrime,gamma);
  }
  return -w_i*sum/c_sqrdPrime;
}

//---------------EQUILIBRIUM FUNCTIONS-------------
double LatticeBoltzmann::feq(vector3D & E,vector3D & B,double c_sqrdPrime,
                               double w_i,int i,int alpha,int beta){
  return w_i/c_sqrdPrime*(Lambda(E,alpha,beta) + V[i][alpha]*B.entry(beta)-V[i][beta]*B.entry(alpha));
}


//-------------------SIMULATION FUNCTIONS ----------------------------
void LatticeBoltzmann::Start(void){
  int ix,iy,iz,i,r,id;
  vector3D E0,B0,D0;
  double sigma0,mur0,epsilonr0,speedCell,c_sqrdPrime,w_i;
  double alpha0=Lz*0.5/(10*sqrt(2)),iz0= Lz/2.0 - Lz/6.0;
  for(ix=0;ix<Lx;ix++) //para cada celda
    for(iy=0;iy<Ly;iy++)
      for(iz=0;iz<Lz;iz++){
        //Compute the constants
        sigma0=sigma(ix,iy,iz); mur0=mur(ix,iy,iz); epsilonr0=epsilonr(ix,iy,iz);
        //speedCell = Ccell(epsilonr0,mur0);
        c_sqrdPrime = c2Prime(ix,iy,iz);
        speedCell = sqrt(c2(ix,iy,iz));
        B0.cargue(0,0,0);
        E0.cargue(0,0,0);

        for(int alpha=0; alpha < 3; alpha++)
          for(int beta=0; beta < 3; beta++)
            for(i=0;i<Qi;i++){
              if(i==0)
                w_i = weight_0(ix,iy,iz);
              else
                w_i = weight_i(ix,iy,iz);

              id = index(ix,iy,iz,i,alpha,beta);
              fnew[id]=f[id]=feq(E0,B0,c_sqrdPrime,w_i,i,alpha,beta);
            }
      }
}

void LatticeBoltzmann::Collision(void){
  int ix,iy,iz,i,r,id;
  double mur0,epsilonr0,sigma0,c_sqrdPrime,w_i;;
  vector3D E0,B0,J;
  for(ix=0;ix<Lx;ix++) //para cada celda
    for(iy=0;iy<Ly;iy++)
      for(iz=0;iz<Lz;iz++){
        //Compute the constants
        sigma0=sigma(ix,iy,iz); mur0=mur(ix,iy,iz); epsilonr0=epsilonr(ix,iy,iz);
        c_sqrdPrime = c2Prime(ix,iy,iz);
        //Compute the fields
        E0=E(ix,iy,iz,false); B0=B(ix,iy,iz,false);
        J = sigma0*E0/Xfactor;
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
              fnew[id]=2*feq(E0,B0,c_sqrdPrime,w_i,i,alpha,beta)-f[id] +
                T(J,B0,c_sqrdPrime,w_i,i,alpha,beta);
            }

      }
}

void LatticeBoltzmann::ImposeFields(int t){
  int ix,iy,iz,i,id;
  vector3D E0,B0,J;
  double sigma0,speedCell,c_sqrdPrime,w_i;
  //Impose fields in iz=0
  iz = 0;
  for(ix=0;ix<Lx;ix++) //para cada celda
    for(iy=0;iy<Ly;iy++){
        //Compute the constants
        sigma0=sigma(ix,iy,iz);
        //speedCell = Ccell(epsilonr0,mur0);
        c_sqrdPrime = c2Prime(ix,iy,iz);
        speedCell = sqrt(c2(ix,iy,iz));
        B0.cargue(0,(E00/speedCell/Xfactor)*sin(omega*t),0);
        E0.cargue(E00*sin(omega*t),0,0);

        for(int alpha=0; alpha < 3; alpha++)
          for(int beta=0; beta < 3; beta++)
            for(i=0;i<Qi;i++){
              if(i==0)
                w_i = weight_0(ix,iy,iz);
              else
                w_i = weight_i(ix,iy,iz);

              id = index(ix,iy,iz,i,alpha,beta);
              fnew[id]=f[id]=feq(E0,B0,c_sqrdPrime,w_i,i,alpha,beta);
            }
      }

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


void LatticeBoltzmann::Print(vector<double>& EAmplit){
  int ix=0,iy=0,iz; double sigma0,mur0,epsilonr0;
  vector3D B0,E0; double E2,B2,eps,mus;
  ofstream Myfile("data.dat");
  int iAmpl = 0;
  double delta, mu,aux;

  for(iz=Lz/4.0;iz<Lz/2.0;iz++){
    //Compute the electromagnetic constants
    sigma0=sigma(ix,iy,iz); mur0=mur(ix,iy,iz); epsilonr0=epsilonr(ix,iy,iz);
    mu = mur0*Mu0;
    aux = omega*epsilonr0*Epsilon0/sigma0;
    delta = sqrt(2/(sigma0*omega*mu)*(sqrt(1+aux*aux)+aux));
    //Compute the Fields
    E0=E(ix,iy,iz,true); B0=B(ix,iy,iz,true)*Xfactor;
    //Print
    E2=norma2(E0); B2=norma2(B0);
    eps = Epsilon0*epsilonr0;
    mus = Mu0*mur0;
    //cout<<iz<<" "<<0.5*(eps*E2+B2/mus)<<endl;
    Myfile<<double(iz)/Lz<<"\t"<<EAmplit[iAmpl]<<"\t"<<exp(-(iz-Lz/4.0)/delta)
          << "\t" << E0.x()/E00<< endl;
    iAmpl ++;
  }
  Myfile.close();
}

void LatticeBoltzmann::MeasureAmplitude(vector<double> &EAmplit){
  int ix=0,iy=0,iz; double sigma0,mur0,epsilonr0;
  vector3D B0,E0; double E2,B2,eps,mus;
  int iAmpl = 0;
  for(iz=Lz/4.0;iz<Lz/2.0;iz++){
        //Compute the electromagnetic constants
    sigma0=sigma(ix,iy,iz); mur0=mur(ix,iy,iz); epsilonr0=epsilonr(ix,iy,iz);

    //Compute the Fields
    E0=E(ix,iy,iz,true); B0=B(ix,iy,iz,true)*Xfactor;

    EAmplit[iAmpl] = abs(E0.x()/E00);
    iAmpl ++;
  }
}


int main(){
  LatticeBoltzmann SkinWave;
  int t, tmax=Lz/(2*C)*Xfactor;

  int numE = Lz/2 - Lz/4;
  vector<double> EAmplit(numE,0);
  vector<double> OldEAmplit(numE,0);

  SkinWave.Start();
  SkinWave.ImposeFields(0);

  for(t=0;t<tmax+T;t++){
    SkinWave.Collision();
    SkinWave.ImposeFields(t);
    SkinWave.Advection();
    if(t == tmax) SkinWave.MeasureAmplitude(OldEAmplit);
    if(t > tmax){
      SkinWave.MeasureAmplitude(EAmplit);
      Amplitud(EAmplit,OldEAmplit);
      OldEAmplit = EAmplit;
    }
  }

  SkinWave.Print(EAmplit);

  return 0;
}

void Amplitud(vector<double>& Amplitud,  vector<double>& oldAmplitud) {
  transform(oldAmplitud.begin(), oldAmplitud.end(), Amplitud.begin(), Amplitud.begin(), [](double a, double b) {
    double a_abs = abs(a), b_abs = abs(b);
    return (a_abs > b_abs) ? a_abs : b_abs;
  });
}
