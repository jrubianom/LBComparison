#include <fstream>
#include <iostream>
#include <cmath>
#include <string>
#include <time.h>
#include "fstream"
#include "Vector.h"
#include "Parameters.h"
#include <algorithm>
#include <vector>

using namespace std;

//------------------------CONSTANTS-------------------------------
const int Qr = 2, Qi = 6;

//---------Auxiliar funcion
void Amplitud(vector<double>& Amplitud, vector<double>& oldAmplitud);

double levi(int i,int j,int k){
  return ((i-j)*(j-k)*(k-i))/2.0;
}


//--------------------- class LatticeBoltzmann ------------
class LatticeBoltzmann{
  private:
    int V[Qi][3], V0[3];   vector3D v[Qi];
    vector3D *f=nullptr,*fnew=nullptr;//f[ix][iy][iz][r][i]

    int Lx,Ly,Lz;
    int scale = 1;
    //------------
    double Epsilon0=3, Mu0=3;
    double C=1.0/sqrt(Epsilon0*Mu0);
    double Sigma0=0.0125/3.0;
    //-----------------
    double E00=1,B00=E00/C;
    //-----------------
    double Period, omega;
  public:
    LatticeBoltzmann(Parameter P);
    ~LatticeBoltzmann();
    double Ccell(double epsilonr0,double mur0){return C/sqrt(epsilonr0*mur0);};
    int index(int ix,int iy,int iz,int r,int i);
    int index0(int ix,int iy,int iz);
    double getPeriod(){return Period;}
    //Fields from direct sums
    vector3D D(int ix,int iy,int iz,bool UseNew);
    vector3D B(int ix,int iy,int iz,bool UseNew);
    //Fields deduced from the first ones through electromagnetic constants
    vector3D E(vector3D & D0,double Epsilonr);
    vector3D H(vector3D & B0,double Mur);
    //Force terms
    vector3D T(vector3D JPrime,double Ccell0, int r,int i);
    //Equilibrium Functions
    vector3D feq(vector3D & D,vector3D & B,int r,int i,
                 double epsr,double mur0);
    double mur(int ix,int iy,int iz){return 1.0;}
    double epsilonr(int ix,int iy, int iz){return 1.0;}
    double sigma(int ix,int iy,int iz){return Sigma0*(1+tanh(iz-Lz/4.0));}
    //Simulation Functions
    void Start(void);
    void Collision(void);
    void ImposeFields(int t);
    void Advection(void);
    void Print(vector<double>& EAmplit, int iter);
    void MeasureAmplitude(vector<double> &EAmplit);
};

LatticeBoltzmann::LatticeBoltzmann(Parameter P){
  scale = P.scale;
  Lx = P.Lx;
  Ly = P.Ly;
  Lz = P.Lz;
  //------------Characteristic period-----------------------
  Period = 17.68/100*Lz/C*sqrt(epsilonr(0,0,0)*mur(0,0,0)); 
  omega = 2*M_PI/Period;
  //-------------------
  Epsilon0 = P.Epsilon0; Mu0 = P.Mu0;
  Sigma0 = 0.1*Epsilon0/Period;
  C = P.C;
  E00 = P.E00; B00 = P.B00;
  //-------------------
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
//---------------Force Terms----------------------
vector3D LatticeBoltzmann::T(vector3D JPrime,double Ccell0, int r,int i){
  vector3D aux;
  aux.cargue(0,0,0);
  if(r == 0)
    aux = v[i]*(-1/2.0*JPrime*v[i]);
    //aux = v[i]*(-1/(6.0*Ccell0*Ccell0)*JPrime*v[i]);
  return aux;
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
        D0 = epsilonr0*Epsilon0*E0;

        for(r=0;r<Qr;r++)
          for(i=0;i<Qi;i++){
            id = index(ix,iy,iz,r,i);
            fnew[id]=f[id]=feq(D0,B0,r,i,epsilonr0,mur0);
          }
      }
}

void LatticeBoltzmann::Collision(void){
  int ix,iy,iz,i,r,id;
  double mur0,epsilonr0,eps,sigma0,Ccell0;
  vector3D D0,B0,J0;
  for(ix=0;ix<Lx;ix++) //para cada celda
    for(iy=0;iy<Ly;iy++)
      for(iz=0;iz<Lz;iz++){
        //Compute the constants
        sigma0=sigma(ix,iy,iz); mur0=mur(ix,iy,iz); epsilonr0=epsilonr(ix,iy,iz);
        eps = epsilonr0*Epsilon0;
        Ccell0 = Ccell(epsilonr0,mur0);
        //Compute the fields
        D0=D(ix,iy,iz,false);
        B0=B(ix,iy,iz,false);
        //J0 = sigma0*D0/pow(epsilonr0*Epsilon0,2)/mur0/1.0;
        J0 = sigma0*D0/pow(eps,1);
        //E0 = E(D0,epsilonr0); H0 = H(B0,mur0);
        //BGK evolution rule
        for(r = 0; r < Qr; r++)
          for(i=0; i < Qi; i++){
            id = index(ix,iy,iz,r,i);
            fnew[id]=2*feq(D0,B0,r,i,epsilonr0,mur0)-f[id]
              + T(J0,Ccell0,r,i);
              //+(r-1)*v[i]*(1/(6.0*C*C)*sigma(ix,iy,iz)*D0/pow(epsilonr0*Epsilon0,2))*v[i];
          }

      }
}

void LatticeBoltzmann::ImposeFields(int t){
  int ix,iy,iz,i,r,id;
  //Impose file, plane wave to the right in z=0
  iz = 0;
  vector3D E0,B0,D0,H0,Jprima0,P0;
  double sigma0,mur0,epsilonr0,speedCell;
  for(ix=0;ix<Lx;ix++) //para cada celda
    for(iy=0;iy<Ly;iy++){
      //Compute the constants
      sigma0=sigma(ix,iy,iz); mur0=mur(ix,iy,iz); epsilonr0=epsilonr(ix,iy,iz);
      speedCell = Ccell(epsilonr0,mur0);
      //Impose the fields
      //rhoc0=0; Jprima0.cargue(0,0,0);
      B0.cargue(0,(E00/speedCell)*sin(omega*t),0);
      E0.cargue(E00*sin(omega*t),0,0);
      D0 = epsilonr0*Epsilon0*E0;
      //H0.cargue(0,sin(iz*2*M_PI/Lz),0);
      //E0.cargue(sin(iz*2*M_PI/Lz),0,0);
      //Impose f=fnew=feq with the desired fields
      for(r=0; r < Qr; r++)
        for(i=0;i < Qi; i++){
          id = index(ix,iy,iz,r,i);
          fnew[id]=f[id]=feq(D0,B0,r,i,epsilonr0,mur0);
        }
    }
  //Impose field of free conditions in iz = -1
  iz = Lz-1;
  for(ix=0;ix<Lx;ix++) //para cada celda
    for(iy=0;iy<Ly;iy++){
      //Compute the constants
      sigma0=sigma(ix,iy,iz); mur0=mur(ix,iy,iz); epsilonr0=epsilonr(ix,iy,iz);
      speedCell = Ccell(epsilonr0,mur0);
      //Impose the fields
      //rhoc0=0; Jprima0.cargue(0,0,0);
      B0.cargue(0,0,0);
      E0.cargue(0,0,0);
      D0 = epsilonr0*Epsilon0*E0;
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

void LatticeBoltzmann::Print(vector<double>& EAmplit, int iter){
  int ix=0,iy=0,iz; double sigma0,mur0,epsilonr0;
  vector3D D0,B0,E0; double E2,B2,eps,mus;

  string filename = "data_" + to_string(iter) + ".dat";
  ofstream Myfile(filename);

  int iAmpl = 0;
  double delta, mu,aux;

  for(iz=0;iz<Lz;iz++){
  //for(iz=Lz/4.0;iz<Lz/2.0;iz++){
    //Compute the electromagnetic constants
    sigma0=sigma(ix,iy,iz); mur0=mur(ix,iy,iz); epsilonr0=epsilonr(ix,iy,iz);
    mu = mur0*Mu0;
    aux = omega*epsilonr0*Epsilon0/sigma0;
    delta = sqrt(2/(sigma0*omega*mu)*(sqrt(1+aux*aux)+aux));
//Compute the Fields
    D0=D(ix,iy,iz,true); B0=B(ix,iy,iz,true);
    E0=E(D0,epsilonr0);
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
  vector3D D0,B0,E0; double E2,B2,eps,mus;
  int iAmpl = 0;
  for(iz=0;iz<Lz;iz++){
  //for(iz=Lz/4.0;iz<Lz/2.0;iz++){
        //Compute the electromagnetic constants
    sigma0=sigma(ix,iy,iz); mur0=mur(ix,iy,iz); epsilonr0=epsilonr(ix,iy,iz);

    //Compute the Fields
    D0=D(ix,iy,iz,true); B0=B(ix,iy,iz,true);
    E0=E(D0,epsilonr0);

    EAmplit[iAmpl] = abs(E0.x()/E00);
    iAmpl ++;
  }
}


int main(){
  // Data container to initialize the LB
  Parameter P;
  P.Update(1, 1, 100, 1, 3, 3, 0.0125/3.0, 1);
  // Simulation variables
  int Lz;
  double C = P.C;
  int t,tmax,period;

  // Output file for time data
  ofstream MyTime("times.dat");
  MyTime << "Refinement\tTime\n";

  // Time measurement
  struct timespec begin, end;

  for(int i=1; i<10; i++){
    // Refine the mesh
    Lz = 100*i;
    P.Lz = Lz;
    P.scale = i;

    // Measure the time several times to average the results
    for (int iteration=0; iteration<10; iteration++){
      // Initialize LB for this iteration
      LatticeBoltzmann SkinWave(P);
      SkinWave.Start();
      SkinWave.ImposeFields(0);

      // Max simulation time
      tmax=Lz/C;
      period = SkinWave.getPeriod();

      int numE = Lz;//Lz/2 - Lz/4;
      vector<double> EAmplit(numE,0);
      vector<double> OldEAmplit(numE,0);
      
      // Time measurement
      struct timespec begin, end;

      clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &begin);
      for(t=0; t<tmax+period; t++){
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
      clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
      long seconds = end.tv_sec - begin.tv_sec;
      long nanoseconds = end.tv_nsec - begin.tv_nsec;
      double ellapsed = seconds + nanoseconds*1e-9;

      // Print fields
      SkinWave.Print(EAmplit, i);
      // Print time
      MyTime << i << "\t" << ellapsed << "\n";
  }
    }
  return 0;
}

void Amplitud(vector<double>& Amplitud,  vector<double>& oldAmplitud) {
  transform(oldAmplitud.begin(), oldAmplitud.end(), Amplitud.begin(), Amplitud.begin(), [](double a, double b) {
    double a_abs = abs(a), b_abs = abs(b);
    return (a_abs > b_abs) ? a_abs : b_abs;
  });
}
