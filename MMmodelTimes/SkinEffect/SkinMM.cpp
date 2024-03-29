#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <algorithm>
#include <vector>
#include "Vector.h"
#include "Parameters.h"

using namespace std;

void Amplitud(vector<double>& Amplitud, vector<double>& oldAmplitud);

//------------------------CONSTANTS-------------------------------
const int Qr = 2, Qp = 3, Qi = 4, Qj = 2;
//-------------------
const double Tau = 0.5;
const double UTau = 1/Tau;
const double UmUTau=1-1/Tau;

//--------------------- class LatticeBoltzmann ------------
class LatticeBoltzmann{
  private:
    int V[3][3][4], V0[3]; /*V[xyz][p][i]*/  vector3D v[3][4],v0; //v[p][i]
    vector3D e[3][4][2], e0; //e[p][i][j]
    vector3D b[3][4][2], b0; //b[p][i][j]
    double *f=nullptr,*fnew=nullptr;//f[ix][iy][iz][r][p][i][j]
    double *f0 = nullptr,*f0new = nullptr;//f0[ix][iy][iz] (r=0)

    int Lx,Ly,Lz;
    int scale = 1;
    //------------
    double Epsilon0, Mu0;
    double C;
    double Sigma0;
    //-----------------
    double E00,B00;
    //-----------------
    double T, omega;

  public:
    LatticeBoltzmann(Parameter P);
    ~LatticeBoltzmann();
    int index(int ix,int iy,int iz,int r,int p,int i,int j);
    int index0(int ix,int iy,int iz);
    //Auxiliary variables
    double Ccell(double epsilonr0,double mur0){return C/sqrt(epsilonr0*mur0);};
    double prefactor(double epsilonr0,double sigma0){return sigma0/(1+(sigma0*Mu0)/(4*epsilonr0));};
    //Constants for the medium
    double mur(int ix,int iy,int iz){return 1.0;};
    double epsilonr(int ix,int iy,int iz){return 1.0;};
    double sigma(int ix,int iy,int iz){return Sigma0*(1+tanh(iz-Lz/4.0));};
    double getPeriod(){return T;}
    //Fields from direct sums
    double rhoc(int ix,int iy,int iz,bool UseNew);
    vector3D D(int ix,int iy,int iz,bool UseNew);
    vector3D B(int ix,int iy,int iz,bool UseNew);
    //Fields deduced from the first ones through electromagnetic constants
    vector3D E(vector3D & D0,double Epsilonr);
    vector3D H(vector3D & B0,double Mur);
    //Forced (Actual) Macroscopic Fields (inline)
    vector3D Jprima(vector3D & E0,double prefactor0){return E0*prefactor0;};
    vector3D Eprima(vector3D & E0,vector3D & Jprima0,double epsilonr0){return E0-Jprima0*(Mu0/(4*epsilonr0));};
    vector3D Dprima(vector3D & Eprima0,double epsilonr0){return Eprima0/epsilonr0;};
    //Equilibrium Functions
    double feq(vector3D & Jprima0,vector3D & Eprima0,vector3D & B0,
               double Epsilonr,double Mur,
               int r,int p,int i,int j);
    double feq0(double rhoc0);
    //Simulation Functions
    void Start(void);
    void Collision(void);
    void ImposeFields(int t);
    void Advection(void);
    void Print(vector<double>& EAmplit, int iter);
    void MeasureAmplitude(vector<double> &EAmplit);
};

LatticeBoltzmann::LatticeBoltzmann(Parameter P){
  Lx = P.Lx;
  Ly = P.Ly;
  Lz = P.Lz;
  //------------Characteristic period-----------------------
  T = 17.68/1000*Lz/C*sqrt(epsilonr(0,0,0)*mur(0,0,0)); 
  omega = 2*M_PI/T;
  //-------------------
  Epsilon0 = P.Epsilon0; Mu0 = P.Mu0;
  Sigma0 = 0.1*Epsilon0/T;
  C = P.C;
  E00 = P.E00; B00 = P.B00;
  //-------------------
  int ix,iy,iz,alpha,r,p,i,j;
  //Velocity vectors V[p][i]=V^p_i (in components)
  V0[0]=V0[1]=V0[2]=0;

  V[0][0][0]=V[0][1][0]=V[1][2][0]=1;
  V[1][0][0]=V[2][1][0]=V[2][2][0]=1;
  V[2][0][0]=V[1][1][0]=V[0][2][0]=0;

  V[0][0][1]=V[0][1][1]=V[1][2][1]=-1;
  V[1][0][1]=V[2][1][1]=V[2][2][1]=1;
  V[2][0][1]=V[1][1][1]=V[0][2][1]=0;

  V[0][0][2]=V[0][1][2]=V[1][2][2]=-1;
  V[1][0][2]=V[2][1][2]=V[2][2][2]=-1;
  V[2][0][2]=V[1][1][2]=V[0][2][2]=0;

  V[0][0][3]=V[0][1][3]=V[1][2][3]=1;
  V[1][0][3]=V[2][1][3]=V[2][2][3]=-1;
  V[2][0][3]=V[1][1][3]=V[0][2][3]=0;
  //Velocity vectors V[p][i]=V^p_i (as vectors)
  v0.cargue(V0[0],V0[1],V0[2]); //cargue= load (in Spanish)
  for(p=0;p<3;p++)
    for(i=0;i<4;i++){
      v[p][i].cargue(V[0][p][i],V[1][p][i],V[2][p][i]);
    }
  //Electric vectors e[p][i][j]=e^p_{ij}
  e0.cargue(0,0,0);
  for(p=0;p<3;p++)
    for(i=0;i<4;i++){
      e[p][i][0]=v[p][(i+1)%4]*0.5;
      e[p][i][1]=v[p][(i+3)%4]*0.5;
    }
  //Magnetic vectors b[p][i][j]=b^p_{ij}=v^p_i x e^p_{ij}
  b0.cargue(0,0,0);  
  for(p=0;p<3;p++)
    for(i=0;i<4;i++)
      for(j=0;j<2;j++)
        b[p][i][j]=(v[p][i]^e[p][i][j]);

  f = new double[Lx*Ly*Lz*Qr*Qp*Qi*Qj]; fnew = new double[Lx*Ly*Lz*Qr*Qp*Qi*Qj];
  f0 = new double[Lx*Ly*Lz]; f0new=new double[Lx*Ly*Lz];
}

LatticeBoltzmann::~LatticeBoltzmann(void){
  delete[] f;  delete[] fnew;
  delete[] f0; delete[] f0new;
}

int LatticeBoltzmann::index(int ix,int iy,int iz,int r,int p,int i,int j){
  return (iz*Lx*Ly+iy*Lx+ix)*Qr*Qp*Qi*Qj + (r*Qp*Qi*Qj + p*Qi*Qj + i*Qj +j);
}

int LatticeBoltzmann::index0(int ix,int iy,int iz){
  return (iz*Lx*Ly+iy*Lx+ix);
}

//-----------------MACROSCOPIC FIELDS------------------
//Fields from direct sums
double LatticeBoltzmann::rhoc(int ix,int iy,int iz,bool UseNew){
  int p,i,j; double sum;
  int id0,id;
  id0 = index0(ix,iy,iz);
  //Start for the distribution for the central (zero) vector
  if(UseNew) 
    sum=f0new[id0];
  else 
    sum=f0[id0];
  //Add all the others
  for(p=0;p<2;p++)
    for(i=0;i<4;i++)
      for(j=0;j<2;j++){
        id = index(ix,iy,iz,0,p,i,j);
        if(UseNew)
          sum+=fnew[id];
        else
          sum+=f[id];
      }
  return sum;
}
vector3D LatticeBoltzmann::D(int ix,int iy,int iz,bool UseNew){
  int p,i,j,id; vector3D sum; sum.cargue(0,0,0);
  for(p=0;p<3;p++)
    for(i=0;i<4;i++)
      for(j=0;j<2;j++){
        id = index(ix,iy,iz,0,p,i,j);
        if(UseNew)
          sum+=e[p][i][j]*fnew[id];
        else
          sum+=e[p][i][j]*f[id];
      }
  return sum;
}
vector3D LatticeBoltzmann::B(int ix,int iy,int iz,bool UseNew){
  int p,i,j,id; vector3D sum; sum.cargue(0,0,0);
  for(p=0;p<3;p++)
    for(i=0;i<4;i++)
      for(j=0;j<2;j++){
        id = index(ix,iy,iz,1,p,i,j);
        if(UseNew)
          sum+=b[p][i][j]*fnew[id];
        else
          sum+=b[p][i][j]*f[id];
      }
  return sum;
}
//Fields deduced from the first ones through electromagnetic constants
vector3D LatticeBoltzmann::E(vector3D & D0,double Epsilonr){
  return D0*(1.0/Epsilonr);
}
vector3D LatticeBoltzmann::H(vector3D & B0,double Mur){
  return B0*(1.0/Mur);
}
//---------------EQUILIBRIUM FUNCTIONS-------------
double LatticeBoltzmann::feq(vector3D & Jprima0,vector3D & Eprima0,vector3D & B0,
                             double epsilonr0,double mur0,
                             int r,int p,int i,int j){
  double VdotJp=(v[p][i]*Jprima0),Epdote=(e[p][i][j]*Eprima0),Bdotb=(b[p][i][j]*B0),aux;
  if(r==0)
    aux=0.25*(0.25*VdotJp+epsilonr0*Epdote+0.5/mur0*Bdotb);
  if(r==1)
    aux=0.25*(0.25*VdotJp+Epdote+0.5*Bdotb);
  return aux;
}
double LatticeBoltzmann::feq0(double rhoc0){
  return rhoc0;
}

//-------------------SIMULATION FUNCTIONS ----------------------------
void LatticeBoltzmann::Start(void){
  int ix,iy,iz,r,p,i,j; double sigma0,mur0,epsilonr0,prefactor0;
  int id0,id;
  double rhoc0; vector3D D0,B0,E0,H0,Jprima0,Eprima0;
  for(ix=0;ix<Lx;ix++) //para cada celda
    for(iy=0;iy<Ly;iy++)
      for(iz=0;iz<Lz;iz++){
        //Compute the constants
        sigma0=sigma(ix,iy,iz); mur0=mur(ix,iy,iz); epsilonr0=epsilonr(ix,iy,iz);
        prefactor0=prefactor(epsilonr0,sigma0);
        //Impose the fields
        rhoc0=0; D0.cargue(0,0,0); B0.cargue(0,0,0);
        E0=E(D0,epsilonr0); H0=H(B0,mur0);
        Jprima0=Jprima(E0,prefactor0); Eprima0=Eprima(E0,Jprima0,epsilonr0);
        //Impose f=fnew=feq with the desired fields
        id0 = index0(ix,iy,iz);
        f0new[id0]=f0[id0]=feq0(rhoc0);
        for(r=0;r<2;r++)
          for(p=0;p<3;p++)
            for(i=0;i<4;i++)
              for(j=0;j<2;j++){
                id = index(ix,iy,iz,r,p,i,j);
                fnew[id]=f[id]=feq(Jprima0,Eprima0,B0,epsilonr0,mur0,r,p,i,j);
              }
      }
}

void LatticeBoltzmann::Collision(void){
  int ix,iy,iz,r,p,i,j; double sigma0,mur0,epsilonr0,prefactor0;
  int id0,id;
  double rhoc0; vector3D D0,B0,E0,H0,Jprima0,Eprima0;
  for(ix=0;ix<Lx;ix++) //para cada celda
    for(iy=0;iy<Ly;iy++)
      for(iz=0;iz<Lz;iz++){
        //Compute the constants
        sigma0=sigma(ix,iy,iz); mur0=mur(ix,iy,iz); epsilonr0=epsilonr(ix,iy,iz);
        prefactor0=prefactor(epsilonr0,sigma0);
        //Compute the fields
        rhoc0=rhoc(ix,iy,iz,false); D0=D(ix,iy,iz,false); B0=B(ix,iy,iz,false);
        E0=E(D0,epsilonr0); H0=H(B0,mur0);
        Jprima0=Jprima(E0,prefactor0); Eprima0=Eprima(E0,Jprima0,epsilonr0);
        //BGK evolution rule
        id0 = index0(ix,iy,iz);
        f0new[id0]=UmUTau*f0[id0]+UTau*feq0(rhoc0);
        for(r=0;r<2;r++)
          for(p=0;p<3;p++)
            for(i=0;i<4;i++)
              for(j=0;j<2;j++){
                id = index(ix,iy,iz,r,p,i,j);
                fnew[id]=UmUTau*f[id]+UTau*feq(Jprima0,Eprima0,B0,epsilonr0,mur0,r,p,i,j);
              }
      }
}
void LatticeBoltzmann::ImposeFields(int t){
  int ix,iy,iz,r,p,i,j; double sigma0,mur0,epsilonr0,prefactor0,speedCell;
  int id0,id;
  double rhoc0; vector3D D0,B0,E0,H0,Jprima0,Eprima0;
  //25<T<50
  iz=0; //On the whole incident plane
  for(ix=0;ix<Lx;ix++) //for each cell
    for(iy=0;iy<Ly;iy++){
      //Compute the constants
      sigma0=sigma(ix,iy,iz); mur0=mur(ix,iy,iz); epsilonr0=epsilonr(ix,iy,iz);
      prefactor0=prefactor(epsilonr0,sigma0);
      speedCell = Ccell(epsilonr0,mur0);
      //Compute the fields
      //Primary fields (i.e. those from the sums)
      rhoc0=rhoc(ix,iy,iz,false); 
      D0.cargue(E00*sin(omega*t)*epsilonr0,0,0);
      B0.cargue(0,(E00/speedCell)*sin(omega*t),0);
      //Secundary fields (i.e. computed from the primary fields)
      E0=E(D0,epsilonr0); H0=H(B0,mur0);
      Jprima0=Jprima(E0,prefactor0); Eprima0=Eprima(E0,Jprima0,epsilonr0); 
      //Impose fnew=feq with the desired fields
      for(r=0;r<2;r++)
        for(p=0;p<3;p++)
          for(i=0;i<4;i++)
            for(j=0;j<2;j++){
              id0 = index0(ix,iy,iz); id = index(ix,iy,iz,r,p,i,j);
              fnew[id]=feq(Jprima0,Eprima0,B0,epsilonr0,mur0,r,p,i,j);
              f0new[id0]=feq0(rhoc0);
            }
    }
  iz=Lz-1; //Free condition
  for(ix=0;ix<Lx;ix++) //for each cell
    for(iy=0;iy<Ly;iy++){
      //Compute the constants
      sigma0=sigma(ix,iy,iz); mur0=mur(ix,iy,iz); epsilonr0=epsilonr(ix,iy,iz);
      prefactor0=prefactor(epsilonr0,sigma0);
      speedCell = Ccell(epsilonr0,mur0);
      //Compute the fields
      //Primary fields (i.e. those from the sums)
      rhoc0=rhoc(ix,iy,iz,false);
      D0.cargue(0,0,0);
      B0.cargue(0,0,0);
      //Secundary fields (i.e. computed from the primary fields)
      E0=E(D0,epsilonr0); H0=H(B0,mur0);
      Jprima0=Jprima(E0,prefactor0); Eprima0=Eprima(E0,Jprima0,epsilonr0);
      //Impose fnew=feq with the desired fields
      for(r=0;r<2;r++)
        for(p=0;p<3;p++)
          for(i=0;i<4;i++)
            for(j=0;j<2;j++){
              id0 = index0(ix,iy,iz); id = index(ix,iy,iz,r,p,i,j);
              fnew[id]=feq(Jprima0,Eprima0,B0,epsilonr0,mur0,r,p,i,j);
              f0new[id0]=feq0(rhoc0);
            }
    }
}

void LatticeBoltzmann::Advection(void){
  int ix,iy,iz,r,p,i,j,ixnew,iynew,iznew;
  int id0,id,idnew;
  for(ix=0;ix<Lx;ix++) 
    for(iy=0;iy<Ly;iy++)
      for(iz=0;iz<Lz;iz++){//for each cell
        for(r=0;r<2;r++)
          for(p=0;p<3;p++)
            for(i=0;i<4;i++)
              for(j=0;j<2;j++){
                ixnew=(ix+V[0][p][i]+Lx)%Lx; iynew=(iy+V[1][p][i]+Ly)%Ly; iznew=(iz+V[2][p][i]+Lz)%Lz;
                id0=index0(ix,iy,iz);
                id = index(ix,iy,iz,r,p,i,j);idnew=index(ixnew,iynew,iznew,r,p,i,j);
                f[idnew]=fnew[id];
                f0new[id0]=f0[id0];
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
  int ix=0,iy=0,iz,r,p,i,j; double sigma0,mur0,epsilonr0,prefactor0;
  double rhoc0; vector3D D0,B0,E0,H0,Jprima0,Eprima0; double E2,B2;
  double iAmpl = 0;
  for(iz=0;iz<Lz;iz++){
    //Compute the electromagnetic constants
    sigma0=sigma(ix,iy,iz); mur0=mur(ix,iy,iz); epsilonr0=epsilonr(ix,iy,iz);
    prefactor0=prefactor(epsilonr0,sigma0);
    //Compute the Fields
    rhoc0=rhoc(ix,iy,iz,true); D0=D(ix,iy,iz,true); B0=B(ix,iy,iz,true);
    E0=E(D0,epsilonr0); H0=H(B0,mur0);
    Jprima0=Jprima(E0,prefactor0); Eprima0=Eprima(E0,Jprima0,epsilonr0);
    //Print
    EAmplit[iAmpl] = abs(Eprima0.x()/E00);
    iAmpl ++;
  }
}


int main(){
  // Data container to initialize the LB
  Parameter P;
  P.Update(1, 1, 100, 1, 1, 2, 0.0125/3.0, 1);
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
    Lz = 1000+100*(i-1);
    P.Lz = Lz;
    P.scale = i;

    // Measure the time several times to average the results
    for (int iteration=0; iteration<10; iteration++){
      // Initialize LB for this iteration
      LatticeBoltzmann SkinWave(P);
      SkinWave.Start();
      SkinWave.ImposeFields(0);

      // Max simulation time
      period = SkinWave.getPeriod();
      tmax=Lz/C;

      int numE = Lz;
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
