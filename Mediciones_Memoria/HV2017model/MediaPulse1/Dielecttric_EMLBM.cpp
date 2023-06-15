#include <fstream>
#include <iostream>
#include <cmath>
#include <string>
#include "fstream"
#include "Vector.h"
#include "Parameters.h"
#include <sys/resource.h>

using namespace std;

//------------------------CONSTANTS-------------------------------

const int Qr = 2, Qi = 6;

//------------------Electromagnetic Constants for the Media------------------------------


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
    double Sigma0=0.0;
    //-----------------
    double E00=0.001,B00=E00/C;
    double epsr1 = 1, epsr2 = 2.5;

  public:
    LatticeBoltzmann(Parameter P, ofstream &fileMem, int Ref);
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

    double mur(int ix,int iy,int iz){return 1.0;}

    double epsilonr(int ix,int iy, int iz){
      double epsr = epsr1;
      if(iz > Lz/2)
        epsr = epsr2;
      //epsr = (epsr1 + epsr2)/2.0 + (epsr2 - epsr1)/2*tanh(iz-Lz/2.0);
      return epsr;
    }
    double sigma(int ix,int iy,int iz){return 0.0;}
    //Simulation Functions
    void Start(void);
    void Collision(void);
    void ImposeFields(int t);
    void Advection(void);
    void Print(ofstream &fileAmpl, ofstream &fileError,int i);
    void PrintFrame(string CurrentFrame);
    //Memory Usage
    long get_mem_usage(){
      struct rusage my_usage;
      getrusage(RUSAGE_SELF, &my_usage);
      return my_usage.ru_maxrss;
  }
};

LatticeBoltzmann::LatticeBoltzmann(Parameter P, ofstream &fileMem, int Ref){
  scale = P.scale;
  Lx = P.Lx;
  Ly = P.Ly;
  Lz = P.Lz;
//-------------------
  Epsilon0 = P.Epsilon0; Mu0 = P.Mu0;
  Sigma0 = P.Sigma0;
  C = P.C;
  E00 = P.E00; B00 = P.B00;
  epsr1 = P.epsr1; epsr2 = P.epsr2;
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
  fileMem << Ref << "\t" << sizeof(vector3D)*2*Lx*Ly*Lz*Qr*Qi/1000 <<endl; 

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

        B0.cargue(0,(E00/speedCell)*exp(-pow(iz-iz0,2)/(2*pow(alpha0,2))),0);
        E0.cargue(E00*exp(-pow(iz-iz0,2)/(2*pow(alpha0,2))),0,0);
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
  double mur0,epsilonr0,sigma0;
  vector3D D0,B0;
  for(ix=0;ix<Lx;ix++) //para cada celda
    for(iy=0;iy<Ly;iy++)
      for(iz=0;iz<Lz;iz++){
        //Compute the constants
        sigma0=sigma(ix,iy,iz); mur0=mur(ix,iy,iz); epsilonr0=epsilonr(ix,iy,iz);
        //Compute the fields
        D0=D(ix,iy,iz,false); B0=B(ix,iy,iz,false);
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

void LatticeBoltzmann::Print(ofstream &fileAmpl, ofstream &fileError,int i){
  int ix=0,iy=0,iz; double sigma0,mur0,epsilonr0;
  vector3D D0,B0,E0; double E2,B2,eps,mus;
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
    fileAmpl<<double(iz)/Lz<<" "<<E0.x()/E00<<endl;


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
  Er = Er/E00; Et = Et/E00;
  double Ertheo,Ettheo,ratio = sqrt(epsr2/epsr1);
  Ertheo = abs((ratio -1 )/(ratio + 1));
  Ettheo = abs(2/(ratio + 1));
  fileError << i << "\t" << Er << "\t"  << Ertheo << "\t" <<  abs(100*(Ertheo - Er)/Ertheo) << "\t" <<
    Et << "\t" << Ettheo << "\t" << abs(100*(Ettheo - Et)/Ettheo) << endl;

}

void LatticeBoltzmann::PrintFrame(string CurrentFrame){
  int ix=0,iy=0,iz; double sigma0,mur0,epsilonr0;
  vector3D D0,B0,E0; double E2,B2,eps,mus;
  double Et, Er,EAmp;
  ofstream fileAmpl(CurrentFrame);
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
    fileAmpl<<double(iz)/Lz<<" "<<E0.x()/E00<<endl;
  }
  fileAmpl.close();
}

//----------- Animation functions -----------
void OpenFrame(ofstream &file,string frame){
    file<<"plot '" << frame <<"'"<<" u 1:2 w l"<<endl;
}

void CloseFrame(ofstream &file){
    file<<endl;
}

void StartAnimation(ofstream &file,double xmax,double ymax){
  file<<"set terminal gif animate"<<endl;
  file<<"set output 'Animation.gif'"<<endl;
  file<<"set xlabel sprintf(\"z/L_z position\")"<<endl;
  file<<"set ylabel sprintf(\"E/E_0 Electric field\")"<<endl;
  file<<"set g"<<endl;
  file<<"unset key"<<endl;
  file<<"set xrange[0:"<<xmax<<"]"<<endl;
  file<<"set yrange["<<-ymax*0.3<<":"<<ymax<<"]"<<endl;
}
//______________________________

int main(){
  Parameter P;
  int Lz;
  double C = P.C;
  int t,tmax;
  string AmplFileName,ErrorFileName, MemFileName;
  double xmax = 1, ymax = 1.2;
  ErrorFileName = "Errors.dat";
  MemFileName = "Memory_Pulse_HV.dat";
  ofstream MyError("Errors.dat");
  ofstream MyMem(MemFileName);
  MyError << "Refinment\tSimulated Er \tTheorical Er \tRelative Error\t"<<
    "Simulated Et \tTheorical Et\tRelative Error\n";
  MyMem << "Refinement\tMemory Used\n"; 

  string current_frame = "Animation/frame.dat";
  ofstream AnimFile("animation.gp");
  StartAnimation(AnimFile,xmax,ymax);
  for (int i = 1; i < 6; i++){
    //Refine the mesh
    Lz = 100*i;
    P.Lz = Lz;

    LatticeBoltzmann DielectricPulse(P,MyMem,i);
    DielectricPulse.Start();
    DielectricPulse.ImposeFields(0);

    tmax=Lz/(2*C);
    for(t=0;t<tmax;t++){
      DielectricPulse.Collision();
      DielectricPulse.ImposeFields(t);
      DielectricPulse.Advection();
      if(i == 5){
        DielectricPulse.PrintFrame(current_frame);
        OpenFrame(AnimFile,current_frame);
        current_frame = "Animation/frame"+to_string(t)+".dat";
      }
    }
    AmplFileName = "data"+to_string(i)+".dat";
    ofstream Myfile(AmplFileName);
    DielectricPulse.Print(Myfile,MyError,i);
    Myfile.close();
  }
  MyMem.close();
  MyError.close();
  AnimFile.close();
  return 0;
}
