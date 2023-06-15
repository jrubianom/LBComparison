#include<iostream>
#include<cmath>
#include "Vector.h"
#include <fstream>
#include <vector>
#include <algorithm>
#include "Theory.h"
#include <string>

using namespace std;

double Interpolate(double q1,double q2,double q3,double q4,
                   double q5,double q6,double q7,double q8,
                   double x,double y,double z);
double Interpolatebi(double q1,double q2,double q3,double q4,
                     double x,double y);

void SphericalCoordinates(double R,double theta,double phi,double &x,double &y, double &z);
double smooth(double a,double b,double beta,double z); //smooth the current
void compare(vector<double> &V1,vector<double> &V2,int N);
void Set(vector<vector<double>> &Ss,int Ntheta,int Nphi);

//--------------------- class LatticeBoltzmann ------------
class LatticeBoltzmann{
  private:
    int V[6][3], V0[3];   vector3D v[6];
    vector3D *f=nullptr,*fnew=nullptr;//f[ix][iy][iz][r][i]
    int Qr = 2, Qi = 6; //Numbers of total r,i
    //------//Dimensions
    int Lx,Ly,Lz;
    //------------
    double Epsilon0,Mu0;
    double C,Sigma0,Z0;
    //-----------------
    double E00,B00,J0;
    double alpha; //alpha = parameter related with width of the Gaussian
    double T,lambda,k,omega;
    Theo formulas;
    int t;

  public:
    LatticeBoltzmann(Parameter Params0);
    ~LatticeBoltzmann();
    void UpdateTime(int t0){t=t0;};
    double GetTime(){return t;};
    double mur(int ix,int iy,int iz){ return 1.0;}
    double epsilonr(int ix,int iy,int iz){return 1.0;}
    double sigma(int ix,int iy,int iz){return Sigma0;}
    double Ccell(double epsilonr0,double mur0){return C/sqrt(epsilonr0*mur0);};
    int index(int ix,int iy,int iz,int r,int i);
    int index0(int ix,int iy,int iz);
    //Fields from direct sums
    vector3D D(int ix,int iy,int iz,bool UseNew);
    vector3D B(int ix,int iy,int iz,bool UseNew);
    //Fields deduced from the first ones through electromagnetic constants
    vector3D E(vector3D & D0,double Epsilonr);
    vector3D H(vector3D & B0,double Mur);
    vector3D JprimaSource(int ix,int iy,int iz);
    //Force terms
    vector3D TSource(vector3D JPrime,double Ccell0, int r,int i);
    //Equilibrium Functions
    vector3D feq(vector3D & D,vector3D & B,int r,int i,
                 double epsr,double mur0);
    //Simulation Functions
    void Start(void);
    void Collision(void);
    void ImposeFields(int t);
    void Advection(void);

    void MacroscopicFields(int ix,int iy, int iz,
                           vector3D &E0,vector3D &H0, vector3D &B0);
    double Poynting(int ix,int iy,int iz);
    double PowerAtPoint(double x,double y, double z);
    double PowerAtPointTheta(double x,double y);
    double PowerAtPointPhi(double x,double z);
    void PowerPlanePhi(double R,int Ntheta,vector<double> &V);
    void PowerPlaneTheta(double R,int Nphi,vector<double> &V);
    double PowerAtPointPlane(double x,double y,double phi0);
    void PowerPlane(double R,int Nphi,double phi0,vector<double> &V);
    void PowerByPlanes(double R,int Ntheta,int Nphi,
                       vector<vector<double>> &Ss,vector<vector<double>> &SsCurrent,
                       vector<vector<double>> &SsPlus);
    void PrintPowerPlanes(vector<vector<double>> &Ss,int Ntheta,int Nphi,double R);
    void Print(ofstream &fileB, ofstream &fileE, ofstream &contour);
};


LatticeBoltzmann::LatticeBoltzmann(Parameter Params0){

  Lx = Params0.Lx; Ly = Params0.Ly; Lz = Params0.Lz;
  Epsilon0 = Params0.Epsilon0; Mu0 = Params0.Mu0; C=Params0.C; Sigma0 = Params0.Sigma0;
  Z0 = sqrt(Mu0/Epsilon0);
  E00=Params0.E00; B00=Params0.B00;
  J0=Params0.J0;
  alpha=Params0.alpha;
  T = Params0.T;
  omega = 2*M_PI/T; lambda = C*T; k = omega/C;
  formulas.Init(Params0);
  t=0;

  int ix,iy,iz,alpha,i;
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

vector3D LatticeBoltzmann::JprimaSource(int ix,int iy,int iz){
  vector3D Jprima0;
  double Jz;
  Jz=J0*exp(-alpha*( pow((ix-Lx/2),2) + pow((iy-Ly/2),2) + pow(iz-Lz/2,2)) ) *sin(omega*t);
  Jprima0.cargue(0,0,Jz);
  return Jprima0;
}

//---------------Force Terms----------------------
vector3D LatticeBoltzmann::TSource(vector3D JPrime,double Ccell0, int r,int i){
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
        J0 = JprimaSource(ix,iy,iz);
        //E0 = E(D0,epsilonr0); H0 = H(B0,mur0);
        //BGK evolution rule
        for(r = 0; r < Qr; r++)
          for(i=0; i < Qi; i++){
            id = index(ix,iy,iz,r,i);
            fnew[id]=2*feq(D0,B0,r,i,epsilonr0,mur0)-f[id]
              + TSource(J0,Ccell0,r,i);
              //+(r-1)*v[i]*(1/(6.0*C*C)*sigma(ix,iy,iz)*D0/pow(epsilonr0*Epsilon0,2))*v[i];
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


void LatticeBoltzmann::Print(ofstream &fileB, ofstream &fileE, ofstream &contour){
  int ix,iy=Ly/2,iz;
  vector3D B0,E0,H0; double E2,B2,Eteo,Bteo,r,rxp;

  for(ix=0;ix < Lx; ix++){
    for(iz=0;iz<Lz;iz++){
      //Compute the electromagnetic constants
      MacroscopicFields(ix,iy,iz,E0,H0,B0);
      //Print
      B2 = B0.y()/J0;
      r = sqrt((ix-Lx/2)*(ix-Lx/2) + (iz-Lz/2)*(iz-Lz/2)); 
      Bteo = formulas.B(ix-Lx/2, iz-Lz/2, t)/J0;
      contour << ix <<" " << iz << " " << B2 << " " << Bteo << endl;

      if(iz==Lz/2){
        E2 = -E0.z()/(J0*Z0);
        Eteo = formulas.E_z(abs(ix-Lx/2),t)/(J0*Z0);
        //Bteo = formulas.B_y(abs(ix-Lx/2),t)/J0;
        fileB << ix << " " << B2 << " " << Bteo <<endl;
        fileE << ix << " " << E2 << " " << Eteo << endl;
      }
    }

  }
  contour.close();
  fileB.close();
  fileE.close();
}



//Return the Macroscopic Fields of E,H and B evaluated in a lattice point
void LatticeBoltzmann::MacroscopicFields(int ix,int iy, int iz,
                                         vector3D &E0,vector3D &H0, vector3D &B0){
  double mur0,epsilonr0;
  vector3D D0;

  //Compute the electromagnetic constants
  mur0=mur(ix,iy,iz); epsilonr0=epsilonr(ix,iy,iz);
  //Compute the Fields
  D0=D(ix,iy,iz,true); B0=B(ix,iy,iz,true);
  E0=E(D0,epsilonr0); H0=H(B0,mur0);
}


void LatticeBoltzmann::PowerPlaneTheta(double R,int Nphi,vector<double> &V){
  double phi, dphi = 2*M_PI/(Nphi-1);
  double x,y,z;
  for(int i=0; i < Nphi; i++){
    phi = i*dphi;
    SphericalCoordinates(R,M_PI/2.0,phi,x,y,z);
    x += Lx/2.0; y += Ly/2.0;
    V[i] = PowerAtPointTheta(x,y);
  }
}

void LatticeBoltzmann::PowerPlanePhi(double R,int Ntheta,vector<double> &V){
  double phi, dphi = 2*M_PI/(Ntheta-1);
  double x,y,z;
  for(int i=0; i < Ntheta; i++){
    phi = i*dphi;
    SphericalCoordinates(R,M_PI/2.0,phi,x,z,y);
    x += Lx/2.0; z += Lz/2.0;
    V[i] = PowerAtPointPhi(x,z);
  }
}



//Returns the interpolated value of S* r_hat of a point x,y,z
double LatticeBoltzmann::PowerAtPoint(double x,double y, double z){
  int x1 = (int) x, y1 = (int) y, z1 = (int) z;
  double s1,s2,s3,s4,s5,s6,s7,s8;
  s1 = Poynting(x1,y1,z1);
  s2 = Poynting(x1+1,y1,z1);
  s3 = Poynting(x1+1,y1+1,z1);
  s4 = Poynting(x1,y1+1,z1);
  s5 = Poynting(x1,y1,z1+1);
  s6 = Poynting(x1+1,y1,z1+1);
  s7 = Poynting(x1+1,y1+1,z1+1);
  s8 = Poynting(x1,y1+1,z1+1);
  return Interpolate(s1,s2,s3,s4,s5,s6,s7,s8,x,y,z);
}

double LatticeBoltzmann::PowerAtPointTheta(double x,double y){
  int x1 = (int) x, y1 = (int) y;
  double s1,s2,s3,s4,s5,s6,s7,s8;
  s1 = Poynting(x1,y1,Lz/2);
  s2 = Poynting(x1+1,y1,Lz/2);
  s3 = Poynting(x1,y1+1,Lz/2);
  s4 = Poynting(x1+1,y1+1,Lz/2);
  return Interpolatebi(s1,s2,s3,s4,x,y);
}

double LatticeBoltzmann::PowerAtPointPhi(double x,double z){
  int x1 = (int) x, z1 = (int) z;
  double s1,s2,s3,s4,s5,s6,s7,s8;
  s1 = Poynting(x1,Ly/2,z1);
  s2 = Poynting(x1+1,Ly/2,z1);
  s3 = Poynting(x1,Ly/2,z1+1);
  s4 = Poynting(x1+1,Ly/2,z1+1);
  return Interpolatebi(s1,s2,s3,s4,x,z);
}


void LatticeBoltzmann::PowerPlane(double R,int Ntheta,double phi0,vector<double> &V){
  double theta, dtheta = 2*M_PI/(Ntheta-1);
  double x,y,z;
  for(int i=0; i < Ntheta/2; i++){
    theta = i*dtheta;
    SphericalCoordinates(R,theta,phi0,x,y,z);
    x += Lx/2; y += Ly/2; z+=Lz/2;
    V[i] = PowerAtPoint(x,y,z);
  }
  for(int i=0; i < Ntheta/2; i++){
    theta = M_PI-i*dtheta;
    SphericalCoordinates(R,theta,phi0+M_PI,x,y,z);
    x += Lx/2; y += Ly/2; z+=Lz/2;
    V[i+Ntheta/2] = PowerAtPoint(x,y,z);
  }
}

void LatticeBoltzmann::PowerByPlanes(double R,int Ntheta,int Nphi,
                   vector<vector<double>> &Ss,vector<vector<double>> &SsCurrent,
                   vector<vector<double>> &SsPlus){
  double dphi = 2*M_PI/(Nphi-1),phi;
  for(int i=0;i < Nphi; i++){
    phi = i*dphi;
    PowerPlane(R,Ntheta,phi,SsCurrent[i]);
    compare(Ss[i],SsCurrent[i],Ntheta);
  }
}


void LatticeBoltzmann::PrintPowerPlanes(vector<vector<double>> &Ss,int Ntheta,int Nphi,double R){

  ofstream file("Datos/PowerByPlanes.txt");
  ofstream fileE("Datos/EPlane.txt");
  ofstream fileETeo("Datos/TeoEPlane.txt");
  ofstream fileB("Datos/BPlane.txt");
  ofstream fileBTeo("Datos/TeoBPlane.txt");

  double dphi,phi,dtheta,theta;
  double Pteo;
  dtheta = 2*M_PI/(Ntheta-1); dphi = 2*M_PI/(Nphi-1);
  for(int i=0; i < Nphi; i++){
    phi = i*dphi;
    for(int j=0; j < Ntheta/2; j++){
      theta = j*dtheta;
      file << theta << " " << phi << " " << Ss[i][j]/(Z0*J0*J0)*0.5 << endl;
      if(i==0){
        Pteo = formulas.Power(theta,0,R);
        fileETeo << M_PI/2-theta  << " " << Pteo/(Z0*J0*J0) << endl;
        fileE << M_PI/2-theta << " " << Ss[i][j]/(Z0*J0*J0)*0.5 << endl;
      }
    }
  }

  for(int j=0; j < Ntheta/2; j++){
    theta = (j+Ntheta/2.0)*dtheta;
    Pteo = formulas.Power(theta,0,R);
    fileETeo << M_PI/2-theta  << " " << Pteo/(Z0*J0*J0) << endl;
    fileE << M_PI/2-theta << " " << Ss[0][j+Ntheta/2]/(Z0*J0*J0)*0.5 << endl;
  }


  file.close();
  fileETeo.close();
  fileE.close();

  ofstream file2("Datos/PowerByPlanesOrder.txt");
  dtheta = 2*M_PI/(Ntheta-1); dphi = 2*M_PI/(Nphi-1);
  for(int j=0; j < Ntheta/2.0; j++){
    for(int i=0; i < Nphi; i++){
      phi = i*dphi;
      theta = j*dtheta;
      file2 << theta << " " << phi << " " << Ss[i][j]/(Z0*J0*J0)*0.5 << endl;
      if(j == (int) Ntheta/4){
        Pteo = formulas.Power(M_PI/2,phi,R);
        fileBTeo << phi << " " << Pteo/(Z0*J0*J0) << endl;
        fileB << phi << " " << Ss[i][j]/(Z0*J0*J0)*0.5<< endl;
      }
    }

  }
  file2.close();

}




//////////////////
//////////////////
//////////////////
//////////////////
//////////////////
//////////////////
//--------------Funciones auxuliares
//given a lattice point x,y,z gives the dot product of the poynting vector S* r_hat
double LatticeBoltzmann::Poynting(int ix,int iy,int iz){
  vector3D E,H,B;
  MacroscopicFields(ix,iy,iz,E,H,B);
  vector3D S = E^H;
  vector3D r; r.cargue(ix-Lx/2.0,iy-Ly/2.0,iz-Lz/2.0);
  return (S*r)*norma(r);
}

//Makes an trilinear interpolation at point x,y,z.
//q_i is the field evaluated in the i-th point of the cube
double Interpolate(double q1,double q2,double q3,double q4,
                   double q5,double q6,double q7,double q8,
                   double x,double y,double z){
  /*_              8-------------7
    _             -|            -|
    _           -  |           - |
    _         -    |          -  |
    _       5----------------6   |
    _       |     4---------|----3
    _       |    -          |    -
    _       |   -           |   -
    _       |  -            |  -
    _       |-              | -
    _      1----------------2

  */


  int x0 = (int) x, y0 = (int) y, z0 = (int) z;
  double xd = (x-x0), yd = (y-y0), zd = (z-z0);
  double c00,c01,c10,c11;
  c00 = q1*(1-xd)+q2*xd;
  c01 = q5*(1-xd)+q6*xd;
  c10 = q4*(1-xd)+q3*xd;
  c11 = q8*(1-xd)+q7*xd;

  double c0,c1,c;
  c0 = c00*(1-yd)+c10*yd;
  c1 = c01*(1-yd)+c11*yd;
  c = c0*(1-zd)+c1*zd;

  return c;
}

//Transforms from spherical cooridnates to cartesian ones
void SphericalCoordinates(double R,double theta,double phi,double &x,double &y, double &z){
  x = R*sin(theta)*cos(phi);
  y = R*sin(theta)*sin(phi);
  z = R*cos(theta);
}


double Interpolatebi(double q1,double q2,double q3,double q4,
                     double x,double y){
  int ix = (int) x, iy = (int) y;
  double u = x-ix, v = y-iy;
  double qxy = q1*(1-u)*(1-v)+q2*u*(1-v)+q3*(1-u)*v+q4*u*v;
  return qxy;
}


double smooth(double a,double b,double beta, double z){
  return tanh(beta*(z-a)) - tanh(beta*(z-b));
}


void compare(vector<double> &V1,vector<double> &V2,int N){
  for(int i=0; i < N; i++)
    if(V1[i] < V2[i])
      V1[i] = V2[i];
}

void Set(vector<vector<double>> &Ss,int Ntheta,int Nphi){
  for(int i=0;i < Nphi; i++)
    Ss.push_back(vector<double>(Ntheta,0));
}
