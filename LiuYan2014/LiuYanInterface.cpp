//Dielectric Interface in 1D
#include<iostream>
#include<cmath>
#include "Vector.h"
using namespace std;

//------------------------CONSTANTS-------------------------------
const int Lx = 1;   //
const int Ly = 1;   //
const int Lz = 200; //
//-------------------
const double Tau = 0.5;
const double UTau = 1/Tau;
const double UmUTau=1-1/Tau;
//-------------------
const double Epsilon0=1, Mu0=1;
const double Sigma0=0.0125;
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
  return 1.0;
}

//--------------------- class LatticeBoltzmann ------------
class LatticeBoltzmann{
private:
  vector3D v[13]; //v[alpha]
  double f[Lx][Ly][Lz][13][3],fnew[Lx][Ly][Lz][13][3];//f[ix][iy][iz][alpha][i]
  double g[Lx][Ly][Lz][13][3],gnew[Lx][Ly][Lz][13][3];//f[ix][iy][iz][alpha][i]
public:
  LatticeBoltzmann(void);
  //Auxiliary variables
  double Ccell(double epsilonr0,double mur0){return C/sqrt(epsilonr0*mur0);};
  //Fields from direct sums
  vector3D E(int ix,int iy,int iz,bool UseNew);
  vector3D H(int ix,int iy,int iz,bool UseNew);
  //Factores auxiliares para feq y geq
  double A0(vector3D & E0, double epsilonr0, double mur0);
  double B0(vector3D & E0, double epsilonr0, double mur0);
  double A2(vector3D & E0, double epsilonr0, double mur0);
  double B2(vector3D & E0, double epsilonr0, double mur0);
  //Equilibrium Functions
  double feq(vector3D & H0,vector3D & E0,
	     double Epsilonr,double Mur, int alpha, int i);
  double feq0(vector3D & H0, int i);
  double geq(vector3D & H0,vector3D & E0,
	     double Epsilonr,double Mur, int alpha, int i);
  double geq0(vector3D & E0, int i);
  //Simulation Functions
  void Start(void);
  void Collision(void);
  void ImposeFields(int t);
  void Advection(void);
  void Print(void);
};

LatticeBoltzmann::LatticeBoltzmann(void){
  int ix,iy,iz,alpha,r,p,i,j;
  //Velocity vectors 
  v[0].cargue(0,0,0);
  v[1].cargue(1,0,0);
  v[2].cargue(-1,0,0);
  v[3].cargue(0,1,0);
  v[4].cargue(0,-1,0);
  v[5].cargue(0,0,1);
  v[6].cargue(0,0,-1);
  v[7].cargue(2,0,0);
  v[8].cargue(-2,0,0);
  v[9].cargue(0,2,0);
  v[10].cargue(0,-2,0);
  v[11].cargue(0,0,2);
  v[12].cargue(0,0,-2);
}
//-----------------MACROSCOPIC FIELDS------------------
//Fields from direct sums

vector3D LatticeBoltzmann::E(int ix,int iy,int iz,bool UseNew){
  int alpha; vector3D sum,x,y,z;
  sum.cargue(0,0,0);
  x.cargue(1,0,0), y.cargue(0,1,0), z.cargue(0,0,1);
  
  //Cargar f
  for(alpha=0;alpha<13;alpha++){
	if(UseNew) 
	    sum+=fnew[ix][iy][iz][alpha][0]*x +fnew[ix][iy][iz][alpha][1]*y+fnew[ix][iy][iz][alpha][2]*z;
	else
    	sum+=f[ix][iy][iz][alpha][0]*x +f[ix][iy][iz][alpha][1]*y+f[ix][iy][iz][alpha][2]*z;
  }
  return sum;
}
  
vector3D LatticeBoltzmann::H(int ix,int iy,int iz,bool UseNew){
  int alpha; vector3D sum,x,y,z;
  sum.cargue(0,0,0);
  x.cargue(1,0,0), y.cargue(0,1,0), z.cargue(0,0,1);

  //Cargar otras g
  for(alpha=0;alpha<13;alpha++){
	if(UseNew) 
	    sum+=gnew[ix][iy][iz][alpha][0]*x +gnew[ix][iy][iz][alpha][1]*y+gnew[ix][iy][iz][alpha][2]*z;
	else
    	sum+=g[ix][iy][iz][alpha][0]*x +g[ix][iy][iz][alpha][1]*y+g[ix][iy][iz][alpha][2]*z;
  }
  return sum;
}

//Factores auxiliares funciones de equilibrio
double LatticeBoltzmann::A0(vector3D & E0, double epsilonr0, double mur0){
  double b = 6;
  double c = Ccell(epsilonr0,mur0);
  vector3D D0 =  (epsilonr0*Epsilon0)*E0;
  double D = norma(D0);

  double aux1 = 4*D/(mur0*Mu0*epsilonr0*Epsilon0*3*b*pow(c,2));
  double aux2 = D/(pow(mur0*Mu0*epsilonr0*Epsilon0,2)*3*b*pow(c,4));
  return aux1 - aux2;
};

double LatticeBoltzmann::B0(vector3D & E0, double epsilonr0, double mur0){
  double b = 6;
  double c = Ccell(epsilonr0,mur0);
  vector3D D0 =  (epsilonr0*Epsilon0)*E0;
  double D = norma(D0);

  double aux1 = D/(pow(mur0*Mu0*epsilonr0*Epsilon0,2)*12*b*pow(c,4));
  double aux2 = D/(mur0*Mu0*epsilonr0*Epsilon0*12*b*pow(c,2));
  
  return aux1 - aux2;
};

double LatticeBoltzmann::A2(vector3D & E0, double epsilonr0, double mur0){
  double b = 6;
  double c = Ccell(epsilonr0,mur0);
  vector3D D0 =  (epsilonr0*Epsilon0)*E0;
  double D = norma(D0);

  double aux1 = 4*D/(3*b*pow(c,2));
  double aux2 = D/(mur0*Mu0*epsilonr0*Epsilon0*3*b*pow(c,4));
  
  return aux1 - aux2;
};

double LatticeBoltzmann::B2(vector3D & E0, double epsilonr0, double mur0){
  double b = 6;
  double c = Ccell(epsilonr0,mur0);
  vector3D D0 =  (epsilonr0*Epsilon0)*E0;
  double D = norma(D0);
  
  double aux1 = D/(mur0*Mu0*epsilonr0*Epsilon0*12*b*pow(c,4));
  double aux2 = D/(12*b*pow(c,2));

  return aux1 - aux2;
};

//---------------EQUILIBRIUM FUNCTIONS-------------
double LatticeBoltzmann::feq(vector3D & E0,vector3D & H0, double epsilonr0,double mur0, int alpha, int i){
  double A00 = A0(E0,epsilonr0,mur0);
  double B00 = B0(E0,epsilonr0,mur0);
  double A20 = A2(E0,epsilonr0,mur0);
  double B20 = B2(E0,epsilonr0,mur0);
  double D0 = 1-6*(A00+B00);

  vector3D EcrossV = E0^v[alpha];

  double aux;
  if(alpha == 0)
    aux = D0*H0[i];
  if(alpha <= 6)
    aux = A00*H0[i]+ (A20*EcrossV[i])/(mur0*Mu0);
  if(alpha > 6)
    aux = B00*H0[i]+ (B20*EcrossV[i])/(mur0*Mu0);
  
  return aux;
}

double LatticeBoltzmann::geq(vector3D & E0,vector3D & H0, double epsilonr0,double mur0, int alpha, int i){
  double a00 = A0(E0,epsilonr0,mur0);
  double b00 = B0(E0,epsilonr0,mur0);
  double a20 = A2(E0,epsilonr0,mur0);
  double b20 = B2(E0,epsilonr0,mur0);
  double d0 = 1-6*(a00+b00);

  vector3D HcrossV = H0^v[alpha];

  double aux;
  if(alpha == 0)
    aux = d0*E0[i];
  if(alpha <= 6)
    aux = a00*E0[i]+ (a20*HcrossV[i])/(mur0*Mu0);
  if(alpha > 6)
    aux = b00*E0[i]+ (b20*HcrossV[i])/(mur0*Mu0);
  return aux;
}

//-------------------SIMULATION FUNCTIONS ----------------------------
void LatticeBoltzmann::Start(void){
  int ix,iy,iz,alpha,i; double mur0,epsilonr0;
  vector3D E0,H0;
  double amp = 0.01, iz0 = 40;
  for(ix=0;ix<Lx;ix++) //para cada celda
    for(iy=0;iy<Ly;iy++)
      for(iz=0;iz<Lz;iz++){
	//Compute the constants
  mur0=mur(ix,iy,iz); epsilonr0=epsilonr(ix,iy,iz);
	//Impose the fields
  H0.cargue(0,B00*exp(-amp*(iz -iz0)*(iz - iz0)),0); //pulso
  E0.cargue(E00*exp(-amp*(iz -iz0)*(iz - iz0)),0,0); //pulso
	//Impose f=fnew=feq with the desired fields
	for(alpha=0;alpha<13;alpha++)
    for(i=0;i<3;i++){
      fnew[ix][iy][iz][alpha][i]=f[ix][iy][iz][alpha][i]=feq(E0,H0,epsilonr0,mur0,alpha,i);
      gnew[ix][iy][iz][alpha][i]=g[ix][iy][iz][alpha][i]=geq(E0,H0,epsilonr0,mur0,alpha,i);
    }
      }
}

void LatticeBoltzmann::Collision(void){
  int ix,iy,iz,alpha,i; double mur0,epsilonr0;
  vector3D E0,H0;
  for(ix=0;ix<Lx;ix++) //para cada celda
    for(iy=0;iy<Ly;iy++)
      for(iz=0;iz<Lz;iz++){
	//Compute the constants
  mur0=mur(ix,iy,iz); epsilonr0=epsilonr(ix,iy,iz);
	//Compute the fields
  E0=E(ix,iy,iz,false); H0=H(ix,iy,iz,false);
	//BGK evolution rule
	for(alpha=0;alpha<13;alpha++)
    for(i=0;i<3;i++){
      fnew[ix][iy][iz][alpha][i]=UmUTau*f[ix][iy][iz][alpha][i]
		  +UTau*feq(E0,H0,epsilonr0,mur0,alpha,i); 
      gnew[ix][iy][iz][alpha][i]=UmUTau*g[ix][iy][iz][alpha][i]
		  +UTau*geq(E0,H0,epsilonr0,mur0,alpha,i);
    }
      }
}
/* 
void LatticeBoltzmann::ImposeFields(int t){
  int ix,iy,iz,r,p,i,j; double sigma0,mur0,epsilonr0,prefactor0;
   double rhoc0; vector3D D0,B0,E0,H0,Jprima0,Eprima0;
  double T=25.0,omega=2*M_PI/T;//25<T<50
  double alpha = 0.01, z0 = 40;
  iz=0; //On the whole incident plane
  for(ix=0;ix<Lx;ix++) //for each cell
    for(iy=0;iy<Ly;iy++){
      //Compute the constants
      sigma0=sigma(ix,iy,iz); mur0=mur(ix,iy,iz); epsilonr0=epsilonr(ix,iy,iz);
      prefactor0=prefactor(epsilonr0,sigma0);
      //Compute the fields
      //Primary fields (i.e. those from the sums)
      rhoc0=rhoc(ix,iy,iz,false); 
      D0.cargue(E00*sin(omega*t)*epsilonr0,0,0);
      B0.cargue(0,B00*sin(omega*t),0);
      //Secundary fields (i.e. computed from the primary fields)
      E0=E(D0,epsilonr0); H0=H(B0,mur0);
      Jprima0=Jprima(E0,prefactor0); Eprima0=Eprima(E0,Jprima0,epsilonr0); 
      //Impose fnew=feq with the desired fields
      for(r=0;r<2;r++)
	for(p=0;p<3;p++)
	  for(i=0;i<4;i++)
	    for(j=0;j<2;j++){
	      fnew[ix][iy][iz][r][p][i][j]=feq(Jprima0,Eprima0,B0,epsilonr0,mur0,r,p,i,j);
	      f0new[ix][iy][iz]=feq0(rhoc0);
	    }
    }
}

*/
void LatticeBoltzmann::Advection(void){
  int ix,iy,iz,alpha,i,ixnew,iynew,iznew;
  for(ix=0;ix<Lx;ix++) 
    for(iy=0;iy<Ly;iy++)
      for(iz=0;iz<Lz;iz++){//for each cell
        for(alpha=0;alpha<13;alpha++)
            ixnew=(int(ix+v[alpha][0]+Lx))%Lx; iynew=(int(iy+v[alpha][1]+Ly))%Ly; iznew=(int(iz+v[alpha][2]+Lz))%Lz;
        for(i=0;i<3;i++){
          f[ixnew][iynew][iznew][alpha][i]=fnew[ix][iy][iz][alpha][i];
          g[ixnew][iynew][iznew][alpha][i]=gnew[ix][iy][iz][alpha][i];
        }
      }
}

void LatticeBoltzmann::Print(void){
  int ix=0,iy=0,iz,alpha,i; double mur0,epsilonr0;
  vector3D E0,H0; double E2,H2;
  for(iz=0;iz<Lz;iz++){
    //Compute the electromagnetic constante
    mur0=mur(ix,iy,iz); 
    epsilonr0=epsilonr(ix,iy,iz);
    //Compute the Fields
    E0=E(ix,iy,iz,true); H0=H(ix,iy,iz,true);
    //Print
    cout<<iz<<" "<<E0.x()/E00<<endl;
  }
}


int main(){
  LatticeBoltzmann OndaInterface;
  int t, tmax=10;
  
  OndaInterface.Start();
  //OndaSkin.ImposeFields(0);
  
  for(t=0;t<tmax;t++){
    OndaInterface.Collision();
    //OndaSkin.ImposeFields(t);
    OndaInterface.Advection();
  }
  
  OndaInterface.Print();
  
  return 0;
}
