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
const double Epsilon0=1.0, Mu0=1.0;
const double Sigma0=0.0125;
const double C=2.0; //Velocidad de la partícula(no necesariamente vel. de la luz)
const double EA = 1.3*C*C;
const double EB = -1.0*C*C;
const double ED = 0.5*C*C;

const double E00=1.0,B00=E00/C;

const double epsr1 = 1, epsr2 = 2.5;
//------------------Electromagnetic Constants for the Media------------------------------
double mur(int ix,int iy,int iz){
  return 1.0;
}

double epsilonr(int ix,int iy, int iz){
  //Vacio
  return 1.0;
  //Interface dielectrica
  //return (epsr1 + epsr2)/2.0 + (epsr2 - epsr1)/2*tanh(iz-Lz/2.0);
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
  double wH(int ix,int iy,int iz,bool UseNew);
  double wE(int ix,int it,int iz,bool UseNew);
  //Equilibrium Functions
  double feq(vector3D & H0,vector3D & E0, double wH0,double Epsilonr,double Mur, int alpha, int i);
  double geq(vector3D & H0,vector3D & E0, double wH0,double Epsilonr,double Mur, int alpha, int i);
  //Simulation Functions
  void Start(void);
  void Collision(void);
  void ImposeFields(int t);
  void Advection(void);
  void Print(void);
};

LatticeBoltzmann::LatticeBoltzmann(void){
  //Velocity vectors 
  v[0].cargue(0,0,0);
  v[1].cargue(C,0,0);
  v[2].cargue(-C,0,0);
  v[3].cargue(0,C,0);
  v[4].cargue(0,-C,0);
  v[5].cargue(0,0,C);
  v[6].cargue(0,0,-C);
  v[7].cargue(2*C,0,0);
  v[8].cargue(-2*C,0,0);
  v[9].cargue(0,2*C,0);
  v[10].cargue(0,-2*C,0);
  v[11].cargue(0,0,2*C);
  v[12].cargue(0,0,-2*C);
}
//-----------------MACROSCOPIC FIELDS------------------
//Fields from direct sums

vector3D LatticeBoltzmann::H(int ix,int iy,int iz,bool UseNew){
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
  
vector3D LatticeBoltzmann::E(int ix,int iy,int iz,bool UseNew){
  int alpha; vector3D sum,x,y,z;
  sum.cargue(0,0,0);
  x.cargue(1,0,0), y.cargue(0,1,0), z.cargue(0,0,1);

  //Cargar g
  for(alpha=0;alpha<13;alpha++){
	if(UseNew)
	  sum+=gnew[ix][iy][iz][alpha][0]*x +gnew[ix][iy][iz][alpha][1]*y+gnew[ix][iy][iz][alpha][2]*z;
	else
	  sum+=g[ix][iy][iz][alpha][0]*x +g[ix][iy][iz][alpha][1]*y+g[ix][iy][iz][alpha][2]*z;
  }
  return sum;
}

double LatticeBoltzmann::wE(int ix,int iy,int iz, bool UseNew){
  int alpha; double Ex, sum = 0.0;

  //Arbitrariamente se usa la componente x
  for(alpha=0;alpha<13;alpha++){
    if(UseNew){
      Ex+=gnew[ix][iy][iz][alpha][0];
      if(alpha == 0)
	sum+=gnew[ix][iy][iz][alpha][0]*ED;
      else if(alpha <=6)
	sum+=gnew[ix][iy][iz][alpha][0]*EA;
      else
	sum+=gnew[ix][iy][iz][alpha][0]*EB;
    }else{
      Ex+=g[ix][iy][iz][alpha][0];
      if(alpha == 0)
	sum+=g[ix][iy][iz][alpha][0]*ED;
      else if(alpha <=6)
	sum+=g[ix][iy][iz][alpha][0]*EA;
      else
	sum+=g[ix][iy][iz][alpha][0]*EB;
    }
  }

  return sum/Ex;
  
}

double LatticeBoltzmann::wH(int ix,int iy,int iz, bool UseNew){
  int alpha; double Hx, sum = 0.0;

  //Arbitrariamente se usa la componente y
  for(alpha=0;alpha<13;alpha++){
    if(UseNew){
      Hx+=fnew[ix][iy][iz][alpha][1];
      if(alpha == 0)
	sum+=fnew[ix][iy][iz][alpha][1]*ED;
      else if(alpha <=6)
	sum+=fnew[ix][iy][iz][alpha][1]*EA;
      else
	sum+=fnew[ix][iy][iz][alpha][1]*EB;
    }else{
      Hx+=f[ix][iy][iz][alpha][1];
      if(alpha == 0)
	sum+=f[ix][iy][iz][alpha][1]*ED;
      else if(alpha <=6)
	sum+=f[ix][iy][iz][alpha][1]*EA;
      else
	sum+=f[ix][iy][iz][alpha][1]*EB;
    }
  }
  return sum/Hx;
  
}
//---------------EQUILIBRIUM FUNCTIONS-------------
double LatticeBoltzmann::feq(vector3D & E0,vector3D & H0, double wH0, double epsilonr0,double mur0, int alpha, int i){
  double b = 6;
  double c = Ccell(epsilonr0,mur0);
  double D = 3;

  //Parameters
  double A0 = (D*(ED-EB) - mur0*Mu0*epsilonr0*Epsilon0*pow(c,2)*(ED-wH0))/(mur0*Mu0*epsilonr0*Epsilon0*b*pow(c,2)*(4*EA-3*ED-EB));

  double B0 = (D*(EA-ED) + mur0*Mu0*epsilonr0*Epsilon0*pow(c,2)*(ED-wH0))/(mur0*Mu0*epsilonr0*Epsilon0*b*pow(c,2)*(4*EA-3*ED-EB));

  double A2 = (D*(wH0 - EB))/(b*pow(c,2)*(EA-EB));

  double B2 = (D*(EA - wH0))/(4*b*pow(c,2)*(EA-EB));
 
  double D0 = 1-b*(A0+B0);

  vector3D EcrossV = E0^v[alpha];

  double f_eq;
  if(alpha == 0)
    f_eq = D0*H0[i];
  else if(alpha <= 6)
    f_eq = A0*H0[i]+ (A2*EcrossV[i])/(mur0*Mu0);
  else
    f_eq = B0*H0[i]+ (B2*EcrossV[i])/(mur0*Mu0);
  
  return f_eq;
}

double LatticeBoltzmann::geq(vector3D & E0,vector3D & H0,double wH0, double epsilonr0,double mur0, int alpha, int i){
  double b = 6;
  double c = Ccell(epsilonr0,mur0);
  double D = 3;

  //Parameters
  double a0 = (D*(ED-EB) - mur0*Mu0*epsilonr0*Epsilon0*pow(c,2)*(ED-wH0))/(mur0*Mu0*epsilonr0*Epsilon0*b*pow(c,2)*(4*EA-3*ED-EB));

  double b0 = (D*(EA-ED) + mur0*Mu0*epsilonr0*Epsilon0*pow(c,2)*(ED-wH0))/(mur0*Mu0*epsilonr0*Epsilon0*b*pow(c,2)*(4*EA-3*ED-EB));

  double a2 = (D*(wH0 - EB))/(b*pow(c,2)*(EA-EB));

  double b2 = (D*(EA - wH0))/(4*b*pow(c,2)*(EA-EB));
 
  double d0 = 1-b*(a0+b0);

  vector3D HcrossV = H0^v[alpha];

  double g_eq;
  if(alpha == 0)
    g_eq = d0*E0[i];
  else if(alpha <= 6)
    g_eq = a0*E0[i]+ (a2*HcrossV[i])/(mur0*Mu0);
  else 
    g_eq = b0*E0[i]+ (b2*HcrossV[i])/(mur0*Mu0);
  return g_eq;
}

//-------------------SIMULATION FUNCTIONS ----------------------------
void LatticeBoltzmann::Start(void){
  int ix,iy,iz,alpha,i; double mur0,epsilonr0;
  vector3D E0,H0; double wH0, H2;
  double amp = 0.01, iz0 = 40;
  for(ix=0;ix<Lx;ix++) //para cada celda
    for(iy=0;iy<Ly;iy++)
      for(iz=0;iz<Lz;iz++){
	//Compute the constants
	mur0=mur(ix,iy,iz); epsilonr0=epsilonr(ix,iy,iz);
	//Impose the fields
	H0.cargue(0,B00*exp(-amp*(iz -iz0)*(iz - iz0)),0); //pulso
	E0.cargue(E00*exp(-amp*(iz -iz0)*(iz - iz0)),0,0); //pulso
	H2 = norma2(H0); wH0 = (mur0*mur0)*Mu0*H2/2.0;
	//Impose f=fnew=feq with the desired fields
	for(alpha=0;alpha<13;alpha++)
	  for(i=0;i<3;i++){
	    fnew[ix][iy][iz][alpha][i]=f[ix][iy][iz][alpha][i]=feq(E0,H0,wH0,epsilonr0,mur0,alpha,i);
	    gnew[ix][iy][iz][alpha][i]=g[ix][iy][iz][alpha][i]=geq(E0,H0,wH0,epsilonr0,mur0,alpha,i);
	  }
      }
}

void LatticeBoltzmann::Collision(void){
  int ix,iy,iz,alpha,i; double mur0,epsilonr0;
  vector3D E0,H0; double H2, wH0;
  for(ix=0;ix<Lx;ix++) //para cada celda
    for(iy=0;iy<Ly;iy++)
      for(iz=0;iz<Lz;iz++){
	//Compute the constants
	mur0=mur(ix,iy,iz); epsilonr0=epsilonr(ix,iy,iz); 
	//Compute the fields
	E0=E(ix,iy,iz,false); H0=H(ix,iy,iz,false);
	//Computer energy density
	//wH0 = wH(ix,iy,iz,false);
     	H2 = norma2(H0); wH0 = (mur0*mur0)*Mu0*H2/2.0;
	//BGK evolution rule
	for(alpha=0;alpha<13;alpha++)
	  for(i=0;i<3;i++){
	    fnew[ix][iy][iz][alpha][i]=UmUTau*f[ix][iy][iz][alpha][i]+UTau*feq(E0,H0,wH0,epsilonr0,mur0,alpha,i);
	    gnew[ix][iy][iz][alpha][i]=UmUTau*g[ix][iy][iz][alpha][i]+UTau*geq(E0,H0,wH0,epsilonr0,mur0,alpha,i);
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
        for(alpha=0;alpha<13;alpha++){
            ixnew=(int(ix+v[alpha][0]+Lx))%Lx; iynew=(int(iy+v[alpha][1]+Ly))%Ly; iznew=(int(iz+v[alpha][2]+Lz))%Lz;
	    for(i=0;i<3;i++){
	      f[ixnew][iynew][iznew][alpha][i]=fnew[ix][iy][iz][alpha][i];
	      g[ixnew][iynew][iznew][alpha][i]=gnew[ix][iy][iz][alpha][i];
	    }
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
  int t, tmax=20;
  
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
