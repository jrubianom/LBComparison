#include<iostream>
#include<cmath>

using namespace std;
class Parameter{
    private:

    public:
        //------//Dimensions
        int Lx = 1,Ly = 1,Lz = 50;
        int scale = 1;
        //------------
        double Epsilon0=1, Mu0=2;
        double C=1.0/sqrt(Epsilon0*Mu0);
        double Sigma0=0.0;
        //-----------------
        double E00=0.001,B00=E00/C;
        double epsr1 = 1, epsr2 = 2.5;

        void Update(int Lx0,int Ly0,int Lz0,int scale0,double Epsilon00,
                        double Mu00,double Sigma00,double Ei);

        friend class LatticeBoltzmann;
};

void Parameter::Update(int Lx0,int Ly0,int Lz0,int scale0,double Epsilon00,
                     double Mu00,double Sigma00,double Ei){
    Lx=Lx0; Ly=Ly0; Lz=Lz0;
    scale = scale0;
    Epsilon0=Epsilon00; Mu0=Mu00;  C=1.0/sqrt(Epsilon0*Mu0); Sigma0=Sigma00;
    E00=Ei;
}
