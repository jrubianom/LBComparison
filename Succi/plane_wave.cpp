#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

//----------- Constants -----------

// Spatial Dimensions
const int D=2;

// Domain size
const int Lx=200;
const int Ly=1;
const int Lz=1;

// Time step
const double Dt=1.0;

// Relaxation time
 const double tau=0.5;

// Physical constants
const double Epsilon0=1.0, Mu0=1.0;

// Properties of the medium
double mur(int ix, int iy, int iz){
    return 1.0;
}
double epsilonr(int ix, int iy, int iz){
    return 1.0;
}
double c2(int ix, int iy, int iz){
    double mur0=mur(ix, iy, iz), epsilonr0=epsilonr(ix, iy, iz);
    return 1.0/(Mu0*mur0*Epsilon0*epsilonr0);
}
//----------- LB Class -----------
class LB{
    private:
        // Velocity vectors
        double v[2*D+1][3];
        
        //Vacuum distribution functions:
        double g[Lx][Ly][Lz][2*D+1][3][3];
        double gnew[Lx][Ly][Lz][2*D+1][3][3];
        //std::vector<std::vector<std::vector<double>>> g(2*D+1, std::vector<std::vector<double>>(3, std::vector<double>(3, 0)));
        //std::vector<double> gnew;
    
    public:
        LB(void); // Constructor
        int index(int ix, int iy, int iz, int a, int b, int i);
        // Lattice weights (defined at the end of page 2 of the paper)
        double weight_0(int ix, int iy, int iz);
        double weight_i(int ix, int iy, int iz);
        // Macroscopic tensors
        double Lambda(int ix, int iy, int iz, int a, int b, bool UseNew);
        double B(int ix, int iy, int iz, int b, bool UseNew);
        // Sources
        double T(int ix, int iy, int iz, int a, int b, int i);
        // Equilibrium functions
        double geq(int ix, int iy, int iz, int a, int b, int i, double lambda_12, double B_1, double B_2);
        // Simulation functions
        void Collision(void);
        void ImposeFields(int t);
        void Advection(void);
        void Print(int t);
        void Print_geq(void);
        void Print_debug(int t);
};

//----------- LB Functions -----------
LB::LB(void){
    // Initialize containers
    // (2*D+1)*9 functions per cell
    // Lx*Ly*Lz cells in the domain
    for(int ix=0; ix<Lx; ix++)
        for(int iy=0; iy<Ly; iy++)
            for(int iz=0; iz<Lz; iz++)
                for(int i=0; i<2*D+1; i++)
                    for(int a=0; a<3; a++)
                        for(int b=0; b<3; b++){
                            g[ix][iy][iz][i][a][b] = 0;
                            gnew[ix][iy][iz][i][a][b] = 0;
            }

    // Initialize velocity vectors
    v[0][0] = 0.0; v[0][1] = 0.0; v[0][2] = 0.0;
    v[1][0] = 1.0; v[1][1] = 0.0; v[1][2] = 0.0;
    v[2][0] = 0.0; v[2][1] = 1.0; v[2][2] = 0.0;
    v[3][0] =-1.0; v[3][1] = 0.0; v[3][2] = 0.0;
    v[4][0] = 0.0; v[4][1] =-1.0; v[4][2] = 0.0;
    //v[5][0] = 0.0; v[5][1] = 0.0; v[5][2] = 1.0;
    //v[6][0] = 0.0; v[6][1] = 0.0; v[6][2] =-1.0;
}
// Functions to index the arrays
int LB::index(int ix, int iy, int iz, int a, int b, int i){
    /*
    Returns the position of the (a,b) component
    of the ith distribution function in the cell (ix,iy,iz)

    Parameters
    ----------
    ix : x index of the cell
    iy : y index of the cell
    iz : z index of the cell
    i : Distribution function number
    a : First component
    b : Second component
    */
    return (iz*Lx*Ly+iy*Lx+ix)*(2*D+1)*9 + a*(2*D+1)*3 + b*(2*D+1) + i;
}
// Velocity weights
double LB::weight_0(int ix, int iy, int iz){
    /*
    Returns the weight of the rest velocity at cell (ix,iy,iz).
    It is computed from the speed of light, as described on page 3
    of the paper. Note that I removed the factor of 3 accompanying c2;
    the reason is that the weights summation used to define c2
    does not actually yields this factor of 3.

    Parameters
    ----------
    ix : x index of the cell
    iy : y index of the cell
    iz : z index of the cell
    */
    return 1.0 - c2(ix, iy, iz);
}
double LB::weight_i(int ix, int iy, int iz){
    /*
    Returns the weight of the streaming velocities at cell (ix,iy,iz).
    They're computed from w_0, as described on page 3 of the paper.

    Parameters
    ----------
    ix : x index of the cell
    iy : y index of the cell
    iz : z index of the cell
    */
    return (1.0 - weight_0(ix, iy, iz))/(2.0*D);
}
// Macroscopic tensors
double LB::Lambda(int ix, int iy, int iz, int a, int b, bool UseNew){
    /*
    Returns the component (a,b) of the antisymmetric tensor Lambda,
    at cell (ix,iy,iz).
    It is computed from equation (9) of the paper.

    Parameters
    ----------
    ix : x index of the cell
    iy : y index of the cell
    iz : z index of the cell
    a : First component
    b : Second component
    */
    double sum = 0;
    for(int i=0; i<2*D+1; i++){
        if(UseNew)
            sum += gnew[ix][iy][iz][i][a][b];
        else
            sum += g[ix][iy][iz][i][a][b];
    }
    return c2(ix, iy, iz)*sum;
}
double LB::B(int ix, int iy, int iz, int b, bool UseNew){
    /*
    Returns the component (b) of the magnetic field,
    at cell (ix,iy,iz).
    It is computed from equation (10) of the paper, with alpha = gamma
    and gamma != beta. Note that alpha is completely arbitrary, here I
    choose it to be the component following b: x -> y 
                                               y -> x
                                               z -> x

    Parameters
    ----------
    ix : x index of the cell
    iy : y index of the cell
    iz : z index of the cell
    b : Magnetic field component to compute
    */
    double sum = 0;
    int a = (b+1)%3; // alpha
    for(int i=0; i<2*D+1; i++){
        if(UseNew)
            sum += gnew[ix][iy][iz][i][a][b]*v[i][a];
        else
            sum += g[ix][iy][iz][i][a][b]*v[i][a];
    }
    return sum;
}
// Sources
double LB::T(int ix, int iy, int iz, int a, int b, int i){
    return 0;    
}
// Equilibrium functions
double LB::geq(int ix, int iy, int iz, int a, int b, int i, double lambda_12, double B_1, double B_2){
    /*
    Returns the component (i,a,b) of the equillibrium distribution function tensor,
    at cell (ix,iy,iz).
    It is computed from equation (8) of the paper.
    
    Parameters
    ----------
    ix : x index of the cell
    iy : y index of the cell
    iz : z index of the cell
    a : First spatial component
    b : Second spatial component
    i : Streaming index
    lambda_12 : Component (a,b) of the tensor Lambda at the current cell
    B_1 : Component (a) of the magnetic field at the current cell
    B_2 : Component (b) of the magnetic field at the current cell
    */
    double v_i1 = v[i][a];
    double v_i2 = v[i][b];
    double w_i = weight_i(ix, iy, iz);
    double c_sqrd = c2(ix, iy, iz);

    double sum = lambda_12 + v_i1*B_2 - v_i2*B_1;

    return w_i*sum/c_sqrd;
}
//----------- Simulation Functions -----------
void LB::Collision(void){
    /*
    Performs the collision step, as described in equation (17) of the paper.
    The post-collisional state is stored in the array gnew.
    */
    // For each cell
    for(int ix=0; ix<Lx; ix++)
        for(int iy=0; iy<Ly; iy++)
            for(int iz=0; iz<Lz; iz++)
                // For each distribution function
                for(int i=0; i<2*D+1; i++)
                    for(int a=0; a<3; a++)
                        for(int b=0; b<3; b++)
                        {
                            // Compute fields
                            double lambda_12 = Lambda(ix, iy, iz, a, b, false);
                            double B_2 = B(ix, iy, iz, b, false);
                            double B_1 = B(ix, iy, iz, a, false);
                            // Distribution at current position
                            double g0 = g[ix][iy][iz][i][a][b];
                            // Equilibrium
                            double g_equi = geq(ix, iy, iz, a, b, i, lambda_12, B_1, B_2);
                            // Sources
                            double T0 = T(ix, iy, iz, a, b, i);
                            // BGK Evolution
                            gnew[ix][iy][iz][i][a][b]=(1.0-1.0/tau)*g0 + (1.0/tau)*g_equi + T0*Dt;
                            }
}
void LB::ImposeFields(int t){
    /*
    Imposes a gaussian pulse at the beggining of the simulation.
    It is defined from equations (24) of Miller´s paper, but for 
    propagation in the x direction.
    */
    // Impose gnew=geq with the desired fields
    double amp = 0.1, omega=M_PI/100;
    double lambda_12, B_2, B_1, v_i1, v_i2;
    // For each cell
    for(int ix=0; ix<Lx; ix++)
        for(int iy=0; iy<Ly; iy++)
            for(int iz=0; iz<Lz; iz++)
                // For each distribution function
                for(int i=0; i<2*D+1; i++)
                    for(int a=0; a<3; a++)
                        for(int b=0; b<3; b++){
                            // Electric field
                            // Only y component
                            double Ey = amp*std::sin(omega*t);
                            if(ix==0 && a==2 && b==0){
                                lambda_12 = -Ey;
                            } else if(ix==0 && a==0 && b==2){
                                lambda_12 = Ey;
                            } else{
                                lambda_12 = 0;
                            }
                            // Magnetic field
                            // Only z component
                            double c = std::sqrt(c2(ix,iy,iz));
                            double Bz = Ey/c;
                            if(ix==0 && b==2){
                                B_2 = Bz;
                            } else{
                                B_2 = 0;
                            }
                            if(ix==0 && a==2){
                                B_1 = Bz;
                            } else{
                                B_1 = 0;
                            }
                            // Impose field
                            gnew[ix][iy][iz][i][a][b] = geq(ix, iy, iz, a, b, i, lambda_12, B_1, B_2);
                        }

}
void LB::Advection(void){
    // For each cell
    for(int ix=0; ix<Lx; ix++)
        for(int iy=0; iy<Ly; iy++)
            for(int iz=0; iz<Lz; iz++)
                // For each distribution function
                for(int i=0; i<2*D+1; i++)
                    for(int a=0; a<3; a++)
                        for(int b=0; b<3; b++){
                            // Periodic boundary conditions
                            int ixnew=(ix+int(v[i][0])+Lx); ixnew = ixnew%Lx;
                            int iynew=(iy+int(v[i][1])+Ly); iynew = iynew%Ly;
                            int iznew=(iz+int(v[i][2])+Lz); iznew = iznew%Lz;
                            
                            g[ixnew][iynew][iznew][i][a][b]=gnew[ix][iy][iz][i][a][b];
                }
}
void LB::Print(int t){
    // Names
    std::string E_string("E_" + std::to_string(t) + ".dat");
    std::string B_string("B_" + std::to_string(t) + ".dat");
    // Output files
    std::fstream E_file;
    E_file.open(E_string, std::ios_base::out);
    std::fstream B_file;
    B_file.open(B_string, std::ios_base::out);
    // For each cell
    for(int iz=0; iz<Lz; iz++)
        for(int iy=0; iy<Ly; iy++)
            for(int ix=0; ix<Lx; ix++){
                // Electric field
                double Ex = Lambda(ix, iy, iz, 2, 1, true);
                double Ey = Lambda(ix, iy, iz, 0, 2, true);
                double Ez = Lambda(ix, iy, iz, 1, 0, true);
                // Magnetic field
                double Bx = B(ix, iy, iz, 0, true);
                double By = B(ix, iy, iz, 1, true);
                double Bz = B(ix, iy, iz, 2, true);
                // Print
                E_file << Ex << ' ' << Ey << ' ' << Ez << '\n';
                B_file << Bx << ' ' << By << ' ' << Bz << '\n';
            }
    E_file.close();
    B_file.close();
}
int main(void){
    LB Pulse;
    int t, tmax=10;

    Pulse.ImposeFields(0);

    for(t=0; t<=tmax; t++){
        Pulse.Collision();
        Pulse.ImposeFields(t);
        Pulse.Advection();
        if(t%2==0) Pulse.Print(t);
    }

    return 0;
}