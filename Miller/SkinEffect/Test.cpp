#include <iostream>
#include <algorithm>
#include <vector>
#include <math.h>

using namespace std;
const int Lz=12;

void Amplitud(vector<double>& Amplitud, vector<double>& oldAmplitud);

void MeasureAmplitude(vector<double> &EAmplit){
    for(int iz=0; iz < Lz ; iz++)
        EAmplit[iz] = iz;
}


int main() {
    /*
    int numE = Lz;//Lz/2 - Lz/4;
    vector<double> EAmplit(numE,0);
    vector<double> OldEAmplit(numE,0);
    MeasureAmplitude(OldEAmplit);
    */

    vector<double> E = {1.0, -2.0, 3.0, 4.0,-4,9};
    vector<double> E_old = {5.0, 4.0, 3.0, 2.0,-9,-4};
    E_old = E;
    E[2]=100;
    E_old[0] = -90;
    for (auto e : E) {
        std::cout << e << " ";
    }
    cout << endl;
    for (auto e : E_old) {
        std::cout << e << " ";
    }



    return 0;
}

void Amplitud(vector<double>& Amplitud,  vector<double>& oldAmplitud) {
  transform(oldAmplitud.begin(), oldAmplitud.end(), Amplitud.begin(), Amplitud.begin(), [](double a, double b) {
    double a_abs = abs(a), b_abs = abs(b);
    return (a_abs > b_abs) ? a_abs : b_abs;
  });
}
