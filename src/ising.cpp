/*
 * Exercise 28.7 from Concepts in Thermal Physics, 2nd edition (Blundell & Blundell)
 */

#include <bits/stdc++.h>
#include "../inc/ising_model.h"
#include "../inc/constants.h"

using namespace std;


int stabilizationDuration = 100000;
int averagingDuration = 100;
int M = 50;
vector<double> E_t(M);
vector<double> E2_t(M);
vector<double> m_t(M);
vector<double> C_t(M);
vector<double> T_t(M);

void print() {
    cout << endl;
}

template<typename T>
void print(T t) {
    cout << t << endl;
}

template<typename T, typename... Args>
void print(T t, Args... args) {
    cout << t << " ";
    print(args...);
}

void iteration(int j) {
    ModelOptions options = ModelOptions();
    options.J = 1e-22;
    options.setGridSize(20, 20);
    int T_start = 30;
    int T_end = 100;
    double deltaT = (T_end - T_start)/(double)M;
    double temp = T_start + j*deltaT;
    options.setTemperature(temp);
    IsingModel model = IsingModel(options);
    for(int i = 0; i < stabilizationDuration; i++) {
        model.Step();
    }
    double E;
    double E2;
    double m;
    double C;
    for(int i = 0; i < averagingDuration; i++) {
        model.Step();
        double energy = model.SystemEnergy();
        E += energy;
        E2 += pow(energy, 2);
        m += model.SystemMagnetization();
    }
    E /= averagingDuration;
    E2 /= averagingDuration;
    m /= averagingDuration;
    C = kB * pow(model.options.beta, 2) * (E2 - pow(E, 2));
    E_t[j] = E;
    E2_t[j] = E2;
    m_t[j] = m;
    C_t[j] = C;
    T_t[j] = temp;
    model.PrintGrid();
    print("Expected energy:", E);
    print("Expected magnetization:", m);
    print("Heat capacity:", C);
    print();
}

int main() {
    srand(time(0));
    vector<thread> threads;

    for(int j = 0; j < M; j++)
        threads.push_back(thread(iteration, j));

    for(int j = 0; j < M; j++)
        threads[j].join();

    ofstream f;
    f.open("results");
    for(int j = 0; j < M; j++)
        f << E_t[j] << "\t" << m_t[j] << "\t" << C_t[j] << "\t" << T_t[j] << endl;
    f.close();
    return 0;
}
