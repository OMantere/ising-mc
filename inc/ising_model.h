#include <bits/stdc++.h>
#include "constants.h"

using namespace std;

enum class Initialization { PositiveConstant, NegativeConstant, Random };

struct ModelOptions {
    ModelOptions() : T(273.15), SIZE_X(10), SIZE_Y(10) {
        beta = 1/T/kB;
        N = SIZE_X * SIZE_Y;
        J = 1e-20 / N;
        initialization = Initialization::PositiveConstant;
    }
    void setTemperature(double temp) {
        T = temp;
        beta = 1/T/kB;
    }
    void setGridSize(int x_size, int y_size) {
        SIZE_Y = y_size;
        SIZE_X = x_size;
        N = SIZE_Y * SIZE_X;
    }
    double beta;
    double T;
    double J;
    int SIZE_X;
    int SIZE_Y;
    int N;
    Initialization initialization;
};

class IsingModel {
    public:
        ModelOptions options;
        IsingModel(ModelOptions);
        IsingModel();
        void Step();
        void PrintGrid();
        double deltaEnergy();
        double SystemEnergy();
        double SystemMagnetization();
    private:
        int g(int, int);
        template<typename T>
        T GridLoop(function<T(int x, int y)> f);
        vector<int> grid;
        double previousSystemEnergy;
};