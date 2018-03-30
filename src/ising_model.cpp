#include "../inc/ising_model.h"

IsingModel::IsingModel(ModelOptions options) : options(options), grid(options.N) {
    for (auto& spin : grid) {
        switch(options.initialization) {
            case Initialization::Random:
                if(rand() % 2)
                    spin = -1;
                else
                    spin = 1;
                break;
            case Initialization::PositiveConstant:
                spin = 1;
                break;
            case Initialization::NegativeConstant:
                spin = -1;
                break;
            default:
                throw invalid_argument("Unknown enum");
        }
    }
    previousSystemEnergy = SystemEnergy();
}

IsingModel::IsingModel() : IsingModel(ModelOptions()) {}

void IsingModel::PrintGrid() {
    int i = 0;
    while(i < options.N) {
        if(grid[i] == -1)
            cout << "-1" << " ";
        else
            cout << " 1" << " ";
        i++;
        if(!(i % options.SIZE_X))
            cout << "\n";
    }
}

int IsingModel::g(int x, int y) {
    return grid[y * options.SIZE_X + x];
}

template<typename T>
T IsingModel::GridLoop(function<T(int x, int y)> f) {
    int x = 0;
    int y = 0;
    T acc = 0;
    while(y < options.SIZE_Y) {
        acc += f(x, y);
        x++;
        if(x == options.SIZE_X) {
            x = 0;
            y++;
        }
    }
    return acc;
}

double IsingModel::SystemEnergy() {
    return GridLoop<double>([&](int x, int y) {
        int s = 0;
        int si = g(x, y);
        if(y > 0)
            s += g(x, y - 1) * si;
        if(y < options.SIZE_Y - 1)
            s += g(x, y + 1) * si;
        if(x > 0) {
            s += g(x - 1, y) * si;
            if(y > 0)
                s += g(x - 1, y - 1) * si;
            if(y < options.SIZE_Y - 1)
                s += g(x - 1, y + 1) * si;
        }
        if(x < options.SIZE_X - 1) {
            s += g(x + 1, y) * si;
            if(y > 0)
                s += g(x + 1, y - 1) * si;
            if(y < options.SIZE_Y - 1)
                s += g(x + 1, y + 1) * si;
        }
        return -options.J * s;
    });
}

double IsingModel::SystemMagnetization() {
    return GridLoop<double>([&](int x, int y) {
        return (double)g(x, y) / options.N;
    });
}

void IsingModel::Step() {
    int spin = rand() % options.N;
    grid[spin] *= -1;
    double systemEnergy = SystemEnergy();
    double deltaE = systemEnergy - previousSystemEnergy;
    double r = (double)rand()/RAND_MAX;
    if(deltaE >= 0 && 1 - exp(-options.beta * deltaE) > r) {
        grid[spin] *= -1;
    } else {
        previousSystemEnergy = systemEnergy;
    }
}
