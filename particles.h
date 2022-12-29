#ifndef _particles_h
#define _particles_h
#include<vector>
#include<string>
using namespace std;

class PhaseSpace
{
  public:
    double x;
    double y;
    double vx;
    double vy;

    PhaseSpace(double _x = 0, double _y = 0.0, double _vx = 0.0, double _vy = 0.0);
    ~PhaseSpace(void);
};
class Particles
{
  public:
    const int num;//alpha-particle numbers
    const double q;//alpha-particle charge
    const double m;//alpha-particle mass
    const string name;//alpha-particle name
    vector<PhaseSpace> rv;

    Particles(double N, double Q, double M, string name);
    void InitializeXV_Random(double (*Distribution)(double, double, double, double), double v_max, double Lx, double Ly);
    //void InitializeXV_Quiet(double (*DistributionX)(double), double (*DistributionV)(double), double L, double dx, int nx_grids);
};

#endif
