#ifndef _plasma_h
#define _plasma_h

#include<eigen3/Eigen/Core>

#include"particles.h"
#include"poisson_solver_2d.h"
#include"bfield.h"
#include"input.h"

using namespace std;
using namespace Eigen;

class PlasmaSystem: public Input
{
  public:
    //PoissonSolver2D_DirichletBC poisson_solver;
    PoissonSolver2D_PeriodicBC poisson_solver;
    ConstMagneticField B;
    MatrixXd charge;

    vector<double> Ek;
    vector<double> Ep;
    vector<double> Et;
    vector<double> T;
    vector<double> num_in_grid;

    PlasmaSystem();

    //print important parameters in simulation
    void PrintParameters() const;

    //trace energy during simulation
    void CalculateE();

    //Setup Charge
    void SetupSpeciesChargeOnGrids();
    void SetupBackgroundChargeOnGrids();

    void Run();

};
#endif
