#ifndef _poisson_solver_2d_h
#define _poisson_solver_2d_h
#include<eigen3/Eigen/Core>
using namespace Eigen;

class PoissonSolver2D_DirichletBC
{
  public:
    const int nx;
    const int ny;
    const double dx;
    const double dy;
    const double Lx;
    const double Ly;
    const double dx2;
    const double dy2;

    //BC
    VectorXd phi_0k;
    VectorXd phi_j0;
    VectorXd phi_n1k;
    VectorXd phi_jn1;

    //oprator
    MatrixXd A;
    MatrixXd B;
    MatrixXd C;

    //electric field
    MatrixXd Ex;
    MatrixXd Ey;
    MatrixXd phi;

    PoissonSolver2D_DirichletBC(int _nx, double _dx, int _ny, double _dy, VectorXd _phi_0k, VectorXd _phi_j0, VectorXd _phi_n1k, VectorXd _phi_jn1);
    void Solve(MatrixXd _charge);

    double GetEx(int x_idx, int y_idx)
    {
        return Ex(x_idx, y_idx);
    }
    double GetEy(int x_idx, int y_idx)
    {
        return Ey(x_idx, y_idx);
    }
    double GetPhi(int x_idx, int y_idx)
    {
        return phi(x_idx, y_idx);
    }
};

class PoissonSolver2D_PeriodicBC
{
  public:
    const int nx;
    const int ny;
    const double dx;
    const double dy;
    const double Lx;
    const double Ly;
    const double dx2;
    const double dy2;

    //oprator
    MatrixXd A;
    MatrixXd B;
    MatrixXd C;

    //electric field
    MatrixXd Ex;
    MatrixXd Ey;
    MatrixXd phi;

    PoissonSolver2D_PeriodicBC(int _nx, double _dx, int _ny, double _dy);
    void Solve(MatrixXd _charge);

    double GetEx(int x_idx, int y_idx)
    {
        return Ex(x_idx, y_idx);
    }
    double GetEy(int x_idx, int y_idx)
    {
        return Ey(x_idx, y_idx);
    }
    double GetPhi(int x_idx, int y_idx)
    {
        return phi(x_idx, y_idx);
    }
};

class PoissonSolver2D_XPeriodic_YDirichletBC
{
  public:
    const int nx;
    const int ny;
    const double dx;
    const double dy;
    const double Lx;
    const double Ly;
    const double dx2;
    const double dy2;

    //BC
    VectorXd phi_j0;
    VectorXd phi_jn1;

    //oprator
    MatrixXd A;
    MatrixXd B;
    MatrixXd C;

    //electric field
    MatrixXd Ex;
    MatrixXd Ey;
    MatrixXd phi;

    PoissonSolver2D_XPeriodic_YDirichletBC(int _nx, double _dx, int _ny, double _dy, VectorXd _phi_j0, VectorXd _phi_jn1);
    void Solve(MatrixXd _charge);

    double GetEx(int x_idx, int y_idx)
    {
        return Ex(x_idx, y_idx);
    }
    double GetEy(int x_idx, int y_idx)
    {
        return Ey(x_idx, y_idx);
    }
    double GetPhi(int x_idx, int y_idx)
    {
        return phi(x_idx, y_idx);
    }
};

#endif
