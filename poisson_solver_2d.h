#ifndef _poisson_solver_2d_h
#define _poisson_solver_2d_h
#include<eigen3/Eigen/Core>
using namespace Eigen;

class PoissonSolver2D_DirichletBC
{
    public:
        const int nx_grids;
        const int ny_grids;
        const double dx;
        const double dy;
        const double Lx;
        const double Ly;

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

        PoissonSolver2D_DirichletBC(int _nx_grids, double _dx, int _ny_grids, double _dy, VectorXd _phi_0k, VectorXd _phi_j0, VectorXd _phi_n1k, VectorXd _phi_jn1);
        void Solve(MatrixXd _charge);
        double GetEx(int x_idx, int y_idx);
        double GetEy(int x_idx, int y_idx);
        double GetPhi(int x_idx, int y_idx);
};

#endif
