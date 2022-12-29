/******************************************
 * 2D Solver of Poisson equation
 * Author: Ran Guo
 * Data: 2022-12-12
 * Note: Only for Dirichlet Boundary Condition
******************************************/

#include<eigen3/Eigen/Core>
#include<eigen3/Eigen/Dense>
#include<eigen3/Eigen/QR>
#include"poisson_solver_2d.h"

using namespace Eigen;

PoissonSolver2D_DirichletBC::PoissonSolver2D_DirichletBC(int _nx_grids, double _dx, int _ny_grids, double _dy, VectorXd _phi_0k, VectorXd _phi_j0, VectorXd _phi_n1k, VectorXd _phi_jn1):
    nx_grids(_nx_grids), dx(_dx), ny_grids(_ny_grids), dy(_dy), Lx(_nx_grids * dx), Ly(_ny_grids * dy), phi_0k(_phi_0k), phi_j0(_phi_j0), phi_n1k(_phi_n1k), phi_jn1(_phi_jn1)
{
    phi.resize(nx_grids, ny_grids);
    Ex.resize(nx_grids, ny_grids);
    Ey.resize(nx_grids, ny_grids);
    phi.row(0) = _phi_0k;
    phi.col(0) = _phi_j0;
    phi.row(nx_grids - 1) = _phi_n1k;
    phi.col(ny_grids - 1) = _phi_jn1;

    //construct oprator A
    double dx2 = dx * dx;
    double dy2 = dy * dy;

    B.resize(nx_grids - 2, ny_grids - 2);
    B.setZero();
    B.diagonal() = VectorXd::Constant(nx_grids - 2, 2.0 / dx2 + 2.0 / dy2);
    B.diagonal(1) = VectorXd::Constant(nx_grids - 3,  -1.0 / dx2);
    B.diagonal(-1) = VectorXd::Constant(nx_grids - 3, -1.0 / dx2);

    C.resize(nx_grids - 2, ny_grids - 2);
    C = MatrixXd::Identity(nx_grids - 2, ny_grids - 2) / dy2;

    A.resize((nx_grids - 2) * (ny_grids - 2), (nx_grids - 2) * (ny_grids - 2));
    for(int i = 0; i < nx_grids - 2; i++)
    {
        int idx = i * (nx_grids - 2);
        A.block(idx, idx, nx_grids - 2, nx_grids - 2) = B;
        if(i < nx_grids - 3)
        {
            A.block(idx + nx_grids - 2, idx, nx_grids - 2, nx_grids - 2) = -C;
            A.block(idx, idx + nx_grids - 2, nx_grids - 2, nx_grids - 2) = -C;
        }
    }
}

void PoissonSolver2D_DirichletBC::Solve(MatrixXd charge)
{
    double dx2 = dx * dx;
    double dy2 = dy * dy;

    MatrixXd charge_modified = charge.block(1, 1, nx_grids - 2, ny_grids - 2);
    charge_modified.row(0) += phi.row(0).segment(1, nx_grids - 2) / dx2;
    charge_modified.row(nx_grids - 3) += phi.row(nx_grids - 1).segment(1, nx_grids - 2) / dx2;
    charge_modified.col(0) += phi.col(0).segment(1, nx_grids - 2) / dy2;
    charge_modified.col(ny_grids - 3) += phi.col(ny_grids - 1).segment(1, nx_grids - 2) / dy2;

    VectorXd phi_vec = A.fullPivLu().solve(charge_modified.reshaped());
    phi.block(1, 1, nx_grids - 2, ny_grids - 2) = phi_vec.reshaped(nx_grids - 2, ny_grids - 2).eval();
}

double PoissonSolver2D_DirichletBC::GetEx(int x_idx, int y_idx)
{
    int x_idx_p = (x_idx + 1 > nx_grids - 1) ? x_idx : x_idx + 1;
    int x_idx_m = (x_idx - 1 < 0) ? x_idx : x_idx - 1;
    double Ex = (phi(x_idx_p, y_idx) - phi(x_idx_m, y_idx)) / 2.0 / dx;
    return Ex;
}

double PoissonSolver2D_DirichletBC::GetEy(int x_idx, int y_idx)
{
    int y_idx_p = (y_idx + 1 > ny_grids - 1) ? y_idx : y_idx + 1;
    int y_idx_m = (y_idx - 1 < 0) ? y_idx : y_idx - 1;
    double Ex = (phi(x_idx, y_idx_p) - phi(x_idx, y_idx_m)) / 2.0 / dx;
    return Ex;
}
