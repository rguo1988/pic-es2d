/******************************************
 * 2D Solver of Poisson equation
 * Author: Ran Guo
 * Data: 2022-12-12
 * Note: Dirichlet, Periodic, and Complex Boundary Condition
******************************************/

#include<eigen3/Eigen/Core>
#include<eigen3/Eigen/Dense>
#include<eigen3/Eigen/QR>
#include<eigen3/Eigen/SparseCholesky>
#include"poisson_solver_2d.h"
#include<vector>

using namespace Eigen;

//Dirichlet Boundary Condition
PoissonSolver2D_DirichletBC::PoissonSolver2D_DirichletBC(int _nx, double _dx, int _ny, double _dy, VectorXd _phi_0k, VectorXd _phi_j0, VectorXd _phi_n1k, VectorXd _phi_jn1):
    nx(_nx), ny(_ny), dx(_dx), dy(_dy), Lx(_nx * _dx), Ly(_ny * _dy),
    dx2(_dx * _dx), dy2(_dy * _dy),
    phi_0k(_phi_0k), phi_j0(_phi_j0), phi_n1k(_phi_n1k), phi_jn1(_phi_jn1)
{
    phi.resize(nx, ny);
    Ex.resize(nx, ny);
    Ey.resize(nx, ny);
    phi.row(0) = _phi_0k;
    phi.col(0) = _phi_j0;
    phi.row(nx - 1) = _phi_n1k;
    phi.col(ny - 1) = _phi_jn1;

    /*********** construct oprator A (dx=dy=1) ************
    B =  4 -1  0  0     C = -1  0  0  0     A =  B  C  0  0
        -1  4 -1  0          0 -1  0  0          C  B  C  0
         0 -1  4 -1          0  0 -1  0          0  C  B  C
         0  0 -1  4          0  0  0 -1          0  0  C  B
    *******************************************************/
    /*
    B.resize(nx - 2, nx - 2);
    B.setZero();
    B.diagonal() = VectorXd::Constant(nx - 2, 2.0 / dx2 + 2.0 / dy2);
    B.diagonal(1) = VectorXd::Constant(nx - 3,  -1.0 / dx2);
    B.diagonal(-1) = VectorXd::Constant(nx - 3, -1.0 / dx2);

    C.resize(nx - 2, nx - 2);
    C = -1.0 * MatrixXd::Identity(nx - 2, nx - 2) / dy2;

    A.resize((nx - 2) * (ny - 2), (nx - 2) * (ny - 2));
    A.setZero();
    for(int i = 0; i < ny - 2; i++)
    {
        int idx = i * (nx - 2);
        A.block(idx, idx, nx - 2, nx - 2) = B;
        if(i != ny - 3)
        {
            A.block(idx + nx - 2, idx, nx - 2, nx - 2) = C;
            A.block(idx, idx + nx - 2, nx - 2, nx - 2) = C;
        }
    }
    */
    //construct Sparse SpA
    SpA.resize((nx - 2) * (ny - 2), (nx - 2) * (ny - 2));
    //SpA = A.sparseView();
    //construct Sparse SpA by setFromTriplets
    std::vector<Triplet<double>>triplet_SpA;
    for(int i = 0; i < (nx - 2) * (ny - 2); i++)
    {
        //B diagonal
        triplet_SpA.push_back(Triplet<double>(i, i, 2.0 / dx2 + 2.0 / dy2));
        //B sub-diagonal; Noting that it is B's sub-diagnoal rather than A's sub-diagonal
        if( (i + 1) % (nx - 2) != 0 )
        {
            triplet_SpA.push_back(Triplet<double>(i, i + 1, -1.0 / dx2));
            triplet_SpA.push_back(Triplet<double>(i + 1, i, -1.0 / dx2));
        }
        //C
        if( i + nx - 2 < (nx - 2) * (ny - 2) )
        {
            triplet_SpA.push_back(Triplet<double>(i, i + nx - 2, -1.0 / dy2));
            triplet_SpA.push_back(Triplet<double>(i + nx - 2, i, -1.0 / dy2));
        }
    }
    SpA.setFromTriplets(triplet_SpA.begin(), triplet_SpA.end());
}

void PoissonSolver2D_DirichletBC::Solve(MatrixXd charge)
{
    MatrixXd charge_modified = charge.block(1, 1, nx - 2, ny - 2);
    charge_modified.row(0) += phi.row(0).segment(1, ny - 2) / dx2;
    charge_modified.row(nx - 3) += phi.row(nx - 1).segment(1, ny - 2) / dx2;
    charge_modified.col(0) += phi.col(0).segment(1, nx - 2) / dy2;
    charge_modified.col(ny - 3) += phi.col(ny - 1).segment(1, nx - 2) / dy2;

    //solving by DenseMatrix Method
    //VectorXd phi_vec = A.llt().solve(charge_modified.reshaped());

    //solving SparseMatrix Method
    SimplicialLDLT<SparseMatrix<double>> SpA_solver;
    SpA_solver.compute(SpA);
    VectorXd phi_vec = SpA_solver.solve(charge_modified.reshaped());

    phi.block(1, 1, nx - 2, ny - 2) = phi_vec.reshaped(nx - 2, ny - 2).eval();

    //Calculate E
    #pragma omp parallel for
    for(int x_idx = 0; x_idx < nx; x_idx++)
    {
        for(int y_idx = 0; y_idx < ny; y_idx++)
        {
            int x_idx_p = (x_idx == nx - 1) ? x_idx : x_idx + 1;
            int x_idx_m = (x_idx == 0) ? x_idx : x_idx - 1;
            int y_idx_p = (y_idx == ny - 1) ? y_idx : y_idx + 1;
            int y_idx_m = (y_idx == 0) ? y_idx : y_idx - 1;
            Ex(x_idx, y_idx) = (phi(x_idx_m, y_idx) - phi(x_idx_p, y_idx)) / 2.0 / dx;
            Ey(x_idx, y_idx) = (phi(x_idx, y_idx_m) - phi(x_idx, y_idx_p)) / 2.0 / dy;
        }
    }
}

//Periodic Boundary Condition
PoissonSolver2D_PeriodicBC::PoissonSolver2D_PeriodicBC(int _nx, double _dx, int _ny, double _dy):
    nx(_nx), ny(_ny), dx(_dx), dy(_dy), Lx(_nx * _dx), Ly(_ny * _dy),
    dx2(_dx * _dx), dy2(_dy * _dy)
{
    phi.resize(nx, ny);
    Ex.resize(nx, ny);
    Ey.resize(nx, ny);

    /*********** construct oprator A (dx=dy=1) ************
    B =  4 -1  0 -1     C = -1  0  0  0     A =  B  C  0  C
        -1  4 -1  0          0 -1  0  0          C  B  C  0
         0 -1  4 -1          0  0 -1  0          0  C  B  C
        -1  0 -1  4          0  0  0 -1          C  0  C  B
    *******************************************************/

    /*
    B.resize(nx - 1, nx - 1);
    B.setZero();
    B.diagonal() = VectorXd::Constant(nx - 1, 2.0 / dx2 + 2.0 / dy2);
    B.diagonal(1) = VectorXd::Constant(nx - 2,  -1.0 / dx2);
    B.diagonal(-1) = VectorXd::Constant(nx - 2, -1.0 / dx2);
    B(0, nx - 2) = -1.0 / dx2;
    B(nx - 2, 0) = -1.0 / dx2;

    C.resize(nx - 1, nx - 1);
    C = -1.0 * MatrixXd::Identity(nx - 1, nx - 1) / dy2;

    A.resize((nx - 1) * (ny - 1), (nx - 1) * (ny - 1));
    A.setZero();
    for(int i = 0; i < ny - 1; i++)
    {
        int idx = i * (nx - 1);
        A.block(idx, idx, nx - 1, nx - 1) = B;
        if(i != ny - 2)
        {
            A.block(idx + nx - 1, idx, nx - 1, nx - 1) = C;
            A.block(idx, idx + nx - 1, nx - 1, nx - 1) = C;
        }
    }
    A.block(0, (nx - 1) * (ny - 2), nx - 1, nx - 1) = C;
    A.block((nx - 1) * (ny - 2), 0, nx - 1, nx - 1) = C;
    */

    ////construct Sparse SpA
    SpA.resize((nx - 1) * (ny - 1), (nx - 1) * (ny - 1));
    //SpA = A.sparseView();
    //construct Sparse SpA by setFromTriplets
    std::vector<Triplet<double>>triplet_SpA;
    for(int i = 0; i < (nx - 1) * (ny - 1); i++)
    {
        //B diagonal
        triplet_SpA.push_back(Triplet<double>(i, i, 2.0 / dx2 + 2.0 / dy2));
        //B sub-diagonal; Noting that it is B's sub-diagnoal rather than A's sub-diagonal
        if( (i + 1) % (nx - 1) != 0 )
        {
            triplet_SpA.push_back(Triplet<double>(i, i + 1, -1.0 / dx2));
            triplet_SpA.push_back(Triplet<double>(i + 1, i, -1.0 / dx2));
        }
        //C
        if( i + nx - 1 < (nx - 1) * (ny - 1) )
        {
            triplet_SpA.push_back(Triplet<double>(i, i + nx - 1, -1.0 / dy2));
            triplet_SpA.push_back(Triplet<double>(i + nx - 1, i, -1.0 / dy2));
        }
        if( i + (nx - 1) * (ny - 2) < (nx - 1) * (ny - 1) )
        {
            triplet_SpA.push_back(Triplet<double>(i, i + (nx - 1) * (ny - 2), -1.0 / dy2));
            triplet_SpA.push_back(Triplet<double>(i + (nx - 1) * (ny - 2), i, -1.0 / dy2));
        }
    }
    //B(0,nx-2) and B(nx-2,0)
    for(int i = 0; i < ny - 1; i++)
    {
        int idx = i * (nx - 1);
        triplet_SpA.push_back(Triplet<double>(idx, idx + nx - 2, -1.0 / dx2));
        triplet_SpA.push_back(Triplet<double>(idx + nx - 2, idx, -1.0 / dx2));
    }
    SpA.setFromTriplets(triplet_SpA.begin(), triplet_SpA.end());
}

void PoissonSolver2D_PeriodicBC::Solve(MatrixXd charge)
{
    MatrixXd charge_modified = charge.block(0, 0, nx - 1, ny - 1);

    //solving by DenseMatrix Method
    //VectorXd phi_vec = A.llt().solve(charge_modified.reshaped());

    //solving SparseMatrix Method
    SimplicialLDLT<SparseMatrix<double>> SpA_solver;
    SpA_solver.compute(SpA);
    VectorXd phi_vec = SpA_solver.solve(charge_modified.reshaped());

    phi.block(0, 0, nx - 1, ny - 1) = phi_vec.reshaped(nx - 1, ny - 1).eval();
    phi.col(ny - 1) = phi.col(0);
    phi.row(nx - 1) = phi.row(0);

    //Calculate E
    #pragma omp parallel for
    for(int x_idx = 0; x_idx < nx; x_idx++)
    {
        for(int y_idx = 0; y_idx < ny; y_idx++)
        {
            int x_idx_p = (x_idx == nx - 1) ? 1 : x_idx + 1;
            int x_idx_m = (x_idx == 0) ? nx - 2 : x_idx - 1;
            int y_idx_p = (y_idx == ny - 1) ? 1 : y_idx + 1;
            int y_idx_m = (y_idx == 0) ? ny - 2 : y_idx - 1;
            Ex(x_idx, y_idx) = (phi(x_idx_m, y_idx) - phi(x_idx_p, y_idx)) / 2.0 / dx;
            Ey(x_idx, y_idx) = (phi(x_idx, y_idx_m) - phi(x_idx, y_idx_p)) / 2.0 / dy;
        }
    }
}

//X Periodic + Y Dirichlet Boundary Condition
PoissonSolver2D_XPeriodic_YDirichletBC::PoissonSolver2D_XPeriodic_YDirichletBC(int _nx, double _dx, int _ny, double _dy, VectorXd _phi_j0, VectorXd _phi_jn1):
    nx(_nx), ny(_ny), dx(_dx), dy(_dy), Lx(_nx * _dx), Ly(_ny * _dy),
    dx2(_dx * _dx), dy2(_dy * _dy),
    phi_j0(_phi_j0), phi_jn1(_phi_jn1)
{
    phi.resize(nx, ny);
    Ex.resize(nx, ny);
    Ey.resize(nx, ny);
    phi.col(0) = _phi_j0;
    phi.col(ny - 1) = _phi_jn1;

    /*********** construct oprator A (dx=dy=1) ************
    B =  4 -1  0 -1     C = -1  0  0  0     A =  B  C  0  0
        -1  4 -1  0          0 -1  0  0          C  B  C  0
         0 -1  4 -1          0  0 -1  0          0  C  B  C
        -1  0 -1  4          0  0  0 -1          0  0  C  B
    *******************************************************/
    /*
    B.resize(nx - 1, nx - 1);
    B.setZero();
    B.diagonal() = VectorXd::Constant(nx - 1, 2.0 / dx2 + 2.0 / dy2);
    B.diagonal(1) = VectorXd::Constant(nx - 2,  -1.0 / dx2);
    B.diagonal(-1) = VectorXd::Constant(nx - 2, -1.0 / dx2);
    B(0, nx - 2) = -1.0 / dx2;
    B(nx - 2, 0) = -1.0 / dx2;

    C.resize(nx - 1, nx - 1);
    C = -1.0 * MatrixXd::Identity(nx - 1, nx - 1) / dy2;

    A.resize((nx - 1) * (ny - 2), (nx - 1) * (ny - 2));
    A.setZero();
    for(int i = 0; i < ny - 2; i++)
    {
        int idx = i * (nx - 1);
        A.block(idx, idx, nx - 1, nx - 1) = B;
        if(i != ny - 3)
        {
            A.block(idx + nx - 1, idx, nx - 1, nx - 1) = C;
            A.block(idx, idx + nx - 1, nx - 1, nx - 1) = C;
        }
    }
    */
    //construct Sparse SpA
    SpA.resize((nx - 1) * (ny - 2), (nx - 1) * (ny - 2));
    //SpA = A.sparseView();
    //construct Sparse SpA by setFromTriplets
    std::vector<Triplet<double>>triplet_SpA;
    for(int i = 0; i < (nx - 1) * (ny - 2); i++)
    {
        //B diagonal
        triplet_SpA.push_back(Triplet<double>(i, i, 2.0 / dx2 + 2.0 / dy2));
        //B sub-diagonal; Noting that it is B's sub-diagnoal rather than A's sub-diagonal
        if( (i + 1) % (nx - 1) != 0 )
        {
            triplet_SpA.push_back(Triplet<double>(i, i + 1, -1.0 / dx2));
            triplet_SpA.push_back(Triplet<double>(i + 1, i, -1.0 / dx2));
        }
        //C
        if( i + nx - 1 < (nx - 1) * (ny - 2) )
        {
            triplet_SpA.push_back(Triplet<double>(i, i + nx - 1, -1.0 / dy2));
            triplet_SpA.push_back(Triplet<double>(i + nx - 1, i, -1.0 / dy2));
        }
    }
    //B(0,nx-2) and B(nx-2,0)
    for(int i = 0; i < ny - 2; i++)
    {
        int idx = i * (nx - 1);
        triplet_SpA.push_back(Triplet<double>(idx, idx + nx - 2, -1.0 / dx2));
        triplet_SpA.push_back(Triplet<double>(idx + nx - 2, idx, -1.0 / dx2));
    }
    SpA.setFromTriplets(triplet_SpA.begin(), triplet_SpA.end());
}

void PoissonSolver2D_XPeriodic_YDirichletBC::Solve(MatrixXd charge)
{
    MatrixXd charge_modified = charge.block(0, 1, nx - 1, ny - 2);
    charge_modified.col(0) += phi.col(0).segment(0, nx - 1) / dy2;
    charge_modified.col(ny - 3) += phi.col(ny - 1).segment(0, nx - 1) / dy2;

    //solving by DenseMatrix Method
    //VectorXd phi_vec = A.llt().solve(charge_modified.reshaped());

    //solving by SparseMatrix Method
    SimplicialLDLT<SparseMatrix<double>> SpA_solver;
    SpA_solver.compute(SpA);
    VectorXd phi_vec = SpA_solver.solve(charge_modified.reshaped());

    phi.block(0, 1, nx - 1, ny - 2) = phi_vec.reshaped(nx - 1, ny - 2).eval();
    phi.row(nx - 1) = phi.row(0);

    //Calculate E
    #pragma omp parallel for
    for(int x_idx = 0; x_idx < nx; x_idx++)
    {
        for(int y_idx = 0; y_idx < ny; y_idx++)
        {
            int x_idx_p = (x_idx == nx - 1) ? 1 : x_idx + 1;
            int x_idx_m = (x_idx == 0) ? nx - 2 : x_idx - 1;
            int y_idx_p = (y_idx == ny - 1) ? y_idx : y_idx + 1;
            int y_idx_m = (y_idx == 0) ? y_idx : y_idx - 1;
            Ex(x_idx, y_idx) = (phi(x_idx_m, y_idx) - phi(x_idx_p, y_idx)) / 2.0 / dx;
            Ey(x_idx, y_idx) = (phi(x_idx, y_idx_m) - phi(x_idx, y_idx_p)) / 2.0 / dy;
        }
    }
}
