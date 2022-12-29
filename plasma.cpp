//class plasma realization
#include"plasma.h"
#include"input.h"
#include"diagnose.h"
#include"partition2d.h"
#include<iostream>
#include<iomanip>
#include<fstream>

using namespace std;

void PlasmaSystem::CalculateE()
{
    double tempEk = 0.0;
    int particles_tot_num = 0;
    for (auto particles_a : species)
    {
        for(auto prv : particles_a.rv)
        {
            tempEk += 0.5 * particles_a.m * (prv.vx * prv.vx);
        }
        particles_tot_num += particles_a.num;
    }
    tempEk /= particles_tot_num;
    Ek.push_back(tempEk);

    double tempEp = 0.0;
    for(int i = 0; i < nx_grids; i++)
    {
        for(int j = 0; j < ny_grids; j++)
            tempEp += 0.5 * pow(poisson_solver.GetEx(i, j), 2) * dx;
    }
    tempEp /= particles_tot_num;
    Ep.push_back(tempEp);
    Et.push_back(tempEk + tempEp);
}
void PlasmaSystem::CalculateT()
{
    T.clear();
    T.resize(nx_grids, 0.0);
    num_in_grid.clear();
    num_in_grid.resize(nx_grids, 0.0);

    for (auto particles_a : species)
    {
        for(auto prv : particles_a.rv)
        {
            int idx_grid = floor(prv.x / dx);
            if (idx_grid == nx_grids)
                idx_grid = 0;
            T[idx_grid] += particles_a.m * (prv.vx * prv.vx);
            num_in_grid[idx_grid]++;
        }
    }
    for(int i = 0; i < nx_grids; i++)
    {
        T[i] /= num_in_grid[i];
    }
}
PlasmaSystem::PlasmaSystem():
    B(0, 0, 0), poisson_solver(nx_grids, dx, ny_grids, dy, VectorXd::Zero(ny_grids), VectorXd::Zero(nx_grids), VectorXd::Zero(ny_grids), VectorXd::Zero(nx_grids))
{
    Ek.clear();
    Ep.clear();
    Et.clear();
    charge.resize(nx_grids, ny_grids);
}

void PlasmaSystem::PushOneStep(int if_init)
{
    for(auto &particles_a : species)
    {
        #pragma omp parallel for
        for(int j = 0; j < particles_a.num; j++)
        {
            PartitionToGrids partition(dx, dy, particles_a.rv[j].x, particles_a.rv[j].y); //linear interpolation //ec scheme
            double fex = 0.0;
            double fey = 0.0;
            for(int i = 0; i < 4; i++)
            {
                fex += particles_a.q * poisson_solver.GetEx(partition.x_idx[i], partition.y_idx[i]) * partition.contrib[i];
                fey += particles_a.q * poisson_solver.GetEy(partition.x_idx[i], partition.y_idx[i]) * partition.contrib[i];
            }

            if(if_init == 0)
            {
                particles_a.rv[j].vx += 0.5 * fex * dt / particles_a.m;
                particles_a.rv[j].x += particles_a.rv[j].vx * dt;
                particles_a.rv[j].vy += 0.5 * fey * dt / particles_a.m;
                particles_a.rv[j].y += particles_a.rv[j].vy * dt;
            }
            else
            {
                particles_a.rv[j].vx += fex * dt / particles_a.m;
                particles_a.rv[j].x += particles_a.rv[j].vx * dt;
                particles_a.rv[j].vy += fey * dt / particles_a.m;
                particles_a.rv[j].y += particles_a.rv[j].vy * dt;
            }
            //period condition
            while(particles_a.rv[j].x < 0.0)
            {
                particles_a.rv[j].x += Lx;
            }
            while(particles_a.rv[j].x >= Lx)
            {
                particles_a.rv[j].x -= Lx;
            }
            while(particles_a.rv[j].y < 0.0)
            {
                particles_a.rv[j].y += Ly;
            }
            while(particles_a.rv[j].y >= Ly)
            {
                particles_a.rv[j].y -= Ly;
            }
        }
    }
}


void PlasmaSystem::Run()
{
    Initialize();
    PrintParameters();
    PrintSpecialInformation();

    //main loop
    for(int n = 0; n < maxsteps + 1; n++)
    {
        CalculateT();
        int percent = 100 * n / (maxsteps);
        //print running process
        if(percent % 5 == 0)
        {
            cout << "\r" << "  Calculation Process:" << percent << "%" << flush;
        }

        //diagnose plasma every timestep
        if(n % data_steps == 0)
        {
            string p = to_string(n / data_steps);

            for(auto particles_a : species)
            {
                string filenameV = data_path + particles_a.name + "v_data" + p;
                string filenameX = data_path + particles_a.name + "x_data" + p;
                //output particles x v
                OutputData(filenameV, GetParticlesVX(particles_a));
                OutputData(filenameX,  GetParticlesX(particles_a));
                //output temperature distribution
                //OutputData(data_path + "temperature" + p, T);
                //OutputData(data_path + "num_in_grid" + p, num_in_grid);
            }
        }

        //calculate E
        charge.resize(nx_grids,ny_grids);
        charge.setZero();
        SetupSpeciesChargeOnGrids();
        //SetupBackgroundChargeOnGrids();
        poisson_solver.Solve(charge);
        PushOneStep(n);
        CalculateE();
    }
    //output energy evolution
    OutputData(data_path + "/tot_energy", Et);
    OutputData(data_path + "/kin_energy", Ek);
    OutputData(data_path + "/pot_energy", Ep);

    cout << endl << "  Complete!" << endl;
}

void PlasmaSystem::PrintParameters() const
{
    cout << "--------------------------------------------" << endl;
    cout << "  PIC Simulation Start!" << endl;
    cout << "--------------------------------------------" << endl;
    cout << "  Simulation Parameters:" << endl;
    cout << setw(13) << "  Length"
         << setw(13) << "       k"
         << setw(13) << "nx_grids"
         << setw(13) << "Lambda_D"
         << setw(13) << setprecision(6) << "      dx" << endl;

    cout << setw(13) << Lx
         << setw(13) << kx
         << setw(13) << nx_grids
         << setw(13) << lambda_D
         << setw(13) << dx << endl;

    cout << setw(13) << "MaxSteps"
         << setw(13) << "      dt"
         << setw(13) << "    Time" << endl;

    cout << setw(13) << maxsteps
         << setw(13) << dt
         << setw(13) << maxsteps*dt << endl;

    cout << "--------------------------------------------" << endl;
    if(dx > lambda_D)
    {
        cout << "WARNING: DebyeL is NOT satisfied!" << endl;
        cout << "--------------------------------------------" << endl;
    }
    for(auto particles_a : species)
    {
        cout << "  Species: " << particles_a.name << endl;
        cout  << setw(8) << "N = "  << particles_a.num
              << setw(15) << "N/Cell = " << particles_a.num / nx_grids
              << setw(8) << "q = " << setprecision(6) << particles_a.q
              << setw(8) << "m = " << particles_a.m << endl;
    }
    cout << "--------------------------------------------" << endl;
    cout << "  Data: " << endl;
    if(!if_continued)
        cout << "  Evolving from a NEW initial state!" << endl;
    else
        cout << "  Evolving from a CONTINUED state!" << endl;
    cout << "  data_num = " << data_num << endl;
    cout << "--------------------------------------------" << endl;

}

void PlasmaSystem::SetupBackgroundChargeOnGrids()
{
    double net_charge = 0.0;
    for(auto p : species)
    {
        net_charge += p.num * p.q;
    }
    double normalization = 0.0;
    MatrixXd rho(nx_grids, ny_grids);
    for(int i = 0; i < nx_grids; i++)
    {
        double x = i * dx;
        for(int j = 0; j < ny_grids; j++)
        {
            double y = j * dy;
            rho(i, j) = GetBackgroundDensity(x, y);
            normalization += rho(i, j);
        }
    }
    charge -= net_charge * rho / dx / normalization;
    //for(int i = 0; i < nx_grids; i++)
    //{
    //for(int j = 0; j < ny_grids; j++)
    //{
    //rho(i, j) /= normalization;
    //charge(i, j) -= net_charge * rho(i, j) / dx;
    //}
    //}
}

void PlasmaSystem::SetupSpeciesChargeOnGrids()
{
    //CIC 2-order interpolation
    for(auto p : species)
    {
        for(int i = 0; i < p.num; i++)
        {
            //map<int, double> partition_contrib = PartitionToGrid(dx, nx_grids, p.rv[i].x, 2); //interpolation
            PartitionToGrids partition(dx, dy, p.rv[i].x, p.rv[i].y); //linear interpolation //ec scheme
            for(int i = 0; i < 4; i++)
            {
                charge(partition.x_idx[i], partition.y_idx[i]) += p.q * partition.contrib[i] / dx / dy; // / pow(gridWidth, 1);
            }
        }
    }
}
