//class plasma realization
#include"plasma.h"
#include"input.h"
#include"diagnose.h"
#include"partition2d.h"
#include<iostream>
#include<iomanip>

using namespace std;

void PlasmaSystem::CalculateE()
{
    double tempEk = 0.0;
    int particles_tot_num = 0;
    for (auto particles_a : species)
    {
        #pragma omp parallel for reduction(+:tempEk)
        for(auto prv : particles_a.rv)
        {
            tempEk += 0.5 * particles_a.m * (prv.vx * prv.vx + prv.vy * prv.vy);
        }
        particles_tot_num += particles_a.num;
    }
    tempEk /= particles_tot_num;
    Ek.push_back(tempEk);

    double tempEp = 0.0;
    #pragma omp parallel for reduction(+:tempEp) collapse(2)
    for(int i = 0; i < nx_grids; i++)
    {
        for(int j = 0; j < ny_grids; j++)
            tempEp += 0.5 * (pow(poisson_solver.GetEx(i, j), 2) + pow(poisson_solver.GetEy(i, j), 2)) * dx * dy;
    }
    tempEp /= particles_tot_num;
    Ep.push_back(tempEp);
    Et.push_back(tempEk + tempEp);
}
PlasmaSystem::PlasmaSystem():
    poisson_solver(nx, dx, ny, dy),
    B(0, 0, 0)
{
    Ek.clear();
    Ep.clear();
    Et.clear();
    charge.resize(nx, ny);
}

void PlasmaSystem::Run()
{
    Initialize();
    PrintParameters();
    PrintSpecialInformation();

    //main loop
    for(int n = 0; n < maxsteps + 1; n++)
    {
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
                string filenameVX = data_path + particles_a.name + "_vx_data" + p;
                string filenameVY = data_path + particles_a.name + "_vy_data" + p;
                string filenameX = data_path + particles_a.name + "_x_data" + p;
                string filenameY = data_path + particles_a.name + "_y_data" + p;
                //output particles x v
                //OutputData(filenameVX, GetParticlesVX(particles_a));
                //OutputData(filenameVY, GetParticlesVY(particles_a));
                //OutputData(filenameX,  GetParticlesX(particles_a));
                //OutputData(filenameY,  GetParticlesY(particles_a));
            }
        }

        //calculate E
        charge.setZero();
        SetupSpeciesChargeOnGrids();
        SetupBackgroundChargeOnGrids();
        poisson_solver.Solve(charge);

        //Push One Step
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
                    fex += particles_a.q * partition.contrib[i] * poisson_solver.GetEx(partition.x_idx[i], partition.y_idx[i]);
                    fey += particles_a.q * partition.contrib[i] * poisson_solver.GetEy(partition.x_idx[i], partition.y_idx[i]);
                }

                if(n == 0)
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

        //diagnose energy
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
    cout << "  2D ES PIC Simulation Start!" << endl;
    cout << "--------------------------------------------" << endl;
    cout << "  Simulation Parameters:" << endl;
    cout << "    Lx = " << left << setw(7) << setprecision(4)  << Lx
         << "    kx = " << left << setw(7) << kx
         << "    nx = " << left << setw(7) << nx
         << "    dx = " << left << setw(7) << setprecision(4) << dx << endl;

    cout << "    Ly = " << left << setw(7) << setprecision(4) << Ly
         << "    ky = " << left << setw(7) << ky
         << "    ny = " << left << setw(7) << ny
         << "    dy = " << left << setw(7) << setprecision(4) << dy << endl;

    cout << " Steps = " << left << setw(7) << maxsteps
         << "    dt = " << left << setw(7)  << dt
         << "  Time = " << left << setw(7)  << maxsteps*dt << endl;

    cout << "--------------------------------------------" << endl;
    if(dx > lambda_D)
    {
        cout << "WARNING: DebyeL is NOT satisfied!" << endl;
        cout << "--------------------------------------------" << endl;
    }
    for(auto particles_a : species)
    {
        cout << "  Species: " << particles_a.name << endl;
        cout << "     m = " << left << setw(7) << particles_a.m
             << "N/Cell = " << left << setw(7) << particles_a.num / nx_grids / ny_grids
             << "     N = " << left << setw(7) << particles_a.num
             << "     q = " << left << setw(7) << setprecision(4) << particles_a.q << endl;
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
    MatrixXd temp_rho(nx, ny);
    for(int i = 0; i < nx; i++)
    {
        double x = i * dx;
        for(int j = 0; j < ny; j++)
        {
            double y = j * dy;
            temp_rho(i, j) = GetBackgroundDensity(x, y);
            normalization += temp_rho(i, j);
        }
    }
    charge -= net_charge * temp_rho / dx / dy / normalization;
}

void PlasmaSystem::SetupSpeciesChargeOnGrids()
{
    //CIC 2-order interpolation
    for(auto p : species)
    {
        for(int i = 0; i < p.num; i++)
        {
            PartitionToGrids partition(dx, dy, p.rv[i].x, p.rv[i].y); //linear interpolation //ec scheme
            for(int i = 0; i < 4; i++)
            {
                charge(partition.x_idx[i], partition.y_idx[i]) += p.q * partition.contrib[i] / dx / dy; // charge is density of charge
            }
        }
    }
    charge.col(0) += charge.col(ny - 1);
    charge.col(ny - 1) = charge.col(0);
    charge.row(0) += charge.row(nx - 1);
    charge.row(nx - 1) = charge.row(0);
}
