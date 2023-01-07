#include "particles.h"
#include<chrono>
#include<random>
#include<iostream>
#include<omp.h>

using namespace std;

PhaseSpace::PhaseSpace(double _x, double _y, double _vx, double _vy)
{
    x = _x;
    y = _y;
    vx = _vx;
    vy = _vy;
}
PhaseSpace::~PhaseSpace(void)
{}

Particles::Particles(double _n, double _q, double _m, string _name):
    num(_n), q(_q),  m(_m), name(_name)
{
    rv.resize(num, 0.0);
}
void Particles::InitializeXV_Random(double (*Distribution)(double, double, double, double), double v_max, double Lx, double Ly)
{
    cout << "--------------------------------------------" << endl;
    cout << "  Initializing " << name << "..." << endl;

    //set random number
    //gsl_rng_default_seed = ((unsigned long)(time(NULL)));
    //gsl_rng *r;
    //r = gsl_rng_alloc(gsl_rng_default);
    default_random_engine e;
    uniform_real_distribution<double> uniform_dist(0.0, 1.0);

    //initialize particles
    double max_probability_density = Distribution(0.0, 0.0, 0.0, 0.0);
    #pragma omp parallel firstprivate(e)
    {
        e.seed(chrono::system_clock::now().time_since_epoch().count() + omp_get_thread_num());

        #pragma omp for
        for(int i = 0; i < num; i++)
        {
            //double temp_vx = gsl_rng_uniform(r) * 2 * v_max - v_max;
            //double temp_vy = gsl_rng_uniform(r) * 2 * v_max - v_max;
            //double temp_x = gsl_rng_uniform(r) * Lx;
            //double temp_y = gsl_rng_uniform(r) * Ly;

            //while(gsl_rng_uniform(r) * max_probability_density > Distribution(temp_x, temp_y, temp_vx, temp_vy))
            //{
            //temp_x = gsl_rng_uniform(r) * Lx;
            //temp_y = gsl_rng_uniform(r) * Ly;
            //temp_vx = gsl_rng_uniform(r) * 2 * v_max - v_max;
            //temp_vy = gsl_rng_uniform(r) * 2 * v_max - v_max;
            //}

            double temp_vx = uniform_dist(e) * 2 * v_max - v_max;
            double temp_vy = uniform_dist(e) * 2 * v_max - v_max;
            double temp_x = uniform_dist(e) * Lx;
            double temp_y = uniform_dist(e) * Ly;

            while(uniform_dist(e) * max_probability_density > Distribution(temp_x, temp_y, temp_vx, temp_vy))
            {
                temp_x = uniform_dist(e) * Lx;
                temp_y = uniform_dist(e) * Ly;
                temp_vx = uniform_dist(e) * 2 * v_max - v_max;
                temp_vy = uniform_dist(e) * 2 * v_max - v_max;
            }

            rv[i].x = temp_x;
            rv[i].y = temp_y;
            rv[i].vx = temp_vx;
            rv[i].vy = temp_vy;
        }
    }
    //gsl_rng_free(r);
    cout << "  Finish!" << endl;
}
