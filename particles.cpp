#include "particles.h"
#include<gsl/gsl_rng.h>
#include<sys/timeb.h>
#include<cmath>
#include<iostream>

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
    //set random number
    struct timeb time_seed;
    ftime(&time_seed);
    gsl_rng_default_seed = (time_seed.time * 1000 + time_seed.millitm);
    gsl_rng *r;
    r = gsl_rng_alloc(gsl_rng_default);

    //initialize particles
    double max_probability_density = Distribution(0.0, 0.0, 0.0, 0.0);
    for(int i = 0; i < num; i++)
    {
        double temp_vx = gsl_rng_uniform(r) * 2 * v_max - v_max;
        double temp_vy = gsl_rng_uniform(r) * 2 * v_max - v_max;
        double temp_x = gsl_rng_uniform(r) * Lx;
        double temp_y = gsl_rng_uniform(r) * Ly;

        while(gsl_rng_uniform(r) * max_probability_density > Distribution(temp_x, temp_y, temp_vx, temp_vy))
        {
            temp_x = gsl_rng_uniform(r) * Lx;
            temp_y = gsl_rng_uniform(r) * Ly;
            temp_vx = gsl_rng_uniform(r) * 2 * v_max - v_max;
            temp_vy = gsl_rng_uniform(r) * 2 * v_max - v_max;
        }
        rv[i].x = temp_x;
        rv[i].y = temp_y;
        rv[i].vx = temp_vx;
        rv[i].vy = temp_vy;
    }
}
