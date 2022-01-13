
#include "particle_transport.h"

void solve_particles_decoupled_c(double *initial_concentration, double *inlet_concentration, double *particle_size);

void solve_particles_decoupled(double initial_concentration, double inlet_concentration, double particle_size)
{
  solve_particles_decoupled_c(&initial_concentration, &inlet_concentration, &particle_size);
}

