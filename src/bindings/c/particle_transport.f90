module particle_transport_c
  implicit none
  private

contains

!!!###################################################################################

  subroutine solve_particles_decoupled_c(initial_concentration,inlet_concentration,particle_size) &
    bind(C, name="solve_particles_decoupled_c")
    use arrays, only: dp
    use particle_transport, only: solve_particles_decoupled
    implicit none
    
    real(dp),intent(in) :: initial_concentration, inlet_concentration, particle_size

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_solve_particles_decoupled(initial_concentration,inlet_concentration,particle_size)
#else
    call solve_particles_decoupled(initial_concentration, inlet_concentration,particle_size)
#endif

  end subroutine solve_particles_decoupled_c


!###################################################################################
end module particle_transport_c
