program main

  use iso_fortran_env, only : dp => real64, i32 => int32
  use parameters
  use arrays
  use functions
  use statistics
  use measurements
  implicit none

  call read_input()
  call cpu_time(starting)
  call init_vecs()
  !call thermalize(-2._dp,1)
  !call acceptance_rate(-3._dp,-0.5_dp,11)
  call vary_m0(-4._dp,2._dp,13)

  call cpu_time(ending)
  write(*,*) "Elapsed time: ", (ending-starting), " s"

end program main
