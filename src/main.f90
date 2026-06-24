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
  !call acceptance_rate(-1.0_dp,0.0_dp,21)
  !call autocorrelation(m02,300,4)
  call vary_m04(-1.4_dp,-0.6_dp,11)


  call cpu_time(ending)
  write(*,*) "Elapsed time: ", (ending-starting), " s"

end program main
