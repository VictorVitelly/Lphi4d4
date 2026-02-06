module parameters

    use iso_fortran_env, only : dp => real64, i32 => int32
    implicit none
    integer(i32) :: L,thermalization,eachsweep
    real(dp) :: lamb0

    integer(i32) :: Nmsrs=30, Nmsrs2=120
    real(dp) :: dphi=0.4_dp, hotphi=1._dp
    integer(i32), parameter :: Mbin(5)=(/4,5,10,15,20/)

    character(100) :: input_file
    real :: starting,ending

    namelist /input_parameters/ L, thermalization, eachsweep, lamb0

  contains

  subroutine read_input
    integer(i32) :: unit
    write(*,*) "Dame input file"
    read(*,*) input_file
    open(newunit = unit, file = input_file)
    read(unit, nml = input_parameters)
    close(unit)
    write(*, nml = input_parameters)
  end subroutine read_input

end module parameters
