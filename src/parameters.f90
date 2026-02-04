module parameters

    use iso_fortran_env, only : dp => real64, i32 => int32
    implicit none
    integer(i32), parameter :: L=6,thermalization=5000,eachsweep=100
    real(dp), parameter :: lamb0=1.0_dp
    integer(i32) :: Nmsrs=100, Nmsrs2=120
    real(dp) :: dphi=0.4_dp, hotphi=1._dp
    integer(i32), parameter :: Mbin(5)=(/4,5,10,15,20/)

    real :: starting,ending

end module parameters
