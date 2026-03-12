module measurements
  use iso_fortran_env, only : dp => real64, i32 => int32
  use parameters
  use arrays
  use functions
  use statistics
  implicit none

contains

  subroutine thermalize(m0,montecarlos)
  real(dp), intent(in) :: m0
  integer(i32),intent(in) :: montecarlos
  integer(i32) :: i
  open(10, file = 'data/history.dat', status = 'replace')
  !call hot_start(phi,hotphi)
  call cold_start(phi)
  do i=1,thermalization
    !write(10,*) i, action(m0,phi)/real(L**4,dp)
    call cycles(m0,phi,montecarlos)
  end do
  close(10)
  end subroutine thermalize

  subroutine acceptance_rate(mi,mf,Nps)
  real(dp), intent(in) :: mi,mf
  integer(i32), intent(in) :: Nps
  integer(i32) :: i0,i,j,k
  real(dp) :: m0,AR,ARR(Nmsrs2),AR_ave,AR_err
  real(dp) :: magnet(Nmsrs2),magnet_ave,magnet_err
  open(20, file = 'data/magnet.dat', status = 'replace')
  do i0=1,Nps
    m0=mi+(mf-mi)*real(i0-1,dp)/real(Nps-1)
    call thermalize(m0,1)
    ARR=0._dp
    magnet=0._dp
    do i=1,Nmsrs2
      do j=1,Nmsrs
        do k=1,eachsweep
          call montecarlo(m0,dphi,phi,AR)
        end do
        ARR(i)=ARR(i)+AR
        magnet(i)=magnet(i)+abs(mean(phi))
      end do
    end do
    ARR(:)=ARR(:)/real(Nmsrs,dp)
    magnet(:)=magnet(:)/real(Nmsrs,dp)
    call mean_scalar(ARR,AR_ave,AR_err)
    call mean_scalar(magnet,magnet_ave,magnet_err)
    write(*,*) m0, AR_ave, AR_err
    write(20,*) m0, magnet_ave/real(L**4,dp), magnet_err/real(L**4,dp)
  end do
  close(20)
  end subroutine acceptance_rate

  subroutine vary_m0(mi,mf,Nps)
  real(dp), intent(in) :: mi,mf
  integer(i32), intent(in) :: Nps
  integer(i32) :: i0,i,j,k
  real(dp) :: m0,measurements
  real(dp) :: magnet(Nmsrs2),magnet_ave,magnet_err
  real(dp) :: action(Nmsrs2),action_ave,action_err
  open(20, file = 'data/magnet.dat', status = 'replace')
  open(30, file = 'data/action.dat', status = 'replace')
  measurements=real(Nmsrs,dp)
  do i0=1,Nps
    write(*,*) i0, 'de', Nps
    m0=mi+(mf-mi)*real(i0-1,dp)/real(Nps-1)
    call thermalize(m0,4)
    magnet=0._dp
    action=0._dp
    do i=1,Nmsrs2
      do j=1,Nmsrs
        do k=1,eachsweep
          call cycles(m0,phi,4)
        end do
        magnet(i)=magnet(i)+abs(mean(phi))
        action(i)=action(i)+S(m0,phi)
      end do
    end do
    magnet(:)=magnet(:)/measurements
    action(:)=action(:)/measurements
    call mean_scalar(magnet,magnet_ave,magnet_err)
    call mean_scalar(action,action_ave,action_err)
    write(20,*) m0, magnet_ave/real(L**4,dp), magnet_err/real(L**4,dp)
    write(30,*) m0, action_ave/real(L**4,dp), action_err/real(L**4,dp)
  end do
  close(20)
  close(30)
  end subroutine vary_m0

end module measurements
