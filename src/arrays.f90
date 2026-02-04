module arrays
    use iso_fortran_env, only : dp => real64, i32 => int32
    use parameters, only : L
    implicit none

    real(dp), allocatable :: phi(:,:,:,:)
    integer(i32), allocatable :: ip(:), im(:)

contains

  subroutine init_vecs()
  integer(i32) :: i
  allocate(phi(L,L,L,L))
  allocate(ip(L),im(L))
  do i=1,L-1
    ip(i)=i+1
  end do
  ip(L)=1
  do i=2,L
    im(i)=i-1
  end do
  im(1)=L
  end subroutine init_vecs


end module arrays
