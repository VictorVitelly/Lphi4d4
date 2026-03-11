module statistics
  use iso_fortran_env, only : dp => real64, i32 => int32
  use parameters
  use arrays
  use functions
  implicit none

contains

  subroutine random_phi(x,bound)
    implicit none
    real(dp),intent(out) :: x
    real(dp), intent(in) :: bound
    real(dp) :: y
    call random_number(y)
    x = 2._dp*bound*y -bound
  end subroutine random_phi

  subroutine cold_start(phi)
  real(dp) :: phi(:,:,:,:)
    phi=0.0_dp
  end subroutine cold_start

  subroutine hot_start(phi,hotphi)
    real(dp) :: phi(:,:,:,:)
    real(dp), intent(in) :: hotphi
    integer(i32) :: i1,i2,i3,i4
    do i1=1,L
      do i2=1,L
        do i3=1,L
          do i4=1,L
            call random_phi(phi(i1,i2,i3,i4),hotphi)
          end do
        end do
      end do
    end do
  end subroutine hot_start

  subroutine montecarlo(m02,dphi,phi,AR)
    real(dp), intent(in) :: m02,dphi
    real(dp), dimension(:,:,:,:), intent(inout) :: phi
    real(dp), intent(out) :: AR
    real(dp) :: deltaphi,phi2,DS,r,p
    integer(i32) :: i1,i2,i3,i4
    AR=0._dp
    do i1=1,L
      do i2=1,L
        do i3=1,L
          do i4=1,L
            call random_phi(deltaphi,dphi)
            phi2=phi(i1,i2,i3,i4)+deltaphi
            DS=DeltaS(m02,phi,i1,i2,i3,i4,phi2)
            if(DS .le. 0._dp) then
                phi(i1,i2,i3,i4)=phi2
                AR=AR+1._dp
            else
                call random_number(r)
                p=Exp(-DS)
                AR=AR+p
                if(r < p ) then
                    phi(i1,i2,i3,i4)=phi2
                end if
            end if
          end do
        end do
      end do
    end do
    AR=AR/real(L**4,dp)
  end subroutine montecarlo

  subroutine metropolis(m02,phi)
    real(dp), intent(in) :: m02
    real(dp), dimension(:,:,:,:), intent(inout) :: phi
    real(dp) :: deltaphi,phi2,DS,r,p
    integer(i32) :: i1,i2,i3,i4
    do i1=1,L
      do i2=1,L
        do i3=1,L
          do i4=1,L
            call random_phi(deltaphi,dphi)
            phi2=phi(i1,i2,i3,i4)+deltaphi
            DS=DeltaS(m02,phi,i1,i2,i3,i4,phi2)
            if(DS .le. 0._dp) then
                phi(i1,i2,i3,i4)=phi2
            else
                call random_number(r)
                p=Exp(-DS)
                if(r < p ) then
                    phi(i1,i2,i3,i4)=phi2
                end if
            end if
          end do
        end do
      end do
    end do
  end subroutine metropolis

  subroutine cycles(m02,phi,montecarlos)
  real(dp), intent(in) :: m02
  integer(i32), intent(in) :: montecarlos
  real(dp), dimension(:,:,:,:), intent(inout) :: phi
  integer(i32) :: i
  do i=1,montecarlos
    call metropolis(m02,phi)
  end do
  end subroutine cycles

  !ERRORS

  subroutine jackknife(x,y,deltay)
    real(dp), dimension(:), intent(in) :: x
    real(dp), intent(in) :: y
    real(dp), intent(out) :: deltay
    real(dp) :: jackk
    real(dp), allocatable :: xmean(:), delta_y(:)
    integer(i32) :: k,Narr,i,j
      Narr=size(x)
      allocate(delta_y(size(Mbin)))
      do j=1,size(Mbin)
        allocate(xmean(Mbin(j)))
        jackk=0._dp
        xmean=0._dp
        do i=1,Mbin(j)
          do k=1,Narr
            if(k .le. (i-1)*Narr/Mbin(j)) then
              xmean(i)=xmean(i)+x(k)
            else if(k > i*Narr/Mbin(j)) then
              xmean(i)=xmean(i)+x(k)
            end if
          end do
          xmean(i)=xmean(i)/(real(Narr,dp) -real(Narr/Mbin(j),dp))
        end do
        do k=1,Mbin(j)
          jackk=jackk+(xmean(k)-y )**2
        end do
        delta_y(j)=Sqrt(real(Mbin(j)-1,dp)*jackk/real(Mbin(j),dp))
        deallocate(xmean)
      end do
      deltay=maxval(delta_y)
  end subroutine jackknife

  subroutine standard_error(x,y,deltay)
    real(dp), dimension(:), intent(in) :: x
    real(dp), intent(in) :: y
    real(dp), intent(out) :: deltay
    real(dp) :: variance
    integer(i32) :: k,Narr
    Narr=size(x)
    deltay=0._dp
    variance=0._dp
    do k=1,Narr
      variance=variance+(x(k) -y)**2
    end do
    variance=variance/real(Narr-1,dp)
    deltay=Sqrt(variance/real(Narr,dp))
  end subroutine standard_error

  subroutine mean_0(x,y)
    real(dp), dimension(:), intent(in) :: x
    real(dp), intent(out) :: y
    integer(i32) :: k,Narr
    Narr=size(x)
    y=0._dp
    do k=1,Narr
      y=y+x(k)
    end do
    y=y/real(Narr,dp)
  end subroutine mean_0

  subroutine mean_scalar(x,y,deltay)
    real(dp), dimension(:), intent(in) :: x
    real(dp), intent(out) :: y,deltay
    call mean_0(x,y)
    call standard_error(x,y,deltay)
    !call jackknife(x,y,deltay)
  end subroutine mean_scalar

end module statistics
