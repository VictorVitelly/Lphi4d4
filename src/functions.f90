module functions
    use iso_fortran_env, only : dp => real64, i32 => int32
    use parameters
    use arrays
    implicit none

    contains


  function lagrangian(m02,phi,i1,i2,i3,i4)
    real(dp), intent(in) :: m02
    real(dp), dimension(:,:,:,:), intent(in) :: phi
    integer(i32), intent(in) :: i1,i2,i3,i4
    real(dp) :: lag1,lag2,lag3,lag4,lag5
    real(dp) :: lagrangian
    lagrangian=0._dp
    lag1=(phi(ip(i1),i2,i3,i4)-phi(i1,i2,i3,i4))**2
    lag2=(phi(i1,ip(i2),i3,i4)-phi(i1,i2,i3,i4))**2
    lag3=(phi(i1,i2,ip(i3),i4)-phi(i1,i2,i3,i4))**2
    lag4=(phi(i1,i2,i3,ip(i4))-phi(i1,i2,i3,i4))**2
    lag5=(m02+lamb0*0.5_dp*phi(i1,i2,i3,i4)**2 )*phi(i1,i2,i3,i4)**2
    lagrangian=lag1+lag2+lag3+lag4+lag5
  end function lagrangian

  function S(m02,phi)
    real(dp), intent(in) :: m02
    real(dp), dimension(:,:,:,:), intent(in) :: phi
    integer(i32) :: i1,i2,i3,i4
    real(dp) :: S
    S=0._dp
    do i1=1,L
      do i2=1,L
        do i3=1,L
          do i4=1,L
            S=S+lagrangian(m02,phi,i1,i2,i3,i4)
          end do
        end do
      end do
    end do
    S=0.5_dp*S
  end function S

  function DeltaS(m02,phi,i1,i2,i3,i4,phi2)
    real(dp), intent(in) :: m02
    real(dp), dimension(:,:,:,:), intent(in) :: phi
    integer(i32), intent(in) :: i1,i2,i3,i4
    real(dp), intent(in) :: phi2
    real(dp) :: ds01,ds02,ds03,ds1,ds2,ds3,ds4
    real(dp) :: DeltaS
    ds01=(4._dp+0.5_dp*m02)*(phi2**2-phi(i1,i2,i3,i4)**2)
    ds02=0.25_dp*lamb0*(phi2**4-phi(i1,i2,i3,i4)**4)
    ds1=phi(ip(i1),i2,i3,i4)+phi(im(i1),i2,i3,i4)
    ds2=phi(i1,ip(i2),i3,i4)+phi(i1,im(i2),i3,i4)
    ds3=phi(i1,i2,ip(i3),i4)+phi(i1,i2,im(i3),i4)
    ds4=phi(i1,i2,i3,ip(i4))+phi(i1,i2,i3,im(i4))
    ds03=-(phi2-phi(i1,i2,i3,i4))*(ds1+ds2+ds3+ds4)
    DeltaS=ds01+ds02+ds03
  end function DeltaS

  function mean(phi)
  integer(i32) :: i1,i2,i3,i4
  real(dp) :: phi(:,:,:,:)
  real(dp) :: mean
  mean=0._dp
    do i1=1,L
      do i2=1,L
        do i3=1,L
          do i4=1,L
            mean=mean+phi(i1,i2,i3,i4)
          end do
        end do
      end do
    end do
  end function

  recursive function find(x,parent) result(out)
    integer(i32), intent(in) :: x
    integer(i32), intent(inout) :: parent(:)
    integer(i32) :: out
    if(parent(x) /= x) then
      parent(x)=find(parent(x),parent )
    end if
    out=parent(x)
  end function find

  subroutine union(x,y,parent)
    integer(i32),intent(in) :: x,y
    integer(i32),intent(inout) :: parent(:)
    integer :: root_x, root_y
    root_x=find(x,parent)
    root_y=find(y,parent)
    if (root_x /= root_y) then
      parent(root_y)=root_x
    end if
  end subroutine union


end module functions
