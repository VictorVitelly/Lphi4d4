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

  subroutine cluster(phi)
  real(dp), dimension(L,L,L,L), intent(inout) :: phi
  integer, allocatable :: spin(:,:,:,:)
  logical, allocatable :: bond_x(:,:,:,:),bond_y(:,:,:,:),bond_z(:,:,:,:),bond_w(:,:,:,:)
  integer(i32), allocatable :: label(:,:,:,:),parent(:)
  logical, allocatable :: flip_cluster(:)
  integer(i32) :: i1,i2,i3,i4,nb,next_label,min_label,nb_label(4)
  real(dp) :: beta,r,p
  allocate(spin(L,L,L,L))
  allocate(bond_x(L,L,L,L), bond_y(L,L,L,L))
  allocate(bond_z(L,L,L,L), bond_w(L,L,L,L))
  allocate(label(L,L,L,L))
  allocate(parent(L*L*L*L))

    spin(:,:,:,:) = nint(sign(1._dp, phi(:,:,:,:)), i32)
    bond_x(:,:,:,:)=(.false.)
    bond_y(:,:,:,:)=(.false.)
    bond_z(:,:,:,:)=(.false.)
    bond_w(:,:,:,:)=(.false.)

  !--- Bond formation in all 4 directions (with periodic BCs via mod) ---
    do i1=1,L
      do i2=1,L
        do i3=1,L
          do i4=1,L
            if(spin(i1,i2,i3,i4)==spin(mod(i1,L)+1,i2,i3,i4) ) then
              beta=abs(phi(i1,i2,i3,i4))*abs(phi(mod(i1,L)+1,i2,i3,i4))
              p=1._dp-exp(-2._dp*beta )
              call random_number(r)
              bond_x(i1,i2,i3,i4)=(r<p)
            end if
            if(spin(i1,i2,i3,i4)==spin(i1,mod(i2,L)+1,i3,i4) ) then
              beta=abs(phi(i1,i2,i3,i4))*abs(phi(i1,mod(i2,L)+1,i3,i4))
              p=1._dp-exp(-2._dp*beta )
              call random_number(r)
              bond_y(i1,i2,i3,i4)=(r<p)
            end if
            if(spin(i1,i2,i3,i4)==spin(i1,i2,mod(i3,L)+1,i4) ) then
              beta=abs(phi(i1,i2,i3,i4))*abs(phi(i1,i2,mod(i3,L)+1,i4))
              p=1._dp-exp(-2._dp*beta )
              call random_number(r)
              bond_z(i1,i2,i3,i4)=(r<p)
            end if
            if(spin(i1,i2,i3,i4)==spin(i1,i2,i3,mod(i4,L)+1) ) then
              beta=abs(phi(i1,i2,i3,i4))*abs(phi(i1,i2,i3,mod(i4,L)+1))
              p=1._dp-exp(-2._dp*beta )
              call random_number(r)
              bond_w(i1,i2,i3,i4)=(r<p)
            end if
          end do
        end do
      end do
    end do

  !--- Initialise union-find structure ---
    label(:,:,:,:) = 0
    do i1=1,L*L*L*L
      parent(i1) = i1
    end do
    next_label=1

  !--- Hoshen-Kopelman pass: check the 4 back-neighbours ---
    do i1=1,L
      do i2=1,L
        do i3=1,L
          do i4=1,L
            nb_label(:) = 0
            if(i1>1 .and. bond_x(i1-1,i2,i3,i4) ) then
              nb_label(1)=label(i1-1,i2,i3,i4)
            end if
            if(i2>1 .and. bond_y(i1,i2-1,i3,i4) ) then
              nb_label(2)=label(i1,i2-1,i3,i4)
            end if
            if(i3>1 .and. bond_z(i1,i2,i3-1,i4) ) then
              nb_label(3)=label(i1,i2,i3-1,i4)
            end if
            if(i4>1 .and. bond_w(i1,i2,i3,i4-1) ) then
              nb_label(4)=label(i1,i2,i3,i4-1)
            end if

          ! Find minimum non-zero neighbour label
            min_label = 0
            do nb=1,4
              if (nb_label(nb) /= 0) then
                if (min_label==0) then
                  min_label = nb_label(nb)
                else
                  min_label = min(min_label,nb_label(nb))
                end if
              end if
            end do

            if (min_label==0) then
              ! No bonded back-neighbours: assign a new label
              label(i1,i2,i3,i4)=next_label
              next_label=next_label+1
            else
            ! Assign minimum and union all differing back-neighbour labels
              label(i1,i2,i3,i4)=min_label
              do nb=1,4
                if (nb_label(nb) /= 0 .and. nb_label(nb) /= min_label) then
                  call union(min_label,nb_label(nb),parent)
                end if
              end do
            end if
          end do
        end do
      end do
    end do

  !--- Periodic boundary unions (wrap-around edges) ---
  ! x wrap: face i=L bonded to face i=1
    do i2=1,L
      do i3=1,L
        do i4=1,L
          if (bond_x(L,i2,i3,i4)) then
            call union(label(1,i2,i3,i4),label(L,i2,i3,i4),parent)
          end if
        end do
      end do
    end do

  ! y wrap
    do i1=1,L
      do i3=1,L
        do i4=1,L
          if (bond_y(i1,L,i3,i4)) then
            call union(label(i1,1,i3,i4),label(i1,L,i3,i4),parent)
          end if
        end do
      end do
    end do

  ! z wrap
    do i1=1,L
      do i2=1,L
        do i4=1,L
          if (bond_z(i1,i2,L,i4)) then
            call union(label(i1,i2,1,i4),label(i1,i2,L,i4),parent)
          end if
        end do
      end do
    end do

  ! w wrap
    do i1=1,L
      do i2=1,L
        do i3=1,L
          if (bond_w(i1,i2,i3,L)) then
            call union(label(i1,i2,i3,1),label(i1,i2,i3,L),parent)
          end if
        end do
      end do
    end do

  !--- Flatten all labels to their root ---
    do i1=1,L
      do i2=1,L
        do i3=1,L
          do i4=1,L
            label(i1,i2,i3,i4)=find(label(i1,i2,i3,i4),parent)
          end do
        end do
      end do
    end do

  !--- Randomly decide whether to flip each cluster ---
    allocate(flip_cluster(next_label))
    flip_cluster(:) = .false.
    do i1 = 1, next_label - 1
      call random_number(r)
      flip_cluster(i1) = (r < 0.5_dp)
    end do

  !--- Apply flips ---
    do i1=1,L
      do i2=1,L
        do i3=1,L
          do i4=1,L
            if (flip_cluster(label(i1,i2,i3,i4))) then
              phi(i1,i2,i3,i4) = -phi(i1,i2,i3,i4)
            end if
          end do
        end do
      end do
    end do

    deallocate(flip_cluster)

  end subroutine cluster


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
