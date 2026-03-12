program test
  use iso_fortran_env, only : dp => real64, i32 => int32
  implicit none
    integer(i32), parameter :: L=3
    real(dp) :: phi(L,L,L,L)
    
    call hot_start(phi,0.5_dp)
    write(*,*) phi(:,:,:,:)
    call cluster(phi)
    write(*,*) phi(:,:,:,:)
    

contains

  subroutine random_phi(x,bound)
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
              p=1._dp!-exp(-2._dp*beta )
              call random_number(r)
              bond_x(i1,i2,i3,i4)=(r<p)
            end if
            if(spin(i1,i2,i3,i4)==spin(i1,mod(i2,L)+1,i3,i4) ) then
              beta=abs(phi(i1,i2,i3,i4))*abs(phi(i1,mod(i2,L)+1,i3,i4))
              p=1._dp!-exp(-2._dp*beta )
              call random_number(r)
              bond_y(i1,i2,i3,i4)=(r<p)
            end if
            if(spin(i1,i2,i3,i4)==spin(i1,i2,mod(i3,L)+1,i4) ) then
              beta=abs(phi(i1,i2,i3,i4))*abs(phi(i1,i2,mod(i3,L)+1,i4))
              p=1._dp!-exp(-2._dp*beta )
              call random_number(r)
              bond_z(i1,i2,i3,i4)=(r<p)
            end if
            if(spin(i1,i2,i3,i4)==spin(i1,i2,i3,mod(i4,L)+1) ) then
              beta=abs(phi(i1,i2,i3,i4))*abs(phi(i1,i2,i3,mod(i4,L)+1))
              p=1._dp!-exp(-2._dp*beta )
              call random_number(r)
              bond_w(i1,i2,i3,i4)=(r<p)
            end if
          end do
        end do
      end do
    end do

    do i1=1,L
      do i2=1,L
        do i3=1,L
          do i4=1,L
            write(*,*) 'El sitio',i1,i2,i3,i4, &
            &'tiene enlaces', bond_x(i1,i2,i3,i4),bond_y(i1,i2,i3,i4), &
            &bond_z(i1,i2,i3,i4),bond_w(i1,i2,i3,i4)
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
      !call random_number(r)
      r=0.1_dp
      flip_cluster(i1) = (r < 0.5_dp)
    end do

  !--- Apply flips ---
    do i1=1,L
      do i2=1,L
        do i3=1,L
          do i4=1,L
            write(*,*) 'El sitio',i1,i2,i3,i4, 'tiene label', label(i1,i2,i3,i4)
            if (flip_cluster(label(i1,i2,i3,i4))) then
              phi(i1,i2,i3,i4) = -phi(i1,i2,i3,i4)
            end if
          end do
        end do
      end do
    end do

    deallocate(flip_cluster)

end subroutine cluster
  

end program test
