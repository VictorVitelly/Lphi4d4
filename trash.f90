subroutine cluster(phi)
  real(dp), dimension(N,N,N,N), intent(inout) :: phi
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
    bond_t(:,:,:,:)=(.false.)

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
      parent(i) = i
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
      end do
    end do
  end do

  !--- Flatten all labels to their root ---
  do i1=1,L  
    do i2=1,L  
      do i3=1,L  
        do i4=1,L
          label(i1,i2,i3,i4)=find(label(i1,i2,i3,i4),parent)
        end do;  
      end do;  
    end do;  
  end do

  !--- Randomly decide whether to flip each cluster ---
  allocate(flip_cluster(next_label))
  flip_cluster(:) = .false.
  do i1 = 1, next_label - 1
    call random_number(r)
    flip_cluster(i) = (r < 0.5_dp)
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

  subroutine clusterold(phi) 
    real(dp), dimension(L,L,L,L),intent(inout) :: phi
    integer(i32), dimension(L,L,L,L) :: spin
    logical, dimension(L,L,L,L) :: bond_x,bond_y,bond_z,bond_t
    integer(i32) :: i1,i2,i3,i4,label(L,L,L,L),parent(L**4)
    integer(i32) :: next_label,x_label,y_label,z_label,t_label
    logical, allocatable :: flip_cluster(:)
    real(dp) :: beta,r,p
    
    spin(:,:,:,:)=nint(sign(1._dp,phi(:,:,:,:)),i32)
    bond_x(:,:,:,:)=(.false.)
    bond_y(:,:,:,:)=(.false.)
    bond_z(:,:,:,:)=(.false.)
    bond_t(:,:,:,:)=(.false.)
    
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
              bond_y(i1,i2)=(r<p)
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
              bond_t(i1,i2,i3,i4)=(r<p)
            end if
          end do
        end do
      end do
    end do

    label(:,:,:,:)=0
    do i=1,L**4
      parent(i)=i
    end do
    next_label=1

    do i1=1,L
      do i2=1,L
        do i3=1,L
          do i4=1,L
            x_label=0
            y_label=0
            z_label=0
            t_label=0
            if(i1>1 .and. bond_x(i1-1,i2,i3,i4) ) then
              x_label=label(i1-1,i2,i3,i4)
            end if
            if(i2>1 .and. bond_y(i1,i2-1,i3,i4) ) then
              y_label=label(i1,i2-1,i3,i4)
            end if
            if(i3>1 .and. bond_z(i1,i2,i3-1,i4) ) then
              z_label=label(i1,i2,i3-1,i4)
            end if
            if(i4>1 .and. bond_t(i1,i2,i3,i4-1) ) then
              t_label=label(i1,i2,i3,i4-1)
            end if
            
            if(left_label==0 .and. up_label==0 .and. front_label=0 .and. future_label=0) then
              label(i1,i2,i3,i4)=next_label
              next_label=next_label+1  
            else if(left_label/=0 .and. up_label==0 .and. front_label==0 .and. future_label==0) then
              label(i1,i2,i3,i4)=left_label
            else if(left_label==0 .and. up_label/=0 .and. front_label==0 .and. future_label==0) then
              label(i1,i2,i3,i4)=up_label
            else if(left_label==0 .and. up_label==0 .and. front_label/=0 .and. future_label==0) then
              label(i1,i2,i3,i4)=front_label
            else if(left_label==0 .and. up_label==0 .and. front_label==0 .and. future_label/=0) then
              label(i1,i2,i3,i4)=future_label
            else
              label(i1,i2)=min(left_label,up_label)
              call union(left_label,up_label,parent)
            end if
            
      end do
    end do
    
    do i2=1,L
      if(bond_x(L,i2) ) then
        call union(label(1,i2),label(L,i2),parent )
      end if
    end do
    do i1=1,L
      if(bond_y(i1,L) ) then
        call union(label(i1,1),label(i1,L),parent )
      end if 
    end do
    do i1=1,L
      do i2=1,L
        label(i1,i2)=find(label(i1,i2),parent)
      end do
    end do

    allocate(flip_cluster(next_label) )
    flip_cluster(:)=.false.
    do i=1,next_label-1
      call random_number(r)
      flip_cluster(i)=(r<0.5_dp)
    end do
    
    do i1=1,L
      do i2=1,L
        if(flip_cluster(label(i1,i2))) then
          phi(i1,i2)=-phi(i1,i2)
        end if
      end do
    end do
    deallocate(flip_cluster)
    
  end subroutine clusterold
