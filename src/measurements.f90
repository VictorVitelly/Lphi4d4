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
  !open(10, file = 'data/history.dat', status = 'replace')
  !call hot_start(phi,hotphi)
  call cold_start(phi)
  do i=1,thermalization
    !write(10,*) i, action(m0,phi)/real(L**4,dp)
    call cycles(m0,phi,montecarlos)
  end do
  !close(10)
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
  real(dp) :: m0,norm,vol
  real(dp) :: magnet(Nmsrs2),magnet_ave,magnet_err
  real(dp) :: action(Nmsrs2),action_ave,action_err
  real(dp), allocatable :: corr2(:),CF(:,:),CF_ave(:,:),CF_err(:,:)
  real(dp) :: xi2_ave,xi2_err
  open(20, file = 'data/magnet.dat', status = 'replace')
  open(30, file = 'data/action.dat', status = 'replace')
  open(60, file = 'data/corrfunc.dat', status = 'replace')
  open(70, file = 'data/corrlen2.dat', status = 'replace')
  allocate(corr2(L))
  allocate(CF(L,Nmsrs2))
  allocate(CF_ave(L,Nps))
  allocate(CF_err(L,Nps))
  vol=real(L**4,dp)
  norm=real(Nmsrs,dp)
  do i0=1,Nps
    CF(:,:)=0._dp
    write(*,*) i0, 'de', Nps
    m0=mi+(mf-mi)*real(i0-1,dp)/real(Nps-1)
    write(*,*) m0
    call thermalize(m0,4)
    magnet=0._dp
    action=0._dp
    do i=1,Nmsrs2
      corr2(:)=0._dp
      do j=1,Nmsrs
        do k=1,eachsweep
          call cycles(m0,phi,4)
        end do
        magnet(i)=magnet(i)+abs(mean(phi))
        action(i)=action(i)+S(m0,phi)
        call correlation(phi,corr2)
      end do
      CF(:,i)=corr2(:)/(norm*vol)
    end do
    magnet(:)=magnet(:)/norm
    action(:)=action(:)/norm
    call mean_scalar(magnet,magnet_ave,magnet_err)
    call mean_scalar(action,action_ave,action_err)
    call secondmomentum(CF,xi2_ave,xi2_err)
    do i=1,L
      call mean_scalar(CF(i,:)-magnet_ave/vol,CF_ave(i,i0),CF_err(i,i0))
    end do
    write(20,*) m0, magnet_ave/vol, magnet_err/vol
    write(30,*) m0, action_ave/vol, action_err/vol
    write(70,*) m0, xi2_ave, xi2_err
  end do
  do k=1,L
    write(60,*) abs(k-1), CF_ave(k,:), CF_err(k,:)
  end do
  close(20)
  close(30)
  close(60)
  close(70)
  deallocate(corr2,CF,CF_ave,CF_err)
  end subroutine vary_m0
  
  subroutine fixed_m0(mi,mf,Nps,kkkk)
  real(dp), intent(in) :: mi,mf
  integer(i32), intent(in) :: Nps,kkkk
  integer(i32) :: i0,i,j,k
  real(dp) :: m0,norm,vol,MM,SS,M2,M4,S2
  real(dp),dimension(Nmsrs2) :: magnet,action,sus,hea,U4 
  real(dp) :: magnet_ave,magnet_err,action_ave,action_err
  real(dp) :: sus_ave,sus_err,hea_ave,hea_err,U4_ave,U4_err
  real(dp), allocatable :: corr2(:),CF(:,:),CF_ave(:),CF_err(:)
  real(dp) :: xi2_ave,xi2_err
  character(len=32) f1,f2,f3,f4,f5,f6,f7
  write(f1, '("data/mag", I0, ".dat")') kkkk
  open(20, file = f1, status = 'replace')
  write(f2, '("data/ene", I0, ".dat")') kkkk
  open(30, file = f2, status = 'replace')
  write(f3, '("data/sus", I0, ".dat")') kkkk
  open(40, file = f3, status = 'replace')
  write(f4, '("data/hea", I0, ".dat")') kkkk
  open(50, file = f4, status = 'replace')
  write(f5, '("data/bin", I0, ".dat")') kkkk
  open(60, file = f5, status = 'replace')  
  write(f6, '("data/cfs", I0, ".dat")') kkkk
  open(70, file = f6, status = 'replace')
  write(f7, '("data/xi2", I0, ".dat")') kkkk
  open(80, file = f7, status = 'replace')
  allocate(corr2(L))
  allocate(CF(L,Nmsrs2))
  allocate(CF_ave(L))
  allocate(CF_err(L))
  vol=real(L**4,dp)
  norm=real(Nmsrs,dp)

  CF(:,:)=0._dp
  m0=mi+(mf-mi)*real(kkkk-1,dp)/real(Nps-1)
  call thermalize(m0,4)
  magnet=0._dp
  action=0._dp
  do i=1,Nmsrs2
    corr2(:)=0._dp
    do j=1,Nmsrs
      do k=1,eachsweep
        call cycles(m0,phi,4)
      end do
      MM=mean(phi)
      SS=S(m0,phi)
      magnet(i)=magnet(i)+abs(MM)
      action(i)=action(i)+SS
      S2=S2+SS**2
      M2=M2+MM**2
      M4=M4+MM**4
      call correlation2(phi,corr2)
    end do
    magnet(i)=magnet(i)/norm
    action(i)=action(i)/norm
    CF(:,i)=corr2(:)/(norm*vol)
    sus(i)=M2-magnet(i)**2 
    hea(i)=S2-action(i)**2
    U4(i)=1._dp-M4/(3._dp*M2**2)
  end do
  call mean_scalar(magnet,magnet_ave,magnet_err)
  call mean_scalar(action,action_ave,action_err)
  call mean_scalar(sus,sus_ave,sus_err)
  call mean_scalar(hea,hea_ave,hea_err)
  call mean_scalar(U4,U4_ave,U4_err)
  call secondmomentum2(CF,xi2_ave,xi2_err)
  do i=1,L
    call mean_scalar(CF(i,:),CF_ave(i),CF_err(i))
  end do
  
  write(20,*) m0, magnet_ave/vol, magnet_err/vol
  write(30,*) m0, action_ave/vol, action_err/vol
  write(40,*) m0, sus_ave/vol, sus_err/vol
  write(50,*) m0, hea_ave/vol, hea_err/vol
  write(60,*) m0, U4_ave/vol, U4_err/vol
  do k=1,L
    write(70,*) abs(k-1), CF_ave(k), CF_err(k)
  end do
  write(80,*) m0, xi2_ave, xi2_err
  
  close(20)
  close(30)
  close(40)
  close(50)
  close(60)
  close(70)
  close(80)
  deallocate(corr2,CF,CF_ave,CF_err)
  end subroutine fixed_m0
  
  subroutine fixed_m0_configs(m0)
  real(dp), intent(in) :: m0
  integer(i32) :: i1,i2,i3
  character(len=8)  :: cL, clambda, cm0
  character(len=64) :: filename
  write(cL,      '(I0)')    L
  write(clambda, '(F6.3)')  lamb0
  write(cm0,     '(F7.3)')  m0

  filename = 'data/chain_' // trim(adjustl(cL))      // &
                      '_' // trim(adjustl(clambda)) // &
                      '_' // trim(adjustl(cm0))     // '.dat'
  open(10, file = filename, form='unformatted', action='write')
    call thermalize(m0,4)
    do i1=1,Nmsrs
      do i2=1,Nmsrs2
        do i3=1,eachsweep
          call cycles(m0,phi,4)
        end do
        write(10,*) phi
      end do
    end do  
  close(10)  
  end subroutine fixed_m0_configs
  
  subroutine autocorrelation(m0,tmax,montecarlos)
    integer(i32), intent(in) :: tmax,montecarlos
    real(dp), intent(in) :: m0
    real(dp), dimension(tmax+1) :: auto,auto_delta,autob,autob_delta
    real(dp) :: vol,norm1,norm2
    real(dp) :: M(Nmsrs+tmax),E(Nmsrs+tmax), auto1(Nmsrs),auto1b(Nmsrs)
    real(dp) :: M_ave,E_ave,auto1_ave,auto1b_ave
    real(dp) :: autoj(tmax+1,Nmsrs2),autojb(tmax+1,Nmsrs2)
    integer(i32) :: i,j,tt
    character(len=8)  :: cL, cm0
    character(len=64) :: filename1,filename2
    integer(i32) :: f1,f2
    write(cL,      '(I0)')    L
    write(cm0,     '(F7.2)')  m0

    filename1 = 'data/chainM_' // trim(adjustl(cL))      // &
                      '_' // trim(adjustl(cm0))     // '.dat'
    open(newunit=f1, file=filename1, form='formatted', action='write')
    
    filename2 = 'data/chainE_' // trim(adjustl(cL))      // &
                      '_' // trim(adjustl(cm0))     // '.dat'
    open(newunit=f2, file=filename2, form='formatted', action='write')
    print*, 'Files created at:', filename1,filename2
    
    vol=real(L**4,dp)
    call thermalize(m0,montecarlos)
    do j=1,Nmsrs2
      do i =1,thermalization
        call cycles(m0,phi,montecarlos)
      end do
      do i=1,Nmsrs+tmax
        call cycles(m0,phi,montecarlos)
        E(i)=S(m0,phi)/vol
        M(i)=abs(mean(phi))/vol
      end do
      call mean_0(M,M_ave )
      call mean_0(E,E_ave )
      
      do tt=0,tmax
        do i=1,Nmsrs
          auto1(i)=M(i)*M(i+tt)
          auto1b(i)=E(i)*E(i+tt)
        end do
        call mean_0(auto1,auto1_ave)
        call mean_0(auto1b,auto1b_ave)
        autoj(tt+1,j)=auto1_ave-(M_ave**2)
        autojb(tt+1,j)=auto1b_ave-(E_ave**2)
      end do
    end do
    do tt=0,tmax
      call mean_scalar(autoj(tt+1,:),auto(tt+1),auto_delta(tt+1))
      call mean_scalar(autojb(tt+1,:),autob(tt+1),autob_delta(tt+1))
    end do
    norm1=abs(auto(1))
    auto=auto/norm1
    auto_delta=auto_delta/norm1
    norm2=abs(autob(1))
    autob=autob/norm2
    autob_delta=autob_delta/norm2
    do tt=0,100
      write(f2,*) tt,autob(tt+1),autob_delta(tt+1)
      write(f1,*) tt, auto(tt+1), auto_delta(tt+1)
    end do
    close(f1)
    close(f2)
  end subroutine autocorrelation
  
  
  subroutine time_test(m0,montecarlos)
  integer(i32), intent(in) :: montecarlos
  real(dp), intent(in) :: m0
  integer(i32) :: i,j
  real(dp) :: ti,tf,time_ave,time_delta
  real(dp) :: time(120)
  call cold_start(phi)
  call thermalize(m0,montecarlos)
  do j=1,120
    call cpu_time(ti)
    do i=1,5000
      call cycles(m0,phi,montecarlos)
    end do
    call cpu_time(tf)
    time(j)=tf-ti
  end do
  call mean_scalar(time,time_ave,time_delta)
  write(*,*) montecarlos, m0, time_ave, time_delta
  end subroutine time_test
  
  subroutine vary_m02(mi,mf,Nps)
  real(dp), intent(in) :: mi,mf
  integer(i32), intent(in) :: Nps
  integer(i32) :: i0,i,j,k
  real(dp) :: m0,norm,vol
  real(dp) :: magnet(Nmsrs2),magnet_ave,magnet_err
  real(dp) :: action(Nmsrs2),action_ave,action_err
  real(dp), allocatable :: corr2(:,:),CF(:,:,:),CF_ave(:,:),CF_err(:,:)
  real(dp) :: xi2_ave,xi2_err
  open(20, file = 'data/magnet.dat', status = 'replace')
  open(30, file = 'data/action.dat', status = 'replace')
  open(60, file = 'data/corrfunc.dat', status = 'replace')
  open(70, file = 'data/corrlen2.dat', status = 'replace')
  allocate(corr2(L,L))
  allocate(CF(L,L,Nmsrs2))
  allocate(CF_ave(L,Nps))
  allocate(CF_err(L,Nps))
  vol=real(L**4,dp)
  norm=real(Nmsrs,dp)
  do i0=1,Nps
    CF=0._dp
    write(*,*) i0, 'de', Nps
    m0=mi+(mf-mi)*real(i0-1,dp)/real(Nps-1)
    write(*,*) m0
    call thermalize(m0,4)
    magnet=0._dp
    action=0._dp
    do i=1,Nmsrs2
      corr2=0._dp
      do j=1,Nmsrs
        do k=1,eachsweep
          call cycles(m0,phi,4)
        end do
        magnet(i)=magnet(i)+abs(mean(phi))
        action(i)=action(i)+S(m0,phi)
        call correlation2(phi,corr2)
      end do
      CF(:,:,i)=corr2(:,:)/(norm*vol)
    end do
    magnet(:)=magnet(:)/norm
    action(:)=action(:)/norm
    call mean_scalar(magnet,magnet_ave,magnet_err)
    call mean_scalar(action,action_ave,action_err)
    call secondmomentum2(CF,xi2_ave,xi2_err)
    do i=1,L
      call mean_scalar(CF(i,1,:),CF_ave(i,i0),CF_err(i,i0))
    end do
    write(20,*) m0, magnet_ave/vol, magnet_err/vol
    write(30,*) m0, action_ave/vol, action_err/vol
    write(70,*) m0, xi2_ave, xi2_err
  end do
  do k=1,L
    write(60,*) abs(k-1), CF_ave(k,:), CF_err(k,:)
  end do
  close(20)
  close(30)
  close(60)
  close(70)
  deallocate(corr2,CF,CF_ave,CF_err)
  end subroutine vary_m02

  subroutine vary_m04(mi,mf,Nps)
  real(dp), intent(in) :: mi,mf
  integer(i32), intent(in) :: Nps
  integer(i32) :: i0,i,j,k
  real(dp) :: m0,norm,vol
  real(dp) :: magnet(Nmsrs2),magnet_ave,magnet_err
  real(dp) :: action(Nmsrs2),action_ave,action_err
  real(dp), allocatable :: corr2(:,:,:,:),CF(:,:,:,:,:),CF_ave(:,:),CF_err(:,:)
  real(dp) :: xi2_ave,xi2_err
  open(20, file = 'data/magnet.dat', status = 'replace')
  open(30, file = 'data/action.dat', status = 'replace')
  open(60, file = 'data/corrfunc.dat', status = 'replace')
  open(70, file = 'data/corrlen2.dat', status = 'replace')
  allocate(corr2(L,L,L,L))
  allocate(CF(L,L,L,L,Nmsrs2))
  allocate(CF_ave(L,Nps))
  allocate(CF_err(L,Nps))
  vol=real(L**4,dp)
  norm=real(Nmsrs,dp)
  do i0=1,Nps
    CF=0._dp
    write(*,*) i0, 'de', Nps
    m0=mi+(mf-mi)*real(i0-1,dp)/real(Nps-1)
    write(*,*) m0
    call thermalize(m0,4)
    magnet=0._dp
    action=0._dp
    do i=1,Nmsrs2
      corr2=0._dp
      do j=1,Nmsrs
        do k=1,eachsweep
          call cycles(m0,phi,4)
        end do
        magnet(i)=magnet(i)+abs(mean(phi))
        action(i)=action(i)+S(m0,phi)
        call correlation4(phi,corr2)
      end do
      CF(:,:,:,:,i)=corr2(:,:,:,:)/(norm*vol)
    end do
    magnet(:)=magnet(:)/norm
    action(:)=action(:)/norm
    call mean_scalar(magnet,magnet_ave,magnet_err)
    call mean_scalar(action,action_ave,action_err)
    call secondmomentum4(CF,xi2_ave,xi2_err)
    do i=1,L
      call mean_scalar(CF(i,1,1,1,:),CF_ave(i,i0),CF_err(i,i0))
    end do
    write(20,*) m0, magnet_ave/vol, magnet_err/vol
    write(30,*) m0, action_ave/vol, action_err/vol
    write(70,*) m0, xi2_ave, xi2_err
  end do
  do k=1,L
    write(60,*) abs(k-1), CF_ave(k,:), CF_err(k,:)
  end do
  close(20)
  close(30)
  close(60)
  close(70)
  deallocate(corr2,CF,CF_ave,CF_err)
  end subroutine vary_m04
end module measurements
