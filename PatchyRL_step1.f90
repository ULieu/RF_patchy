!****************************************************************************
!  PURPOSE: Dynamics of patchy particles described by spherical harmonic Y_55. 
!           The particles freely rotate in 3d while their position is confined to 2d flat space.!           
!****************************************************************************
!  For RL program:
!       - in PatchyPlaneM_RF_irun0.f90, at irun=0, the initial position is random.  
!       - in PatchyPlaneM_RF.f90, the initial position is from last snapshot of previous irun.  

!------------------------------------------------------------------------
    module para
    implicit double precision (a-h,o-z)
    parameter (nmax=5000) !CHANGEABLE as long as nmax > npm
    parameter (maxnab=nmax*25)
    !parameter (pi=3.1415927d0)
    common /pi/pi
    common /npm/npm /dt/dt
    common /dens/dens /bigE/bigE /temp/temp /sigma/sigma /rcut/rcut /radius/radius /bl/bl
    common /Pe/Pe,pbcx,Pe0,omg   
    common /istep/istep /istepmax/istepmax /ioutstep/ioutstep /ioutstep2/ioutstep2 
    
    dimension ::    xp(nmax,3), fp(nmax,3), tp(nmax,3)     
    dimension :: xic(nmax*3)    !thermal fluctuation -translation 
    dimension :: xicr(nmax*3)   !thermal fluctuation -rotation
    
    common /lrank/lrank /mrank/mrank /coe/coe /cmd/cmd /cmr/cmr 
    dimension ::     P0(NMAX,3),P1(NMAX,3),P2(NMAX,3)
    common/amax/amax 
    common/ncode/ncode 
    common /coef_har/coef_har /comp_rate/comp_rate
    common /volfra/volfra /volmax/volmax
    
    common /tinterval/tinterval /tmax/tmax
    common /bige_program/ ibigemax,istep_interval,itime(30),timeMD
    
    common /energyI/energyI(nmax) /aveEnergy/aveEnergyI   
    common /nseed/nseed(2) !generated seed 
    
    !variable for verlet scheme:---
    common /rlist/rlist
    integer point(nmax),list(maxnab)
    logical update
    dimension xp0(nmax,3) !the position in the last update
        
    end module para
    
    !................................    
    module crossproduct
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    contains 
    FUNCTION CROSS(A,B)
    DIMENSION :: CROSS(3)
    DIMENSION :: A(3),B(3)
    cross(1) = a(2) * b(3) - a(3) * b(2)
    cross(2) = a(3) * b(1) - a(1) * b(3)
    cross(3) = a(1) * b(2) - a(2) * b(1)
    END FUNCTION cross 
    end module crossproduct
!-----------------------------------------------------------------------       
    
!   MAIN PROGRAM
    program patchyonplaneTau
    use para
    implicit double precision (a-h,o-z) 
    common /kcr/kcr 
    
    !initialise system_clock
    call SYSTEM_CLOCK(count_rate=kcr)
    call SYSTEM_CLOCK(count_max =kcm)
    call SYSTEM_CLOCK(kcbeg)
    
!---set parameters either by matlab or by hand
    call defprm_file     !if input by matlab   
    
!---set other condition    
    call condition    
    
!---set initial position and direction
    call initia   
    !call initia_file 
    
!---output particle   
        open (53,file='00pos.dat'   ,status='unknown',access='append')          
        open (54,file='00istep.dat' ,status='unknown',access='append')          
                                                                                
        open (30,file='00orient0.dat',status='unknown',access='append')         
        open (31,file='00orient1.dat',status='unknown',access='append')         
        open (32,file='00orient2.dat',status='unknown',access='append')         
        open (33,file='00aveEnergyI.dat',status='unknown',access='append')     
        open (34,file='00EnergyI.dat',status='unknown',access='append')   
        open (35,file='00PBC.dat',status='unknown',access='append')  
        open (36,file='00Pe_t.dat',status='unknown',access='append')   
        
    !just write the INITIAL condition of step
        if (mod(istep,ioutstep)==0) then
        write(*,'(i8,3f10.5)')istep,xp(1,1:3)
        do i=1,npm
          write(53,'(2e13.5)')xp(i,1:2)
          write(54,'(i9)')istep
      
          write(30,'(3e13.5)')p0(i,1:3)  
          write(31,'(3e13.5)')p1(i,1:3)  
          write(32,'(3e13.5)')p2(i,1:3)       
        enddo
         write(35,'(e13.5)')pbcx
         Pe=Pe0*sin(omg*istep*dt)
         write(36,'(i9,e13.5)')istep,Pe
        endif 
    
    call CPU_TIME(t1)
    call SYSTEM_CLOCK(kcini)
    istep=0    
100 istep=istep+1
    !call bigetemp !increase e/kt with timestep
!---calculate force and torque   
    call check         
    call force_verlet  
    
!---calculate correlated normal variates (noise)
    call covar     
              
!---calculate new configuration 
    call calpos   !translation 
    call calorient !orientation   
     
    
!---output for position 
        if (mod(istep,ioutstep)==0) then
        write(*,'(i8,3f10.5)')istep,xp(1,1:3)
          
        do 200 i=1,npm    
          write(53,'(2e13.5)')xp(i,1:2)
          write(54,'(i9)')istep
      
          write(30,'(3e13.5)')p0(i,1:3)  
          write(31,'(3e13.5)')p1(i,1:3)  
          write(32,'(3e13.5)')p2(i,1:3)       
    200 continue
        write(35,'(e13.5)')pbcx
        write(36,'(i9,e13.5)')istep,Pe
        if (xp(1,1) /= xp(1,1)) stop !stop if there is NaN    
        endif    
    
!---output of energy 
        if (mod(istep-1,ioutstep)==0) then   
        write(33,'(i9,x,e13.5)')istep-1,aveEnergyI 
        do ipm=1,npm   
            write(34,'(i9,x,e13.5)')istep-1,EnergyI(ipm) 
        enddo
        endif
     
    if (istep < istepmax) goto 100
    
    call CPU_TIME(t2)
    call SYSTEM_CLOCK(kcend)
    open(60,file='00time.dat',status='unknown',access='append')
    write(60,'(i9,f10.3,i10,2e13.5)')npm,rlist,istepmax,real(kcini-kcbeg)/real(kcr)/3600.d0,real(kcend-kcini)/real(kcr)/3600.d0
    close(60)
    
    
!---final output of energy
        call force
        write(33,'(i9,x,e13.5)')istep,aveEnergyI
        do ipm=1,npm
            write(34,'(i9,x,e13.5)')istep-1,EnergyI(ipm) !201904
        enddo
    
        close(53); close(54); 
        close(30); close(31); close(32); 
        close(33); close(34)
        close(35); close(36)    

        
    stop
    end program patchyonplaneTau
    
!****************************************************************************
 
    
!--------------------------------------    
    subroutine defprm_file
    use para
    implicit double precision (a-h,o-z)
    !set pi value
    pi=acos(-1.d0)
    
    Write(*,*)'READ lrank, mrank, ncode, npm, volfra, bige, dt, istepmax, ioutstep  FROM MATLAB'
    open(10,file='00defprm.dat',status='old')
    read(10,*)lrank, mrank, ncode, npm, volfra, ioutstep,timeMD, ibigemax, istep_interval, tinterval, tmax, rlist
    close(10)  
    write(*,*)lrank, mrank, ncode, npm, volfra, ioutstep,          ibigemax, istep_interval, tinterval, tmax, rlist
    
    open(12,file='00Pe_osci.dat',status='old')
    read(12,*)Pe0,omg
    close(12) 
    
    !timestep
    istep=0
        
    !energy scale initial
     temp= tmax
     bige=1.d0/temp
    
    !IF use fixed dt (case a)
    istepmax = ibigemax * istep_interval    
    dt=1.d-4  ! Now the time of each simulation is dt*istepmax=10
        
    end subroutine defprm_file    
    
!------------------------------------------------    
    subroutine condition
    use para
    implicit double precision (a-h,o-z)
    
    !normalized diameter and radius of particle
      sigma = 2.0d0             
      radius = sigma * 0.5d0 
      rcut = 6.d0*radius  !use in subroutine pair
        pbcx=0.d0
      !volume fraction    
    bL = (2.d0/3.d0*pi*dble(npm)/volfra)**(1.0d0/2.0d0)*radius  !calculation box size  ! ON PLANE
    if (rcut > bl*0.5d0) write(*,*)'caution rcut > box size/2' 
    areafraction=pi*dble(npm)*radius*radius/bl/bl  
    
    !condition for ncode on potential and rcut
          
    if (ncode ==4) then
            coe = 1.d0 !h-h or t-t
            cmd=2.294d0 !wide 
            cmr=1.0d0
            rcut = 6.d0*radius
    endif 
    
    !normalize value for Aaniso    
    if (lrank==5.and. mrank==5)    then  
        amax= 11563.d0*pi/352.d0                                        
    endif
    
    open(85,file='00condition.dat',status='unknown')
    write(85,'(''Number_of_particle     '', i8 )')npm
    write(85,'(''istepmax               '', i8 )')istepmax
    write(85,'(''ioutstep               '', i8 )')ioutstep
    write(85,'(''bige                   '', f10.4 )')bige
    write(85,'(''temp                   '', f10.4 )')temp
  
    write(85,'('' Potential cutoff      '',f10.4 )') rcut
    write(85,'('' Radius                '',f10.4 )') radius
    write(85,'('' Periodic box size     '',f10.4 )') bl
    write(85,'('' Volume volfration     '',f10.5 )') volfra
    write(85,'(''particle_lrank         '', i8 )')lrank
    write(85,'(''particle_mrank         '', i8 )')mrank
    write(85,'(''ncode                  '', i8 )')ncode
    
    write(85,'(''tmax                   '', f10.5 )')tmax
    write(85,'(''tinterval              '', f10.5 )')tinterval
    write(85,'(''ibigemax               '', i8 )')ibigemax
    write(85,'(''istep_interval         '', i8 )')istep_interval
    write(85,'(''timeMD         '', f10.2 )')timeMD
        
    write(85,'(''verlet_rlist(00time)   '', f10.4 )')rlist  
    write(85,'(''Pe                     '', f10.4 )')Pe
    
    close(85)  
    
    end subroutine condition
    
 
!----------------------------------------------------
    subroutine initia
    use para
    use crossproduct
    implicit double precision (a-h,o-z)    
    DIMENSION:: AA(3),BB(3),AXB(3),ei(3)
    common /dxx/dxx
    common /kcr/kcr !use in initia + main routine if the time is too long
    dimension :: ntime(8)

    nseed=0
    !---if turn on/off random_seed auto generated/set (max int32 = 2147483647)
        !call RANDOM_SEED ()  !tunr off for fix chain of numbers
        !call RANDOM_SEED (get=nseed) ! get system generated seed
        open(486,file='00seedparfor.dat',status='old') !matlab  parfor rep2
        read(486,*)nseedparfor
        close(486)
    
        call DATE_AND_TIME()
        call DATE_AND_TIME(values=ntime)
        nseed = ntime(5)*1d7 + ntime(6)*1d5 + ntime(7)*1d3 + ntime(8)*1d0 + & !   hhmmssmmm 
                               ntime(1)*1d5 + ntime(2)*1d3 + ntime(3)*1d1 + & ! + yyyymmdd0  
                               nseedparfor*1d6                                ! +  nn000000 !rep2 in parfor     
        !nseed = 335710487
        call random_seed(put=nseed) ! set current seed
        call random_seed(get=nseed) ! get current seed
    !-----------------------------
    write(*,*)'system generated seed',    nseed   ! writes system generated seed
    open(482,file='00seed.dat',status='unknown')
    write(482,'(11i13)')ntime, nseedparfor, nseed
    close(482)
    
    
!---assign position of particle if volume volfration < 0.35
    rrsq=2.d0*radius*2.d0*radius
    if (volfra < 0.35) then   
       
       do 100 i=1,npm 
   300   call RANDOM_NUMBER (rx)
         call RANDOM_NUMBER (ry)
         !call RANDOM_NUMBER (rz)
         
         xx1=(rx-0.5d0)*bl ! (-bl/2,bl/2) MODIFIED
         yy1=(ry-0.5d0)*bl ! (-bl/2,bl/2) MODIFIED
         !zz1=(rz-0.5d0)*bl ! (-bl/2,bl/2) MODIFIED
         zz1=0.d0
         
         do 200 j=1,i-1
            !do 250 kk=-1,1
            do 250 jj=-1,1
            do 250 ii=-1,1
                xx2=xp(j,1)+ bl*dble(ii)
                yy2=xp(j,2)+ bl*dble(jj)
                !zz2=xp(j,3)+ bl*dble(kk)
                
                delx=xx1-xx2
                dely=yy1-yy2
                delz=0.d0
                
                rr = delx*delx + dely*dely !rr^2
                if (rr < rrsq) goto 300 
   250      continue 
   200   continue 
         xp(i,1)=xx1
         xp(i,2)=yy1
         xp(i,3)=zz1
100    continue  
    endif 
    
!---assign position of particle if volume fraction > 0.35       
!---assign position of particle if volume fraction > 0.35       
    if (volfra >=  0.35d0) then  
            write(*,*)'NOTE volume fraction > 0.35', volfra
            nout1=0
            nout2=0
            !inflate condition
            iratemax=100
            gap=0.d0 !mote this number
            
            !distribution of point particle
            do 400 ipm=1,npm                
              call RANDOM_NUMBER (xp(ipm,1))
              call RANDOM_NUMBER (xp(ipm,2))
              !call RANDOM_NUMBER (xp(ipm,3))
              xp(ipm,1)=(xp(ipm,1)-0.5d0)*bl ! box is [-L/2,L/2]
              xp(ipm,2)=(xp(ipm,2)-0.5d0)*bl ! box is [-L/2,L/2]
              !xp(ipm,3)=(xp(ipm,3)-0.5d0)*bl ! box is [-L/2,L/2]           
              xp(ipm,3)=0.d0
400         continue
            write(*,*)'initia',xp(1,1:3)
            !inflate the particles
            irate=0                     
            r=0.0d0 
410         irate=irate+1
            r=dble(irate)*1.0d0/(iratemax) ! radius during the inflation, final radius is 1
            dis=2.0d0*r + gap
            !write(*,*)'irate,r',irate,r
            
            !check
            do 420 ipm=1,npm       
            do 420 jpm=1,npm
               if (ipm >= jpm) goto 420               ! or only == case...  
               do k=1,3
                   ei(k)=xp(jpm,k)-xp(ipm,k)
                   ei(k)=ei(k)-bl*anint(ei(k)/bl) ! distance for minimum image
               enddo
                   ei(3)=0.d0 
               rr=sqrt(ei(1)*ei(1)+ei(2)*ei(2)+ei(3)*ei(3))
               if (rr < dis) then    !if overlap then we apply artificial repulsion...
                   xx1=xp(ipm,1)
                   yy1=xp(ipm,2)
                   zz1=xp(ipm,3) 
                      
                   xx2=xp(ipm,1)+ei(1)
                   yy2=xp(ipm,2)+ei(2)
                   zz2=xp(ipm,3)+ei(3)    
                          
                   430 delx=xx2-xx1   !from i to j
                       dely=yy2-yy1   !from i to j
                       delz=zz2-zz1   !from i to j            
                       call artificial_force(rr,r)
                       xx1=xx1-dxx*delx/rr ! we can simplify the code here
                       xx2=xx2+dxx*delx/rr
                       yy1=yy1-dxx*dely/rr
                       yy2=yy2+dxx*dely/rr
                       zz1=zz1-dxx*delz/rr
                       zz2=zz2+dxx*delz/rr  
                       rr=sqrt((xx2-xx1)*(xx2-xx1) +(yy2-xy1)*(yy2-yy1)+(zz2-zz1)*(zz2-zz1))
                       if (rr < dis) goto 430 !...until there is no overlap
                  
                   xx1=xx1-bl*anint(xx1/bl)     !if (xx1 < 0.0d0) xx1=xx1+side
                   yy1=yy1-bl*anint(yy1/bl)     !if (yy1 < 0.0d0) yy1=yy1+side
                   zz1=zz1-bl*anint(zz1/bl)     !if (zz1 < 0.0d0) zz1=zz1+side
                                                !if (xx1 > side)  xx1=xx1-side
                                                !if (yy1 > side)  yy1=yy1-side
                                                !if (zz1 > side)  zz1=zz1-side
                   xx2=xx2-bl*anint(xx2/bl)     !if (xx2 < 0.0d0) xx2=xx2+side
                   yy2=yy2-bl*anint(yy2/bl)     !if (yy2 < 0.0d0) yy2=yy2+side
                   zz2=zz2-bl*anint(zz2/bl)     !if (zz2 < 0.0d0) zz2=zz2+side
                                                !if (xx2 > side)  xx2=xx2-side
                                                !if (yy2 > side)  yy2=yy2-side
                                                !if (zz2 > side)  zz2=zz2-side  
                   xp(ipm,1)=xx1
                   xp(ipm,2)=yy1
                   xp(ipm,3)=zz1
                   xp(jpm,1)=xx2
                   xp(jpm,2)=yy2
                   xp(jpm,3)=zz2
                   do k=1,3
                     if (abs(xp(ipm,k)) > bl*0.5d0) then
                         write(*,*)'i out of range 1',ipm,k,xp(ipm,k)
                         nout1=nout1+1
                     endif 
                     if (abs(xp(jpm,k)) > bl*0.5d0) then
                         write(*,*)'j out of range 1',jpm,k,xp(jpm,k)
                         nout1=nout1+1
                     endif
                   enddo
                endif
        420 continue    
            
            !check the overlap again    
        600 ncheck=0
            dis=2.0d0 *r +gap 
            do 610 ipm=1,npm       
            do 610 jpm=1,npm
               if (ipm >= jpm) goto 610     
               do k=1,3
                   ei(k)=xp(jpm,k)-xp(ipm,k)
                   ei(k)=ei(k)-bl*anint(ei(k)/bl) !image
               enddo 
               rr=sqrt(ei(1)*ei(1)+ei(2)*ei(2)+ei(3)*ei(3))
               if (rr < dis) then 
                  ncheck=ncheck+1 
                  xx1=xp(ipm,1)
                  yy1=xp(ipm,2)
                  zz1=xp(ipm,3) 
                      
                  xx2=xp(ipm,1)+ei(1)
                  yy2=xp(ipm,2)+ei(2)
                  zz2=xp(ipm,3)+ei(3)    
                            
                  620  delx=xx2-xx1   !from i to j
                       dely=yy2-yy1   !from i to j
                       delz=zz2-zz1   !from i to j            
                       call artificial_force(rr,r)
                       xx1=xx1-dxx*delx/rr
                       xx2=xx2+dxx*delx/rr
                       yy1=yy1-dxx*dely/rr
                       yy2=yy2+dxx*dely/rr
                       zz1=zz1-dxx*delz/rr
                       zz2=zz2+dxx*delz/rr  
                       rr=sqrt((xx2-xx1)*(xx2-xx1) +(yy2-yy1)*(yy2-yy1)+(zz2-zz1)*(zz2-zz1))
                       if (rr < dis) goto 620
                       
                  xx1=xx1-bl*anint(xx1/bl)     !if (xx1 < 0.0d0) xx1=xx1+side
                  yy1=yy1-bl*anint(yy1/bl)     !if (yy1 < 0.0d0) yy1=yy1+side
                  zz1=zz1-bl*anint(zz1/bl)     !if (zz1 < 0.0d0) zz1=zz1+side
                                               !if (xx1 > side)  xx1=xx1-side
                                               !if (yy1 > side)  yy1=yy1-side
                                               !if (zz1 > side)  zz1=zz1-side
                  xx2=xx2-bl*anint(xx2/bl)     !if (xx2 < 0.0d0) xx2=xx2+side
                  yy2=yy2-bl*anint(yy2/bl)     !if (yy2 < 0.0d0) yy2=yy2+side
                  zz2=zz2-bl*anint(zz2/bl)     !if (zz2 < 0.0d0) zz2=zz2+side
                                               !if (xx2 > side)  xx2=xx2-side
                                               !if (yy2 > side)  yy2=yy2-side
                                               !if (zz2 > side)  zz2=zz2-side  
                  
                  xp(ipm,1)=xx1
                  xp(ipm,2)=yy1
                  xp(ipm,3)=zz1
                  xp(jpm,1)=xx2
                  xp(jpm,2)=yy2
                  xp(jpm,3)=zz2
                  
                  do k=1,3
                    if (abs(xp(ipm,k)) > bl*0.5d0) then
                        write(*,*)'i out of range 2',ipm,k,xp(ipm,k)
                        nout2=nout2+1
                    endif
                    if (abs(xp(jpm,k)) > bl*0.5d0) then
                        write(*,*)'j out of range 2',jpm,k,xp(jpm,k)
                        nout2=nout2+1
                    endif 
                  enddo
                endif        
            610  continue    
            !write(*,*)'ncheck',ncheck
            if (ncheck > 0) goto 600  
    if (irate < iratemax) goto 410
     !---final check for distance
         do 290 ipm=1,npm
         do 290 jpm=1,npm
            if(ipm >= jpm) go to 290          
            xx = sqrt( (xp(ipm,1)-xp(jpm,1))**2 + (xp(ipm,2)-xp(jpm,2))**2 + (xp(ipm,3)-xp(jpm,3))**2 )
            if (xx < 2.d0) write(*,*)'warning xx < 2',xx
        290   continue
  endif 
       
        
        
!---assign 3 random orthogonal vector of particle i.e. p0, p1, p2, where...
    !the relation to basis of SH is p0 = POLAR (Z), P1 AZIMUTHAL X, P2 AZIMUTHAL Y 
    AXB=0.D0;    AA=0.D0;    BB=0.D0
    do 500 i=1,npm 
    call RANDOM_NUMBER (pz)                       ! (0,1)    
    call RANDOM_NUMBER (phi)                      ! (0,1)          
    pz=(pz-0.5d0)*2.d0        ! (-1,1)
    phi=phi*2.d0*pi           ! (0, 2pi)
    p0(i,1)=cos(phi)*sqrt(1.d0-pz*pz)
    p0(i,2)=sin(phi)*sqrt(1.d0-pz*pz)
    p0(i,3)=pz       !at this point, p0 is truly random and |p0|=1   
    
    !choose another random point BB on sphere:
    call RANDOM_NUMBER (pz)                       ! (0,1)    
    call RANDOM_NUMBER (phi)                      ! (0,1)          
    pz=(pz-0.5d0)*2.d0        ! (-1,1)
    phi=phi*2.d0*pi           ! (0, 2pi)
    BB(1)=cos(phi)*sqrt(1.d0-pz*pz)
    BB(2)=sin(phi)*sqrt(1.d0-pz*pz)
    BB(3)=pz       !at this point, BB is truly random, |bb|=1
    ! next, calculate vector p1 = p0 x BB = AAxBB
    AA(1:3)=p0(i,1:3)
    AxB = cross(AA,BB)
    rr = sqrt( AxB(1)*AxB(1) + AxB(2)*AxB(2) + AxB(3)*AxB(3))
    p1(i,1:3)=AxB(1:3)/rr !normalize so that |p1|=1
            
    !CAL. P2 = P0 X P1         ! P2, P1 and P0 are Orthogonal 
    DO J=1,3
        AA(J)=P0(I,J)
        BB(J)=P1(I,J)
    ENDDO
     
    AxB=CROSS(AA,BB)      
    P2(I,1:3)=AxB(1:3) !|P2|=1 because p0=1, p1=1 and p0.p1=0
            
500 continue     
    end subroutine initia  
!:::::::::::::::::::::::::
    subroutine artificial_force(rr,r)
    implicit double precision (a-h,o-z)
    common /dxx/dxx
    
    
    !write(*,*)'force',rr,r
    rr0=rr/(2.0d0*r)
    dxx=(1.0d0/rr0)**1
    dxx=dxx*0.0001d0
    if (dxx > 1.0d0) write(*,*)'warning dxx too large, do not choose this set',dxx
    !write(*,*)'rr0,dxx',rr0,dxx
    end subroutine artificial_force
    
!----------------------------------------------------
    subroutine force
    use para
    use crossproduct !!!!!!!!!
    implicit double precision (a-h,o-z)
    
    dimension:: r(3),rhat(3)!,u(3),v(3) !xu(3),xv(3),
    
    dimension:: a0per(3),b0per(3),a1per(3),b1per(3), a2per(3),b2per(3)
    dimension:: a0xr(3),b0xr(3),a1xr(3),b1xr(3),a2xr(3),b2xr(3)
    dimension:: a0xb0(3),a0xb1(3),a0xb2(3)
    dimension:: a1xb0(3),a1xb1(3),a1xb2(3)
    dimension:: a2xb0(3),a2xb1(3),a2xb2(3)
               !per = perpendicular, x = cross product 
    dimension:: rtestipm(3),rtestjpm(3)     !test  ????
    common /ftab/ fa(3),fb(3),ta(3),tb(3)   !output f,t for pair
    
    !common variables for Aaniso
    common /ab/a0(3),b0(3),a1(3),b1(3),a2(3),b2(3)  !orientation vector  
    common /angle/a0r,b0r, a1r,b1r, a2r,b2r, a0b0,a0b1,a0b2, a1b0,a1b1,a1b2, a2b0,a2b1,a2b2     
    common /Aaniso/Aaniso,dAaniso_da0r,dAaniso_da1r,dAaniso_da2r,dAaniso_db0r,dAaniso_db1r,dAaniso_db2r, &
                          dAaniso_da0b0, dAaniso_da0b1, dAaniso_da0b2, &
                          dAaniso_da1b0, dAaniso_da1b1, dAaniso_da1b2, &
                          dAaniso_da2b0, dAaniso_da2b1, dAaniso_da2b2
    
    fp= 0.d0 ; tp=0.0d0
    energyI=0.d0
    do 100 ipm=1,npm-1      
    do 100 jpm=ipm+1,npm    
        
        fa = 0.0d0; ta= 0.0d0
        fb = 0.0d0; tb= 0.0d0           
        !distance vector from jpm* to ipm
        
        
        do k=1,2                            ! only x and y direction ON PLANE                
            r(k)=xp(jpm,k)-xp(ipm,k)        ! from ipm to jpm !Feb                   
        enddo
        r(1)=r(1)-anint(r(1)/bl)*bl !if cross x
        r(1)=r(1)-anint(r(2)/bl)*bl*pbcx; r(2)=r(2)-anint(r(2)/bl)*bl !if cross y
      
        rr = sqrt( r(1)*r(1) + r(2)*r(2) + r(3)*r(3) )  !distance
       
                        
        if (rr > rcut ) go to 200                 
        rhat=r/rr           !unit distance vector            
        a0(1:3)=p0(ipm,1:3);    b0(1:3)=p0(jpm,1:3)
        a1(1:3)=p1(ipm,1:3);    b1(1:3)=p1(jpm,1:3)
        a2(1:3)=p2(ipm,1:3);    b2(1:3)=p2(jpm,1:3)
               
        !calculate cosine i.e. dot product of basis vector an.r, bn.r
        a0r=0.d0; b0r=0.d0 
        a1r=0.d0; b1r=0.d0
        a2r=0.d0; b2r=0.d0
        a0b0=0.d0; a0b1=0.d0; a0b2=0.d0    !Add on Jan3 2018
        a1b0=0.d0; a1b1=0.d0; a1b2=0.d0    !Add on Jan3 2018
        a2b0=0.d0; a2b1=0.d0; a2b2=0.d0    !Add on Jan3 2018
             
        do i=1,3 !cosa=a.rhat
            a0r= a0r + a0(i)*rhat(i);    b0r= b0r + b0(i)*rhat(i)     
            a1r= a1r + a1(i)*rhat(i);    b1r= b1r + b1(i)*rhat(i)
            a2r= a2r + a2(i)*rhat(i);    b2r= b2r + b2(i)*rhat(i)             
            a0b0=a0b0 + a0(i)*b0(i)    !cosan.bn = an.bn
            a0b1=a0b1 + a0(i)*b1(i)
            a0b2=a0b2 + a0(i)*b2(i)
            a1b0=a1b0 + a1(i)*b0(i)
            a1b1=a1b1 + a1(i)*b1(i)
            a1b2=a1b2 + a1(i)*b2(i)
            a2b0=a2b0 + a2(i)*b0(i)
            a2b1=a2b1 + a2(i)*b1(i)
            a2b2=a2b2 + a2(i)*b2(i)
        enddo  
                 
        !---Calculate F on u 
            !calculate prefactor
            !vector !a perpendicular = an - (an.rhat)rhat
            a0per=a0-a0r*rhat;  b0per=b0-b0r*rhat   
            a1per=a1-a1r*rhat;  b1per=b1-b1r*rhat 
            a2per=a2-a2r*rhat;  b2per=b2-b2r*rhat         
            !an x r, bn x r, an x bm
            a0xr = cross(a0,rhat);  b0xr = cross(b0,rhat)
            a1xr = cross(a1,rhat);  b1xr = cross(b1,rhat)
            a2xr = cross(a2,rhat);  b2xr = cross(b2,rhat)                
            a0xb0=cross(a0,b0); a0xb1=cross(a0,b1); a0xb2=cross(a0,b2)
            a1xb0=cross(a1,b0); a1xb1=cross(a1,b1); a1xb2=cross(a1,b2)
            a2xb0=cross(a2,b0); a2xb1=cross(a2,b1); a2xb2=cross(a2,b2)
               
            !Potential 
            !Uiso = U_WCA  
            sigma12 =(sigma/rr)**12 ! rr is length distance
            sigma6  =(sigma/rr)**6        
            Uiso    = 4.0d0*(sigma12-sigma6+0.25d0)            !FOR rr <= 2^(1/6)   ! +0.25 due to WCA   
            dUiso_dr= -4.0d0*6.0d0*(2.0d0*sigma12-sigma6)/rr   !FOR rr <= 2^(1/6)
            
            rrlim=2.D0*RADIUS*2.d0**(1.d0/6.d0)                     !FOR rr > 2^(1/6)
            if (rr > rrlim) then                              
                Uiso = 0.d0
                dUiso_dr = 0.d0
            endif 
                   
            !Uaniso = U_Morse (refer to our paper)             
            req = 1.878d0   !pair equilibrium distance
            Uaniso = cmd*( (1.d0-exp(-(rr-req)/cmr))**2 - 1.d0 )
            dUaniso_dr = cmd*2.d0/cmr*(1.d0-exp((-rr+req)/cmr))*exp((-rr+req)/cmr)                                                                                
            
                 
            !Calculate Anisotropic interaction function and its derivatives              
            call calAaniso2 
            energyI(ipm)=energyI(ipm)+ Uiso + Aaniso*Uaniso
            energyI(jpm)=energyI(jpm)+ Uiso + Aaniso*Uaniso 
            
            do 10 i=1,3
            fa(i)=   ( (dUiso_dr + Aaniso* dUaniso_dr)*rhat(i) + &
                       1/rr*Uaniso*(dAaniso_da0r* a0per(i) + dAaniso_da1r* a1per(i) + dAaniso_da2r* a2per(i)+ &
                                    dAaniso_db0r* b0per(i) + dAaniso_db1r* b1per(i) + dAaniso_db2r* b2per(i)   ) )
            fb(i)= -fa(i)
             
            ta(i)= Uaniso* ( (-1.d0)*(dAaniso_da0r * a0xr(i) +  dAaniso_da1r * a1xr(i) + dAaniso_da2r * a2xr(i)) + &
                              (-1.d0)*(dAaniso_da0b0 * a0xb0(i) +  dAaniso_da0b1 * a0xb1(i) + dAaniso_da0b2 * a0xb2(i) +&
                                      dAaniso_da1b0 * a1xb0(i) +  dAaniso_da1b1 * a1xb1(i) + dAaniso_da1b2 * a1xb2(i) +&
                                      dAaniso_da2b0 * a2xb0(i) +  dAaniso_da2b1 * a2xb1(i) + dAaniso_da2b2 * a2xb2(i)) )
            tb(i)= Uaniso* ( (-1.d0)*(dAaniso_db0r * b0xr(i) +  dAaniso_db1r * b1xr(i) + dAaniso_db2r * b2xr(i)) + &
                                    (dAaniso_da0b0 * a0xb0(i) +  dAaniso_da0b1 * a0xb1(i) + dAaniso_da0b2 * a0xb2(i) +&
                                    dAaniso_da1b0 * a1xb0(i) +  dAaniso_da1b1 * a1xb1(i) + dAaniso_da1b2 * a1xb2(i) +&
                                    dAaniso_da2b0 * a2xb0(i) +  dAaniso_da2b1 * a2xb1(i) + dAaniso_da2b2 * a2xb2(i)) )
      10    continue
          
       !give back value from subroutine pair(ipm,jpm)
200    fp(ipm,1:3)=fp(ipm,1:3) + fa(1:3)  !!!!
       fp(jpm,1:3)=fp(jpm,1:3) + fb(1:3)
       tp(ipm,1:3)=tp(ipm,1:3) + ta(1:3)
       tp(jpm,1:3)=tp(jpm,1:3) + tb(1:3)
100 continue
  
    aveEnergyI=0.d0
    do ipm=1,npm
    aveEnergyI=aveEnergyI+energyI(ipm)
    enddo
    aveEnergyI=aveEnergyI/dble(npm)
    
    rlist=0.0 !if force is called then set rlist=0 (caution: force is called after write rlist data to file) 
    end subroutine force    
    
!----------------------------------------------------- 
    subroutine calAaniso2 
    USE para
    implicit double precision (a-h,o-z)
    
    !common variables for Aaniso
    common /ab/a0(3),b0(3),a1(3),b1(3),a2(3),b2(3)  !orientation vector an and bn  
    common /angle/a0r,b0r, a1r,b1r, a2r,b2r, a0b0,a0b1,a0b2, a1b0,a1b1,a1b2, a2b0,a2b1,a2b2 !dot product     
    common /Aaniso/Aaniso,dAaniso_da0r,dAaniso_da1r,dAaniso_da2r,dAaniso_db0r,dAaniso_db1r,dAaniso_db2r, &
                          dAaniso_da0b0, dAaniso_da0b1, dAaniso_da0b2, &
                          dAaniso_da1b0, dAaniso_da1b1, dAaniso_da1b2, &
                          dAaniso_da2b0, dAaniso_da2b1, dAaniso_da2b2
    
  
    
    
    if (lrank == 5 .and. mrank==5) then 
    Aaniso = ((8*a1b1**5 - 440*a1b1**4*a1r*b1r - 440*a1b2**4*a1r*b1r - 40*a1b1**3*(2*a1b2**2 + 2*a2b1**2 - 2*a2b2**2 - 44*a2b1*a2r*b1r - 143*a1r**2*b1r**2 + 143*a2r**2*b1r**2 - 44*a1b2*a1r*b2r + 44*a2b2*a2r*b2r + 143*a1r**2*b2r**2 - 143*a2r**2*b2r**2) - 80*a1b2**3*(2*a2b1*(a2b2 - 11*a2r*b2r) - 11*b1r*(2*a2b2*a2r + 13*(a1r**2 - a2r**2)*b2r)) + 440*a1b2**2*a1r*(-6*a2b1**2*b1r + 6*a2b2**2*b1r + 12*a2b1*a2b2*b2r - 156*a2b2*a2r*b1r*b2r + 65*(a1r**2 - 3*a2r**2)*b1r*(b1r**2 - 3*b2r**2) + 78*a2b1*a2r*(b1r**2 - b2r**2)) + 40*a1b1**2*(66*a1b2**2*a1r*b1r + 6*a1b2*(2*a2b1*(a2b2 - 11*a2r*b2r) - 11*b1r*(2*a2b2*a2r + 13*(a1r**2 - a2r**2)*b2r)) - 11*a1r*(-6*a2b1**2*b1r + 6*a2b2**2*b1r + 12*a2b1*a2b2*b2r - 156*a2b2*a2r*b1r*b2r + 65*(a1r**2 - 3*a2r**2)*b1r*(b1r**2 - 3*b2r**2) + 78*a2b1*a2r*(b1r**2 - b2r**2))) + 5*a1b1*(8*a1b2**4 + 8*a2b1**4 + 8*a2b2**4 - 352*a2b1**3*a2r*b1r + 3432*a1r**2*a2b2**2*b1r**2 - 3432*a2b2**2*a2r**2*b1r**2 + 12155*a1r**4*b1r**4 - 72930*a1r**2*a2r**2*b1r**4 + 12155*a2r**4*b1r**4 - 352*a1b2**3*a1r*b2r - 352*a2b2**3*a2r*b2r - 102960*a1r**2*a2b2*a2r*b1r**2*b2r + 34320*a2b2*a2r**3*b1r**2*b2r - 3432*a1r**2*a2b2**2*b2r**2 + 3432*a2b2**2*a2r**2*b2r**2 - 72930*a1r**4*b1r**2*b2r**2 + 437580*a1r**2*a2r**2*b1r**2*b2r**2 - 72930*a2r**4*b1r**2*b2r**2 + 34320*a1r**2*a2b2*a2r*b2r**3 - 11440*a2b2*a2r**3*b2r**3 + 12155*a1r**4*b2r**4 - 72930*a1r**2*a2r**2*b2r**4 + 12155*a2r**4*b2r**4 + 176*a2b1*b1r*(6*a2b2**2*a2r + 78*a2b2*(a1r**2 - a2r**2)*b2r - 65*a2r*(-3*a1r**2 + a2r**2)*(b1r**2 - 3*b2r**2)) + 24*a1b2**2*(2*a2b1**2 - 2*a2b2**2 - 44*a2b1*a2r*b1r + 44*a2b2*a2r*b2r - 143*(a1r**2 - a2r**2)*(b1r**2 - b2r**2)) - 24*a2b1**2*(2*a2b2**2 - 44*a2b2*a2r*b2r + 143*(a1r**2 - a2r**2)*(b1r**2 - b2r**2)) + 176*a1b2*a1r*(-6*a2b1**2*b2r + 6*a2b2**2*b2r - 12*a2b1*b1r*(a2b2 - 13*a2r*b2r) + 78*a2b2*a2r*(b1r**2 - b2r**2) - 65*(a1r**2 - 3*a2r**2)*b2r*(-3*b1r**2 + b2r**2))) - 20*a1b2*(8*a2b1**3*(a2b2 - 11*a2r*b2r) - 132*a2b1**2*b1r*(2*a2b2*a2r + 13*(a1r**2 - a2r**2)*b2r) + 11*b1r*(8*a2b2**3*a2r + 156*a2b2**2*(a1r**2 - a2r**2)*b2r - 260*a2b2*a2r*(-3*a1r**2 + a2r**2)*(b1r**2 - 3*b2r**2) + 1105*(a1r**4 - 6*a1r**2*a2r**2 + a2r**4)*b2r*(b1r**2 - b2r**2)) - 4*a2b1*(2*a2b2**3 - 66*a2b2**2*a2r*b2r + 429*a2b2*(a1r**2 - a2r**2)*(b1r**2 - b2r**2) - 715*a2r*(-3*a1r**2 + a2r**2)*b2r*(-3*b1r**2 + b2r**2))) - 11*a1r*(40*a2b1**4*b1r - 40*a2b1**2*b1r*(6*a2b2**2 - 156*a2b2*a2r*b2r + 65*(a1r**2 - 3*a2r**2)*(b1r**2 - 3*b2r**2)) - 80*a2b1**3*(2*a2b2*b2r + 13*a2r*(b1r**2 - b2r**2)) + 20*a2b1*(8*a2b2**3*b2r + 156*a2b2**2*a2r*(b1r**2 - b2r**2) - 260*a2b2*(a1r**2 - 3*a2r**2)*b2r*(-3*b1r**2 + b2r**2) + 1105*a2r*(a1r**2 - a2r**2)*(b1r**4 - 6*b1r**2*b2r**2 + b2r**4)) + b1r*(40*a2b2**4 - 2080*a2b2**3*a2r*b2r + 2600*a2b2**2*(a1r**2 - 3*a2r**2)*(b1r**2 - 3*b2r**2) - 88400*a2b2*a2r*(-a1r**2 + a2r**2)*b2r*(-b1r**2 + b2r**2) + 4199*(a1r**4 - 10*a1r**2*a2r**2 + 5*a2r**4)*(b1r**4 - 10*b1r**2*b2r**2 + 5*b2r**4))))*Pi)/704.
    dAaniso_da0r = 0.d0
    dAaniso_da1r = ((-440*a1b1**4*b1r - 440*a1b2**4*b1r + 22880*a1b2**3*a1r*b1r*b2r + 57200*a1b2**2*a1r**2*b1r*(b1r**2 - 3*b2r**2) - 40*a1b1**3*(-286*a1r*b1r**2 - 44*a1b2*b2r + 286*a1r*b2r**2) + 440*a1b2**2*(-6*a2b1**2*b1r + 6*a2b2**2*b1r + 12*a2b1*a2b2*b2r - 156*a2b2*a2r*b1r*b2r + 65*(a1r**2 - 3*a2r**2)*b1r*(b1r**2 - 3*b2r**2) + 78*a2b1*a2r*(b1r**2 - b2r**2)) + 40*a1b1**2*(66*a1b2**2*b1r - 1716*a1b2*a1r*b1r*b2r - 1430*a1r**2*b1r*(b1r**2 - 3*b2r**2) - 11*(-6*a2b1**2*b1r + 6*a2b2**2*b1r + 12*a2b1*a2b2*b2r - 156*a2b2*a2r*b1r*b2r + 65*(a1r**2 - 3*a2r**2)*b1r*(b1r**2 - 3*b2r**2) + 78*a2b1*a2r*(b1r**2 - b2r**2))) - 20*a1b2*(-3432*a1r*a2b1**2*b1r*b2r + 11*b1r*(312*a1r*a2b2**2*b2r + 1560*a1r*a2b2*a2r*(b1r**2 - 3*b2r**2) + 1105*(4*a1r**3 - 12*a1r*a2r**2)*b2r*(b1r**2 - b2r**2)) - 4*a2b1*(858*a1r*a2b2*(b1r**2 - b2r**2) + 4290*a1r*a2r*b2r*(-3*b1r**2 + b2r**2))) + 5*a1b1*(6864*a1r*a2b2**2*b1r**2 + 48620*a1r**3*b1r**4 - 145860*a1r*a2r**2*b1r**4 - 352*a1b2**3*b2r - 205920*a1r*a2b2*a2r*b1r**2*b2r - 6864*a1r*a2b2**2*b2r**2 - 291720*a1r**3*b1r**2*b2r**2 + 875160*a1r*a2r**2*b1r**2*b2r**2 + 68640*a1r*a2b2*a2r*b2r**3 + 48620*a1r**3*b2r**4 - 145860*a1r*a2r**2*b2r**4 - 6864*a1b2**2*a1r*(b1r**2 - b2r**2) - 6864*a1r*a2b1**2*(b1r**2 - b2r**2) - 22880*a1b2*a1r**2*b2r*(-3*b1r**2 + b2r**2) + 176*a2b1*b1r*(156*a1r*a2b2*b2r + 390*a1r*a2r*(b1r**2 - 3*b2r**2)) + 176*a1b2*(-6*a2b1**2*b2r + 6*a2b2**2*b2r - 12*a2b1*b1r*(a2b2 - 13*a2r*b2r) + 78*a2b2*a2r*(b1r**2 - b2r**2) - 65*(a1r**2 - 3*a2r**2)*b2r*(-3*b1r**2 + b2r**2))) - 11*a1r*(-5200*a1r*a2b1**2*b1r*(b1r**2 - 3*b2r**2) + 20*a2b1*(-520*a1r*a2b2*b2r*(-3*b1r**2 + b2r**2) + 2210*a1r*a2r*(b1r**4 - 6*b1r**2*b2r**2 + b2r**4)) + b1r*(5200*a1r*a2b2**2*(b1r**2 - 3*b2r**2) + 176800*a1r*a2b2*a2r*b2r*(-b1r**2 + b2r**2) + 4199*(4*a1r**3 - 20*a1r*a2r**2)*(b1r**4 - 10*b1r**2*b2r**2 + 5*b2r**4))) - 11*(40*a2b1**4*b1r - 40*a2b1**2*b1r*(6*a2b2**2 - 156*a2b2*a2r*b2r + 65*(a1r**2 - 3*a2r**2)*(b1r**2 - 3*b2r**2)) - 80*a2b1**3*(2*a2b2*b2r + 13*a2r*(b1r**2 - b2r**2)) + 20*a2b1*(8*a2b2**3*b2r + 156*a2b2**2*a2r*(b1r**2 - b2r**2) - 260*a2b2*(a1r**2 - 3*a2r**2)*b2r*(-3*b1r**2 + b2r**2) + 1105*a2r*(a1r**2 - a2r**2)*(b1r**4 - 6*b1r**2*b2r**2 + b2r**4)) + b1r*(40*a2b2**4 - 2080*a2b2**3*a2r*b2r + 2600*a2b2**2*(a1r**2 - 3*a2r**2)*(b1r**2 - 3*b2r**2) - 88400*a2b2*a2r*(-a1r**2 + a2r**2)*b2r*(-b1r**2 + b2r**2) + 4199*(a1r**4 - 10*a1r**2*a2r**2 + 5*a2r**4)*(b1r**4 - 10*b1r**2*b2r**2 + 5*b2r**4))))*Pi)/704.
    dAaniso_da2r = ((-40*a1b1**3*(-44*a2b1*b1r + 286*a2r*b1r**2 + 44*a2b2*b2r - 286*a2r*b2r**2) - 80*a1b2**3*(-22*a2b1*b2r - 11*b1r*(2*a2b2 - 26*a2r*b2r)) + 440*a1b2**2*a1r*(-156*a2b2*b1r*b2r - 390*a2r*b1r*(b1r**2 - 3*b2r**2) + 78*a2b1*(b1r**2 - b2r**2)) + 40*a1b1**2*(6*a1b2*(-22*a2b1*b2r - 11*b1r*(2*a2b2 - 26*a2r*b2r)) - 11*a1r*(-156*a2b2*b1r*b2r - 390*a2r*b1r*(b1r**2 - 3*b2r**2) + 78*a2b1*(b1r**2 - b2r**2))) + 5*a1b1*(-352*a2b1**3*b1r - 6864*a2b2**2*a2r*b1r**2 - 145860*a1r**2*a2r*b1r**4 + 48620*a2r**3*b1r**4 - 352*a2b2**3*b2r - 102960*a1r**2*a2b2*b1r**2*b2r + 102960*a2b2*a2r**2*b1r**2*b2r + 6864*a2b2**2*a2r*b2r**2 + 875160*a1r**2*a2r*b1r**2*b2r**2 - 291720*a2r**3*b1r**2*b2r**2 + 34320*a1r**2*a2b2*b2r**3 - 34320*a2b2*a2r**2*b2r**3 - 145860*a1r**2*a2r*b2r**4 + 48620*a2r**3*b2r**4 + 176*a2b1*b1r*(6*a2b2**2 - 156*a2b2*a2r*b2r - 130*a2r**2*(b1r**2 - 3*b2r**2) - 65*(-3*a1r**2 + a2r**2)*(b1r**2 - 3*b2r**2)) - 24*a2b1**2*(-44*a2b2*b2r - 286*a2r*(b1r**2 - b2r**2)) + 24*a1b2**2*(-44*a2b1*b1r + 44*a2b2*b2r + 286*a2r*(b1r**2 - b2r**2)) + 176*a1b2*a1r*(156*a2b1*b1r*b2r + 78*a2b2*(b1r**2 - b2r**2) + 390*a2r*b2r*(-3*b1r**2 + b2r**2))) - 20*a1b2*(-88*a2b1**3*b2r - 132*a2b1**2*b1r*(2*a2b2 - 26*a2r*b2r) + 11*b1r*(8*a2b2**3 - 312*a2b2**2*a2r*b2r - 520*a2b2*a2r**2*(b1r**2 - 3*b2r**2) - 260*a2b2*(-3*a1r**2 + a2r**2)*(b1r**2 - 3*b2r**2) + 1105*(-12*a1r**2*a2r + 4*a2r**3)*b2r*(b1r**2 - b2r**2)) - 4*a2b1*(-66*a2b2**2*b2r - 858*a2b2*a2r*(b1r**2 - b2r**2) - 1430*a2r**2*b2r*(-3*b1r**2 + b2r**2) - 715*(-3*a1r**2 + a2r**2)*b2r*(-3*b1r**2 + b2r**2))) - 11*a1r*(-1040*a2b1**3*(b1r**2 - b2r**2) - 40*a2b1**2*b1r*(-156*a2b2*b2r - 390*a2r*(b1r**2 - 3*b2r**2)) + 20*a2b1*(156*a2b2**2*(b1r**2 - b2r**2) + 1560*a2b2*a2r*b2r*(-3*b1r**2 + b2r**2) - 2210*a2r**2*(b1r**4 - 6*b1r**2*b2r**2 + b2r**4) + 1105*(a1r**2 - a2r**2)*(b1r**4 - 6*b1r**2*b2r**2 + b2r**4)) + b1r*(-2080*a2b2**3*b2r - 15600*a2b2**2*a2r*(b1r**2 - 3*b2r**2) - 176800*a2b2*a2r**2*b2r*(-b1r**2 + b2r**2) - 88400*a2b2*(-a1r**2 + a2r**2)*b2r*(-b1r**2 + b2r**2) + 4199*(-20*a1r**2*a2r + 20*a2r**3)*(b1r**4 - 10*b1r**2*b2r**2 + 5*b2r**4))))*Pi)/704.
    dAaniso_db0r = 0.d0
    dAaniso_db1r = ((-440*a1b1**4*a1r - 440*a1b2**4*a1r - 40*a1b1**3*(-44*a2b1*a2r - 286*a1r**2*b1r + 286*a2r**2*b1r) + 880*a1b2**3*(2*a2b2*a2r + 13*(a1r**2 - a2r**2)*b2r) + 440*a1b2**2*a1r*(-6*a2b1**2 + 6*a2b2**2 + 156*a2b1*a2r*b1r + 130*(a1r**2 - 3*a2r**2)*b1r**2 - 156*a2b2*a2r*b2r + 65*(a1r**2 - 3*a2r**2)*(b1r**2 - 3*b2r**2)) + 40*a1b1**2*(66*a1b2**2*a1r - 66*a1b2*(2*a2b2*a2r + 13*(a1r**2 - a2r**2)*b2r) - 11*a1r*(-6*a2b1**2 + 6*a2b2**2 + 156*a2b1*a2r*b1r + 130*(a1r**2 - 3*a2r**2)*b1r**2 - 156*a2b2*a2r*b2r + 65*(a1r**2 - 3*a2r**2)*(b1r**2 - 3*b2r**2))) + 5*a1b1*(-352*a2b1**3*a2r + 6864*a1r**2*a2b2**2*b1r - 6864*a2b2**2*a2r**2*b1r - 6864*a2b1**2*(a1r**2 - a2r**2)*b1r - 22880*a2b1*a2r*(-3*a1r**2 + a2r**2)*b1r**2 + 48620*a1r**4*b1r**3 - 291720*a1r**2*a2r**2*b1r**3 + 48620*a2r**4*b1r**3 + 24*a1b2**2*(-44*a2b1*a2r - 286*(a1r**2 - a2r**2)*b1r) - 205920*a1r**2*a2b2*a2r*b1r*b2r + 68640*a2b2*a2r**3*b1r*b2r - 145860*a1r**4*b1r*b2r**2 + 875160*a1r**2*a2r**2*b1r*b2r**2 - 145860*a2r**4*b1r*b2r**2 + 176*a1b2*a1r*(156*a2b2*a2r*b1r + 390*(a1r**2 - 3*a2r**2)*b1r*b2r - 12*a2b1*(a2b2 - 13*a2r*b2r)) + 176*a2b1*(6*a2b2**2*a2r + 78*a2b2*(a1r**2 - a2r**2)*b2r - 65*a2r*(-3*a1r**2 + a2r**2)*(b1r**2 - 3*b2r**2))) - 20*a1b2*(-132*a2b1**2*(2*a2b2*a2r + 13*(a1r**2 - a2r**2)*b2r) - 4*a2b1*(858*a2b2*(a1r**2 - a2r**2)*b1r + 4290*a2r*(-3*a1r**2 + a2r**2)*b1r*b2r) + 11*b1r*(-520*a2b2*a2r*(-3*a1r**2 + a2r**2)*b1r + 2210*(a1r**4 - 6*a1r**2*a2r**2 + a2r**4)*b1r*b2r) + 11*(8*a2b2**3*a2r + 156*a2b2**2*(a1r**2 - a2r**2)*b2r - 260*a2b2*a2r*(-3*a1r**2 + a2r**2)*(b1r**2 - 3*b2r**2) + 1105*(a1r**4 - 6*a1r**2*a2r**2 + a2r**4)*b2r*(b1r**2 - b2r**2))) - 11*a1r*(40*a2b1**4 + 40*a2b2**4 - 2080*a2b1**3*a2r*b1r - 5200*a2b1**2*(a1r**2 - 3*a2r**2)*b1r**2 - 2080*a2b2**3*a2r*b2r + 2600*a2b2**2*(a1r**2 - 3*a2r**2)*(b1r**2 - 3*b2r**2) - 88400*a2b2*a2r*(-a1r**2 + a2r**2)*b2r*(-b1r**2 + b2r**2) + 4199*(a1r**4 - 10*a1r**2*a2r**2 + 5*a2r**4)*(b1r**4 - 10*b1r**2*b2r**2 + 5*b2r**4) - 40*a2b1**2*(6*a2b2**2 - 156*a2b2*a2r*b2r + 65*(a1r**2 - 3*a2r**2)*(b1r**2 - 3*b2r**2)) + b1r*(5200*a2b2**2*(a1r**2 - 3*a2r**2)*b1r + 176800*a2b2*a2r*(-a1r**2 + a2r**2)*b1r*b2r + 4199*(a1r**4 - 10*a1r**2*a2r**2 + 5*a2r**4)*(4*b1r**3 - 20*b1r*b2r**2)) + 20*a2b1*(312*a2b2**2*a2r*b1r + 1560*a2b2*(a1r**2 - 3*a2r**2)*b1r*b2r + 1105*a2r*(a1r**2 - a2r**2)*(4*b1r**3 - 12*b1r*b2r**2))))*Pi)/704.
    dAaniso_db2r = ((-80*a1b2**3*(-22*a2b1*a2r - 143*(a1r**2 - a2r**2)*b1r) - 40*a1b1**3*(-44*a1b2*a1r + 44*a2b2*a2r + 286*a1r**2*b2r - 286*a2r**2*b2r) + 440*a1b2**2*a1r*(12*a2b1*a2b2 - 156*a2b2*a2r*b1r - 156*a2b1*a2r*b2r - 390*(a1r**2 - 3*a2r**2)*b1r*b2r) + 40*a1b1**2*(6*a1b2*(-22*a2b1*a2r - 143*(a1r**2 - a2r**2)*b1r) - 11*a1r*(12*a2b1*a2b2 - 156*a2b2*a2r*b1r - 156*a2b1*a2r*b2r - 390*(a1r**2 - 3*a2r**2)*b1r*b2r)) + 5*a1b1*(-352*a1b2**3*a1r - 352*a2b2**3*a2r - 102960*a1r**2*a2b2*a2r*b1r**2 + 34320*a2b2*a2r**3*b1r**2 - 6864*a1r**2*a2b2**2*b2r + 6864*a2b2**2*a2r**2*b2r - 145860*a1r**4*b1r**2*b2r + 875160*a1r**2*a2r**2*b1r**2*b2r - 145860*a2r**4*b1r**2*b2r + 102960*a1r**2*a2b2*a2r*b2r**2 - 34320*a2b2*a2r**3*b2r**2 + 48620*a1r**4*b2r**3 - 291720*a1r**2*a2r**2*b2r**3 + 48620*a2r**4*b2r**3 - 24*a2b1**2*(-44*a2b2*a2r - 286*(a1r**2 - a2r**2)*b2r) + 24*a1b2**2*(44*a2b2*a2r + 286*(a1r**2 - a2r**2)*b2r) + 176*a2b1*b1r*(78*a2b2*(a1r**2 - a2r**2) + 390*a2r*(-3*a1r**2 + a2r**2)*b2r) + 176*a1b2*a1r*(-6*a2b1**2 + 6*a2b2**2 + 156*a2b1*a2r*b1r - 156*a2b2*a2r*b2r - 130*(a1r**2 - 3*a2r**2)*b2r**2 - 65*(a1r**2 - 3*a2r**2)*(-3*b1r**2 + b2r**2))) - 20*a1b2*(-88*a2b1**3*a2r - 1716*a2b1**2*(a1r**2 - a2r**2)*b1r + 11*b1r*(156*a2b2**2*(a1r**2 - a2r**2) + 1560*a2b2*a2r*(-3*a1r**2 + a2r**2)*b2r - 2210*(a1r**4 - 6*a1r**2*a2r**2 + a2r**4)*b2r**2 + 1105*(a1r**4 - 6*a1r**2*a2r**2 + a2r**4)*(b1r**2 - b2r**2)) - 4*a2b1*(-66*a2b2**2*a2r - 858*a2b2*(a1r**2 - a2r**2)*b2r - 1430*a2r*(-3*a1r**2 + a2r**2)*b2r**2 - 715*a2r*(-3*a1r**2 + a2r**2)*(-3*b1r**2 + b2r**2))) - 11*a1r*(-80*a2b1**3*(2*a2b2 - 26*a2r*b2r) - 40*a2b1**2*b1r*(-156*a2b2*a2r - 390*(a1r**2 - 3*a2r**2)*b2r) + 20*a2b1*(8*a2b2**3 - 312*a2b2**2*a2r*b2r - 520*a2b2*(a1r**2 - 3*a2r**2)*b2r**2 - 260*a2b2*(a1r**2 - 3*a2r**2)*(-3*b1r**2 + b2r**2) + 1105*a2r*(a1r**2 - a2r**2)*(-12*b1r**2*b2r + 4*b2r**3)) + b1r*(-2080*a2b2**3*a2r - 15600*a2b2**2*(a1r**2 - 3*a2r**2)*b2r - 176800*a2b2*a2r*(-a1r**2 + a2r**2)*b2r**2 - 88400*a2b2*a2r*(-a1r**2 + a2r**2)*(-b1r**2 + b2r**2) + 4199*(a1r**4 - 10*a1r**2*a2r**2 + 5*a2r**4)*(-20*b1r**2*b2r + 20*b2r**3))))*Pi)/704.
    dAaniso_da0b0 = 0.d0
    dAaniso_da0b1 = 0.d0
    dAaniso_da0b2 = 0.d0
    dAaniso_da1b0 = 0.d0
    dAaniso_da1b1 = ((40*a1b1**4 - 1760*a1b1**3*a1r*b1r - 120*a1b1**2*(2*a1b2**2 + 2*a2b1**2 - 2*a2b2**2 - 44*a2b1*a2r*b1r - 143*a1r**2*b1r**2 + 143*a2r**2*b1r**2 - 44*a1b2*a1r*b2r + 44*a2b2*a2r*b2r + 143*a1r**2*b2r**2 - 143*a2r**2*b2r**2) + 80*a1b1*(66*a1b2**2*a1r*b1r + 6*a1b2*(2*a2b1*(a2b2 - 11*a2r*b2r) - 11*b1r*(2*a2b2*a2r + 13*(a1r**2 - a2r**2)*b2r)) - 11*a1r*(-6*a2b1**2*b1r + 6*a2b2**2*b1r + 12*a2b1*a2b2*b2r - 156*a2b2*a2r*b1r*b2r + 65*(a1r**2 - 3*a2r**2)*b1r*(b1r**2 - 3*b2r**2) + 78*a2b1*a2r*(b1r**2 - b2r**2))) + 5*(8*a1b2**4 + 8*a2b1**4 + 8*a2b2**4 - 352*a2b1**3*a2r*b1r + 3432*a1r**2*a2b2**2*b1r**2 - 3432*a2b2**2*a2r**2*b1r**2 + 12155*a1r**4*b1r**4 - 72930*a1r**2*a2r**2*b1r**4 + 12155*a2r**4*b1r**4 - 352*a1b2**3*a1r*b2r - 352*a2b2**3*a2r*b2r - 102960*a1r**2*a2b2*a2r*b1r**2*b2r + 34320*a2b2*a2r**3*b1r**2*b2r - 3432*a1r**2*a2b2**2*b2r**2 + 3432*a2b2**2*a2r**2*b2r**2 - 72930*a1r**4*b1r**2*b2r**2 + 437580*a1r**2*a2r**2*b1r**2*b2r**2 - 72930*a2r**4*b1r**2*b2r**2 + 34320*a1r**2*a2b2*a2r*b2r**3 - 11440*a2b2*a2r**3*b2r**3 + 12155*a1r**4*b2r**4 - 72930*a1r**2*a2r**2*b2r**4 + 12155*a2r**4*b2r**4 + 176*a2b1*b1r*(6*a2b2**2*a2r + 78*a2b2*(a1r**2 - a2r**2)*b2r - 65*a2r*(-3*a1r**2 + a2r**2)*(b1r**2 - 3*b2r**2)) + 24*a1b2**2*(2*a2b1**2 - 2*a2b2**2 - 44*a2b1*a2r*b1r + 44*a2b2*a2r*b2r - 143*(a1r**2 - a2r**2)*(b1r**2 - b2r**2)) - 24*a2b1**2*(2*a2b2**2 - 44*a2b2*a2r*b2r + 143*(a1r**2 - a2r**2)*(b1r**2 - b2r**2)) + 176*a1b2*a1r*(-6*a2b1**2*b2r + 6*a2b2**2*b2r - 12*a2b1*b1r*(a2b2 - 13*a2r*b2r) + 78*a2b2*a2r*(b1r**2 - b2r**2) - 65*(a1r**2 - 3*a2r**2)*b2r*(-3*b1r**2 + b2r**2))))*Pi)/704.
    dAaniso_da1b2 = ((-1760*a1b2**3*a1r*b1r - 40*a1b1**3*(4*a1b2 - 44*a1r*b2r) - 240*a1b2**2*(2*a2b1*(a2b2 - 11*a2r*b2r) - 11*b1r*(2*a2b2*a2r + 13*(a1r**2 - a2r**2)*b2r)) + 880*a1b2*a1r*(-6*a2b1**2*b1r + 6*a2b2**2*b1r + 12*a2b1*a2b2*b2r - 156*a2b2*a2r*b1r*b2r + 65*(a1r**2 - 3*a2r**2)*b1r*(b1r**2 - 3*b2r**2) + 78*a2b1*a2r*(b1r**2 - b2r**2)) + 40*a1b1**2*(132*a1b2*a1r*b1r + 6*(2*a2b1*(a2b2 - 11*a2r*b2r) - 11*b1r*(2*a2b2*a2r + 13*(a1r**2 - a2r**2)*b2r))) + 5*a1b1*(32*a1b2**3 - 1056*a1b2**2*a1r*b2r + 48*a1b2*(2*a2b1**2 - 2*a2b2**2 - 44*a2b1*a2r*b1r + 44*a2b2*a2r*b2r - 143*(a1r**2 - a2r**2)*(b1r**2 - b2r**2)) + 176*a1r*(-6*a2b1**2*b2r + 6*a2b2**2*b2r - 12*a2b1*b1r*(a2b2 - 13*a2r*b2r) + 78*a2b2*a2r*(b1r**2 - b2r**2) - 65*(a1r**2 - 3*a2r**2)*b2r*(-3*b1r**2 + b2r**2))) - 20*(8*a2b1**3*(a2b2 - 11*a2r*b2r) - 132*a2b1**2*b1r*(2*a2b2*a2r + 13*(a1r**2 - a2r**2)*b2r) + 11*b1r*(8*a2b2**3*a2r + 156*a2b2**2*(a1r**2 - a2r**2)*b2r - 260*a2b2*a2r*(-3*a1r**2 + a2r**2)*(b1r**2 - 3*b2r**2) + 1105*(a1r**4 - 6*a1r**2*a2r**2 + a2r**4)*b2r*(b1r**2 - b2r**2)) - 4*a2b1*(2*a2b2**3 - 66*a2b2**2*a2r*b2r + 429*a2b2*(a1r**2 - a2r**2)*(b1r**2 - b2r**2) - 715*a2r*(-3*a1r**2 + a2r**2)*b2r*(-3*b1r**2 + b2r**2))))*Pi)/704.
    dAaniso_da2b0 = 0.d0
    dAaniso_da2b1 = ((-40*a1b1**3*(4*a2b1 - 44*a2r*b1r) - 160*a1b2**3*(a2b2 - 11*a2r*b2r) + 440*a1b2**2*a1r*(-12*a2b1*b1r + 12*a2b2*b2r + 78*a2r*(b1r**2 - b2r**2)) + 40*a1b1**2*(12*a1b2*(a2b2 - 11*a2r*b2r) - 11*a1r*(-12*a2b1*b1r + 12*a2b2*b2r + 78*a2r*(b1r**2 - b2r**2))) + 5*a1b1*(32*a2b1**3 - 1056*a2b1**2*a2r*b1r + 24*a1b2**2*(4*a2b1 - 44*a2r*b1r) + 176*a1b2*a1r*(-12*a2b1*b2r - 12*b1r*(a2b2 - 13*a2r*b2r)) + 176*b1r*(6*a2b2**2*a2r + 78*a2b2*(a1r**2 - a2r**2)*b2r - 65*a2r*(-3*a1r**2 + a2r**2)*(b1r**2 - 3*b2r**2)) - 48*a2b1*(2*a2b2**2 - 44*a2b2*a2r*b2r + 143*(a1r**2 - a2r**2)*(b1r**2 - b2r**2))) - 20*a1b2*(24*a2b1**2*(a2b2 - 11*a2r*b2r) - 264*a2b1*b1r*(2*a2b2*a2r + 13*(a1r**2 - a2r**2)*b2r) - 4*(2*a2b2**3 - 66*a2b2**2*a2r*b2r + 429*a2b2*(a1r**2 - a2r**2)*(b1r**2 - b2r**2) - 715*a2r*(-3*a1r**2 + a2r**2)*b2r*(-3*b1r**2 + b2r**2))) - 11*a1r*(160*a2b1**3*b1r - 80*a2b1*b1r*(6*a2b2**2 - 156*a2b2*a2r*b2r + 65*(a1r**2 - 3*a2r**2)*(b1r**2 - 3*b2r**2)) - 240*a2b1**2*(2*a2b2*b2r + 13*a2r*(b1r**2 - b2r**2)) + 20*(8*a2b2**3*b2r + 156*a2b2**2*a2r*(b1r**2 - b2r**2) - 260*a2b2*(a1r**2 - 3*a2r**2)*b2r*(-3*b1r**2 + b2r**2) + 1105*a2r*(a1r**2 - a2r**2)*(b1r**4 - 6*b1r**2*b2r**2 + b2r**4))))*Pi)/704.
    dAaniso_da2b2 = ((-80*a1b2**3*(2*a2b1 - 22*a2r*b1r) - 40*a1b1**3*(-4*a2b2 + 44*a2r*b2r) + 440*a1b2**2*a1r*(12*a2b2*b1r + 12*a2b1*b2r - 156*a2r*b1r*b2r) + 40*a1b1**2*(6*a1b2*(2*a2b1 - 22*a2r*b1r) - 11*a1r*(12*a2b2*b1r + 12*a2b1*b2r - 156*a2r*b1r*b2r)) + 5*a1b1*(32*a2b2**3 + 6864*a1r**2*a2b2*b1r**2 - 6864*a2b2*a2r**2*b1r**2 - 1056*a2b2**2*a2r*b2r - 102960*a1r**2*a2r*b1r**2*b2r + 34320*a2r**3*b1r**2*b2r - 6864*a1r**2*a2b2*b2r**2 + 6864*a2b2*a2r**2*b2r**2 + 34320*a1r**2*a2r*b2r**3 - 11440*a2r**3*b2r**3 - 24*a2b1**2*(4*a2b2 - 44*a2r*b2r) + 24*a1b2**2*(-4*a2b2 + 44*a2r*b2r) + 176*a2b1*b1r*(12*a2b2*a2r + 78*(a1r**2 - a2r**2)*b2r) + 176*a1b2*a1r*(-12*a2b1*b1r + 12*a2b2*b2r + 78*a2r*(b1r**2 - b2r**2))) - 20*a1b2*(8*a2b1**3 - 264*a2b1**2*a2r*b1r + 11*b1r*(24*a2b2**2*a2r + 312*a2b2*(a1r**2 - a2r**2)*b2r - 260*a2r*(-3*a1r**2 + a2r**2)*(b1r**2 - 3*b2r**2)) - 4*a2b1*(6*a2b2**2 - 132*a2b2*a2r*b2r + 429*(a1r**2 - a2r**2)*(b1r**2 - b2r**2))) - 11*a1r*(-160*a2b1**3*b2r - 40*a2b1**2*b1r*(12*a2b2 - 156*a2r*b2r) + 20*a2b1*(24*a2b2**2*b2r + 312*a2b2*a2r*(b1r**2 - b2r**2) - 260*(a1r**2 - 3*a2r**2)*b2r*(-3*b1r**2 + b2r**2)) + b1r*(160*a2b2**3 - 6240*a2b2**2*a2r*b2r + 5200*a2b2*(a1r**2 - 3*a2r**2)*(b1r**2 - 3*b2r**2) - 88400*a2r*(-a1r**2 + a2r**2)*b2r*(-b1r**2 + b2r**2))))*Pi)/704.    
    
    endif  
    
   !Normalisation
    Aaniso          =coe* Aaniso/amax       
    dAaniso_da0r    =coe* dAaniso_da0r/amax
    dAaniso_da1r    =coe* dAaniso_da1r/amax
    dAaniso_da2r    =coe* dAaniso_da2r/amax
    dAaniso_db0r    =coe* dAaniso_db0r/amax
    dAaniso_db1r    =coe* dAaniso_db1r/amax
    dAaniso_db2r    =coe* dAaniso_db2r/amax
                     
    dAaniso_da0b0   =coe* dAaniso_da0b0/amax
    dAaniso_da0b1   =coe* dAaniso_da0b1/amax
    dAaniso_da0b2   =coe* dAaniso_da0b2/amax
    dAaniso_da1b0   =coe* dAaniso_da1b0/amax
    dAaniso_da1b1   =coe* dAaniso_da1b1/amax
    dAaniso_da1b2   =coe* dAaniso_da1b2/amax
    dAaniso_da2b0   =coe* dAaniso_da2b0/amax
    dAaniso_da2b1   =coe* dAaniso_da2b1/amax
    dAaniso_da2b2   =coe* dAaniso_da2b2/amax
  
    end subroutine calAaniso2     
    
!----------------------------------------------------  
    subroutine covar  
    use para
    implicit double precision (a-h,o-z)
            
    !calculate correlated random displacement
    do 10 i=1,npm*3
        xic(i) =gauss(dummy) * sqrt(2.0d0*dt)   !noise
        xicr(i)=gauss(dummy) * sqrt(1.5D0*dt)   !noise
       
10  continue
    
    end subroutine covar   
    
!----------------------------------------------------
    subroutine calpos  
    use para
    implicit double precision (a-h,o-z)
    dimension ::fpA(npm*3),xpA(npm*3)
    
    fpa=0.d0; xpA=0.d0
    !fp(2,2)=100.d0
    jj=0
    do  j=1,3                     
    do  i=1,npm                        
        jj=jj+1
        fpa(jj)=fp(i,j)
        xpa(jj)=xp(i,j)
    enddo
    enddo   
        
    !cal. Dij*Fj *dt/tem +ri    
    do i=1,npm*3
        xpa(i)=xpa(i)+fpa(i)*dt*BIGE + xic(i)            
    enddo
    
    jj=0
    do  j=1,3                    
    do  i=1,npm                      
        jj=jj+1
        xp(i,j)=xpa(jj)      
    enddo
    enddo 
     
    do ipm=1,npm            
        xp(ipm,3)=0.d0      
    enddo
    
    !ADD the shear displacement
    do ipm=1,npm
        xp(ipm,1)=xp(ipm,1)+Pe*dt*xp(ipm,2)    
    enddo
    
    !Treat the PBC issue if sinusoidal shear is applied
    Pe=Pe0*sin(omg*dt*istep) ! the sinusoidal shear rate gamadot 
    pbcx=pbcx+ Pe*dt     
    do ipm=1,npm       
       xp(ipm,1)= xp(ipm,1) - anint((xp(ipm,1)-pbcx*xp(ipm,2))/bl)*bl       
       xp(ipm,1)= xp(ipm,1) - anint((xp(ipm,2))/bl)*bl *pbcx                           
       xp(ipm,2)= xp(ipm,2) - anint((xp(ipm,2))/bl)*bl                 
    enddo
             
    end subroutine calpos 
    
!----------------------------------------------------
    subroutine calorient  !if d is constant ROTATION only
    use para
    implicit double precision (a-h,o-z)
    dimension :: ps0(3),ps1(3),ps2(3),prot0(3),prot1(3),prot2(3),rx(3,3),ry(3,3),rz(3,3),rzry(3,3),rzryrx(3,3) 
    
   
    !for particle i--
    do 100 i=1,npm      
    ps0(1:3)=p0(i,1:3); ps1(1:3)=p1(i,1:3); ps2(1:3)=p2(i,1:3)
    
    !-------------------------------------------INCLUDE ROTATION FROM TORQUE
    ii=(i-1)*3
    ax=0.75D0*tp(i,1)*dt*BIGE + xicr(ii+1)         ! noise + torque
    ay=0.75D0*tp(i,2)*dt*BIGE + xicr(ii+2)         ! noise + torque           
    az=0.75D0*tp(i,3)*dt*BIGE + xicr(ii+3) - Pe*0.5d0*dt        ! noise + torque   + shear (plus or minus?)
    !-------------------------------------------INCLUDE ROTATION FROM TORQUE
    
    !Consider the (orientation vector of) particle rotates about x,y,z angle ax,ay,az 
    !                                                   [Refer to Rodrigue formula]    
    !rotation vector
    rx(1,1)=1.d0;        rx(1,2)=0.d0;       rx(1,3)= 0.d0
    rx(2,1)=0.d0;        rx(2,2)=cos(ax);    rx(2,3)=-sin(ax)
    rx(3,1)=0.d0;        rx(3,2)=sin(ax);    rx(3,3)= cos(ax)
    
    ry(1,1)= cos(ay);   ry(1,2)=0.d0;        ry(1,3)= sin(ay)
    ry(2,1)= 0.d0;      ry(2,2)=1.d0;        ry(2,3)= 0.d0
    ry(3,1)=-sin(ay);   ry(3,2)=0.d0;        ry(3,3)= cos(ay) 
    
    rz(1,1)=cos(az);    rz(1,2)=-sin(az);   rz(1,3)= 0.d0
    rz(2,1)=sin(az);    rz(2,2)= cos(az);   rz(2,3)= 0.d0
    rz(3,1)=0.d0;       rz(3,2)=0.d0;       rz(3,3)= 1.d0
    
    rzry    = 0.0d0;    rzryrz  = 0.0d0;    prot    = 0.0d0
    rzry    = matmul(rz,ry)
    rzryrx  = matmul(rzry,rx)
    prot0    = matmul(rzryrx,ps0) !ROTATED orientation vector
    prot1    = matmul(rzryrx,ps1)
    prot2    = matmul(rzryrx,ps2)
    p0(i,1:3)=prot0(1:3)
    p1(i,1:3)=prot1(1:3)
    p2(i,1:3)=prot2(1:3)
   

100 continue 

    end subroutine calorient 

!----------------------------------------------------     
    function gauss (dummy)
   !return a uniform random normal variate from distribution zero mean and unit variance
   !Ref:Knuth D, The art of computer programming 2ed                                           
   
   implicit double precision (a-h,o-z)
   A1 = 3.949846138; A3 = 0.252408784; A5 = 0.076542912 
   A7 = 0.008355968; A9 = 0.02989977

   SUM = 0.0d0
   !call RANDOM_SEED() 
   do i = 1, 12
 !     SUM = SUM + RANF ( DUMMY )
      
      call RANDOM_NUMBER(a)
      sum =sum +  a
   enddo

   R  = ( SUM - 6.0 ) / 4.0
   R2 = R * R

   gauss=((((A9*R2 + A7 ) * R2 + A5 ) * R2 + A3 ) * R2 + A1 )* R

    end function gauss
    
!---------------------------------------------------- 
    subroutine force_verlet
    use para
    use crossproduct !!!!!!!!!
    implicit double precision (a-h,o-z)
    
    dimension:: r(3),rhat(3)!,u(3),v(3) !xu(3),xv(3),
    
    dimension:: a0per(3),b0per(3),a1per(3),b1per(3), a2per(3),b2per(3)
    dimension:: a0xr(3),b0xr(3),a1xr(3),b1xr(3),a2xr(3),b2xr(3)
    dimension:: a0xb0(3),a0xb1(3),a0xb2(3)
    dimension:: a1xb0(3),a1xb1(3),a1xb2(3)
    dimension:: a2xb0(3),a2xb1(3),a2xb2(3)
               !per = perpendicular, x = cross product 
    dimension:: rtestipm(3),rtestjpm(3)     !test  ????
    common /ftab/ fa(3),fb(3),ta(3),tb(3)   !output f,t for pair
    
    !common variables for Aaniso
    common /ab/a0(3),b0(3),a1(3),b1(3),a2(3),b2(3)  !orientation vector  
    common /angle/a0r,b0r, a1r,b1r, a2r,b2r, a0b0,a0b1,a0b2, a1b0,a1b1,a1b2, a2b0,a2b1,a2b2     
    common /Aaniso/Aaniso,dAaniso_da0r,dAaniso_da1r,dAaniso_da2r,dAaniso_db0r,dAaniso_db1r,dAaniso_db2r, &
                          dAaniso_da0b0, dAaniso_da0b1, dAaniso_da0b2, &
                          dAaniso_da1b0, dAaniso_da1b1, dAaniso_da1b2, &
                          dAaniso_da2b0, dAaniso_da2b1, dAaniso_da2b2
    
    rlstsq=rlist*rlist
    rcutsq=rcut*rcut    
    
    fp= 0.d0 ; tp=0.0d0
    energyI=0.d0
    
    if (update) then !save current configuration , construct neighbour list and calculate force)
                
        xp0=xp ! = (call save_verlet)
        nlist=0
        do 100 ipm=1,npm-1
            point(ipm)=nlist+1               
        do 100 jpm=ipm+1,npm           
            do k=1,2                            ! only x and y direction ON PLANE
            r(k)=xp(jpm,k)-xp(ipm,k)        ! from ipm to jpm !Feb
            enddo
            r(1)=r(1)-anint(r(1)/bl)*bl !if cross x
            r(1)=r(1)-anint(r(2)/bl)*bl*pbcx; r(2)=r(2)-anint(r(2)/bl)*bl !if cross y
            
            rr = r(1)*r(1) + r(2)*r(2)   !distance^2
            
            
            
            if (rr < rlstsq) then
                nlist=nlist+1
                list(nlist)=jpm
                !remove this check if maxnab is appropriate)
                if (nlist == maxnab) stop 'list too small'
                if (rr < rcutsq) then
                
                !CALCUALTE F AND T:------------------------
                    rr=sqrt(rr) !now rr is distance   #########
                                        
                    rhat=r/rr           !unit distance vector                     
                    a0(1:3)=p0(ipm,1:3);    b0(1:3)=p0(jpm,1:3)
                    a1(1:3)=p1(ipm,1:3);    b1(1:3)=p1(jpm,1:3)
                    a2(1:3)=p2(ipm,1:3);    b2(1:3)=p2(jpm,1:3)
                   !calculate cosine i.e. dot product of basis vector an.r, bn.r
                    a0r=0.d0; b0r=0.d0 
                    a1r=0.d0; b1r=0.d0
                    a2r=0.d0; b2r=0.d0
                    a0b0=0.d0; a0b1=0.d0; a0b2=0.d0    !Add on Jan3 2018
                    a1b0=0.d0; a1b1=0.d0; a1b2=0.d0    !Add on Jan3 2018
                    a2b0=0.d0; a2b1=0.d0; a2b2=0.d0    !Add on Jan3 2018
                     
                    do i=1,3 !cosa=a.rhat
                         a0r= a0r + a0(i)*rhat(i);  b0r= b0r + b0(i)*rhat(i)     
                         a1r= a1r + a1(i)*rhat(i);  b1r= b1r + b1(i)*rhat(i)
                         a2r= a2r + a2(i)*rhat(i);  b2r= b2r + b2(i)*rhat(i)                        
                         a0b0=a0b0 + a0(i)*b0(i)    !cosan.bn = an.bn
                         a0b1=a0b1 + a0(i)*b1(i)
                         a0b2=a0b2 + a0(i)*b2(i)
                         a1b0=a1b0 + a1(i)*b0(i)
                         a1b1=a1b1 + a1(i)*b1(i)
                         a1b2=a1b2 + a1(i)*b2(i)
                         a2b0=a2b0 + a2(i)*b0(i)
                         a2b1=a2b1 + a2(i)*b1(i)
                         a2b2=a2b2 + a2(i)*b2(i)
                    enddo  
                   
                    !---Calculate F on u 
                        !calculate prefactor
                        !vector !a perpendicular = an - (an.rhat)rhat
                        a0per=a0-a0r*rhat;  b0per=b0-b0r*rhat     
                        a1per=a1-a1r*rhat;  b1per=b1-b1r*rhat 
                        a2per=a2-a2r*rhat;  b2per=b2-b2r*rhat                                          
                        !an x r, bn x r, an x bm
                        a0xr = cross(a0,rhat);  b0xr = cross(b0,rhat)
                        a1xr = cross(a1,rhat);  b1xr = cross(b1,rhat)
                        a2xr = cross(a2,rhat);  b2xr = cross(b2,rhat)                             
                        a0xb0=cross(a0,b0); a0xb1=cross(a0,b1); a0xb2=cross(a0,b2)
                        a1xb0=cross(a1,b0); a1xb1=cross(a1,b1); a1xb2=cross(a1,b2)
                        a2xb0=cross(a2,b0); a2xb1=cross(a2,b1); a2xb2=cross(a2,b2)
                             
                        !Potential 
                        !Uiso = U_WCA  
                        sigma12 =(sigma/rr)**12 ! rr is length distance
                        sigma6  =(sigma/rr)**6        
                        Uiso    = 4.0d0*(sigma12-sigma6+0.25d0)            !FOR rr <= 2^(1/6)   ! +0.25 due to WCA   
                        dUiso_dr= -4.0d0*6.0d0*(2.0d0*sigma12-sigma6)/rr   !FOR rr <= 2^(1/6)
                        
                        rrlim=2.D0*RADIUS*2.d0**(1.d0/6.d0)                     !FOR rr > 2^(1/6)
                        if (rr > rrlim) then                              
                            Uiso = 0.d0
                            dUiso_dr = 0.d0
                        endif 
                               
                        !Uaniso = U_Morse (refer to our paper)             
                        req = 1.878d0   !pair equilibrium distance
                        Uaniso = cmd*( (1.d0-exp(-(rr-req)/cmr))**2 - 1.d0 )
                        dUaniso_dr = cmd*2.d0/cmr*(1.d0-exp((-rr+req)/cmr))*exp((-rr+req)/cmr)                                                                                
                                                       
                        !Calculate Anisotropic interaction function and its derivatives              
                        call calAaniso2 
                        energyI(ipm)=energyI(ipm)+ Uiso + Aaniso*Uaniso
                        energyI(jpm)=energyI(jpm)+ Uiso + Aaniso*Uaniso !add 201904
                        
                        do 10 i=1,3
                        fa(i)=   ( (dUiso_dr + Aaniso* dUaniso_dr)*rhat(i) + &
                                   1/rr*Uaniso*(dAaniso_da0r* a0per(i) + dAaniso_da1r* a1per(i) + dAaniso_da2r* a2per(i)+ &
                                                dAaniso_db0r* b0per(i) + dAaniso_db1r* b1per(i) + dAaniso_db2r* b2per(i)   ) )
                        fb(i)= -fa(i)
                         
                        ta(i)= Uaniso* ( (-1.d0)*(dAaniso_da0r * a0xr(i) +  dAaniso_da1r * a1xr(i) + dAaniso_da2r * a2xr(i)) + &
                                          (-1.d0)*(dAaniso_da0b0 * a0xb0(i) +  dAaniso_da0b1 * a0xb1(i) + dAaniso_da0b2 * a0xb2(i) +&
                                                  dAaniso_da1b0 * a1xb0(i) +  dAaniso_da1b1 * a1xb1(i) + dAaniso_da1b2 * a1xb2(i) +&
                                                  dAaniso_da2b0 * a2xb0(i) +  dAaniso_da2b1 * a2xb1(i) + dAaniso_da2b2 * a2xb2(i)) )
                        tb(i)= Uaniso* ( (-1.d0)*(dAaniso_db0r * b0xr(i) +  dAaniso_db1r * b1xr(i) + dAaniso_db2r * b2xr(i)) + &
                                                  (dAaniso_da0b0 * a0xb0(i) +  dAaniso_da0b1 * a0xb1(i) + dAaniso_da0b2 * a0xb2(i) +&
                                                  dAaniso_da1b0 * a1xb0(i) +  dAaniso_da1b1 * a1xb1(i) + dAaniso_da1b2 * a1xb2(i) +&
                                                  dAaniso_da2b0 * a2xb0(i) +  dAaniso_da2b1 * a2xb1(i) + dAaniso_da2b2 * a2xb2(i)) )
                    10  continue   
        
                     !give back value from subroutine pair(ipm,jpm)
                     fp(ipm,1:3)=fp(ipm,1:3) + fa(1:3)  !!!!
                     fp(jpm,1:3)=fp(jpm,1:3) + fb(1:3)
                     tp(ipm,1:3)=tp(ipm,1:3) + ta(1:3)
                     tp(jpm,1:3)=tp(jpm,1:3) + tb(1:3)
                endif
            endif
                    
100     continue 

        aveEnergyI=0.d0
        do ipm=1,npm
        aveEnergyI=aveEnergyI+energyI(ipm)
        enddo
        aveEnergyI=aveEnergyI/dble(npm)
                  
        point(npm)=nlist+1
            
    else !use the list to find the neighbours:
        do 200 ipm=1,npm-1
        jbeg=point(ipm)
        jend=point(ipm+1)-1
            
        !check that ipm has neighbours:
        if (jbeg <= jend) then
        do 210 jnab=jbeg,jend
            jpm=list(jnab)
                
            do k=1,2                            ! only x and y direction ON PLANE
            r(k)=xp(jpm,k)-xp(ipm,k)        ! from ipm to jpm !Feb
            enddo
            r(1)=r(1)-anint(r(1)/bl)*bl !if cross x
            r(1)=r(1)-anint(r(2)/bl)*bl*pbcx; r(2)=r(2)-anint(r(2)/bl)*bl !if cross y
            
            rr = r(1)*r(1) + r(2)*r(2)          !distance^2
            if (rr < rcutsq) then
                
            !CALCUALTE F AND T:------------------------
                rr=sqrt(rr) !now rr is distance   #########
                    
                rhat=r/rr           !unit distance vector                
                a0(1:3)=p0(ipm,1:3);    b0(1:3)=p0(jpm,1:3)
                a1(1:3)=p1(ipm,1:3);    b1(1:3)=p1(jpm,1:3)
                a2(1:3)=p2(ipm,1:3);    b2(1:3)=p2(jpm,1:3)                                    
                !calculate cosine i.e. dot product of basis vector an.r, bn.r
                a0r=0.d0; b0r=0.d0 
                a1r=0.d0; b1r=0.d0
                a2r=0.d0; b2r=0.d0
                a0b0=0.d0; a0b1=0.d0; a0b2=0.d0    !Add on Jan3 2018
                a1b0=0.d0; a1b1=0.d0; a1b2=0.d0    !Add on Jan3 2018
                a2b0=0.d0; a2b1=0.d0; a2b2=0.d0    !Add on Jan3 2018
                 
                do  i=1,3 !cosa=a.rhat
                    a0r= a0r + a0(i)*rhat(i);  b0r= b0r + b0(i)*rhat(i)    
                    a1r= a1r + a1(i)*rhat(i);  b1r= b1r + b1(i)*rhat(i)
                    a2r= a2r + a2(i)*rhat(i);  b2r= b2r + b2(i)*rhat(i)                   
                    a0b0=a0b0 + a0(i)*b0(i)    !cosan.bn = an.bn
                    a0b1=a0b1 + a0(i)*b1(i)
                    a0b2=a0b2 + a0(i)*b2(i)
                    a1b0=a1b0 + a1(i)*b0(i)
                    a1b1=a1b1 + a1(i)*b1(i)
                    a1b2=a1b2 + a1(i)*b2(i)
                    a2b0=a2b0 + a2(i)*b0(i)
                    a2b1=a2b1 + a2(i)*b1(i)
                    a2b2=a2b2 + a2(i)*b2(i)
                enddo  
                    
                !---Calculate F on u 
                !calculate prefactor
                    !vector !a perpendicular = an - (an.rhat)rhat
                    a0per=a0-a0r*rhat;  b0per=b0-b0r*rhat    
                    a1per=a1-a1r*rhat;  b1per=b1-b1r*rhat 
                    a2per=a2-a2r*rhat;  b2per=b2-b2r*rhat                                
                    !an x r, bn x r, an x bm
                    a0xr = cross(a0,rhat);  b0xr = cross(b0,rhat)
                    a1xr = cross(a1,rhat);  b1xr = cross(b1,rhat)
                    a2xr = cross(a2,rhat);  b2xr = cross(b2,rhat)                         
                    a0xb0=cross(a0,b0); a0xb1=cross(a0,b1); a0xb2=cross(a0,b2)
                    a1xb0=cross(a1,b0); a1xb1=cross(a1,b1); a1xb2=cross(a1,b2)
                    a2xb0=cross(a2,b0); a2xb1=cross(a2,b1); a2xb2=cross(a2,b2)
                     
                !Potential 
                    !Uiso = U_WCA  
                    sigma12 =(sigma/rr)**12 ! rr is length distance
                    sigma6  =(sigma/rr)**6        
                    Uiso    = 4.0d0*(sigma12-sigma6+0.25d0)            !FOR rr <= 2^(1/6)   ! +0.25 due to WCA   
                    dUiso_dr= -4.0d0*6.0d0*(2.0d0*sigma12-sigma6)/rr   !FOR rr <= 2^(1/6)
                    
                    rrlim=2.D0*RADIUS*2.d0**(1.d0/6.d0)                     !FOR rr > 2^(1/6)
                    if (rr > rrlim) then                              
                        Uiso = 0.d0
                        dUiso_dr = 0.d0
                    endif 
                           
                    !Uaniso = U_Morse (refer to Bianchi's paper)             
                    req = 1.878d0   !pair equilibrium distance
                    Uaniso = cmd*( (1.d0-exp(-(rr-req)/cmr))**2 - 1.d0 )
                    dUaniso_dr = cmd*2.d0/cmr*(1.d0-exp((-rr+req)/cmr))*exp((-rr+req)/cmr)                                                                                
                       
                !Calculate Anisotropic interaction function and its derivatives              
                    call calAaniso2 
                    energyI(ipm)=energyI(ipm)+ Uiso + Aaniso*Uaniso
                    energyI(jpm)=energyI(jpm)+ Uiso + Aaniso*Uaniso !add 201904
                    
                    do 20 i=1,3
                    fa(i)=   ( (dUiso_dr + Aaniso* dUaniso_dr)*rhat(i) + &
                               1/rr*Uaniso*(dAaniso_da0r* a0per(i) + dAaniso_da1r* a1per(i) + dAaniso_da2r* a2per(i)+ &
                                            dAaniso_db0r* b0per(i) + dAaniso_db1r* b1per(i) + dAaniso_db2r* b2per(i)   ) )
                    fb(i)= -fa(i)
                     
                    ta(i)= Uaniso* ( (-1.d0)*(dAaniso_da0r * a0xr(i) +  dAaniso_da1r * a1xr(i) + dAaniso_da2r * a2xr(i)) + &
                                      (-1.d0)*(dAaniso_da0b0 * a0xb0(i) +  dAaniso_da0b1 * a0xb1(i) + dAaniso_da0b2 * a0xb2(i) +&
                                              dAaniso_da1b0 * a1xb0(i) +  dAaniso_da1b1 * a1xb1(i) + dAaniso_da1b2 * a1xb2(i) +&
                                              dAaniso_da2b0 * a2xb0(i) +  dAaniso_da2b1 * a2xb1(i) + dAaniso_da2b2 * a2xb2(i)) )
                    tb(i)= Uaniso* ( (-1.d0)*(dAaniso_db0r * b0xr(i) +  dAaniso_db1r * b1xr(i) + dAaniso_db2r * b2xr(i)) + &
                                              (dAaniso_da0b0 * a0xb0(i) +  dAaniso_da0b1 * a0xb1(i) + dAaniso_da0b2 * a0xb2(i) +&
                                              dAaniso_da1b0 * a1xb0(i) +  dAaniso_da1b1 * a1xb1(i) + dAaniso_da1b2 * a1xb2(i) +&
                                              dAaniso_da2b0 * a2xb0(i) +  dAaniso_da2b1 * a2xb1(i) + dAaniso_da2b2 * a2xb2(i)) )
                20  continue   
        
                !give back value from subroutine pair(ipm,jpm)
                fp(ipm,1:3)=fp(ipm,1:3) + fa(1:3)  !!!!
                fp(jpm,1:3)=fp(jpm,1:3) + fb(1:3)
                tp(ipm,1:3)=tp(ipm,1:3) + ta(1:3)
                tp(jpm,1:3)=tp(jpm,1:3) + tb(1:3)
            endif
210     continue
        endif
                                 
200     continue 
        
        aveEnergyI=0.d0
        do ipm=1,npm
        aveEnergyI=aveEnergyI+energyI(ipm)
        enddo
        aveEnergyI=aveEnergyI/dble(npm)           
    
    endif
    end subroutine force_verlet
    
!---------------------------------------------------
    subroutine check
    use para
    
    implicit double precision (a-h,o-z)
    
    if (istep ==1) then
        dispmx=100.d0 !a large value
    else     
        !calculate maximum displacement since last update
        dispmx=0.d0
        do 10 j=1,2 ! for in 2d
        do 10 ipm=1,npm
        dispmx=max(abs(xp(ipm,j)-xp0(ipm,j)),dispmx)
        10 continue
        
        !a conservative test of the list skin crossing (?)
        dispmx=2.d0*sqrt(2.d0*dispmx*dispmx) !the 1st 2.d0 is by definitaion 2.d0 in sqrt is for 2d case)
   
    endif
  
    update= ( dispmx > (rlist-rcut))
    end subroutine check
    
!-----------------------------------------------    
    

        
        



