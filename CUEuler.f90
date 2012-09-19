! Author of the program Vladimir Fuka
! Modified by
! Donato Cecere
! 1 July 2007 

module CUeuler
    implicit none

    integer,parameter:: KND=KIND(1.0D0)

    integer RKtype   !(1 - 1. order ... 3 - 3. order, 4. fully discrete)
    integer FLUXtype !(1 - KURGTADMOR, 2 - CENTUPW, 3 - LAXFRIED, 4 - LAXWENDROFF,5 - MUSTA)
    integer MUSTAstages!(Number of stages for MUSTA flux, for 1 it reduces to FORCE flux)
    integer RECtype  !(1 - piecewise constant,  2 - piecewise linear)
    integer limitertype !(1 - minmod, 2 - modified minmod (MC), 3 - modified superbee)
    real(KND) limparam !(parameter for limiter - only if needed)
    integer extype   !(1 - shock tube test, 2 - smooth test case with periodic BC)
    integer outputtype !(1 - single frame textfiles for gnuplot, 2- one Tecplot file) Tecplot output by Donato Cecere, beta
    integer nx,nt
    
    real(KND) time,dx,xl,xr,dt,tend,CFL,rhoLEFT,rhouLEFT,eLEFT,rhoRIGHT,rhouRIGHT,eRIGHT

    real(KND),parameter:: gamma=1.4_KND, pi=3.1415926535_KND
    
    interface MINMOD
         module procedure MINMOD2, MINMOD3
    end interface
    
    contains
    
    subroutine TMARCH(rho,rhou,e,x)
     integer i
     real(KND),dimension(-2:):: rho,rhou,e,x
	 !real(KND),dimension(-2:nx+3):: rho,rhou,e
	 real(KND),dimension(LBOUND(rho,1):UBOUND(rho,1)):: rho2,rhou2,e2,drho,drhou,de

     if (RKtype==1) then
      call DRK(drho,drhou,de,rho,rhou,e,x)
      rho=rho + dt*drho
      rhou=rhou + dt*drhou
      e=e + dt*de
      
     elseif (RKtype==2) then
      call DRK(drho,drhou,de,rho,rhou,e,x)
      rho2=rho + dt*drho
      rhou2=rhou + dt*drhou
      e2=e + dt*de

      call DRK(drho,drhou,de,rho2,rhou2,e2,x)
      rho=(1._KND/2._KND)*rho + (1._KND/2._KND)*rho2 + (1._KND/2._KND)*dt*drho
      rhou=(1._KND/2._KND)*rhou + (1._KND/2._KND)*rhou2 + (1._KND/2._KND)*dt*drhou
      e=(1._KND/2._KND)*e + (1._KND/2._KND)*e2 + (1._KND/2._KND)*dt*de

     elseif (RKtype==3) then

      call DRK(drho,drhou,de,rho,rhou,e,x) !1
      rho2=rho + dt*drho
      rhou2=rhou + dt*drhou
      e2=e + dt*de
     
      call DRK(drho,drhou,de,rho2,rhou2,e2,x) !2
      rho2=(3._KND/4._KND)*rho + (1._KND/4._KND)*rho2 + (1._KND/4._KND)*dt*drho
      rhou2=(3._KND/4._KND)*rhou + (1._KND/4._KND)*rhou2 + (1._KND/4._KND)*dt*drhou
      e2=(3._KND/4._KND)*e + (1._KND/4._KND)*e2 + (1._KND/4._KND)*dt*de
     
      call DRK(drho,drhou,de,rho2,rhou2,e2,x) !3
      rho=(1._KND/3._KND)*rho + (2._KND/3._KND)*rho2 + (2._KND/3._KND)*dt*drho
      rhou=(1._KND/3._KND)*rhou + (2._KND/3._KND)*rhou2 + (2._KND/3._KND)*dt*drhou
      e=(1._KND/3._KND)*e + (2._KND/3._KND)*e2 + (2._KND/3._KND)*dt*de
     
     elseif (RKtype==4) then
      call DFD(rho,rhou,e)
     endif    
    endsubroutine TMARCH    

    
    subroutine DRK(drho,drhou,de,rho,rhou,e,x)
     real(KND),dimension(-2:):: rho,rhou,e,drho,drhou,de,x
     real(KND),allocatable,dimension(:):: rhoL,rhouL,eL,rhoR,rhouR,eR
     drho=0
     drhou=0
     de=0

     allocate(rhoL(-2:nx+2))
     allocate(rhouL(-2:nx+2))
     allocate(eL(-2:nx+2))
    
     allocate(rhoR(-2:nx+2))
     allocate(rhouR(-2:nx+2))
     allocate(eR(-2:nx+2))
     
     
     call BOUND(rho,rhou,e) !boundary condition
     if (RECtype==1) then
         call PCONST(rhoL,rhouL,eL,rhoR,rhouR,eR,rho,rhou,e)
     elseif (RECtype==2) then
         call MUSCL(rhoL,rhouL,eL,rhoR,rhouR,eR,rho,rhou,e,x)
     endif

     if (FLUXtype==1) then
         call FLUXESKT(drho,drhou,de,rhoL,rhouL,eL,rhoR,rhouR,eR,x)
     elseif (FLUXtype==2) then
         call FLUXESCU(drho,drhou,de,rhoL,rhouL,eL,rhoR,rhouR,eR)
     elseif (FLUXtype==3) then
         call FLUXESLF(drho,drhou,de,rhoL,rhouL,eL,rhoR,rhouR,eR)
     elseif (FLUXtype==4) then
         call FLUXESLW(drho,drhou,de,rhoL,rhouL,eL,rhoR,rhouR,eR)
     elseif (FLUXtype==5) then
         call FLUXESMUSTA(drho,drhou,de,rhoL,rhouL,eL,rhoR,rhouR,eR)
     elseif (FLUXtype==6) then
         call FLUXESCULD(drho,drhou,de,rhoL,rhouL,eL,rhoR,rhouR,eR)
     endif

     deallocate(rhoL,rhoR,rhouL,rhouR,eL,eR)
    endsubroutine DRK

    subroutine DFD(rho,rhou,e)
     real(KND),dimension(-2:):: rho,rhou,e

     if (FLUXtype==3) then
         call FDLF(rho,rhou,e)
     elseif (FLUXtype==4) then
         call FDLW(rho,rhou,e)
     endif
    endsubroutine DFD
    

! This procedure calculates  extrapolated values  of the dependent variables using 
! a central difference scheme: k = 1, see Tannehill-Anderson-Pletcher, pp. 204-205.
! The values for  the  extrapolated variables  on both sides  of the interface are 
! calculated on  the left of i point. 

    subroutine MUSCL(rhoL,rhouL,eL,rhoR,rhouR,eR,rho,rhou,e,x)
     real(KND),dimension(-2:)::rhoL,rhouL,eL,rhoR,rhouR,eR,rho,rhou,e,x
     real(KND) sL,sC,sR
     integer i

        do i=0,nx+1
!         sL=(rho(i)-rho(i-1))
         sL=(rho(i)-rho(i-1))/(x(i)-x(i-1))
!         sR=(rho(i+1)-rho(i))
         sR=(rho(i+1)-rho(i))/(x(i+1)-x(i))
         sL=sL
         sR=sR
         sC=LIMITER(sL,sR)
         rhoR(i)=rho(i)-sC*(0.5_KND)*(x(i)-x(i-1))
         rhoL(i+1)=rho(i)+sC*(0.5_KND)*(x(i+1)-x(i))
        enddo
        do i=0,nx+1
!         sL=(rhou(i)-rhou(i-1))
         sL=(rhou(i)-rhou(i-1))/(x(i)-x(i-1))         
!         sR=(rhou(i+1)-rhou(i))
         sR=(rhou(i+1)-rhou(i))/(x(i+1)-x(i))
         sL=sL
         sR=sR
         sC=LIMITER(sL,sR)
         rhouR(i)=rhou(i)-sC*(0.5_KND)*(x(i)-x(i-1))
         rhouL(i+1)=rhou(i)+sC*(0.5_KND)*(x(i+1)-x(i))
        enddo
        do i=0,nx+1
!         sL=(e(i)-e(i-1))
         sL=(e(i)-e(i-1))/(x(i)-x(i-1))
!         sR=(e(i+1)-e(i))
         sR=(e(i+1)-e(i))/(x(i+1)-x(i))
         sL=sL
         sR=sR
         sC=LIMITER(sL,sR)
         eR(i)=e(i)-sC*(0.5_KND)*(x(i)-x(i-1))
         eL(i+1)=e(i)+sC*(0.5_KND)*(x(i+1)-x(i))
        enddo
   endsubroutine MUSCL
! Linear Piecewise extrapolation   
   subroutine PCONST(rhoL,rhouL,eL,rhoR,rhouR,eR,rho,rhou,e)
     real(KND),dimension(-2:)::rhoL,rhouL,eL,rhoR,rhouR,eR,rho,rhou,e

     integer i

        do i=0,nx+1
         rhoR(i)=rho(i)
         rhoL(i+1)=rho(i)
        enddo
        do i=0,nx+1
         rhouR(i)=rhou(i)
         rhouL(i+1)=rhou(i)
        enddo
        do i=0,nx+1
         eR(i)=e(i)
         eL(i+1)=e(i)
        enddo
   endsubroutine PCONST
   

   
   
   subroutine FLUXESLF(drho,drhou,de,rhoL,rhouL,eL,rhoR,rhouR,eR)
    real(KND),dimension(-2:)::rhoL,rhouL,eL,rhoR,rhouR,eR
    real(KND),dimension(LBOUND(rhoL,1):UBOUND(rhoL,1))::drho,drhou,de,pL,pR
    real(KND) p,A
    integer i

    do i=1,nx+1
        pL(i)=(gamma-1._KND)*(eL(i)-(1._KND/2._KND)*(rhouL(i)*rhouL(i))/rhoL(i))
        pR(i)=(gamma-1._KND)*(eR(i)-(1._KND/2._KND)*(rhouR(i)*rhouR(i))/rhoR(i))
        if ((pL(i)<0).or.(pR(i)<0)) write(*,*) "ERROR! pressure < 0",pL(i),pR(i)
    enddo

   A=dx/(dt)
    !
    !flux of rho
    !
   drho(1)=-(frho(rhoL(1),rhouL(1),eL(1),pL(1))+frho(rhoR(1),rhouR(1),eR(1),pR(1))&
             -A*(rhoR(1)-rhoL(1)))/dx

     do i=2,nx
      p=frho(rhoL(i),rhouL(i),eL(i),pL(i))+frho(rhoR(i),rhouR(i),eR(i),pR(i))&
         -A*(rhoR(i)-rhoL(i))
      drho(i-1)=drho(i-1)+p/dx
      drho(i)=drho(i)-p/dx
     enddo
    drho(nx)=drho(nx)+(frho(rhoL(nx+1),rhouL(nx+1),eL(nx+1),pL(1))&
                     +frho(rhoR(nx+1),rhouR(nx+1),eR(nx+1),pR(1))&
         -A*(rhoR(nx+1)-rhoL(nx+1)))/dx

    drho=-drho/2._KND

    !
    !flux of rho*u
    !
    drhou(1)=-(frhou(rhoL(1),rhouL(1),eL(1),pL(1))+frhou(rhoR(1),rhouR(1),eR(1),pR(1))&
         -A*(rhouR(1)-rhouL(1)))/dx
    do i=2,nx
      p=frhou(rhoL(i),rhouL(i),eL(i),pL(i))+frhou(rhoR(i),rhouR(i),eR(i),pR(i))&
         -A*(rhouR(i)-rhouL(i))
      drhou(i-1)=drhou(i-1)+p/dx
      drhou(i)=drhou(i)-p/dx
    enddo
   drhou(nx)=drhou(nx)+(frhou(rhoL(nx+1),rhouL(nx+1),eL(nx+1),pL(nx+1))&
                       +frhou(rhoR(nx+1),rhouR(nx+1),eR(nx+1),pR(nx+1))&
         -A*(rhouR(nx+1)-rhouL(nx+1)))/dx

    drhou=-drhou/2._KND


    !
    !flux ef energy
    !
    de(1)=-(fe(rhoL(1),rhouL(1),eL(1),pL(1))+fe(rhoR(1),rhouR(1),eR(1),pR(1))&
         -A*(eR(1)-eL(1)))/dx
    do i=2,nx
      p=fe(rhoL(i),rhouL(i),eL(i),pL(i))+fe(rhoR(i),rhouR(i),eR(i),pR(i))&
         -A*(eR(i)-eL(i))
      de(i-1)=de(i-1)+p/dx
      de(i)=de(i)-p/dx
    enddo
    de(nx)=de(nx)+(fe(rhoL(nx+1),rhouL(nx+1),eL(nx+1),pL(nx+1))&
                  +fe(rhoR(nx+1),rhouR(nx+1),eR(nx+1),pR(nx+1))&
         -A*(eR(nx+1)-eL(nx+1)))/dx

    de=-de/2._KND
   endsubroutine FLUXESLF

   subroutine FLUXESLW(drho,drhou,de,rhoL,rhouL,eL,rhoR,rhouR,eR)
    real(KND),dimension(-2:)::rhoL,rhouL,eL,rhoR,rhouR,eR
    real(KND),dimension(LBOUND(rhoL,1):UBOUND(rhoL,1))::drho,drhou,de,flrho,flrhou,fle
    real(KND) A,pL,pR,pM
    real(KND),DIMENSION(1:3):: qL,qR,qM,fl,fr
    integer i
    real(KND),PARAMETER:: eps=0.0001

    A=dx/dt
    do i=1,nx+1
      qL=(/ rhoL(i),rhouL(i),eL(i) /)
      qR=(/ rhoR(i),rhouR(i),eR(i) /)
      pL=(gamma-1._KND)*(eL(i)-0.5_KND*(rhouL(i)*rhouL(i))/rhoL(i))
      pR=(gamma-1._KND)*(eR(i)-0.5_KND*(rhouR(i)*rhouR(i))/rhoR(i))
      if (pL<0)  pL=0
      if (pR<0)  pR=0

      fl(1)=frho(qL(1),qL(2),qL(3),pL)
      fl(2)=frhou(qL(1),qL(2),qL(3),pL)
      fl(3)=fe(qL(1),qL(2),qL(3),pL)
      fr(1)=frho(qR(1),qR(2),qR(3),pR)
      fr(2)=frhou(qR(1),qR(2),qR(3),pR)
      fr(3)=fe(qR(1),qR(2),qR(3),pR)

      qM=0.5_KND*(qL+qR)-(0.5_KND/A)*(fr-fl)
      if (qM(1)<eps)  qM(1)=eps
      if (qM(3)<eps)  qM(3)=eps
      pM=(gamma-1._KND)*(qM(3)-0.5_KND*(qM(2)*qM(2))/qM(1))
      if (pM<0)  pM=0
      
      flrho(i)=frho(qM(1),qM(2),qM(3),pM)
      flrhou(i)=frhou(qM(1),qM(2),qM(3),pM)
      fle(i)=fe(qM(1),qM(2),qM(3),pM)
    enddo
    
    
    drho(1)=flrho(1)/dx
    drhou(1)=flrhou(1)/dx
    de(1)=fle(1)/dx
    do i=2,nx
     drho(i-1)=drho(i-1)-flrho(i)/dx
     drhou(i-1)=drhou(i-1)-flrhou(i)/dx
     de(i-1)=de(i-1)-fle(i)/dx
     drho(i)=flrho(i)/dx
     drhou(i)=flrhou(i)/dx
     de(i)=fle(i)/dx
    enddo
    drho(nx)=drho(nx)-flrho(nx+1)/dx
    drhou(nx)=drhou(nx)-flrhou(nx+1)/dx
    de(nx)=de(nx)-fle(nx+1)/dx
    
   endsubroutine FLUXESLW


   subroutine FLUXESKT(drho,drhou,de,rhoL,rhouL,eL,rhoR,rhouR,eR,x) !!!!!!!!!!!!!!!!!!!!KURGANOV TADMOR SCHEME !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(KND),dimension(-2:)::rhoL,rhouL,eL,rhoR,rhouR,eR,x
    real(KND),dimension(LBOUND(rhoL,1):UBOUND(rhoL,1))::drho,drhou,de,A,pL,pR
    real(KND) p,cL,cR,dxnotu
    integer i

    do i=1,nx+1
        pL(i)=(gamma-1._KND)*(eL(i)-(1._KND/2._KND)*(rhouL(i)*rhouL(i))/rhoL(i))
        pR(i)=(gamma-1._KND)*(eR(i)-(1._KND/2._KND)*(rhouR(i)*rhouR(i))/rhoR(i))
        if ((pL(i)<0).or.(pR(i)<0)) write(*,*) "ERROR! pressure < 0",pL(i),pR(i), &
         eL(i),(1._KND/2._KND)*(rhouL(i)*rhouL(i))/rhoL(i),eR(i),(1._KND/2._KND)*(rhouR(i)*rhouR(i))/rhoR(i)

        cL=SQRT(gamma*pL(i)/rhoL(i))
        cR=SQRT(gamma*pR(i)/rhoR(i))
        A(i)=MAX(ABS(rhouL(1)/rhoL(1)-cL),ABS(rhouR(1)/rhoR(1)-cR),&
              ABS(rhouL(1)/rhoL(1)+cL),ABS(rhouR(1)/rhoR(1)+cR))
    enddo
    
	drho  = 0.d0 !PER CHIAREZZA
    drhou = 0.d0
	de    = 0.d0
    
	
	!    
    !flux of rho
    ! 
!     drho(1)=-(frho(rhoL(1),rhouL(1),eL(1),pL(1))+frho(rhoR(1),rhouR(1),eR(1),pR(1))&
!             -A(1)*(rhoR(1)-rhoL(1)))/dx
     drho(1)=-(frho(rhoL(1),rhouL(1),eL(1),pL(1))+frho(rhoR(1),rhouR(1),eR(1),pR(1))&
             -A(1)*(rhoR(1)-rhoL(1)))



     do i=2,nx
      p=frho(rhoL(i),rhouL(i),eL(i),pL(i))+frho(rhoR(i),rhouR(i),eR(i),pR(i))&
         -A(i)*(rhoR(i)-rhoL(i))
!      drho(i-1)=drho(i-1)+p/dx
      drho(i-1)=drho(i-1)+p
!      drho(i)=drho(i)-p/dx
      drho(i)=drho(i)-p
     enddo
!    drho(nx)=drho(nx)+(frho(rhoL(nx+1),rhouL(nx+1),eL(nx+1),pL(1))&
!                     +frho(rhoR(nx+1),rhouR(nx+1),eR(nx+1),pR(1))&
!         -A(nx+1)*(rhoR(nx+1)-rhoL(nx+1)))/dx
!    drho=-drho/2._KND
    drho(nx)=drho(nx)+(frho(rhoL(nx+1),rhouL(nx+1),eL(nx+1),pL(1))&
                     +frho(rhoR(nx+1),rhouR(nx+1),eR(nx+1),pR(1))&
         -A(nx+1)*(rhoR(nx+1)-rhoL(nx+1)))

     do i=1,nx  
	  dxnotu = (x(i+1)-x(i-1))
      drho(i) = -drho(i)/dxnotu
     enddo



    !
    !flux of rho*u 
    !
!    drhou(1)=-(frhou(rhoL(1),rhouL(1),eL(1),pL(1))+frhou(rhoR(1),rhouR(1),eR(1),pR(1))&
!         -A(1)*(rhouR(1)-rhouL(1)))/dx
    drhou(1)=-(frhou(rhoL(1),rhouL(1),eL(1),pL(1))+frhou(rhoR(1),rhouR(1),eR(1),pR(1))&
         -A(1)*(rhouR(1)-rhouL(1)))
    do i=2,nx
      p=frhou(rhoL(i),rhouL(i),eL(i),pL(i))+frhou(rhoR(i),rhouR(i),eR(i),pR(i))&
         -A(i)*(rhouR(i)-rhouL(i))
!      drhou(i-1)=drhou(i-1)+p/dx
      drhou(i-1)=drhou(i-1)+p
!      drhou(i)=drhou(i)-p/dx
      drhou(i)=drhou(i)-p
    enddo
!    drhou(nx)=drhou(nx)+(frhou(rhoL(nx+1),rhouL(nx+1),eL(nx+1),pL(nx+1))&
!                       +frhou(rhoR(nx+1),rhouR(nx+1),eR(nx+1),pR(nx+1))&
!         -A(nx+1)*(rhouR(nx+1)-rhouL(nx+1)))/dx
!    drhou=-drhou/2._KND   
    drhou(nx)=drhou(nx)+(frhou(rhoL(nx+1),rhouL(nx+1),eL(nx+1),pL(nx+1))&
                       +frhou(rhoR(nx+1),rhouR(nx+1),eR(nx+1),pR(nx+1))&
         -A(nx+1)*(rhouR(nx+1)-rhouL(nx+1)))
    
	 do i=1,nx  
	  dxnotu = (x(i+1)-x(i-1))
      drhou(i) = -drhou(i)/dxnotu
     enddo

    



    
    !
    !flux ef energy  
    !
!    de(1)=-(fe(rhoL(1),rhouL(1),eL(1),pL(1))+fe(rhoR(1),rhouR(1),eR(1),pR(1))&
!         -A(1)*(eR(1)-eL(1)))/dx
    de(1)=-(fe(rhoL(1),rhouL(1),eL(1),pL(1))+fe(rhoR(1),rhouR(1),eR(1),pR(1))&
         -A(1)*(eR(1)-eL(1)))
    do i=2,nx
      p=fe(rhoL(i),rhouL(i),eL(i),pL(i))+fe(rhoR(i),rhouR(i),eR(i),pR(i))&
         -A(i)*(eR(i)-eL(i))
!      de(i-1)=de(i-1)+p/dx
!      de(i)=de(i)-p/dx
      de(i-1)=de(i-1)+p
      de(i)=de(i)-p
    enddo
!    de(nx)=de(nx)+(fe(rhoL(nx+1),rhouL(nx+1),eL(nx+1),pL(nx+1))&
!                  +fe(rhoR(nx+1),rhouR(nx+1),eR(nx+1),pR(nx+1))&
!         -A(nx+1)*(eR(nx+1)-eL(nx+1)))/dx      
!    de=-de/2._KND
    de(nx)=de(nx)+(fe(rhoL(nx+1),rhouL(nx+1),eL(nx+1),pL(nx+1))&
                  +fe(rhoR(nx+1),rhouR(nx+1),eR(nx+1),pR(nx+1))&
         -A(nx+1)*(eR(nx+1)-eL(nx+1)))     
	 do i=1,nx  
	  dxnotu = (x(i+1)-x(i-1))
      de(i) = -de(i)/dxnotu
     enddo
   endsubroutine FLUXESKT
   
   
   subroutine FLUXESCU(drho,drhou,de,rhoL,rhouL,eL,rhoR,rhouR,eR)
    real(KND),dimension(-2:)::rhoL,rhouL,eL,rhoR,rhouR,eR
    real(KND),dimension(LBOUND(rhoL,1):UBOUND(rhoL,1))::drho,drhou,de,aL,aR,pL,pR
    real(KND) cL,cR,p
    integer i
    
    do i=1,nx+1
        pL(i)=(gamma-1._KND)*(eL(i)-(1._KND/2._KND)*(rhouL(i)*rhouL(i))/rhoL(i))
        pR(i)=(gamma-1._KND)*(eR(i)-(1._KND/2._KND)*(rhouR(i)*rhouR(i))/rhoR(i))
        if ((pL(i)<0).or.(pR(i)<0)) then
          write(*,*) "ERROR! pressure < 0",pL(i),pR(i)
        endif
        cL=SQRT(gamma*pL(i)/rhoL(i))
        cR=SQRT(gamma*pR(i)/rhoR(i))
        aL(i)=MIN(rhouL(i)/rhoL(i)-cL,rhouR(i)/rhoR(i)-cR,0._KND)
        aR(i)=MAX(rhouL(i)/rhoL(i)+cL,rhouR(i)/rhoR(i)+cR,0._KND)        
    enddo

    !    
    !flux of rho
    ! 

   drho(1)=-((aR(1)*frho(rhoL(1),rhouL(1),eL(1),pL(1))&
             -aL(1)*frho(rhoR(1),rhouR(1),eR(1),pR(1)))/(aR(1)-aL(1))&
             +(aL(1)*aR(1)/(aR(1)-aL(1)))*(rhoR(1)-rhoL(1)))/dx

     do i=2,nx
      p=(aR(i)*frho(rhoL(i),rhouL(i),eL(i),pL(i))&
            -aL(i)*frho(rhoR(i),rhouR(i),eR(i),pR(i)))/(aR(i)-aL(i))&
            +(aL(i)*aR(i)/(aR(i)-aL(i)))*(rhoR(i)-rhoL(i))
      drho(i-1)=drho(i-1)+p/dx
      drho(i)=drho(i)-p/dx
     enddo
    drho(nx)=drho(nx)+((aR(nx+1)*frho(rhoL(nx+1),rhouL(nx+1),eL(nx+1),pL(nx+1))&
                       -aL(nx+1)*frho(rhoR(nx+1),rhouR(nx+1),eR(nx+1),pR(nx+1)))/(aR(nx+1)-aL(nx+1))&
                       +(aL(nx+1)*aR(nx+1)/(aR(nx+1)-aL(nx+1)))*(rhoR(nx+1)-rhoL(nx+1)))/dx
      
    drho=-drho

    !
    !flux of rho*u 
    !
   drhou(1)=-((aR(1)*frhou(rhoL(1),rhouL(1),eL(1),pL(1))&
             -aL(1)*frhou(rhoR(1),rhouR(1),eR(1),pR(1)))/(aR(1)-aL(1))&
             +(aL(1)*aR(1)/(aR(1)-aL(1)))*(rhouR(1)-rhouL(1)))/dx

     do i=2,nx
      p=(aR(i)*frhou(rhoL(i),rhouL(i),eL(i),pL(i))&
            -aL(i)*frhou(rhoR(i),rhouR(i),eR(i),pR(i)))/(aR(i)-aL(i))&
            +(aL(i)*aR(i)/(aR(i)-aL(i)))*(rhouR(i)-rhouL(i))
      drhou(i-1)=drhou(i-1)+p/dx
      drhou(i)=drhou(i)-p/dx
     enddo
    drhou(nx)=drhou(nx)+((aR(nx+1)*frhou(rhoL(nx+1),rhouL(nx+1),eL(nx+1),pL(nx+1))&
                       -aL(nx+1)*frhou(rhoR(nx+1),rhouR(nx+1),eR(nx+1),pR(nx+1)))/(aR(nx+1)-aL(nx+1))&
                       +(aL(nx+1)*aR(nx+1)/(aR(nx+1)-aL(nx+1)))*(rhouR(nx+1)-rhouL(nx+1)))/dx
      
    drhou=-drhou

    
    
    !
    !flux ef energy
    !
   de(1)=-((aR(1)*fe(rhoL(1),rhouL(1),eL(1),pL(1))&
             -aL(1)*fe(rhoR(1),rhouR(1),eR(1),pR(1)))/(aR(1)-aL(1))&
             +(aL(1)*aR(1)/(aR(1)-aL(1)))*(eR(1)-eL(1)))/dx

     do i=2,nx
      p=(aR(i)*fe(rhoL(i),rhouL(i),eL(i),pL(i))&
            -aL(i)*fe(rhoR(i),rhouR(i),eR(i),pR(i)))/(aR(i)-aL(i))&
            +(aL(i)*aR(i)/(aR(i)-aL(i)))*(eR(i)-eL(i))
      de(i-1)=de(i-1)+p/dx
      de(i)=de(i)-p/dx
     enddo
    de(nx)=de(nx)+((aR(nx+1)*fe(rhoL(nx+1),rhouL(nx+1),eL(nx+1),pL(nx+1))&
                       -aL(nx+1)*fe(rhoR(nx+1),rhouR(nx+1),eR(nx+1),pR(nx+1)))/(aR(nx+1)-aL(nx+1))&
                       +(aL(nx+1)*aR(nx+1)/(aR(nx+1)-aL(nx+1)))*(eR(nx+1)-eL(nx+1)))/dx
      
    de=-de

   endsubroutine FLUXESCU



   subroutine FLUXESCULD(drho,drhou,de,rhoL,rhouL,eL,rhoR,rhouR,eR) !by Donato Cecere
    real(KND),dimension(-2:)::rhoL,rhouL,eL,rhoR,rhouR,eR
    real(KND),dimension(LBOUND(rhoL,1):UBOUND(rhoL,1))::drho,drhou,de,aL,aR,pL,pR
    real(KND) cL,cR,p,sL,sR,wint,adt
    integer i
    
    do i=1,nx+1
        pL(i)=(gamma-1._KND)*(eL(i)-(1._KND/2._KND)*(rhouL(i)*rhouL(i))/rhoL(i))
        pR(i)=(gamma-1._KND)*(eR(i)-(1._KND/2._KND)*(rhouR(i)*rhouR(i))/rhoR(i))
        if ((pL(i)<0).or.(pR(i)<0)) then
          write(*,*) "ERROR! pressure < 0",pL(i),pR(i)
        endif
        cL=SQRT(gamma*pL(i)/rhoL(i))
        cR=SQRT(gamma*pR(i)/rhoR(i))
        aL(i)=MIN(rhouL(i)/rhoL(i)-cL,rhouR(i)/rhoR(i)-cR,0._KND)
        aR(i)=MAX(rhouL(i)/rhoL(i)+cL,rhouR(i)/rhoR(i)+cR,0._KND)        
    enddo

    !    
    !flux of rho
    ! 
     ! 1
     wint   = (aR(1)*rhoR(1) -&
      aL(1)*rhoL(1)-(frho(rhoR(1),rhouR(1),eR(1),pR(1)) -&
       frho(rhoL(1),rhouL(1),eL(1),pL(1))))/(aR(1)-aL(1))
     sR     = (rhoR(1) - wint)/(aR(1)-aL(1))
	 sL     = (wint - rhoL(1))/(aR(1)-aL(1))
	 adt    = LIMITER(sL,sR)     
     drho(1)=-((aR(1)*frho(rhoL(1),rhouL(1),eL(1),pL(1))&
             -aL(1)*frho(rhoR(1),rhouR(1),eR(1),pR(1)))/(aR(1)-aL(1))&
             +(aL(1)*aR(1))*((rhoR(1)-rhoL(1))/(aR(1)-aL(1)) - adt))/dx

     do i=2,nx
      wint  = (aR(i)*rhoR(i) -&
       aL(i)*rhoL(i)-(frho(rhoR(i),rhouR(i),eR(i),pR(i)) -&
        frho(rhoL(i),rhouL(i),eL(i),pL(i))))/(aR(i)-aL(i))
      sR    = (rhoR(i) - wint)/(aR(i)-aL(i))
          sL    = (wint - rhoL(i))/(aR(i)-aL(i))
          adt   = LIMITER(sL,sR)
         p  =(aR(i)*frho(rhoL(i),rhouL(i),eL(i),pL(i))&
            -aL(i)*frho(rhoR(i),rhouR(i),eR(i),pR(i)))/(aR(i)-aL(i))&
            +(aL(i)*aR(i)*(   (rhoR(i)-rhoL(i))/(aR(i)-aL(i)) - adt ))
      drho(i-1)=drho(i-1)+p/dx
      drho(i)=drho(i)-p/dx
     enddo
     
	 ! nx
     wint = (aR(nx)*rhoR(nx) -&
      aL(nx)*rhoL(nx)-(frho(rhoR(nx),rhouR(nx),eR(nx),pR(nx)) -&
       frho(rhoL(nx),rhouL(nx),eL(nx),pL(nx))))/(aR(nx)-aL(nx))
     sR   = (rhoR(nx) - wint)/(aR(nx)-aL(nx))
     sL   = (wint - rhoL(nx))/(aR(nx)-aL(nx))
     adt  = LIMITER(sL,sR)
     drho(nx)=drho(nx)+((aR(nx+1)*frho(rhoL(nx+1),rhouL(nx+1),eL(nx+1),pL(nx+1))&
      -aL(nx+1)*frho(rhoR(nx+1),rhouR(nx+1),eR(nx+1),pR(nx+1)))/(aR(nx+1)-aL(nx+1))&
      +(aL(nx+1)*aR(nx+1))*((rhoR(nx+1)-rhoL(nx+1))/(aR(nx+1)-aL(nx+1)) - adt))/dx
      
     drho=-drho

     !
     !flux of rho*u 
     !
	 ! 1
     wint = (aR(1)*rhouR(1) -&
      aL(1)*rhouL(1)-(frhou(rhoR(1),rhouR(1),eR(1),pR(1)) -&
       frhou(rhoL(1),rhouL(1),eL(1),pL(1))))/(aR(1)-aL(1))
     sR   = (rhouR(1) - wint)/(aR(1)-aL(1))
     sL   = (wint - rhouL(1))/(aR(1)-aL(1))
     adt  = LIMITER(sL,sR)
     drhou(1)=-((aR(1)*frhou(rhoL(1),rhouL(1),eL(1),pL(1))&
             -aL(1)*frhou(rhoR(1),rhouR(1),eR(1),pR(1)))/(aR(1)-aL(1))&
             +(aL(1)*aR(1))*((rhouR(1)-rhouL(1))/(aR(1)-aL(1)) - adt))/dx

     do i=2,nx
     wint = (aR(i)*rhouR(i) -&
      aL(i)*rhouL(i)-(frhou(rhoR(i),rhouR(i),eR(i),pR(i)) -&
       frhou(rhoL(i),rhouL(i),eL(i),pL(i))))/(aR(i)-aL(i))
     sR   = (rhouR(i) - wint)/(aR(i)-aL(i))
     sL   = (wint - rhouL(i))/(aR(i)-aL(i))
     adt  = LIMITER(sL,sR)
     p=(aR(i)*frhou(rhoL(i),rhouL(i),eL(i),pL(i))&
            -aL(i)*frhou(rhoR(i),rhouR(i),eR(i),pR(i)))/(aR(i)-aL(i))&
            +(aL(i)*aR(i))*((rhouR(i)-rhouL(i))/(aR(i)-aL(i)) - adt)
      drhou(i-1)=drhou(i-1)+p/dx
      drhou(i)=drhou(i)-p/dx
     enddo
	 ! nx
     wint = (aR(nx)*rhouR(nx) -&
      aL(nx)*rhouL(nx)-(frhou(rhoR(nx),rhouR(nx),eR(nx),pR(nx)) -&
       frhou(rhoL(nx),rhouL(nx),eL(nx),pL(nx))))/(aR(nx)-aL(nx))
     sR   = (rhouR(nx) - wint)/(aR(nx)-aL(nx))
     sL   = (wint - rhouL(nx))/(aR(i)-aL(nx))
     adt  = LIMITER(sL,sR)
     drhou(nx)=drhou(nx)+((aR(nx+1)*frhou(rhoL(nx+1),rhouL(nx+1),eL(nx+1),pL(nx+1))&
                       -aL(nx+1)*frhou(rhoR(nx+1),rhouR(nx+1),eR(nx+1),pR(nx+1)))/(aR(nx+1)-aL(nx+1))&
                       +(aL(nx+1)*aR(nx+1))*((rhouR(nx+1)-rhouL(nx+1))/(aR(nx+1)-aL(nx+1))))/dx
      
     drhou=-drhou

    
    
     !
     !flux ef energy
     !
	 ! 1
     wint = (aR(1)*rhouR(1) - aL(1)*rhouL(1)-(fe(rhoR(1),rhouR(1),eR(1),pR(1)) - fe(rhoL(1),rhouL(1),eL(1),pL(1))))/(aR(1)-aL(1))
     sR   = (rhouR(1) - wint)/(aR(1)-aL(1))
     sL   = (wint - rhouL(1))/(aR(1)-aL(1))
     adt  = LIMITER(sL,sR)
     de(1)=-((aR(1)*fe(rhoL(1),rhouL(1),eL(1),pL(1))&
             -aL(1)*fe(rhoR(1),rhouR(1),eR(1),pR(1)))/(aR(1)-aL(1))&
             +(aL(1)*aR(1))*((eR(1)-eL(1))/(aR(1)-aL(1)) - adt))/dx

     do i=2,nx
     wint = (aR(i)*rhouR(i) - aL(i)*rhouL(i)-(fe(rhoR(i),rhouR(i),eR(i),pR(i)) - fe(rhoL(i),rhouL(i),eL(i),pL(i))))/(aR(i)-aL(i))
     sR   = (rhouR(i) - wint)/(aR(i)-aL(i))
     sL   = (wint - rhouL(i))/(aR(i)-aL(i))
     adt  = LIMITER(sL,sR)
      p=(aR(i)*fe(rhoL(i),rhouL(i),eL(i),pL(i))&
            -aL(i)*fe(rhoR(i),rhouR(i),eR(i),pR(i)))/(aR(i)-aL(i))&
            +(aL(i)*aR(i))*((eR(i)-eL(i))/(aR(i)-aL(i)) - adt)
      de(i-1)=de(i-1)+p/dx
      de(i)=de(i)-p/dx
     enddo
	 ! nx
     wint = (aR(nx)*rhouR(nx) -&
      aL(nx)*rhouL(nx)-(fe(rhoR(nx),rhouR(nx),eR(nx),pR(nx)) -&
       fe(rhoL(nx),rhouL(nx),eL(nx),pL(nx))))/(aR(nx)-aL(nx))
     sR   = (rhouR(nx) - wint)/(aR(nx)-aL(nx))
     sL   = (wint - rhouL(nx))/(aR(nx)-aL(nx))
     adt  = LIMITER(sL,sR)
     de(nx)=de(nx)+((aR(nx+1)*fe(rhoL(nx+1),rhouL(nx+1),eL(nx+1),pL(nx+1))&
                       -aL(nx+1)*fe(rhoR(nx+1),rhouR(nx+1),eR(nx+1),pR(nx+1)))/(aR(nx+1)-aL(nx+1))&
                       +(aL(nx+1)*aR(nx+1))*(eR(nx+1)-eL(nx+1))/(aR(nx+1)-aL(nx+1)) - adt)/dx
      
    de=-de

   endsubroutine FLUXESCULD
   
    subroutine FLUXESMUSTA(drho,drhou,de,rhoL,rhouL,eL,rhoR,rhouR,eR) !MUSTA
    real(KND),dimension(-2:)::rhoL,rhouL,eL,rhoR,rhouR,eR
    real(KND),dimension(LBOUND(rhoL,1):UBOUND(rhoL,1))::drho,drhou,de,&
                                                      flrho,flrhou,fle
    real(KND) p,A,pL,pR,pM
    real(KND),dimension(1:3):: qL,qR,qM,fL,fR,fM
    integer i,l

    A=dx/dt
    do i=1,nx+1
stg: do l=1,MUSTAstages
      if (l==1) then
       qL=(/ rhoL(i),rhouL(i),eL(i) /)
       qR=(/ rhoR(i),rhouR(i),eR(i) /)
       pL=(gamma-1._KND)*(eL(i)-0.5_KND*(rhouL(i)*rhouL(i))/rhoL(i))
       pR=(gamma-1._KND)*(eR(i)-0.5_KND*(rhouR(i)*rhouR(i))/rhoR(i))
      endif
      !predictor step
      fl(1)=frho(qL(1),qL(2),qL(3),pL)
      fl(2)=frhou(qL(1),qL(2),qL(3),pL)
      fl(3)=fe(qL(1),qL(2),qL(3),pL)
      fr(1)=frho(qR(1),qR(2),qR(3),pR)
      fr(2)=frhou(qR(1),qR(2),qR(3),pR)
      fr(3)=fe(qR(1),qR(2),qR(3),pR)

      qM=0.5_KND*(qL+qR)-(0.5_KND/A)*(fr-fl)
      pM=(gamma-1._KND)*(qM(3)-0.5_KND*(qM(2)*qM(2))/qM(1))
      
      fM(1)=frho(qM(1),qM(2),qM(3),pM)
      fM(2)=frhou(qM(1),qM(2),qM(3),pM)
      fM(3)=fe(qM(1),qM(2),qM(3),pM)
      fM=0.25_KND*(fL+2*fM+fR-A*(qR-qL))
      if (l==MUSTAstages) then
       flrho(i)=fM(1)
       flrhou(i)=fM(2)
       fle(i)=fM(3)
       EXIT stg
      else
       !corrector step - opens Riemann fan
       qL=qL-(fM-fL)/A
       qR=qR-(fR-fM)/A
       pL=(gamma-1._KND)*(qL(3)-0.5_KND*(qL(2)*qL(2))/qL(1))
       pR=(gamma-1._KND)*(qR(3)-0.5_KND*(qR(2)*qR(2))/qR(1))
      endif
     enddo stg
    enddo
    
    !update differences
    drho(1)=flrho(1)/dx
    drhou(1)=flrhou(1)/dx
    de(1)=fle(1)/dx
    do i=2,nx
     drho(i-1)=drho(i-1)-flrho(i)/dx
     drhou(i-1)=drhou(i-1)-flrhou(i)/dx
     de(i-1)=de(i-1)-fle(i)/dx
     drho(i)=flrho(i)/dx
     drhou(i)=flrhou(i)/dx
     de(i)=fle(i)/dx
    enddo
    drho(nx)=drho(nx)-flrho(nx+1)/dx
    drhou(nx)=drhou(nx)-flrhou(nx+1)/dx
    de(nx)=de(nx)-fle(nx+1)/dx
   endsubroutine FLUXESMUSTA
 

   subroutine FDLF(rho,rhou,e)
    real(KND),dimension(-2:)::rho,rhou,e
    real(KND),dimension(LBOUND(rho,1):UBOUND(rho,1))::p,rho2,rhou2,e2
    real(KND) A
    integer i
    call BOUND(rho,rhou,e)
    do i=0,nx+1
        p(i)=(gamma-1._KND)*(e(i)-(1._KND/2._KND)*(rhou(i)*rhou(i))/rho(i))
        if ((p(i)<0).or.(p(i)<0)) write(*,*) "ERROR! pressure < 0",p(i)
    enddo

   A=(dt/dx)
   do i=1,nx
    rho2(i)=(rho(i-1)+rho(i+1))/2._KND&
            -A*(frho(rho(i+1),rhou(i+1),e(i+1),p(i+1))-frho(rho(i-1),rhou(i-1),e(i-1),p(i-1)))/2._KND
   rhou2(i)=(rhou(i-1)+rhou(i+1))/2._KND&
            -A*(frhou(rho(i+1),rhou(i+1),e(i+1),p(i+1))-frhou(rho(i-1),rhou(i-1),e(i-1),p(i-1)))/2._KND
   e2(i)=(e(i-1)+e(i+1))/2._KND&
            -A*(fe(rho(i+1),rhou(i+1),e(i+1),p(i+1))-fe(rho(i-1),rhou(i-1),e(i-1),p(i-1)))/2._KND
   enddo
   rho=rho2
   rhou=rhou2
   e=e2
   endsubroutine FDLF

   subroutine FDLW(rho,rhou,e)
    real(KND),dimension(-2:)::rho,rhou,e
    real(KND),dimension(LBOUND(rho,1):UBOUND(rho,1))::p,rho2,rhou2,e2
    real(KND) A
    integer i
    real(KND),PARAMETER:: eps=0.001
    call BOUND(rho,rhou,e)
    do i=0,nx+1
        p(i)=(gamma-1._KND)*(e(i)-(1._KND/2._KND)*(rhou(i)*rhou(i))/rho(i))
        if (p(i)<0) p=0
    enddo

   A=dt/dx
   do i=1,nx
    rho2(i)=(rho(i)+rho(i+1))/2._KND&
            -A*(frho(rho(i+1),rhou(i+1),e(i+1),p(i+1))-frho(rho(i),rhou(i),e(i),p(i)))/2._KND
    if (rho2(i)<eps) rho2(i)=eps
    rhou2(i)=(rhou(i)+rhou(i+1))/2._KND&
            -A*(frhou(rho(i+1),rhou(i+1),e(i+1),p(i+1))-frhou(rho(i),rhou(i),e(i),p(i)))/2._KND
    e2(i)=(e(i)+e(i+1))/2._KND&
            -A*(fe(rho(i+1),rhou(i+1),e(i+1),p(i+1))-fe(rho(i),rhou(i),e(i),p(i)))/2._KND
    if (e2(i)<eps) e2(i)=eps        
   enddo
   
  
   call BOUND(rho2,rhou2,e2)
    do i=0,nx+1
        p(i)=(gamma-1._KND)*(e2(i)-(1._KND/2._KND)*(rhou2(i)*rhou2(i))/rho2(i))
        if (p(i)<0) p=0
    enddo

   do i=1,nx
    rho(i)=rho(i)-A*(frho(rho2(i),rhou2(i),e2(i),p(i))-frho(rho2(i-1),rhou2(i-1),e2(i-1),p(i-1)))
    if (rho(i)<eps) rho(i)=eps
    rhou(i)=rhou(i)-A*(frhou(rho2(i),rhou2(i),e2(i),p(i))-frhou(rho2(i-1),rhou2(i-1),e2(i-1),p(i-1)))
    e(i)=e(i)-A*(fe(rho2(i),rhou2(i),e2(i),p(i))-fe(rho2(i-1),rhou2(i-1),e2(i-1),p(i-1)))
    if (e(i)<eps) e(i)=eps
   enddo
   
   endsubroutine FDLW

   real(KND) function frho(rho,rhou,e,p)
    real(KND) rho,rhou,e,p
    frho=rhou
   endfunction frho
   
   real(KND) function frhou(rho,rhou,e,p)
    real(KND) rho,rhou,e,p
    frhou=rhou*rhou/rho+p
   endfunction frhou

   real(KND) function fe(rho,rhou,e,p)
    real(KND) rho,rhou,e,p
    fe=(rhou/rho)*(e+p)
    endfunction fe


  real(KND) function LIMITER(a,b)
  real(KND) a,b
  if (limitertype==1) then
     LIMITER=MINMOD(a,b)
   elseif (limitertype==2) then
     LIMITER=MINMOD(limparam*a,(a+b)/2._KND,limparam*b)
   elseif (limitertype==3) then
     LIMITER=SUPERBEE(a,b)
   endif      
   
  endfunction LIMITER
    
  real(KND) function MINMOD2(a,b)
  real(KND) a,b
      MINMOD2=(SIGN(1._KND,a)+SIGN(1._KND,b))*MIN(ABS(a),ABS(b))/2._KND
  endfunction MINMOD2

  real(KND) function MINMOD3(a,b,c)
  real(KND) a,b,c
    if ((a>0._KND).and.(b>0._KND).and.(c>0._KND)) then  
        MINMOD3=MIN(a,b,c)
    elseif ((a<0._KND).and.(b<0._KND).and.(c<0._KND)) then
        MINMOD3=MAX(a,b,c)
    else
        MINMOD3=0._KND
    endif
  endfunction MINMOD3

  real(KND) function SUPERBEE(a,b)
  real(KND) a,b,ab,b2

  ab=a*b
  b2=b*b
  if ((a*b)<=0._KND) then
   SUPERBEE=0._KND
  else if (limparam*ab<b2) then
   SUPERBEE=a*limparam
  else if (ab<b2) then
   SUPERBEE=b
  else if (ab<limparam*b2) then
   SUPERBEE=a
  else
   SUPERBEE=b*limparam
  endif
  endfunction SUPERBEE
    
  subroutine BOUND(rho,rhou,e)
   real(KND),dimension(-2:)::rho,rhou,e  

   if (extype==1) then
    rho(-2:0)=rho(1)
    rho(nx+1:nx+3)=rho(nx)

    rhou(-2:0)=rhou(1)
    rhou(nx+1:nx+3)=rhou(nx)


    e(-2:0)=e(1)
    e(nx+1:nx+3)=e(nx)
   else
    rho(-2:0)=rho(nx-2:nx)
    rho(nx+1:nx+3)=rho(1:3)

    rhou(-2:0)=rhou(nx-2:nx)
    rhou(nx+1:nx+3)=rhou(1:3)

    e(-2:0)=e(nx-2:nx)
    e(nx+1:nx+3)=e(1:3)
   endif
  endsubroutine BOUND

 
integer function lowerIndex(xvalue, x, nx)
    integer nx,i
    real(KND) x(nx)
    real(KND) xvalue
     !print*,' '
    do i=2,nx-1
        if (xvalue < x(i)) exit
    end do
    lowerIndex = i-1
    return
endfunction lowerIndex

subroutine linearMesh(yp, r, n, y, ny)
!
! yp  locations of segment matching points (returns xp)
! r   relative stretching to be applied at matching points
! n   number of matching points (including end points)
! y   mesh coordinates (output)
! ny  number of mesh points
    integer n,ny,i,j
    real(KND) yp(n), r(n)
    real(KND) y(ny)
    real(KND)xp(n)
    real(KND)a,b,x

    xp(1) = 0.
    do i=2,n
        a = r(i-1)
        b = (r(i) - r(i-1)) / (yp(i) - yp(i-1))
        if (b == 0) then
            xp(i) = xp(i-1) + (1./a)*(yp(i) - yp(i-1))
        else
            xp(i) = xp(i-1) + (1./b)*log((a + b*(yp(i) - yp(i-1))) / a)
        end if
    end do


    do j=1,ny
        x = (j-1.)/(ny-1.)*xp(n)
        i = lowerIndex(x, xp, n)
        a = r(i)
        b = (r(i+1) - r(i)) / (yp(i+1) - yp(i))
        if (b == 0) then
            y(j) = yp(i) + a*(x - xp(i))
        else
            y(j) = yp(i) + (a/b)*(exp(b*(x - xp(i))) - 1.)
        end if
    end do

    ! Normalize xp-coordinates
    xp = xp / xp(n)

    ! Return xp-coordinates in yp
    yp = xp

    return

 
 endsubroutine linearMesh
endmodule CUEuler