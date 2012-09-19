! Author: Vladimir Fuka
! Solution of Euler equations for compressible fluid
! in 1D - Riemann problem
! added subroutine by Donato Cecere :FLUXESCULD, and a subroutine for not uniform grid.

program Euler1D
    use CUeuler
    implicit none
    
    real(KND),allocatable,dimension(:):: rho,rhou,e,x
    integer i,fnum,frames
    real(KND) timefram1,timefram2

 !Not uniform grid
    real(KND) yp(3), r(3)
    real(KND) xp(3)
    yp(1) = -1.0
    yp(2) =  0.6
    yp(3) =  1.0
    r(1)  =  0.5
    r(2)  =  0.3
    r(3)  =  1.0
!Not uniform grid	
    
    open(11,file="euler.conf")
    
    
    read (11,FMT='(/)')
    read (11,*) FLUXtype
    write(*,*) 'Flux type:',FLUXtype
    if ((FLUXtype<1).or.(FLUXtype>5)) then
     write (*,*) "Wrong flux type, set to 5...MUSTA"
     FLUXtype=5
    endif 

    read (11,FMT='(/)')
    read (11,*) MUSTAstages
    write (*,*) 'Number of MUSTA stages:',MUSTAstages
    if (MUSTAstages<1) then
     write (*,*) "Number of MUSTA stages must be positive, set to 4."
     MUSTAstages=4
    endif
    
    read (11,FMT='(/)')
    read (11,*) RECtype
    write (*,*) 'Reconstruction type:',RECtype
    if ((RECtype<1).or.(RECtype>2)) then
     write (*,*) "Wrong reconstruction type, set to 2...p. linear."
     FLUXtype=5
    endif 
    
    read (11,FMT='(/)')
    read (11,*) limitertype
    write (*,*) 'Limiter type:',limitertype
    if ((limitertype<1).or.(limitertype>3)) then
     write (*,*) "Wrong limiter type, set to 1...minmod."
     FLUXtype=1
    endif 
    
    read (11,FMT='(/)')
    read (11,*) limparam
    read (11,FMT='(/)')
    read (11,*) RKtype
    write (*,*) 'Runge-Kutta method:',RKtype
    if ((RKtype<1).or.(RKtype>4)) then
     write(*,*) "Wrong Runge-Kutta method type, set to 2...2. order."
     RKtype=2
    elseif ((RKtype==4).and.((FLUXtype<3).or.(FLUXtype>4))) then
     write(*,*) "Full discrete scheme can be only Lax-Friedrichs or Lax-Wendroff. Runge-Kutta type set to 2...2. order."
     RKtype=2
    endif
     
    read (11,FMT='(/)')
    read (11,*) nt
    if (nt<1) then
     write (*,*) "Number of timesteps must be positive, set to 1000."
     nt=1000
    endif 

    read (11,FMT='(/)')
    read (11,*) tend
    write (*,*) 'End time:',tend
    if (tend<=0) then
     write (*,*) "End time must be positive, set to 0.1"
     tend=0.1
    endif 
    
    read (11,FMT='(/)')
    read (11,*) xl
    read (11,FMT='(/)')
    read (11,*) xr
    if (xl>=xr) then
     write (*,*) "Right margin must be greater then the left one, set to xl+2."
     xr=xl+2
    endif 
    
    read (11,FMT='(/)')
    read (11,*) nx
    if (nt<1) then
     write (*,*) "Number of cells must be positive, set to 100."
     nx=100
    endif 
    
    read (11,FMT='(/)')
    read (11,*) CFL
    if (CFL<=0) then
     write (*,*) "CFL must be positive, set to 0.5."
     CFL=0.5
    endif 

    read (11,FMT='(/)')
    read (11,*) extype
    write (*,*) 'Problem type: ',extype
    if ((extype<1).or.(extype>2)) then
     write (*,*) "Wrong type of problem to solution, sett to 1...Riemann problem."
     extype=1
    endif 

    read (11,FMT='(/)')
    read (11,*) rhoLEFT
    if (rhoLEFT<=0) then
     write (*,*) "Left state density must be positive, set to 0.1."
     rhoLEFT=0.1
    endif 
    
    read (11,FMT='(/)')
    read (11,*) rhouLEFT
   
    read (11,FMT='(/)')
    read (11,*) eLEFT
    if (eLEFT<=0) then
     write (*,*) "Left state energy must be positive, set to 0.5."
     eLEFT=0.5
    endif 
    
    read (11,FMT='(/)')
    read (11,*) rhoRIGHT
    if (rhoRIGHT<=0) then
     write (*,*) "RIGHT state density must be positive, set to 0.1."
     rhoRIGHT=0.1
    endif 
    
    read (11,FMT='(/)')
    read (11,*) rhouRIGHT
    
    read (11,FMT='(/)')
    read (11,*) eRIGHT
    if (eRIGHT<=0) then
     write (*,*) "Right state energy must be positive, set to 0.5."
     eRIGHT=0.5
    endif 
    
    read (11,FMT='(/)')
    read (11,*) frames
    read (11,FMT='(/)')
    read (11,*) timefram1
    if ((frames>0).and.(timefram1<0)) then
     write (*,*) "Time of first frame must be nonnegative, set to 0.001."
     timefram1=0.001
    endif 
    

    read (11,FMT='(/)')
    read (11,*) timefram2
    if ((frames>0).and.(timefram2<timefram1)) then
     write (*,*) "Time of last frame must be nonnegative, set to tend-0.001."
     timefram2=tend-0.001
    elseif ((frames>0).and.(timefram2>tend)) then
     write (*,*) "Time of last frame can't be greater then tend, set to tend-0.001.",frames,timefram2,tend
     timefram2=tend-0.001
    endif 

    read (11,FMT='(/)')
    read (11,*) outputtype
    if ((frames>0).and.(outputtype==1)) then
     write (*,*) "Frames for animation will be stored as muliple *.dat files."
    elseif ((frames>0).and.(outputtype==2)) then
    elseif ((frames>0).and.((outputtype<1).or.(outputtype>2))) then
     write (*,*) "Warning! No filetype set for animation frames, defaulting to text files."
     outputtype=1
    endif
    close(11)




    dx=(xr-xl)/nx

    allocate(x(-2:nx+3))



!   UNIFORM GRID
    forall(i=-2:nx+3)
     x(i)=xl+dx*(i-1._KND/2._KND)
    endforall

!   NOT UNIFORM GRID
!    call linearMesh(yp, r, 3, x(1:nx), nx)     
!	 dx = x(nx)-x(nx-1)
!	do i=1,3
!     x(nx+i) = x(nx+i-1) + dx 
!	enddo
!     dx = x(2)-x(1)
!	do i=0,-2,-1
!     x(i) = x(i+1) - dx 
!	enddo



    allocate(rho(-2:nx+3))
    allocate(rhou(-2:nx+3))
    allocate(e(-2:nx+3))


    if (extype==1) then
     forall(i=-2:nx+3,x(i)<xl+(xr-xl)/2._KND) !forall(i=-1:nx+2,x(i)<xl+(xr-xl)/2._KND)
      rho(i)=rhoLEFT
      rhou(i)=rhouLEFT
      e(i)=eLEFT
     endforall
     forall(i=-2:nx+3,x(i)>=xl+(xr-xl)/2._KND) !forall(i=-1:nx+2,x(i)>=xl+(xr-xl)/2._KND)
      rho(i)=rhoRIGHT
      rhou(i)=rhouRIGHT
      e(i)=eRIGHT
     endforall
    elseif (extype==2) then
     forall(i=1:nx)
      rho(i)=2 - (cos(pi*x(i)+dx*2._KND)-cos(pi*x(i)-dx*2._KND))/(dx*pi)
      rhou(i)=rho(i)
      e(i)=rho(i)
     endforall
    endif

   if (outputtype==2)    call open_time_evol !Create the history file
   
    time=0
    do i=1,nt  
     write(*,*) "i=",i
     write(*,*) "time=",time
     dt=DTIME(rho,rhou,e)

     if (dt>(tend-time)) dt=(tend-time)
     write (*,*) 'Initial dt ', dt
     call TMARCH(rho,rhou,e,x)
     time=time+dt
     if (frames>0)then
       if ((time>=timefram1).and.(time<=timefram2+(timefram2-timefram1)/(frames-1))&
       .and.(time>=timefram1+fnum*(timefram2-timefram1)/(frames-1))) then
        fnum=fnum+1
       if (outputtype==1) then
        call FRAME(x,rho,rhou,e,fnum)
       elseif (outputtype==2) then
        call time_evol(x,rho,rhou,e,fnum)
       endif  
   endif
     endif
     if (time>=tend) EXIT
    enddo
    write(*,*) "endtime=",time

    call OUTPUT(x,rho,rhou,e)
    if (outputtype==1) call close_time_evol
    deallocate(x,rho,rhou,e)


 contains

  real(KND) function DTIME(rho,rhou,e)
  real(KND),dimension(-2:):: rho,rhou,e
  real(KND) p,c,S,pom
  integer i
       
     S=0
     do i=1,nx
         p=(gamma-1._KND)*(e(i)-(1._KND/2._KND)*(rhou(i)*rhou(i))/rho(i))
	 if (p<0) then
	           write (*,*) "p<0!!!",e(i),rho(i),rhou(i)/rho(i)
                   p=0
         endif
         c=SQRT(gamma*p/rho(i))
         pom=MAX(ABS(rhou(i)/rho(i)-c),ABS(rhou(i)/rho(i)+c))
         if (S<p) S=pom
     enddo
    DTIME=CFL*dx/S
  endfunction DTIME

  subroutine OUTPUT(x,rho,rhou,e)
  real(KND),dimension(-2:):: x,rho,rhou,e
  integer i  
    write (*,"(a)",advance="no") "Saving....."
    open(11,file="rho.txt")
    do i=1,nx
        write(11,*) x(i),rho(i)
    enddo
    close(11)
    open(11,file="u.txt")
    do i=1,nx
        write(11,*) x(i),rhou(i)/rho(i)
    enddo
    close(11)
    open(11,file="rhou.txt")
    do i=1,nx
        write(11,*) x(i),rhou(i)
    enddo
    close(11)
    open(11,file="e.txt")
    do i=1,nx
        write(11,*) x(i),e(i)
    enddo
    close(11)
    open(11,file="p.txt")
    do i=1,nx
        write(11,*) x(i),(gamma-1._KND)*(e(i)-(1._KND/2._KND)*(rhou(i))*(rhou(i))/rho(i))
    enddo
    close(11)
    open(11,file="m.txt")
    do i=1,nx
        write(11,*) x(i),rhou(i)/(rho(i)*SQRT(ABS(gamma*(gamma-1._KND)&
                          *(e(i)-(1._KND/2._KND)*(rhou(i))*(rhou(i))/rho(i))/rho(i))))
    enddo
    close(11)
    write(*,"(a)") "saved"
  endsubroutine OUTPUT

  subroutine FRAME(x,rho,rhou,e,n)
  real(KND),dimension(-2:):: x,rho,rhou,e
  integer i,n
  character fname*10,str*70


   write(fname(1:3),"(I3.3)") n
   fname(4:7)=".dat"
   write(*,*) "Saving frame:",fname(1:3),"   time:",time
  
   open(11,file=fname(1:7))

   do i=1,nx
       write(11,*) x(i),rho(i),rhou(i)/rho(i), &
	   (gamma-1._KND)*(e(i)-(1._KND/2._KND)*(rhou(i))*(rhou(i))/rho(i))
   enddo
   close(11)
      
   write (*,*) "Saved frame",n
  endsubroutine FRAME

     SUBROUTINE open_time_evol

      INTEGER s

      CHARACTER VAR_FLUID*298, VAR_Yi*298, VAR_CK*298, VAR_PROP*298, Tec_string*500


      VAR_FLUID = "variables = ""z [m]"", ""rho [kg/m^3]"", ""Uz [m/s]"", ""P [Pa]""" 

      Tec_string = VAR_FLUID 

      open (30,file='HISTORY.plt', POSITION='APPEND', STATUS='UNKNOWN', ACTION='WRITE')
        write(30,'(500(x,a))') Tec_string

      
      END subroutine open_time_evol


      SUBROUTINE time_evol(x,rho,rhou,e,n)
      real(KND),dimension(-2:):: x,rho,rhou,e
      integer i, n, s

      CHARACTER*15 FORCHR,nvar
      write(Nvar,45) 4
   45 format(I1)
   !43 format('(F20.12,1X,F20.12,1X,F20.12,1X,F20.12,1X)')

      write(30,*) 'ZONE T="time=',time,'"'
      do i= 1,nx
        write(30,*) x(i),rho(i),rhou(i)/rho(i), &
	   (gamma-1.D0)*(e(i)-(1.D0/2.D0)*(rhou(i))*(rhou(i))/rho(i))
      enddo

      
      END subroutine time_evol


      SUBROUTINE close_time_evol

      close(30)
      END subroutine close_time_evol





end


