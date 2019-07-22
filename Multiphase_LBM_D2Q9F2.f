	PROGRAM LBM
	IMPLICIT NONE  

	Integer, parameter:: Xmax=101,Ymax=101
	Integer, parameter:: ph=2

	LOGICAL wall(Xmax,Ymax)	
	INTEGER MaxStep,NUMAX,NUMB,istep,start,i,j,k,l
	REAL*8 tau(ph),mm(ph),f(ph,0:8,Xmax,Ymax)
     &     ,ftemp(ph,0:8,Xmax,Ymax),ftemp_1(ph,0:8,Xmax,Ymax)
     &     ,nu(ph),residual,G12,Gw(ph),small
     &     ,Fx(ph,Xmax,Ymax),Fy(ph,Xmax,Ymax)

	MaxStep=1000000
	NUMAX =10000
	start=0
	NUMB=1
	residual=0.0
	ftemp_1=1.d0
	small=0.0
c*****************************************************
	tau(1)=1.0
	tau(2)=1.0
	
	mm(1)=1.0
	mm(2)=1.0

	G12=5.0

	Gw(1)= 0.060
	Gw(2)= -0.060
	
	do i=1,ph
	nu(i)=(tau(i)-0.5)/3.0
	enddo
c*****************************************************
	CALL Wall_coordinate(wall,Xmax,Ymax)
	CALL Init_density(wall,mm,ph,Xmax,Ymax,f,start)
c*****************************************************
	ITERATIONS: DO istep = 1, MaxStep
	IF (NUMB.LE.NUMAX) THEN
	NUMB=NUMB+1
	END IF

	CALL Stream(ph,Xmax,Ymax,f)
	CALL Velocity_BC(ph,Xmax,Ymax,mm,wall,f)
	CALL Interphase_Interactions(wall,ph,Xmax,Ymax,f
     &     ,Fx,Fy,G12,Gw)
	CALL Velocity_Collision(tau,mm,ph,Xmax,Ymax,f,wall
     &     ,Fx,Fy,small)
	
	IF (NUMB.GT.NUMAX) THEN
	CALL Write_results(ph,Xmax,Ymax,mm,f,istep,NUMAX
     &     ,wall,small,Fx,Fy,G12)
	NUMB=1
	END IF
	
	CALL convergence_check(f,ftemp_1,residual,ph,Xmax,Ymax)
	write(*,*)istep,residual
	residual=0.0
	ftemp_1=f
	END DO ITERATIONS

	END
c****************************************************
	SUBROUTINE Wall_coordinate(wall,Xmax,Ymax)
	IMPLICIT NONE  
	
	INTEGER  Xmax,Ymax,x,y
	LOGICAL  wall(Xmax,Ymax)

	DO  y = 1, Ymax
	DO  x = 1, Xmax
	wall(x,y) = .false.
	END DO
	END DO

c	DO  x = 1, Xmax
c	wall(x,1) = .true.
c	wall(x,Ymax) = .true.
c	END DO
c	DO  y = 1, Ymax
c	wall(1,y) = .true.
c	wall(Xmax,y) = .true.
c	END DO
		
	END
c*****************************************************
	SUBROUTINE Init_density(wall,mm,ph,Xmax,Ymax,f,start)
	IMPLICIT NONE
	
	LOGICAL  wall(Xmax,Ymax)
	INTEGER  ph,Xmax,Ymax,start,x,y,a,b,s
	REAL*8   mm(ph),f(ph,0:8,Xmax,Ymax)
     &      ,w0,w1,w2,d(ph),t1,t2
     &     ,rhonwpin,rhonwpout,rhowpin,rhowpout


	IF(start .eq. 0)THEN
	w0 =  4.d0 / 9.d0
	w1 = 1.0 /  9.d0
	w2 = 1.0 / 36.d0

c	rhonwpin=150.0
c	rhonwpout=1.e-5
c	rhowpin=1.e-5
c	rhowpout=100.0

	DO x = 1, Xmax
	DO y = 1, Ymax
	IF(.not. wall(x,y))THEN
c	t1=(x-101)**2+(y-101)**2
c	IF(t1 .le. 100)THEN
c	d(1)=2.093
c	d(2)=0.0524
c	ELSE
c	d(1)=0.0638
c	d(2)=2.011
c	ENDIF	
c

c	t1=(x-(Xmax+1)/2)**2+(y-(Ymax+1)/2)**2
	t1=(x-(Xmax+1)/2)**2+(y-(Ymax+1)/2)**2
	IF(t1 .le. 144)THEN
	d(1)=1.0
	d(2)=1.e-5
	ELSE
	d(1)=1.e-5
	d(2)=1.0 
	ENDIF	

c	IF(y .le. 21 .or. y .ge. 81)THEN
c	d(1)=1.0
c	d(2)=0.000000001
c	ELSE
c	d(1)=0.000000001
c	d(2)=1.0
c	ENDIF	

c	if(.not. wall(x,y)) then
c	if(x .eq. 1) then
c	d(1)=rhowpin
c	d(2)=rhonwpin
c	endif
c	if(x .eq. 2) then
c	d(1)=rhowpout
c	d(2)=rhonwpin
c	endif
c	if(x .eq. Xmax-1) then
c	d(1)=rhowpout
c	d(2)=rhonwpout
c	endif
c	if(x .eq. Xmax) then
c	d(1)=rhowpout
c	d(2)=rhonwpout
c	endif
c	if((x .gt. 2) .and. (x .lt. Xmax-1)) then
c	d(1)=rhowpout
c	d(2)=rhonwpout
c	endif
	else
	d(1)=0.0d0
	d(2)=0.0d0
	endif

	DO s = 1, ph
	f(s,0,x,y) = d(s)*w0/mm(s)
	f(s,1,x,y) = d(s)*w1/mm(s)
	f(s,2,x,y) = d(s)*w1/mm(s)
	f(s,3,x,y) = d(s)*w1/mm(s)
	f(s,4,x,y) = d(s)*w1/mm(s)
	f(s,5,x,y) = d(s)*w2/mm(s)
	f(s,6,x,y) = d(s)*w2/mm(s)
	f(s,7,x,y) = d(s)*w2/mm(s)
	f(s,8,x,y) = d(s)*w2/mm(s)
	END DO
	END DO
	END DO

	ELSE
	open(100,file='f1_values.dat')
	do y=1,Ymax
	do x=1,Xmax
	READ(100,*) a,b,f(1,0,x,y),f(1,1,x,y),f(1,2,x,y),f(1,3,x,y)
     &     ,f(1,4,x,y),f(1,5,x,y),f(1,6,x,y),f(1,7,x,y),f(1,8,x,y)
	end do
	end do
	close(100) 
	
	if(ph .eq. 2)then
	open(101,file='f2_values.dat')
	do y=1,Ymax
	do x=1,Xmax
	READ(101,*) a,b,f(2,0,x,y),f(2,1,x,y),f(2,2,x,y),f(2,3,x,y)
     &     ,f(2,4,x,y),f(2,5,x,y),f(2,6,x,y),f(2,7,x,y),f(2,8,x,y)
	end do
	end do
	close(101) 
	endif

	END IF

	END
c***********************************************************
	SUBROUTINE Stream(ph,Xmax,Ymax,f)
	IMPLICIT NONE

	INTEGER ph,Xmax,Ymax,x,y,xr,xl,yt,yb,s,xrr,xll,ytt,ybb
	REAL*8 f(ph,0:8,Xmax,Ymax),ftemp(ph,0:8,Xmax,Ymax)

	DO x = 1, Xmax
	DO y = 1, Ymax
	
	yt = MOD(y,Ymax) + 1
	xr = MOD(x,Xmax) + 1
	yb = Ymax - MOD(Ymax + 1 - y, Ymax)
	xl = Xmax - MOD(Xmax + 1 - x, Xmax)
	
	DO s= 1, ph
	ftemp(s,0,x ,y ) = f(s,0,x,y)
	ftemp(s,1,xr,y ) = f(s,1,x,y)
	ftemp(s,2,x ,yt) = f(s,2,x,y)
	ftemp(s,3,xl,y ) = f(s,3,x,y)
	ftemp(s,4,x ,yb) = f(s,4,x,y)
	ftemp(s,5,xr,yt) = f(s,5,x,y)
	ftemp(s,6,xl,yt) = f(s,6,x,y)
	ftemp(s,7,xl,yb) = f(s,7,x,y)
	ftemp(s,8,xr,yb) = f(s,8,x,y)
	ENDDO
	ENDDO
	ENDDO

	f=ftemp
	
	END
c**************************************************************
	SUBROUTINE Velocity_BC(ph,Xmax,Ymax,mm,wall,f)
	IMPLICIT NONE
	
	INTEGER ph,Xmax,Ymax,x,y,i,s
	LOGICAL wall(Xmax,Ymax)
	REAL*8 f(ph,0:8,Xmax,Ymax),C,mm(ph),d(ph)
     &     ,density1(ph),density2(ph),w0,w1,w2
     &     ,ftemp(ph,0:8,Xmax,Ymax)
	
	w0 =  4.d0 / 9.d0
	w1 = 1.0 /  9.d0
	w2 = 1.0 / 36.d0

	DO x = 1, Xmax
	DO y = 1, Ymax
	IF (wall(x,y)) THEN
	ftemp(:,0,x,y) = f(:,0,x,y)
	ftemp(:,1,x,y) = f(:,3,x,y)
	ftemp(:,2,x,y) = f(:,4,x,y)
	ftemp(:,3,x,y) = f(:,1,x,y)
	ftemp(:,4,x,y) = f(:,2,x,y)
	ftemp(:,5,x,y) = f(:,7,x,y)
	ftemp(:,6,x,y) = f(:,8,x,y)
	ftemp(:,7,x,y) = f(:,5,x,y)
	ftemp(:,8,x,y) = f(:,6,x,y)

	f(:,:,x,y)=ftemp(:,:,x,y)
	END IF     
	END DO
	END DO
c
c	DO s=1,ph
c	DO y=2,Ymax-1
c	f(s,1,1,y)=f(s,1,Xmax,y)
c	f(s,5,1,y)=f(s,5,Xmax,y)
c	f(s,8,1,y)=f(s,8,Xmax,y)
c
c	f(s,3,Xmax,y)=f(s,3,1,y)
c	f(s,6,Xmax,y)=f(s,6,1,y)
c	f(s,7,Xmax,y)=f(s,7,1,y)
c	ENDDO
c	DO x=1,Xmax
c	f(s,2,x,1)=f(s,2,x,Ymax)
c	f(s,5,x,1)=f(s,5,x,Ymax)
c	f(s,6,x,1)=f(s,6,x,Ymax)
c
c	f(s,4,x,Ymax)=f(s,4,x,1)
c	f(s,7,x,Ymax)=f(s,7,x,1)
c	f(s,8,x,Ymax)=f(s,8,x,1)
c	ENDDO
cc	goto 100
c	f(s,1,1,1)=f(s,1,Xmax,1)
c	f(s,2,1,1)=f(s,2,1,Ymax)
c	f(s,5,1,1)=f(s,5,Xmax,Ymax)
cc	d(s)=f(s,0,2,2)+f(s,1,2,2)+f(s,2,2,2)+f(s,3,2,2)
cc     &+f(s,4,2,2)+f(s,5,2,2)+f(s,6,2,2)+f(s,7,2,2)+f(s,8,2,2)
cc	f(s,6,1,1)=0.5*(d(s)-f(s,0,1,1)-f(s,1,1,1)-f(s,2,1,1)
cc     &     -f(s,3,1,1)-f(s,4,1,1)-f(s,5,1,1)-f(s,7,1,1))
c	f(s,8,1,1)=f(s,6,1,1)
c
c	f(s,2,Xmax,1)=f(s,2,Xmax,Ymax)
c	f(s,3,Xmax,1)=f(s,3,1,1)
c	f(s,6,Xmax,1)=f(s,6,1,Ymax)
cc	d(s)=f(s,0,Xmax-1,2)+f(s,1,Xmax-1,2)+f(s,2,Xmax-1,2)
cc     &+f(s,3,Xmax-1,2)+f(s,4,Xmax-1,2)+f(s,5,Xmax-1,2)
cc     &+f(s,6,Xmax-1,2)+f(s,7,Xmax-1,2)+f(s,8,Xmax-1,2)
cc	f(s,5,Xmax,1)=0.5*(d(s)-f(s,0,Xmax,1)-f(s,1,Xmax,1)
cc     &-f(s,2,Xmax,1)-f(s,3,Xmax,1)-f(s,4,Xmax,1)
cc     &-f(s,6,Xmax,1)-f(s,8,Xmax,1))
c	f(s,7,Xmax,1)=f(s,5,Xmax,1)
c
c	f(s,1,1,Ymax)=f(s,1,Xmax,Ymax)
c	f(s,4,1,Ymax)=f(s,4,1,1)
c	f(s,8,1,Ymax)=f(s,8,Xmax,1)
c	d(s)=f(s,0,2,Ymax-1)+f(s,1,2,Ymax-1)+f(s,2,2,Ymax-1)
cc     &+f(s,3,2,Ymax-1)+f(s,4,2,Ymax-1)+f(s,5,2,Ymax-1)
cc     &+f(s,6,2,Ymax-1)+f(s,7,2,Ymax-1)+f(s,8,2,Ymax-1)
cc	f(s,5,1,Ymax)=0.5*(d(s)-f(s,0,1,Ymax)-f(s,1,1,Ymax)
cc     &-f(s,2,1,Ymax)-f(s,3,1,Ymax)-f(s,4,1,Ymax)
cc     &-f(s,6,1,Ymax)-f(s,8,1,Ymax))
c	f(s,7,1,Ymax)=f(s,5,1,Ymax)
cc
c	f(s,3,Xmax,Ymax)=f(s,3,1,Ymax)
c	f(s,4,Xmax,Ymax)=f(s,4,Xmax,1)
c	f(s,7,Xmax,Ymax)=f(s,7,1,1)
cc	d(s)=f(s,0,Xmax-1,Ymax-1)+f(s,1,Xmax-1,Ymax-1)
cc     &+f(s,2,Xmax-1,Ymax-1)+f(s,3,Xmax-1,Ymax-1)
cc     &+f(s,4,Xmax-1,Ymax-1)+f(s,5,Xmax-1,Ymax-1)
cc     &+f(s,6,Xmax-1,Ymax-1)+f(s,7,Xmax-1,Ymax-1)
cc     &+f(s,8,Xmax-1,Ymax-1)
cc	f(s,6,Xmax,Ymax)=0.5*(d(s)-f(s,0,Xmax,Ymax)
cc     &-f(s,1,Xmax,Ymax)-f(s,2,Xmax,Ymax)
cc     &-f(s,3,Xmax,Ymax)-f(s,4,Xmax,Ymax)
cc     &-f(s,5,Xmax,Ymax)-f(s,7,Xmax,Ymax))
c	f(s,8,Xmax,Ymax)=f(s,6,Xmax,Ymax)
c 100          continue
c	ENDDO

	END
c***********************************************************
	SUBROUTINE Interphase_Interactions(wall,ph,Xmax,Ymax,f
     &     ,Fx,Fy,G12,Gw)
	IMPLICIT NONE

	INTEGER ph,Xmax,Ymax,x,y,s,yt,xr,yb,xl,xll,xrr,ytt,ybb
	LOGICAL wall(Xmax,Ymax)
	REAL*8 f(ph,0:8,Xmax,Ymax),Fx(ph,Xmax,Ymax),Fy(ph,Xmax,Ymax)
     &     ,G12,rho(ph,Xmax,Ymax),Gw(ph),D0,w0,w1,w2

	w0 = 4.d0 /  9.d0
	w1 = 1.d0 /  9.d0
	w2 = 1.d0 / 36.d0

	Fx=0.0
	Fy=0.0

	DO x=1,Xmax
	DO y=1,Ymax

	yt = MOD(y,Ymax) + 1
	xr = MOD(x,Xmax) + 1
	yb = Ymax - MOD(Ymax + 1 - y, Ymax)
	xl = Xmax - MOD(Xmax + 1 - x, Xmax)
	
	DO s=1,ph
	IF(.not. wall(x,y))THEN
	rho(s,x,y)=f(s,0,x,y)+f(s,1,x,y)+f(s,2,x,y)+f(s,3,x,y)
     &     +f(s,4,x,y)+f(s,5,x,y)+f(s,6,x,y)+f(s,7,x,y)+f(s,8,x,y)
c	rho(s,x,y)=1.0-exp(-rho(s,x,y))
	ELSE
	rho(s,x,y)=0.0
	ENDIF
	ENDDO
	ENDDO
	ENDDO

	DO x=1,Xmax
	DO y=1,Ymax

	yt = MOD(y,Ymax) + 1
	xr = MOD(x,Xmax) + 1
	yb = Ymax - MOD(Ymax+1-y,Ymax)
	xl = Xmax - MOD(Xmax+1-x,Xmax)
        xrr = MOD(xr,Xmax) + 1
        xll = Xmax - MOD(Xmax+1-xl,Xmax)
        ytt = MOD(yt,Ymax)+1
        ybb = Ymax - MOD(Ymax+1-yb,Ymax)

	

	IF(.not. wall(x,y))THEN
c	Fx(1,x,y)=-rho(1,x,y)*G12*((rho(2,xr,y)+rho(2,xr,yt)/2.0
c     &     +rho(2,xr,yb)/2.0-(rho(2,xl,y)+rho(2,xl,yt)/2.0
c     &     +rho(2,xl,yb)/2.0)))
c	Fx(2,x,y)=-rho(2,x,y)*G12*((rho(1,xr,y)+rho(1,xr,yt)/2.0
c     &     +rho(1,xr,yb)/2.0-(rho(1,xl,y)+rho(1,xl,yt)/2.0
c     &     +rho(1,xl,yb)/2.0)))
c	Fy(1,x,y)=-rho(1,x,y)*G12*(rho(2,x,yt)+rho(2,xr,yt)/2.0
c     &     +rho(2,xl,yt)/2.0-(rho(2,x,yb)+rho(2,xr,yb)/2.0
c     &     +rho(2,xl,yb)/2.0))
c	Fy(2,x,y)=-rho(2,x,y)*G12*((rho(1,x,yt)+rho(1,xr,yt)/2.0
c     &     +rho(1,xl,yt)/2.0-(rho(1,x,yb)+rho(1,xr,yb)/2.0
c     &     +rho(1,xl,yb)/2.0)))

c	Fx(1,x,y)=-rho(1,x,y)*G12*((rho(2,xr,y)+rho(2,xr,yt)/1.414
c     &     +rho(2,xr,yb)/1.414-(rho(2,xl,y)+rho(2,xl,yt)/1.414
c     &     +rho(2,xl,yb)/1.414)))
c	Fx(2,x,y)=-rho(2,x,y)*G12*((rho(1,xr,y)+rho(1,xr,yt)/1.414
c     &     +rho(1,xr,yb)/1.414-(rho(1,xl,y)+rho(1,xl,yt)/1.414
c     &     +rho(1,xl,yb)/1.414)))
c	Fy(1,x,y)=-rho(1,x,y)*G12*(rho(2,x,yt)+rho(2,xr,yt)/1.414
c     &     +rho(2,xl,yt)/1.414-(rho(2,x,yb)+rho(2,xr,yb)/1.414
c     &     +rho(2,xl,yb)/1.414))
c	Fy(2,x,y)=-rho(2,x,y)*G12*((rho(1,x,yt)+rho(1,xr,yt)/1.414
c     &     +rho(1,xl,yt)/1.414-(rho(1,x,yb)+rho(1,xr,yb)/1.414
c     &     +rho(1,xl,yb)/1.414)))

	Fx(1,x,y)=-rho(1,x,y)*G12*(1.0/3.0)*(4.0*(rho(2,xr,y)
     & -rho(2,xl,y))/21.0 + 4.0*(rho(2,xr,yt)-rho(2,xl,yt))/45.0 
     & + 4.0*(rho(2,xr,yb)-
     & rho(2,xl,yb))/45.0 + (rho(2,xrr,y)-rho(2,xll,y))/30.0 
     & + 4.0*(rho(2,xrr,yt)-rho(2,xll,yt))/315.0 + 4.0*(rho(2,xrr,yb)
     & -rho(2,xll,yb))/315.0 + 2.0*(rho(2,xr,ytt)-rho(2,xl,ytt))
     & /315.0 + 2.0*(rho(2,xr,ybb)-rho(2,xl,ybb))/315.0 
     & + (rho(2,xrr,ytt)
     & -rho(2,xll,ytt))/2520.0 + (rho(2,xrr,ybb)-rho(2,xll,ybb))
     & /2520.0)

	Fx(2,x,y)=-rho(2,x,y)*G12*(1.0/3.0)*(4.0*(rho(1,xr,y)
     & -rho(1,xl,y))/21.0 + 4.0*(rho(1,xr,yt)-rho(1,xl,yt))/45.0 
     & + 4.0*(rho(1,xr,yb)-
     & rho(1,xl,yb))/45.0 + (rho(1,xrr,y)-rho(1,xll,y))/30.0 
     & + 4.0*(rho(1,xrr,yt)-rho(1,xll,yt))/315.0 + 4.0*(rho(1,xrr,yb)
     & -rho(1,xll,yb))/315.0 + 2.0*(rho(1,xr,ytt)-rho(1,xl,ytt))
     & /315.0
     & + 2.0*(rho(1,xr,ybb)-rho(1,xl,ybb))/315.0 + (rho(1,xrr,ytt)
     & -rho(1,xll,ytt))/2520.0 + (rho(1,xrr,ybb)-rho(1,xll,ybb))
     & /2520.0)

	Fy(1,x,y)=-rho(1,x,y)*G12*(1.0/3.0)*(4.0*(rho(2,x,yt)
     & -rho(2,x,yb))/21.0+
     & 4.0*(rho(2,xr,yt)-rho(2,xr,yb))/45.0 + 4.0*(rho(2,xl,yt)-
     & rho(2,xl,yb))/45.0 + (rho(2,x,ytt)-rho(2,x,ybb))/30.0 + 
     & 4.0*(rho(2,xr,ytt)-rho(2,xr,ybb))/315.0 + 4.0*(rho(2,xl,ytt)-
     & rho(2,xl,ybb))/315.0 + 2.0*(rho(2,xrr,yt)-rho(2,xrr,yb))/315.0
     & + 2.0*(rho(2,xll,yt)-rho(2,xll,yb))/315.0 + (rho(2,xll,ytt)-
     & rho(2,xll,ybb))/2520.0 + (rho(2,xrr,ytt)-rho(2,xrr,ybb))/
     & 2520.0)

	Fy(2,x,y)=-rho(2,x,y)*G12*(1.0/3.0)*(4.0*(rho(1,x,yt)
     & -rho(1,x,yb))/21.0+
     & 4.0*(rho(1,xr,yt)-rho(1,xr,yb))/45.0 + 4.0*(rho(1,xl,yt)-
     & rho(1,xl,yb))/45.0 + (rho(1,x,ytt)-rho(1,x,ybb))/30.0 + 
     & 4.0*(rho(1,xr,ytt)-rho(1,xr,ybb))/315.0 + 4.0*(rho(1,xl,ytt)-
     & rho(1,xl,ybb))/315.0 + 2.0*(rho(1,xrr,yt)-rho(1,xrr,yb))/315.0
     & + 2.0*(rho(1,xll,yt)-rho(1,xll,yb))/315.0 + (rho(1,xll,ytt)-
     & rho(1,xll,ybb))/2520.0 + (rho(1,xrr,ytt)-rho(1,xrr,ybb))/
     & 2520.0)

c	IF(wall(x+1,y) .or. wall(x-1,y))THEN
c	Fx(1,x,y)=0.0
c	Fx(2,x,y)=0.0
c	ENDIF
c	IF(wall(x,y+1) .or. wall(x,y-1))THEN
c	Fy(1,x,y)=0.0
c	Fy(2,x,y)=0.0
c	ENDIF

	IF(wall(xr,y))THEN
	Fx(:,x,y)=Fx(:,x,y)-rho(:,x,y)*Gw(:)
	ENDIF
	IF(wall(xl,y))THEN
	Fx(:,x,y)=Fx(:,x,y)+rho(:,x,y)*Gw(:)
	ENDIF
	IF(wall(x,yt))THEN
	Fy(:,x,y)=Fy(:,x,y)-rho(:,x,y)*Gw(:)
	ENDIF
	IF(wall(x,yb))THEN
	Fy(:,x,y)=Fy(:,x,y)+rho(:,x,y)*Gw(:)
	ENDIF
	ENDIF
	ENDDO
	ENDDO

c	open(101,file='temp.dat')
c	do y=1,Ymax
c	do x=1,Xmax
c	WRITE(101,*) x,y,Fx(1,x,y),Fy(1,x,y),Fx(2,x,y),Fy(2,x,y)
c	end do
c	end do
c	close(101) 

	END
c***********************************************************
	SUBROUTINE Velocity_Collision(tau,mm,ph,Xmax,Ymax,f,wall
     &     ,Fx,Fy,small)	
	IMPLICIT NONE

	INTEGER  ph,Xmax,Ymax,x,y,i,j,k,l,s
	LOGICAL  wall(Xmax,Ymax)
	REAL*8   tau(ph),mm(ph),f(ph,0:8,Xmax,Ymax)
     &     ,C1,C2,C3,w0,w1,w2,cs_sq,small
     &     ,Dxy(ph),U(ph),V(ph),Ut,Vt,uv(0:8),fq(ph,0:8)
     &     ,Fx(ph,Xmax,Ymax),Fy(ph,Xmax,Ymax)
     &     ,rhonwpin,rhonwpout,rhowpin,rhowpout


	rhonwpin=150.0
	rhonwpout=1.e-5
	rhowpin=1.e-5
	rhowpout=100.0

c--------SRT Collision--------------
	w0 = 4.d0 /  9.d0
	w1 = 1.d0 /  9.d0
	w2 = 1.d0 / 36.d0
	
	C1 = 1.d0 / 3.d0
	C2 = 2.d0 / 9.d0
	C3 = 2.d0 / 3.d0

	DO x = 1, Xmax
	DO y = 1, Ymax
	IF (.not. wall(x,y)) THEN
	DO s= 1,ph
	Dxy(s) =f(s,0,x,y)+f(s,1,x,y)+f(s,2,x,y)
     &     +f(s,3,x,y)+f(s,4,x,y)+f(s,5,x,y)
     &     +f(s,6,x,y)+f(s,7,x,y)+f(s,8,x,y)
	ENDDO

	Ut = (((f(1,1,x,y) + f(1,5,x,y) + f(1,8,x,y)
     &-(f(1,3,x,y) + f(1,6,x,y) + f(1,7,x,y)))*mm(1)/tau(1))
     &+((f(2,1,x,y) + f(2,5,x,y) + f(2,8,x,y)
     &-(f(2,3,x,y) + f(2,6,x,y) + f(2,7,x,y)))*mm(2)/tau(2)))
     &/((mm(1)*Dxy(1)/tau(1))+(mm(2)*Dxy(2)/tau(2)))
	Vt = (((f(1,2,x,y) + f(1,5,x,y) + f(1,6,x,y)
     &-(f(1,4,x,y) + f(1,7,x,y) + f(1,8,x,y)))*mm(1)/tau(1))
     &+((f(2,2,x,y) + f(2,5,x,y) + f(2,6,x,y)
     &-(f(2,4,x,y) + f(2,7,x,y) + f(2,8,x,y)))*mm(2)/tau(2)))
     &/((mm(1)*Dxy(1)/tau(1))+(mm(2)*Dxy(2)/tau(2)))

	DO s= 1,ph
	U(s)=Ut+tau(s)*(Fx(s,x,y))/(Dxy(s)+small)
	V(s)=Vt+tau(s)*Fy(s,x,y)/(Dxy(s)+small)
c	ELSE
c	Dxy(s)=0.0
c	U(s)=0.0
c	V(s)=0.0
c	ENDIF
c	IF(x .eq. 1)THEN
c	Dxy(1)=rhowpin
c	Dxy(2)=rhonwpin
c	U(s)=0.0
c	V(s)=0.0
c	ELSEIF(x.eq.Xmax)THEN
c	Dxy(1)=rhowpout
c	Dxy(2)=rhonwpout
c	U(s)=0.0
c	V(s)=0.0
c	ENDIF	
	
	uv(0) =   U(s) * U(s) + V(s) * V(s)
	uv(1) =   U(s)
	uv(2) =   V(s)
	uv(3) = - U(s)
	uv(4) = - V(s)
	uv(5) =   U(s) + V(s)
	uv(6) = - U(s) + V(s)
	uv(7) = - U(s) - V(s)
	uv(8) =   U(s) - V(s)

	fq(s,0)=w0*Dxy(s)*(1.d0-uv(0)/C3)
	fq(s,1)=w1*Dxy(s)*(1.d0+uv(1)/C1+uv(1)**2.d0/C2-uv(0)/C3)
	fq(s,2)=w1*Dxy(s)*(1.d0+uv(2)/C1+uv(2)**2.d0/C2-uv(0)/C3)
	fq(s,3)=w1*Dxy(s)*(1.d0+uv(3)/C1+uv(3)**2.d0/C2-uv(0)/C3)
	fq(s,4)=w1*Dxy(s)*(1.d0+uv(4)/C1+uv(4)**2.d0/C2-uv(0)/C3)
	fq(s,5)=w2*Dxy(s)*(1.d0+uv(5)/C1+uv(5)**2.d0/C2-uv(0)/C3)
	fq(s,6)=w2*Dxy(s)*(1.d0+uv(6)/C1+uv(6)**2.d0/C2-uv(0)/C3)
	fq(s,7)=w2*Dxy(s)*(1.d0+uv(7)/C1+uv(7)**2.d0/C2-uv(0)/C3)
	fq(s,8)=w2*Dxy(s)*(1.d0+uv(8)/C1+uv(8)**2.d0/C2-uv(0)/C3)
	
	DO i = 0, 8
	f(s,i,x,y)=(1.0-(1.0/tau(s)))*f(s,i,x,y)+(1.0/tau(s))*fq(s,i) 
	END DO

c	if((x .eq. 1) .or. (x .eq. Xmax)) then
c	if(wall(x,y)) then
c	do i=0,8
c	f(s,i,x,y)=fq(s,i)
c	enddo
c	endif
c	endif

	END DO
	END IF
	END DO
	END DO
	END
c*************************************************************
	SUBROUTINE Write_results(ph,Xmax,Ymax,mm,f,istep,NUMAX
     &     ,wall,small,Fx,Fy,G12)
	IMPLICIT NONE

	INTEGER  ph,Xmax,Ymax,istep,NUMAX,x,y,i,s,ttt
     &     ,a,b
	LOGICAL  wall(Xmax,Ymax)
	REAL*8  f(ph,0:8,Xmax,Ymax),mm(ph),small
     &     ,U(ph,Xmax,Ymax),V(ph,Xmax,Ymax)
     &     ,Dxy(ph,Xmax,Ymax),cs,Vmag(ph,Xmax,Ymax)
     &     ,Density(Xmax,Ymax),Utotal(Xmax,Ymax)
     &     ,Vtotal(Xmax,Ymax),Pressure(Xmax,Ymax)
     &     ,Davg,Fx(ph,Xmax,Ymax),Fy(ph,Xmax,Ymax)
     &     ,G12,R,tan_theta
 	CHARACTER *6 rr

	PRINT*,'istep= ', istep
	cs = 1.d0 / 3.d0   

	WRITE (rr,200) 100000+istep/NUMAX
200     FORMAT (I6) 

	DO x=1,Xmax
	DO y=1,Ymax
	DO s= 1,ph
	Dxy(s,x,y) =mm(s)*(f(s,0,x,y)+f(s,1,x,y)+f(s,2,x,y)
     &     +f(s,3,x,y)+f(s,4,x,y)+f(s,5,x,y)
     &     +f(s,6,x,y)+f(s,7,x,y)+f(s,8,x,y))
	U(s,x,y) =(f(s,1,x,y) + f(s,5,x,y) + f(s,8,x,y)
     &-(f(s,3,x,y) + f(s,6,x,y) + f(s,7,x,y)))*mm(s)
     &/(Dxy(s,x,y)+small)
	V(s,x,y) = (f(s,2,x,y) + f(s,5,x,y) + f(s,6,x,y)
     &-(f(s,4,x,y) + f(s,7,x,y) + f(s,8,x,y)))*mm(s)
     &/(Dxy(s,x,y)+small)
	Vmag(s,x,y)=sqrt(U(s,x,y)**2.0+V(s,x,y)**2.0)
	END DO

	Density(x,y)=Dxy(1,x,y)+Dxy(2,x,y)
	Utotal(x,y)=(Dxy(1,x,y)*U(1,x,y)+Dxy(2,x,y)
     &  *U(2,x,y)+0.5*(Fx(1,x,y)+Fx(2,x,y)))/Density(x,y)
	Vtotal(x,y)=(Dxy(1,x,y)*V(1,x,y)+Dxy(2,x,y)
     &  *V(2,x,y)+0.5*(Fy(1,x,y)+Fy(2,x,y)))/Density(x,y)
c	Pressure(x,y)=(Density(x,y)/3.0)+3.0*G12*
c     &     Dxy(1,x,y)*Dxy(2,x,y)
	Pressure(x,y)=(Density(x,y)/3.0)+G12*
     &     Dxy(1,x,y)*Dxy(2,x,y)*6.6188*
     &     (2.0-0.5*((1.0/1.0)+(1.0/1.42)))	
c	Pressure(x,y)=(Density(x,y)/3.0)+6.0*G12*
c     &     Dxy(1,x,y)*Dxy(2,x,y)

	IF(wall(x,y))THEN
	Dxy(:,x,y) =0.0
	U(:,x,y) =0.0
	V(:,x,y) = 0.0
	Vmag(:,x,y)=0.0

	Density(x,y)=0.0
	Utotal(x,y)=0.0
	Vtotal(x,y)=0.0
	Pressure(x,y)=0.0
	ENDIF

	END DO
	END DO

	Davg=0.0
	ttt=0
	DO x=1,Xmax
	DO y=1,Ymax
	IF(.not. wall(x,y))THEN
	ttt=ttt+1
	Davg=Davg+Dxy(1,x,y)
	ENDIF
	ENDDO
	ENDDO
	Davg=Davg/ttt

99            format(2I8,4E15.6)

c	open(11,file='T1_'//TRIM(rr)//'.dat')
	open(11,file='T1.dat')
	WRITE(11,*) 'VARIABLES = X, Y, D1, U1, V1, Vmag1' 
	WRITE(11,*) 'ZONE I=', Xmax, ', J=', Ymax, ', F=POINT'
	DO y=1,Ymax
	DO x=1,Xmax
	WRITE(11,99) x,y,Dxy(1,x,y),U(1,x,y),V(1,x,y),Vmag(1,x,y)
	END DO
	END DO
	CLOSE(11)


c	open(12,file='T2_'//TRIM(rr)//'.dat')
	open(12,file='T2.dat')
	WRITE(12,*) 'VARIABLES = X, Y, D2, U2, V2, Vmag2 ' 
	WRITE(12,*) 'ZONE I=', Xmax, ', J=', Ymax, ', F=POINT'
	DO y=1,Ymax
	DO x=1,Xmax
	WRITE(12,99) x,y,Dxy(2,x,y),U(2,x,y),V(2,x,y),Vmag(2,x,y)
	END DO
	END DO
	CLOSE(12)

c	open(11,file='T_'//TRIM(rr)//'.dat')
	open(11,file='T.dat')
	WRITE(11,*) 'VARIABLES = X, Y, D, U, V, P' 
	WRITE(11,*) 'ZONE I=', Xmax, ', J=', Ymax, ', F=POINT'
	DO y=1,Ymax
	DO x=1,Xmax
	WRITE(11,99) x,y,Density(x,y),Utotal(x,y)
     &     ,Vtotal(x,y),Pressure(x,y)
	END DO
	END DO
	CLOSE(11)

	open(100,file='f1_values.dat')
	do y=1,Ymax
	do x=1,Xmax
	WRITE(100,*) x,y,f(1,0,x,y),f(1,1,x,y),f(1,2,x,y),f(1,3,x,y)
     &     ,f(1,4,x,y),f(1,5,x,y),f(1,6,x,y),f(1,7,x,y),f(1,8,x,y)
	end do
	end do
	close(100) 

	
	open(101,file='f2_values.dat')
	do y=1,Ymax
	do x=1,Xmax
	WRITE(101,*) x,y,f(2,0,x,y),f(2,1,x,y),f(2,2,x,y),f(2,3,x,y)
     &     ,f(2,4,x,y),f(2,5,x,y),f(2,6,x,y),f(2,7,x,y),f(2,8,x,y)
	end do
	end do
	close(101) 

	open(11,file='D_avg.dat')
	WRITE(11,*) 'Davg(1)=',Davg
	CLOSE(11)

	a=0
	DO y=2,Ymax
	x=(Xmax+1)/2
	IF(Dxy(1,x,y) .ge. 0.9*Density(x,y))THEN
	a=a+1
	ENDIF
	ENDDO

	b=0
	DO x=1,(Xmax+1)/2
	y=2
	IF(Dxy(1,x,y) .ge. 0.9*Density(x,y))THEN
	b=b+1
	ENDIF
	ENDDO
	b=b*2

	R=(a/2.0)+(b**2)/(8.0*a)
	tan_theta=b/(2.0*(R-a))

	open(11,file='R_vs_theta.dat')
	WRITE(11,*) a,b,R, tan_theta
	CLOSE(11)
	
	END
	
c********************************************************************
	SUBROUTINE convergence_check(f,ftemp_1,residual,ph,Xmax,Ymax)
	IMPLICIT NONE  

	INTEGER i,x,y,Xmax,Ymax,ph,s
	REAL*8  f(ph,0:8,Xmax,Ymax),ftemp_1(ph,0:8,Xmax,Ymax)
     &     ,residual

	DO x=1,Xmax
	DO y=1,Ymax
	DO i=0,8
	DO s=1,ph
	residual=residual+(f(s,i,x,y)-ftemp_1(s,i,x,y))**2
	ENDDO
	ENDDO
	ENDDO
	ENDDO
	residual=residual/(ph*Xmax*Ymax*9)
	END
c*********************************************************************
