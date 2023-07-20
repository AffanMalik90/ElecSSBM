PROGRAM ES_1D_NMC_T03
IMPLICIT NONE

INTEGER :: nx, ife, ifp, iter, itlp, nlp, lp, cnt
REAL(KIND=8) :: dx, dt
REAL(KIND=8) :: rho, Frd, Cst1, De, t_minus, BvP, BvE, Dmp, Kpl, Xrf
REAL(KIND=8) :: CvgP, CvgL, Vsr, tm

REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: CnP, DvC, DfC

REAL(KIND=8) :: RMS, CnE0, CnE1
REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: CnE, CEo, RHS, RES

REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: phP, ppO

REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: phE, RHE, peO

REAL(KIND=8) :: rxn, ioC, OCV, Kfw, Kbw, dPh, rxf, trgI, Cr, Scrnt
REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: Kap

nx = 200
ife = 150
ifp = 151

dx = 1.0d-5
dt = 0.01 *10
tm = 0.0

BvP =  2.5
BvE = -1.36

rho = 0.0501		!! mole/cm^3
Frd = 96485.3365
Cst1 = 1.6021766d-19/(1.3806488d-23*300.0)
De = 2.0d0*(0.73d-6*1.5d-6)/(0.73d-6+1.5d-6)
Dmp = 1.5d-6 - 0.73d-6
Kpl = Cst1*(0.73d-6+1.5d-6)
t_minus = 1.5d-6/(0.73d-6+1.5d-6)

Cr = 0.5
Vsr = 1.0d-2 *0.1

rxn = 1.0d-9
rxf = rxn/Frd

ALLOCATE(CnP(0:nx+1), DfC(0:nx+1), DvC(1:nx) )
CnP(0:nx+1) = 0.2

ALLOCATE(CnE(0:nx+1), CEo(1:nx), RHS(1:nx), RES(1:nx) )
CnE(0:nx+1) = 1.0d-3


ALLOCATE(Kap(0:nx+1), phP(0:nx+1), phE(0:nx+1), RHE(1:nx), ppO(1:nx), peO(1:nx))
phP = BvP
phE = BvE

! trgI = 0.75*rho*Frd/(3600/Cr)
trGI = 0.75*rho*Frd/(3600/Cr)*dx*(50)

cnt = 1
DO iter = 1, 1000001

	
	!! Particle Conc - ifx+1:nx+1
	DfC(ife:nx+1) = (0.0277-0.0840*CnP(ife:nx+1)+0.1003*CnP(ife:nx+1)**2)*1.0d-8
	DvC(ifp+1:nx) = 0.5*(DfC(ifp+2:nx+1)+DfC(ifp+1:nx))*(CnP(ifp+2:nx+1)-CnP(ifp+1:nx)) - &
				    0.5*(DfC(ifp+1:nx)  +DfC(ifp:nx-1))*(CnP(ifp+1:nx)-  CnP(ifp:nx-1))
	DvC(ifp+1:nx) = DvC(ifp+1:nx)/dx**2
	
	CnP(ifp+1:nx) = CnP(ifp+1:nx) + dt*DvC(ifp+1:nx)
	
	!! BC
	CnP(nx+1) = CnP(nx)
	CnP(ifp) = CnP(ifp) + dt/dx*( 0.5*(DfC(ifp+1)+Dfc(ifp))*(CnP(ifp+1)-CnP(ifp))/dx + rxf/rho )

	Xrf = SUM(CnP(ifp:nx))/(nx-ifp+1)

! 	if ( mod(iter,100) == 1 ) Print*, iter, 'CnP', CnP(ifp), Xrf

	!! Electrolyte concentration
	RHS(1:ife) = CnE(1:ife)

	RMS = 1.0
	itlp = 1

	DO WHILE (RMS > 1.0d-7)
! 	DO itlp = 1, 10001
	
		CnE(1) = (RHS(1) + dt/dx*(CnE(2)*De/dx + rxf*t_minus))/(1.0d0+dt*De/dx**2)	
		CnE(2:ife-1) = (RHS(2:ife-1) + dt*De/dx**2*(CnE(3:ife)+CnE(1:ife-2)))/(1.0d0+2*dt*De/dx**2)
		CnE(ife) = (RHS(ife) - dt/dx*(rxf*t_minus - CnE(ife-1)*De/dx))/(1.0d0+dt*De/dx**2)

		!! residual
		RES(1) = CnE(1) - dt/dx*( De*(CnE(2)-CnE(1))/dx + rxf*t_minus ) - RHS(1)
		RES(2:ife-1) = CnE(2:ife-1) - dt*De/dx**2*(CnE(3:ife)-2*CnE(2:ife-1)+CnE(1:ife-2)) - RHS(2:ife-1)
		RES(ife) = CnE(ife) + dt/dx*(rxf*t_minus + De/dx*(CnE(ife)-CnE(ife-1))) - RHS(ife)
		
		RMS = SQRT(SUM(RES(1:ife)**2))

! 		if ( mod(itlp,20) == 1 ) Print*, 'CnE', iter, itlp, RMS
		IF ( itlp > 500000 ) EXIT
! 
		itlp = itlp + 1
	ENDDO
	
	CnE(0) = CnE(1) + rxf*t_minus*dx/De
	CnE(ife+1) = CnE(ife) - rxf*t_minus*dx/De

	!! potential
	!! conductivity
	Kap(ife:nx+1) = 0.01929 + 0.7045*TANH(2.399*CnP(ife:nx+1)) - 0.7238*TANH(2.412*CnP(ife:nx+1)) - 4.2106d-6
! 
! 	!! Maybe extrapolate CnP etc. to the boundary point instead of CnP(1)
	ioC = 10.0**(-0.2*(CnP(ifp)-0.37)-1.559-0.9376*TANH(8.961*CnP(ifp)-3.195))  
	OCV = 1.095*CnP(ifp)**2-8.324d-7*EXP(14.31*CnP(ifp))+4.692*EXP(-0.5389*CnP(ifp))

	Kfw = ioC/(Frd*0.001 )*  EXP( 0.5*Cst1*OCV)
	Kbw = ioC/(Frd*CnP(ifp))*EXP(-0.5*Cst1*OCV)
 
	!! For Electrolyte potential
	RHE(1:ife) = (Dmp/dx**2)*(CnE(2:ife+1)-2*(CnE(1:ife))+CnE(0:ife-1))

	CvgP = 1.0
	CvgL = 1.0
	nlp = 1
! 
	DO WHILE ( CvgP > 1.0d-8 .AND. CvgL > 1.0d-8 ) 
! 	DO nlp = 1, 11
		dPh = phP(ifp) - phE(ife) !might change

		rxn = Kfw*CnE(ife)*EXP(-0.5*Cst1*dPh) - &
			  Kbw*CnP(ifp)*EXP( 0.5*Cst1*dPh)
			  
		rxn = rxn * 1.0d-3		!! to the unit of Amp
		rxf = rxn/Frd			!! to the unit of flux
		
! 		! Print*, rxn, Kfw, Kbw, OCV, ioC, dPh
! 
		!! particle potential
		ppO(ifp:nx) = phP(ifp:nx)

		RMS = 1.0
		lp = 1
		DO WHILE ( RMS > 1.0d-6 )
! 		DO lp = 1,100001
		
		
			phP(ifp+1:nx) = (-0.5*(Kap(ifp+2:nx+1)+Kap(ifp+1:nx))*phP(ifp+2:nx+1)  &
							 -0.5*(Kap(ifp+1:nx)+  Kap(ifp:nx-2))*phP(ifp:nx-2)  )/ &
							(-0.5*(Kap(ifp+2:nx+1)+2*Kap(ifp+1:nx)+Kap(ifp:nx-2)))

			phP(ifp) = phP(ifp+1) + rxn*2*dx/(Kap(ifp+1)+Kap(ifp))
			phP(nx+1) = BvP
				
			!! residual
			RES(ifp) = (0.5*(Kap(ifp+1)+Kap(ifp))*(phP(ifp+1)-phP(ifp))/dx + rxn)/dx
			RES(ifp+1:nx) = ( 0.5*(Kap(ifp+2:nx+1)+Kap(ifp+1:nx))*(phP(ifp+2:nx+1) - phP(ifp+1:nx)) - &
							  0.5*(Kap(ifp+1:nx)+  Kap(ifp:nx-1))*(phP(ifp+1:nx) -   phP(ifp:nx-1)))/dx**2
! 
			RMS = SQRT(SUM(RES(ifp:nx)**2))
! 
! 			if ( mod(lp,2000) == 1 ) Print*, 'phP', iter, nlp, lp, RMS	
			IF ( lp > 300000 ) EXIT

			lp = lp + 1
		ENDDO
		CvgP = SQRT(SUM((phP(ifp:nx)-ppO(ifp:nx))**2))
! 
		peO(1:ife) = phE(1:ife) 
		RHS(1:ife) = RHE(1:ife)
! 
		RMS = 1.0
		lp = 1
		DO WHILE ( RMS > 1.0d-7 )
! 		DO lp = 1,100001
			phE(1:ife-1) = ( -0.5*(CnE(2:ife)  +CnE(1:ife-1))*phE(2:ife) &
							 -0.5*(CnE(1:ife-1)+CnE(0:ife-2))*phE(0:ife-2) + RHS(1:ife-1)*dx**2/Kpl)/ &
						   ( -0.5*(CnE(2:ife)+2*CnE(1:ife-1)+CnE(0:ife-2)))

			phE(0) = BvE						   
			phE(ife) = phE(ife-1) + (RHS(ife)*dx-rxf)*2*dx/Kpl/(CnE(ife)+CnE(ife-1))
			
			!! residual
			RES(1:ife-1) = Kpl*(0.5*(CnE(2:ife)+  CnE(1:ife-1))*(phE(2:ife)-phE(1:ife-1)) - &
								0.5*(CnE(1:ife-1)+CnE(0:ife-2))*(phE(1:ife-1)-phE(0:ife-2)))/dx**2 - RHS(1:ife-1)
			RES(ife) = (-rxf - Kpl*0.5*(CnE(ife)+CnE(ife-1))*(phE(ife)-phE(ife-1))/dx)/dx - RHS(ife)
							
			RMS = SQRT(SUM(RES(1:ife)**2))
! 
! 			if ( mod(lp,2000) == 1 ) Print*, 'phE', iter, nlp, lp, RMS	
			IF ( lp > 300000 ) EXIT
! 		
			lp = lp + 1
		ENDDO 
		CvgL = SQRT(SUM((phE(1:ife)-peO(1:ife))**2))
		
		IF (MOD(nlp,500) == 1 .AND. MOD(iter,500) == 1) THEN
			Print*, 'Converg', iter, nlp, CvgP, CvgL
		ENDIF

		nlp = nlp + 1
	ENDDO

	! Scrnt = rxn/dx
	Scrnt = rxn
	BvE = BvE + DSIGN(Vsr,trgI-sCrnt)*dt
	phE(0:ife) = phE(0:ife) + DSIGN(Vsr,trgI-sCrnt)*dt
	! phE(0) = BvE
	! phE(ife) = phE(ife-1) + (RHS(ife)*dx-rxf)*2*dx/Kpl/(CnE(ife)+CnE(ife-1))

	IF (MOD(iter,500) == 1) THEN
		Print*, 'Current', iter, trGI, Scrnt, BvE, Xrf
	ENDIF
	! BvE = BvE + DSIGN(Vsr,trgI-sCrnt)*dt

	tm = tm + dt
	! IF ( MOD(iter,500) == 1 ) THEN
	IF ( Xrf > 0.2+(cnt-1)*0.005 ) THEN
		CALL OutPut1D_A(cnt, nx, 'CONPRT', 1, CnP(1:nx))
		CALL OutPut1D_A(cnt, nx, 'CONELY', 1, CnE(1:nx))
		CALL OutPut1D_A(cnt, nx, 'PHIPRT', 1, phP(1:nx))
		CALL OutPut1D_A(cnt, nx, 'PHIELY', 1, phE(1:nx))
		CALL ReC_TmFX(cnt, 'ES1DSH', iter, tm, Xrf, sCrnt, trgI, BvE)
		cnt = cnt + 1
	ENDIF

	! CALL CALL ReC_TmFX(fr, 'ES1D', iter, tm, Xrf, sCrnt, trgI, BvE)

ENDDO

END PROGRAM ES_1D_NMC_T03

!! ============================================================================================== !!
!! SUBROUTINES
!! ============================================================================================== !!
!! Data IOs
!! ============================================================================================== !!
SUBROUTINE OutPut1D_A(fr,px,prfn,AB,Cn1)
IMPLICIT NONE

INTEGER,INTENT(IN) :: fr, px, AB
REAL(KIND=8),INTENT(IN),DIMENSION(1:px) :: Cn1
CHARACTER(LEN=6),INTENT(IN) :: prfn
CHARACTER(LEN=4) :: fs
CHARACTER(LEN=14) :: flnm

WRITE(fs,"(I4)") 1000+fr

!! case 1 ==> ascii; case 2 ==> binary
SELECT CASE (AB)
	CASE (1)
		flnm = prfn//fs//'.txt'
		OPEN(UNIT=20,FILE=flnm,FORM='FORMATTED',STATUS='REPLACE',ACTION='WRITE')
			WRITE(20,*) Cn1
		CLOSE(20)
	CASE (2)
		flnm = prfn//fs//'.dat'
		OPEN(UNIT=20,FILE=flnm,FORM='UNFORMATTED',STATUS='REPLACE',ACTION='WRITE')
! 			WRITE(20) SNGL(Cn1)
			WRITE(20) Cn1
		CLOSE(20)
END SELECT

END SUBROUTINE OutPut1D_A
!! ============================================================================================== !!

SUBROUTINE ReC_TmFX(fr, prfn, iter, tm, frx, crnt, Tg, BvE)
IMPLICIT NONE

INTEGER, INTENT(IN) :: fr, iter
REAL(KIND=8), INTENT(In) :: tm, frx, crnt, Tg, BvE
CHARACTER(LEN=6),INTENT(In) :: prfn
CHARACTER(LEN=14) :: flname

flname = prfn//'TmFx.txt'
IF ( fr == 1 ) THEN 
	OPEN(UNIT=5,FILE=flname,STATUS='REPLACE',ACTION='WRITE')
		WRITE(5,*) iter, tm, frx, crnt, Tg, BvE
	CLOSE(5)	
ELSE
	OPEN(UNIT=6,FILE=flname,STATUS='OLD',ACTION='WRITE',POSITION='APPEND')
		WRITE(6,*) iter, tm, frx, crnt, Tg, BvE
	CLOSE(6)	
ENDIF


END SUBROUTINE ReC_TmFX
!! ==================== ==================== ==================== ==================== !!
