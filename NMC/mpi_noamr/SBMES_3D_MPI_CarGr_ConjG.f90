!! SMOOTHED BOUNDARY METHOD ELECTROCHEMICAL SIMULATION 3D

PROGRAM SBMES_3D_MPI_CarGr_ConjG_v01
IMPLICIT NONE

INCLUDE 'mpif.h'

INTEGER :: nx, ny, nz, i, j, k, iter, itlp, nlp, cnt, itle
REAL(KIND=8) :: dx, dt, rad, tm

REAL(KIND=8) :: rho, Xrf
REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: CnP, DvC, DfC

REAL(KIND=8) :: zeta, tPP
REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: psP, AvP, AvPx
REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:,:) :: ppN
! ppN - 4th dim is average of psP in each degree of freedom

REAL(KIND=8) :: BvP, CrtP, CvgP, AvrgP
REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: Kap, phP, ppO, tmp
REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:,:) :: kpN

REAL(KIND=8) :: De, Dmp, Kpl, t_minus, infx, tPE, CrtE
REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: psE
REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:,:) :: peN
REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: CnE, CEo
REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:,:) :: ceN

REAL(KIND=8) :: BvE, CrtL, CvgL, AvrgL, Vsr
REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: phE, RHE, peO

REAL(KIND=8) :: Cst1, RMS, Frd, TrgI, Cr, sCrnt
REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: Rxn, RHS, RES, ioC, OCV, Kfw, RxC
REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: Kbw, dPh

INTEGER :: gy, gx, gz, rnb, cnb, ddR, ddC, lwy, upy, lwz, upz
REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: GLP

REAL(KIND=8) :: Lsm
REAL(KIND=8),ALLOCATABLE :: VaER(:)

CHARACTER(LEN=4) :: RkNb

!! variables for material parameters
INTEGER :: cse
REAL(KIND=8),ALLOCATABLE :: Vrb(:)
REAL(KIND=8),ALLOCATABLE :: RdTb(:,:)
CHARACTER(LEN=6) :: TmNmRt

!! for MPI functions
INTEGER :: errcode, rank, np
INTEGER,DIMENSION(MPI_STATUS_SIZE) :: nstatus

REAL(KIND=8) :: DMY, alp, bet, BTM
REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: PRJ, WMX, TOR
REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:,:) :: pkN, pcN


!! -- initiate MPI
CALL MPI_INIT(errcode)
CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,errcode)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD,np,errcode)

gy = 40
gx = 65
gz = 25

!! ============================================ !!
!! -- domain decomposition, only in Y and Z directions
rnb =  8		!! rows of ranks
cnb =  5		!! columns of ranks
CALL MyDmnCmpn_YZ(rank, gy, gz, rnb, cnb, ddR, ddC, ny, lwy, upy, nz, lwz, upz)
!! y index along R; z index along C
nx = gx
!! ============================================ !!

! print*, rank, ny, lwy, upy, nz, lwz, upz

dx = 2.0d-5		!! cm	double precision
dt = 2.0d-3		!! s

zeta = 2.0

!! optimize  convergence criteria for speed
CrtE = 1.0d-7	!! convergency criterion for CnE, concentration in electrolyte
CrtP = 1.0d-10	!! convergency criterion for phP, potential in particle
CrtL = 1.0d-10	!! convergency criterion for phE, potential in electrolyte

rho = 0.0501		!! mole/cm^3; lattice density of particle
Cst1 = 1.6021766d-19/(1.3806488d-23*300.0)		!! constant F/(RT), e/(kT)
Frd = 96485.3365		!! Faraday constant

!!!!!!!!!!!!!!!!! tune these params !!!!!!!!!!!!!!!!!

BvP =  2.5				!! Boundary condition for phP
BvE = -1.7447			!! Boundary condition for phE

Cr = 1.0 					!! charge rate (0.5 -> 2 hours to charge)
Vsr = 4.0d-2      !! voltage scanning rate

!!!!!!!!!!!!!!!!! tune these params !!!!!!!!!!!!!!!!!

!! -- LiPF6 electrolyte used. D_PF6 = 1.5d-6 and D_Li = 0.73d-6 . Newman p 284
! De = 2.0d0*(0.73d-6 * 1.5d-6)/(0.73d-6 + 1.5d-6)	!! ambipolar diffusivity
! Dmp = 1.5d-6 - 0.73d-6								!! D_m - D_p
! Kpl = Cst1 * (0.73d-6 + 1.5d-6)						!! conductivity of electrolyte
! t_minus = 1.5d-6/(0.73d-6 + 1.5d-6)					!! transference number

!! D_PF6 = 4.0d-6 and D_Li = 1.25d-6 in Bernardo p 46
De = 2.0d0*(1.25d-6 * 4.0d-6)/(1.25d-6 + 4.0d-6)	!! ambipolar diffusivity
Dmp = 4.0d-6 - 1.25d-6								!! D_m - D_p
Kpl = Cst1 * (1.25d-6 + 4.0d-6)						!! conductivity of electrolyte
t_minus = 4.0d-6/(1.25d-6 + 4.0d-6)					!! transference number


!! reading variables
ALLOCATE(Vrb(5), RdTb(101,5))
OPEN(UNIT=5,FILE='newData_101x5.txt',STATUS='OLD',ACTION='READ')
		READ(5,*) RdTB(1:101,1:5)
CLOSE(5)
!! switch cases
cse = 1
TmNmRt = 'NMC101'
Vrb(1:5) = RdTb(cse,1:5)


if ( rank == 0 ) print*, cse, Vrb
		

ALLOCATE(psP(0:ny+1,0:nx+1,0:nz+1), AvP(ny,nx,nz), ppN(ny,nx,nz,6), &
	AvPx(ny,nx,nz), psE(0:ny+1,0:nx+1,0:nz+1), peN(ny,nx,nz,6))

!! domain parameter of particle (global)
ALLOCATE(GLP(0:gy+1,0:gx+1,0:gz+1) )
! OPEN(UNIT=20,FILE='/scratch/01629/hcy/NMC_2019/2019_1127_A/Geom/PS_360x118x320.dat', &
! 	FORM='UNFORMATTED',STATUS='OLD',ACTION='READ')
! 	READ(20) GLP(1:gy,1:gx,1:gz)
! CLOSE(20)
OPEN(UNIT=20,FILE='/mnt/scratch/hcy/SBMES_2020/2020_0922_A/Data0/'//'DSTF_40x65x25.txt', &
	FORM='FORMATTED',STATUS='OLD',ACTION='READ')
	READ(20,*) GLP(1:gy,1:gx,1:gz)
CLOSE(20)
GLP(1:gy,1:gx,1:gz) = 0.5*(1.0+TANH(GLP(1:gy,1:gx,1:gz)/zeta))
!! spherical particle
! DO k = 0, gz+1
! 	DO j = 0, gx+1
! 		DO i = 0, gy+1
! 			rad = SQRT( (DBLE(i)-30.0)**2 + (DBLE(j)-82.0)**2 + (DBLE(k)-24.5)**2)
! 			GLP(i,j,k) = 0.5*(1.0-TANH((rad-20.0)/zeta))
! 		ENDDO
! 	ENDDO
! ENDDO
!! pseudo 1D
! DO k = 1, gz; DO j = 1, gx; DO i = 1, gy
! 	GLP(i,j,k) = DBLE(j) - 150.5
! ENDDO; ENDDO; ENDDO
! GLP(1:gy,1:gx,1:gz) = 0.5*(1.0+TANH(GLP(1:gy,1:gx,1:gz)/zeta))


!! local psi
psP(1:ny,1:nx,1:nz) = GLP(lwy:upy,1:gx,lwz:upz)
DEALLOCATE(GLP)

CALL BcMpiYZ(rank, rnb, cnb, ddR, ddC, ny, nx, nz, psP(0:ny+1,0:nx+1,0:nz+1) )

psP(0:ny+1,0,0:nz+1) = psP(0:ny+1,1,0:nz+1)
psP(0:ny+1,nx+1,0:nz+1) = psP(0:ny+1,nx,0:nz+1)

!! domain parameter of electrolyte
psE(0:ny+1,0:nx+1,0:nz+1) = 1.0 - psP(0:ny+1,0:nx+1,0:nz+1)

!! added 1.0d-7 to avoid numerical instability
psP(0:ny+1,0:nx+1,0:nz+1) = psP(0:ny+1,0:nx+1,0:nz+1) + 1.0d-7
Lsm = SUM(psP(1:ny,1:nx,1:nz))							!! total amount of local psi_P

ALLOCATE(VaER(0:np-1))
CALL MPI_ALLGATHER(Lsm,1,MPI_DOUBLE_PRECISION,VaER(0:np-1),1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,errcode)
tPP = SUM(VaER(0:np-1))

! print*, rank, Lsm, tPP

ppN(1:ny,1:nx,1:nz,1) = (psP(0:ny-1,1:nx,1:nz)+psP(1:ny,1:nx,1:nz))*0.5d0
ppN(1:ny,1:nx,1:nz,2) = (psP(2:ny+1,1:nx,1:nz)+psP(1:ny,1:nx,1:nz))*0.5d0
ppN(1:ny,1:nx,1:nz,3) = (psP(1:ny,0:nx-1,1:nz)+psP(1:ny,1:nx,1:nz))*0.5d0
ppN(1:ny,1:nx,1:nz,4) = (psP(1:ny,2:nx+1,1:nz)+psP(1:ny,1:nx,1:nz))*0.5d0
!! add delta r in +z, -z directions
ppN(1:ny,1:nx,1:nz,5) = (psP(1:ny,1:nx,0:nz-1)+psP(1:ny,1:nx,1:nz))*0.5d0
ppN(1:ny,1:nx,1:nz,6) = (psP(1:ny,1:nx,2:nz+1)+psP(1:ny,1:nx,1:nz))*0.5d0

!! added 1.0d-7 to avoid numerical instability
psE(0:ny+1,0:nx+1,0:nz+1) = psE(0:ny+1,0:nx+1,0:nz+1) + 1.0d-9

Lsm = SUM(psE(1:ny,1:nx,1:nz))
CALL MPI_ALLGATHER(Lsm,1,MPI_DOUBLE_PRECISION,VaER(0:np-1),1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,errcode)
tPE = SUM(VaER(0:np-1))			!! total amount of psi_E
print*, rank, Lsm, tPE

! similar to above
peN(1:ny,1:nx,1:nz,1) = (psE(0:ny-1,1:nx,1:nz)+psE(1:ny,1:nx,1:nz))*0.5d0
peN(1:ny,1:nx,1:nz,2) = (psE(2:ny+1,1:nx,1:nz)+psE(1:ny,1:nx,1:nz))*0.5d0
peN(1:ny,1:nx,1:nz,3) = (psE(1:ny,0:nx-1,1:nz)+psE(1:ny,1:nx,1:nz))*0.5d0
peN(1:ny,1:nx,1:nz,4) = (psE(1:ny,2:nx+1,1:nz)+psE(1:ny,1:nx,1:nz))*0.5d0
peN(1:ny,1:nx,1:nz,5) = (psE(1:ny,1:nx,0:nz-1)+psE(1:ny,1:nx,1:nz))*0.5d0
peN(1:ny,1:nx,1:nz,6) = (psE(1:ny,1:nx,2:nz+1)+psE(1:ny,1:nx,1:nz))*0.5d0

!! particle surface  !! absolute value of gradient psi
AvP(1:ny,1:nx,1:nz) = SQRT(((psP(2:ny+1,1:nx,1:nz)-psP(0:ny-1,1:nx,1:nz))/2)**2 + &
	  				       ((psP(1:ny,2:nx+1,1:nz)-psP(1:ny,0:nx-1,1:nz))/2)**2 + &
						   ((psP(1:ny,1:nx,2:nz+1)-psP(1:ny,1:nx,0:nz-1))/2)**2)

WHERE ( AvP(1:ny,1:nx,1:nz) < 1.0d-4 )
	AvP(1:ny,1:nx,1:nz) = 0.0d0
END WHERE

AvPx(1:ny,1:nx,1:nz) = AvP(1:ny,1:nx,1:nz)/dx

!! Li fraction, Divergence of CnP, diffusivity of Li in particle
ALLOCATE(CnP(0:ny+1,0:nx+1,0:nz+1), DvC(ny,nx,nz), DfC(0:ny+1,0:nx+1,0:nz+1))
CnP(0:ny+1,0:nx+1,0:nz+1) = 0.2	!! initial value of concentration

!! average concentration in the particle
Lsm = SUM(CnP(1:ny,1:nx,1:nz)*psP(1:ny,1:nx,1:nz))
CALL MPI_ALLGATHER(Lsm,1,MPI_DOUBLE_PRECISION,VaER(0:np-1),1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,errcode)
Xrf = SUM(VaER(0:np-1))/tPP

!! salt concentration, old value of CnE,
ALLOCATE(CnE(0:ny+1,0:nx+1,0:nz+1), CEo(ny,nx,nz), ceN(ny,nx,nz,6))
CnE(0:ny+1,0:nx+1,0:nz+1) = 1.0d-3 		!! mol/cm^3

!! reaction rate
ALLOCATE(Rxn(ny,nx,nz), RHS(ny,nx,nz), RES(ny,nx,nz), ioC(ny,nx,nz), &
	OCV(ny,nx,nz), Kfw(ny,nx,nz), Kbw(ny,nx,nz), dPh(ny,nx,nz), RxC(ny,nx,nz))
! Rxn(1:ny,1:nx,1:nz) = 1.0d-7
Rxn(1:ny,1:nx,1:nz) = 0.0d0

!! potential in particle, Kap conductivity
ALLOCATE(Kap(0:ny+1,0:nx+1,0:nz+1), phP(0:ny+1,0:nx+1,0:nz+1), kpN(ny,nx,nz,6), &
	ppO(ny,nx,nz))
phP(0:ny+1,0:nx+1,0:nz+1) = BvP

!! potential in electrolyte
ALLOCATE(RHE(ny,nx,nz), phE(0:ny+1,0:nx+1,0:nz+1), peO(ny,nx,nz))
phE(0:ny+1,0:nx+1,0:nz+1) = BvE

ALLOCATE(tmp(ny,nx,nz), PRJ(0:ny+1,0:nx+1,0:nz+1), WMX(ny,nx,nz), TOR(ny,nx,nz), pkN(ny,nx,nz,6), &
	pcN(ny,nx,nz,6) )

trgI = 0.75*tPP*rho/(3600/Cr) !! target current

tm = 0.0d0

!! -- file names for output data
WRITE(RkNb,"(I4)") 1000+rank

CALL OutPut3D_A(1, ny, nx, nz, 'PS'//RkNb, 2, psP(1:ny,1:nx,1:nz))


!! time evolution
cnt = 1
iter = 1
! DO iter = 1, 1!!0001
DO WHILE ( Xrf < 0.995 )


! !  ===============================================
! !    _____      _   _               _      
! !   / ____|    | | | |             | |     
! !  | |     __ _| |_| |__   ___   __| | ___ 
! !  | |    / _` | __| '_ \ / _ \ / _` |/ _ \
! !  | |___| (_| | |_| | | | (_) | (_| |  __/
! !   \_____\__,_|\__|_| |_|\___/ \__,_|\___|
! ! ===============================================

	!! cathode particle
	!! diffusivity as a function of concentration (variables Xi1 & Xi2)
	DfC(0:ny+1,0:nx+1,0:nz+1) = (Vrb(1)*1.0775d-2 + 1.6300d-2*(CnP(0:ny+1,0:nx+1,0:nz+1)-0.5) + &
								 Vrb(2)*1.0030d-1*(CnP(0:ny+1,0:nx+1,0:nz+1)-0.5)**2)*1.0d-8		!! cm^2/s

	!! calculate the divergency of grad (psi D grad C)
	DvC(1:ny,1:nx,1:nz) = &
					 0.5*ppN(1:ny,1:nx,1:nz,2)*(DfC(2:ny+1,1:nx,1:nz)+DfC(1:ny,1:nx,1:nz))*(CnP(2:ny+1,1:nx,1:nz)-CnP(1:ny,1:nx,1:nz)) - &
					 0.5*ppN(1:ny,1:nx,1:nz,1)*(DfC(1:ny,1:nx,1:nz)+DfC(0:ny-1,1:nx,1:nz))*(CnP(1:ny,1:nx,1:nz)-CnP(0:ny-1,1:nx,1:nz)) + &
					 0.5*ppN(1:ny,1:nx,1:nz,4)*(DfC(1:ny,2:nx+1,1:nz)+DfC(1:ny,1:nx,1:nz))*(CnP(1:ny,2:nx+1,1:nz)-CnP(1:ny,1:nx,1:nz)) - &
					 0.5*ppN(1:ny,1:nx,1:nz,3)*(DfC(1:ny,1:nx,1:nz)+DfC(1:ny,0:nx-1,1:nz))*(CnP(1:ny,1:nx,1:nz)-CnP(1:ny,0:nx-1,1:nz)) + &
					 0.5*ppN(1:ny,1:nx,1:nz,6)*(DfC(1:ny,1:nx,2:nz+1)+DfC(1:ny,1:nx,1:nz))*(CnP(1:ny,1:nx,2:nz+1)-CnP(1:ny,1:nx,1:nz)) - &
					 0.5*ppN(1:ny,1:nx,1:nz,5)*(DfC(1:ny,1:nx,1:nz)+DfC(1:ny,1:nx,0:nz-1))*(CnP(1:ny,1:nx,1:nz)-CnP(1:ny,1:nx,0:nz-1))

	DvC(1:ny,1:nx,1:nz) = DvC(1:ny,1:nx,1:nz)/dx**2

	WHERE (psP(1:ny,1:nx,1:nz) > 1.0d-5)
		CnP(1:ny,1:nx,1:nz) = CnP(1:ny,1:nx,1:nz) + dt*(DvC(1:ny,1:nx,1:nz) + &
			(Rxn(1:ny,1:nx,1:nz)/rho)*AvPx(1:ny,1:nx,1:nz))/psP(1:ny,1:nx,1:nz)
	END WHERE

	!! no flux BC on X sides
	CALL BcMpiYZ(rank, rnb, cnb, ddR, ddC, ny, nx, nz, CnP(0:ny+1,0:nx+1,0:nz+1) )
	CnP(0:ny+1,0,0:nz+1) = CnP(0:ny+1,1,0:nz+1)
	CnP(0:ny+1,nx+1,0:nz+1) = CnP(0:ny+1,nx,0:nz+1)

	!! average concentration in the particle
	Lsm = SUM(CnP(1:ny,1:nx,1:nz)*psP(1:ny,1:nx,1:nz))
	CALL MPI_ALLGATHER(Lsm,1,MPI_DOUBLE_PRECISION,VaER(0:np-1),1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,errcode)
	Xrf = SUM(VaER(0:np-1))/tPP

! 	IF ( MOD(iter,20) == 1 .AND. rank == 1 ) print*, 'Xrf = ', Xrf


 
! ! ===============================================
! !   ______ _           _             _       _       
! !  |  ____| |         | |           | |     | |      
! !  | |__  | | ___  ___| |_ _ __ ___ | |_   _| |_ ___ 
! !  |  __| | |/ _ \/ __| __| '__/ _ \| | | | | __/ _ \
! !  | |____| |  __/ (__| |_| | | (_) | | |_| | ||  __/
! !  |______|_|\___|\___|\__|_|  \___/|_|\__, |\__\___|
! !                                       __/ |        
! !                                      |___/      
! ! ===============================================  
 
	!! electrolyte concentration
	CEo(1:ny,1:nx,1:nz) = CnE(1:ny,1:nx,1:nz)
	RHS(1:ny,1:nx,1:nz) = psE(1:ny,1:nx,1:nz)*CEo(1:ny,1:nx,1:nz) - dt*(Rxn(1:ny,1:nx,1:nz))*AvPx(1:ny,1:nx,1:nz)*t_minus

	!! influx of Li
	Lsm = SUM(Rxn(1:ny,1:nx,1:nz)*AvPx(1:ny,1:nx,1:nz))*(dx**3)   !!/(ny*nz*dx**2)
	CALL MPI_ALLGATHER(Lsm,1,MPI_DOUBLE_PRECISION,VaER(0:np-1),1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,errcode)
	infx = SUM(VaER(0:np-1))/(gy*gz*dx**2)

! ! 	IF ( MOD(iter,20) == 1 .AND. rank == 1 ) print*, 'infx = ', infx

	CnE(0:ny+1,0,0:nz+1) = infx*t_minus*dx/De + CnE(0:ny+1,1,0:nz+1)

	!! conjugate gradient method
	RES(1:ny,1:nx,1:nz) = RHS(1:ny,1:nx,1:nz) - &
		(psE(1:ny,1:nx,1:nz)*CnE(1:ny,1:nx,1:nz) - dt*De/dx**2 * &
		(peN(1:ny,1:nx,1:nz,2)*(CnE(2:ny+1,1:nx,1:nz)-CnE(1:ny,1:nx,1:nz)) - &
		 peN(1:ny,1:nx,1:nz,1)*(CnE(1:ny,1:nx,1:nz)-CnE(0:ny-1,1:nx,1:nz)) + &
		 peN(1:ny,1:nx,1:nz,4)*(CnE(1:ny,2:nx+1,1:nz)-CnE(1:ny,1:nx,1:nz)) - &
		 peN(1:ny,1:nx,1:nz,3)*(CnE(1:ny,1:nx,1:nz)-CnE(1:ny,0:nx-1,1:nz)) + &
		 peN(1:ny,1:nx,1:nz,6)*(CnE(1:ny,1:nx,2:nz+1)-CnE(1:ny,1:nx,1:nz)) - &
		 peN(1:ny,1:nx,1:nz,5)*(CnE(1:ny,1:nx,1:nz)-CnE(1:ny,1:nx,0:nz-1))))

	PRJ(1:ny,1:nx,1:nz) = RES(1:ny,1:nx,1:nz)
		
	RMS = 1.0
	itlp = 1
	IF ( iter > 1 ) THEN
	
		!! internal loop
! ! 	DO itlp = 1, 1001	!!0001
		DO WHILE ( RMS > CrtE )

			CALL BcMpiYZ(rank, rnb, cnb, ddR, ddC, ny, nx, nz, PRJ(0:ny+1,0:nx+1,0:nz+1) )
			PRJ(0:ny+1,nx+1,0:nz+1) = PRJ(0:ny+1,nx,0:nz+1)
			PRJ(0:ny+1,0,0:nz+1) = PRJ(0:ny+1,1,0:nz+1)
		
			WMX(1:ny,1:nx,1:nz) = psE(1:ny,1:nx,1:nz)*PRJ(1:ny,1:nx,1:nz) - dt*De/dx**2 * &
				(peN(1:ny,1:nx,1:nz,2)*(PRJ(2:ny+1,1:nx,1:nz)-PRJ(1:ny,1:nx,1:nz)) - &
				 peN(1:ny,1:nx,1:nz,1)*(PRJ(1:ny,1:nx,1:nz)-PRJ(0:ny-1,1:nx,1:nz)) + &
				 peN(1:ny,1:nx,1:nz,4)*(PRJ(1:ny,2:nx+1,1:nz)-PRJ(1:ny,1:nx,1:nz)) - &
				 peN(1:ny,1:nx,1:nz,3)*(PRJ(1:ny,1:nx,1:nz)-PRJ(1:ny,0:nx-1,1:nz)) + &
				 peN(1:ny,1:nx,1:nz,6)*(PRJ(1:ny,1:nx,2:nz+1)-PRJ(1:ny,1:nx,1:nz)) - &
				 peN(1:ny,1:nx,1:nz,5)*(PRJ(1:ny,1:nx,1:nz)-PRJ(1:ny,1:nx,0:nz-1))) 

			RMS = SUM(RES(1:ny,1:nx,1:nz)**2)
			CALL MPI_ALLGATHER(RMS,1,MPI_DOUBLE_PRECISION,VaER(0:np-1),1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,errcode)			
			DMY = SUM(VaER(0:np-1))
	
			RMS = SUM(PRJ(1:ny,1:nx,1:nz)*WMX(1:ny,1:nx,1:nz))
			CALL MPI_ALLGATHER(RMS,1,MPI_DOUBLE_PRECISION,VaER(0:np-1),1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,errcode)			
			BTM = SUM(VaER(0:np-1))
								
			alp = DMY/BTM
		
			CnE(1:ny,1:nx,1:nz) = CnE(1:ny,1:nx,1:nz) + alp*PRJ(1:ny,1:nx,1:nz)
		
			RES(1:ny,1:nx,1:nz) = RES(1:ny,1:nx,1:nz) - alp*WMX(1:ny,1:nx,1:nz)

			RMS = SUM(RES(1:ny,1:nx,1:nz)**2)
			CALL MPI_ALLGATHER(RMS,1,MPI_DOUBLE_PRECISION,VaER(0:np-1),1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,errcode)			
			BTM = SUM(VaER(0:np-1))
		
			bet = BTM/DMY
			PRJ(1:ny,1:nx,1:nz) = RES(1:ny,1:nx,1:nz) + bet*PRJ(1:ny,1:nx,1:nz)
			
		
			Lsm = SUM(RES(1:ny,1:nx,1:nz)**2*psE(1:ny,1:nx,1:nz))
			CALL MPI_ALLGATHER(Lsm,1,MPI_DOUBLE_PRECISION,VaER(0:np-1),1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,errcode)
			RMS = SQRT(SUM(VaER(0:np-1))/tPE)			


			if ( mod(itlp,200) == 1 .AND. mod(iter,100) == 1 .AND. rank == 0 ) Print*, 'CnE', iter, itlp, RMS
			IF ( itlp > 50000 ) EXIT

			itlp = itlp + 1
		ENDDO
  	ENDIF
	!! boundary condition
	CALL BcMpiYZ(rank, rnb, cnb, ddR, ddC, ny, nx, nz, CnE(0:ny+1,0:nx+1,0:nz+1) )
	CnE(0:ny+1,0,0:nz+1) = infx*t_minus*dx/De + CnE(0:ny+1,1,0:nz+1)
	CnE(0:ny+1,nx+1,0:nz+1) = CnE(0:ny+1,nx,0:nz+1)

	!! electrolyte potential !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	RHE(1:ny,1:nx,1:nz) = Dmp/dx**2*(peN(1:ny,1:nx,1:nz,2)*(CnE(2:ny+1,1:nx,1:nz)-CnE(1:ny,1:nx,1:nz)) - &
									 peN(1:ny,1:nx,1:nz,1)*(CnE(1:ny,1:nx,1:nz)-CnE(0:ny-1,1:nx,1:nz)) + &
									 peN(1:ny,1:nx,1:nz,4)*(CnE(1:ny,2:nx+1,1:nz)-CnE(1:ny,1:nx,1:nz)) - &
									 peN(1:ny,1:nx,1:nz,3)*(CnE(1:ny,1:nx,1:nz)-CnE(1:ny,0:nx-1,1:nz)) + &
									 peN(1:ny,1:nx,1:nz,6)*(CnE(1:ny,1:nx,2:nz+1)-CnE(1:ny,1:nx,1:nz)) - &
									 peN(1:ny,1:nx,1:nz,5)*(CnE(1:ny,1:nx,1:nz)-CnE(1:ny,1:nx,0:nz-1)))

	ceN(1:ny,1:nx,1:nz,1) = (CnE(0:ny-1,1:nx,1:nz)+CnE(1:ny,1:nx,1:nz))*0.5d0
	ceN(1:ny,1:nx,1:nz,2) = (CnE(2:ny+1,1:nx,1:nz)+CnE(1:ny,1:nx,1:nz))*0.5d0
	ceN(1:ny,1:nx,1:nz,3) = (CnE(1:ny,0:nx-1,1:nz)+CnE(1:ny,1:nx,1:nz))*0.5d0
	ceN(1:ny,1:nx,1:nz,4) = (CnE(1:ny,2:nx+1,1:nz)+CnE(1:ny,1:nx,1:nz))*0.5d0
	ceN(1:ny,1:nx,1:nz,5) = (CnE(1:ny,1:nx,0:nz-1)+CnE(1:ny,1:nx,1:nz))*0.5d0
	ceN(1:ny,1:nx,1:nz,6) = (CnE(1:ny,1:nx,2:nz+1)+CnE(1:ny,1:nx,1:nz))*0.5d0
	
	pcN(1:ny,1:nx,1:nz,1:6) = peN(1:ny,1:nx,1:nz,1:6)*ceN(1:ny,1:nx,1:nz,1:6)

	!! particle conductivity (variable Xi3)
	Kap(1:ny,1:nx,1:nz) = Vrb(3)*(0.01929 - 4.2106d-6) + 0.7045*TANH(2.399*CnP(1:ny,1:nx,1:nz)) & 
													   - 0.7238*TANH(2.412*CnP(1:ny,1:nx,1:nz)) 

	CALL BcMpiYZ(rank, rnb, cnb, ddR, ddC, ny, nx, nz, Kap(0:ny+1,0:nx+1,0:nz+1) )
	Kap(0:ny+1,0,0:nz+1) = Kap(0:ny+1,1,0:nz+1)
	Kap(0:ny+1,nx+1,0:nz+1) = Kap(0:ny+1,nx,0:nz+1)

	kpN(1:ny,1:nx,1:nz,1) = (Kap(0:ny-1,1:nx,1:nz)+Kap(1:ny,1:nx,1:nz))*0.5d0
	kpN(1:ny,1:nx,1:nz,2) = (Kap(2:ny+1,1:nx,1:nz)+Kap(1:ny,1:nx,1:nz))*0.5d0
	kpN(1:ny,1:nx,1:nz,3) = (Kap(1:ny,0:nx-1,1:nz)+Kap(1:ny,1:nx,1:nz))*0.5d0
	kpN(1:ny,1:nx,1:nz,4) = (Kap(1:ny,2:nx+1,1:nz)+Kap(1:ny,1:nx,1:nz))*0.5d0
	kpN(1:ny,1:nx,1:nz,5) = (Kap(1:ny,1:nx,0:nz-1)+Kap(1:ny,1:nx,1:nz))*0.5d0
	kpN(1:ny,1:nx,1:nz,6) = (Kap(1:ny,1:nx,2:nz+1)+Kap(1:ny,1:nx,1:nz))*0.5d0
	
	pkN(1:ny,1:nx,1:nz,1:6) = ppN(1:ny,1:nx,1:nz,1:6)*kpN(1:ny,1:nx,1:nz,1:6)

	!! exchange current density (variables Xi4 & Xi5)
	WHERE ( AvP(1:ny,1:nx,1:nz) > 1.0d-4 )
		ioC(1:ny,1:nx,1:nz) = 10.0**(-0.2*(CnP(1:ny,1:nx,1:nz)-0.37) - 1.559 - &
									 Vrb(4)*0.9376*TANH(Vrb(5)*8.961*CnP(1:ny,1:nx,1:nz)-3.195))
		ioC = ioC * 1.0d-3 !! in unit of Amp

		OCV(1:ny,1:nx,1:nz) = 1.095*CnP(1:ny,1:nx,1:nz)**2-8.324d-7*EXP(14.31*CnP(1:ny,1:nx,1:nz))+4.692*EXP(-0.5389*CnP(1:ny,1:nx,1:nz))

		Kfw(1:ny,1:nx,1:nz) = ioC(1:ny,1:nx,1:nz)/(Frd*0.001              )*EXP( 0.5*Cst1*OCV(1:ny,1:nx,1:nz))
		Kbw(1:ny,1:nx,1:nz) = ioC(1:ny,1:nx,1:nz)/(Frd*CnP(1:ny,1:nx,1:nz))*EXP(-0.5*Cst1*OCV(1:ny,1:nx,1:nz))
	END WHERE

	CvgP = 1.0
	CvgL = 1.0
	nlp = 1
! ! 	! DO WHILE ( CvgP > 1.0d-11 .AND. CvgL > 1.0d-11 )
	DO WHILE ( CvgP > 1.0d-8 .AND. CvgL > 1.0d-8 )
! ! 	DO nlp = 1, 101

		WHERE ( AvP(1:ny,1:nx,1:nz) > 1.0d-4 )
			!! in the unit of Volt
			dPh(1:ny,1:nx,1:nz) = phP(1:ny,1:nx,1:nz) - phE(1:ny,1:nx,1:nz)

			!! in the unit of flux
			Rxn(1:ny,1:nx,1:nz) = Kfw(1:ny,1:nx,1:nz)*CnE(1:ny,1:nx,1:nz)*EXP(-0.5*Cst1*dPh(1:ny,1:nx,1:nz)) - &
						          Kbw(1:ny,1:nx,1:nz)*CnP(1:ny,1:nx,1:nz)*EXP( 0.5*Cst1*dPh(1:ny,1:nx,1:nz))

			!! in the unit of Amp
			RxC(1:ny,1:nx,1:nz) = Rxn(1:ny,1:nx,1:nz)*Frd
		END WHERE

		!! particle potential
		ppO(1:ny,1:nx,1:nz) = phP(1:ny,1:nx,1:nz)
		RHS(1:ny,1:nx,1:nz) = -RxC(1:ny,1:nx,1:nz)*AvPx(1:ny,1:nx,1:nz)	!! in unit of Amp

		phP(0:ny+1,nx+1,0:nz+1) = BvP
		!! conjugate gradient method
		RES(1:ny,1:nx,1:nz) = RHS(1:ny,1:nx,1:nz) - 1.0/dx**2* &
			(pkN(1:ny,1:nx,1:nz,2)*(phP(2:ny+1,1:nx,1:nz)-phP(1:ny,1:nx,1:nz)) - &
			 pkN(1:ny,1:nx,1:nz,1)*(phP(1:ny,1:nx,1:nz)-phP(0:ny-1,1:nx,1:nz)) + &
			 pkN(1:ny,1:nx,1:nz,4)*(phP(1:ny,2:nx+1,1:nz)-phP(1:ny,1:nx,1:nz)) - &
			 pkN(1:ny,1:nx,1:nz,3)*(phP(1:ny,1:nx,1:nz)-phP(1:ny,0:nx-1,1:nz)) + &
			 pkN(1:ny,1:nx,1:nz,6)*(phP(1:ny,1:nx,2:nz+1)-phP(1:ny,1:nx,1:nz)) - &
			 pkN(1:ny,1:nx,1:nz,5)*(phP(1:ny,1:nx,1:nz)-phP(1:ny,1:nx,0:nz-1)))
			 
		PRJ(1:ny,1:nx,1:nz) = RES(1:ny,1:nx,1:nz)
		
		RMS = 1.0
		itlp = 1
		DO WHILE ( RMS > CrtP )
! 		DO itlp = 1, 10001
		
			tmp(1:ny,1:nx,1:nz) = phP(1:ny,1:nx,1:nz)
			
			CALL BcMpiYZ(rank, rnb, cnb, ddR, ddC, ny, nx, nz, PRJ(0:ny+1,0:nx+1,0:nz+1) )
			PRJ(0:ny+1,nx+1,0:nz+1) = 0.0d0
			PRJ(0:ny+1,0,0:nz+1) = PRJ(0:ny+1,1,0:nz+1)
			
			WMX(1:ny,1:nx,1:nz) = 1.0/dx**2* &
				(pkN(1:ny,1:nx,1:nz,2)*(PRJ(2:ny+1,1:nx,1:nz)-PRJ(1:ny,1:nx,1:nz)) - &
				 pkN(1:ny,1:nx,1:nz,1)*(PRJ(1:ny,1:nx,1:nz)-PRJ(0:ny-1,1:nx,1:nz)) + &
				 pkN(1:ny,1:nx,1:nz,4)*(PRJ(1:ny,2:nx+1,1:nz)-PRJ(1:ny,1:nx,1:nz)) - &
				 pkN(1:ny,1:nx,1:nz,3)*(PRJ(1:ny,1:nx,1:nz)-PRJ(1:ny,0:nx-1,1:nz)) + &
				 pkN(1:ny,1:nx,1:nz,6)*(PRJ(1:ny,1:nx,2:nz+1)-PRJ(1:ny,1:nx,1:nz)) - &
				 pkN(1:ny,1:nx,1:nz,5)*(PRJ(1:ny,1:nx,1:nz)-PRJ(1:ny,1:nx,0:nz-1))) 
			
			RMS = SUM(RES(1:ny,1:nx,1:nz)**2)
			CALL MPI_ALLGATHER(RMS,1,MPI_DOUBLE_PRECISION,VaER(0:np-1),1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,errcode)			
			DMY = SUM(VaER(0:np-1))
		
			RMS = SUM(PRJ(1:ny,1:nx,1:nz)*WMX(1:ny,1:nx,1:nz))
			CALL MPI_ALLGATHER(RMS,1,MPI_DOUBLE_PRECISION,VaER(0:np-1),1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,errcode)			
			BTM = SUM(VaER(0:np-1))
									
			alp = DMY/BTM
			
			phP(1:ny,1:nx,1:nz) = phP(1:ny,1:nx,1:nz) + alp*PRJ(1:ny,1:nx,1:nz)
			
			RES(1:ny,1:nx,1:nz) = RES(1:ny,1:nx,1:nz) - alp*WMX(1:ny,1:nx,1:nz)

			RMS = SUM(RES(1:ny,1:nx,1:nz)**2)
			CALL MPI_ALLGATHER(RMS,1,MPI_DOUBLE_PRECISION,VaER(0:np-1),1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,errcode)			
			BTM = SUM(VaER(0:np-1))
			
			bet = BTM/DMY
			PRJ(1:ny,1:nx,1:nz) = RES(1:ny,1:nx,1:nz) + bet*PRJ(1:ny,1:nx,1:nz)
			
			TOR(1:ny,1:nx,1:nz) = phP(1:ny,1:nx,1:nz) - tmp(1:ny,1:nx,1:nz)

! 			Lsm = SUM(RES(1:ny,1:nx,1:nz)**2*psP(1:ny,1:nx,1:nz))
			Lsm = SUM(TOR(1:ny,1:nx,1:nz)**2*psP(1:ny,1:nx,1:nz))
			
			CALL MPI_ALLGATHER(Lsm,1,MPI_DOUBLE_PRECISION,VaER(0:np-1),1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,errcode)
			RMS = SQRT(SUM(VaER(0:np-1))/tPP)

			if ( mod(itlp,200) == 1 .AND. mod(iter,500) == 1 .AND. rank == 0 ) Print*, 'phP', iter, nlp, itlp, RMS
			
			IF ( itlp > 30000 ) EXIT

			itlp = itlp + 1

		ENDDO
		
		!! BC
		CALL BcMpiYZ(rank, rnb, cnb, ddR, ddC, ny, nx, nz, phP(0:ny+1,0:nx+1,0:nz+1) )	
		phP(0:ny+1,nx+1,0:nz+1) = BvP
		phP(0:ny+1,0,0:nz+1) = phP(0:ny+1,1,0:nz+1)			

		Lsm = SUM((phP(1:ny,1:nx,1:nz)-ppO(1:ny,1:nx,1:nz))**2*psP(1:ny,1:nx,1:nz))
		CALL MPI_ALLGATHER(Lsm,1,MPI_DOUBLE_PRECISION,VaER(0:np-1),1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,errcode)
		CvgP = SQRT(SUM(VaER(0:np-1))/tPP)
! 
		!! electrolyte potential
		peO(1:ny,1:nx,1:nz) = phE(1:ny,1:nx,1:nz)
		RHS(1:ny,1:nx,1:nz) = (Rxn(1:ny,1:nx,1:nz))*AvPx(1:ny,1:nx,1:nz) + RHE(1:ny,1:nx,1:nz) !! flux

		RES(1:ny,1:nx,1:nz) = RHS(1:ny,1:nx,1:nz) - Kpl/dx**2* &
			(pcN(1:ny,1:nx,1:nz,2)*(phE(2:ny+1,1:nx,1:nz)-phE(1:ny,1:nx,1:nz)) - &
			 pcN(1:ny,1:nx,1:nz,1)*(phE(1:ny,1:nx,1:nz)-phE(0:ny-1,1:nx,1:nz)) + &
			 pcN(1:ny,1:nx,1:nz,4)*(phE(1:ny,2:nx+1,1:nz)-phE(1:ny,1:nx,1:nz)) - &
			 pcN(1:ny,1:nx,1:nz,3)*(phE(1:ny,1:nx,1:nz)-phE(1:ny,0:nx-1,1:nz)) + &
			 pcN(1:ny,1:nx,1:nz,6)*(phE(1:ny,1:nx,2:nz+1)-phE(1:ny,1:nx,1:nz)) - &
			 pcN(1:ny,1:nx,1:nz,5)*(phE(1:ny,1:nx,1:nz)-phE(1:ny,1:nx,0:nz-1)))

		PRJ(1:ny,1:nx,1:nz) = RES(1:ny,1:nx,1:nz)
		PRJ(0:ny+1,0,0:nz+1) = 0.0d0

		RMS = 1.0
		itle = 1
		DO WHILE ( RMS > CrtL )
! 		DO itlp = 1, 10001
! 
			tmp(1:ny,1:nx,1:nz) = phE(1:ny,1:nx,1:nz)

			CALL BcMpiYZ(rank, rnb, cnb, ddR, ddC, ny, nx, nz, PRJ(0:ny+1,0:nx+1,0:nz+1) )
			PRJ(0:ny+1,nx+1,0:nz+1) = PRJ(0:ny+1,nx,0:nz+1)
			PRJ(0:ny+1,0,0:nz+1) = 0.0d0

			WMX(1:ny,1:nx,1:nz) = Kpl/dx**2* &
				(pcN(1:ny,1:nx,1:nz,2)*(PRJ(2:ny+1,1:nx,1:nz)-PRJ(1:ny,1:nx,1:nz)) - &
				 pcN(1:ny,1:nx,1:nz,1)*(PRJ(1:ny,1:nx,1:nz)-PRJ(0:ny-1,1:nx,1:nz)) + &
				 pcN(1:ny,1:nx,1:nz,4)*(PRJ(1:ny,2:nx+1,1:nz)-PRJ(1:ny,1:nx,1:nz)) - &
				 pcN(1:ny,1:nx,1:nz,3)*(PRJ(1:ny,1:nx,1:nz)-PRJ(1:ny,0:nx-1,1:nz)) + &
				 pcN(1:ny,1:nx,1:nz,6)*(PRJ(1:ny,1:nx,2:nz+1)-PRJ(1:ny,1:nx,1:nz)) - &
				 pcN(1:ny,1:nx,1:nz,5)*(PRJ(1:ny,1:nx,1:nz)-PRJ(1:ny,1:nx,0:nz-1))) 
			
			RMS = SUM(RES(1:ny,1:nx,1:nz)**2)
			CALL MPI_ALLGATHER(RMS,1,MPI_DOUBLE_PRECISION,VaER(0:np-1),1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,errcode)			
			DMY = SUM(VaER(0:np-1))
			
			RMS = SUM(PRJ(1:ny,1:nx,1:nz)*WMX(1:ny,1:nx,1:nz))
			CALL MPI_ALLGATHER(RMS,1,MPI_DOUBLE_PRECISION,VaER(0:np-1),1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,errcode)			
			BTM = SUM(VaER(0:np-1))
									
			alp = DMY/BTM

			phE(1:ny,1:nx,1:nz) = phE(1:ny,1:nx,1:nz) + alp*PRJ(1:ny,1:nx,1:nz)
			
			RES(1:ny,1:nx,1:nz) = RES(1:ny,1:nx,1:nz) - alp*WMX(1:ny,1:nx,1:nz)

			RMS = SUM(RES(1:ny,1:nx,1:nz)**2)
			CALL MPI_ALLGATHER(RMS,1,MPI_DOUBLE_PRECISION,VaER(0:np-1),1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,errcode)			
			BTM = SUM(VaER(0:np-1))
			
			bet = BTM/DMY
			PRJ(1:ny,1:nx,1:nz) = RES(1:ny,1:nx,1:nz) + bet*PRJ(1:ny,1:nx,1:nz)

			TOR(1:ny,1:nx,1:nz) = phE(1:ny,1:nx,1:nz) - tmp(1:ny,1:nx,1:nz)

! 			Lsm = SUM(RES(1:ny,1:nx,1:nz)**2*psE(1:ny,1:nx,1:nz))
			Lsm = SUM(TOR(1:ny,1:nx,1:nz)**2*psE(1:ny,1:nx,1:nz))
			
			CALL MPI_ALLGATHER(Lsm,1,MPI_DOUBLE_PRECISION,VaER(0:np-1),1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,errcode)
			RMS = SQRT(SUM(VaER(0:np-1))/tPE)

			if ( mod(itle,200) == 1 .AND. mod(iter,500) == 1 .AND. rank == 1 ) Print*, 'phE', iter, nlp, itle, RMS
			IF ( itle > 60000 ) EXIT

			itle = itle + 1

		ENDDO

		!! BC
		CALL BcMpiYZ(rank, rnb, cnb, ddR, ddC, ny, nx, nz, phE(0:ny+1,0:nx+1,0:nz+1) )
		phE(0:ny+1,0,0:nz+1) = BvE
		phE(0:ny+1,nx+1,0:nz+1) = phE(0:ny+1,nx,0:nz+1)

		Lsm = SUM((phE(1:ny,1:nx,1:nz)-peO(1:ny,1:nx,1:nz))**2*psE(1:ny,1:nx,1:nz))
		CALL MPI_ALLGATHER(Lsm,1,MPI_DOUBLE_PRECISION,VaER(0:np-1),1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,errcode)
		CvgL = SQRT(SUM(VaER(0:np-1))/tPE)


		if ( rank == 3 .AND. mod(iter,50) == 1 ) print*, 'Converg', iter, nlp, itlp, itle, CvgP, CvgL

		nlp = nlp + 1

	ENDDO


! ! =============================================== 
! !               _ _           _                                    _   
! !      /\      | (_)         | |                                  | |  
! !     /  \   __| |_ _   _ ___| |_    ___ _   _ _ __ _ __ ___ _ __ | |_ 
! !    / /\ \ / _` | | | | / __| __|  / __| | | | '__| '__/ _ \ '_ \| __|
! !   / ____ \ (_| | | |_| \__ \ |_  | (__| |_| | |  | | |  __/ | | | |_ 
! !  /_/    \_\__,_| |\__,_|___/\__|  \___|\__,_|_|  |_|  \___|_| |_|\__|
! !               _/ |                                                   
! !              |__/                                                    
! ! ===============================================             
             
	Lsm = SUM(Rxn(1:ny,1:nx,1:nz)*AvPx(1:ny,1:nx,1:nz))
	CALL MPI_ALLGATHER(Lsm,1,MPI_DOUBLE_PRECISION,VaER(0:np-1),1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,errcode)
	sCrnt = SUM(VaER(0:np-1))

	BvE = BvE + DSIGN(Vsr,trgI-sCrnt)*dt
	phE(1:ny,1:nx,1:nz) = phE(1:ny,1:nx,1:nz) + DSIGN(Vsr,trgI-sCrnt)*dt

	CALL BcMpiYZ(rank, rnb, cnb, ddR, ddC, ny, nx, nz, phE(0:ny+1,0:nx+1,0:nz+1) )
	phE(0:ny+1,0,1:nz) = BvE
	phE(0:ny+1,nx+1,1:nz) = phE(0:ny+1,nx,1:nz)

	if ( rank == 3 .AND. mod(iter,100) == 1 ) print*, 'Current', iter, trgI, sCrnt, BvE, Xrf

	tm = tm + dt

	IF ( Xrf > 0.2 + (cnt-1)*0.005 ) THEN
! ! 	IF ( MOD(iter,30000) == 1 ) THEN


		CALL OutPut3D_A(cnt, ny, nx, nz, 'CP'//RkNb, 2, CnP(1:ny,1:nx,1:nz))
		CALL OutPut3D_A(cnt, ny, nx, nz, 'CE'//RkNb, 2, CnE(1:ny,1:nx,1:nz))
		CALL OutPut3D_A(cnt, ny, nx, nz, 'PP'//RkNb, 2, phP(1:ny,1:nx,1:nz))
		CALL OutPut3D_A(cnt, ny, nx, nz, 'PE'//RkNb, 2, phE(1:ny,1:nx,1:nz))
	
		IF ( rank == 0 ) THEN
			CALL ReC_TmFX(cnt, TmNmRt, iter, tm, Xrf, sCrnt, trgI, BvE)
		ENDIF

		cnt = cnt + 1
	ENDIF
	
	IF ( BvE >= 0.0 .OR. Xrf > 0.9877 ) EXIT
	
	iter = iter + 1

ENDDO


CALL MPI_BARRIER(MPI_COMM_WORLD,errcode)
CALL MPI_FINALIZE(errcode)



END PROGRAM SBMES_3D_MPI_CarGr_ConjG_v01
!! ============================================================================================== !!
!! SUBROUTINES
!! ============================================================================================== !!
!! Data IOs
!! ============================================================================================== !!

SUBROUTINE OutPut3D_A(fr,py,px,pz,prfn,AB,Cn1)
IMPLICIT NONE

INTEGER,INTENT(IN) :: fr, px, py, pz, AB
REAL(KIND=8),INTENT(IN),DIMENSION(1:py,1:px,1:pz) :: Cn1
CHARACTER(LEN=6),INTENT(IN) :: prfn
CHARACTER(LEN=4) :: fs
CHARACTER(LEN=14) :: flnm

WRITE(fs,"(I4)") 1000+fr

!! case 1 ==> ascii; case 2 ==> binary
SELECT CASE (AB)
	CASE (1)
		flnm = prfn//fs//'.txt'
		OPEN(UNIT=20,FILE=flnm,FORM='FORMATTED',STATUS='REPLACE',ACTION='WRITE')
			WRITE(20,*) SNGL(Cn1)
		CLOSE(20)
	CASE (2)
		flnm = prfn//fs//'.dat'
		OPEN(UNIT=20,FILE='/mnt/scratch/hcy/SBMES_2020/2020_0922_D/Data_R1/'//flnm, &
! 		OPEN(UNIT=20,FILE='../Data_R1/'//flnm, &
			FORM='UNFORMATTED',STATUS='REPLACE',ACTION='WRITE')
! 			WRITE(20) SNGL(Cn1)
			WRITE(20) Cn1
		CLOSE(20)
END SELECT

END SUBROUTINE OutPut3D_A

!! ==================== ==================== ==================== ==================== !!
SUBROUTINE ReC_TmFX(fr, prfn, iter, tm, frx, crnt, Tg, BvE)

IMPLICIT NONE

INTEGER, INTENT(IN) :: fr, iter
REAL(KIND=8), INTENT(In) :: tm, frx, crnt, Tg, BvE
CHARACTER(LEN=6),INTENT(In) :: prfn
CHARACTER(LEN=14) :: flname

flname = prfn//'TmFx.txt'
IF ( fr == 1 ) THEN
	OPEN(UNIT=5,FILE=flname,STATUS='REPLACE',ACTION='WRITE')
		WRITE(5,*) fr, iter, tm, frx, crnt, Tg, BvE
	CLOSE(5)
ELSE
	OPEN(UNIT=6,FILE=flname,STATUS='OLD',ACTION='WRITE',POSITION='APPEND')
		WRITE(6,*) fr, iter, tm, frx, crnt, Tg, BvE
	CLOSE(6)
ENDIF


END SUBROUTINE ReC_TmFX
!! ==================== ==================== ==================== ==================== !!
!! Domain decomposition
!! ==================== ==================== ==================== ==================== !!
SUBROUTINE MyDmnCmpn_YZ(rank,tlyg,tlzg,rnb,cnb,ddR,ddC,py,lwy,upy,pz,lwz,upz)
IMPLICIT NONE
INTEGER,INTENT(IN) :: rank, tlyg, tlzg, rnb, cnb
INTEGER,INTENT(OUT) :: ddR, ddC, py,lwy, upy, pz, lwz, upz
INTEGER :: bivy, rmvy, bivz, rmvz

ddC = INT(rank/rnb)+1
ddR = MOD(rank,rnb)+1

!! determine subdomain size along Z-axis
bivz = INT(tlzg/rnb)
rmvz = MOD(tlzg,rnb)
IF ( (ddR-1) <= rmvz-1 ) THEN
	pz = bivz+1
ELSE
	pz = bivz
ENDIF
!! global index of each subdomain
lwz = (ddR-1)*pz+1
IF ( (ddR-1) >= rmvz ) lwz = lwz+rmvz
upz = lwz+(pz-1)

!! determine subdomain size along X-axis
bivy = INT(tlyg/cnb)
rmvy = MOD(tlyg,cnb)
IF ( (ddC-1) <= rmvy-1 ) THEN
	py = bivy+1
ELSE
	py = bivy
ENDIF
!! global index of each subdomain
lwy = (ddC-1)*py+1
IF ( (ddC-1) >= rmvy ) lwy = lwy+rmvy
upy = lwy+(py-1)

END SUBROUTINE MyDmnCmpn_YZ
!! ==================== ==================== ==================== ==================== !!

!! ==================== ==================== ==================== ==================== !!
SUBROUTINE BcMpiYZ(rank, rnb, cnb, ddR, ddC, py, px, pz, Udis)
IMPLICIT NONE
INCLUDE 'mpif.h'

INTEGER,INTENT(IN) :: rank, rnb, cnb, ddR, ddC, py, px, pz
REAL(KIND=8),INTENT(INOUT) :: Udis(0:py+1,0:px+1,0:pz+1)
REAL(KIND=8),DIMENSION(py, px) :: BcMpiA, BcMpiB
REAL(KIND=8),DIMENSION(px,0:pz+1) :: BcMpiC, BcMpiD
INTEGER :: errcode
INTEGER,DIMENSION(MPI_STATUS_SIZE) :: nstatus

!! apply BCs along Z-axis using MPI communication (sending up)
BcMpiA(1:py,1:px) = Udis(1:py,1:px,pz)
IF ( ddR /= rnb ) THEN
	CALL MPI_SEND(BcMpiA(1:py,1:px),py*px,MPI_DOUBLE_PRECISION,rank+1,99,MPI_COMM_WORLD,errcode)
ENDIF
IF ( ddR /= 1 ) THEN
	CALL MPI_RECV(BcMpiB(1:py,1:px),py*px,MPI_DOUBLE_PRECISION,rank-1,99,MPI_COMM_WORLD,nstatus,errcode)
ELSE
	BcMpiB(1:py,1:px) = Udis(1:py,1:px,1)
ENDIF
Udis(1:py,1:px,0) = BcMpiB(1:py,1:px)

!! apply BCs along Z-axis using MPI communication (sending down)
BcMpiA(1:py,1:px) = Udis(1:py,1:px,1)
IF ( ddR /= 1 ) THEN
	CALL MPI_SEND(BcMpiA(1:py,1:px),py*px,MPI_DOUBLE_PRECISION,rank-1,98,MPI_COMM_WORLD,errcode)
ENDIF
IF ( ddR /= rnb ) THEN
	CALL MPI_RECV(BcMpiB(1:py,1:px),py*px,MPI_DOUBLE_PRECISION,rank+1,98,MPI_COMM_WORLD,nstatus,errcode)
ELSE
	BcMpiB(1:py,1:px) = Udis(1:py,1:px,pz)
ENDIF
Udis(1:py,1:px,pz+1) = BcMpiB(1:py,1:px)

!! apply BCs along Y-axis using MPI communication (sending right)
BcMpiC(1:px,0:pz+1) = Udis(py,1:px,0:pz+1)
IF ( ddC /= cnb ) THEN
	CALL MPI_SEND(BcMpiC(1:px,0:pz+1),px*(pz+2),MPI_DOUBLE_PRECISION,rank+rnb,97,MPI_COMM_WORLD,errcode)
ENDIF
IF ( ddC /= 1 ) THEN
	CALL MPI_RECV(BcMpiD(1:px,0:pz+1),px*(pz+2),MPI_DOUBLE_PRECISION,rank-rnb,97,MPI_COMM_WORLD,nstatus,errcode)
ELSE
	BcMpiD(1:px,0:pz+1) = Udis(1,1:px,0:pz+1)
ENDIF
Udis(0,1:px,0:pz+1) = BcMpiD(1:px,0:pz+1)

!! apply BCs along X-axis using MPI communication (sending left)
BcMpiC(1:px,0:pz+1) = Udis(1,1:px,0:pz+1)
IF ( ddC /= 1 ) THEN
	CALL MPI_SEND(BcMpiC(1:px,0:pz+1),px*(pz+2),MPI_DOUBLE_PRECISION,rank-rnb,96,MPI_COMM_WORLD,errcode)
ENDIF
IF ( ddC /= cnb ) THEN
	CALL MPI_RECV(BcMpiD(1:px,0:pz+1),px*(pz+2),MPI_DOUBLE_PRECISION,rank+rnb,96,MPI_COMM_WORLD,nstatus,errcode)
ELSE
	BcMpiD(1:px,0:pz+1) = Udis(py,1:px,0:pz+1)
ENDIF
Udis(py+1,1:px,0:pz+1) = BcMpiD(1:px,0:pz+1)

END SUBROUTINE BcMpiYZ
!! ==================== ==================== ==================== ==================== !!
