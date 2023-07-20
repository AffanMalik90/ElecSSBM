!! SMOOTHED BOUNDARY METHOD ELECTROCHEMICAL SIMULATION 3D

PROGRAM SBMES_3D_AS_MPI_v1
IMPLICIT NONE

INCLUDE 'mpif.h'

INTEGER :: nx, ny, nz, i, j, k, iter, itlp, nlp, cnt
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

!! for MPI functions
INTEGER :: errcode, rank, np
INTEGER,DIMENSION(MPI_STATUS_SIZE) :: nstatus


!! -- initiate MPI
CALL MPI_INIT(errcode)
CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,errcode)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD,np,errcode)

gy = 360
gx = 118
gz = 320

!! ============================================ !!
!! -- domain decomposition, only in Y and Z directions
rnb =  36		!! rows of ranks
cnb =  24		!! columns of ranks
CALL MyDmnCmpn_YZ(rank, gy, gz, rnb, cnb, ddR, ddC, ny, lwy, upy, nz, lwz, upz)
!! y index along R; z index along C
nx = gx
!! ============================================ !!

print*, rank, ny, lwy, upy, nz, lwz, upz

dx = 1.0d-5		!! cm	double precision
dt = 0.01		!! s

zeta = 1.0
!! optimize  convergence criteria for speed
CrtE = 1.0d-6	!! convergency criterion for CnE, concentration in electrolyte
CrtP = 1.0d-8	!! convergency criterion for phP, potential in particle
CrtL = 2.0d-10	!! convergency criterion for phE, potential in electrolyte

!! optimize  convergence criteria for accuracy
! CrtE = 1.0d-11	!! convergency criterion for CnE, concentration in electrolyte
! CrtP = 1.0d-8	!! convergency criterion for phP, potential in particle
! CrtL = 1.0d-8	!! convergency criterion for phE, potential in electrolyte

rho = 0.0501		!! mole/cm^3; lattice density of particle
Cst1 = 1.6021766d-19/(1.3806488d-23*300.0)		!! constant F/(RT), e/(kT)
Frd = 96485.3365		!! Faraday constant

!!!!!!!!!!!!!!!!! tune these params !!!!!!!!!!!!!!!!!

BvP =  2.5				!! Boundary condition for phP
BvE = -0.325			!! Boundary condition for phE

Cr = 0.5 					!! charge rate (0.5 -> 2 hours to charge)
Vsr = 2.0d-2      !! voltage scanning rate

!!!!!!!!!!!!!!!!! tune these params !!!!!!!!!!!!!!!!!

!! -- LiPF6 electrolyte used. D_PF6 = 1.5d-6 and D_Li = 0.73d-6 . Newman p 284
De = 2.0d0*(0.73d-6*1.5d-6)/(0.73d-6+1.5d-6)		!! ambipolar diffusivity
Dmp = 1.5d-6 - 0.73d-6								!! D_m - D_p
Kpl = Cst1 * (0.73d-6+1.5d-6)							!! conductivity of electrolyte
t_minus = 1.5d-6/(0.73d-6+1.5d-6)					!! transference number

ALLOCATE(psP(0:ny+1,0:nx+1,0:nz+1), AvP(ny,nx,nz), ppN(ny,nx,nz,6), &
	AvPx(ny,nx,nz), psE(0:ny+1,0:nx+1,0:nz+1), peN(ny,nx,nz,6))

!! domain parameter of particle (global)
ALLOCATE(GLP(0:gy+1,0:gx+1,0:gz+1) )
OPEN(UNIT=20,FILE='/scratch/01629/hcy/NMC_2019/2019_1127_A/Geom/PS_360x118x320.dat', &
	FORM='UNFORMATTED',STATUS='OLD',ACTION='READ')
	READ(20) GLP(1:gy,1:gx,1:gz)
CLOSE(20)
! GLP(1:gy,1:gx,1:gz) = 0.5*(1.0+TANH(GLP(1:gy,1:gx,1:gz)/zeta))

!! local psi
psP(1:ny,1:nx,1:nz) = GLP(lwy:upy,1:gx,lwz:upz)

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
psE(0:ny+1,0:nx+1,0:nz+1) = psE(0:ny+1,0:nx+1,0:nz+1) + 1.0d-7

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

WHERE ( AvP(1:ny,1:nx,1:nz) < 1.0d-2 )
	AvP(1:ny,1:nx,1:nz) = 0.0d0
END WHERE

AvPx(1:ny,1:nx,1:nz) = AvP(1:ny,1:nx,1:nz)/dx

!! Li fraction, Divergence of CnP, diffusivity of Li in particle
ALLOCATE(CnP(0:ny+1,0:nx+1,0:nz+1), DvC(ny,nx,nz), DfC(0:ny+1,0:nx+1,0:nz+1))
CnP(0:ny+1,0:nx+1,0:nz+1) = 0.99	!! initial value of concentration

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
! rxn(1:ny,1:nx,1:nz) = 1.0d-7
Rxn(1:ny,1:nx,1:nz) = 0.0d0

!! potential in particle, Kap conductivity
ALLOCATE(Kap(0:ny+1,0:nx+1,0:nz+1), phP(0:ny+1,0:nx+1,0:nz+1), kpN(ny,nx,nz,6), &
	ppO(ny,nx,nz))
phP(0:ny+1,0:nx+1,0:nz+1) = BvP

!! potential in electrolyte
ALLOCATE(RHE(ny,nx,nz), phE(0:ny+1,0:nx+1,0:nz+1), peO(ny,nx,nz))
phE(0:ny+1,0:nx+1,0:nz+1) = BvE

ALLOCATE(tmp(ny,nx,nz) )

trgI = -0.75*tPP*rho/(3600/Cr) !! target current

tm = 0.0d0

!! -- file names for output data
WRITE(RkNb,"(I4)") 1000+rank

! CALL OutPut3D_A(1, ny, nx, nz, 'PS'//RkNb, 2, psP(1:ny,1:nx,1:nz))

if ( rank == 3) print*, ny, nx, nz



!! time evolution
cnt = 1
iter = 1
! DO iter = 1, 180001
DO WHILE ( Xrf > 0.2 )


	!! cathode particle
	!! diffusivity as a function of concentration
	DfC(0:ny+1,0:nx+1,0:nz+1) = (0.0277-0.0840*CnP(0:ny+1,0:nx+1,0:nz+1) + &
		0.1003*CnP(0:ny+1,0:nx+1,0:nz+1)**2)*1.0d-8		!! cm^2/s

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
!
	!! no flux BC on X sides
	CALL BcMpiYZ(rank, rnb, cnb, ddR, ddC, ny, nx, nz, CnP(0:ny+1,0:nx+1,0:nz+1) )
	CnP(0:ny+1,0,0:nz+1) = CnP(0:ny+1,1,0:nz+1)
	CnP(0:ny+1,nx+1,0:nz+1) = CnP(0:ny+1,nx,0:nz+1)

	!! average concentration in the particle
	Lsm = SUM(CnP(1:ny,1:nx,1:nz)*psP(1:ny,1:nx,1:nz))
	CALL MPI_ALLGATHER(Lsm,1,MPI_DOUBLE_PRECISION,VaER(0:np-1),1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,errcode)
	Xrf = SUM(VaER(0:np-1))/tPP

! 	! CALL OutPut3D_A(iter, ny, nx, nz, 'CONPRT', 2, CnP(1:ny,1:nx,1:nz))

	!! electrolyte concentration
	CEo(1:ny,1:nx,1:nz) = CnE(1:ny,1:nx,1:nz)
	RHS(1:ny,1:nx,1:nz) = psE(1:ny,1:nx,1:nz)*CEo(1:ny,1:nx,1:nz) - dt*(Rxn(1:ny,1:nx,1:nz))*AvPx(1:ny,1:nx,1:nz)*t_minus

	! influx of Li
	Lsm = SUM(Rxn(1:ny,1:nx,1:nz)*AvPx(1:ny,1:nx,1:nz))*(dx**3)   !!/(ny*nz*dx**2)
	CALL MPI_ALLGATHER(Lsm,1,MPI_DOUBLE_PRECISION,VaER(0:np-1),1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,errcode)
	infx = SUM(VaER(0:np-1))/(gy*gz*dx**2)


	IF ( MOD(iter,20) == 1 .AND. rank == 1 ) print*, 'infx = ', infx

	RMS = 1.0
	itlp = 1
	!! internal loop
! 	DO itlp = 1, 10001
	DO WHILE ( RMS > CrtE )
!
		!! point-wise relaxation
		CnE(1:ny,1:nx,1:nz) = (RHS(1:ny,1:nx,1:nz) + &
		   dt*De/dx**2*(peN(1:ny,1:nx,1:nz,2)*CnE(2:ny+1,1:nx,1:nz) + peN(1:ny,1:nx,1:nz,1)*CnE(0:ny-1,1:nx,1:nz) + &
						peN(1:ny,1:nx,1:nz,4)*CnE(1:ny,2:nx+1,1:nz) + peN(1:ny,1:nx,1:nz,3)*CnE(1:ny,0:nx-1,1:nz) + &
						peN(1:ny,1:nx,1:nz,6)*CnE(1:ny,1:nx,2:nz+1) + peN(1:ny,1:nx,1:nz,5)*CnE(1:ny,1:nx,0:nz-1)))/ &
			(psE(1:ny,1:nx,1:nz) + dt*De/dx**2*(peN(1:ny,1:nx,1:nz,2) + peN(1:ny,1:nx,1:nz,1) + &
												peN(1:ny,1:nx,1:nz,4) + peN(1:ny,1:nx,1:nz,3) + &
												peN(1:ny,1:nx,1:nz,6) + peN(1:ny,1:nx,1:nz,5)))

		!! boundary condition
		CALL BcMpiYZ(rank, rnb, cnb, ddR, ddC, ny, nx, nz, CnE(0:ny+1,0:nx+1,0:nz+1) )
		CnE(0:ny+1,0,0:nz+1) = infx*t_minus*dx/De + CnE(0:ny+1,1,0:nz+1)
		CnE(0:ny+1,nx+1,0:nz+1) = CnE(0:ny+1,nx,0:nz+1)


		!! residual
		RES(1:ny,1:nx,1:nz) = psE(1:ny,1:nx,1:nz)*CnE(1:ny,1:nx,1:nz) - &
			dt*De/dx**2*(peN(1:ny,1:nx,1:nz,2)*(CnE(2:ny+1,1:nx,1:nz)-CnE(1:ny,1:nx,1:nz)) - &
						 peN(1:ny,1:nx,1:nz,1)*(CnE(1:ny,1:nx,1:nz)-CnE(0:ny-1,1:nx,1:nz))	+ &
						 peN(1:ny,1:nx,1:nz,4)*(CnE(1:ny,2:nx+1,1:nz)-CnE(1:ny,1:nx,1:nz)) - &
						 peN(1:ny,1:nx,1:nz,3)*(CnE(1:ny,1:nx,1:nz)-CnE(1:ny,0:nx-1,1:nz)) + &
						 peN(1:ny,1:nx,1:nz,6)*(CnE(1:ny,1:nx,2:nz+1)-CnE(1:ny,1:nx,1:nz)) - &
						 peN(1:ny,1:nx,1:nz,5)*(CnE(1:ny,1:nx,1:nz)-CnE(1:ny,1:nx,0:nz-1))) - &
						 RHS(1:ny,1:nx,1:nz)

		Lsm = SUM(RES(1:ny,1:nx,1:nz)**2*psE(1:ny,1:nx,1:nz))
		CALL MPI_ALLGATHER(Lsm,1,MPI_DOUBLE_PRECISION,VaER(0:np-1),1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,errcode)
		RMS = SQRT(SUM(VaER(0:np-1))/tPE)

! 		! CALL OutPut3D_A(iter, ny, nx, nz, 'CONELY', 2, CnE(1:ny,1:nx,1:nz))
! 		! CALL OutPut3D_A(iter, ny, nx, nz, 'PSIELY', 2, psE(1:ny,1:nx,1:nz))

		if ( mod(itlp,2000) == 1 .AND. mod(iter,50) == 1 .AND. rank == 2 ) Print*, 'CnE', iter, itlp, RMS
		IF ( itlp > 50000 ) EXIT
!
		itlp = itlp + 1
  	ENDDO


	! electrolyte potential !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

	!! particle conductivity
	Kap(1:ny,1:nx,1:nz) = 0.01929 + 0.7045*TANH(2.399*CnP(1:ny,1:nx,1:nz)) - 0.7238*TANH(2.412*CnP(1:ny,1:nx,1:nz)) - 4.2106d-6

	CALL BcMpiYZ(rank, rnb, cnb, ddR, ddC, ny, nx, nz, Kap(0:ny+1,0:nx+1,0:nz+1) )
	Kap(0:ny+1,0,0:nz+1) = Kap(0:ny+1,1,0:nz+1)
	Kap(0:ny+1,nx+1,0:nz+1) = Kap(0:ny+1,nx,0:nz+1)

	kpN(1:ny,1:nx,1:nz,1) = (Kap(0:ny-1,1:nx,1:nz)+Kap(1:ny,1:nx,1:nz))*0.5d0
	kpN(1:ny,1:nx,1:nz,2) = (Kap(2:ny+1,1:nx,1:nz)+Kap(1:ny,1:nx,1:nz))*0.5d0
	kpN(1:ny,1:nx,1:nz,3) = (Kap(1:ny,0:nx-1,1:nz)+Kap(1:ny,1:nx,1:nz))*0.5d0
	kpN(1:ny,1:nx,1:nz,4) = (Kap(1:ny,2:nx+1,1:nz)+Kap(1:ny,1:nx,1:nz))*0.5d0
	kpN(1:ny,1:nx,1:nz,5) = (Kap(1:ny,1:nx,0:nz-1)+Kap(1:ny,1:nx,1:nz))*0.5d0
	kpN(1:ny,1:nx,1:nz,6) = (Kap(1:ny,1:nx,2:nz+1)+Kap(1:ny,1:nx,1:nz))*0.5d0

	!! exchange current density
	WHERE ( AvP(1:ny,1:nx,1:nz) > 1.0d-2 )
		ioC(1:ny,1:nx,1:nz) = 10.0**(-0.2*(CnP(1:ny,1:nx,1:nz)-0.37)-1.559-0.9376*TANH(8.961*CnP(1:ny,1:nx,1:nz)-3.195))
		ioC = ioC * 1.0d-3 !! in unit of Amp

		OCV(1:ny,1:nx,1:nz) = 1.095*CnP(1:ny,1:nx,1:nz)**2-8.324d-7*EXP(14.31*CnP(1:ny,1:nx,1:nz))+4.692*EXP(-0.5389*CnP(1:ny,1:nx,1:nz))

		Kfw(1:ny,1:nx,1:nz) = ioC(1:ny,1:nx,1:nz)/(Frd*0.001              )*EXP( 0.5*Cst1*OCV(1:ny,1:nx,1:nz))
		Kbw(1:ny,1:nx,1:nz) = ioC(1:ny,1:nx,1:nz)/(Frd*CnP(1:ny,1:nx,1:nz))*EXP(-0.5*Cst1*OCV(1:ny,1:nx,1:nz))
	END WHERE

	CvgP = 1.0
	CvgL = 1.0
	nlp = 1
! 	! DO WHILE ( CvgP > 1.0d-11 .AND. CvgL > 1.0d-11 )
	DO WHILE ( CvgP > 1.0d-8 .AND. CvgL > 1.0d-8 )
! 	DO nlp = 1, 101

		WHERE ( AvP(1:ny,1:nx,1:nz) > 1.0d-2 )
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

		RMS = 1.0
		itlp = 1
		DO WHILE ( RMS > CrtP )
! 		DO itlp = 1, 101
			tmp(1:ny,1:nx,1:nz) = phP(1:ny,1:nx,1:nz)

			phP(1:ny,1:nx,1:nz) = (RHS(1:ny,1:nx,1:nz) - &
				1.0/dx**2*(ppN(1:ny,1:nx,1:nz,2)*kpN(1:ny,1:nx,1:nz,2)*phP(2:ny+1,1:nx,1:nz) + &
						   		 ppN(1:ny,1:nx,1:nz,1)*kpN(1:ny,1:nx,1:nz,1)*phP(0:ny-1,1:nx,1:nz) + &
						   		 ppN(1:ny,1:nx,1:nz,4)*kpN(1:ny,1:nx,1:nz,4)*phP(1:ny,2:nx+1,1:nz) + &
						   		 ppN(1:ny,1:nx,1:nz,3)*kpN(1:ny,1:nx,1:nz,3)*phP(1:ny,0:nx-1,1:nz) + &
								 ppN(1:ny,1:nx,1:nz,6)*kpN(1:ny,1:nx,1:nz,6)*phP(1:ny,1:nx,2:nz+1) + &
								 ppN(1:ny,1:nx,1:nz,5)*kpN(1:ny,1:nx,1:nz,5)*phP(1:ny,1:nx,0:nz-1)))/ &
			(-1.0/dx**2*(ppN(1:ny,1:nx,1:nz,2)*kpN(1:ny,1:nx,1:nz,2) + ppN(1:ny,1:nx,1:nz,1)*kpN(1:ny,1:nx,1:nz,1) + &
						 ppN(1:ny,1:nx,1:nz,4)*kpN(1:ny,1:nx,1:nz,4) + ppN(1:ny,1:nx,1:nz,3)*kpN(1:ny,1:nx,1:nz,3) + &
						 ppN(1:ny,1:nx,1:nz,6)*kpN(1:ny,1:nx,1:nz,6) + ppN(1:ny,1:nx,1:nz,5)*kpN(1:ny,1:nx,1:nz,5)))
!
!
			!! BC
			CALL BcMpiYZ(rank, rnb, cnb, ddR, ddC, ny, nx, nz, phP(0:ny+1,0:nx+1,0:nz+1) )
			phP(0:ny+1,0,0:nz+1) = phP(0:ny+1,1,0:nz+1)
			phP(0:ny+1,nx+1,0:nz+1) = BvP*psP(0:ny+1,nx+1,0:nz+1)+phP(0:ny+1,nx,0:nz+1)*(1.0-psP(0:ny+1,nx+1,0:nz+1))

			!! residual
! 			RES(1:ny,1:nx,1:nz) = 1.0/dx**2* &
! 				(ppN(1:ny,1:nx,1:nz,2)*kpN(1:ny,1:nx,1:nz,2)*(phP(2:ny+1,1:nx,1:nz)-phP(1:ny,1:nx,1:nz)) - &
! 				 ppN(1:ny,1:nx,1:nz,1)*kpN(1:ny,1:nx,1:nz,1)*(phP(1:ny,1:nx,1:nz)-phP(0:ny-1,1:nx,1:nz)) + &
! 				 ppN(1:ny,1:nx,1:nz,4)*kpN(1:ny,1:nx,1:nz,4)*(phP(1:ny,2:nx+1,1:nz)-phP(1:ny,1:nx,1:nz)) - &
! 				 ppN(1:ny,1:nx,1:nz,3)*kpN(1:ny,1:nx,1:nz,3)*(phP(1:ny,1:nx,1:nz)-phP(1:ny,0:nx-1,1:nz)) + &
! 				 ppN(1:ny,1:nx,1:nz,6)*kpN(1:ny,1:nx,1:nz,6)*(phP(1:ny,1:nx,2:nz+1)-phP(1:ny,1:nx,1:nz)) - &
! 				 ppN(1:ny,1:nx,1:nz,5)*kpN(1:ny,1:nx,1:nz,5)*(phP(1:ny,1:nx,1:nz)-phP(1:ny,1:nx,0:nz-1))) - &
! 				 RHS(1:ny,1:nx,1:nz)



			RES(1:ny,1:nx,1:nz) = phP(1:ny,1:nx,1:nz) - tmp(1:ny,1:nx,1:nz)

			Lsm = SUM(RES(1:ny,1:nx,1:nz)**2*psP(1:ny,1:nx,1:nz))
			CALL MPI_ALLGATHER(Lsm,1,MPI_DOUBLE_PRECISION,VaER(0:np-1),1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,errcode)
			RMS = SQRT(SUM(VaER(0:np-1))/tPP)/BvP

			if ( mod(itlp,2000) == 1 .AND. mod(iter,50) == 1 .AND. rank == 1 ) Print*, 'phP', iter, nlp, itlp, RMS
			IF ( itlp > 30000 ) EXIT

			itlp = itlp + 1

		ENDDO

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


		Lsm = SUM((phP(1:ny,1:nx,1:nz)-ppO(1:ny,1:nx,1:nz))**2*psP(1:ny,1:nx,1:nz))
		CALL MPI_ALLGATHER(Lsm,1,MPI_DOUBLE_PRECISION,VaER(0:np-1),1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,errcode)
		CvgP = SQRT(SUM(VaER(0:np-1))/tPP)

		!! electrolyte potential
		peO(1:ny,1:nx,1:nz) = phE(1:ny,1:nx,1:nz)
		RHS(1:ny,1:nx,1:nz) = (Rxn(1:ny,1:nx,1:nz))*AvPx(1:ny,1:nx,1:nz) + RHE(1:ny,1:nx,1:nz) !! flux

		RMS = 1.0
		itlp = 1
		DO WHILE ( RMS > CrtL )
! 		DO itlp = 1, 101

			tmp(1:ny,1:nx,1:nz) = phE(1:ny,1:nx,1:nz)

			phE(1:ny,1:nx,1:nz) = (RHS(1:ny,1:nx,1:nz) - &
				Kpl/dx**2*(peN(1:ny,1:nx,1:nz,2)*ceN(1:ny,1:nx,1:nz,2)*phE(2:ny+1,1:nx,1:nz) + &
						   peN(1:ny,1:nx,1:nz,1)*ceN(1:ny,1:nx,1:nz,1)*phE(0:ny-1,1:nx,1:nz) + &
						   peN(1:ny,1:nx,1:nz,4)*ceN(1:ny,1:nx,1:nz,4)*phE(1:ny,2:nx+1,1:nz) + &
						   peN(1:ny,1:nx,1:nz,3)*ceN(1:ny,1:nx,1:nz,3)*phE(1:ny,0:nx-1,1:nz) + &
						   peN(1:ny,1:nx,1:nz,6)*ceN(1:ny,1:nx,1:nz,6)*phE(1:ny,1:nx,2:nz+1) + &
						   peN(1:ny,1:nx,1:nz,5)*ceN(1:ny,1:nx,1:nz,5)*phE(1:ny,1:nx,0:nz-1)))/ &
			  (-Kpl/dx**2*(peN(1:ny,1:nx,1:nz,2)*ceN(1:ny,1:nx,1:nz,2) + peN(1:ny,1:nx,1:nz,1)*ceN(1:ny,1:nx,1:nz,1) + &
						   peN(1:ny,1:nx,1:nz,4)*ceN(1:ny,1:nx,1:nz,4) + peN(1:ny,1:nx,1:nz,3)*ceN(1:ny,1:nx,1:nz,3) + &
						   peN(1:ny,1:nx,1:nz,6)*ceN(1:ny,1:nx,1:nz,6) + peN(1:ny,1:nx,1:nz,5)*ceN(1:ny,1:nx,1:nz,5)))

			!! BC
			CALL BcMpiYZ(rank, rnb, cnb, ddR, ddC, ny, nx, nz, phE(0:ny+1,0:nx+1,0:nz+1) )
			phE(0:ny+1,0,0:nz+1) = BvE
			phE(0:ny+1,nx+1,0:nz+1) = phE(0:ny+1,nx,0:nz+1)


			!! residual
! 			RES(1:ny,1:nx,1:nz) = Kpl/dx**2*( &
! 				peN(1:ny,1:nx,1:nz,2)*ceN(1:ny,1:nx,1:nz,2)*(phE(2:ny+1,1:nx,1:nz)-phE(1:ny,1:nx,1:nz)) - &
! 				peN(1:ny,1:nx,1:nz,1)*ceN(1:ny,1:nx,1:nz,1)*(phE(1:ny,1:nx,1:nz)-phE(0:ny-1,1:nx,1:nz)) + &
! 				peN(1:ny,1:nx,1:nz,4)*ceN(1:ny,1:nx,1:nz,4)*(phE(1:ny,2:nx+1,1:nz)-phE(1:ny,1:nx,1:nz)) - &
! 				peN(1:ny,1:nx,1:nz,3)*ceN(1:ny,1:nx,1:nz,3)*(phE(1:ny,1:nx,1:nz)-phE(1:ny,0:nx-1,1:nz)) + &
! 				peN(1:ny,1:nx,1:nz,6)*ceN(1:ny,1:nx,1:nz,6)*(phE(1:ny,1:nx,2:nz+1)-phE(1:ny,1:nx,1:nz)) - &
! 				peN(1:ny,1:nx,1:nz,5)*ceN(1:ny,1:nx,1:nz,5)*(phE(1:ny,1:nx,1:nz)-phE(1:ny,1:nx,0:nz-1))) - &
! 				RHS(1:ny,1:nx,1:nz)

			RES(1:ny,1:nx,1:nz) = phE(1:ny,1:nx,1:nz) - tmp(1:ny,1:nx,1:nz)

			Lsm = SUM(RES(1:ny,1:nx,1:nz)**2*psE(1:ny,1:nx,1:nz))
			CALL MPI_ALLGATHER(Lsm,1,MPI_DOUBLE_PRECISION,VaER(0:np-1),1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,errcode)
			RMS = SQRT(SUM(VaER(0:np-1))/tPE)


! 			RMS = SQRT(SUM(RES(1:ny,1:nx,1:nz)**2*psE(1:ny,1:nx,1:nz))/tPE)

			if ( mod(itlp,2000) == 1 .AND. mod(iter,50) == 1 .AND. rank == 0 ) Print*, 'phE', iter, nlp, itlp, RMS
			IF ( itlp > 60000 ) EXIT

			itlp = itlp + 1

		ENDDO

		Lsm = SUM((phE(1:ny,1:nx,1:nz)-peO(1:ny,1:nx,1:nz))**2*psE(1:ny,1:nx,1:nz))
		CALL MPI_ALLGATHER(Lsm,1,MPI_DOUBLE_PRECISION,VaER(0:np-1),1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,errcode)
		CvgL = SQRT(SUM(VaER(0:np-1))/tPE)

! 		CvgL = SQRT(SUM((phE(1:ny,1:nx,1:nz)-peO(1:ny,1:nx,1:nz))**2*psE(1:ny,1:nx,1:nz))/tPE)
! 		CvgL = CvgL
!
		if ( rank == 0 .AND. mod(iter,50) == 1 ) print*, 'Converg', iter, nlp, CvgP, CvgL

		nlp = nlp + 1

	ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! check this !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	Lsm = SUM(Rxn(1:ny,1:nx,1:nz)*AvPx(1:ny,1:nx,1:nz))
	CALL MPI_ALLGATHER(Lsm,1,MPI_DOUBLE_PRECISION,VaER(0:np-1),1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,errcode)
	sCrnt = SUM(VaER(0:np-1))

	BvE = BvE + DSIGN(Vsr,trgI-sCrnt)*dt
	phE(1:ny,1:nx,1:nz) = phE(1:ny,1:nx,1:nz) + DSIGN(Vsr,trgI-sCrnt)*dt

	CALL BcMpiYZ(rank, rnb, cnb, ddR, ddC, ny, nx, nz, phE(0:ny+1,0:nx+1,0:nz+1) )
	phE(0:ny+1,0,1:nz) = BvE
	phE(0:ny+1,nx+1,1:nz) = phE(0:ny+1,nx,1:nz)

!
! 	phE(0,1:nx,1:nz) = phE(1,1:nx,1:nz)
! 	phE(ny+1,1:nx,1:nz) = phE(ny,1:nx,1:nz)
! 	phE(0:ny+1,0,1:nz) = BvE
! 	phE(0:ny+1,nx+1,1:nz) = phE(0:ny+1,nx,1:nz)
! 	phE(0:ny+1,0:nx+1,0) = phE(0:ny+1,0:nx+1,1)
! 	phE(0:ny+1,0:nx+1,nz+1) = phE(0:ny+1,0:nx+1,nz)

	if ( rank == 3 .AND. mod(iter,100) == 1 ) print*, 'Current', iter, trgI, sCrnt, BvE, Xrf

	tm = tm + dt

	IF ( Xrf < 0.99 - (cnt-1)*0.005 ) THEN
! 	IF ( MOD(iter,30000) == 1 ) THEN


		CALL OutPut3D_A(cnt, ny, nx, nz, 'CP'//RkNb, 2, CnP(1:ny,1:nx,1:nz))
		CALL OutPut3D_A(cnt, ny, nx, nz, 'CE'//RkNb, 2, CnE(1:ny,1:nx,1:nz))
		CALL OutPut3D_A(cnt, ny, nx, nz, 'PP'//RkNb, 2, phP(1:ny,1:nx,1:nz))
		CALL OutPut3D_A(cnt, ny, nx, nz, 'PE'//RkNb, 2, phE(1:ny,1:nx,1:nz))
	
		IF ( rank == 0 ) THEN
			CALL ReC_TmFX(cnt, 'NMC3DB', iter, tm, Xrf, sCrnt, trgI, BvE)
		ENDIF

		cnt = cnt + 1
	ENDIF
	
	IF ( BvE <= -2.0 ) EXIT
	
	iter = iter + 1

ENDDO


CALL MPI_BARRIER(MPI_COMM_WORLD,errcode)
CALL MPI_FINALIZE(errcode)



END PROGRAM SBMES_3D_AS_MPI_v1
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
		OPEN(UNIT=20,FILE='/scratch/01629/hcy/NMC_2019/2019_1127_A/Data_R2/'//flnm, &
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
