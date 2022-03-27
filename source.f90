 !************************************************************
!*															*
!*    2-Dimensional Coupling Bingham and Fredrickson Model  *
!*                                                          *
!*    Modified Bingham Fluid and Fredrickson Model          *
!*    Finite Element Method                                 *
!*    Triangular Linear Mesh                                *
!*    Velocity Correction Method                            *
!*    SUPG Version                                          *
!*    Runge-Kutta method (2nd order)                        *
!*                                                          *
!*    Dimensional form                                      *
!*                                                          *
!*    Coded by Takuro TOKUNAGA, Feb,06,2011                 *
!*                                                          *
!************************************************************

!*********************************
!*  Main Routine                 *
!*  Last Modified : DEC 21,2010  *
!*********************************


!*     use portlib

!*     [CAUTION]
!*        The statement "use portlib" is for Compaq (Visual) 
!*        Fortran. If your system works on UNIX, remove this 
!*        statement.

!* rewrite common to module 10/13
!* finished rewrite common to module 10/19

MODULE PARAMETERS

IMPLICIT NONE

	  integer,PARAMETER :: COUNT=100,FNUMBER=130
	  integer,PARAMETER :: EMAX=5000,NMAX=5000
      integer,PARAMETER :: BCMAX1=500,BCMAX2=500,BCMAX3=200
      integer,PARAMETER :: BCMAX4=200,BCMAX5=200,BCMAX6=200
      integer,PARAMETER :: BCMAXN=50
	  integer NE,NP
	  integer :: NF=0,COUNT2,FEM=1
	  integer NOP(EMAX,3)
	  double precision  CORD(NMAX,2)
	  double precision  TY,MM,N0,DENS				  
      double precision  EHX(EMAX,3,3),EHY(EMAX,3,3)		  !/MTR1/
      double precision  EKX(EMAX,3,3,3),EKY(EMAX,3,3,3)   !/MTR2/
      double precision  EMM(EMAX,3)						  !/MTR3/
      double precision  ESS(EMAX,3,3)					  !/MTR4/
      double precision  LMM(NMAX),ILMM(NMAX)			  !/MTR5/
      double precision  EB(EMAX,3),EC(EMAX,3),AR(EMAX)	  !/MTR6/
      double precision  VX0(NMAX),VY0(NMAX)			  	  !/V0/
      double precision  VX1(NMAX),VY1(NMAX)			  	  !/V1/
      double precision  PR0(NMAX),PR1(NMAX)			      !/PR/
	  double precision  TXX0(NMAX),TXY0(NMAX),TYY0(NMAX)
	  double precision  TXX1(NMAX),TXY1(NMAX),TYY1(NMAX),FNSD(NMAX)
	  double precision  TXXBing(NMAX),TXYBing(NMAX),TYYBing(NMAX)
	  double precision  VTXX(NMAX),VTXY(NMAX),VTYY(NMAX)
	  
	  double precision  PHI0(NMAX),PHI1(NMAX)
	  double precision  PHI(NMAX),VPHI(NMAX),ETA(NMAX)
      double precision  VVX(NMAX),VVY(NMAX)				  !/VEC1/
      double precision  VBC1(BCMAX1),NBC1                 !/BC1/
      double precision  VBC2(BCMAX2),NBC2                 !/BC2/
      double precision  VBC3(BCMAX3),NBC3                 !/BC3/
	  double precision  VBC4(BCMAX4),NBC4                 !/BC4/
      double precision  VBC5(BCMAX5),NBC5                 !/BC5/
      double precision  VBC6(BCMAX6),NBC6                 !/BC6/
      double precision  VBCN(BCMAXN),NBCN                 !/BCN/
	  double precision  VBCF(BCMAX4),NBCF              	  !/BCF/
	  integer NNBC1(BCMAX1)
	  integer NNBC2(BCMAX2)
	  integer NNBC3(BCMAX3)
	  integer NNBC4(BCMAX4)
	  integer NNBC5(BCMAX5)
	  integer NNBC6(BCMAX6)
	  integer NNBCN(BCMAXN,2)
	  integer NNBCF(BCMAX4)

      double precision  XX(NMAX)	  !/DP/
	  double precision  INTVAL        !the interval to output data
		
      !* model parameters *
	  double precision F1,F0,LAMBDA,K0

END MODULE PARAMETERS


      PROGRAM MAIN
      use parameters
!*
      IMPLICIT NONE
!*
!*     [CAUTION]
!*        The option "NONE" is available in Fortran 90, so you
!*        have to comment out this line  if you want to use  a
!*        FORTRAN 77 compiler.

!*	   [postscript]
!*        this program is rewrited to f90


       INTEGER I	
	   INTEGER NCOUNT
       DOUBLE PRECISION TIME
       DOUBLE PRECISION STD,RES,TEND,DT,VMEAN
       CHARACTER FNAME*20,FNAME2*20,ONAME*20,CMNT*25

	   REAL time0,time1,cpu,cpu0,cpu1,dtime

       dimension time0(2),time1(2)

	
!*  Title Page  *
      CALL TITLE

      STD=1.0D-08

!*  Model Parameters for PTT Model  *
      CALL INPUT1

!*  Mesh Data  *
      CALL INPUT2(DT)

!*  Initial Data  *
      CALL INPUT3(TIME,NP)

!*  Boundary Conditions  *
      CALL INPUT4(VMEAN)

      WRITE(*,*) ' Enter Final Time: '
      READ(5,*) TEND

      WRITE(*,*) ' Enter Output Data File Name: '
	  WRITE(*,*) ' (uo to 20 characters)          : '
      READ(5,*) ONAME

	  WRITE(*,*) ' Enter Output Result Data File Name: '
	  WRITE(*,*) ' (uo to 20 characters)          : '
      READ(5,*) FNAME

	  WRITE(*,*) ' Enter Output Result Data File Name: '
	  WRITE(*,*) ' (uo to 20 characters)          : '
      READ(5,*) FNAME2

	  WRITE(*,*) ' Enter Output Interval (step): '
	  read(5,*) INTVAL

      WRITE(*,*) ' Enter Comment within 25 Letters: '
      READ(5,*) CMNT


!* Open File *
		OPEN(UNIT=30, FILE=FNAME,STATUS='UNKNOWN')
	
!*  Information  *
        CALL INFO(NE,NP,DT,TIME,TEND,ONAME,FNAME,CMNT)

!*  Beginning of Measurement of CPU Time
      cpu0=dtime(time0)

!*  Calculate Element Matrices  *
      DO 10 I=1,NP
        LMM(I) = 0.D0
   10 CONTINUE

      DO 20 I=1,NE
        CALL MKMTR(I)
   20 CONTINUE

!*  Inverse Matrix of LMM
      DO 30 I=1,NP
        ILMM(I) = 1.D0 / LMM(I)
   30 CONTINUE

!*  Time Loop Begins  *

      NCOUNT = 0
  999 TIME = TIME + DT
      NCOUNT = NCOUNT + 1
	  COUNT2 = COUNT2 + 1

      IF (NCOUNT.EQ.COUNT) THEN
        WRITE(*,*) TIME,RES
        NCOUNT = 0
      ENDIF

!***** STEP 1a/1b *****
      CALL STEP1(DT)

!***** STEP 2 *****
      CALL STEP2(DT)

!***** STEP 3 *****
      CALL STEP3(DT,DENS)

!*  Convergence Check  *
      CALL CONV(RES,NP,VMEAN)

!*      IF ((TIME.LT.TEND).AND.(RES.GT.STD)) THEN *
	    IF (TIME.LT.TEND) THEN

!*	Update Calculated Data Per Step (added 10/10/20) *
		CALL UPDATE_STEP(NE,NP)

		if (COUNT2.EQ.INTVAL) THEN
			NF=NF+1 
			COUNT2=0 
			CALL OUTPUT_FEM(NP,NE,FEM,TIME,TEND,DT,INTVAL,FNAME)
			CALL OUTPUT_ALL(FNAME2,NE,NP,TIME,NF,DT)
			FEM=FEM+1
		END IF 
			GOTO 999
		ELSE

!*  Output Results (Output interval states) *
			CALL OUTPUT_FEM(NP,NE,FEM,TIME,TEND,DT,INTVAL,FNAME)		
			CALL OUTPUT_ALL(FNAME2,NE,NP,TIME,NF,DT)	
		END IF
	
!*  Output CPU Time  *
      cpu1= dtime(time1)
      cpu = cpu1
      write(*,*) 'CPU Time: ',cpu,' [s]'

	  write(*,*) 'ALL DONE'

	  CLOSE(30)

      STOP
      END 


!*********************************
!*  Display Title Page           *
!*  Last Modified : NOV 12,2010  *
!*********************************

      SUBROUTINE TITLE

      WRITE(*,*) '***************************************************'
      WRITE(*,*) '*                                                 *'
      WRITE(*,*) '*  2-Dimensional Modified Bingham and Fredrickson *'
      WRITE(*,*) '*                                                 *'
      WRITE(*,*) '*  Source : source.f90	                        *'
      WRITE(*,*) '*  Version: 1.00 (Nov, 2010)                      *'
      WRITE(*,*) '*                                                 *'
      WRITE(*,*) '*  Coded by                                       *'
      WRITE(*,*) '*     Takuro TOKUNAGA,      ###########           *'
      WRITE(*,*) '*     Yamamoto Laboratory,       ######           *'
      WRITE(*,*) '*     Osaka University              ###           *'
      WRITE(*,*) '*                                                 *'
      WRITE(*,*) '***************************************************'

      RETURN
      END


!*********************************
!*  Input Model Parameter        *
!*  Last Modified : DEC 21,2010  *
!*********************************

     SUBROUTINE INPUT1
      use parameters
      IMPLICIT NONE
  
      WRITE(*,100)
  100 FORMAT('* Enter Parameter:'/)

      WRITE(*,*) ' YIELD STRESS: '
      READ(*,*) TY

      WRITE(*,*) ' POWER INDEX 1: '
      READ(*,*) MM
	
	  WRITE(*,*) ' FLUIDITY ZERO S-R: '
      READ(*,*) F0

      WRITE(*,*) ' FLUIDITY HIGH S-R: '
      READ(*,*) F1
      
      WRITE(*,*) ' STRUCTURAL RELAXATION TIEM: '
      READ(*,*) LAMBDA
		
	  WRITE(*,*) ' KINETIC CONSTANT: '
      READ(*,*) K0
	
	  WRITE(*,*) ' DENSITY: '
      READ(*,*) DENS

	  RETURN

      STOP
      END


!*********************************
!*  Input Element Information    *
!*  Last Modified : AUG 9,1999   *
!*********************************

      SUBROUTINE INPUT2(DT) 
      use parameters
      IMPLICIT NONE

      INTEGER I,N
      DOUBLE PRECISION DT
      CHARACTER CMNT*25,MNAME*20

      WRITE(*,100)
      WRITE(*,110)
  100 FORMAT(/'[ Mesh Data ]')
  110 FORMAT('* Enter Mesh File Name:'/)

      READ(5,*) MNAME

      OPEN(UNIT=15, ERR=90, FILE=MNAME, STATUS='OLD')

      READ(15,'(A25)') CMNT
      WRITE(*,120) CMNT
  120 FORMAT('>>>> Comment: ',A)

      READ(15,*) NP
      READ(15,*) NE
      READ(15,*) DT

!C      dt=dt*0.5d-01

      READ(15,'(A80)')
      READ(15,*) (N,CORD(N,1),CORD(N,2),I=1,NP)
      READ(15,'(A80)')

      READ(15,*) (N,NOP(N,1),NOP(N,2),NOP(N,3),I=1,NE)

!C  100 FORMAT(I4,2E15.7/)
!C  110 FORMAT(4I4/)

      CLOSE(15)
      RETURN

   90 WRITE(*,130) MNAME
  130 FORMAT(/'>>>>>> Error (1): Cannot find ',A/)

      STOP
      END

!*********************************
!*  Input Initial Conditions     *
!*  Last Modified : NOV 05,2010  *
!*********************************

      SUBROUTINE INPUT3(TIME,NPSUB)
      
      use parameters
IMPLICIT NONE
      INTEGER I,N,NPSUB
      DOUBLE PRECISION TIME
      CHARACTER CMNT*25,INAME*20
      
	  WRITE(*,100)      
	  WRITE(*,110)
  100 FORMAT(/'[ Initial Data ]')
  110 FORMAT('* Enter Initial Data File Name:'/)

      READ(5,*) INAME

      IF ((INAME.EQ.'REST').OR.(INAME.EQ.'rest')) THEN

!*       When Fluid is at Rest at Initial Time

        WRITE(*,120) 'Fluid is at rest.'
        TIME = 0.D0

        DO 10 I=1,NPSUB
           VX0(I) = 0.D0
           VY0(I) = 0.D0
           PR0(I) = 0.D0
           VX1(I) = 0.D0
           VY1(I) = 0.D0
           PR1(I) = 0.D0
	      TXX0(I) = 0.D0
		  TXY0(I) = 0.D0
		  TYY0(I) = 0.D0
		  TXX1(I) = 0.D0
		  TXY1(I) = 0.D0
		  TYY1(I) = 0.D0
		  PHI0(I) = F0
		  PHI1(I) = F1
   10   CONTINUE

      ELSE

        OPEN(UNIT=20, ERR=90, FILE=INAME, STATUS='OLD')

!*  Initial Time  *

        READ(20,'(A25)') CMNT
        WRITE(*,120) CMNT
  120   FORMAT('>>>> Comment: ',A)

        READ(20,*) TIME

!*  Initial Data for U,V,P  *
        READ(20,'(A80)')
        READ(20,*) (N,VX0(N),VY0(N),PR0(N),I=1,NP)
 
!*  Initial Data for Tij  *
		READ(20,'(A80)')
        READ(20,*) (N,TXX0(N),TXY0(N),TYY0(N),I=1,NP)

!*  Initial Data for PHI  *
		READ(20,'(A80)')
        READ(20,*) (N,PHI0(N),I=1,NP)
        
       DO 20 I=1,NPSUB
           VX1(I) = VX0(I)
           VY1(I) = VY0(I)
           PR1(I) = PR0(I)
		  TXX1(I) = TXX0(I)
          TXY1(I) = TXY0(I)
		  TYY1(I) = TYY0(I)
          PHI1(I) = PHI0(I)
   20   CONTINUE

      ENDIF

!C  200 FORMAT(I5,3E15.7)

      CLOSE(20)
      RETURN

   90 WRITE(*,130) INAME
  130 FORMAT(/'>>>>>> Error (1): Cannot find ',A/)

      STOP
      END


!*********************************
!*  Input Boundary Conditions    *
!*  Last Modified : NOV 16,2010  *
!*********************************

      SUBROUTINE INPUT4(VMEAN)
     

      use parameters 

 IMPLICIT NONE
      INTEGER I,M,N,N1,N2
      
      DOUBLE PRECISION LL,X1,X2,Y1,Y2,VMEAN,PG

      
      CHARACTER CMNT*25,BCNAME*20



!*     *********************************
!*     *  NBC# : Number of B.C		   *
!*     *  NNBC#: Node Number		   *
!*     *  VBC# : Value				   *
!*     *                               *
!*     *  # | 1 | 2 | 3 | 4 | 5 | 6 |  *
!*     * ---+---+---+---+---+---+---+  *
!*     * var| u | v | p |Txx|Txy|Tyy|  *
!*     *							   *
!*     *********************************

      WRITE(*,100)
      WRITE(*,110)
  100 FORMAT(/'[ Boundary Condition ]')
  110 FORMAT('* Enter B.C. File Name:'/)

      READ(5,*) BCNAME

      OPEN(UNIT=25, ERR=90, FILE=BCNAME, STATUS='OLD')

      READ(25,'(A25)') CMNT
      WRITE(*,120) CMNT
  120 FORMAT('>>>> Comment: ',A)

      READ(25,*) NBC1
      READ(25,*) NBC2
      READ(25,*) NBC3
! 	  READ(25,*) NBC4
!     READ(25,*) NBC5
!     READ(25,*) NBC6
      READ(25,*) NBCN
!  	  READ(25,*) NBCF
      READ(25,*) VMEAN

!*  B.C. DATA FOR U  *
      READ(25,'(A80)')
      READ(25,*) (NNBC1(N),VBC1(N),N=1,NBC1)

!*  B.C. DATA FOR V  *
      READ(25,'(A80)')
      READ(25,*) (NNBC2(N),VBC2(N),N=1,NBC2)

!*  B.C. DATA FOR P  *
      READ(25,'(A80)')
      READ(25,*) (NNBC3(N),VBC3(N),N=1,NBC3)

!*  NATURAL B.C.  *
      READ(25,'(A80)')
      DO 10 I=1,NBCN
        READ(25,*) N,M,PG

        IF (M.EQ.1) THEN
          N1 = NOP(N,2)
          N2 = NOP(N,3)

        ELSE IF (M.EQ.2) THEN
          N1 = NOP(N,1)
          N2 = NOP(N,3)

        ELSE IF (M.EQ.3) THEN
          N1 = NOP(N,1)
          N2 = NOP(N,2)

        ENDIF

        NNBCN(I,1) = N1
        NNBCN(I,2) = N2

        X1 = CORD(N1,1)
        Y1 = CORD(N1,2)
        X2 = CORD(N2,1)
        Y2 = CORD(N2,2)

        LL = SQRT( (X1-X2)**2 + (Y1-Y2)**2 )
        VBCN(I) = PG*LL

   10 CONTINUE

  
	!*  B.C. Data for Txx  *
		
!          READ(25,'(A80)')		
!          READ(25,*) (NNBC4(N),VBC4(N),N=1,NBC4)
		

	!*  B.C. Data for Txy  *
        
!  		READ(25,'(A80)')
!  		READ(25,*) (NNBC5(N),VBC5(N),N=1,NBC5)
		

	!*  B.C. Data for Tyy  *
        
!  		READ(25,'(A80)')
!  		READ(25,*) (NNBC6(N),VBC6(N),N=1,NBC6)
		

	!*  B.C. Data for Phi  *
        
!  		READ(25,'(A80)')
!  		READ(25,*) (NNBCF(N),VBCF(N),N=1,NBCF)
		
      CLOSE(25)
      RETURN

   90 WRITE(*,130) BCNAME
  130 FORMAT(/,'>>>>>> Error (1): Cannot find ',A/)

      STOP
      END

!*************************************
!*  Output Results per interval      *
!*  Output ALL Nodes DATUM           *
!*  Last Modified : OCT 20,2010      *
!*************************************
	subroutine OUTPUT_ALL(FNAME,NESB,NPSB,TIME,NFSB,DTSB)
	use parameters

	IMPLICIT NONE
      INTEGER N,NESB,NPSB
      DOUBLE PRECISION TIME,DTSB
	  INTEGER l1,l2,l3,NFSB
      CHARACTER FNAME*5,ONAME*20

	!Update File Name
		l1=NFSB/100
		l2=(NFSB-100*l1)/10
		l3=NFSB-100*l1-10*l2
		ONAME=FNAME//CHAR(48+l1)//CHAR(48+l2)//CHAR(48+l3)//'.dat' !Making File Name
	
	!Open File
      OPEN(UNIT=50, FILE=ONAME, STATUS='UNKNOWN')

!C  Comment, Time
!*      WRITE(30,*) TIME

!C  Data for U,V,P
        WRITE(50,110)
        WRITE(50,200) (N,TIME,VX1(N),VY1(N),PR1(N),N=1,NP)
		
!*  Data for Tij
        WRITE(50,120)
        WRITE(50,210) (N,TIME,TXX1(N),TXY1(N),TYY1(N),N=1,NP)
		
!*  Data for phi
        WRITE(50,130)
        WRITE(50,220) (N,TIME,PHI1(N),ETA(N),FNSD(N),N=1,NP)


		100 FORMAT(E15.7/)
		110 FORMAT(' ',' NODE    TIME  [s]     Vx  [m/s]      Vy  [m/s]        p  [Pa]')
		120 FORMAT(' ',' NODE    TIME  [s]    Txx  [Pa]      Txy  [Pa]       Tyy  [Pa]')	
		130 FORMAT(' ',' NODE    TIME  [s]    PHI  [1/Pa*s]  ETA  [1/Pa*s]    N1  [Pa]')
		200 FORMAT(' ',I5,4E15.7E2)
		210 FORMAT(' ',I5,4E15.7E2)
		220 FORMAT(' ',I5,4E15.7E2)
		
		CLOSE(50)

		RETURN
		END
		
		
!*************************************
!*  Output Results                   *
!*  For Micro AVS FEM                *
!*  Last Modified : NOV 16,2010      *
!*************************************
   SUBROUTINE OUTPUT_FEM(NPSB,NESB,FEMCOUNT,TIME,TENDSB,DTSB,INTVALSB,FNAMESB)
  
   USE PARAMETERS 

 IMPLICIT NONE
      INTEGER N,NPSB,NESB,FEMCOUNT
	  INTEGER :: W,X=0,Y=1,Z=2
	  
	  INTEGER :: COUNTER=1
	  DOUBLE PRECISION TIME
	  DOUBLE PRECISION TENDSB,DTSB,INTVALSB
      CHARACTER FNAMESB*20

!* WRITE OUTPUT DATA FILE *
		W = TENDSB/(DTSB*INTVALSB)
	IF(COUNTER==1) THEN
		
		WRITE(30,100) 
		WRITE(30,110) W !* W is Integer type. Might wrong type.
		WRITE(30,120) 
		WRITE(30,130) 'step',COUNTER
		WRITE(30,140) NPSB,NESB
		WRITE(30,150) (N,CORD(N,1),CORD(N,2),0.d0,N=1,NPSB)
		WRITE(30,160) (N,0,'tri',NOP(N,1),NOP(N,2),NOP(N,3),N=1,NESB)
		WRITE(30,170) Z,X
		WRITE(30,180) Z,Y,Y
		WRITE(30,190) 'VX,'
		WRITE(30,200) 'VY,'
		WRITE(30,210) (N,VX1(N),VY1(N),N=1,NPSB)
		
		100 FORMAT('# mbm data')
		110 FORMAT(I0) !*('WRITE HERE TOTAL STEP')
		120 FORMAT('data')
		130 FORMAT(A4,I0)
		140 FORMAT(I0,1x,I4)
		150 FORMAT(I0,3E14.7E2)
		160 FORMAT(I0,1X,I1,1X,A3,1x,I0,1X,I0,1X,I0)
		170 FORMAT(I0,I2)
		180 FORMAT(I0,I2,I2)
		190 FORMAT(A3)
		200 FORMAT(A3)
		210 FORMAT(I0,2E15.7E2)
		
		COUNTER=COUNTER+1

		ELSE

			WRITE(30,220) 'step',FEMCOUNT
			WRITE(30,230) Z,X
			WRITE(30,240) Z,Y,Y
			WRITE(30,250) 'VX,'
			WRITE(30,260) 'VY,'
			WRITE(30,270) (N,VX1(N),VY1(N),N=1,NPSB)							
	
			220 FORMAT(A4,I0)
			230 FORMAT(I0,I2)
			240 FORMAT(I0,I2,I2)
			250 FORMAT(A3)
			260 FORMAT(A3)
			270 FORMAT(I0,2E15.7E2)		
			
		END IF
		

      RETURN
      END

!*********************************
!*  Display Information          *
!*  Last Modified : DEC 21,2010  *
!*********************************

      SUBROUTINE INFO(NESB,NPSB,DT,TIME,TEND,ONAME,FNAME,CMNT)
	  USE PARAMETERS
      IMPLICIT NONE

      INTEGER NESB,NPSB
      DOUBLE PRECISION DT,TIME,TEND
      
      CHARACTER ONAME*20,FNAME*20,CMNT*25

      WRITE(*,50)
      WRITE(*,100) NPSB,NESB
      WRITE(*,110) TIME,TEND,DT
      WRITE(*,60)
      
      WRITE(*,200) TY
      WRITE(*,210) MM
	  WRITE(*,220) F0
	  WRITE(*,230) F1
	  WRITE(*,240) LAMBDA
	  WRITE(*,250) K0
	  WRITE(*,260) DENS
      WRITE(*,270) CMNT
      WRITE(*,280) ONAME
	  WRITE(*,290) FNAME
	  WRITE(*,300) 

   50 FORMAT(' '/'*** Condition for Calculation ***'/)
   60 FORMAT(' ','*** Test Fluid ***'/)
  100 FORMAT(' ',I4,'nodes, ',I4,'element'/)
  110 FORMAT(' ','Start Time:',E15.7, ', Final Time:',E15.7/,' ','Incriment :',E15.7/)
  
  200 FORMAT(' ','Dimensional Yield Stress: ',F10.3)
  210 FORMAT(' ','Power-Index: ',F10.3)
  220 FORMAT(' ','FLUIDITY ZERO S-R: ',F10.3)
  230 FORMAT(' ','FLUIDITY HIGH S-R: ',F10.3)
  240 FORMAT(' ','STRUCTURAL RELAXATION TIME: ',F10.3)
  250 FORMAT(' ','KINETIC CONSTANT: ',F10.3)
  260 FORMAT(' ','DENSITY: ',F10.3)
  270 FORMAT(' ','Comment: ',A)
  280 FORMAT(' ','Output File: ',A//)
  290 FORMAT(' ','Result File: ',A//)
  300 FORMAT(' ',' Time    Res')	

      RETURN
      END

!*********************************
!*  Convergence Check            *
!*  Last Modified : NOV 11,1999  *
!*********************************

      SUBROUTINE CONV(RESSUB,NPSUB,VMEANSUB)
    
      use parameters

  IMPLICIT NONE
      INTEGER I,NPSUB
      DOUBLE PRECISION RESSUB,TMP,V0,V1,VMEANSUB

      RESSUB = 0.0
      DO 10 I=1,NPSUB
        V0 = SQRT(VX0(I)**2 + VY0(I)**2)
        V1 = SQRT(VX1(I)**2 + VY1(I)**2)
        TMP = DABS(V1-V0)/VMEANSUB
        RESSUB = MAX(RESSUB,TMP)
   10 CONTINUE

      RETURN
      END


!*********************************
!*  Calculate Matrices           *
!*  Last Modified : AUG 9,1999   *
!*********************************

      SUBROUTINE MKMTR(ELEM)
    
      use parameters 

  IMPLICIT NONE
      INTEGER I,J,K,N,N1,N2,N3
	  INTEGER ELEM
      DOUBLE PRECISION AREA,AREA02,AREA3,AREA12,A2
      DOUBLE PRECISION B(3),C(3),DD(3,3)

      
     
      DOUBLE PRECISION X1,X2,X3,Y1,Y2,Y3      

!*  Initialization  *
      DO 10 I=1,3
        EMM(ELEM,I) = 0.D0
        DO 20 J=1,3
          EHX(ELEM,I,J) = 0.D0
          EHY(ELEM,I,J) = 0.D0
          ESS(ELEM,I,J) = 0.D0
          DO 30 K=1,3
            EKX(ELEM,I,J,K) = 0.D0
            EKY(ELEM,I,J,K) = 0.D0
   30     CONTINUE
   20   CONTINUE
   10 CONTINUE

      DD(1,1) = 2.0D0
      DD(1,2) = 1.0D0
      DD(1,3) = 1.0D0

      DD(2,1) = 1.0D0
      DD(2,2) = 2.0D0
      DD(2,3) = 1.0D0

      DD(3,1) = 1.0D0
      DD(3,2) = 1.0D0
      DD(3,3) = 2.0D0

!*  Calculate Area of Element  *
!*     n1 -- number of 1st node
!*     n2 -- number of 2nd node
!*     n3 -- number of 3rd node
!*     area -- area of element

      N1 = NOP(ELEM,1)
      N2 = NOP(ELEM,2)
      N3 = NOP(ELEM,3)

      X1 = CORD(N1,1)
      Y1 = CORD(N1,2)

      X2 = CORD(N2,1)
      Y2 = CORD(N2,2)

      X3 = CORD(N3,1)
      Y3 = CORD(N3,2)

      AREA02 = X1*(Y2-Y3) + X2*(Y3-Y1) + X3*(Y1-Y2)
      AREA = 0.5D0*AREA02

      IF (AREA.LE.0.D0) THEN
        WRITE(*,*) 'ERROR: Area is less than zero (',ELEM, ')'
        STOP
      END IF

      A2      = 1.0D0 / AREA02
      AREA3  = AREA / 3.0D0
      AREA12 = AREA / 12.0D0

!*  Differential of Shape Function  *
      B(1) = (Y2-Y3) * A2
      B(2) = (Y3-Y1) * A2
      B(3) = (Y1-Y2) * A2
      C(1) = (X3-X2) * A2
      C(2) = (X1-X3) * A2
      C(3) = (X2-X1) * A2

!C
      AR(ELEM) = AREA
      EB(ELEM,1) = B(1)
      EB(ELEM,2) = B(2)
      EB(ELEM,3) = B(3)
      EC(ELEM,1) = C(1)
      EC(ELEM,2) = C(2)
      EC(ELEM,3) = C(3)
!C

!*  Calculate Matrices  *
!*     [Hx], [Hy]

      DO 100 I=1,3
        DO 110 J=1,3
          EHX(ELEM,I,J) = AREA3 * B(J)
          EHY(ELEM,I,J) = AREA3 * C(J)
  110   CONTINUE
  100 CONTINUE

!*     [Kxx], [Kyy]

      DO 200 I=1,3
        DO 210 J=1,3
          DO 220 K=1,3
            EKX(ELEM,I,J,K) = AREA12*DD(I,J)*B(K)
            EKY(ELEM,I,J,K) = AREA12*DD(I,J)*C(K)
  220     CONTINUE
  210   CONTINUE
  200 CONTINUE

!*      [S]

       DO 300 I=1,3
         DO 310 J=1,3
           ESS(ELEM,I,J) = AREA*(B(I)*B(J)+C(I)*C(J))
  310    CONTINUE
  300  CONTINUE

!*       _
!*      [M] (lumped mass matrix)
!*
       DO 400 I=1,3
         EMM(ELEM,I) = AREA3
         N = NOP(ELEM,I)
         LMM(N) = LMM(N) + AREA3
  400  CONTINUE
!*
       RETURN
       END


!*********************************
!*  Step 1                       *
!*  Last Modified : AUG 05,2000  *
!*********************************

      SUBROUTINE STEP1(DT)
   
	  use parameters

   IMPLICIT NONE
      
      INTEGER I,N
      DOUBLE PRECISION C1,DT,DDT


!*  Step 1(a)  *
!  Calculate Velocity  *
      C1 = 0.5D0*DT/DENS
      DDT = 0.5D0*DT


!*  Calculate R.H.S. Vector (U)
      CALL RVVX()

!*  Calculate R.H.S. Vector (V)
      CALL RVVY()

!*  B.C. (u)  *
      DO 100 I=1,NBC1
        N = NNBC1(I)
        VX0(N) = VBC1(I)
  100 CONTINUE

!*  Explicit Method
      DO 10 I=1,NP
        VX1(I)  = VX0(I) + C1*VVX(I)*ILMM(I)
   10 CONTINUE

!*  B.C. (v)  *
      DO 110 I=1,NBC2
        N = NNBC2(I)
        VY0(N) = VBC2(I)
  110 CONTINUE

!*  Explicit Method
      DO 20 I=1,NP
        VY1(I)  = VY0(I) + C1*VVY(I)*ILMM(I)
   20 CONTINUE

!*  Calculate FLUIDITY (Thixotoropic Fluids/MBM model) *

	CALL RHSV_FLUID(DT)
	
	DO I=1,NP
		PHI1(I)=PHI0(I)+DDT*VPHI(I)*ILMM(I)
		END DO
	
!* B.C. (FLUIDITY)
!* 	DO I=1,NBCF
!* 		N=NNBCF(I)
!* 		PHI1(N) = VBCF(I)
!* 		END DO

!*  Calculate Viscosity
	DO I=1,NP
		IF(PHI1(I) .NE. 0) THEN
			ETA(I) = 1.D0/PHI1(I)
		END IF
	END DO

!*  Calculate STRESS (Thixotoropic Fluids/MBM model) *

	CALL RHSV_T(DT)
	DO I=1,NP
		TXX1(I) = TXX0(I)+TXXBing(I)
		TXY1(I) = TXY0(I)+TXYBing(I)
		TYY1(I) = TYY0(I)+TYYBing(I)
			
	END DO

!* COMMENT OUT B.C. APPLY FOR CIRCLE
!* B.C. (Tij) 
!*   	DO I=1,NBC4
!*   	N = NNBC4(I)
!*   	TXX1(N) = VBC4(I)
!*   	END DO

!*   	DO I=1,NBC5
!*   	N = NNBC5(I)
!*   	TXY1(N) = VBC5(I)
!*   	END DO
	
!*   	DO I=1,NBC6
!*   	N = NNBC6(I)
!*   	TYY1(N) = VBC6(I)

!*	END DO

	!*  Calculate First Normal Stress Difference *
	DO I=1,NP
		FNSD(I) = TXX1(I) - TYY1(I)
	END DO	
	

!*  Step 1(b)  *
!*  Calculate Velocity  *
       C1 = DT/DENS
	   DDT = 0.5D0*DT
      
!*  Calculate R.H.S. Vector
      CALL RVVX()

!*  Calculate R.H.S. Vector
      CALL RVVY()

!*  B.C. (u)  *
      DO 300 I=1,NBC1
        N = NNBC1(I)
        VX0(N) = VBC1(I)
  300 CONTINUE

!*  Explicit Method
      DO 30 I=1,NP
        VX1(I) = VX0(I) + C1*VVX(I)*ILMM(I)
   30 CONTINUE

!*  B.C. (v)  *
      DO 310 I=1,NBC2
        N = NNBC2(I)
        VY0(N) = VBC2(I)
  310 CONTINUE

!*  Explicit Method
      DO 40 I=1,NP
        VY1(I) = VY0(I) + C1*VVY(I)*ILMM(I)
   40 CONTINUE

   	CALL RHSV_FLUID(DT)

!*  Calculate Stress (Thixotoropic Fluids/MBM model) *

!*  Calculate Viscosity
!*	DO I=1,NP
!*		PHI1(I)=PHI0(I)+DDT*VPHI(I)*ILMM(I)
!*		END DO

			
!* B.C. (FLUIDITY)
!* 	DO I=1,NBCF
!* 		N=NNBCF(I)
!* 		PHI1(N) = VBCF(I)
!* 		END DO

!*  Calculate Viscosity
	DO I=1,NP
		IF(PHI1(I) .NE. 0) THEN
			ETA(I) = 1.D0/PHI1(I)
		END IF
	END DO

!*  Calculate STRESS (Thixotoropic Fluids/MBM model) *

	CALL RHSV_T(DT)
	DO I=1,NP	
		TXX1(I) = TXX0(I)+TXXBing(I)
		TXY1(I) = TXY0(I)+TXYBing(I)
		TYY1(I) = TYY0(I)+TYYBing(I)
	END DO

!* B.C. (Tij)
!*   	DO I=1,NBC4
!*   	N = NNBC4(I)
!*   	TXX1(N) = VBC4(I)
!*   	END DO

!*   	DO I=1,NBC5
!*   	N = NNBC5(I)
!*   	TXY1(N) = VBC5(I)
!*   	END DO
	
!*   	DO I=1,NBC6
!*   	N = NNBC6(I)
!*   	TYY1(N) = VBC6(I)
!*   	END DO

!*  Calculate First Normal Stress Difference *
	DO I=1,NP
		FNSD(I) = TXX1(I) - TYY1(I)
	END DO		

	RETURN
	END


!************************************
!*  Step 2: Pressure (Poisson Eq.)  *
!*  Last Modified : AUG 19,1999     *
!************************************

      SUBROUTINE STEP2(DT)
  

      use parameters
    IMPLICIT NONE
      DOUBLE PRECISION DT
      INTEGER I


!**********************************************
!*  CG Method Routine                         *
!*     EPS : Condition for Convergence Check  *
!*     LMAX: Maximum Number of Iteration      *
!*     XX  : Solution                         *
!**********************************************

!C     LMAX = 1000
      CALL CGPOIS(DT)

!C  Update Pressure
      DO 20 I=1,NP
        PR1(I) = XX(I)
   20 CONTINUE

      RETURN
      END


!***********************************
!*  Step 3: Velocity Correction    *
!*  Last Modified : DEC 21, 2010   *
!***********************************

      SUBROUTINE STEP3(DT,DENSSB)
     
      use parameters
				
 IMPLICIT NONE

      INTEGER I,J,N,N1,N2 
      DOUBLE PRECISION C1,DT,DENSSB

      
	  C1 = 1.D0/DENSSB

!*  Initialization  *
      DO 10 I=1,NP
        VVX(I) = 0.D0
        VVY(I) = 0.D0
   10 CONTINUE

!*  Calculate Vector  *
      DO 100 N=1,NE
        DO 110 I=1,3
          N1 = NOP(N,I)
          DO 120 J=1,3
            N2 = NOP(N,J)
            VVX(N1) = VVX(N1) + EHX(N,I,J)*PR1(N2)*C1
            VVY(N1) = VVY(N1) + EHY(N,I,J)*PR1(N2)*C1
  120     CONTINUE
  110   CONTINUE
  100 CONTINUE

      DO 300 I=1,NP
        VX1(I) = VX1(I) - DT*VVX(I)*ILMM(I)
        VY1(I) = VY1(I) - DT*VVY(I)*ILMM(I)
  300 CONTINUE

!*  Boundary Conditions  *
      DO 200 I=1,NBC1
        VX1(NNBC1(I)) = VBC1(I)
  200 CONTINUE

      DO 210 I=1,NBC2
        VY1(NNBC2(I)) = VBC2(I)
  210 CONTINUE

      RETURN
      END


!************************************************
!*  Calculate R.H.S. Vector for Eq. Motion (u)  *
!*  Step 1                                      *
!*  Last Modified : NOV 17, 2010                *
!************************************************

      SUBROUTINE RVVX()
 
      use parameters
				
     IMPLICIT NONE      
      INTEGER I,J,K,N,N1,N2,N3
      DOUBLE PRECISION KK,II,FF
      DOUBLE PRECISION DUDX,DVDY,DUYVX
	  


!*  Initialization
      DO 500 I=1,NP
        VVX(I) = 0.D0
  500 CONTINUE

      	  
!*  Calculate R.H.S. Vector
	DO N=1,NE
		DO I=1,3
			N1=NOP(N,I)
			DO J=1,3
				N2=NOP(N,J)
				VVX(N1) = VVX(N1) + EHX(N,I,J)*TXX1(N2)&
								  + EHY(N,I,J)*TXY1(N2)
								  !* ADD SOMETHING HERE ? 

				DO K=1,3
					N3=NOP(N,K)
					VVX(N1) = VVX(N1) - DENS*EKX(N,I,J,K)*VX1(N2)*VX1(N3)&
									  - DENS*EKY(N,I,J,K)*VY1(N2)*VX1(N3)
				END DO
			END DO
		END DO
	END DO

!*  B.C.  *
      DO 100 I=1,NBC1
        VVX(NNBC1(I)) = 0.D0
  100 CONTINUE

      RETURN
      END


!************************************************
!*  Calculate R.H.S. Vector for Eq. Motion (v)  *
!*  Step 1                                      *
!*  Last Modified : NOV 17, 2010                *
!************************************************

      SUBROUTINE RVVY()


      use parameters 
      IMPLICIT NONE
      INTEGER I,J,K,N,N1,N2,N3
      DOUBLE PRECISION KK,II,REI,FF
      DOUBLE PRECISION DUDX,DVDY,DUYVX
	  


!*  Initialization
      DO 500 I=1,NP
        VVY(I) = 0.D0
  500 CONTINUE

!*  Calculate R.H.S. Vector
	DO N=1,NE
		DO I=1,3
			N1=NOP(N,I)
			DO J=1,3
				N2=NOP(N,J)
				VVY(N1) = VVY(N1) + EHX(N,I,J)*TXY1(N2)&
								  + EHY(N,I,J)*TYY1(N2)
								  !* ADD SOMETHING HERE ? 

				DO K=1,3
					N3=NOP(N,K)
					VVY(N1) = VVY(N1) - DENS*EKX(N,I,J,K)*VX1(N2)*VY1(N3)&
									  - DENS*EKY(N,I,J,K)*VY1(N2)*VY1(N3)
				END DO
			END DO
		END DO
	END DO


!*  B.C.  *
      DO 100 I=1,NBC2
        VVY(NNBC2(I)) = 0.D0
  100 CONTINUE

      RETURN
      END

!********************************************************************
!*  Calculate R.H.S. Vector for Fredrickson Model >> Stress tensor  *
!*  Step 1a                                                         *
!*  Last Modified : NOV 13, 2010                                    *
!*  Added		  : NOV 22, 2010                                    *
!********************************************************************

SUBROUTINE RHSV_T(DT)
	USE PARAMETERS
	implicit none

	integer I,J,K,N,N1,N2,N3
	DOUBLE PRECISION C1,C2,C3,GG,DT
	DOUBLE PRECISION AREA,AREA02,DX(3),DY(3)
	DOUBLE PRECISION FF,KK,II
	DOUBLE PRECISION M,MMM,M5,M6
	DOUBLE PRECISION TT,X1,X2,X3,Y1,Y2,Y3
	DOUBLE PRECISION DUDX,DVDY,DUYVX

	!* Initialization *
	DO I=1,NP
		TXXBing(I) = 0.D0
		TXYBing(I) = 0.D0
		TYYBing(I) = 0.D0
		END DO

	!* Calculate Vector *
	DO N=1,NE
	!* Calculate Area, Dx,Dy *
	N1 = NOP(N,1)
	N2 = NOP(N,2)
	N3 = NOP(N,3)
	
		!*  Calculate R.H.S. Vector Second Invariation

	    DUDX = EB(N,1)*VX1(N1)&
			 + EB(N,2)*VX1(N2)&
			 + EB(N,3)*VX1(N3)

        DVDY = EC(N,1)*VY1(N1)&
			 + EC(N,2)*VY1(N2)&
			 + EC(N,3)*VY1(N3)

        DUYVX = EC(N,1)*VX1(N1) + EC(N,2)*VX1(N2)&
			  + EC(N,3)*VX1(N3) + EB(N,1)*VY1(N1)&
			  + EB(N,2)*VY1(N2) + EB(N,3)*VY1(N3)


        II = DSQRT(2.D0*DUDX*DUDX&
				+ 2.D0*DVDY*DVDY&
				+ DUYVX*DUYVX)

		!* Fluidity *
		FF = (PHI1(N1)+PHI1(N2)+PHI1(N3))/3.D0


!*		if II=0 then TAU_{ij} = 0

	IF(II .NE. 0) THEN
		
		!* KK = FF*II**(N0-1.d0) + TY*(1.d0-DEXP(-MM*II))/II
		KK = 1.d0/FF + TY*(1.d0-DEXP(-MM*II))/II

	!* Calculation of Stress *
	
	C1 = 2.D0*KK
	
	DO I=1,3
		N1 = NOP(N,I)
		DO J = 1,3
		N2 = NOP(N,J)
			TXXBing(N1) = TXXBing(N1) + C1*EHX(N,I,J)*VX1(N2)
			TXYBing(N1) = TXYBing(N1) + KK*(EHY(N,I,J)*VX1(N2)+EHX(N,I,J)*VY1(N2))
			TYYBing(N1) = TYYBing(N1) + C1*EHY(N,I,J)*VY1(N2)
		END DO
	END DO
   
   END IF
	
   END DO

      RETURN
      END


!*********************************************************
!*  Calculate R.H.S. Vector for Bautista Model>>Fluidity *
!*  Step 1a                                              *
!*  Last Modified : DEC 14, 2010                         *
!*  Added		  : OCT 25, 2010                         *
!*********************************************************
	subroutine RHSV_FLUID(DT)
	use parameters
	implicit none

	integer I,J,K,N
	integer N1,N2,N3
	DOUBLE PRECISION AREA,AREA02,DX(3),DY(3),B(3),C(3)
	double precision DUDX1,DUDY1,DVDX1,DVDY1
	DOUBLE PRECISION X1,X2,X3,Y1,Y2,Y3
	DOUBLE PRECISION C1,C2
	DOUBLE PRECISION FF
	DOUBLE PRECISION TT,M,MMM
	double precision DT

!* initialization *
	do i=1,NP
	VPHI(I)=0.d0
	end do
	
!* Calculate Vector *
	do N=1,NE

	!* Calculate Area, Dx, Dy
	N1=NOP(N,1)
	N2=NOP(N,2)
	N3=NOP(N,3)

	X1=CORD(N1,1)
	Y1=CORD(N1,2)
	X2=CORD(N2,1)
	Y2=CORD(N2,2)
	X3=CORD(N3,1)
	Y3=CORD(N3,2)

	DX(1) = X3-X2
	DX(2) = X1-X3
	DX(3) = X2-X1
	DY(1) = Y2-Y3
	DY(2) = Y3-Y1
	DY(3) = Y1-Y2

	AREA02=X1*DY(1)+X2*DY(2)+X3*DY(3)
	AREA = 0.5D0*AREA02

	B(1) = DY(1) * AREA02
	B(2) = DY(2) * AREA02
	B(3) = DY(3) * AREA02
	
	C(1) = DX(1) * AREA02
	C(2) = DX(2) * AREA02
	C(3) = DX(3) * AREA02

	FF = (PHI1(N1)+PHI1(N2)+PHI(N3))/3.D0
	C1 = (F0-FF)/LAMBDA
	C2 = K0*(F1-FF)
	
	DO I=1,3
		N1 = NOP(N,I)
!* SUPG
	TT = 0.5D0*DT
	M = VX1(N1)*DY(I)+VY1(N1)*DX(I)
	MMM = TT*M/AREA/12.D0
!* HERE WRONG? 10/11/05
!* ADDED "*PHI1(N1) 10/12/14"
!* OMIT "*PHI1(N1) 10/12/19"
	VPHI(N1) = VPHI(N1)+C1*EMM(N,I)
	VPHI(N1) = VPHI(N1)+C1*EMM(N,I)
	VPHI(N1) = VPHI(N1)+C1*EMM(N,I)

		DO J=1,3
			N2=NOP(N,J)
			DO K=1,3
			N3=NOP(N,K)
			VPHI(N1) = VPHI(N1)+C2*(&
							EKX(N,I,J,K)*VX1(N3)*TXX1(N2)&
						   +EKY(N,I,J,K)*VX1(N3)*TXY1(N2)&
						   +EKX(N,I,J,K)*VY1(N3)*TXY1(N2)&
						   +EKY(N,I,J,K)*VY1(N3)*TYY1(N2)&
						     	   ) 
			END DO
		END DO
	END DO
	END DO

	RETURN
	END

!*********************************
!*  CG Method                    *
!*  Last Modified : AUG 19,1999  *
!*********************************

      SUBROUTINE CGPOIS(DT)


      use parameters
      IMPLICIT NONE

      INTEGER I,J,K,N,N1,N2,ICOUNT,KEND

      DOUBLE PRECISION ALPHA,BETA,DELTA,EPS,EPS1,EPS2
      DOUBLE PRECISION B2,DT,RUR0,RUR1,PAP,RES2

	  DOUBLE PRECISION R(NMAX),P(NMAX)
      DOUBLE PRECISION AP(NMAX),B(NMAX)

    
!C     EPS1 = 1.D-20
      EPS1 = 1.D-04
      EPS2 = EPS1*EPS1
      KEND=1000
      DELTA=-1.D0/DT

!*  Initialization  *
      DO 100 I=1,NP
         R(I)  = 0.D0
         B(I)  = 0.D0
         XX(I) = 0.D0
  100 CONTINUE

!*     r0 = b-Ax0 (x0=0)
      DO 110 N=1,NE
        DO 120 I=1,3
          N1 = NOP(N,I)
          DO 130 J=1,3
            N2 = NOP(N,J)
            B(N1) = B(N1)+DELTA*(EHX(N,I,J)*VX1(N2)&
								+EHY(N,I,J)*VY1(N2))
  130     CONTINUE
  120   CONTINUE
  110 CONTINUE

!  NATURAL B.C.  *
      DO 600 I=1,NBCN
        N1 = NNBCN(I,1)
        N2 = NNBCN(I,2)
        B(N1) = B(N1) + 0.5D0*VBCN(I)
        B(N2) = B(N2) + 0.5D0*VBCN(I)
  600 CONTINUE

!*  B.C. (DP=0)  *
      DO 700 I=1,NBC3
        B(NNBC3(I)) = 0.D0
  700 CONTINUE

!*     p0=r0

      DO 510 I=1,NP
        R(I) = B(I)
        P(I) = B(I)
  510 CONTINUE

      B2 = 0.0
      DO 530 I=1,NP
        B2 = B2 + B(I)**2
  530 CONTINUE
 
      IF(B2 .LT. EPS2) THEN
        EPS = 0.D0
        ICOUNT = 1
        GOTO 999
      END IF

!*     Iterative Loop

      ICOUNT = 0
      DO 10 K=1,KEND
        ICOUNT = ICOUNT + 1

        DO 300 I=1,NP
          AP(I) = 0.D0
  300   CONTINUE

!*     AP : Apk      P: pk
!*     PAP: (pk,Apk) 

        DO 310 N=1,NE
          DO 320 I=1,3
            N1 = NOP(N,I)
            DO 330 J=1,3
              N2 = NOP(N,J)
              AP(N1) = AP(N1) + ESS(N,I,J)*P(N2)
  330       CONTINUE
  320     CONTINUE
  310   CONTINUE

!*  B.C. (DP=0)  *
      DO 800 I=1,NBC3
        AP(NNBC3(I)) = 0.D0
  800 CONTINUE


        RUR0 = 0.D0
        DO 400 I=1,NP
          RUR0 = RUR0 + P(I)*R(I)
  400   CONTINUE

        PAP=0.D0
        DO 540 I=1,NP
          PAP = PAP + P(I)*AP(I)
  540   CONTINUE

        IF (PAP .EQ. 0.0) THEN
          ALPHA = 0.D0
        ELSE
          ALPHA = RUR0/PAP
        END IF

        DO 550 I=1,NP
           XX(I) = XX(I) + ALPHA*P(I)
           R(I)  = R(I)  - ALPHA*AP(I)
  550   CONTINUE

!*  Boundary Condition  *
!*        DO 410 I=1,NBC3
!*           XX(NNBC3(I)) = 0.D0
!*            R(NNBC3(I)) = 0.D0
!*  410   CONTINUE

        RUR1=0.D0
        RES2=0.D0
        DO 560 I=1,NP
          RUR1 = RUR1 + R(I)*AP(I)
          RES2 = RES2 + R(I)**2
  560   CONTINUE

        EPS = RES2/B2
        IF(EPS .LE. EPS2)  THEN
          GO TO 999
        END IF

!*  NOT CONVERGED  *
        BETA = -RUR1/PAP
        DO 570 I=1,NP

!*     p(k+1) = r(k+1)+beta*p(k)

          P(I) = R(I) + BETA*P(I)

  570   CONTINUE

   10 CONTINUE

  999 KEND = ICOUNT
      RETURN
      END


!**********************************
!*  Update Results                *
!*  Last Modified : OCT 10, 2010  *
!**********************************

      SUBROUTINE UPDATE(NPSUB)

      use parameters

      IMPLICIT NONE
      INTEGER I,NPSUB 


      DO 10 I=1,NPSUB
        VX0(I) = VX1(I)
        VY0(I) = VY1(I)
        PR0(I) = PR1(I)
   10 CONTINUE

      RETURN
      END

!*  Added 10/10/20

!**********************************
!*  Update Results per step       *
!*  Last Modified : OCT 10, 2010  *
!**********************************

      SUBROUTINE UPDATE_STEP(NESB,NPSB)



      use parameters

      IMPLICIT NONE
      INTEGER I,NESB,NPSB  


      DO  I=1,NPSB
        VX0(I) = VX1(I)
        VY0(I) = VY1(I)
        PR0(I) = PR1(I)
   END DO

      RETURN
      END


!********************************************
!*   HISTORY                                *
!*                                          *
!*   November 12, 2010: Modifying was begun *
!*   November 30, 2010: ver.1.00 Complete   *
!*   December 02, 2010: CHECKING HAS BEGAN  * 
!********************************************