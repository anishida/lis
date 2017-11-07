! (This is a demo program for the alternative workflow PSD, or
! Preconditioner and Solver Decoupled, contributed by Neil Hodge.
! To execute this program, it must be built with an MPI library,
! and the Fortran 90 interface and the SA-AMG preconditioner must 
! be enabled.)
    
PROGRAM lis_driver

    IMPLICIT NONE

! per the Intel fortran user forums, including a *.h file
! via the preprocessor mechanism works properly, versus using
! the fortran include mechanism

#include "lisf.h"

    INTERFACE
        SUBROUTINE plotxy(plotdata,filename,datalength, &
                          grid,limits,title1,xlabel,ylabel,series_labels, &
                          legend,lines,points,colors,widths,output_animate)
            REAL(kind(0.d0)),INTENT(IN) :: plotdata(:,:)
            CHARACTER(LEN=*),INTENT(IN) :: filename
            INTEGER,INTENT(IN),OPTIONAL :: datalength(:)
            INTEGER,INTENT(IN),OPTIONAL :: grid
            REAL(kind(0.d0)),INTENT(IN),OPTIONAL :: limits(:,:)
            CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: title1
            CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: xlabel,ylabel
            CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: series_labels(:)
            INTEGER,INTENT(IN),OPTIONAL :: legend
            INTEGER,INTENT(IN),OPTIONAL :: lines(:)
            INTEGER,INTENT(IN),OPTIONAL :: points(:)
            INTEGER,INTENT(IN),OPTIONAL :: colors(:)
            INTEGER,INTENT(IN),OPTIONAL :: widths(:)
            INTEGER,INTENT(IN),OPTIONAL :: output_animate
        END SUBROUTINE
    END INTERFACE

    !======================================================================
    ! parameters
    !======================================================================
    INTEGER,PARAMETER :: singI=kind(0)
    INTEGER,PARAMETER :: fullR=kind(0.d0)

    REAL(fullR),PARAMETER :: ZERO      = 0.0_fullR 
    REAL(fullR),PARAMETER :: ROOT_EPS  = 0.0000000596046448_fullR ! 2^(-24) = 5.96046448E-8
    REAL(fullR),PARAMETER :: ONE       = 1.0_fullR 
    REAL(fullR),PARAMETER :: TWO       = 2.0_fullR 
    REAL(fullR),PARAMETER :: THREE     = 3.0_fullR 
    REAL(fullR),PARAMETER :: FOUR      = 4.0_fullR 
    REAL(fullR),PARAMETER :: FIVE      = 5.0_fullR 
    REAL(fullR),PARAMETER :: SIX       = 6.0_fullR
    REAL(fullR),PARAMETER :: SEVEN     = 7.0_fullR
    REAL(fullR),PARAMETER :: EIGHT     = 8.0_fullR 
    REAL(fullR),PARAMETER :: NINE      = 9.0_fullR 
    REAL(fullR),PARAMETER :: TEN       = 10.0_fullR 
    REAL(fullR),PARAMETER :: HALF      = 0.500000000000000_fullR 
    REAL(fullR),PARAMETER :: THIRD     = 0.333333333333333_fullR 
    REAL(fullR),PARAMETER :: FOURTH    = 0.250000000000000_fullR 
    REAL(fullR),PARAMETER :: TENTH     = 0.100000000000000_fullR 
    REAL(fullR),PARAMETER :: PI        = 3.141592653589793_fullR

    !======================================================================
    ! regular old variables
    !======================================================================
    INTEGER(singI) :: ierr,MPIsize,rank
    LIS_INTEGER :: i,j,gn,n,is,ie,iter
    LIS_INTEGER,POINTER :: proc_n(:)
    LIS_VECTOR :: b,x
    LIS_MATRIX :: A
    LIS_PRECON :: precon
    LIS_SOLVER :: solver
    CHARACTER(LEN=1024) :: options
    REAL(fullR),POINTER :: RHS(:),LHS(:,:)
    LOGICAL,POINTER :: LHSassem(:,:)
    INTEGER(singI),POINTER :: HomeProc(:)
    LOGICAL :: LinSysDefFlag
    REAL(fullR) :: ResNorm0,ResNorm,AbsResTol,RelResTol
    REAL(fullR) :: IncNorm0,IncNorm,AbsIncTol,RelIncTol
    INTEGER(singI) :: IterCount,MaxIterCount
    INTEGER(singI) :: nel,nnp,inp,DOF1,DOF2,iProc
    REAL(fullR) :: h
    REAL(fullR),POINTER :: xhat(:),That(:),dT(:),plotdata(:,:)
    INTEGER(singI) :: PlotCount
    INTEGER(singI),POINTER :: lps(:),widths(:)
    CHARACTER(LEN=1024) :: title,PlotFilename
    INTEGER(singI) :: NonlinLHSUpdateFreq,NumberLHSUpdates

    proc_n=>NULL()
    LHS=>NULL()
    RHS=>NULL()
    LHSassem=>NULL()
    HomeProc=>NULL()
    xhat=>NULL()
    That=>NULL()
    dT=>NULL()
    plotdata=>NULL()
    lps=>NULL()
    widths=>NULL()


    !======================================================================
    ! proceed
    !======================================================================
    CALL MPI_Init(ierr)
!    CALL ErrorCheck("MPI_Init",ierr)

    CALL MPI_Comm_size(MPI_COMM_WORLD,MPIsize,ierr)
!    CALL ErrorCheck("MPI_Comm_rank",ierr)
    CALL MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
!    CALL ErrorCheck("MPI_Comm_rank",ierr)

    IF ( MPIsize.LT.4 ) THEN
       IF ( rank.EQ.0 ) THEN
          WRITE(UNIT=*,FMT='(a)') "number of processes must be 4 or more"
          CALL lis_finalize(ierr)
       ENDIF
       STOP
    ENDIF

    ! do it!
    CALL lis_initialize(ierr)
!    CALL ErrorCheck("lis_initialize",ierr)

    ! mesh definition
    nel=100
    nnp=nel+1

    ! allocate data
    ALLOCATE(xhat(nnp))
    ALLOCATE(That(nnp))
    ALLOCATE(RHS(nnp))
    ALLOCATE(LHS(nnp,nnp))
    ALLOCATE(LHSassem(nnp,nnp))
    IF (rank==0) THEN
        ALLOCATE(plotdata(nnp,2))
        ALLOCATE(lps(1))
        lps=1
        ALLOCATE(widths(1))
        widths=3
    END IF

    ! parallelism data
    ALLOCATE(proc_n(0:MPIsize-1))
    proc_n=0
    !======================================
    ! DOFs on all processes
    !======================================
!    proc_n(0)=INT(REAL(nnp,fullR)/TEN,singI)
!    proc_n(1)=2*proc_n(0)
!    proc_n(2)=3*proc_n(0)
!    proc_n(3)=nnp-SUM(proc_n(0:2))
    !======================================
    ! alternate cases, testing out various versions of n_local=0
    !======================================
!    proc_n(0)=0
!    proc_n(1)=INT(nnp/4,singI)
!    proc_n(2)=proc_n(1)
!    proc_n(3)=nnp-SUM(proc_n(0:2))
    !======================================
    proc_n(0)=0
    proc_n(1)=INT(REAL(nnp,fullR)/THREE,singI)
    proc_n(2)=0
    proc_n(3)=nnp-SUM(proc_n(0:2))
    ! output DOFs per process
    IF (rank==0) THEN
        DO i=0,MPIsize-1
            WRITE(UNIT=*,FMT='(a,i0,a,i0)') "proc_n(",i,")=",proc_n(i)
        END DO
        WRITE(UNIT=*,FMT='(a)') ""
    END IF

    ! initial condition
    h=ONE/REAL(nel,fullR)
    DO inp=1,nnp
        xhat(inp)=REAL((inp-1),fullR)*h
        That(inp)=202.0+202.0*COS(FOUR*THIRD*SEVEN*xhat(inp))
        IF (rank==0) plotdata(inp,1)=xhat(inp)
    END DO
    ! Dirichlet BC
    That(1)=101.0
    That(nnp)=808.0
    ! initial state plot
    IF (rank==0) THEN
        DO inp=1,nnp
            plotdata(inp,2)=That(inp)
        END DO
        WRITE(UNIT=title,FMT='(a,i0)') "nonlinear iteration ",0
        PlotCount=1
        WRITE(UNIT=PlotFilename,FMT='(a,i0)') "blah",PlotCount
        CALL plotxy(plotdata=plotdata,filename=TRIM(PlotFilename), &
                    lines=lps,points=lps,widths=widths, &
                    title1=TRIM(title),output_animate=1)
        ! write initial state twice to deal with ffmpeg bug . . .
        PlotCount=PlotCount+1
        WRITE(UNIT=PlotFilename,FMT='(a,i0)') "blah",PlotCount
        CALL plotxy(plotdata=plotdata,filename=TRIM(PlotFilename), &
                    lines=lps,points=lps,widths=widths, &
                    title1=TRIM(title),output_animate=1)
    END IF

    ! convergence and iteration-related criteria
    AbsResTol=1.0d-06
    RelResTol=1.0d-06
    AbsIncTol=1.0d-05
    RelIncTol=1.0d-05
    MaxIterCount=50
    NonlinLHSUpdateFreq=3
    NumberLHSUpdates=0

    ! update linear system
    ! do this outside the nonlinear loop, just for the first iteration
    CALL UpdateLinearSystemThermal(nel,xhat,That,RHS,LHS,LHSassem)
    IncNorm0=ZERO
    ResNorm0=ZERO
    gn=SUM(proc_n)
    DO i=1,gn
        ResNorm0=ResNorm0+RHS(i)**2
    END DO
    ResNorm0=SQRT(ResNorm0)
    IF (rank==0) THEN
        WRITE(UNIT=*,FMT='(a,es22.15)') "initial residual norm:",ResNorm0
        WRITE(UNIT=*,FMT='(a,es22.15)') "absolute residual tolerance:",AbsResTol
        WRITE(UNIT=*,FMT='(a,es22.15)') "relative residual tolerance:",RelResTol
        WRITE(UNIT=*,FMT='(a,es22.15)') "absolute increment tolerance:",AbsIncTol
        WRITE(UNIT=*,FMT='(a,es22.15)') "relative increment tolerance:",RelIncTol
        WRITE(UNIT=*,FMT='(a,i0)') "MaxIterCount: ",MaxIterCount
        WRITE(UNIT=*,FMT='(a,i0)') "NonlinLHSUpdateFreq: ",NonlinLHSUpdateFreq
        WRITE(UNIT=*,FMT='(a)') ""
    END IF
    ResNorm=ResNorm0

    ! HomeProc indicates the sender in Bcast below
    ALLOCATE(HomeProc(gn))
    DOF2=0
    DO iProc=0,(MPIsize-1)
        DOF1=DOF2+1
        DOF2=SUM(proc_n(0:iProc))
        HomeProc(DOF1:DOF2)=iProc
    END DO
    ALLOCATE(dT(gn))

    ! calculate solution
    LinSysDefFlag=.FALSE.
    IterCount=0
    DO
        IterCount=IterCount+1

        CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
        CALL pauser(HALF)
        IF (rank==0) THEN
            WRITE(UNIT=*,FMT='(a)') ""
            WRITE(UNIT=*,FMT='(a)') "======================================================================"
            WRITE(UNIT=*,FMT='(a,i0)') "Start of nonlinear iteration ",IterCount
            WRITE(UNIT=*,FMT='(a)') "======================================================================"
        END IF

        ! define once, and assemble
        IF (.NOT.LinSysDefFlag) THEN
            CALL lis_matrix_create(MPI_COMM_WORLD,A,ierr)
!            CALL ErrorCheck("lis_matrix_set_size",ierr)
            CALL lis_matrix_set_size(A,proc_n(rank),0,ierr)
!            CALL ErrorCheck("lis_matrix_set_size",ierr)
            CALL lis_matrix_get_size(A,n,gn,ierr)
!            CALL ErrorCheck("lis_matrix_get_size",ierr)
            CALL lis_matrix_get_range(A,is,ie,ierr)
!            CALL ErrorCheck("lis_matrix_get_range",ierr)
        END IF

        IF (IterCount==1) THEN
            DO i=is,ie-1
                DO j=1,gn
                    IF (.NOT.LHSassem(i,j)) CYCLE
                    CALL lis_matrix_set_value(LIS_INS_VALUE,i,j,LHS(i,j),A,ierr)
!                    CALL ErrorCheck("lis_matrix_set_value",ierr)
                END DO
            END DO
        ELSE
            IF (NonlinLHSUpdateFreq>0 .AND. MOD(IterCount,NonlinLHSUpdateFreq)==0) THEN
                NumberLHSUpdates=NumberLHSUpdates+1
                DO i=is,ie-1
                    DO j=1,gn
                        IF (.NOT.LHSassem(i,j)) CYCLE
                        CALL lis_matrix_psd_set_value(LIS_INS_VALUE,i,j,LHS(i,j),A,ierr)
!                        CALL ErrorCheck("lis_matrix_psd_set_value",ierr)
                    END DO
                END DO
                CALL lis_matrix_psd_reset_scale(A,ierr)
!                CALL ErrorCheck("lis_matrix_psd_reset_scale",ierr)
            END IF
        END IF

        IF (.NOT.LinSysDefFlag) THEN
            CALL lis_matrix_set_type(A,LIS_MATRIX_CSR,ierr)
!            CALL ErrorCheck("lis_matrix_set_type",ierr)
            CALL lis_matrix_assemble(A,ierr)
!            CALL ErrorCheck("lis_matrix_assemble",ierr)
            CALL lis_vector_duplicate(A,b,ierr)
!            CALL ErrorCheck("lis_vector_duplicate",ierr)
            CALL lis_vector_duplicate(A,x,ierr)
!            CALL ErrorCheck("lis_vector_duplicate",ierr)
        END IF

        ! update the RHS _every_ nonlinear iteration
        DO i=is,ie-1
            CALL lis_vector_set_value(LIS_INS_VALUE,i,RHS(i),b,ierr)
!            CALL ErrorCheck("lis_vector_set_value",ierr)
        END DO
        CALL lis_vector_psd_reset_scale(b,ierr)
!        CALL ErrorCheck("lis_vector_psd_reset_scale",ierr)

        IF (.NOT.LinSysDefFlag) THEN
            CALL lis_solver_create(solver,ierr)
!            CALL ErrorCheck("lis_solver_create",ierr)
!            WRITE(UNIT=options,FMT='(a)') "-p none -i gmres -print out -scale none"
!            WRITE(UNIT=options,FMT='(a)') "-p none -i gmres -print out -scale jacobi"
!            WRITE(UNIT=options,FMT='(a)') "-p ilu -i gmres -print out -scale none"
            WRITE(UNIT=options,FMT='(a)') "-p saamg -saamg_unsym true -i gmres -print out -scale none"
            WRITE(UNIT=options,FMT='(a,a,es22.15)') TRIM(options)," -tol ", &
                 & MAX(MIN(AbsResTol/ResNorm,ONE),MIN(RelResTol*ResNorm0/ResNorm,ONE))/TEN
            IF (rank==0) WRITE(UNIT=*,FMT='(a,a,a)') "solver options: """,TRIM(options),""""
            CALL lis_solver_set_option(TRIM(options),solver,ierr)
!            CALL ErrorCheck("lis_solver_set_option",ierr)
            ! be sure to create and set options for solver before
            ! calling precon_create . . .
            ! A has to be set in "solver" _before_ precon_create is called . . .
            CALL lis_solver_set_matrix(A,solver,ierr)
!            CALL ErrorCheck("lis_solver_set_matrix",ierr)
            CALL lis_precon_psd_create(solver,precon,ierr)
!            CALL ErrorCheck("lis_precon_psd_create",ierr)
            CALL lis_precon_psd_update(solver,precon,ierr)
!            CALL ErrorCheck("lis_precon_psd_update",ierr)
            LinSysDefFlag=.TRUE.
        ELSE
            WRITE(UNIT=options,FMT='(a,es22.15)') "-tol ",MAX(MIN(AbsResTol/ResNorm,ONE),MIN(RelResTol*ResNorm0/ResNorm,ONE))/TEN
            IF (rank==0) WRITE(UNIT=*,FMT='(a,a,a)') "solver options: """,TRIM(options),""""
            CALL lis_solver_set_option(TRIM(options),solver,ierr)
!            CALL ErrorCheck("lis_solver_set_option",ierr)
        END IF
 
        IF (NonlinLHSUpdateFreq>0 .AND. MOD(IterCount,NonlinLHSUpdateFreq)==0) THEN
            IF (rank==0) WRITE(UNIT=*,FMT='(a)') "updating preconditioner now . . ."
            CALL lis_precon_psd_update(solver,precon,ierr)
!            CALL ErrorCheck("lis_precon_psd_update",ierr)
        END IF
        CALL lis_solve_kernel(A,b,x,solver,precon,ierr)
!        CALL ErrorCheck("lis_solve_kernel",ierr)
        CALL lis_solver_get_iter(solver,iter,ierr)
!        CALL ErrorCheck("lis_solver_get_iter",ierr)
        IF (rank==0) WRITE(UNIT=*,FMT='(a,i0)') "number of iterations = ",iter
        
        ResNorm=ZERO
        DO i=is,ie-1
            CALL lis_vector_get_value(x,i,dT(i),ierr)
!            CALL ErrorCheck("lis_vector_get_value",ierr)
            dT(i)=-dT(i)
        END DO
        CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
        CALL pauser(HALF)
        IF (rank==0) WRITE(UNIT=*,FMT='(a)') ""

        ! make sure increment is consistent everywhere, since
        ! (1) system is nonlinear in the independent variables, and
        ! (2) to ensure calculation of convergence criterion is correct
        ! use 0-based MPI rank values
        IF (rank==0) WRITE(UNIT=*,FMT='(a)') "calling MPI_Bcast now . . ."
        DO i=1,gn
            CALL MPI_Bcast(dT(i),1,MPI_DOUBLE_PRECISION,HomeProc(i),MPI_COMM_WORLD,ierr)
        END DO

        ! update solution vector
        That=That+dT

        ! plot solution
        IF (rank==0) THEN
            DO inp=1,nnp
                plotdata(inp,2)=That(inp)
            END DO
            WRITE(UNIT=title,FMT='(a,i0)') "nonlinear iteration ",IterCount
            PlotCount=PlotCount+1
            WRITE(UNIT=PlotFilename,FMT='(a,i0)') "blah",PlotCount
            CALL plotxy(plotdata=plotdata,filename=TRIM(PlotFilename), &
                        lines=lps,points=lps,widths=widths, &
                        title1=TRIM(title),output_animate=1)
        END IF

        ! update linear system
        CALL UpdateLinearSystemThermal(nel,xhat,That,RHS,LHS,LHSassem)

        ! check convergence
        ResNorm=ZERO
        IncNorm=ZERO
        DO i=1,gn
            ResNorm=ResNorm+RHS(i)**2
            IncNorm=IncNorm+dT(i)**2
        END DO
        ResNorm=SQRT(ResNorm)
        IncNorm=SQRT(IncNorm)
        IF (IncNorm0==ZERO) IncNorm0=IncNorm
        IF (rank==0) WRITE(UNIT=*,FMT='(a)') ""
        IF (rank==0) WRITE(UNIT=*,FMT='(a,es22.15)') "absolute residual 2-norm=",ResNorm
        IF (rank==0) WRITE(UNIT=*,FMT='(a,es22.15)') "relative residual 2-norm=",ResNorm/ResNorm0
        IF (rank==0) WRITE(UNIT=*,FMT='(a,es22.15)') "absolute increment 2-norm=",IncNorm
        IF (rank==0) WRITE(UNIT=*,FMT='(a,es22.15)') "relative increment 2-norm=",IncNorm/IncNorm0
        IF (rank==0) WRITE(UNIT=*,FMT='(a,i0)') "number of LHS updates: ",NumberLHSUpdates
        IF (ResNorm<MAX(AbsResTol,RelResTol*ResNorm0) .AND. IncNorm<MAX(AbsIncTol,RelIncTol*IncNorm0)) THEN
            IF (rank==0) WRITE(UNIT=*,FMT='(a)') ""
            IF (rank==0) WRITE(UNIT=*,FMT='(a)') "nonlinear termination criterion ""convergence"" satisfied, stopping"
            EXIT
        END IF

        IF (IterCount==MaxIterCount) THEN
            IF (rank==0) WRITE(UNIT=*,FMT='(a)') ""
            IF (rank==0) WRITE(UNIT=*,FMT='(a,i0,a)') &
                 & "nonlinear termination criterion num_iter=", &
                 & MaxIterCount," satisfied, stopping"
            EXIT
        END IF
    END DO

    IF (rank==0) WRITE(UNIT=*,FMT='(a)') ""

    CALL lis_precon_destroy(precon,ierr)
    CALL lis_solver_destroy(solver,ierr)
    CALL lis_vector_destroy(x,ierr)
    CALL lis_vector_destroy(b,ierr)
    CALL lis_matrix_destroy(A,ierr)

    DEALLOCATE(dT)
    DEALLOCATE(HomeProc)

    IF (rank==0) THEN
        DEALLOCATE(widths)
        DEALLOCATE(lps)
        DEALLOCATE(plotdata)
    END IF
    DEALLOCATE(LHSassem)
    DEALLOCATE(LHS)
    DEALLOCATE(RHS)
    DEALLOCATE(That)
    DEALLOCATE(xhat)

    DEALLOCATE(proc_n)

    CALL lis_finalize(ierr)
!    CALL ErrorCheck("lis_finalize",ierr)
    
    CALL MPI_Finalize(ierr)
!    CALL ErrorCheck("MPI_Finalize",ierr)

END PROGRAM lis_driver


!==========================================================================
SUBROUTINE pauser(pause_time)

    REAL(kind(0.d0)),INTENT(IN) :: pause_time

    INTEGER,PARAMETER :: singI=kind(0)
    INTEGER,PARAMETER :: fullR=kind(0.d0)
    INTEGER(singI) :: time_data(8)
    REAL(fullR) :: start_time,end_time,dt

    CALL DATE_AND_TIME(values=time_data)
    start_time=time_data(5)*3600 + time_data(6)*60 + time_data(7)
    DO
        CALL DATE_AND_TIME(values=time_data)
        end_time=time_data(5)*3600 + time_data(6)*60 + time_data(7)
        dt=end_time-start_time
        IF (dt>pause_time) EXIT
    END DO

END SUBROUTINE pauser


!==========================================================================
!SUBROUTINE ErrorCheck(FuncName,ierr)

!    CHARACTER(LEN=*) :: FuncName
!    INTEGER(kind(0)),INTENT(IN) :: ierr

!    IF (ierr/=0) THEN
!        WRITE(UNIT=*,FMT='(a,a,a,i0,a)') "ERROR: ",TRIM(FuncName)," returned error ",ierr,", quitting . . ."
!        STOP
!    END IF

!END SUBROUTINE ErrorCheck


!==========================================================================
SUBROUTINE UpdateLinearSystemThermal(nel,xhat,That,RHS,LHS,LHSassem)

    IMPLICIT NONE

    INTEGER(kind(0)),INTENT(IN) :: nel
    REAL(kind(0.d0)),INTENT(IN) :: xhat(nel+1),That(nel+1)
    REAL(kind(0.d0)),INTENT(OUT) :: RHS(nel+1),LHS(nel+1,nel+1)
    LOGICAL,INTENT(OUT) :: LHSassem(nel+1,nel+1)

    INTEGER,PARAMETER :: singI=kind(0)
    INTEGER,PARAMETER :: fullR=kind(0.d0)

    REAL(fullR),PARAMETER :: ZERO      = 0.0_fullR 
    REAL(fullR),PARAMETER :: ROOT_EPS  = 0.0000000596046448_fullR ! 2^(-24) = 5.96046448E-8
    REAL(fullR),PARAMETER :: ONE       = 1.0_fullR 
    REAL(fullR),PARAMETER :: TWO       = 2.0_fullR 
    REAL(fullR),PARAMETER :: THREE     = 3.0_fullR 
    REAL(fullR),PARAMETER :: FOUR      = 4.0_fullR 
    REAL(fullR),PARAMETER :: FIVE      = 5.0_fullR 
    REAL(fullR),PARAMETER :: SIX       = 6.0_fullR
    REAL(fullR),PARAMETER :: SEVEN     = 7.0_fullR
    REAL(fullR),PARAMETER :: EIGHT     = 8.0_fullR 
    REAL(fullR),PARAMETER :: NINE      = 9.0_fullR 
    REAL(fullR),PARAMETER :: TEN       = 10.0_fullR 
    REAL(fullR),PARAMETER :: HALF      = 0.500000000000000_fullR 
    REAL(fullR),PARAMETER :: THIRD     = 0.333333333333333_fullR 
    REAL(fullR),PARAMETER :: FOURTH    = 0.250000000000000_fullR 
    REAL(fullR),PARAMETER :: TENTH     = 0.100000000000000_fullR 
    REAL(fullR),PARAMETER :: PI        = 3.141592653589793_fullR

    LOGICAL :: BuildAssemFlag
    REAL(fullR) :: xi(2),N(2,2)
    INTEGER(singI) :: nnp,inp,jnp,iel,iip,DOF1,DOF2
    REAL(fullR) :: xhat_elem(2),xbar_elem,That_elem(2)
    REAL(fullR) :: dN_dxi(2,2),dN_dx(2)
    REAL(fulLR) :: Tip,k,dk_dT,dT_dx
    REAL(fullR) :: Relem(2),Kelem(2,2)


    BuildAssemFlag=.FALSE.
    IF (.NOT.LHSassem(1,1)) BuildAssemFlag=.TRUE.

    ! initialization
    nnp=nel+1
    DO inp=1,nnp
        RHS(inp)=ZERO
        DO jnp=1,nnp
            IF (BuildAssemFlag) LHSassem(inp,jnp)=.FALSE.
            LHS(inp,jnp)=ZERO
        END DO
    END DO

    xi(1)=-SQRT(THIRD)
    xi(2)=+SQRT(THIRD)

    CALL CalcN(xi,N)
    CALL CalcdNdxi(dN_dxi)

    DO iel=1,nel
        ! local element data
        DOF1=iel
        DOF2=iel+1

        xhat_elem(1)=xhat(DOF1)
        xhat_elem(2)=xhat(DOF2)
        xbar_elem=(xhat_elem(1)+xhat_elem(2))/TWO
        That_elem(1)=That(DOF1)
        That_elem(2)=That(DOF2)
        Relem=ZERO
        Kelem=ZERO
        DO iip=1,2
            dN_dx(1)=ONE/(dN_dxi(1,iip)*xhat_elem(1)+dN_dxi(2,iip)*xhat_elem(2))*dN_dxi(1,iip)
            dN_dx(2)=ONE/(dN_dxi(1,iip)*xhat_elem(1)+dN_dxi(2,iip)*xhat_elem(2))*dN_dxi(2,iip)
            ! Temperature at the current integration point
            Tip=N(1,iip)*That_elem(1)+N(2,iip)*That_elem(2)
            ! Calculate temperature-dependent material properties
            CALL CalcCond(Tip,xbar_elem,k,dk_dT)
            dT_dx=dN_dx(1)*That_elem(1)+dN_dx(2)*That_elem(2)

            Relem=Relem+dN_dx*k*dT_dx

            Kelem(1,1)=Kelem(1,1)+dN_dx(1)*k*dN_dx(1)+dN_dx(1)*dT_dx*dk_dT*N(1,iip)
            Kelem(1,2)=Kelem(1,2)+dN_dx(1)*k*dN_dx(2)+dN_dx(1)*dT_dx*dk_dT*N(2,iip)
            Kelem(2,1)=Kelem(2,1)+dN_dx(2)*k*dN_dx(1)+dN_dx(2)*dT_dx*dk_dT*N(1,iip)
            Kelem(2,2)=Kelem(2,2)+dN_dx(2)*k*dN_dx(2)+dN_dx(2)*dT_dx*dk_dT*N(2,iip)
        END DO

        ! assemble
        RHS(DOF1)=RHS(DOF1)+Relem(1)
        RHS(DOF2)=RHS(DOF2)+Relem(2)

        IF (BuildAssemFlag) THEN
            LHSassem(DOF1,DOF1)=.TRUE.
            LHSassem(DOF1,DOF2)=.TRUE.
            LHSassem(DOF2,DOF1)=.TRUE.
            LHSassem(DOF2,DOF2)=.TRUE.
        END IF
        LHS(DOF1,DOF1)=LHS(DOF1,DOF1)+Kelem(1,1)
        LHS(DOF1,DOF2)=LHS(DOF1,DOF2)+Kelem(1,2)
        LHS(DOF2,DOF1)=LHS(DOF2,DOF1)+Kelem(2,1)
        LHS(DOF2,DOF2)=LHS(DOF2,DOF2)+Kelem(2,2)
    END DO

    ! apply dirichlet BC
    RHS(1)=ZERO
    RHS(nnp)=ZERO
    DO inp=1,nnp
        LHS(1,inp)=ZERO
        LHS(inp,1)=ZERO
        LHS(nnp,inp)=ZERO
        LHS(inp,nnp)=ZERO
    END DO
    LHS(1,1)=ONE
    LHS(nnp,nnp)=ONE

END SUBROUTINE UpdateLinearSystemThermal


!==========================================================================
SUBROUTINE CalcN(xi,N)

    IMPLICIT NONE

    REAL(kind(0.d0)),INTENT(IN) :: xi(2)
    REAL(kind(0.d0)),INTENT(OUT) :: N(2,2)

    INTEGER,PARAMETER :: singI=kind(0)
    INTEGER,PARAMETER :: fullR=kind(0.d0)

    REAL(fullR),PARAMETER :: ZERO      = 0.0_fullR 
    REAL(fullR),PARAMETER :: ROOT_EPS  = 0.0000000596046448_fullR ! 2^(-24) = 5.96046448E-8
    REAL(fullR),PARAMETER :: ONE       = 1.0_fullR 
    REAL(fullR),PARAMETER :: TWO       = 2.0_fullR 
    REAL(fullR),PARAMETER :: THREE     = 3.0_fullR 
    REAL(fullR),PARAMETER :: FOUR      = 4.0_fullR 
    REAL(fullR),PARAMETER :: FIVE      = 5.0_fullR 
    REAL(fullR),PARAMETER :: SIX       = 6.0_fullR
    REAL(fullR),PARAMETER :: SEVEN     = 7.0_fullR
    REAL(fullR),PARAMETER :: EIGHT     = 8.0_fullR 
    REAL(fullR),PARAMETER :: NINE      = 9.0_fullR 
    REAL(fullR),PARAMETER :: TEN       = 10.0_fullR 
    REAL(fullR),PARAMETER :: HALF      = 0.500000000000000_fullR 
    REAL(fullR),PARAMETER :: THIRD     = 0.333333333333333_fullR 
    REAL(fullR),PARAMETER :: FOURTH    = 0.250000000000000_fullR 
    REAL(fullR),PARAMETER :: TENTH     = 0.100000000000000_fullR 
    REAL(fullR),PARAMETER :: PI        = 3.141592653589793_fullR

    ! ip=1
    N(1,1)=HALF*(ONE-xi(1))
    N(2,1)=HALF*(ONE+xi(1))

    ! ip=2
    N(1,2)=HALF*(ONE-xi(2))
    N(2,2)=HALF*(ONE+xi(2))

END SUBROUTINE CalcN


!==========================================================================
! only need the integration points themselves here for higher order elements
!==========================================================================
SUBROUTINE CalcdNdxi(dN_dxi)

    IMPLICIT NONE

    REAL(kind(0.d0)),INTENT(OUT) :: dN_dxi(2,2)

    INTEGER,PARAMETER :: singI=kind(0)
    INTEGER,PARAMETER :: fullR=kind(0.d0)

    REAL(fullR),PARAMETER :: ZERO      = 0.0_fullR 
    REAL(fullR),PARAMETER :: ROOT_EPS  = 0.0000000596046448_fullR ! 2^(-24) = 5.96046448E-8
    REAL(fullR),PARAMETER :: ONE       = 1.0_fullR 
    REAL(fullR),PARAMETER :: TWO       = 2.0_fullR 
    REAL(fullR),PARAMETER :: THREE     = 3.0_fullR 
    REAL(fullR),PARAMETER :: FOUR      = 4.0_fullR 
    REAL(fullR),PARAMETER :: FIVE      = 5.0_fullR 
    REAL(fullR),PARAMETER :: SIX       = 6.0_fullR
    REAL(fullR),PARAMETER :: SEVEN     = 7.0_fullR
    REAL(fullR),PARAMETER :: EIGHT     = 8.0_fullR 
    REAL(fullR),PARAMETER :: NINE      = 9.0_fullR 
    REAL(fullR),PARAMETER :: TEN       = 10.0_fullR 
    REAL(fullR),PARAMETER :: HALF      = 0.500000000000000_fullR 
    REAL(fullR),PARAMETER :: THIRD     = 0.333333333333333_fullR 
    REAL(fullR),PARAMETER :: FOURTH    = 0.250000000000000_fullR 
    REAL(fullR),PARAMETER :: TENTH     = 0.100000000000000_fullR 
    REAL(fullR),PARAMETER :: PI        = 3.141592653589793_fullR

    ! ip=1
    dN_dxi(1,1)=-HALF
    dN_dxi(2,1)=+HALF

    ! ip=2
    dN_dxi(1,2)=-HALF
    dN_dxi(2,2)=+HALF

END SUBROUTINE CalcdNdxi


!==========================================================================
SUBROUTINE CalcCond(T,xbar,k,dk_dT)

    IMPLICIT NONE

    REAL(kind(0.d0)),INTENT(IN) :: T,xbar
    REAL(kind(0.d0)),INTENT(OUT) :: k,dk_dT

    INTEGER,PARAMETER :: singI=kind(0)
    INTEGER,PARAMETER :: fullR=kind(0.d0)

    REAL(fullR),PARAMETER :: ZERO      = 0.0_fullR 
    REAL(fullR),PARAMETER :: ROOT_EPS  = 0.0000000596046448_fullR ! 2^(-24) = 5.96046448E-8
    REAL(fullR),PARAMETER :: ONE       = 1.0_fullR 
    REAL(fullR),PARAMETER :: TWO       = 2.0_fullR 
    REAL(fullR),PARAMETER :: THREE     = 3.0_fullR 
    REAL(fullR),PARAMETER :: FOUR      = 4.0_fullR 
    REAL(fullR),PARAMETER :: FIVE      = 5.0_fullR 
    REAL(fullR),PARAMETER :: SIX       = 6.0_fullR
    REAL(fullR),PARAMETER :: SEVEN     = 7.0_fullR
    REAL(fullR),PARAMETER :: EIGHT     = 8.0_fullR 
    REAL(fullR),PARAMETER :: NINE      = 9.0_fullR 
    REAL(fullR),PARAMETER :: TEN       = 10.0_fullR 
    REAL(fullR),PARAMETER :: HALF      = 0.500000000000000_fullR 
    REAL(fullR),PARAMETER :: THIRD     = 0.333333333333333_fullR 
    REAL(fullR),PARAMETER :: FOURTH    = 0.250000000000000_fullR 
    REAL(fullR),PARAMETER :: TENTH     = 0.100000000000000_fullR 
    REAL(fullR),PARAMETER :: PI        = 3.141592653589793_fullR

    REAL(fullR) :: TabT(16),T1,T2
    REAL(fullR) :: Tabk(16),k1,k2
    INTEGER(singI) :: NumInts,iInt
    LOGICAL :: FoundFlag

    TabT(1)=-1.0e+10
    TabT(2)=273.0
    TabT(3)=431.6
    TabT(4)=590.0
    TabT(5)=748.7
    TabT(6)=907.2
    TabT(7)=1065.8
    TabT(8)=1224.3
    TabT(9)=1382.9
    TabT(10)=1541.4
    TabT(11)=1650.0
    TabT(12)=1660.0
    TabT(13)=1690.0
    TabT(14)=1700.0
    TabT(15)=3200.0
    TabT(16)=1.0e+10

    Tabk(1)=3.0e-02
    Tabk(2)=3.0e-02
    Tabk(3)=3.0e-02
    Tabk(4)=3.0e-02
    Tabk(5)=3.0e-02
    Tabk(6)=3.0e-02
    Tabk(7)=3.0e-02
    Tabk(8)=3.0e-02
    Tabk(9)=3.0e-02
    Tabk(10)=3.0e-02
    Tabk(11)=3.0e-02
    Tabk(12)=3.0e-02
    Tabk(13)=3.0e-02
    Tabk(14)=3.0e-02
    Tabk(15)=3.0e-02
    Tabk(16)=3.0e-02
    !=============================================
    ! scaling to induce ill-conditioning . . .
    !=============================================
    IF (xbar<=HALF) THEN
!        Tabk=Tabk*(TENTH**2)
    ELSE IF (HALF<xbar) THEN
!        Tabk=Tabk*(TEN**2)
    END IF
    !=============================================
    ! scaling to induce nonlinearity . . .
    !=============================================
    Tabk(1)=Tabk(1)*ONE
    Tabk(2)=Tabk(2)*TWO
    Tabk(3)=Tabk(3)*FIVE
    Tabk(4)=Tabk(4)*TEN
    Tabk(5)=Tabk(5)*TEN*TWO
    Tabk(6)=Tabk(6)*TEN*FIVE
    Tabk(7)=Tabk(7)*TEN*EIGHT
    Tabk(8)=Tabk(8)*TEN*TEN
    Tabk(9)=Tabk(9)*TEN*NINE
    Tabk(10)=Tabk(10)*TEN*EIGHT
    Tabk(11)=Tabk(11)*TEN*SEVEN
    Tabk(12)=Tabk(12)*TEN*SIX
    Tabk(13)=Tabk(13)*TEN*FIVE
    Tabk(14)=Tabk(14)*TEN*FOUR
    Tabk(15)=Tabk(15)*TEN*THREE
    Tabk(16)=Tabk(16)*TEN*TWO


    NumInts=15
    FoundFlag=.FALSE.
    DO iInt=1,NumInts
        T1=TabT(iInt)
        T2=TabT(iInt+1)
        IF (T<T1) CYCLE
        IF (T>T2) CYCLE

        FoundFlag=.TRUE.
        k1=Tabk(iInt)
        k2=Tabk(iInt+1)
        k=k1 + (T-T1)/(T2-T1)*(k2-k1)
        dk_dT=(k2-k1)/(T2-T1)
        EXIT
    END DO

    IF (.NOT.FoundFlag) THEN
        WRITE(UNIT=*,FMT='(a)') "ERROR: conductivities not calculated, quitting . . ."
        STOP
    END IF

END SUBROUTINE CalcCond


!==========================================================================
! Plot (x,y) data using gnuplot
! * legend removal
! * series names
! * axes names
! * title
!
! NOTE: There is a "bug" (probably) in gnuplot.  Basically, "term wxt", which is
! required if one wants to do fancy zooming in and out of plots, does not work
! for the way I originally had my plot files structured.  Thus, I had to re-work
! them, so that the data is in a seperate file.  Lame . . .
!
! Some of the behavior of this subroutine, as well as the options shown below
! require some explanation . . .
!
! * the data is assumed to be in one of two formats, the first being
!
!   x  y1  y2  . . .   yn
!  =======================
!   #   #   #           #
!   #   #   #           #
!   #   #   #           #
!   . . . . . . . . . . .
!
!   For this type of data, it is assumed that each data series y1, y2, .
!   . ., yn has a value corresponding to each value of x.
!
!   Alternately, for data series that do not share a common domain, the data
!   should be formatted like
!
!   x1  y1
!  ========
!    #  #
!    #  #
!    #  #
!   x2  y2
!  ========
!    #  #
!    #  #
!    #  #
!   x3  y3
!  ========
!    #  #
!    #  #
!    #  #
!      .
!      .
!      .
!   xn  yn
!  ========
!    #  #
!    #  #
!    #  #
!
! In this case, the variable "datalength" should also be passed in,
! which defines the length of each data series.
!
! * color codes:
!   1=black
!   2=blue
!   3=red
!   4=green
!
!==========================================================================
SUBROUTINE plotxy(plotdata,filename,datalength, &
                  grid,limits,title1,xlabel,ylabel,series_labels, &
                  legend,lines,points,colors,widths,output_animate)

    REAL(kind(0.d0)),INTENT(IN) :: plotdata(:,:)
    CHARACTER(*),INTENT(IN) :: filename
    INTEGER,INTENT(IN),OPTIONAL :: datalength(:)
    INTEGER,INTENT(IN),OPTIONAL :: grid
    REAL(kind(0.d0)),INTENT(IN),OPTIONAL :: limits(:,:)
    CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: title1
    CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: xlabel,ylabel
    CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: series_labels(:)
    INTEGER,INTENT(IN),OPTIONAL :: legend
    INTEGER,INTENT(IN),OPTIONAL :: lines(:)
    INTEGER,INTENT(IN),OPTIONAL :: points(:)
    INTEGER,INTENT(IN),OPTIONAL :: colors(:)
    INTEGER,INTENT(IN),OPTIONAL :: widths(:)
    INTEGER,INTENT(IN),OPTIONAL :: output_animate

    INTEGER,PARAMETER :: singI=kind(0)
    INTEGER,PARAMETER :: fullR=kind(0.d0)
    INTEGER(singI),PARAMETER :: PLOT_FID=3

    REAL(fullR),PARAMETER :: ZERO      = 0.0_fullR 
    REAL(fullR),PARAMETER :: ROOT_EPS  = 0.0000000596046448_fullR ! 2^(-24) = 5.96046448E-8
    REAL(fullR),PARAMETER :: ONE       = 1.0_fullR 
    REAL(fullR),PARAMETER :: TWO       = 2.0_fullR 
    REAL(fullR),PARAMETER :: THREE     = 3.0_fullR 
    REAL(fullR),PARAMETER :: FOUR      = 4.0_fullR 
    REAL(fullR),PARAMETER :: FIVE      = 5.0_fullR 
    REAL(fullR),PARAMETER :: SIX       = 6.0_fullR
    REAL(fullR),PARAMETER :: SEVEN     = 7.0_fullR
    REAL(fullR),PARAMETER :: EIGHT     = 8.0_fullR 
    REAL(fullR),PARAMETER :: NINE      = 9.0_fullR 
    REAL(fullR),PARAMETER :: TEN       = 10.0_fullR 
    REAL(fullR),PARAMETER :: HALF      = 0.500000000000000_fullR 
    REAL(fullR),PARAMETER :: THIRD     = 0.333333333333333_fullR 
    REAL(fullR),PARAMETER :: FOURTH    = 0.250000000000000_fullR 
    REAL(fullR),PARAMETER :: TENTH     = 0.100000000000000_fullR 
    REAL(fullR),PARAMETER :: PI        = 3.141592653589793_fullR

    INTEGER(singI) :: length,nseries,i
    REAL(fullR) :: minx,maxx
    CHARACTER(LEN=1024) :: textline
    INTEGER(singI) :: llegend,loutput_animate
    INTEGER(singI),POINTER :: llines(:),lpoints(:),lcolors(:),lwidths(:)
    INTEGER(singI) :: linestat,pointstat,colorstat,widthstat
    INTEGER(singI) :: iseries,ndatalength,idatalength,startrow,endrow
    CHARACTER(LEN=1024),POINTER :: filename_complete(:)

    INTEGER :: debug=0

    IF (debug>=1) THEN
        WRITE(UNIT=*,FMT='(a)') "status of optional parameters:"
        IF (PRESENT(datalength)) THEN
            WRITE(UNIT=*,FMT='(a)') "  ""datalength"" exists"
        ELSE
            WRITE(UNIT=*,FMT='(a)') "  ""datalength"" does not exist"
        END IF
        IF (PRESENT(grid)) THEN
            WRITE(UNIT=*,FMT='(a)') "  ""grid"" exists"
        ELSE
            WRITE(UNIT=*,FMT='(a)') "  ""grid"" does not exist"
        END IF
        IF (PRESENT(limits)) THEN
            WRITE(UNIT=*,FMT='(a)') "  ""limits"" exists"
        ELSE
            WRITE(UNIT=*,FMT='(a)') "  ""limits"" does not exist"
        END IF
        IF (PRESENT(title1)) THEN
            WRITE(UNIT=*,FMT='(a)') "  ""title"" exists"
        ELSE
            WRITE(UNIT=*,FMT='(a)') "  ""title"" does not exist"
        END IF
        IF (PRESENT(xlabel)) THEN
            WRITE(UNIT=*,FMT='(a)') "  ""xlabel"" exists"
        ELSE
            WRITE(UNIT=*,FMT='(a)') "  ""xlabel"" does not exist"
        END IF
        IF (PRESENT(ylabel)) THEN
            WRITE(UNIT=*,FMT='(a)') "  ""ylabel"" exists"
        ELSE
            WRITE(UNIT=*,FMT='(a)') "  ""ylabel"" does not exist"
        END IF
        IF (PRESENT(series_labels)) THEN
            WRITE(UNIT=*,FMT='(a)') "  ""series_labels"" exists"
        ELSE
            WRITE(UNIT=*,FMT='(a)') "  ""series_labels"" does not exist"
        END IF
        IF (PRESENT(legend)) THEN
            WRITE(UNIT=*,FMT='(a)') "  ""legend"" exists"
        ELSE
            WRITE(UNIT=*,FMT='(a)') "  ""legend"" does not exist"
        END IF
        IF (PRESENT(lines)) THEN
            WRITE(UNIT=*,FMT='(a)') "  ""lines"" exists"
        ELSE
            WRITE(UNIT=*,FMT='(a)') "  ""lines"" does not exist"
        END IF
        IF (PRESENT(points)) THEN
            WRITE(UNIT=*,FMT='(a)') "  ""points"" exists"
        ELSE
            WRITE(UNIT=*,FMT='(a)') "  ""points"" does not exist"
        END IF
        IF (PRESENT(colors)) THEN
            WRITE(UNIT=*,FMT='(a)') "  ""colors"" exists"
        ELSE
            WRITE(UNIT=*,FMT='(a)') "  ""colors"" does not exist"
        END IF
        IF (PRESENT(widths)) THEN
            WRITE(UNIT=*,FMT='(a)') "  ""widths"" exists"
        ELSE
            WRITE(UNIT=*,FMT='(a)') "  ""widths"" does not exist"
        END IF
        IF (PRESENT(output_animate)) THEN
            WRITE(UNIT=*,FMT='(a)') "  ""output_animate"" exists"
        ELSE
            WRITE(UNIT=*,FMT='(a)') "  ""output_animate"" does not exist"
        END IF
    END IF

    length=SIZE(plotdata,1)
    nseries=SIZE(plotdata,2)-1
    IF (PRESENT(datalength)) THEN
        ndatalength=SIZE(datalength)
    ELSE
        ndatalength=1
    END IF 

    OPEN(UNIT=PLOT_FID,FILE=TRIM(filename)//".plot",ACTION="WRITE")
    IF (debug>=1) THEN
        WRITE(UNIT=*,FMT='(a)') "writing plot commands"
    END IF

    loutput_animate=0
    IF (PRESENT(output_animate)) THEN
        loutput_animate=output_animate
    END IF

    !
    ! write to a file, as opposed to the screen
    !
    ! note that it is not guaranteed that all plot types will be available for
    ! all gnuplot installations
    !
    IF (loutput_animate==1) THEN
        WRITE(UNIT=PLOT_FID,FMT='(a,a,a)') "set terminal png truecolor size 1600,1200; set output '",TRIM(filename),".png'"
    ELSE
        WRITE(UNIT=PLOT_FID,FMT='(a)') "set terminal wxt size 1600,1200"
!        WRITE(UNIT=PLOT_FID,FMT='(a)') "set terminal X11"
    END IF

    !
    ! set min and max value for x and y axes
    !
    IF (.not.PRESENT(limits)) THEN
        minx=MINVAL(plotdata(:,1))
        maxx=MAXVAL(plotdata(:,1))
        WRITE(UNIT=PLOT_FID,FMT='(a,es22.15,a,es22.15,a)') "set xrange [",(minx-(maxx-minx)/TEN),":",(maxx+(maxx-minx)/TEN),"]"
        minx = HUGE(minx)
        maxx = -HUGE(maxx)
        DO i = 1,nseries
            minx=MIN(minx,MINVAL(plotdata(:,i+1)))
            maxx=MAX(maxx,MAXVAL(plotdata(:,i+1)))
        END DO
        WRITE(UNIT=PLOT_FID,FMT='(a,es22.15,a,es22.15,a)') "set yrange [",(minx-(maxx-minx)/TEN),":",(maxx+(maxx-minx)/TEN),"]"
    ELSE
        minx=limits(1,1)
        maxx=limits(1,2)
        WRITE(UNIT=PLOT_FID,FMT='(a,es22.15,a,es22.15,a)') "set xrange [",(minx),":",(maxx),"]"
        minx=limits(2,1)
        maxx=limits(2,2)
        WRITE(UNIT=PLOT_FID,FMT='(a,es22.15,a,es22.15,a)') "set yrange [",(minx),":",(maxx),"]"
    END IF

    !
    ! set the aspect ratio
    !
    WRITE(UNIT=PLOT_FID,FMT='(a)') "set size square"

    !
    ! turn on the grid
    !
    IF (PRESENT(grid)) THEN
        IF (grid==1) WRITE(UNIT=PLOT_FID,FMT='(a)') "set grid"
    ELSE
        WRITE(UNIT=PLOT_FID,FMT='(a)') "set grid"
    END IF

    !
    ! add a title
    !
    IF (PRESENT(title1)) THEN
        WRITE(UNIT=PLOT_FID,FMT='(a,a,a)') "set title """,TRIM(title1),""""
    END IF

    !
    ! add axis labels
    !
    IF (PRESENT(xlabel)) WRITE(UNIT=PLOT_FID,FMT='(a,a,a)') "set xlabel """,TRIM(xlabel),""""
    IF (PRESENT(ylabel)) WRITE(UNIT=PLOT_FID,FMT='(a,a,a)') "set ylabel """,TRIM(ylabel),""""

    !
    ! add various other stuff
    !
    llegend=0
    IF (PRESENT(legend)) THEN
        llegend=legend
    END IF
    ALLOCATE(llines(1))
    llines=0
    IF (PRESENT(lines)) THEN
        DEALLOCATE(llines)
        ALLOCATE(llines(SIZE(lines,1)))
        llines=lines
    END IF
    ALLOCATE(lpoints(1))
    lpoints=0
    IF (PRESENT(points)) THEN
        DEALLOCATE(lpoints)
        ALLOCATE(lpoints(SIZE(points,1)))
        lpoints=points
    END IF
    ALLOCATE(lcolors(1))
    lcolors=1
    IF (PRESENT(colors)) THEN
        DEALLOCATE(lcolors)
        ALLOCATE(lcolors(SIZE(colors,1)))
        lcolors=colors
    END IF
    ALLOCATE(lwidths(1))
    lwidths=1
    IF (PRESENT(widths)) THEN
        DEALLOCATE(lwidths)
        ALLOCATE(lwidths(SIZE(widths,1)))
        lwidths=widths
    END IF

    IF (llegend==0) THEN
        WRITE(UNIT=PLOT_FID,FMT='(a,a,a)') "set nokey"
    END IF

    ! create the names of the data files and store them
    IF (nseries>1) THEN
        ALLOCATE(filename_complete(nseries))
    ELSE IF (ndatalength>=1) THEN
        ALLOCATE(filename_complete(ndatalength))
    END IF
    IF (nseries>1) THEN
        DO iseries=1,nseries
            WRITE(UNIT=filename_complete(iseries),FMT='(a,a,i0.1,a)') TRIM(filename),"_data",iseries,".plot"
            IF (debug>=1) WRITE(UNIT=*,FMT='(a,a,a)') "creating plot data file named """,TRIM(filename_complete(iseries)),""""
        END DO
    ELSE IF (nseries==1) THEN
        DO idatalength=1,ndatalength
            WRITE(UNIT=filename_complete(idatalength),FMT='(a,a,i0.1,a)') TRIM(filename),"_data",idatalength,".plot"
            IF (debug>=1) WRITE(UNIT=*,FMT='(a,a,a)') "creating plot data file named """,TRIM(filename_complete(idatalength)),""""
        END DO
    END IF

    IF (nseries>1) THEN
        DO iseries=1,nseries
            IF (iseries==1) THEN
                WRITE(UNIT=textline,FMT='(a,a,a)') "plot """,TRIM(filename_complete(iseries)),""""
            ELSE
                WRITE(UNIT=textline,FMT='(a,a,a)') "     """,TRIM(filename_complete(iseries)),""""
            END IF
            IF (PRESENT(series_labels)) THEN
                WRITE(UNIT=textline,FMT='(a,a,a,a)') TRIM(textline), " title """,TRIM(series_labels(iseries)),""""
            END IF
            IF (SIZE(llines,1)==1) THEN
                linestat=llines(1)
            ELSE
                linestat=llines(iseries)
            END IF
            IF (SIZE(lpoints,1)==1) THEN
                pointstat=lpoints(1)
            ELSE
                pointstat=lpoints(iseries)
            END IF
            IF ((linestat==1).AND.(pointstat==0)) THEN
                WRITE(UNIT=textline,FMT='(a,a)') TRIM(textline), " with lines"
            ELSE IF ((linestat==0).AND.(pointstat==1)) THEN
                WRITE(UNIT=textline,FMT='(a,a)') TRIM(textline), " with points"
            ELSE IF ((linestat==1).AND.(pointstat==1)) THEN
                WRITE(UNIT=textline,FMT='(a,a)') TRIM(textline), " with linespoints"
            END IF
            colorstat=0
            IF (SIZE(lcolors,1)==1) THEN
                colorstat=lcolors(1)
            ELSE
                colorstat=lcolors(iseries)
            END IF
            IF (colorstat==1) THEN
                WRITE(UNIT=textline,FMT='(a,a)') TRIM(textline), " linecolor -1"
            ELSE IF (colorstat==2) THEN
                WRITE(UNIT=textline,FMT='(a,a)') TRIM(textline), " linecolor 3"
            ELSE IF (colorstat==3) THEN
                WRITE(UNIT=textline,FMT='(a,a)') TRIM(textline), " linecolor 1"
            ELSE IF (colorstat==4) THEN
                WRITE(UNIT=textline,FMT='(a,a)') TRIM(textline), " linecolor 2"
            END IF
            widthstat=0
            IF (SIZE(lwidths,1)==1) THEN
                widthstat=lwidths(1)
            ELSE
                widthstat=lwidths(iseries)
            END IF
            IF (widthstat>1) THEN
                WRITE(UNIT=textline,FMT='(a,a,i0.1)') TRIM(textline), " linewidth ",widthstat
            END IF
            IF ((nseries>1) .AND. (iseries<nseries)) THEN
                WRITE(UNIT=textline,FMT='(a,a)') TRIM(textline), ", \"
            END IF
            WRITE(UNIT=PLOT_FID,FMT='(a)') TRIM(textline)
        END DO
    ELSE IF (nseries==1) THEN
        DO idatalength=1,ndatalength
            IF (idatalength==1) THEN
                WRITE(UNIT=textline,FMT='(a,a,a)') "plot """,TRIM(filename_complete(idatalength)),""""
            ELSE
                WRITE(UNIT=textline,FMT='(a,a,a)') "     """,TRIM(filename_complete(idatalength)),""""
            END IF
            IF (PRESENT(series_labels)) THEN
                WRITE(UNIT=textline,FMT='(a,a,a,a)') TRIM(textline), " title """,TRIM(series_labels(idatalength)),""""
            END IF
            IF (SIZE(llines,1)==1) THEN
                linestat=llines(1)
            ELSE
                linestat=llines(idatalength)
            END IF
            IF (SIZE(lpoints,1)==1) THEN
                pointstat=lpoints(1)
            ELSE
                pointstat=lpoints(idatalength)
            END IF
            IF ((linestat==1).AND.(pointstat==0)) THEN
                WRITE(UNIT=textline,FMT='(a,a)') TRIM(textline), " with lines"
            ELSE IF ((linestat==0).AND.(pointstat==1)) THEN
                WRITE(UNIT=textline,FMT='(a,a)') TRIM(textline), " with points"
            ELSE IF ((linestat==1).AND.(pointstat==1)) THEN
                WRITE(UNIT=textline,FMT='(a,a)') TRIM(textline), " with linespoints"
            END IF
            colorstat=0
            IF (SIZE(lcolors,1)==1) THEN
                colorstat=lcolors(1)
            ELSE
                colorstat=lcolors(idatalength)
            END IF
            IF (colorstat==1) THEN
                WRITE(UNIT=textline,FMT='(a,a)') TRIM(textline), " linecolor -1"
            ELSE IF (colorstat==2) THEN
                WRITE(UNIT=textline,FMT='(a,a)') TRIM(textline), " linecolor 3"
            ELSE IF (colorstat==3) THEN
                WRITE(UNIT=textline,FMT='(a,a)') TRIM(textline), " linecolor 1"
            ELSE IF (colorstat==4) THEN
                WRITE(UNIT=textline,FMT='(a,a)') TRIM(textline), " linecolor 2"
            END IF
            widthstat=0
            IF (SIZE(lwidths,1)==1) THEN
                widthstat=lwidths(1)
            ELSE
                widthstat=lwidths(idatalength)
            END IF
            IF (widthstat>1) THEN
                WRITE(UNIT=textline,FMT='(a,a,i0.1)') TRIM(textline), " linewidth ",widthstat
            END IF
            IF ((ndatalength>1) .AND. (idatalength<ndatalength)) THEN
                WRITE(UNIT=textline,FMT='(a,a)') TRIM(textline), ", \"
            END IF
            WRITE(UNIT=PLOT_FID,FMT='(a)') TRIM(textline)
        END DO
    END IF


    IF (loutput_animate==0) THEN
        WRITE(UNIT=PLOT_FID,FMT='(a)') "pause -1"

        ! do it all over again, to create figures for journal articles
        WRITE(UNIT=PLOT_FID,FMT='(a)') ""
        WRITE(UNIT=PLOT_FID,FMT='(a)') "# comment the following line to create a journal-formatted .eps file"
        WRITE(UNIT=PLOT_FID,FMT='(a)') "exit"
        WRITE(UNIT=PLOT_FID,FMT='(a)') ""

        WRITE(UNIT=PLOT_FID,FMT='(a,a,a)') &
            "set terminal postscript enhanced color; set output '",TRIM(filename),".eps'"
        WRITE(UNIT=PLOT_FID,FMT='(a)') &
            "# use the following command to convert the .eps file to a .pdf file:"
        WRITE(UNIT=PLOT_FID,FMT='(a,a,a,a,a)') &
            "# ps2pdf -dEPSCrop ",TRIM(filename),".eps ",TRIM(filename),".pdf"

        !
        ! set min and max value for x and y axes
        !
        IF (.NOT.PRESENT(limits)) THEN
            minx=MINVAL(plotdata(:,1))
            maxx=MAXVAL(plotdata(:,1))
            WRITE(UNIT=PLOT_FID,FMT='(a,es22.15,a,es22.15,a)') "set xrange [",(minx-(maxx-minx)/TEN),":",(maxx+(maxx-minx)/TEN),"]"
            minx = HUGE(minx)
            maxx = -HUGE(maxx)
            DO i = 1,nseries
                minx=MIN(minx,MINVAL(plotdata(:,i+1)))
                maxx=MAX(maxx,MAXVAL(plotdata(:,i+1)))
            END DO
            WRITE(UNIT=PLOT_FID,FMT='(a,es22.15,a,es22.15,a)') "set yrange [",(minx-(maxx-minx)/TEN),":",(maxx+(maxx-minx)/TEN),"]"
        ELSE
            minx=limits(1,1)
            maxx=limits(1,2)
            WRITE(UNIT=PLOT_FID,FMT='(a,es22.15,a,es22.15,a)') "set xrange [",(minx),":",(maxx),"]"
            minx=limits(2,1)
            maxx=limits(2,2)
            WRITE(UNIT=PLOT_FID,FMT='(a,es22.15,a,es22.15,a)') "set yrange [",(minx),":",(maxx),"]"
        END IF

        !
        ! set the aspect ratio
        !
        WRITE(UNIT=PLOT_FID,FMT='(a)') "set size square"

        !
        ! turn on the grid
        !
        IF (PRESENT(grid)) THEN
            IF (grid==1) THEN
                WRITE(UNIT=PLOT_FID,FMT='(a)') "set grid"
            END IF
        ELSE
            WRITE(UNIT=PLOT_FID,FMT='(a)') "set grid"
        END IF

        !
        ! add a title
        !
        IF (PRESENT(title1)) THEN
            WRITE(UNIT=PLOT_FID,FMT='(a)') "set title"
        END IF

        !
        ! add axis labels
        !
        WRITE(UNIT=PLOT_FID,FMT='(a)') "set xlabel"
        WRITE(UNIT=PLOT_FID,FMT='(a)') "set ylabel"

        !
        ! add various other stuff
        !
        llegend=0
        IF (PRESENT(legend)) THEN
            llegend=legend
        END IF
        ALLOCATE(llines(1))
        llines=0
        IF (PRESENT(lines)) THEN
            DEALLOCATE(llines)
            ALLOCATE(llines(SIZE(lines,1)))
            llines=lines
        END IF
        ALLOCATE(lpoints(1))
        lpoints=0
        IF (PRESENT(points)) THEN
            DEALLOCATE(lpoints)
            ALLOCATE(lpoints(SIZE(points,1)))
            lpoints=points
        END IF
        ALLOCATE(lcolors(1))
        lcolors=1
        IF (PRESENT(colors)) THEN
            DEALLOCATE(lcolors)
            ALLOCATE(lcolors(SIZE(colors,1)))
            lcolors=colors
        END IF
        ALLOCATE(lwidths(1))
        lwidths=1
        IF (PRESENT(widths)) THEN
            DEALLOCATE(lwidths)
            ALLOCATE(lwidths(SIZE(widths,1)))
            lwidths=widths
        END IF

        IF (llegend==0) THEN
            WRITE(UNIT=PLOT_FID,FMT='(a,a,a)') "set nokey"
        END IF

        ! create the names of the data files and store them
        IF (nseries>1) THEN
            ALLOCATE(filename_complete(nseries))
        ELSE IF (ndatalength>=1) THEN
            ALLOCATE(filename_complete(ndatalength))
        END IF
        IF (nseries>1) THEN
            DO iseries=1,nseries
                WRITE(UNIT=filename_complete(iseries),FMT='(a,a,i0.1,a)') TRIM(filename),"_data",iseries,".plot"
                IF (debug>=1) THEN
                    WRITE(UNIT=*,FMT='(a,a,a)') "creating plot data file named """,TRIM(filename_complete(iseries)),""""
                END IF
            END DO
        ELSE IF (nseries==1) THEN
            DO idatalength=1,ndatalength
                WRITE(UNIT=filename_complete(idatalength),FMT='(a,a,i0.1,a)') TRIM(filename),"_data",idatalength,".plot"
                IF (debug>=1) THEN
                    WRITE(UNIT=*,FMT='(a,a,a)') "creating plot data file named """,TRIM(filename_complete(idatalength)),""""
                END IF
            END DO
        END IF

        IF (nseries>1) THEN
            DO iseries=1,nseries
                IF (iseries==1) THEN
                    WRITE(UNIT=textline,FMT='(a,a,a)') "plot """,TRIM(filename_complete(iseries)),""""
                ELSE
                    WRITE(UNIT=textline,FMT='(a,a,a)') "     """,TRIM(filename_complete(iseries)),""""
                END IF
                IF (PRESENT(series_labels)) THEN
                    WRITE(UNIT=textline,FMT='(a,a,a,a)') TRIM(textline), " title """,TRIM(series_labels(iseries)),""""
                END IF
                IF (SIZE(llines,1)==1) THEN
                    linestat=llines(1)
                ELSE
                    linestat=llines(iseries)
                END IF
                IF (SIZE(lpoints,1)==1) THEN
                    pointstat=lpoints(1)
                ELSE
                    pointstat=lpoints(iseries)
                END IF
                WRITE(UNIT=textline,FMT='(a,a)') TRIM(textline), " with lines"
                WRITE(UNIT=textline,FMT='(a,a)') TRIM(textline), " linecolor 3"
                WRITE(UNIT=textline,FMT='(a,a)') TRIM(textline), " linewidth 4"
                IF ((nseries>1) .AND. (iseries<nseries)) THEN
                    WRITE(UNIT=textline,FMT='(a,a)') TRIM(textline), ", \"
                END IF
                WRITE(UNIT=PLOT_FID,FMT='(a)') TRIM(textline)
            END DO
        ELSE IF (nseries==1) THEN
            DO idatalength=1,ndatalength
                IF (idatalength==1) THEN
                    WRITE(UNIT=textline,FMT='(a,a,a)') "plot """,TRIM(filename_complete(idatalength)),""""
                ELSE
                    WRITE(UNIT=textline,FMT='(a,a,a)') "     """,TRIM(filename_complete(idatalength)),""""
                END IF
                IF (PRESENT(series_labels)) THEN
                    WRITE(UNIT=textline,FMT='(a,a,a,a)') TRIM(textline), " title """,TRIM(series_labels(idatalength)),""""
                END IF
                IF (SIZE(llines,1)==1) THEN
                    linestat=llines(1)
                ELSE
                    linestat=llines(idatalength)
                END IF
                IF (SIZE(lpoints,1)==1) THEN
                    pointstat=lpoints(1)
                ELSE
                    pointstat=lpoints(idatalength)
                END IF
                WRITE(UNIT=textline,FMT='(a,a)') TRIM(textline), " with lines"
                WRITE(UNIT=textline,FMT='(a,a)') TRIM(textline), " linecolor 3"
                WRITE(UNIT=textline,FMT='(a,a)') TRIM(textline), " linewidth 4"
                WRITE(UNIT=PLOT_FID,FMT='(a)') TRIM(textline)
            END DO
        END IF

        WRITE(UNIT=PLOT_FID,FMT='(a)') ""
        CLOSE(UNIT=PLOT_FID)

    END IF ! loutput_animate

    !
    ! plot the data itself, taking into account the continuity type
    !
    IF (nseries>1) THEN
        DO iseries=1,nseries
            IF (debug>=1) THEN
                WRITE(UNIT=*,FMT='(a,a,a)') "writing plot data to file """,TRIM(filename_complete(iseries)),""""
            END IF
            OPEN(UNIT=PLOT_FID,FILE=TRIM(filename_complete(iseries)))

            WRITE(UNIT=PLOT_FID,FMT='(a)') "#   X1      X2"
            linestat=0
            DO i = 1,length
                WRITE(UNIT=PLOT_FID,FMT='(e22.15,a,e22.15)') plotdata(i,1)," ",plotdata(i,iseries+1)
            END DO

        END DO

    ELSE IF (nseries==1) THEN
        DO idatalength=1,ndatalength
            IF (debug>=1) THEN
                WRITE(UNIT=*,FMT='(a,a,a)') "writing plot data to file """,TRIM(filename_complete(idatalength)),""""
            END IF
            OPEN(UNIT=PLOT_FID,FILE=TRIM(filename_complete(idatalength)))

            IF (idatalength==1) THEN
                WRITE(UNIT=PLOT_FID,FMT='(a)') "#   X1      X2"
            ELSE
                WRITE(UNIT=PLOT_FID,FMT='(a)') ""
            END IF
            IF (PRESENT(datalength)) THEN
                startrow=SUM(datalength(1:idatalength-1))+1
                endrow=startrow+datalength(idatalength)-1
            ELSE
                startrow=1
                endrow=length
            END IF
            linestat=0
            DO i = startrow,endrow
                WRITE(UNIT=PLOT_FID,FMT='(e22.15,a,e22.15)') plotdata(i,1)," ",plotdata(i,2)
            END DO

        END DO

    END IF

    CLOSE(UNIT=PLOT_FID)

    DEALLOCATE(llines)
    DEALLOCATE(lpoints)
    DEALLOCATE(filename_complete)

END SUBROUTINE plotxy


