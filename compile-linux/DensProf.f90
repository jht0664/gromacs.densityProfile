PROGRAM DENSITYPROFILE
  INTEGER, PARAMETER :: FILE_LEN = 20

INTERFACE
  SUBROUTINE READ_XTC_R(INPUT_NAME,FILE_LEN,AR_DATA,BOX,N_ATOM,N_DATA)
    USE xtc_mod, only: xtcfile ! Use the xdr interface
    TYPE(xtcfile) :: xtc ! Declare a variable of type xtcfile
    INTEGER, INTENT(IN) :: FILE_LEN
    CHARACTER(LEN=FILE_LEN), INTENT(IN) :: INPUT_NAME
    REAL, DIMENSION(:,:,:), ALLOCATABLE, INTENT(OUT) :: AR_DATA
    REAL, DIMENSION(:,:,:), ALLOCATABLE :: TEMP
    INTEGER :: MAX_DATA = 3000, ICR = 1000, IERR, I, J, K, L
    INTEGER, INTENT(OUT) :: N_DATA, N_ATOM ! Total time steps & total number of atoms
    REAL, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: BOX
    REAL, DIMENSION(:,:), ALLOCATABLE :: TEMP_BOX
    REAL, DIMENSION(1:3) :: BOX_INIT
  END SUBROUTINE READ_XTC_R
  SUBROUTINE POS_NOJUMP_PBC_3DR(POSI,BOX,N_ATOM,N_DATA)
    INTEGER, INTENT(IN) :: N_ATOM, N_DATA
    REAL, DIMENSION(1:3,1:N_ATOM,1:N_DATA), INTENT(INOUT) :: POSI
    REAL, DIMENSION(1:3,1:N_DATA), INTENT(IN) :: BOX
    INTEGER :: I,J,K
  END SUBROUTINE POS_NOJUMP_PBC_3DR
  SUBROUTINE SLICE_BOX_3DR(D_XYZ,BOX,N_DATA,T_SLICE,N_SLICE)
    REAL, DIMENSION(1:3), INTENT(IN) :: D_XYZ
    REAL, DIMENSION(1:3,1:N_DATA), INTENT(IN) :: BOX
    INTEGER, DIMENSION(1:3), INTENT(OUT) :: N_SLICE
    INTEGER, INTENT(IN) :: N_DATA
    INTEGER :: I,J,K, IERR
    REAL, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: T_SLICE
    REAL, DIMENSION(1:3) :: T_SLICE_AVG, T_SLICE_STD
  END SUBROUTINE SLICE_BOX_3DR
  SUBROUTINE DENS_3DR(POSI,T_SLICE,N_SLICE,N_ATOM,N_DATA,CONC_3D)
    INTEGER, INTENT(IN) :: N_ATOM, N_DATA
    REAL, DIMENSION(1:3,1:N_DATA), INTENT(IN) :: T_SLICE
    REAL, DIMENSION(1:3,1:N_ATOM,1:N_DATA), INTENT(IN) :: POSI
    INTEGER :: I,J,K,L,IERR,N_CUBE,I_CUBE
    REAL, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: CONC_3D
    INTEGER, DIMENSION(1:3), INTENT(IN) :: N_SLICE
    INTEGER, DIMENSION(1:4) :: T_INDEX
  END SUBROUTINE DENS_3DR
  SUBROUTINE BLOCK_AVG_2DR(MATX_2D,N_VAL,N_TIME,SIZE_BLOCK,BLOCK,N_BLOCK)
    INTEGER, INTENT(IN) :: N_VAL, N_TIME, SIZE_BLOCK
    REAL, DIMENSION(1:N_VAL,1:N_TIME), INTENT(IN) :: MATX_2D
    REAL, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: BLOCK
    INTEGER :: I,J,K, I_BLOCK
    INTEGER, INTENT(OUT) :: N_BLOCK 
  END SUBROUTINE BLOCK_AVG_2DR
  SUBROUTINE DENS_AXIS_1DR(AXIS,CONC_3D,N_SLICE,N_TIME,DENS_AXIS)
    INTEGER, INTENT(IN) :: AXIS ! 1, 2, or 3 menas average on x, y, or z axis.
    INTEGER, DIMENSION(1:3), INTENT(IN) :: N_SLICE
    INTEGER, INTENT(IN) :: N_TIME
    REAL, DIMENSION(1:PRODUCT(N_SLICE),1:N_TIME), INTENT(IN) :: CONC_3D
    INTEGER, DIMENSION(1:N_SLICE(1),1:N_SLICE(2),1:N_SLICE(3)) :: TEMP_3D
    INTEGER :: I,J,K
    REAL, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: DENS_AXIS
  END SUBROUTINE DENS_AXIS_1DR
  SUBROUTINE ANIMATION_1DR(MATX_2D,N_X,N_T,SAVE_FILE,FILE_LEN)
    INTEGER, INTENT(IN) :: N_X, N_T, FILE_LEN
    REAL, DIMENSION(1:N_X,1:N_T), INTENT(IN) :: MATX_2D
    CHARACTER(LEN=FILE_LEN), INTENT(IN) :: SAVE_FILE
    CHARACTER(LEN=FILE_LEN) :: myfmt
    INTEGER :: I,J
  END SUBROUTINE ANIMATION_1DR
  SUBROUTINE SAVE_3DR(MATX_2D,N_X,N_T,SAVE_FILE,FILE_LEN)
    INTEGER, INTENT(IN) :: FILE_LEN
    INTEGER, DIMENSION(1:3), INTENT(IN) :: N_X
    INTEGER, INTENT(IN) :: N_T
    REAL, DIMENSION(1:PRODUCT(N_X),1:N_T), INTENT(IN) :: MATX_2D
    CHARACTER(LEN=FILE_LEN), INTENT(IN) :: SAVE_FILE
    INTEGER :: I,J,K,L
  END SUBROUTINE SAVE_3DR
  SUBROUTINE SET_BASIS_MATRIX_B_UNI(MATX_2D,N_X,N_T,B_SIZE,N_SAMP)
  INTEGER, INTENT(IN) :: N_X, N_T, B_SIZE, N_SAMP
  INTEGER :: DATA_SIZE
  INTEGER :: I, J, K, K_HIGH, LEFT
  CHARACTER :: MARK
  REAL(KIND=8), DIMENSION(1:B_SIZE,1:B_SIZE) :: BASIS
  REAL, DIMENSION(1:N_X,1:N_T), INTENT(INOUT) :: MATX_2D
  REAL(KIND=8), DIMENSION(1:N_X) :: DATA_X, DATA_Y ! MATX_2D is divided to DATA_X and DATA_Y
  REAL(KIND=8) :: X_HIGH, X_LOW, X_VAL, Y_VAL
  END SUBROUTINE SET_BASIS_MATRIX_B_UNI
  SUBROUTINE FIT_RMSD_1DR(MATX_2D,N_X,N_T,SAVE_FILE,FILE_LEN)
    INTEGER, INTENT(IN) :: N_X, N_T, FILE_LEN
    CHARACTER(LEN=FILE_LEN), INTENT(IN) :: SAVE_FILE
    REAL, DIMENSION(1:N_X,1:N_T), INTENT(INOUT) :: MATX_2D
    REAL, DIMENSION(1:N_X) :: TEMP_2D
    REAL, DIMENSION(0:N_X-1) :: RMSD
    INTEGER :: I,J,K,L,M
    INTEGER :: TRANS_MIN_RMSD, TRY_TRANS, NEW_TRANS
  END SUBROUTINE FIT_RMSD_1DR
  SUBROUTINE TIME_AVG_1DR(MATX_2D,N_X,N_T,SAVE_FILE,FILE_LEN)
    INTEGER, INTENT(IN) :: N_X, N_T, FILE_LEN
    CHARACTER(LEN=FILE_LEN), INTENT(IN) :: SAVE_FILE
    REAL, DIMENSION(1:N_X,1:N_T), INTENT(IN) :: MATX_2D
    REAL :: X_AVG, X_STD
    INTEGER :: I,J
  END SUBROUTINE TIME_AVG_1DR
END INTERFACE

! command
  INTEGER :: NUM_OF_ARGS, I_ARG
  CHARACTER(LEN=FILE_LEN) :: FILE_NAME, INPUT_NAME, OUTPUT_NAME1,OUTPUT_NAME2,OUTPUT_NAME3, PREFIX_NAME, AXIS 
  LOGICAL :: LookForAXIS=.FALSE.,LookForSLICE=.FALSE.,LookForBLOCK=.FALSE.,LookForPREFIX=.FALSE.
  LOGICAL :: FileExist
  INTEGER :: SIZE_BLOCK, SELECT_AXIS
  REAL :: SLICE, CPU_START, CPU_FINISH
! READ_XTC_R
  REAL, DIMENSION(:,:,:), ALLOCATABLE :: AR_DATA
  REAL, DIMENSION(:,:), ALLOCATABLE :: BOX
  INTEGER :: N_ATOM, N_DATA ! Total time steps & total number of atoms
! SLICE_BOX_3DR
  REAL, DIMENSION(1:3) :: D_XYZ
  INTEGER, DIMENSION(1:3) :: N_SLICE
  REAL, DIMENSION(:,:), ALLOCATABLE :: T_SLICE
! DENS_3DR
  REAL, DIMENSION(:,:), ALLOCATABLE :: CONC_3D
! BLOCK_AVG_2DR
  REAL, DIMENSION(:,:), ALLOCATABLE :: BLOCK
! DENS_AXIS_1DR
  REAL, DIMENSION(:,:), ALLOCATABLE :: DENS_AXIS
! Set parameter
  PREFIX_NAME = "D-prof"
  AXIS = "Z"
  SLICE = 0.12D0
  SIZE_BLOCK = 5

! =========== START(Update 03/31/2016) =============
! Command line
  NUM_OF_ARGS = command_argument_count()
  IF (NUM_OF_ARGS == 0) THEN
    WRITE(*,*) "USAGE : ./D-prof.x -pf [PREFIX] -slice [real] -block [integer]" 
    WRITE(*,*) "-pf   D-prof              Prefix of save files"
    WRITE(*,*) "-axis Z                   Density profiling axis (e.g. X, Y, Z, or XY)"
    WRITE(*,*) "-slice 0.12               Slab thickness (nm)" !(usually cut-off electrostatic interaction is 1.4)
    WRITE(*,*) "-block 5                  Size of block (5 = every 5 steps (time) average) for cal. density"
    WRITE(*,*) "PF.xtc              Input xtc file"
    WRITE(*,*) "PF.AXIS.dens.ani    Time change of density w.r.t position (animation) "
    WRITE(*,*) "PF.AXIS.dens.align  Aligned time change of density w.r.t position"
    WRITE(*,*) "PF.AXIS.dens.avg    Block average density from aligned density-position"
    WRITE(*,*) "It calcualtes density profiles."
    WRITE(*,*) "Assume that it has only one interface"
    WRITE(*,*) "xtc file should be created by g_traj -jump option"
    WRITE(*,*) "EXAMPLE: make_ndx_mpi -f conf.gro (select only one kind of particles);"
    WRITE(*,*) "         g_traj_mpi -n index.ndx -nonojump -oxt D-in.xtc; ./D-prof.x -pf D-in"
    STOP
    ELSE IF(NUM_OF_ARGS > 0 ) THEN
    !loop across options
      DO I_ARG=1,NUM_OF_ARGS
        CALL get_command_argument(I_ARG,FILE_NAME)
        FILE_NAME = adjustl(FILE_NAME)
        SELECT CASE(FILE_NAME)
          CASE("-pf")
            LookForPREFIX=.TRUE. 
            CYCLE
          CASE("-axis")
            LookForAXIS=.TRUE. 
            CYCLE
          CASE("-slice")
            LookForSLICE=.TRUE.
            CYCLE
          CASE("-block")
            LookForBLOCK=.TRUE.
            CYCLE
          CASE default
            IF(LookForPREFIX) THEN
                READ(FILE_NAME,*) PREFIX_NAME
                LookForPREFIX = .FALSE.
                INPUT_NAME = TRIM(PREFIX_NAME)//'.xtc'
                inquire(file=INPUT_NAME,exist=FileExist)
                IF(.NOT.FileExist) THEN
                  WRITE(*,*) 'Input File ',FILE_NAME,' is not found'
                  STOP
                ENDIF
                CYCLE
              ELSE IF(LookForAXIS) THEN
                READ(FILE_NAME,*) AXIS
                LookForAXIS = .FALSE.
                CYCLE
              ELSE IF(LookForSLICE) THEN
                READ(FILE_NAME,*) SLICE
                LookForSLICE = .FALSE.
                CYCLE
              ELSE IF(LookForBLOCK) THEN
                READ(FILE_NAME,*) SIZE_BLOCK
                LookForBLOCK = .FALSE.
                CYCLE
              ELSE
                PRINT*, "Option ", FILE_NAME, "is unknown"
                STOP
            ENDIF
         END SELECT
      ENDDO
  ENDIF
  SELECT CASE (AXIS)
    CASE ("X")
      SELECT_AXIS = 1
    CASE ("Y")
      SELECT_AXIS = 2
    CASE ("Z")
      SELECT_AXIS = 3
    CASE ("XY")
      SELECT_AXIS = 4
    CASE default
      PRINT*, "You put the Wrong axis,",AXIS
      STOP
  END SELECT
  OUTPUT_NAME1 = TRIM(PREFIX_NAME)//'.'//TRIM(AXIS)//'.dens.ani' ! raw save file
  OUTPUT_NAME2 = TRIM(PREFIX_NAME)//'.'//TRIM(AXIS)//'.dens.align' ! interface fluctuation 
  OUTPUT_NAME3 = TRIM(PREFIX_NAME)//'.'//TRIM(AXIS)//'.dens.avg' ! set x,y,or z = 0 has interface (assume only one interface)
  inquire(file=TRIM(OUTPUT_NAME1),exist=FileExist)
  IF(FileExist) WRITE(*,*) OUTPUT_NAME1,' already exists'
  inquire(file=TRIM(OUTPUT_NAME2),exist=FileExist)
  IF(FileExist) WRITE(*,*) OUTPUT_NAME2,' already exists'
  inquire(file=TRIM(OUTPUT_NAME3),exist=FileExist)
  IF(FileExist) WRITE(*,*) OUTPUT_NAME3,' already exists'
  PRINT*, "Your command line: D-prof.x -pf ",TRIM(PREFIX_NAME)," -axis ",TRIM(AXIS),&
          " -slice ",SLICE," -block ",SIZE_BLOCK
! =========== END(Update 03/31/2016) =============
  CALL CPU_TIME(CPU_START)
  CALL READ_XTC_R(INPUT_NAME,FILE_LEN,AR_DATA,BOX,N_ATOM,N_DATA) ! read xtc file to save position and box information
  CALL POS_NOJUMP_PBC_3DR(AR_DATA,BOX,N_ATOM,N_DATA)! position goes inside PBC in 3D
  D_XYZ = SLICE ! we make small cubes of box 
  CALL SLICE_BOX_3DR(D_XYZ,BOX,N_DATA,T_SLICE,N_SLICE) ! slice the box in 3D and get thickness of cubes
  CALL DENS_3DR(AR_DATA,T_SLICE,N_SLICE,N_ATOM,N_DATA,CONC_3D) ! density of cubes
  IF(SELECT_AXIS < 4) THEN ! density profiling along 1D
    IF(SIZE_BLOCK /= 1) THEN
      CALL BLOCK_AVG_2DR(CONC_3D,PRODUCT(N_SLICE),N_DATA,SIZE_BLOCK,BLOCK,N_BLOCK) ! Block average
      CALL DENS_AXIS_1DR(SELECT_AXIS,BLOCK,N_SLICE,N_BLOCK,DENS_AXIS) ! average density on z axis
      CALL ANIMATION_1DR(DENS_AXIS,N_SLICE(SELECT_AXIS),N_BLOCK,OUTPUT_NAME1,FILE_LEN)
    ELSE
      CALL DENS_AXIS_1DR(SELECT_AXIS,CONC_3D,N_SLICE,N_DATA,DENS_AXIS)    
      CALL ANIMATION_1DR(DENS_AXIS,N_SLICE(SELECT_AXIS),N_DATA,OUTPUT_NAME1,FILE_LEN)
      N_BLOCK = N_DATA
    ENDIF
!    CALL SET_BASIS_MATRIX_B_UNI(DENS_AXIS,N_SLICE(SELECT_AXIS),N_BLOCK,4,4) ! not completed - Curve Fitting for spline fuction
    CALL FIT_RMSD_1DR(DENS_AXIS,N_SLICE(SELECT_AXIS),N_BLOCK,OUTPUT_NAME2,FILE_LEN) ! find minimum RMSD plots by testing translating
    CALL TIME_AVG_1DR(DENS_AXIS,N_SLICE(SELECT_AXIS),N_BLOCK,OUTPUT_NAME3,FILE_LEN) ! time average of functions
  ELSE ! density profiling along XY
    CALL BLOCK_AVG_2DR(CONC_3D,PRODUCT(N_SLICE),N_DATA,SIZE_BLOCK,BLOCK,N_BLOCK) ! Block average
    CALL SAVE_3DR(BLOCK,N_SLICE,N_BLOCK,OUTPUT_NAME1,FILE_LEN)
  ENDIF
  CALL CPU_TIME(CPU_FINISH)
  WRITE(*,'("TIME =",F6.3," seconds.")') CPU_FINISH-CPU_START
END PROGRAM DENSITYPROFILE

SUBROUTINE TIME_AVG_1DR(MATX_2D,N_X,N_T,SAVE_FILE,FILE_LEN)
!! 02/23/2016 - coding
!! calculate AVG and STD of 1D data
  INTEGER, INTENT(IN) :: N_X, N_T, FILE_LEN
  CHARACTER(LEN=FILE_LEN), INTENT(IN) :: SAVE_FILE
  REAL, DIMENSION(1:N_X,1:N_T), INTENT(IN) :: MATX_2D
  REAL :: X_AVG, X_STD
  INTEGER :: I,J
  PRINT*, "==============================="
  PRINT*, "Time average 1D data"
  OPEN(UNIT=11,FILE=TRIM(SAVE_FILE),STATUS='unknown')
  WRITE(11,*) "# I, X_AVG, X_STD"
  DO I=1,N_X ! select the graph which will translate 
    X_AVG = 0.0D0
    X_STD = 0.0D0
    DO J=1,N_T
      X_AVG = X_AVG + MATX_2D(I,J)
      X_STD = X_STD + MATX_2D(I,J)**2
    ENDDO
    X_AVG = X_AVG/REAL(N_T)
    X_STD = X_STD/REAL(N_T) - X_AVG**2 
    WRITE(11,*) I, X_AVG, SQRT(X_STD)
  ENDDO
  CLOSE(12)
END SUBROUTINE TIME_AVG_1DR

SUBROUTINE FIT_RMSD_1DR(MATX_2D,N_X,N_T,SAVE_FILE,FILE_LEN)
!! 02/23/2016 - coding
!! calculate RMSD of translated plots and find optimum plot with minimum RMSDs 
!! assume the MATX_2D is periodic
  INTEGER, INTENT(IN) :: N_X, N_T, FILE_LEN
  CHARACTER(LEN=FILE_LEN), INTENT(IN) :: SAVE_FILE
  REAL, DIMENSION(1:N_X,1:N_T), INTENT(INOUT) :: MATX_2D
  REAL, DIMENSION(1:N_X) :: TEMP_2D
  REAL, DIMENSION(0:N_X-1) :: RMSD
  INTEGER :: I,J,K,L,M
  INTEGER :: TRANS_MIN_RMSD, TRY_TRANS, NEW_TRANS
  PRINT*, "==============================="
  PRINT*, "Fitting graphs with RMSD"
  DO I=2,N_T ! select the graph which will translate 
    RMSD = 0.0D0
    DO J=I-1,1,-1 ! select the graph(s) which is(are) reference.
!      PRINT*,I,J
      DO L=0,N_X-1 ! amount of translation which will be tested
        DO K=1,N_X ! select the 1D position of reference
          IF(L+K > N_X) THEN ! periodic boundary condition 
            TRY_TRANS = L+K-N_X
          ELSE 
            TRY_TRANS = L+K
          ENDIF
!          PRINT*, L,"Try movement of ",TRY_TRANS,K,I,J
          RMSD(L) = RMSD(L) + SQRT((MATX_2D(TRY_TRANS,I) - MATX_2D(K,J))**2)
        ENDDO
      ENDDO
    ENDDO
!    DO L=0,N_X-1
!      PRINT*, L, RMSD(L)
!    ENDDO
    TRANS_MIN_RMSD = SUM(MINLOC(RMSD))-1 ! finding translation with minimum RMSD 
!    IF(I > 1000) THEN
!      PRINT*, I,TRANS_MIN_RMSD, MINLOC(RMSD)
!      PRINT*, RMSD(0:N_X-1)
!    ENDIF
!    PRINT*, "Optimum translation = ",TRANS_MIN_RMSD," for time step ",I,"RMSD = ",MINVAL(RMSD)
    TEMP_2D = 0.0D0
    DO M=1,N_X ! translate and save temporary matrix
      IF((M-TRANS_MIN_RMSD) <= 0) THEN ! periodic boundary condition
        NEW_TRANS = M-TRANS_MIN_RMSD + N_X
      ELSE
        NEW_TRANS = M-TRANS_MIN_RMSD
      ENDIF
      TEMP_2D(NEW_TRANS) = MATX_2D(M,I)
    ENDDO
    DO M=1,N_X ! assign new translated matrix
      MATX_2D(M,I) = TEMP_2D(M)
!      PRINT*, "New assignment", I, MATX_2D(M,I)
    ENDDO
  ENDDO
  OPEN(UNIT=12,FILE=TRIM(SAVE_FILE),STATUS='unknown')
  WRITE(12,*) "# Time step, X-values, Y-values"
  DO I=1,N_T
    DO J=1,N_X 
      WRITE(12,*) I, J, MATX_2D(J,I)
    ENDDO
    WRITE(12,*)
    WRITE(12,*)
  ENDDO
  CLOSE(12)
END SUBROUTINE FIT_RMSD_1DR

SUBROUTINE SET_BASIS_MATRIX_B_UNI(MATX_2D,N_X,N_T,B_SIZE,N_SAMP)
!! 02/21/2016 - coding from test03 subroutine of spline_prb.f90  
!! Refer: Numerical Recipes in Fortran 90 (http://people.sc.fsu.edu/~jburkardt/f_src/spline/spline.html)  
!! Should be edited later!
  INTEGER, INTENT(IN) :: N_X, N_T, B_SIZE, N_SAMP
  INTEGER :: DATA_SIZE
  INTEGER :: I, J, K, K_HIGH, LEFT
  CHARACTER :: MARK
  REAL(KIND=8) :: BASIS(B_SIZE,B_SIZE)
  REAL, DIMENSION(1:N_X,1:N_T), INTENT(INOUT) :: MATX_2D
  REAL(KIND=8), DIMENSION(0:N_X) :: DATA_X, DATA_Y ! MATX_2D is divided to DATA_X and DATA_Y
  REAL(KIND=8) :: X_HIGH, X_LOW, X_VAL, Y_VAL
  ! CALL BASIS_MATRIX_B_UNI(BASIS) ! initialize
  DO I=0,B_SIZE
    DO J=0,B_SIZE
  !    PRINT*, I,J,BASIS(I,J)
    ENDDO
  ENDDO
  DATA_SIZE = N_X
  DO J=1,N_X !make matrix of data for x
    DATA_X(J) = DBLE(J)
  ENDDO
  DO I=1,N_T
    DO J=1,N_X !make matrix of data for y in a certain time, I
      DATA_Y(J) = DBLE(MATX_2D(J,I))
    ENDDO
    LEFT = 2
    DO J=0,DATA_SIZE
      IF(J == 0) THEN
        X_LOW = DATA_X(1) - 0.5D0*(DATA_X(2)-DATA_X(1))
        X_HIGH = DATA_X(1)
      ELSE IF(J < DATA_SIZE) THEN
        X_LOW = DATA_X(J)
        X_HIGH = DATA_X(J+1)
      ELSE IF(J >= DATA_SIZE) THEN
        X_LOW = DATA_X(DATA_SIZE)
        X_HIGH = DATA_X(DATA_SIZE) + 0.5D0*(DATA_X(DATA_SIZE)-DATA_X(DATA_SIZE-1))
      ELSE
        STOP "something wrong!"
      ENDIF
      IF(J < DATA_SIZE) THEN
        K_HIGH = N_SAMP - 1
      ELSE
        K_HIGH = N_SAMP
      ENDIF
      DO K=0,K_HIGH
        X_VAL = (DBLE(N_SAMP-K)*X_LOW + DBLE(K)*X_HIGH)/DBLE(N_SAMP)
        PRINT*, "XVAL",X_VAL,K
        PRINT*, LEFT,DATA_SIZE,BASIS,DATA_X,DATA_Y,X_VAL,Y_VAL
   !     CALL BASIS_MATRIX_TMP(LEFT,DATA_SIZE,BASIS,DATA_X,DATA_Y,X_VAL,Y_VAL)
        IF(0 < J .and. K == 0) THEN
          MARK = '*'
        ELSE
          MARK = ' '
        ENDIF
        WRITE(*, '(2x,a1,2g14.6)' ) mark, tval, yval
      ENDDO
    ENDDO
  ENDDO
END SUBROUTINE SET_BASIS_MATRIX_B_UNI

SUBROUTINE SAVE_3DR(MATX_2D,N_X,N_T,SAVE_FILE,FILE_LEN)
!! 03/31/2016 gnuplot pm3d plot for X/Y/Z/color 
!! 
  INTEGER, INTENT(IN) :: FILE_LEN
  INTEGER, DIMENSION(1:3), INTENT(IN) :: N_X
  INTEGER, INTENT(IN) :: N_T
  REAL, DIMENSION(1:PRODUCT(N_X),1:N_T), INTENT(IN) :: MATX_2D
  CHARACTER(LEN=FILE_LEN), INTENT(IN) :: SAVE_FILE
!  CHARACTER(LEN=FILE_LEN) :: myfmt
  INTEGER :: I,J,K,L
  PRINT*, "==============================="
  PRINT*, "Making data file for 3D matrix plotting (pm3d)"
  OPEN(UNIT=11,FILE=TRIM(SAVE_FILE),STATUS='unknown')
  WRITE(11,*) "# Data File generated by SAVE_3DR"
  WRITE(11,*) "# MAX value", MAXVAL(MATX_2D)
!  WRITE(myfmt, '("(",I0,"(E20.14,2X))")') N_T 
!  PRINT*, myfmt
  DO L=50,50!N_T
    DO K=1,N_X(3)
    DO I=1,N_X(1)
    DO J=1,N_X(2)
      WRITE(11,"(3(I5,1X),F8.4)") I,J,K,MATX_2D(I+(J-1)*N_X(1)+(K-1)*N_X(1)*N_X(2),L)
    ENDDO
      WRITE(11,*)
    ENDDO
    WRITE(11,*)
    WRITE(11,*)
    ENDDO
  ENDDO
  CLOSE(11)
END SUBROUTINE SAVE_3DR

SUBROUTINE ANIMATION_1DR(MATX_2D,N_X,N_T,SAVE_FILE,FILE_LEN)
!! 02/21/2016 - Ref. http://www.gnuplotting.org/tag/animation/ and http://www.oc.nps.edu/~bird/oc3030_online/fortran/io/io.html
!! Making animations with respect to time
  INTEGER, INTENT(IN) :: N_X, N_T, FILE_LEN
  REAL, DIMENSION(1:N_X,1:N_T), INTENT(IN) :: MATX_2D
  CHARACTER(LEN=FILE_LEN), INTENT(IN) :: SAVE_FILE
  CHARACTER(LEN=FILE_LEN) :: myfmt
  INTEGER :: I,J
  PRINT*, "==============================="
  PRINT*, "Making data file for 1D animation"
  OPEN(UNIT=11,FILE=TRIM(SAVE_FILE),STATUS='unknown')
  WRITE(11,*) "# Data File generated by ANIMATION_1DR"
  WRITE(11,*) "# "
  WRITE(myfmt, '("(",I0,"(E20.14,2X))")') N_T 
!  PRINT*, myfmt
  DO I=1,N_X
    WRITE(11,fmt=myfmt) (MATX_2D(I,J),J=1,N_T)
!    PRINT*, (MATX_2D(I,J),J=1,N_T)
  ENDDO
  CLOSE(11)
END SUBROUTINE ANIMATION_1DR

SUBROUTINE DENS_AXIS_1DR(AXIS,CONC_3D,N_SLICE,N_TIME,DENS_AXIS)
!! 02/21/2016 - subroutine version
!! 01/29/2016 - coding
!! 1D density profiling with respect to principal axis, AXIS, x=1, y=2, z=3
  INTEGER, INTENT(IN) :: AXIS ! 1, 2, or 3 menas average on x, y, or z axis.
  INTEGER, DIMENSION(1:3), INTENT(IN) :: N_SLICE
  INTEGER, INTENT(IN) :: N_TIME
  REAL, DIMENSION(1:PRODUCT(N_SLICE),1:N_TIME), INTENT(IN) :: CONC_3D
  INTEGER, DIMENSION(1:N_SLICE(1),1:N_SLICE(2),1:N_SLICE(3)) :: TEMP_3D
  INTEGER :: I,J,K, IERR
  REAL, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: DENS_AXIS
  PRINT*, "==============================="
  PRINT*, "Starting 1D-Density AVG w.r.t principal axis"
  TEMP_3D = 0
  DO K=1,N_SLICE(3) ! make list
    DO J=1,N_SLICE(2)
      DO I=1,N_SLICE(1)
        TEMP_3D(I,J,K) = I+(J-1)*N_SLICE(1)+(K-1)*N_SLICE(1)*N_SLICE(2)
!        IF(K == 9) PRINT*, I,J,K,TEMP_3D(I,J,K),CONC_3D(TEMP_3D(I,J,K),1)
      ENDDO
    ENDDO
  ENDDO
  ALLOCATE(DENS_AXIS(1:N_SLICE(AXIS),1:N_TIME),stat=IERR) 
  IF (IERR /= 0) STOP "DENS_AXIS in DENS_AXIS_1DR: problem in attempt to allocate memory"
  DENS_AXIS = 0.0D0
  DO L=1,N_TIME ! calculate average density on principal axis
    DO K=1,N_SLICE(3)
      DO J=1,N_SLICE(2)
        DO I=1,N_SLICE(1)
          IF(AXIS == 3) THEN
            DENS_AXIS(K,L) = DENS_AXIS(K,L) + CONC_3D(TEMP_3D(I,J,K),L)
          ELSE IF(AXIS == 2) THEN
            DENS_AXIS(J,L) = DENS_AXIS(J,L) + CONC_3D(TEMP_3D(I,J,K),L)
          ELSE IF(AXIS == 1) THEN
            DENS_AXIS(I,L) = DENS_AXIS(I,L) + CONC_3D(TEMP_3D(I,J,K),L)
          ELSE
            STOP "Wrong AXIS value"
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    DO K=1,N_SLICE(AXIS)
      DENS_AXIS(K,L) = DENS_AXIS(K,L)/REAL(PRODUCT(N_SLICE)/N_SLICE(AXIS)) ! this is for real density
    ENDDO
  ENDDO
  ! plotting text
!  MAX_DENS = MAXVAL(DENS_AXIS)
!  MIN_DENS = MINVAL(DENS_AXIS)
!  OPEN(UNIT=11,FILE=TRIM(OUTPUT_NAME1),STATUS='unknown')  
!  CLOSE(11)
!    WRITE(11,'(I8,1X,I3,1X,F13.6)') I, J, REAL(SUM_Z)/(T_SLICE(3,I)*BOX(2,I)*BOX(1,I))
!    WRITE(11,*)
!    WRITe(11,*)
END SUBROUTINE DENS_AXIS_1DR
       
SUBROUTINE BLOCK_AVG_2DR(MATX_2D,N_VAL,N_TIME,SIZE_BLOCK,BLOCK,N_BLOCK)
!! 02/21/2016 - Subroutine version
!! Basic Block Average (no checking) 
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: N_VAL, N_TIME, SIZE_BLOCK
  REAL, DIMENSION(1:N_VAL,1:N_TIME), INTENT(IN) :: MATX_2D
  REAL, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: BLOCK
  INTEGER :: I,J,K, IERR, I_BLOCK
  INTEGER, INTENT(OUT) :: N_BLOCK
  PRINT*, "==============================="
  PRINT*, "Starting Block averaging"
  IF(MOD(N_TIME,SIZE_BLOCK) /= 0) THEN
    PRINT*, "Ignore first ",MOD(N_TIME,SIZE_BLOCK)," steps"
  ENDIF
  N_BLOCK = (N_TIME-MOD(N_TIME,SIZE_BLOCK))/SIZE_BLOCK
  PRINT*, "Total time",N_TIME,", Size of Block, ",SIZE_BLOCK,", Number of Blocks,",N_BLOCK
  ALLOCATE(BLOCK(1:N_VAL,1:N_BLOCK),stat=IERR)
  IF (IERR /= 0) STOP "BLOCK in BLOCK_AVG_2DR: problem in attempt to allocate memory"
  BLOCK = 0.0D0
  DO K=1,N_VAL ! each index
    I_BLOCK = 0
    DO I=MOD(N_TIME,SIZE_BLOCK)+1,N_TIME-SIZE_BLOCK+1,SIZE_BLOCK! select only first time steps of each block except initial steps
      I_BLOCK = I_BLOCK + 1
      DO J=I,I+SIZE_BLOCK-1 ! do within each block
!        PRINT*,K,J,I,I_BLOCK
        BLOCK(K,I_BLOCK) = BLOCK(K,I_BLOCK) + MATX_2D(K,J)
      ENDDO
!      IF(K==22) PRINT*, "TEST?",BLOCK(K,I_BLOCK)
      BLOCK(K,I_BLOCK) = BLOCK(K,I_BLOCK)/REAL(SIZE_BLOCK)
    ENDDO
  ENDDO
  IF(N_BLOCK /= I_BLOCK) STOP "Memory size is not the same as loop times"
END SUBROUTINE BLOCK_AVG_2DR

SUBROUTINE DENS_3DR(POSI,T_SLICE,N_SLICE,N_ATOM,N_DATA,CONC_3D)
!! 02/20/2016 - subroutine version
!! 01/29/2016 - 3D version
!! calculate density in 3D which is the number of particle in volume of a cube
!! CONC_3D array should be 2D because I reduced the rank of POSI, 3D, to 2D.
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: N_ATOM, N_DATA
  REAL, DIMENSION(1:3,1:N_DATA), INTENT(IN) :: T_SLICE
  REAL, DIMENSION(1:3,1:N_ATOM,1:N_DATA), INTENT(IN) :: POSI
  INTEGER :: I,J,K,L,IERR,N_CUBE,I_CUBE
  REAL, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: CONC_3D
  INTEGER, DIMENSION(1:3), INTENT(IN) :: N_SLICE
  INTEGER, DIMENSION(1:4) :: T_INDEX
  N_CUBE = N_SLICE(1)*N_SLICE(2)*N_SLICE(3) ! total cubes
  ALLOCATE(CONC_3D(1:N_CUBE,1:N_DATA),stat=IERR) ! save memory using 2D arrary
  IF (IERR /= 0) STOP "CONC_3D in DENS_3DR: problem in attempt to allocate memory"
  CONC_3D = 0.0D0
  DO I=1,N_DATA
    DO J=1,N_ATOM
      DO K=1,3
        T_INDEX(K) = CEILING(POSI(K,J,I)/T_SLICE(K,I))
        IF(POSI(K,J,I) == 0.0D0) T_INDEX(K) = N_SLICE(K) ! zero position is the same as boundary of box
!        PRINT*, INDEX, COM(3,J,I)/T_SLICE, CEILING(COM(3,J,I)/T_SLICE), COM(3,J,I)/BOX(3,I),FLOOR(COM(3,J,I)/BOX(3,I))
        IF(T_INDEX(K) <= 0) THEN ! check something weird
          PRINT*, POSI(K,J,I),I,J
          STOP "index of slabs is less than 0 or zero"
          ELSE IF(T_INDEX(K) > N_SLICE(K)) THEN
            PRINT*, "WARNING: index of slabs is greater than MAX!"
            PRINT*, "In a",I,"th time step,",J,"th particle is placed at",POSI(K,J,I),"(nm) along",K,"-axis."
            PRINT*, "Because",K,"-axis edge of a cube is",T_SLICE(K,I),"at",I,"th time step,"
            PRINT*, "the cube index should be",T_INDEX(K),"but is beyond max.",N_SLICE(K)
            PRINT*, "If the particle could be at boundary, but not a problem."
            PRINT*, "If the same particle occurs errors during long time, DO DEBUGGING!"
            T_INDEX(K) = N_SLICE(K) !If POSI(K,J,I) is on boundary, precision matters and makes wrong T_INDEX.
        ENDIF
      ENDDO
      I_CUBE = (T_INDEX(3)-1)*N_SLICE(1)*N_SLICE(2) + (T_INDEX(2)-1)*N_SLICE(1) + T_INDEX(1)
!      IF(I_CUBE > 34*N_SLICE(1)*N_SLICE(2)) THEN
!        IF(I_CUBE < 35*N_SLICE(1)*N_SLICE(2)) THEN
!          PRINT*, "I_CUBE",I_CUBE
!        ENDIF
!      ENDIF
      IF(I_CUBE > N_CUBE) STOP "Wrong index of cube"
      CONC_3D(I_CUBE,I) = CONC_3D(I_CUBE,I) + 1.0D0
!        PRINT*, INDEX,I,CONC(INDEX,I),INT(CONC(INDEX,I))
    ENDDO
    DO L=1,N_CUBE
!      PRINT*,"Before",CONC_3D(L,I),(T_SLICE(1,I)*T_SLICE(2,I)*T_SLICE(3,I))
      CONC_3D(L,I) = CONC_3D(L,I)/(T_SLICE(1,I)*T_SLICE(2,I)*T_SLICE(3,I)) ! this is for real density. [number/(nm^3)] 
!      IF(L == 22) PRINT*, "TEST?",CONC_3D(L,I)
!      PRINT*, "After",CONC_3D(L,I)
    ENDDO
  ENDDO
END SUBROUTINE DENS_3DR

SUBROUTINE SLICE_BOX_3DR(D_XYZ,BOX,N_DATA,T_SLICE,N_SLICE)
!! 02/20/2016 - subroutine version
!! slice the box in 3D and get thickness of cubes
  IMPLICIT NONE
  REAL, DIMENSION(1:3), INTENT(IN) :: D_XYZ
  REAL, DIMENSION(1:3,1:N_DATA), INTENT(IN) :: BOX
  INTEGER, DIMENSION(1:3), INTENT(OUT) :: N_SLICE
  INTEGER, INTENT(IN) :: N_DATA
  INTEGER :: I,J,K, IERR
  REAL, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: T_SLICE
  REAL, DIMENSION(1:3) :: T_SLICE_AVG, T_SLICE_STD
  PRINT*, "==============================="
  PRINT*, "Slicing the simulation box"
  IF(D_XYZ(3) <= 0.1D0) THEN ! # of Slices from initial box_z length
    PRINT*, "WARNING: Thickness of slices, ",D_XYZ(3),"(nm) may be too small, < 1 A"
    PRINT*, "Because of single precision simulation, very small thickness is not recommended."
  ENDIF
  DO I=1,3 ! initially, set number of slices on axes from zero step
    N_SLICE(I) = CEILING(BOX(I,1)/D_XYZ(I))
  ENDDO
  PRINT*, "set ",N_SLICE(1:3)," = # of slices (x,y,z)"
  PRINT*, "set ",PRODUCT(N_SLICE)," = # of infinitesimal cubes"
  ALLOCATE(T_SLICE(1:3,1:N_DATA),stat=IERR)
  IF (IERR /= 0) STOP "T_SLICE in SLICE_BOX_3DR: problem in attempt to allocate memory"
  T_SLICE(I,J) = 0
  T_SLICE_AVG = 0.0D0
  T_SLICE_STD = 0.0D0
  DO I=1,N_DATA ! calculate thickness of slices
    DO J=1,3
      T_SLICE(J,I) = BOX(J,I)/REAL(N_SLICE(J))
      T_SLICE_AVG(J) = T_SLICE_AVG(J) + T_SLICE(J,I) 
      T_SLICE_STD(J) = T_SLICE_STD(J) + T_SLICE(J,I)**2 
!      PRINT*, I,J,T_SLICE(J,I),T_SLICE(J,I)**2,T_SLICE_STD(J)
    ENDDO
  ENDDO
  DO J=1,3
    T_SLICE_AVG(J) = T_SLICE_AVG(J)/N_DATA
    T_SLICE_STD(J) = T_SLICE_STD(J)/N_DATA - T_SLICE_AVG(J)**2
!    PRINT*,J,"TSLICE_STD",T_SLICE_STD(J),T_SLICE_AVG(J)
  ENDDO
  PRINT*, "Time AVG thickness of slices on x,y,z (nm),",T_SLICE_AVG(1:3)
  PRINT*, "Time STD thickness of slices on x,y,z (nm),",(SQRT(T_SLICE_STD(J)),J=1,3)
  PRINT*, "NaN means the length of the box is fixed in simulation."
  PRINT*, " But for a large box, it would be not Nan even though the lengh is fixed, but still very small."
END SUBROUTINE SLICE_BOX_3DR

SUBROUTINE POS_NOJUMP_PBC_3DR(POSI,BOX,N_ATOM,N_DATA)
!! 02/20/2016 - subroutine version
!! I use floor function to assign inside PBC
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: N_ATOM, N_DATA
  REAL, DIMENSION(1:3,1:N_ATOM,1:N_DATA), INTENT(INOUT) :: POSI
  REAL, DIMENSION(1:3,1:N_DATA), INTENT(IN) :: BOX
  INTEGER :: I,J,K
  DO I=1,N_DATA
    DO J=1,N_ATOM
      DO K=1,3
        POSI(K,J,I)=POSI(K,J,I)-BOX(K,I)*FLOOR(POSI(K,J,I)/BOX(K,I))
      ENDDO
    ENDDO
  ENDDO
END SUBROUTINE POS_NOJUMP_PBC_3DR

SUBROUTINE READ_XTC_R(INPUT_NAME,FILE_LEN,AR_DATA,BOX,N_ATOM,N_DATA)
!! 02/20/2016 - Subroutine version
!! 12/22/2015 - Update for semi-isotropic box(a,b,c)
!! I follow the instruction of James' code! More detail you can see:  
!! XDR Fortran Interface XTC Example Program with Wrapper
!! 2014 (c) James W. Barnett <jbarnet4@tulane.edu>
!! https://github.com/wesbarnett/
  USE xtc_mod, only: xtcfile ! Use the xdr interface
  IMPLICIT NONE ! IMPLICIT should be next to USE module
  TYPE(xtcfile) :: xtc ! Declare a variable of type xtcfile
  INTEGER, INTENT(IN) :: FILE_LEN
  CHARACTER(LEN=FILE_LEN), INTENT(IN) :: INPUT_NAME
  REAL, DIMENSION(:,:,:), ALLOCATABLE, INTENT(OUT) :: AR_DATA
  REAL, DIMENSION(:,:,:), ALLOCATABLE :: TEMP
  INTEGER :: MAX_DATA = 3000, ICR = 1000, IERR, I, J, K, L
  INTEGER, INTENT(OUT) :: N_DATA, N_ATOM ! Total time steps & total number of atoms
  REAL, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: BOX
  REAL, DIMENSION(:,:), ALLOCATABLE :: TEMP_BOX
  REAL, DIMENSION(1:3) :: BOX_INIT
! Initially, read Input file
  PRINT*, "==============================="
  PRINT*, "Reading Input File, ",TRIM(INPUT_NAME)
  CALL xtc%init(TRIM(INPUT_NAME)) ! Initialize it with the name of xtc file you want to read in.
  CALL xtc%read ! Read first information
  IF ( xtc%STAT == 0 ) THEN
    IF( xtc%prec /= 1000.00D0 ) THEN ! check precision of xtc file
      PRINT*, "xtc file precision is not 1000.0 (real)!, but ",xtc%prec
      PRINT*, "this program does only support real precision."
      STOP
    ENDIF
    PRINT*, "This program only support a CUBIC system."
    ! save original position 
    N_ATOM = xtc%NATOMS
    DO J=1,3
      BOX_INIT(J) = xtc%box(J,J) ! assume simulation box is rectangular parallelepiped.
    ENDDO
    PRINT*, "This initial BOX size (t=0) (",(BOX_INIT(J),J=1,3),")" 
    ELSE
      PRINT*, TRIM(INPUT_NAME)," has no trajectory. Check again!"
      STOP
  ENDIF
  ALLOCATE(AR_DATA(1:3,1:N_ATOM,1:MAX_DATA),stat=IERR)
  IF (IERR /= 0) STOP "DATA: problem in attempt to allocate memory"
  AR_DATA = 0.0D0
  ALLOCATE(BOX(1:3,1:MAX_DATA),stat=IERR)
  BOX = 0.0D0
  I = 1    
  DO ! save position and pbc cell
    IF(MOD(I,500) == 0) PRINT*,"Reading ",I,"th time"
    call xtc%read ! read data
    IF ( xtc%STAT /= 0 ) THEN !checking end of the xtc file
      PRINT*, "Reading COMPLETE"
      EXIT
    ENDIF
    IF(I > MAX_DATA) THEN ! check array memory for coordinates (AR_DATA) and box (BOX)
      PRINT*, "Memory increases from size of 3 x ",N_ATOM," x ",SIZE(AR_DATA)/(3*N_ATOM)
      ALLOCATE(TEMP(1:3,1:N_ATOM,1:MAX_DATA+ICR),stat=IERR)
      ALLOCATE(TEMP_BOX(1:3,1:MAX_DATA+ICR),stat=IERR)
      IF (IERR /= 0) STOP "reallocate_temp: problem in attempt to allocate memory"
      TEMP = 0.0D0
      TEMP_BOX = 0.0D0
      DO J=1,3
        DO K=1,N_ATOM
          DO L=1,MAX_DATA
            TEMP(J,K,L) = AR_DATA(J,K,L)
          ENDDO
        ENDDO
      ENDDO
      DO J=1,3
        DO L=1,MAX_DATA
          TEMP_BOX(J,L) = BOX(J,L)
        ENDDO
      ENDDO
      DEALLOCATE(AR_DATA)
      DEALLOCATE(BOX)
      ALLOCATE(AR_DATA(1:3,1:N_ATOM,1:MAX_DATA+ICR),stat=IERR)
      ALLOCATE(BOX(1:3,1:MAX_DATA+ICR),stat=IERR)
      IF (IERR /= 0) STOP "reallocate_array: problem in attempt to allocate memory"
      DO J=1,3
        DO K=1,N_ATOM
          DO L=1,MAX_DATA+ICR
            AR_DATA(J,K,L) = TEMP(J,K,L)
          ENDDO
        ENDDO
      ENDDO
      DO J=1,3
        DO L=1,MAX_DATA+ICR
          BOX(J,L) = TEMP_BOX(J,L)
        ENDDO
      ENDDO
      DEALLOCATE(TEMP)
      DEALLOCATE(TEMP_BOX)
      MAX_DATA = MAX_DATA+ICR
    ENDIF
    AR_DATA(:,:,I) = xtc%pos(:,:) ! save position into newly expanded matrix
    DO J=1,3
      BOX(J,I) = xtc%box(J,J)
    ENDDO
    I = I + 1
  ENDDO
  N_DATA = I-1
  PRINT*, "Total frame is ",N_DATA
  call xtc % close ! Close the file
! Just an example to show what was read in
! write(*,'(a,f12.6,a,i0)') " Time (ps): ", xtc % time, "  Step: ", xtc % STEP
! write(*,'(a,f12.6,a,i0)') " Precision: ", xtc % prec, "  No. Atoms: ", xtc % NATOMS
! write(*,'(3f9.3)') xtc % pos
! This is the same order as found in the GRO format fyi
! write(*,'(11f9.5)') xtc % box(1,1), xtc % box(2,2), xtc % box(3,3), &
!                     xtc % box(1,2), xtc % box(1,3), & 
!                     xtc % box(2,1), xtc % box(2,3), &
!                     xtc % box(3,1), xtc % box(3,2) 
END SUBROUTINE READ_XTC_R
