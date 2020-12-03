module geom_parameters_mod

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! DESCRIPTION
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use iso_fortran_env, only : REAL64, INT32

 ! DATA TYPES 
    type geom_pars_t
        real(kind=REAL64), allocatable :: R_0(:,:)
        real(kind=REAL64), allocatable :: f_R(:,:)
        real(kind=REAL64), allocatable :: ang0_p(:,:)
        real(kind=REAL64), allocatable :: ang0_m(:,:)
        real(kind=REAL64), allocatable :: f_ang(:,:)
        real(kind=REAL64), allocatable :: f_ang_m(:,:)
        real(kind=REAL64), allocatable :: dih_0(:,:)
        real(kind=REAL64), allocatable :: f_dih(:,:)
    end type geom_pars_t


 ! PARAMETERS
    real(kind=REAL64), parameter, private :: pi = dacos(-1.d0)
    real(kind=REAL64), parameter, private :: R_constants(3) = (/1.479,1.458,1.401/)
    real(kind=REAL64), parameter, private :: bond_angles(2) = (/108.0*pi/180.0,120.0*pi/180.0/)
    real(kind=REAL64), parameter, private :: dih_constants(2,2,2) = reshape(dcos((/0.652358d0,0.509674d0,0.509674d0, &
                                             0.345123d0,0.615841d0,0.417884d0,0.417884d0,0.0d0/)), shape=(/2,2,2/))
    real(kind=REAL64), parameter, private :: f_const1(2) = (/100.0,100.0/), f_const2(3) = (/260.0,390.0,450.0/),     &
                                             f_const3(4) = (/35.0,65.0,85.0,270.0/)
    real(kind=REAL64), parameter, private :: f_const1_seminario(2) = (/207.924,216.787/),                            &
                                             f_const2_seminario(3) = (/260.0,353.377,518.992/),                      &
                                             f_const3_seminario(4) = (/35.0,65.0,3.772,270.0/) 


contains

    subroutine geom_parameters(face_right,N_atoms,neighbours,n_neigh,seminario,geom_pars)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! DESCRIPTION
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        implicit none

    ! FORMAL ARGUMENTS
        integer(kind=INT32), intent(in)   :: N_atoms
        integer(kind=INT32), intent(in)   :: n_neigh
        integer(kind=INT32), intent(in)   :: face_right(N_atoms,n_neigh)
        integer(kind=INT32), intent(in)   :: neighbours(N_atoms,n_neigh)
        logical,             intent(in)   :: seminario
        type(geom_pars_t),   intent(out)  :: geom_pars

    ! VARIABLES
        integer(kind=INT32)    :: i, j
        integer(kind=INT32)    :: left1, right1, left2
        integer(kind=INT32)    :: left(N_atoms), right(N_atoms)
        !integer(kind=INT32), allocatable :: face_left(:,:)

        !allocate (face_left, mold=face_right)
        allocate (geom_pars%R_0(N_atoms,n_neigh))
        allocate (geom_pars%f_R(N_atoms,n_neigh))
        allocate (geom_pars%ang0_p(N_atoms,n_neigh))
        allocate (geom_pars%ang0_m(N_atoms,n_neigh))
        allocate (geom_pars%f_ang(N_atoms,n_neigh))
        allocate (geom_pars%f_ang_m(N_atoms,n_neigh))
        allocate (geom_pars%dih_0(N_atoms,n_neigh))
        allocate (geom_pars%f_dih(N_atoms,n_neigh))

!        ! cshift version
!        face_left = cshift(face_right, -1, 2)
!        DO i = 1, n_neigh
!            geom_pars%R_0(:,i) = R_constants(face_right(:,i)+face_left(:,i)-10 +1) 
!        END DO
!        ! double loop version
!        DO i = 1, n_neigh
!            DO j = 1, N_atoms
!                right = face_right(j,i) - 5
!                left  = face_right(j,modulo(i-2,n_neigh)+1) - 5
!                geom_pars%R_0(j,i) = R_constants(left+right+1)
!            END DO
!        END DO

        ! fill equil. bond distances and angles into array
        DO i = 1, n_neigh
            right = face_right(:,i) - 5
            left  = face_right(:,modulo(i-2,n_neigh)+1) - 5        !py roll +1
            geom_pars%R_0(:,i)     = R_constants(left+right+1)
            geom_pars%f_R(:,i)     = f_const2(left+right+1)
            geom_pars%ang0_p(:,i)  = dcos(bond_angles(right+1))
            geom_pars%f_ang(:,i)   = f_const1(right+1)
            !todo: fang_p is equal to fang
            geom_pars%f_ang_m(:,i) = f_const1(left+1)
            geom_pars%ang0_m(:,i)  = dcos(bond_angles(left+1))
        END DO

        ! fill equil. dih. constants into array (todo: look into getting rid of the inner loop)
        DO i = 1, n_neigh
            DO j = 1, N_atoms
                right1 = face_right(j,i) - 5
                left2  = face_right(j,mod(i+n_neigh,n_neigh)+1) - 5    !py roll -1
                left1  = face_right(j,modulo(i-2,n_neigh)+1) - 5       !py roll -2 (same as +1)
                geom_pars%dih_0(j,i) = dih_constants(left1+1,left2+1,right1+1)
                geom_pars%f_dih(j,i) = f_const3(right1+left1+left2+1)
            END DO
        END DO
       
!        write(*,*) 'dih constants order as in python'
!        do i = 1,2
!            write(*,*) 'i= ', i
!            do j = 1, 2
!                write(*,*) dih_constants(:,j,i)
!            end do
!        end do

        write(*,*) 'f dih'
        do i = 1, N_atoms
            write(*,*) geom_pars%f_dih(i,:)
        end do

    end subroutine
end module
