module geom_parameters_mod

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! DESCRIPTION
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use iso_fortran_env, only : REAL64, INT32

    type geom_pars_t
        real(kind=REAL64), allocatable :: R_0(:,:)
        real(kind=REAL64), allocatable :: ang_0(:,:)
    end type geom_pars_t


 ! PARAMETERS
    real(kind=REAL64), parameter, private :: R_constants(3) = (/1.479,1.458,1.401/)

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
        !integer(kind=INT32)    :: left, right
        integer(kind=INT32)    :: left(N_atoms), right(N_atoms)
        !integer(kind=INT32), allocatable :: face_left(:,:)

        allocate (geom_pars%R_0(N_atoms,n_neigh))
        !allocate (face_left, mold=face_right)

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
        DO i = 1, n_neigh
            right = face_right(:,i) - 5
            left  = face_right(:,modulo(i-2,n_neigh)+1) - 5
            geom_pars%R_0(:,i) = R_constants(left+right+1)
        END DO

    end subroutine
end module
