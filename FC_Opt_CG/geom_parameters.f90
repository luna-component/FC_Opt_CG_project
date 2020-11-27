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
    real(kind=REAL64), parameter, private :: R_constants(3) = (/1.479,1.458, 1.401/)

contains

    subroutine geom_parameters(face_right,N_atoms,neighbours,n_neigh,seminario, geom_pars)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! DESCRIPTION
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        implicit none

    ! FORMAL ARGUMENTS
        integer(kind=INT32), intent(in)   :: N_atoms
        integer(kind=INT32), intent(in)   :: n_neigh
        integer(kind=INT32), intent(in)   :: face_right(n_neigh,N_atoms)
        integer(kind=INT32), intent(in)   :: neighbours(n_neigh,N_atoms)
        logical,             intent(in)   :: seminario
        type(geom_pars_t),   intent(out)  :: geom_pars

    ! VARIABLES
        integer(kind=INT32)    :: i, j
        integer(kind=INT32)    :: left, right
        integer(kind=INT32), allocatable :: face_left(:,:)

        allocate (geom_pars%R_0(n_neigh,N_atoms))
        allocate (face_left, mold=face_right)

        face_left = cshift(face_right, -1, 1)
        DO i = 1, n_neigh
            geom_pars%R_0(i,:) = R_constants(face_right(i,:)+face_left(i,:)-10 +1) 
        END DO
!        ! fill the array of R_0 constants
!        DO j = 1, N_atoms
!            DO i = 1, n_neigh
!                right = face_right(i,j) - 5
!                left  = face_right(modulo(i-2,n_neigh)+1,j) - 5
!                geom_pars%R_0(i,j) = R_constants(left+right+1)
!            END DO
!        END DO

    end subroutine
end module
