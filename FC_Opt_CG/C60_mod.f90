module C60_mod

    use iso_fortran_env, only : REAL64, INT32
    implicit none

    integer(kind=INT32), parameter, private :: n = 60, m = 3
    integer(kind=INT32), parameter :: N_atoms = 60, n_neigh = 3

    integer(kind=INT32) :: cubic_neighbours(m,n), face_right(m,n), & 
                           next_on_face(m,n), prev_on_face(m,n),   &
                           dual_neighbours(m,n), next_on_tri(m,n), &
                           rspi(12,1)
    integer(kind=INT32) :: pentagons(5,12), hexagons(6,20)
    integer(kind=INT32) :: triangles(m,n)

    real(kind=REAL64)    :: points_start(m,n), points_opt(m,n), tutte_layout(2,n)

end module
