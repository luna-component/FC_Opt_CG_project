module C60_mod

    use iso_fortran_env, only : REAL64, INT32
    implicit none

    integer(kind=INT32), parameter, private :: n = 60, m = 3
    integer(kind=INT32), parameter :: N_atoms = 60, n_neigh = 3

    integer(kind=INT32) :: cubic_neighbours(n,m), face_right(n,m), & 
                           next_on_face(n,m), prev_on_face(n,m),   &
                           dual_neighbours(32,6), next_on_tri(32,6), &
                           rspi(1,12)
    integer(kind=INT32) :: pentagons(12,5), hexagons(20,6)
    integer(kind=INT32) :: triangles(n,m)

    real(kind=REAL64)    :: points_start(n,m), points_opt(n,m), tutte_layout(n,2)

end module
