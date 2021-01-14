program homolumogap
  use huckelmethod
  use gnufor2
  implicit none

  integer :: j
  integer, parameter :: n = 15
  real(kind=8), dimension(2:n) :: chain_hlgap, num_atoms
  real(kind=8) :: tempgap
  real(kind=8), dimension(:), allocatable :: eval
  real(kind=8), dimension(:,:), allocatable :: evec

  do j = 2, 15
    call chain_eigenval_and_eigenvector(j,eval,evec,tempgap)
    num_atoms(j) = j
    chain_hlgap(j) = tempgap
  end do

call plot_1(num_atoms, chain_hlgap, style = ' 7.')

end program homolumogap
