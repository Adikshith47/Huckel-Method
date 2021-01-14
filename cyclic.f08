program CHuckelMethod
  use gnufor2
  use huckelmethod
  implicit none

integer :: num = 5
real(kind=8) :: ring_hlgap

call ring_eigenval_and_eigenvector(num, ring_hlgap)
write(*,*) ring_hlgap

end program CHuckelMethod
