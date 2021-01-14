program HuckelMethod
  implicit none

  integer :: i, ldz, n, info, j
  integer, dimension(4) :: allst
  double precision, dimension(:), allocatable :: diag, subdiag, work
  double precision, dimension(:,:), allocatable :: z
  real :: a, b

  n = 6
  ldz = n

! Allocating memory to store the diagonal and subdiagonal elements of the matrix
  allocate(diag(n), stat=allst(1))
  allocate(subdiag(n-1), stat=allst(2))
  allocate(z(ldz, n), stat=allst(3))
  allocate(work(2*n-2), stat=allst(4))

! Assigning the corresponding energy values for the diagonal elements
  a = -11.40000
  b = -3.630000

  do i = 1, n
    diag(i) = a
  end do

  do i = 1, n-1
    subdiag(i) = b
  end do

!Implementing the sub routine to calculate the Eigenvalues and Eigenvectors of the molecule
  call dstev('V', n, diag, subdiag, z, ldz, work, info)

!Verifying the calculation was done properly
  if(info==0) then
    write(*,*) "The Eigenvalues are"
    write(*,*) (diag(i), i = 1,n)
    write(*,*) "The Eigenvectors are"
    do i = 1, n
      write(*,*) i,".", (z(j,i), j = 1,n)
    end do
  end if

end program HuckelMethod
