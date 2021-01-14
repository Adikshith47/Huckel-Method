module HuckelMethod
  implicit none

  real(kind=8) :: a = -11.40000000, b = -3.63000000

contains
  subroutine chain_eigenval_and_eigenvector(num, diag, z, gap)
    integer, intent(in) :: num
    integer :: i, ldz, info
    integer, dimension(4) :: allst
    double precision, dimension(:), allocatable :: subdiag, work
    double precision, dimension(:), allocatable, intent(out) :: diag
    double precision, dimension(:,:), allocatable, intent(out) ::  z
    real(kind=8), intent(out) :: gap

    ldz = num

    allocate(diag(num), stat=allst(1))
    allocate(subdiag(num-1), stat=allst(2))
    allocate(work(2*num-2), stat=allst(3))
    allocate(z(ldz, num), stat=allst(4))

    diag = (/(a, i=1,num)/)
    subdiag = (/(b, i=1,num-1)/)

    call dstev('V', num, diag, subdiag, z, ldz, work, info)

    if(MOD(num,2)==0) then
      gap = -1.0*diag(num/2) + diag((num/2)+1)
    else
      gap = -1.0*diag((num+1)/2) + diag(((num+1)/2)+1)
    end if

    return
  end subroutine chain_eigenval_and_eigenvector

  subroutine ring_eigenval_and_eigenvector(num, gap)
    integer, intent(in) :: num
    real(kind=8), intent(out) :: gap
    double precision, dimension(:,:), allocatable :: mat
    double precision, dimension(:), allocatable :: w, work
    integer, dimension(3) :: allst
    integer :: i, j, lwork, lda, info

    lda = num
    lwork = -1

    allocate(w(num), stat=allst(2))
    allocate(work(MAX(1, lwork)), stat=allst(3))
    allocate(mat(lda, num), stat=allst(1))

    do i = 1, num
      do j = 1, num
        if(i==j) then
          mat(i,j) = a
        else if(i-j==1.or.j-i==1.or.i-j==(num-1).or.j-i==(num-1)) then
          mat(i,j) = b
        end if
      end do
    end do


    call dsyev('V', 'U', num, mat, lda, w, work, lwork, info)

    if(info==0) then
      lwork = nint(work(1))
    end if

    call dsyev('V', 'U', num, mat, lda, w, work, lwork, info)

    if(MOD(num,2)==0) then
      gap = -1.0*w(num/2) + w((num/2)+1)
    else
      gap = -1.0*w((num+1)/2) + w(((num+1)/2)+1)
    end if

    deallocate(mat)
    deallocate(w)
    deallocate(work)

    return
  end subroutine ring_eigenval_and_eigenvector

  subroutine chain_laplacian_matrix(num, diag, z)

    integer, intent(in) :: num
    integer :: i, ldz, info
    integer, dimension(4) :: allst
    double precision, dimension(:), allocatable :: subdiag, work
    double precision, dimension(:), allocatable, intent(out) :: diag
    double precision, dimension(:,:), allocatable, intent(out) ::  z

    ldz = num

    allocate(diag(num), stat=allst(1))
    allocate(subdiag(num-1), stat=allst(2))
    allocate(work(2*num-2), stat=allst(3))
    allocate(z(ldz, num), stat=allst(4))

    do i = 2, num-1
      diag(i) = 2
    end do
    diag(1) = 1
    diag(num) = 1
    subdiag = (/(-1, i=1,num-1)/)

    call dstev('V', num, diag, subdiag, z, ldz, work, info)
    return
  end subroutine chain_laplacian_matrix

  subroutine ring_laplacian_matrix(num, mat)
    integer, intent(in) :: num
    double precision, dimension(:,:), allocatable, intent(out) :: mat
    double precision, dimension(:), allocatable :: w, work
    integer, dimension(3) :: allst
    integer :: i, j, lwork, lda, info

    lda = num
    lwork = -1

    allocate(w(num), stat=allst(2))
    allocate(work(MAX(1, lwork)), stat=allst(3))
    allocate(mat(lda, num), stat=allst(1))

    do i = 1, num
      do j = 1, num
        if(i==j) then
          mat(i,j) = 2
        else if(i-j==1.or.j-i==1.or.i-j==(num-1).or.j-i==(num-1)) then
          mat(i,j) = -1
        end if
      end do
    end do


    call dsyev('V', 'U', num, mat, lda, w, work, lwork, info)

    if(info==0) then
      lwork = nint(work(1))
    end if

    call dsyev('V', 'U', num, mat, lda, w, work, lwork, info)

    return
  end subroutine ring_laplacian_matrix


end module HuckelMethod
