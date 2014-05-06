program main
  implicit none


  integer*4 :: g_duration
  complex*16 :: g_sigmax(2,2), g_sigmaz(2,2), g_a, g_b, g_twoi
  complex*16 :: g_psi(2), g_HAM(2,2), g_MTOP(2,2), g_MBOT(2,2), g_I(2,2)
  real*8 :: g_gamma = 1d0, g_alpha = 1d0, g_B0 = 1d0
  real*8 :: g_omega = 1d0, g_xi = 1d0, g_dt = 0.01 , g_t0 = 0d0, g_t = 0d0
    g_a = dcmplx(1d0,0d0) 
    g_b = dcmplx(1d0,0d0)
    g_twoi = dcmplx(0d0, 2d0)
    g_I(1,1) = dcmplx(1d0,0d0)
    g_I(1,2) = dcmplx(0d0,0d0)
    g_I(2,1) = dcmplx(0d0,0d0)
    g_I(2,2) = dcmplx(1d0,0d0)
    g_sigmax(1,1) = dcmplx(0d0,0d0)
    g_sigmax(1,2) = dcmplx(1d0,0d0)
    g_sigmax(2,1) = dcmplx(1d0,0d0)
    g_sigmax(2,2) = dcmplx(0d0,0d0)
    g_sigmaz(1,1) = dcmplx(1d0,0d0)
    g_sigmaz(1,2) = dcmplx(0d0,0d0)
    g_sigmaz(2,1) = dcmplx(0d0,0d0)
    g_sigmaz(2,2) = dcmplx(-1d0,0d0)
    g_duration = ceiling(1d0 / g_dt)
  g_HAM = -g_gamma * g_B0 * g_sigmaz + g_alpha * exp(-g_xi * (g_t - g_t0)**2) * sin(g_omega * g_t) * g_sigmax
  g_MTOP = g_twoi*g_I + g_HAM
  g_MBOT = g_twoi*g_I - g_HAM
  
  call invertComplex(g_MBOT)  

  call calcPsi()

Contains

subroutine invertComplex(A)
  complex*16, intent(inout) :: A(:,:)

  integer,parameter::M=300
  complex*16,allocatable,dimension(:,:)::A
  complex*16,allocatable,dimension(:)::WORK
  integer,allocatable,dimension(:)::IPIV
  integer i,j,info,error

  allocate(A(M,M),WORK(M),IPIV(M),stat=error)
  if (error.ne.0)then
    print *,"error:not enough memory"
    stop
  end if

  allocate(work(lwork)) 
  print*, "qae"
  call CGETRI(n, A, lda, ipiv, work, lwork, info)
  print*, info
  if(info /= 0) stop "error in call to CGETRI"

  print*, "wae"
  deallocate(work)

end subroutine invertComplex

subroutine calcPsi()
  complex*16 :: M(2,2)
  integer*4 :: i
  
  g_psi(1) = g_a / sqrt(abs(g_a)**2 + abs(g_b)**2)
  g_psi(2) = g_b / sqrt(abs(g_a)**2 + abs(g_b)**2)

  do i = 1, g_duration

    g_t = g_t + g_dt

    g_HAM = -g_gamma * g_B0 * g_sigmaz + g_alpha * exp(-g_xi * (g_t - g_t0)**2) * sin(g_omega * g_t) * g_sigmax
    g_MBOT = g_twoi - g_HAM
    call invertComplex(g_MBOT)

    M = MATMUL(g_MBOT, g_MTOP)
    g_psi = MATMUL(M,g_psi)
    g_MTOP = g_twoi + g_HAM
 
  end do
  
end subroutine calcPsi

end program