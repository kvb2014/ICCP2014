program main
  implicit none

  integer*4  :: i,k, loop !counters
  integer*4 :: g_duration !Duration of time evolution
  complex*16 :: g_sigmax(2,2), g_sigmaz(2,2) !x and z Pauli matrices
  complex*16 :: g_a, g_b !a, b = linear combination coefficients of initial condition
  complex*16 :: g_twoi !shortkey for 2*i
  complex*16 :: g_psi(2) !wave function vector
  complex*16 :: g_HAM(2,2) !2x2 Hamiltonian matrix
  complex*16 :: g_MTOP(2,2), g_MBOT(2,2) !MTOP = (2i+HAM), MBOT = inverse(2i-HAM)
  complex*16 :: g_I(2,2) !complex-valued identity matrix
  real*8 :: g_gamma !constant for effect of magnetic field on energy
  real*8 :: g_B0 !constant magnetic field strength 
  real*8 :: g_alpha !Driving force coupling constant
  real*8 :: g_omega !Driving frequency
  real*8 :: g_xi !pulse width constant  
  real*8 :: g_dt, g_t0, g_t !g_dt = timestep, g_t0 = time of pulse, g_t = current time
  integer :: CL_narg
  character(len=256) :: CL_input, g_constants_file

  call check_commandline_arguments

  call setup_variables

  call calcPsi()

Contains

subroutine check_commandline_arguments
  CL_narg = command_argument_count() !count number of input arguments

  if (CL_narg.gt.0) then
    do loop = 1,CL_narg
      call get_command_argument(loop,CL_input) !get the 'loop'th input argument
      select case(adjustl(CL_input))
        case ("-i","-input")
          call get_command_argument(loop+1,CL_input) !makes CL_input the name of the constants file
          g_constants_file = CL_input
      end select
    end do 
  end if  

end subroutine check_commandline_arguments

subroutine setup_variables
  
  !Read in the constants
  open(11, file=g_constants_file)
  read(11,*) g_gamma
  read(11,*) g_B0
  read(11,*) g_alpha
  read(11,*) g_omega
  read(11,*) g_xi
  read(11,*) g_t0
  read(11,*) g_dt
  read(11,*) g_t
  g_duration = ceiling(1d0 / g_dt)

  !Define initial condition linear coefficients
    g_a = dcmplx(1d0,0d0) 
    g_b = dcmplx(1d0,0d0)

  !Define the matrices
    g_twoi = dcmplx(0d0, 2d0) !shortcut for 2*i

    do i=1,2
      do k=1,2
        if (i.eq.k) then
        g_I(i,k) = dcmplx(1d0,0d0)
        g_sigmax(i,k) = dcmplx(0d0,0d0)
        g_sigmaz(i,k) = dcmplx((-1d0)**(i+1),0d0)
        else 
        g_I(i,k) = dcmplx(0d0,0d0)
        g_sigmax = dcmplx(1d0,0d0)
        end if
      end do
    end do

  g_HAM = -g_gamma * g_B0 * g_sigmaz & 
  + g_alpha * exp(-g_xi * (g_t - g_t0)**2) * sin(g_omega * g_t) * g_sigmax
  
  g_MTOP = g_twoi*g_I + g_HAM

  g_MBOT = g_twoi*g_I - g_HAM
  call invertComplex(g_MBOT)
  close(11)

end subroutine setup_variables

subroutine invertComplex(A)
  complex*16, intent(inout) :: A(:,:)

  complex*16, allocatable, dimension(:) :: WORK
  integer, allocatable, dimension(:) :: IPIV
  integer info, error, M
  
  M = size(A,1)

  allocate(WORK(M), IPIV(M), stat=error)
  if (error.ne.0)then
    print *,"error:not enough memory"
    stop
  end if

  call ZGETRF(M,M,A,M,IPIV,info)
  if(info .eq. 0) then
    write(*,*)"succeded"
  else
   write(*,*)"failed"
  end if

  call ZGETRI(M,A,M,IPIV,WORK,M,info)
  if(info .eq. 0) then
    write(*,*)"succeded"
  else
   write(*,*)"failed"
  end if
  deallocate(IPIV, WORK, stat=error)
  if (error.ne.0)then
    print *,"error:fail to release"
    stop
  end if

end subroutine invertComplex

subroutine calcPsi()
  complex*16 :: M(2,2)
  integer*4 :: i

  g_psi(1) = g_a / sqrt(abs(g_a)**2 + abs(g_b)**2)
  g_psi(2) = g_b / sqrt(abs(g_a)**2 + abs(g_b)**2)
  
  open(unit = 20,file = 'data.txt')
  !open(unit = 21, file = 'inverses.txt')
  !open(unit = 22, file = 'M.txt')
  !open(unit = 23, file = 'g_psi')
  !open(unit = 24, file = 'g_MTOP')

  do i = 1, g_duration
    write(20, *) abs(dot_product(g_psi,g_psi))
    g_t = g_t + g_dt

    !write(21, *) g_MBOT(1,1), ' ', g_MBOT(1,2), ' ', g_MBOT(2,1), ' ', g_MBOT(2,2)
    !write(22, *) M(1,1), ' ', M(1,2), ' ', M(2,1), ' ', M(2,2)
    !write(23, *) g_psi(1), ' ', g_psi(2)
    !write(24, *) g_MTOP(1,1), ' ', g_MTOP(1,2), ' ', g_MTOP(2,1), ' ', g_MTOP(2,2)

    g_HAM = -g_gamma * g_B0 * g_sigmaz + g_alpha * exp(-g_xi * (g_t - g_t0)**2) * sin(g_omega * g_t) * g_sigmax
    g_MBOT = g_twoi * g_I - g_HAM
    call invertComplex(g_MBOT)

    M = MATMUL(g_MBOT, g_MTOP)
    g_psi = MATMUL(M,g_psi)
    g_MTOP = g_twoi * g_I + g_HAM
 
  end do

  close(20)
  !close(21)
  !close(22)
  !close(23)
  !close(24)

end subroutine calcPsi

end program
