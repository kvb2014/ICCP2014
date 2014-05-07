program twolevel

   implicit none
   
   !- Constants
   
   !-- Simulation Variables
   complex*16 :: g_omega0
   complex*16 :: g_freq
   complex*16 :: g_period1
   complex*16 :: g_period2
   
   !-- Quantum Variables
   complex*16 :: g_sigmax(3,3), g_sigmaz(3,3), g_sigmay(3,3) !Pauli matrices
   complex*16 :: g_unity(3,3) ! Complex unity matrix
   
   complex*16 :: g_psi(3) ! The wavefunction
   
   complex*16 :: g_Ham(3,3) ! The Hamiltonian
   
   !-- Program Variables   
   integer*4  :: i,k, loop !counters
   real*8 :: g_dt, g_t0, g_t !g_dt = timestep, g_t0 = time of pulse, g_t = current time
   integer*4 :: g_duration !Duration of time evolution
   complex*16 :: g_HamNumer(3,3) ! Numerator of the Hamiltonian
   complex*16 :: g_HamDenom(3,3) ! Denominator of the Hamiltonian
   
   call setUpVariables()
   do i = 0,g_duration
      call calcPsi()
   end do
   
contains

   subroutine calcPsi
      complex*16 :: up(2) = (/dcmplx(0d0,0d0), dcmplx(1d0,0d0)/)
      complex*16 :: down(2) = (/dcmplx(1d0,0d0), dcmplx(0d0,0d0)/), ii=cmplx(0._8,1._8)

      g_Ham = g_omega0/2d0 * g_sigmaz + g_period1 * sin(g_freq * g_t) * g_sigmay &
         + g_period2 * cos(g_freq * g_t) * g_sigmax   
   
      g_HamNumer = 1*g_unity - g_Ham * g_dt * II * 0.5_8

      g_t = g_t + g_dt

      g_Ham = g_omega0/2d0 * g_sigmaz + g_period1 * sin(g_freq * g_t) * g_sigmay &
         + g_period2 * cos(g_freq * g_t) * g_sigmax

      g_HamDenom = 1*g_unity + g_Ham * g_dt * II * 0.5_8
      call invertComplex(g_HamDenom)
      
      g_Ham = matmul(g_HamDenom,g_HamNumer)
      
      g_psi = matmul(g_Ham,g_psi)
      
      print *, dble(g_psi(1)*conjg(g_psi(1)))
      
   end subroutine calcPsi
   
   subroutine setUpVariables
      
      !- Define the Pauli Spin Matrices
      
      g_sigmax(1,:)  = (/ dcmplx(0d0,0d0), dcmplx(1d0,0d0), dcmplx(0d0,0d0)/)
      g_sigmax(2,:)  = (/ dcmplx(1d0,0d0), dcmplx(0d0,0d0), dcmplx(1d0,0d0)/)
      g_sigmax(3,:)  = (/ dcmplx(0d0,0d0), dcmplx(1d0,0d0), dcmplx(0d0,0d0)/)
      
      g_sigmay(1,:)  = (/ dcmplx(0d0,0d0), dcmplx(0d0,-1d0), dcmplx(0d0,0d0)/)
      g_sigmay(2,:)  = (/ dcmplx(0d0,1d0), dcmplx(0d0,0d0), dcmplx(0d0,-1d0)/)
      g_sigmay(3,:)  = (/ dcmplx(0d0,0d0), dcmplx(0d0,1d0), dcmplx(0d0,0d0)/)
      
      g_sigmaz(1,:)  = (/ dcmplx(0d0,0d0), dcmplx(0d0,0d0), dcmplx(0d0,0d0)/)
      g_sigmaz(2,:)  = (/ dcmplx(0d0,0d0), dcmplx(1d0,0d0), dcmplx(0d0,0d0)/)
      g_sigmaz(3,:)  = (/ dcmplx(0d0,0d0), dcmplx(0d0,0d0), dcmplx(-2d0,0d0)/)
     
      g_unity(1,:)   = (/ dcmplx(1d0,0d0), dcmplx(0d0,0d0), dcmplx(0d0,0d0)/)
      g_unity(2,:)   = (/ dcmplx(0d0,0d0), dcmplx(1d0,0d0), dcmplx(0d0,0d0)/)
      g_unity(2,:)   = (/ dcmplx(0d0,0d0), dcmplx(0d0,0d0), dcmplx(1d0,0d0)/)
      
      !- Define the simulation constants
      g_omega0       = dcmplx(1d0,0d0);
      g_freq         = dcmplx(1d0,0d0);
      g_period1      = dcmplx(1d0,0d0);
      g_period2      = dcmplx(1d0,0d0);
      
      !- Define the time-related variables
      g_dt  = .01d0
      g_t   = 0d0
      g_duration = 1000
      
      !- Set up the inital wavefunction (the up state)
      g_psi = (/ dcmplx(1d0,0d0), dcmplx(0d0,0d0), dcmplx(0d0,0d0)/)
      !- Set up the inital wavefunction (the down state)
      !g_psi = (/ dcmplx(1d0,0d0), dcmplx(0d0,0d0) /)
      
   end subroutine setUpVariables
   
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

     call ZGETRI(M,A,M,IPIV,WORK,M,info)

     deallocate(IPIV, WORK, stat=error)
     if (error.ne.0)then
       print *,"error:fail to release"
       stop
     end if

   end subroutine invertComplex
   

end program twolevel
