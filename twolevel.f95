program twolevel

   implicit none
   
   !- Constants
   
   !-- Simulation Variables
   complex*16 :: g_omega0
   complex*16 :: g_freq
   complex*16 :: g_period1
   complex*16 :: g_period2
   
   !-- Quantum Variables
   complex*16 :: g_sigmax(2,2), g_sigmaz(2,2), g_sigmay(2,2) !Pauli matrices
   complex*16 :: g_unity(2,2) ! Complex unity matrix
   
   complex*16 :: g_psi(2) ! The wavefunction
   
   complex*16 :: g_Ham(2,2) ! The Hamiltonian
   complex*16 :: sigma_plus(2,2), sigma_min(2,2) !jump matrices
   complex*16 :: jump_minplus(2,2), jump_plusmin(2,2) !matrix products plus-min & min-plus
   
   !-- Program Variables   
   integer*4  :: i,k, loop !counters
   real*8 :: g_dt, g_t0, g_t !g_dt = timestep, g_t0 = time of pulse, g_t = current time
   integer*4 :: g_duration !Duration of time evolution
   complex*16 :: g_HamNumer(2,2) ! Numerator of the Hamiltonian
   complex*16 :: g_HamDenom(2,2) ! Denominator of the Hamiltonian
   real*8 :: dp(2), dp_tot
   
   call setUpVariables()

   call initialize_files
   
   do i = 0,g_duration
      call calcPsi()
      call quantum_mc
   end do

   call close_files

contains
    subroutine quantum_mc
      real(8) :: eps(2) !dp = transition probability, eps = random number

      call random_number(eps)

      !Jump from ground to excited
      dp(1) = g_dt*realpart(dot_product(g_psi,matmul(jump_minplus,g_psi)))
      
      !Jump from excited to ground
      dp(2) = g_dt*realpart(dot_product(g_psi,matmul(jump_plusmin,g_psi)))
     
      !Total jump probability
      dp_tot = sum(dp)
      

      if (eps(1).lt.dp_tot) then
        if (eps(2).lt.dp(1)/dp_tot) then
          !Jump from ground state to excited state
          g_psi = matmul(sigma_plus,g_psi)/sqrt(dp(1)/g_dt)
          !DEBUG 
          !print *, abs(dot_product(g_psi,g_psi))**2
        else
          !Jump from excited state to ground state
          g_psi = matmul(sigma_plus,g_psi)/sqrt(dp(2)/g_dt)
        end if
      end if

      !Normalize wavefunction
      g_psi = g_psi/(sqrt(dot_product(g_psi,g_psi)))

      !Problem: probabilities are not complementary (don't add up to 1)
      !Solution?: normalize wave function after each iteration? (see 'ch. 4 Quantum Trajectories')
      write(11,*) dble(g_psi(1)*conjg(g_psi(1)))
      write(12,*) dble(g_psi(2)*conjg(g_psi(2)))

    end subroutine quantum_mc


   subroutine calcPsi
      complex*16 :: up(2) = (/dcmplx(0d0,0d0), dcmplx(1d0,0d0)/)
      complex*16 :: down(2) = (/dcmplx(1d0,0d0), dcmplx(0d0,0d0)/), ii=cmplx(0._8,1._8)

      !Define Hamiltonian
      g_Ham = g_omega0/2d0 * g_sigmaz + g_period1 * sin(g_freq * g_t) * g_sigmay &
         + g_period2 * cos(g_freq * g_t) * g_sigmax   
      
      !Define Crank-Nicholson 'numerator operator' at time t
      g_HamNumer = 1*g_unity - g_Ham * g_dt *II*0.5_8

      g_t = g_t + g_dt

      !Define Crank-Nicholson 'denominator operator' at time t+dt
      g_Ham = g_omega0/2d0 * g_sigmaz + g_period1 * sin(g_freq * g_t) * g_sigmay &
         + g_period2 * cos(g_freq * g_t) * g_sigmax

      g_HamDenom = 1*g_unity + g_Ham * g_dt * II*0.5_8
      call invertComplex(g_HamDenom)
      
      !Define time-evolution operator from t to t+dt
      g_Ham = matmul(g_HamDenom,g_HamNumer)
      
      !Update wave function
      g_psi = matmul(g_Ham,g_psi)
      
      !print *, dble(g_psi(1)*conjg(g_psi(1)))
      
   end subroutine calcPsi
   
   subroutine setUpVariables
      
      !- Define the Pauli Spin Matrices
      
      g_sigmax(1,:)  = (/ dcmplx(0d0,0d0), dcmplx(1d0,0d0) /)
      g_sigmax(2,:)  = (/ dcmplx(1d0,0d0), dcmplx(0d0,0d0) /)
      
      g_sigmay(1,:)  = (/ dcmplx(0d0,0d0), dcmplx(0d0,-1d0) /)
      g_sigmay(2,:)  = (/ dcmplx(0d0,1d0), dcmplx(0d0,0d0) /)
      
      g_sigmaz(1,:)  = (/ dcmplx(1d0,0d0), dcmplx(0d0,0d0) /)
      g_sigmaz(2,:)  = (/ dcmplx(0d0,0d0), dcmplx(-1d0,0d0) /)
      
      g_unity(1,:)   = (/ dcmplx(1d0,0d0), dcmplx(0d0,0d0) /)
      g_unity(2,:)   = (/ dcmplx(0d0,0d0), dcmplx(1d0,0d0) /)

      !- Define the simulation constants
      g_omega0       = dcmplx(1d0,0d0);
      g_freq         = dcmplx(1d0,0d0);
      g_period1      = dcmplx(1d0,0d0);
      g_period2      = dcmplx(1d0,0d0);
      
      !- Define the time-related variables
      g_dt  = 0.01d0
      g_t   = 0d0
      g_duration = 1000

      !Define Monte-carlo parameters
      sigma_min(1,:) = (/ dcmplx(0d0,0d0), dcmplx(1d0,0d0) /)
      sigma_min(2,:) = (/ dcmplx(0d0,0d0), dcmplx(0d0,0d0) /)
          
      sigma_plus(1,:) = (/ dcmplx(0d0,0d0), dcmplx(0d0,0d0) /)
      sigma_plus(2,:) = (/ dcmplx(1d0,0d0), dcmplx(0d0,0d0) /)

      jump_minplus = matmul(sigma_min,sigma_plus)
      jump_plusmin = matmul(sigma_plus,sigma_min)



      !- Set up the inital wavefunction (the up state)
      g_psi = (/ dcmplx(1d0,0d0), dcmplx(0d0,0d0) /)
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

   subroutine initialize_files

     open(11, file='prob_1.dat')
     open(12, file='prob_2.dat')

   end subroutine initialize_files
   
   subroutine close_files

     close(11)
     close(12)

   end subroutine close_files

end program twolevel
