program threelevel

   implicit none
   
   !- Constants
   
   !-- Simulation Variables
   complex*16 :: g_omega0
   complex*16 :: g_freq
   complex*16 :: g_period1
   complex*16 :: g_period2
   complex*16 :: g_delta
   
   !-- Quantum Variables
   complex*16 :: g_sigmax(3,3), g_sigmaz(3,3), g_sigmay(3,3) !Pauli matrices
   complex*16 :: g_unity(3,3) ! Complex unity matrix
   complex*16 :: g_proj2(3,3) ! The projection matrix voor level 2
   complex*16 :: sigma_0to1(3,3), sigma_1to0(3,3), sigma_1to2(3,3), sigma_2to1(3,3) !Jump matrices
   complex*16 :: jump_0to1(3,3), jump_1to0(3,3), jump_1to2(3,3), jump_2to1(3,3) ! conjugate matrix products of jump matrices
   complex*16, allocatable :: g_psi(:,:) ! The wavefunction
   
   complex*16 :: g_Ham(3,3) ! The Hamiltonian
   
   !-- Program Variables   
   integer*4  :: i,k, loop !counters
   real*8 :: g_dt, g_t0, g_t !g_dt = timestep, g_t0 = time of pulse, g_t = current time
   integer*4 :: g_duration !Duration of time evolution
   complex*16 :: g_HamNumer(3,3) ! Numerator of the Hamiltonian
   complex*16 :: g_HamDenom(3,3) ! Denominator of the Hamiltonian
   integer*4 :: g_ensemblesize !Size of the MC ensemble
   complex*16, allocatable :: g_psiAvg(:,:) !Ensemble average of psi
   
   call initialize_files

   call setUpVariables()

   do i = 0,g_duration
      call calcPsi()
      call quantum_mc(i)
   end do

   call write_files

   call close_files
   
contains

  subroutine quantum_mc(iteration)
    integer*4, intent(in) :: iteration !Iteration point where we're at
    real(8) :: eps(2), dp(4),dp_tot !dp = transition probability, eps = random number

     do k = 1, g_ensemblesize 
          call random_number(eps)

          !Compute transition probabilities
          dp(1) = g_dt * realpart(dot_product(g_psi(:,k),matmul(jump_0to1,g_psi(:,k))))
          
          dp(2) = g_dt * realpart(dot_product(g_psi(:,k),matmul(jump_1to0,g_psi(:,k))))

          dp(3) = g_dt * realpart(dot_product(g_psi(:,k),matmul(jump_1to2,g_psi(:,k))))

          dp(4) = g_dt * realpart(dot_product(g_psi(:,k),matmul(jump_2to1,g_psi(:,k))))

          dp_tot = sum(dp)

          !print *, dp_tot

          if  (eps(1).lt.dp_tot) then
              if ( eps(2).lt.(dp(1)/dp_tot) ) then
              g_psi(:,k) = matmul(sigma_0to1,g_psi(:,k))
              !print *, '0 to 1'
              elseif ( (eps(2).lt.(dp(2)/dp_tot) ) ) then
              g_psi(:,k) = matmul(sigma_1to0,g_psi(:,k))
              !print *, '1 to 0'
              elseif ( (eps(2).lt.(dp(3)/dp_tot) ) ) then
              g_psi(:,k) = matmul(sigma_1to2,g_psi(:,k))
              !print *, '1 to 2'              
              elseif ( (eps(2).lt.1d0) ) then              
              g_psi (:,k) = matmul(sigma_2to1,g_psi(:,k))
              !print *, '2 to 1'
            end if
          end if

          !Normalize wavefunction
          g_psi(:,k) = g_psi(:,k)/(sqrt(dot_product(g_psi(:,k),g_psi(:,k))))   
      end do

      !Average the wavefunction evolution
      do k=1,3
        g_psiAvg(k,iteration) = sum(g_psi(k,:))
      end do

      end subroutine quantum_mc

   subroutine calcPsi
      complex*16 :: ground(3) = (/dcmplx(1d0,0d0), dcmplx(0d0,0d0), dcmplx(0d0,0d0)/)
      complex*16 :: excited(3) = (/dcmplx(0d0,0d0), dcmplx(1d0,0d0), dcmplx(0d0,0d0)/)
      complex*16 :: leak(3) = (/dcmplx(0d0,0d0), dcmplx(0d0,0d0), dcmplx(1d0,0d0)/)
      complex*16 :: ii=dcmplx(0._8,1._8)

      call calcHamiltonian(g_Ham) 
   
      g_HamNumer = 1*g_unity - g_Ham * g_dt * II * 0.5_8

      g_t = g_t + g_dt

      call calcHamiltonian(g_Ham)  

      g_HamDenom = 1*g_unity + g_Ham * g_dt * II * 0.5_8
      call invertComplex(g_HamDenom)
      
      g_Ham = matmul(g_HamDenom,g_HamNumer)
      
      do k=1, g_ensemblesize
        g_psi(:,k) = matmul(g_Ham,g_psi(:,k))
      end do
      
      ! Calculate the Hamiltonian for the Energy
      call calcHamiltonian(g_Ham) 
      
      !print *, abs(dot_product(ground,g_psi))**2d0, &
      !         abs(dot_product(excited,g_psi))**2d0, &
      !         abs(dot_product(leak,g_psi))**2d0, &
      !         dot_product(g_psi,matmul(g_Ham,g_psi))
      
   end subroutine calcPsi
   
   subroutine calcHamiltonian(H)
      complex*16, intent(inout) :: H(:,:)
      complex*16 :: H0(3,3), H_perturb(3,3)
      
      H0 = g_omega0/2d0 * g_sigmaz &
               + g_proj2 * g_delta &
               + g_period1 * sin(g_freq * g_t) * g_sigmay &
               + g_period2 * cos(g_freq * g_t) * g_sigmax

      !Perturbing jump elements
      H_perturb = -0.5d0*dcmplx(0d0,1d0)*(jump_0to1+jump_1to0+jump_1to2+jump_2to1)

      !Effective Hamiltonian
      H = H0 + H_perturb
          
   end subroutine calcHamiltonian
   
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
      g_unity(3,:)   = (/ dcmplx(0d0,0d0), dcmplx(0d0,0d0), dcmplx(1d0,0d0)/)
      
      !- Define the projection matrices
      g_proj2(1,:)  = (/ dcmplx(0d0,0d0), dcmplx(0d0,0d0), dcmplx(0d0,0d0)/)
      g_proj2(2,:)  = (/ dcmplx(0d0,0d0), dcmplx(0d0,0d0), dcmplx(0d0,0d0)/)
      g_proj2(3,:)  = (/ dcmplx(0d0,0d0), dcmplx(0d0,0d0), dcmplx(1d0,0d0)/)
      
      !- Define the Monte Carlo parameters
      sigma_0to1(1,:) = (/ dcmplx(0d0,0d0), dcmplx(0d0,0d0), dcmplx(0d0,0d0) /)
      sigma_0to1(2,:) = (/ dcmplx(1d0,0d0), dcmplx(0d0,0d0), dcmplx(0d0,0d0) /)
      sigma_0to1(3,:) = (/ dcmplx(0d0,0d0), dcmplx(0d0,0d0), dcmplx(0d0,0d0) /)

      sigma_1to0(1,:) = (/ dcmplx(0d0,0d0), dcmplx(1d0,0d0), dcmplx(0d0,0d0) /)
      sigma_1to0(2,:) = (/ dcmplx(0d0,0d0), dcmplx(0d0,0d0), dcmplx(0d0,0d0) /)
      sigma_1to0(3,:) = (/ dcmplx(0d0,0d0), dcmplx(0d0,0d0), dcmplx(0d0,0d0) /)

      sigma_1to2(1,:) = (/ dcmplx(0d0,0d0), dcmplx(0d0,0d0), dcmplx(0d0,0d0) /)
      sigma_1to2(2,:) = (/ dcmplx(0d0,0d0), dcmplx(0d0,0d0), dcmplx(0d0,0d0) /)
      sigma_1to2(3,:) = (/ dcmplx(0d0,0d0), dcmplx(1d0,0d0), dcmplx(0d0,0d0) /)

      sigma_2to1(1,:) = (/ dcmplx(0d0,0d0), dcmplx(0d0,0d0), dcmplx(0d0,0d0) /)
      sigma_2to1(2,:) = (/ dcmplx(0d0,0d0), dcmplx(0d0,0d0), dcmplx(1d0,0d0) /)
      sigma_2to1(3,:) = (/ dcmplx(0d0,0d0), dcmplx(0d0,0d0), dcmplx(0d0,0d0) /)      

      jump_0to1 = matmul(sigma_1to0,sigma_0to1)
      jump_1to0 = matmul(sigma_0to1, sigma_1to0)
      jump_1to2 = matmul(sigma_2to1, sigma_1to2)
      jump_2to1 = matmul(sigma_1to2, sigma_2to1)

      g_ensemblesize = 5000

      !- Define the simulation constants
      g_omega0       = dcmplx(1d0,0d0);
      g_freq         = dcmplx(1d0,0d0);
      g_period1      = dcmplx(1d0,0d0);
      g_period2      = dcmplx(1d0,0d0);
      g_delta        = dcmplx(-0.1d0,0d0);

      !- Define the time-related variables
      g_dt  = .01d0
      g_t   = 0d0
      g_duration = 1000
      
      allocate(g_psi(3,g_ensemblesize))
      allocate(g_psiAvg(3,g_duration))

      !- Set up the inital wavefunction (the down state)
      do k=1,g_ensemblesize
        g_psi(:,k) = (/ dcmplx(1d0,0d0), dcmplx(0d0,0d0), dcmplx(0d0,0d0)/)
      end do
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

     open(11, file='prob1.dat')
     open(12, file='prob2.dat')
     open(13, file='prob3.dat')

   end subroutine initialize_files

   subroutine write_files

      g_psiAvg = 1/real(g_ensemblesize)*g_psiAvg

      do k=1,g_duration
        !Write out average wave function
         write(11,*) dble(g_psiAvg(1,k)*conjg(g_psiAvg(1,k)))
         write(12,*) dble(g_psiAvg(2,k)*conjg(g_psiAvg(2,k)))
         write(13,*) dble(g_psiAvg(3,k)*conjg(g_psiAvg(3,k)))

      end do

   end subroutine write_files

  subroutine close_files

     close(11)
     close(12)
     close(13)

   end subroutine close_files
   

end program threelevel
