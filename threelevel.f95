program threelevel

   implicit none
   
   !- Constants
   real*8, parameter :: pi = 4.0*atan(1d0)
   
   !-- Simulation Variables
   complex*16 :: g_omega1
   complex*16 :: g_omega2
   complex*16 :: g_freq
   complex*16 :: g_omegax
   complex*16 :: g_omegay
   complex*16 :: g_delta, g_sdelta
   complex*16 :: g_lamb1, g_lamb2
   
   !-- Quantum Variables
   complex*16 :: g_sigmax(3,3), g_sigmaz(3,3), g_sigmay(3,3) !Pauli matrices
   complex*16 :: g_unity(3,3) ! Complex unity matrix
   complex*16 :: g_E1(3,3) ! The corrected energy matrix for qubit 1
   complex*16 :: g_E2(3,3) ! The corrected energy matrix for qubit 2
   
   complex*16 :: g_psi(3,2) ! Both wavefunctions
   
   complex*16 :: g_Ham(3,3) ! The Hamiltonian
   
   !-- Program Variables   
   integer*4  :: i !counters
   real*8 :: g_dt,g_t !g_dt = timestep,g_t = current time
   integer*4 :: g_duration !Duration of time evolution
   
   real*8 :: g_tg, g_Am, g_Api ! (See eq. (23) from Schutjens)
    
   complex*16 :: g_HamNumer(3,3) ! Numerator of the Hamiltonian
   complex*16 :: g_HamDenom(3,3) ! Denominator of the Hamiltonian
   
   complex*16 :: ii=dcmplx(0._8,1._8) ! Complex 'i'
   
   !-- Dummy testing stuff
   complex*16 :: g_zero_one(9), g_zero_two(9) 
   complex*16 :: g_one_one(9), g_one_two(9)
   real*8     :: expect_rho(4)
   
   open(11,file='psi1.dat')
   open(12,file='psi2.dat')
   open(13,file='omega.dat')
   open(14,file='rho.dat')
   
   call setUpVariables()
      
   do i = 1,g_duration
      
      if(modulo(i,100)==0) then
         call expectDensityMatrix()
      end if
      call calcPsi()
!       call printOut()
      call writeOut()
      
!       print *, i
      
   end do
   
   close(11)
   close(12)
   close(13)
   
contains

   subroutine calcPsi

      call calcPsiSub(1, g_t)
      call calcPsiSub(2, g_t)
      g_t = g_t + g_dt
      
      write(13, *) realpart(g_omegax), realpart(g_omegay)
      
      ! Normalise the result
      g_psi(:,1) = g_psi(:,1) / sqrt(dot_product(g_psi(:,1),g_psi(:,1)))
      g_psi(:,2) = g_psi(:,2) / sqrt(dot_product(g_psi(:,2),g_psi(:,2)))
      
      ! Calculate the Hamiltonian for the Energy
      call calcHamiltonian(g_Ham,1,g_t) 
      
   end subroutine calcPsi
   
   subroutine calcPsiSub(qubit,t)
      integer*4, intent(in) :: qubit
      real*8, intent(in) :: t
      real*8 :: localT
      
      localT = t
      
      call calcHamiltonian(g_Ham,qubit,localT)
      
      g_HamNumer = 1*g_unity - g_Ham * g_dt * II * 0.5_8
      
      localT = t + g_dt
      
      call calcHamiltonian(g_Ham,qubit, localT)
      
      g_HamDenom = 1*g_unity + g_Ham * g_dt * II * 0.5_8
      
      call invertComplex(g_HamDenom)
      
      g_Ham = matmul(g_HamDenom,g_HamNumer)  
      g_psi(:,qubit) = matmul(g_Ham,g_psi(:,qubit))      
      
   end subroutine calcPsiSub
   
   subroutine calcHamiltonian(H,qubit,t)
      complex*16, intent(inout) :: H(:,:)
      integer*4, intent(in) :: qubit
      real*8, intent(in) :: t
      
      real*8 :: sig 
      real*8 :: wx
      real*8 :: beta
      
      real*8 :: A,B,C,D ! Holding variables for OmegaY
      
      sig = g_tg / 6d0
      wx = g_sdelta * pi!/ 2d0
      beta = - 1d0 / (2d0 * g_delta)
      
      A = (g_Api * exp( -(t - g_tg/2d0) ** 2d0 / (2d0 * sig**2d0) ) )
      B = t - g_tg/2d0
      C = 1 - g_Am * cos( (t - g_tg/2d0) * wx)
      D = A * g_Am * wx * sin( B * wx )

      g_omegax = g_Api * exp(-(t-g_tg/2)**2d0 / (2d0 * sig**2d0) ) &
         * (1 - g_Am * cos( wx * (t-g_tg/2d0) ))
      
      ! There is a reason this variable is called as is...   
      g_omegay  = - beta * ( A * B * C / sig**2d0 + D )
      
!       g_omegax       = dcmplx(100d0,0d0);
!       g_omegay       = dcmplx(100d0,0d0);
                  
      if(qubit==1) then      
         H  = g_E1
      else
         H  = g_E2
      end if
      
      H  = H  & 
         + g_omegax/2d0 * sin(g_omega1 * t) * g_sigmay &
         + g_omegay/2d0 * cos(g_omega1 * t) * g_sigmax
               
   end subroutine calcHamiltonian
   
   subroutine printOut
      complex*16 :: ground(3) = (/dcmplx(1d0,0d0), dcmplx(0d0,0d0), dcmplx(0d0,0d0)/)
      complex*16 :: excited(3) = (/dcmplx(0d0,0d0), dcmplx(1d0,0d0), dcmplx(0d0,0d0)/)
      complex*16 :: leak(3) = (/dcmplx(0d0,0d0), dcmplx(0d0,0d0), dcmplx(1d0,0d0)/)
      
      print *, abs(dot_product(g_psi(:,1),g_psi(:,1)))**2d0, &
               abs(dot_product(g_psi(:,2),g_psi(:,2)))**2d0!, &
!                dot_product(g_psi(:,1),matmul(g_Ham,g_psi(:,1)))
      
   end subroutine printOut
   
   subroutine writeOut
      complex*16 :: ground(3) = (/dcmplx(1d0,0d0), dcmplx(0d0,0d0), dcmplx(0d0,0d0)/)
      complex*16 :: excited(3) = (/dcmplx(0d0,0d0), dcmplx(1d0,0d0), dcmplx(0d0,0d0)/)
      complex*16 :: leak(3) = (/dcmplx(0d0,0d0), dcmplx(0d0,0d0), dcmplx(1d0,0d0)/)
      
      write(11, *) abs(dot_product(ground,g_psi(:,1)))**2d0, &
               abs(dot_product(excited,   g_psi(:,1)))**2d0, &
               abs(dot_product(leak,      g_psi(:,1)))**2d0
               
      write(12, *) abs(dot_product(ground,g_psi(:,2)))**2d0, &
               abs(dot_product(excited,   g_psi(:,2)))**2d0, &
               abs(dot_product(leak,      g_psi(:,2)))**2d0
      
   end subroutine writeOut
   
   subroutine setUpVariables
      complex*16 :: zero(3), one(3), two(3) ! The Basis States
      
      !- Define the simulation constants
      ! All times are defines in ns; Therefore all frequencies are in GHz
      g_omega1       = dcmplx(5.508d0,0d0)*2d0*pi    ! ⍵1     = 5.508 GHz
      g_omega2       = dcmplx(5.902d0,0d0)*2d0*pi    ! ⍵2     = 5.902 GHz
      g_delta        = dcmplx(-350E-3,0d0)*2d0*pi    ! ∆/2π   = -350 MHz
      g_sdelta       = dcmplx(45E-3,0d0)*2d0*pi      ! ∂/2π   = 45 MHz
      g_lamb1        = dcmplx(1d0,0d0)               ! λ1     = 1
      g_lamb2        = dcmplx(sqrt(2d0),0d0)         ! λ2     = √2
      g_tg           = 24d0                          ! gt     = 17 ns
      
      g_Am           = 0.9d0
      g_Api          = 1d0
      
      g_freq         = dcmplx(1d0,0d0);
      
      !- Define the Corrected energies
      ! Harmonic plus/ minus the corrections (∂ and ∆) on the energy
      g_E1(1,:)  = (/ dcmplx(0d0,0d0), dcmplx(0d0,0d0), dcmplx(0d0,0d0)/)
      g_E1(2,:)  = (/ dcmplx(0d0,0d0), g_omega1, dcmplx(0d0,0d0)/)
      g_E1(3,:)  = (/ dcmplx(0d0,0d0), dcmplx(0d0,0d0), 2d0*g_omega1 + g_delta/)
        
      g_E2(1,:)  = (/ dcmplx(0d0,0d0), dcmplx(0d0,0d0), dcmplx(0d0,0d0)/)
      g_E2(2,:)  = (/ dcmplx(0d0,0d0), g_omega2, dcmplx(0d0,0d0)/)
      g_E2(3,:)  = (/ dcmplx(0d0,0d0), dcmplx(0d0,0d0), 2d0*g_omega2 + g_delta/)
      
      !- Define the Pauli Spin Matrices (with resp. lambda's)
      g_sigmax(1,:)  = (/ dcmplx(0d0,0d0), g_lamb1, dcmplx(0d0,0d0)/)
      g_sigmax(2,:)  = (/ g_lamb1, dcmplx(0d0,0d0), g_lamb2/)
      g_sigmax(3,:)  = (/ dcmplx(0d0,0d0), g_lamb2, dcmplx(0d0,0d0)/)
      
      g_sigmay(1,:)  = (/ dcmplx(0d0,0d0), - II * g_lamb1, dcmplx(0d0,0d0)/)
      g_sigmay(2,:)  = (/ II * g_lamb1, dcmplx(0d0,0d0), - II * g_lamb2/)
      g_sigmay(3,:)  = (/ dcmplx(0d0,0d0), II * g_lamb2, dcmplx(0d0,0d0)/)
      
      g_sigmaz(1,:)  = (/ dcmplx(0d0,0d0), dcmplx(0d0,0d0), dcmplx(0d0,0d0)/)
      g_sigmaz(2,:)  = (/ dcmplx(0d0,0d0), dcmplx(1d0,0d0), dcmplx(0d0,0d0)/)
      g_sigmaz(3,:)  = (/ dcmplx(0d0,0d0), dcmplx(0d0,0d0), dcmplx(-2d0,0d0)/)
     
      g_unity(1,:)   = (/ dcmplx(1d0,0d0), dcmplx(0d0,0d0), dcmplx(0d0,0d0)/)
      g_unity(2,:)   = (/ dcmplx(0d0,0d0), dcmplx(1d0,0d0), dcmplx(0d0,0d0)/)
      g_unity(3,:)   = (/ dcmplx(0d0,0d0), dcmplx(0d0,0d0), dcmplx(1d0,0d0)/)  
      
      !- Define the time-related variables
      g_dt  = .01d0 ! Timestep in nanoseconds
      g_t   = 0d0
      g_duration = 5000
      
      !- Set up the inital wavefunction (the up state)
      g_psi(:,1) = (/ dcmplx(0d0,0d0), dcmplx(1d0,0d0), dcmplx(0d0,0d0)/)
      g_psi(:,2) = (/ dcmplx(1d0,0d0), dcmplx(0d0,0d0), dcmplx(0d0,0d0)/)
      !- Set up the inital wavefunction (the down state)
      !g_psi = (/ dcmplx(1d0,0d0), dcmplx(0d0,0d0) /)
      
      !- Set all the possible basis states
      zero  = (/ dcmplx(1d0,0d0), dcmplx(0d0,0d0), dcmplx(0d0,0d0)/)
      one   = (/ dcmplx(0d0,0d0), dcmplx(1d0,0d0), dcmplx(0d0,0d0)/)
      two   = (/ dcmplx(0d0,0d0), dcmplx(0d0,0d0), dcmplx(1d0,0d0)/)
      
      call tensorProduct(zero,one,g_zero_one)
      call tensorProduct(zero,two,g_zero_two)
      call tensorProduct(one,one,g_one_one)
      call tensorProduct(one,two,g_one_two)
      
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
   
   subroutine tensorProduct(a,b,c)
      complex*16, intent(in) :: a(:),b(:)
      complex*16, intent(out) :: c(:)
      integer*4 :: k
      integer*4 :: sz
      
      sz = size(a,dim=1)
      
      do k=0,sz-1
         c(k*sz+1:(k+1)* sz) = a(k+1) * b
      end do
      
   end subroutine tensorProduct
   
   subroutine calcDensityMatrix(a,b,rho)
      complex*16, intent(in) :: a(:),b(:)
      complex*16, intent(out) :: rho(:,:)
      integer*4 :: k
      integer*4 :: sz
      
      sz = size(a,dim=1)
      
      do k=1,sz
         rho(k,:) = a(k) * conjg(b)
      end do
      
   end subroutine calcDensityMatrix
   
   subroutine expectDensityMatrix
      complex*16 :: c(9), rho(9,9)
      complex*16 :: tr
      integer*4 :: k,l
      
      call tensorProduct(g_psi(:,2),g_psi(:,1),c)
      call calcDensityMatrix(c,c,rho)
      
      do k=1,9
         do l=1,9
            write(14, *) realpart(rho(k,l))
         end do
      end do
      
   end subroutine expectDensityMatrix
   
   subroutine trace(A,tr)
      complex*16, intent(in) :: A(:,:)
      complex*16, intent(out) :: tr
      
      tr = dcmplx(0d0,0d0)
      
      do i=0,size(A,dim=1)
         tr = tr + A(i,i)
      end do
      
   end subroutine trace

end program threelevel
