! In this program, we apply diffusion filtering to the "Sheff Logan" test filtering
! This is part of a course for parallel computing. We have been provided a serial version of the code and our supposed to parallelize it
! The first step is to write code to 

!The first step is to distribute the image over the processors. Generate the noise-free image on the main processor. Define a block-wise domain decomposition, with dx blocks in x-direction and dy 
!blocks in y direction. Assign a domain to each processor and send the corresponding part of the image to this processor. Every processor should add noise to its part of the image. At the end of 
!the program, all the processors have to send their part of the (noisy and filtered) image to the main processor. This processor should assemble all these parts of the image to one image and write 
!this to file in gif format. Add this to the code. Do not apply filtering yet, test if it works correctly for the noisy image.


program phantom

   use gif_module

   implicit none

! Input parameters that define the problem:
   character(len=80)           :: arg, value

   integer                     :: n, nx, ny, np_root
   integer                     :: i, k, l
   real(kind=8)                :: stddev, err, lambda
   real(kind=8), allocatable   :: f_phantom(:,:), f_noisy(:,:), f_filtered(:,:)
   real(kind=8), allocatable   :: f_local(:,:)[:]
   real(kind = 8) :: f_test(512,512)

! For timing:
   integer                     :: t0, t1, tb, te, clock_rate, clock_max

! For parallelizing the program
   integer :: np, me

   me = this_image()
   np = num_images()
   np_root = int(sqrt( real(np) ))

! Defaults:
   n = 512
   stddev = 0.5
   lambda = 0.

! Read parameters from command line:
   do i = 1, command_argument_count()
      call get_command_argument(i, arg)

      select case(arg)
      case('-npx')
         call get_command_argument(i+1, value)
         read(value,'(i6)') n
      case('-sigma')
         call get_command_argument(i+1, value)
         read(value,'(e9.4)') stddev
      case('-lambda')
         call get_command_argument(i+1, value)
         read(value,'(e9.4)') lambda
      end select
   end do

   ! Check if the square root of np is an integer
   if ( abs(sqrt(real(np)) - np_root) > 1e-10)  then
      if(me == 1) then
         write(*,*) abs( sqrt(real(np))-int(sqrt(real(np))) )>1e-10
         write(*,*) 'Error: np must be a perfect square'
      end if
      stop
   end if

   ! Check if the square root of np divides n
   if ( mod(n, np_root) /= 0) then
      if(me == 1) then
         write(*,*) 'Error: n must be divisible by sqrt(np)'
      end if
      stop
   end if

   if (me == 1) then
      write(*,*) 
      call system_clock ( t0, clock_rate, clock_max )
   end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Input is set up correctly, we now start with the actual program.

   ! Block sizes
   nx = n/np_root
   ny = n/np_root
   
   ! Allocate the memory for all arrays
   if( me == 1) then
      allocate( f_phantom(n,n), f_noisy(n,n), f_filtered(n,n) )
   end if

   allocate( f_local(nx,ny)[*] )

! Generate phantom image
   if (me == 1) then
      call shepp_logan( f_phantom, n)

      ! Timing for adding noise and filtering (part of the program that will be parallelised)
      call system_clock ( tb, clock_rate, clock_max )
   end if 

   sync all

! Distribute the image
   call distribute_image( f_phantom, f_local, nx, ny, me, np, np_root)

! Add noise
   call add_noise( f_local, stddev, nx, ny )

   ! Reassemble the image in the main processor
   call gather_image(f_noisy, f_local, nx, ny, me, np, np_root)

   f_local = cg(f_local, lambda, nx, ny, me, np_root)

   call gather_image(f_filtered, f_local, nx, ny, me, np, np_root)

! Write the image to file
   if (me == 1) then
      call gif_images(f_phantom, f_noisy, f_filtered)
      call system_clock ( t1, clock_rate, clock_max )
      write(*,'(a,e8.2)') 'Time: ',real(t1-t0)/real(clock_rate)
   end if
   
! Timing for adding noise and filtering (part of the program that will be parallelised)
   call system_clock ( tb, clock_rate, clock_max )

contains

subroutine distribute_image( f_phantom, f_local, nx, ny, me, np, np_root)
   real(kind=8)              :: f_phantom(:,:), f_local(:,:)[*]
   integer, intent(in)       ::  nx, ny, me, np, np_root

   integer :: i, x_start, y_start

   if (me == 1) then
      do i = 1, np
         x_start = nx*(mod(i-1, np_root))
         y_start = ny*((i-1)/np_root)
         f_local(:,:)[i] = f_phantom(x_start+1 : x_start+nx, y_start + 1 : y_start+ny)
      end do
   end if
   sync all
end subroutine distribute_image

subroutine add_noise( f_local, stddev, nx, ny )
   real(kind=8), intent(in)  :: stddev
   integer, intent(in)       :: nx, ny
   real(kind=8)              :: f_local(:,:)[*]

   real(kind=8), allocatable :: fn(:,:), stddev_rand


   allocate(fn(nx,ny))

   call RANDOM_SEED
   call RANDOM_NUMBER(fn)

   fn = (fn - 0.5)
   stddev_rand = sqrt(sum( fn**2)/(nx*ny))

   f_local(:,:)[me] = f_local(:,:)[me] + fn
   sync all
end subroutine add_noise

subroutine gather_image(gather_array, local_array, nx, ny, me, np, np_root )
   real(kind=8)             :: gather_array(:,:), local_array(:,:)[*]
   integer, intent(in)      :: nx, ny, me, np, np_root
   integer                  :: i, k, l

   if (me == 1) then
      do i = 1, np
         k = nx*(mod(i-1, np_root))
         l = ny*((i-1)/np_root)

         gather_array(k+1 : k+nx, l + 1 : l+ny) = local_array(:, :)[i]
      end do
   end if
   sync all
end subroutine gather_image

subroutine gif_images( model_image, noisy_image, filtered_image )  
   !
   ! Create gif-images and write to file
   !
   implicit none

   real(kind=8), allocatable     :: model_image(:,:), noisy_image(:,:), filtered_image(:,:)
   integer, parameter            :: n_colours = 256
   integer, allocatable          :: gif_image(:,:,:)
   integer                       :: map(1:3,0:n_colours-1)
   character(len=23)             :: gif_name
   integer                       :: i,j,nx,ny

   nx = size(noisy_image,1)
   ny = size(noisy_image,2)
   !
   ! Colour map grey
   do i = 0,n_colours-1
      map(:,i ) = i
   end do

   ! Model image
   allocate(gif_image(1,nx,ny))
   gif_name = 'phantom.gif'
   do j = 1, ny
      do i = 1, nx
         gif_image(1,i,j) = int(model_image(i,n-j+1)*n_colours)
      end do
   end do
   where ( gif_image > n_colours-1 ) gif_image = n_colours-1
   where ( gif_image < 0 ) gif_image = 0
   call write_animated_gif( gif_name, gif_image, map )

   ! Noisy image
   gif_name = 'noisy.gif'
   do j = 1, ny
      do i = 1, nx
         gif_image(1,i,j) = int(noisy_image(i,n-j+1)*n_colours)
      end do
   end do
   where ( gif_image > n_colours-1 ) gif_image = n_colours-1
   where ( gif_image < 0 ) gif_image = 0
   call write_animated_gif( gif_name, gif_image, map )

   ! Filtered image
   gif_name = 'filtered.gif'
   do j = 1, ny
      do i = 1, nx
         gif_image(1,i,j) = int(filtered_image(i,n-j+1)*n_colours)
      end do
   end do
   where ( gif_image > n_colours-1 ) gif_image = n_colours-1
   where ( gif_image < 0 ) gif_image = 0
   call write_animated_gif( gif_name, gif_image, map )

end subroutine gif_images

function cg( f_local, lambda, nx, ny, me, np_root ) result(f_local_filtered)
   !
   ! CG algorithm 
   !
   implicit none

   real(kind=8), intent(in)    :: lambda
   real(kind=8), intent(in)    :: f_local(:,:)[*]
   integer, intent(in)         :: nx, ny, me, np_root
   real(kind=8)                :: f_local_filtered(size(f_local,1),size(f_local,2))

   real(kind=8)                :: b(size(f_local,1),size(f_local,2))
   real(kind=8)                :: x(size(f_local,1),size(f_local,2))
   real(kind=8)                :: r(size(f_local,1),size(f_local,2))
   real(kind=8), allocatable   :: p(:,:)[:]
   real(kind=8)                :: Ap(size(f_local,1),size(f_local,2))
   real(kind=8)                :: alpha, beta, gamma, xi
   real(kind=8)                :: tol, normr, normb
   integer                     :: it, maxit


   allocate(p(nx,ny)[*])

   if ( lambda == 0. ) then
      f_local_filtered = f_local(:,:)[me]
      return
   end if

   maxit = 100
   tol   = 1.e-8

   b = lambda*f_local(:,:)[me]

   r = b - mv( f_local, lambda, nx, ny, me, np_root )

   p(:,:)[me] = r

   sync all

   Ap = mv( p, lambda, nx, ny, me, np_root )

   gamma = dot_butterfly_2d(r,r)
   normr = sqrt( gamma )
   normb = sqrt( dot_butterfly_2d(b,b) )

   it = 0
   do while ( normr/normb > tol .and. it < maxit )

      xi = dot_butterfly_2d(p(:,:)[me],Ap)
      alpha = gamma/xi

      x = x + alpha*p(:,:)[me]
      r = r - alpha*Ap

      beta  = 1/gamma
      gamma = dot_butterfly_2d(r,r)
      normr = sqrt(gamma)
      beta  = gamma*beta
      p(:,:)[me] = r + beta*p(:,:)[me]
      Ap    = mv( p, lambda, nx, ny, me, np_root )

      it = it+1
   end do

   f_local_filtered = x
   write(*,'(a,i4,1x,a,e9.3)') 'CG terminated after ',it, &
           'iterations with relative residual norm ', normr/normb

   deallocate(p)
end function cg

function mv( v, lambda, nx, ny, me, np_root) result(w)
   ! Matrix vector multiplication: w = lambda v + A v
   ! in which A represents -Laplacian with Neumann boundary conditions
   !
   implicit none

   real(kind=8), intent(in)    :: lambda
   real(kind=8), intent(in)    :: v(:,:)[*]
   integer, intent(in)         :: np_root

   real(kind=8)                :: w(size(v,1),size(v,2))
   real(kind=8)                :: v_halo(0:size(v,1)+1,0:size(v,2)+1)
   integer                     :: nx, ny
   integer                     :: i, j
   integer                     :: me, p

   ! Add extra row/column around the domain to account for boundary conditions
   ! or data on other processor. Note that these points are outside the domain.
   ! They can be eliminated by the homogeneous Neumann condition. This condition
   ! says that the points on a row/column outside the boundary are equal to the 
   ! neighbouring points on the row/column inside the boundary.

   v_halo = 0.
   v_halo(1 : nx, 1: ny) = v  ! Points in the interior

   ! If you are a left processor, set a left halo. Else, retreive the boundary from your neighbour
   if (mod(me-1,p) == 0) then
      v_halo(0,1:ny) = v_halo(1,1:ny)
   else if (mod(me,p) /= 0) then
      v_halo(0, 1:ny) = v(nx, 1:ny)[me - 1]
   end if

   ! If you are the right processor, set a right halo
   if (mod(me,p) == 0) then
      v_halo(nx+1 , 1:ny) = v_halo(nx,1:ny)
   else if (mod(me,p) /= 0) then
      v_halo(nx+1, 1:ny) = v(1, 1:ny)[me + 1]
   end if

   ! If you are the top processor, set a top halo
   if (me > np - p) then
      v_halo(1:nx, ny+1) = v_halo(1:nx, ny)
   else if (me <= np - p) then
      v_halo(1: nx, ny+1 ) = v( 1: nx , 1)[me + p]
   end if

   ! If you are the bottom processor, set a bottom halo
   if (me <= p) then
      v_halo(1:nx,0) = v_halo(1:nx,1)
   else if (me > p) then
      v_halo(1 : nx , 0) = v(1 :  nx , ny)[me - p]
   end if

   sync all

   ! Now we can just perform the stencil operation:
   w = lambda*v(:,:)[me]
   do j = 1, ny
      do i = 1, nx
         w(i,j) = w(i,j) &
             + 4.*v_halo(i, j) -v_halo( i-1,  j) -v_halo( i+1, j) -v_halo( i, j-1) -v_halo( i, j+1)
      end do
   end do
end function mv

real(8) function dot_butterfly_2d(x,y)
   ! Flatten the 2D arrays and call the 1D butterfly reduction
   implicit none
   real(kind=8), dimension(:,:), intent(in) :: x, y
   real(kind=8), dimension(:), allocatable :: x_flat, y_flat
   real(kind=8) :: s

   ! Allocate the 1D arrays to hold flattened data
   allocate(x_flat(size(x)))
   allocate(y_flat(size(y)))

   ! Flatten the 2D arrays
   x_flat = reshape(x, [size(x)])
   y_flat = reshape(y, [size(y)])

   s = dot_butterfly(x_flat, y_flat)

   dot_butterfly_2d = s

   deallocate(x_flat, y_flat)
end function dot_butterfly_2d

real(8) function dot_butterfly(x, y)
   implicit none
   real(kind=8), dimension(:), intent(in) :: x, y
   ! note: the standard requires coarrays inside functions to be
   ! persistent in one way or the other so that one doesn't accidently
   ! recreate the coarray many times (which causes an implicit sync)
   real(kind=8), save :: s[*], t[*]
   integer :: me, np
   integer :: nlev, l, other
   
   np = num_images()
   me = this_image()
   

   if (np /= 2**(log(dble(np))/log(2.d0))) then
       stop 'The butterfly_dot variant only works for num_images() a power of 2'
   end if
   
   s = dot_product(x, y)
   
   !
   ! Butterfly reduction: E.g., for np=4 we get
   !
   !  s1  s2  s3  s4
   !    \/      \/
   !    /\      /\
   ! s12 s12 s34  s34
   !    \___/
   !    /   \ (and similar for p2 <-> P4
   nlev = log(dble(np))/log(2.d0)
   
   do l=0, nlev-1
       other = ieor(me-1, 2**l)+1
       t[other] = s
       sync all
       !write(*,'(A,I1,A,I1,A,I1,A,G8.2,A,G8.2)') 'level ',l,' ',me,'->',other, ', s=',s, ', t=',t
       s=s+t
       sync all
   end do
   
   dot_butterfly = s
end function dot_butterfly

subroutine shepp_logan( f, n )

!
!    Larry Shepp, Ben Logan,
!    The Fourier reconstruction of a head section,
!    IEEE Transactions on Nuclear Science,
!    Volume  NS-21, June 1974, pages 21-43.
!
!  Parameters:
!
!    Input, integer n, the number of points in each direction.
!
!    Output, real f(n,n), the image values.
!
!  Local parameters:
!
!    Local, integer CHOICE:
!    1, use Archibald's (and Shepp and Logan's) level values;
!    2, use Matlab's level values;
!    3, use Matlab's enhanced contrast level values.
!

   implicit none

   integer ( kind = 4 ) i,j,n
   real ( kind = 8 ) c(4)
   real ( kind = 8 ), dimension ( 4 ) :: c1 = (/ 2.0E+00, -0.98E+00, -0.02E+00, +0.01E+00 /) 
   real ( kind = 8 ), dimension ( 4 ) :: c2 = (/ 1.0E+00, -0.98E+00, -0.02E+00, +0.01E+00 /) 
   real ( kind = 8 ), dimension ( 4 ) :: c3 = (/ 1.0E+00, -0.8E+00,  -0.2E+00,  +0.1E+00  /) 
   integer ( kind = 4 ) choice
   real ( kind = 8 ) eta1
   real ( kind = 8 ) eta2
   real ( kind = 8 ) f(n,n)
   real ( kind = 8 ), parameter :: pi = 3.141593E+00
   real ( kind = 8 ) xi1
   real ( kind = 8 ) xi2
   real ( kind = 8 ) x, x_max, x_min, dx
   real ( kind = 8 ) y, y_max, y_min, dy

   x_min = -1.0E+00
   x_max = +1.0E+00
   dx = (x_max-x_min)/n
   y_min = -1.0E+00
   y_max = +1.0E+00
   dy = (y_max-y_min)/n
 
   choice = 3

   if ( choice == 1 ) then
     c = c1
   else if ( choice == 2 ) then
     c = c2
   else
     c = c3
   end if


   y = y_min + dy/2
   do j = 1, n
      x = x_min + dx/2
      do i = 1, n

         if ( ( x / 0.69E+00 )**2 + ( y / 0.92E+00 )**2 <= 1.0E+00 ) then
            f(i,j) = f(i,j) + c(1)
         end if
 
         if ( ( x / 0.6624E+00 )**2 + ( ( y + 0.0184E+00 ) / 0.874E+00 )**2 <= 1.0E+00 ) then
            f(i,j) = f(i,j) + c(2)
         end if

         xi1  =   ( x - 0.22E+00 ) * cos ( 0.4E+00 * pi ) &
                  + y              * sin ( 0.4E+00 * pi )
         eta1 = - ( x - 0.22E+00 ) * sin ( 0.4E+00 * pi ) &
                  + y              * cos ( 0.4E+00 * pi )

         xi2  =   ( x + 0.22E+00 ) * cos ( 0.6E+00 * pi ) &
                  + y              * sin ( 0.6E+00 * pi )
         eta2 = - ( x + 0.22E+00 ) * sin ( 0.6E+00 * pi ) &
                  + y              * cos ( 0.6E+00 * pi )
 
         if ( ( xi1 / 0.31E+00 )**2 + ( eta1 / 0.11E+00 )**2 <= 1.0E+00 .or. &
              ( xi2 / 0.41E+00 )**2 + ( eta2 / 0.16E+00 )**2 <= 1.0E+00 ) then
            f(i,j) = f(i,j) + c(3)
         end if

         if ( (( x              /0.21E+00  )**2 + (( y - 0.35E+00 )/0.25E+00  )**2 <= 1.0E+00 ) .or. &
              (( x              /0.046E+00 )**2 + (( y - 0.10E+00 )/0.046E+00 )**2 <= 1.0E+00 ) .or. &
              (( x              /0.046E+00 )**2 + (( y + 0.1E0+00 )/0.046E+00 )**2 <= 1.0E+00 ) .or. &
              (((x + 0.080E+00) /0.046E+00 )**2 + (( y + 0.605E+00)/0.023E+00 )**2 <= 1.0E+00 ) .or. &
              (( x              /0.023E+00 )**2 + (( y + 0.605E+00)/0.023E+00 )**2 <= 1.0E+00 ) .or. &
              (((x - 0.06E+00 ) /0.023E+00 )**2 + (( y + 0.605E+00)/0.023E+00 )**2 <= 1.0E+00 ) ) then
            f(i,j) = f(i,j) + c(4)
         end if
         x = x + dx
      end do
      y = y + dy
   end do

   return
end subroutine shepp_logan

end program phantom