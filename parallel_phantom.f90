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

   integer                     :: n, nx, ny, p, x_start, y_start
   integer                     :: i, k, l
   real(kind=8)                :: stddev, err, lambda
   real(kind=8), allocatable   :: f_phantom(:,:)[:], f_noisy(:,:)[:], f_filtered(:,:)[:]
   real(kind = 8) :: f_test(512,512)

! For timing:
   integer                     :: t0, t1, tb, te, clock_rate, clock_max

! For parallelizing the program
   integer :: np, me

   me = this_image()
   np = num_images()
   p = int(sqrt( real(np) ))

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
   if ( abs(sqrt(real(np)) - p) > 1e-10)  then
      if(me == 1) then
         write(*,*) abs( sqrt(real(np))-int(sqrt(real(np))) )>1e-10
         write(*,*) 'Error: np must be a perfect square'
      end if
      stop
   end if

   ! Check if the square root of np divides n
   if ( mod(n, p) /= 0) then
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
   nx = n/p
   ny = n/p

   ! Starting places
   x_start = nx*(mod(me-1, p))
   y_start = ny*((me-1)/p)
   
   ! Allocate the memory for all arrays
   allocate( f_phantom(n,n)[*], f_noisy(n,n)[*], f_filtered(n,n)[*] )


! Generate phantom image
   if (me == 1) then
      call shepp_logan( f_phantom, n)

      ! Timing for adding noise and filtering (part of the program that will be parallelised)
      call system_clock ( tb, clock_rate, clock_max )
   end if 

   sync all

! Distribute the image
   call distribute_image( f_phantom, x_start, y_start, nx, ny, me)
   sync all

! Add noise
   call add_noise( f_phantom, f_noisy, stddev, x_start, y_start, nx, ny )
   sync all


! Reassemble the image in the main processor
   if (me == 1) then
      do i = 2, np
         k = nx*(mod(i-1, p))
         l = ny*((i-1)/p)

         f_noisy(k+1 : k+nx, l + 1 : l+ny) = f_noisy(k+1 : k+nx, l + 1 : l+ny)[i]
      end do
   end if
   sync all

! Write the image to file
   if (me == 1) then
      call gif_images(f_phantom, f_noisy)
      call system_clock ( t1, clock_rate, clock_max )
      write(*,'(a,e8.2)') 'Time: ',real(t1-t0)/real(clock_rate)
   end if
   
! Timing for adding noise and filtering (part of the program that will be parallelised)
   call system_clock ( tb, clock_rate, clock_max )

contains

subroutine distribute_image( f_phantom, x_start, y_start, nx, ny, me)

   real(kind=8)              :: f_phantom(:,:)[*]
   integer, intent(in)       ::  x_start, y_start, nx, ny, me

   integer :: i


   f_phantom(x_start+1 : x_start+nx, y_start + 1 : y_start+ny)[me] = f_phantom(x_start+1 : x_start+nx, y_start + 1 : y_start+ny)[1]

end subroutine distribute_image

subroutine add_noise( f_phantom, f_noisy, stddev, x_start, y_start, nx, ny )

   real(kind=8), intent(in)  :: f_phantom(:,:)[*], stddev
   real(kind=8), intent(out) :: f_noisy(:,:)[*]
   integer, intent(in)       :: x_start, y_start, nx, ny

   real(kind=8), allocatable :: fn(:,:), stddev_rand

   allocate(fn(nx,ny))

   call RANDOM_SEED
   call RANDOM_NUMBER(fn)

   fn = (fn - 0.5)
   stddev_rand = sqrt(sum( fn**2)/(nx*ny))

   f_noisy(x_start+1 : x_start+nx, y_start + 1 : y_start+ny) = f_phantom(x_start+1 : x_start+nx, y_start + 1 : y_start+ny) + fn( 1:nx, 1:ny)

end subroutine add_noise


subroutine gif_images( model_image, noisy_image )
!
! Create gif-images and write to file
!
   implicit none

   real(kind=8)                  :: model_image(:,:), noisy_image(:,:)
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

end subroutine gif_images

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