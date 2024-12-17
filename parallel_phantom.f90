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

   integer                     :: n, nx, ny, np_row, np_col
   integer                     :: i, k, l
   real(kind=8)                :: stddev, err, lambda
   real(kind=8), allocatable   :: f_phantom(:,:), f_noisy(:,:), f_filtered(:,:)
   real(kind=8), allocatable   :: f_local(:,:)[:]

! For timing:
   integer                     :: t0, t_gen, t_distr, t_noise, t_gather, t_comp, t1, clock_rate, clock_max

! For parallelizing the program
   integer :: np, me

   ! Debugging
   logical :: debug = .false.


   OPEN(UNIT=10, FILE='output.out', STATUS='old', POSITION='APPEND')

   me = this_image()
   np = num_images()

   ! Defaults:
   n = 512
   stddev = 0.5
   lambda = 0.1

   ! amount of processors in a row and in a column respectively
   np_row = 4
   np_col = 1

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
         ! case('-debug')
         !    call get_command_argument(i+1, value)
         !    read(value, '(l1)') debug
         ! case ('-np_row')
         !    call get_command_argument(i+1, value)
         !    read(value, '(i6') np_row
         ! case ('-np_col')
         !    call get_command_argument(i+1, value)
         !    read(value, '(i6') np_col
      end select
   end do

   ! Check if np_row times np_col is equal to np
   if ( abs(np - np_row*np_col) > 1e-10)  then
      if(me == 1) then
         write(*,*) 'Error: np_row*np_col must be equal to np'
      end if
      stop
   end if

   ! Check if number of pixels is neatly divisible by the amount of processors
   if ( mod(n, np) /= 0) then
      if(me == 1) then
         write(*,*) 'Error: n must be divisible by np'
         ! Add a small amount of pixels to make n divisible by np
         write (*,*) "Added ", np - mod(n, np), " pixels to n"
      end if
      n = n + np - mod(n, np)
   end if

   ! Start clock on the program
   if (me == 1) then
      write(*,*) 
      call system_clock ( t0, clock_rate, clock_max )
   end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Input is set up correctly, we now start with the actual program.

   ! Block sizes
   nx = n/np_row
   ny = n/np_col
   
   ! Allocate the memory for all arrays
   if( me == 1) then
      allocate( f_phantom(n,n), f_noisy(n,n), f_filtered(n,n) )
   end if

   allocate( f_local(nx,ny)[*] )

! Generate phantom image
   if (me == 1) then
      call shepp_logan( f_phantom, n)

      call system_clock(t_gen, clock_rate, clock_max )
      if (debug) then
         write(*,'(a, e8.2)') 'Time for generating the phantom image: ',real(t_gen-t0)/real(clock_rate)
      end if
   end if 

   sync all

   ! Distribute the image
   call distribute_image( f_phantom, f_local, me, np_row, np_col)

   if (me == 1) then
      call system_clock(t_distr, clock_rate, clock_max )
      if (debug) then
         write(*,'(a, e8.2)') 'Time for distributing the image: ',real(t_distr-t_gen)/real(clock_rate)
      end if
   end if

   ! Add noise
   call add_noise( f_local, stddev, nx, ny )

   ! Debug timing for noise adding
   if (me == 1) then
      call system_clock(t_noise, clock_rate, clock_max )
      if (debug) then
         write(*,'(a, e8.2)') 'Time for adding noise: ',real(t_noise-t_distr)/real(clock_rate)
      end if
   end if

   ! Reassemble the image in the main processor
   call gather_image(f_noisy, f_local, me, np_row, np_col)

   ! Debugging timing for gathering image
   if (me == 1) then
      call system_clock(t_gather, clock_rate, clock_max )
      if (debug) then
         write(*,'(a, e8.2)') 'Time for gathering the noisy image: ',real(t_gather-t_noise)/real(clock_rate)
         write(*,*)
      end if
   end if

   f_local = cg(f_local, lambda, me, np_row, np_col)

   ! Debugging timing for conjugate gradient
   if (me == 1) then
      call system_clock(t_comp, clock_rate, clock_max )
      if (debug) then
         write(*,'(a, e8.2)') 'Time for computation: ',real(t_comp-t_gather)/real(clock_rate)
      end if
   end if

   call gather_image(f_filtered, f_local, me, np_row, np_col)

   ! Write the image to file
   ! If debug is on, report timing on gathering the image
   ! Always report the total time of the program
   if (me == 1) then
      call gif_images(f_phantom, f_noisy, f_filtered)
      call system_clock ( t1, clock_rate, clock_max )
      if (debug) then
         write(*,'(a, e8.2)') 'Time for gathering and writing to file ',real(t1-t_comp)/real(clock_rate)
         write(*,*)
      end if
      WRITE(10, '(A1, E15.7, A1, I0, A1, I0)') ',',  real(t1-t0)/real(clock_rate), ',', np, ',', n
   end if

   ! Close the output file
   CLOSE(10)

contains

subroutine distribute_image( f_phantom, f_local, me, np_row, np_col)
   real(kind=8)                :: f_phantom(:,:), f_local(:,:)[*]
   real(kind = 8), allocatable :: f_flattened(:)[:]
   integer, intent(in)         ::  me, np_row, np_col

   integer :: i, x_start, y_start, nx, ny

   nx = size(f_local,1)
   ny = size(f_local,2)

   allocate(f_flattened(nx*ny)[*])
   if (me == 1) then
      do i = 1, np_row*np_col
         ! Calculate the starting position of the left lowest pixel of the starting array according to x_start = nx * (i-1) mod np_row
         x_start = nx*(mod(i-1, np_row))
         ! Calculate the starting y position according to y_start = ny * floor((i-1)/np_col)
         y_start = ny*((i-1)/np_row)

         ! Fortran starts indexing at 1, so add 1 to the starting index.
         f_flattened(:)[i] = reshape(f_phantom(x_start+1 : x_start+nx, y_start + 1 : y_start+ny), [nx*ny])
      end do
   end if
   sync all

   f_local(:,:)[me] = reshape(f_flattened(:)[me], [nx,ny])
   sync all
end subroutine distribute_image

subroutine add_noise( f_local, stddev, nx, ny )
   real(kind=8), intent(in)  :: stddev
   integer, intent(in)       :: nx, ny
   real(kind=8)              :: f_local(:,:)[*]

   real(kind=8), allocatable :: fn(:,:), stddev_rand


   allocate(fn(nx,ny))

   call RANDOM_SEED()
   call RANDOM_NUMBER(fn)

   fn = (fn - 0.5)
   stddev_rand = sqrt(sum( fn**2)/(nx*ny))

   f_local(:,:)[me] = f_local(:,:)[me] + (stddev/stddev_rand)*fn
   sync all
end subroutine add_noise

subroutine gather_image(gather_array, local_array, me, np_row, np_col )
   real(kind=8)              :: gather_array(:,:), local_array(:,:)[*] 
   real(kind=8), allocatable :: local_flattened(:)[:]
   integer, intent(in)       :: me, np_row, np_col
   integer                   :: i, x_start, y_start, nx, ny

   ! Setting up some dummy variables.
   nx = size(local_array,1)
   ny = size(local_array,2)


   ! Setup a flattened version of the coarray due to 2D coarrays being very slow.
   allocate(local_flattened(nx*ny)[*])
   local_flattened(:)[me] = reshape(local_array(:,:)[me], [nx*ny])


   if (me == 1) then
      do i = 1, np_row*np_col
         ! Calculate the starting position of the left lowest pixel of the starting array according to x_start = nx * (i-1) mod np_row
         x_start = nx*(mod(i-1, np_row))
         ! Calculate the starting y position according to y_start = ny * floor((i-1)/np_col)
         y_start = ny*((i-1)/np_row)

         ! Fortran starts indexing at 1, so add 1 to the starting index.
         gather_array(x_start+1 : x_start+nx, y_start + 1 : y_start+ny) = reshape(local_flattened(:)[i], [nx,ny])
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

function cg( f_local, lambda, me, np_row, np_col ) result(f_local_filtered)
   !
   ! Conjugate gradient algorithm with parallel MV products and inproducts
   !
   implicit none

   real(kind=8), intent(in)    :: lambda
   integer, intent(in)         :: me, np_row, np_col
   real(kind=8)                :: f_local(:,:)[*]
   real(kind=8)                :: f_local_filtered(size(f_local,1),size(f_local,2))

   real(kind=8)                :: b(size(f_local,1),size(f_local,2))
   real(kind=8)                :: r(size(f_local,1),size(f_local,2))
   real(kind=8), allocatable   :: p(:,:)[:]
   real(kind=8)                :: Ap(size(f_local,1),size(f_local,2))
   real(kind=8)                :: alpha, beta, gamma, xi
   real(kind=8)                :: tol, normr, normb
   integer                     :: it, maxit, nx, ny


   nx = size(f_local,1)
   ny = size(f_local,2)

   allocate(p(nx,ny)[*])


   ! If lambda = 0 then no filtering is applied
   if ( lambda == 0. ) then
      f_local_filtered = f_local(:,:)
      return
   end if

   ! Maximum amount of allowed tolerance for relative residual.
   ! Stop once the relative residual is smaller than tol or the amount of iterations is bigger than maxit.

   maxit = 100
   tol   = 1.e-8

   ! First steps of the algorithm
   b = lambda*f_local(:,:)
   r = b - mv( f_local, lambda, me, np_row, np_col )

   p(:,:) = r
   sync all
   Ap = mv( p, lambda, me, np_row, np_col )
   gamma = dot_butterfly_2d(r,r)
   normr = sqrt( gamma )
   normb = sqrt( dot_butterfly_2d(b,b) )

   ! Start main iteration loop of CG
   it = 0
   do while ( normr/normb > tol .and. it < maxit )

      xi = dot_butterfly_2d(p(:,:),Ap)
      alpha = gamma/xi

      f_local(:,:) = f_local(:,:) + alpha*p(:,:)
      r = r - alpha*Ap

      beta  = 1/gamma
      gamma = dot_butterfly_2d(r,r)
      normr = sqrt(gamma)
      beta  = gamma*beta
      p(:,:) = r + beta*p(:,:)
      Ap    = mv( p, lambda, me, np_row, np_col)

      it = it+1
   end do

   f_local_filtered = f_local(:,:)
   if (me == 1) then
      WRITE(10, '(I0, A1, E15.7)', ADVANCE='NO') it, ',', normr/normb
   end if

   deallocate(p)
end function cg

function mv( v, lambda, me, np_row, np_col) result(w)
   ! Matrix vector multiplication: w = lambda v + A v
   ! in which A represents -Laplacian with Neumann boundary conditions
   !
   implicit none

   real(kind=8), intent(in)    :: lambda
   real(kind=8), intent(in)    :: v(:,:)[*]
   integer, intent(in)         :: np_row, np_col, me

   real(kind=8), allocatable   :: v_flattened(:)[:]
   real(kind=8)                :: w(size(v,1),size(v,2))
   real(kind=8)                :: v_halo(0:size(v,1)+1,0:size(v,2)+1)
   integer                     :: i, j, nx, ny
   ! Add extra row/column around the domain to account for boundary conditions
   ! or data on other processor. Note that these points are outside the domain.
   ! They can be eliminated by the homogeneous Neumann condition. This condition
   ! says that the points on a row/column outside the boundary are equal to the 
   ! neighbouring points on the row/column inside the boundary.

   nx = size(v,1)
   ny = size(v,2)

   v_halo = 0.
   v_halo(1 : nx, 1: ny) = v  ! Points in the interior


   ! Because Fortran is slow with using 2D matrices, we will cast our coarray to a 1D array
   allocate(v_flattened(nx*ny)[*])
   v_flattened = 0.
   v_flattened = reshape(v, [nx*ny])

   ! Set the neumann boundary conditions
   ! Left boundary
   if (mod(me-1, np_row) == 0) then
      ! Set imaginary boundary if you have no left neighbour
      v_halo(0,1:ny) = v_halo(1,1:ny)
   end if

   ! Right boundary
   if (mod(me, np_row) == 0) then
      v_halo(nx+1, 1:ny) = v_halo(nx,1:ny)
   end if

   ! Top boundary
   if (me > np_row*np_col - np_row) then
      v_halo(1:nx, ny+1) = v_halo(1:nx, ny)
   end if
   
   ! Bottom boundary
   if (me <= np_row) then
      v_halo(1:nx,0) = v_halo(1:nx,1)
   end if 
   sync all

   ! Go across each processor and add a virtual boundary or retreive data from a neighbouring processor.
   ! Left boundary is given by v_flattened(1 : (nx-1)*ny + 1 : nx)
   ! Right boundary is given by v_flattened(nx:nx*ny:nx)
   ! Top boundary is given by v_flattened(size(v_flattened) - nx + 1 :)
   ! Bottom boundary is given by v_flattened(1 : nx )

   ! set LEFT boundary
   if (mod(me-1, np_row) /= 0) then
      v_halo(0, 1:ny) = v_flattened(nx:nx*ny:nx)[me-1]
   end if


   ! set RIGHT boundary
   if (mod(me, np_row) /= 0) then      
      v_halo(nx+1, 1:ny) = v_flattened(1 : (nx-1)*ny + 1 : nx)[me+1]
   end if
   
   ! set TOP boundary
      ! Remember, np_row*np_col = number of processors
   if (me <= np_row*np_col - np_row) then
      v_halo(1 : nx, ny+1) = v_flattened(1 : nx )[me + np_row]
   end if

   ! set BOTTOM boundary
   if (me > np_row) then
      v_halo( 1: nx, 0 ) = v_flattened(size(v_flattened) - nx + 1 :)[me - np_row]
   end if
   sync all

   ! Now we can just perform the stencil operation:
   w = lambda*v(:,:)
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
   integer                                  :: np
   real(kind=8), dimension(:), allocatable :: x_flat, y_flat
   real(kind=8) :: s

   ! Allocate the 1D arrays to hold flattened data
   allocate(x_flat(size(x)))
   allocate(y_flat(size(y)))

   ! Flatten the 2D arrays
   x_flat = reshape(x, [size(x)])
   y_flat = reshape(y, [size(y)])

   ! If the amount of processors is a power of 2 use butterfly,
   ! Else, use allgather
   np = num_images()

   if (np /= 2**(log(dble(np))/log(2.d0))) then
      s = dot_allgather(x_flat, y_flat)
   else
      s = dot_butterfly(x_flat, y_flat)
   end if

   dot_butterfly_2d = s

   deallocate(x_flat, y_flat)
end function dot_butterfly_2d

real(8) function dot_allgather(x, y)

implicit none
real(kind=8), dimension(:), intent(in) :: x, y
real(kind=8), save, allocatable :: local_dots(:)[:]
integer :: me, np, i

np = num_images()
me = this_image()

if (.not. allocated(local_dots)) then
    allocate(local_dots(np)[*])
end if

if (size(local_dots) /= np) then
    allocate(local_dots(np)[*])
end if

local_dots(me) = dot_product(x, y)

do i=1,np
    local_dots(me)[i] = local_dots(me)
end do

sync all

dot_allgather = sum(local_dots)

end function dot_allgather

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