!*****************************************************************************************
!>
!  Conversion of raster data to GIF89 format.
!
!# See also
!   * The original code (Licence: public domain) was from
!     [here](http://fortranwiki.org/fortran/show/writegif)
!
!# History
!   * Version 1.01, August 1999, Written by Jos Bergervoet
!   * 2008 Jan 28: Modified by Clive Page to use stream I/O, array as colourmap.
!   * Jacob Williams, 7/27/2014. Refactored, updated, added ability to export animated gifs.
!
!*****************************************************************************************

    module gif_module

    implicit none

    private
    
    public :: write_animated_gif

    contains
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!  date: 7/27/2014
!
!  Writes gif89 image, given pixel array and color map
!  This version can create an animated gif:
!
!  * The pixel matrix is rank 3: image i is pixel(i,:,:)
!  * If size(pixel,1) is 1, then a regular gif is produced.
!
!# See also
!   1. [writegif](http://fortranwiki.org/fortran/show/writegif)
!   2. [GIF format](http://www.onicos.com/staff/iz/formats/gif.html#aeb)
!   3. [GIF File Format Summary](http://www.fileformat.info/format/gif/egff.htm)

    subroutine write_animated_gif(filename,pixel,colormap,transparent,delay)
    
    character(len=*),intent(in)         :: filename    !! file to create or replace
    integer,intent(in),dimension(:,:,:) :: pixel       !! pixel values [0 to ncol]
    integer,intent(in),dimension(:,0:)  :: colormap    !! [r,g,b (0:255)] , [0:ncol] colors
    integer,intent(in),optional         :: transparent !! transparent color
    integer,intent(in),optional         :: delay       !! delay time [1/100 of seconds]

    integer,parameter     :: bufend=260
    character(len=bufend) :: buf
    integer               :: ibuf ! output buffer vars
    integer,parameter     :: maxcode = 4095 
    integer,parameter     :: nocode = maxcode+1 !! definitions for lzw

    ! define lzw code tables for hashing:
    !
    ! for any code p, which codes for a sequence af pixel-values, endbyte(p)
    ! is the last pixel-value, follow(p) points to another code (if it exists)
    ! which codes for this same sequence, but with one more pixel-value
    ! appended.
    !   for each code p, next(p) points to another code which codes for a
    ! similar sequence with only the endbyte different. this is a hashing
    ! pointer, for fast look-up.
    !   all pointers are 'nocode' if they point to nothing
    !
    character(len=1),dimension(0:maxcode+1) :: endbyte
    integer,dimension(0:maxcode) :: follow, next
    integer :: ncod, curmaxcode, eoi, cc, p, k, child, &
                    maxbase, skip, slen, blen, accum, nout

    integer :: infobyte,nx,ny,cblen,hasmap,maxincol,istat,&
                maxgifcol,background,i,iunit,iframe,n,dt
    character(len=1),dimension(2) :: t
    
    !delay time:
    if (present(delay)) then
        dt = delay
    else
        dt = 1
    end if
    
    !transparency info:
    if (present(transparent)) then
        t(1) = char(1) !Reserved+Disposal Method+User Input Flag+Transparent Color Flag
        t(2) = char(transparent)
    else
        t(1) = char(0)
        t(2) = char(0)
    end if
                
    open(    newunit=iunit,&
            file=trim(filename),&
            access='STREAM',&
            status='REPLACE',&
            iostat=istat)
    
    if (istat==0) then

        n = size(pixel,1) !number of images
        nx = ubound(pixel, 2)
        ny = ubound(pixel, 3)
        maxincol = ubound(colormap,2) 

        do i=1,8                            ! find the bitsize, blen, for pixels
            blen = i
            maxgifcol = 2**blen - 1         ! Number of colors has to be power of 2 
            if (maxgifcol>=maxincol) exit   ! now blen and maxgifcol are correct
                                            ! [only up to 256 colors]
        end do
    
        !------------
        ! GIF Header
        !------------
    
        write(iunit) 'GIF89a'
    
        ! create information for screen descriptor
        if (present(transparent)) then
            background = transparent
        else
            background = 0
        end if
        hasmap = 1
        cblen = blen
        infobyte = hasmap * 128 + (cblen-1) * 16 + blen-1
    
        ! write the screen descriptor
        write(iunit)   char2(nx),&         ! logical screen width
                        char2(ny),&         ! logical screen height
                        char(infobyte),&    ! screen and color map information
                        char(background),&  ! background color index
                        char(0)             ! pixel aspect ratio
    
        ! write global colormap
        do i=0,maxgifcol
            write(iunit)  char(colormap(1,min(i,maxincol))), &
                          char(colormap(2,min(i,maxincol))), &
                          char(colormap(3,min(i,maxincol)))
        end do
    
        if (n>1) then    !it is an animated gif
    
            !-----------------------------
            ! Application Extension Block
            !-----------------------------
    
            ! See: http://odur.let.rug.nl/kleiweg/gif/netscape.html
    
            write(iunit)   '!',&           ! Extension Introducer (0x21)
                            char(255),&     ! GIF Extension code
                            char(11),&      ! Length of Application Block
                            'NETSCAPE',&    ! Application Identifier
                            '2.0',&         ! Application Authentication Code
                            char(3),&       ! Length of Data Sub-Block
                            char(1),&        ! 1 (0x01)
                            char(0),&       ! number of loop iterations
                            char(0),&       ! Data Sub-Block Terminator
                            char(0)         ! Block Terminator (0x00)
        
        end if
    
        !each frame of the animated gif:
        do iframe = 1,n
                
            !---------------------------------
            ! Graphic Control Extension Block
            !---------------------------------
        
            write(iunit)   '!',&           ! Extension Introducer (0x21)
                            char(249),&    ! Graphic Control Label (0xF9)
                            char(4),&      ! Block Size (0x04)
                            t(1),&         !
                            char2(dt),&    ! Delay Time
                            t(2),&         ! Transparent Color Index
                            char(0)        ! Block Terminator (0x00)

            !-------------
            ! Image Block
            !-------------

            ! now create and write image descriptor
            hasmap = 0
            infobyte = hasmap * 128 + blen-1    ! add 64, if interlaced
 
            write(iunit)   ',',&            ! Image Separator (0x2C)
                            char2(0),&      ! Image Left Position
                            char2(0),&      ! Image Top Position
                            char2(nx),&     ! Image Width
                            char2(ny),&     ! Image Height
                            char(infobyte)  ! Image and Color Table Data Information
                
            call giflzw(iunit,pixel(iframe,:,:))  ! now the raster data

            write(iunit) char(0)    ! Block Terminator (0x00)
            
        end do
        
        !---------
        ! Trailer
        !---------

        write(iunit) ';'
    
    else
        write(*,*) 'Error opening :'//trim(filename)
    end if
    
    !close the gif file:
    close(unit=iunit,iostat=istat)
    
    contains
!*****************************************************************************************
    
    !*************************************************************************************
    !>
    !  Convert the two least sig bytes of an integer to a 2-character string

        function char2(ival) result(c)
    
        integer, intent(in) :: ival
        character(len=2)    :: c
    
        c = achar(mod(ival,256)) // achar(mod(ival/256,256))
    
        end function char2
    !*************************************************************************************

    !*************************************************************************************
    !>
    !  Flushes up to 255 bytes to output file if buffer contains data, keeping
    !  rest of data in buffer. If skip>0 there is a partially filled last byte
    !  in buf[ibuf]. This byte will be written only if ibuf<256. That should be
    !  the last call to flushbuffer.

        subroutine flushbuffer(iunit)

        integer, intent(in) :: iunit   !! i/o unit to use
    
        integer :: bl   !! number of bytes to write (to be determined)

        if (ibuf > 255) then     ! we will write buf[1..255]
            bl = 255
        else if (skip /= 0) then ! buf[ibuf] is partially used, write buf[1..ibuf]
            bl = ibuf
        else if (ibuf > 1) then  ! write buf[1..ibuf-1], there is no partial byte
            bl = ibuf-1
        else                     ! nothing to write
            return
        end if

        write(iunit) char(bl)
        write(iunit) buf(1:bl)
        buf(1:ibuf-bl) = buf(bl+1:ibuf) ! shift down remaining data
        ibuf = ibuf - bl

        end subroutine flushbuffer
    !*************************************************************************************

    !*************************************************************************************
    !>
    !  routine for LZW coding

        subroutine giflzw(iunit, pixel)          

        integer, intent(in)                 :: iunit
        integer, intent(in), dimension(:,:) :: pixel
        integer                             :: i
        integer                             :: j

        nout=0                        ! for counting the codes going out
        if (blen<2) then
            blen=2                    ! pixel code-length, 2 is minimum for gif
        end if
        write(iunit) char(blen)
        maxbase = 2**blen - 1
        call inittable()
        call slicewrite(iunit, cc)

        do j=1, ubound(pixel,2)
            do i=1, ubound(pixel,1)
                k = modulo(pixel(i,j), maxbase+1)    ! take next byte, prevent overflow
                if (i==1 .and. j==1) then
                    p = k                       ! first raster byte has one-byte code p
                    cycle                       ! for the first byte no further action
                end if
                                            ! now see if code exists for sequence [.p.]k
                child = follow(p)           ! [.p.]k is "string coded by p" followed by k
                childloop: do
                    if ((child == nocode) .or. (ichar(endbyte(child)) == k)) then
                        exit childloop
                    end if
                    child = next(child)
                end do childloop

                if (child /= nocode) then    ! if code for [.p.]k was found, store it in p
                    p = child
                else                         ! if not: output p and create code for [.p.]k
                    call slicewrite(iunit, p)
                    if (ncod > maxcode) then        ! check if a new code can be added
                        call slicewrite(iunit, cc)  ! if not: tell listener to clear table
                        call inittable()            ! and clear our own table
                    else
                        if (ncod > curmaxcode) then
                            slen = slen+1               ! new codes will be one bit longer
                            curmaxcode = curmaxcode * 2 + 1 ! and more codes are possible
                        end if
                        endbyte(ncod) = char(k)   ! ncod is the new code for [.p.]k
                        follow(ncod) = nocode
                        next(ncod) = follow(p)    ! include ncod in the hashing list
                        follow(p) = ncod          !     of codes with same start-sequence
                        ncod = ncod+1
                    end if
                    p = k
                end if
            end do
        end do
    
        call slicewrite(iunit, p)      ! send the last code to buffer
        call slicewrite(iunit, eoi)    ! send 'end of image' to buffer
        call flushbuffer(iunit)        ! extra flush, including partial last byte

        end subroutine giflzw
    !*************************************************************************************

    !*************************************************************************************
    !>
    !  Initialize table.

        subroutine inittable()

        integer :: i

        do i=0,maxbase                  ! start with defining the codes 0..maxbase
            endbyte(i) = char(i)        ! for one-pixel sequences (code=pixelvalue)
        end do                          ! initially no multi-pixel codes exist
        follow(0:maxbase) = nocode
        next(0:maxbase) = nocode
        cc = maxbase+1                  ! 'clear code-tabel', a control code
        eoi = maxbase+2                 ! 'end of image', another control code
        ncod = cc + 2                   ! ncod = number of currently defined codes
        slen = blen + 1                 ! current number of bits to write one code
        curmaxcode = 2**slen - 1        ! currently the highest, until slen increases

        end subroutine inittable
    !*************************************************************************************

    !*************************************************************************************
    !>
    !  add some bits (a 'slice') to output buffer

        subroutine slicewrite(iunit, code) 
    
        integer, intent(in)  :: iunit
        integer, intent(in)  :: code

        if (nout == 0) then             ! initiate output buffer
            ibuf = 1
            skip = 0
            accum = 0
        end if
        nout = nout+1

        accum = accum + code * 2**skip  ! add bits at correct position in accum
        skip = skip + slen              ! slen is current slice length, in bits

        shiftout: do
            buf(ibuf:ibuf) = char(modulo(accum,256))
            if (skip<8) exit shiftout
            ibuf = ibuf+1               ! last written buffer-byte is now permanent
            accum = accum/256           ! remove that byte from accum
            skip = skip-8               ! skip points to next bit to write in accum
        end do shiftout

        if (ibuf>255) then
            call flushbuffer(iunit)    ! won't write unfinished byte in buf[ibuf]
        end if
    
        ! at most 255 bytes will be left in buffer
    
        end subroutine slicewrite
    !*************************************************************************************

    end subroutine write_animated_gif
!*****************************************************************************************

!*****************************************************************************************
    end module gif_module
!*****************************************************************************************
