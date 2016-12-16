module FILTER
    implicit none
    real(kind=8), dimension(3,3), protected :: lpf_filter_6
    real(kind=8), dimension(3,3), protected :: lpf_filter_9
    real(kind=8), dimension(3,3), protected :: lpf_filter_10
    real(kind=8), dimension(3,3), protected :: lpf_filter_16
    real(kind=8), dimension(3,3), protected :: lpf_filter_32
    real(kind=8), dimension(3,3), protected :: hpf_filter_1
    real(kind=8), dimension(3,3), protected :: hpf_filter_2
    real(kind=8), dimension(3,3), protected :: hpf_filter_3

contains
    subroutine filter_image(the_image, out_image, nrows, ncols, &
        bits_per_pixel, type, low_high)
        implicit none
        integer, intent(in) :: nrows, ncols, bits_per_pixel, type
        real(kind=8), intent(in), dimension(nrows, ncols) :: the_image
        real(kind=8), intent(out), dimension(nrows, ncols) :: out_image
        real(kind=8), dimension(3,3) :: filter
        character(len=*), intent(in) :: low_high

        integer :: a, b, d, i, j, k, length, width, imax, isum
        call setup_filters(type, low_high, filter)
        d = type
        imax = 255
        if(type == 2 .or. type==3) then
            d = 1
        endif
        if(bits_per_pixel == 4) then
            imax= 16
        endif
        !         /* Do convolution over image array */
        print*, ""
        do j=2, ncols
            do i=2, nrows
                isum = 0.0d0
                do a=0,2
                    do b=0,2
                        isum = isum + the_image(i+a,j+b)*filter(a+1, b+1)
                    enddo
                enddo
                isum = isum/d
                if (isum < 0) then
                    isum = 0
                elseif (isum > imax) then
                    isum = imax
                endif
                out_image(i, j) = isum
            enddo
        enddo
    end subroutine filter_image

    subroutine adaptive_median_filter(the_image, out_image, nrows,&
        ncols, mask_size, noise_variance, cutoff)
       implicit none
        integer, intent(in) :: nrows, ncols, mask_size
        real(kind=8), intent(in), dimension(nrows, ncols) :: the_image
        real(kind=8), intent(out), dimension(nrows, ncols) :: out_image
        real(kind=8), intent(in) :: noise_variance, cutoff
        integer :: a, b, icount, i, j, k, length, sd2, sd2p1, ss, width, avg, tmp, var
        real(kind=8), dimension(mask_size*mask_size) :: elements
!        real(kind=8), parameter :: cutoff = 3.0D4
        sd2   = mask_size/2
!        sd2p1 = sd2 + 1
        ss = mask_size*mask_size
        !   Loop over image array
        do j=sd2+1, ncols-sd2
            do i=sd2+1, nrows-sd2
                icount=0
                do a=-sd2, sd2
                    do b=-sd2, sd2
                        icount = icount + 1
                        elements(icount) = the_image(i+a,j+b)
                    enddo
                enddo
                avg = sum(elements)/dble(icount)
                var = sum((elements - avg)**2)/dble(icount)
                tmp = noise_variance/var
                if (tmp > 1.0d0) then
                    tmp = 1.0d0
                endif
                if (avg < the_image(i, j) .and. the_image(i, j) < cutoff) then
                   out_image(i, j) = 0.0d0
                else
                    out_image(i, j) = the_image(i, j) - tmp*(the_image(i, j) - avg)
                endif
            enddo
        enddo
   end subroutine adaptive_median_filter

   subroutine median_filter(the_image, out_image,nrows, ncols, mask_size)
        implicit none
        integer, intent(in) :: nrows, ncols, mask_size
        real(kind=8), intent(in), dimension(nrows, ncols) :: the_image
        real(kind=8), intent(out), dimension(nrows, ncols) :: out_image
        integer :: a, b, icount, i, j, k, length, sd2, sd2p1, ss, width
        real(kind=8), dimension(mask_size*mask_size) :: elements
        sd2   = mask_size/2
!!!        sd2p1 = sd2 + 1
        ss = mask_size*mask_size
!!!        Loop over image array
        do j=sd2+1, ncols-sd2
            do i=sd2+1, nrows-sd2
                icount=0
                do a=-sd2, sd2
                    do b=-sd2, sd2
                        icount = icount + 1
                        elements(icount) = the_image(i+a,j+b)
                    enddo
                enddo
                !!out_image(i, j) = sum(the_image)/dble(icount)
                out_image(i, j) = getMedian(elements, ss)
            enddo
        enddo
    end subroutine median_filter

    function getMedian(arr, arrsize) result(median)
        use MERGE_RANK, only: mrgrnk
        implicit none
        integer, intent(in) :: arrsize
        real(kind=8), intent(inout), dimension(arrsize) :: arr
        real(kind=8) :: median
        integer, dimension(arrsize) :: arrindx
        call mrgrnk(arr, arrindx)
        arr = arr(arrindx)
        if(mod(arrsize,2)==1) then
            median = arr(arrsize/2)
        else
            median = (arr(arrsize/2) + arr(arrsize/2+1))/2
        endif
    end function getMedian

    subroutine setup_filters(filter_type, low_high, filter)
        implicit none
        integer, intent(in) :: filter_type
        real(kind=8), intent(out), dimension(3,3) :: filter
        character(len=*), intent(in) :: low_high
        character(len=1) :: lowHigh
        integer, save :: i,j
        logical :: init=.True.
        if(init) then
            call setup_masks()
            init = .false.
        endif
        lowHigh = Adjustl(trim(low_high))
        if(lowHigh == "l" .or. lowHigh=="L") then
            print*, "Low pass"
            if (filter_type == 6) then
                filter = lpf_filter_6
            elseif (filter_type == 9) then
                filter = lpf_filter_9
            elseif (filter_type == 10) then
                filter = lpf_filter_10
            elseif (filter_type == 16) then
                filter = lpf_filter_16
            elseif (filter_type == 32) then
                filter = lpf_filter_32
            else
                print*, "Error! Filter type is not recognized."
                stop
            endif
        elseif(lowHigh == "h" .or. lowHigh=="H") then
            print*, "High pass"
            if (filter_type == 1) then
                filter = hpf_filter_1
            elseif(filter_type == 2) then
                filter = hpf_filter_2
            elseif(filter_type == 3) then
                filter = hpf_filter_3
            else
                print*, "Filter type is not recognized."
            endif
        else
            print*, "Filter type is not recognized."
        endif
    end  subroutine setup_filters

    subroutine setup_masks()
        implicit none
        lpf_filter_6  = reshape((/ 0, 1, 0, 1, 2,  1, 0, 1, 0 /),(/3,3/))
        lpf_filter_9  = reshape((/ 1, 1, 1, 1, 1,  1, 1, 1, 1 /),(/3,3/))
        lpf_filter_10 = reshape((/ 1, 1, 1, 1, 2,  1, 1, 1, 1 /),(/3,3/))
        lpf_filter_16 = reshape((/ 1, 2, 1, 2, 4,  2, 1, 2, 1 /),(/3,3/))
        lpf_filter_32 = reshape((/ 1, 4, 1, 4, 12, 4, 1, 4, 1 /),(/3,3/))
        hpf_filter_1  = reshape((/ 0,-1, 0,-1, 5, -1, 0,-1, 0 /),(/3,3/))
        hpf_filter_2  = reshape((/-1,-1,-1,-1, 9, -1,-1,-1,-1 /),(/3,3/))
        hpf_filter_3  = reshape((/ 1,-2, 1,-2, 5, -2, 1,-2, 1 /),(/3,3/))
    end subroutine setup_masks

    subroutine print_mask(mask, myunit)
        implicit none
        integer :: i, n
        integer, intent(in) :: myunit
        real(kind=8), intent(in), dimension(:,:) :: mask
        n=size(mask(:,1))
        write(myunit, *) ""
        write(myunit, "(A, I0)") "Mask : ", nint(sum(mask))
        if(all(mask==0)) then
            write(myunit, *) "WARNING: Mask not initialized. Please initialize mask before using."
        endif
        do i=1, n
            write(myunit, "(3(f3.0, 2x))") mask(i, :)
        enddo
    end subroutine print_mask

end module FILTER
