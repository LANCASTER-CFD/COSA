!-----------------------------------------------------------------------
      subroutine copy_array(d4l_v,d4u_v,d5l_v,d5u_v,off_v,nl,v1, v2, &
                            d4l_l,d4u_l,d5l_l,d5u_l,off_l)
!-----------------------------------------------------------------------

      use cosa_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,nl,iblock,iv,dim4a,dim5a
      integer(kind=cosa_int) &
        d1l_v,d1u_v,d2l_v,d2u_v,d3l_v,d3u_v,d4l_v,d4u_v,d5l_v,d5u_v, &
        d1l_l,d1u_l,d2l_l,d2u_l,d3l_l,d3u_l,d4l_l,d4u_l,d5l_l,d5u_l, &
        off_v, off_l
      real (kind=cosa_real) v1(*),v2(*)

      dim4a  = d4u_v - d4l_v + 1
      dim5a  = d5u_v - d5l_v + 1

      do iblock = 1,mynblocks
        imax = i_imax  (iblock,nl)
        jmax = j_jmax  (iblock,nl)
        kmax = k_kmax  (iblock,nl)

        if (off_v.eq. 3) then
          iv   = 1 + off_p3  (iblock,nl) * dim4a * dim5a
          d1l_v = -1
          d1u_v = imax+1
          d2l_v = -1
          d2u_v = jmax+1
          d3l_v = -1
          d3u_v = kmax+1
        else if (off_v.eq. 2) then
          iv   = 1 + off_p2  (iblock,nl) * dim4a * dim5a
          d1l_v =  0
          d1u_v = imax+1
          d2l_v =  0
          d2u_v = jmax+1
          d3l_v =  0
          d3u_v = kmax+1
        else if (off_v.eq. 1) then
          iv   = 1 + off_p1  (iblock,nl) * dim4a * dim5a
          d1l_v =  0
          d1u_v = imax
          d2l_v =  0
          d2u_v = jmax
          d3l_v =  0
          d3u_v = kmax
        else
          write(*,*) 'Error in add_array. Aborting!'
          stop
        end if

        if (off_l.eq.-1) then
          d1l_l =  1
          d1u_l = imax-1
          d2l_l =  1
          d2u_l = jmax-1
          d3l_l =  1
          d3u_l = kmax-1
        else if (off_l.eq. 0) then
          d1l_l =  1
          d1u_l = imax
          d2l_l =  1
          d2u_l = jmax
          d3l_l =  1
          d3u_l = kmax
        else if (off_l.eq. 1) then
          d1l_l =  0
          d1u_l = imax
          d2l_l =  0
          d2u_l = jmax
          d3l_l =  0
          d3u_l = kmax
        else if (off_l.eq. 2) then
          d1l_l =  0
          d1u_l = imax+1
          d2l_l =  0
          d2u_l = jmax+1
          d3l_l =  0
          d3u_l = kmax+1
        else if (off_l.eq. 3) then
          d1l_l =  -1
          d1u_l = imax+1
          d2l_l =  -1
          d2u_l = jmax+1
          d3l_l =  -1
          d3u_l = kmax+1
        else
          write(*,*) 'Error in add_array. Aborting!'
          stop
        end if
   
        call copy_barray(v1(iv),v2(iv), &
           d1l_v,d1u_v,d2l_v,d2u_v,d3l_v,d3u_v,d4l_v,d4u_v,d5l_v,d5u_v, &
           d1l_l,d1u_l,d2l_l,d2u_l,d3l_l,d3u_l,d4l_l,d4u_l,d5l_l,d5u_l)
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine copy_barray(vn,vo, &
           d1l_v,d1u_v,d2l_v,d2u_v,d3l_v,d3u_v,d4l_v,d4u_v,d5l_v,d5u_v, &
           d1l_l,d1u_l,d2l_l,d2u_l,d3l_l,d3u_l,d4l_l,d4u_l,d5l_l,d5u_l)
!---------------------------------------------------------------------       
       
      use cosa_precision

      implicit none

      integer(kind=cosa_int) &
        d1l_v,d1u_v,d2l_v,d2u_v,d3l_v,d3u_v,d4l_v,d4u_v,d5l_v,d5u_v, &
        d1l_l,d1u_l,d2l_l,d2u_l,d3l_l,d3u_l,d4l_l,d4u_l,d5l_l,d5u_l
      integer(kind=cosa_int) i1,i2,i3,i4,i5
      real (kind=cosa_real) &
        vn(d1l_v:d1u_v,d2l_v:d2u_v,d3l_v:d3u_v,d4l_v:d4u_v,d5l_v:d5u_v), &
        vo(d1l_v:d1u_v,d2l_v:d2u_v,d3l_v:d3u_v,d4l_v:d4u_v,d5l_v:d5u_v)

      do i5 = d5l_l,d5u_l
        do i4 = d4l_l,d4u_l
          do i3 = d3l_l,d3u_l
            do i2 = d2l_l,d2u_l
              do i1 = d1l_l,d1u_l
                vn(i1,i2,i3,i4,i5) = vo(i1,i2,i3,i4,i5)
              end do
            end do
          end do
        end do
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine add_array(d4l_v,d4u_v,d5l_v,d5u_v,off_v,nl,v1, v2, &
                           d4l_l,d4u_l,d5l_l,d5u_l,off_l,fact)
!-----------------------------------------------------------------------

      use cosa_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,nl,iblock,iv,dim4a,dim5a
      integer(kind=cosa_int) &
        d1l_v,d1u_v,d2l_v,d2u_v,d3l_v,d3u_v,d4l_v,d4u_v,d5l_v,d5u_v, &
        d1l_l,d1u_l,d2l_l,d2u_l,d3l_l,d3u_l,d4l_l,d4u_l,d5l_l,d5u_l, &
        off_v, off_l
      real (kind=cosa_real) v1(*),v2(*)
      real (kind=cosa_real) fact

      dim4a  = d4u_v - d4l_v + 1
      dim5a  = d5u_v - d5l_v + 1

      do iblock = 1,mynblocks
        imax = i_imax  (iblock,nl)
        jmax = j_jmax  (iblock,nl)
        kmax = k_kmax  (iblock,nl)

        if (off_v.eq. 3) then
          iv   = 1 + off_p3  (iblock,nl) * dim4a * dim5a
          d1l_v = -1
          d1u_v = imax+1
          d2l_v = -1
          d2u_v = jmax+1
          d3l_v = -1
          d3u_v = kmax+1
        else
          write(*,*) 'Error in add_array. Aborting!'
          stop
        end if

        if (off_l.eq.-1) then
          d1l_l =  1
          d1u_l = imax-1
          d2l_l =  1
          d2u_l = jmax-1
          d3l_l =  1
          d3u_l = kmax-1
        else if (off_l.eq. 1) then
          d1l_l =  0
          d1u_l = imax
          d2l_l =  0
          d2u_l = jmax
          d3l_l =  0
          d3u_l = kmax
        else
          write(*,*) 'Error in add_array. Aborting!'
          stop
        end if
   
        call add_barray(v1(iv),v2(iv),fact, &
           d1l_v,d1u_v,d2l_v,d2u_v,d3l_v,d3u_v,d4l_v,d4u_v,d5l_v,d5u_v, &
           d1l_l,d1u_l,d2l_l,d2u_l,d3l_l,d3u_l,d4l_l,d4u_l,d5l_l,d5u_l)
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine add_barray(v1,v2,fact, &
           d1l_v,d1u_v,d2l_v,d2u_v,d3l_v,d3u_v,d4l_v,d4u_v,d5l_v,d5u_v, &
           d1l_l,d1u_l,d2l_l,d2u_l,d3l_l,d3u_l,d4l_l,d4u_l,d5l_l,d5u_l)
!---------------------------------------------------------------------       
       
      use cosa_precision

      implicit none

      integer(kind=cosa_int) &
        d1l_v,d1u_v,d2l_v,d2u_v,d3l_v,d3u_v,d4l_v,d4u_v,d5l_v,d5u_v, &
        d1l_l,d1u_l,d2l_l,d2u_l,d3l_l,d3u_l,d4l_l,d4u_l,d5l_l,d5u_l
      integer(kind=cosa_int) i1,i2,i3,i4,i5
      real (kind=cosa_real) &
        v1(d1l_v:d1u_v,d2l_v:d2u_v,d3l_v:d3u_v,d4l_v:d4u_v,d5l_v:d5u_v), &
        v2(d1l_v:d1u_v,d2l_v:d2u_v,d3l_v:d3u_v,d4l_v:d4u_v,d5l_v:d5u_v)
      real (kind=cosa_real) fact

      do i5 = d5l_l,d5u_l
        do i4 = d4l_l,d4u_l
          do i3 = d3l_l,d3u_l
            do i2 = d2l_l,d2u_l
              do i1 = d1l_l,d1u_l
                v1(i1,i2,i3,i4,i5) = v1(i1,i2,i3,i4,i5) + &
                                  v2(i1,i2,i3,i4,i5) * fact
              end do
            end do
          end do
        end do
      end do

      return
      end

!-------------------------------------------------------------------------
      subroutine zero(d4l_v,d4u_v,d5l_v,d5u_v,off_v,nl,v1, &
                      d4l_l,d4u_l,d5l_l,d5u_l,off_l)
!-----------------------------------------------------------------------

      use cosa_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,nl,iblock,iv,dim4a,dim5a
      integer(kind=cosa_int) &
        d1l_v,d1u_v,d2l_v,d2u_v,d3l_v,d3u_v,d4l_v,d4u_v,d5l_v,d5u_v, &
        d1l_l,d1u_l,d2l_l,d2u_l,d3l_l,d3u_l,d4l_l,d4u_l,d5l_l,d5u_l, &
        off_v, off_l
      real (kind=cosa_real) v1(*)

      dim4a  = d4u_v - d4l_v + 1
      dim5a  = d5u_v - d5l_v + 1

      do iblock = 1,mynblocks
        imax = i_imax  (iblock,nl)
        jmax = j_jmax  (iblock,nl)
        kmax = k_kmax  (iblock,nl)

        if (off_v.eq. 3) then
          iv   = 1 + off_p3  (iblock,nl) * dim4a * dim5a
          d1l_v = -1
          d1u_v = imax+1
          d2l_v = -1
          d2u_v = jmax+1
          d3l_v = -1
          d3u_v = kmax+1
        else if (off_v.eq. 2) then
          iv   = 1 + off_p2  (iblock,nl) * dim4a * dim5a
          d1l_v = 0
          d1u_v = imax+1
          d2l_v = 0
          d2u_v = jmax+1
          d3l_v = 0
          d3u_v = kmax+1
        else if (off_v.eq. 1) then
          iv   = 1 + off_p1  (iblock,nl) * dim4a * dim5a
          d1l_v =  0
          d1u_v = imax
          d2l_v =  0
          d2u_v = jmax
          d3l_v =  0
          d3u_v = kmax
        else
          write(*,*) 'routine zero: variable index error. Aborting!'
          stop
        end if

        if (off_l.eq.-1) then
          d1l_l =  1
          d1u_l = imax-1
          d2l_l =  1
          d2u_l = jmax-1
          d3l_l =  1
          d3u_l = kmax-1
        else if (off_l.eq. 0) then
          d1l_l =  1
          d1u_l = imax
          d2l_l =  1
          d2u_l = jmax
          d3l_l =  1
          d3u_l = kmax
        else if (off_l.eq. 1) then
          d1l_l =  0
          d1u_l = imax
          d2l_l =  0
          d2u_l = jmax
          d3l_l =  0
          d3u_l = kmax
        else if (off_l.eq. 2) then
          d1l_l =  0
          d1u_l = imax+1
          d2l_l =  0
          d2u_l = jmax+1
          d3l_l =  0
          d3u_l = kmax+1
        else if (off_l.eq. 3) then
          d1l_l = -1
          d1u_l = imax+1
          d2l_l = -1
          d2u_l = jmax+1
          d3l_l = -1
          d3u_l = kmax+1
        else
          write(*,*) 'routine zero: loop index error. Aborting!'
          stop
        end if
   
        call zerob(v1(iv), &
           d1l_v,d1u_v,d2l_v,d2u_v,d3l_v,d3u_v,d4l_v,d4u_v,d5l_v,d5u_v, &
           d1l_l,d1u_l,d2l_l,d2u_l,d3l_l,d3u_l,d4l_l,d4u_l,d5l_l,d5u_l)
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine zerob(v, &
           d1l_v,d1u_v,d2l_v,d2u_v,d3l_v,d3u_v,d4l_v,d4u_v,d5l_v,d5u_v, &
           d1l_l,d1u_l,d2l_l,d2u_l,d3l_l,d3u_l,d4l_l,d4u_l,d5l_l,d5u_l)
!-----------------------------------------------------------------------
       
      use cosa_precision

      implicit none

      integer(kind=cosa_int) &
        d1l_v,d1u_v,d2l_v,d2u_v,d3l_v,d3u_v,d4l_v,d4u_v,d5l_v,d5u_v, &
        d1l_l,d1u_l,d2l_l,d2u_l,d3l_l,d3u_l,d4l_l,d4u_l,d5l_l,d5u_l
      integer(kind=cosa_int) i1,i2,i3,i4,i5
      real (kind=cosa_real) &
        v(d1l_v:d1u_v,d2l_v:d2u_v,d3l_v:d3u_v,d4l_v:d4u_v,d5l_v:d5u_v)

      do i5 = d5l_l,d5u_l
        do i4 = d4l_l,d4u_l
          do i3 = d3l_l,d3u_l
            do i2 = d2l_l,d2u_l
              do i1 = d1l_l,d1u_l
                v(i1,i2,i3,i4,i5) = 0.d0
              end do
            end do
          end do
        end do
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine zeroall(dim,v)
!---------------------------------------------------------------------       
       
      use cosa_precision

      implicit none

      integer*8 i,dim
      real (kind=cosa_real) v(dim)

      do i = 1,dim
        v(i) = 0.d0
      end do

      return
      end
