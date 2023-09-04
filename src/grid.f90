!-----------------------------------------------------------------------
      subroutine coordex(nl,x,y,z)
!-----------------------------------------------------------------------
!     extend grid to enable calculation of inviscid metrix of auxiliary
!     cells, required to compute viscous metrix.
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,iblock,imax,jmax,kmax,ixyz
      real(kind=cosa_real) x(*),y(*),z(*)

      do iblock = 1,mynblocks
        imax  = i_imax  (iblock,nl)
        jmax  = j_jmax  (iblock,nl)
        kmax  = k_kmax  (iblock,nl)
        ixyz  = 1 + off_p2(iblock,nl)         * dim5h
        call cordex_b(x(ixyz),y(ixyz),z(ixyz),imax,jmax,kmax,nharms)
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine cordex_b(x,y,z,imax,jmax,kmax,nharms)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,lmet,nharms
      integer(kind=cosa_int) i,j,k,n
      real(kind=cosa_real) &
          x(0:imax+1,0:jmax+1,0:kmax+1,0:2*nharms*hbmove), &
          y(0:imax+1,0:jmax+1,0:kmax+1,0:2*nharms*hbmove), &
          z(0:imax+1,0:jmax+1,0:kmax+1,0:2*nharms*hbmove)
      real(kind=cosa_real) voldot

      do n=0,2*nharms*hbmove

        do k=1,kmax
          do i = 1,imax
            x(i,0     ,k,n) = x(i,1   ,k,n)
            y(i,0     ,k,n) = y(i,1   ,k,n)
            z(i,0     ,k,n) = z(i,1   ,k,n)
            x(i,jmax+1,k,n) = x(i,jmax,k,n)
            y(i,jmax+1,k,n) = y(i,jmax,k,n)
            z(i,jmax+1,k,n) = z(i,jmax,k,n)
          end do
          do j = 1,jmax
            x(0     ,j,k,n) = x(1   ,j,k,n)
            y(0     ,j,k,n) = y(1   ,j,k,n)
            z(0     ,j,k,n) = z(1   ,j,k,n)
            x(imax+1,j,k,n) = x(imax,j,k,n)
            y(imax+1,j,k,n) = y(imax,j,k,n)
            z(imax+1,j,k,n) = z(imax,j,k,n)
          end do
        end do

        do j=1,jmax
          do i = 1,imax
            x(i,j,0     ,n) = x(i,j,1   ,n)
            y(i,j,0     ,n) = y(i,j,1   ,n)
            z(i,j,0     ,n) = z(i,j,1   ,n)
            x(i,j,kmax+1,n) = x(i,j,kmax,n)
            y(i,j,kmax+1,n) = y(i,j,kmax,n)
            z(i,j,kmax+1,n) = z(i,j,kmax,n)
          end do
        end do

      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine volex(nl,vol)
!-----------------------------------------------------------------------
!     extend volume array to auxiliary cells  to enable calculation of 
!     viscous metrix, and/or correct use of MG prolongation
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,iblock,imax,jmax,kmax,ivol
      real(kind=cosa_real) vol(*)

      do iblock = 1,mynblocks
        imax  = i_imax  (iblock,nl)
        jmax  = j_jmax  (iblock,nl)
        kmax  = k_kmax  (iblock,nl)
        ivol  = 1 + off_p1(iblock,nl)
        call volex_b(vol(ivol),imax,jmax,kmax)
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine volex_b(vol,imax,jmax,kmax)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax
      integer(kind=cosa_int) iblock,i,j,k,n
      real (kind=cosa_real) vol(0:imax,0:jmax,0:kmax)

      do k=1,kmax-1
        do i = 1,imax-1
          vol(i,0,k) = vol(i,1,k)
          vol(i,jmax,k) = vol(i,jmax-1,k)
        end do
        do j = 1,jmax-1
          vol(0,j,k) = vol(1,j,k)
          vol(imax,j,k) = vol(imax-1,j,k)
        end do
      end do

      do j = 1,jmax-1
        do i = 1,imax-1
          vol(i,j,0) = vol(i,j,1)
          vol(i,j,kmax) = vol(i,j,kmax-1)
        end do
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine cutman_x(nl,cutopo,x,y,z)
!-----------------------------------------------------------------------

      use parallelutils
      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,iblk1,icut,iblk2,imax1,jmax1,kmax1,imax2,jmax2,kmax2, &
                ixyz1,ixyz2
      integer(kind=cosa_int) cutopo(21,myncuts),a,b,c,i,reqnum,n
      logical myblock1, myblock2
      real(kind=cosa_real) x(*),y(*),z(*)
      integer(kind=cosa_int) sendrequestnum,receiverequestnum
      integer(kind=cosa_int) i_id,i_in_id
      integer(kind=cosa_int) sendid,recvid

      numrecv = 0  ! array syntax
      numsend = 0  ! array syntax
      sendrecv_datasize = dim5h

      allocate(sendarray(sendrecv_datasize, mrequest_perproc, &
                         maxproc_cutman))
      allocate(receivearray(sendrecv_datasize, mrequest_perproc, &
                            maxproc_cutman))

      do icut = 1,myncuts
        iblk1 = cutopo( 1,icut)
        iblk2 = cutopo(10,icut)
        call myblock(iblk1,myblock1,.true.)
        call myblock(iblk2,myblock2,.true.)
        if(myblock1) then
          imax1  = i_imax     (iblk1,nl)
          jmax1  = j_jmax     (iblk1,nl)
          kmax1  = k_kmax     (iblk1,nl)
          ixyz1  = 1 + off_p2 (iblk1,nl) * dim5h
        else
          ixyz1  = 1
        end if
        if(myblock2) then
          imax2  = i_imax     (iblk2,nl)
          jmax2  = j_jmax     (iblk2,nl)
          kmax2  = k_kmax     (iblk2,nl)
          ixyz2  = 1 + off_p2 (iblk2,nl) * dim5h
        else
          ixyz2  = 1
        end if
        call cut_x(imax1,jmax1,kmax1,x(ixyz1),y(ixyz1),z(ixyz1),imax2, &
             jmax2,kmax2,x(ixyz2),y(ixyz2),z(ixyz2),nharms,dim5h, &
             cutopo(1,icut),receivearray,receiverequests,numrecv, &
             receiveindices,sendarray,sendrequests,numsend,ixyz1, &
             mrequest_perproc,maxproc_cutman,maxproc,myrecvid_num, &
             mysendid_num)
      end do

      do i_id=1,maxrecvid_num
        recvid = myrecvid(i_id)
        call receiveblockdata_fromid(receivearray(1,1,i_id),recvid, &
             numrecv(i_id)*sendrecv_datasize,receiverequests(i_id))
      enddo
      do i_id=1,maxsendid_num
        sendid = mysendid(i_id)
        call sendblockdata_toid(sendarray(1,1,i_id),sendid, &
             numsend(i_id)*sendrecv_datasize,sendrequests(i_id))
      enddo
      do i_id=1,maxrecvid_num
        call waitanymessages(receiverequests,maxrecvid_num,reqnum)
        do i_in_id=1,numrecv(reqnum)
          a     = receiveindices(i_in_id,1,reqnum)
          b     = receiveindices(i_in_id,2,reqnum)
          c     = receiveindices(i_in_id,3,reqnum)
          ixyz1 = receiveindices(i_in_id,4,reqnum)
          imax1 = receiveindices(i_in_id,6,reqnum)
          jmax1 = receiveindices(i_in_id,7,reqnum)
          kmax1 = receiveindices(i_in_id,8,reqnum)
          if (receiveindices(i_in_id,5,reqnum) .eq. 0) then
            call copyxreceivedata(imax1,jmax1,kmax1,a,b,c,x(ixyz1), &
                 receivearray,i_in_id,reqnum)
          else if (receiveindices(i_in_id,5,reqnum) .eq. 1) then
            call copyxreceivedata(imax1,jmax1,kmax1,a,b,c,y(ixyz1), &
                 receivearray,i_in_id,reqnum)
          else if (receiveindices(i_in_id,5,reqnum) .eq. 2) then
            call copyxreceivedata(imax1,jmax1,kmax1,a,b,c,z(ixyz1), &
                 receivearray,i_in_id,reqnum)
          endif
        end do
      end do

      call waitallmessages(sendrequests,maxsendid_num)

      deallocate(sendarray, receivearray)

      return

      end
!-----------------------------------------------------------------------
      subroutine cut_x(imax1,jmax1,kmax1,x1,y1,z1,imax2,jmax2,kmax2,x2, &
        y2,z2,nharms,dim5h,cutopo,receivearray,receiverequests,numrecv, &
        receiveindices,sendarray,sendrequests,numsend,ixyz1, &
        mrequest_perproc,maxproc_cutman,maxproc,myrecvid_num, &
        mysendid_num)
!-----------------------------------------------------------------------

!     Routine to do cut boundary condition.
!     Flow data stored in Q1 are updated from data in the interior of Q2.
!      * sendrequests/receivrequests have indexes which is the process I am communicating to
!      * sendarray/recvarray have three components:
!        - index spanning q or mut component to send/recv
!        - index of point to send/recv
!        - index of process I have to send/recv to/from
!      * receiveindices have three components:
!        - index of point to send/recv
!        - index of information to trasmit (1-6)
!        - index of process to send/recv to/from

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax1,jmax1,kmax1,imax2,jmax2,kmax2,nharms,dim5h,iblk1, &
                iblk2
      integer(kind=cosa_int) cutopo(21),ixyz1,mrequest
      integer(kind=cosa_int) ijkmax1(3),ijkmax2(3),istr1(3),iend1(3),istr2(3), &
           iend2(3),isgn1(3),isgn2(3),iord(3),len(3),idir1,idir2,inout1, &
           inout2,ibcpt,ibcinc,ibcpt2,inr,inrinc,inr2,l,n,ic1,ic2,ic3, &
           jc1,jc2,jc3,i2,i3,ii1,jj1,kk1,in1,jn1,kn1
      real(kind=cosa_real) &
           x1(0:imax1+1,0:jmax1+1,0:kmax1+1,0:2*nharms*hbmove), &
           y1(0:imax1+1,0:jmax1+1,0:kmax1+1,0:2*nharms*hbmove), &
           z1(0:imax1+1,0:jmax1+1,0:kmax1+1,0:2*nharms*hbmove), &
           x2(0:imax2+1,0:jmax2+1,0:kmax2+1,0:2*nharms*hbmove), &
           y2(0:imax2+1,0:jmax2+1,0:kmax2+1,0:2*nharms*hbmove), &
           z2(0:imax2+1,0:jmax2+1,0:kmax2+1,0:2*nharms*hbmove)
      integer(kind=cosa_int) tempindex, datasize
      integer(kind=cosa_int) receiverequestnum,sendrequestnum
      logical myblock1, myblock2
      integer(kind=cosa_int) mrequest_perproc,maxproc_cutman
      integer(kind=cosa_int) sendid,recvid,numrecvid,numsendid,maxproc
      integer(kind=cosa_int) sendrequests(maxproc_cutman), &
                receiverequests(maxproc_cutman)
      integer(kind=cosa_int) receiveindices(mrequest_perproc,8,maxproc_cutman)
      real(kind=cosa_real) &
!del    sendarray   (2*nharms+1,mrequest_perproc,maxproc_cutman), &
!del    receivearray(2*nharms+1,mrequest_perproc,maxproc_cutman) &
        sendarray   (dim5h,mrequest_perproc,maxproc_cutman), &
        receivearray(dim5h,mrequest_perproc,maxproc_cutman)
      integer(kind=cosa_int) numsend(maxproc_cutman),numrecv(maxproc_cutman)
      integer(kind=cosa_int) myrecvid_num(0:maxproc-1),mysendid_num(0:maxproc-1)
!--------------------------------------------------------------------------

!     Store imax, jmax in ijmax for looping
      ijkmax1(1) = imax1
      ijkmax1(2) = jmax1
      ijkmax1(3) = kmax1
      ijkmax2(1) = imax2
      ijkmax2(2) = jmax2
      ijkmax2(3) = kmax2

!     Store boundary condition data in mnemonic names

      iblk1    = cutopo( 1)
      idir1    = cutopo( 2)
      inout1   = cutopo( 3)
      istr1(1) = cutopo( 4)
      iend1(1) = cutopo( 5)
      istr1(2) = cutopo( 6)
      iend1(2) = cutopo( 7)
      istr1(3) = cutopo( 8)
      iend1(3) = cutopo( 9)

      iblk2    = cutopo(10)
      idir2    = cutopo(11)
      inout2   = cutopo(12)
      istr2(1) = cutopo(13)
      iend2(1) = cutopo(14)
      istr2(2) = cutopo(15)
      iend2(2) = cutopo(16)
      istr2(3) = cutopo(17)
      iend2(3) = cutopo(18)

      iord(1)  = cutopo(19)
      iord(2)  = cutopo(20)
      iord(3)  = cutopo(21)

      call myblock(iblk1,myblock1,.false.)
      call myblock(iblk2,myblock2,.false.)

!     Set needed variables depending on whether the boundary is
!     the inner boundary (inout1 = 1) or the outer boundary (inout1 > 1)
!         ibcpt  = boundary point of block 1
!         ibcinc = increment to second boundary point of block 1
!         ibcpt2 = ibcpt + ibcinc
!         inr    = interior point of block 2.
!         inrinc = increment to second interior point of block 2
!                  inr2 = inr + inrinc

      if (inout1 .eq. 1) then
        ibcpt  = 0
      else
        ibcpt  = ijkmax1(idir1) + 1
      end if
!     
      if (inout2 .eq. 1) then
        inr    = 2
      else
        inr    = ijkmax2(idir2) - 1
      end if

      do l = 1, 3
        if (l .ne. idir1) then
          if (iend1(l) .gt. istr1(l)) then
            iend1(l)  = iend1(l) + 1
          else
            istr1(l)  = istr1(l) + 1
          end if
        end if
      end do

      do l = 1, 3
        if (l .ne. idir2) then
          if (iend2(l) .gt. istr2(l)) then
            iend2(l)  = iend2(l) + 1
          else
            istr2(l)  = istr2(l) + 1
          end if
        end if
      end do

!     Find the length of the two outer loops and loop over these using
!     offsets and delta function to set the two cut data points to the
!     two interior data points of block 2.

      do l = 1, 3

        len(l) = abs ( iend1(l) - istr1(l) )

!       Increment/Decrement 

        if ( iend1(l) .gt. istr1(l) ) then
          isgn1(l) =   1
        else
          isgn1(l) = - 1
        end if

!       Increment/Decrement 

        if ( iend2(l) .gt. istr2(l) ) then
          isgn2(l) =   1
        else
          isgn2(l) = - 1
        end if

      end do

!     ii1 first  boundary point of block 1
!     ii2 second boundary point of block 1
!     in1 first  interior point of block 2
!     in2 second interior point of block 2

      ic1 = cyc (idir1, 1)
      ic2 = cyc (idir1, 2)
      ic3 = cyc (idir1, 3)

      jc1 = iord (ic1)
      jc2 = iord (ic2)
      jc3 = iord (ic3)

      do i3 = 0, len(ic3)
        do i2 = 0, len(ic2)

          ii1 = ibcpt                        * krd (ic1, 1) + &
                (istr1(ic2) + isgn1(ic2)*i2) * krd (ic2, 1) + &
                (istr1(ic3) + isgn1(ic3)*i3) * krd (ic3, 1)
          jj1 = ibcpt                        * krd (ic1, 2) + &
                (istr1(ic2) + isgn1(ic2)*i2) * krd (ic2, 2) + &
                (istr1(ic3) + isgn1(ic3)*i3) * krd (ic3, 2)
          kk1 = ibcpt                        * krd (ic1, 3) + &
                (istr1(ic2) + isgn1(ic2)*i2) * krd (ic2, 3) + &
                (istr1(ic3) + isgn1(ic3)*i3) * krd (ic3, 3)

          in1 = inr                          * krd (jc1, 1) + &
                (istr2(jc2) + isgn2(jc2)*i2) * krd (jc2, 1) + &
                (istr2(jc3) + isgn2(jc3)*i3) * krd (jc3, 1)
          jn1 = inr                          * krd (jc1, 2) + &
                (istr2(jc2) + isgn2(jc2)*i2) * krd (jc2, 2) + &
                (istr2(jc3) + isgn2(jc3)*i3) * krd (jc3, 2)
          kn1 = inr                          * krd (jc1, 3) + &
                (istr2(jc2) + isgn2(jc2)*i2) * krd (jc2, 3) + &
                (istr2(jc3) + isgn2(jc3)*i3) * krd (jc3, 3)

          if(myblock1 .and. myblock2) then

            do n = 0, 2*nharms*hbmove
              x1(ii1,jj1,kk1,n) = x2(in1,jn1,kn1,n)
              y1(ii1,jj1,kk1,n) = y2(in1,jn1,kn1,n)
              z1(ii1,jj1,kk1,n) = z2(in1,jn1,kn1,n)
            end do

          else if(myblock1) then

            call getownerid(iblk2,recvid)
            numrecvid = myrecvid_num(recvid)
            numrecv(numrecvid) = numrecv(numrecvid)+1 
            receiveindices(numrecv(numrecvid),1,numrecvid) = ii1
            receiveindices(numrecv(numrecvid),2,numrecvid) = jj1
            receiveindices(numrecv(numrecvid),3,numrecvid) = kk1
            receiveindices(numrecv(numrecvid),4,numrecvid) = ixyz1
            receiveindices(numrecv(numrecvid),5,numrecvid) = 0
            receiveindices(numrecv(numrecvid),6,numrecvid) = imax1
            receiveindices(numrecv(numrecvid),7,numrecvid) = jmax1
            receiveindices(numrecv(numrecvid),8,numrecvid) = kmax1
            numrecv(numrecvid) = numrecv(numrecvid)+1
            receiveindices(numrecv(numrecvid),1,numrecvid) = ii1
            receiveindices(numrecv(numrecvid),2,numrecvid) = jj1
            receiveindices(numrecv(numrecvid),3,numrecvid) = kk1
            receiveindices(numrecv(numrecvid),4,numrecvid) = ixyz1
            receiveindices(numrecv(numrecvid),5,numrecvid) = 1
            receiveindices(numrecv(numrecvid),6,numrecvid) = imax1
            receiveindices(numrecv(numrecvid),7,numrecvid) = jmax1
            receiveindices(numrecv(numrecvid),8,numrecvid) = kmax1
            numrecv(numrecvid) = numrecv(numrecvid)+1
            receiveindices(numrecv(numrecvid),1,numrecvid) = ii1
            receiveindices(numrecv(numrecvid),2,numrecvid) = jj1
            receiveindices(numrecv(numrecvid),3,numrecvid) = kk1
            receiveindices(numrecv(numrecvid),4,numrecvid) = ixyz1
            receiveindices(numrecv(numrecvid),5,numrecvid) = 2
            receiveindices(numrecv(numrecvid),6,numrecvid) = imax1
            receiveindices(numrecv(numrecvid),7,numrecvid) = jmax1
            receiveindices(numrecv(numrecvid),8,numrecvid) = kmax1
            if (numrecv(numrecvid).gt. mrequest_perproc) then
              write(*,*) 'numrecv (cut_x)--> INCREASE MREQUEST_PERPROC'
              call abortmpi()
            end if
            if (numrecvid.gt. maxproc_cutman) then
              write(*,*) 'numrecv (cut_x)--> INCREASE MAXPROC_CUTMAN'
              call abortmpi()
            end if
          else if (myblock2) then
            call getownerid(iblk1,sendid)
            numsendid = mysendid_num(sendid)
            numsend(numsendid) = numsend(numsendid)+1 
            tempindex = 1
            do n = 0, 2*nharms*hbmove
              sendarray(tempindex,numsend(numsendid),numsendid) = &
                x2(in1,jn1,kn1,n)
              tempindex = tempindex + 1
            end do
            numsend(numsendid) = numsend(numsendid)+1
            tempindex  = 1
            do n = 0, 2*nharms*hbmove
              sendarray(tempindex,numsend(numsendid),numsendid) = &
                y2(in1,jn1,kn1,n)
              tempindex = tempindex + 1
            end do
            numsend(numsendid) = numsend(numsendid)+1
            tempindex  = 1
            do n = 0, 2*nharms*hbmove
              sendarray(tempindex,numsend(numsendid),numsendid) = &
                z2(in1,jn1,kn1,n)
              tempindex = tempindex + 1
            end do
            if (numsend(numsendid).gt. mrequest_perproc) then
              write(*,*) 'numsend (cut_x)--> INCREASE MREQUEST_PERPROC'
              call abortmpi()
            end if
            if (numsendid.gt. maxproc_cutman) then
              write(*,*) 'numsend (cut_x)--> INCREASE MAXPROC_CUTMAN'
              call abortmpi()
            end if
          end if

        end do
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine cutman_v(nl,cutopo,vol)
!-----------------------------------------------------------------------

      use parallelutils
      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,iblk1,icut,iblk2,imax1,jmax1,kmax1,imax2,jmax2,kmax2, &
                iv1,iv2
      integer(kind=cosa_int) cutopo(21,myncuts),a,b,c,i,reqnum
      logical myblock1, myblock2
      real(kind=cosa_real) vol(*)
      integer(kind=cosa_int) sendrequestnum,receiverequestnum
      integer(kind=cosa_int) i_id,i_in_id
      real(kind=cosa_real), allocatable :: lsendarray(:,:),lreceivearray(:,:)
      integer(kind=cosa_int) sendid,recvid

      numrecv = 0  ! array syntax
      numsend = 0  ! array syntax

      allocate(lsendarray(mrequest_perproc, maxproc_cutman))
      allocate(lreceivearray(mrequest_perproc, maxproc_cutman))

      do icut = 1,myncuts
        iblk1 = cutopo( 1,icut)
        iblk2 = cutopo(10,icut)
        call myblock(iblk1,myblock1,.true.)
        call myblock(iblk2,myblock2,.true.)
        if (myblock1) then
          imax1  = i_imax     (iblk1,nl)
          jmax1  = j_jmax     (iblk1,nl)
          kmax1  = k_kmax     (iblk1,nl)
          iv1    = 1 + off_p1 (iblk1,nl)
        else
          iv1    = 1
        end if
        if (myblock2) then
          imax2  = i_imax     (iblk2,nl)
          jmax2  = j_jmax     (iblk2,nl)
          kmax2  = k_kmax     (iblk2,nl)
          iv2    = 1 + off_p1 (iblk2,nl)
        else
          iv2    = 1
        end if
        call cut_v(imax1,jmax1,kmax1,vol(iv1),imax2,jmax2,kmax2, &
             vol(iv2),cutopo(1,icut),lreceivearray,receiverequests, &
             numrecv,receiveindices,lsendarray,sendrequests,numsend,iv1, &
             mrequest_perproc,maxproc_cutman,maxproc,myrecvid_num, &
             mysendid_num)
      end do

      do i_id=1,maxrecvid_num
        recvid = myrecvid(i_id)
        call receiveblockdata_fromid(lreceivearray(1,i_id),recvid, &
             numrecv(i_id),receiverequests(i_id))
      end do
      do i_id=1,maxsendid_num
        sendid = mysendid(i_id)
        call sendblockdata_toid(lsendarray(1,i_id),sendid, &
             numsend(i_id),sendrequests(i_id))
      end do
      do i_id=1,maxrecvid_num
        call waitanymessages(receiverequests,maxrecvid_num,reqnum)
        do i_in_id=1,numrecv(reqnum)
          a     = receiveindices(i_in_id,1,reqnum)
          b     = receiveindices(i_in_id,2,reqnum)
          c     = receiveindices(i_in_id,3,reqnum)
          iv1   = receiveindices(i_in_id,4,reqnum)
          imax1 = receiveindices(i_in_id,5,reqnum)
          jmax1 = receiveindices(i_in_id,6,reqnum)
          kmax1 = receiveindices(i_in_id,7,reqnum)
          call copyvreceivedata(imax1,jmax1,kmax1,a,b,c,vol(iv1), &
               lreceivearray,i_in_id,reqnum)
        end do
      enddo

      call waitallmessages(sendrequests,maxsendid_num)

      deallocate(lsendarray, lreceivearray)

      return

      end
!-----------------------------------------------------------------------
      subroutine cut_v(imax1,jmax1,kmax1,vol1,imax2,jmax2,kmax2,vol2, &
           cutopo,receivearray,receiverequests,numrecv, &
           receiveindices,sendarray,sendrequests,numsend,iv1, &
           mrequest_perproc,maxproc_cutman,maxproc,myrecvid_num, &
           mysendid_num)
!-----------------------------------------------------------------------

!     Routine to do cut boundary condition.
!     Flow data stored in Q1 are updated from data in the interior of Q2.
!      * sendrequests/receivrequests have indexes which is the process I am communicating to
!      * sendarray/recvarray have two components:
!        - index of point to send/recv
!        - index of process I have to send/recv to/from
!      * receiveindices have three components:
!        - index of point to send/recv
!        - index of information to trasmit (1-6)
!        - index of process to send/recv to/from

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax1,jmax1,kmax1,imax2,jmax2,kmax2,iblk1,iblk2
      integer(kind=cosa_int) cutopo(21),mrequest
      integer(kind=cosa_int) ijkmax1(3),ijkmax2(3),istr1(3),iend1(3),istr2(3), &
           iend2(3),isgn1(3),isgn2(3),iord(3),len(3),idir1,idir2,inout1, &
           inout2,ibcpt,ibcpt2,inr,inr2,l,ic1,ic2,ic3,jc1,jc2,jc3,i2,i3, &
           ii1,jj1,kk1,in1,jn1,kn1
      real(kind=cosa_real) &
        vol1(0:imax1,0:jmax1,0:kmax1), &
        vol2(0:imax2,0:jmax2,0:kmax2)
      integer(kind=cosa_int) tempindex, datasize
      integer(kind=cosa_int) receiverequestnum,iv1,sendrequestnum
      logical myblock1, myblock2
      integer(kind=cosa_int) mrequest_perproc,maxproc_cutman
      integer(kind=cosa_int) sendid,recvid,numrecvid,numsendid,maxproc
      integer(kind=cosa_int) sendrequests(maxproc_cutman), &
                receiverequests(maxproc_cutman)
      integer(kind=cosa_int) receiveindices(mrequest_perproc,8,maxproc_cutman)
      real(kind=cosa_real) &
        sendarray   (mrequest_perproc,maxproc_cutman), &
        receivearray(mrequest_perproc,maxproc_cutman)
      integer(kind=cosa_int) numsend(maxproc_cutman),numrecv(maxproc_cutman)
      integer(kind=cosa_int) myrecvid_num(0:maxproc-1),mysendid_num(0:maxproc-1)
!--------------------------------------------------------------------------

!     Store imax, jmax in ijmax for looping
      ijkmax1(1) = imax1
      ijkmax1(2) = jmax1
      ijkmax1(3) = kmax1
      ijkmax2(1) = imax2
      ijkmax2(2) = jmax2
      ijkmax2(3) = kmax2

!     Store boundary condition data in mnemonic names

      iblk1    = cutopo( 1)
      idir1    = cutopo( 2)
      inout1   = cutopo( 3)
      istr1(1) = cutopo( 4)
      iend1(1) = cutopo( 5)
      istr1(2) = cutopo( 6)
      iend1(2) = cutopo( 7)
      istr1(3) = cutopo( 8)
      iend1(3) = cutopo( 9)

      iblk2    = cutopo(10)
      idir2    = cutopo(11)
      inout2   = cutopo(12)
      istr2(1) = cutopo(13)
      iend2(1) = cutopo(14)
      istr2(2) = cutopo(15)
      iend2(2) = cutopo(16)
      istr2(3) = cutopo(17)
      iend2(3) = cutopo(18)

      iord(1)  = cutopo(19)
      iord(2)  = cutopo(20)
      iord(3)  = cutopo(21)

      call myblock(iblk1,myblock1,.false.)
      call myblock(iblk2,myblock2,.false.)

!     Set needed variables depending on whether the boundary is
!     the inner boundary (inout1 = 1) or the outer boundary (inout1 > 1)
!         ibcpt  = boundary point of block 1
!         ibcinc = increment to second boundary point of block 1
!         ibcpt2 = ibcpt + ibcinc
!         inr    = interior point of block 2.
!         inrinc = increment to second interior point of block 2
!                  inr2 = inr + inrinc

      if (inout1 .eq. 1) then
        ibcpt  =   0
      else
        ibcpt  = ijkmax1(idir1)
      end if
!     
      if (inout2 .eq. 1) then
        inr    =   1
      else
        inr    =   ijkmax2(idir2) - 1
      end if

!     Find the length of the two outer loops and loop over these using
!     offsets and delta function to set the two cut data points to the
!     two interior data points of block 2.

      do l = 1, 3

         len(l) = abs ( iend1(l) - istr1(l) )

!     Increment/Decrement 

         if ( iend1(l) .gt. istr1(l) ) then
            isgn1(l) =   1
         else
            isgn1(l) = - 1
         end if

!     Increment/Decrement 

         if ( iend2(l) .gt. istr2(l) ) then
            isgn2(l) =   1
         else
            isgn2(l) = - 1
         end if

      end do

!     ii1 first  boundary point of block 1
!     ii2 second boundary point of block 1
!     in1 first  interior point of block 2
!     in2 second interior point of block 2

      ic1 = cyc (idir1, 1)
      ic2 = cyc (idir1, 2)
      ic3 = cyc (idir1, 3)

      jc1 = iord (ic1)
      jc2 = iord (ic2)
      jc3 = iord (ic3)

      do i3 = 0, len(ic3)
        do i2 = 0, len(ic2)

          ii1 = ibcpt                        * krd (ic1, 1) + &
                (istr1(ic2) + isgn1(ic2)*i2) * krd (ic2, 1) + &
                (istr1(ic3) + isgn1(ic3)*i3) * krd (ic3, 1)
          jj1 = ibcpt                        * krd (ic1, 2) + &
                (istr1(ic2) + isgn1(ic2)*i2) * krd (ic2, 2) + &
                (istr1(ic3) + isgn1(ic3)*i3) * krd (ic3, 2)
          kk1 = ibcpt                        * krd (ic1, 3) + &
                (istr1(ic2) + isgn1(ic2)*i2) * krd (ic2, 3) + &
                (istr1(ic3) + isgn1(ic3)*i3) * krd (ic3, 3)

          in1 = inr                          * krd (jc1, 1) + &
                (istr2(jc2) + isgn2(jc2)*i2) * krd (jc2, 1) + &
                (istr2(jc3) + isgn2(jc3)*i3) * krd (jc3, 1)
          jn1 = inr                          * krd (jc1, 2) + &
                (istr2(jc2) + isgn2(jc2)*i2) * krd (jc2, 2) + &
                (istr2(jc3) + isgn2(jc3)*i3) * krd (jc3, 2)
          kn1 = inr                          * krd (jc1, 3) + &
                (istr2(jc2) + isgn2(jc2)*i2) * krd (jc2, 3) + &
                (istr2(jc3) + isgn2(jc3)*i3) * krd (jc3, 3)

         if(myblock1 .and. myblock2) then

           vol1(ii1,jj1,kk1) = vol2(in1,jn1,kn1)

         else if(myblock1) then

           call getownerid(iblk2,recvid)
           numrecvid = myrecvid_num(recvid)
           numrecv(numrecvid) = numrecv(numrecvid)+1 
           receiveindices(numrecv(numrecvid),1,numrecvid) = ii1
           receiveindices(numrecv(numrecvid),2,numrecvid) = jj1
           receiveindices(numrecv(numrecvid),3,numrecvid) = kk1
           receiveindices(numrecv(numrecvid),4,numrecvid) = iv1
           receiveindices(numrecv(numrecvid),5,numrecvid) = imax1
           receiveindices(numrecv(numrecvid),6,numrecvid) = jmax1
           receiveindices(numrecv(numrecvid),7,numrecvid) = kmax1
           if (numrecv(numrecvid).gt. mrequest_perproc) then
             write(*,*) 'numrecv (cut_v)--> INCREASE MREQUEST_PERPROC'
             call abortmpi()
           end if

         else if(myblock2) then

           call getownerid(iblk1,sendid)
           numsendid = mysendid_num(sendid)
           numsend(numsendid) = numsend(numsendid)+1 
           sendarray(numsend(numsendid),numsendid) = vol2(in1,jn1,kn1)
           if (numsend(numsendid).gt. mrequest_perproc) then
             write(*,*) 'numsend (cut_v)--> INCREASE MREQUEST_PERPROC'
             call abortmpi()
           end if

         end if

        end do
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine imetrix(nl,x,y,z,xdot,ydot,zdot,si,sj,sk)
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,iblock,imax,jmax,kmax,ixyz,iimt
      real(kind=cosa_real) x(*),y(*),z(*),xdot(*),ydot(*),zdot(*),si(*),sj(*), &
          sk(*)

      do iblock = 1,mynblocks
        imax  = i_imax  (iblock,nl)
        jmax  = j_jmax  (iblock,nl)
        kmax  = k_kmax  (iblock,nl)
        ixyz  = 1 + off_p2(iblock,nl)         * dim5h
        iimt  = 1 + off_p2 (iblock,nl) * lmet * dim5h
        call imetrix_b(x(ixyz),y(ixyz),z(ixyz),xdot(ixyz),ydot(ixyz), &
                       zdot(ixyz),si(iimt),sj(iimt),sk(iimt),iblock, &
                       imax,jmax,kmax,lmet,nharms)
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine imetrix_b(x,y,z,xdot,ydot,zdot,si,sj,sk,iblock,imax, &
                           jmax,kmax,lmet,nharms)
!-----------------------------------------------------------------------

!     si(1,i,j,k) = x-component of unit normal of xi-face
!     si(2,i,j,k) = y-component of unit normal of xi-face
!     si(3,i,j,k) = z-component of unit normal of xi-face
!     si(4,i,j,k) = area of xi-face
!     si(5,i,j,k) = x-component of velocity of xi-face midpoint
!     si(6,i,j,k) = y-component of velocity of xi-face midpoint
!     si(7,i,j,k) = z-component of velocity of xi-face midpoint

!     sj(1,i,j,k) = x-component of unit normal of eta-face
!     sj(2,i,j,k) = y-component of unit normal of eta-face
!     sj(3,i,j,k) = z-component of unit normal of eta-face
!     sj(4,i,j,k) = area of eta-face
!     sj(5,i,j,k) = x-component of velocity of eta-face midpoint
!     sj(6,i,j,k) = y-component of velocity of eta-face midpoint
!     sj(7,i,j,k) = y-component of velocity of eta-face midpoint

!     sk(1,i,j,k) = x-component of unit normal of zeta-face
!     sk(2,i,j,k) = y-component of unit normal of zeta-face
!     sk(3,i,j,k) = z-component of unit normal of zeta-face
!     sk(4,i,j,k) = area of zeta-face
!     sk(5,i,j,k) = x-component of velocity of zeta-face midpoint
!     sk(6,i,j,k) = y-component of velocity of zeta-face midpoint
!     sk(7,i,j,k) = z-component of velocity of zeta-face midpoint

!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,lmet,nharms
      integer(kind=cosa_int) iblock,i,j,k,idir,n
      real(kind=cosa_real) &
          x   (         0:imax+1,0:jmax+1,0:kmax+1,0:2*nharms*hbmove), &
          y   (         0:imax+1,0:jmax+1,0:kmax+1,0:2*nharms*hbmove), &
          z   (         0:imax+1,0:jmax+1,0:kmax+1,0:2*nharms*hbmove), &
          xdot(         0:imax+1,0:jmax+1,0:kmax+1,0:2*nharms*hbmove), &
          ydot(         0:imax+1,0:jmax+1,0:kmax+1,0:2*nharms*hbmove), &
          zdot(         0:imax+1,0:jmax+1,0:kmax+1,0:2*nharms*hbmove), &
          si  (lmet    ,0:imax+1,0:jmax+1,0:kmax+1,0:2*nharms*hbmove), &
          sj  (lmet    ,0:imax+1,0:jmax+1,0:kmax+1,0:2*nharms*hbmove), &
          sk  (lmet    ,0:imax+1,0:jmax+1,0:kmax+1,0:2*nharms*hbmove)
      real(kind=cosa_real) r(3), s(3), m(3)

      do n=0,2*nharms*hbmove

!------ xi-face normals
        do k = 1,kmax-1
          do j = 1,jmax-1
            do i = 1,imax
              si(1,i,j,k,n) = ( (y(i,j+1,k+1,n) - y(i,j  ,k,n)) * &
                                (z(i,j  ,k+1,n) - z(i,j+1,k,n)) - &
                                (z(i,j+1,k+1,n) - z(i,j  ,k,n)) * &
                                (y(i,j  ,k+1,n) - y(i,j+1,k,n)) ) / 2
              si(2,i,j,k,n) = ( (z(i,j+1,k+1,n) - z(i,j  ,k,n)) * &
                                (x(i,j  ,k+1,n) - x(i,j+1,k,n)) - &
                                (x(i,j+1,k+1,n) - x(i,j  ,k,n)) * &
                                (z(i,j  ,k+1,n) - z(i,j+1,k,n)) ) / 2
              si(3,i,j,k,n) = ( (x(i,j+1,k+1,n) - x(i,j  ,k,n)) * &
                                (y(i,j  ,k+1,n) - y(i,j+1,k,n)) - &
                                (y(i,j+1,k+1,n) - y(i,j  ,k,n)) * &
                                (x(i,j  ,k+1,n) - x(i,j+1,k,n)) ) / 2
              si(4,i,j,k,n) = sqrt ( si(1,i,j,k,n) * si(1,i,j,k,n) + &
                                     si(2,i,j,k,n) * si(2,i,j,k,n) + &
                                     si(3,i,j,k,n) * si(3,i,j,k,n) )

              if (si(4,i,j,k,n).le.0) then
                write(*,*) 'xi-area is le 0 at n,i,j,k,iblock: ', &
                           n,i,j,k,iblock
                if (.not.dgn_cell) stop
              end if

            end do
          end do
        end do

!------ normalize xi-face normals
        do k = 1,kmax-1
          do j = 1,jmax-1
            do i = 1,imax
              do idir = 1,3
                if ((si(4,i,j,k,n).eq.0).and.dgn_cell) cycle
                si(idir,i,j,k,n) = si(idir,i,j,k,n) / si(4,i,j,k,n)
              end do
            end do
          end do
        end do

!------ eta-face normals
        do k = 1,kmax-1
          do i = 1,imax-1
            do j = 1,jmax
              sj(1,i,j,k,n) = ( (z(i+1,j,k+1,n) - z(i  ,j,k,n)) * &
                                (y(i  ,j,k+1,n) - y(i+1,j,k,n)) - &
                                (y(i+1,j,k+1,n) - y(i  ,j,k,n)) * &
                                (z(i  ,j,k+1,n) - z(i+1,j,k,n)) ) / 2
              sj(2,i,j,k,n) = ( (x(i+1,j,k+1,n) - x(i  ,j,k,n)) * &
                                (z(i  ,j,k+1,n) - z(i+1,j,k,n)) - &
                                (z(i+1,j,k+1,n) - z(i  ,j,k,n)) * &
                                (x(i  ,j,k+1,n) - x(i+1,j,k,n)) ) / 2
              sj(3,i,j,k,n) = ( (y(i+1,j,k+1,n) - y(i  ,j,k,n)) * &
                                (x(i  ,j,k+1,n) - x(i+1,j,k,n)) - &
                                (x(i+1,j,k+1,n) - x(i  ,j,k,n)) * &
                                (y(i  ,j,k+1,n) - y(i+1,j,k,n)) ) / 2
              sj(4,i,j,k,n) = sqrt ( sj(1,i,j,k,n) * sj(1,i,j,k,n) + &
                                     sj(2,i,j,k,n) * sj(2,i,j,k,n) + &
                                     sj(3,i,j,k,n) * sj(3,i,j,k,n) )

              if (sj(4,i,j,k,n).le.0) then
                write(*,*) 'eta-area is le 0 at n,i,j,k,iblock: ', &
                           n,i,j,k,iblock
                if (.not.dgn_cell) stop
              end if

            end do
          end do
        end do

!------ normalize eta-face normals
        do k = 1,kmax-1
          do i = 1,imax-1
            do j = 1,jmax
              do idir = 1,3
                if ((sj(4,i,j,k,n).eq.0).and.dgn_cell) cycle
                sj(idir,i,j,k,n) = sj(idir,i,j,k,n) / sj(4,i,j,k,n)
              end do
            end do
          end do
        end do

!------ zeta-face normals
        do j = 1,jmax-1
          do i = 1,imax-1
            do k = 1,kmax
              sk(1,i,j,k,n) = ( (y(i+1,j+1,k,n) - y(i  ,j,k,n)) * &
                                (z(i  ,j+1,k,n) - z(i+1,j,k,n)) - &
                                (z(i+1,j+1,k,n) - z(i  ,j,k,n)) * &
                                (y(i  ,j+1,k,n) - y(i+1,j,k,n)) ) / 2
              sk(2,i,j,k,n) = ( (z(i+1,j+1,k,n) - z(i  ,j,k,n)) * &
                                (x(i  ,j+1,k,n) - x(i+1,j,k,n)) - &
                                (x(i+1,j+1,k,n) - x(i  ,j,k,n)) * &
                                (z(i  ,j+1,k,n) - z(i+1,j,k,n)) ) / 2
              sk(3,i,j,k,n) = ( (x(i+1,j+1,k,n) - x(i  ,j,k,n)) * &
                                (y(i  ,j+1,k,n) - y(i+1,j,k,n)) - &
                                (y(i+1,j+1,k,n) - y(i  ,j,k,n)) * &
                                (x(i  ,j+1,k,n) - x(i+1,j,k,n)) ) / 2
              sk(4,i,j,k,n) = sqrt ( sk(1,i,j,k,n)*sk(1,i,j,k,n) + &
                                     sk(2,i,j,k,n)*sk(2,i,j,k,n) + &
                                     sk(3,i,j,k,n)*sk(3,i,j,k,n) )

              if (sk(4,i,j,k,n).le.0) then
                write(*,*) 'zeta-area is le 0 at n,i,j,k,iblock: ', &
                           n,i,j,k,iblock
                if (.not.dgn_cell) stop
              end if

            end do
          end do
        end do

!------ normalize zeta-face normals
        do j = 1,jmax-1
          do i = 1,imax-1
            do k = 1,kmax
              do idir = 1,3
                if ((sk(4,i,j,k,n).eq.0).and.dgn_cell) cycle
                sk(idir,i,j,k,n) = sk(idir,i,j,k,n) / sk(4,i,j,k,n)
              end do
            end do
          end do
        end do

        if (moving) then

!-------- xi-face velocities
          do k = 1,kmax-1
            do j = 1,jmax-1
              do i = 1,imax
                si(5,i,j,k,n)= ( xdot(i,j+1,k+1,n) + xdot(i,j+1,k,n) + &
                                 xdot(i,j  ,k+1,n) + xdot(i,j  ,k,n) )/4
                si(6,i,j,k,n)= ( ydot(i,j+1,k+1,n) + ydot(i,j+1,k,n) + &
                                 ydot(i,j  ,k+1,n) + ydot(i,j  ,k,n) )/4
                si(7,i,j,k,n)= ( zdot(i,j+1,k+1,n) + zdot(i,j+1,k,n) + &
                                 zdot(i,j  ,k+1,n) + zdot(i,j  ,k,n) )/4
              end do
            end do
          end do
!-------- eta-face velocities
          do k = 1,kmax-1
            do i = 1,imax-1
              do j = 1,jmax
                sj(5,i,j,k,n)= ( xdot(i+1,j,k+1,n) + xdot(i+1,j,k,n) + &
                                 xdot(i  ,j,k+1,n) + xdot(i  ,j,k,n) )/4
                sj(6,i,j,k,n)= ( ydot(i+1,j,k+1,n) + ydot(i+1,j,k,n) + &
                                 ydot(i  ,j,k+1,n) + ydot(i  ,j,k,n) )/4
                sj(7,i,j,k,n)= ( zdot(i+1,j,k+1,n) + zdot(i+1,j,k,n) + &
                                 zdot(i  ,j,k+1,n) + zdot(i  ,j,k,n) )/4
              end do
            end do
          end do
!-------- zeta-face velocities
          do j = 1,jmax-1
            do i = 1,imax-1
              do k = 1,kmax
                sk(5,i,j,k,n)= ( xdot(i+1,j+1,k,n) + xdot(i+1,j,k,n) + &
                                 xdot(i  ,j+1,k,n) + xdot(i  ,j,k,n) )/4
                sk(6,i,j,k,n)= ( ydot(i+1,j+1,k,n) + ydot(i+1,j,k,n) + &
                                 ydot(i  ,j+1,k,n) + ydot(i  ,j,k,n) )/4
                sk(7,i,j,k,n)= ( zdot(i+1,j+1,k,n) + zdot(i+1,j,k,n) + &
                                 zdot(i  ,j+1,k,n) + zdot(i  ,j,k,n) )/4
              end do
            end do
          end do

        end if

!------ normal component of face velocity computed as by 
!       Allen, IJNMF 2002 and Obayashi, AIAA J. 1992.

        if (hawt.and.rotating) then

!------   xi-face area moments
          do k = 1,kmax-1
            do j = 1,jmax-1
              do i = 1,imax
                r(1) = ( x(i,j,k,n) + x(i,j+1,k,n) + x(i,j+1,k+1,n) )/3- &
                       xrotc
                r(2) = ( y(i,j,k,n) + y(i,j+1,k,n) + y(i,j+1,k+1,n) )/3- &
                       yrotc
                r(3) = ( z(i,j,k,n) + z(i,j+1,k,n) + z(i,j+1,k+1,n) )/3

                s(1) = ( (y(i,j+1,k  ,n) - y(i,j  ,k  ,n)) * &
                         (z(i,j+1,k+1,n) - z(i,j  ,k  ,n)) - &
                         (z(i,j+1,k  ,n) - z(i,j  ,k  ,n)) * &
                         (y(i,j+1,k+1,n) - y(i,j  ,k  ,n)) ) / 2
                s(2) = ( (z(i,j+1,k  ,n) - z(i,j  ,k  ,n)) * &
                         (x(i,j+1,k+1,n) - x(i,j  ,k  ,n)) - &
                         (x(i,j+1,k  ,n) - x(i,j  ,k  ,n)) * &
                         (z(i,j+1,k+1,n) - z(i,j  ,k  ,n)) ) / 2
                s(3) = ( (x(i,j+1,k  ,n) - x(i,j  ,k  ,n)) * &
                         (y(i,j+1,k+1,n) - y(i,j  ,k  ,n)) - &
                         (y(i,j+1,k  ,n) - y(i,j  ,k  ,n)) * &
                         (x(i,j+1,k+1,n) - x(i,j  ,k  ,n)) ) / 2

                m(1) = r(2) * s(3) - r(3) * s(2)
                m(2) = r(3) * s(1) - r(1) * s(3)
                m(3) = r(1) * s(2) - r(2) * s(1)

                r(1) = ( x(i,j,k,n) + x(i,j+1,k+1,n) + x(i,j,k+1,n) )/3- &
                       xrotc
                r(2) = ( y(i,j,k,n) + y(i,j+1,k+1,n) + y(i,j,k+1,n) )/3- &
                       yrotc
                r(3) = ( z(i,j,k,n) + z(i,j+1,k+1,n) + z(i,j,k+1,n) )/3

                s(1) = ( (y(i,j+1,k+1,n) - y(i,j  ,k  ,n)) * &
                         (z(i,j  ,k+1,n) - z(i,j  ,k  ,n)) - &
                         (z(i,j+1,k+1,n) - z(i,j  ,k  ,n)) * &
                         (y(i,j  ,k+1,n) - y(i,j  ,k  ,n)) ) / 2
                s(2) = ( (z(i,j+1,k+1,n) - z(i,j  ,k  ,n)) * &
                         (x(i,j  ,k+1,n) - x(i,j  ,k  ,n)) - &
                         (x(i,j+1,k+1,n) - x(i,j  ,k  ,n)) * &
                         (z(i,j  ,k+1,n) - z(i,j  ,k  ,n)) ) / 2
                s(3) = ( (x(i,j+1,k+1,n) - x(i,j  ,k  ,n)) * &
                         (y(i,j  ,k+1,n) - y(i,j  ,k  ,n)) - &
                         (y(i,j+1,k+1,n) - y(i,j  ,k  ,n)) * &
                         (x(i,j  ,k+1,n) - x(i,j  ,k  ,n)) ) / 2
 
                m(1) = m(1) + r(2) * s(3) - r(3) * s(2)
                m(2) = m(2) + r(3) * s(1) - r(1) * s(3)
                m(3) = m(3) + r(1) * s(2) - r(2) * s(1)

!               At present rotational axis is z, hence only Mz is required.
                if ((si(4,i,j,k,n).eq.0).and.dgn_cell) cycle
                si(8,i,j,k,n) = ( m(3) * omegas ) / si(4,i,j,k,n)

              end do
            end do
          end do

!-------- eta-face area moments
          do k = 1,kmax-1
            do i = 1,imax-1
              do j = 1,jmax
                r(1) = ( x(i,j,k,n) + x(i,j,k+1,n) + x(i+1,j,k+1,n) )/3- &
                       xrotc
                r(2) = ( y(i,j,k,n) + y(i,j,k+1,n) + y(i+1,j,k+1,n) )/3- &
                       yrotc
                r(3) = ( z(i,j,k,n) + z(i,j,k+1,n) + z(i+1,j,k+1,n) )/3

                s(1) = ( (y(i  ,j,k+1,n) - y(i  ,j,k  ,n)) * &
                         (z(i+1,j,k+1,n) - z(i  ,j,k  ,n)) - &
                         (z(i  ,j,k+1,n) - z(i  ,j,k  ,n)) * &
                         (y(i+1,j,k+1,n) - y(i  ,j,k  ,n)) ) / 2
                s(2) = ( (z(i  ,j,k+1,n) - z(i  ,j,k  ,n)) * &
                         (x(i+1,j,k+1,n) - x(i  ,j,k  ,n)) - &
                         (x(i  ,j,k+1,n) - x(i  ,j,k  ,n)) * &
                         (z(i+1,j,k+1,n) - z(i  ,j,k  ,n)) ) / 2
                s(3) = ( (x(i  ,j,k+1,n) - x(i  ,j,k  ,n)) * &
                         (y(i+1,j,k+1,n) - y(i  ,j,k  ,n)) - &
                         (y(i  ,j,k+1,n) - y(i  ,j,k  ,n)) * &
                         (x(i+1,j,k+1,n) - x(i  ,j,k  ,n)) ) / 2

                m(1) = r(2) * s(3) - r(3) * s(2)
                m(2) = r(3) * s(1) - r(1) * s(3)
                m(3) = r(1) * s(2) - r(2) * s(1)

                r(1) = ( x(i,j,k,n) + x(i+1,j,k+1,n) + x(i+1,j,k,n) )/3- &
                       xrotc
                r(2) = ( y(i,j,k,n) + y(i+1,j,k+1,n) + y(i+1,j,k,n) )/3- &
                       yrotc
                r(3) = ( z(i,j,k,n) + z(i+1,j,k+1,n) + z(i+1,j,k,n) )/3

                s(1) = ( (y(i+1,j,k+1,n) - y(i  ,j,k  ,n)) * &
                         (z(i+1,j,k  ,n) - z(i  ,j,k  ,n)) - &
                         (z(i+1,j,k+1,n) - z(i  ,j,k  ,n)) * &
                         (y(i+1,j,k  ,n) - y(i  ,j,k  ,n)) ) / 2
                s(2) = ( (z(i+1,j,k+1,n) - z(i  ,j,k  ,n)) * &
                         (x(i+1,j,k  ,n) - x(i  ,j,k  ,n)) - &
                         (x(i+1,j,k+1,n) - x(i  ,j,k  ,n)) * &
                         (z(i+1,j,k  ,n) - z(i  ,j,k  ,n)) ) / 2
                s(3) = ( (x(i+1,j,k+1,n) - x(i  ,j,k  ,n)) * &
                         (y(i+1,j,k  ,n) - y(i  ,j,k  ,n)) - &
                         (y(i+1,j,k+1,n) - y(i  ,j,k  ,n)) * &
                         (x(i+1,j,k  ,n) - x(i  ,j,k  ,n)) ) / 2
 
                m(1) = m(1) + r(2) * s(3) - r(3) * s(2)
                m(2) = m(2) + r(3) * s(1) - r(1) * s(3)
                m(3) = m(3) + r(1) * s(2) - r(2) * s(1)

!               At present rotational axis is z, hence only Mz is required.
                if ((sj(4,i,j,k,n).eq.0).and.dgn_cell) cycle
                sj(8,i,j,k,n) = ( m(3) * omegas ) / sj(4,i,j,k,n)

              end do
            end do
          end do

!-------- zeta-face area moments
          do j = 1,jmax-1
            do i = 1,imax-1
              do k = 1,kmax
                r(1) = ( x(i,j,k,n) + x(i+1,j,k,n) + x(i+1,j+1,k,n) )/3- &
                       xrotc
                r(2) = ( y(i,j,k,n) + y(i+1,j,k,n) + y(i+1,j+1,k,n) )/3- &
                       yrotc
                r(3) = ( z(i,j,k,n) + z(i+1,j,k,n) + z(i+1,j+1,k,n) )/3

                s(1) = ( (y(i+1,j  ,k,n) - y(i  ,j  ,k,n)) * &
                         (z(i+1,j+1,k,n) - z(i  ,j  ,k,n)) - &
                         (z(i+1,j  ,k,n) - z(i  ,j  ,k,n)) * &
                         (y(i+1,j+1,k,n) - y(i  ,j  ,k,n)) ) / 2
                s(2) = ( (z(i+1,j  ,k,n) - z(i  ,j  ,k,n)) * &
                         (x(i+1,j+1,k,n) - x(i  ,j  ,k,n)) - &
                         (x(i+1,j  ,k,n) - x(i  ,j  ,k,n)) * &
                         (z(i+1,j+1,k,n) - z(i  ,j  ,k,n)) ) / 2
                s(3) = ( (x(i+1,j  ,k,n) - x(i  ,j  ,k,n)) * &
                         (y(i+1,j+1,k,n) - y(i  ,j  ,k,n)) - &
                         (y(i+1,j  ,k,n) - y(i  ,j  ,k,n)) * &
                         (x(i+1,j+1,k,n) - x(i  ,j  ,k,n)) ) / 2

                m(1) = r(2) * s(3) - r(3) * s(2)
                m(2) = r(3) * s(1) - r(1) * s(3)
                m(3) = r(1) * s(2) - r(2) * s(1)

                r(1) = ( x(i,j,k,n) + x(i+1,j+1,k,n) + x(i,j+1,k,n) )/3- &
                       xrotc
                r(2) = ( y(i,j,k,n) + y(i+1,j+1,k,n) + y(i,j+1,k,n) )/3- &
                       yrotc
                r(3) = ( z(i,j,k,n) + z(i+1,j+1,k,n) + z(i,j+1,k,n) )/3

                s(1) = ( (y(i+1,j+1,k,n) - y(i  ,j  ,k,n)) * &
                         (z(i  ,j+1,k,n) - z(i  ,j  ,k,n)) - &
                         (z(i+1,j+1,k,n) - z(i  ,j  ,k,n)) * &
                         (y(i  ,j+1,k,n) - y(i  ,j  ,k,n)) ) / 2
                s(2) = ( (z(i+1,j+1,k,n) - z(i  ,j  ,k,n)) * &
                         (x(i  ,j+1,k,n) - x(i  ,j  ,k,n)) - &
                         (x(i+1,j+1,k,n) - x(i  ,j  ,k,n)) * &
                         (z(i  ,j+1,k,n) - z(i  ,j  ,k,n)) ) / 2
                s(3) = ( (x(i+1,j+1,k,n) - x(i  ,j  ,k,n)) * &
                         (y(i  ,j+1,k,n) - y(i  ,j  ,k,n)) - &
                         (y(i+1,j+1,k,n) - y(i  ,j  ,k,n)) * &
                         (x(i  ,j+1,k,n) - x(i  ,j  ,k,n)) ) / 2
 
                m(1) = m(1) + r(2) * s(3) - r(3) * s(2)
                m(2) = m(2) + r(3) * s(1) - r(1) * s(3)
                m(3) = m(3) + r(1) * s(2) - r(2) * s(1)

!               At present rotational axis is z, hence only Mz is required.
                if ((sk(4,i,j,k,n).eq.0).and.dgn_cell) cycle
                sk(8,i,j,k,n) = ( m(3) * omegas ) / sk(4,i,j,k,n)

              end do
            end do
          end do

        end if

      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine vmetrix_old(nl,x,y,z,si,sj,sk,xideri,xiderj,xiderk,etaderi, &
                         etaderj,etaderk,zetaderi,zetaderj,zetaderk,vol)
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,iblock,imax,jmax,kmax,ivol,ixyz,iimt,ivmt
      real (kind=cosa_real) x(*),y(*),z(*),si(*),sj(*),sk(*),xideri(*), &
           xiderj(*),xiderk(*),etaderi(*),etaderj(*),etaderk(*), &
           zetaderi(*),zetaderj(*),zetaderk(*),vol(*)

      do iblock = 1,mynblocks
        imax  = i_imax  (iblock,nl)
        jmax  = j_jmax  (iblock,nl)
        kmax  = k_kmax  (iblock,nl)
        ixyz  = 1 + off_p2(iblock,nl)         * dim5h
        ivol  = 1 + off_p1(iblock,nl)
        iimt  = 1 + off_p2 (iblock,nl) * lmet * dim5h
        ivmt  = 1 + off_0  (iblock,nl) * 3    * dim5h
        call vmetrix_old_b(x(ixyz),y(ixyz),z(ixyz),si(iimt),sj(iimt), &
                       sk(iimt),xideri(ivmt),xiderj(ivmt),xiderk(ivmt), &
                       etaderi(ivmt),etaderj(ivmt),etaderk(ivmt), &
                       zetaderi(ivmt),zetaderj(ivmt),zetaderk(ivmt), &
                       vol(ivol),iblock,imax,jmax,kmax,lmet,nharms)
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine vmetrix_old_b(x,y,z,si,sj,sk,xideri,xiderj,xiderk,etaderi, &
                           etaderj,etaderk,zetaderi,zetaderj,zetaderk, &
                           vol,iblock,imax,jmax,kmax,lmet,nharms)
!-----------------------------------------------------------------------
!     Oriented xi-surface   has components: vol*[  xi_x   xi_y   xi_z]^T
!     Oriented eta-surface  has components: vol*[ eta_x  eta_y  eta_z]^T
!     Oriented zeta-surface has components: vol*[zeta_x zeta_y zeta_z]^T
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,lmet,nharms
      integer(kind=cosa_int) iblock,i,j,k,idir,n
      real(kind=cosa_real) &
          x       (     0:imax+1,0:jmax+1,0:kmax+1,0:2*nharms*hbmove), &
          y       (     0:imax+1,0:jmax+1,0:kmax+1,0:2*nharms*hbmove), &
          z       (     0:imax+1,0:jmax+1,0:kmax+1,0:2*nharms*hbmove), &
          si      (lmet,0:imax+1,0:jmax+1,0:kmax+1,0:2*nharms*hbmove), &
          sj      (lmet,0:imax+1,0:jmax+1,0:kmax+1,0:2*nharms*hbmove), &
          sk      (lmet,0:imax+1,0:jmax+1,0:kmax+1,0:2*nharms*hbmove), &
          xideri  (3   ,imax    ,jmax    ,kmax    ,0:2*nharms*hbmove), &
          xiderj  (3   ,imax    ,jmax    ,kmax    ,0:2*nharms*hbmove), &
          xiderk  (3   ,imax    ,jmax    ,kmax    ,0:2*nharms*hbmove), &
          etaderi (3   ,imax    ,jmax    ,kmax    ,0:2*nharms*hbmove), &
          etaderj (3   ,imax    ,jmax    ,kmax    ,0:2*nharms*hbmove), &
          etaderk (3   ,imax    ,jmax    ,kmax    ,0:2*nharms*hbmove), &
          zetaderi(3   ,imax    ,jmax    ,kmax    ,0:2*nharms*hbmove), &
          zetaderj(3   ,imax    ,jmax    ,kmax    ,0:2*nharms*hbmove), &
          zetaderk(3   ,imax    ,jmax    ,kmax    ,0:2*nharms*hbmove), &
          vol     (0:imax,0:jmax,0:kmax)
      real(kind=cosa_real) xix, etax, zetax, xiy, etay, zetay, xiz, etaz, zetaz, &
                   volvisc

      do n=0,2*nharms*hbmove

!------ calculate face normals of auxiliary cells: xi-faces

        do k = 1,kmax-1
          do j = 0,jmax,jmax
            do i = 1,imax
              si(1,i,j,k,n) = ( (y(i,j+1,k+1,n) - y(i,j  ,k,n)) * &
                                (z(i,j  ,k+1,n) - z(i,j+1,k,n)) - &
                                (z(i,j+1,k+1,n) - z(i,j  ,k,n)) * &
                                (y(i,j  ,k+1,n) - y(i,j+1,k,n)) ) / 2
              si(2,i,j,k,n) = ( (z(i,j+1,k+1,n) - z(i,j  ,k,n)) * &
                                (x(i,j  ,k+1,n) - x(i,j+1,k,n)) - &
                                (x(i,j+1,k+1,n) - x(i,j  ,k,n)) * &
                                (z(i,j  ,k+1,n) - z(i,j+1,k,n)) ) / 2
              si(3,i,j,k,n) = ( (x(i,j+1,k+1,n) - x(i,j  ,k,n)) * &
                                (y(i,j  ,k+1,n) - y(i,j+1,k,n)) - &
                                (y(i,j+1,k+1,n) - y(i,j  ,k,n)) * &
                                (x(i,j  ,k+1,n) - x(i,j+1,k,n)) ) / 2
              si(4,i,j,k,n) = sqrt ( si(1,i,j,k,n) * si(1,i,j,k,n) + &
                                     si(2,i,j,k,n) * si(2,i,j,k,n) + &
                                     si(3,i,j,k,n) * si(3,i,j,k,n) )
              if (si(4,i,j,k,n).le.0) then
!-------------- MSC: areas are zero ONLY for auxiliary cells adjacent
!                    to grid boundaries which are NOT cuts. For auxiliary
!                    cells ajacent to cuts areas CANNOT be zero due to
!                    due to cuman_x
                write(*,*) &
                  'vmetrix warning: zero xi-area at n,i,j,k,iblock: ', &
                           n,i,j,k,iblock
                if (j.eq.0) then
                  do idir = 1,4
                    si(idir,i,j,k,n) = si(idir,i,j+1,k,n)
                  end do
                else if (j.eq.jmax) then
                  do idir = 1,4
                    si(idir,i,j,k,n) = si(idir,i,j-1,k,n)
                  end do
                end if
              end if

              do idir = 1,3
                si(idir,i,j,k,n) = si(idir,i,j,k,n) / si(4,i,j,k,n)
              end do

            end do
          end do
        end do

        do k = 0,kmax,kmax
          do j = 1,jmax-1
            do i = 1,imax
              si(1,i,j,k,n) = ( (y(i,j+1,k+1,n) - y(i,j  ,k,n)) * &
                                (z(i,j  ,k+1,n) - z(i,j+1,k,n)) - &
                                (z(i,j+1,k+1,n) - z(i,j  ,k,n)) * &
                                (y(i,j  ,k+1,n) - y(i,j+1,k,n)) ) / 2
              si(2,i,j,k,n) = ( (z(i,j+1,k+1,n) - z(i,j  ,k,n)) * &
                                (x(i,j  ,k+1,n) - x(i,j+1,k,n)) - &
                                (x(i,j+1,k+1,n) - x(i,j  ,k,n)) * &
                                (z(i,j  ,k+1,n) - z(i,j+1,k,n)) ) / 2
              si(3,i,j,k,n) = ( (x(i,j+1,k+1,n) - x(i,j  ,k,n)) * &
                                (y(i,j  ,k+1,n) - y(i,j+1,k,n)) - &
                                (y(i,j+1,k+1,n) - y(i,j  ,k,n)) * &
                                (x(i,j  ,k+1,n) - x(i,j+1,k,n)) ) / 2
              si(4,i,j,k,n) = sqrt ( si(1,i,j,k,n) * si(1,i,j,k,n) + &
                                     si(2,i,j,k,n) * si(2,i,j,k,n) + &
                                     si(3,i,j,k,n) * si(3,i,j,k,n) )
              if (si(4,i,j,k,n).le.0) then
!-------------- MSC: areas are zero ONLY for auxiliary cells adjacent
!                    to grid boundaries which are NOT cuts. For auxiliary
!                    cells ajacent to cuts areas CANNOT be zero due to
!                    due to cuman_x
                write(*,*) &
                  'vmetrix warning: zero xi-area at n,i,j,k,iblock: ', &
                           n,i,j,k,iblock
                if (k.eq.0) then
                  do idir = 1,3
                    si(idir,i,j,k,n) = si(idir,i,j,k+1,n)
                  end do
                else if (k.eq.kmax) then
                  do idir = 1,3
                    si(idir,i,j,k,n) = si(idir,i,j,k-1,n)
                  end do
                end if
              end if

              do idir = 1,3
                si(idir,i,j,k,n) = si(idir,i,j,k,n) / si(4,i,j,k,n)
              end do

            end do
          end do
        end do

!------ calculate face normals of auxiliary cells: eta-faces

        do k = 0,kmax,kmax
          do i = 1,imax-1
            do j = 1,jmax
              sj(1,i,j,k,n) = ( (z(i+1,j,k+1,n) - z(i  ,j,k,n)) * &
                                (y(i  ,j,k+1,n) - y(i+1,j,k,n)) - &
                                (y(i+1,j,k+1,n) - y(i  ,j,k,n)) * &
                                (z(i  ,j,k+1,n) - z(i+1,j,k,n)) ) / 2
              sj(2,i,j,k,n) = ( (x(i+1,j,k+1,n) - x(i  ,j,k,n)) * &
                                (z(i  ,j,k+1,n) - z(i+1,j,k,n)) - &
                                (z(i+1,j,k+1,n) - z(i  ,j,k,n)) * &
                                (x(i  ,j,k+1,n) - x(i+1,j,k,n)) ) / 2
              sj(3,i,j,k,n) = ( (y(i+1,j,k+1,n) - y(i  ,j,k,n)) * &
                                (x(i  ,j,k+1,n) - x(i+1,j,k,n)) - &
                                (x(i+1,j,k+1,n) - x(i  ,j,k,n)) * &
                                (y(i  ,j,k+1,n) - y(i+1,j,k,n)) ) / 2
              sj(4,i,j,k,n) = sqrt ( sj(1,i,j,k,n) * sj(1,i,j,k,n) + &
                                     sj(2,i,j,k,n) * sj(2,i,j,k,n) + &
                                     sj(3,i,j,k,n) * sj(3,i,j,k,n) )
              if (sj(4,i,j,k,n).le.0) then
!-------------- MSC: areas are zero ONLY for auxiliary cells adjacent
!                    to grid boundaries which are NOT cuts. For auxiliary
!                    cells ajacent to cuts areas CANNOT be zero due to
!                    due to cuman_x
                write(*,*) &
                  'vmetrix warning: zero eta-area at n,i,j,k,iblock: ', &
                           n,i,j,k,iblock
                if (k.eq.0) then
                  do idir = 1,4
                    sj(idir,i,j,k,n) = sj(idir,i,j,k+1,n)
                  end do
                else if (k.eq.kmax) then
                  do idir = 1,4
                    sj(idir,i,j,k,n) = sj(idir,i,j,k-1,n)
                  end do
                end if
              end if

              do idir = 1,3
                sj(idir,i,j,k,n) = sj(idir,i,j,k,n) / sj(4,i,j,k,n)
              end do

            end do
          end do
        end do

        do k = 1,kmax-1
          do i = 0,imax,imax
            do j = 1,jmax
              sj(1,i,j,k,n) = ( (z(i+1,j,k+1,n) - z(i  ,j,k,n)) * &
                                (y(i  ,j,k+1,n) - y(i+1,j,k,n)) - &
                                (y(i+1,j,k+1,n) - y(i  ,j,k,n)) * &
                                (z(i  ,j,k+1,n) - z(i+1,j,k,n)) ) / 2
              sj(2,i,j,k,n) = ( (x(i+1,j,k+1,n) - x(i  ,j,k,n)) * &
                                (z(i  ,j,k+1,n) - z(i+1,j,k,n)) - &
                                (z(i+1,j,k+1,n) - z(i  ,j,k,n)) * &
                                (x(i  ,j,k+1,n) - x(i+1,j,k,n)) ) / 2
              sj(3,i,j,k,n) = ( (y(i+1,j,k+1,n) - y(i  ,j,k,n)) * &
                                (x(i  ,j,k+1,n) - x(i+1,j,k,n)) - &
                                (x(i+1,j,k+1,n) - x(i  ,j,k,n)) * &
                                (y(i  ,j,k+1,n) - y(i+1,j,k,n)) ) / 2
              sj(4,i,j,k,n) = sqrt ( sj(1,i,j,k,n) * sj(1,i,j,k,n) + &
                                     sj(2,i,j,k,n) * sj(2,i,j,k,n) + &
                                     sj(3,i,j,k,n) * sj(3,i,j,k,n) )
              if (sj(4,i,j,k,n).le.0) then
!-------------- MSC: areas are zero ONLY for auxiliary cells adjacent
!                    to grid boundaries which are NOT cuts. For auxiliary
!                    cells ajacent to cuts areas CANNOT be zero due to
!                    due to cuman_x
                write(*,*) &
                  'vmetrix warning: zero eta-area at n,i,j,k,iblock: ', &
                           n,i,j,k,iblock
                if (i.eq.0) then
                  do idir = 1,4
                    sj(idir,i,j,k,n) = sj(idir,i+1,j,k,n)
                  end do
                else if (i.eq.imax) then
                  do idir = 1,4
                    sj(idir,i,j,k,n) = sj(idir,i-1,j,k,n)
                  end do
                end if
              end if

              do idir = 1,3
                sj(idir,i,j,k,n) = sj(idir,i,j,k,n) / sj(4,i,j,k,n)
              end do

            end do
          end do
        end do

!------ calculate face normals of auxiliary cells: zeta-faces

        do j = 0,jmax,jmax
          do i = 1,imax-1
            do k = 1,kmax
              sk(1,i,j,k,n) = ( (y(i+1,j+1,k,n) - y(i  ,j,k,n)) * &
                                (z(i  ,j+1,k,n) - z(i+1,j,k,n)) - &
                                (z(i+1,j+1,k,n) - z(i  ,j,k,n)) * &
                                (y(i  ,j+1,k,n) - y(i+1,j,k,n)) ) / 2
              sk(2,i,j,k,n) = ( (z(i+1,j+1,k,n) - z(i  ,j,k,n)) * &
                                (x(i  ,j+1,k,n) - x(i+1,j,k,n)) - &
                                (x(i+1,j+1,k,n) - x(i  ,j,k,n)) * &
                                (z(i  ,j+1,k,n) - z(i+1,j,k,n)) ) / 2
              sk(3,i,j,k,n) = ( (x(i+1,j+1,k,n) - x(i  ,j,k,n)) * &
                                (y(i  ,j+1,k,n) - y(i+1,j,k,n)) - &
                                (y(i+1,j+1,k,n) - y(i  ,j,k,n)) * &
                                (x(i  ,j+1,k,n) - x(i+1,j,k,n)) ) / 2
              sk(4,i,j,k,n) = sqrt ( sk(1,i,j,k,n)*sk(1,i,j,k,n) + &
                                     sk(2,i,j,k,n)*sk(2,i,j,k,n) + &
                                     sk(3,i,j,k,n)*sk(3,i,j,k,n) )

              if (sk(4,i,j,k,n).le.0) then
!-------------- MSC: areas are zero ONLY for auxiliary cells adjacent
!                    to grid boundaries which are NOT cuts. For auxiliary
!                    cells ajacent to cuts areas CANNOT be zero due to
!                    due to cuman_x
                write(*,*) &
                  'vmetrix warning: zero zeta-area at n,i,j,k,iblock: ', &
                           n,i,j,k,iblock
                if (j.eq.0) then
                  do idir = 1,4
                    sk(idir,i,j,k,n) = sk(idir,i,j+1,k,n)
                  end do
                else if (i.eq.imax) then
                  do idir = 1,4
                    sk(idir,i,j,k,n) = sk(idir,i,j-1,k,n)
                  end do
                end if
              end if

              do idir = 1,3
                sk(idir,i,j,k,n) = sk(idir,i,j,k,n) / sk(4,i,j,k,n)
              end do

            end do
          end do
        end do

        do j = 1,jmax-1
          do i = 0,imax,imax
            do k = 1,kmax
              sk(1,i,j,k,n) = ( (y(i+1,j+1,k,n) - y(i  ,j,k,n)) * &
                                (z(i  ,j+1,k,n) - z(i+1,j,k,n)) - &
                                (z(i+1,j+1,k,n) - z(i  ,j,k,n)) * &
                                (y(i  ,j+1,k,n) - y(i+1,j,k,n)) ) / 2
              sk(2,i,j,k,n) = ( (z(i+1,j+1,k,n) - z(i  ,j,k,n)) * &
                                (x(i  ,j+1,k,n) - x(i+1,j,k,n)) - &
                                (x(i+1,j+1,k,n) - x(i  ,j,k,n)) * &
                                (z(i  ,j+1,k,n) - z(i+1,j,k,n)) ) / 2
              sk(3,i,j,k,n) = ( (x(i+1,j+1,k,n) - x(i  ,j,k,n)) * &
                                (y(i  ,j+1,k,n) - y(i+1,j,k,n)) - &
                                (y(i+1,j+1,k,n) - y(i  ,j,k,n)) * &
                                (x(i  ,j+1,k,n) - x(i+1,j,k,n)) ) / 2
              sk(4,i,j,k,n) = sqrt ( sk(1,i,j,k,n)*sk(1,i,j,k,n) + &
                                     sk(2,i,j,k,n)*sk(2,i,j,k,n) + &
                                     sk(3,i,j,k,n)*sk(3,i,j,k,n) )

              if (sk(4,i,j,k,n).le.0) then
!-------------- MSC: areas are zero ONLY for auxiliary cells adjacent
!                    to grid boundaries which are NOT cuts. For auxiliary
!                    cells ajacent to cuts areas CANNOT be zero due to
!                    due to cuman_x
                write(*,*) &
                  'vmetrix warning: zero zeta-area at n,i,j,k,iblock: ', &
                           n,i,j,k,iblock
                if (i.eq.0) then
                  do idir = 1,4
                    sk(idir,i,j,k,n) = sk(idir,i+1,j,k,n)
                  end do
                else if (i.eq.imax) then
                  do idir = 1,4
                    sk(idir,i,j,k,n) = sk(idir,i-1,j,k,n)
                  end do
                end if
              end if

              do idir = 1,3
                sk(idir,i,j,k,n) = sk(idir,i,j,k,n) / sk(4,i,j,k,n)
              end do

            end do
          end do
        end do

!------ xi-face metrics
        do k = 1,kmax-1
          do j = 1,jmax-1
            do i = 1,imax

              volvisc   = (vol(i-1,j,k)+vol(i,j,k)) / 2

              xix       =   si(1,i,j,k,n)*si(4,i,j,k,n) / volvisc
              xiy       =   si(2,i,j,k,n)*si(4,i,j,k,n) / volvisc
              xiz       =   si(3,i,j,k,n)*si(4,i,j,k,n) / volvisc

              etax      = ( sj(1,i-1,j  ,k,n)*sj(4,i-1,j  ,k,n) + &
                            sj(1,i  ,j  ,k,n)*sj(4,i  ,j  ,k,n) + &
                            sj(1,i-1,j+1,k,n)*sj(4,i-1,j+1,k,n) + &
                            sj(1,i  ,j+1,k,n)*sj(4,i  ,j+1,k,n) ) / 4 / &
                          volvisc
              etay      = ( sj(2,i-1,j  ,k,n)*sj(4,i-1,j  ,k,n) + &
                            sj(2,i  ,j  ,k,n)*sj(4,i  ,j  ,k,n) + &
                            sj(2,i-1,j+1,k,n)*sj(4,i-1,j+1,k,n) + &
                            sj(2,i  ,j+1,k,n)*sj(4,i  ,j+1,k,n) ) / 4 / &
                          volvisc
              etaz      = ( sj(3,i-1,j  ,k,n)*sj(4,i-1,j  ,k,n) + &
                            sj(3,i  ,j  ,k,n)*sj(4,i  ,j  ,k,n) + &
                            sj(3,i-1,j+1,k,n)*sj(4,i-1,j+1,k,n) + &
                            sj(3,i  ,j+1,k,n)*sj(4,i  ,j+1,k,n) ) / 4 / &
                          volvisc

              zetax     = ( sk(1,i-1,j,k  ,n)*sk(4,i-1,j,k  ,n) + &
                            sk(1,i  ,j,k  ,n)*sk(4,i  ,j,k  ,n) + &
                            sk(1,i-1,j,k+1,n)*sk(4,i-1,j,k+1,n) + &
                            sk(1,i  ,j,k+1,n)*sk(4,i  ,j,k+1,n) ) / 4 / &
                          volvisc
              zetay     = ( sk(2,i-1,j,k  ,n)*sk(4,i-1,j,k  ,n) + &
                            sk(2,i  ,j,k  ,n)*sk(4,i  ,j,k  ,n) + &
                            sk(2,i-1,j,k+1,n)*sk(4,i-1,j,k+1,n) + &
                            sk(2,i  ,j,k+1,n)*sk(4,i  ,j,k+1,n) ) / 4 / &
                          volvisc
              zetaz     = ( sk(3,i-1,j,k  ,n)*sk(4,i-1,j,k  ,n) + &
                            sk(3,i  ,j,k  ,n)*sk(4,i  ,j,k  ,n) + &
                            sk(3,i-1,j,k+1,n)*sk(4,i-1,j,k+1,n) + &
                            sk(3,i  ,j,k+1,n)*sk(4,i  ,j,k+1,n) ) / 4 / &
                          volvisc

              xideri  (1,i,j,k,n) = xix
              xideri  (2,i,j,k,n) = xiy
              xideri  (3,i,j,k,n) = xiz
              etaderi (1,i,j,k,n) = etax
              etaderi (2,i,j,k,n) = etay
              etaderi (3,i,j,k,n) = etaz
              zetaderi(1,i,j,k,n) = zetax
              zetaderi(2,i,j,k,n) = zetay
              zetaderi(3,i,j,k,n) = zetaz

            end do
          end do
        end do

!------ eta-face metrics
        do k = 1,kmax-1
          do j = 1,jmax
            do i = 1,imax-1

              volvisc   = (vol(i,j-1,k)+vol(i,j,k)) / 2

              xix       = ( si(1,i  ,j-1,k,n)*si(4,i  ,j-1,k,n) + &
                            si(1,i  ,j  ,k,n)*si(4,i  ,j  ,k,n) + &
                            si(1,i+1,j-1,k,n)*si(4,i+1,j-1,k,n) + &
                            si(1,i+1,j  ,k,n)*si(4,i+1,j  ,k,n) ) / 4 / &
                          volvisc
              xiy       = ( si(2,i  ,j-1,k,n)*si(4,i  ,j-1,k,n) + &
                            si(2,i  ,j  ,k,n)*si(4,i  ,j  ,k,n) + &
                            si(2,i+1,j-1,k,n)*si(4,i+1,j-1,k,n) + &
                            si(2,i+1,j  ,k,n)*si(4,i+1,j  ,k,n) ) / 4 / &
                          volvisc
              xiz       = ( si(3,i  ,j-1,k,n)*si(4,i  ,j-1,k,n) + &
                            si(3,i  ,j  ,k,n)*si(4,i  ,j  ,k,n) + &
                            si(3,i+1,j-1,k,n)*si(4,i+1,j-1,k,n) + &
                            si(3,i+1,j  ,k,n)*si(4,i+1,j  ,k,n) ) / 4 / &
                          volvisc

              etax      =   sj(1,i,j,k,n)*sj(4,i,j,k,n) / volvisc
              etay      =   sj(2,i,j,k,n)*sj(4,i,j,k,n) / volvisc
              etaz      =   sj(3,i,j,k,n)*sj(4,i,j,k,n) / volvisc

              zetax     = ( sk(1,i,j-1,k  ,n)*sk(4,i,j-1,k  ,n) + &
                            sk(1,i,j  ,k  ,n)*sk(4,i,j  ,k  ,n) + &
                            sk(1,i,j-1,k+1,n)*sk(4,i,j-1,k+1,n) + &
                            sk(1,i,j  ,k+1,n)*sk(4,i,j  ,k+1,n) ) / 4 / &
                          volvisc
              zetay     = ( sk(2,i,j-1,k  ,n)*sk(4,i,j-1,k  ,n) + &
                            sk(2,i,j  ,k  ,n)*sk(4,i,j  ,k  ,n) + &
                            sk(2,i,j-1,k+1,n)*sk(4,i,j-1,k+1,n) + &
                            sk(2,i,j  ,k+1,n)*sk(4,i,j  ,k+1,n) ) / 4 / &
                          volvisc
              zetaz     = ( sk(3,i,j-1,k  ,n)*sk(4,i,j-1,k  ,n) + &
                            sk(3,i,j  ,k  ,n)*sk(4,i,j  ,k  ,n) + &
                            sk(3,i,j-1,k+1,n)*sk(4,i,j-1,k+1,n) + &
                            sk(3,i,j  ,k+1,n)*sk(4,i,j  ,k+1,n) ) / 4 / &
                          volvisc

              xiderj  (1,i,j,k,n) = xix
              xiderj  (2,i,j,k,n) = xiy
              xiderj  (3,i,j,k,n) = xiz
              etaderj (1,i,j,k,n) = etax
              etaderj (2,i,j,k,n) = etay
              etaderj (3,i,j,k,n) = etaz
              zetaderj(1,i,j,k,n) = zetax
              zetaderj(2,i,j,k,n) = zetay
              zetaderj(3,i,j,k,n) = zetaz

            end do
          end do
        end do

!------ zeta-face metrics
        do k = 1,kmax
          do j = 1,jmax-1
            do i = 1,imax-1

              volvisc   = (vol(i,j,k-1)+vol(i,j,k)) / 2

              xix       = ( si(1,i  ,j,k-1,n)*si(4,i  ,j,k-1,n) + &
                            si(1,i  ,j,k  ,n)*si(4,i  ,j,k  ,n) + &
                            si(1,i+1,j,k-1,n)*si(4,i+1,j,k-1,n) + &
                            si(1,i+1,j,k  ,n)*si(4,i+1,j,k  ,n) ) / 4 / &
                          volvisc
              xiy       = ( si(2,i  ,j,k-1,n)*si(4,i  ,j,k-1,n) + &
                            si(2,i  ,j,k  ,n)*si(4,i  ,j,k  ,n) + &
                            si(2,i+1,j,k-1,n)*si(4,i+1,j,k-1,n) + &
                            si(2,i+1,j,k  ,n)*si(4,i+1,j,k  ,n) ) / 4 / &
                          volvisc
              xiz       = ( si(3,i  ,j,k-1,n)*si(4,i  ,j,k-1,n) + &
                            si(3,i  ,j,k  ,n)*si(4,i  ,j,k  ,n) + &
                            si(3,i+1,j,k-1,n)*si(4,i+1,j,k-1,n) + &
                            si(3,i+1,j,k  ,n)*si(4,i+1,j,k  ,n) ) / 4 / &
                          volvisc

              etax      = ( sj(1,i,j  ,k-1,n)*sj(4,i,j  ,k-1,n) + &
                            sj(1,i,j  ,k  ,n)*sj(4,i,j  ,k  ,n) + &
                            sj(1,i,j+1,k-1,n)*sj(4,i,j+1,k-1,n) + &
                            sj(1,i,j+1,k  ,n)*sj(4,i,j+1,k  ,n) ) / 4 / &
                          volvisc
              etay      = ( sj(2,i,j  ,k-1,n)*sj(4,i,j  ,k-1,n) + &
                            sj(2,i,j  ,k  ,n)*sj(4,i,j  ,k  ,n) + &
                            sj(2,i,j+1,k-1,n)*sj(4,i,j+1,k-1,n) + &
                            sj(2,i,j+1,k  ,n)*sj(4,i,j+1,k  ,n) ) / 4 / &
                          volvisc
              etaz      = ( sj(3,i,j  ,k-1,n)*sj(4,i,j  ,k-1,n) + &
                            sj(3,i,j  ,k  ,n)*sj(4,i,j  ,k  ,n) + &
                            sj(3,i,j+1,k-1,n)*sj(4,i,j+1,k-1,n) + &
                            sj(3,i,j+1,k  ,n)*sj(4,i,j+1,k  ,n) ) / 4 / &
                          volvisc

              zetax     =   sk(1,i,j,k,n)*sk(4,i,j,k,n) / volvisc
              zetay     =   sk(2,i,j,k,n)*sk(4,i,j,k,n) / volvisc
              zetaz     =   sk(3,i,j,k,n)*sk(4,i,j,k,n) / volvisc

              xiderk  (1,i,j,k,n) = xix
              xiderk  (2,i,j,k,n) = xiy
              xiderk  (3,i,j,k,n) = xiz
              etaderk (1,i,j,k,n) = etax
              etaderk (2,i,j,k,n) = etay
              etaderk (3,i,j,k,n) = etaz
              zetaderk(1,i,j,k,n) = zetax
              zetaderk(2,i,j,k,n) = zetay
              zetaderk(3,i,j,k,n) = zetaz

            end do
          end do
        end do

      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine movegrid(nl,x,y,z,xdot,ydot,zdot,dx,dy,dz,x0,y0,z0)
!-----------------------------------------------------------------------

      use cosa_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,iblock,imax,jmax,kmax
      integer(kind=cosa_int) ixyz,ixyz0,idxyz
      real(kind=cosa_real) x(*),y(*),z(*),x0(*),y0(*),z0(*),dx(*),dy(*),dz(*), &
                   xdot(*),ydot(*),zdot(*)

      do iblock = 1,mynblocks
        imax  = i_imax  (iblock,nl)
        jmax  = j_jmax  (iblock,nl)
        kmax  = k_kmax  (iblock,nl)
        ixyz  = 1 + off_p2(iblock,nl)         * dim5h
        ixyz0 = 1 + off_p2(iblock,nl)
        idxyz = 1 + off_p1(iblock,nl)
        call movegrid_b(x(ixyz),y(ixyz),z(ixyz),dx(idxyz),dy(idxyz), &
          dz(idxyz),x0(ixyz0),y0(ixyz0),z0(ixyz0),xdot(ixyz),ydot(ixyz), &
          zdot(ixyz),imax,jmax,kmax,nharms)
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine movegrid_b(x,y,z,dx,dy,dz,x0,y0,z0,xdot,ydot,zdot,imax, &
                            jmax,kmax,nharms)
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,nharms
      integer(kind=cosa_int) n
      real(kind=cosa_real) &
          x    (0:imax+1,0:jmax+1,0:kmax+1,0:2*nharms*hbmove), &
          y    (0:imax+1,0:jmax+1,0:kmax+1,0:2*nharms*hbmove), &
          z    (0:imax+1,0:jmax+1,0:kmax+1,0:2*nharms*hbmove), &
          xdot (0:imax+1,0:jmax+1,0:kmax+1,0:2*nharms*hbmove), &
          ydot (0:imax+1,0:jmax+1,0:kmax+1,0:2*nharms*hbmove), &
          zdot (0:imax+1,0:jmax+1,0:kmax+1,0:2*nharms*hbmove), &
          x0   (0:imax+1,0:jmax+1,0:kmax+1                   ), &
          y0   (0:imax+1,0:jmax+1,0:kmax+1                   ), &
          z0   (0:imax+1,0:jmax+1,0:kmax+1                   ), &
          dx   (0:imax  ,0:jmax  ,0:kmax                     ), &
          dy   (0:imax  ,0:jmax  ,0:kmax                     ), &
          dz   (0:imax  ,0:jmax  ,0:kmax                     )

!msc- 03 January 2015: from version 05_1 onwards, the following
!                      statement has been moved to the dualt main loop
!                      in main of cosa.f. This was done for need of
!                      zeit(0) in bc_far0 when using dualt for hawt
!                      rotors and sectors in relative frame.

      do n=0,2*nharms*hbmove

        if (pitching) then
          call pitch(zeit(n),dtheta(n),dtheta0,omegas,xh,yh, &
            x0,y0,z0,dx,dy,dz,x(0,0,0,n),y(0,0,0,n),z(0,0,0,n), &
            xdot(0,0,0,n),ydot(0,0,0,n),zdot(0,0,0,n),imax,jmax,kmax)
        else if (plunging) then
          call plunge(zeit(n),dh0x,dh0y,omegas,x0,y0,z0,x(0,0,0,n), &
            y(0,0,0,n),z(0,0,0,n),xdot(0,0,0,n),ydot(0,0,0,n), &
            zdot(0,0,0,n),imax,jmax,kmax)
        else if (plupitching) then
          call plupitch(zeit(n),dtheta(n),dtheta0,omegas,xh,yh, &
            dh0x,dh0y,phipp,x0,y0,z0,x(0,0,0,n),y(0,0,0,n),z(0,0,0,n), &
            xdot(0,0,0,n),ydot(0,0,0,n),zdot(0,0,0,n),imax,jmax,kmax)
        else if (rotating) then
          call rotate(zeit(n),dtheta(n),thetatp(n),theta0tp,xrotc, &
            yrotc,yhtp,zhtp,omegas,omegatp,betaw,phitp,x0,y0,z0, &
            x(0,0,0,n),y(0,0,0,n),z(0,0,0,n),xdot(0,0,0,n), &
            ydot(0,0,0,n),zdot(0,0,0,n),imax,jmax,kmax,relframe,tpitch)
        end if

      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine pitch(zeit,dtheta,dtheta0,omegas,xh,yh,x0,y0,z0,dx,dy, &
        dz,x,y,z,xdot,ydot,zdot,imax,jmax,kmax)
!-----------------------------------------------------------------------
!     calculates time-dependent grid position and velocity when the whole
!     grid undergoes a Rigid Body Pitching Motion.
!
!     REMARK 1: it may appear silly to divide and multiply by dx or dy.
!     This is done in view of the deforming grid implementation, whereby 
!     disp_tfact and vel_tfact will ONLY depend on time, whereas dx and 
!     dy will come from the mass-spring analogy.
!
!     REMARK 2: for validation I consider the case of a pitching box in 
!     a uniform flow. In this case, the centre of rotation belongs to
!     the domain. The dx or dy close to the hinge will be close to zero.
!     In order to avoid division by zero, I put the tests on the size of
!     the displacements.
!-----------------------------------------------------------------------

      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax
      integer(kind=cosa_int) i,j,k
      real(kind=cosa_real) &
          x    (0:imax+1,0:jmax+1,0:kmax+1), &
          y    (0:imax+1,0:jmax+1,0:kmax+1), &
          z    (0:imax+1,0:jmax+1,0:kmax+1), &
          xdot (0:imax+1,0:jmax+1,0:kmax+1), &
          ydot (0:imax+1,0:jmax+1,0:kmax+1), &
          zdot (0:imax+1,0:jmax+1,0:kmax+1), &
          x0   (0:imax+1,0:jmax+1,0:kmax+1), &
          y0   (0:imax+1,0:jmax+1,0:kmax+1), &
          z0   (0:imax+1,0:jmax+1,0:kmax+1), &
          dx   (0:imax  ,0:jmax  ,0:kmax  ), &
          dy   (0:imax  ,0:jmax  ,0:kmax  ), &
          dz   (0:imax  ,0:jmax  ,0:kmax  )
      real(kind=cosa_real) zeit,dtheta,dtheta0,omegas,xh,yh,disp_tfac(2), &
          vel_tfac(2)

      do k = 1,kmax
        do i = 1,imax
          do j = 1,jmax

            z   (i,j,k) = z0(i,j,k)
            zdot(i,j,k) = 0.d0

            if ((abs(dx(i,j,k)).gt.1.d-4).and.(abs(dy(i,j,k)).gt.1.d-4)) &
                                                                    then

              disp_tfac(1) = ( (x0(i,j,k) - xh) * (dcos(dtheta) - 1) - &
                               (y0(i,j,k) - yh) *  dsin(dtheta)      ) / &
                             dx(i,j,k)
              disp_tfac(2) = ( (x0(i,j,k) - xh) *  dsin(dtheta)      + &
                               (y0(i,j,k) - yh) * (dcos(dtheta) - 1) ) / &
                             dy(i,j,k)
              x(i,j,k) = x0(i,j,k) + disp_tfac(1) * dx(i,j,k)
              y(i,j,k) = y0(i,j,k) + disp_tfac(2) * dy(i,j,k)

              vel_tfac(1)  = (-(x0(i,j,k) - xh) * dsin(dtheta) - &
                               (y0(i,j,k) - yh) * dcos(dtheta) ) * &
                             dtheta0*dcos(omegas*zeit) / dx(i,j,k)
              vel_tfac(2)  = ( (x0(i,j,k) - xh) * dcos(dtheta) - &
                               (y0(i,j,k) - yh) * dsin(dtheta) ) * &
                             dtheta0*dcos(omegas*zeit) / dy(i,j,k)
              xdot(i,j,k) = omegas * vel_tfac(1) * dx(i,j,k)
              ydot(i,j,k) = omegas * vel_tfac(2) * dy(i,j,k)

            else

              disp_tfac(1) = ( (x0(i,j,k) -xh) * (dcos(dtheta) - 1.d0) - &
                          y0(i,j,k) * dsin(dtheta) )
              disp_tfac(2) = ( (x0(i,j,k) -xh) * dsin(dtheta) + &
                          y0(i,j,k) * (dcos(dtheta) - 1.d0) )
              x(i,j,k) = x0(i,j,k) + disp_tfac(1)
              y(i,j,k) = y0(i,j,k) + disp_tfac(2)
              vel_tfac(1) = dtheta0*dcos(omegas*zeit)* &
                           ( -(x0(i,j,k) -xh) * dsin(dtheta) - &
                             y0(i,j,k) * dcos(dtheta) )
              vel_tfac(2) = dtheta0*dcos(omegas*zeit)* &
                           ( (x0(i,j,k) -xh) * dcos(dtheta) - &
                             y0(i,j,k) * dsin(dtheta) )
              xdot(i,j,k) = omegas * vel_tfac(1)
              ydot(i,j,k) = omegas * vel_tfac(2)

            end if

          end do
        end do
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine plunge(zeit,dh0x,dh0y,omegas,x0,y0,z0,x,y,z,xdot,ydot, &
        zdot,imax,jmax,kmax)
!-----------------------------------------------------------------------
!     calculates time-dependent grid position and velocity when the whole
!     grid undergoes a Rigid Body Plunging Motion.
!-----------------------------------------------------------------------

      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax
      integer(kind=cosa_int) i,j,k
      real(kind=cosa_real) &
          x    (0:imax+1,0:jmax+1,0:kmax+1), &
          y    (0:imax+1,0:jmax+1,0:kmax+1), &
          z    (0:imax+1,0:jmax+1,0:kmax+1), &
          xdot (0:imax+1,0:jmax+1,0:kmax+1), &
          ydot (0:imax+1,0:jmax+1,0:kmax+1), &
          zdot (0:imax+1,0:jmax+1,0:kmax+1), &
          x0   (0:imax+1,0:jmax+1,0:kmax+1), &
          y0   (0:imax+1,0:jmax+1,0:kmax+1), &
          z0   (0:imax+1,0:jmax+1,0:kmax+1)
      real(kind=cosa_real) zeit,dh0x,dh0y,omegas,disp_tfac(2),vel_tfac(2)

      disp_tfac(1) = dsin(omegas*zeit)
      disp_tfac(2) = dsin(omegas*zeit)
      vel_tfac(1)  = dcos(omegas*zeit)
      vel_tfac(2)  = dcos(omegas*zeit)

      do k = 1,kmax
        do j = 1,jmax
          do i = 1,imax
            x   (i,j,k) = x0(i,j,k) + dh0x*disp_tfac(1)
            y   (i,j,k) = y0(i,j,k) + dh0y*disp_tfac(2)
            z   (i,j,k) = z0(i,j,k)
            xdot(i,j,k) = omegas*vel_tfac(1)*dh0x
            ydot(i,j,k) = omegas*vel_tfac(2)*dh0y
            zdot(i,j,k) = 0.d0
          end do
        end do
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine rotate(zeit,theta,thetatp,theta0tp,xrotc,yrotc,yhtp, &
        zhtp,omegas,omegatp,betaw,phitp,x0,y0,z0,x,y,z,xdot, &
        ydot,zdot,imax,jmax,kmax,relframe,tpitch)
!-----------------------------------------------------------------------
!     calculates time-dependent grid position and velocity when the whole
!     grid undergoes a Rigid Body Plunging Motion.
!-----------------------------------------------------------------------

      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax
      integer(kind=cosa_int) i,j,k
      real (kind=cosa_real) &
          x    (0:imax+1,0:jmax+1,0:kmax+1), &
          y    (0:imax+1,0:jmax+1,0:kmax+1), &
          z    (0:imax+1,0:jmax+1,0:kmax+1), &
          xdot (0:imax+1,0:jmax+1,0:kmax+1), &
          ydot (0:imax+1,0:jmax+1,0:kmax+1), &
          zdot (0:imax+1,0:jmax+1,0:kmax+1), &
          x0   (0:imax+1,0:jmax+1,0:kmax+1), &
          y0   (0:imax+1,0:jmax+1,0:kmax+1), &
          z0   (0:imax+1,0:jmax+1,0:kmax+1)
      real (kind=cosa_real) zeit,theta,thetatp,theta0tp,xrotc,yrotc,yhtp,zhtp, &
        xrtw,yrtw,zrtw,xrtwdot,yrtwdot,zrtwdot,omegas,omegatp,betaw, &
        phitp,disp_tfac(3),vel_tfac(3),cw,sw,ctp,stp
      logical relframe,tpitch

      if (.not.relframe) then

        do k = 1,kmax
          do j = 1,jmax
            do i = 1,imax

              disp_tfac(1) = ( (x0(i,j,k) - xrotc) * (dcos(theta) - 1) - &
                               (y0(i,j,k) - yrotc) *  dsin(theta)      )
              disp_tfac(2) = ( (x0(i,j,k) - xrotc) *  dsin(theta)      + &
                               (y0(i,j,k) - yrotc) * (dcos(theta) - 1) )

              x(i,j,k) = x0(i,j,k) + disp_tfac(1)
              y(i,j,k) = y0(i,j,k) + disp_tfac(2)
              z(i,j,k) = z0(i,j,k)

              vel_tfac(1)  = (-(x0(i,j,k) - xrotc) * dsin(theta) - &
                               (y0(i,j,k) - yrotc) * dcos(theta) )
              vel_tfac(2)  = ( (x0(i,j,k) - xrotc) * dcos(theta) - &
                               (y0(i,j,k) - yrotc) * dsin(theta) )
              xdot(i,j,k) = omegas * vel_tfac(1)
              ydot(i,j,k) = omegas * vel_tfac(2)
              zdot(i,j,k) = 0.d0

              if (tpitch) then

                cw      = dcos(betaw)
                sw      = dsin(betaw)
                xrtw    = cw*x(i,j,k) - sw*z(i,j,k)
                yrtw    =    y(i,j,k)
                zrtw    = sw*x(i,j,k) + cw*z(i,j,k)
                xrtwdot = cw*xdot(i,j,k)
                yrtwdot =    ydot(i,j,k)
                zrtwdot = sw*xdot(i,j,k)

                ctp     = dcos(thetatp)
                stp     = dsin(thetatp)
                disp_tfac(1) = xrtw
                disp_tfac(2) = yrtw + &
                  (-(zrtw - zhtp) * stp  + (yrtw - yhtp) * (ctp - 1) )
                disp_tfac(3) = zrtw + &
                  ( (zrtw - zhtp) * (ctp - 1) + (yrtw - yhtp) * stp )
                vel_tfac(1) = xrtwdot
                vel_tfac(2)  = &
                  (-(zrtw - zhtp) * ctp - (yrtw - yhtp) * stp ) * &
                    omegatp * theta0tp * dcos(omegatp*zeit+phitp) + &
                  yrtwdot * ctp - zrtwdot * stp
                vel_tfac(3)  = &
                  (-(zrtw - zhtp) * stp + (yrtw - yhtp) * ctp ) * &
                    omegatp * theta0tp * dcos(omegatp*zeit+phitp) + &
                  yrtwdot * stp + zrtwdot * ctp

                x(i,j,k)    = sw*disp_tfac(3) + cw*disp_tfac(1)
                y(i,j,k)    = disp_tfac(2)
                z(i,j,k)    = cw*disp_tfac(3) - sw*disp_tfac(1)
                xdot(i,j,k) = sw*vel_tfac(3) + cw*vel_tfac(1)
                ydot(i,j,k) = vel_tfac(2)
                zdot(i,j,k) = cw*vel_tfac(3) - sw*vel_tfac(1)

              end if

            end do
          end do
        end do

      else

        do k = 1,kmax
          do j = 1,jmax
            do i = 1,imax

              x(i,j,k) = x0(i,j,k)
              y(i,j,k) = y0(i,j,k)
              z(i,j,k) = z0(i,j,k)

              vel_tfac(1)  = -(y0(i,j,k) - yrotc)
              vel_tfac(2)  =  (x0(i,j,k) - xrotc)
              xdot(i,j,k) = omegas * vel_tfac(1)
              ydot(i,j,k) = omegas * vel_tfac(2)
              zdot(i,j,k) = 0.d0

            end do
          end do
        end do

      end if

      return
      end

!-----------------------------------------------------------------------
      subroutine vols(nl,x,y,z,xdot,ydot,zdot,vol,rad)
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,imax,jmax,kmax
      integer(kind=cosa_int) iblock,ixyz,ivol,irad
      real (kind=cosa_real) x(*),y(*),z(*),xdot(*),ydot(*),zdot(*),vol(*),rad(*)

      do iblock = 1,mynblocks
        imax  = i_imax  (iblock,nl)
        jmax  = j_jmax  (iblock,nl)
        kmax  = k_kmax  (iblock,nl)
        ixyz  = 1 + off_p2(iblock,nl) * dim5h
        ivol  = 1 + off_p1(iblock,nl)
        irad  = 1 + off_m1(iblock,nl) * dim5h
        call bvols(iblock,x(ixyz),y(ixyz),z(ixyz),xdot(ixyz),ydot(ixyz), &
                   zdot(ixyz),vol(ivol),rad(irad),imax,jmax,kmax,nharms, &
                   lowernblock)
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine bvols(iblock,x,y,z,xdot,ydot,zdot,vol,rad,imax,jmax, &
                       kmax,nharms,lowernblock)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,nharms
      integer(kind=cosa_int) iblock,i,j,k,n,lowernblock,liblock
      real (kind=cosa_real) &
          x   ( 0:imax+1,0:jmax+1,0:kmax+1,0:2*nharms*hbmove), &
          y   ( 0:imax+1,0:jmax+1,0:kmax+1,0:2*nharms*hbmove), &
          z   ( 0:imax+1,0:jmax+1,0:kmax+1,0:2*nharms*hbmove), &
          xdot( 0:imax+1,0:jmax+1,0:kmax+1,0:2*nharms*hbmove), &
          ydot( 0:imax+1,0:jmax+1,0:kmax+1,0:2*nharms*hbmove), &
          zdot( 0:imax+1,0:jmax+1,0:kmax+1,0:2*nharms*hbmove), &
          vol ( 0:imax  ,0:jmax  ,0:kmax), &
          rad (   imax-1,  jmax-1,  kmax-1,0:2*nharms*hbmove)
      real (kind=cosa_real) si(3), sj(3), sk(3), dx, dy, dz, sidot(3), sjdot(3), &
           skdot(3), dxdot, dydot, dzdot, voldot, radl(8)

      liblock = lowernblock+iblock-1

!---- Calculate the cell volume from eqn. 2c of AIAA J. Vol. 21
!     No. 6 pp. 917-918 by Kordulla and Vinokur

!---- MSC: for now grid does not deform, so volumes do not depend on n
      n = 0

      do k = 1, kmax-1
        do j = 1, jmax-1
          do i = 1, imax-1
            si(1) = ( (y(i,j+1,k+1,n) - y(i,j  ,k,n)) * &
                      (z(i,j  ,k+1,n) - z(i,j+1,k,n)) - &
                      (z(i,j+1,k+1,n) - z(i,j  ,k,n)) * &
                      (y(i,j  ,k+1,n) - y(i,j+1,k,n)) ) / 2
            sj(1) = ( (z(i+1,j,k+1,n) - z(i  ,j,k,n)) * &
                      (y(i  ,j,k+1,n) - y(i+1,j,k,n)) - &
                      (y(i+1,j,k+1,n) - y(i  ,j,k,n)) * &
                      (z(i  ,j,k+1,n) - z(i+1,j,k,n)) ) / 2
            sk(1) = ( (y(i+1,j+1,k,n) - y(i  ,j,k,n)) * &
                      (z(i  ,j+1,k,n) - z(i+1,j,k,n)) - &
                      (z(i+1,j+1,k,n) - z(i  ,j,k,n)) * &
                      (y(i  ,j+1,k,n) - y(i+1,j,k,n)) ) / 2
            dx = ( si(1) + sj(1) + sk(1) ) * &
                 ( x(i+1,j+1,k+1,n) - x(i  ,j  ,k  ,n) )

            si(2) = ( (z(i,j+1,k+1,n) - z(i,j  ,k,n)) * &
                      (x(i,j  ,k+1,n) - x(i,j+1,k,n)) - &
                      (x(i,j+1,k+1,n) - x(i,j  ,k,n)) * &
                      (z(i,j  ,k+1,n) - z(i,j+1,k,n)) ) / 2
            sj(2) = ( (x(i+1,j,k+1,n) - x(i  ,j,k,n)) * &
                      (z(i  ,j,k+1,n) - z(i+1,j,k,n)) - &
                      (z(i+1,j,k+1,n) - z(i  ,j,k,n)) * &
                      (x(i  ,j,k+1,n) - x(i+1,j,k,n)) ) / 2
            sk(2) = ( (z(i+1,j+1,k,n) - z(i  ,j,k,n)) * &
                      (x(i  ,j+1,k,n) - x(i+1,j,k,n)) - &
                      (x(i+1,j+1,k,n) - x(i  ,j,k,n)) * &
                      (z(i  ,j+1,k,n) - z(i+1,j,k,n)) ) / 2
            dy = ( si(2) + sj(2) + sk(2) ) * &
                 ( y(i+1,j+1,k+1,n) - y(i  ,j  ,k  ,n) )

            si(3) = ( (x(i,j+1,k+1,n) - x(i,j  ,k,n)) * &
                      (y(i,j  ,k+1,n) - y(i,j+1,k,n)) - &
                      (y(i,j+1,k+1,n) - y(i,j  ,k,n)) * &
                      (x(i,j  ,k+1,n) - x(i,j+1,k,n)) ) / 2
            sj(3) = ( (y(i+1,j,k+1,n) - y(i  ,j,k,n)) * &
                      (x(i  ,j,k+1,n) - x(i+1,j,k,n)) - &
                      (x(i+1,j,k+1,n) - x(i  ,j,k,n)) * &
                      (y(i  ,j,k+1,n) - y(i+1,j,k,n)) ) / 2
            sk(3) = ( (x(i+1,j+1,k,n) - x(i  ,j,k,n)) * &
                      (y(i  ,j+1,k,n) - y(i+1,j,k,n)) - &
                      (y(i+1,j+1,k,n) - y(i  ,j,k,n)) * &
                      (x(i  ,j+1,k,n) - x(i+1,j,k,n)) ) / 2
            dz = ( si(3) + sj(3) + sk(3) ) * &
                 ( z(i+1,j+1,k+1,n) - z(i  ,j  ,k  ,n) )

            vol(i,j,k) = (dx + dy + dz) / 3

            if (vol(i,j,k).le.0.d0) then
!AJ              write(2,100)
!AJ              write(2,*) liblock,i,j,k
!AJ              write(*,100)
!AJ              write(*,*) liblock,i,j,k
!AJ              write(*,*) dx,dy,dz
#if MPI
!AJ              call abortmpi()
#else
!AJ              stop
#endif
            end if

          end do
        end do
      end do

      if (debug.and.moving) then
!------ make sure that rate of volume change is zero at all times
        do n=0,2*nharms*hbmove
          do k = 1, kmax-1
            do j = 1, jmax-1
              do i = 1, imax-1
                si(1) = ( (y(i,j+1,k+1,n) - y(i,j  ,k,n)) * &
                          (z(i,j  ,k+1,n) - z(i,j+1,k,n)) - &
                          (z(i,j+1,k+1,n) - z(i,j  ,k,n)) * &
                          (y(i,j  ,k+1,n) - y(i,j+1,k,n)) ) / 2
                sj(1) = ( (z(i+1,j,k+1,n) - z(i  ,j,k,n)) * &
                          (y(i  ,j,k+1,n) - y(i+1,j,k,n)) - &
                          (y(i+1,j,k+1,n) - y(i  ,j,k,n)) * &
                          (z(i  ,j,k+1,n) - z(i+1,j,k,n)) ) / 2
                sk(1) = ( (y(i+1,j+1,k,n) - y(i  ,j,k,n)) * &
                          (z(i  ,j+1,k,n) - z(i+1,j,k,n)) - &
                          (z(i+1,j+1,k,n) - z(i  ,j,k,n)) * &
                          (y(i  ,j+1,k,n) - y(i+1,j,k,n)) ) / 2
                sidot(1) = ( ((ydot(i,j+1,k+1,n) - ydot(i,j  ,k,n)) * &
                              (z(i,j  ,k+1,n) - z(i,j+1,k,n)) + &
                              (y(i,j+1,k+1,n) - y(i,j  ,k,n)) * &
                              (zdot(i,j  ,k+1,n) - zdot(i,j+1,k,n)))- &
                             ((zdot(i,j+1,k+1,n) - zdot(i,j  ,k,n)) * &
                              (y(i,j  ,k+1,n) - y(i,j+1,k,n)) + &
                              (z(i,j+1,k+1,n) - z(i,j  ,k,n)) * &
                              (ydot(i,j  ,k+1,n) - ydot(i,j+1,k,n))) )/2
                sjdot(1) = ( ((zdot(i+1,j,k+1,n) - zdot(i  ,j,k,n)) * &
                              (y(i  ,j,k+1,n) - y(i+1,j,k,n)) + &
                              (z(i+1,j,k+1,n) - z(i  ,j,k,n)) * &
                              (ydot(i  ,j,k+1,n) - ydot(i+1,j,k,n))) - &
                             ((ydot(i+1,j,k+1,n) - ydot(i  ,j,k,n)) * &
                              (z(i  ,j,k+1,n) - z(i+1,j,k,n)) + &
                              (y(i+1,j,k+1,n) - y(i  ,j,k,n)) * &
                              (zdot(i  ,j,k+1,n) - zdot(i+1,j,k,n))) )/2
                skdot(1) = ( ((ydot(i+1,j+1,k,n) - ydot(i  ,j,k,n)) * &
                              (z(i  ,j+1,k,n) - z(i+1,j,k,n)) + &
                              (y(i+1,j+1,k,n) - y(i  ,j,k,n)) * &
                              (zdot(i  ,j+1,k,n) - zdot(i+1,j,k,n))) - &
                             ((zdot(i+1,j+1,k,n) - zdot(i  ,j,k,n)) * &
                              (y(i  ,j+1,k,n) - y(i+1,j,k,n)) + &
                              (z(i+1,j+1,k,n) - z(i  ,j,k,n)) * &
                              (ydot(i  ,j+1,k,n) - ydot(i+1,j,k,n))) )/2
            dxdot = ( sidot(1) + sjdot(1) + skdot(1) ) * &
                    ( x   (i+1,j+1,k+1,n) - x   (i  ,j  ,k  ,n) ) + &
                    ( si   (1) + sj   (1) + sk   (1) ) * &
                    ( xdot(i+1,j+1,k+1,n) - xdot(i  ,j  ,k  ,n) )

                si(2) = ( (z(i,j+1,k+1,n) - z(i,j  ,k,n)) * &
                          (x(i,j  ,k+1,n) - x(i,j+1,k,n)) - &
                          (x(i,j+1,k+1,n) - x(i,j  ,k,n)) * &
                          (z(i,j  ,k+1,n) - z(i,j+1,k,n)) ) / 2
                sj(2) = ( (x(i+1,j,k+1,n) - x(i  ,j,k,n)) * &
                          (z(i  ,j,k+1,n) - z(i+1,j,k,n)) - &
                          (z(i+1,j,k+1,n) - z(i  ,j,k,n)) * &
                          (x(i  ,j,k+1,n) - x(i+1,j,k,n)) ) / 2
                sk(2) = ( (z(i+1,j+1,k,n) - z(i  ,j,k,n)) * &
                          (x(i  ,j+1,k,n) - x(i+1,j,k,n)) - &
                          (x(i+1,j+1,k,n) - x(i  ,j,k,n)) * &
                          (z(i  ,j+1,k,n) - z(i+1,j,k,n)) ) / 2
                sidot(2) = ( ((zdot(i,j+1,k+1,n) - zdot(i,j  ,k,n)) * &
                              (x(i,j  ,k+1,n) - x(i,j+1,k,n)) + &
                              (z(i,j+1,k+1,n) - z(i,j  ,k,n)) * &
                              (xdot(i,j  ,k+1,n) - xdot(i,j+1,k,n)))- &
                             ((xdot(i,j+1,k+1,n) - xdot(i,j  ,k,n)) * &
                              (z(i,j  ,k+1,n) - z(i,j+1,k,n)) + &
                              (x(i,j+1,k+1,n) - x(i,j  ,k,n)) * &
                              (zdot(i,j  ,k+1,n) - zdot(i,j+1,k,n))) )/2
                sjdot(2) = ( ((xdot(i+1,j,k+1,n) - xdot(i  ,j,k,n)) * &
                              (z(i  ,j,k+1,n) - z(i+1,j,k,n)) + &
                              (x(i+1,j,k+1,n) - x(i  ,j,k,n)) * &
                              (zdot(i  ,j,k+1,n) - zdot(i+1,j,k,n))) - &
                             ((zdot(i+1,j,k+1,n) - zdot(i  ,j,k,n)) * &
                              (x(i  ,j,k+1,n) - x(i+1,j,k,n)) + &
                              (z(i+1,j,k+1,n) - z(i  ,j,k,n)) * &
                              (xdot(i  ,j,k+1,n) - xdot(i+1,j,k,n))) )/2
                skdot(2) = ( ((zdot(i+1,j+1,k,n) - zdot(i  ,j,k,n)) * &
                              (x(i  ,j+1,k,n) - x(i+1,j,k,n)) + &
                              (z(i+1,j+1,k,n) - z(i  ,j,k,n)) * &
                              (xdot(i  ,j+1,k,n) - xdot(i+1,j,k,n))) - &
                             ((xdot(i+1,j+1,k,n) - xdot(i  ,j,k,n)) * &
                              (z(i  ,j+1,k,n) - z(i+1,j,k,n)) + &
                              (x(i+1,j+1,k,n) - x(i  ,j,k,n)) * &
                              (zdot(i  ,j+1,k,n) - zdot(i+1,j,k,n))) )/2
            dydot = ( sidot(2) + sjdot(2) + skdot(2) ) * &
                    ( y   (i+1,j+1,k+1,n) - y   (i  ,j  ,k  ,n) ) + &
                    ( si   (2) + sj   (2) + sk   (2) ) * &
                    ( ydot(i+1,j+1,k+1,n) - ydot(i  ,j  ,k  ,n) )

                si(3) = ( (x(i,j+1,k+1,n) - x(i,j  ,k,n)) * &
                          (y(i,j  ,k+1,n) - y(i,j+1,k,n)) - &
                          (y(i,j+1,k+1,n) - y(i,j  ,k,n)) * &
                          (x(i,j  ,k+1,n) - x(i,j+1,k,n)) ) / 2
                sj(3) = ( (y(i+1,j,k+1,n) - y(i  ,j,k,n)) * &
                          (x(i  ,j,k+1,n) - x(i+1,j,k,n)) - &
                          (x(i+1,j,k+1,n) - x(i  ,j,k,n)) * &
                          (y(i  ,j,k+1,n) - y(i+1,j,k,n)) ) / 2
                sk(3) = ( (x(i+1,j+1,k,n) - x(i  ,j,k,n)) * &
                          (y(i  ,j+1,k,n) - y(i+1,j,k,n)) - &
                          (y(i+1,j+1,k,n) - y(i  ,j,k,n)) * &
                          (x(i  ,j+1,k,n) - x(i+1,j,k,n)) ) / 2
                sidot(3) = ( ((xdot(i,j+1,k+1,n) - xdot(i,j  ,k,n)) * &
                              (y(i,j  ,k+1,n) - y(i,j+1,k,n)) + &
                              (x(i,j+1,k+1,n) - x(i,j  ,k,n)) * &
                              (ydot(i,j  ,k+1,n) - ydot(i,j+1,k,n))) - &
                             ((ydot(i,j+1,k+1,n) - ydot(i,j  ,k,n)) * &
                              (x(i,j  ,k+1,n) - x(i,j+1,k,n)) + &
                              (y(i,j+1,k+1,n) - y(i,j  ,k,n)) * &
                              (xdot(i,j  ,k+1,n) - xdot(i,j+1,k,n))) )/2
                sjdot(3) = ( ((ydot(i+1,j,k+1,n) - ydot(i  ,j,k,n)) * &
                              (x(i  ,j,k+1,n) - x(i+1,j,k,n)) + &
                              (y(i+1,j,k+1,n) - y(i  ,j,k,n)) * &
                              (xdot(i  ,j,k+1,n) - xdot(i+1,j,k,n))) - &
                             ((xdot(i+1,j,k+1,n) - xdot(i  ,j,k,n)) * &
                              (y(i  ,j,k+1,n) - y(i+1,j,k,n)) + &
                              (x(i+1,j,k+1,n) - x(i  ,j,k,n)) * &
                              (ydot(i  ,j,k+1,n) - ydot(i+1,j,k,n))) )/2
                skdot(3) = ( ((xdot(i+1,j+1,k,n) - xdot(i  ,j,k,n)) * &
                              (y(i  ,j+1,k,n) - y(i+1,j,k,n)) + &
                              (x(i+1,j+1,k,n) - x(i  ,j,k,n)) * &
                              (ydot(i  ,j+1,k,n) - ydot(i+1,j,k,n))) - &
                             ((ydot(i+1,j+1,k,n) - ydot(i  ,j,k,n)) * &
                              (x(i  ,j+1,k,n) - x(i+1,j,k,n)) + &
                              (y(i+1,j+1,k,n) - y(i  ,j,k,n)) * &
                              (xdot(i  ,j+1,k,n) - xdot(i+1,j,k,n))) )/2
            dzdot = ( sidot(3) + sjdot(3) + skdot(3) ) * &
                    ( z   (i+1,j+1,k+1,n) - z   (i  ,j  ,k  ,n) ) + &
                    ( si   (3) + sj   (3) + sk   (3) ) * &
                    ( zdot(i+1,j+1,k+1,n) - zdot(i  ,j  ,k  ,n) )

                voldot = (dxdot + dydot + dzdot ) / 3

                if (abs(voldot).gt.1.d-14) then
!AJ                  write(*,*) 'voldot is nonzero at (itime, i, j, k): ',
!AJ     &                       itime,i,j,k,dxdot,dydot,dzdot
!AJ                  write(*,*) ' voldot = ', voldot
#if MPI
!AJ                  call abortmpi()
#else
                  stop
#endif
                end if
              end do
            end do
          end do
        end do
        write(*,*) 'voldot check is OK at itime ', itime
      end if

      if (1.lt.0) then
        do n=0,2*nharms*hbmove
          do k = 1,kmax-1
            do j = 1,jmax-1
              do i = 1,imax-1
                radl(1) = &
                  sqrt( y(i+1,j+1,k+1,n)**2 + z(i+1,j+1,k+1,n)**2 )
                radl(2) = &
                  sqrt( y(i+1,j  ,k+1,n)**2 + z(i+1,j  ,k+1,n)**2 )
                radl(3) = &
                  sqrt( y(i  ,j+1,k+1,n)**2 + z(i  ,j+1,k+1,n)**2 )
                radl(4) = &
                  sqrt( y(i  ,j  ,k+1,n)**2 + z(i  ,j  ,k+1,n)**2 )
                radl(5) = &
                  sqrt( y(i+1,j+1,k  ,n)**2 + z(i+1,j+1,k  ,n)**2 )
                radl(6) = &
                  sqrt( y(i+1,j  ,k  ,n)**2 + z(i+1,j  ,k  ,n)**2 )
                radl(7) = &
                  sqrt( y(i  ,j+1,k  ,n)**2 + z(i  ,j+1,k  ,n)**2 )
                radl(8) = &
                  sqrt( y(i  ,j  ,k  ,n)**2 + z(i  ,j  ,k  ,n)**2 )
                rad(i,j,k,n) = (radl(1)+radl(2)+radl(3)+radl(4)+ &
                                radl(5)+radl(6)+radl(7)+radl(8)) / 8
              end do
            end do
          end do
        end do
      end if

 100  format(1x,'Negative volumes. (iblock, i, j, k):',i5,1x,3i6)

      return
      end

!-----------------------------------------------------------------------
      subroutine plupitch(zeit,dtheta,dtheta0,omegas,xh,yh,dh0x,dh0y, &
        phipp,x0,y0,z0,x,y,z,xdot,ydot,zdot,imax,jmax,kmax)
!-----------------------------------------------------------------------
!     calculates time-dependent grid position and velocity when the whole
!     grid undergoes a harmonic rigid body plunging/pitching motion.
!-----------------------------------------------------------------------

      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax
      integer(kind=cosa_int) i,j,k
      real(kind=cosa_real) &
          x    (0:imax+1,0:jmax+1,0:kmax+1), &
          y    (0:imax+1,0:jmax+1,0:kmax+1), &
          z    (0:imax+1,0:jmax+1,0:kmax+1), &
          xdot (0:imax+1,0:jmax+1,0:kmax+1), &
          ydot (0:imax+1,0:jmax+1,0:kmax+1), &
          zdot (0:imax+1,0:jmax+1,0:kmax+1), &
          x0   (0:imax+1,0:jmax+1,0:kmax+1), &
          y0   (0:imax+1,0:jmax+1,0:kmax+1), &
          z0   (0:imax+1,0:jmax+1,0:kmax+1)
      real(kind=cosa_real) zeit,dtheta,dtheta0,omegas,xh,yh,dh0x,dh0y,phipp, &
          disp_tfac(2),vel_tfac(2),disph_tfac(2),velh_tfac(2)

      disph_tfac(1) = dh0x * dsin(omegas*zeit)
      disph_tfac(2) = dh0y * dsin(omegas*zeit)
      velh_tfac(1)  = dh0x * dcos(omegas*zeit)
      velh_tfac(2)  = dh0y * dcos(omegas*zeit)

      do k = 1,kmax
        do j = 1,jmax
          do i = 1,imax

            disp_tfac(1) = ( (x0(i,j,k) - xh) * (dcos(dtheta) - 1) - &
                             (y0(i,j,k) - yh) *  dsin(dtheta)      )
            disp_tfac(2) = ( (x0(i,j,k) - xh) *  dsin(dtheta)      + &
                             (y0(i,j,k) - yh) * (dcos(dtheta) - 1) )

            vel_tfac(1)  = (-(x0(i,j,k) - xh) * dsin(dtheta) - &
                             (y0(i,j,k) - yh) * dcos(dtheta) ) * &
                           dtheta0*dcos(omegas*zeit+phipp)
            vel_tfac(2)  = ( (x0(i,j,k) - xh) * dcos(dtheta) - &
                             (y0(i,j,k) - yh) * dsin(dtheta) ) * &
                           dtheta0*dcos(omegas*zeit+phipp)


            x   (i,j,k) = disp_tfac(1) + disph_tfac(1) + x0(i,j,k)
            y   (i,j,k) = disp_tfac(2) + disph_tfac(2) + y0(i,j,k)
            z   (i,j,k) = z0(i,j,k)
            xdot(i,j,k) = (vel_tfac(1) + velh_tfac(1)) * omegas
            ydot(i,j,k) = (vel_tfac(2) + velh_tfac(2)) * omegas
            zdot(i,j,k) = 0.d0

          end do
        end do
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine dist2wall(nl,x,y,z,xgwall,ygwall,zgwall,dist,work, &
                           bctopo)
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,iblock,imax,jmax,kmax,iblk,ixyz,iw,idist
      integer(kind=cosa_int) bctopo(*)
      real (kind=cosa_real) x(*),y(*),z(*),xgwall(ng_wall(nl)), &
           ygwall(ng_wall(nl)),zgwall(ng_wall(nl)),dist(*),work(*)

      do iblock = 1,mynblocks
        imax   = i_imax     (iblock,nl)
        jmax   = j_jmax     (iblock,nl)
        kmax   = k_kmax     (iblock,nl)
        idist  = 1 + off_p1 (iblock,nl)
        call init_dist(dist(idist),imax,jmax,kmax)
      end do

      do iw = 1,ng_wall(nl)
        do iblock = 1,mynblocks
          imax   = i_imax     (iblock,nl)
          jmax   = j_jmax     (iblock,nl)
          kmax   = k_kmax     (iblock,nl)
          ixyz   = 1 + off_p2 (iblock,nl)
          idist  = 1 + off_p1 (iblock,nl)
          call dist2wall_b(imax,jmax,kmax,x(ixyz),y(ixyz),z(ixyz), &
                           xgwall(iw),ygwall(iw),zgwall(iw),dist(idist))
        end do
      end do

      do iblock = 1,mynblocks
        imax   = i_imax     (iblock,nl)
        jmax   = j_jmax     (iblock,nl)
        kmax   = k_kmax     (iblock,nl)
        idist  = 1 + off_p1 (iblock,nl)
        call sqrt_dist(dist(idist),imax,jmax,kmax)
        if (debug) then
          ixyz   = 1 + off_p2 (iblock,nl)
          if (iblock.le.9) then
            write(filename,'(''dist_nl'',i1,''_blk0'',i1,''.dat'')') &
              nl,iblock
          else if (iblock.le.99) then
            write(filename,'(''dist_nl'',i1,''_blk'',i2,''.dat'')') &
              nl,iblock
          end if
          call wrdist(imax,jmax,kmax,x(ixyz),y(ixyz),z(ixyz), &
                      dist(idist),nl)
        end if
      end do

      call copy_array(1,1,1,1, 1,nl,work,dist, &
                      1,1,1,1, 1)

      do iblock = 1,mynblocks
        imax   = i_imax     (iblock,nl)
        jmax   = j_jmax     (iblock,nl)
        kmax   = k_kmax     (iblock,nl)
        idist  = 1 + off_p1 (iblock,nl)
        call dist_vert2cell(dist(idist),work(idist),imax,jmax,kmax,nl, &
                            iblock)
      end do

      call zero(1,1,1,1, 1,nl,work, &
                1,1,1,1, 1)

      return
      end

!-----------------------------------------------------------------------
      subroutine dist2wall_b(imax,jmax,kmax,x,y,z,xgwall,ygwall,zgwall, &
                             dist)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax
      real (kind=cosa_real) &
           x   (0:imax+1 ,0:jmax+1 ,0:kmax+1), &
           y   (0:imax+1 ,0:jmax+1 ,0:kmax+1), &
           z   (0:imax+1 ,0:jmax+1 ,0:kmax+1), &
           dist(0:imax   ,0:jmax   ,0:kmax  ), &
           xgwall, ygwall, zgwall

      call distance(xgwall,ygwall,zgwall,imax,jmax,kmax,x,y,z,dist)

      return
      end

!-----------------------------------------------------------------------
      subroutine init_dist(dist,imax,jmax,kmax)
!-----------------------------------------------------------------------

      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax
      integer(kind=cosa_int) i,j,k
      real (kind=cosa_real) fact
      real (kind=cosa_real) dist(0:imax,0:jmax,0:kmax)

      fact = 1.0d20

      do k = 1,kmax
        do j = 1,jmax
          do i = 1,imax
            dist(i,j,k) = fact
          end do
        end do
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine distance(xw,yw,zw,imax,jmax,kmax,x,y,z,dist)
!-----------------------------------------------------------------------

      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax
      integer(kind=cosa_int) i,j,k
      real (kind=cosa_real) xw,yw,zw,dx,dy,dz,d2
      real (kind=cosa_real) &
           x   (0:imax+1,0:jmax+1,0:kmax+1), &
           y   (0:imax+1,0:jmax+1,0:kmax+1), &
           z   (0:imax+1,0:jmax+1,0:kmax+1), &
           dist(0:imax  ,0:jmax  ,0:kmax  )

      do k=1,kmax
        do j=1,jmax
          do i=1,imax
            dx          = x(i,j,k)-xw
            dy          = y(i,j,k)-yw
            dz          = z(i,j,k)-zw
            d2          = dx*dx + dy*dy + dz*dz
            dist(i,j,k) = min(dist(i,j,k),d2)
          end do
        end do
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine sqrt_dist(dist,imax,jmax,kmax)
!-----------------------------------------------------------------------

      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax
      integer(kind=cosa_int) i,j,k
      real (kind=cosa_real) dist(0:imax,0:jmax,0:kmax)

      do k = 1,kmax
        do j = 1,jmax
          do i = 1,imax
            dist(i,j,k) = dsqrt(dist(i,j,k))
          end do
        end do
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine dist_vert2cell(dist,work,imax,jmax,kmax,nl,iblock)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax
      integer(kind=cosa_int) i,j,k,nl,iblock
      real (kind=cosa_real) &
        dist(0:imax,0:jmax,0:kmax), &
        work(0:imax,0:jmax,0:kmax)

!---- corners
      dist(0   ,0   ,0   ) = work(1   ,1   ,1   )
      dist(imax,0   ,0   ) = work(imax,1   ,1   )
      dist(0   ,jmax,0   ) = work(1   ,jmax,1   )
      dist(imax,jmax,0   ) = work(imax,jmax,1   )
      dist(0   ,0   ,kmax) = work(1   ,1   ,kmax)
      dist(imax,0   ,kmax) = work(imax,1   ,kmax)
      dist(0   ,jmax,kmax) = work(1   ,jmax,kmax)
      dist(imax,jmax,kmax) = work(imax,jmax,kmax)

!---- xi-junctions
      do i = 1,imax-1
        dist(i,0   ,0   ) = &
          (work(i+1,1   ,1   ) + work(i  ,1   ,1   )) /2
        dist(i,jmax,0   ) = &
          (work(i+1,jmax,1   ) + work(i  ,jmax,1   )) /2
        dist(i,0   ,kmax) = &
          (work(i+1,1   ,kmax) + work(i  ,1   ,kmax)) /2
        dist(i,jmax,kmax) = &
          (work(i+1,jmax,kmax) + work(i  ,jmax,kmax)) /2
      end do

!---- eta-junctions
      do j = 1,jmax-1
        dist(0   ,j,0   ) = &
          (work(1   ,j+1,1   ) + work(1   ,j  ,1   )) /2
        dist(imax,j,0   ) = &
          (work(imax,j+1,1   ) + work(imax,j  ,1   )) /2
        dist(0   ,j,kmax) = &
          (work(1   ,j+1,kmax) + work(1   ,j  ,kmax)) /2
        dist(imax,j,kmax) = &
          (work(imax,j+1,kmax) + work(imax,j  ,kmax)) /2
      end do

!---- zeta-junctions
      do k = 1,kmax-1
        dist(0   ,0   ,k) = &
          (work(1   ,1   ,k+1) + work(1   ,1   ,k  )) /2
        dist(imax,0   ,k) = &
          (work(imax,1   ,k+1) + work(imax,1   ,k  )) /2
        dist(0   ,jmax,k) = &
          (work(1   ,jmax,k+1) + work(1   ,jmax,k  )) /2
        dist(imax,jmax,k) = &
          (work(imax,jmax,k+1) + work(imax,jmax,k  )) /2
      end do

!---- xi-eta boundaries
      do j = 1,jmax-1
        do i = 1,imax-1
          dist(i,j,0   ) = (work(i  ,j,1   ) + work(i  ,j+1,1   ) + &
                            work(i+1,j,1   ) + work(i+1,j+1,1   )) / 4
          dist(i,j,kmax) = (work(i  ,j,kmax) + work(i  ,j+1,kmax) + &
                            work(i+1,j,kmax) + work(i+1,j+1,kmax)) / 4
        end do
      end do

!---- xi-zeta boundaries
      do k = 1,kmax-1
        do i = 1,imax-1
          dist(i,0   ,k) = (work(i,1   ,k  )    + work(i+1,1,k) + &
                            work(i,1   ,k+1)    + work(i+1,1,k+1)) / 4
          dist(i,jmax,k) = (work(i,jmax,k  ) + work(i+1,jmax,k  ) + &
                            work(i,jmax,k+1) + work(i+1,jmax,k+1)) / 4
        end do
      end do

!---- eta-zeta boundaries
      do k = 1,kmax-1
        do j = 1,jmax-1
          dist(0   ,j,k) = (work(1   ,j,k  ) + work(1   ,j+1,k  ) + &
                            work(1   ,j,k+1) + work(1   ,j+1,k+1)) / 4
          dist(imax,j,k) = (work(imax,j,k  ) + work(imax,j+1,k  ) + &
                            work(imax,j,k+1) + work(imax,j+1,k+1)) / 4
        end do
      end do

!---- internal nodes
      do k = 1,kmax-1
        do j = 1,jmax-1
          do i = 1,imax-1
            dist(i,j,k) = ( work(i  ,j  ,k  ) + work(i+1,j  ,k  ) + &
                            work(i+1,j+1,k  ) + work(i  ,j+1,k  ) + &
                            work(i  ,j  ,k+1) + work(i+1,j  ,k+1) + &
                            work(i+1,j+1,k+1) + work(i  ,j+1,k+1) ) / 8
          end do
        end do
      end do

      if (debug) then
        if (iblock.le.9) then
          write(filename,'(''celldist_nl'',i1,''_blk0'',i1,''.dat'')') &
            nl,iblock
        else if (iblock.le.99) then
          write(filename,'(''celldist_nl'',i1,''_blk'',i2,''.dat'')') &
            nl,iblock
        end if
        open(1,file=filename,status='replace')
        do k = 0,kmax
          do j = 0,jmax
            do i = 0,imax
              write(1,*) i,j,k,dist(i,j,k)
            end do
          end do
        end do
        close(unit=1)
      end if

      return
      end

!-----------------------------------------------------------------------
      subroutine facedist(nl,dist,fidist,fjdist,fkdist)
!-----------------------------------------------------------------------

      use cosa_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,iblock,imax,jmax,kmax,idist,ifdist
      real (kind=cosa_real) dist(*),fidist(*),fjdist(*),fkdist(*)

      do iblock = 1,mynblocks
        imax   = i_imax     (iblock,nl)
        jmax   = j_jmax     (iblock,nl)
        kmax   = k_kmax     (iblock,nl)
        idist  = 1 + off_p1 (iblock,nl)
        ifdist = 1 + off_0  (iblock,nl)
        call facedist_b(dist(idist),fidist(ifdist),fjdist(ifdist), &
                        fkdist(ifdist),imax,jmax,kmax)
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine facedist_b(dist,fidist,fjdist,fkdist,imax,jmax,kmax)
!-----------------------------------------------------------------------

      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax
      integer(kind=cosa_int) i,j,k
      real (kind=cosa_real) &
           dist  (0:imax,0:jmax,0:kmax), &
           fidist(  imax,  jmax,  kmax), &
           fjdist(  imax,  jmax,  kmax), &
           fkdist(  imax,  jmax,  kmax)

!---- distance of eta-zeta faces
      do k = 1,kmax-1
        do j = 1,jmax-1
          fidist(1   ,j,k) = dist(0   ,j,k)
          do i = 2,imax-1
            fidist(i,j,k) = (dist(i-1,j,k)+dist(i,j,k)) / 2
          end do
          fidist(imax,j,k) = dist(imax,j,k)
        end do
      end do

!---- distance of xi-zeta faces
      do k = 1,kmax-1
        do i = 1,imax-1
          fjdist(i,1   ,k) = dist(i,0   ,k)
          do j = 2,jmax-1
            fjdist(i,j,k) = (dist(i,j-1,k)+dist(i,j,k)) / 2
          end do
          fjdist(i,jmax,k) = dist(i,jmax,k)
        end do
      end do

!---- distance of xi-eta faces
      do j = 1,jmax-1
        do i = 1,imax-1
          fkdist(i,j,1   ) = dist(i,j,0   )
          do k = 2,kmax-1
            fkdist(i,j,k) = (dist(i,j,k-1)+dist(i,j,k)) / 2
          end do
          fkdist(i,j,kmax) = dist(i,j,kmax)
        end do
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine extract_wall(nl,x,y,z,xwall,ywall,zwall,bctopo)
!-----------------------------------------------------------------------
!---- msc, 03. March 2012: when COSA will use deforming grids with the
!          K-\omega SST model, wall distances will have to be computed
!          at each physical time.
!          To date, the wall distances are computed only once, using the
!          steady grid coordinates. Both COSA and POSCOSA call present 
!          routine extract_wall soon after reading the steady coordinates
!          from mesh.dat. The parameter dim5h is not included in the 
!          definition of ixy of this routine (see below) for consistency
!          with the construction of the x and y arrays in readgrid.
!          When deforming grid problems will be used with K-\omega SST
!          model, however, distances will have to be computed for each
!          physical time, and construction of xgwall and ygwall will 
!          require some code changes.
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,iblock,iblk,imax,jmax,kmax,ixyz,ixyzw,nw
      integer(kind=cosa_int) bctopo(*)
      real (kind=cosa_real) x(*),y(*),z(*),xwall(*),ywall(*),zwall(*)

      do iblock = 1,mynblocks
        imax   = i_imax     (iblock,nl)
        jmax   = j_jmax     (iblock,nl)
        kmax   = k_kmax     (iblock,nl)
        iblk   = 1 + off_bct(iblock,nl)
        ixyz   = 1 + off_p2 (iblock,nl)
        nw     = n_wall     (iblock,nl)
        ixyzw  = 1 + off_1dw(iblock,nl)
!dbg        write(*,*) 'iblock,nw,ixyw,ng_wall(1)',iblock,nw,ixyzw,ng_wall(1)
        call extract_bwall(imax,jmax,kmax,x(ixyz),y(ixyz),z(ixyz), &
                           nbcs(iblock),bctopo(iblk),nw,xwall(ixyzw), &
                           ywall(ixyzw),zwall(ixyzw))
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine extract_bwall(imax,jmax,kmax,x,y,z,nbcs,bctopo,nw, &
                               xwall,ywall,zwall)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,nbcs,nw
      integer(kind=cosa_int) bctopo(10,nbcs)
      real (kind=cosa_real) &
           x   (0:imax+1 ,0:jmax+1 ,0:kmax+1), &
           y   (0:imax+1 ,0:jmax+1 ,0:kmax+1), &
           z   (0:imax+1 ,0:jmax+1 ,0:kmax+1), &
           xwall(nw), ywall(nw), zwall(nw)

      call extract_bhwall(imax,jmax,kmax,x,y,z,nbcs,bctopo,nw,xwall, &
                          ywall,zwall)

      return
      end

!-----------------------------------------------------------------------
      subroutine extract_bhwall(imax,jmax,kmax,x,y,z,nbcs,bctopo,nw, &
                          xwall,ywall,zwall)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,nbcs,nw
      integer(kind=cosa_int) ibc,ibctyp,i2,i3,ic1,ic2,ic3,ijkmax(3),idir,istrt(3), &
                iend(3),bctyp,inrout,ibcpt,ibcw,jbcw,kbcw,iw
      integer(kind=cosa_int) bctopo(10,nbcs)
      real (kind=cosa_real) &
           x   (0:imax+1 ,0:jmax+1 ,0:kmax+1), &
           y   (0:imax+1 ,0:jmax+1 ,0:kmax+1), &
           z   (0:imax+1 ,0:jmax+1 ,0:kmax+1), &
           xwall(nw), ywall(nw), zwall(nw)

      ijkmax(1) = imax
      ijkmax(2) = jmax
      ijkmax(3) = kmax

      iw = 0

      do ibc = 1,nbcs

        ibctyp = bctopo(1,ibc)

        if (ibctyp/100.eq.15) then

          idir     = bctopo(2,ibc)
          inrout   = bctopo(3,ibc)
          istrt(1) = bctopo(4,ibc)
          iend (1) = bctopo(5,ibc)
          istrt(2) = bctopo(6,ibc)
          iend (2) = bctopo(7,ibc)
          istrt(3) = bctopo(8,ibc)
          iend (3) = bctopo(9,ibc)

          do i2 = 1,3
            if (i2.ne.idir) then
              iend(i2) = iend(i2) + 1
            end if
          end do

!         Set needed variables depending on whether the boundary is
!         the inner boundary (INROUT = 1) or 
!         the outer boundary (INROUT > 1).
!         IBCPT  is the grid boundary location

          if (inrout.eq.1) then
            ibcpt  = 1
          else
            ibcpt  = ijkmax(idir)
 
         end if
!
          ic1 = cyc (idir, 1)
          ic2 = cyc (idir, 2)
          ic3 = cyc (idir, 3)

!         Loop over all grid points on the current block boundary
!
          do i3 = istrt(ic3),iend(ic3)
            do i2 = istrt(ic2),iend(ic2)

              ibcw = ibcpt *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
              jbcw = ibcpt *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
              kbcw = ibcpt *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

              iw        = iw + 1
              xwall(iw) = x(ibcw,jbcw,kbcw)
              ywall(iw) = y(ibcw,jbcw,kbcw)
              zwall(iw) = z(ibcw,jbcw,kbcw)

            end do
          end do

        end if

      end do

      if (iw.ne.nw) then
        write(*,*) 'mismatch in extraction of wall node coordinates. Abo &
     &rting!'
        stop
      end if

      return
      end

!-----------------------------------------------------------------------
      subroutine merge_wall(nl,xwall,ywall,zwall,xgwall,ygwall,zgwall)
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,iblock,ixyz,ixyzw,nw
      real (kind=cosa_real) xwall(*),ywall(*),zwall(*),xgwall(ng_wall(nl)), &
           ygwall(ng_wall(nl)),zgwall(ng_wall(nl))

      logical parallel                                      
                      
      if(debug) then
        write(*,*) 'merge_wall'                               
      end if
                                             
      call usingparallel(parallel)                          
                                                            
      if (parallel) then                                     
                                                            
        call parallelmergewall(nl,xwall,ywall,zwall,xgwall,ygwall, &
                               zgwall)
                                                            
      else                                                  
                                                            
        do iblock = 1,mynblocks                            
           nw     = n_wall     (iblock,nl)                 
           ixyzw  = 1 + off_1dw(iblock,nl)                 
           if(debug) then
              write(*,*) 'iblock1,nw,ixyzw1,ng_wall1(1)', &
                   iblock,nw,ixyzw,ng_wall(nl)
           end if
           call merge_bwall(nw,xwall(ixyzw),ywall(ixyzw),zwall(ixyzw), &
                            ixyzw,ng_wall(nl),xgwall,ygwall,zgwall)
        end do                                             
                                                            
      end if                                                

      return
      end

!-----------------------------------------------------------------------
      subroutine merge_bwall(nw,xwall,ywall,zwall,ixyzw,ng_wall,xgwall, &
                             ygwall,zgwall)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nw,ixyzw,ng_wall
      integer(kind=cosa_int) iw
      real (kind=cosa_real) &
           xwall(nw),       ywall(nw),       zwall(nw), &
           xgwall(ng_wall), ygwall(ng_wall), zgwall(ng_wall)

      do iw = 1,nw
        xgwall(ixyzw-1+iw) = xwall(iw)
        ygwall(ixyzw-1+iw) = ywall(iw)
        zgwall(ixyzw-1+iw) = zwall(iw)
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine vmetrix(nl,x,y,z,xideri,xiderj,xiderk,etaderi,etaderj, &
        etaderk,zetaderi,zetaderj,zetaderk,vol)
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,iblock,imax,jmax,kmax,ivol,ixyz,ivmt
      real (kind=cosa_real) x(*),y(*),z(*),xideri(*), &
           xiderj(*),xiderk(*),etaderi(*),etaderj(*),etaderk(*), &
           zetaderi(*),zetaderj(*),zetaderk(*),vol(*)

      do iblock = 1,mynblocks
        imax  = i_imax  (iblock,nl)
        jmax  = j_jmax  (iblock,nl)
        kmax  = k_kmax  (iblock,nl)
        ixyz  = 1 + off_p2(iblock,nl)         * dim5h
        ivol  = 1 + off_p1(iblock,nl)
        ivmt  = 1 + off_0  (iblock,nl) * 3    * dim5h
        call vmetrix_b(x(ixyz),y(ixyz),z(ixyz),xideri(ivmt), &
          xiderj(ivmt),xiderk(ivmt),etaderi(ivmt),etaderj(ivmt), &
          etaderk(ivmt),zetaderi(ivmt),zetaderj(ivmt),zetaderk(ivmt), &
          vol(ivol),iblock,imax,jmax,kmax,nharms)
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine vmetrix_b(x,y,z,xideri,xiderj,xiderk,etaderi,etaderj, &
        etaderk,zetaderi,zetaderj,zetaderk,vol,iblock,imax,jmax,kmax, &
        nharms)
!-----------------------------------------------------------------------

!     xix   =  (z_zeta * y_eta  - z_eta  * y_zeta ) / inv(|jac|)
!     xiy   =  (x_zeta * z_eta  - x_eta  * z_zeta ) / inv(|jac|)
!     xiz   =  (y_zeta * x_eta  - y_eta  * x_zeta ) / inv(|jac|)
!     etax  =  (z_xi   * y_zeta - z_zeta * y_xi   ) / inv(|jac|)
!     etay  =  (x_xi   * z_zeta - x_zeta * z_xi   ) / inv(|jac|)
!     etaz  =  (y_xi   * x_zeta - y_zeta * x_xi   ) / inv(|jac|)
!     zetax =  (z_eta  * y_xi   - z_xi   * y_eta  ) / inv(|jac|)
!     zetay =  (x_eta  * z_xi   - x_xi   * z_eta  ) / inv(|jac|)
!     zetaz =  (y_eta  * x_xi   - y_xi   * x_eta  ) / inv(|jac|)

!     xideri  (1,i,j,n) = xix   on i-faces
!     xideri  (2,i,j,n) = xiy   on i-faces
!     xideri  (3,i,j,n) = xiz   on i-faces
!     etaderi (1,i,j,n) = etax  on i-faces
!     etaderi (2,i,j,n) = etay  on i-faces
!     etaderi (3,i,j,n) = etaz  on i-faces
!     zetaderi(1,i,j,n) = zetax on i-faces
!     zetaderi(2,i,j,n) = zetay on i-faces
!     zetaderi(3,i,j,n) = zetaz on i-faces

!     xiderj,etaderj,zetaderj are the corresponding values for j-faces
!     xiderk,etaderk,zetaderk are the corresponding values for k-faces
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,nharms
      integer(kind=cosa_int) iblock,i,j,k,n,i1
      real(kind=cosa_real) &
          x       (     0:imax+1,0:jmax+1,0:kmax+1,0:2*nharms*hbmove), &
          y       (     0:imax+1,0:jmax+1,0:kmax+1,0:2*nharms*hbmove), &
          z       (     0:imax+1,0:jmax+1,0:kmax+1,0:2*nharms*hbmove), &
          xideri  (3   ,imax    ,jmax    ,kmax    ,0:2*nharms*hbmove), &
          xiderj  (3   ,imax    ,jmax    ,kmax    ,0:2*nharms*hbmove), &
          xiderk  (3   ,imax    ,jmax    ,kmax    ,0:2*nharms*hbmove), &
          etaderi (3   ,imax    ,jmax    ,kmax    ,0:2*nharms*hbmove), &
          etaderj (3   ,imax    ,jmax    ,kmax    ,0:2*nharms*hbmove), &
          etaderk (3   ,imax    ,jmax    ,kmax    ,0:2*nharms*hbmove), &
          zetaderi(3   ,imax    ,jmax    ,kmax    ,0:2*nharms*hbmove), &
          zetaderj(3   ,imax    ,jmax    ,kmax    ,0:2*nharms*hbmove), &
          zetaderk(3   ,imax    ,jmax    ,kmax    ,0:2*nharms*hbmove), &
          vol     (0:imax,0:jmax,0:kmax)
      real(kind=cosa_real) xix,xiy,xiz,etax,etay,etaz,zetax,zetay,zetaz,volvisc, &
        x_xi(0:4),y_xi(0:4),z_xi(0:4),x_eta(0:4),y_eta(0:4),z_eta(0:4), &
       x_zeta(0:4),y_zeta(0:4),z_zeta(0:4),dist(4)

      do n=0,2*nharms*hbmove

!------ metrix on xi-faces
        do k = 1,kmax-1
          do j = 1,jmax-1
            do i = 1,imax

              x_xi(1)   = ( x(i+1,j  ,k  ,n) - x(i-1,j  ,k  ,n) ) / 2
              y_xi(1)   = ( y(i+1,j  ,k  ,n) - y(i-1,j  ,k  ,n) ) / 2
              z_xi(1)   = ( z(i+1,j  ,k  ,n) - z(i-1,j  ,k  ,n) ) / 2
              x_xi(2)   = ( x(i+1,j+1,k  ,n) - x(i-1,j+1,k  ,n) ) / 2
              y_xi(2)   = ( y(i+1,j+1,k  ,n) - y(i-1,j+1,k  ,n) ) / 2
              z_xi(2)   = ( z(i+1,j+1,k  ,n) - z(i-1,j+1,k  ,n) ) / 2
              x_xi(3)   = ( x(i+1,j+1,k+1,n) - x(i-1,j+1,k+1,n) ) / 2
              y_xi(3)   = ( y(i+1,j+1,k+1,n) - y(i-1,j+1,k+1,n) ) / 2
              z_xi(3)   = ( z(i+1,j+1,k+1,n) - z(i-1,j+1,k+1,n) ) / 2
              x_xi(4)   = ( x(i+1,j  ,k+1,n) - x(i-1,j  ,k+1,n) ) / 2
              y_xi(4)   = ( y(i+1,j  ,k+1,n) - y(i-1,j  ,k+1,n) ) / 2
              z_xi(4)   = ( z(i+1,j  ,k+1,n) - z(i-1,j  ,k+1,n) ) / 2

              x_eta(1)  = ( x(i  ,j+1,k  ,n) - x(i  ,j  ,k  ,n) )
              y_eta(1)  = ( y(i  ,j+1,k  ,n) - y(i  ,j  ,k  ,n) )
              z_eta(1)  = ( z(i  ,j+1,k  ,n) - z(i  ,j  ,k  ,n) )
              x_eta(2)  = ( x(i  ,j+1,k+1,n) - x(i  ,j  ,k+1,n) )
              y_eta(2)  = ( y(i  ,j+1,k+1,n) - y(i  ,j  ,k+1,n) )
              z_eta(2)  = ( z(i  ,j+1,k+1,n) - z(i  ,j  ,k+1,n) )

              x_zeta(1) = ( x(i  ,j  ,k+1,n) - x(i  ,j  ,k  ,n) )
              y_zeta(1) = ( y(i  ,j  ,k+1,n) - y(i  ,j  ,k  ,n) )
              z_zeta(1) = ( z(i  ,j  ,k+1,n) - z(i  ,j  ,k  ,n) )
              x_zeta(2) = ( x(i  ,j+1,k+1,n) - x(i  ,j+1,k  ,n) )
              y_zeta(2) = ( y(i  ,j+1,k+1,n) - y(i  ,j+1,k  ,n) )
              z_zeta(2) = ( z(i  ,j+1,k+1,n) - z(i  ,j+1,k  ,n) )

              volvisc   = (vol(i-1,j,k)+vol(i,j,k)) / 2

              if (i.eq.1) then
                dist(1)=sqrt( (x(i-1,j  ,k  ,n) - x(i  ,j  ,k  ,n))**2 + &
                              (y(i-1,j  ,k  ,n) - y(i  ,j  ,k  ,n))**2 + &
                              (z(i-1,j  ,k  ,n) - z(i  ,j  ,k  ,n))**2 )
                dist(2)=sqrt( (x(i-1,j+1,k  ,n) - x(i  ,j+1,k  ,n))**2 + &
                              (y(i-1,j+1,k  ,n) - y(i  ,j+1,k  ,n))**2 + &
                              (z(i-1,j+1,k  ,n) - z(i  ,j+1,k  ,n))**2 )
                dist(3)=sqrt( (x(i-1,j+1,k+1,n) - x(i  ,j+1,k+1,n))**2 + &
                              (y(i-1,j+1,k+1,n) - y(i  ,j+1,k+1,n))**2 + &
                              (z(i-1,j+1,k+1,n) - z(i  ,j+1,k+1,n))**2 )
                dist(4)=sqrt( (x(i-1,j  ,k+1,n) - x(i  ,j  ,k+1,n))**2 + &
                              (y(i-1,j  ,k+1,n) - y(i  ,j  ,k+1,n))**2 + &
                              (z(i-1,j  ,k+1,n) - z(i  ,j  ,k+1,n))**2 )
                do i1=1,4
                  if (dist(i1).lt.1.d-14) then
                    x_xi(i1) = 2 * x_xi(i1)
                    y_xi(i1) = 2 * y_xi(i1)
                    z_xi(i1) = 2 * z_xi(i1)
                  end if
                end do
              else if (i.eq.imax) then
                dist(1)=sqrt( (x(i+1,j  ,k  ,n) - x(i  ,j  ,k  ,n))**2 + &
                              (y(i+1,j  ,k  ,n) - y(i  ,j  ,k  ,n))**2 + &
                              (z(i+1,j  ,k  ,n) - z(i  ,j  ,k  ,n))**2 )
                dist(2)=sqrt( (x(i+1,j+1,k  ,n) - x(i  ,j+1,k  ,n))**2 + &
                              (y(i+1,j+1,k  ,n) - y(i  ,j+1,k  ,n))**2 + &
                              (z(i+1,j+1,k  ,n) - z(i  ,j+1,k  ,n))**2 )
                dist(3)=sqrt( (x(i+1,j+1,k+1,n) - x(i  ,j+1,k+1,n))**2 + &
                              (y(i+1,j+1,k+1,n) - y(i  ,j+1,k+1,n))**2 + &
                              (z(i+1,j+1,k+1,n) - z(i  ,j+1,k+1,n))**2 )
                dist(4)=sqrt( (x(i+1,j  ,k+1,n) - x(i  ,j  ,k+1,n))**2 + &
                              (y(i+1,j  ,k+1,n) - y(i  ,j  ,k+1,n))**2 + &
                              (z(i+1,j  ,k+1,n) - z(i  ,j  ,k+1,n))**2 )
                do i1=1,4
                  if (dist(i1).lt.1.d-14) then
                    x_xi(i1) = 2 * x_xi(i1)
                    y_xi(i1) = 2 * y_xi(i1)
                    z_xi(i1) = 2 * z_xi(i1)
                  end if
                end do
              end if

              x_xi(0)   = &
                ( x_xi(1)   + x_xi(2)   + x_xi(3)   + x_xi(4)   ) / 4
              y_xi(0)   = &
                ( y_xi(1)   + y_xi(2)   + y_xi(3)   + y_xi(4)   ) / 4
              z_xi(0)   = &
                ( z_xi(1)   + z_xi(2)   + z_xi(3)   + z_xi(4)   ) / 4

              x_eta(0)  = ( x_eta(1)  + x_eta(2)  ) / 2
              y_eta(0)  = ( y_eta(1)  + y_eta(2)  ) / 2
              z_eta(0)  = ( z_eta(1)  + z_eta(2)  ) / 2

              x_zeta(0) = ( x_zeta(1) + x_zeta(2) ) / 2
              y_zeta(0) = ( y_zeta(1) + y_zeta(2) ) / 2
              z_zeta(0) = ( z_zeta(1) + z_zeta(2) ) / 2

              xix   =  (z_zeta(0) * y_eta(0)  - z_eta(0)  * y_zeta(0) )/ &
                       volvisc
              xiy   =  (x_zeta(0) * z_eta(0)  - x_eta(0)  * z_zeta(0) )/ &
                       volvisc
              xiz   =  (y_zeta(0) * x_eta(0)  - y_eta(0)  * x_zeta(0) )/ &
                       volvisc
              etax  =  (z_xi(0)   * y_zeta(0) - z_zeta(0) * y_xi(0)   )/ &
                       volvisc
              etay  =  (x_xi(0)   * z_zeta(0) - x_zeta(0) * z_xi(0)   )/ &
                       volvisc
              etaz  =  (y_xi(0)   * x_zeta(0) - y_zeta(0) * x_xi(0)   )/ &
                       volvisc
              zetax =  (z_eta(0)  * y_xi(0)   - z_xi(0)   * y_eta(0)  )/ &
                       volvisc
              zetay =  (x_eta(0)  * z_xi(0)   - x_xi(0)   * z_eta(0)  )/ &
                       volvisc
              zetaz =  (y_eta(0)  * x_xi(0)   - y_xi(0)   * x_eta(0)  )/ &
                       volvisc

              xideri  (1,i,j,k,n) = xix
              xideri  (2,i,j,k,n) = xiy
              xideri  (3,i,j,k,n) = xiz
              etaderi (1,i,j,k,n) = etax
              etaderi (2,i,j,k,n) = etay
              etaderi (3,i,j,k,n) = etaz
              zetaderi(1,i,j,k,n) = zetax
              zetaderi(2,i,j,k,n) = zetay
              zetaderi(3,i,j,k,n) = zetaz

            end do
          end do
        end do

!------ metrix on eta-faces
        do k = 1,kmax-1
          do j = 1,jmax
            do i = 1,imax-1

              x_xi(1)   = ( x(i+1,j  ,k  ,n) - x(i  ,j  ,k  ,n) )
              y_xi(1)   = ( y(i+1,j  ,k  ,n) - y(i  ,j  ,k  ,n) )
              z_xi(1)   = ( z(i+1,j  ,k  ,n) - z(i  ,j  ,k  ,n) )
              x_xi(2)   = ( x(i+1,j  ,k+1,n) - x(i  ,j  ,k+1,n) )
              y_xi(2)   = ( y(i+1,j  ,k+1,n) - y(i  ,j  ,k+1,n) )
              z_xi(2)   = ( z(i+1,j  ,k+1,n) - z(i  ,j  ,k+1,n) )

              x_eta(1)  = ( x(i  ,j+1,k  ,n) - x(i  ,j-1,k  ,n) ) / 2
              y_eta(1)  = ( y(i  ,j+1,k  ,n) - y(i  ,j-1,k  ,n) ) / 2
              z_eta(1)  = ( z(i  ,j+1,k  ,n) - z(i  ,j-1,k  ,n) ) / 2
              x_eta(2)  = ( x(i  ,j+1,k+1,n) - x(i  ,j-1,k+1,n) ) / 2
              y_eta(2)  = ( y(i  ,j+1,k+1,n) - y(i  ,j-1,k+1,n) ) / 2
              z_eta(2)  = ( z(i  ,j+1,k+1,n) - z(i  ,j-1,k+1,n) ) / 2
              x_eta(3)  = ( x(i+1,j+1,k+1,n) - x(i+1,j-1,k+1,n) ) / 2
              y_eta(3)  = ( y(i+1,j+1,k+1,n) - y(i+1,j-1,k+1,n) ) / 2
              z_eta(3)  = ( z(i+1,j+1,k+1,n) - z(i+1,j-1,k+1,n) ) / 2
              x_eta(4)  = ( x(i+1,j+1,k  ,n) - x(i+1,j-1,k  ,n) ) / 2
              y_eta(4)  = ( y(i+1,j+1,k  ,n) - y(i+1,j-1,k  ,n) ) / 2
              z_eta(4)  = ( z(i+1,j+1,k  ,n) - z(i+1,j-1,k  ,n) ) / 2

              x_zeta(1) = ( x(i  ,j  ,k+1,n) - x(i  ,j  ,k  ,n) )
              y_zeta(1) = ( y(i  ,j  ,k+1,n) - y(i  ,j  ,k  ,n) )
              z_zeta(1) = ( z(i  ,j  ,k+1,n) - z(i  ,j  ,k  ,n) )
              x_zeta(2) = ( x(i+1,j  ,k+1,n) - x(i+1,j  ,k  ,n) )
              y_zeta(2) = ( y(i+1,j  ,k+1,n) - y(i+1,j  ,k  ,n) )
              z_zeta(2) = ( z(i+1,j  ,k+1,n) - z(i+1,j  ,k  ,n) )

              volvisc   = (vol(i,j,k)+vol(i,j-1,k)) / 2

              if (j.eq.1) then
                dist(1)=sqrt( (x(i  ,j-1,k  ,n) - x(i  ,j  ,k  ,n))**2 + &
                              (y(i  ,j-1,k  ,n) - y(i  ,j  ,k  ,n))**2 + &
                              (z(i  ,j-1,k  ,n) - z(i  ,j  ,k  ,n))**2 )
                dist(2)=sqrt( (x(i  ,j-1,k+1,n) - x(i  ,j  ,k+1,n))**2 + &
                              (y(i  ,j-1,k+1,n) - y(i  ,j  ,k+1,n))**2 + &
                              (z(i  ,j-1,k+1,n) - z(i  ,j  ,k+1,n))**2 )
                dist(3)=sqrt( (x(i+1,j-1,k+1,n) - x(i+1,j  ,k+1,n))**2 + &
                              (y(i+1,j-1,k+1,n) - y(i+1,j  ,k+1,n))**2 + &
                              (z(i+1,j-1,k+1,n) - z(i+1,j  ,k+1,n))**2 )
                dist(4)=sqrt( (x(i+1,j-1,k  ,n) - x(i+1,j  ,k  ,n))**2 + &
                              (y(i+1,j-1,k  ,n) - y(i+1,j  ,k  ,n))**2 + &
                              (z(i+1,j-1,k  ,n) - z(i+1,j  ,k  ,n))**2 )
                do i1=1,4
                  if (dist(i1).lt.1.d-14) then
                    x_eta(i1) = 2 * x_eta(i1)
                    y_eta(i1) = 2 * y_eta(i1)
                    z_eta(i1) = 2 * z_eta(i1)
                  end if
                end do
              else if (j.eq.jmax) then
                dist(1)=sqrt( (x(i  ,j+1,k  ,n) - x(i  ,j  ,k  ,n))**2 + &
                              (y(i  ,j+1,k  ,n) - y(i  ,j  ,k  ,n))**2 + &
                              (z(i  ,j+1,k  ,n) - z(i  ,j  ,k  ,n))**2 )
                dist(2)=sqrt( (x(i  ,j+1,k+1,n) - x(i  ,j  ,k+1,n))**2 + &
                              (y(i  ,j+1,k+1,n) - y(i  ,j  ,k+1,n))**2 + &
                              (z(i  ,j+1,k+1,n) - z(i  ,j  ,k+1,n))**2 )
                dist(3)=sqrt( (x(i+1,j+1,k+1,n) - x(i+1,j  ,k+1,n))**2 + &
                              (y(i+1,j+1,k+1,n) - y(i+1,j  ,k+1,n))**2 + &
                              (z(i+1,j+1,k+1,n) - z(i+1,j  ,k+1,n))**2 )
                dist(4)=sqrt( (x(i+1,j+1,k  ,n) - x(i+1,j  ,k  ,n))**2 + &
                              (y(i+1,j+1,k  ,n) - y(i+1,j  ,k  ,n))**2 + &
                              (z(i+1,j+1,k  ,n) - z(i+1,j  ,k  ,n))**2 )
                do i1=1,4
                  if (dist(i1).lt.1.d-14) then
                    x_eta(i1) = 2 * x_eta(i1)
                    y_eta(i1) = 2 * y_eta(i1)
                    z_eta(i1) = 2 * z_eta(i1)
                  end if
                end do
              end if

              x_xi(0)   = ( x_xi(1)   + x_xi(2)   ) / 2
              y_xi(0)   = ( y_xi(1)   + y_xi(2)   ) / 2
              z_xi(0)   = ( z_xi(1)   + z_xi(2)   ) / 2

              x_eta(0)  = &
                ( x_eta(1)  + x_eta(2)  + x_eta(3)  + x_eta(4)  ) / 4
              y_eta(0)  = &
                ( y_eta(1)  + y_eta(2)  + y_eta(3)  + y_eta(4)  ) / 4
              z_eta(0)  = &
                ( z_eta(1)  + z_eta(2)  + z_eta(3)  + z_eta(4)  ) / 4

              x_zeta(0) = ( x_zeta(1) + x_zeta(2) ) / 2
              y_zeta(0) = ( y_zeta(1) + y_zeta(2) ) / 2
              z_zeta(0) = ( z_zeta(1) + z_zeta(2) ) / 2

              xix   =  (z_zeta(0) * y_eta(0)  - z_eta(0)  * y_zeta(0) )/ &
                       volvisc
              xiy   =  (x_zeta(0) * z_eta(0)  - x_eta(0)  * z_zeta(0) )/ &
                       volvisc
              xiz   =  (y_zeta(0) * x_eta(0)  - y_eta(0)  * x_zeta(0) )/ &
                       volvisc
              etax  =  (z_xi(0)   * y_zeta(0) - z_zeta(0) * y_xi(0)   )/ &
                       volvisc
              etay  =  (x_xi(0)   * z_zeta(0) - x_zeta(0) * z_xi(0)   )/ &
                       volvisc
              etaz  =  (y_xi(0)   * x_zeta(0) - y_zeta(0) * x_xi(0)   )/ &
                       volvisc
              zetax =  (z_eta(0)  * y_xi(0)   - z_xi(0)   * y_eta(0)  )/ &
                       volvisc
              zetay =  (x_eta(0)  * z_xi(0)   - x_xi(0)   * z_eta(0)  )/ &
                       volvisc
              zetaz =  (y_eta(0)  * x_xi(0)   - y_xi(0)   * x_eta(0)  )/ &
                       volvisc

              xiderj  (1,i,j,k,n) = xix
              xiderj  (2,i,j,k,n) = xiy
              xiderj  (3,i,j,k,n) = xiz
              etaderj (1,i,j,k,n) = etax
              etaderj (2,i,j,k,n) = etay
              etaderj (3,i,j,k,n) = etaz
              zetaderj(1,i,j,k,n) = zetax
              zetaderj(2,i,j,k,n) = zetay
              zetaderj(3,i,j,k,n) = zetaz

            end do
          end do
        end do

!------ metrix on zeta-faces
        do k = 1,kmax
          do j = 1,jmax-1
            do i = 1,imax-1

              x_xi(1)   = ( x(i+1,j  ,k  ,n) - x(i  ,j  ,k  ,n) )
              y_xi(1)   = ( y(i+1,j  ,k  ,n) - y(i  ,j  ,k  ,n) )
              z_xi(1)   = ( z(i+1,j  ,k  ,n) - z(i  ,j  ,k  ,n) )
              x_xi(2)   = ( x(i+1,j+1,k  ,n) - x(i  ,j+1,k  ,n) )
              y_xi(2)   = ( y(i+1,j+1,k  ,n) - y(i  ,j+1,k  ,n) )
              z_xi(2)   = ( z(i+1,j+1,k  ,n) - z(i  ,j+1,k  ,n) )

              x_eta(1)  = ( x(i  ,j+1,k  ,n) - x(i  ,j  ,k  ,n) )
              y_eta(1)  = ( y(i  ,j+1,k  ,n) - y(i  ,j  ,k  ,n) )
              z_eta(1)  = ( z(i  ,j+1,k  ,n) - z(i  ,j  ,k  ,n) )
              x_eta(2)  = ( x(i+1,j+1,k  ,n) - x(i+1,j  ,k  ,n) )
              y_eta(2)  = ( y(i+1,j+1,k  ,n) - y(i+1,j  ,k  ,n) )
              z_eta(2)  = ( z(i+1,j+1,k  ,n) - z(i+1,j  ,k  ,n) )

              x_zeta(1) = ( x(i  ,j  ,k+1,n) - x(i  ,j  ,k-1,n) ) / 2
              y_zeta(1) = ( y(i  ,j  ,k+1,n) - y(i  ,j  ,k-1,n) ) / 2
              z_zeta(1) = ( z(i  ,j  ,k+1,n) - z(i  ,j  ,k-1,n) ) / 2
              x_zeta(2) = ( x(i+1,j  ,k+1,n) - x(i+1,j  ,k-1,n) ) / 2
              y_zeta(2) = ( y(i+1,j  ,k+1,n) - y(i+1,j  ,k-1,n) ) / 2
              z_zeta(2) = ( z(i+1,j  ,k+1,n) - z(i+1,j  ,k-1,n) ) / 2
              x_zeta(3) = ( x(i+1,j+1,k+1,n) - x(i+1,j+1,k-1,n) ) / 2
              y_zeta(3) = ( y(i+1,j+1,k+1,n) - y(i+1,j+1,k-1,n) ) / 2
              z_zeta(3) = ( z(i+1,j+1,k+1,n) - z(i+1,j+1,k-1,n) ) / 2
              x_zeta(4) = ( x(i  ,j+1,k+1,n) - x(i  ,j+1,k-1,n) ) / 2
              y_zeta(4) = ( y(i  ,j+1,k+1,n) - y(i  ,j+1,k-1,n) ) / 2
              z_zeta(4) = ( z(i  ,j+1,k+1,n) - z(i  ,j+1,k-1,n) ) / 2

              volvisc   = (vol(i,j,k)+vol(i,j,k-1)) / 2

              if (k.eq.1) then
                dist(1)=sqrt( (x(i  ,j  ,k-1,n) - x(i  ,j  ,k  ,n))**2 + &
                              (y(i  ,j  ,k-1,n) - y(i  ,j  ,k  ,n))**2 + &
                              (z(i  ,j  ,k-1,n) - z(i  ,j  ,k  ,n))**2 )
                dist(2)=sqrt( (x(i+1,j  ,k-1,n) - x(i+1,j  ,k  ,n))**2 + &
                              (y(i+1,j  ,k-1,n) - y(i+1,j  ,k  ,n))**2 + &
                              (z(i+1,j  ,k-1,n) - z(i+1,j  ,k  ,n))**2 )
                dist(3)=sqrt( (x(i+1,j+1,k-1,n) - x(i+1,j+1,k  ,n))**2 + &
                              (y(i+1,j+1,k-1,n) - y(i+1,j+1,k  ,n))**2 + &
                              (z(i+1,j+1,k-1,n) - z(i+1,j+1,k  ,n))**2 )
                dist(4)=sqrt( (x(i  ,j+1,k-1,n) - x(i  ,j+1,k  ,n))**2 + &
                              (y(i  ,j+1,k-1,n) - y(i  ,j+1,k  ,n))**2 + &
                              (z(i  ,j+1,k-1,n) - z(i  ,j+1,k  ,n))**2 )
                do i1=1,4
                  if (dist(i1).lt.1.d-14) then
                    x_zeta(i1) = 2 * x_zeta(i1)
                    y_zeta(i1) = 2 * y_zeta(i1)
                    z_zeta(i1) = 2 * z_zeta(i1)
                  end if
                end do
              else if (k.eq.kmax) then
                dist(1)=sqrt( (x(i  ,j  ,k+1,n) - x(i  ,j  ,k  ,n))**2 + &
                              (y(i  ,j  ,k+1,n) - y(i  ,j  ,k  ,n))**2 + &
                              (z(i  ,j  ,k+1,n) - z(i  ,j  ,k  ,n))**2 )
                dist(2)=sqrt( (x(i+1,j  ,k+1,n) - x(i+1,j  ,k  ,n))**2 + &
                              (y(i+1,j  ,k+1,n) - y(i+1,j  ,k  ,n))**2 + &
                              (z(i+1,j  ,k+1,n) - z(i+1,j  ,k  ,n))**2 )
                dist(3)=sqrt( (x(i+1,j+1,k+1,n) - x(i+1,j+1,k  ,n))**2 + &
                              (y(i+1,j+1,k+1,n) - y(i+1,j+1,k  ,n))**2 + &
                              (z(i+1,j+1,k+1,n) - z(i+1,j+1,k  ,n))**2 )
                dist(4)=sqrt( (x(i  ,j+1,k+1,n) - x(i  ,j+1,k  ,n))**2 + &
                              (y(i  ,j+1,k+1,n) - y(i  ,j+1,k  ,n))**2 + &
                              (z(i  ,j+1,k+1,n) - z(i  ,j+1,k  ,n))**2 )
                do i1=1,4
                  if (dist(i1).lt.1.d-14) then
                    x_zeta(i1) = 2 * x_zeta(i1)
                    y_zeta(i1) = 2 * y_zeta(i1)
                    z_zeta(i1) = 2 * z_zeta(i1)
                  end if
                end do
              end if

              x_xi(0)   = ( x_xi(1)   + x_xi(2)   ) / 2
              y_xi(0)   = ( y_xi(1)   + y_xi(2)   ) / 2
              z_xi(0)   = ( z_xi(1)   + z_xi(2)   ) / 2

              x_eta(0)  = ( x_eta(1)  + x_eta(2)  ) / 2
              y_eta(0)  = ( y_eta(1)  + y_eta(2)  ) / 2
              z_eta(0)  = ( z_eta(1)  + z_eta(2)  ) / 2

              x_zeta(0)  = &
                ( x_zeta(1)  + x_zeta(2)  + x_zeta(3)  + x_zeta(4)  ) /4
              y_zeta(0)  = &
                ( y_zeta(1)  + y_zeta(2)  + y_zeta(3)  + y_zeta(4)  ) /4
              z_zeta(0)  = &
                ( z_zeta(1)  + z_zeta(2)  + z_zeta(3)  + z_zeta(4)  ) /4

              xix   =  (z_zeta(0) * y_eta(0)  - z_eta(0)  * y_zeta(0) )/ &
                       volvisc
              xiy   =  (x_zeta(0) * z_eta(0)  - x_eta(0)  * z_zeta(0) )/ &
                       volvisc
              xiz   =  (y_zeta(0) * x_eta(0)  - y_eta(0)  * x_zeta(0) )/ &
                       volvisc
              etax  =  (z_xi(0)   * y_zeta(0) - z_zeta(0) * y_xi(0)   )/ &
                       volvisc
              etay  =  (x_xi(0)   * z_zeta(0) - x_zeta(0) * z_xi(0)   )/ &
                       volvisc
              etaz  =  (y_xi(0)   * x_zeta(0) - y_zeta(0) * x_xi(0)   )/ &
                       volvisc
              zetax =  (z_eta(0)  * y_xi(0)   - z_xi(0)   * y_eta(0)  )/ &
                       volvisc
              zetay =  (x_eta(0)  * z_xi(0)   - x_xi(0)   * z_eta(0)  )/ &
                       volvisc
              zetaz =  (y_eta(0)  * x_xi(0)   - y_xi(0)   * x_eta(0)  )/ &
                       volvisc

              xiderk  (1,i,j,k,n) = xix
              xiderk  (2,i,j,k,n) = xiy
              xiderk  (3,i,j,k,n) = xiz
              etaderk (1,i,j,k,n) = etax
              etaderk (2,i,j,k,n) = etay
              etaderk (3,i,j,k,n) = etaz
              zetaderk(1,i,j,k,n) = zetax
              zetaderk(2,i,j,k,n) = zetay
              zetaderk(3,i,j,k,n) = zetaz

            end do
          end do
        end do

      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine vol_edges(nl,vol)
!-----------------------------------------------------------------------

      use cosa_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,iblock,imax,jmax,kmax,ivol
      real (kind=cosa_real) vol(*)

      do iblock = 1,mynblocks
        imax   = i_imax     (iblock,nl)
        jmax   = j_jmax     (iblock,nl)
        kmax   = k_kmax     (iblock,nl)
        ivol   = 1 + off_p1(iblock,nl)
        call vol_bedges(vol(ivol),imax,jmax,kmax)
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine vol_bedges(vol,imax,jmax,kmax)
!-----------------------------------------------------------------------
!     This routine has to be called before the prol_lin routine.
!     This is because the prol. routines use the volumes of the auxiliary
!     cells of the block edges. The volume values of these cells do not
!     affect the solution but they must be nonzero, otherwise floating
!     point exception occur in the prol_lin_b routine.
!-----------------------------------------------------------------------
!     ijkc: array of size 2 x 4 x 12 . First component is (i,j) or (i,k)
!           or (j,k); second
!           component is cell index of particulr cell of onsidered edge
!           region; third component denotes edge region.
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax
      integer(kind=cosa_int) ic,ijkc(2,0:3,12),i,j,k
      real (kind=cosa_real) vol(0:imax,0:jmax,0:kmax)
      real (kind=cosa_real) a,b,c

!---- region 1
      ijkc(1,0,1) = 0
      ijkc(2,0,1) = 0
      ijkc(1,1,1) = 0
      ijkc(2,1,1) = 1
      ijkc(1,2,1) = 1
      ijkc(2,2,1) = 0
      ijkc(1,3,1) = 1
      ijkc(2,3,1) = 1

!---- region 2
      ijkc(1,0,2) = imax
      ijkc(2,0,2) = 0
      ijkc(1,1,2) = imax
      ijkc(2,1,2) = 1
      ijkc(1,2,2) = imax-1
      ijkc(2,2,2) = 0
      ijkc(1,3,2) = imax-1
      ijkc(2,3,2) = 1

!---- region 3
      ijkc(1,0,3) = imax
      ijkc(2,0,3) = jmax
      ijkc(1,1,3) = imax
      ijkc(2,1,3) = jmax-1
      ijkc(1,2,3) = imax-1
      ijkc(2,2,3) = jmax
      ijkc(1,3,3) = imax-1
      ijkc(2,3,3) = jmax-1

!---- region 4
      ijkc(1,0,4) = 0
      ijkc(2,0,4) = jmax
      ijkc(1,1,4) = 0
      ijkc(2,1,4) = jmax-1
      ijkc(1,2,4) = 1
      ijkc(2,2,4) = jmax
      ijkc(1,3,4) = 1
      ijkc(2,3,4) = jmax-1

!---- region 5
      ijkc(1,0,5) = 0
      ijkc(2,0,5) = 0
      ijkc(1,1,5) = 0
      ijkc(2,1,5) = 1
      ijkc(1,2,5) = 1
      ijkc(2,2,5) = 0
      ijkc(1,3,5) = 1
      ijkc(2,3,5) = 1

!---- region 6
      ijkc(1,0,6) = imax
      ijkc(2,0,6) = 0
      ijkc(1,1,6) = imax
      ijkc(2,1,6) = 1
      ijkc(1,2,6) = imax-1
      ijkc(2,2,6) = 0
      ijkc(1,3,6) = imax-1
      ijkc(2,3,6) = 1

!---- region 7
      ijkc(1,0,7) = imax
      ijkc(2,0,7) = kmax
      ijkc(1,1,7) = imax
      ijkc(2,1,7) = kmax-1
      ijkc(1,2,7) = imax-1
      ijkc(2,2,7) = kmax
      ijkc(1,3,7) = imax-1
      ijkc(2,3,7) = kmax-1

!---- region 8
      ijkc(1,0,8) = 0
      ijkc(2,0,8) = kmax
      ijkc(1,1,8) = 0
      ijkc(2,1,8) = kmax-1
      ijkc(1,2,8) = 1
      ijkc(2,2,8) = kmax
      ijkc(1,3,8) = 1
      ijkc(2,3,8) = kmax-1

!---- region 9
      ijkc(1,0,9) = 0
      ijkc(2,0,9) = 0
      ijkc(1,1,9) = 0
      ijkc(2,1,9) = 1
      ijkc(1,2,9) = 1
      ijkc(2,2,9) = 0
      ijkc(1,3,9) = 1
      ijkc(2,3,9) = 1

!---- region 10
      ijkc(1,0,10) = jmax
      ijkc(2,0,10) = 0
      ijkc(1,1,10) = jmax
      ijkc(2,1,10) = 1
      ijkc(1,2,10) = jmax-1
      ijkc(2,2,10) = 0
      ijkc(1,3,10) = jmax-1
      ijkc(2,3,10) = 1

!---- region 11
      ijkc(1,0,11) = jmax
      ijkc(2,0,11) = kmax
      ijkc(1,1,11) = jmax
      ijkc(2,1,11) = kmax-1
      ijkc(1,2,11) = jmax-1
      ijkc(2,2,11) = kmax
      ijkc(1,3,11) = jmax-1
      ijkc(2,3,11) = kmax-1

!---- region 12
      ijkc(1,0,12) = 0
      ijkc(2,0,12) = kmax
      ijkc(1,1,12) = 0
      ijkc(2,1,12) = kmax-1
      ijkc(1,2,12) = 1
      ijkc(2,2,12) = kmax
      ijkc(1,3,12) = 1
      ijkc(2,3,12) = kmax-1

!---- regions 1 to 4
      do k=1,kmax-1
        do ic = 1,4

          b = ( (vol(ijkc(1,3,ic),ijkc(2,3,ic),k) - &
                 vol(ijkc(1,1,ic),ijkc(2,1,ic),k)) * &
                (ijkc(1,2,ic) - ijkc(1,1,ic)) - &
                (vol(ijkc(1,2,ic),ijkc(2,2,ic),k) - &
                 vol(ijkc(1,1,ic),ijkc(2,1,ic),k)) * &
                (ijkc(1,3,ic) - ijkc(1,1,ic)) ) / &
              ( (ijkc(2,3,ic) - ijkc(2,1,ic)) * &
                (ijkc(1,2,ic) - ijkc(1,1,ic)) - &
                (ijkc(2,2,ic) - ijkc(2,1,ic)) * &
                (ijkc(1,3,ic) - ijkc(1,1,ic)) )

          a = ( vol(ijkc(1,2,ic),ijkc(2,2,ic),k) - &
                vol(ijkc(1,1,ic),ijkc(2,1,ic),k) - &
                b * (ijkc(2,2,ic)-ijkc(2,1,ic)) ) / &
              ( ijkc(1,2,ic) - ijkc(1,1,ic) )

          c = vol(ijkc(1,1,ic),ijkc(2,1,ic),k) - &
              a*ijkc(1,1,ic) -b*ijkc(2,1,ic)

          vol(ijkc(1,0,ic),ijkc(2,0,ic),k) = &
            a*ijkc(1,0,ic) + b*ijkc(2,0,ic) + c

        end do
      end do

!---- regions 5 to 8
      do j=1,jmax-1
        do ic = 5,8

          b = ( (vol(ijkc(1,3,ic),j,ijkc(2,3,ic)) - &
                 vol(ijkc(1,1,ic),j,ijkc(2,1,ic))) * &
                (ijkc(1,2,ic) - ijkc(1,1,ic)) - &
                (vol(ijkc(1,2,ic),j,ijkc(2,2,ic)) - &
                 vol(ijkc(1,1,ic),j,ijkc(2,1,ic))) * &
                (ijkc(1,3,ic) - ijkc(1,1,ic)) ) / &
              ( (ijkc(2,3,ic) - ijkc(2,1,ic)) * &
                (ijkc(1,2,ic) - ijkc(1,1,ic)) - &
                (ijkc(2,2,ic) - ijkc(2,1,ic)) * &
                (ijkc(1,3,ic) - ijkc(1,1,ic)) )

          a = ( vol(ijkc(1,2,ic),j,ijkc(2,2,ic)) - &
                vol(ijkc(1,1,ic),j,ijkc(2,1,ic)) - &
                b * (ijkc(2,2,ic)-ijkc(2,1,ic)) ) / &
              ( ijkc(1,2,ic) - ijkc(1,1,ic) )

          c = vol(ijkc(1,1,ic),j,ijkc(2,1,ic)) - &
              a*ijkc(1,1,ic) -b*ijkc(2,1,ic)

          vol(ijkc(1,0,ic),j,ijkc(2,0,ic)) = &
            a*ijkc(1,0,ic) + b*ijkc(2,0,ic) + c

        end do
      end do

!---- regions 9 to 12
      do i=1,imax-1
        do ic = 9,12

          b = ( (vol(i,ijkc(1,3,ic),ijkc(2,3,ic)) - &
                 vol(i,ijkc(1,1,ic),ijkc(2,1,ic))) * &
                (ijkc(1,2,ic) - ijkc(1,1,ic)) - &
                (vol(i,ijkc(1,2,ic),ijkc(2,2,ic)) - &
                 vol(i,ijkc(1,1,ic),ijkc(2,1,ic))) * &
                (ijkc(1,3,ic) - ijkc(1,1,ic)) ) / &
              ( (ijkc(2,3,ic) - ijkc(2,1,ic)) * &
                (ijkc(1,2,ic) - ijkc(1,1,ic)) - &
                (ijkc(2,2,ic) - ijkc(2,1,ic)) * &
                (ijkc(1,3,ic) - ijkc(1,1,ic)) )

          a = ( vol(i,ijkc(1,2,ic),ijkc(2,2,ic)) - &
                vol(i,ijkc(1,1,ic),ijkc(2,1,ic)) - &
                b * (ijkc(2,2,ic)-ijkc(2,1,ic)) ) / &
              ( ijkc(1,2,ic) - ijkc(1,1,ic) )

          c = vol(i,ijkc(1,1,ic),ijkc(2,1,ic)) - &
              a*ijkc(1,1,ic) -b*ijkc(2,1,ic)

          vol(i,ijkc(1,0,ic),ijkc(2,0,ic)) = &
            a*ijkc(1,0,ic) + b*ijkc(2,0,ic) + c

        end do
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine pcutman_v(nl,percutopo,vol)
!-----------------------------------------------------------------------

      use parallelutils
      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,iblk1,icut,iblk2,imax1,jmax1,kmax1,imax2,jmax2,kmax2, &
                iv1,iv2
      integer(kind=cosa_int) percutopo(22,mynpcuts),a,b,c,i,reqnum
      logical myblock1, myblock2
      real(kind=cosa_real) vol(*)
      integer(kind=cosa_int) sendrequestnum,receiverequestnum
      integer(kind=cosa_int) i_id,i_in_id
      real(kind=cosa_real), allocatable :: lsendarray(:,:),lreceivearray(:,:)
      integer(kind=cosa_int) sendid,recvid

      pnumrecv = 0  ! array syntax
      pnumsend = 0  ! array syntax

      allocate(lsendarray(mrequest_pperproc, maxproc_pcutman))
      allocate(lreceivearray(mrequest_pperproc, maxproc_pcutman))

      do icut = 1,mynpcuts
        iblk1 = percutopo( 1,icut)
        iblk2 = percutopo(10,icut)
        call myblock(iblk1,myblock1,.true.)
        call myblock(iblk2,myblock2,.true.)
        if (myblock1) then
          imax1  = i_imax     (iblk1,nl)
          jmax1  = j_jmax     (iblk1,nl)
          kmax1  = k_kmax     (iblk1,nl)
          iv1    = 1 + off_p1 (iblk1,nl)
        else
          iv1    = 1
        end if
        if (myblock2) then
          imax2  = i_imax     (iblk2,nl)
          jmax2  = j_jmax     (iblk2,nl)
          kmax2  = k_kmax     (iblk2,nl)
          iv2    = 1 + off_p1 (iblk2,nl)
        else
          iv2    = 1
        end if
        call pcut_v(imax1,jmax1,kmax1,vol(iv1),imax2,jmax2,kmax2, &
             vol(iv2),percutopo(1,icut),lreceivearray,preceiverequests, &
             pnumrecv,preceiveindices,lsendarray,psendrequests,pnumsend,iv1, &
             mrequest_pperproc,maxproc_pcutman,maxproc,myprecvid_num, &
             mypsendid_num)
      end do

      do i_id=1,maxprecvid_num
        recvid = myprecvid(i_id)
        call receiveblockdata_fromid(lreceivearray(1,i_id),recvid, &
             pnumrecv(i_id),preceiverequests(i_id))
      end do
      do i_id=1,maxpsendid_num
        sendid = mypsendid(i_id)
        call sendblockdata_toid(lsendarray(1,i_id),sendid, &
             pnumsend(i_id),psendrequests(i_id))
      end do
      do i_id=1,maxprecvid_num
        call waitanymessages(preceiverequests,maxprecvid_num,reqnum)
        do i_in_id=1,pnumrecv(reqnum)
          a     = preceiveindices(i_in_id,1,reqnum)
          b     = preceiveindices(i_in_id,2,reqnum)
          c     = preceiveindices(i_in_id,3,reqnum)
          iv1   = preceiveindices(i_in_id,4,reqnum)
          imax1 = preceiveindices(i_in_id,5,reqnum)
          jmax1 = preceiveindices(i_in_id,6,reqnum)
          kmax1 = preceiveindices(i_in_id,7,reqnum)
          call copypvreceivedata(imax1,jmax1,kmax1,a,b,c,vol(iv1), &
               lreceivearray,i_in_id,reqnum)
        end do
      enddo

      call waitallmessages(psendrequests,maxpsendid_num)

      deallocate(lsendarray, lreceivearray)

      return

      end

!-----------------------------------------------------------------------
      subroutine pcut_v(imax1,jmax1,kmax1,vol1,imax2,jmax2,kmax2,vol2, &
           percutopo,receivearray,receiverequests,numrecv, &
           receiveindices,sendarray,sendrequests,numsend,iv1, &
           mrequest_perproc,maxproc_cutman,maxproc,myrecvid_num, &
           mysendid_num)
!-----------------------------------------------------------------------

!     Routine to do cut boundary condition.
!     Flow data stored in Q1 are updated from data in the interior of Q2.
!      * sendrequests/receivrequests have indexes which is the process I am communicating to
!      * sendarray/recvarray have two components:
!        - index of point to send/recv
!        - index of process I have to send/recv to/from
!      * receiveindices have three components:
!        - index of point to send/recv
!        - index of information to trasmit (1-6)
!        - index of process to send/recv to/from

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax1,jmax1,kmax1,imax2,jmax2,kmax2,iblk1,iblk2
      integer(kind=cosa_int) percutopo(22),mrequest
      integer(kind=cosa_int) ijkmax1(3),ijkmax2(3),istr1(3),iend1(3),istr2(3), &
           iend2(3),isgn1(3),isgn2(3),iord(3),len(3),idir1,idir2,inout1, &
           inout2,ibcpt,ibcpt2,inr,inr2,l,ic1,ic2,ic3,jc1,jc2,jc3,i2,i3, &
           ii1,jj1,kk1,in1,jn1,kn1
      real(kind=cosa_real) &
        vol1(0:imax1,0:jmax1,0:kmax1), &
        vol2(0:imax2,0:jmax2,0:kmax2)
      integer(kind=cosa_int) tempindex, datasize
      integer(kind=cosa_int) receiverequestnum,iv1,sendrequestnum
      logical myblock1, myblock2
      integer(kind=cosa_int) mrequest_perproc,maxproc_cutman
      integer(kind=cosa_int) sendid,recvid,numrecvid,numsendid,maxproc
      integer(kind=cosa_int) sendrequests(maxproc_cutman), &
                receiverequests(maxproc_cutman)
      integer(kind=cosa_int) receiveindices(mrequest_perproc,8,maxproc_cutman)
      real(kind=cosa_real) &
        sendarray   (mrequest_perproc,maxproc_cutman), &
        receivearray(mrequest_perproc,maxproc_cutman)
      integer(kind=cosa_int) numsend(maxproc_cutman),numrecv(maxproc_cutman)
      integer(kind=cosa_int) myrecvid_num(0:maxproc-1),mysendid_num(0:maxproc-1)
!--------------------------------------------------------------------------

!     Store imax, jmax in ijmax for looping
      ijkmax1(1) = imax1
      ijkmax1(2) = jmax1
      ijkmax1(3) = kmax1
      ijkmax2(1) = imax2
      ijkmax2(2) = jmax2
      ijkmax2(3) = kmax2

!     Store boundary condition data in mnemonic names

      iblk1    = percutopo( 1)
      idir1    = percutopo( 2)
      inout1   = percutopo( 3)
      istr1(1) = percutopo( 4)
      iend1(1) = percutopo( 5)
      istr1(2) = percutopo( 6)
      iend1(2) = percutopo( 7)
      istr1(3) = percutopo( 8)
      iend1(3) = percutopo( 9)

      iblk2    = percutopo(10)
      idir2    = percutopo(11)
      inout2   = percutopo(12)
      istr2(1) = percutopo(13)
      iend2(1) = percutopo(14)
      istr2(2) = percutopo(15)
      iend2(2) = percutopo(16)
      istr2(3) = percutopo(17)
      iend2(3) = percutopo(18)

      iord(1)  = percutopo(19)
      iord(2)  = percutopo(20)
      iord(3)  = percutopo(21)

      call myblock(iblk1,myblock1,.false.)
      call myblock(iblk2,myblock2,.false.)

!     Set needed variables depending on whether the boundary is
!     the inner boundary (inout1 = 1) or the outer boundary (inout1 > 1)
!         ibcpt  = boundary point of block 1
!         ibcinc = increment to second boundary point of block 1
!         ibcpt2 = ibcpt + ibcinc
!         inr    = interior point of block 2.
!         inrinc = increment to second interior point of block 2
!                  inr2 = inr + inrinc

      if (inout1 .eq. 1) then
        ibcpt  =   0
      else
        ibcpt  = ijkmax1(idir1)
      end if
!     
      if (inout2 .eq. 1) then
        inr    =   1
      else
        inr    =   ijkmax2(idir2) - 1
      end if

!     Find the length of the two outer loops and loop over these using
!     offsets and delta function to set the two cut data points to the
!     two interior data points of block 2.

      do l = 1, 3

         len(l) = abs ( iend1(l) - istr1(l) )

!     Increment/Decrement 

         if ( iend1(l) .gt. istr1(l) ) then
            isgn1(l) =   1
         else
            isgn1(l) = - 1
         end if

!     Increment/Decrement 

         if ( iend2(l) .gt. istr2(l) ) then
            isgn2(l) =   1
         else
            isgn2(l) = - 1
         end if

      end do

!     ii1 first  boundary point of block 1
!     ii2 second boundary point of block 1
!     in1 first  interior point of block 2
!     in2 second interior point of block 2

      ic1 = cyc (idir1, 1)
      ic2 = cyc (idir1, 2)
      ic3 = cyc (idir1, 3)

      jc1 = iord (ic1)
      jc2 = iord (ic2)
      jc3 = iord (ic3)

      do i3 = 0, len(ic3)
        do i2 = 0, len(ic2)

          ii1 = ibcpt                        * krd (ic1, 1) + &
                (istr1(ic2) + isgn1(ic2)*i2) * krd (ic2, 1) + &
                (istr1(ic3) + isgn1(ic3)*i3) * krd (ic3, 1)
          jj1 = ibcpt                        * krd (ic1, 2) + &
                (istr1(ic2) + isgn1(ic2)*i2) * krd (ic2, 2) + &
                (istr1(ic3) + isgn1(ic3)*i3) * krd (ic3, 2)
          kk1 = ibcpt                        * krd (ic1, 3) + &
                (istr1(ic2) + isgn1(ic2)*i2) * krd (ic2, 3) + &
                (istr1(ic3) + isgn1(ic3)*i3) * krd (ic3, 3)

          in1 = inr                          * krd (jc1, 1) + &
                (istr2(jc2) + isgn2(jc2)*i2) * krd (jc2, 1) + &
                (istr2(jc3) + isgn2(jc3)*i3) * krd (jc3, 1)
          jn1 = inr                          * krd (jc1, 2) + &
                (istr2(jc2) + isgn2(jc2)*i2) * krd (jc2, 2) + &
                (istr2(jc3) + isgn2(jc3)*i3) * krd (jc3, 2)
          kn1 = inr                          * krd (jc1, 3) + &
                (istr2(jc2) + isgn2(jc2)*i2) * krd (jc2, 3) + &
                (istr2(jc3) + isgn2(jc3)*i3) * krd (jc3, 3)

         if(myblock1 .and. myblock2) then

           vol1(ii1,jj1,kk1) = vol2(in1,jn1,kn1)

         else if(myblock1) then

           call getownerid(iblk2,recvid)
           numrecvid = myrecvid_num(recvid)
           numrecv(numrecvid) = numrecv(numrecvid)+1 
           receiveindices(numrecv(numrecvid),1,numrecvid) = ii1
           receiveindices(numrecv(numrecvid),2,numrecvid) = jj1
           receiveindices(numrecv(numrecvid),3,numrecvid) = kk1
           receiveindices(numrecv(numrecvid),4,numrecvid) = iv1
           receiveindices(numrecv(numrecvid),5,numrecvid) = imax1
           receiveindices(numrecv(numrecvid),6,numrecvid) = jmax1
           receiveindices(numrecv(numrecvid),7,numrecvid) = kmax1
           if (numrecv(numrecvid).gt. mrequest_perproc) then
             write(*,*) 'numrecv (pcut_v)--> INCREASE MREQUEST_PERPROC'
             call abortmpi()
           end if

         else if(myblock2) then

           call getownerid(iblk1,sendid)
           numsendid = mysendid_num(sendid)
           numsend(numsendid) = numsend(numsendid)+1 
           sendarray(numsend(numsendid),numsendid) = vol2(in1,jn1,kn1)
           if (numsend(numsendid).gt. mrequest_perproc) then
             write(*,*) 'numsend (pcut_v)--> INCREASE MREQUEST_PERPROC'
             call abortmpi()
           end if

         end if

        end do
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine pcutman_x(nl,percutopo,x,y,z)
!-----------------------------------------------------------------------

      use parallelutils
      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,iblk1,icut,iblk2,imax1,jmax1,kmax1,imax2,jmax2,kmax2, &
                ixyz1,ixyz2
      integer(kind=cosa_int) percutopo(22,mynpcuts),a,b,c,i,reqnum,n
      logical myblock1, myblock2
      real(kind=cosa_real) x(*),y(*),z(*)
      integer(kind=cosa_int) sendrequestnum,receiverequestnum
      integer(kind=cosa_int) i_id,i_in_id
      integer(kind=cosa_int) sendid,recvid

      pnumrecv = 0  ! array syntax
      pnumsend = 0  ! array syntax
      sendrecv_datasize = dim5h

      allocate(sendarray(sendrecv_datasize, mrequest_pperproc, &
                         maxproc_pcutman))
      allocate(receivearray(sendrecv_datasize, mrequest_pperproc, &
                            maxproc_pcutman))

      do icut = 1,mynpcuts
        iblk1 = percutopo( 1,icut)
        iblk2 = percutopo(10,icut)
        call myblock(iblk1,myblock1,.true.)
        call myblock(iblk2,myblock2,.true.)
        if(myblock1) then
          imax1  = i_imax     (iblk1,nl)
          jmax1  = j_jmax     (iblk1,nl)
          kmax1  = k_kmax     (iblk1,nl)
          ixyz1  = 1 + off_p2 (iblk1,nl) * dim5h
        else
          ixyz1  = 1
        end if
        if(myblock2) then
          imax2  = i_imax     (iblk2,nl)
          jmax2  = j_jmax     (iblk2,nl)
          kmax2  = k_kmax     (iblk2,nl)
          ixyz2  = 1 + off_p2 (iblk2,nl) * dim5h
        else
          ixyz2  = 1
        end if
        call pcut_x(imax1,jmax1,kmax1,x(ixyz1),y(ixyz1),z(ixyz1),imax2, &
             jmax2,kmax2,x(ixyz2),y(ixyz2),z(ixyz2),nharms,dim5h, &
             percutopo(1,icut),receivearray,preceiverequests,pnumrecv, &
             preceiveindices,sendarray,psendrequests,pnumsend,ixyz1, &
             mrequest_pperproc,maxproc_pcutman,maxproc,myprecvid_num, &
             mypsendid_num)
      end do

      do i_id=1,maxprecvid_num
        recvid = myprecvid(i_id)
        call receiveblockdata_fromid(receivearray(1,1,i_id),recvid, &
             pnumrecv(i_id)*sendrecv_datasize,preceiverequests(i_id))
      enddo
      do i_id=1,maxpsendid_num
        sendid = mypsendid(i_id)
        call sendblockdata_toid(sendarray(1,1,i_id),sendid, &
             pnumsend(i_id)*sendrecv_datasize,psendrequests(i_id))
      enddo
      do i_id=1,maxprecvid_num
        call waitanymessages(preceiverequests,maxprecvid_num,reqnum)
        do i_in_id=1,pnumrecv(reqnum)
          a     = preceiveindices(i_in_id,1,reqnum)
          b     = preceiveindices(i_in_id,2,reqnum)
          c     = preceiveindices(i_in_id,3,reqnum)
          ixyz1 = preceiveindices(i_in_id,4,reqnum)
          imax1 = preceiveindices(i_in_id,6,reqnum)
          jmax1 = preceiveindices(i_in_id,7,reqnum)
          kmax1 = preceiveindices(i_in_id,8,reqnum)
          if (preceiveindices(i_in_id,5,reqnum) .eq. 0) then
            call copypxreceivedata(imax1,jmax1,kmax1,a,b,c,x(ixyz1), &
                 receivearray,i_in_id,reqnum)
          else if (preceiveindices(i_in_id,5,reqnum) .eq. 1) then
            call copypxreceivedata(imax1,jmax1,kmax1,a,b,c,y(ixyz1), &
                 receivearray,i_in_id,reqnum)
          else if (preceiveindices(i_in_id,5,reqnum) .eq. 2) then
            call copypxreceivedata(imax1,jmax1,kmax1,a,b,c,z(ixyz1), &
                 receivearray,i_in_id,reqnum)
          endif
        end do
      end do

      call waitallmessages(psendrequests,maxpsendid_num)

      deallocate(sendarray, receivearray)

      return

      end

!-----------------------------------------------------------------------
      subroutine pcut_x(imax1,jmax1,kmax1,x1,y1,z1,imax2,jmax2,kmax2,x2, &
        y2,z2,nharms,dim5h,percutopo,receivearray,receiverequests, &
        numrecv,receiveindices,sendarray,sendrequests,numsend,ixyz1, &
        mrequest_perproc,maxproc_cutman,maxproc,myrecvid_num, &
        mysendid_num)
!-----------------------------------------------------------------------

!     Routine to do cut boundary condition.
!     Flow data stored in Q1 are updated from data in the interior of Q2.
!      * sendrequests/receivrequests have indexes which is the process I am communicating to
!      * sendarray/recvarray have three components:
!        - index spanning q or mut component to send/recv
!        - index of point to send/recv
!        - index of process I have to send/recv to/from
!      * receiveindices have three components:
!        - index of point to send/recv
!        - index of information to trasmit (1-6)
!        - index of process to send/recv to/from

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax1,jmax1,kmax1,imax2,jmax2,kmax2,nharms,dim5h,iblk1, &
                iblk2
      integer(kind=cosa_int) percutopo(22),ixyz1,mrequest
      integer(kind=cosa_int) ijkmax1(3),ijkmax2(3),istr1(3),iend1(3),istr2(3), &
           iend2(3),isgn1(3),isgn2(3),iord(3),irot,len(3),idir1,idir2, &
           inout1,inout2,ibcpt,ibcinc,ibcpt2,inr,inrinc,inr2,l,n,ic1, &
           ic2,ic3,jc1,jc2,jc3,i2,i3,ii1,jj1,kk1,in1,jn1,kn1
      real(kind=cosa_real) &
           x1(0:imax1+1,0:jmax1+1,0:kmax1+1,0:2*nharms*hbmove), &
           y1(0:imax1+1,0:jmax1+1,0:kmax1+1,0:2*nharms*hbmove), &
           z1(0:imax1+1,0:jmax1+1,0:kmax1+1,0:2*nharms*hbmove), &
           x2(0:imax2+1,0:jmax2+1,0:kmax2+1,0:2*nharms*hbmove), &
           y2(0:imax2+1,0:jmax2+1,0:kmax2+1,0:2*nharms*hbmove), &
           z2(0:imax2+1,0:jmax2+1,0:kmax2+1,0:2*nharms*hbmove)
      integer(kind=cosa_int) tempindex, datasize
      integer(kind=cosa_int) receiverequestnum,sendrequestnum
      logical myblock1, myblock2
      integer(kind=cosa_int) mrequest_perproc,maxproc_cutman
      integer(kind=cosa_int) sendid,recvid,numrecvid,numsendid,maxproc
      integer(kind=cosa_int) sendrequests(maxproc_cutman), &
                receiverequests(maxproc_cutman)
      integer(kind=cosa_int) receiveindices(mrequest_perproc,8,maxproc_cutman)
      real(kind=cosa_real) &
!del    sendarray   (2*nharms+1,mrequest_perproc,maxproc_cutman), &
!del    receivearray(2*nharms+1,mrequest_perproc,maxproc_cutman) &
        sendarray   (dim5h,mrequest_perproc,maxproc_cutman), &
        receivearray(dim5h,mrequest_perproc,maxproc_cutman)
      real(kind=cosa_real) tmp(2)
      integer(kind=cosa_int) numsend(maxproc_cutman),numrecv(maxproc_cutman)
      integer(kind=cosa_int) myrecvid_num(0:maxproc-1),mysendid_num(0:maxproc-1)
!--------------------------------------------------------------------------

!     Store imax, jmax in ijmax for looping
      ijkmax1(1) = imax1
      ijkmax1(2) = jmax1
      ijkmax1(3) = kmax1
      ijkmax2(1) = imax2
      ijkmax2(2) = jmax2
      ijkmax2(3) = kmax2

!     Store boundary condition data in mnemonic names

      iblk1    = percutopo( 1)
      idir1    = percutopo( 2)
      inout1   = percutopo( 3)
      istr1(1) = percutopo( 4)
      iend1(1) = percutopo( 5)
      istr1(2) = percutopo( 6)
      iend1(2) = percutopo( 7)
      istr1(3) = percutopo( 8)
      iend1(3) = percutopo( 9)

      iblk2    = percutopo(10)
      idir2    = percutopo(11)
      inout2   = percutopo(12)
      istr2(1) = percutopo(13)
      iend2(1) = percutopo(14)
      istr2(2) = percutopo(15)
      iend2(2) = percutopo(16)
      istr2(3) = percutopo(17)
      iend2(3) = percutopo(18)

      iord(1)  = percutopo(19)
      iord(2)  = percutopo(20)
      iord(3)  = percutopo(21)

!---- msc, 09/05/2021: irot definition was missing and I added it now
      irot     = percutopo(22)

      call myblock(iblk1,myblock1,.false.)
      call myblock(iblk2,myblock2,.false.)

!     Set needed variables depending on whether the boundary is
!     the inner boundary (inout1 = 1) or the outer boundary (inout1 > 1)
!         ibcpt  = boundary point of block 1
!         ibcinc = increment to second boundary point of block 1
!         ibcpt2 = ibcpt + ibcinc
!         inr    = interior point of block 2.
!         inrinc = increment to second interior point of block 2
!                  inr2 = inr + inrinc

      if (inout1 .eq. 1) then
        ibcpt  = 0
      else
        ibcpt  = ijkmax1(idir1) + 1
      end if
!     
      if (inout2 .eq. 1) then
        inr    = 2
      else
        inr    = ijkmax2(idir2) - 1
      end if

      do l = 1, 3
        if (l .ne. idir1) then
          if (iend1(l) .gt. istr1(l)) then
            iend1(l)  = iend1(l) + 1
          else
            istr1(l)  = istr1(l) + 1
          end if
        end if
      end do

      do l = 1, 3
        if (l .ne. idir2) then
          if (iend2(l) .gt. istr2(l)) then
            iend2(l)  = iend2(l) + 1
          else
            istr2(l)  = istr2(l) + 1
          end if
        end if
      end do

!     Find the length of the two outer loops and loop over these using
!     offsets and delta function to set the two cut data points to the
!     two interior data points of block 2.

      do l = 1, 3

        len(l) = abs ( iend1(l) - istr1(l) )

!       Increment/Decrement 

        if ( iend1(l) .gt. istr1(l) ) then
          isgn1(l) =   1
        else
          isgn1(l) = - 1
        end if

!       Increment/Decrement 

        if ( iend2(l) .gt. istr2(l) ) then
          isgn2(l) =   1
        else
          isgn2(l) = - 1
        end if

      end do

!     ii1 first  boundary point of block 1
!     ii2 second boundary point of block 1
!     in1 first  interior point of block 2
!     in2 second interior point of block 2

      ic1 = cyc (idir1, 1)
      ic2 = cyc (idir1, 2)
      ic3 = cyc (idir1, 3)

      jc1 = iord (ic1)
      jc2 = iord (ic2)
      jc3 = iord (ic3)

      do i3 = 0, len(ic3)
        do i2 = 0, len(ic2)

          ii1 = ibcpt                        * krd (ic1, 1) + &
                (istr1(ic2) + isgn1(ic2)*i2) * krd (ic2, 1) + &
                (istr1(ic3) + isgn1(ic3)*i3) * krd (ic3, 1)
          jj1 = ibcpt                        * krd (ic1, 2) + &
                (istr1(ic2) + isgn1(ic2)*i2) * krd (ic2, 2) + &
                (istr1(ic3) + isgn1(ic3)*i3) * krd (ic3, 2)
          kk1 = ibcpt                        * krd (ic1, 3) + &
                (istr1(ic2) + isgn1(ic2)*i2) * krd (ic2, 3) + &
                (istr1(ic3) + isgn1(ic3)*i3) * krd (ic3, 3)

          in1 = inr                          * krd (jc1, 1) + &
                (istr2(jc2) + isgn2(jc2)*i2) * krd (jc2, 1) + &
                (istr2(jc3) + isgn2(jc3)*i3) * krd (jc3, 1)
          jn1 = inr                          * krd (jc1, 2) + &
                (istr2(jc2) + isgn2(jc2)*i2) * krd (jc2, 2) + &
                (istr2(jc3) + isgn2(jc3)*i3) * krd (jc3, 2)
          kn1 = inr                          * krd (jc1, 3) + &
                (istr2(jc2) + isgn2(jc2)*i2) * krd (jc2, 3) + &
                (istr2(jc3) + isgn2(jc3)*i3) * krd (jc3, 3)

          if(myblock1 .and. myblock2) then

            do n = 0, 2*nharms*hbmove
              tmp(1) = x2(in1,jn1,kn1,n) - xrotc
              tmp(2) = y2(in1,jn1,kn1,n) - yrotc
              x1(ii1,jj1,kk1,n) = xrotc &
                + rotmat(1,1) * tmp(1) - rotmat(1,2) * tmp(2) * irot
              y1(ii1,jj1,kk1,n) = yrotc &
                - rotmat(2,1) * tmp(1) * irot + rotmat(2,2) * tmp(2)
              z1(ii1,jj1,kk1,n) = z2(in1,jn1,kn1,n)
            end do

          else if(myblock1) then

            call getownerid(iblk2,recvid)
            numrecvid = myrecvid_num(recvid)
            numrecv(numrecvid) = numrecv(numrecvid)+1 
            receiveindices(numrecv(numrecvid),1,numrecvid) = ii1
            receiveindices(numrecv(numrecvid),2,numrecvid) = jj1
            receiveindices(numrecv(numrecvid),3,numrecvid) = kk1
            receiveindices(numrecv(numrecvid),4,numrecvid) = ixyz1
            receiveindices(numrecv(numrecvid),5,numrecvid) = 0
            receiveindices(numrecv(numrecvid),6,numrecvid) = imax1
            receiveindices(numrecv(numrecvid),7,numrecvid) = jmax1
            receiveindices(numrecv(numrecvid),8,numrecvid) = kmax1
            numrecv(numrecvid) = numrecv(numrecvid)+1
            receiveindices(numrecv(numrecvid),1,numrecvid) = ii1
            receiveindices(numrecv(numrecvid),2,numrecvid) = jj1
            receiveindices(numrecv(numrecvid),3,numrecvid) = kk1
            receiveindices(numrecv(numrecvid),4,numrecvid) = ixyz1
            receiveindices(numrecv(numrecvid),5,numrecvid) = 1
            receiveindices(numrecv(numrecvid),6,numrecvid) = imax1
            receiveindices(numrecv(numrecvid),7,numrecvid) = jmax1
            receiveindices(numrecv(numrecvid),8,numrecvid) = kmax1
            numrecv(numrecvid) = numrecv(numrecvid)+1
            receiveindices(numrecv(numrecvid),1,numrecvid) = ii1
            receiveindices(numrecv(numrecvid),2,numrecvid) = jj1
            receiveindices(numrecv(numrecvid),3,numrecvid) = kk1
            receiveindices(numrecv(numrecvid),4,numrecvid) = ixyz1
            receiveindices(numrecv(numrecvid),5,numrecvid) = 2
            receiveindices(numrecv(numrecvid),6,numrecvid) = imax1
            receiveindices(numrecv(numrecvid),7,numrecvid) = jmax1
            receiveindices(numrecv(numrecvid),8,numrecvid) = kmax1
            if (numrecv(numrecvid).gt. mrequest_perproc) then
              write(*,*) 'numrecv (pcut_x)--> INCREASE MREQUEST_PERPROC'
              call abortmpi()
            end if
          else if (myblock2) then
            call getownerid(iblk1,sendid)
            numsendid = mysendid_num(sendid)
            numsend(numsendid) = numsend(numsendid)+1 
            tempindex = 1
            do n = 0, 2*nharms*hbmove
              tmp(1) = x2(in1,jn1,kn1,n) - xrotc
              tmp(2) = y2(in1,jn1,kn1,n) - yrotc
              sendarray(tempindex,numsend(numsendid),numsendid) = xrotc &
                + rotmat(1,1) * tmp(1) - rotmat(1,2) * tmp(2) * irot
              tempindex = tempindex + 1
            end do
            numsend(numsendid) = numsend(numsendid)+1
            tempindex  = 1
            do n = 0, 2*nharms*hbmove
              tmp(1) = x2(in1,jn1,kn1,n) - xrotc
              tmp(2) = y2(in1,jn1,kn1,n) - yrotc
              sendarray(tempindex,numsend(numsendid),numsendid) = yrotc &
                - rotmat(2,1) * tmp(1) * irot + rotmat(2,2) * tmp(2)
              tempindex = tempindex + 1
            end do
            numsend(numsendid) = numsend(numsendid)+1
            tempindex  = 1
            do n = 0, 2*nharms*hbmove
              sendarray(tempindex,numsend(numsendid),numsendid) = &
                z2(in1,jn1,kn1,n)
              tempindex = tempindex + 1
            end do
            if (numsend(numsendid).gt. mrequest_perproc) then
              write(*,*) 'numsend (pcut_x)--> INCREASE MREQUEST_PERPROC'
              call abortmpi()
            end if
          end if

        end do
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine bodypos()
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) n

      do n=0,2*nharms

        if (pitching) then
          dtheta(n) = dtheta0 * dsin(omegas*zeit(n))
        else if (plunging) then
          continue
        else if (plupitching) then
          dtheta(n) = dtheta0 * dsin(omegas*zeit(n)+phipp)
        else if (rotating) then
          dtheta(n) = omegas*zeit(n)
          if (tpitch) then
            thetatp(n) = theta0tp * dsin(omegatp*zeit(n)+phitp)
          end if
        end if

      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine checkcut_init(cutopo)
!-----------------------------------------------------------------------
!     The subroutine initializes data exchange along cuts for check cuts
!     Each rank of process I have to perform send/receive is associated
!     to a number ranging from 1 to maxsendid_num/maxrecvid_num.
!     The inverse association is provided by another array.
!     More precisely, the results are:
!     * maxsendid_num = number of processes which I have to send to 
!     * maxrecvid_num = number of processes which I have to receive from
!     * mysendid(1:maxsendid_num) ---> from number to rank
!     * myrecvid(1:maxrecvid_num) ---> from number to rank
!     * mysendid_num(1:maxproc) ---> from rank to number
!     * myrecvid_num(1:maxproc) ---> from rank to number
!     Maybe send and receive arrays are always identical, I do not know...
!     * maxproc ----> maximum number of MPI processes
!     * maxproc_check_cutman ---> maximum number of MPI processes one process exchange data
!     * mrequest_check_perproc ---> maximum number of exchanged points for each send/recv operation
!-----------------------------------------------------------------------

      use parallelutils
      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) iblk1,icut,iblk2
      integer(kind=cosa_int) cutopo(21,myncuts)
      logical myblock1, myblock2
      integer(kind=cosa_int) length(3),istr1(3),iend1(3),idir1,l,ic2,ic3
      integer(kind=cosa_int) recvfromid,sendtoid,i_id,rank,ierr,msize
      logical found
      integer(kind=cosa_int) maxsend(0:maxproc-1), maxrecv(0:maxproc-1)

      maxsend = 0
      maxrecv = 0

      call allocate_check_init_memory()

      check_maxsendid_num = 0
      check_maxrecvid_num = 0
      check_mysendid = -1 ! array syntax
      check_myrecvid = -1 ! array syntax 
      check_myrecvid_num = -1 ! array syntax
      check_mysendid_num = -1 ! array syntax

      do icut = 1,myncuts
         iblk1 = cutopo( 1,icut)
         iblk2 = cutopo(10,icut)
         idir1    = cutopo(2,icut)
         istr1(1) = cutopo(4,icut)
         iend1(1) = cutopo(5,icut)
         istr1(2) = cutopo(6,icut)
         iend1(2) = cutopo(7,icut)
         istr1(3) = cutopo(8,icut)
         iend1(3) = cutopo(9,icut)
         do l = 1, 3
            if (l .ne. idir1) then
               if (iend1(l) .gt. istr1(l)) then
                  iend1(l)  = iend1(l) + 1
               else
                  istr1(l)  = istr1(l) + 1
              end if
            end if
         end do         
         do l = 1, 3
            length(l) = abs ( iend1(l) - istr1(l) )
         end do
         ic2 = cyc (idir1, 2)
         ic3 = cyc (idir1, 3)
         call myblock(iblk1,myblock1,.false.)
         call myblock(iblk2,myblock2,.false.)
         if(myblock1 .and. myblock2) then
         else if(myblock1) then
!           print*,'iblk2: ',iblk2
            call getownerid(iblk2,recvfromid)
            found = .false.
            do i_id=1,check_maxrecvid_num
               if(check_myrecvid(i_id) .eq. recvfromid) found = .true. 
            enddo
            if(.not.found) then
              check_maxrecvid_num = check_maxrecvid_num + 1
              check_myrecvid(check_maxrecvid_num) = recvfromid
              check_myrecvid_num(recvfromid) = check_maxrecvid_num
            endif
            maxrecv(recvfromid) = maxrecv(recvfromid) + &
                                  (length(ic3)+1)*(length(ic2)+1)*3
            found = .false.
            do i_id=1,check_maxsendid_num
               if(check_mysendid(i_id) .eq. recvfromid) found = .true. 
            enddo
            if(.not.found) then
              check_maxsendid_num = check_maxsendid_num + 1
              check_mysendid(check_maxsendid_num) = recvfromid
              check_mysendid_num(recvfromid) = check_maxsendid_num
            endif
            maxsend(recvfromid) = maxsend(recvfromid) + &
                                  (length(ic3)+1)*(length(ic2)+1)*3
         else if(myblock2) then
            call getownerid(iblk1,sendtoid)
!            print*,'iblk1: ',iblk1
            found = .false.
            do i_id=1,check_maxsendid_num
               if(check_mysendid(i_id) .eq. sendtoid) found = .true. 
            enddo
            if(.not.found) then
              check_maxsendid_num = check_maxsendid_num + 1
              check_mysendid(check_maxsendid_num) = sendtoid
              check_mysendid_num(sendtoid) = check_maxsendid_num
            endif
            maxsend(sendtoid) = maxsend(sendtoid) + &
                                (length(ic3)+1)*(length(ic2)+1)*3
            found = .false.
            do i_id=1,check_maxrecvid_num
               if(check_myrecvid(i_id) .eq. sendtoid) found = .true. 
            enddo
            if(.not.found) then
              check_maxrecvid_num = check_maxrecvid_num + 1
              check_myrecvid(check_maxrecvid_num) = sendtoid
              check_myrecvid_num(sendtoid) = check_maxrecvid_num
            endif
            maxrecv(sendtoid) = maxrecv(sendtoid) + &
                                (length(ic3)+1)*(length(ic2)+1)*3
         end if
      end do

      mrequest_check_perproc = max(maxval(maxsend), maxval(maxrecv))
      maxproc_check_cutman = max(check_maxsendid_num,check_maxrecvid_num)

      call allocate_check_memory()

! debug      call getmpiid(rank)
! debug      call getmpisize(msize)
! debug      write(100+rank,*) 'Rank:',rank,'send to: ',mysendid(1:maxsendid_num)
! debug      write(100+rank,*) 'Rank:',rank,'recv fr: ',myrecvid(1:maxrecvid_num)
! debug      write(100+rank,*) 'Rank:',rank,'to proc: ',mysendid_num(0:msize)
! debug      write(100+rank,*) 'Rank:',rank,'fr proc: ',myrecvid_num(0:msize)
! debug      CALL MPI_FINALIZE(ierr)
! debug      STOP

      return
      end

!-----------------------------------------------------------------------
      subroutine check_cuts(nl,cutopo,x,y,z)
!-----------------------------------------------------------------------

      use parallelutils
      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,iblk1,icut,iblk2,imax1,jmax1,kmax1,imax2,jmax2,kmax2, &
                ixyz1,ixyz2,ii1,jj1,kk1,in1,jn1,kn1
      integer(kind=cosa_int) cutopo(21,myncuts),i,reqnum,n
      logical myblock1,myblock2
      integer(kind=cosa_int) mxd_ijk(2,3,0:2*nharms*hbmove,myncuts),temp_ijk(2,3), &
                mxd_loc(0:2*nharms*hbmove)
      real(kind=cosa_real) avd(0:2*nharms*hbmove), &
                   maxharmdiff(0:2*nharms*hbmove), &
                   mxd(0:2*nharms*hbmove,myncuts), &
                   mxd_xyz(2,3,0:2*nharms*hbmove,myncuts),temp_xyz(2,3), &
                   tempdata(0:2*nharms*hbmove,3)
      real(kind=cosa_real) x(*),y(*),z(*)
      integer(kind=cosa_int) sendrequestnum,receiverequestnum
      integer(kind=cosa_int) requests(4)
      integer(kind=cosa_int) i_id,i_in_id
      integer(kind=cosa_int) sendid,recvid,blocknum
      integer maxdiff_loc(0:2*nharms*hbmove)
      integer gcut_num(0:2*nharms*hbmove)
      integer myid
      logical amcontrol

      call checkcut_init(cutopo)

      check_numrecv = 0  ! array syntax
      check_numsend = 0  ! array syntax
      sendrecv_datasize = dim5h

      allocate(sendarray(sendrecv_datasize, mrequest_check_perproc, &
                         maxproc_check_cutman))
      allocate(receivearray(sendrecv_datasize, mrequest_check_perproc, &
                            maxproc_check_cutman))

      sendarray = 0.d0
      receivearray = 0.d0
      mxd = 1.d-28

      call amcontroller(amcontrol)
      call getmpiid(myid)

      do icut = 1,myncuts
        iblk1 = cutopo( 1,icut)
        iblk2 = cutopo(10,icut)
        call myblock(iblk1,myblock1,.true.)
        call myblock(iblk2,myblock2,.true.)
        if(myblock1) then
          imax1  = i_imax     (iblk1,nl)
          jmax1  = j_jmax     (iblk1,nl)
          kmax1  = k_kmax     (iblk1,nl)
          ixyz1  = 1 + off_p2 (iblk1,nl) * dim5h
        else
          ixyz1  = 1
        end if
        if(myblock2) then
          imax2  = i_imax     (iblk2,nl)
          jmax2  = j_jmax     (iblk2,nl)
          kmax2  = k_kmax     (iblk2,nl)
          ixyz2  = 1 + off_p2 (iblk2,nl) * dim5h
        else
          ixyz2  = 1
        end if
        call check_cut(imax1,jmax1,kmax1,x(ixyz1),y(ixyz1),z(ixyz1), &
             imax2,jmax2,kmax2,x(ixyz2),y(ixyz2),z(ixyz2),nharms,dim5h, &
             cutopo(1,icut),mxd_ijk(:,:,:,icut),mxd(:,icut), &
             mxd_xyz(:,:,:,icut),check_receiverequests, &
             check_receiveindices,sendarray,check_sendrequests, &
             ixyz1,ixyz2,mrequest_check_perproc, &
             maxproc_check_cutman,maxproc,myblock1,myblock2, &
             check_myrecvid_num, check_mysendid_num, check_numrecv, &
             check_numsend, icut)
      end do

      do i_id=1,check_maxrecvid_num
         recvid = check_myrecvid(i_id)
         call receiveblockdata_fromid(receivearray(1,1,i_id),recvid, &
              check_numrecv(i_id)*sendrecv_datasize, &
              check_receiverequests(i_id))
      enddo
      do i_id=1,check_maxsendid_num
         sendid = check_mysendid(i_id)
         call sendblockdata_toid(sendarray(1,1,i_id),sendid, &
              check_numsend(i_id)*sendrecv_datasize, &
              check_sendrequests(i_id))
      enddo

      call waitallmessages(check_sendrequests,check_maxsendid_num)

      do i_id=1,check_maxrecvid_num
         call waitanymessages(check_receiverequests,check_maxrecvid_num, &
                              reqnum)
         do i_in_id=1,check_numrecv(reqnum),3
            ii1 = check_receiveindices(i_in_id,1,reqnum)
            jj1 = check_receiveindices(i_in_id,2,reqnum)
            kk1 = check_receiveindices(i_in_id,3,reqnum)
            ixyz1 = check_receiveindices(i_in_id,4,reqnum)
            imax1 = check_receiveindices(i_in_id,5,reqnum)
            jmax1 = check_receiveindices(i_in_id,6,reqnum)
            kmax1 = check_receiveindices(i_in_id,7,reqnum)
            in1 = check_receiveindices(i_in_id,8,reqnum)
            jn1 = check_receiveindices(i_in_id,9,reqnum)
            kn1 = check_receiveindices(i_in_id,10,reqnum)
            ixyz2 = check_receiveindices(i_in_id,11,reqnum)
            imax2 = check_receiveindices(i_in_id,12,reqnum)
            jmax2 = check_receiveindices(i_in_id,13,reqnum)
            kmax2 = check_receiveindices(i_in_id,14,reqnum)            
!     This entry is set to 2 if I am receiving block 2 and 1 if I am receiving block 1
            blocknum =  check_receiveindices(i_in_id,15,reqnum)
            icut =  check_receiveindices(i_in_id,16,reqnum)
            if(blocknum .eq. 2) then
!     copy block1 data into tempdata for the next calculatio
               call getblockdata(imax1,jmax1,kmax1,x(ixyz1),y(ixyz1), &
                    z(ixyz1),nharms,ii1,jj1,kk1,tempdata)
               call check_received_cut( &
                    tempdata(:,1),tempdata(:,2),tempdata(:,3), &
                    receivearray(:,i_in_id,reqnum), &
                    receivearray(:,i_in_id+1,reqnum), &
                    receivearray(:,i_in_id+2,reqnum), &
                    nharms,mxd_ijk(:,:,:,icut),mxd(:,icut), &
                    mxd_xyz(:,:,:,icut),ii1,jj1,kk1,in1,jn1,kn1)
            else
              call getblockdata(imax2,jmax2,kmax2,x(ixyz2),y(ixyz2), &
                    z(ixyz2),nharms,in1,jn1,kn1,tempdata)
               call check_received_cut(receivearray(:,i_in_id,reqnum), &
                    receivearray(:,i_in_id+1,reqnum), &
                    receivearray(:,i_in_id+2,reqnum), &
                    tempdata(:,1),tempdata(:,2),tempdata(:,3), &
                    nharms,mxd_ijk(:,:,:,icut),mxd(:,icut), &
                    mxd_xyz(:,:,:,icut),ii1,jj1,kk1,in1,jn1,kn1)
            end if
         end do
      end do

      do n = 0,2*nharms*hbmove
         maxharmdiff(n) = 1.d-28 
         mxd_loc(n) = 0
         maxdiff_loc(n) = 0
      end do

      do n = 0,2*nharms*hbmove
!------- msc, 09/05/2021: avd initialisation was missing and I added it now
         avd(n) = 0.d0
         do icut = 1,myncuts
            iblk1 = cutopo( 1,icut)
            call myblock(iblk1,myblock1,.false.)
!           This is required to ensure that cuts aren't double 
!           counted in parallel
            if (mxd(n,icut).gt.maxharmdiff(n)) then
               maxharmdiff(n) = mxd(n,icut)
               mxd_loc(n) = icut
            end if
            if (myblock1) then
               avd(n) = avd(n) + mxd(n,icut)
            end if
         end do     
      end do

      call realsumreduce(avd,(2*nharms*hbmove)+1)
      call realmaxlocallreduce(maxharmdiff,(2*nharms*hbmove)+1,myid, &
           maxdiff_loc)
      if (amcontrol) then
        write(*,*) '     *** Outcome of cut test ***'
      end if
      do n = 0,2*nharms*hbmove         
         if(amcontrol.or.(myid.eq.maxdiff_loc(n))) then
!           If the max location is on the root/control process then no data 
!           needs to be sent, we already have it and can just copy it internally
            if(amcontrol.and.(myid.eq.maxdiff_loc(n))) then
               if(mxd_loc(n).ne.0) then
                  temp_ijk = mxd_ijk(:,:,n,mxd_loc(n))
                  temp_xyz = mxd_xyz(:,:,n,mxd_loc(n)) 
                  gcut_num(n) = mycuts(mxd_loc(n))
               end if
            else
               if(amcontrol) then
                  call mpireceivewt(maxdiff_loc(n),temp_ijk,6, &
                       requests(1),n+1)
                  call mpireceivewt(maxdiff_loc(n),temp_xyz,6, &
                       requests(2),n+2)
                  call mpireceivewt(maxdiff_loc(n),mxd_loc(n),1, &
                       requests(3),n+3)
                  call mpireceivewt(maxdiff_loc(n),gcut_num(n),1, &
                       requests(4),n+4)
                  call waitallmessages(requests,4)
               else
!                If the max location is not on the root/control process then 
!                we need to send the data to that process
                  call mpisendwt(0,mxd_ijk(:,:,n,mxd_loc(n)),6, &
                       requests(1),n+1)
                  call mpisendwt(0,mxd_xyz(:,:,n,mxd_loc(n)),6, &
                       requests(2),n+2)
                  call mpisendwt(0,mxd_loc(n),1,requests(3),n+3)
                  call mpisendwt(0,mycuts(mxd_loc(n)),1,requests(4),n+4)
                  call waitallmessages(requests,4)
               end if
            end if
!           Print out the data on the control process
            if(amcontrol.and.(maxharmdiff(n).gt.1.d-28)) then
               maxharmdiff(n) = dsqrt(maxharmdiff(n))
               avd(n)         = dsqrt(avd(n) / ncuts)
               write(*,*) 'harm., cut no.', n, gcut_num(n)
               write(*,*) 'max diff, average diff', maxharmdiff(n), &
                             avd(n)
               write(*,*) 'i1, j1, k1', temp_ijk(1,1:3)
               write(*,*) 'x1, y1, z1', temp_xyz(1,1:3)
               write(*,*) 'i2, j2, k2', temp_ijk(2,1:3)
               write(*,*) 'x2, y2, z2', temp_xyz(2,1:3)
               write(*,*)
            else if (amcontrol) then
               write(*,*) 'CUT TEST SUCCESSFUL --> harm.', n
            end if
         end if
      end do

      deallocate(sendarray, receivearray)

      call deallocate_check_memory()
      
      return
      end

!-----------------------------------------------------------------------
      subroutine check_cut(imax1,jmax1,kmax1,x1,y1,z1,imax2,jmax2,kmax2, &
        x2,y2,z2,nharms,dim5h,cutopo,mxd_ijk,mxd,mxd_xyz, &
        receiverequests, &
        receiveindices,sendarray,sendrequests,ixyz1, &
        ixyz2,mrequest_perproc,maxproc_cutman,maxproc, &
        myblock1,myblock2, myrecvid_num, mysendid_num, &
        numrecv, numsend, icut)
!-----------------------------------------------------------------------

!     Routine to do check on cuts.
!     Flow data stored in Q1 are updated from data in the interior of Q2.
!      * sendrequests/receivrequests have indexes which is the process I am communicating to
!      * sendarray/recvarray have three components:
!        - index spanning q or mut component to send/recv
!        - index of point to send/recv
!        - index of process I have to send/recv to/from
!      * receiveindices have three components:
!        - index of point to send/recv
!        - index of information to trasmit (1-6)
!        - index of process to send/recv to/from

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax1,jmax1,kmax1,imax2,jmax2,kmax2,nharms,dim5h,iblk1, &
                iblk2,icut
      integer(kind=cosa_int) cutopo(21),ixyz1,ixyz2,mrequest
      integer(kind=cosa_int) ijkmax1(3),ijkmax2(3),istr1(3),iend1(3),istr2(3), &
           iend2(3),isgn1(3),isgn2(3),iord(3),len(3),idir1,idir2,inout1, &
           inout2,ibcpt,ibcinc,ibcpt2,inr,inrinc,inr2,l,n,ic1,ic2,ic3, &
           jc1,jc2,jc3,i2,i3,ii1,jj1,kk1,in1,jn1,kn1
      integer(kind=cosa_int) mxd_ijk(2,3,0:2*nharms*hbmove)
      real(kind=cosa_real) distsq,mxd(0:2*nharms*hbmove), &
           mxd_xyz(2,3,0:2*nharms*hbmove)
      real(kind=cosa_real) &
           x1(0:imax1+1,0:jmax1+1,0:kmax1+1,0:2*nharms*hbmove), &
           y1(0:imax1+1,0:jmax1+1,0:kmax1+1,0:2*nharms*hbmove), &
           z1(0:imax1+1,0:jmax1+1,0:kmax1+1,0:2*nharms*hbmove), &
           x2(0:imax2+1,0:jmax2+1,0:kmax2+1,0:2*nharms*hbmove), &
           y2(0:imax2+1,0:jmax2+1,0:kmax2+1,0:2*nharms*hbmove), &
           z2(0:imax2+1,0:jmax2+1,0:kmax2+1,0:2*nharms*hbmove)
      integer(kind=cosa_int) tempindex, datasize
      integer(kind=cosa_int) receiverequestnum,sendrequestnum
      logical myblock1, myblock2
      integer(kind=cosa_int) mrequest_perproc,maxproc_cutman
      integer(kind=cosa_int) sendid,recvid,maxproc,numsendid,numrecvid
      integer(kind=cosa_int) sendrequests(maxproc_cutman), &
                receiverequests(maxproc_cutman)
      integer(kind=cosa_int) receiveindices(mrequest_perproc,16,maxproc_cutman)
      real(kind=cosa_real) &
           sendarray   (dim5h,mrequest_perproc,maxproc_cutman)
      integer(kind=cosa_int) numsend(mrequest_perproc),numrecv(mrequest_perproc)
      integer(kind=cosa_int) myrecvid_num(0:maxproc-1),mysendid_num(0:maxproc-1)
      integer(kind=cosa_int) inneri
!--------------------------------------------------------------------------

      integer myid

      call getmpiid(myid)

!     Store imax, jmax in ijmax for looping
      ijkmax1(1) = imax1
      ijkmax1(2) = jmax1
      ijkmax1(3) = kmax1
      ijkmax2(1) = imax2
      ijkmax2(2) = jmax2
      ijkmax2(3) = kmax2

!     Store boundary condition data in mnemonic names

      iblk1    = cutopo( 1)
      idir1    = cutopo( 2)
      inout1   = cutopo( 3)
      istr1(1) = cutopo( 4)
      iend1(1) = cutopo( 5)
      istr1(2) = cutopo( 6)
      iend1(2) = cutopo( 7)
      istr1(3) = cutopo( 8)
      iend1(3) = cutopo( 9)

      iblk2    = cutopo(10)
      idir2    = cutopo(11)
      inout2   = cutopo(12)
      istr2(1) = cutopo(13)
      iend2(1) = cutopo(14)
      istr2(2) = cutopo(15)
      iend2(2) = cutopo(16)
      istr2(3) = cutopo(17)
      iend2(3) = cutopo(18)

      iord(1)  = cutopo(19)
      iord(2)  = cutopo(20)
      iord(3)  = cutopo(21)


!     Set needed variables depending on whether the boundary is
!     the inner boundary (inout1 = 1) or the outer boundary (inout1 > 1)
!         ibcpt  = boundary point of block 1
!         ibcinc = increment to second boundary point of block 1
!         ibcpt2 = ibcpt + ibcinc
!         inr    = interior point of block 2.
!         inrinc = increment to second interior point of block 2
!                  inr2 = inr + inrinc

      if (inout1 .eq. 1) then
        ibcpt  = 1
      else
        ibcpt  = ijkmax1(idir1)
      end if
!     
      if (inout2 .eq. 1) then
        inr    = 1
      else
        inr    = ijkmax2(idir2)
      end if

      do l = 1, 3
        if (l .ne. idir1) then
          if (iend1(l) .gt. istr1(l)) then
            iend1(l)  = iend1(l) + 1
          else
            istr1(l)  = istr1(l) + 1
          end if
        end if
      end do

      do l = 1, 3
        if (l .ne. idir2) then
          if (iend2(l) .gt. istr2(l)) then
            iend2(l)  = iend2(l) + 1
          else
            istr2(l)  = istr2(l) + 1
          end if
        end if
      end do

!     Find the length of the two outer loops and loop over these using
!     offsets and delta function to set the two cut data points to the
!     two interior data points of block 2.

      do l = 1, 3

        len(l) = abs ( iend1(l) - istr1(l) )

!       Increment/Decrement 

        if ( iend1(l) .gt. istr1(l) ) then
          isgn1(l) =   1
        else
          isgn1(l) = - 1
        end if

!       Increment/Decrement 

        if ( iend2(l) .gt. istr2(l) ) then
          isgn2(l) =   1
        else
          isgn2(l) = - 1
        end if

      end do

!     ii1 first  boundary point of block 1
!     ii2 second boundary point of block 1
!     in1 first  interior point of block 2
!     in2 second interior point of block 2

      ic1 = cyc (idir1, 1)
      ic2 = cyc (idir1, 2)
      ic3 = cyc (idir1, 3)

      jc1 = iord (ic1)
      jc2 = iord (ic2)
      jc3 = iord (ic3)

      do i3 = 0, len(ic3)
         do i2 = 0, len(ic2)
            
            ii1 = ibcpt                        * krd (ic1, 1) + &
                 (istr1(ic2) + isgn1(ic2)*i2) * krd (ic2, 1) + &
                 (istr1(ic3) + isgn1(ic3)*i3) * krd (ic3, 1)
            jj1 = ibcpt                        * krd (ic1, 2) + &
                 (istr1(ic2) + isgn1(ic2)*i2) * krd (ic2, 2) + &
                 (istr1(ic3) + isgn1(ic3)*i3) * krd (ic3, 2)
            kk1 = ibcpt                        * krd (ic1, 3) + &
                 (istr1(ic2) + isgn1(ic2)*i2) * krd (ic2, 3) + &
                 (istr1(ic3) + isgn1(ic3)*i3) * krd (ic3, 3)
            
            in1 = inr                          * krd (jc1, 1) + &
                 (istr2(jc2) + isgn2(jc2)*i2) * krd (jc2, 1) + &
                 (istr2(jc3) + isgn2(jc3)*i3) * krd (jc3, 1)
            jn1 = inr                          * krd (jc1, 2) + &
                 (istr2(jc2) + isgn2(jc2)*i2) * krd (jc2, 2) + &
                 (istr2(jc3) + isgn2(jc3)*i3) * krd (jc3, 2)
            kn1 = inr                          * krd (jc1, 3) + &
                 (istr2(jc2) + isgn2(jc2)*i2) * krd (jc2, 3) + &
                 (istr2(jc3) + isgn2(jc3)*i3) * krd (jc3, 3)
            
            if(myblock1 .and. myblock2) then
               
               do n = 0, 2*nharms*hbmove
                  distsq = (x2(in1,jn1,kn1,n)-x1(ii1,jj1,kk1,n))**2 + &
                           (y2(in1,jn1,kn1,n)-y1(ii1,jj1,kk1,n))**2 + &
                           (z2(in1,jn1,kn1,n)-z1(ii1,jj1,kk1,n))**2
                  if (distsq.gt.mxd(n)) then
                     mxd(n)      = distsq
                     mxd_xyz(1,1,n) = x1(ii1,jj1,kk1,n) 
                     mxd_xyz(1,2,n) = y1(ii1,jj1,kk1,n)
                     mxd_xyz(1,3,n) = z1(ii1,jj1,kk1,n)
                     mxd_xyz(2,1,n) = x2(in1,jn1,kn1,n)
                     mxd_xyz(2,2,n) = y2(in1,jn1,kn1,n)
                     mxd_xyz(2,3,n) = z2(in1,jn1,kn1,n)
                     mxd_ijk(1,1,n) = ii1 
                     mxd_ijk(1,2,n) = jj1
                     mxd_ijk(1,3,n) = kk1
                     mxd_ijk(2,1,n) = in1
                     mxd_ijk(2,2,n) = jn1
                     mxd_ijk(2,3,n) = kn1
                  end if
               end do
               
            else if(myblock1) then
               
               call getownerid(iblk2,recvid)
               numrecvid = myrecvid_num(recvid)

               do inneri=1,3
                  numrecv(numrecvid) = numrecv(numrecvid)+1 

                  if (numrecv(numrecvid).gt. mrequest_perproc) then
                     write(*,*) &
                     'numrecv (check_cut)--> INCREASE MREQUEST_PERPROC', &
                     numrecv(numrecvid),mrequest_perproc
                     call abortmpi()
                  end if


                  if (numrecvid.gt. maxproc_cutman) then
                     write(*,*) &
                     'numrecv (check_cut)--> INCREASE MAXPROC_CUTMAN', &
                     numrecvid,maxproc_cutman
                     call abortmpi()
                  end if


                  receiveindices(numrecv(numrecvid),1,numrecvid) = ii1
                  receiveindices(numrecv(numrecvid),2,numrecvid) = jj1
                  receiveindices(numrecv(numrecvid),3,numrecvid) = kk1
                  receiveindices(numrecv(numrecvid),4,numrecvid) = ixyz1
                  receiveindices(numrecv(numrecvid),5,numrecvid) = imax1
                  receiveindices(numrecv(numrecvid),6,numrecvid) = jmax1
                  receiveindices(numrecv(numrecvid),7,numrecvid) = kmax1
                  receiveindices(numrecv(numrecvid),8,numrecvid) = in1
                  receiveindices(numrecv(numrecvid),9,numrecvid) = jn1
                  receiveindices(numrecv(numrecvid),10,numrecvid) = kn1
                  receiveindices(numrecv(numrecvid),11,numrecvid) = ixyz2
                  receiveindices(numrecv(numrecvid),12,numrecvid) = imax2
                  receiveindices(numrecv(numrecvid),13,numrecvid) = jmax2
                  receiveindices(numrecv(numrecvid),14,numrecvid) = kmax2
!     This entries is set to 2 if I am receiving block 2 and 1 if I am receiving block 1
                  receiveindices(numrecv(numrecvid),15,numrecvid) = 2
                  receiveindices(numrecv(numrecvid),16,numrecvid) = icut
               end do

               call getownerid(iblk2,sendid)
               numsendid = mysendid_num(sendid)
               numsend(numsendid) = numsend(numsendid)+1 

               if (numsend(numsendid).gt.mrequest_perproc) then
                  write(*,*) &
                  'numsend (check_cut) --> INCREASE MREQUEST_PERPROC', &
                  numsend(numsendid),mrequest_perproc
                  call abortmpi()
               end if

                  if (numsendid.gt. maxproc_cutman) then
                     write(*,*) &
                     'numsendid (check_cut)--> INCREASE MAXPROC_CUTMAN', &
                     numsendid,maxproc_cutman
                     call abortmpi()
                  end if


               tempindex = 1
               do n = 0, 2*nharms*hbmove
                  sendarray(tempindex,numsend(numsendid),numsendid) = &
                       x1(ii1,jj1,kk1,n)
                  tempindex = tempindex + 1
               end do
               numsend(numsendid) = numsend(numsendid)+1 
               if (numsend(numsendid).gt. mrequest_perproc) then
                  write(*,*) &
                  'numsend (check_cut) --> INCREASE MREQUEST_PERPROC', &
                  numsend(numsendid),mrequest_perproc
                  call abortmpi()
               end if


                  if (numsendid.gt. maxproc_cutman) then
                     write(*,*) &
                     'numsendid (check_cut)--> INCREASE MAXPROC_CUTMAN', &
                     numsendid,maxproc_cutman
                     call abortmpi()
                  end if


               tempindex = 1
               do n = 0, 2*nharms*hbmove
                  sendarray(tempindex,numsend(numsendid),numsendid) = &
                       y1(ii1,jj1,kk1,n)
                  tempindex = tempindex + 1
               end do
               numsend(numsendid) = numsend(numsendid)+1 
               if (numsend(numsendid).gt. mrequest_perproc) then
                  write(*,*) &
                  'numsend (check_cut) --> INCREASE MREQUEST_PERPROC', &
                  numsend(numsendid),mrequest_perproc
                  call abortmpi()
               end if


                  if (numsendid.gt. maxproc_cutman) then
                     write(*,*) &
                     'numsendid (check_cut)--> INCREASE MAXPROC_CUTMAN', &
                     numsendid,maxproc_cutman
                     call abortmpi()
                  end if



               tempindex = 1
               do n = 0, 2*nharms*hbmove
                  sendarray(tempindex,numsend(numsendid),numsendid) = &
                       z1(ii1,jj1,kk1,n)
                  tempindex = tempindex + 1
               end do
           
          
            else if (myblock2) then
               
               call getownerid(iblk1,recvid)
               numrecvid = myrecvid_num(recvid)
            
               do inneri=1,3
                  numrecv(numrecvid) = numrecv(numrecvid)+1 

                  if (numrecv(numrecvid).gt. mrequest_perproc) then
                     write(*,*) &
                     'numrecv (check_cut) --> INCREASE MREQUEST_PERPROC', &
                     numrecv(numrecvid),mrequest_perproc
                     call abortmpi()
                  end if

                  if (numrecvid.gt. maxproc_cutman) then
                     write(*,*) &
                     'numrecv (check_cut)--> INCREASE MAXPROC_CUTMAN', &
                     numrecvid,maxproc_cutman
                     call abortmpi()
                  end if


                  receiveindices(numrecv(numrecvid),1,numrecvid) = ii1
                  receiveindices(numrecv(numrecvid),2,numrecvid) = jj1
                  receiveindices(numrecv(numrecvid),3,numrecvid) = kk1
                  receiveindices(numrecv(numrecvid),4,numrecvid) = ixyz1
                  receiveindices(numrecv(numrecvid),5,numrecvid) = imax1
                  receiveindices(numrecv(numrecvid),6,numrecvid) = jmax1
                  receiveindices(numrecv(numrecvid),7,numrecvid) = kmax1
                  receiveindices(numrecv(numrecvid),8,numrecvid) = in1
                  receiveindices(numrecv(numrecvid),9,numrecvid) = jn1
                  receiveindices(numrecv(numrecvid),10,numrecvid) = kn1
                  receiveindices(numrecv(numrecvid),11,numrecvid) = ixyz2
                  receiveindices(numrecv(numrecvid),12,numrecvid) = imax2
                  receiveindices(numrecv(numrecvid),13,numrecvid) = jmax2
                  receiveindices(numrecv(numrecvid),14,numrecvid) = kmax2
!     This entries is set to 2 if I am receiving block 2 and 1 if I am receiving block 1
                  receiveindices(numrecv(numrecvid),15,numrecvid) = 1
                  receiveindices(numrecv(numrecvid),16,numrecvid) = icut
              
               end do

               call getownerid(iblk1,sendid)
               numsendid = mysendid_num(sendid)
               numsend(numsendid) = numsend(numsendid)+1 

               if (numsend(numsendid).gt. mrequest_perproc) then
                  write(*,*) &
                  'numsend (check_cut) --> INCREASE MREQUEST_PERPROC', &
                  numsend(numsendid),mrequest_perproc
                  call abortmpi()
               end if


                  if (numsendid.gt. maxproc_cutman) then
                     write(*,*) &
                     'numsendid (check_cut)--> INCREASE MAXPROC_CUTMAN', &
                     numsendid,maxproc_cutman
                     call abortmpi()
                  end if


               tempindex = 1
               do n = 0, 2*nharms*hbmove
                  sendarray(tempindex,numsend(numsendid),numsendid) = &
                       x2(in1,jn1,kn1,n)
                  tempindex = tempindex + 1
               end do
               numsend(numsendid) = numsend(numsendid)+1 
               if (numsend(numsendid).gt. mrequest_perproc) then
                  write(*,*) &
                  'numsend (check_cut) --> INCREASE MREQUEST_PERPROC', &
                  numsend(numsendid),mrequest_perproc
                  call abortmpi()
               end if


                  if (numsendid.gt. maxproc_cutman) then
                     write(*,*) &
                     'numsendid (check_cut)--> INCREASE MAXPROC_CUTMAN', &
                     numsendid,maxproc_cutman
                     call abortmpi()
                  end if


               tempindex = 1
               do n = 0, 2*nharms*hbmove
                  sendarray(tempindex,numsend(numsendid),numsendid) = &
                       y2(in1,jn1,kn1,n)
                  tempindex = tempindex + 1
               end do
               numsend(numsendid) = numsend(numsendid)+1 
               if (numsend(numsendid).gt. mrequest_perproc) then
                  write(*,*) &
                  'numsend (check_cut) --> INCREASE MREQUEST_PERPROC', &
                  numsend(numsendid),mrequest_perproc
                  call abortmpi()
               end if


                  if (numsendid.gt. maxproc_cutman) then
                     write(*,*) &
                     'numsendid (check_cut)--> INCREASE MAXPROC_CUTMAN', &
                     numsendid,maxproc_cutman
                     call abortmpi()
                  end if



               tempindex = 1
               do n = 0, 2*nharms*hbmove
                  sendarray(tempindex,numsend(numsendid),numsendid) = &
                       z2(in1,jn1,kn1,n)
                  tempindex = tempindex + 1
               end do
               
            end if
            
         end do
      end do
      
      return
      end

!-----------------------------------------------------------------------
      subroutine check_received_cut(x1,y1,z1,x2,y2,z2,nharms, &
        mxd_ijk,mxd,mxd_xyz,ii1,jj1,kk1,in1,jn1,kn1)
!-----------------------------------------------------------------------

!     Routine to do check on cuts.
!     Flow data stored in Q1 are updated from data in the interior of Q2.
!      * sendrequests/receivrequests have indexes which is the process I am communicating to
!      * sendarray/recvarray have three components:
!        - index spanning q or mut component to send/recv
!        - index of point to send/recv
!        - index of process I have to send/recv to/from
!      * receiveindices have three components:
!        - index of point to send/recv
!        - index of information to trasmit (1-6)
!        - index of process to send/recv to/from

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nharms,n,ii1,jj1,kk1,in1,jn1,kn1
      integer(kind=cosa_int) mxd_ijk(2,3,0:2*nharms*hbmove)
      real(kind=cosa_real) distsq,mxd(0:2*nharms*hbmove), &
           mxd_xyz(2,3,0:2*nharms*hbmove)
      real(kind=cosa_real) &
           x1(0:2*nharms*hbmove), &
           y1(0:2*nharms*hbmove), &
           z1(0:2*nharms*hbmove), &
           x2(0:2*nharms*hbmove), &
           y2(0:2*nharms*hbmove), &
           z2(0:2*nharms*hbmove)
!--------------------------------------------------------------------------


      
      do n = 0, 2*nharms*hbmove
         distsq = (x2(n)-x1(n))**2 + &
                  (y2(n)-y1(n))**2 + &
                  (z2(n)-z1(n))**2
         if (distsq.gt.mxd(n)) then
            mxd(n)      = distsq
            mxd_xyz(1,1,n) = x1(n) 
            mxd_xyz(1,2,n) = y1(n)
            mxd_xyz(1,3,n) = z1(n)
            mxd_xyz(2,1,n) = x2(n)
            mxd_xyz(2,2,n) = y2(n)
            mxd_xyz(2,3,n) = z2(n)
            mxd_ijk(1,1,n) = ii1 
            mxd_ijk(1,2,n) = jj1
            mxd_ijk(1,3,n) = kk1
            mxd_ijk(2,1,n) = in1
            mxd_ijk(2,2,n) = jn1
            mxd_ijk(2,3,n) = kn1
         end if
      end do
      
      
      return
      end

!-----------------------------------------------------------------------
      subroutine getblockdata(imax,jmax,kmax,x,y,z,nharms,ii1,jj1,kk1, &
                              tempdata)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,nharms,ii1,jj1,kk1,n
      real(kind=cosa_real) &
           x(0:imax+1,0:jmax+1,0:kmax+1,0:2*nharms*hbmove), &
           y(0:imax+1,0:jmax+1,0:kmax+1,0:2*nharms*hbmove), &
           z(0:imax+1,0:jmax+1,0:kmax+1,0:2*nharms*hbmove), &
           tempdata(0:2*nharms*hbmove,3)

      do n = 0, 2*nharms*hbmove
         tempdata(n,1) = x(ii1,jj1,kk1,n)
         tempdata(n,2) = y(ii1,jj1,kk1,n)
         tempdata(n,3) = z(ii1,jj1,kk1,n)
      end do

      return
      end
