
!4adrian: new routines pcutman_def and pcutman_def and cutman_def are needed.
!         their definition is sketched below (next c4adrian).
!         percutopo and cutopo will have to be passed to these new
!         routines.
!         I am not sure if it best to put pointerd here in defgrid and not having
!         them in argument of subroutine, or to use existing pointer setting in main COSA.
!-----------------------------------------------------------------------
      subroutine defgrid(nl,x0,y0,z0,workx,worky,workz,work1,dx,dy,dz,bctopo, &
        cazz,percutopo)
!-----------------------------------------------------------------------

      implicit none
      include 'cosa.inc'

      integer*4 nl,iblock,imax,jmax,kmax
      integer*4 ixyz,ixyz0,idxyz,iblk
      integer*4 bctopo(*)
      real(kind=8) x0(*),y0(*),z0(*),workx(*),worky(*),workz(*), &
        work1(*),dx(*),dy(*),dz(*)

      call zero(1,   1,0,2*nharms, 1,nl,    dx, &
                1,   1,0,2*nharms, 1)
      call zero(1,   1,0,2*nharms, 1,nl,    dy, &
                1,   1,0,2*nharms, 1)
      call zero(1,   1,0,2*nharms, 1,nl,    dzm &
                1,   1,0,2*nharms, 1)

      do iblock = 1,mynblocks
        imax  = i_imax  (iblock,nl)
        jmax  = j_jmax  (iblock,nl)
        kmax  = k_kmax  (iblock,nl)
        ixyz0 = 1 + off_p2(iblock,nl)
        idxyz = 1 + off_p1(iblock,nl)
        iblk  = 1 + off_bct(iblock,nl)
        call ass_boundaries_b(x0(ixyz0),y0(ixyz0),z0(ixyz0),dx(idxyz),dy(idxyz), &
          dz(idxyz),nbcs(iblock),bctopo(iblk),imax,jmax,kmax)
      end do

      do iter_def=1,3000

!---- warning: in HB case, work arrays have size (....0:2*nharms..).
!              however here we need only componnent 0 and NOT all
!              (2*nharms+1) components.

      call zero(1,   1,0,       0, 2,nl,workx, &
                1,   1,0,       0, 2)
      call zero(1,   1,0,       0, 2,nl,worky, &
                1,   1,0,       0, 2)
      call zero(1,   1,0,       0, 2,nl,workz, &
                1,   1,0,       0, 2)

      do iblock = 1,mynblocks
        imax  = i_imax  (iblock,nl)
        jmax  = j_jmax  (iblock,nl)
        kmax  = k_kmax  (iblock,nl)
        ixyz  = 1 + off_p3(iblock,nl)
        ixyz0 = 1 + off_p2(iblock,nl)
        idxyz = 1 + off_p1(iblock,nl)
        iblk  = 1 + off_bct(iblock,nl)
        call dxyz_increment_b(x0(ixyz0),y0(ixyz0),z0(ixyz0),dx(idxyz), &
          dy(idxyz),dz(idxyz),workx(ixyz0),worky(ixyz0),workz(ixyz0), &
          work1(ixyz),imax,jmax,kmax)
      end do

!4adrian:
!def- routines cutman_def and pcutman_def can be written starting from
!       cutman_x and pcutman_x respectively. They differ from the
!       indicated source for 2 reasons:
!         1) routines *_def act only on component 0 of work arrays
!            (if they act on components 0:2*nharms it should only be
!            waste of computing.
!         2) equality is exnforced exactly on cut, and not one or two
!            interior grid verteices.

      call  cutman_def(nl,   cutopo,workx,worky,workz)
      call pcutman_def(nl,percutopo,workx,worky,workz)

      do iblock = 1,mynblocks
        imax  = i_imax  (iblock,nl)
        jmax  = j_jmax  (iblock,nl)
        kmax  = k_kmax  (iblock,nl)
        ixyz  = 1 + off_p3(iblock,nl)
        idxyz = 1 + off_p1(iblock,nl)
        ixyz0 = 1 + off_p2(iblock,nl)
        call zero_boundaries_b(workx(ixyz0),worky(ixyz0), &
          workz(ixyz0),nbcs(iblock),bctopo(iblk),imax,jmax,kmax)
        call update_dxyz(dx(idxyz),dy(idxyz),dz(idxyz), &
          workx(ixyz0),worky(ixyz0),workz(ixyz0),work1(ixyz), &
          imax,jmax,kmax)
      end do

      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine dxyz_increment_b(x0,y0,z0,dx,dy,dz,workx,worky,workz,work1, &
        imax,jmax,kmax)
!-----------------------------------------------------------------------

      implicit none

      include 'common.block'

      integer*4 imax,jmax,kmax,i,j,k
      real(kind=8) &
        x0    ( 0:imax+1, 0:jmax+1, 0:kmax+1), &
        y0    ( 0:imax+1, 0:jmax+1, 0:kmax+1), &
        z0    ( 0:imax+1, 0:jmax+1, 0:kmax+1), &
        workx ( 0:imax+1, 0:jmax+1, 0:kmax+1), &
        worky ( 0:imax+1, 0:jmax+1, 0:kmax+1), &
        workz ( 0:imax+1, 0:jmax+1, 0:kmax+1), &
        work1 (-1:imax+1,-1:jmax+1, -1:kmax+1), &
        dx    ( 0:imax  , 0:jmax  , 0:kmax  ), &
        dy    ( 0:imax  , 0:jmax  , 0:kmax  ), &
        dz    ( 0:imax  , 0:jmax  , 0:kmax  )
      real(kind=8) rds

!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------

!def    WARNING: implementation below is done for
!         for 2D pitching case only if z coordinate always has index k.
!         It may work in more general cases of axis orientation, I need 
!         to check. However this is fine, in theory, for test case
!         we will use now.
!         
!------ j=const. edges
        do k = 1, kmax
        do j = 1, jmax
          do i = 1, imax-1

            rds = ( (x0(i+1,j,k) - x0(i,j,k))**2 + &
                    (y0(i+1,j,k) - y0(i,j,k))**2 + &
                    (z0(i+1,j,k) - z0(i,j,k))**2 )**(-0.5)

            workx(i+1,j,k) = workx(i+1,j,k) + &
                             rds*(dx(i,j,k)-dx(i+1,j,k))
            workx(i  ,j,k) = workx(i  ,j,k) - &
                             rds*(dx(i,j,k)-dx(i+1,j,k))
            worky(i+1,j,k) = worky(i+1,j,k) + &
                             rds*(dy(i,j,k)-dy(i+1,j,k))
            worky(i  ,j,k) = worky(i  ,j,k) - &
                             rds*(dy(i,j,k)-dy(i+1,j,k))
            workz(i+1,j,k) = workz(i+1,j,k) + &
                             rds*(dz(i,j,k)-dz(i+1,j,k))
            workz(i  ,j,k) = workz(i  ,j,k) - &
                             rds*(dz(i,j,k)-dz(i+1,j,k))

            work1(i+1,j,k) = work1(i+1,j,k) + rds
            work1(i  ,j,k) = work1(i  ,j,k) + rds
          end do
        end do
        end do

!msc--- WARNING: in loop below I swapped order of i and j loops w.r.t.
!                original DEFGRID. I think it should be mathematically 
!                equivalent.

!------ j=const. edges
        do k = 1, kmax
        do j = 1, jmax-1
          do i = 1, imax

            rds = ( (x0(i,j,k) - x0(i,j+1,k))**2 + &
                    (y0(i,j,k) - y0(i,j+1,k))**2 + &
                    (z0(i,j,k) - z0(i,j+1,k))**2 )**(-0.5)

            workx(i,j  ,k) = workx(i,j  ,k) + &
                             rds*(dx(i,j+1,k)-dx(i,j,k))
            workx(i,j+1,k) = workx(i,j+1,k) - &
                             rds*(dx(i,j+1,k)-dx(i,j,k))
            worky(i,j  ,k) = worky(i,j  ,k) + &
                             rds*(dy(i,j+1,k)-dy(i,j,k))
            worky(i,j+1,k) = worky(i,j+1,k) - &
                             rds*(dy(i,j+1,k)-dy(i,j,k))
            workz(i,j  ,k) = workz(i,j  ,k) + &
                             rds*(dz(i,j+1,k)-dz(i,j,k))
            workz(i,j+1,k) = workz(i,j+1,k) - &
                             rds*(dz(i,j+1,k)-dz(,i,j,k))

            work1(i,j  ,k) = work1(i,j  ,k) + rds
            work1(i,j+1,k) = work1(i,j+1,k) + rds
          end do
        end do

!------ z=const. edges
        do k = 1, kmax-1
        do j = 1, jmax
          do i = 1, imax

            rds = ( (x0(i,j,k) - x0(i,j,k+1))**2 + &
                    (y0(i,j,k) - y0(i,j,k+1))**2 + &
                    (z0(i,j,k) - z0(i,j,k+1))**2 )**(-0.5)

            workx(i,j,k  ) = workx(i,j,l  ) + &
                             rds*(dx(i,j,k+1)-dx(i,j,k))
            workx(i,j,k+1) = workx(i,j,k+1) - &
                             rds*(dx(i,j,k+1)-dx(i,j,k))

            worky(i,j,k  ) = worky(i,j,l  ) + &
                             rds*(dy(i,j,k+1)-dy(i,j,k))
            worky(i,j,k+1) = worky(i,j,k+1) - &
                             rds*(dy(i,j,k+1)-dy(i,j,k))

            workz(i,j,k  ) = workz(i,j,l  ) + &
                             rds*(dz(i,j,k+1)-dz(i,j,k))
            workz(i,j,k+1) = workz(i,j,k+1) - &
                             rds*(dz(i,j,k+1)-dz(i,j,k))

            work1(i,j,k  ) = work1(i,j,k  ) + rds
            work1(i,j,k+1) = work1(i,j,k+1) + rds
          end do
        end do

        return
        end

!-----------------------------------------------------------------------
      subroutine update_dxyz(dx,dy,dz,workx,worky,workz,work1, &
        imax,jmax,kmax)
!-----------------------------------------------------------------------

      implicit none

      include 'common.block'

      integer*4 imax,jmax,kmax,i,j,k
      real(kind=8) &
        workx ( 0:imax+1, 0:jmax+1, 0:kmax+1), &
        worky ( 0:imax+1, 0:jmax+1, 0:kmax+1), &
        workz ( 0:imax+1, 0:jmax+1, 0:kmax+1), &
        work1 (-1:imax+1,-1:jmax+1,-1:kmax+1), &
        dx    ( 0:imax  , 0:jmax  , 0:kmax  ), &
        dy    ( 0:imax  , 0:jmax  , 0:kmax  ), &
        dz (0:imax  ,0:jmax  ,0:kmax  )
      real(kind=8) fac

!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------

!...... add correction and reset work array
        do k = 1,kmax
        do j = 1,jmax
          do i = 1,imax

            fac = 0.5d0 / work1(i,j,k)

!def------- update array of nodal displacements
            do m = 1, 2
              dx(i,j,k) = dx(i,j,k) + fac*workx(i,j,k)
              dy(i,j,k) = dy(i,j,k) + fac*worky(i,j,k)
              dz(i,j,k) = dz(i,j,k) + fac*workz(i,j,k)
            end do

          end do
        enddo
        enddo

      return
      end

!-----------------------------------------------------------------------
      subroutine ass_boundaries_b(x0,y0,z0,dx,dy,dz,bctopo,nbcs,imax,jmax,kmax)
!-----------------------------------------------------------------------

      implicit none

      include 'common.block'

      integer*4 imax,jmax,kmax,nbcs
      integer*4 nl,ibc
      integer*4 bctopo(10,nbcs)
      real(kind=8) &
        x0 (0:imax+1,0:jmax+1,0:kmax+1), &
        y0 (0:imax+1,0:jmax+1,0:kmax+1), &
        z0 (0:imax+1,0:jmax+1,0:kmax+1), &
        dx (0:imax  ,0:jmax  ,0:kmax  ), &
        dy (0:imax  ,0:jmax  ,0:kmax  ), &
        dz (0:imax  ,0:jmax  ,0:kmax  )

!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
      do ibc = 1,nbcs

!------ viscous wall (jbctyp=15)
        if (bctopo(1,ibc)/100.eq.15) then

          call set_wall_disp(x0,y0,z0,dx,gy,dz,bctopo(1,ibc),imax,jmax, &
                             kmax)

        else

          call set_zero_disp(x0,y0,z0,dx,gy,dz,bctopo(1,ibc),imax,jmax, &
                             kmax)
        end if

      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine set_wall_disp(x0,y0,z0,dx,gy,dz,bctopo,imax,jmax,kmax)
!-----------------------------------------------------------------------

      implicit none

      include 'common.block'

      integer*4 imax,jmax,kmax,nharms
      integer*4 bctopo(10)
      integer*4 n,nh,i1,i2,i3,ijkmax(3),idir,istrt(3),iend(3),bctyp, &
                inrout,ibcpt,ibcpt2,ibcn,ibcn2,ibcm,ibc,jbc,kbc,ibc2, &
                jbc2,kbc2,in,jn,kn,in2,jn2,kn2,im,jm,km,ic1,ic2,ic3
      real (kind=8) &
        x0 (0:imax+1,0:jmax+1,0:kmax+1), &
        y0 (0:imax+1,0:jmax+1,0:kmax+1), &
        z0 (0:imax+1,0:jmax+1,0:kmax+1), &
        dx (0:imax  ,0:jmax  ,0:kmax  ), &
        dy (0:imax  ,0:jmax  ,0:kmax  ), &
        dz (0:imax  ,0:jmax  ,0:kmax  )
      real (kind=8) rr11,rr12,rr21,rr22.xtemp.ytemp

!     store imax, jmax in ijmax to enable boundary data location

      ijkmax(1) = imax
      ijkmax(2) = jmax
      ijkmax(3) = kmax

!     store boundary topology in mnemonic names

      bctyp    = bctopo(1)
      idir     = bctopo(2)
      inrout   = bctopo(3)
      istrt(1) = bctopo(4)
      iend(1)  = bctopo(5)
      istrt(2) = bctopo(6)
      iend(2)  = bctopo(7)
      istrt(3) = bctopo(8)
      iend(3)  = bctopo(9)

      rr11 =  dcos(dtheta0)
      rr12 = -dsin(dtheta0)
      rr21 =  dsin(dtheta0)
      rr22 =  dcos(dtheta0)

!     modify beginning, ending indices to extend boundary condition to
!     edge/corner

!tmp  do i1 = 1,2
!tmp    if (i1.ne.idir) then
!tmp      if (istrt(i1) .eq. 1          ) istrt(i1) = 0
!tmp      if (iend (i1) .eq. ijmax(i1)-1) iend (i1) = ijmax(i1)
!tmp    end if
!tmp  end do

!     set needed variables depending on whether the boundary is the
!     the inner boundary (inrout = 1) or the outer boundary (inrout > 1)
!         ibcpt : boundary condition location (first aux. cell)
!         ibcpt2: boundary condition location outside the block from
!                 ibcpt (second aux. cell)
!         ibcn  : point to the inside of the block from ibcpt
!         ibcm  : location of the metrics

      if (inrout.eq.1) then
!del    ibcpt  =  0
!del    ibcpt2 = -1
!del    ibcn   =  1
!del    ibcn2  =  2
        ibcm   =  1
!del    sgnm   = - 1.d0
      else
!del    ibcpt  = ijkmax(idir)
!del    ibcpt2 = ijkmax(idir) + 1
!del    ibcn   = ijkmax(idir) - 1
!del    ibcn2  = ijkmax(idir) - 2
        ibcm   = ijkmax(idir)
!del    sgnm   = 1.d0
      end if

      ic1 = cyc (idir, 1)
      ic2 = cyc (idir, 2)
      ic3 = cyc (idir, 3)

        do i3 = istrt(ic3),iend(ic3)
          do i2 = istrt(ic2),iend(ic2)

!del        ibc  = ibcpt *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
!del        jbc  = ibcpt *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
!del        kbc  = ibcpt *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

!del        ibc2 = ibcpt2*krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
!del        jbc2 = ibcpt2*krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
!del        kbc2 = ibcpt2*krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

!del        in   = ibcn  *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
!del        jn   = ibcn  *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
!del        kn   = ibcn  *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

!del        in2  = ibcn2 *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
!del        jn2  = ibcn2 *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
!del        kn2  = ibcn2 *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

            im   = ibcm  *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
            jm   = ibcm  *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
            km   = ibcm  *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

            xtemp = x0(im,jm,km) - xh
            ytemp = y0(im,jm,km) - yh
            dx(im,jm,km) = rr11*xtemp + rr12*ytemp
            dy(im,jm,km) = rr21*xtemp + rr22*ytemp
            dx(im,jm,km) = dx(im,jm,km) - xtemp
            dy(im,jm,km) = dy(im,jm,km) - ytemp

          end do
        end do

      return
      end

!-----------------------------------------------------------------------
      subroutine set_zero_disp(x0,y0,z0,dx,gy,dz,bctopo,imax,jmax,kmax)
!-----------------------------------------------------------------------

      implicit none

      include 'common.block'

      integer*4 imax,jmax,kmax,nharms
      integer*4 bctopo(10)
      integer*4 n,nh,i1,i2,i3,ijkmax(3),idir,istrt(3),iend(3),bctyp, &
                inrout,ibcpt,ibcpt2,ibcn,ibcn2,ibcm,ibc,jbc,kbc,ibc2, &
                jbc2,kbc2,in,jn,kn,in2,jn2,kn2,im,jm,km,ic1,ic2,ic3
      real (kind=8) &
        x0 (0:imax+1,0:jmax+1,0:kmax+1), &
        y0 (0:imax+1,0:jmax+1,0:kmax+1), &
        z0 (0:imax+1,0:jmax+1,0:kmax+1), &
        dx (0:imax  ,0:jmax  ,0:kmax  ), &
        dy (0:imax  ,0:jmax  ,0:kmax  ), &
        dz (0:imax  ,0:jmax  ,0:kmax  )
      real (kind=8) rr11,rr12,rr21,rr22.xtemp.ytemp

!     store imax, jmax in ijmax to enable boundary data location

      ijkmax(1) = imax
      ijkmax(2) = jmax
      ijkmax(3) = kmax

!     store boundary topology in mnemonic names

      bctyp    = bctopo(1)
      idir     = bctopo(2)
      inrout   = bctopo(3)
      istrt(1) = bctopo(4)
      iend(1)  = bctopo(5)
      istrt(2) = bctopo(6)
      iend(2)  = bctopo(7)
      istrt(3) = bctopo(8)
      iend(3)  = bctopo(9)

      rr11 =  dcos(dtheta0)
      rr12 = -dsin(dtheta0)
      rr21 =  dsin(dtheta0)
      rr22 =  dcos(dtheta0)

!     modify beginning, ending indices to extend boundary condition to
!     edge/corner

!tmp  do i1 = 1,2
!tmp    if (i1.ne.idir) then
!tmp      if (istrt(i1) .eq. 1          ) istrt(i1) = 0
!tmp      if (iend (i1) .eq. ijmax(i1)-1) iend (i1) = ijmax(i1)
!tmp    end if
!tmp  end do

!     set needed variables depending on whether the boundary is the
!     the inner boundary (inrout = 1) or the outer boundary (inrout > 1)
!         ibcpt : boundary condition location (first aux. cell)
!         ibcpt2: boundary condition location outside the block from
!                 ibcpt (second aux. cell)
!         ibcn  : point to the inside of the block from ibcpt
!         ibcm  : location of the metrics

      if (inrout.eq.1) then
!del    ibcpt  =  0
!del    ibcpt2 = -1
!del    ibcn   =  1
!del    ibcn2  =  2
        ibcm   =  1
!del    sgnm   = - 1.d0
      else
!del    ibcpt  = ijkmax(idir)
!del    ibcpt2 = ijkmax(idir) + 1
!del    ibcn   = ijkmax(idir) - 1
!del    ibcn2  = ijkmax(idir) - 2
        ibcm   = ijkmax(idir)
!del    sgnm   = 1.d0
      end if

      ic1 = cyc (idir, 1)
      ic2 = cyc (idir, 2)
      ic3 = cyc (idir, 3)

        do i3 = istrt(ic3),iend(ic3)
          do i2 = istrt(ic2),iend(ic2)

!del        ibc  = ibcpt *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
!del        jbc  = ibcpt *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
!del        kbc  = ibcpt *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

!del        ibc2 = ibcpt2*krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
!del        jbc2 = ibcpt2*krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
!del        kbc2 = ibcpt2*krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

!del        in   = ibcn  *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
!del        jn   = ibcn  *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
!del        kn   = ibcn  *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

!del        in2  = ibcn2 *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
!del        jn2  = ibcn2 *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
!del        kn2  = ibcn2 *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

            im   = ibcm  *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
            jm   = ibcm  *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
            km   = ibcm  *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

            dx(im,jm,km) = 0.d0
            dy(im,jm,km) = 0.d0

          end do
        end do

      return
      end

!-----------------------------------------------------------------------
      subroutine zero_boundaries_b(workx,worky,workz,bctopo,nbcs,imax,jmax,kmax)
!-----------------------------------------------------------------------

      implicit none

      include 'common.block'

      integer*4 imax,jmax,kmax,nbcs
      integer*4 nl,ibc
      integer*4 bctopo(10,nbcs)
      real(kind=8) &
        workx (0:imax+1,0:jmax+1,0:kmax+1), &
        worky (0:imax+1,0:jmax+1,0:kmax+1), &
        workz (0:imax+1,0:jmax+1,0:kmax+1)

!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
      do ibc = 1,nbcs

          call set_zero_disp1(workx,worky,workz,bctopo(1,ibc),imax,jmax, &
                             kmax)

      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine set_zero_disp1(workx,worky,workz,bctopo,imax,jmax,kmax)
!-----------------------------------------------------------------------

      implicit none

      include 'common.block'

      integer*4 imax,jmax,kmax,nharms
      integer*4 bctopo(10)
      integer*4 n,nh,i1,i2,i3,ijkmax(3),idir,istrt(3),iend(3),bctyp, &
                inrout,ibcpt,ibcpt2,ibcn,ibcn2,ibcm,ibc,jbc,kbc,ibc2, &
                jbc2,kbc2,in,jn,kn,in2,jn2,kn2,im,jm,km,ic1,ic2,ic3
      real (kind=8) &
        workx (0:imax+1,0:jmax+1,0:kmax+1), &
        worky (0:imax+1,0:jmax+1,0:kmax+1), &
        workz (0:imax+1,0:jmax+1,0:kmax+1)
      real (kind=8) rr11,rr12,rr21,rr22.xtemp.ytemp

!     store imax, jmax in ijmax to enable boundary data location

      ijkmax(1) = imax
      ijkmax(2) = jmax
      ijkmax(3) = kmax

!     store boundary topology in mnemonic names

      bctyp    = bctopo(1)
      idir     = bctopo(2)
      inrout   = bctopo(3)
      istrt(1) = bctopo(4)
      iend(1)  = bctopo(5)
      istrt(2) = bctopo(6)
      iend(2)  = bctopo(7)
      istrt(3) = bctopo(8)
      iend(3)  = bctopo(9)

      rr11 =  dcos(dtheta0)
      rr12 = -dsin(dtheta0)
      rr21 =  dsin(dtheta0)
      rr22 =  dcos(dtheta0)

!     modify beginning, ending indices to extend boundary condition to
!     edge/corner

!tmp  do i1 = 1,2
!tmp    if (i1.ne.idir) then
!tmp      if (istrt(i1) .eq. 1          ) istrt(i1) = 0
!tmp      if (iend (i1) .eq. ijmax(i1)-1) iend (i1) = ijmax(i1)
!tmp    end if
!tmp  end do

!     set needed variables depending on whether the boundary is the
!     the inner boundary (inrout = 1) or the outer boundary (inrout > 1)
!         ibcpt : boundary condition location (first aux. cell)
!         ibcpt2: boundary condition location outside the block from
!                 ibcpt (second aux. cell)
!         ibcn  : point to the inside of the block from ibcpt
!         ibcm  : location of the metrics

      if (inrout.eq.1) then
!del    ibcpt  =  0
!del    ibcpt2 = -1
!del    ibcn   =  1
!del    ibcn2  =  2
        ibcm   =  1
!del    sgnm   = - 1.d0
      else
!del    ibcpt  = ijkmax(idir)
!del    ibcpt2 = ijkmax(idir) + 1
!del    ibcn   = ijkmax(idir) - 1
!del    ibcn2  = ijkmax(idir) - 2
        ibcm   = ijkmax(idir)
!del    sgnm   = 1.d0
      end if

      ic1 = cyc (idir, 1)
      ic2 = cyc (idir, 2)
      ic3 = cyc (idir, 3)

        do i3 = istrt(ic3),iend(ic3)
          do i2 = istrt(ic2),iend(ic2)

!del        ibc  = ibcpt *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
!del        jbc  = ibcpt *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
!del        kbc  = ibcpt *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

!del        ibc2 = ibcpt2*krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
!del        jbc2 = ibcpt2*krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
!del        kbc2 = ibcpt2*krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

!del        in   = ibcn  *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
!del        jn   = ibcn  *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
!del        kn   = ibcn  *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

!del        in2  = ibcn2 *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
!del        jn2  = ibcn2 *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
!del        kn2  = ibcn2 *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

            im   = ibcm  *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
            jm   = ibcm  *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
            km   = ibcm  *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

            workx(im,jm,km) = 0
            worky(im,jm,km) = 0
            workz(im,jm,km) = 0

          end do
        end do

      return
      end
