
      module postprocessing

      use cosa_precision

      implicit none

      private

      save

      contains
      

!-----------------------------------------------------------------------
      subroutine sfht(nl,itime1,ijkmaxg,nsurf,msurf,surf_dim,cf,yplus, &
                      xc,yc,zc,xyzcdot,dist,si,sj,sk,xideri,xiderj, &
                      xiderk,etaderi,etaderj,etaderk,zetaderi,zetaderj, &
                      zetaderk,q1,q2,bctopo)
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,itime1,iblock,iblk_pos,imax,jmax,kmax,ijkmaxg,ixyzc, &
                ixyzc1,idist,iimt,ivmt,iq,iblk,nsurf,msurf, &
                surf_dim(msurf),nsurf1
      real (kind=cosa_real) xc(*),yc(*),zc(*),xyzcdot(*),dist(*),si(*),sj(*), &
           sk(*),xideri(*),xiderj(*),xiderk(*),etaderi(*),etaderj(*), &
           etaderk(*),zetaderi(*),zetaderj(*),zetaderk(*),q1(*),q2(*), &
           cf(ijkmaxg-1,msurf),yplus(ijkmaxg-1,msurf)
      integer(kind=cosa_int) bctopo(*)

      nsurf1 = 0
      do iblk_pos = 1,nblock_pos
        
        iblock   = iblock_pos(iblk_pos)
        imax     = i_imax(iblock,nl)
        jmax     = j_jmax(iblock,nl)
        kmax     = k_kmax(iblock,nl)
        ixyzc    = 1 + off_p2 (iblock,nl)
        ixyzc1   = 1 + off_p3 (iblock,nl) * npde * dim5
        idist    = 1 + off_p1 (iblock,nl)
        iimt     = 1 + off_p2 (iblock,nl) * lmet * dim5h
        ivmt     = 1 + off_0  (iblock,nl) * 3    * dim5h
        iq       = 1 + off_p3 (iblock,nl) * npde
        iblk     = 1 + off_bct(iblock,nl)

        call sfht_b(imax,jmax,kmax,ijkmaxg,npde,lmet,nharms,itime1, &
                    nsurf1,msurf,surf_dim,cf,yplus,nbcs(iblock), &
                    bctopo(iblk),xc(ixyzc),yc(ixyzc),zc(ixyzc), &
                    xyzcdot(ixyzc1),dist(idist),si(iimt),sj(iimt), &
                    sk(iimt),xideri(ivmt),xiderj(ivmt),xiderk(ivmt), &
                    etaderi(ivmt),etaderj(ivmt),etaderk(ivmt), &
                    zetaderi(ivmt),zetaderj(ivmt),zetaderk(ivmt),q1(iq), &
                    q2(iq))

      end do

      if (nsurf1.ne.nsurf) then
        write(*,*) 'Incongruency in sfht. Aborting!'
        write(*,*) 'nsurf1 = ', nsurf1
        write(*,*) 'nsurf  = ', nsurf
        stop
      end if

      return
      end

!-----------------------------------------------------------------------
      subroutine sfht_b(imax,jmax,kmax,ijkmaxg,npde,lmet,nharms,itime1, &
                        nsurf1,msurf,surf_dim,cf,yplus,nbcs,bctopo,xc, &
                        yc,zc,xyzcdot,dist,si,sj,sk,xideri,xiderj, &
                        xiderk,etaderi,etaderj,etaderk,zetaderi, &
                        zetaderj,zetaderk,q1,q2)
!-----------------------------------------------------------------------
!     Routine to output the skin friction, heat transfer, and boundary 
!     layer quantities.
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,ijkmaxg,npde,nharms,nbcs,lmet
      integer(kind=cosa_int) bctopo(10,nbcs)
      integer(kind=cosa_int) itime1,n,nh,msurf,surf_dim(msurf),nsurf1
      real (kind=cosa_real) &
           xc(0:imax+1,0:jmax+1,0:kmax+1), &
           yc(0:imax+1,0:jmax+1,0:kmax+1), &
           zc(0:imax+1,0:jmax+1,0:kmax+1), &
           xyzcdot (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms    ), &
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
           q1(-1:imax+1,-1:jmax+1,-1:kmax+1,npde), &
           q2(-1:imax+1,-1:jmax+1,-1:kmax+1,npde), &
           dist(0:imax,0:jmax,0:kmax)
      real (kind=cosa_real) cf(ijkmaxg-1,msurf),yplus(ijkmaxg-1,msurf)

      if ((.not.unsteady).or.dualt) then
        nh = 0
      else if (harbal) then
        n  = itime1 - 1
        nh = n*hbmove
      end if
      call sfht_b_h(imax,jmax,kmax,ijkmaxg,npde,lmet,msurf,surf_dim, &
                    nsurf1,cf,yplus,xc,yc,zc,xyzcdot(-1,-1,-1,1,nh), &
                    dist,si(1,0,0,0,nh),sj(1,0,0,0,nh),sk(1,0,0,0,nh), &
                    xideri(1,1,1,1,nh),xiderj(1,1,1,1,nh), &
                    xiderk(1,1,1,1,nh),etaderi(1,1,1,1,nh), &
                    etaderj(1,1,1,1,nh),etaderk(1,1,1,1,nh), &
                    zetaderi(1,1,1,1,nh),zetaderj(1,1,1,1,nh), &
                    zetaderk(1,1,1,1,nh),q1,q2,bctopo,nbcs)

      return
      end

!-----------------------------------------------------------------------
      subroutine sfht_b_h(imax,jmax,kmax,ijkmaxg,npde,lmet,msurf, &
                          surf_dim,nsurf1,cf,yplus,xc,yc,zc,xyzcdot, &
                          dist,si,sj,sk,xideri,xiderj,xiderk,etaderi, &
                          etaderj,etaderk,zetaderi,zetaderj,zetaderk,q1, &
                          q2,bctopo,nbcs)
!-----------------------------------------------------------------------
!     Routine to output the skin friction, heat transfer, and boundary 
!     layer quantities.
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax, jmax, kmax, npde, nbcs, lmet
      integer(kind=cosa_int) ijkmaxg,msurf,surf_dim(msurf),nsurf1,sur_dim,i1,i2,i3, &
                ijkmax(3),idir,istrt(3),iend(3),bctyp,inrout, &
                ibcpt,ibcn,ibcn2,ibcm,ibc,jbc,kbc,in,jn,kn,in2,jn2,kn2, &
                im,jm,km,ic1,ic2,ic3,ibbc,ioff1,ioff2,joff1,joff2,koff1, &
                koff2
      real (kind=cosa_real) &
           xc(0:imax+1,0:jmax+1,0:kmax+1), &
           yc(0:imax+1,0:jmax+1,0:kmax+1), &
           zc(0:imax+1,0:jmax+1,0:kmax+1), &
           xyzcdot (-1:imax+1,-1:jmax+1,-1:kmax+1,npde), &
           si      (lmet,0:imax+1,0:jmax+1,0:kmax+1), &
           sj      (lmet,0:imax+1,0:jmax+1,0:kmax+1), &
           sk      (lmet,0:imax+1,0:jmax+1,0:kmax+1), &
           xideri  (3   ,imax    ,jmax    ,kmax    ), &
           xiderj  (3   ,imax    ,jmax    ,kmax    ), &
           xiderk  (3   ,imax    ,jmax    ,kmax    ), &
           etaderi (3   ,imax    ,jmax    ,kmax    ), &
           etaderj (3   ,imax    ,jmax    ,kmax    ), &
           etaderk (3   ,imax    ,jmax    ,kmax    ), &
           zetaderi(3   ,imax    ,jmax    ,kmax    ), &
           zetaderj(3   ,imax    ,jmax    ,kmax    ), &
           zetaderk(3   ,imax    ,jmax    ,kmax    ), &
           q1(-1:imax+1,-1:jmax+1,-1:kmax+1,npde), &
           q2(-1:imax+1,-1:jmax+1,-1:kmax+1,npde), &
           dist(0:imax,0:jmax,0:kmax)
      real (kind=cosa_real) cf(ijkmaxg-1,msurf), yplus(ijkmaxg-1,msurf), sc, &
           muw,rho1, rho2, rhow, p1, &
           p2, pw, t1, t2, tw, nx, ny, nz, ds, kx, ky, kz, uc, vc, wc, &
           dn1, tauw, tauwx, tauwy, tauwz, tauwpx, tauwpy, tauwpz, &
           u1, u2, dudxi, dudeta, dudzeta, v1, v2, dvdxi, dvdeta, &
           dvdzeta, w1, w2, dwdxi, dwdeta, dwdzeta, dudx, dudy, dudz, &
           dvdx, dvdy, dvdz, dwdx, dwdy, dwdz, xix, xiy, xiz, etax, &
           etay, etaz, zetax, zetay, zetaz, divv, txx, txy, txz, tyy, &
           tyz, tzz, utau
      integer(kind=cosa_int) bctopo(10,nbcs)

      do ibbc=1,nbcs

        if (bctopo(1,ibbc)/100.eq.wallmask) then

          sur_dim = 0
          nsurf1 = nsurf1 + 1

!         store imax, jmax in ijmax to enable boundary data location

          ijkmax(1) = imax
          ijkmax(2) = jmax
          ijkmax(3) = kmax

!         store boundary topology in mnemonic names

          bctyp    = bctopo(1,ibbc)
          idir     = bctopo(2,ibbc)
          inrout   = bctopo(3,ibbc)
          istrt(1) = bctopo(4,ibbc)
          iend(1)  = bctopo(5,ibbc)
          istrt(2) = bctopo(6,ibbc)
          iend(2)  = bctopo(7,ibbc)
          istrt(3) = bctopo(8,ibbc)
          iend(3)  = bctopo(9,ibbc)


!          print*,bctyp,idir,inrout, istrt(1),iend(1),istrt(2),iend(2)
!         modify beginning, ending indices to extend boundary condition
!         to edge/corner

!tmp      do i1 = 1,2
!tmp        if (i1.ne.idir) then
!tmp          if (istrt(i1) .eq. 1          ) istrt(i1) = 0
!tmp          if (iend (i1) .eq. ijmax(i1)-1) iend (i1) = ijmax(i1)
!tmp        end if
!tmp      end do

!         set needed variables depending on whether the boundary is the
!         the inner boundary (inrout = 1) or 
!         the outer boundary (inrout > 1)
!             ibcpt : boundary condition location (first aux. cell)
!             ibcpt2: boundary condition location outside the block from
!                     ibcpt (second aux. cell)
!             ibcn  : point to the inside of the block from ibcpt
!             ibcm  : location of the metrics

          if (inrout.eq.1) then
            ibcpt  =  0
            ibcn   =  1
            ibcn2  =  2
            ibcm   =  1
          else
            ibcpt  = ijkmax(idir)
            ibcn   = ijkmax(idir) - 1
            ibcn2  = ijkmax(idir) - 2
            ibcm   = ijkmax(idir)
          end if

          ic1 = cyc (idir, 1)
          ic2 = cyc (idir, 2)
          ic3 = cyc (idir, 3)

          do i3 = istrt(ic3),iend(ic3)
            do i2 = istrt(ic2),iend(ic2)
      
              sur_dim = sur_dim + 1

              ibc  = ibcpt *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
              jbc  = ibcpt *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
              kbc  = ibcpt *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

              in   = ibcn  *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
              jn   = ibcn  *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
              kn   = ibcn  *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

              in2  = ibcn2 *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
              jn2  = ibcn2 *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
              kn2  = ibcn2 *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

              im   = ibcm  *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
              jm   = ibcm  *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
              km   = ibcm  *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

              if (idir.eq.1) then
                nx    = si(1,im,jm,km)
                ny    = si(2,im,jm,km)
                nz    = si(3,im,jm,km)
                ds    = si(4,im,jm,km)
              else if (idir.eq.2) then
                nx    = sj(1,im,jm,km)
                ny    = sj(2,im,jm,km)
                nz    = sj(3,im,jm,km)
                ds    = sj(4,im,jm,km)
              else if (idir.eq.3) then
                nx    = sk(1,im,jm,km)
                ny    = sk(2,im,jm,km)
                nz    = sk(3,im,jm,km)
                ds    = sk(4,im,jm,km)
              end if

              kx   = (-1)**(1+1/inrout) * nx
              ky   = (-1)**(1+1/inrout) * ny
              kz   = (-1)**(1+1/inrout) * nz

              rhow = q1(ibc,jbc,kbc,1)
              rho1 = q1(in ,jn ,kn ,1)
              rho2 = q1(in2,jn2,kn2,1)

              pw   = q1(ibc,jbc,kbc,5)
              p1   = q1(in ,jn ,kn ,5)
              p2   = q1(in2,jn2,kn2,5)

              tw   = q2(ibc,jbc,kbc,1)
              t1   = q2(in ,jn ,kn ,1)
              t2   = q2(in2,jn2,kn2,1)

              muw  = (stemp+1.d0)/(tw+stemp) * (tw**1.5d0)

              dn1     = dist(in ,jn ,kn )

              if (idir.eq.1) then
                ioff1 = (-1)**(1+1/inrout)
                ioff2 = 2*ioff1
                joff1 = 0
                joff2 = 0
                koff1 = 0
                koff2 = 0
              else if (idir.eq.2) then
                ioff1 = 0
                ioff2 = 0
                joff1 = (-1)**(1+1/inrout)
                joff2 = 2*joff1
                koff1 = 0
                koff2 = 0
              else if (idir.eq.3) then
                ioff1 = 0
                ioff2 = 0
                joff1 = 0
                joff2 = 0
                koff1 = (-1)**(1+1/inrout)
                koff2 = 2*koff1
              end if
 
              if (idir.eq.1) then
   
                dudxi = (-1)**(1+1/inrout) * &
                        (-8*q1(ibc,jbc,kbc,2) + &
                          9*q1(ibc+ioff1,jbc,kbc,2) - &
                            q1(ibc+ioff2,jbc,kbc,2) ) / 3
                dvdxi = (-1)**(1+1/inrout) * &
                        (-8*q1(ibc,jbc,kbc,3) + &
                          9*q1(ibc+ioff1,jbc,kbc,3) - &
                            q1(ibc+ioff2,jbc,kbc,3) ) / 3
                dwdxi = (-1)**(1+1/inrout) * &
                        (-8*q1(ibc,jbc,kbc,4) + &
                          9*q1(ibc+ioff1,jbc,kbc,4) - &
                            q1(ibc+ioff2,jbc,kbc,4) ) / 3

                if (i2.eq.istrt(ic2)) then
                  u2 = (  q1(ibc,jbc,kbc,2) + q1(ibc,jbc+1,kbc,2)) / 2
                  u1 = (3*q1(ibc,jbc,kbc,2) - q1(ibc,jbc+1,kbc,2)) / 2
                  v2 = (  q1(ibc,jbc,kbc,3) + q1(ibc,jbc+1,kbc,3)) / 2
                  v1 = (3*q1(ibc,jbc,kbc,3) - q1(ibc,jbc+1,kbc,3)) / 2
                  w2 = (  q1(ibc,jbc,kbc,4) + q1(ibc,jbc+1,kbc,4)) / 2
                  w1 = (3*q1(ibc,jbc,kbc,4) - q1(ibc,jbc+1,kbc,4)) / 2
                else if (i2.eq.iend(ic2)) then
                  u2 = (3*q1(ibc,jbc,kbc,2) - q1(ibc,jbc-1,kbc,2)) / 2
                  u1 = (  q1(ibc,jbc,kbc,2) + q1(ibc,jbc-1,kbc,2)) / 2
                  v2 = (3*q1(ibc,jbc,kbc,3) - q1(ibc,jbc-1,kbc,3)) / 2
                  v1 = (  q1(ibc,jbc,kbc,3) + q1(ibc,jbc-1,kbc,3)) / 2
                  w2 = (3*q1(ibc,jbc,kbc,4) - q1(ibc,jbc-1,kbc,4)) / 2
                  w1 = (  q1(ibc,jbc,kbc,4) + q1(ibc,jbc-1,kbc,4)) / 2
                else
                  u2 = (q1(ibc,jbc  ,kbc,2) + q1(ibc,jbc+1,kbc,2)) / 2
                  u1 = (q1(ibc,jbc-1,kbc,2) + q1(ibc,jbc  ,kbc,2)) / 2
                  v2 = (q1(ibc,jbc  ,kbc,3) + q1(ibc,jbc+1,kbc,3)) / 2
                  v1 = (q1(ibc,jbc-1,kbc,3) + q1(ibc,jbc  ,kbc,3)) / 2
                  w2 = (q1(ibc,jbc  ,kbc,4) + q1(ibc,jbc+1,kbc,4)) / 2
                  w1 = (q1(ibc,jbc-1,kbc,4) + q1(ibc,jbc  ,kbc,4)) / 2
                end if
                dudeta = u2 - u1
                dvdeta = v2 - v1
                dwdeta = w2 - w1

                if (i3.eq.istrt(ic3)) then
                  u2 = (  q1(ibc,jbc,kbc,2) + q1(ibc,jbc,kbc+1,2)) / 2
                  u1 = (3*q1(ibc,jbc,kbc,2) - q1(ibc,jbc,kbc+1,2)) / 2
                  v2 = (  q1(ibc,jbc,kbc,3) + q1(ibc,jbc,kbc+1,3)) / 2
                  v1 = (3*q1(ibc,jbc,kbc,3) - q1(ibc,jbc,kbc+1,3)) / 2
                  w2 = (  q1(ibc,jbc,kbc,4) + q1(ibc,jbc,kbc+1,4)) / 2
                  w1 = (3*q1(ibc,jbc,kbc,4) - q1(ibc,jbc,kbc+1,4)) / 2
                else if (i3.eq.iend(ic3)) then
                  u2 = (3*q1(ibc,jbc,kbc,2) - q1(ibc,jbc,kbc-1,2)) / 2
                  u1 = (  q1(ibc,jbc,kbc,2) + q1(ibc,jbc,kbc-1,2)) / 2
                  v2 = (3*q1(ibc,jbc,kbc,3) - q1(ibc,jbc,kbc-1,3)) / 2
                  v1 = (  q1(ibc,jbc,kbc,3) + q1(ibc,jbc,kbc-1,3)) / 2
                  w2 = (3*q1(ibc,jbc,kbc,4) - q1(ibc,jbc,kbc-1,4)) / 2
                  w1 = (  q1(ibc,jbc,kbc,4) + q1(ibc,jbc,kbc-1,4)) / 2
                else
                  u2 = (q1(ibc,jbc,kbc  ,2) + q1(ibc,jbc,kbc+1,2)) / 2
                  u1 = (q1(ibc,jbc,kbc-1,2) + q1(ibc,jbc,kbc  ,2)) / 2
                  v2 = (q1(ibc,jbc,kbc  ,3) + q1(ibc,jbc,kbc+1,3)) / 2
                  v1 = (q1(ibc,jbc,kbc-1,3) + q1(ibc,jbc,kbc  ,3)) / 2
                  w2 = (q1(ibc,jbc,kbc  ,4) + q1(ibc,jbc,kbc+1,4)) / 2
                  w1 = (q1(ibc,jbc,kbc-1,4) + q1(ibc,jbc,kbc  ,4)) / 2
                end if
                dudzeta = u2 - u1
                dvdzeta = v2 - v1
                dwdzeta = w2 - w1

                xix   = xideri  (1,im,jm,km)
                xiy   = xideri  (2,im,jm,km)
                xiz   = xideri  (3,im,jm,km)
                etax  = etaderi (1,im,jm,km)
                etay  = etaderi (2,im,jm,km)
                etaz  = etaderi (3,im,jm,km)
                zetax = zetaderi(1,im,jm,km)
                zetay = zetaderi(2,im,jm,km)
                zetaz = zetaderi(3,im,jm,km)

              else if (idir.eq.2) then

                if (i3.eq.istrt(ic3)) then
                  u2 = (  q1(ibc+1,jbc,kbc,2) + q1(ibc  ,jbc,kbc,2))/2
                  u1 = (3*q1(ibc  ,jbc,kbc,2) - q1(ibc+1,jbc,kbc,2))/2
                  v2 = (  q1(ibc+1,jbc,kbc,3) + q1(ibc  ,jbc,kbc,3))/2
                  v1 = (3*q1(ibc  ,jbc,kbc,3) - q1(ibc+1,jbc,kbc,3))/2
                  w2 = (  q1(ibc+1,jbc,kbc,4) + q1(ibc  ,jbc,kbc,4))/2
                  w1 = (3*q1(ibc  ,jbc,kbc,4) - q1(ibc+1,jbc,kbc,4))/2
                else if (i3.eq.iend(ic3)) then
                  u2 = (3*q1(ibc,jbc,kbc,2) - q1(ibc-1,jbc,kbc,2)) / 2
                  u1 = (  q1(ibc,jbc,kbc,2) + q1(ibc-1,jbc,kbc,2)) / 2
                  v2 = (3*q1(ibc,jbc,kbc,3) - q1(ibc-1,jbc,kbc,3)) / 2
                  v1 = (  q1(ibc,jbc,kbc,3) + q1(ibc-1,jbc,kbc,3)) / 2
                  w2 = (3*q1(ibc,jbc,kbc,4) - q1(ibc-1,jbc,kbc,4)) / 2
                  w1 = (  q1(ibc,jbc,kbc,4) + q1(ibc-1,jbc,kbc,4)) / 2
                else
                  u2 = (q1(ibc+1,jbc,kbc,2) + q1(ibc  ,jbc,kbc,2)) / 2
                  u1 = (q1(ibc  ,jbc,kbc,2) + q1(ibc-1,jbc,kbc,2)) / 2
                  v2 = (q1(ibc+1,jbc,kbc,3) + q1(ibc  ,jbc,kbc,3)) / 2
                  v1 = (q1(ibc  ,jbc,kbc,3) + q1(ibc-1,jbc,kbc,3)) / 2
                  w2 = (q1(ibc+1,jbc,kbc,4) + q1(ibc  ,jbc,kbc,4)) / 2
                  w1 = (q1(ibc  ,jbc,kbc,4) + q1(ibc-1,jbc,kbc,4)) / 2
                end if
                dudxi = u2 - u1
                dvdxi = v2 - v1
                dwdxi = w2 - w1

                dudeta = (-1)**(1+1/inrout) * &
                         (-8*q1(ibc,jbc,kbc,2) + &
                           9*q1(ibc,jbc+joff1,kbc,2) - &
                             q1(ibc,jbc+joff2,kbc,2) ) / 3
                dvdeta = (-1)**(1+1/inrout) * &
                         (-8*q1(ibc,jbc,kbc,3) + &
                           9*q1(ibc,jbc+joff1,kbc,3) - &
                             q1(ibc,jbc+joff2,kbc,3) ) / 3
                dwdeta = (-1)**(1+1/inrout) * &
                         (-8*q1(ibc,jbc,kbc,4) + &
                           9*q1(ibc,jbc+joff1,kbc,4) - &
                             q1(ibc,jbc+joff2,kbc,4) ) / 3

                if (i2.eq.istrt(ic2)) then
                  u2 = (  q1(ibc,jbc,kbc+1,2) + q1(ibc,jbc,kbc  ,2))/2
                  u1 = (3*q1(ibc,jbc,kbc  ,2) - q1(ibc,jbc,kbc+1,2))/2
                  v2 = (  q1(ibc,jbc,kbc+1,3) + q1(ibc,jbc,kbc  ,3))/2
                  v1 = (3*q1(ibc,jbc,kbc  ,3) - q1(ibc,jbc,kbc+1,3))/2
                  w2 = (  q1(ibc,jbc,kbc+1,4) + q1(ibc,jbc,kbc  ,4))/2
                  w1 = (3*q1(ibc,jbc,kbc  ,4) - q1(ibc,jbc,kbc+1,4))/2
                else if (i2.eq.iend(ic2)) then
                  u2 = (3*q1(ibc,jbc,kbc,2) - q1(ibc,jbc,kbc-1,2)) / 2
                  u1 = (  q1(ibc,jbc,kbc,2) + q1(ibc,jbc,kbc-1,2)) / 2
                  v2 = (3*q1(ibc,jbc,kbc,3) - q1(ibc,jbc,kbc-1,3)) / 2
                  v1 = (  q1(ibc,jbc,kbc,3) + q1(ibc,jbc,kbc-1,3)) / 2
                  w2 = (3*q1(ibc,jbc,kbc,4) - q1(ibc,jbc,kbc-1,4)) / 2
                  w1 = (  q1(ibc,jbc,kbc,4) + q1(ibc,jbc,kbc-1,4)) / 2
                else
                  u2 = (q1(ibc,jbc,kbc+1,2) + q1(ibc,jbc,kbc  ,2)) / 2
                  u1 = (q1(ibc,jbc,kbc  ,2) + q1(ibc,jbc,kbc-1,2)) / 2
                  v2 = (q1(ibc,jbc,kbc+1,3) + q1(ibc,jbc,kbc  ,3)) / 2
                  v1 = (q1(ibc,jbc,kbc  ,3) + q1(ibc,jbc,kbc-1,3)) / 2
                  w2 = (q1(ibc,jbc,kbc+1,4) + q1(ibc,jbc,kbc  ,4)) / 2
                  w1 = (q1(ibc,jbc,kbc  ,4) + q1(ibc,jbc,kbc-1,4)) / 2
                end if
                dudzeta = u2 - u1
                dvdzeta = v2 - v1
                dwdzeta = w2 - w1
   
                xix   = xiderj  (1,im,jm,km)
                xiy   = xiderj  (2,im,jm,km)
                xiz   = xiderj  (3,im,jm,km)
                etax  = etaderj (1,im,jm,km)
                etay  = etaderj (2,im,jm,km)
                etaz  = etaderj (3,im,jm,km)
                zetax = zetaderj(1,im,jm,km)
                zetay = zetaderj(2,im,jm,km)
                zetaz = zetaderj(3,im,jm,km)

              else if (idir.eq.3) then
                if (i2.eq.istrt(ic2)) then
                  u2 = (  q1(ibc,jbc,kbc,2) + q1(ibc+1,jbc,kbc,2)) / 2
                  u1 = (3*q1(ibc,jbc,kbc,2) - q1(ibc+1,jbc,kbc,2)) / 2
                  v2 = (  q1(ibc,jbc,kbc,3) + q1(ibc+1,jbc,kbc,3)) / 2
                  v1 = (3*q1(ibc,jbc,kbc,3) - q1(ibc+1,jbc,kbc,3)) / 2
                  w2 = (  q1(ibc,jbc,kbc,4) + q1(ibc+1,jbc,kbc,4)) / 2
                  w1 = (3*q1(ibc,jbc,kbc,4) - q1(ibc+1,jbc,kbc,4)) / 2
                else if (i2.eq.iend(ic2)) then
                  u2 = (3*q1(ibc,jbc,kbc,2) - q1(ibc-1,jbc,kbc,2)) / 2
                  u1 = (  q1(ibc,jbc,kbc,2) + q1(ibc-1,jbc,kbc,2)) / 2
                  v2 = (3*q1(ibc,jbc,kbc,3) - q1(ibc-1,jbc,kbc,3)) / 2
                  v1 = (  q1(ibc,jbc,kbc,3) + q1(ibc-1,jbc,kbc,3)) / 2
                  w2 = (3*q1(ibc,jbc,kbc,4) - q1(ibc-1,jbc,kbc,4)) / 2
                  w1 = (  q1(ibc,jbc,kbc,4) + q1(ibc-1,jbc,kbc,4)) / 2
                else
                  u2 = (q1(ibc  ,jbc,kbc,2) + q1(ibc+1,jbc,kbc,2)) / 2
                  u1 = (q1(ibc-1,jbc,kbc,2) + q1(ibc  ,jbc,kbc,2)) / 2
                  v2 = (q1(ibc  ,jbc,kbc,3) + q1(ibc+1,jbc,kbc,3)) / 2
                  v1 = (q1(ibc-1,jbc,kbc,3) + q1(ibc  ,jbc,kbc,3)) / 2
                  w2 = (q1(ibc  ,jbc,kbc,4) + q1(ibc+1,jbc,kbc,4)) / 2
                  w1 = (q1(ibc-1,jbc,kbc,4) + q1(ibc  ,jbc,kbc,4)) / 2
                end if
                dudxi = u2 - u1
                dvdxi = v2 - v1
                dwdxi = w2 - w1
 
                if (i3.eq.istrt(ic3)) then
                  u2 = (  q1(ibc,jbc,kbc,2) + q1(ibc,jbc+1,kbc,2)) / 2
                  u1 = (3*q1(ibc,jbc,kbc,2) - q1(ibc,jbc+1,kbc,2)) / 2
                  v2 = (  q1(ibc,jbc,kbc,3) + q1(ibc,jbc+1,kbc,3)) / 2
                  v1 = (3*q1(ibc,jbc,kbc,3) - q1(ibc,jbc+1,kbc,3)) / 2
                  w2 = (  q1(ibc,jbc,kbc,4) + q1(ibc,jbc+1,kbc,4)) / 2
                  w1 = (3*q1(ibc,jbc,kbc,4) - q1(ibc,jbc+1,kbc,4)) / 2
                else if (i3.eq.iend(ic3)) then
                  u2 = (3*q1(ibc,jbc,kbc,2) - q1(ibc,jbc-1,kbc,2)) / 2
                  u1 = (  q1(ibc,jbc,kbc,2) + q1(ibc,jbc-1,kbc,2)) / 2
                  v2 = (3*q1(ibc,jbc,kbc,3) - q1(ibc,jbc-1,kbc,3)) / 2
                  v1 = (  q1(ibc,jbc,kbc,3) + q1(ibc,jbc-1,kbc,3)) / 2
                  w2 = (3*q1(ibc,jbc,kbc,4) - q1(ibc,jbc-1,kbc,4)) / 2
                  w1 = (  q1(ibc,jbc,kbc,4) + q1(ibc,jbc-1,kbc,4)) / 2
                else
                  u2 = (q1(ibc,jbc  ,kbc,2) + q1(ibc,jbc+1,kbc,2)) / 2
                  u1 = (q1(ibc,jbc-1,kbc,2) + q1(ibc,jbc  ,kbc,2)) / 2
                  v2 = (q1(ibc,jbc  ,kbc,3) + q1(ibc,jbc+1,kbc,3)) / 2
                  v1 = (q1(ibc,jbc-1,kbc,3) + q1(ibc,jbc  ,kbc,3)) / 2
                  w2 = (q1(ibc,jbc  ,kbc,4) + q1(ibc,jbc+1,kbc,4)) / 2
                  w1 = (q1(ibc,jbc-1,kbc,4) + q1(ibc,jbc  ,kbc,4)) / 2
                end if
                dudeta = u2 - u1
                dvdeta = v2 - v1
                dwdeta = w2 - w1
 
                dudzeta = (-1)**(1+1/inrout) * &
                          (-8*q1(ibc,jbc,kbc,2) + &
                            9*q1(ibc,jbc,kbc+koff1,2) - &
                              q1(ibc,jbc,kbc+koff2,2) ) / 3
                dvdzeta = (-1)**(1+1/inrout) * &
                          (-8*q1(ibc,jbc,kbc,3) + &
                            9*q1(ibc,jbc,kbc+koff1,3) - &
                              q1(ibc,jbc,kbc+koff2,3) ) / 3
                dwdzeta = (-1)**(1+1/inrout) * &
                          (-8*q1(ibc,jbc,kbc,4) + &
                            9*q1(ibc,jbc,kbc+koff1,4) - &
                              q1(ibc,jbc,kbc+koff2,4) ) / 3
 
                xix   = xiderk  (1,im,jm,km)
                xiy   = xiderk  (2,im,jm,km)
                xiz   = xiderk  (3,im,jm,km)
                etax  = etaderk (1,im,jm,km)
                etay  = etaderk (2,im,jm,km)
                etaz  = etaderk (3,im,jm,km)
                zetax = zetaderk(1,im,jm,km)
                zetay = zetaderk(2,im,jm,km)
                zetaz = zetaderk(3,im,jm,km)
   
              end if

              dudx = dudxi*xix + dudeta*etax + dudzeta*zetax
              dudy = dudxi*xiy + dudeta*etay + dudzeta*zetay
              dudz = dudxi*xiz + dudeta*etaz + dudzeta*zetaz
              dvdx = dvdxi*xix + dvdeta*etax + dvdzeta*zetax
              dvdy = dvdxi*xiy + dvdeta*etay + dvdzeta*zetay
              dvdz = dvdxi*xiz + dvdeta*etaz + dvdzeta*zetaz
              dwdx = dwdxi*xix + dwdeta*etax + dwdzeta*zetax
              dwdy = dwdxi*xiy + dwdeta*etay + dwdzeta*zetay
              dwdz = dwdxi*xiz + dwdeta*etaz + dwdzeta*zetaz

              divv = dudx + dvdy + dwdz

              txx = 2*muw * (dudx - divv/3)
              txy =   muw * (dudy + dvdx)
              txz =   muw * (dudz + dwdx)
              tyy = 2*muw * (dvdy - divv/3)
              tyz =   muw * (dvdz + dwdy)
              tzz = 2*muw * (dwdz - divv/3)


              tauwx = txx*kx + txy*ky + txz*kz
              tauwy = txy*kx + tyy*ky + tyz*kz
              tauwz = txz*kx + tyz*ky + tzz*kz

              tauwpx = tauwx - (tauwx*kx+tauwy*ky+tauwz*kz)*kx
              tauwpy = tauwy - (tauwx*kx+tauwy*ky+tauwz*kz)*ky
              tauwpz = tauwz - (tauwx*kx+tauwy*ky+tauwz*kz)*kz

              tauw   = dsqrt(tauwpx**2+tauwpy**2+tauwpz**2)
              utau   = dsqrt(tauw/rhow)

              cf   (sur_dim,nsurf1) = 2.d0 * tauw / ( reyno * machfs )
              yplus(sur_dim,nsurf1) = &
                dsqrt(reyno/machfs) * rhow * dn1 * utau / muw

            end do
          end do

          if (sur_dim.ne.surf_dim(nsurf1)) then
              write(*,*) 'Incongruency in sfht_b_h. Aborting!'
            stop
          end if

        end if

      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine blprofil(nl,itime1,xc,yc,xycdot,si,sj,q1,q2,dist,bctopo)
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,itime1,iblock,iblk_pos,imax,jmax,ijmaxg,ixyc,ixyc1, &
                idist,iimt,iq,iblk
      real (kind=cosa_real) xc(*),yc(*),xycdot(*),dist(*),si(*),sj(*),q1(*), &
           q2(*)
      integer(kind=cosa_int) bctopo(*)

      do iblock = 1,nblocks

        if (iblock.eq.profpar(1)) then

          imax     = i_imax(iblock,nl)
          jmax     = j_jmax(iblock,nl)
          ixyc     = 1 + off_p2 (iblock,nl)
          ixyc1    = 1 + off_p3 (iblock,nl) * npde * dim5
          idist    = 1 + off_p1 (iblock,nl)
          iimt     = 1 + off_p2 (iblock,nl) * lmet * dim5h
          iq       = 1 + off_p3 (iblock,nl) * npde
          iblk     = 1 + off_bct(iblock,nl)

                 
          call blprofil_b(imax,jmax,npde,lmet,nharms,itime1, &
                          nbcs(iblock),bctopo(iblk),xc(ixyc),yc(ixyc), &
                          xycdot(ixyc1),dist(idist),si(iimt),sj(iimt), &
                          q1(iq),q2(iq))

        end if

      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine blprofil_b(imax,jmax,npde,lmet,nharms,itime1,nbcs, &
                            bctopo,xc,yc,xycdot,dist,si,sj,q1,q2)
!-----------------------------------------------------------------------
!     Routine to output a boundary layer profile at user-given location.
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,npde,nharms,nbcs,lmet
      integer(kind=cosa_int) bctopo(7,nbcs)
      integer(kind=cosa_int) itime1,n,nh
      real (kind=cosa_real) &
           xc(0:imax+1,0:jmax+1), &
           yc(0:imax+1,0:jmax+1), &
           xycdot (-1:imax+1,-1:jmax+1,npde,0:2*nharms       ), &
           si     (lmet,0:imax+1 ,0:jmax+1 ,0:2*nharms*hbmove), &
           sj     (lmet,0:imax+1 ,0:jmax+1 ,0:2*nharms*hbmove), &
           q1(-1:imax+1,-1:jmax+1,npde), &
           q2(-1:imax+1,-1:jmax+1,npde), &
           dist(0:imax,0:jmax)

      n  = itime1 - 1
      nh = n*hbmove
      call blprofil_b_h(imax,jmax,npde,lmet,xc,yc,xycdot(-1,-1,1,n), &
                        dist,si(1,0,0,n),sj(1,0,0,n),q1,q2,bctopo,nbcs)

      return
      end

!-----------------------------------------------------------------------
      subroutine blprofil_b_h(imax,jmax,npde,lmet,xc,yc,xycdot,dist,si, &
                              sj,q1,q2,bctopo,nbcs)
!-----------------------------------------------------------------------
!     Routine to output the boundary layer profile at an i=const. line 
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax, jmax, npde, nbcs, lmet
      integer(kind=cosa_int) iblock, i, j, jw, jpstart, jpend, jpstep, jb, &
                iprof, off1, off2, offbot, offtop, ibc
      real (kind=cosa_real) &
           xc(0:imax+1,0:jmax+1), &
           yc(0:imax+1,0:jmax+1), &
           xycdot (-1:imax+1,-1:jmax+1,npde), &
           si(lmet,0:imax+1,0:jmax+1), sj(lmet,0:imax+1,0:jmax+1), &
           q1(-1:imax+1,-1:jmax+1,npde),q2(-1:imax+1,-1:jmax+1,npde), &
           dist(0:imax,0:jmax)
      real (kind=cosa_real) mup, rhop, pp, tp, up, vp, uc, vc, uorth, upar, &
           udotn1, udotn2, upar1, upar2, dn1, dn2, dudn, dtdn, cf, &
           rhow, pw, tw, muw, uw, vw, uworth, uwpar, tauw, utau, yplus, &
           uplus, dupardn, fact
      integer(kind=cosa_int) bctopo(7,nbcs)
      character*72 name1,name2

      fact = 1.d-12

      iblock = profpar(1)
      if (.not.unsteady) then

        if (iblock.le.9) then
          write(filename,'(''bl_profil_b0'',i1,''.dat'')') iblock
        else
          write(filename,'(''bl_profil_b'',i2,''.dat'')') iblock
        end if

      else

        if (iblock.le.9) then
          write(name1,'(''b0'',i1,''_nt_'')') iblock
        else
          write(name1,'(''b'',i2,''_nt_'')') iblock
        end if

        if (itime.le.9) then
          write(name2,'(''000'',i1)') itime
        elseif (itime.le.99) then
          write(name2,'(''00'',i2)') itime
        elseif (itime.le.999) then
          write(name2,'(''0'',i3)') itime
        else
          write(name2,'(i4)') itime
        end if
        filename = 'bl_profil_'//name1
        filename = trim(filename)//name2
        filename = trim(filename)//'.dat'
      end if

      open(1,file=filename,status='replace')

      if (kom.or.kom_bsl.or.kom_sst) then
        write (1,1006) 
      else
        write (1,1005) 
      end if

      iprof   = profpar(2)
      jpstart = profpar(3)
      jpend   = profpar(4)
      jpstep  = sign(1,jpend-jpstart)

!     check if wall is on j=1 or on j=jmax boundary
      if (jpstep.eq.1) then
        jb = 1
      else if (jpstep.eq.-1) then
        jb = 2
      end if

      if ((jpstart.lt.1).or.(jpstart.gt.jmax-1)) then
        write(*,*) 'jpstart outside admissable bounds. Aborting!'
        stop
      else if ((jpend.lt.1).or.(jpend.gt.jmax-1)) then
        write(*,*) 'jpend outside admissable bounds. Aborting!'
        stop
      end if

      off2  = (jb-1)*(jmax-2)
      off1  = (jb-1)*(jmax-1)

! MSC 22 February 2018. Loop below works only if j dir is 
!                       normal wall. Otherwise generalisation of 
!                       implementation is needed.
      do ibc=1,nbcs

        if ((bctopo(2,ibc).eq.2).and.(bctopo(1,ibc)/100.eq.wallmask)) then

          if (((bctopo(3,ibc).eq.1).and.(jb.eq.1)).or. &
              ((bctopo(3,ibc).eq.jmax).and.(jb.eq.2))) then
            continue
          else 
            write(*,*) 'Inconsistency in blprof. Aborting!'
            stop
          end if

          do i=bctopo(4,ibc),bctopo(5,ibc)

!---------- check if required profile is actually on a wall
            if (i.eq.iprof) then
!------------ check if required profile actually starts on wall
              if (jpstart.eq.1+off2) then
                goto 10
              else
                write(*,*) 'Requested profile does not start on wall. Ab &
     &orting!'
                stop
              end if
            end if
          end do
          write(*,*) 'Requested profile does not lie on a wall. Aborting &
     &!'
          stop
        end if
      end do
      write(*,*) 'No wall boundary under this profile has been found. Ab &
     &orting!'
      stop

 10   offbot = (-1)**jb
      offtop = (-1)**(jb-1)

        do j = jpstart+offbot,jpend+offtop,jpstep

          rhop = q1(iprof,j,1)
          up   = q1(iprof,j,2)
          vp   = q1(iprof,j,3)
          pp   = q1(iprof,j,4)
          tp   = q2(iprof,j,1)
          mup  = (stemp+1.d0)/(tp+stemp) * (tp**1.5d0)

!-------- turbulent case: determine profile of y+ vs. u+
          if (kom.or.kom_bsl.or.kom_sst) then
!           determine tauw
            if (j.eq.jpstart+offbot) then

              jw   = 1 + off1
              rhow = rhop
              pw   = pp
              tw   = tp
              muw  = mup

!             calculate the velocity parallel to the wall at the cell
!             centre off the wall.  
!
!                     U(parallel) = U - U(orthogonal) 
!
!             where U is the total velocity vector
!             and U(orthogonal) is the velocity normal to the wall.  
!             For fixed-nody problems, velocity on the wall is zero.
!             For moving-nody problems, velocity on the wall is nonzero.
!             uorth is the velocity component (with sign) normal to the
!             wall at the first cell centre off the wall.
!             upar is the modulus of the velocity component parallel to 
!             the wall at the first cell centre off the wall.

              uc = q1(iprof,j+jpstep,2)
              vc = q1(iprof,j+jpstep,3)
              if (moving) then
                uc = uc - xycdot(iprof,j+jpstep,1)
                vc = vc - xycdot(iprof,j+jpstep,2)
              end if
              uorth = sj(1,iprof,jw) * uc + &
                      sj(2,iprof,jw) * vc
              upar  = sqrt( (uc-uorth*sj(1,iprof,jw))**2 + &
                            (vc-uorth*sj(2,iprof,jw))**2 )


!cc           if (moving) then
!cc             uw =  sj(iprof,jw,6)
!cc             vw =  sj(iprof,jw,7)
!cc             uworth = sj(iprof,jw,1) * uw + &
!cc                      sj(iprof,jw,2) * vw
!cc             uwpar  = sqrt( (uw-uorth*sj(iprof,jw,1))**2 + &
!cc                            (vw-uorth*sj(iprof,jw,2))**2 )
!cc           else
!cc             uw     = 0.d0
!cc             vw     = 0.d0
!cc             uworth = 0.d0
!cc             uwpar  = 0.d0
!cc           endif 

!             Tau(wall) = mu * du/dn

!cc           dupardn = (upar-uwpar) / (dist(iprof,jpstart) + fact)
              dupardn = (upar      ) / (dist(iprof,jpstart) + fact)
              tauw    = muw * dupardn
              utau    = dsqrt(tauw/rhow)

            end if

            yplus = dsqrt(reyno/machfs) * rhow*dist(iprof,j)*utau/muw

            uc = q1(iprof,j,2)
            vc = q1(iprof,j,3)
            if (moving) then
              uc = uc - xycdot(iprof,j,1)
              vc = vc - xycdot(iprof,j,2)
            end if
            uorth = sj(1,iprof,jw) * uc + &
                    sj(2,iprof,jw) * vc
            upar  = sqrt( (uc-uorth*sj(1,iprof,jw))**2 + &
                          (vc-uorth*sj(2,iprof,jw))**2 )

            uplus = dsqrt(reyno/machfs) * upar / utau

            write(1,1011) iprof, j, xc(iprof,j), yc(iprof,j), rhop, up, &
                          vp, pp, tp, mup, yplus, uplus, q2(iprof,j,3), &
                          q2(iprof,j,4), q2(iprof,j,5)

          else

            write(1,1010) iprof, j, xc(iprof,j), yc(iprof,j), rhop, up, &
                          vp, pp, tp, mup

          end if

        end do

      close(unit=1)

 1005 format (/' ',3x,'i',4x,'j',2x,'x',12x,'y',12x,'rho',10x,'u',12x, &
         'v',12x,'p',12x,'t',12x,'mu')
 1006 format (/' ',3x,'i',4x,'j',2x,'x',12x,'y',12x,'rho',10x,'u',12x, &
         'v',12x,'p',12x,'t',12x,'mul',12x,'yplus',12x,'uplus',12x, &
         'tke',12x,'omega',12x,'mut')
 1010 format (' ',i4,1x,i4,8(e13.5))
 1011 format (' ',i4,1x,i4,13(e13.5))

      return
      end

!-----------------------------------------------------------------------
      subroutine rdtec(nl,xc,yc,zc,q1,q2)
!-----------------------------------------------------------------------
!     q1(:,:,1) : rho
!     q1(:,:,2) : u
!     q1(:,:,3) : v
!     q1(:,:,4) : w
!     q1(:,:,5) : p
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,iblock,imax,jmax,kmax,iq,ixyzc,fid,dummy
      real(kind=cosa_real) xc(*),yc(*),zc(*),q1(*),q2(*)
      integer, allocatable :: blockindexes(:,:)

      character*72 line

      fid = 2

      if (tecbin) then
!-------Open unformatted file
        open(unit=fid,file=flowtec,form='unformatted',access='stream', &
             status='old')
!-------reading logicals and no. of blocks
        read(fid) dummy
        read(fid) dummy
        read(fid) dummy
        read(fid) dummy
        read(fid) dummy
        read(fid) dummy
        write(*,*) 'No of blocks=', dummy
        if (unsteady) then
          read(fid) simtime
          write(*,*) 'simtime=', simtime
        end if
        read(fid) dummy
        write(*,*) 'No of harms.=', dummy

        allocate(blockindexes(3,nblocks))
        read(fid) blockindexes
        deallocate(blockindexes)
      else
!-------Open formatted file
        open(fid,file=flowtec,status='old')
        read(fid,*)
        read(fid,*)
      end if

      do iblock = 1,nblocks
        imax   = i_imax     (iblock,nl)
        jmax   = j_jmax     (iblock,nl)
        kmax   = k_kmax     (iblock,nl)
        iq     = 1 + off_p3 (iblock,nl) * npde
        ixyzc  = 1 + off_p2 (iblock,nl)
        call rdtec_b(imax,jmax,kmax,npde,xc(ixyzc),yc(ixyzc),zc(ixyzc), &
                     q1(iq),q2(iq),fid)
      end do

      if (unsteady.and.(.not.tecbin)) then
        read(fid,15) simtime
        write(*,*) 'simtime=', simtime
      else
        simtime = 0.d0
      end if

      write(line,'(''Read in file '',a40)') flowtec
      write(*,'(a)') line

      close(unit=fid)

 15   format(e21.14)

      return
      end

!-----------------------------------------------------------------------
      subroutine rdtec_b(imax,jmax,kmax,npde,xc,yc,zc,q1,q2,fid)
!-----------------------------------------------------------------------
!     q1(:,:,1) : rho
!     q1(:,:,2) : u
!     q1(:,:,3) : v
!     q1(:,:,4) : w
!     q1(:,:,5) : p
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,npde,fid
      integer(kind=cosa_int) i,j,k,ipde,idummy
      real(kind=cosa_real) xc(0:imax+1,0:jmax+1,0:kmax+1), &
                   yc(0:imax+1,0:jmax+1,0:kmax+1), &
                   zc(0:imax+1,0:jmax+1,0:kmax+1), &
                   q1(-1:imax+1,-1:jmax+1,-1:kmax+1,npde), &
                   q2(-1:imax+1,-1:jmax+1,-1:kmax+1,npde)
      real(kind=cosa_real) dummy

      if (tecbin) then
        if (kom.or.kom_bsl.or.kom_sst) then
          read (fid) (((xc(i,j,k),i=0,imax),j=0,jmax),k=0,kmax), &
                     (((yc(i,j,k),i=0,imax),j=0,jmax),k=0,kmax), &
                     (((zc(i,j,k),i=0,imax),j=0,jmax),k=0,kmax), &
                     ((((q1(i,j,k,ipde),i=0,imax),j=0,jmax), &
                       k=0,kmax),ipde=1,5), &
                     ((((q2(i,j,k,ipde),i=0,imax),j=0,jmax), &
                       k=0,kmax),ipde=1,5), &
                     (((dummy,i=0,imax),j=0,jmax),k=0,kmax) 
        else
          read (fid) (((xc(i,j,k),i=0,imax),j=0,jmax),k=0,kmax), &
                     (((yc(i,j,k),i=0,imax),j=0,jmax),k=0,kmax), &
                     (((zc(i,j,k),i=0,imax),j=0,jmax),k=0,kmax), &
                     ((((q1(i,j,k,ipde),i=0,imax),j=0,jmax), &
                       k=0,kmax),ipde=1,5), &
                     ((((q2(i,j,k,ipde),i=0,imax),j=0,jmax), &
                       k=0,kmax),ipde=1,2)
        end if
      else
        read(fid,*)
        if (kom.or.kom_bsl.or.kom_sst) then
          read (fid,14) (((xc(i,j,k),i=0,imax),j=0,jmax),k=0,kmax), &
                        (((yc(i,j,k),i=0,imax),j=0,jmax),k=0,kmax), &
                        (((zc(i,j,k),i=0,imax),j=0,jmax),k=0,kmax), &
                        ((((q1(i,j,k,ipde),i=0,imax),j=0,jmax), &
                          k=0,kmax),ipde=1,5), &
                        ((((q2(i,j,k,ipde),i=0,imax),j=0,jmax), &
                          k=0,kmax),ipde=1,5), &
                        (((dummy,i=0,imax),j=0,jmax),k=0,kmax) 
        else
          read (fid,10) (((xc(i,j,k),i=0,imax),j=0,jmax),k=0,kmax), &
                        (((yc(i,j,k),i=0,imax),j=0,jmax),k=0,kmax), &
                        (((zc(i,j,k),i=0,imax),j=0,jmax),k=0,kmax), &
                        ((((q1(i,j,k,ipde),i=0,imax),j=0,jmax), &
                          k=0,kmax),ipde=1,5), &
                        ((((q2(i,j,k,ipde),i=0,imax),j=0,jmax), &
                          k=0,kmax),ipde=1,2)
        end if
      end if

 14   format(3e22.14)
 10   format(3e22.14)

      return
      end

!-----------------------------------------------------------------------
      subroutine cp_wall(nl,ijkmaxg,nsurf,msurf,p,q1,bctopo,surf_dim)
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,iblock,iblk_pos,imax,jmax,kmax,ijkmaxg,iq,iblk,nsurf, &
                msurf,surf_dim(msurf),nsurf1
      real(kind=cosa_real) q1(*),p(ijkmaxg-1,msurf)
      integer(kind=cosa_int) bctopo(*)

      nsurf1 = 0
      do iblk_pos = 1,nblock_pos

        iblock   = iblock_pos(iblk_pos)
        imax     = i_imax(iblock,nl)
        jmax     = j_jmax(iblock,nl)
        kmax     = k_kmax(iblock,nl)
        iq       = 1 + off_p3 (iblock,nl) * npde
        iblk     = 1 + off_bct(iblock,nl)

        call cp_wall_b(imax,jmax,kmax,ijkmaxg,npde,nsurf1,msurf, &
                       nbcs(iblock),p,surf_dim,q1(iq),bctopo(iblk))

      end do

      if (nsurf1.ne.nsurf) then
        write(*,*) 'Incongruency in cp_wall. Aborting!'
        stop
      end if

      return
      end

!-----------------------------------------------------------------------
      subroutine cp_wall_b(imax,jmax,kmax,ijkmaxg,npde,nsurf1,msurf, &
                           nbcs,p,surf_dim,q1,bctopo)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,ijkmaxg,nbcs,npde,nsurf1,msurf, &
                surf_dim(msurf)
      integer(kind=cosa_int) i2,i3,ibbc,sur_dim,ijkmax(3),bctyp,idir,inrout,istrt(3), &
                iend(3),ibcpt,ic1,ic2,ic3,ibc,jbc,kbc,bctopo1
      real(kind=cosa_real) q1(-1:imax+1,-1:jmax+1,-1:kmax+1,npde), &
                   p(ijkmaxg-1,msurf)
      integer(kind=cosa_int) bctopo(10,nbcs)

      do ibbc=1,nbcs

        bctopo1 = bctopo(1,ibbc)
        if (bctopo1.eq.14) bctopo1=bctopo1*100
        if (bctopo1/100.eq.wallmask) then
!       if (bctopo(1,ibbc)/100.eq.wallmask) then

          sur_dim = 0
          nsurf1 = nsurf1 + 1

!         store imax, jmax in ijmax to enable boundary data location

          ijkmax(1) = imax
          ijkmax(2) = jmax
          ijkmax(3) = kmax

!         store boundary topology in mnemonic names

          bctyp    = bctopo(1,ibbc)
          idir     = bctopo(2,ibbc)
          inrout   = bctopo(3,ibbc)
          istrt(1) = bctopo(4,ibbc)
          iend(1)  = bctopo(5,ibbc)
          istrt(2) = bctopo(6,ibbc)
          iend(2)  = bctopo(7,ibbc)
          istrt(3) = bctopo(8,ibbc)
          iend(3)  = bctopo(9,ibbc)

!         modify beginning, ending indices to extend boundary condition
!         to edge/corner

!tmp      do i1 = 1,2
!tmp        if (i1.ne.idir) then
!tmp          if (istrt(i1) .eq. 1          ) istrt(i1) = 0
!tmp          if (iend (i1) .eq. ijmax(i1)-1) iend (i1) = ijmax(i1)
!tmp        end if
!tmp      end do

!         set needed variables depending on whether the boundary is the
!         the inner boundary (inrout = 1) or 
!         the outer boundary (inrout > 1)
!             ibcpt : boundary condition location (first aux. cell)
!             ibcpt2: boundary condition location outside the block from
!                     ibcpt (second aux. cell)
!             ibcn  : point to the inside of the block from ibcpt
!             ibcm  : location of the metrics

          if (inrout.eq.1) then
            ibcpt  =  0
          else
            ibcpt  = ijkmax(idir)
          end if

          ic1 = cyc (idir, 1)
          ic2 = cyc (idir, 2)
          ic3 = cyc (idir, 3)

          do i3 = istrt(ic3),iend(ic3)
            do i2 = istrt(ic2),iend(ic2)

              ibc  = ibcpt *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
              jbc  = ibcpt *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
              kbc  = ibcpt *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

              sur_dim = sur_dim + 1
              if (machfs.gt.1d-14) then
                p(sur_dim,nsurf1) = 2 * (q1(ibc,jbc,kbc,5) - 1/gamma) / &
                                    machfs**2
              else
                p(sur_dim,nsurf1) = (q1(ibc,jbc,kbc,5) - 1/gamma)
              end if

            end do
          end do

          if (sur_dim.ne.surf_dim(nsurf1)) then
            write(*,*) 'Incongruency in cp_wall_b. Aborting!'
            stop
          end if

        end if

      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine wr_cp(ijkmaxg,p,xsurf,ysurf,zsurf,nsurf,msurf,surf_dim, &
                       ijkmax_surf)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) ijkmaxg,nsurf,msurf,surf_dim(msurf)
      integer(kind=cosa_int) i,itime1,isurf,ijkmax_surf(msurf,3)
      real(kind=cosa_real) p(ijkmaxg-1,msurf,ntime),xsurf(ijkmaxg-1,msurf), &
                   ysurf(ijkmaxg-1,msurf),zsurf(ijkmaxg-1,msurf)

      if (dualt) then
        if (itime.le.9) then
          write(flowtec,'(''000'',i1)') itime
        elseif (itime.le.99) then
          write(flowtec,'(''00'',i2)') itime
        elseif (itime.le.999) then
          write(flowtec,'(''0'',i3)') itime
        else
          write(flowtec,'(i4)') itime
        end if
      else if (harbal) then
        if (itime.le.9) then
          write(flowtec,'(''00'',i1)') itime
        elseif (itime.le.99) then
          write(flowtec,'(''0'',i2)') itime
        else
          write(flowtec,'(i3)') itime
        end if
      else
        flowtec = 'steady'
      end if

      flowtec = 'flow_tec_cp_'//flowtec
      flowtec = trim(flowtec)//'.dat'

!---- write header of cp tecplot file
      open(200,file=flowtec,status='replace')
      write(200,*) 'TITLE="Static pressure coefficient"'
      write(200,*) 'VARIABLES=x,y,z,cp'

      do isurf = 1,nsurf
        write(200,11) ijkmax_surf(isurf,1), ijkmax_surf(isurf,2), &
                      ijkmax_surf(isurf,3)
        do i=1,surf_dim(isurf)
          write(200,10) xsurf(i,isurf),ysurf(i,isurf),zsurf(i,isurf), &
                        p(i,isurf,1)
        end do
      end do
      close(unit=200)    

!----- The lines below are for writing out the ascii file for each patch 
!      separately, as it is done in 2D code.

!     do isurf = 1,nsurf
!       if (unsteady) then
!
!         if (dualt) then
!           if (itime.le.9) then
!             write(flowtec,'(''000'',i1)') itime
!           elseif (itime.le.99) then
!             write(flowtec,'(''00'',i2)') itime
!           elseif (itime.le.999) then
!             write(flowtec,'(''0'',i3)') itime
!           else
!             write(flowtec,'(i4)') itime
!           end if
!         else if (harbal) then
!           if (itime.le.9) then
!             write(flowtec,'(''00'',i1)') itime
!           elseif (itime.le.99) then
!             write(flowtec,'(''0'',i2)') itime
!           else
!             write(flowtec,'(i3)') itime
!           end if
!         end if
!
!         if (isurf.le.9) then
!           write(filename,'(''000'',i1)') isurf
!         else if (isurf.le.99) then
!           write(filename,'(''00'',i2)') isurf
!         else if (isurf.le.999) then
!           write(filename,'(''0'',i3)') isurf
!         end if
!         flowtec = 'unst_cp_TS_'//flowtec
!         filename = trim(flowtec)//'_s_'//filename
!         filename = trim(filename)//'.dat'
!         open(7,file=filename,status='replace')
!         write(7,10) simtime
!         do i=1,surf_dim(isurf)
!           write(7,10) xsurf(i,isurf),ysurf(i,isurf),zsurf(i,isurf), &
!                       p(i,isurf,1)
!         end do
!         close(unit=7)    
!
!       else
!
!         if (isurf.le.9) then
!           write(filename,'(''000'',i1)') isurf
!         else if (isurf.le.99) then
!           write(filename,'(''00'',i2)') isurf
!         else if (isurf.le.999) then
!           write(filename,'(''0'',i3)') isurf
!         end if
!         filename = 'stea_cp_s'//filename
!         filename = trim(filename)//'.dat'
!         open(7,file=filename,status='replace')
!         do i=1,surf_dim(isurf)
!           write(7,10) xsurf(i,isurf),ysurf(i,isurf),zsurf(i,isurf), &
!                       p(i,isurf,1)
!         end do
!         close(unit=7)    
!
!       end if
!
!     end do


 10   format(300e14.6)
 11   format('ZONE T="arturo",I=',i4,',J=',i4,',K=',i4,', F=POINT, DT=(S &
     &INGLE SINGLE SINGLE SINGLE)')

      return
      end

!-----------------------------------------------------------------------
      subroutine surf_coord(imax,jmax,kmax,npde,ijkmaxg,nsurf,msurf, &
                            nbcs,xyz0,bctopo,xsurf,ysurf,zsurf,surf_dim, &
                            ijkmax_surf)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,npde,ijkmaxg,nsurf,msurf,nbcs, &
                surf_dim(msurf),ijkmax_surf(msurf,3)
      integer(kind=cosa_int) i2,i3,ibbc,sur_dim,ijkmax(3),bctyp,idir,inrout,istrt(3), &
                iend(3),ibcpt,ic1,ic2,ic3,ibc,jbc,kbc,bctopo1
      real(kind=cosa_real) xsurf(ijkmaxg-1,msurf),ysurf(ijkmaxg-1,msurf), &
         zsurf(ijkmaxg-1,msurf),xyz0(-1:imax+1,-1:jmax+1,-1:kmax+1,npde)
      integer(kind=cosa_int) bctopo(10,nbcs)

      do ibbc=1,nbcs

        bctopo1 = bctopo(1,ibbc)
        if (bctopo1.eq.14) bctopo1=bctopo1*100
        if (bctopo1/100.eq.wallmask) then
!       if (bctopo(1,ibbc)/100.eq.wallmask) then

          sur_dim = 0
          nsurf = nsurf + 1
          if (nsurf.eq.msurf) then
            write(*,*) 'Exceeded maximum number of surfaces. Aborting!'
            stop
          end if

!         store imax, jmax in ijmax to enable boundary data location

          ijkmax(1) = imax
          ijkmax(2) = jmax
          ijkmax(3) = kmax

!         store boundary topology in mnemonic names

          bctyp    = bctopo(1,ibbc)
          idir     = bctopo(2,ibbc)
          inrout   = bctopo(3,ibbc)
          istrt(1) = bctopo(4,ibbc)
          iend(1)  = bctopo(5,ibbc)
          istrt(2) = bctopo(6,ibbc)
          iend(2)  = bctopo(7,ibbc)
          istrt(3) = bctopo(8,ibbc)
          iend(3)  = bctopo(9,ibbc)

!         modify beginning, ending indices to extend boundary condition
!         to edge/corner

!tmp      do i1 = 1,2
!tmp        if (i1.ne.idir) then
!tmp          if (istrt(i1) .eq. 1          ) istrt(i1) = 0
!tmp          if (iend (i1) .eq. ijmax(i1)-1) iend (i1) = ijmax(i1)
!tmp        end if
!tmp      end do

!         set needed variables depending on whether the boundary is the
!         the inner boundary (inrout = 1) or
!         the outer boundary (inrout > 1)
!             ibcpt : boundary condition location (first aux. cell)
!             ibcpt2: boundary condition location outside the block from
!                     ibcpt (second aux. cell)
!             ibcn  : point to the inside of the block from ibcpt
!             ibcm  : location of the metrics

          if (inrout.eq.1) then
            ibcpt  =  0
          else
            ibcpt  = ijkmax(idir)
          end if

          ic1 = cyc (idir, 1)
          ic2 = cyc (idir, 2)
          ic3 = cyc (idir, 3)

          do i3 = istrt(ic3),iend(ic3)
            do i2 = istrt(ic2),iend(ic2)

              ibc  = ibcpt *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
              jbc  = ibcpt *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
              kbc  = ibcpt *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

              sur_dim = sur_dim + 1
              xsurf(sur_dim,nsurf) = xyz0(ibc,jbc,kbc,1)
              ysurf(sur_dim,nsurf) = xyz0(ibc,jbc,kbc,2)
              zsurf(sur_dim,nsurf) = xyz0(ibc,jbc,kbc,3)

            end do
          end do

          ijkmax_surf(nsurf,1) = iend(ic2) - istrt(ic2) + 1
          ijkmax_surf(nsurf,2) = 1
          ijkmax_surf(nsurf,3) = iend(ic3) - istrt(ic3) + 1
          surf_dim(nsurf) = sur_dim

        end if

      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine vorticity(q1,q2,imax,jmax,npde,xc,yc)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,npde
      integer(kind=cosa_int) i,j,ipde

      real(kind=cosa_real) xc(0:imax,0:jmax),yc(0:imax,0:jmax), &
                   q1(-1:imax+1,-1:jmax+1,npde), &
                   q2(-1:imax+1,-1:jmax+1,npde)
      real(kind=cosa_real) dxi, dxj, dyi, dyj, dui, duj, dvi, dvj 

      do i=0,imax
        do j=0,jmax
          if (i.eq.0) then
            dxi = -xc(i+2,j)   + 4*xc(i+1,j)   - 3*xc(i,j) 
            dyi = -yc(i+2,j)   + 4*yc(i+1,j)   - 3*yc(i,j)
            dui = -q1(i+2,j,2) + 4*q1(i+1,j,2) - 3*q1(i,j,2)
            dvi = -q1(i+2,j,3) + 4*q1(i+1,j,3) - 3*q1(i,j,3)
          end if
          if (i.eq.imax) then
            dxi =  xc(i-2,j)   - 4*xc(i-1,j)   + 3*xc(i,j)
            dyi =  yc(i-2,j)   - 4*yc(i-1,j)   + 3*yc(i,j)
            dui =  q1(i-2,j,2) - 4*q1(i-1,j,2) + 3*q1(i,j,2)
            dvi =  q1(i-2,j,3) - 4*q1(i-1,j,3) + 3*q1(i,j,3)
          end if
          if (j.eq.0) then
            dxj = -xc(i,j+2)   + 4*xc(i,j+1)   - 3*xc(i,j)
            dyj = -yc(i,j+2)   + 4*yc(i,j+1)   - 3*yc(i,j)
            duj = -q1(i,j+2,2) + 4*q1(i,j+1,2) - 3*q1(i,j,2)
            dvj = -q1(i,j+2,3) + 4*q1(i,j+1,3) - 3*q1(i,j,3)
          end if
          if (j.eq.jmax) then
            dxj =  xc(i,j-2)   - 4*xc(i,j-1)   + 3*xc(i,j)
            dyj =  yc(i,j-2)   - 4*yc(i,j-1)   + 3*yc(i,j)
            duj =  q1(i,j-2,2) - 4*q1(i,j-1,2) + 3*q1(i,j,2)
            dvj =  q1(i,j-2,3) - 4*q1(i,j-1,3) + 3*q1(i,j,3)
          end if
          if ((i.ne.0).and.(i.ne.imax)) then
            dxi =  xc(i+1,j)   - xc(i-1,j) 
            dyi =  yc(i+1,j)   - yc(i-1,j)
            dui =  q1(i+1,j,2) - q1(i-1,j,2)
            dvi =  q1(i+1,j,3) - q1(i-1,j,3)
          end if
          if ((j.ne.0).and.(j.ne.jmax)) then
            dxj =  xc(i,j+1)   - xc(i,j-1) 
            dyj =  yc(i,j+1)   - yc(i,j-1)
            duj =  q1(i,j+1,2) - q1(i,j-1,2)
            dvj =  q1(i,j+1,3) - q1(i,j-1,3)
          end if
          q2(i,j,3)=(dyj*dvi-dyi*dvj+dxj*dui-dxi*duj)/(dxi*dyj-dyi*dxj)
        end do
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine cp_glob(q1,q2,imax,jmax,npde)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,npde
      integer(kind=cosa_int) i,j,ipde

      real(kind=cosa_real) q1(-1:imax+1,-1:jmax+1,npde), &
                   q2(-1:imax+1,-1:jmax+1,npde)

      if (calvort) then
        ipde = 4
      else
        ipde = 3
      end if

      do i=0,imax
        do j=0,jmax
          q2(i,j,ipde) = 2 * (q1(i,j,4) - 1/gamma) / machfs**2
        end do
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine wr_cfyp(ijkmaxg,cf,yplus,xsurf,ysurf,zsurf,nsurf,msurf, &
                         surf_dim,ijkmax_surf)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) ijkmaxg,nsurf,msurf,surf_dim(msurf)
      integer(kind=cosa_int) i,itime1,isurf,ijkmax_surf(msurf,3)
      real(kind=cosa_real) cf(ijkmaxg-1,msurf,ntime), &
        yplus(ijkmaxg-1,msurf,ntime),xsurf(ijkmaxg-1,msurf), &
        ysurf(ijkmaxg-1,msurf),zsurf(ijkmaxg-1,msurf)

!---- write cf files

      if (dualt) then
        if (itime.le.9) then
          write(flowtec,'(''000'',i1)') itime
        elseif (itime.le.99) then
          write(flowtec,'(''00'',i2)') itime
        elseif (itime.le.999) then
          write(flowtec,'(''0'',i3)') itime
        else
          write(flowtec,'(i4)') itime
        end if
      else if (harbal) then
        if (itime.le.9) then
          write(flowtec,'(''00'',i1)') itime
        elseif (itime.le.99) then
          write(flowtec,'(''0'',i2)') itime
        else
          write(flowtec,'(i3)') itime
        end if
      else
        flowtec = 'steady'
      end if

      flowtec = 'flow_tec_cf_'//flowtec
      flowtec = trim(flowtec)//'.dat'

!---- write header of cf tecplot file
      open(200,file=flowtec,status='replace')
      write(200,*) 'TITLE="Skin friction coefficient"'
      write(200,*) 'VARIABLES=x,y,z,cf'

      do isurf = 1,nsurf

        write(200,11) ijkmax_surf(isurf,1), ijkmax_surf(isurf,2), &
                      ijkmax_surf(isurf,3)
        do i=1,surf_dim(isurf)
          write(200,10) xsurf(i,isurf),ysurf(i,isurf),zsurf(i,isurf), &
                        cf(i,isurf,1)
        end do

      end do

      close(unit=200)    


!----- The lines below are for writing out the ascii file for each patch 
!      separately, as it is done in 2D code.
!     do isurf = 1,nsurf
!
!       if (unsteady) then
!
!         if (dualt) then
!           if (itime.le.9) then
!             write(flowtec,'(''000'',i1)') itime
!           elseif (itime.le.99) then
!             write(flowtec,'(''00'',i2)') itime
!           elseif (itime.le.999) then
!             write(flowtec,'(''0'',i3)') itime
!           else
!             write(flowtec,'(i4)') itime
!           end if
!         else if (harbal) then
!           if (itime.le.9) then
!             write(flowtec,'(''00'',i1)') itime
!           elseif (itime.le.99) then
!             write(flowtec,'(''0'',i2)') itime
!           else
!             write(flowtec,'(i3)') itime
!           end if
!         end if
!
!         if (isurf.le.9) then
!           write(filename,'(''000'',i1)') isurf
!         else if (isurf.le.99) then
!           write(filename,'(''00'',i2)') isurf
!         else if (isurf.le.9) then
!           write(filename,'(''0'',i3)') isurf
!         end if
!         flowtec = 'unst_cf_TS_'//flowtec
!         filename = trim(flowtec)//'_s_'//filename
!         filename = trim(filename)//'.dat'
!         open(8,file=filename,status='replace')
!         write(8,10) simtime
!         do i=1,surf_dim(isurf)
!           write(8,10) xsurf(i,isurf),ysurf(i,isurf),zsurf(i,isurf), &
!                       cf(i,isurf,1)
!         end do
!         close(unit=8)    
!
!       else
!
!         if (isurf.le.9) then
!           write(filename,'(''000'',i1)') isurf
!         else if (isurf.le.99) then
!           write(filename,'(''00'',i2)') isurf
!         else if (isurf.le.9) then
!           write(filename,'(''0'',i3)') isurf
!         end if
!         filename = 'stea_cf_s'//filename
!         filename = trim(filename)//'.dat'
!         open(8,file=filename,status='replace')
!         do i=1,surf_dim(isurf)
!           write(8,10) xsurf(i,isurf),ysurf(i,isurf),zsurf(i,isurf), &
!                       cf(i,isurf,1)
!         end do
!         close(unit=8)    
!
!       end if
!
!     end do
 
!---- write yplus files

      if (dualt) then
        if (itime.le.9) then
          write(flowtec,'(''000'',i1)') itime
        elseif (itime.le.99) then
          write(flowtec,'(''00'',i2)') itime
        elseif (itime.le.999) then
          write(flowtec,'(''0'',i3)') itime
        else
          write(flowtec,'(i4)') itime
        end if
      else if (harbal) then
        if (itime.le.9) then
          write(flowtec,'(''00'',i1)') itime
        elseif (itime.le.99) then
          write(flowtec,'(''0'',i2)') itime
        else
          write(flowtec,'(i3)') itime
        end if
      else
        flowtec = 'steady'
      end if

      flowtec = 'flow_tec_yp_'//flowtec
      flowtec = trim(flowtec)//'.dat'

!---- write header of yp tecplot file
      open(200,file=flowtec,status='replace')
      write(200,*) 'TITLE="Skin friction coefficient"'
      write(200,*) 'VARIABLES=x,y,z,y<sup>+'

      do isurf = 1,nsurf

        write(200,11) ijkmax_surf(isurf,1), ijkmax_surf(isurf,2), &
                      ijkmax_surf(isurf,3)
        do i=1,surf_dim(isurf)
          write(200,10) xsurf(i,isurf),ysurf(i,isurf),zsurf(i,isurf), &
                        yplus(i,isurf,1)
        end do

      end do

      close(unit=200)    


!----- The lines below are for writing out the ascii file for each patch 
!      separately, as it is done in 2D code.
!     do isurf = 1,nsurf
!
!       if (unsteady) then
!
!         if (dualt) then
!           if (itime.le.9) then
!             write(flowtec,'(''000'',i1)') itime
!           elseif (itime.le.99) then
!             write(flowtec,'(''00'',i2)') itime
!           elseif (itime.le.999) then
!             write(flowtec,'(''0'',i3)') itime
!           else
!             write(flowtec,'(i4)') itime
!           end if
!         else if (harbal) then
!           if (itime.le.9) then
!             write(flowtec,'(''00'',i1)') itime
!           elseif (itime.le.99) then
!             write(flowtec,'(''0'',i2)') itime
!           else
!             write(flowtec,'(i3)') itime
!           end if
!         end if
!
!         if (isurf.le.9) then
!           write(filename,'(''000'',i1)') isurf
!         else if (isurf.le.99) then
!           write(filename,'(''00'',i2)') isurf
!         end if
!         flowtec = 'unst_yp_TS_'//flowtec
!         filename = trim(flowtec)//'_s_'//filename
!         filename = trim(filename)//'.dat'
!         open(8,file=filename,status='replace')
!         write(8,10) simtime
!         do i=1,surf_dim(isurf)
!           write(8,10) xsurf(i,isurf),ysurf(i,isurf),zsurf(i,isurf), &
!                       yplus(i,isurf,1)
!         end do
!         close(unit=8)    
!
!       else
!
!         if (isurf.le.9) then
!           write(filename,'(''000'',i1)') isurf
!         else if (isurf.le.99) then
!           write(filename,'(''00'',i2)') isurf
!         end if
!         filename = 'stea_yp_s'//filename
!         filename = trim(filename)//'.dat'
!         open(8,file=filename,status='replace')
!         do i=1,surf_dim(isurf)
!           write(8,10) xsurf(i,isurf),ysurf(i,isurf),zsurf(i,isurf), &
!                       yplus(i,isurf,1)
!         end do
!         close(unit=8)    
!
!       end if
!
!     end do

 10   format(300e14.6)
 11   format('ZONE T="arturo",I=',i4,',J=',i4,',K=',i4,', F=POINT, DT=(S &
     &INGLE SINGLE SINGLE SINGLE)')
!ccc  1005 format (/' ',3x,'i',3x,'j',6x,'x',12x,'y',9x,'cf')
!ccc  1010 format (' ',2i4,2(e13.5),1x,e16.8)

      return
      end

!-----------------------------------------------------------------------
      subroutine vert2cell_gb_c(x,y,z,xyzc,imax,jmax,kmax,npde,nharms)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,npde,nharms,nl
      integer(kind=cosa_int) i,im,i1,j,jm,j1,k,km,k1,n,nh
      real (kind=cosa_real) &
           x     (0:imax+1 ,0:jmax+1 ,0:kmax+1      ,0:2*nharms*hbmove), &
           y     (0:imax+1 ,0:jmax+1 ,0:kmax+1      ,0:2*nharms*hbmove), &
           z     (0:imax+1 ,0:jmax+1 ,0:kmax+1      ,0:2*nharms*hbmove), &
           xyzc  (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms)

      im = imax-1
      jm = jmax-1
      km = kmax-1

      do n = 0,2*nharms
        nh = n*hbmove

!------ block vertices

        do k=0,kmax,kmax
          if (k.eq.0) then
            k1 = 1
          else if (k.eq.kmax) then
            k1 = kmax
          end if
          xyzc(0   ,0   ,k,1,n) = x(1   ,1   ,k1   ,nh)
          xyzc(0   ,0   ,k,2,n) = y(1   ,1   ,k1   ,nh)
          xyzc(0   ,0   ,k,3,n) = z(1   ,1   ,k1   ,nh)
          xyzc(imax,0   ,k,1,n) = x(imax,1   ,k1   ,nh)    
          xyzc(imax,0   ,k,2,n) = y(imax,1   ,k1   ,nh)
          xyzc(imax,0   ,k,3,n) = z(imax,1   ,k1   ,nh)
          xyzc(0   ,jmax,k,1,n) = x(1   ,jmax,k1   ,nh)
          xyzc(0   ,jmax,k,2,n) = y(1   ,jmax,k1   ,nh)
          xyzc(0   ,jmax,k,3,n) = z(1   ,jmax,k1   ,nh)
          xyzc(imax,jmax,k,1,n) = x(imax,jmax,k1   ,nh)
          xyzc(imax,jmax,k,2,n) = y(imax,jmax,k1   ,nh)
          xyzc(imax,jmax,k,3,n) = z(imax,jmax,k1   ,nh)
        end do

!------ block edges

        do k=0,kmax,kmax
          if (k.eq.0) then
            k1 = 1
          else if (k.eq.kmax) then
            k1 = kmax
          end if
          do j=0,jmax,jmax
            if (j.eq.0) then
              j1 = 1
            else if (j.eq.jmax) then
              j1 = jmax
            end if
            do i = 1,im
              xyzc(i,j,k,1,n) = (x(i,j1,k1,nh) + x(i+1,j1,k1,nh)) / 2
              xyzc(i,j,k,2,n) = (y(i,j1,k1,nh) + y(i+1,j1,k1,nh)) / 2
              xyzc(i,j,k,3,n) = (z(i,j1,k1,nh) + z(i+1,j1,k1,nh)) / 2
            end do
          end do
        end do

        do k=0,kmax,kmax
          if (k.eq.0) then
            k1 = 1
          else if (k.eq.kmax) then
            k1 = kmax
          end if
          do i=0,imax,imax
            if (i.eq.0) then
              i1 = 1
            else if (i.eq.imax) then
              i1 = imax
            end if
            do j = 1,jm
              xyzc(i,j,k,1,n) = (x(i1,j,k1,nh) + x(i1,j+1,k1,nh)) / 2
              xyzc(i,j,k,2,n) = (y(i1,j,k1,nh) + y(i1,j+1,k1,nh)) / 2
              xyzc(i,j,k,3,n) = (z(i1,j,k1,nh) + z(i1,j+1,k1,nh)) / 2
            end do
          end do
        end do

        do j=0,jmax,jmax
          if (j.eq.0) then
            j1 = 1
          else if (j.eq.jmax) then
            j1 = jmax
          end if
          do i=0,imax,imax
            if (i.eq.0) then
              i1 = 1
            else if (i.eq.imax) then
              i1 = imax
            end if
            do k = 1,km
              xyzc(i,j,k,1,n) = (x(i1,j1,k,nh) + x(i1,j1,k+1,nh)) / 2
              xyzc(i,j,k,2,n) = (y(i1,j1,k,nh) + y(i1,j1,k+1,nh)) / 2
              xyzc(i,j,k,3,n) = (z(i1,j1,k,nh) + z(i1,j1,k+1,nh)) / 2
            end do
          end do
        end do

!------ block bounding surfaces

        do k=0,kmax,kmax
          if (k.eq.0) then
            k1 = 1
          else if (k.eq.kmax) then
            k1 = kmax
          end if
          do j = 1,jm
            do i = 1,im
              xyzc(i,j,k,1,n) = &
                (x(i  ,j  ,k1,nh) + x(i+1,j  ,k1,nh) + &
                 x(i  ,j+1,k1,nh) + x(i+1,j+1,k1,nh)) / 4
              xyzc(i,j,k,2,n) = &
                (y(i  ,j  ,k1,nh) + y(i+1,j  ,k1,nh) + &
                 y(i  ,j+1,k1,nh) + y(i+1,j+1,k1,nh)) / 4
              xyzc(i,j,k,3,n) = &
                (z(i  ,j  ,k1,nh) + z(i+1,j  ,k1,nh) + &
                 z(i  ,j+1,k1,nh) + z(i+1,j+1,k1,nh)) / 4
            end do
          end do
        end do

        do j=0,jmax,jmax
          if (j.eq.0) then
            j1 = 1
          else if (j.eq.jmax) then
            j1 = jmax
          end if
          do k = 1,km
            do i = 1,im
              xyzc(i,j,k,1,n) = &
                (x(i  ,j1,k  ,nh) + x(i+1,j1,k  ,nh) + &
                 x(i  ,j1,k+1,nh) + x(i+1,j1,k+1,nh)) / 4
              xyzc(i,j,k,2,n) = &
                (y(i  ,j1,k  ,nh) + y(i+1,j1,k  ,nh) + &
                 y(i  ,j1,k+1,nh) + y(i+1,j1,k+1,nh)) / 4
              xyzc(i,j,k,3,n) = &
                (z(i  ,j1,k  ,nh) + z(i+1,j1,k  ,nh) + &
                 z(i  ,j1,k+1,nh) + z(i+1,j1,k+1,nh)) / 4
            end do
          end do
        end do

        do i=0,imax,imax
          if (i.eq.0) then
            i1 = 1
          else if (i.eq.imax) then
            i1 = imax
          end if
          do k = 1,km
            do j = 1,jm
              xyzc(i,j,k,1,n) = &
                (x(i1,j  ,k  ,nh) + x(i1,j  ,k+1,nh) + &
                 x(i1,j+1,k  ,nh) + x(i1,j+1,k+1,nh)) / 4
              xyzc(i,j,k,2,n) = &
                (y(i1,j  ,k  ,nh) + y(i1,j  ,k+1,nh) + &
                 y(i1,j+1,k  ,nh) + y(i1,j+1,k+1,nh)) / 4
              xyzc(i,j,k,3,n) = &
                (z(i1,j  ,k  ,nh) + z(i1,j  ,k+1,nh) + &
                 z(i1,j+1,k  ,nh) + z(i1,j+1,k+1,nh)) / 4
            end do
          end do
        end do

!------ internal nodes

        do k = 1,km
          do j = 1,jm
            do i = 1,im
              xyzc(i,j,k,1,n) = &
                (x(i  ,j  ,k  ,nh) + x(i+1,j  ,k  ,nh) + &
                 x(i  ,j+1,k  ,nh) + x(i+1,j+1,k  ,nh) + &
                 x(i  ,j  ,k+1,nh) + x(i+1,j  ,k+1,nh) + &
                 x(i  ,j+1,k+1,nh) + x(i+1,j+1,k+1,nh) ) / 8
              xyzc(i,j,k,2,n) = &
                (y(i  ,j  ,k  ,nh) + y(i+1,j  ,k  ,nh) + &
                 y(i  ,j+1,k  ,nh) + y(i+1,j+1,k  ,nh) + &
                 y(i  ,j  ,k+1,nh) + y(i+1,j  ,k+1,nh) + &
                 y(i  ,j+1,k+1,nh) + y(i+1,j+1,k+1,nh) ) / 8
              xyzc(i,j,k,3,n) = &
                (z(i  ,j  ,k  ,nh) + z(i+1,j  ,k  ,nh) + &
                 z(i  ,j+1,k  ,nh) + z(i+1,j+1,k  ,nh) + &
                 z(i  ,j  ,k+1,nh) + z(i+1,j  ,k+1,nh) + &
                 z(i  ,j+1,k+1,nh) + z(i+1,j+1,k+1,nh) ) / 8
            end do
          end do
        end do

      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine vert2cell_gb_v(xdot,ydot,zdot,xyzcdot,imax,jmax,kmax, &
                                npde,nharms)
!-----------------------------------------------------------------------
       
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,npde,nharms,nl
      integer(kind=cosa_int) i,im,i1,j,jm,j1,k,km,k1,n,nh
      real (kind=cosa_real) &
          xdot   (0:imax+1 ,0:jmax+1 ,0:kmax+1      ,0:2*nharms*hbmove), &
          ydot   (0:imax+1 ,0:jmax+1 ,0:kmax+1      ,0:2*nharms*hbmove), &
          zdot   (0:imax+1 ,0:jmax+1 ,0:kmax+1      ,0:2*nharms*hbmove), &
          xyzcdot(-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms)

      im = imax-1
      jm = jmax-1
      km = kmax-1

      do n = 0,2*nharms
        nh = n*hbmove

!------ block vertices

        do k=0,kmax,kmax
          if (k.eq.0) then
            k1 = 1
          else if (k.eq.kmax) then
            k1 = kmax
          end if
          xyzcdot(0   ,0   ,k,1,n) = xdot(1   ,1   ,k1   ,nh)
          xyzcdot(0   ,0   ,k,2,n) = ydot(1   ,1   ,k1   ,nh)
          xyzcdot(0   ,0   ,k,3,n) = zdot(1   ,1   ,k1   ,nh)
          xyzcdot(imax,0   ,k,1,n) = xdot(imax,1   ,k1   ,nh)    
          xyzcdot(imax,0   ,k,2,n) = ydot(imax,1   ,k1   ,nh)
          xyzcdot(imax,0   ,k,3,n) = zdot(imax,1   ,k1   ,nh)
          xyzcdot(0   ,jmax,k,1,n) = xdot(1   ,jmax,k1   ,nh)
          xyzcdot(0   ,jmax,k,2,n) = ydot(1   ,jmax,k1   ,nh)
          xyzcdot(0   ,jmax,k,3,n) = zdot(1   ,jmax,k1   ,nh)
          xyzcdot(imax,jmax,k,1,n) = xdot(imax,jmax,k1   ,nh)
          xyzcdot(imax,jmax,k,2,n) = ydot(imax,jmax,k1   ,nh)
          xyzcdot(imax,jmax,k,3,n) = zdot(imax,jmax,k1   ,nh)
        end do

!------ block edges

        do k=0,kmax,kmax
          if (k.eq.0) then
            k1 = 1
          else if (k.eq.kmax) then
            k1 = kmax
          end if
          do j=0,jmax,jmax
            if (j.eq.0) then
              j1 = 1
            else if (j.eq.jmax) then
              j1 = jmax
            end if
            do i = 1,im
              xyzcdot(i,j,k,1,n) = &
                     (xdot(i,j1,k1,nh) + xdot(i+1,j1,k1,nh)) / 2
              xyzcdot(i,j,k,2,n) = &
                     (ydot(i,j1,k1,nh) + ydot(i+1,j1,k1,nh)) / 2
              xyzcdot(i,j,k,3,n) = &
                     (zdot(i,j1,k1,nh) + zdot(i+1,j1,k1,nh)) / 2
            end do
          end do
        end do

        do k=0,kmax,kmax
          if (k.eq.0) then
            k1 = 1
          else if (k.eq.kmax) then
            k1 = kmax
          end if
          do i=0,imax,imax
            if (i.eq.0) then
              i1 = 1
            else if (i.eq.imax) then
              i1 = imax
            end if
            do j = 1,jm
              xyzcdot(i,j,k,1,n) = &
                     (xdot(i1,j,k1,nh) + xdot(i1,j+1,k1,nh)) / 2
              xyzcdot(i,j,k,2,n) = &
                     (ydot(i1,j,k1,nh) + ydot(i1,j+1,k1,nh)) / 2
              xyzcdot(i,j,k,3,n) = &
                     (zdot(i1,j,k1,nh) + zdot(i1,j+1,k1,nh)) / 2
            end do
          end do
        end do

        do j=0,jmax,jmax
          if (j.eq.0) then
            j1 = 1
          else if (j.eq.jmax) then
            j1 = jmax
          end if
          do i=0,imax,imax
            if (i.eq.0) then
              i1 = 1
            else if (i.eq.imax) then
              i1 = imax
            end if
            do k = 1,km
              xyzcdot(i,j,k,1,n) = &
                     (xdot(i1,j1,k,nh) + xdot(i1,j1,k+1,nh)) / 2
              xyzcdot(i,j,k,2,n) = &
                     (ydot(i1,j1,k,nh) + ydot(i1,j1,k+1,nh)) / 2
              xyzcdot(i,j,k,3,n) = &
                     (zdot(i1,j1,k,nh) + zdot(i1,j1,k+1,nh)) / 2
            end do
          end do
        end do

!------ block bounding surfaces

        do k=0,kmax,kmax
          if (k.eq.0) then
            k1 = 1
          else if (k.eq.kmax) then
            k1 = kmax
          end if
          do j = 1,jm
            do i = 1,im
              xyzcdot(i,j,k,1,n) = &
                (xdot(i  ,j  ,k1,nh) + xdot(i+1,j  ,k1,nh) + &
                 xdot(i  ,j+1,k1,nh) + xdot(i+1,j+1,k1,nh)) / 4
              xyzcdot(i,j,k,2,n) = &
                (ydot(i  ,j  ,k1,nh) + ydot(i+1,j  ,k1,nh) + &
                 ydot(i  ,j+1,k1,nh) + ydot(i+1,j+1,k1,nh)) / 4
              xyzcdot(i,j,k,3,n) = &
                (zdot(i  ,j  ,k1,nh) + zdot(i+1,j  ,k1,nh) + &
                 zdot(i  ,j+1,k1,nh) + zdot(i+1,j+1,k1,nh)) / 4
            end do
          end do
        end do

        do j=0,jmax,jmax
          if (j.eq.0) then
            j1 = 1
          else if (j.eq.jmax) then
            j1 = jmax
          end if
          do k = 1,km
            do i = 1,im
              xyzcdot(i,j,k,1,n) = &
                (xdot(i  ,j1,k  ,nh) + xdot(i+1,j1,k  ,nh) + &
                 xdot(i  ,j1,k+1,nh) + xdot(i+1,j1,k+1,nh)) / 4
              xyzcdot(i,j,k,2,n) = &
                (ydot(i  ,j1,k  ,nh) + ydot(i+1,j1,k  ,nh) + &
                 ydot(i  ,j1,k+1,nh) + ydot(i+1,j1,k+1,nh)) / 4
              xyzcdot(i,j,k,3,n) = &
                (zdot(i  ,j1,k  ,nh) + zdot(i+1,j1,k  ,nh) + &
                 zdot(i  ,j1,k+1,nh) + zdot(i+1,j1,k+1,nh)) / 4
            end do
          end do
        end do

        do i=0,imax,imax
          if (i.eq.0) then
            i1 = 1
          else if (i.eq.imax) then
            i1 = imax
          end if
          do k = 1,km
            do j = 1,jm
              xyzcdot(i,j,k,1,n) = &
                (xdot(i1,j  ,k  ,nh) + xdot(i1,j  ,k+1,nh) + &
                 xdot(i1,j+1,k  ,nh) + xdot(i1,j+1,k+1,nh)) / 4
              xyzcdot(i,j,k,2,n) = &
                (ydot(i1,j  ,k  ,nh) + ydot(i1,j  ,k+1,nh) + &
                 ydot(i1,j+1,k  ,nh) + ydot(i1,j+1,k+1,nh)) / 4
              xyzcdot(i,j,k,3,n) = &
                (zdot(i1,j  ,k  ,nh) + zdot(i1,j  ,k+1,nh) + &
                 zdot(i1,j+1,k  ,nh) + zdot(i1,j+1,k+1,nh)) / 4
            end do
          end do
        end do

!------ internal nodes

        do k = 1,km
          do j = 1,jm
            do i = 1,im
              xyzcdot(i,j,k,1,n) = &
                (xdot(i  ,j  ,k  ,nh) + xdot(i+1,j  ,k  ,nh) + &
                 xdot(i  ,j+1,k  ,nh) + xdot(i+1,j+1,k  ,nh) + &
                 xdot(i  ,j  ,k+1,nh) + xdot(i+1,j  ,k+1,nh) + &
                 xdot(i  ,j+1,k+1,nh) + xdot(i+1,j+1,k+1,nh) ) / 8
              xyzcdot(i,j,k,2,n) = &
                (ydot(i  ,j  ,k  ,nh) + ydot(i+1,j  ,k  ,nh) + &
                 ydot(i  ,j+1,k  ,nh) + ydot(i+1,j+1,k  ,nh) + &
                 ydot(i  ,j  ,k+1,nh) + ydot(i+1,j  ,k+1,nh) + &
                 ydot(i  ,j+1,k+1,nh) + ydot(i+1,j+1,k+1,nh) ) / 8
              xyzcdot(i,j,k,3,n) = &
                (zdot(i  ,j  ,k  ,nh) + zdot(i+1,j  ,k  ,nh) + &
                 zdot(i  ,j+1,k  ,nh) + zdot(i+1,j+1,k  ,nh) + &
                 zdot(i  ,j  ,k+1,nh) + zdot(i+1,j  ,k+1,nh) + &
                 zdot(i  ,j+1,k+1,nh) + zdot(i+1,j+1,k+1,nh) ) / 8
            end do
          end do
        end do

      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine writec(imax,jmax,npde,xc,yc,q1,q2,tempus)
!-----------------------------------------------------------------------

!mpi THIS WILL NOT WORK WITH MPI BUT ONLY SEEMS TO BE CALLED FROM POSCOSA
!mpi AT THE MOMENT SO SHOULD NOT BE A PROBLEM

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,npde
      integer(kind=cosa_int) i,j,ipde,npde2
      real(kind=cosa_real) xc(0:imax,0:jmax),yc(0:imax,0:jmax), &
          q1(-1:imax+1,-1:jmax+1,npde),q2(-1:imax+1,-1:jmax+1,npde)
      real(kind=cosa_real) tempus

!---- write tecplot file  
      open(60,file=flowtec,status='replace')

      write(60,*) 'TITLE="meshflow"'

      if (calvort.and.(.not.calcp)) then
        write(60,*) 'VARIABLES=xc,yc,rho,u,v,p,t,mach,vort'
        npde2 = 3
      else if ((.not.calvort).and.(calcp)) then
        write(60,*) 'VARIABLES=xc,yc,rho,u,v,p,t,mach,cp'
        npde2 = 3
      else if ((calvort).and.(calcp)) then
        write(60,*) 'VARIABLES=xc,yc,rho,u,v,p,t,mach,vort,cp'
        npde2 = 4
      end if

      write(60,11) imax+1,jmax+1
!     write(line,'(''ZONE T="arturo",I='',i4,'', J='',i4,'',F=POINT, DT= &
!    &(SINGLE SINGLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)'')') &
!           imax1,jmax1
!     write(60,'(a)') line
      do j=0,jmax
        do i=0,imax
          write(60,8) xc(i,j), yc(i,j), (q1(i,j,ipde),ipde=1,4), &
               (q2(i,j,ipde),ipde=1,npde2)
        end do
      end do

      if (unsteady) then
        write(60,15) tempus 
      end if

 8    format(2e16.8,6e22.14)
 11   format('ZONE T="arturo",I=',i4,', J=',i4,',F=POINT, DT=(SINGLE SIN &
     &GLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)')
 15   format(e21.14)

      close(unit=60)

      return
      end

!-----------------------------------------------------------------------
      subroutine forces(nl,ijkmaxg,nsurf,msurf,q1,q2,bctopo,surf_dim, &
                        si,sj,sk,xideri,xiderj,xiderk,etaderi,etaderj, &
                        etaderk,zetaderi,zetaderj,zetaderk,x,y,z, &
                        xyzcdot,cl_p,cl_v,cd_p,cd_v,cz_p,cz_v,cmx_p, &
                        cmx_v,cmy_p,cmy_v,cm_p,cm_v)
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,iblock,iblk_pos,imax,jmax,kmax,ijkmaxg,iq,iblk,nsurf, &
                msurf,surf_dim(msurf),nsurf1,iimt,ivmt,ixyz,ixyzc1
      real (kind=cosa_real) x(*),y(*),z(*),xyzcdot(*),si(*),sj(*),sk(*), &
                    xideri(*),xiderj(*),xiderk(*),etaderi(*),etaderj(*), &
                    etaderk(*),zetaderi(*),zetaderj(*),zetaderk(*), &
                    q1(*),q2(*)
      real(kind=cosa_real) &
                   cl_p (ijkmaxg-1,msurf),cl_v (ijkmaxg-1,msurf), &
                   cd_p (ijkmaxg-1,msurf),cd_v (ijkmaxg-1,msurf), &
                   cz_p (ijkmaxg-1,msurf),cz_v (ijkmaxg-1,msurf), &
                   cmx_p(ijkmaxg-1,msurf),cmx_v(ijkmaxg-1,msurf), &
                   cmy_p(ijkmaxg-1,msurf),cmy_v(ijkmaxg-1,msurf), &
                   cm_p (ijkmaxg-1,msurf),cm_v (ijkmaxg-1,msurf)
      integer(kind=cosa_int) bctopo(*)

      nsurf1 = 0
      do iblk_pos = 1,nblock_pos
        
        iblock   = iblock_pos(iblk_pos)
        imax     = i_imax(iblock,nl)
        jmax     = j_jmax(iblock,nl)
        kmax     = k_kmax(iblock,nl)
        iq       = 1 + off_p3 (iblock,nl) * npde
        iblk     = 1 + off_bct(iblock,nl)
        iimt     = 1 + off_p2 (iblock,nl) * lmet * dim5h
        ivmt     = 1 + off_0  (iblock,nl) * 3    * dim5h
        ixyz     = 1 + off_p2 (iblock,nl)        * dim5h
        ixyzc1   = 1 + off_p3 (iblock,nl) * npde * dim5

        call forces_b(imax,jmax,kmax,ijkmaxg,npde,lmet,nharms,nsurf1, &
                      msurf,nbcs(iblock),surf_dim,q1(iq),q2(iq), &
                      bctopo(iblk), &
                      si(iimt),sj(iimt),sk(iimt),xideri(ivmt), &
                      xiderj(ivmt),xiderk(ivmt),etaderi(ivmt), &
                      etaderj(ivmt),etaderk(ivmt),zetaderi(ivmt), &
                      zetaderj(ivmt),zetaderk(ivmt),x(ixyz),y(ixyz), &
                      z(ixyz),xyzcdot(ixyzc1), &
                      cl_p,cl_v,cd_p,cd_v,cz_p,cz_v,cmx_p,cmx_v,cmy_p, &
                      cmy_v,cm_p,cm_v)

      end do

      if (nsurf1.ne.nsurf) then
        write(*,*) 'Incongruency in cp_wall. Aborting!'
        stop
      end if

      return
      end
!-----------------------------------------------------------------------
      subroutine forces_b(imax,jmax,kmax,ijkmaxg,npde,lmet,nharms, &
                          nsurf1,msurf,nbcs,surf_dim,q1,q2,bctopo,si,sj, &
                          sk,xideri, &
                          xiderj,xiderk,etaderi,etaderj,etaderk, &
                          zetaderi,zetaderj,zetaderk,x,y,z,xyzcdot, &
                          cl_p,cl_v,cd_p,cd_v,cz_p,cz_v,cmx_p,cmx_v, &
                          cmy_p,cmy_v,cm_p,cm_v)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,ijkmaxg,npde,nharms,nbcs,lmet
      integer(kind=cosa_int) bctopo(10,nbcs)
      integer(kind=cosa_int) n,nh,msurf,surf_dim(msurf),nsurf1
      real (kind=cosa_real) &
           x       (0:imax+1 ,0:jmax+1 ,0:kmax+1   ,0:2*nharms*hbmove), &
           y       (0:imax+1 ,0:jmax+1 ,0:kmax+1   ,0:2*nharms*hbmove), &
           z       (0:imax+1 ,0:jmax+1 ,0:kmax+1   ,0:2*nharms*hbmove), &
           xyzcdot (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms    ), &
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
           q1(-1:imax+1,-1:jmax+1,-1:kmax+1,npde), &
           q2(-1:imax+1,-1:jmax+1,-1:kmax+1,npde)

      real(kind=cosa_real) itime1, &
                   cl_p (ijkmaxg-1,msurf),cl_v (ijkmaxg-1,msurf), &
                   cd_p (ijkmaxg-1,msurf),cd_v (ijkmaxg-1,msurf), &
                   cz_p (ijkmaxg-1,msurf),cz_v (ijkmaxg-1,msurf), &
                   cmx_p(ijkmaxg-1,msurf),cmx_v(ijkmaxg-1,msurf), &
                   cmy_p(ijkmaxg-1,msurf),cmy_v(ijkmaxg-1,msurf), &
                   cm_p (ijkmaxg-1,msurf),cm_v (ijkmaxg-1,msurf)

      if ((.not.unsteady).or.dualt) then
        nh = 0
      else if (harbal) then
        n  = itime1 - 1
        nh = n*hbmove
      end if
      call forces_b_h(imax,jmax,kmax,ijkmaxg,npde,lmet,nsurf1, &
                      msurf,nbcs,surf_dim,q1,q2,bctopo,si(1,0,0,0,nh), &
                      sj(1,0,0,0,nh),sk(1,0,0,0,nh), &
                      xideri(1,1,1,1,nh),xiderj(1,1,1,1,nh), &
                      xiderk(1,1,1,1,nh),etaderi(1,1,1,1,nh), &
                      etaderj(1,1,1,1,nh),etaderk(1,1,1,1,nh), &
                      zetaderi(1,1,1,1,nh),zetaderj(1,1,1,1,nh), &
                      zetaderk(1,1,1,1,nh),x(0,0,0,nh),y(0,0,0,nh), &
                      z(0,0,0,nh),xyzcdot(-1,-1,-1,1,nh), &
                      cl_p,cl_v,cd_p,cd_v,cz_p,cz_v,cmx_p,cmx_v,cmy_p, &
                      cmy_v,cm_p,cm_v,nh)

      return
      end

!-----------------------------------------------------------------------
      subroutine forces_b_h(imax,jmax,kmax,ijkmaxg,npde,lmet,nsurf1, &
                            msurf,nbcs,surf_dim,q1,q2,bctopo,si,sj,sk, &
                            xideri,xiderj,xiderk,etaderi,etaderj,etaderk, &
                            zetaderi,zetaderj,zetaderk,x,y,z,xyzcdot, &
                            cl_p,cl_v,cd_p,cd_v,cz_p,cz_v,cmx_p,cmx_v, &
                            cmy_p,cmy_v,cm_p,cm_v,n)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,ijkmaxg,nbcs,npde,nsurf1,msurf, &
                surf_dim(msurf),lmet,n
!     integer(kind=cosa_int) i1,i2,i3,ibbc,sur_dim,ijkmax(3),bctyp,idir,inrout, &
!               istrt(3),iend(3),ibcpt,ic1,ic2,ic3,ibc,jbc,kbc
      integer(kind=cosa_int) i1,i2,i3,ibbc,sur_dim,ijkmax(3),idir,istrt(3),iend(3), &
                bctyp,bctopo1, &
                inrout,ibcpt,ibcpt2,ibcn,ibcn2,ibcm,ibcw,ibc,jbc,kbc, &
                ibc2,jbc2,kbc2,in,jn,kn,in2,jn2,kn2,im,jm,km,ic1,ic2, &
                ic3,sysize,ioff1,ioff2,joff1,joff2,koff1,koff2,iw,jw,kw
      real(kind=cosa_real) &
           q1(-1:imax+1,-1:jmax+1,-1:kmax+1,npde), &
           q2(-1:imax+1,-1:jmax+1,-1:kmax+1,npde), &
           x(0:imax+1,0:jmax+1,0:kmax+1), &
           y(0:imax+1,0:jmax+1,0:kmax+1), &
           z(0:imax+1,0:jmax+1,0:kmax+1), &
           xyzcdot (-1:imax+1,-1:jmax+1,-1:kmax+1,npde), &
           si      (lmet,0:imax+1,0:jmax+1,0:kmax+1), &
           sj      (lmet,0:imax+1,0:jmax+1,0:kmax+1), &
           sk      (lmet,0:imax+1,0:jmax+1,0:kmax+1), &
           xideri  (3   ,imax    ,jmax    ,kmax    ), &
           xiderj  (3   ,imax    ,jmax    ,kmax    ), &
           xiderk  (3   ,imax    ,jmax    ,kmax    ), &
           etaderi (3   ,imax    ,jmax    ,kmax    ), &
           etaderj (3   ,imax    ,jmax    ,kmax    ), &
           etaderk (3   ,imax    ,jmax    ,kmax    ), &
           zetaderi(3   ,imax    ,jmax    ,kmax    ), &
           zetaderj(3   ,imax    ,jmax    ,kmax    ), &
           zetaderk(3   ,imax    ,jmax    ,kmax    )
      integer(kind=cosa_int) bctopo(10,nbcs)

      real(kind=cosa_real) fp(6),fv(6), &
           cl_p (ijkmaxg-1,msurf), cd_p (ijkmaxg-1,msurf), &
           cz_p (ijkmaxg-1,msurf), cmx_p(ijkmaxg-1,msurf), &
           cmy_p(ijkmaxg-1,msurf), cm_p (ijkmaxg-1,msurf), &
           cl_v (ijkmaxg-1,msurf), cd_v (ijkmaxg-1,msurf), &
           cz_v (ijkmaxg-1,msurf), cmx_v(ijkmaxg-1,msurf), &
           cmy_v(ijkmaxg-1,msurf), cm_v (ijkmaxg-1,msurf), &
           den,ssx,ssy,ssz,aax,aay,aaz,bbx,bby,bbz, &
           nx,ny,nz,kx,ky,kz,ds,sgnm,xarm,yarm,zarm,xmomr,ymomr,zmomr, &
           dr1(2),dr2(2),rhow,pw,tw,muw,u1,u2,uw,u1r,u2r,uwr,v1,v2,vw, &
           v1r,v2r,vwr,w1,w2,ww,w1r,w2r,wwr,ueta,veta,weta, &
           dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,divv, &
           txx,txy,txz,tyy,tyz,tzz,mach,dudxi,dudeta,dudzeta,dvdxi, &
           dvdeta, dvdzeta, dwdxi,dwdeta,dwdzeta,xix,xiy,xiz,etax, &
           etay,etaz,zetax,zetay,zetaz 

!---- mach is defined for force coefficient normalization. For viscous
      if (machfs.gt.1d-14) then
        mach = machfs
      else
        mach = 1.d0
      end if

      if (.not.harbal) zeit(0)=simtime

      if (aircraft) then
        ssx = dcos(alpha) * dcos(beta)
        ssy = dsin(alpha) * dcos(beta)
        ssz = dsin(beta)
        den = dsqrt(dsin(beta)**2 + dcos(alpha)**2 * dcos(beta)**2)
        aax = - dsin(alpha) * dcos(alpha) * dcos(beta)**2 / den
        aay = den
        aaz = dsin(alpha) * dsin(beta) * dcos(beta) / den
        bbx = -dsin(beta) / den
        bby = 0.d0
        bbz = dcos(alpha) * dcos(beta) / den
      else if (hawt.or.vawt) then
        ssx = 0.d0
        ssy = 1.d0
        ssz = 0.d0
        aax = 1.d0
        aay = 0.d0
        aaz = 0.d0
        bbx = 0.d0
        bby = 0.d0
        bbz = 1.d0
      end if

      do i1=1,6
        fp(i1)=0.d0
        fv(i1)=0.d0
      end do

      if (moving) then
        if (pitching) then
          dr2(1) =  (xmom - xh) * (dcos(dtheta(n)) - 1) - &
                    (ymom - yh) *  dsin(dtheta(n))
          dr2(2) =  (xmom - xh) *  dsin(dtheta(n))      + &
                    (ymom - yh) * (dcos(dtheta(n)) - 1)
          xmomr  =  xmom + dr2(1)
          ymomr  =  ymom + dr2(2)
          zmomr  =  zmom
        else if (plunging) then
          dr1(1) = dh0x * dsin(omega*zeit(n))
          dr1(2) = dh0y * dsin(omega*zeit(n))
          xmomr  = xmom + dr1(1)
          ymomr  = ymom + dr1(2)
          zmomr  =  zmom
        else if (plupitching) then
          dr1(1) = dh0x * dsin(omega*zeit(n))
          dr1(2) = dh0y * dsin(omega*zeit(n))
          dr2(1) =  (xmom - xh) * (dcos(dtheta(n)) - 1) - &
                    (ymom - yh) *  dsin(dtheta(n))
          dr2(2) =  (xmom - xh) *  dsin(dtheta(n))      + &
                    (ymom - yh) * (dcos(dtheta(n)) - 1)
          xmomr  = xmom + dr1(1) + dr2(1)
          ymomr  = ymom + dr1(2) + dr2(2)
          zmomr  =  zmom
        else if (rotating) then
          dr1(1) = ( (xmom - xrotc) * (dcos(dtheta(n)) - 1) - &
                     (ymom - yrotc) *  dsin(dtheta(n))      )
          dr1(2) = ( (xmom - xrotc) *  dsin(dtheta(n))      + &
                     (ymom - yrotc) * (dcos(dtheta(n)) - 1) )
          xmomr = xmom + dr1(1)
          ymomr = ymom + dr1(2)
          zmomr  =  zmom
        end if
      else
        xmomr = xmom
        ymomr = ymom
        zmomr = zmom
      end if

      do ibbc=1,nbcs

        bctopo1 = bctopo(1,ibbc)
        if (bctopo1.eq.14) bctopo1=bctopo1*100
        if (bctopo1/100.eq.wallmask) then
!       if (bctopo(1,ibbc)/100.eq.wallmask) then

          sur_dim = 0
          nsurf1 = nsurf1 + 1

!         store imax, jmax in ijmax to enable boundary data location

          ijkmax(1) = imax
          ijkmax(2) = jmax
          ijkmax(3) = kmax

!         store boundary topology in mnemonic names

          bctyp    = bctopo(1,ibbc)
          idir     = bctopo(2,ibbc)
          inrout   = bctopo(3,ibbc)
          istrt(1) = bctopo(4,ibbc)
          iend(1)  = bctopo(5,ibbc)
          istrt(2) = bctopo(6,ibbc)
          iend(2)  = bctopo(7,ibbc)
          istrt(3) = bctopo(8,ibbc)
          iend(3)  = bctopo(9,ibbc)

!         modify beginning, ending indices to extend boundary condition
!         to edge/corner

!tmp      do i1 = 1,2
!tmp        if (i1.ne.idir) then
!tmp          if (istrt(i1) .eq. 1          ) istrt(i1) = 0
!tmp          if (iend (i1) .eq. ijmax(i1)-1) iend (i1) = ijmax(i1)
!tmp        end if
!tmp      end do

!         set needed variables depending on whether the boundary is the
!         the inner boundary (inrout = 1) or 
!         the outer boundary (inrout > 1)
!             ibcpt : boundary condition location (first aux. cell)
!             ibcpt2: boundary condition location outside the block from
!                     ibcpt (second aux. cell)
!             ibcn  : point to the inside of the block from ibcpt
!             ibcm  : location of the metrics
          if (inrout.eq.1) then
            ibcpt  =  0
            ibcpt2 = -1
            ibcn   =  1
            ibcn2  =  2
            ibcm   =  1
            ibcw   =  0
            sgnm   = - 1.d0
          else
            ibcpt  = ijkmax(idir)
            ibcpt2 = ijkmax(idir) + 1
            ibcn   = ijkmax(idir) - 1
            ibcn2  = ijkmax(idir) - 2
            ibcm   = ijkmax(idir)
            ibcw   = ijkmax(idir)
            sgnm   = 1.d0
          end if

          if (idir.eq.1) then
            ioff1 = (-1)**(1+1/inrout)
            ioff2 = 2*ioff1
            joff1 = 0
            joff2 = 0
            koff1 = 0
            koff2 = 0
          else if (idir.eq.2) then
            ioff1 = 0
            ioff2 = 0
            joff1 = (-1)**(1+1/inrout)
            joff2 = 2*joff1
            koff1 = 0
            koff2 = 0
          else if (idir.eq.3) then
            ioff1 = 0
            ioff2 = 0
            joff1 = 0
            joff2 = 0
            koff1 = (-1)**(1+1/inrout)
            koff2 = 2*koff1
          end if

          ic1 = cyc (idir, 1)
          ic2 = cyc (idir, 2)
          ic3 = cyc (idir, 3)

!-------- calculate pressure force
          do i3 = istrt(ic3),iend(ic3)
            do i2 = istrt(ic2),iend(ic2)

              ibc  = ibcpt *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
              jbc  = ibcpt *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
              kbc  = ibcpt *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

              ibc2 = ibcpt2*krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
              jbc2 = ibcpt2*krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
              kbc2 = ibcpt2*krd(Ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

              in   = ibcn  *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
              jn   = ibcn  *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
              kn   = ibcn  *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

              in2  = ibcn2 *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
              jn2  = ibcn2 *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
              kn2  = ibcn2 *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

              im   = ibcm  *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
              jm   = ibcm  *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
              km   = ibcm  *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

              iw   = ibcw  *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
              jw   = ibcw  *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
              kw   = ibcw  *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

              if (idir.eq.1) then
                nx    = si(1,im,jm,km)
                ny    = si(2,im,jm,km)
                nz    = si(3,im,jm,km)
                ds    = si(4,im,jm,km)
              else if (idir.eq.2) then
                nx    = sj(1,im,jm,km)
                ny    = sj(2,im,jm,km)
                nz    = sj(3,im,jm,km)
                ds    = sj(4,im,jm,km)
              else if (idir.eq.3) then
                nx    = sk(1,im,jm,km)
                ny    = sk(2,im,jm,km)
                nz    = sk(3,im,jm,km)
                ds    = sk(4,im,jm,km)
              end if

              kx    = (-1)**(1+1/inrout) * nx
              ky    = (-1)**(1+1/inrout) * ny
              kz    = (-1)**(1+1/inrout) * nz
              pw    = q1(ibc,jbc,kbc,5)
              fp(1) = - pw * kx * ds 
              fp(2) = - pw * ky * ds
              fp(3) = - pw * kz * ds
              xarm  = x(iw,jw,kw) - xmomr
              yarm  = y(iw,jw,kw) - ymomr
              zarm  = z(iw,jw,kw) - zmomr
              fp(4) = (ky * zarm - kz *yarm) * pw * ds
              fp(5) = (kz * xarm - kx *zarm) * pw * ds
              fp(6) = (kx * yarm - ky *xarm) * pw * ds

              sur_dim = sur_dim + 1

              cl_p(sur_dim,nsurf1) = (fp(1)*aax+fp(2)*aay+fp(3)*aaz) * &
                                     2/(lngth(1)*lngth(2)*mach**2)
              cd_p(sur_dim,nsurf1) = (fp(1)*ssx+fp(2)*ssy+fp(3)*ssz) * &
                                     2/(lngth(1)*lngth(2)*mach**2)
              cz_p(sur_dim,nsurf1) = (fp(1)*bbx+fp(2)*bby+fp(3)*bbz) * &
                                     2/(lngth(1)*lngth(2)*mach**2)
              cmx_p(sur_dim,nsurf1)= fp(4) * &
                                     2/(lngth(1)*lngth(2)*lngth(3)*mach**2)
              cmy_p(sur_dim,nsurf1)= fp(5) * &
                                     2/(lngth(1)*lngth(2)*lngth(3)*mach**2)
              cm_p(sur_dim,nsurf1) = fp(6) * &
                                     2/(lngth(1)*lngth(2)*lngth(3)*mach**2)

            end do
          end do

          if ((viscous).and.(bctopo(1,ibbc)/100.eq.15)) then
!---------- calculate viscous force
          sur_dim = 0

            do i3 = istrt(ic3),iend(ic3)
              do i2 = istrt(ic2),iend(ic2)

                ibc  = ibcpt *krd(ic1,1) +i2*krd(ic2,1) +i3*krd(ic3,1)
                jbc  = ibcpt *krd(ic1,2) +i2*krd(ic2,2) +i3*krd(ic3,2)
                kbc  = ibcpt *krd(ic1,3) +i2*krd(ic2,3) +i3*krd(ic3,3)

                ibc2 = ibcpt2*krd(ic1,1) +i2*krd(ic2,1) +i3*krd(ic3,1)
                jbc2 = ibcpt2*krd(ic1,2) +i2*krd(ic2,2) +i3*krd(ic3,2)
                kbc2 = ibcpt2*krd(Ic1,3) +i2*krd(ic2,3) +i3*krd(ic3,3)

                in   = ibcn  *krd(ic1,1) +i2*krd(ic2,1) +i3*krd(ic3,1)
                jn   = ibcn  *krd(ic1,2) +i2*krd(ic2,2) +i3*krd(ic3,2)
                kn   = ibcn  *krd(ic1,3) +i2*krd(ic2,3) +i3*krd(ic3,3)

                in2  = ibcn2 *krd(ic1,1) +i2*krd(ic2,1) +i3*krd(ic3,1)
                jn2  = ibcn2 *krd(ic1,2) +i2*krd(ic2,2) +i3*krd(ic3,2)
                kn2  = ibcn2 *krd(ic1,3) +i2*krd(ic2,3) +i3*krd(ic3,3)

                im   = ibcm  *krd(ic1,1) +i2*krd(ic2,1) +i3*krd(ic3,1)
                jm   = ibcm  *krd(ic1,2) +i2*krd(ic2,2) +i3*krd(ic3,2)
                km   = ibcm  *krd(ic1,3) +i2*krd(ic2,3) +i3*krd(ic3,3)

                iw   = ibcw  *krd(ic1,1) +i2*krd(ic2,1) +i3*krd(ic3,1)
                jw   = ibcw  *krd(ic1,2) +i2*krd(ic2,2) +i3*krd(ic3,2)
                kw   = ibcw  *krd(ic1,3) +i2*krd(ic2,3) +i3*krd(ic3,3)

                if (idir.eq.1) then
                  nx    = si(1,im,jm,km)
                  ny    = si(2,im,jm,km)
                  nz    = si(3,im,jm,km)
                  ds    = si(4,im,jm,km)
                else if (idir.eq.2) then
                  nx    = sj(1,im,jm,km)
                  ny    = sj(2,im,jm,km)
                  nz    = sj(3,im,jm,km)
                  ds    = sj(4,im,jm,km)
                else if (idir.eq.3) then
                  nx    = sk(1,im,jm,km)
                  ny    = sk(2,im,jm,km)
                  nz    = sk(3,im,jm,km)
                  ds    = sk(4,im,jm,km)
                end if

                kx   = (-1)**(1+1/inrout) * nx
                ky   = (-1)**(1+1/inrout) * ny
                kz   = (-1)**(1+1/inrout) * nz
                rhow = q1(ibc,jbc,kbc,1)
                pw   = q1(ibc,jbc,kbc,5)
                tw   = q2(ibc,jbc,kbc,1)
                muw  = (stemp+1)/(tw+stemp) * (tw**1.5d0)

                if (idir.eq.1) then
                  ioff1 = (-1)**(1+1/inrout)
                  ioff2 = 2*ioff1
                  joff1 = 0
                  joff2 = 0
                  koff1 = 0
                  koff2 = 0
                else if (idir.eq.2) then
                  ioff1 = 0
                  ioff2 = 0
                  joff1 = (-1)**(1+1/inrout)
                  joff2 = 2*joff1
                  koff1 = 0
                  koff2 = 0
                else if (idir.eq.3) then
                  ioff1 = 0
                  ioff2 = 0
                  joff1 = 0
                  joff2 = 0
                  koff1 = (-1)**(1+1/inrout)
                  koff2 = 2*koff1
                end if
   
                if (idir.eq.1) then
   
                  dudxi = (-1)**(1+1/inrout) * &
                          (-8*q1(ibc,jbc,kbc,2) + &
                            9*q1(ibc+ioff1,jbc,kbc,2) - &
                              q1(ibc+ioff2,jbc,kbc,2) ) / 3
                  dvdxi = (-1)**(1+1/inrout) * &
                          (-8*q1(ibc,jbc,kbc,3) + &
                            9*q1(ibc+ioff1,jbc,kbc,3) - &
                              q1(ibc+ioff2,jbc,kbc,3) ) / 3
                  dwdxi = (-1)**(1+1/inrout) * &
                          (-8*q1(ibc,jbc,kbc,4) + &
                            9*q1(ibc+ioff1,jbc,kbc,4) - &
                              q1(ibc+ioff2,jbc,kbc,4) ) / 3

                  if (i2.eq.istrt(ic2)) then
                    u2 = (  q1(ibc,jbc,kbc,2) + q1(ibc,jbc+1,kbc,2)) / 2
                    u1 = (3*q1(ibc,jbc,kbc,2) - q1(ibc,jbc+1,kbc,2)) / 2
                    v2 = (  q1(ibc,jbc,kbc,3) + q1(ibc,jbc+1,kbc,3)) / 2
                    v1 = (3*q1(ibc,jbc,kbc,3) - q1(ibc,jbc+1,kbc,3)) / 2
                    w2 = (  q1(ibc,jbc,kbc,4) + q1(ibc,jbc+1,kbc,4)) / 2
                    w1 = (3*q1(ibc,jbc,kbc,4) - q1(ibc,jbc+1,kbc,4)) / 2
                  else if (i2.eq.iend(ic2)) then
                    u2 = (3*q1(ibc,jbc,kbc,2) - q1(ibc,jbc-1,kbc,2)) / 2
                    u1 = (  q1(ibc,jbc,kbc,2) + q1(ibc,jbc-1,kbc,2)) / 2
                    v2 = (3*q1(ibc,jbc,kbc,3) - q1(ibc,jbc-1,kbc,3)) / 2
                    v1 = (  q1(ibc,jbc,kbc,3) + q1(ibc,jbc-1,kbc,3)) / 2
                    w2 = (3*q1(ibc,jbc,kbc,4) - q1(ibc,jbc-1,kbc,4)) / 2
                    w1 = (  q1(ibc,jbc,kbc,4) + q1(ibc,jbc-1,kbc,4)) / 2
                  else
                    u2 = (q1(ibc,jbc  ,kbc,2) + q1(ibc,jbc+1,kbc,2)) / 2
                    u1 = (q1(ibc,jbc-1,kbc,2) + q1(ibc,jbc  ,kbc,2)) / 2
                    v2 = (q1(ibc,jbc  ,kbc,3) + q1(ibc,jbc+1,kbc,3)) / 2
                    v1 = (q1(ibc,jbc-1,kbc,3) + q1(ibc,jbc  ,kbc,3)) / 2
                    w2 = (q1(ibc,jbc  ,kbc,4) + q1(ibc,jbc+1,kbc,4)) / 2
                    w1 = (q1(ibc,jbc-1,kbc,4) + q1(ibc,jbc  ,kbc,4)) / 2
                  end if
                  dudeta = u2 - u1
                  dvdeta = v2 - v1
                  dwdeta = w2 - w1

                  if (i3.eq.istrt(ic3)) then
                    u2 = (  q1(ibc,jbc,kbc,2) + q1(ibc,jbc,kbc+1,2)) / 2
                    u1 = (3*q1(ibc,jbc,kbc,2) - q1(ibc,jbc,kbc+1,2)) / 2
                    v2 = (  q1(ibc,jbc,kbc,3) + q1(ibc,jbc,kbc+1,3)) / 2
                    v1 = (3*q1(ibc,jbc,kbc,3) - q1(ibc,jbc,kbc+1,3)) / 2
                    w2 = (  q1(ibc,jbc,kbc,4) + q1(ibc,jbc,kbc+1,4)) / 2
                    w1 = (3*q1(ibc,jbc,kbc,4) - q1(ibc,jbc,kbc+1,4)) / 2
                  else if (i3.eq.iend(ic3)) then
                    u2 = (3*q1(ibc,jbc,kbc,2) - q1(ibc,jbc,kbc-1,2)) / 2
                    u1 = (  q1(ibc,jbc,kbc,2) + q1(ibc,jbc,kbc-1,2)) / 2
                    v2 = (3*q1(ibc,jbc,kbc,3) - q1(ibc,jbc,kbc-1,3)) / 2
                    v1 = (  q1(ibc,jbc,kbc,3) + q1(ibc,jbc,kbc-1,3)) / 2
                    w2 = (3*q1(ibc,jbc,kbc,4) - q1(ibc,jbc,kbc-1,4)) / 2
                    w1 = (  q1(ibc,jbc,kbc,4) + q1(ibc,jbc,kbc-1,4)) / 2
                  else
                    u2 = (q1(ibc,jbc,kbc  ,2) + q1(ibc,jbc,kbc+1,2)) / 2
                    u1 = (q1(ibc,jbc,kbc-1,2) + q1(ibc,jbc,kbc  ,2)) / 2
                    v2 = (q1(ibc,jbc,kbc  ,3) + q1(ibc,jbc,kbc+1,3)) / 2
                    v1 = (q1(ibc,jbc,kbc-1,3) + q1(ibc,jbc,kbc  ,3)) / 2
                    w2 = (q1(ibc,jbc,kbc  ,4) + q1(ibc,jbc,kbc+1,4)) / 2
                    w1 = (q1(ibc,jbc,kbc-1,4) + q1(ibc,jbc,kbc  ,4)) / 2
                  end if
                  dudzeta = u2 - u1
                  dvdzeta = v2 - v1
                  dwdzeta = w2 - w1

                  xix   = xideri  (1,im,jm,km)
                  xiy   = xideri  (2,im,jm,km)
                  xiz   = xideri  (3,im,jm,km)
                  etax  = etaderi (1,im,jm,km)
                  etay  = etaderi (2,im,jm,km)
                  etaz  = etaderi (3,im,jm,km)
                  zetax = zetaderi(1,im,jm,km)
                  zetay = zetaderi(2,im,jm,km)
                  zetaz = zetaderi(3,im,jm,km)
   
                else if (idir.eq.2) then

                  if (i3.eq.istrt(ic3)) then
                    u2 = (  q1(ibc+1,jbc,kbc,2) + q1(ibc  ,jbc,kbc,2))/2
                    u1 = (3*q1(ibc  ,jbc,kbc,2) - q1(ibc+1,jbc,kbc,2))/2
                    v2 = (  q1(ibc+1,jbc,kbc,3) + q1(ibc  ,jbc,kbc,3))/2
                    v1 = (3*q1(ibc  ,jbc,kbc,3) - q1(ibc+1,jbc,kbc,3))/2
                    w2 = (  q1(ibc+1,jbc,kbc,4) + q1(ibc  ,jbc,kbc,4))/2
                    w1 = (3*q1(ibc  ,jbc,kbc,4) - q1(ibc+1,jbc,kbc,4))/2
                  else if (i3.eq.iend(ic3)) then
                    u2 = (3*q1(ibc,jbc,kbc,2) - q1(ibc-1,jbc,kbc,2)) / 2
                    u1 = (  q1(ibc,jbc,kbc,2) + q1(ibc-1,jbc,kbc,2)) / 2
                    v2 = (3*q1(ibc,jbc,kbc,3) - q1(ibc-1,jbc,kbc,3)) / 2
                    v1 = (  q1(ibc,jbc,kbc,3) + q1(ibc-1,jbc,kbc,3)) / 2
                    w2 = (3*q1(ibc,jbc,kbc,4) - q1(ibc-1,jbc,kbc,4)) / 2
                    w1 = (  q1(ibc,jbc,kbc,4) + q1(ibc-1,jbc,kbc,4)) / 2
                  else
                    u2 = (q1(ibc+1,jbc,kbc,2) + q1(ibc  ,jbc,kbc,2)) / 2
                    u1 = (q1(ibc  ,jbc,kbc,2) + q1(ibc-1,jbc,kbc,2)) / 2
                    v2 = (q1(ibc+1,jbc,kbc,3) + q1(ibc  ,jbc,kbc,3)) / 2
                    v1 = (q1(ibc  ,jbc,kbc,3) + q1(ibc-1,jbc,kbc,3)) / 2
                    w2 = (q1(ibc+1,jbc,kbc,4) + q1(ibc  ,jbc,kbc,4)) / 2
                    w1 = (q1(ibc  ,jbc,kbc,4) + q1(ibc-1,jbc,kbc,4)) / 2
                  end if
                  dudxi = u2 - u1
                  dvdxi = v2 - v1
                  dwdxi = w2 - w1

                  dudeta = (-1)**(1+1/inrout) * &
                           (-8*q1(ibc,jbc,kbc,2) + &
                             9*q1(ibc,jbc+joff1,kbc,2) - &
                               q1(ibc,jbc+joff2,kbc,2) ) / 3
                  dvdeta = (-1)**(1+1/inrout) * &
                           (-8*q1(ibc,jbc,kbc,3) + &
                             9*q1(ibc,jbc+joff1,kbc,3) - &
                               q1(ibc,jbc+joff2,kbc,3) ) / 3
                  dwdeta = (-1)**(1+1/inrout) * &
                           (-8*q1(ibc,jbc,kbc,4) + &
                             9*q1(ibc,jbc+joff1,kbc,4) - &
                               q1(ibc,jbc+joff2,kbc,4) ) / 3

                  if (i2.eq.istrt(ic2)) then
                    u2 = (  q1(ibc,jbc,kbc+1,2) + q1(ibc,jbc,kbc  ,2))/2
                    u1 = (3*q1(ibc,jbc,kbc  ,2) - q1(ibc,jbc,kbc+1,2))/2
                    v2 = (  q1(ibc,jbc,kbc+1,3) + q1(ibc,jbc,kbc  ,3))/2
                    v1 = (3*q1(ibc,jbc,kbc  ,3) - q1(ibc,jbc,kbc+1,3))/2
                    w2 = (  q1(ibc,jbc,kbc+1,4) + q1(ibc,jbc,kbc  ,4))/2
                    w1 = (3*q1(ibc,jbc,kbc  ,4) - q1(ibc,jbc,kbc+1,4))/2
                  else if (i2.eq.iend(ic2)) then
                    u2 = (3*q1(ibc,jbc,kbc,2) - q1(ibc,jbc,kbc-1,2)) / 2
                    u1 = (  q1(ibc,jbc,kbc,2) + q1(ibc,jbc,kbc-1,2)) / 2
                    v2 = (3*q1(ibc,jbc,kbc,3) - q1(ibc,jbc,kbc-1,3)) / 2
                    v1 = (  q1(ibc,jbc,kbc,3) + q1(ibc,jbc,kbc-1,3)) / 2
                    w2 = (3*q1(ibc,jbc,kbc,4) - q1(ibc,jbc,kbc-1,4)) / 2
                    w1 = (  q1(ibc,jbc,kbc,4) + q1(ibc,jbc,kbc-1,4)) / 2
                  else
                    u2 = (q1(ibc,jbc,kbc+1,2) + q1(ibc,jbc,kbc  ,2)) / 2
                    u1 = (q1(ibc,jbc,kbc  ,2) + q1(ibc,jbc,kbc-1,2)) / 2
                    v2 = (q1(ibc,jbc,kbc+1,3) + q1(ibc,jbc,kbc  ,3)) / 2
                    v1 = (q1(ibc,jbc,kbc  ,3) + q1(ibc,jbc,kbc-1,3)) / 2
                    w2 = (q1(ibc,jbc,kbc+1,4) + q1(ibc,jbc,kbc  ,4)) / 2
                    w1 = (q1(ibc,jbc,kbc  ,4) + q1(ibc,jbc,kbc-1,4)) / 2
                  end if
                  dudzeta = u2 - u1
                  dvdzeta = v2 - v1
                  dwdzeta = w2 - w1
   
                  xix   = xiderj  (1,im,jm,km)
                  xiy   = xiderj  (2,im,jm,km)
                  xiz   = xiderj  (3,im,jm,km)
                  etax  = etaderj (1,im,jm,km)
                  etay  = etaderj (2,im,jm,km)
                  etaz  = etaderj (3,im,jm,km)
                  zetax = zetaderj(1,im,jm,km)
                  zetay = zetaderj(2,im,jm,km)
                  zetaz = zetaderj(3,im,jm,km)

                else if (idir.eq.3) then
                  if (i2.eq.istrt(ic2)) then
                    u2 = (  q1(ibc,jbc,kbc,2) + q1(ibc+1,jbc,kbc,2)) / 2
                    u1 = (3*q1(ibc,jbc,kbc,2) - q1(ibc+1,jbc,kbc,2)) / 2
                    v2 = (  q1(ibc,jbc,kbc,3) + q1(ibc+1,jbc,kbc,3)) / 2
                    v1 = (3*q1(ibc,jbc,kbc,3) - q1(ibc+1,jbc,kbc,3)) / 2
                    w2 = (  q1(ibc,jbc,kbc,4) + q1(ibc+1,jbc,kbc,4)) / 2
                    w1 = (3*q1(ibc,jbc,kbc,4) - q1(ibc+1,jbc,kbc,4)) / 2
                  else if (i2.eq.iend(ic2)) then
                    u2 = (3*q1(ibc,jbc,kbc,2) - q1(ibc-1,jbc,kbc,2)) / 2
                    u1 = (  q1(ibc,jbc,kbc,2) + q1(ibc-1,jbc,kbc,2)) / 2
                    v2 = (3*q1(ibc,jbc,kbc,3) - q1(ibc-1,jbc,kbc,3)) / 2
                    v1 = (  q1(ibc,jbc,kbc,3) + q1(ibc-1,jbc,kbc,3)) / 2
                    w2 = (3*q1(ibc,jbc,kbc,4) - q1(ibc-1,jbc,kbc,4)) / 2
                    w1 = (  q1(ibc,jbc,kbc,4) + q1(ibc-1,jbc,kbc,4)) / 2
                  else
                    u2 = (q1(ibc  ,jbc,kbc,2) + q1(ibc+1,jbc,kbc,2)) / 2
                    u1 = (q1(ibc-1,jbc,kbc,2) + q1(ibc  ,jbc,kbc,2)) / 2
                    v2 = (q1(ibc  ,jbc,kbc,3) + q1(ibc+1,jbc,kbc,3)) / 2
                    v1 = (q1(ibc-1,jbc,kbc,3) + q1(ibc  ,jbc,kbc,3)) / 2
                    w2 = (q1(ibc  ,jbc,kbc,4) + q1(ibc+1,jbc,kbc,4)) / 2
                    w1 = (q1(ibc-1,jbc,kbc,4) + q1(ibc  ,jbc,kbc,4)) / 2
                  end if
                  dudxi = u2 - u1
                  dvdxi = v2 - v1
                  dwdxi = w2 - w1

                  if (i3.eq.istrt(ic3)) then
                    u2 = (  q1(ibc,jbc,kbc,2) + q1(ibc,jbc+1,kbc,2)) / 2
                    u1 = (3*q1(ibc,jbc,kbc,2) - q1(ibc,jbc+1,kbc,2)) / 2
                    v2 = (  q1(ibc,jbc,kbc,3) + q1(ibc,jbc+1,kbc,3)) / 2
                    v1 = (3*q1(ibc,jbc,kbc,3) - q1(ibc,jbc+1,kbc,3)) / 2
                    w2 = (  q1(ibc,jbc,kbc,4) + q1(ibc,jbc+1,kbc,4)) / 2
                    w1 = (3*q1(ibc,jbc,kbc,4) - q1(ibc,jbc+1,kbc,4)) / 2
                  else if (i3.eq.iend(ic3)) then
                    u2 = (3*q1(ibc,jbc,kbc,2) - q1(ibc,jbc-1,kbc,2)) / 2
                    u1 = (  q1(ibc,jbc,kbc,2) + q1(ibc,jbc-1,kbc,2)) / 2
                    v2 = (3*q1(ibc,jbc,kbc,3) - q1(ibc,jbc-1,kbc,3)) / 2
                    v1 = (  q1(ibc,jbc,kbc,3) + q1(ibc,jbc-1,kbc,3)) / 2
                    w2 = (3*q1(ibc,jbc,kbc,4) - q1(ibc,jbc-1,kbc,4)) / 2
                    w1 = (  q1(ibc,jbc,kbc,4) + q1(ibc,jbc-1,kbc,4)) / 2
                  else
                    u2 = (q1(ibc,jbc  ,kbc,2) + q1(ibc,jbc+1,kbc,2)) / 2
                    u1 = (q1(ibc,jbc-1,kbc,2) + q1(ibc,jbc  ,kbc,2)) / 2
                    v2 = (q1(ibc,jbc  ,kbc,3) + q1(ibc,jbc+1,kbc,3)) / 2
                    v1 = (q1(ibc,jbc-1,kbc,3) + q1(ibc,jbc  ,kbc,3)) / 2
                    w2 = (q1(ibc,jbc  ,kbc,4) + q1(ibc,jbc+1,kbc,4)) / 2
                    w1 = (q1(ibc,jbc-1,kbc,4) + q1(ibc,jbc  ,kbc,4)) / 2
                  end if
                  dudeta = u2 - u1
                  dvdeta = v2 - v1
                  dwdeta = w2 - w1

                  dudzeta = (-1)**(1+1/inrout) * &
                            (-8*q1(ibc,jbc,kbc,2) + &
                              9*q1(ibc,jbc,kbc+koff1,2) - &
                                q1(ibc,jbc,kbc+koff2,2) ) / 3
                  dvdzeta = (-1)**(1+1/inrout) * &
                            (-8*q1(ibc,jbc,kbc,3) + &
                              9*q1(ibc,jbc,kbc+koff1,3) - &
                                q1(ibc,jbc,kbc+koff2,3) ) / 3
                  dwdzeta = (-1)**(1+1/inrout) * &
                            (-8*q1(ibc,jbc,kbc,4) + &
                              9*q1(ibc,jbc,kbc+koff1,4) - &
                                q1(ibc,jbc,kbc+koff2,4) ) / 3

                  xix   = xiderk  (1,im,jm,km)
                  xiy   = xiderk  (2,im,jm,km)
                  xiz   = xiderk  (3,im,jm,km)
                  etax  = etaderk (1,im,jm,km)
                  etay  = etaderk (2,im,jm,km)
                  etaz  = etaderk (3,im,jm,km)
                  zetax = zetaderk(1,im,jm,km)
                  zetay = zetaderk(2,im,jm,km)
                  zetaz = zetaderk(3,im,jm,km)
   
                end if

                dudx = dudxi*xix + dudeta*etax + dudzeta*zetax
                dudy = dudxi*xiy + dudeta*etay + dudzeta*zetay
                dudz = dudxi*xiz + dudeta*etaz + dudzeta*zetaz
                dvdx = dvdxi*xix + dvdeta*etax + dvdzeta*zetax
                dvdy = dvdxi*xiy + dvdeta*etay + dvdzeta*zetay
                dvdz = dvdxi*xiz + dvdeta*etaz + dvdzeta*zetaz
                dwdx = dwdxi*xix + dwdeta*etax + dwdzeta*zetax
                dwdy = dwdxi*xiy + dwdeta*etay + dwdzeta*zetay
                dwdz = dwdxi*xiz + dwdeta*etaz + dwdzeta*zetaz

                divv = dudx + dvdy + dwdz

                txx = 2*muw * (dudx - divv/3)
                txy =   muw * (dudy + dvdx)
                txz =   muw * (dudz + dwdx)
                tyy = 2*muw * (dvdy - divv/3)
                tyz =   muw * (dvdz + dwdy)
                tzz = 2*muw * (dwdz - divv/3)

                fv(1) = machfs / reyno * ds * &
                        (txx*kx + txy*ky + txz*kz)
                fv(2) = machfs / reyno * ds * &
                        (txy*kx + tyy*ky + tyz*kz)
                fv(3) = machfs / reyno * ds * &
                        (txz*kx + tyz*ky + tzz*kz)
                xarm = x(iw,jw,kw) - xmomr
                yarm = y(iw,jw,kw) - ymomr
                zarm = z(iw,jw,kw) - zmomr
                fv(4) = machfs / reyno * ds * &
                  (-(txy*kx + tyy*ky + tyz*kz) * zarm + &
                    (txz*kx + tyz*ky + tzz*kz) * yarm )
                fv(5) = machfs / reyno * ds * &
                  (-(txz*kx + tyz*ky + tzz*kz) * xarm + &
                    (txx*kx + txy*ky + txz*kz) * zarm )
                fv(6) = machfs / reyno * ds * &
                  (-(txx*kx + txy*ky + txz*kz) * yarm + &
                    (txy*kx + tyy*ky + tyz*kz) * xarm )

                sur_dim = sur_dim + 1

                cl_v(sur_dim,nsurf1) = (fv(1)*aax+fv(2)*aay+fv(3)*aaz) * &
                           2/(lngth(1)*lngth(2)*mach**2)
                cd_v(sur_dim,nsurf1) = (fv(1)*ssx+fv(2)*ssy+fv(3)*ssz) * &
                           2/(lngth(1)*lngth(2)*mach**2)
                cz_v(sur_dim,nsurf1) = (fv(1)*bbx+fv(2)*bby+fv(3)*bbz) * &
                           2/(lngth(1)*lngth(2)*mach**2)
                cmx_v(sur_dim,nsurf1) = fv(4) * &
                           2/(lngth(1)*lngth(2)*lngth(3)*mach**2)
                cmy_v(sur_dim,nsurf1) = fv(5) * &
                           2/(lngth(1)*lngth(2)*lngth(3)*mach**2)
                cm_v (sur_dim,nsurf1) = fv(6) * &
                           2/(lngth(1)*lngth(2)*lngth(3)*mach**2)

              end do
            end do

          end if

          if (sur_dim.ne.surf_dim(nsurf1)) then
            write(*,*) 'Incongruency in forces_b. Aborting!'
            stop
          end if

        end if

      end do
111   format(4e14.7)

      return
      end

!-----------------------------------------------------------------------
      subroutine wr_forces(ijkmaxg,xsurf,ysurf,zsurf,nsurf,msurf, &
                           surf_dim,ijkmax_surf,cl_p,cl_v,cd_p,cd_v, &
                           cz_p,cz_v,cmx_p,cmx_v,cmy_p,cmy_v,cm_p,cm_v)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) ijkmaxg,nsurf,msurf,surf_dim(msurf)
      integer(kind=cosa_int) i,itime1,isurf,ijkmax_surf(msurf,3)
      real(kind=cosa_real) &
           xsurf(ijkmaxg-1,msurf), ysurf(ijkmaxg-1,msurf), &
           zsurf(ijkmaxg-1,msurf), &
           cl_p (ijkmaxg-1,msurf), cd_p (ijkmaxg-1,msurf), &
           cz_p (ijkmaxg-1,msurf), cmx_p(ijkmaxg-1,msurf), &
           cmy_p(ijkmaxg-1,msurf), cm_p (ijkmaxg-1,msurf), &
           cl_v (ijkmaxg-1,msurf), cd_v (ijkmaxg-1,msurf), &
           cz_v (ijkmaxg-1,msurf), cmx_v(ijkmaxg-1,msurf), &
           cmy_v(ijkmaxg-1,msurf), cm_v (ijkmaxg-1,msurf)

      if (dualt) then
        if (itime.le.9) then
          write(flowtec,'(''000'',i1)') itime
        elseif (itime.le.99) then
          write(flowtec,'(''00'',i2)') itime
        elseif (itime.le.999) then
          write(flowtec,'(''0'',i3)') itime
        else
          write(flowtec,'(i4)') itime
        end if
      else if (harbal) then
        if (itime.le.9) then
          write(flowtec,'(''00'',i1)') itime
        elseif (itime.le.99) then
          write(flowtec,'(''0'',i2)') itime
        else
          write(flowtec,'(i3)') itime
        end if
      else
        flowtec = 'steady'
      end if

      flowtec = 'flow_tec_forces_'//flowtec
      flowtec = trim(flowtec)//'.dat'

!---- write header of tecplot file
      open(200,file=flowtec,status='replace')
      write(200,*) 'TITLE="Surface forces and moment coefficients"'
      write(200,*) 'VARIABLES=x,y,z,cl_p,cl_v,cd_p,cd_v,cz_p,cz_v,cmx_p, &
     &cmx_v,cmy_p,cmy_v,cm_p,cm_v'

      do isurf = 1,nsurf

        write(200,11) ijkmax_surf(isurf,1), ijkmax_surf(isurf,2), &
                      ijkmax_surf(isurf,3)

        do i=1,surf_dim(isurf)
          write(200,20) xsurf(i,isurf),ysurf(i,isurf),zsurf(i,isurf), &
                        cl_p(i,isurf),cl_v(i,isurf),cd_p(i,isurf), &
                        cd_v(i,isurf),cz_p(i,isurf),cz_v(i,isurf), &
                        cmx_p(i,isurf),cmx_v(i,isurf),cmy_p(i,isurf), &
                        cmy_v(i,isurf),cm_p(i,isurf),cm_v(i,isurf)
        end do

      end do

      close(unit=200)    

 20   format(300e14.6)
 11   format('ZONE T="arturo",I=',i4,',J=',i4,',K=',i4,', F=POINT, DT=(S &
      INGLE SINGLE SINGLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUB &
     &LE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)')

!-----The lines below are for writing out the ascii file for each patch 
!     separately.

!     do isurf = 1,nsurf
!
!       if (unsteady) then
!
!         if (dualt) then
!           if (itime.le.9) then
!             write(flowtec,'(''000'',i1)') itime
!           elseif (itime.le.99) then
!             write(flowtec,'(''00'',i2)') itime
!           elseif (itime.le.999) then
!             write(flowtec,'(''0'',i3)') itime
!           else
!             write(flowtec,'(i4)') itime
!           end if
!         else if (harbal) then
!           if (itime.le.9) then
!             write(flowtec,'(''00'',i1)') itime
!           elseif (itime.le.99) then
!             write(flowtec,'(''0'',i2)') itime
!           else
!             write(flowtec,'(i3)') itime
!           end if
!         end if
!
!         if (isurf.le.9) then
!           write(filename,'(''000'',i1)') isurf
!         else if (isurf.le.99) then
!           write(filename,'(''00'',i2)') isurf
!         else if (isurf.le.999) then
!           write(filename,'(''0'',i3)') isurf
!         end if
!         flowtec = 'unst_forces_TS_'//flowtec
!         filename = 'unst_cp_s'//filename
!         filename = trim(flowtec)//'_s_'//filename
!         filename = trim(filename)//'.dat'
!         open(7,file=filename,status='replace')
!         write(7,10) (simutime(itime1),itime1=1,ntime)
!         write(7,20) simtime
!         do i=1,surf_dim(isurf)
!           write(7,20) xsurf(i,isurf),ysurf(i,isurf),zsurf(i,isurf), &
!                       cl_p(i,isurf),cl_v(i,isurf),cd_p(i,isurf), &
!                       cd_v(i,isurf),cz_p(i,isurf),cz_v(i,isurf), &
!                       cmx_p(i,isurf),cmx_v(i,isurf),cmy_p(i,isurf), &
!                       cmy_v(i,isurf),cm_p(i,isurf),cm_v(i,isurf)
!         end do
!         close(unit=7)    
!
!       else
!
!         if (isurf.le.9) then
!           write(filename,'(''000'',i1)') isurf
!         else if (isurf.le.99) then
!           write(filename,'(''00'',i2)') isurf
!         else if (isurf.le.999) then
!           write(filename,'(''0'',i3)') isurf
!         end if
!         filename = 'stea_forces_s'//filename
!         filename = trim(filename)//'.dat'
!         open(7,file=filename,status='replace')
!         do i=1,surf_dim(isurf)
!           write(7,20) xsurf(i,isurf),ysurf(i,isurf),zsurf(i,isurf), &
!                       cl_p(i,isurf),cl_v(i,isurf),cd_p(i,isurf), &
!                       cd_v(i,isurf),cz_p(i,isurf),cz_v(i,isurf), &
!                       cmx_p(i,isurf),cmx_v(i,isurf),cmy_p(i,isurf), &
!                       cmy_v(i,isurf),cm_p(i,isurf),cm_v(i,isurf)
!         end do
!         close(unit=7)    
!
!       end if
!
!     end do

      return
      end

!-----------------------------------------------------------------------
      subroutine vert2cell_g(nl,x,y,z,workx,worky,workz)
!-----------------------------------------------------------------------
       
      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,iblock,imax,jmax,kmax,ixyz,dim5l,dim5u
      real (kind=cosa_real) x(*),y(*),z(*),workx(*),worky(*),workz(*)

      dim5l = 0
      dim5u = 2 * nharms * hbmove
      call copy_array(1,1,dim5l,dim5u, 2,nl,workx,x, &
                      1,1,dim5l,dim5u, 0)
      call copy_array(1,1,dim5l,dim5u, 2,nl,worky,y, &
                      1,1,dim5l,dim5u, 0)
      call copy_array(1,1,dim5l,dim5u, 2,nl,workz,z, &
                      1,1,dim5l,dim5u, 0)

      do iblock = 1,mynblocks
        imax   = i_imax     (iblock,nl)
        jmax   = j_jmax     (iblock,nl)
        kmax   = k_kmax     (iblock,nl)
        ixyz   = 1 + off_p2 (iblock,nl)        * dim5h
        call vert2cell_bg(x(ixyz),y(ixyz),z(ixyz),workx(ixyz), &
          worky(ixyz),workz(ixyz),imax,jmax,kmax,npde,nharms)
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine vert2cell_bg(x,y,z,workx,worky,workz,imax,jmax,kmax, &
        npde,nharms)
!-----------------------------------------------------------------------
       
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,npde,nharms,nl
      integer(kind=cosa_int) i,im,i1,j,jm,j1,k,km,k1,n
      real (kind=cosa_real) &
           x     (0:imax+1,0:jmax+1,0:kmax+1,0:2*nharms*hbmove), &
           y     (0:imax+1,0:jmax+1,0:kmax+1,0:2*nharms*hbmove), &
           z     (0:imax+1,0:jmax+1,0:kmax+1,0:2*nharms*hbmove), &
           workx (0:imax+1,0:jmax+1,0:kmax+1,0:2*nharms*hbmove), &
           worky (0:imax+1,0:jmax+1,0:kmax+1,0:2*nharms*hbmove), &
           workz (0:imax+1,0:jmax+1,0:kmax+1,0:2*nharms*hbmove)

      im = imax-1
      jm = jmax-1
      km = kmax-1

      do n = 0,2*nharms*hbmove

!------ block vertices

        do k=0,kmax,kmax
          if (k.eq.0) then
            k1 = 1
          else if (k.eq.kmax) then
            k1 = kmax
          end if
          x(0   ,0   ,k,n) = workx(1   ,1   ,k1   ,n)
          y(0   ,0   ,k,n) = worky(1   ,1   ,k1   ,n)
          z(0   ,0   ,k,n) = workz(1   ,1   ,k1   ,n)
          x(imax,0   ,k,n) = workx(imax,1   ,k1   ,n)    
          y(imax,0   ,k,n) = worky(imax,1   ,k1   ,n)
          z(imax,0   ,k,n) = workz(imax,1   ,k1   ,n)
          x(0   ,jmax,k,n) = workx(1   ,jmax,k1   ,n)
          y(0   ,jmax,k,n) = worky(1   ,jmax,k1   ,n)
          z(0   ,jmax,k,n) = workz(1   ,jmax,k1   ,n)
          x(imax,jmax,k,n) = workx(imax,jmax,k1   ,n)
          y(imax,jmax,k,n) = worky(imax,jmax,k1   ,n)
          z(imax,jmax,k,n) = workz(imax,jmax,k1   ,n)
        end do

!------ block edges

        do k=0,kmax,kmax
          if (k.eq.0) then
            k1 = 1
          else if (k.eq.kmax) then
            k1 = kmax
          end if
          do j=0,jmax,jmax
            if (j.eq.0) then
              j1 = 1
            else if (j.eq.jmax) then
              j1 = jmax
            end if
            do i = 1,im
              x(i,j,k,n) = (workx(i,j1,k1,n) + workx(i+1,j1,k1,n)) / 2
              y(i,j,k,n) = (worky(i,j1,k1,n) + worky(i+1,j1,k1,n)) / 2
              z(i,j,k,n) = (workz(i,j1,k1,n) + workz(i+1,j1,k1,n)) / 2
            end do
          end do
        end do

        do k=0,kmax,kmax
          if (k.eq.0) then
            k1 = 1
          else if (k.eq.kmax) then
            k1 = kmax
          end if
          do i=0,imax,imax
            if (i.eq.0) then
              i1 = 1
            else if (i.eq.imax) then
              i1 = imax
            end if
            do j = 1,jm
              x(i,j,k,n) = (workx(i1,j,k1,n) + workx(i1,j+1,k1,n)) / 2
              y(i,j,k,n) = (worky(i1,j,k1,n) + worky(i1,j+1,k1,n)) / 2
              z(i,j,k,n) = (workz(i1,j,k1,n) + workz(i1,j+1,k1,n)) / 2
            end do
          end do
        end do

        do j=0,jmax,jmax
          if (j.eq.0) then
            j1 = 1
          else if (j.eq.jmax) then
            j1 = jmax
          end if
          do i=0,imax,imax
            if (i.eq.0) then
              i1 = 1
            else if (i.eq.imax) then
              i1 = imax
            end if
            do k = 1,km
              x(i,j,k,n) = (workx(i1,j1,k,n) + workx(i1,j1,k+1,n)) / 2
              y(i,j,k,n) = (worky(i1,j1,k,n) + worky(i1,j1,k+1,n)) / 2
              z(i,j,k,n) = (workz(i1,j1,k,n) + workz(i1,j1,k+1,n)) / 2
            end do
          end do
        end do

!------ block bounding surfaces

        do k=0,kmax,kmax
          if (k.eq.0) then
            k1 = 1
          else if (k.eq.kmax) then
            k1 = kmax
          end if
          do j = 1,jm
            do i = 1,im
              x(i,j,k,n) = &
                (workx(i  ,j  ,k1,n) + workx(i+1,j  ,k1,n) + &
                 workx(i  ,j+1,k1,n) + workx(i+1,j+1,k1,n)) / 4
              y(i,j,k,n) = &
                (worky(i  ,j  ,k1,n) + worky(i+1,j  ,k1,n) + &
                 worky(i  ,j+1,k1,n) + worky(i+1,j+1,k1,n)) / 4
              z(i,j,k,n) = &
                (workz(i  ,j  ,k1,n) + workz(i+1,j  ,k1,n) + &
                 workz(i  ,j+1,k1,n) + workz(i+1,j+1,k1,n)) / 4
            end do
          end do
        end do

        do j=0,jmax,jmax
          if (j.eq.0) then
            j1 = 1
          else if (j.eq.jmax) then
            j1 = jmax
          end if
          do k = 1,km
            do i = 1,im
              x(i,j,k,n) = &
                (workx(i  ,j1,k  ,n) + workx(i+1,j1,k  ,n) + &
                 workx(i  ,j1,k+1,n) + workx(i+1,j1,k+1,n)) / 4
              y(i,j,k,n) = &
                (worky(i  ,j1,k  ,n) + worky(i+1,j1,k  ,n) + &
                 worky(i  ,j1,k+1,n) + worky(i+1,j1,k+1,n)) / 4
              z(i,j,k,n) = &
                (workz(i  ,j1,k  ,n) + workz(i+1,j1,k  ,n) + &
                 workz(i  ,j1,k+1,n) + workz(i+1,j1,k+1,n)) / 4
            end do
          end do
        end do

        do i=0,imax,imax
          if (i.eq.0) then
            i1 = 1
          else if (i.eq.imax) then
            i1 = imax
          end if
          do k = 1,km
            do j = 1,jm
              x(i,j,k,n) = &
                (workx(i1,j  ,k  ,n) + workx(i1,j  ,k+1,n) + &
                 workx(i1,j+1,k  ,n) + workx(i1,j+1,k+1,n)) / 4
              y(i,j,k,n) = &
                (worky(i1,j  ,k  ,n) + worky(i1,j  ,k+1,n) + &
                 worky(i1,j+1,k  ,n) + worky(i1,j+1,k+1,n)) / 4
              z(i,j,k,n) = &
                (workz(i1,j  ,k  ,n) + workz(i1,j  ,k+1,n) + &
                 workz(i1,j+1,k  ,n) + workz(i1,j+1,k+1,n)) / 4
            end do
          end do
        end do

!------ internal nodes

        do k = 1,km
          do j = 1,jm
            do i = 1,im
              x(i,j,k,n) = &
                (workx(i  ,j  ,k  ,n) + workx(i+1,j  ,k  ,n) + &
                 workx(i  ,j+1,k  ,n) + workx(i+1,j+1,k  ,n) + &
                 workx(i  ,j  ,k+1,n) + workx(i+1,j  ,k+1,n) + &
                 workx(i  ,j+1,k+1,n) + workx(i+1,j+1,k+1,n) ) / 8
              y(i,j,k,n) = &
                (worky(i  ,j  ,k  ,n) + worky(i+1,j  ,k  ,n) + &
                 worky(i  ,j+1,k  ,n) + worky(i+1,j+1,k  ,n) + &
                 worky(i  ,j  ,k+1,n) + worky(i+1,j  ,k+1,n) + &
                 worky(i  ,j+1,k+1,n) + worky(i+1,j+1,k+1,n) ) / 8
              z(i,j,k,n) = &
                (workz(i  ,j  ,k  ,n) + workz(i+1,j  ,k  ,n) + &
                 workz(i  ,j+1,k  ,n) + workz(i+1,j+1,k  ,n) + &
                 workz(i  ,j  ,k+1,n) + workz(i+1,j  ,k+1,n) + &
                 workz(i  ,j+1,k+1,n) + workz(i+1,j+1,k+1,n) ) / 8
            end do
          end do
        end do

      end do

      return
      end

      end module postprocessing
 


!-----------------------------------------------------------------------
      subroutine writerest(nl,q,q1,mut)
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,iblock,imax,jmax,kmax,iq,imut,fid
      real (kind=cosa_real) q(*),q1(*),mut(*)
      logical parallel

      call usingparallel(parallel)

      if(parallel) then

        call parallelwriterestart(q,q1,mut,nl)

      else

        fid = 21
        open(fid,file='restart',status='unknown')
        close(fid,status='delete')
        open(fid,file='restart',form='unformatted',status='new')

        do iblock = 1,mynblocks
          imax = i_imax  (iblock,nl)
          jmax = j_jmax  (iblock,nl)
          kmax = k_kmax  (iblock,nl)
          iq   = 1 + off_p3  (iblock,nl) * npde * dim5
          if (kom.or.kom_bsl.or.kom_sst) then 
             imut   = 1 + off_p3 (iblock,nl) *        dim5
          else 
             imut = 1
          end if
          call write_brest(q(iq),q1(imut),mut(iq),imax,jmax,kmax,npde,nharms,fid)
        end do

        close(fid)

      end if

      return
      end

!-----------------------------------------------------------------------
      subroutine write_brest(q,q1,mut,imax,jmax,kmax,npde,nharms,fid)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,npde,nharms
      integer(kind=cosa_int) i,imax1,j,jmax1,k,kmax1,ipde,key1,n,fid
      real (kind=cosa_real) &
           q (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           q1(-1:imax+1,-1:jmax+1,-1:kmax+1,npde), &
           mut(-1:imax+1,-1:jmax+1,-1:kmax+1 ,0:2*nharms)    

      imax1 = imax+1
      jmax1 = jmax+1
      kmax1 = kmax+1

      write(fid) (((((q(i,j,k,ipde,n), &
        i=-1,imax1),j=-1,jmax1),k=-1,kmax1),ipde=1,npde),n=0,2*nharms)

      if (dualt) then
        write(fid) ((((q1(i,j,k,ipde), &
          i=-1,imax1),j=-1,jmax1),k=-1,kmax1),ipde=1,npde)
      end if

      if (kom.or.kom_bsl.or.kom_sst) then
        write(fid) ((((mut(i,j,k,n), &
          i=-1,imax1),j=-1,jmax1),k=-1,kmax1),n=0,2*nharms)
      end if

      return
      end

!-----------------------------------------------------------------------
      subroutine out_tec(nl)
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables       
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl
      integer(kind=cosa_int) bctopo,cutopo,percutopo
      real (kind=cosa_real) q,qtec,mut,x,y,z,xdot,ydot,zdot,dist,si,sj,sk,xideri, &
        xiderj,xiderk,etaderi,etaderj,etaderk,zetaderi,zetaderj, &
        zetaderk,var1,var2,work
      pointer &
        (pq   ,   q), (pqtec,qtec), (pmut , mut), (pdist,dist), &
        (pvar1,var1), (pvar2,var2), (pwork,work), (px   ,   x), &
        (py   ,   y), (pz   ,   z), (pxdot,xdot), (pydot,ydot), &
        (pzdot,zdot), (psi , si), (psj , sj), (psk , sk), &
        (pxideri  ,  xideri), (pxiderj  ,  xiderj), (pxiderk  ,  xiderk), &
        (petaderi , etaderi), (petaderj , etaderj), (petaderk , etaderk), &
        (pzetaderi,zetaderi), (pzetaderj,zetaderj), (pzetaderk,zetaderk), &
        (pbctopo,bctopo), (pcutopo,cutopo), (ppercutopo,percutopo)

      pq        = p_q(nl)
      pqtec     = p_qold(nl)
      px        = p_x(nl)
      py        = p_y(nl)
      pz        = p_z(nl)
      pbctopo   = p_bctopo(nl)
      pcutopo   = p_cutopo(nl)
      ppercutopo= p_percutopo(nl)
      pvar1     = p_qp(nl)
      pvar2     = p_qm(nl)
      pwork     = p_dq(nl)
      psi       = p_si(nl)
      psj       = p_sj(nl)
      psk       = p_sk(nl)
      pxideri   = p_xideri(nl)
      pxiderj   = p_xiderj(nl)
      pxiderk   = p_xiderk(nl)
      petaderi  = p_etaderi(nl)
      petaderj  = p_etaderj(nl)
      petaderk  = p_etaderk(nl)
      pzetaderi = p_zetaderi(nl)
      pzetaderj = p_zetaderj(nl)
      pzetaderk = p_zetaderk(nl)
      if (kom.or.kom_bsl.or.kom_sst) then
        pmut     = p_mut(nl)
        pdist    = p_dist(nl)
        calmut   = .true.
      else
        calmut   = .false.
      end if
      if (moving) then
        pxdot    = p_xdot(nl)
        pydot    = p_ydot(nl)
        pzdot    = p_zdot(nl)
      end if

      call copy_array(1,npde,0,2*nharms, 3,nl,qtec,q, &
                      1,npde,0,2*nharms, 1)
      call cutman_qtec(nl,cutopo,qtec)
      call pcutman_qtec(nl,percutopo,qtec)
      call bc_tec(nl,bctopo,qtec,mut,xdot,ydot,zdot)
      call q_edges(nl,qtec,mut)
      call q_corns(nl,qtec,mut)
      call cons2prim(nl,1,qtec,work)
      if(write_tec) then
         call q_tec(nl,qtec,mut,x,y,z,xdot,ydot,zdot,dist,var1,var2)
         call wr_tec(nl,var1,var2)
         call wr_tec_surf(nl,qtec,x,y,z,xdot,ydot,zdot,dist,si,sj,sk,xideri, &
              xiderj,xiderk,etaderi,etaderj,etaderk,zetaderi,zetaderj,zetaderk, &
              bctopo)
      else if(write_plt) then
         write(*,*) 'Writing flowtec files as plt format'
         write(*,*) 'is not yet implemented'
      else if(write_szplt) then
         write(*,*) 'Writing flowtec files as plt format'
         write(*,*) 'is not yet implemented'
      else if(write_cgns) then
         call wr_cgns(nl,qtec,mut,x,y,z,xdot,ydot,zdot,dist)
         call wr_cgns_surf(nl,qtec,x,y,z,xdot,ydot,zdot,dist,si,sj,sk,xideri, &
              xiderj,xiderk,etaderi,etaderj,etaderk,zetaderi,zetaderj,zetaderk, &
               bctopo)
      end if



      return
      end

!-----------------------------------------------------------------------
      subroutine bc_tec(nl,bctopo,q,mut,xdot,ydot,zdot)
!-----------------------------------------------------------------------
!     completion of tecplot flow state on symmetry boundaries.
!-----------------------------------------------------------------------
       
      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax
      integer(kind=cosa_int) nl,iblock,iblk,iq,imut,ixyz
      integer(kind=cosa_int) bctopo(*)
      real (kind=cosa_real) q(*), mut(*), xdot(*), ydot(*), zdot(*)

      do iblock = 1,mynblocks
        imax   = i_imax     (iblock,nl)
        jmax   = j_jmax     (iblock,nl)
        kmax   = k_kmax     (iblock,nl)
        iblk   = 1 + off_bct(iblock,nl)
        iq     = 1 + off_p3 (iblock,nl) * npde * dim5
        imut   = 1 + off_p3 (iblock,nl) *        dim5
        ixyz   = 1 + off_p2 (iblock,nl)        * dim5h
        call bbc_tec(bctopo(iblk),nbcs(iblock),q(iq),mut(imut), &
          xdot(ixyz),ydot(ixyz),zdot(ixyz),imax,jmax,kmax,npde,nharms)
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine bbc_tec(bctopo,nbcs,q,mut,xdot,ydot,zdot,imax,jmax, &
                         kmax,npde,nharms)
!-----------------------------------------------------------------------
!              
!-----------------------------------------------------------------------
       
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,npde,nbcs,nharms
      integer(kind=cosa_int) ibc
      integer(kind=cosa_int) bctopo(10,nbcs)
      real (kind=cosa_real) &
           q   (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           mut (-1:imax+1,-1:jmax+1,-1:kmax+1     ,0:2*nharms), &
           xdot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms), &
           ydot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms), &
           zdot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms)

      do ibc = 1,nbcs

!------ symmetry (bctyp=81 or bctyp=82)
        if ((bctopo(1,ibc).eq.81).or.(bctopo(1,ibc).eq.82).or. &
            (bctopo(1,ibc).eq.83)) then

          call bc_sym_tec(q,mut,xdot,ydot,zdot,bctopo(1,ibc),imax,jmax, &
                          kmax,npde,nharms)

        end if

      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine bc_sym_tec(q,mut,xdot,ydot,zdot,bctopo,imax,jmax,kmax, &
                            npde,nharms)
!-----------------------------------------------------------------------
!     symmetry wrt xy-plane                                           81
!     symmetry wrt xz-plane                                           82
!     symmetry wrt yz-plane                                           83
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,npde,nharms
      integer(kind=cosa_int) bctopo(10)
      integer(kind=cosa_int) ipde,n,i1,i2,i3,ijkmax(3),idir,istrt(3),iend(3),bctyp, &
                inrout,ibcpt,ibcpt2,ibcn,ibcn2,ibc,jbc,kbc,ibc2,jbc2, &
                kbc2,in,jn,kn,in2,jn2,kn2,ic1,ic2,ic3,isym(7),isymdr
      real (kind=cosa_real) &
           q   (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           mut (-1:imax+1,-1:jmax+1,-1:kmax+1     ,0:2*nharms), &
           xdot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms), &
           ydot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms), &
           zdot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms)

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
        ibcpt  =  0
        ibcpt2 = -1
        ibcn   =  1
        ibcn2  =  2
      else
        ibcpt  = ijkmax(idir)
        ibcpt2 = ijkmax(idir) + 1
        ibcn   = ijkmax(idir) - 1
        ibcn2  = ijkmax(idir) - 2
      end if

      ic1 = cyc (idir, 1)
      ic2 = cyc (idir, 2)
      ic3 = cyc (idir, 3)

!     Set up the symmetry array such that all variables are
!     reflected across the symmetry plane.
!
      if (bctyp.eq.81) then
         isymdr = 3
      else if (bctyp.eq.82) then
         isymdr = 2
      else if (bctyp.eq.83) then
         isymdr = 1
      end if
      do ipde = 1,npde
         isym(ipde) = 1
      end do

!del  isym(isymdr+1) = -1
      isym(isymdr+1) =  0

      do n = 0,2*nharms

        do ipde = 1,npde

          do i3 = istrt(ic3),iend(ic3)
            do i2 = istrt(ic2),iend(ic2)

              ibc  = ibcpt *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
              jbc  = ibcpt *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
              kbc  = ibcpt *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

              ibc2 = ibcpt2*krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
              jbc2 = ibcpt2*krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
              kbc2 = ibcpt2*krd(Ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

              in   = ibcn  *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
              jn   = ibcn  *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
              kn   = ibcn  *krd(ic1,3) + I2*krd(ic2,3) + i3*krd(ic3,3)

              in2  = ibcn2 *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
              jn2  = ibcn2 *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
              kn2  = ibcn2 *krd(ic1,3) + I2*krd(ic2,3) + i3*krd(ic3,3)

              q(ibc,jbc,kbc,ipde,n) = isym(ipde) * &
                (9*q(in ,jn ,kn ,ipde,n) - q(in2,jn2,kn2,ipde,n)) / 8

            end do
          end do

        end do

        if (kom.or.kom_bsl.or.kom_sst) then

          do i3 = istrt(ic3),iend(ic3)
            do i2 = istrt(ic2),iend(ic2)

              ibc  = ibcpt *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
              jbc  = ibcpt *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
              kbc  = ibcpt *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

              ibc2 = ibcpt2*krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
              jbc2 = ibcpt2*krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
              kbc2 = ibcpt2*krd(Ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

              in   = ibcn  *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
              jn   = ibcn  *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
              kn   = ibcn  *krd(ic1,3) + I2*krd(ic2,3) + i3*krd(ic3,3)

              in2  = ibcn2 *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
              jn2  = ibcn2 *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
              kn2  = ibcn2 *krd(ic1,3) + I2*krd(ic2,3) + i3*krd(ic3,3)

              mut(ibc,jbc,kbc,n) = &
                (9*mut(in ,jn ,kn ,n) - mut(in2,jn2,kn2,n)) / 8

            end do
          end do

        end if

      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine q_tec(nl,q,mut,x,y,z,xdot,ydot,zdot,dist,var1,var2)
!-----------------------------------------------------------------------
!     
!-----------------------------------------------------------------------
       
      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax
      integer(kind=cosa_int) nl,iblock,iq,imut,ixyz,idist,iv,npde1,npde2
      real (kind=cosa_real) q(*),mut(*),x(*),y(*),z(*),xdot(*),ydot(*),zdot(*), &
           dist(*),var1(*),var2(*)

      do iblock = 1,mynblocks
        imax   = i_imax     (iblock,nl)
        jmax   = j_jmax     (iblock,nl)
        kmax   = k_kmax     (iblock,nl)
        iq     = 1 + off_p3 (iblock,nl) * npde * dim5
        imut   = 1 + off_p3 (iblock,nl) *        dim5
        ixyz   = 1 + off_p2 (iblock,nl)        * dim5h
        idist  = 1 + off_p1 (iblock,nl)
        iv     = 1 + off_p1 (iblock,nl) * npde * dim5
        call q_tec_b(q(iq),mut(imut),x(ixyz),y(ixyz),z(ixyz),xdot(ixyz), &
          ydot(ixyz),zdot(ixyz),dist(idist),var1(iv),var2(iv),imax,jmax, &
          kmax,npde,nharms)
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine q_tec_b(q,mut,x,y,z,xdot,ydot,zdot,dist,var1,var2,imax, &
        jmax,kmax,npde,nharms)
!-----------------------------------------------------------------------
       
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,npde,nharms
      integer(kind=cosa_int) i,j,k,n,nh
      integer fid(0:2*mharms)
      real (kind=cosa_real) &
           q   (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           var1( 0:imax  , 0:jmax  , 0:kmax  ,npde,0:2*nharms), &
           var2( 0:imax  , 0:jmax  , 0:kmax  ,npde,0:2*nharms), &
           mut (-1:imax+1,-1:jmax+1,-1:kmax+1     ,0:2*nharms), &
           x   ( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           y   ( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           z   ( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           xdot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           ydot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           zdot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           dist( 0:imax  , 0:jmax  , 0:kmax)


!---- assemble var1 and var2 arrays

      do n = 0,2*nharms
        nh = n*hbmove

        if (kom.or.kom_bsl.or.kom_sst) then
          do k=0,kmax
            do j=0,jmax
              do i=0,imax
                var1(i,j,k,1,n) = x(i,j,k,nh)
                var1(i,j,k,2,n) = y(i,j,k,nh)
                var1(i,j,k,3,n) = z(i,j,k,nh)
                var1(i,j,k,4,n) = q(i,j,k,1,n)
                if (moving) then
                  var1(i,j,k,5,n) = q(i,j,k,2,n) - xdot(i,j,k,nh)
                  var1(i,j,k,6,n) = q(i,j,k,3,n) - ydot(i,j,k,nh)
                  var1(i,j,k,7,n) = q(i,j,k,4,n) - zdot(i,j,k,nh)
                else
                  var1(i,j,k,5,n) = q(i,j,k,2,n)
                  var1(i,j,k,6,n) = q(i,j,k,3,n)
                  var1(i,j,k,7,n) = q(i,j,k,4,n)
                end if
                var2(i,j,k,1,n) = q(i,j,k,5,n)
                var2(i,j,k,2,n) = gamma * var2(i,j,k,1,n) / &
                                  var1(i,j,k,4,n)
                var2(i,j,k,3,n) = &
                  sqrt(var1(i,j,k,5,n)**2+var1(i,j,k,6,n)**2+ &
                       var1(i,j,k,7,n)**2) / sqrt(var2(i,j,k,2,n))
                var2(i,j,k,4,n) = q(i,j,k,6,n)
                var2(i,j,k,5,n) = q(i,j,k,7,n)
                var2(i,j,k,6,n) = mut(i,j,k,n)
                var2(i,j,k,7,n) = dist(i,j,k)
              end do
            end do
          end do
        else
          do k=0,kmax
            do j=0,jmax
              do i=0,imax
                var1(i,j,k,1,n) = x(i,j,k,nh)
                var1(i,j,k,2,n) = y(i,j,k,nh)
                var1(i,j,k,3,n) = z(i,j,k,nh)
                var1(i,j,k,4,n) = q(i,j,k,1,n)
                if (moving) then
                  var1(i,j,k,5,n) = q(i,j,k,2,n) - xdot(i,j,k,nh)
                  var2(i,j,k,1,n) = q(i,j,k,3,n) - ydot(i,j,k,nh)
                  var2(i,j,k,2,n) = q(i,j,k,4,n) - zdot(i,j,k,nh)
                else
                  var1(i,j,k,5,n) = q(i,j,k,2,n)
                  var2(i,j,k,1,n) = q(i,j,k,3,n)
                  var2(i,j,k,2,n) = q(i,j,k,4,n)
                end if
                var2(i,j,k,3,n) = q(i,j,k,5,n)
                var2(i,j,k,4,n) = gamma * var2(i,j,k,3,n) / &
                                  var1(i,j,k,4,n)
                var2(i,j,k,5,n) = &
                  sqrt(var1(i,j,k,5,n)**2+var2(i,j,k,1,n)**2+ &
                       var2(i,j,k,2,n)**2) / sqrt(var2(i,j,k,4,n))
              end do
            end do
          end do
        end if

      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine force(nl,q,x,y,z,xdot,ydot,zdot,si,sj,sk,xideri, &
                       etaderj,zetaderk,bctopo,it)
!-----------------------------------------------------------------------
! Conventions:
!
!   AIRCRAFT. Includes stationary and pitching and/or plunging 
!      wings and aircraft configurations.
!      1) x is fuselage axis, y is altitude, wing is in xz plane.
!      2) force coefficients are normalized by product of absolute
!         dynamic head and reference area.
!      3) reference area for normalizing forces: lnght(1)*lnght(2).
!      4) additional length for normalizing moments: lnght(3).
!
!      cl   : lift force coefficient. Lift force has direction
!             aax, aay, aaz (see below).
!      cd   : drag force coefficient. Drag force has direction
!             ssx, ssy, ssz (see below).
!      cz   : lateral force coefficient. Lateral force has direction
!             bbx,bby,bbz.
!      cmx  : normalized x-component of moment (roll)
!      cmy  : normalized y-component of moment (yaw)
!      cmz  : normalized z-component of moment (pitch)
!
!-----------------------------------------------------------------------
!
!       VAWTs:
!       1) z is parallel to rotational axis.
!       2) force coefficients are normalized by product of absolute
!          dynamic head and reference area.
!       3) reference area is lnght(1)*lnght(2).
!       4) additional length at denominator of moment coefficient is
!          lnght(3).
!
!       cl   : normalized x-component of blade force in abs. frame
!       cd   : normalized y-component of blade force in abs. frame
!       cz   : normalized z-component of blade force in abs. frame
!       cmx  : normalized x-component of torque
!       cmy  : normalized y-component of torque
!       cm   : normalized z-component of torque
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,iblock,iblk,imax,jmax,kmax,n,iq,ixyz,iimt,ivmt,isurface, &
                ibodblk
      integer(kind=cosa_int) it
      integer fid(0:2*mharms)
      integer(kind=cosa_int) bctopo(*)
      integer iomode
      real (kind=cosa_real) q(*),x(*),y(*),z(*),xdot(*),ydot(*),zdot(*),si(*), &
           sj(*),sk(*),xideri(*),etaderj(*),zetaderk(*)
      real (kind=cosa_real) alphae(0:2*mharms), &
           cl_p   (0:2*mharms,msurface), cd_p   (0:2*mharms,msurface), &
           cz_p   (0:2*mharms,msurface), cmx_p  (0:2*mharms,msurface), &
           cmy_p  (0:2*mharms,msurface), cm_p   (0:2*mharms,msurface), &
           cl_v   (0:2*mharms,msurface), cd_v   (0:2*mharms,msurface), &
           cz_v   (0:2*mharms,msurface), cmx_v  (0:2*mharms,msurface), &
           cmy_v  (0:2*mharms,msurface), cm_v   (0:2*mharms,msurface), &
           cl_p_b (0:2*mharms,msurface), cd_p_b (0:2*mharms,msurface), &
           cz_p_b (0:2*mharms,msurface), cmx_p_b(0:2*mharms,msurface), &
           cmy_p_b(0:2*mharms,msurface), cm_p_b (0:2*mharms,msurface), &
           cl_v_b (0:2*mharms,msurface), cd_v_b (0:2*mharms,msurface), &
           cz_v_b (0:2*mharms,msurface), cmx_v_b(0:2*mharms,msurface), &
           cmy_v_b(0:2*mharms,msurface), cm_v_b (0:2*mharms,msurface)
      character*72 forcefile,tag,line
      integer(kind=cosa_int) tempblockindex
      logical amcontrol

      call amcontroller(amcontrol)

!---- set file names and create it/them, or open existing file/files

      do n=0,2*nharms

        fid(n)    = 200 + n
        if (aircraft) then
          alphae(n) = alpha - dtheta(n)
        else if (hawt.or.vawt) then
          alphae(n) = dtheta(n)
        end if

        if(amcontrol) then

        if (.not.harbal) then

          if ((unsteady.and.(itime.eq.0)).or. &
              ((.not.unsteady).and.(it.eq.1))) then
            open(fid(n),file='forces_tec.dat',status='replace', &
                 recl=1000)
            write(fid(n),*) 'TITLE="Force coefficients"'
            if (unsteady) then
              if (moving) then
                  if (nsurface.eq.1) then
                    write(fid(n),*) &
              'VARIABLES="itime","it","time","Cl_p","Cl_v","Cd_p","Cd_v","Cz_p","Cz_v","Cmx_p","Cmx_v","Cmy_p","Cmy_v","Cm_p","Cm &
     &_v","alpha"'
                  else if (nsurface.eq.2) then
                    write(fid(n),*) &
              'VARIABLES="itime","it","time","Cl1_p","Cl1_v","Cd1_p","Cd1_v","Cz1_p","Cz1_v","Cmx1_p","Cmx1_v","Cmy1_p","Cmy1_v", &
     &"Cm1_p","Cm1_v","Cl2_p","Cl2_v","Cd2_p","Cd2_v","Cz2_p","Cz2_v","Cmx2_p","Cmx2_v","Cmy2_p","Cmy2_v","Cm2_p","Cm2_v","alpha"'
                  else if (nsurface.eq.3) then
                    write(fid(n),*) &
              'VARIABLES="itime","it","time","Cl1_p","Cl1_v","Cd1_p","Cd1_v","Cz1_p","Cz1_v","Cmx1_p","Cmx1_v","Cmy1_p","Cmy1_v", &
  &"Cm1_p","Cm1_v","Cl2_p","Cl2_v","Cd2_p","Cd2_v","Cz2_p","Cz2_v","Cmx2_p","Cmx2_v","Cmy2_p","Cmy2_v","Cm2_p","Cm2_v","Cl3_p","C &
  &l3_v","Cd3_p","Cd3_v","Cz3_p","Cz3_v","Cmx3_p","Cmx3_v","Cmy3_p","Cmy3_v","Cm3_p","Cm3_v","alpha"'
                  else if (nsurface.eq.4) then
                    write(fid(n),*) &
             'VARIABLES="itime","it","time","Cl1_p","Cl1_v","Cd1_p","Cd1_v","Cz1_p","Cz1_v","Cmx1_p","Cmx1_v","Cmy1_p","Cmy1_v", &
  &"Cm1_p","Cm1_v","Cl2_p","Cl2_v","Cd2_p","Cd2_v","Cz2_p","Cz2_v","Cmx2_p","Cmx2_v","Cmy2_p","Cmy2_v","Cm2_p","Cm2_v","Cl3_p","C &
  &l3_v","Cd3_p","Cd3_v","Cz3_p","Cz3_v","Cmx3_p","Cmx3_v","Cmy3_p","Cmy3_v","Cm3_p","Cm3_v","Cl4_p","Cl4_v","Cd4_p","Cd4_v","Cz4 &
  &_p","Cz4_v","Cmx4_p","Cmx4_v","Cmy4_p","Cmy4_v","Cm4_p","Cm4_v","alpha"'
                  end if
              else
                  if (nsurface.eq.1) then
                    write(fid(n),*) &
              'VARIABLES="itime","it","time","Cl_p","Cl_v","Cd_p","Cd_v","Cz_p","Cz_v","Cmx_p","Cmx_v","Cmy_p","Cmy_v","Cm_p","Cm &
     &_v"'
                  else if (nsurface.eq.2) then
                    write(fid(n),*) &
     'VARIABLES="itime","it","time","it","Cl1_p","Cl1_v","Cd1_p","Cd1_v","Cz1_p","Cz1_v","Cmx1_p","Cmx1_v","Cmy1_p","Cmy1_v", &
     &"Cm1_p","Cm1_v","Cl2_p","Cl2_v","Cd2_p","Cd2_v","Cz2_p","Cz2_v","Cmx2_p","Cmx2_v","Cmy2_p","Cmy2_v","Cm2_p","Cm2_v"'
                  else if (nsurface.eq.3) then
                    write(fid(n),*) &
      'VARIABLES="itime","it","time","Cl1_p","Cl1_v","Cd1_p","Cd1_v","Cz1_p","Cz1_v","Cmx1_p","Cmx1_v","Cmy1_p","Cmy1_v", &
  &"Cm1_p","Cm1_v","Cl2_p","Cl2_v","Cd2_p","Cd2_v","Cz2_p","Cz2_v","Cmx2_p","Cmx2_v","Cmy2_p","Cmy2_v","Cm2_p","Cm2_v","Cl3_p","C &
  &l3_v","Cd3_p","Cd3_v","Cz3_p","Cz3_v","Cmx3_p","Cmx3_v","Cmy3_p","Cmy3_v","Cm3_p","Cm3_v"'
                  else if (nsurface.eq.4) then
                    write(fid(n),*) &
  'VARIABLES="itime","it","time","it","Cl1_p","Cl1_v","Cd1_p","Cd1_v","Cz1_p","Cz1_v","Cmx1_p","Cmx1_v","Cmy1_p","Cmy1_v", &
  &"Cm1_p","Cm1_v","Cl2_p","Cl2_v","Cd2_p","Cd2_v","Cz2_p","Cz2_v","Cmx2_p","Cmx2_v","Cmy2_p","Cmy2_v","Cm2_p","Cm2_v","Cl3_p","C &
  &l3_v","Cd3_p","Cd3_v","Cz3_p","Cz3_v","Cmx3_p","Cmx3_v","Cmy3_p","Cmy3_v","Cm3_p","Cm3_v","Cl4_p","Cl4_v","Cd4_p","Cd4_v","Cz4 &
  &_p","Cz4_v","Cmx4_p","Cmx4_v","Cmy4_p","Cmy4_v","Cm4_p","Cm4_v"'
                  end if
              end if
            else
                if (nsurface.eq.1) then
               write(fid(n),*) 'VARIABLES="it","Cl_p","Cl_v","Cd_p","Cd_v","Cz_p","Cz_v","Cmx_p","Cmx_v","Cmy_p","Cmy_v","Cm_p"," &
     &Cm_v"'
                else if (nsurface.eq.2) then
                  write(fid(n),*) &
                 'VARIABLES="it","Cl1_p","Cl1_v","Cd1_p","Cd1_v","Cz1_p","Cz1_v","Cmx1_p","Cmx1_v","Cmy1_p","Cmy1_v","Cm1_p","C &
     &m1_v","Cl2_p","Cl2_v","Cd2_p","Cd2_v","Cz2_p","Cz2_v","Cmx2_p","Cmx2_v","Cmy2_p","Cmy2_v","Cm2_p","Cm2_v"'
                else if (nsurface.eq.3) then
                  write(fid(n),*) &
                 'VARIABLES="it","Cl1_p","Cl1_v","Cd1_p","Cd1_v","Cz1_p","Cz1_v","Cmx1_p","Cmx1_v","Cmy1_p","Cmy1_v","Cm1_p","C &
  &m2_v","Cl2_p","Cl2_v","Cd2_p","Cd2_v","Cz2_p","Cz2_v","Cmx2_p","Cmx2_v","Cmy2_p","Cmy2_v","Cm2_p","Cm2_v","Cl3_p","Cl3_v","Cd3 &
  &_p","Cd3_v","Cz3_p","Cz3_v","Cmx3_p","Cmx3_v","Cmy3_p","Cmy3_v","Cm3_p","Cm3_v"'
                else if (nsurface.eq.4) then
                  write(fid(n),*) &
                'VARIABLES="it","Cl1_p","Cl1_v","Cd1_p","Cd1_v","Cz1_p","Cz1_v","Cmx1_p","Cmx1_v","Cmy1_p","Cmy1_v","Cm1_p","C &
  &m1_v","Cl2_p","Cl2_v","Cd2_p","Cd2_v","Cz2_p","Cz2_v","Cmx2_p","Cmx2_v","Cmy2_p","Cmy2_v","Cm2_p","Cm2_v","Cl3_p","Cl3_v","Cd3 &
  &_p","Cd3_v","Cz3_p","Cz3_v","Cmx3_p","Cmx3_v","Cmy3_p","Cmy3_v","Cm3_p","Cm3_v","Cl4_p","Cl4_v","Cd4_p","Cd4_v","Cz4_p","Cz4_v &
  &","Cmx4_p","Cmx4_v","Cmy4_p","Cmy4_v","Cm4_p","Cm4_v"'
                end if
            end if
          else if ((unsteady.and.(itime.gt.0)).or. &
              ((.not.unsteady).and.(it.gt.1))) then
            open(fid(n),file='forces_tec.dat',status='old', &
                 position='append')
          end if

        else if (harbal) then

          if (n.le.9) then
            write(tag,'(''00'',i1)') n
          else if (n.le.99) then
            write(tag,'(''0'',i2)') n
          else if (n.le.999) then
            write(tag,'(i3)') n
          end if

          forcefile = 'forces_hb_'//tag
          forcefile = trim(forcefile)//'.dat'
          if (it.eq.1) then
            open(fid(n),file=forcefile,status='replace',recl=1000)
            write(line,'('' TITLE="Force coefficients, Harmonic '',i2, &
     &                   ''"'')') n
            write(fid(n),'(a)') line
            if (moving) then
                if (nsurface.eq.1) then
                  write(fid(n),*) &
    'VARIABLES="it","time","Cl_p","Cl_v","Cd_p","Cd_v","Cz_p","Cz_v","Cmx_p","Cmx_v","Cmy_p","Cmy_v","Cm_p","Cm_v &
     &","alpha"'
                else if (nsurface.eq.2) then
                  write(fid(n),*) &
     'VARIABLES="it","time","Cl1_p","Cl1_v","Cd1_p","Cd1_v","Cz1_p","Cz1_v","Cmx1_p","Cmx1_v","Cmy1_p","Cmy1_v", &
     &"Cm_p","Cm_v","Cl2_p","Cl2_v","Cd2_p","Cd2_v","Cz2_p","Cz2_v","Cmx2_p","Cmx2_v","Cmy2_p","Cmy2_v","Cm2_p","Cm2_v","alpha"'
                else if (nsurface.eq.3) then
                  write(fid(n),*) &
      'VARIABLES="it","time","Cl1_p","Cl1_v","Cd1_p","Cd1_v","Cz1_p","Cz1_v","Cmx1_p","Cmx1_v","Cmy1_p","Cmy1_v", &
  &"Cm1_p","Cm1_v","Cl2_p","Cl2_v","Cd2_p","Cd2_v","Cz2_p","Cz2_v","Cmx2_p","Cmx2_v","Cmy2_p","Cmy2_v","Cm2_p","Cm2_v","Cl3_p","C &
  &l3_v","Cd3_p","Cd3_v","Cz3_p","Cz3_v","Cmx3_p","Cmx3_v","Cmy3_p","Cmy3_v","Cm3_p","Cm3_v","alpha"'
                else if (nsurface.eq.4) then
                  write(fid(n),*) &
    'VARIABLES="it","time","Cl1_p","Cl1_v","Cd1_p","Cd1_v","Cz1_p","Cz1_v","Cmx1_p","Cmx1_v","Cmy1_p","Cmy1_v", &
  &"Cm1_p","Cm1_v","Cl2_p","Cl2_v","Cd2_p","Cd2_v","Cz2_p","Cz2_v","Cmx2_p","Cmx2_v","Cmy2_p","Cmy2_v","Cm2_p","Cm2_v","Cl3_p","C &
  &l3_v","Cd3_p","Cd3_v","Cz3_p","Cz3_v","Cmx3_p","Cmx3_v","Cmy3_p","Cmy3_v","Cm3_p","Cm3_v","Cl4_p","Cl4_v","Cd4_p","Cd4_v","Cz4 &
  &_p","Cz4_v","Cmx4_p","Cmx4_v","Cmy4_p","Cmy4_v","Cm4_p","Cm4_v","alpha"'
                end if
            else
              write(fid(n),*) 'VARIABLES="it","time","Cl","Cd","Cz"'
            end if
          else if (it.gt.1) then
            open(fid(n),file=forcefile,status='old',position='append')
          end if

        end if

        end if

      end do

!---- initialize

      do isurface=1,nsurface
        do n=0,2*nharms
          cl_p (n,isurface) = 0.d0
          cd_p (n,isurface) = 0.d0
          cz_p (n,isurface) = 0.d0
          cmx_p(n,isurface) = 0.d0
          cmy_p(n,isurface) = 0.d0
          cm_p (n,isurface) = 0.d0
          cl_v (n,isurface) = 0.d0
          cd_v (n,isurface) = 0.d0
          cz_v (n,isurface) = 0.d0
          cmx_v(n,isurface) = 0.d0
          cmy_v(n,isurface) = 0.d0
          cm_v (n,isurface) = 0.d0
        end do
      end do

!---- compute forces

      do iblock = 1,mynblocks
        tempblockindex = (lowernblock - 1) + iblock
        imax   = i_imax     (iblock,nl)
        jmax   = j_jmax     (iblock,nl)
        kmax   = k_kmax     (iblock,nl)
        iblk   = 1 + off_bct(iblock,nl)
        iq     = 1 + off_p3 (iblock,nl) * npde * dim5
        ixyz   = 1 + off_p2 (iblock,nl)        * dim5h
        iimt   = 1 + off_p2 (iblock,nl) * lmet * dim5h
        ivmt   = 1 + off_0  (iblock,nl) * 3    * dim5h
        call &
          force_b(q(iq),x(ixyz),y(ixyz),z(ixyz),xdot(ixyz),ydot(ixyz), &
            zdot(ixyz),si(iimt),sj(iimt),sk(iimt),xideri(ivmt), &
            etaderj(ivmt),zetaderk(ivmt),nbcs(iblock),bctopo(iblk), &
            cl_p_b (0:2*nharms,1:nsurface),cd_p_b (0:2*nharms,1:nsurface), &
            cz_p_b (0:2*nharms,1:nsurface),cmx_p_b(0:2*nharms,1:nsurface), &
            cmy_p_b(0:2*nharms,1:nsurface),cm_p_b (0:2*nharms,1:nsurface), &
            cl_v_b (0:2*nharms,1:nsurface),cd_v_b (0:2*nharms,1:nsurface), &
            cz_v_b (0:2*nharms,1:nsurface),cmx_v_b(0:2*nharms,1:nsurface), &
            cmy_v_b(0:2*nharms,1:nsurface),cm_v_b (0:2*nharms,1:nsurface), &
            imax,jmax,kmax,npde,nharms,lmet)

        do isurface=1,nsurface
          do n=0,2*nharms
!
          cl_p (n,isurface) = cl_p (n,isurface) + cl_p_b (n,isurface)
          cd_p (n,isurface) = cd_p (n,isurface) + cd_p_b (n,isurface)
          cz_p (n,isurface) = cz_p (n,isurface) + cz_p_b (n,isurface)
          cmx_p(n,isurface) = cmx_p(n,isurface) + cmx_p_b(n,isurface)
          cmy_p(n,isurface) = cmy_p(n,isurface) + cmy_p_b(n,isurface)
          cm_p (n,isurface) = cm_p (n,isurface) + cm_p_b (n,isurface)
!
          cl_v (n,isurface) = cl_v (n,isurface) + cl_v_b (n,isurface)
          cd_v (n,isurface) = cd_v (n,isurface) + cd_v_b (n,isurface)
          cz_v (n,isurface) = cz_v (n,isurface) + cz_v_b (n,isurface)
          cmx_v(n,isurface) = cmx_v(n,isurface) + cmx_v_b(n,isurface)
          cmy_v(n,isurface) = cmy_v(n,isurface) + cmy_v_b(n,isurface)
          cm_v (n,isurface) = cm_v (n,isurface) + cm_v_b (n,isurface)
!
          end do
        end do

      end do

!---- write line in force file/files

      do n=0,2*nharms
        call combineforces(cl_p,cd_p ,cm_p ,n)
        call combineforces(cz_p,cmx_p,cmy_p,n)
        call combineforces(cl_v,cd_v ,cm_v ,n)
        call combineforces(cz_v,cmx_v,cmy_v,n)

!------ overall force coefficients

        if(amcontrol) then

          if (dualt.or.rgkuns) then
            if (moving) then
              write(fid(n),212) itime,it,simtime, &
                (cl_p(n,isurface),cl_v(n,isurface),cd_p(n,isurface), &
                 cd_v(n,isurface),cz_p(n,isurface),cz_v(n,isurface), &
                 cmx_p(n,isurface),cmx_v(n,isurface),cmy_p(n,isurface), &
                 cmy_v(n,isurface),cm_p(n,isurface),cm_v(n,isurface), &
                 isurface=1,nsurface),alphae(n)*180.d0/pi
            else
              write(fid(n),214) itime,it,simtime, &
                (cl_p(n,isurface),cl_v(n,isurface),cd_p(n,isurface), &
                 cd_v(n,isurface),cz_p(n,isurface),cz_v(n,isurface), &
                 cmx_p(n,isurface),cmx_v(n,isurface),cmy_p(n,isurface), &
                 cmy_v(n,isurface),cm_p(n,isurface),cm_v(n,isurface), &
                 isurface=1,nsurface)
            end if
          else if (.not.unsteady) then
            write(fid(n),216) it, &
              (cl_p(n,isurface),cl_v(n,isurface),cd_p(n,isurface),cd_v(n,isurface), &
               cz_p(n,isurface),cz_v(n,isurface),cmx_p(n,isurface), &
               cmx_v(n,isurface),cmy_p(n,isurface),cmy_v(n,isurface), &
               cm_p(n,isurface),cm_v(n,isurface),isurface=1,nsurface)
          else if (harbal) then
            if (moving) then
              write(fid(n),218) it,zeit(n), &
                (cl_p(n,isurface),cl_v(n,isurface),cd_p(n,isurface), &
                 cd_v(n,isurface),cz_p(n,isurface),cz_v(n,isurface), &
                 cmx_p(n,isurface),cmx_v(n,isurface),cmy_p(n,isurface), &
                 cmy_v(n,isurface),cm_p(n,isurface),cm_v(n,isurface), &
                 isurface=1,nsurface),alphae(n)*180.d0/pi
            else
              write(fid(n),220) it,zeit(n), &
                (cl_p(n,isurface),cl_v(n,isurface),cd_p(n,isurface), &
                 cd_v(n,isurface),cz_p(n,isurface),cz_v(n,isurface), &
                 cmx_p(n,isurface),cmx_v(n,isurface),cmy_p(n,isurface), &
                 cmy_v(n,isurface),cm_p(n,isurface),cm_v(n,isurface), &
                 isurface=1,nsurface)
            end if
          end if

        close(fid(n))

        end if

      end do

!212  format(2(i8,1x),50(1x,f22.14))
!214  format(2(i8,1x),49(1x,f22.14))
!216  format(i8,1x,48(1x,f22.14))
!218  format(i8,1x,50(1x,f22.14))
!220  format(i8,1x,49(1x,f22.14))

!212  format(2(i8,1x),50(1x,d12.5))
!214  format(2(i8,1x),49(1x,d12.5))
!216  format(i8,1x,48(1x,d12.5))
!218  format(i8,1x,50(1x,d12.5))
!220  format(i8,1x,49(1x,d12.5))

 212  format(2(i8,1x),(1x,d12.5),4(/12(1x,d12.5)),/(1x,d12.5))
 214  format(2(i8,1x),(1x,d12.5),4(/12(1x,d12.5)))
 216  format(i8,1x,4(/12(1x,d12.5)))
 218  format(i8,1x,(1x,d12.5),4(/12(1x,d12.5)),/(1x,d12.5))
 220  format(i8,1x,(1x,d12.5),4(/12(1x,d12.5)))

      return
      end

!-----------------------------------------------------------------------
      subroutine force_b(q,x,y,z,xdot,ydot,zdot,si,sj,sk,xideri,etaderj, &
        zetaderk,nbcs,bctopo,cl_p,cd_p,cz_p,cmx_p,cmy_p,cm_p,cl_v,cd_v, &
        cz_v,cmx_v,cmy_v,cm_v,imax,jmax,kmax,npde,nharms,lmet)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,npde,nharms,lmet,nbcs
      integer(kind=cosa_int) it,ibc1,isurface
      integer(kind=cosa_int) bctopo(10,nbcs)
      integer(kind=cosa_int) n,nh,i1,i2,i3,ijkmax(3),idir,istrt(3),iend(3),bctyp, &
        inrout,ibcpt,ibcpt2,ibcn,ibcn2,ibcm,ibcw,ibc,jbc,kbc,ibc2,jbc2, &
        kbc2,in,jn,kn,in2,jn2,kn2,im,jm,km,ic1,ic2,ic3,sysize,ioff1, &
        ioff2,joff1,joff2,koff1,koff2,iw,jw,kw
      real (kind=cosa_real) &
         q     (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
         x     ( 0:imax+1, 0:jmax+1, 0:kmax+1,     0:2*nharms*hbmove), &
         y     ( 0:imax+1, 0:jmax+1, 0:kmax+1,     0:2*nharms*hbmove), &
         z     ( 0:imax+1, 0:jmax+1, 0:kmax+1,     0:2*nharms*hbmove), &
         xdot  ( 0:imax+1, 0:jmax+1, 0:kmax+1,     0:2*nharms*hbmove), &
         ydot  ( 0:imax+1, 0:jmax+1, 0:kmax+1,     0:2*nharms*hbmove), &
         zdot  ( 0:imax+1, 0:jmax+1, 0:kmax+1,     0:2*nharms*hbmove), &
         si      (lmet, 0:imax+1, 0:jmax+1, 0:kmax+1,0:2*nharms*hbmove), &
         sj      (lmet, 0:imax+1, 0:jmax+1, 0:kmax+1,0:2*nharms*hbmove), &
         sk      (lmet, 0:imax+1, 0:jmax+1, 0:kmax+1,0:2*nharms*hbmove), &
         xideri  (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
         etaderj (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
         zetaderk(3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove)
      real (kind=cosa_real) fp(6),fv(6), &
           cl_p (0:2*nharms,nsurface), cd_p (0:2*nharms,nsurface), &
           cz_p (0:2*nharms,nsurface), cmx_p(0:2*nharms,nsurface), &
           cmy_p(0:2*nharms,nsurface), cm_p (0:2*nharms,nsurface), &
           cl_v (0:2*nharms,nsurface), cd_v (0:2*nharms,nsurface), &
           cz_v (0:2*nharms,nsurface), cmx_v(0:2*nharms,nsurface), &
           cmy_v(0:2*nharms,nsurface), cm_v (0:2*nharms,nsurface), &
           den,ssx,ssy,ssz,aax,aay,aaz,bbx,bby,bbz,nx,ny,nz,kx,ky,kz,ds, &
           sgnm,xarm,yarm,zarm,xmomr,ymomr,zmomr,xmomr_tmp,ymomr_tmp, &
           zmomr_tmp,dr1(2),dr2(2),cw,sw,ctp,stp,xrtw,yrtw,zrtw, &
           rhow,pw,tw,muw,u1,u2,uw,u1r,u2r,uwr,v1,v2,vw,v1r,v2r,vwr, &
           w1,w2,ww,w1r,w2r,wwr,ueta,veta,weta, &
           etax,etay,etaz,dudx,dudy,dudz,dvdx,dvdy, &
           dvdz,dwdx,dwdy,dwdz,divv,txx,txy,txz,tyy,tyz,tzz,mach
      character*72 forcefile,tag,line

!     store imax,jmax,kmax in ijkmax to enable boundary data location

      ijkmax(1) = imax
      ijkmax(2) = jmax
      ijkmax(3) = kmax

!---- mach is defined for force coefficient normalization. For viscous
      if (machfs.gt.1d-14) then
        mach = machfs
      else
        mach = 1.d0
      end if

      if (aircraft) then
        ssx = dcos(alpha) * dcos(beta)
        ssy = dsin(alpha) * dcos(beta)
        ssz = dsin(beta)
        den = dsqrt(dsin(beta)**2 + dcos(alpha)**2 * dcos(beta)**2)
        aax = - dsin(alpha) * dcos(alpha) * dcos(beta)**2 / den
        aay = den
        aaz = dsin(alpha) * dsin(beta) * dcos(beta) / den
        bbx = -dsin(beta) / den
        bby = 0.d0
        bbz = dcos(alpha) * dcos(beta) / den
      else if (hawt.or.vawt) then
        ssx = 0.d0
        ssy = 1.d0
        ssz = 0.d0
        aax = 1.d0
        aay = 0.d0
        aaz = 0.d0
        bbx = 0.d0
        bby = 0.d0
        bbz = 1.d0
      end if

      if (.not.harbal) zeit(0)=simtime

      do n=0,2*nharms
        nh = n*hbmove

        if (moving) then
          if (pitching) then
            dr2(1) =  (xmom - xh) * (dcos(dtheta(n)) - 1) - &
                      (ymom - yh) *  dsin(dtheta(n))
            dr2(2) =  (xmom - xh) *  dsin(dtheta(n))      + &
                      (ymom - yh) * (dcos(dtheta(n)) - 1)
            xmomr  =  xmom + dr2(1)
            ymomr  =  ymom + dr2(2)
            zmomr  =  zmom
          else if (plunging) then
            dr1(1) = dh0x * dsin(omegas*zeit(n))
            dr1(2) = dh0y * dsin(omegas*zeit(n))
            xmomr  = xmom + dr1(1)
            ymomr  = ymom + dr1(2)
            zmomr  =  zmom
          else if (plupitching) then
            dr1(1) = dh0x * dsin(omegas*zeit(n))
            dr1(2) = dh0y * dsin(omegas*zeit(n))
            dr2(1) =  (xmom - xh) * (dcos(dtheta(n)) - 1) - &
                      (ymom - yh) *  dsin(dtheta(n))
            dr2(2) =  (xmom - xh) *  dsin(dtheta(n))      + &
                      (ymom - yh) * (dcos(dtheta(n)) - 1)
            xmomr  = xmom + dr1(1) + dr2(1)
            ymomr  = ymom + dr1(2) + dr2(2)
            zmomr  =  zmom
          else if (rotating) then
            dr1(1) = ( (xmom - xrotc) * (dcos(dtheta(n)) - 1) - &
                       (ymom - yrotc) *  dsin(dtheta(n))      )
            dr1(2) = ( (xmom - xrotc) *  dsin(dtheta(n))      + &
                       (ymom - yrotc) * (dcos(dtheta(n)) - 1) )
            xmomr = xmom + dr1(1)
            ymomr = ymom + dr1(2)
            zmomr  =  zmom
            if (tpitch) then
              xmomr_tmp = xmomr
              ymomr_tmp = ymomr
              zmomr_tmp = zmomr
              cw    = dcos(betaw)
              sw    = dsin(betaw)
              xrtw  = -sw*zmomr_tmp  + cw*xmomr_tmp
              yrtw  =     ymomr_tmp
              zrtw  =  cw*zmomr_tmp  + sw*xmomr_tmp
              ctp   = dcos(thetatp(n))
              stp   = dsin(thetatp(n))
              xmomr_tmp = xrtw
              ymomr_tmp = yrtw &
                -(zrtw - zhtp)*stp + (yrtw - yhtp)*(ctp - 1)
              zmomr_tmp = zrtw &
                +(zrtw - zhtp)*(ctp-1) + (yrtw - yhtp)* stp
              xmomr = sw*zmomr_tmp + cw*xmomr_tmp
              ymomr = ymomr_tmp
              zmomr = cw*zmomr_tmp - sw*xmomr_tmp
            end if
          end if
        else
          xmomr = xmom
          ymomr = ymom
          zmomr = zmom
        end if

        do isurface=1,nsurface

          do i1=1,6
            fp(i1)=0.d0
            fv(i1)=0.d0
          end do

          cl_p(n,isurface)  = 0.d0
          cl_v(n,isurface)  = 0.d0
          cd_p(n,isurface)  = 0.d0
          cd_v(n,isurface)  = 0.d0
          cz_p(n,isurface)  = 0.d0
          cz_v(n,isurface)  = 0.d0
          cmx_p(n,isurface) = 0.d0
          cmx_v(n,isurface) = 0.d0
          cmy_p(n,isurface) = 0.d0
          cmy_v(n,isurface) = 0.d0
          cm_p (n,isurface) = 0.d0
          cm_v (n,isurface) = 0.d0

        do ibc1=1,nbcs

!         if ((bctopo(1,ibc1).eq.14).or.(bctopo(1,ibc1)/100.eq.15)) then
!         if                            (bctopo(1,ibc1)/100.eq.15)  then
          if (bctopo(1,ibc1).eq.1500+(isurface-1))  then

!           store boundary topology in mnemonic names
            bctyp    = bctopo(1,ibc1)
            idir     = bctopo(2,ibc1)
            inrout   = bctopo(3,ibc1)
            istrt(1) = bctopo(4,ibc1)
            iend(1)  = bctopo(5,ibc1)
            istrt(2) = bctopo(6,ibc1)
            iend(2)  = bctopo(7,ibc1)
            istrt(3) = bctopo(8,ibc1)
            iend(3)  = bctopo(9,ibc1)

!           set needed variables depending on whether the boundary is the
!           the inner boundary (inrout = 1) or 
!           the outer boundary (inrout > 1)
!               ibcpt : boundary condition location (first aux. cell)
!               ibcpt2: boundary condition location outside the block from
!                       ibcpt (second aux. cell)
!               ibcn  : point to the inside of the block from ibcpt
!               ibcm  : location of the metrics
!               ibcw  : location of the wall (coordinates)

            if (inrout.eq.1) then
              ibcpt  =  0
              ibcpt2 = -1
              ibcn   =  1
              ibcn2  =  2
              ibcm   =  1
              ibcw   =  0
              sgnm   = - 1.d0
            else
              ibcpt  = ijkmax(idir)
              ibcpt2 = ijkmax(idir) + 1
              ibcn   = ijkmax(idir) - 1
              ibcn2  = ijkmax(idir) - 2
              ibcm   = ijkmax(idir)
              ibcw   = ijkmax(idir)
              sgnm   = 1.d0
            end if

            if (idir.eq.1) then
              ioff1 = (-1)**(1+1/inrout)
              ioff2 = 2*ioff1
              joff1 = 0
              joff2 = 0
              koff1 = 0
              koff2 = 0
            else if (idir.eq.2) then
              ioff1 = 0
              ioff2 = 0
              joff1 = (-1)**(1+1/inrout)
              joff2 = 2*joff1
              koff1 = 0
              koff2 = 0
            else if (idir.eq.3) then
              ioff1 = 0
              ioff2 = 0
              joff1 = 0
              joff2 = 0
              koff1 = (-1)**(1+1/inrout)
              koff2 = 2*koff1
            end if

            ic1 = cyc (idir, 1)
            ic2 = cyc (idir, 2)
            ic3 = cyc (idir, 3)

!---------- calculate pressure force
            do i3 = istrt(ic3),iend(ic3)
              do i2 = istrt(ic2),iend(ic2)

                ibc  = ibcpt *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
                jbc  = ibcpt *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
                kbc  = ibcpt *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

                ibc2 = ibcpt2*krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
                jbc2 = ibcpt2*krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
                kbc2 = ibcpt2*krd(Ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

                in   = ibcn  *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
                jn   = ibcn  *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
                kn   = ibcn  *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

                in2  = ibcn2 *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
                jn2  = ibcn2 *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
                kn2  = ibcn2 *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

                im   = ibcm  *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
                jm   = ibcm  *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
                km   = ibcm  *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

                iw   = ibcw  *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
                jw   = ibcw  *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
                kw   = ibcw  *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

                if (idir.eq.1) then
                  nx    = si(1,im,jm,km,nh)
                  ny    = si(2,im,jm,km,nh)
                  nz    = si(3,im,jm,km,nh)
                  ds    = si(4,im,jm,km,nh)
                else if (idir.eq.2) then
                  nx    = sj(1,im,jm,km,nh)
                  ny    = sj(2,im,jm,km,nh)
                  nz    = sj(3,im,jm,km,nh)
                  ds    = sj(4,im,jm,km,nh)
                else if (idir.eq.3) then
                  nx    = sk(1,im,jm,km,nh)
                  ny    = sk(2,im,jm,km,nh)
                  nz    = sk(3,im,jm,km,nh)
                  ds    = sk(4,im,jm,km,nh)
                end if

                kx   = (-1)**(1+1/inrout) * nx
                ky   = (-1)**(1+1/inrout) * ny
                kz   = (-1)**(1+1/inrout) * nz
                pw   = (gamma-1) * (q(ibc,jbc,kbc,5,n) - &
                     q(ibc,jbc,kbc,2,n)**2 / q(ibc,jbc,kbc,1,n) / 2 - &
                     q(ibc,jbc,kbc,3,n)**2 / q(ibc,jbc,kbc,1,n) / 2 - &
                     q(ibc,jbc,kbc,4,n)**2 / q(ibc,jbc,kbc,1,n) / 2 )
                fp(1) = fp(1) - pw * kx * ds
                fp(2) = fp(2) - pw * ky * ds
                fp(3) = fp(3) - pw * kz * ds
                xarm = x(iw,jw,kw,nh) - xmomr
                yarm = y(iw,jw,kw,nh) - ymomr
                zarm = z(iw,jw,kw,nh) - zmomr
                fp(4) = fp(4) + (ky * zarm - kz *yarm) * pw * ds
                fp(5) = fp(5) + (kz * xarm - kx *zarm) * pw * ds
                fp(6) = fp(6) + (kx * yarm - ky *xarm) * pw * ds

              end do
            end do


!old        if ((viscous).and.(bctopo(1,ibc1)/100.eq.15)) then
!old        if (viscous.and.(bctopo(1,ibc1).eq.1500+(isurface-1)))  then
            if (viscous) then
!------------ calculate viscous force

              do i3 = istrt(ic3),iend(ic3)
                do i2 = istrt(ic2),iend(ic2)

                  ibc  = ibcpt *krd(ic1,1) +i2*krd(ic2,1) +i3*krd(ic3,1)
                  jbc  = ibcpt *krd(ic1,2) +i2*krd(ic2,2) +i3*krd(ic3,2)
                  kbc  = ibcpt *krd(ic1,3) +i2*krd(ic2,3) +i3*krd(ic3,3)

                  ibc2 = ibcpt2*krd(ic1,1) +i2*krd(ic2,1) +i3*krd(ic3,1)
                  jbc2 = ibcpt2*krd(ic1,2) +i2*krd(ic2,2) +i3*krd(ic3,2)
                  kbc2 = ibcpt2*krd(Ic1,3) +i2*krd(ic2,3) +i3*krd(ic3,3)

                  in   = ibcn  *krd(ic1,1) +i2*krd(ic2,1) +i3*krd(ic3,1)
                  jn   = ibcn  *krd(ic1,2) +i2*krd(ic2,2) +i3*krd(ic3,2)
                  kn   = ibcn  *krd(ic1,3) +i2*krd(ic2,3) +i3*krd(ic3,3)

                  in2  = ibcn2 *krd(ic1,1) +i2*krd(ic2,1) +i3*krd(ic3,1)
                  jn2  = ibcn2 *krd(ic1,2) +i2*krd(ic2,2) +i3*krd(ic3,2)
                  kn2  = ibcn2 *krd(ic1,3) +i2*krd(ic2,3) +i3*krd(ic3,3)

                  im   = ibcm  *krd(ic1,1) +i2*krd(ic2,1) +i3*krd(ic3,1)
                  jm   = ibcm  *krd(ic1,2) +i2*krd(ic2,2) +i3*krd(ic3,2)
                  km   = ibcm  *krd(ic1,3) +i2*krd(ic2,3) +i3*krd(ic3,3)

                  iw   = ibcw  *krd(ic1,1) +i2*krd(ic2,1) +i3*krd(ic3,1)
                  jw   = ibcw  *krd(ic1,2) +i2*krd(ic2,2) +i3*krd(ic3,2)
                  kw   = ibcw  *krd(ic1,3) +i2*krd(ic2,3) +i3*krd(ic3,3)

                  if (idir.eq.1) then
                    nx    = si(1,im,jm,km,nh)
                    ny    = si(2,im,jm,km,nh)
                    nz    = si(3,im,jm,km,nh)
                    ds    = si(4,im,jm,km,nh)
                    etax  = xideri(1,im,jm,km,nh)
                    etay  = xideri(2,im,jm,km,nh)
                    etaz  = xideri(3,im,jm,km,nh)
                  else if (idir.eq.2) then
                    nx    = sj(1,im,jm,km,nh)
                    ny    = sj(2,im,jm,km,nh)
                    nz    = sj(3,im,jm,km,nh)
                    ds    = sj(4,im,jm,km,nh)
                    etax  = etaderj(1,im,jm,km,nh)
                    etay  = etaderj(2,im,jm,km,nh)
                    etaz  = etaderj(3,im,jm,km,nh)
                  else if (idir.eq.3) then
                    nx    = sk(1,im,jm,km,nh)
                    ny    = sk(2,im,jm,km,nh)
                    nz    = sk(3,im,jm,km,nh)
                    ds    = sk(4,im,jm,km,nh)
                    etax  = zetaderk(1,im,jm,km,nh)
                    etay  = zetaderk(2,im,jm,km,nh)
                    etaz  = zetaderk(3,im,jm,km,nh)
                  end if

                  kx   = (-1)**(1+1/inrout) * nx
                  ky   = (-1)**(1+1/inrout) * ny
                  kz   = (-1)**(1+1/inrout) * nz
                  rhow = q(ibc,jbc,kbc,1,n)
                  pw   = (gamma-1) * (q(ibc,jbc,kbc,5,n) - &
                       q(ibc,jbc,kbc,2,n)**2 / q(ibc,jbc,kbc,1,n) / 2 - &
                       q(ibc,jbc,kbc,3,n)**2 / q(ibc,jbc,kbc,1,n) / 2 - &
                       q(ibc,jbc,kbc,4,n)**2 / q(ibc,jbc,kbc,1,n) / 2 )
                  tw   = pw *gamma / rhow
                  muw  = (stemp+1)/(tw+stemp) * (tw**1.5d0)

                  uw   = q(ibc,jbc,kbc,2,n) / q(ibc,jbc,kbc,1,n)
                  u1   = q(ibc+ioff1,jbc+joff1,kbc+koff1,2,n) / &
                         q(ibc+ioff1,jbc+joff1,kbc+koff1,1,n)
                  u2   = q(ibc+ioff2,jbc+joff2,kbc+koff2,2,n) / &
                         q(ibc+ioff2,jbc+joff2,kbc+koff2,1,n)
                  if (moving) then
                    uwr  = xdot(ibc      ,jbc      ,kbc      ,nh)
                    u1r  = xdot(ibc+ioff1,jbc+joff1,kbc+koff1,nh)
                    u2r  = xdot(ibc+ioff2,jbc+joff2,kbc+koff2,nh)
                    uw   = uw - uwr
                    u1   = u1 - u1r
                    u2   = u2 - u2r
                  end if

                  vw   = q(ibc,jbc,kbc,3,n) / q(ibc,jbc,kbc,1,n)
                  v1   = q(ibc+ioff1,jbc+joff1,kbc+koff1,3,n) / &
                         q(ibc+ioff1,jbc+joff1,kbc+koff1,1,n)
                  v2   = q(ibc+ioff2,jbc+joff2,kbc+koff2,3,n) / &
                         q(ibc+ioff2,jbc+joff2,kbc+koff2,1,n)
                  if (moving) then
                    vwr  = ydot(ibc      ,jbc      ,kbc      ,nh)
                    v1r  = ydot(ibc+ioff1,jbc+joff1,kbc+koff1,nh)
                    v2r  = ydot(ibc+ioff2,jbc+joff2,kbc+koff2,nh)
                    vw   = vw - vwr
                    v1   = v1 - v1r
                    v2   = v2 - v2r
                  end if

                  ww   = q(ibc,jbc,kbc,4,n) / q(ibc,jbc,kbc,1,n)
                  w1   = q(ibc+ioff1,jbc+joff1,kbc+koff1,4,n) / &
                         q(ibc+ioff1,jbc+joff1,kbc+koff1,1,n)
                  w2   = q(ibc+ioff2,jbc+joff2,kbc+koff2,4,n) / &
                         q(ibc+ioff2,jbc+joff2,kbc+koff2,1,n)
                  if (moving) then
                    wwr  = zdot(ibc      ,jbc      ,kbc      ,nh)
                    w1r  = zdot(ibc+ioff1,jbc+joff1,kbc+koff1,nh)
                    w2r  = zdot(ibc+ioff2,jbc+joff2,kbc+koff2,nh)
                    ww   = ww - wwr
                    w1   = w1 - w1r
                    w2   = w2 - w2r
                  end if

                  ueta = (-1)**(1+1/inrout) * (-8*uw + 9*u1 - u2) / 3
                  veta = (-1)**(1+1/inrout) * (-8*vw + 9*v1 - v2) / 3
                  weta = (-1)**(1+1/inrout) * (-8*ww + 9*w1 - w2) / 3
                  dudx = ueta*etax
                  dudy = ueta*etay
                  dudz = ueta*etaz
                  dvdx = veta*etax
                  dvdy = veta*etay
                  dvdz = veta*etaz
                  dwdx = weta*etax
                  dwdy = weta*etay
                  dwdz = weta*etaz
                  divv = dudx + dvdy + dwdz
                  txx  = 2*muw * (dudx - divv/3)
                  txy  = muw*(dudy + dvdx)
                  txz  = muw*(dudz + dwdx)
                  tyy  = 2*muw * (dvdy - divv/3)
                  tyz  = muw*(dvdz + dwdy)
                  tzz  = 2*muw * (dwdz - divv/3)

                  fv(1) = fv(1) + machfs / reyno * ds * &
                          (txx*kx + txy*ky + txz*kz)
                  fv(2) = fv(2) + machfs / reyno * ds * &
                          (txy*kx + tyy*ky + tyz*kz)
                  fv(3) = fv(3) + machfs / reyno * ds * &
                          (txz*kx + tyz*ky + tzz*kz)
                  xarm = x(iw,jw,kw,nh) - xmomr
                  yarm = y(iw,jw,kw,nh) - ymomr
                  zarm = z(iw,jw,kw,nh) - zmomr
                  fv(4) = fv(4) + machfs / reyno * ds * &
                    (-(txy*kx + tyy*ky + tyz*kz) * zarm + &
                      (txz*kx + tyz*ky + tzz*kz) * yarm )
                  fv(5) = fv(5) + machfs / reyno * ds * &
                    (-(txz*kx + tyz*ky + tzz*kz) * xarm + &
                      (txx*kx + txy*ky + txz*kz) * zarm )
                  fv(6) = fv(6) + machfs / reyno * ds * &
                    (-(txx*kx + txy*ky + txz*kz) * yarm + &
                      (txy*kx + tyy*ky + tyz*kz) * xarm )

                end do
              end do

            end if

          end if

!------ end loop (ibc1=1,nbcs)
        end do

        cl_p(n,isurface) = (fp(1)*aax+fp(2)*aay+fp(3)*aaz) * &
                2/(lngth(1)*lngth(2)*mach**2)
        cd_p(n,isurface) = (fp(1)*ssx+fp(2)*ssy+fp(3)*ssz) * &
                2/(lngth(1)*lngth(2)*mach**2)
        cz_p(n,isurface) = (fp(1)*bbx+fp(2)*bby+fp(3)*bbz) * &
                2/(lngth(1)*lngth(2)*mach**2)
        cmx_p(n,isurface) = fp(4) * &
                   2/(lngth(1)*lngth(2)*lngth(3)*mach**2)
        cmy_p(n,isurface) = fp(5) * &
                   2/(lngth(1)*lngth(2)*lngth(3)*mach**2)
        cm_p (n,isurface) = fp(6) * &
                   2/(lngth(1)*lngth(2)*lngth(3)*mach**2)
        cl_v(n,isurface) = (fv(1)*aax+fv(2)*aay+fv(3)*aaz) * &
                2/(lngth(1)*lngth(2)*mach**2)
        cd_v(n,isurface) = (fv(1)*ssx+fv(2)*ssy+fv(3)*ssz) * &
                2/(lngth(1)*lngth(2)*mach**2)
        cz_v(n,isurface) = (fv(1)*bbx+fv(2)*bby+fv(3)*bbz) * &
                2/(lngth(1)*lngth(2)*mach**2)
        cmx_v(n,isurface) = fv(4) * &
                   2/(lngth(1)*lngth(2)*lngth(3)*mach**2)
        cmy_v(n,isurface) = fv(5) * &
                   2/(lngth(1)*lngth(2)*lngth(3)*mach**2)
        cm_v (n,isurface) = fv(6) * &
                   2/(lngth(1)*lngth(2)*lngth(3)*mach**2)

!---- end loop (isurface=1,nsurface)
      end do

!---- end loop (n=0,2*nharms)
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine output(nl,it,done)
!-----------------------------------------------------------------------
       
      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax
      integer(kind=cosa_int) it,nl,idir,ipde,iprintf
      real (kind=cosa_real) q,q1,pr,qp,qm,dq,dqp,dqm,flux,x,y,z,si,sj,sk,xdot, &
           ydot,zdot,xideri,etaderj,zetaderk,dist,mut
      integer(kind=cosa_int) bctopo
      logical done
      logical amcontrol
      pointer &
       (pq  ,  q), (pq1 , q1), (pqp , qp), (pqm , qm), (pdq , dq), &
       (pdqp,dqp), (pdqm,dqm), (ppr , pr), (px  ,  x), (py  ,  y), &
       (pz  ,  z), (psi , si), (psj , sj), (psk , sk), &
       (pflux , flux), (pdist , dist), &
       (pxdot , xdot), (pydot , ydot), (pzdot , zdot), &
       (pbctopo  ,  bctopo), (pxideri  ,  xideri), (petaderj , etaderj), &
       (pzetaderk,zetaderk), (pmut, mut)

      imax      = i_imax(1,nl)
      jmax      = j_jmax(1,nl)
      kmax      = k_kmax(1,nl)
      pq        = p_q(nl)
      pq1       = p_q1(nl)
      pmut       = p_mut(nl)
      pqp       = p_qp(nl)
      pqm       = p_qm(nl)
      pdq       = p_dq(nl)
      pdqp      = p_dqp(nl)
      pdqm      = p_dqm(nl)
      pflux     = p_flux(nl)
      ppr       = p_pr(nl)
      px        = p_x(nl)
      py        = p_y(nl)
      pz        = p_z(nl)
      pxdot     = p_xdot(nl)
      pydot     = p_ydot(nl)
      pzdot     = p_zdot(nl)
      psi       = p_si(nl)
      psj       = p_sj(nl)
      psk       = p_sk(nl)
      pxideri   = p_xideri(nl)
      petaderj  = p_etaderj(nl)
      pzetaderk = p_zetaderk(nl)
      pbctopo   = p_bctopo(nl)
      pdist     = p_dist(nl)

      iprintf = iprint

      call amcontroller(amcontrol)

!---- write line in hist.dat file
      if (((it.eq.1).or.(it/iprint*iprint.eq.it).or.done.or. &
          (it.eq.ncycle)).and. &
          (it.ne.0)) then
        if(amcontrol) then
          write(37,102) itime,it,(resrms(ipde),ipde=1,5)
          flush(37)
          write(* ,102) itime,it,(resrms(ipde),ipde=1,5)
          if (kom.or.kom_bsl.or.kom_sst) then
            write(41,102) itime,it,(resrms(ipde),ipde=6,npde)
            flush(41)
          end if
        end if
      end if

!---- write line in histg.dat file
      if ((done.or.(it.eq.ncycle)).and.(it.ne.0)) then
        if(amcontrol) then
          write(39,102) itime,it,(resrms(ipde),ipde=1,5)
          flush(39)
!tmp      if (kom.or.kom_bsl.or.kom_sst) then
!tmp        write(xx,102) itime,it,(resrms(ipde),ipde=6,npde)
!tmp      end if
        end if
      end if

!---- write line in forces_tec.dat file
      if (calfor) then
        if (((.not.unsteady).or.harbal).and.(it.ne.0).and. &
            ((it.eq.1).or.(it/iprintf*iprintf.eq.it).or.done.or. &
             (it.eq.ncycle)))                                       then
          call force(nl,q,x,y,z,xdot,ydot,zdot,si,sj,sk,xideri,etaderj, &
                     zetaderk,bctopo,it)
        else if ((dualt.or.rgkuns).and. &
                 ((itime.eq.0).or.(it.eq.1).or. &
                  (it/iprint*iprint.eq.it).or.done.or. &
                  (it.eq.ncycle))) then
          call force(nl,q,x,y,z,xdot,ydot,zdot,si,sj,sk,xideri,etaderj, &
                     zetaderk,bctopo,it)
        end if
      end if

!---- write restart file
      if ((wri_res).and.(it.ne.0)) then
        if ((.not.dualt).and. &
            (done.or.(it/srest*srest.eq.it).or.(it.eq.ncycle)).or. &
            (dualt.and. &
             (done.or. &
              ((itime.eq.ntime).and.(it.eq.ncycle)).or. &
              ((itime/nsave*nsave.eq.itime).and.(it.eq.ncycle)).or. &
              (it/srest*srest.eq.it)))) then
        if(amcontrol) then
          call zeitdata_out
        end if
        call writerest(nl,q,q1,mut)
        end if
      end if

 100  format(i8,4(1x,e13.5))
 102  format(i7,2x,i7,5(1x,e13.5))

      return
      end

!-----------------------------------------------------------------------
      subroutine vert2cell_g(nl,x,y,z,workx,worky,workz)
!-----------------------------------------------------------------------
       
      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,iblock,imax,jmax,kmax,ixyz,dim5l,dim5u
      real (kind=cosa_real) x(*),y(*),z(*),workx(*),worky(*),workz(*)

      dim5l = 0
      dim5u = 2 * nharms * hbmove
      call copy_array(1,1,dim5l,dim5u, 2,nl,workx,x, &
                      1,1,dim5l,dim5u, 0)
      call copy_array(1,1,dim5l,dim5u, 2,nl,worky,y, &
                      1,1,dim5l,dim5u, 0)
      call copy_array(1,1,dim5l,dim5u, 2,nl,workz,z, &
                      1,1,dim5l,dim5u, 0)

      do iblock = 1,mynblocks
        imax   = i_imax     (iblock,nl)
        jmax   = j_jmax     (iblock,nl)
        kmax   = k_kmax     (iblock,nl)
        ixyz   = 1 + off_p2 (iblock,nl)        * dim5h
        call vert2cell_bg(x(ixyz),y(ixyz),z(ixyz),workx(ixyz), &
          worky(ixyz),workz(ixyz),imax,jmax,kmax,npde,nharms)
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine vert2cell_bg(x,y,z,workx,worky,workz,imax,jmax,kmax, &
        npde,nharms)
!-----------------------------------------------------------------------
       
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,npde,nharms,nl
      integer(kind=cosa_int) i,im,i1,j,jm,j1,k,km,k1,n
      real (kind=cosa_real) &
           x     (0:imax+1,0:jmax+1,0:kmax+1,0:2*nharms*hbmove), &
           y     (0:imax+1,0:jmax+1,0:kmax+1,0:2*nharms*hbmove), &
           z     (0:imax+1,0:jmax+1,0:kmax+1,0:2*nharms*hbmove), &
           workx (0:imax+1,0:jmax+1,0:kmax+1,0:2*nharms*hbmove), &
           worky (0:imax+1,0:jmax+1,0:kmax+1,0:2*nharms*hbmove), &
           workz (0:imax+1,0:jmax+1,0:kmax+1,0:2*nharms*hbmove)

      im = imax-1
      jm = jmax-1
      km = kmax-1

      do n = 0,2*nharms*hbmove

!------ block vertices

        do k=0,kmax,kmax
          if (k.eq.0) then
            k1 = 1
          else if (k.eq.kmax) then
            k1 = kmax
          end if
          x(0   ,0   ,k,n) = workx(1   ,1   ,k1   ,n)
          y(0   ,0   ,k,n) = worky(1   ,1   ,k1   ,n)
          z(0   ,0   ,k,n) = workz(1   ,1   ,k1   ,n)
          x(imax,0   ,k,n) = workx(imax,1   ,k1   ,n)    
          y(imax,0   ,k,n) = worky(imax,1   ,k1   ,n)
          z(imax,0   ,k,n) = workz(imax,1   ,k1   ,n)
          x(0   ,jmax,k,n) = workx(1   ,jmax,k1   ,n)
          y(0   ,jmax,k,n) = worky(1   ,jmax,k1   ,n)
          z(0   ,jmax,k,n) = workz(1   ,jmax,k1   ,n)
          x(imax,jmax,k,n) = workx(imax,jmax,k1   ,n)
          y(imax,jmax,k,n) = worky(imax,jmax,k1   ,n)
          z(imax,jmax,k,n) = workz(imax,jmax,k1   ,n)
        end do

!------ block edges

        do k=0,kmax,kmax
          if (k.eq.0) then
            k1 = 1
          else if (k.eq.kmax) then
            k1 = kmax
          end if
          do j=0,jmax,jmax
            if (j.eq.0) then
              j1 = 1
            else if (j.eq.jmax) then
              j1 = jmax
            end if
            do i = 1,im
              x(i,j,k,n) = (workx(i,j1,k1,n) + workx(i+1,j1,k1,n)) / 2
              y(i,j,k,n) = (worky(i,j1,k1,n) + worky(i+1,j1,k1,n)) / 2
              z(i,j,k,n) = (workz(i,j1,k1,n) + workz(i+1,j1,k1,n)) / 2
            end do
          end do
        end do

        do k=0,kmax,kmax
          if (k.eq.0) then
            k1 = 1
          else if (k.eq.kmax) then
            k1 = kmax
          end if
          do i=0,imax,imax
            if (i.eq.0) then
              i1 = 1
            else if (i.eq.imax) then
              i1 = imax
            end if
            do j = 1,jm
              x(i,j,k,n) = (workx(i1,j,k1,n) + workx(i1,j+1,k1,n)) / 2
              y(i,j,k,n) = (worky(i1,j,k1,n) + worky(i1,j+1,k1,n)) / 2
              z(i,j,k,n) = (workz(i1,j,k1,n) + workz(i1,j+1,k1,n)) / 2
            end do
          end do
        end do

        do j=0,jmax,jmax
          if (j.eq.0) then
            j1 = 1
          else if (j.eq.jmax) then
            j1 = jmax
          end if
          do i=0,imax,imax
            if (i.eq.0) then
              i1 = 1
            else if (i.eq.imax) then
              i1 = imax
            end if
            do k = 1,km
              x(i,j,k,n) = (workx(i1,j1,k,n) + workx(i1,j1,k+1,n)) / 2
              y(i,j,k,n) = (worky(i1,j1,k,n) + worky(i1,j1,k+1,n)) / 2
              z(i,j,k,n) = (workz(i1,j1,k,n) + workz(i1,j1,k+1,n)) / 2
            end do
          end do
        end do

!------ block bounding surfaces

        do k=0,kmax,kmax
          if (k.eq.0) then
            k1 = 1
          else if (k.eq.kmax) then
            k1 = kmax
          end if
          do j = 1,jm
            do i = 1,im
              x(i,j,k,n) = &
                (workx(i  ,j  ,k1,n) + workx(i+1,j  ,k1,n) + &
                 workx(i  ,j+1,k1,n) + workx(i+1,j+1,k1,n)) / 4
              y(i,j,k,n) = &
                (worky(i  ,j  ,k1,n) + worky(i+1,j  ,k1,n) + &
                 worky(i  ,j+1,k1,n) + worky(i+1,j+1,k1,n)) / 4
              z(i,j,k,n) = &
                (workz(i  ,j  ,k1,n) + workz(i+1,j  ,k1,n) + &
                 workz(i  ,j+1,k1,n) + workz(i+1,j+1,k1,n)) / 4
            end do
          end do
        end do

        do j=0,jmax,jmax
          if (j.eq.0) then
            j1 = 1
          else if (j.eq.jmax) then
            j1 = jmax
          end if
          do k = 1,km
            do i = 1,im
              x(i,j,k,n) = &
                (workx(i  ,j1,k  ,n) + workx(i+1,j1,k  ,n) + &
                 workx(i  ,j1,k+1,n) + workx(i+1,j1,k+1,n)) / 4
              y(i,j,k,n) = &
                (worky(i  ,j1,k  ,n) + worky(i+1,j1,k  ,n) + &
                 worky(i  ,j1,k+1,n) + worky(i+1,j1,k+1,n)) / 4
              z(i,j,k,n) = &
                (workz(i  ,j1,k  ,n) + workz(i+1,j1,k  ,n) + &
                 workz(i  ,j1,k+1,n) + workz(i+1,j1,k+1,n)) / 4
            end do
          end do
        end do

        do i=0,imax,imax
          if (i.eq.0) then
            i1 = 1
          else if (i.eq.imax) then
            i1 = imax
          end if
          do k = 1,km
            do j = 1,jm
              x(i,j,k,n) = &
                (workx(i1,j  ,k  ,n) + workx(i1,j  ,k+1,n) + &
                 workx(i1,j+1,k  ,n) + workx(i1,j+1,k+1,n)) / 4
              y(i,j,k,n) = &
                (worky(i1,j  ,k  ,n) + worky(i1,j  ,k+1,n) + &
                 worky(i1,j+1,k  ,n) + worky(i1,j+1,k+1,n)) / 4
              z(i,j,k,n) = &
                (workz(i1,j  ,k  ,n) + workz(i1,j  ,k+1,n) + &
                 workz(i1,j+1,k  ,n) + workz(i1,j+1,k+1,n)) / 4
            end do
          end do
        end do

!------ internal nodes

        do k = 1,km
          do j = 1,jm
            do i = 1,im
              x(i,j,k,n) = &
                (workx(i  ,j  ,k  ,n) + workx(i+1,j  ,k  ,n) + &
                 workx(i  ,j+1,k  ,n) + workx(i+1,j+1,k  ,n) + &
                 workx(i  ,j  ,k+1,n) + workx(i+1,j  ,k+1,n) + &
                 workx(i  ,j+1,k+1,n) + workx(i+1,j+1,k+1,n) ) / 8
              y(i,j,k,n) = &
                (worky(i  ,j  ,k  ,n) + worky(i+1,j  ,k  ,n) + &
                 worky(i  ,j+1,k  ,n) + worky(i+1,j+1,k  ,n) + &
                 worky(i  ,j  ,k+1,n) + worky(i+1,j  ,k+1,n) + &
                 worky(i  ,j+1,k+1,n) + worky(i+1,j+1,k+1,n) ) / 8
              z(i,j,k,n) = &
                (workz(i  ,j  ,k  ,n) + workz(i+1,j  ,k  ,n) + &
                 workz(i  ,j+1,k  ,n) + workz(i+1,j+1,k  ,n) + &
                 workz(i  ,j  ,k+1,n) + workz(i+1,j  ,k+1,n) + &
                 workz(i  ,j+1,k+1,n) + workz(i+1,j+1,k+1,n) ) / 8
            end do
          end do
        end do

      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine cutman_qtec(nl,cutopo,q)
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,iblock,iblk,icut,imax,jmax,kmax,iq
      integer(kind=cosa_int) cutopo(21,ncuts)
      logical myblockl
      real(kind=cosa_real) q(*)

!     do iblock = 1,mynblockS
        do icut = 1,myncuts
          iblk = cutopo( 1,icut)
          call myblock(iblk,myblockl,.true.)
          if (myblockl) then
            imax   = i_imax     (iblk,nl)
            jmax   = j_jmax     (iblk,nl)
            kmax   = k_kmax     (iblk,nl)
            iq     = 1 + off_p3 (iblk,nl) * npde * dim5
            call cut_qtec(imax,jmax,kmax,q(iq),npde,nharms, &
                          cutopo(1,icut))
          end if
        end do
!     end do

      return
      end
!-----------------------------------------------------------------------
      subroutine pcutman_qtec(nl,percutopo,q)
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,iblock,iblk,icut,imax,jmax,kmax,iq
      integer(kind=cosa_int) percutopo(22,npcuts)
      logical myblockl
      real(kind=cosa_real) q(*)

!     do iblock = 1,mynblockS
        do icut = 1,mynpcuts
          iblk = percutopo( 1,icut)
          call myblock(iblk,myblockl,.true.)
          if (myblockl) then
            imax   = i_imax     (iblk,nl)
            jmax   = j_jmax     (iblk,nl)
            kmax   = k_kmax     (iblk,nl)
            iq     = 1 + off_p3 (iblk,nl) * npde * dim5
            call cut_qtec(imax,jmax,kmax,q(iq),npde,nharms, &
                          percutopo(1,icut))
          end if
        end do
!     end do

      return
      end

!-----------------------------------------------------------------------
      subroutine cut_qtec(imax,jmax,kmax,q,npde,nharms,cutopo)
!-----------------------------------------------------------------------

!     Routine to calculate flow state on bounday cuts as average of
!     interior cell adjacent to cut and boundary cell (this latter has the
!     flow state at the first interior cell of the adjacent block).

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,npde,nharms
      integer(kind=cosa_int) cutopo(21)
      integer(kind=cosa_int) ijkmax(3),istr(3),iend(3),isgn(3),len(3),idir,inout,l, &
                ibcpt,inr,ipde,n,ic1,ic2,ic3,i2,i3,ibc,jbc,kbc,in,jn,kn
      real(kind=cosa_real) q (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms)

!     Store imax, jmax in ijmax for looping

      ijkmax(1) = imax
      ijkmax(2) = jmax
      ijkmax(3) = kmax

!     Store boundary condition data in mnemonic names

      idir     = cutopo( 2)
      inout    = cutopo( 3)
      istr(1)  = cutopo( 4)
      iend(1)  = cutopo( 5)
      istr(2)  = cutopo( 6)
      iend(2)  = cutopo( 7)
      istr(3)  = cutopo( 8)
      iend(3)  = cutopo( 9)

!     Set needed variables depending on whether the boundary is
!     the inner boundary (inout = 1) or the outer boundary (inout > 1)
!          inr    = interior point of block
!          ibcpt  = boundary point of block

      if (inout .eq. 1) then
         ibcpt  = 0
         inr    = 1
      else
         ibcpt  = ijkmax(idir)
         inr    = ijkmax(idir) - 1
      end if
!
!     Find the length of the two outer loops and loop over these using
!     offsets and delta function to set the two cut data points to the
!     two interior data points of block 2.

      do l = 1, 3

         len(l) = abs ( iend(l) - istr(l) )

!        Increment/Decrement 

         if ( iend(l) .gt. istr(l) ) then
            isgn(l) =   1
         else
            isgn(l) = - 1
         end if

      end do

!     ibc boundary point of block
!     in  interior point of block

      ic1 = cyc (idir, 1)
      ic2 = cyc (idir, 2)
      ic3 = cyc (idir, 3)

      do i3 = 0, len(ic3)
        do i2 = 0, len(ic2)

          ibc = ibcpt                      * krd (ic1, 1) + &
                (istr(ic2) + isgn(ic2)*i2) * krd (ic2, 1) + &
                (istr(ic3) + isgn(ic3)*i3) * krd (ic3, 1)
          jbc = ibcpt                      * krd (ic1, 2) + &
                (istr(ic2) + isgn(ic2)*i2) * krd (ic2, 2) + &
                (istr(ic3) + isgn(ic3)*i3) * krd (ic3, 2)
          kbc = ibcpt                      * krd (ic1, 3) + &
                (istr(ic2) + isgn(ic2)*i2) * krd (ic2, 3) + &
                (istr(ic3) + isgn(ic3)*i3) * krd (ic3, 3)

          in  = inr                        * krd (ic1, 1) + &
                (istr(ic2) + isgn(ic2)*i2) * krd (ic2, 1) + &
                (istr(ic3) + isgn(ic3)*i3) * krd (ic3, 1)
          jn  = inr                        * krd (ic1, 2) + &
                (istr(ic2) + isgn(ic2)*i2) * krd (ic2, 2) + &
                (istr(ic3) + isgn(ic3)*i3) * krd (ic3, 2)
          kn  = inr                        * krd (ic1, 3) + &
                (istr(ic2) + isgn(ic2)*i2) * krd (ic2, 3) + &
                (istr(ic3) + isgn(ic3)*i3) * krd (ic3, 3)
!
        do n = 0, 2*nharms
          do ipde = 1, npde
            q(ibc,jbc,kbc,ipde,n) = (q(in ,jn ,kn ,ipde,n) + &
                                     q(ibc,jbc,kbc,ipde,n)) / 2
          end do
!dbg      write(747,*) ibc,jbc,kbc,in,jn,kn
        end do

        end do
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine wr_tec(nl,var1,var2)
!-----------------------------------------------------------------------
!     
!-----------------------------------------------------------------------
       
      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax
      integer(kind=cosa_int) nl,iblock,iv,n,fid(0:2*mharms)
      logical parallel,amcontrol
      real (kind=cosa_real) var1(*), var2(*)
      character*72 line,tag,flowtec1
      double precision :: starttime, endtime, totaltime, maxtime, mintime

      call usingparallel(parallel)
      call amcontroller(amcontrol)

      call getmpitime(starttime)     

      do n = 0,2*nharms

        fid(n)    = 200 + n

!------ open tecplot file(s)

        if (.not.unsteady) then

          flowtec='flow_tec_steady.dat'

        else if (harbal) then

          if (n.le.9) then
            write(tag,'(''00'',i1)') n
          else if (n.le.99) then
            write(tag,'(''0'',i2)') n
          else if (n.le.999) then
            write(tag,'(i3)') n
          end if
          flowtec = 'flow_tec_hb_'//tag
          flowtec = trim(flowtec)//'.dat'

        else if (dualt) then

          if (itime.le.9) then
            write(flowtec1,'(''000'',i1)') itime
          elseif (itime.le.99) then
            write(flowtec1,'(''00'',i2)')  itime
          elseif (itime.le.999) then
            write(flowtec1,'(''0'',i3)')   itime
          else
            write(flowtec1,'(i4)')         itime
          end if
          flowtec1 = 'flow_tec_'//flowtec1
          flowtec = trim(flowtec1)//'.dat'

        else if (rgkuns) then

          flowtec='flow_tec_urk.dat'

        end if

        if(parallel) then

          call writeparallelflowtecheader(fid(n),flowtec,n,nl, &
                harbal,dualt,rgkuns,kom,kom_bsl,kom_sst,unsteady, &
                simtime,zeit,mharms)

        else

          open(fid(n),file=flowtec,status='replace')

          if (harbal) then
            write(line,'(''TITLE="HB meshflow"'')')
          else if (dualt) then
            write(line,'(''TITLE="DTS meshflow"'')')
          else if (rgkuns) then
            write(line,'(''TITLE="URK meshflow"'')')
          else
            write(line,'(''TITLE="STD meshflow"'')')
          end if

          write(fid(n),'(a)') line
          line=''

          if (kom.or.kom_bsl.or.kom_sst) then
            write(line,'(''VARIABLES=xc,yc,zc,rho,u,v,w,p,t,mach,tke,omega,mut,dist'')')
          else
            write(line,'(''VARIABLES=xc,yc,zc,rho,u,v,w,p,t,mach'')')
          end if

          write(fid(n),'(a)') line

        end if

      end do

      do iblock = 1,mynblocks
        imax   = i_imax     (iblock,nl)
        jmax   = j_jmax     (iblock,nl)
        kmax   = k_kmax     (iblock,nl)
        iv     = 1 + off_p1 (iblock,nl) * npde * dim5
        if (parallel) then
          call write_parallel_wr_tec_b(fid,var1(iv),var2(iv),imax,jmax, &
               kmax,npde,nharms,iblock,nl)
        else
          call wr_tec_b(fid,var1(iv),var2(iv),imax,jmax,kmax,npde, &
                        nharms)
        end if
      end do

      do n = 0,2*nharms
         if (unsteady.and.(.not.parallel)) then
            if (dualt) then
              write(fid(n),15) simtime
            else if (harbal) then
              write(fid(n),15) zeit(n)
            end if
         end if
!------close tecplot file(s)
         if(parallel) then
            call closempiiofile(fid(n))
         else
            close(fid(n))
         end if
      end do

      call getmpitime(endtime)

      totaltime = endtime-starttime
      maxtime = totaltime
      mintime =  totaltime

      call realmaxreduce(maxtime,1)
      call realminreduce(mintime,1)
      
      if(amcontrol) then
         write(*,'(A,F8.2,A,2F8.2,A)') 'Writing flow file time: ',totaltime,' (',maxtime,mintime,') seconds'
      end if


 15   format(e21.14)

      return
      end

!-----------------------------------------------------------------------
      subroutine wr_tec_b(fid,var1,var2,imax,jmax,kmax,npde,nharms)
!-----------------------------------------------------------------------
       
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,npde,nharms
      integer(kind=cosa_int) i,imax1,j,jmax1,k,kmax1,ipde,n,nh
      integer fid(0:2*mharms)
      real (kind=cosa_real) &
           var1  ( 0:imax  , 0:jmax  , 0:kmax  ,npde,0:2*nharms), &
           var2  ( 0:imax  , 0:jmax  , 0:kmax  ,npde,0:2*nharms)
      character*300 line1

      imax1 = imax+1
      jmax1 = jmax+1
      kmax1 = kmax+1

      do n = 0,2*nharms
        nh = n*hbmove

!------ write tecplot block

        if (harbal) then
          if (kom.or.kom_bsl.or.kom_sst) then
            write(line1,'(''ZONE T="arturo HB, mode '',i2,''", I='',i4,'', J='',i4,'', &
      K='',i4,'',F=BLOCK, DT=(DOUBLE DOUBLE DOUBLE &
      DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE &
     &DOUBLE DOUBLE DOUBLE)'')') n,imax1,jmax1,kmax1
          else
            write(line1,'(''ZONE T="arturo HB, mode '',i2,''", I='',i4,'', J='',i4,'', &
      K='',i4,'',F=BLOCK, DT=(DOUBLE DOUBLE DOUBLE &
     &DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)'')') n,imax1,jmax1,kmax1
          end if
        else
          if (kom.or.kom_bsl.or.kom_sst) then
            write(line1,'(''ZONE T="arturo",I='',i4,'', J='',i4,'', &
      K='',i4,'',F=BLOCK, DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBL &
     &E DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)'')') imax1,jmax1,kmax1
          else
            write(line1,'(''ZONE T="arturo",I='',i4,'', J='',i4,'', &
      K='',i4,'',F=BLOCK, DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBL &
     &E DOUBLE DOUBLE DOUBLE DOUBLE)'')') imax1,jmax1,kmax1
          end if
        end if

        write(fid(n),'(a)') line1

        if (kom.or.kom_bsl.or.kom_sst) then
          write (fid(n),14) ((((var1(i,j,k,ipde,n),i=0,imax),j=0,jmax), &
                              k=0,kmax),ipde=1,npde)
          write (fid(n),14) ((((var2(i,j,k,ipde,n),i=0,imax),j=0,jmax), &
                              k=0,kmax),ipde=1,npde)
        else
          write (fid(n),10) ((((var1(i,j,k,ipde,n),i=0,imax),j=0,jmax), &
                              k=0,kmax),ipde=1,npde)
          write (fid(n),10) ((((var2(i,j,k,ipde,n),i=0,imax),j=0,jmax), &
                              k=0,kmax),ipde=1,npde)
        end if

      end do

 14   format(3e22.14)
 10   format(3e22.14)

      return
      end

!-----------------------------------------------------------------------
      subroutine wr_cgns(nl,q,mut,x,y,z,xdot,ydot,zdot,dist)
!-----------------------------------------------------------------------
!     
!-----------------------------------------------------------------------

      use cgns
      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax
      integer(kind=cosa_int) nl,iblock,n,fid(0:2*mharms)
      integer(kind=cosa_int) iq,imut,ixyz,idist,zonenum
      integer :: basenum, ierr, rank
      logical parallel,amcontrol
      real (kind=cosa_real) q(*),mut(*),x(*),y(*),z(*),xdot(*),ydot(*),zdot(*), &
           dist(*)
      character*72 line,tag,flowtec1
      double precision :: starttime, endtime, totaltime, maxtime, mintime

      call usingparallel(parallel)
      call amcontroller(amcontrol)

      call getmpitime(starttime)
      
      do n = 0,2*nharms
         fid(n)    = 200 + n

!------open cgns file(s)
         
         if (.not.unsteady) then
            
            flowtec='flow_steady.cgns'
            
         else if (harbal) then
            
            if (n.le.9) then
               write(tag,'(''00'',i1)') n
            else if (n.le.99) then
               write(tag,'(''0'',i2)') n
            else if (n.le.999) then
               write(tag,'(i3)') n
            end if
            flowtec = 'flow_hb_'//tag
            flowtec = trim(flowtec)//'.cgns'
            
         else if (dualt) then
            
            if (itime.le.9) then
               write(flowtec1,'(''000'',i1)') itime
            elseif (itime.le.99) then
               write(flowtec1,'(''00'',i2)')  itime
            elseif (itime.le.999) then
               write(flowtec1,'(''0'',i3)')   itime
            else
               write(flowtec1,'(i4)')         itime
            end if
            flowtec1 = 'flow_'//flowtec1
            flowtec = trim(flowtec1)//'.cgns'
            
         else if (rgkuns) then
            
            flowtec='flow_urk.cgns'
            
         end if
       
         if(parallel) then

            call writeparallelflowteccgnsheader(fid,basenum,nl,n)
            
            call cgp_open_f(flowtec,CG_MODE_MODIFY,fid(n),ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) 'cgp_open_f error wr_cgns'
               call cg_error_print_f()
            end if
            
         else
            
            flowtec = trim(flowtec)
            call cg_open_f(flowtec,CG_MODE_WRITE,fid(n),ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) 'cg_open_f error on file ',flowtec
               call cg_error_print_f()
            end if
            
            call cg_base_write_f(fid(n),'gridbase',3,3,basenum,ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) 'cg_base_write_f error'
               call cg_error_print_f()
            end if
         end if
         
      end do
     
      do iblock = 1,mynblocks
        imax   = i_imax     (iblock,nl)
        jmax   = j_jmax     (iblock,nl)
        kmax   = k_kmax     (iblock,nl)
        iq     = 1 + off_p3 (iblock,nl) * npde * dim5
        imut   = 1 + off_p3 (iblock,nl) *        dim5
        ixyz   = 1 + off_p2 (iblock,nl)        * dim5h
        idist  = 1 + off_p1 (iblock,nl)
        if (parallel) then
          call myblockglobalid(iblock,zonenum)
          call write_parallel_wr_cgns_b(fid,q(iq),mut(imut),x(ixyz),y(ixyz), &
                z(ixyz),xdot(ixyz),ydot(ixyz),zdot(ixyz),dist(idist),imax,jmax, &
                kmax,npde,nharms,basenum,zonenum)
       else
          call wr_cgns_b(fid,q(iq),mut(imut),x(ixyz),y(ixyz), &
                z(ixyz),xdot(ixyz),ydot(ixyz),zdot(ixyz),dist(idist),imax,jmax, &
                kmax,npde,nharms,basenum,iblock)
        end if
      end do

      do n = 0,2*nharms
!------close cgns file(s)
         if(parallel) then
            call closeparallelcgnsfile(fid(n))
         else
            call cg_close_f(fid(n),ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) 'cg_close_f error'
               call cg_error_print_f()
            end if   

         end if
      end do

      call getmpitime(endtime)

      totaltime = endtime-starttime
      maxtime = totaltime
      mintime =  totaltime

      call realmaxreduce(maxtime,1)
      call realminreduce(mintime,1)

      if(amcontrol) then
         write(*,'(A,F8.2,A,2F8.2,A)') 'Writing flow file time: ',totaltime,' (',maxtime,mintime,') seconds'
      end if


      return
      end

!-----------------------------------------------------------------------
      subroutine wr_cgns_b(fid,q,mut,x,y,z,xdot,ydot,zdot,dist, &
           imax,jmax,kmax,npde,nharms,basenum,zonenum)
!-----------------------------------------------------------------------

      use cgns
      use common_variables
      use cosa_precision

      implicit none
      
      integer(kind=cosa_int) imax,jmax,kmax,npde,nharms
      integer(kind=cosa_int) i,imax1,j,jmax1,k,kmax1,n,nh
      integer fid(0:2*mharms)
      integer ierr,zonenum,basenum,flownum,indexnum
      integer xnum,ynum,znum,disnum
      character*10 zonename
      integer(cgsize_t) ::  sizes(3,3)
      real (kind=cosa_real) &
           q   (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           mut (-1:imax+1,-1:jmax+1,-1:kmax+1     ,0:2*nharms), &
           x   ( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           y   ( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           z   ( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           xdot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           ydot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           zdot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           dist( 0:imax  , 0:jmax  , 0:kmax), &
           velocityX(0:imax  , 0:jmax  , 0:kmax), &
           velocityY(0:imax  , 0:jmax  , 0:kmax), &
           velocityZ(0:imax  , 0:jmax  , 0:kmax), &
           temperature(0:imax  , 0:jmax  , 0:kmax), &
           mach(0:imax  , 0:jmax  , 0:kmax)


      imax1 = imax+1
      jmax1 = jmax+1
      kmax1 = kmax+1
      
      do n = 0,2*nharms
         nh = n*hbmove
         
         write(zonename, "(A5,I5)") "block",zonenum
         
         sizes(1,1) = imax1
         sizes(2,1) = jmax1
         sizes(3,1) = kmax1
         
         sizes(1,2) = imax1-1
         sizes(2,2) = jmax1-1
         sizes(3,2) = kmax1-1
         
         sizes(1,3) = 0
         sizes(2,3) = 0
         sizes(3,3) = 0
         
         call cg_zone_write_f(fid(n),basenum,zonename,sizes, &
              Structured,zonenum,ierr)
         if(ierr .ne. CG_OK) then
            write(*,*) 'cg_zone_write_f error'
            call cg_error_print_f()
         end if
         
         if (kom.or.kom_bsl.or.kom_sst) then                     
            
            do k=0,kmax
               do j=0,jmax
                  do i=0,imax
                     if (moving) then
                        velocityX(i,j,k) = q(i,j,k,2,n) - xdot(i,j,k,nh)
                        velocityY(i,j,k) = q(i,j,k,3,n) - ydot(i,j,k,nh)
                        velocityZ(i,j,k) = q(i,j,k,4,n) - zdot(i,j,k,nh)
                     else
                        velocityX(i,j,k) = q(i,j,k,2,n)
                        velocityY(i,j,k) = q(i,j,k,3,n)
                        velocityZ(i,j,k) = q(i,j,k,4,n)
                     end if
                     temperature(i,j,k) = gamma * q(i,j,k,5,n) / &
                          q(i,j,k,1,n)
                     mach(i,j,k) = &
                          sqrt(velocityX(i,j,k)**2+velocityY(i,j,k)**2+ &
                          velocityZ(i,j,k)**2)/ sqrt(temperature(i,j,k))
                  end do
               end do
            end do
            
            
            call cg_coord_write_f(fid(n),basenum,zonenum,RealDouble, &
                 'CoordinateX',x(0:imax,0:jmax,0:kmax,n),xnum,ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) 'X cg_coord_write_f error'
               call cg_error_print_f()
            end if
            call cg_coord_write_f(fid(n),basenum,zonenum,RealDouble, &
                 'CoordinateY',y(0:imax,0:jmax,0:kmax,n),ynum,ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) 'Y cg_coord_write_f error'
               call cg_error_print_f()
            end if
            call cg_coord_write_f(fid(n),basenum,zonenum,RealDouble, &
                 'CoordinateZ',z(0:imax,0:jmax,0:kmax,n),znum,ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) 'Z cg_coord_write_f error'
               call cg_error_print_f()
            end if
            call cg_sol_write_f(fid(n),basenum,zonenum,'solution', &
                 Vertex,flownum,ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) 'cg_sol_write_f error'
               call cg_error_print_f()
            end if
            call cg_field_write_f(fid(n),basenum,zonenum,flownum, &
                 RealDouble,'Density',q(0:imax,0:jmax,0:kmax,1,n),indexnum,ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) 'cg_field_write_f density error'
               call cg_error_print_f()
            end if
            call cg_field_write_f(fid(n),basenum,zonenum,flownum, &
                 RealDouble,'VelocityX',velocityX,indexnum,ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) 'cg_field_write_f velocityX error'
               call cg_error_print_f()
            end if
            call cg_field_write_f(fid(n),basenum,zonenum,flownum, &
                 RealDouble,'VelocityY',velocityY,indexnum,ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) 'cg_field_write_f velocityY error'
               call cg_error_print_f()
            end if
            call cg_field_write_f(fid(n),basenum,zonenum,flownum, &
                 RealDouble,'VelocityZ',velocityZ,indexnum,ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) 'cg_field_write_f velocityZ error'
               call cg_error_print_f()
            end if
            call cg_field_write_f(fid(n),basenum,zonenum,flownum, &
                 RealDouble,'Pressure',q(0:imax,0:jmax,0:kmax,5,n),indexnum,ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) 'cg_field_write_f pressure error'
               call cg_error_print_f()
            end if
            call cg_field_write_f(fid(n),basenum,zonenum,flownum, &
                 RealDouble,'Temperature',temperature,indexnum,ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) 'cg_field_write_f temperature error'
               call cg_error_print_f()
            end if
            call cg_field_write_f(fid(n),basenum,zonenum,flownum, &
                 RealDouble,'Mach',mach,indexnum,ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) 'cg_field_write_f mach error'
               call cg_error_print_f()
            end if
            call cg_field_write_f(fid(n),basenum,zonenum,flownum, &
                 RealDouble,'TurbulentEnergyKinetic',q(0:imax,0:jmax,0:kmax,6,n), &
                 indexnum,ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) 'cg_field_write_f turbulent kinetic error'
               call cg_error_print_f()
            end if
            call cg_field_write_f(fid(n),basenum,zonenum,flownum, &
                 RealDouble,'TurbulentDissipationRate',q(0:imax,0:jmax,0:kmax,7,n),indexnum,ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) 'cg_field_write_f turbulent disrate error'
               call cg_error_print_f()
            end if
            call cg_field_write_f(fid(n),basenum,zonenum,flownum, &
                 RealDouble,'TurbulentDistance',dist,indexnum,ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) 'cg_field_write_f turbulent distance error'
               call cg_error_print_f()
            end if
            call cg_field_write_f(fid(n),basenum,zonenum,flownum, &
                 RealDouble,'ViscosityEddy',mut(0:imax,0:jmax,0:kmax,n),indexnum,ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) 'cg_field_write_f viscosity eddy error'
               call cg_error_print_f()
            end if
!            call cg_goto_f(fid(n),basenum,ierr,'Zone_t',
!     &           zonenum,'end')
!            if(ierr .ne. CG_OK) then
!               write(*,*) 'cg_goto_f error'
!               call cg_error_print_f()
!            end if
!            call cg_state_write_f('ReferenceQuantities',ierr)
!            if(ierr .ne. CG_OK) then
!               write(*,*) 'cg_state_write_f error'
!               call cg_error_print_f()
!            end if
!            call cg_goto_f(fid(n),basenum,ierr,'Zone_t',zonenum,
!     &           'ReferenceState_t',1,'end')
!            if(ierr .ne. CG_OK) then
!               write(*,*) 'cg_goto_f error'
!               call cg_error_print_f()
!            end if
!            call cg_array_write_f('Mach',RealDouble,3,sizes(:,1),mach,ierr)
            
         else

            do k=0,kmax
               do j=0,jmax
                  do i=0,imax
                     if (moving) then
                        velocityX(i,j,k) = q(i,j,k,2,n) - xdot(i,j,k,nh)
                        velocityY(i,j,k) = q(i,j,k,3,n) - ydot(i,j,k,nh)
                        velocityZ(i,j,k) = q(i,j,k,4,n) - zdot(i,j,k,nh)
                     else
                        velocityX(i,j,k) = q(i,j,k,2,n)
                        velocityY(i,j,k) = q(i,j,k,3,n)
                        velocityZ(i,j,k) = q(i,j,k,4,n)
                     end if
                     temperature(i,j,k) = gamma * q(i,j,k,5,n) / &
                          q(i,j,k,1,n)
                     mach(i,j,k) = &
                          sqrt(velocityX(i,j,k)**2+velocityY(i,j,k)**2+ &
                          velocityZ(i,j,k)**2)/ sqrt(temperature(i,j,k))
                  end do
               end do
            end do
            
            
            call cg_coord_write_f(fid(n),basenum,zonenum,RealDouble, &
                 'CoordinateX',x(0:imax,0:jmax,0:kmax,n),xnum,ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) 'X cg_coord_write_f error'
               call cg_error_print_f()
            end if
            call cg_coord_write_f(fid(n),basenum,zonenum,RealDouble, &
                 'CoordinateY',y(0:imax,0:jmax,0:kmax,n),ynum,ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) 'Y cg_coord_write_f error'
               call cg_error_print_f()
            end if
            call cg_coord_write_f(fid(n),basenum,zonenum,RealDouble, &
                 'CoordinateZ',z(0:imax,0:jmax,0:kmax,n),znum,ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) 'Z cg_coord_write_f error'
               call cg_error_print_f()
            end if
            call cg_sol_write_f(fid(n),basenum,zonenum,'solution', &
                 Vertex,flownum,ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) 'cg_sol_write_f error'
               call cg_error_print_f()
            end if
            call cg_field_write_f(fid(n),basenum,zonenum,flownum, &
                 RealDouble,'Density',q(0:imax,0:jmax,0:kmax,4,n),indexnum,ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) 'cg_field_write_f density error'
               call cg_error_print_f()
            end if
            call cg_field_write_f(fid(n),basenum,zonenum,flownum, &
                 RealDouble,'VelocityX',velocityX,indexnum,ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) 'cg_field_write_f velocityX error'
               call cg_error_print_f()
            end if
            call cg_field_write_f(fid(n),basenum,zonenum,flownum, &
                 RealDouble,'VelocityY',velocityY,indexnum,ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) 'cg_field_write_f velocityY error'
               call cg_error_print_f()
            end if
            call cg_field_write_f(fid(n),basenum,zonenum,flownum, &
                 RealDouble,'VelocityZ',velocityZ,indexnum,ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) 'cg_field_write_f velocityZ error'
               call cg_error_print_f()
            end if
            call cg_field_write_f(fid(n),basenum,zonenum,flownum, &
                 RealDouble,'Pressure',q(0:imax,0:jmax,0:kmax,5,n),indexnum,ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) 'cg_field_write_f pressure error'
               call cg_error_print_f()
            end if
            call cg_field_write_f(fid(n),basenum,zonenum,flownum, &
                 RealDouble,'Temperature',temperature,indexnum,ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) 'cg_field_write_f temperature error'
               call cg_error_print_f()
            end if
            call cg_field_write_f(fid(n),basenum,zonenum,flownum, &
                 RealDouble,'Mach',mach,indexnum,ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) 'cg_field_write_f Mach error'
               call cg_error_print_f()
            end if
!            call cg_goto_f(fid(n),basenum,ierr,'Zone_t',
!     &           zonenum,'end')
!            if(ierr .ne. CG_OK) then
!               write(*,*) 'cg_goto_f error'
!               call cg_error_print_f()
!            end if
!            call cg_state_write_f('ReferenceQuantities',ierr)
!            if(ierr .ne. CG_OK) then
!               write(*,*) 'cg_state_write_f error'
!               call cg_error_print_f()
!            end if
!            call cg_goto_f(fid(n),basenum,ierr,'Zone_t',zonenum,
!     &           'ReferenceState_t',1,'end')
!            if(ierr .ne. CG_OK) then
!               write(*,*) 'cg_goto_f error'
!               call cg_error_print_f()
!            end if
!            call cg_array_write_f('Mach',RealDouble,3,sizes(:,1),mach,ierr)

         end if
         
      end do
      
      return
      end 


!-----------------------------------------------------------------------
      subroutine zeitdata_out
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) key1, fid

      fid = 15
      open(fid,file='zeit.dat',form='formatted',status='replace')

      if (dualt) then
        key1 = 9992999
        write(fid,*) key1
        write(fid,*) simtime
      else  if (rgkuns) then
        key1 = 9991999
        write(fid,*) key1
        write(fid,*) simtime
      else
        key1 = 9990999
        write(fid,*) key1
      end if

      close(fid)

      return
      end

!-----------------------------------------------------------------------
      subroutine wr_tec_surf(nl,q,x,y,z,xdot,ydot,zdot,dist,si,sj,sk,xideri, &
        xiderj,xiderk,etaderi,etaderj,etaderk,zetaderi,zetaderj,zetaderk, &
        bctopo)
!-----------------------------------------------------------------------
!     
!-----------------------------------------------------------------------

      use parallelutils, only: surfblockstarts, surfblockends, &
           surfbcindex, surfblockindex
      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,bctopo(*)
      integer(kind=cosa_int) nl,iblock,ibc,iq,iblk,ixyz,iimt,ivmt,idist,n,fid(0:2*mharms)     
      logical parallel,amcontrol
      real (kind=cosa_real) q(*),x(*),y(*),z(*),xdot(*),ydot(*),zdot(*),dist(*), &
        si(*),sj(*),sk(*),xideri(*),xiderj(*),xiderk(*),etaderi(*),etaderj(*), &
        etaderk(*),zetaderi(*),zetaderj(*),zetaderk(*)
      character*72 line,tag
      double precision :: starttime, endtime, totaltime, maxtime, mintime

      call usingparallel(parallel)
      call amcontroller(amcontrol)


      call getmpitime(starttime)     

      do n = 0,2*nharms

         fid(n)    = 200 + n

!------ build tecplot file(s) name

        if (.not.unsteady) then

          surftec='surf_tec_steady.dat'

        else if (harbal) then

          if (n.le.9) then
            write(tag,'(''00'',i1)') n
          else if (n.le.99) then
            write(tag,'(''0'',i2)') n
          else if (n.le.999) then
            write(tag,'(i3)') n
          end if
          surftec = 'surf_tec_hb_'//tag
          surftec = trim(surftec)//'.dat'

        else if (dualt) then

          if (itime.le.9) then
            write(tag,'(''000'',i1)') itime
          elseif (itime.le.99) then
            write(tag,'(''00'',i2)')  itime
          elseif (itime.le.999) then
            write(tag,'(''0'',i3)')   itime
          else
            write(tag,'(i4)')         itime
          end if
          tag = 'surf_tec_t_'//tag
          surftec = trim(tag)//'.dat'

        else if (rgkuns) then

          surftec='surf_tec_urk.dat'

        end if

!------ open tecplot file(s) and write 2-line file header


        if(parallel) then

           call writeparallelsurftecheader(fid(n),n,nl)

        else

          open(fid(n),file=surftec,status='replace')

          if (harbal) then
            write(line,'(''TITLE="HB surfdata"'')')
          else if (dualt) then
            write(line,'(''TITLE="DTS surfdata"'')')
          else if (rgkuns) then
            write(line,'(''TITLE="URK surfdata"'')')
          else
            write(line,'(''TITLE="STD surfdata"'')')
          end if

          write(fid(n),'(a)') line
          line=''

          write(line,'(''VARIABLES=xc,yc,zc,cp,cf,yp'')')
          write(fid(n),'(a)') line

       end if

      end do

      do iblock = 1,mynblocks

         imax   = i_imax     (iblock,nl)
         jmax   = j_jmax     (iblock,nl)
         kmax   = k_kmax     (iblock,nl)
         iq     = 1 + off_p3 (iblock,nl) * npde * dim5
         iblk   = 1 + off_bct(iblock,nl)
         ixyz   = 1 + off_p2 (iblock,nl)        * dim5h
         iimt   = 1 + off_p2 (iblock,nl) * lmet * dim5h
         ivmt   = 1 + off_0  (iblock,nl) * 3    * dim5h
         idist  = 1 + off_p1 (iblock,nl)
         if (parallel) then

            if(1.lt.0) then

               call write_parallel_wr_tec_surf_b_M1(fid,imax,jmax,kmax,npde, &
                    nharms,iblock,nl,lmet,q(iq),x(ixyz),y(ixyz),z(ixyz), &
                    xdot(ixyz),ydot(ixyz),zdot(ixyz),si(iimt),sj(iimt),sk(iimt), &
                    xideri(ivmt),etaderj(ivmt),zetaderk(ivmt),dist(idist), &
                    bctopo(iblk),nbcs(iblock))

            else

               call write_parallel_wr_tec_surf_b_M2(fid,imax,jmax,kmax,npde, &
                    nharms,iblock,nl,lmet,q(iq),x(ixyz),y(ixyz),z(ixyz), &
                    xdot(ixyz),ydot(ixyz),zdot(ixyz),si(iimt),sj(iimt),sk(iimt), &
                    xideri(ivmt),xiderj(ivmt),xiderk(ivmt),etaderi(ivmt), &
                    etaderj(ivmt),etaderk(ivmt),zetaderi(ivmt), &
                    zetaderj(ivmt),zetaderk(ivmt),dist(idist), &
                    bctopo(iblk),nbcs(iblock))
               
            end if

         else

            if (1.lt.0) then
               call wr_tec_surf_b_M1(fid,q(iq),x(ixyz),y(ixyz),z(ixyz), &
                    xdot(ixyz),ydot(ixyz),zdot(ixyz),si(iimt),sj(iimt),sk(iimt), &
                    xideri(ivmt),etaderj(ivmt),zetaderk(ivmt), &
                    dist(idist),bctopo(iblk),imax,jmax,kmax, &
                    lmet,nbcs(iblock),npde,nharms)
            else
               call wr_tec_surf_b_M2(fid,q(iq),x(ixyz),y(ixyz),z(ixyz), &
                    xdot(ixyz),ydot(ixyz),zdot(ixyz),si(iimt),sj(iimt),sk(iimt), &
                    xideri(ivmt),xiderj(ivmt),xiderk(ivmt),etaderi(ivmt), &
                    etaderj(ivmt),etaderk(ivmt),zetaderi(ivmt),zetaderj(ivmt), &
                    zetaderk(ivmt),dist(idist),bctopo(iblk), &
                    imax,jmax,kmax,lmet,nbcs(iblock),npde,nharms)
            end if 
         end if
      end do
      
      do n = 0,2*nharms
!------close tecplot file(s)
         if(parallel) then
            if(allocated(surfblockstarts)) deallocate(surfblockstarts)
            if(allocated(surfblockends)) deallocate(surfblockends)
            if(allocated(surfblockindex)) deallocate(surfblockindex)
            if(allocated(surfbcindex)) deallocate(surfbcindex)
            call closempiiofile(fid(n))
         else
            close(fid(n))
         end if
      end do

      call getmpitime(endtime)

      totaltime = endtime-starttime
      maxtime = totaltime
      mintime =  totaltime

      call realmaxreduce(maxtime,1)
      call realminreduce(mintime,1)

      if(amcontrol) then        
         write(*,'(A,F8.2,A,2F8.2,A)') 'Writing surface file time: ',totaltime,' (',maxtime,mintime,') seconds'
      end if


      return
      end

!-----------------------------------------------------------------------
      subroutine wr_tec_surf_b_M1(fid,q,x,y,z,xdot,ydot,zdot,si,sj,sk,xideri, &
        etaderj,zetaderk,dist,bctopo, &
        imax,jmax,kmax,lmet,nbcs,npde,nharms)
!-----------------------------------------------------------------------
       
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,npde,nharms,lmet,nbcs
      integer(kind=cosa_int) i,imax1,j,jmax1,k,kmax1,ipde,n,nh,ibc_s,jbc_s,kbc_s, &
        ibc_e,jbc_e,kbc_e
      integer(kind=cosa_int) ibc1,isurface
      integer(kind=cosa_int) bctopo(10,nbcs)
      integer(kind=cosa_int) i1,i2,i3,ijkmax(3),idir,istrt(3),iend(3),bctyp, &
        inrout,ibcpt,ibcpt2,ibcn,ibcn2,ibcm,ibc,jbc,kbc,ibc2,jbc2, &
        kbc2,in,jn,kn,in2,jn2,kn2,im,jm,km,ic1,ic2,ic3,sysize,ioff1, &
        ioff2,joff1,joff2,koff1,koff2

      integer fid(0:2*mharms)
      real (kind=cosa_real) sgnm,nx,ny,nz,kx,ky,kz,rhow,tw,muw,uw,u1,u2,uwr,u1r,u2r, &
           vw,v1,v2,vwr,v1r,v2r,ww,w1,w2,wwr,w1r,w2r,ueta,veta,weta, &
           etax,etay,etaz,dudx,dudy,dudz,dvdx,dvdy, &
           dvdz,dwdx,dwdy,dwdz,divv,txx,txy,txz,tyy,tyz,tzz,mach,tauwx, &
           tauwy,tauwz,tauwpx,tauwpy,tauwpz,dn1,tauw,utau
      real (kind=cosa_real) &
           q   (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           x   ( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           y   ( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           z   ( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           xdot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           ydot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           zdot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           si      (lmet, 0:imax+1, 0:jmax+1, 0:kmax+1,0:2*nharms*hbmove), &
           sj      (lmet, 0:imax+1, 0:jmax+1, 0:kmax+1,0:2*nharms*hbmove), &
           sk      (lmet, 0:imax+1, 0:jmax+1, 0:kmax+1,0:2*nharms*hbmove), &
           xideri  (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           etaderj (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           zetaderk(3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           dist( 0:imax  , 0:jmax  , 0:kmax)
      real (kind=cosa_real),allocatable::cp(:,:,:,:),cf(:,:,:,:),yp(:,:,:,:)
      character*200 line1

      ijkmax(1) = imax
      ijkmax(2) = jmax
      ijkmax(3) = kmax

      imax1 = imax+1
      jmax1 = jmax+1
      kmax1 = kmax+1

!---- loop over harmonics - STARTS
      do n = 0,2*nharms
         nh = n*hbmove
         
!------loop over surfaces - START
!------MSC 12/03/23. In future nsurface loop can be removed HERE if we do
!     not wish to distinguish separate surfaces/walls
!     in surftec file. We will just check for BCs
!     starting with 15 or 14.
!     
!     Keeping nsurface is instead essential to get
!     forces and moments in routine force.
!     Getting nsurface from input.f is trivial when
!     all walls are either visc. (RANS) or inviscid
!     (Euler). We just count them. 
!     NOTE: we will need to use for BC 14 the same syntax as 
!     for BC 1500s (14xx).
!     
!     Less trivial is case when
!     a body has both visc. and inv. BCs (e.g. hub).
!     For this, previous code structure with list of
!     blocks with boundaries on given body was better.
         
         do isurface=1,nsurface
            
            do ibc1=1,nbcs
               if ((bctopo(1,ibc1).eq.1500+(isurface-1)).or. &
                    (bctopo(1,ibc1).eq.14)) then
                  
!     store boundary topology in mnemonic names
                  bctyp    = bctopo(1,ibc1)
                  idir     = bctopo(2,ibc1)
                  inrout   = bctopo(3,ibc1)
                  istrt(1) = bctopo(4,ibc1)
                  iend(1)  = bctopo(5,ibc1)
                  istrt(2) = bctopo(6,ibc1)
                  iend(2)  = bctopo(7,ibc1)
                  istrt(3) = bctopo(8,ibc1)
                  iend(3)  = bctopo(9,ibc1)
                  
!     set needed variables depending on whether the boundary is
!     the inner boundary (inrout = 1) or
!     the outer boundary (inrout > 1)
!     ibcpt : boundary condition location (first aux. cell)
!     ibcpt2: boundary condition location outside the block from
!     ibcpt (second aux. cell)
!     ibcn  : point to the inside of the block from ibcpt
!     ibcm  : location of the metrics
                  
                  if (inrout.eq.1) then
                     ibcpt  =  0
                     ibcpt2 = -1
                     ibcn   =  1
                     ibcn2  =  2
                     ibcm   =  1
                     sgnm   = - 1.d0
                  else
                     ibcpt  = ijkmax(idir)
                     ibcpt2 = ijkmax(idir) + 1
                     ibcn   = ijkmax(idir) - 1
                     ibcn2  = ijkmax(idir) - 2
                     ibcm   = ijkmax(idir)
                     sgnm   = 1.d0
                  end if
                  
                  if (idir.eq.1) then
                     ioff1 = (-1)**(1+1/inrout)
                     ioff2 = 2*ioff1
                     joff1 = 0
                     joff2 = 0
                     koff1 = 0
                     koff2 = 0
                  else if (idir.eq.2) then
                     ioff1 = 0
                     ioff2 = 0
                     joff1 = (-1)**(1+1/inrout)
                     joff2 = 2*joff1
                     koff1 = 0
                     koff2 = 0
                  else if (idir.eq.3) then
                     ioff1 = 0
                     ioff2 = 0
                     joff1 = 0
                     joff2 = 0
                     koff1 = (-1)**(1+1/inrout)
                     koff2 = 2*koff1
                  end if
                  
                  ic1 = cyc (idir, 1)
                  ic2 = cyc (idir, 2)
                  ic3 = cyc (idir, 3)
                  
                  ibc_s  = ibcpt *krd(ic1,1) + istrt(ic2)*krd(ic2,1) + &
                       istrt(ic3)*krd(ic3,1)
                  jbc_s  = ibcpt *krd(ic1,2) + istrt(ic2)*krd(ic2,2) + &
                       istrt(ic3)*krd(ic3,2)
                  kbc_s  = ibcpt *krd(ic1,3) + istrt(ic2)*krd(ic2,3) + &
                       istrt(ic3)*krd(ic3,3)
                  
                  ibc_e  = ibcpt *krd(ic1,1) + iend(ic2)*krd(ic2,1) + &
                       iend(ic3)*krd(ic3,1)
                  jbc_e  = ibcpt *krd(ic1,2) + iend(ic2)*krd(ic2,2) + &
                       iend(ic3)*krd(ic3,2)
                  kbc_e  = ibcpt *krd(ic1,3) + iend(ic2)*krd(ic2,3) + &
                       iend(ic3)*krd(ic3,3)
                  
                  allocate(cp(ibc_s:ibc_e,jbc_s:jbc_e,kbc_s:kbc_e,0:nharms))
                  if (viscous) then
                     allocate(cf(ibc_s:ibc_e,jbc_s:jbc_e,kbc_s:kbc_e,0:nharms), &
                          yp(ibc_s:ibc_e,jbc_s:jbc_e,kbc_s:kbc_e,0:nharms))
                  end if
                  
!------------extract/calculate surface coefficients - START
                  do i3 = istrt(ic3),iend(ic3)
                     do i2 = istrt(ic2),iend(ic2)
                        
                        ibc  = ibcpt *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
                        jbc  = ibcpt *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
                        kbc  = ibcpt *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)
                        
                        im   = ibcm  *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
                        jm   = ibcm  *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
                        km   = ibcm  *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)
                        
                        if (idir.eq.1) then
                           nx    = si(1,im,jm,km,nh)
                           ny    = si(2,im,jm,km,nh)
                           nz    = si(3,im,jm,km,nh)
                        else if (idir.eq.2) then
                           nx    = sj(1,im,jm,km,nh)
                           ny    = sj(2,im,jm,km,nh)
                           nz    = sj(3,im,jm,km,nh)
                        else if (idir.eq.3) then
                           nx    = sk(1,im,jm,km,nh)
                           ny    = sk(2,im,jm,km,nh)
                           nz    = sk(3,im,jm,km,nh)
                        end if
                        
                        kx   = (-1)**(1+1/inrout) * nx
                        ky   = (-1)**(1+1/inrout) * ny
                        kz   = (-1)**(1+1/inrout) * nz
                        
                        cp(ibc,jbc,kbc,n) = 2*(q(ibc,jbc,kbc,5,n)-1/gamma)/ &
                             machfs**2
                        
                        if (viscous) then
                           
                           ibc2 = ibcpt2*krd(ic1,1) +i2*krd(ic2,1) +i3*krd(ic3,1)
                           jbc2 = ibcpt2*krd(ic1,2) +i2*krd(ic2,2) +i3*krd(ic3,2)
                           kbc2 = ibcpt2*krd(Ic1,3) +i2*krd(ic2,3) +i3*krd(ic3,3)
                           
                           in   = ibcn  *krd(ic1,1) +i2*krd(ic2,1) +i3*krd(ic3,1)
                           jn   = ibcn  *krd(ic1,2) +i2*krd(ic2,2) +i3*krd(ic3,2)
                           kn   = ibcn  *krd(ic1,3) +i2*krd(ic2,3) +i3*krd(ic3,3)
                           
                           in2  = ibcn2 *krd(ic1,1) +i2*krd(ic2,1) +i3*krd(ic3,1)
                           jn2  = ibcn2 *krd(ic1,2) +i2*krd(ic2,2) +i3*krd(ic3,2)
                           kn2  = ibcn2 *krd(ic1,3) +i2*krd(ic2,3) +i3*krd(ic3,3)
                           
                           if (idir.eq.1) then
                              etax  = xideri(1,im,jm,km,nh)
                              etay  = xideri(2,im,jm,km,nh)
                              etaz  = xideri(3,im,jm,km,nh)
                           else if (idir.eq.2) then
                              etax  = etaderj(1,im,jm,km,nh)
                              etay  = etaderj(2,im,jm,km,nh)
                              etaz  = etaderj(3,im,jm,km,nh)
                           else if (idir.eq.3) then
                              etax  = zetaderk(1,im,jm,km,nh)
                              etay  = zetaderk(2,im,jm,km,nh)
                              etaz  = zetaderk(3,im,jm,km,nh)
                           end if
                           
                           rhow = q(ibc,jbc,kbc,1,n)
                           tw   = q(ibc,jbc,kbc,5,n) *gamma / q(ibc,jbc,kbc,1,n)
                           muw  = (stemp+1)/(tw+stemp) * (tw**1.5d0)
                           
!------------------MSC, 11/03/2023: dist not computed for laminar now
                           if (kom.or.kom_bsl.or.kom_sst) then
                              dn1  = dist(in ,jn ,kn )
                           else
                              dn1  = 0
                           end if
                           
                           uw   = q(ibc,jbc,kbc,2,n)
                           u1   = q(ibc+ioff1,jbc+joff1,kbc+koff1,2,n)
                           u2   = q(ibc+ioff2,jbc+joff2,kbc+koff2,2,n)
                           if (moving) then
                              uwr  = xdot(ibc      ,jbc      ,kbc      ,nh)
                              u1r  = xdot(ibc+ioff1,jbc+joff1,kbc+koff1,nh)
                              u2r  = xdot(ibc+ioff2,jbc+joff2,kbc+koff2,nh)
                              uw   = uw - uwr
                              u1   = u1 - u1r
                              u2   = u2 - u2r
                           end if
                           
                           vw   = q(ibc,jbc,kbc,3,n)
                           v1   = q(ibc+ioff1,jbc+joff1,kbc+koff1,3,n)
                           v2   = q(ibc+ioff2,jbc+joff2,kbc+koff2,3,n)
                           if (moving) then
                              vwr  = ydot(ibc      ,jbc      ,kbc      ,nh)
                              v1r  = ydot(ibc+ioff1,jbc+joff1,kbc+koff1,nh)
                              v2r  = ydot(ibc+ioff2,jbc+joff2,kbc+koff2,nh)
                              vw   = vw - vwr
                              v1   = v1 - v1r
                              v2   = v2 - v2r
                           end if
                           
                           ww   = q(ibc,jbc,kbc,4,n)
                           w1   = q(ibc+ioff1,jbc+joff1,kbc+koff1,4,n)
                           w2   = q(ibc+ioff2,jbc+joff2,kbc+koff2,4,n)
                           if (moving) then
                              wwr  = zdot(ibc      ,jbc      ,kbc      ,nh)
                              w1r  = zdot(ibc+ioff1,jbc+joff1,kbc+koff1,nh)
                              w2r  = zdot(ibc+ioff2,jbc+joff2,kbc+koff2,nh)
                              ww   = ww - wwr
                              w1   = w1 - w1r
                              w2   = w2 - w2r
                           end if
                           
                           ueta = (-1)**(1+1/inrout) * (-8*uw + 9*u1 - u2) / 3
                           veta = (-1)**(1+1/inrout) * (-8*vw + 9*v1 - v2) / 3
                           weta = (-1)**(1+1/inrout) * (-8*ww + 9*w1 - w2) / 3
                           
                           dudx = ueta*etax
                           dudy = ueta*etay
                           dudz = ueta*etaz
                           dvdx = veta*etax
                           dvdy = veta*etay
                           dvdz = veta*etaz
                           dwdx = weta*etax
                           dwdy = weta*etay
                           dwdz = weta*etaz
                           divv = dudx + dvdy + dwdz
                           
                           txx = 2*muw * (dudx - divv/3)
                           txy =   muw * (dudy + dvdx)
                           txz =   muw * (dudz + dwdx)
                           tyy = 2*muw * (dvdy - divv/3)
                           tyz =   muw * (dvdz + dwdy)
                           tzz = 2*muw * (dwdz - divv/3)
                           
                           tauwx = txx*kx + txy*ky + txz*kz
                           tauwy = txy*kx + tyy*ky + tyz*kz
                           tauwz = txz*kx + tyz*ky + tzz*kz
                           
                           tauwpx = tauwx - (tauwx*kx+tauwy*ky+tauwz*kz)*kx
                           tauwpy = tauwy - (tauwx*kx+tauwy*ky+tauwz*kz)*ky
                           tauwpz = tauwz - (tauwx*kx+tauwy*ky+tauwz*kz)*kz
                           
                           tauw   = dsqrt(tauwpx**2+tauwpy**2+tauwpz**2)
                           utau   = dsqrt(tauw/rhow)
                           
                           cf(ibc,jbc,kbc,n) = 2*tauw / (reyno*machfs)
                           yp(ibc,jbc,kbc,n) = dsqrt(reyno/machfs)* &
                                rhow*dn1*utau/muw
                           
                        end if
                        
                     end do
                  end do
!------------extract/calculate surface coefficients - END
                  
!------------write tecplot surface patch
                  
                  if (harbal) then
                     write(line1,'(''ZONE T="arturo HB, mode '',i2,''", I='', &
      i4,'', J='',i4,'', K='',i4,'',F=BLOCK, DT=(DOUBLE DOUBLE DOUBLE DO &
     &UBLE DOUBLE DOUBLE)'')') n,ibc_e-ibc_s+1,jbc_e-jbc_s+1,kbc_e-kbc_s+1
                  else
                     write(line1,'(''ZONE T="arturo",I='',i4,'', J='',i4,'', &
      K='',i4,'',F=BLOCK, DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE) &
     &'')') ibc_e-ibc_s+1,jbc_e-jbc_s+1,kbc_e-kbc_s+1
                  end if
                  
                  write(fid(n),'(a)') line1
                  
                  write (fid(n),10) (((x(i,j,k,n),i=ibc_s,ibc_e),j=jbc_s,jbc_e), &
                       k=kbc_s,kbc_e)
                  write (fid(n),10) (((y(i,j,k,n),i=ibc_s,ibc_e),j=jbc_s,jbc_e), &
                       k=kbc_s,kbc_e)
                  write (fid(n),10) (((z(i,j,k,n),i=ibc_s,ibc_e),j=jbc_s,jbc_e), &
                       k=kbc_s,kbc_e)
                  write (fid(n),10) (((cp(i,j,k,n),i=ibc_s,ibc_e),j=jbc_s,jbc_e), &
                       k=kbc_s,kbc_e)
                  if (viscous) then
                     write (fid(n),10) (((cf(i,j,k,n),i=ibc_s,ibc_e),j=jbc_s,jbc_e), &
                          k=kbc_s,kbc_e)
                     write (fid(n),10) (((yp(i,j,k,n),i=ibc_s,ibc_e),j=jbc_s,jbc_e), &
                          k=kbc_s,kbc_e)
                  else
                     write (fid(n),10) (((0.d0       ,i=ibc_s,ibc_e),j=jbc_s,jbc_e), &
                          k=kbc_s,kbc_e)
                     write (fid(n),10) (((0.d0       ,i=ibc_s,ibc_e),j=jbc_s,jbc_e), &
                          k=kbc_s,kbc_e)
                  end if
                  
                  deallocate(cp)
                  if (viscous) then
                     deallocate(cf,yp)
                  end if
                  
               end if
               
!--------end loop (ibc1=1,nbcs)
            end do
            
!------loop over surfaces - END
         end do
         
!---- loop over harmonics - END
      end do
      

 14   format(3e22.14)
 10   format(3e22.14)

      return
      end

!-----------------------------------------------------------------------
      subroutine wr_tec_surf_b_M2(fid,q,x,y,z,xdot,ydot,zdot,si,sj,sk,xideri, &
        xiderj,xiderk,etaderi,etaderj,etaderk,zetaderi,zetaderj,zetaderk,dist, &
        bctopo,imax,jmax,kmax,lmet,nbcs,npde,nharms)
!-----------------------------------------------------------------------
       
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,npde,nharms,lmet,nbcs
      integer(kind=cosa_int) i,imax1,j,jmax1,k,kmax1,ipde,n,nh,ibc_s,jbc_s,kbc_s, &
        ibc_e,jbc_e,kbc_e
      integer(kind=cosa_int) ibc1,isurface
      integer(kind=cosa_int) bctopo(10,nbcs)
      integer(kind=cosa_int) i1,i2,i3,ijkmax(3),idir,istrt(3),iend(3),bctyp, &
        inrout,ibcpt,ibcpt2,ibcn,ibcn2,ibcm,ibc,jbc,kbc,ibc2,jbc2, &
        kbc2,in,jn,kn,in2,jn2,kn2,im,jm,km,ic1,ic2,ic3,sysize,ioff1, &
        ioff2,joff1,joff2,koff1,koff2

      integer fid(0:2*mharms)
      real (kind=cosa_real) sgnm,nx,ny,nz,kx,ky,kz,rhow,tw,muw,uw,u1,u2,uwr,u1r,u2r, &
        vw,v1,v2,vwr,v1r,v2r,ww,w1,w2,wwr,w1r,w2r, &
        xix,xiy,xiz,etax,etay,etaz,zetax,zetay,zetaz, &
        dudxi,dudeta,dudzeta,dvdxi,dvdeta,dvdzeta,dwdxi,dwdeta,dwdzeta, &
        dudx,dudy,dudz,dvdx,dvdy, &
        dvdz,dwdx,dwdy,dwdz,divv,txx,txy,txz,tyy,tyz,tzz,mach,tauwx, &
        tauwy,tauwz,tauwpx,tauwpy,tauwpz,dn1,tauw,utau
      real (kind=cosa_real) &
           q   (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           x   ( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           y   ( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           z   ( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           xdot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           ydot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           zdot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           si      (lmet, 0:imax+1, 0:jmax+1, 0:kmax+1,0:2*nharms*hbmove), &
           sj      (lmet, 0:imax+1, 0:jmax+1, 0:kmax+1,0:2*nharms*hbmove), &
           sk      (lmet, 0:imax+1, 0:jmax+1, 0:kmax+1,0:2*nharms*hbmove), &
           xideri  (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           xiderj  (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           xiderk  (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           etaderi (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           etaderj (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           etaderk (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           zetaderi(3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           zetaderj(3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           zetaderk(3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           dist( 0:imax  , 0:jmax  , 0:kmax)
      real (kind=cosa_real),allocatable::cp(:,:,:,:),cf(:,:,:,:),yp(:,:,:,:)
      character*200 line1
      
      ijkmax(1) = imax
      ijkmax(2) = jmax
      ijkmax(3) = kmax
      
      imax1 = imax+1
      jmax1 = jmax+1
      kmax1 = kmax+1
      
!---- loop over harmonics - STARTS
      do n = 0,2*nharms
         nh = n*hbmove
         
!------loop over surfaces - START
!------MSC 12/03/23. In future nbofy loop can be removed HERE if we do
!     not wish to distinguish separate surfaces/walls
!     in surftec file. We will just check for BCs
!     starting with 15 or 14.
!     
!     Keeping nbpdy is instead essential to get
!     forces and moments in routine force.
!     Getting nsurface from input.f is trivial when
!     all walls are either visc. (RANS) or inviscid
!     (Euler). We just count them. 
!     NOTE: we will need to use for BC 14 the same syntax as 
!     for BC 1500s (14xx).
!     
!     Less trivial is case when
!     a body has both visc. and inv. BCs (e.g. hub).
!     For this, previous code structure with list of
!     blocks with boundaries on given body was better.
         
         do isurface=1,nsurface
            
            do ibc1=1,nbcs
               if ((bctopo(1,ibc1).eq.1500+(isurface-1)).or. &
                    (bctopo(1,ibc1).eq.14)) then
                  
!     store boundary topology in mnemonic names
                  bctyp    = bctopo(1,ibc1)
                  idir     = bctopo(2,ibc1)
                  inrout   = bctopo(3,ibc1)
                  istrt(1) = bctopo(4,ibc1)
                  iend(1)  = bctopo(5,ibc1)
                  istrt(2) = bctopo(6,ibc1)
                  iend(2)  = bctopo(7,ibc1)
                  istrt(3) = bctopo(8,ibc1)
                  iend(3)  = bctopo(9,ibc1)
                  
!     set needed variables depending on whether the boundary is
!     the inner boundary (inrout = 1) or
!     the outer boundary (inrout > 1)
!     ibcpt : boundary condition location (first aux. cell)
!     ibcpt2: boundary condition location outside the block from
!     ibcpt (second aux. cell)
!     ibcn  : point to the inside of the block from ibcpt
!     ibcm  : location of the metrics
                  
                  if (inrout.eq.1) then
                     ibcpt  =  0
!     del            ibcpt2 = -1
                     ibcn   =  1
                     ibcn2  =  2
                     ibcm   =  1
!     del            sgnm   = - 1.d0
                  else
                     ibcpt  = ijkmax(idir)
!     del            ibcpt2 = ijkmax(idir) + 1
                     ibcn   = ijkmax(idir) - 1
                     ibcn2  = ijkmax(idir) - 2
                     ibcm   = ijkmax(idir)
!     del            sgnm   = 1.d0
                  end if
                  
                  if (idir.eq.1) then
                     ioff1 = (-1)**(1+1/inrout)
                     ioff2 = 2*ioff1
                     joff1 = 0
                     joff2 = 0
                     koff1 = 0
                     koff2 = 0
                  else if (idir.eq.2) then
                     ioff1 = 0
                     ioff2 = 0
                     joff1 = (-1)**(1+1/inrout)
                     joff2 = 2*joff1
                     koff1 = 0
                     koff2 = 0
                  else if (idir.eq.3) then
                     ioff1 = 0
                     ioff2 = 0
                     joff1 = 0
                     joff2 = 0
                     koff1 = (-1)**(1+1/inrout)
                     koff2 = 2*koff1
                  end if
                  
                  ic1 = cyc (idir, 1)
                  ic2 = cyc (idir, 2)
                  ic3 = cyc (idir, 3)

                  ibc_s  = ibcpt *krd(ic1,1) + istrt(ic2)*krd(ic2,1) + &
                       istrt(ic3)*krd(ic3,1)
                  jbc_s  = ibcpt *krd(ic1,2) + istrt(ic2)*krd(ic2,2) + &
                       istrt(ic3)*krd(ic3,2)
                  kbc_s  = ibcpt *krd(ic1,3) + istrt(ic2)*krd(ic2,3) + &
                       istrt(ic3)*krd(ic3,3)
                  
                  ibc_e  = ibcpt *krd(ic1,1) + iend(ic2)*krd(ic2,1) + &
                       iend(ic3)*krd(ic3,1)
                  jbc_e  = ibcpt *krd(ic1,2) + iend(ic2)*krd(ic2,2) + &
                       iend(ic3)*krd(ic3,2)
                  kbc_e  = ibcpt *krd(ic1,3) + iend(ic2)*krd(ic2,3) + &
                       iend(ic3)*krd(ic3,3)
                  
                  allocate(cp(ibc_s:ibc_e,jbc_s:jbc_e,kbc_s:kbc_e,0:nharms))
                  if (viscous) then
                     allocate(cf(ibc_s:ibc_e,jbc_s:jbc_e,kbc_s:kbc_e,0:nharms), &
                          yp(ibc_s:ibc_e,jbc_s:jbc_e,kbc_s:kbc_e,0:nharms))
                  end if
                  
!------------extract/calculate surface coefficients - START
                  do i3 = istrt(ic3),iend(ic3)
                     do i2 = istrt(ic2),iend(ic2)
                        
                        ibc  = ibcpt *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
                        jbc  = ibcpt *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
                        kbc  = ibcpt *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)
                        
                        im   = ibcm  *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
                        jm   = ibcm  *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
                        km   = ibcm  *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)
                        
                        if (idir.eq.1) then
                           nx    = si(1,im,jm,km,nh)
                           ny    = si(2,im,jm,km,nh)
                           nz    = si(3,im,jm,km,nh)
                        else if (idir.eq.2) then
                           nx    = sj(1,im,jm,km,nh)
                           ny    = sj(2,im,jm,km,nh)
                           nz    = sj(3,im,jm,km,nh)
                        else if (idir.eq.3) then
                           nx    = sk(1,im,jm,km,nh)
                           ny    = sk(2,im,jm,km,nh)
                           nz    = sk(3,im,jm,km,nh)
                        end if
                        
                        kx   = (-1)**(1+1/inrout) * nx
                        ky   = (-1)**(1+1/inrout) * ny
                        kz   = (-1)**(1+1/inrout) * nz
                        
                        cp(ibc,jbc,kbc,n) = 2*(q(ibc,jbc,kbc,5,n)-1/gamma)/ &
                             machfs**2
                        
                        if (viscous) then
                           
!     del                ibc2 = ibcpt2*krd(ic1,1) +i2*krd(ic2,1) +i3*krd(ic3,1)
!     del                jbc2 = ibcpt2*krd(ic1,2) +i2*krd(ic2,2) +i3*krd(ic3,2)
!     del                kbc2 = ibcpt2*krd(Ic1,3) +i2*krd(ic2,3) +i3*krd(ic3,3)
                           
                           in   = ibcn  *krd(ic1,1) +i2*krd(ic2,1) +i3*krd(ic3,1)
                           jn   = ibcn  *krd(ic1,2) +i2*krd(ic2,2) +i3*krd(ic3,2)
                           kn   = ibcn  *krd(ic1,3) +i2*krd(ic2,3) +i3*krd(ic3,3)
                           
                           in2  = ibcn2 *krd(ic1,1) +i2*krd(ic2,1) +i3*krd(ic3,1)
                           jn2  = ibcn2 *krd(ic1,2) +i2*krd(ic2,2) +i3*krd(ic3,2)
                           kn2  = ibcn2 *krd(ic1,3) +i2*krd(ic2,3) +i3*krd(ic3,3)
                           
                           rhow = q(ibc,jbc,kbc,1,n)
                           tw   = q(ibc,jbc,kbc,5,n) *gamma / q(ibc,jbc,kbc,1,n)
                           muw  = (stemp+1)/(tw+stemp) * (tw**1.5d0)
                           
!------------------MSC, 11/03/2023: dist not computed for laminar now
                           if (kom.or.kom_bsl.or.kom_sst) then
                              dn1  = dist(in ,jn ,kn )
                           else
                              dn1  = 0
                           end if
                           
                           
                           if (idir.eq.1) then
                              
                              dudxi = (-1)**(1+1/inrout) * &
                                   (-8*q(ibc,jbc,kbc,2,n) + &
                                   9*q(ibc+ioff1,jbc,kbc,2,n) - &
                                   q(ibc+ioff2,jbc,kbc,2,n) ) / 3
                              dvdxi = (-1)**(1+1/inrout) * &
                                   (-8*q(ibc,jbc,kbc,3,n) + &
                                   9*q(ibc+ioff1,jbc,kbc,3,n) - &
                                   q(ibc+ioff2,jbc,kbc,3,n) ) / 3
                              dwdxi = (-1)**(1+1/inrout) * &
                                   (-8*q(ibc,jbc,kbc,4,n) + &
                                   9*q(ibc+ioff1,jbc,kbc,4,n) - &
                                   q(ibc+ioff2,jbc,kbc,4,n) ) / 3
                              
                              if (i2.eq.istrt(ic2)) then
                                 u2 = (  q(ibc,jbc,kbc,2,n) + q(ibc,jbc+1,kbc,2,n)) / 2
                                 u1 = (3*q(ibc,jbc,kbc,2,n) - q(ibc,jbc+1,kbc,2,n)) / 2
                                 v2 = (  q(ibc,jbc,kbc,3,n) + q(ibc,jbc+1,kbc,3,n)) / 2
                                 v1 = (3*q(ibc,jbc,kbc,3,n) - q(ibc,jbc+1,kbc,3,n)) / 2
                                 w2 = (  q(ibc,jbc,kbc,4,n) + q(ibc,jbc+1,kbc,4,n)) / 2
                                 w1 = (3*q(ibc,jbc,kbc,4,n) - q(ibc,jbc+1,kbc,4,n)) / 2
                              else if (i2.eq.iend(ic2)) then
                                 u2 = (3*q(ibc,jbc,kbc,2,n) - q(ibc,jbc-1,kbc,2,n)) / 2
                                 u1 = (  q(ibc,jbc,kbc,2,n) + q(ibc,jbc-1,kbc,2,n)) / 2
                                 v2 = (3*q(ibc,jbc,kbc,3,n) - q(ibc,jbc-1,kbc,3,n)) / 2
                                 v1 = (  q(ibc,jbc,kbc,3,n) + q(ibc,jbc-1,kbc,3,n)) / 2
                                 w2 = (3*q(ibc,jbc,kbc,4,n) - q(ibc,jbc-1,kbc,4,n)) / 2
                                 w1 = (  q(ibc,jbc,kbc,4,n) + q(ibc,jbc-1,kbc,4,n)) / 2
                              else
                                 u2 = (q(ibc,jbc  ,kbc,2,n) + q(ibc,jbc+1,kbc,2,n)) / 2
                                 u1 = (q(ibc,jbc-1,kbc,2,n) + q(ibc,jbc  ,kbc,2,n)) / 2
                                 v2 = (q(ibc,jbc  ,kbc,3,n) + q(ibc,jbc+1,kbc,3,n)) / 2
                                 v1 = (q(ibc,jbc-1,kbc,3,n) + q(ibc,jbc  ,kbc,3,n)) / 2
                                 w2 = (q(ibc,jbc  ,kbc,4,n) + q(ibc,jbc+1,kbc,4,n)) / 2
                                 w1 = (q(ibc,jbc-1,kbc,4,n) + q(ibc,jbc  ,kbc,4,n)) / 2
                              end if
                              dudeta = u2 - u1
                              dvdeta = v2 - v1
                              dwdeta = w2 - w1
                              
                              if (i3.eq.istrt(ic3)) then
                                 u2 = (  q(ibc,jbc,kbc,2,n) + q(ibc,jbc,kbc+1,2,n)) / 2
                                 u1 = (3*q(ibc,jbc,kbc,2,n) - q(ibc,jbc,kbc+1,2,n)) / 2
                                 v2 = (  q(ibc,jbc,kbc,3,n) + q(ibc,jbc,kbc+1,3,n)) / 2
                                 v1 = (3*q(ibc,jbc,kbc,3,n) - q(ibc,jbc,kbc+1,3,n)) / 2
                                 w2 = (  q(ibc,jbc,kbc,4,n) + q(ibc,jbc,kbc+1,4,n)) / 2
                                 w1 = (3*q(ibc,jbc,kbc,4,n) - q(ibc,jbc,kbc+1,4,n)) / 2
                              else if (i3.eq.iend(ic3)) then
                                 u2 = (3*q(ibc,jbc,kbc,2,n) - q(ibc,jbc,kbc-1,2,n)) / 2
                                 u1 = (  q(ibc,jbc,kbc,2,n) + q(ibc,jbc,kbc-1,2,n)) / 2
                                 v2 = (3*q(ibc,jbc,kbc,3,n) - q(ibc,jbc,kbc-1,3,n)) / 2
                                 v1 = (  q(ibc,jbc,kbc,3,n) + q(ibc,jbc,kbc-1,3,n)) / 2
                                 w2 = (3*q(ibc,jbc,kbc,4,n) - q(ibc,jbc,kbc-1,4,n)) / 2
                                 w1 = (  q(ibc,jbc,kbc,4,n) + q(ibc,jbc,kbc-1,4,n)) / 2
                              else
                                 u2 = (q(ibc,jbc,kbc  ,2,n) + q(ibc,jbc,kbc+1,2,n)) / 2
                                 u1 = (q(ibc,jbc,kbc-1,2,n) + q(ibc,jbc,kbc  ,2,n)) / 2
                                 v2 = (q(ibc,jbc,kbc  ,3,n) + q(ibc,jbc,kbc+1,3,n)) / 2
                                 v1 = (q(ibc,jbc,kbc-1,3,n) + q(ibc,jbc,kbc  ,3,n)) / 2
                                 w2 = (q(ibc,jbc,kbc  ,4,n) + q(ibc,jbc,kbc+1,4,n)) / 2
                                 w1 = (q(ibc,jbc,kbc-1,4,n) + q(ibc,jbc,kbc  ,4,n)) / 2
                              end if
                              dudzeta = u2 - u1
                              dvdzeta = v2 - v1
                              dwdzeta = w2 - w1
                
                              xix   = xideri  (1,im,jm,km,n)
                              xiy   = xideri  (2,im,jm,km,n)
                              xiz   = xideri  (3,im,jm,km,n)
                              etax  = etaderi (1,im,jm,km,n)
                              etay  = etaderi (2,im,jm,km,n)
                              etaz  = etaderi (3,im,jm,km,n)
                              zetax = zetaderi(1,im,jm,km,n)
                              zetay = zetaderi(2,im,jm,km,n)
                              zetaz = zetaderi(3,im,jm,km,n)
                              
                           else if (idir.eq.2) then
                              
                              if (i3.eq.istrt(ic3)) then
                                 u2 = (  q(ibc+1,jbc,kbc,2,n) + q(ibc  ,jbc,kbc,2,n))/2
                                 u1 = (3*q(ibc  ,jbc,kbc,2,n) - q(ibc+1,jbc,kbc,2,n))/2
                                 v2 = (  q(ibc+1,jbc,kbc,3,n) + q(ibc  ,jbc,kbc,3,n))/2
                                 v1 = (3*q(ibc  ,jbc,kbc,3,n) - q(ibc+1,jbc,kbc,3,n))/2
                                 w2 = (  q(ibc+1,jbc,kbc,4,n) + q(ibc  ,jbc,kbc,4,n))/2
                                 w1 = (3*q(ibc  ,jbc,kbc,4,n) - q(ibc+1,jbc,kbc,4,n))/2
                              else if (i3.eq.iend(ic3)) then
                                 u2 = (3*q(ibc,jbc,kbc,2,n) - q(ibc-1,jbc,kbc,2,n)) / 2
                                 u1 = (  q(ibc,jbc,kbc,2,n) + q(ibc-1,jbc,kbc,2,n)) / 2
                                 v2 = (3*q(ibc,jbc,kbc,3,n) - q(ibc-1,jbc,kbc,3,n)) / 2
                                 v1 = (  q(ibc,jbc,kbc,3,n) + q(ibc-1,jbc,kbc,3,n)) / 2
                                 w2 = (3*q(ibc,jbc,kbc,4,n) - q(ibc-1,jbc,kbc,4,n)) / 2
                                 w1 = (  q(ibc,jbc,kbc,4,n) + q(ibc-1,jbc,kbc,4,n)) / 2
                              else
                                 u2 = (q(ibc+1,jbc,kbc,2,n) + q(ibc  ,jbc,kbc,2,n)) / 2
                                 u1 = (q(ibc  ,jbc,kbc,2,n) + q(ibc-1,jbc,kbc,2,n)) / 2
                                 v2 = (q(ibc+1,jbc,kbc,3,n) + q(ibc  ,jbc,kbc,3,n)) / 2
                                 v1 = (q(ibc  ,jbc,kbc,3,n) + q(ibc-1,jbc,kbc,3,n)) / 2
                                 w2 = (q(ibc+1,jbc,kbc,4,n) + q(ibc  ,jbc,kbc,4,n)) / 2
                                 w1 = (q(ibc  ,jbc,kbc,4,n) + q(ibc-1,jbc,kbc,4,n)) / 2
                              end if
                              dudxi = u2 - u1
                              dvdxi = v2 - v1
                              dwdxi = w2 - w1
                              
                              dudeta = (-1)**(1+1/inrout) * &
                                   (-8*q(ibc,jbc,kbc,2,n) + &
                                   9*q(ibc,jbc+joff1,kbc,2,n) - &
                                   q(ibc,jbc+joff2,kbc,2,n) ) / 3
                              dvdeta = (-1)**(1+1/inrout) * &
                                   (-8*q(ibc,jbc,kbc,3,n) + &
                                   9*q(ibc,jbc+joff1,kbc,3,n) - &
                                   q(ibc,jbc+joff2,kbc,3,n) ) / 3
                              dwdeta = (-1)**(1+1/inrout) * &
                                   (-8*q(ibc,jbc,kbc,4,n) + &
                                   9*q(ibc,jbc+joff1,kbc,4,n) - &
                                   q(ibc,jbc+joff2,kbc,4,n) ) / 3
                              
                              if (i2.eq.istrt(ic2)) then
                                 u2 = (  q(ibc,jbc,kbc+1,2,n) + q(ibc,jbc,kbc  ,2,n))/2
                                 u1 = (3*q(ibc,jbc,kbc  ,2,n) - q(ibc,jbc,kbc+1,2,n))/2
                                 v2 = (  q(ibc,jbc,kbc+1,3,n) + q(ibc,jbc,kbc  ,3,n))/2
                                 v1 = (3*q(ibc,jbc,kbc  ,3,n) - q(ibc,jbc,kbc+1,3,n))/2
                                 w2 = (  q(ibc,jbc,kbc+1,4,n) + q(ibc,jbc,kbc  ,4,n))/2
                                 w1 = (3*q(ibc,jbc,kbc  ,4,n) - q(ibc,jbc,kbc+1,4,n))/2
                              else if (i2.eq.iend(ic2)) then
                                 u2 = (3*q(ibc,jbc,kbc,2,n) - q(ibc,jbc,kbc-1,2,n)) / 2
                                 u1 = (  q(ibc,jbc,kbc,2,n) + q(ibc,jbc,kbc-1,2,n)) / 2
                                 v2 = (3*q(ibc,jbc,kbc,3,n) - q(ibc,jbc,kbc-1,3,n)) / 2
                                 v1 = (  q(ibc,jbc,kbc,3,n) + q(ibc,jbc,kbc-1,3,n)) / 2
                                 w2 = (3*q(ibc,jbc,kbc,4,n) - q(ibc,jbc,kbc-1,4,n)) / 2
                                 w1 = (  q(ibc,jbc,kbc,4,n) + q(ibc,jbc,kbc-1,4,n)) / 2
                              else
                                 u2 = (q(ibc,jbc,kbc+1,2,n) + q(ibc,jbc,kbc  ,2,n)) / 2
                                 u1 = (q(ibc,jbc,kbc  ,2,n) + q(ibc,jbc,kbc-1,2,n)) / 2
                                 v2 = (q(ibc,jbc,kbc+1,3,n) + q(ibc,jbc,kbc  ,3,n)) / 2
                                 v1 = (q(ibc,jbc,kbc  ,3,n) + q(ibc,jbc,kbc-1,3,n)) / 2
                                 w2 = (q(ibc,jbc,kbc+1,4,n) + q(ibc,jbc,kbc  ,4,n)) / 2
                                 w1 = (q(ibc,jbc,kbc  ,4,n) + q(ibc,jbc,kbc-1,4,n)) / 2
                              end if
                              dudzeta = u2 - u1
                              dvdzeta = v2 - v1
                              dwdzeta = w2 - w1
                              
                              xix   = xiderj  (1,im,jm,km,n)
                              xiy   = xiderj  (2,im,jm,km,n)
                              xiz   = xiderj  (3,im,jm,km,n)
                              etax  = etaderj (1,im,jm,km,n)
                              etay  = etaderj (2,im,jm,km,n)
                              etaz  = etaderj (3,im,jm,km,n)
                              zetax = zetaderj(1,im,jm,km,n)
                              zetay = zetaderj(2,im,jm,km,n)
                              zetaz = zetaderj(3,im,jm,km,n)
                              
                           else if (idir.eq.3) then
                              if (i2.eq.istrt(ic2)) then
                                 u2 = (  q(ibc,jbc,kbc,2,n) + q(ibc+1,jbc,kbc,2,n)) / 2
                                 u1 = (3*q(ibc,jbc,kbc,2,n) - q(ibc+1,jbc,kbc,2,n)) / 2
                                 v2 = (  q(ibc,jbc,kbc,3,n) + q(ibc+1,jbc,kbc,3,n)) / 2
                                 v1 = (3*q(ibc,jbc,kbc,3,n) - q(ibc+1,jbc,kbc,3,n)) / 2
                                 w2 = (  q(ibc,jbc,kbc,4,n) + q(ibc+1,jbc,kbc,4,n)) / 2
                                 w1 = (3*q(ibc,jbc,kbc,4,n) - q(ibc+1,jbc,kbc,4,n)) / 2
                              else if (i2.eq.iend(ic2)) then
                                 u2 = (3*q(ibc,jbc,kbc,2,n) - q(ibc-1,jbc,kbc,2,n)) / 2
                                 u1 = (  q(ibc,jbc,kbc,2,n) + q(ibc-1,jbc,kbc,2,n)) / 2
                                 v2 = (3*q(ibc,jbc,kbc,3,n) - q(ibc-1,jbc,kbc,3,n)) / 2
                                 v1 = (  q(ibc,jbc,kbc,3,n) + q(ibc-1,jbc,kbc,3,n)) / 2
                                 w2 = (3*q(ibc,jbc,kbc,4,n) - q(ibc-1,jbc,kbc,4,n)) / 2
                                 w1 = (  q(ibc,jbc,kbc,4,n) + q(ibc-1,jbc,kbc,4,n)) / 2
                              else
                                 u2 = (q(ibc  ,jbc,kbc,2,n) + q(ibc+1,jbc,kbc,2,n)) / 2
                                 u1 = (q(ibc-1,jbc,kbc,2,n) + q(ibc  ,jbc,kbc,2,n)) / 2
                                 v2 = (q(ibc  ,jbc,kbc,3,n) + q(ibc+1,jbc,kbc,3,n)) / 2
                                 v1 = (q(ibc-1,jbc,kbc,3,n) + q(ibc  ,jbc,kbc,3,n)) / 2
                                 w2 = (q(ibc  ,jbc,kbc,4,n) + q(ibc+1,jbc,kbc,4,n)) / 2
                                 w1 = (q(ibc-1,jbc,kbc,4,n) + q(ibc  ,jbc,kbc,4,n)) / 2
                              end if
                              dudxi = u2 - u1
                              dvdxi = v2 - v1
                              dwdxi = w2 - w1
                              
                              if (i3.eq.istrt(ic3)) then
                                 u2 = (  q(ibc,jbc,kbc,2,n) + q(ibc,jbc+1,kbc,2,n)) / 2
                                 u1 = (3*q(ibc,jbc,kbc,2,n) - q(ibc,jbc+1,kbc,2,n)) / 2
                                 v2 = (  q(ibc,jbc,kbc,3,n) + q(ibc,jbc+1,kbc,3,n)) / 2
                                 v1 = (3*q(ibc,jbc,kbc,3,n) - q(ibc,jbc+1,kbc,3,n)) / 2
                                 w2 = (  q(ibc,jbc,kbc,4,n) + q(ibc,jbc+1,kbc,4,n)) / 2
                                 w1 = (3*q(ibc,jbc,kbc,4,n) - q(ibc,jbc+1,kbc,4,n)) / 2
                              else if (i3.eq.iend(ic3)) then
                                 u2 = (3*q(ibc,jbc,kbc,2,n) - q(ibc,jbc-1,kbc,2,n)) / 2
                                 u1 = (  q(ibc,jbc,kbc,2,n) + q(ibc,jbc-1,kbc,2,n)) / 2
                                 v2 = (3*q(ibc,jbc,kbc,3,n) - q(ibc,jbc-1,kbc,3,n)) / 2
                                 v1 = (  q(ibc,jbc,kbc,3,n) + q(ibc,jbc-1,kbc,3,n)) / 2
                                 w2 = (3*q(ibc,jbc,kbc,4,n) - q(ibc,jbc-1,kbc,4,n)) / 2
                                 w1 = (  q(ibc,jbc,kbc,4,n) + q(ibc,jbc-1,kbc,4,n)) / 2
                              else
                                 u2 = (q(ibc,jbc  ,kbc,2,n) + q(ibc,jbc+1,kbc,2,n)) / 2
                                 u1 = (q(ibc,jbc-1,kbc,2,n) + q(ibc,jbc  ,kbc,2,n)) / 2
                                 v2 = (q(ibc,jbc  ,kbc,3,n) + q(ibc,jbc+1,kbc,3,n)) / 2
                                 v1 = (q(ibc,jbc-1,kbc,3,n) + q(ibc,jbc  ,kbc,3,n)) / 2
                                 w2 = (q(ibc,jbc  ,kbc,4,n) + q(ibc,jbc+1,kbc,4,n)) / 2
                                 w1 = (q(ibc,jbc-1,kbc,4,n) + q(ibc,jbc  ,kbc,4,n)) / 2
                              end if
                              dudeta = u2 - u1
                              dvdeta = v2 - v1
                              dwdeta = w2 - w1
                              
                              dudzeta = (-1)**(1+1/inrout) * &
                                   (-8*q(ibc,jbc,kbc,2,n) + &
                                   9*q(ibc,jbc,kbc+koff1,2,n) - &
                                   q(ibc,jbc,kbc+koff2,2,n) ) / 3
                              dvdzeta = (-1)**(1+1/inrout) * &
                                   (-8*q(ibc,jbc,kbc,3,n) + &
                                   9*q(ibc,jbc,kbc+koff1,3,n) - &
                                   q(ibc,jbc,kbc+koff2,3,n) ) / 3
                              dwdzeta = (-1)**(1+1/inrout) * &
                                   (-8*q(ibc,jbc,kbc,4,n) + &
                                   9*q(ibc,jbc,kbc+koff1,4,n) - &
                                   q(ibc,jbc,kbc+koff2,4,n) ) / 3
                              
                              xix   = xiderk  (1,im,jm,km,n)
                              xiy   = xiderk  (2,im,jm,km,n)
                              xiz   = xiderk  (3,im,jm,km,n)
                              etax  = etaderk (1,im,jm,km,n)
                              etay  = etaderk (2,im,jm,km,n)
                              etaz  = etaderk (3,im,jm,km,n)
                              zetax = zetaderk(1,im,jm,km,n)
                              zetay = zetaderk(2,im,jm,km,n)
                              zetaz = zetaderk(3,im,jm,km,n)
                              
                           end if
                           
                           dudx = dudxi*xix + dudeta*etax + dudzeta*zetax
                           dudy = dudxi*xiy + dudeta*etay + dudzeta*zetay
                           dudz = dudxi*xiz + dudeta*etaz + dudzeta*zetaz
                           dvdx = dvdxi*xix + dvdeta*etax + dvdzeta*zetax
                           dvdy = dvdxi*xiy + dvdeta*etay + dvdzeta*zetay
                           dvdz = dvdxi*xiz + dvdeta*etaz + dvdzeta*zetaz
                           dwdx = dwdxi*xix + dwdeta*etax + dwdzeta*zetax
                           dwdy = dwdxi*xiy + dwdeta*etay + dwdzeta*zetay
                           dwdz = dwdxi*xiz + dwdeta*etaz + dwdzeta*zetaz
                           
                           divv = dudx + dvdy + dwdz
                           
                           txx = 2*muw * (dudx - divv/3)
                           txy =   muw * (dudy + dvdx)
                           txz =   muw * (dudz + dwdx)
                           tyy = 2*muw * (dvdy - divv/3)
                           tyz =   muw * (dvdz + dwdy)
                           tzz = 2*muw * (dwdz - divv/3)
                           
                           tauwx = txx*kx + txy*ky + txz*kz
                           tauwy = txy*kx + tyy*ky + tyz*kz
                           tauwz = txz*kx + tyz*ky + tzz*kz
                           
                           tauwpx = tauwx - (tauwx*kx+tauwy*ky+tauwz*kz)*kx
                           tauwpy = tauwy - (tauwx*kx+tauwy*ky+tauwz*kz)*ky
                           tauwpz = tauwz - (tauwx*kx+tauwy*ky+tauwz*kz)*kz
                           
                           tauw   = dsqrt(tauwpx**2+tauwpy**2+tauwpz**2)
                           utau   = dsqrt(tauw/rhow)
                           
                           cf(ibc,jbc,kbc,n) = 2*tauw / (reyno*machfs)
                           yp(ibc,jbc,kbc,n) = dsqrt(reyno/machfs)* &
                                rhow*dn1*utau/muw
                           
                        end if
                        
                     end do
                  end do
!------------extract/calculate surface coefficients - END
                  
!------------write tecplot surface patch
                  
                  if (harbal) then
     write(line1,'(''ZONE T="arturo HB, mode '',i2,''", I='',i4,'', J='',i4,'', K='',i4,'',F=BLOCK, DT=(DOUBLE DOUBLE DOUBLE &
     &DOUBLE DOUBLE DOUBLE)'')') &
                          n,ibc_e-ibc_s+1,jbc_e-jbc_s+1,kbc_e-kbc_s+1
                  else
     write(line1,'(''ZONE T="arturo",I='',i4,'', J='',i4,'', K='',i4,'',F=BLOCK, &
     &DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)'')') &
                          ibc_e-ibc_s+1,jbc_e-jbc_s+1,kbc_e-kbc_s+1
                  end if
                  
                  write(fid(n),'(a)') line1
                  
                  write (fid(n),10) (((x(i,j,k,n),i=ibc_s,ibc_e),j=jbc_s,jbc_e), &
                       k=kbc_s,kbc_e)
                  write (fid(n),10) (((y(i,j,k,n),i=ibc_s,ibc_e),j=jbc_s,jbc_e), &
                       k=kbc_s,kbc_e)
                  write (fid(n),10) (((z(i,j,k,n),i=ibc_s,ibc_e),j=jbc_s,jbc_e), &
                       k=kbc_s,kbc_e)
                  write (fid(n),10) (((cp(i,j,k,n),i=ibc_s,ibc_e),j=jbc_s,jbc_e), &
                       k=kbc_s,kbc_e)
                  if (viscous) then
                     write (fid(n),10) (((cf(i,j,k,n),i=ibc_s,ibc_e),j=jbc_s,jbc_e), &
                          k=kbc_s,kbc_e)
                     write (fid(n),10) (((yp(i,j,k,n),i=ibc_s,ibc_e),j=jbc_s,jbc_e), &
                          k=kbc_s,kbc_e)
                  else
                     write (fid(n),10) (((0.d0       ,i=ibc_s,ibc_e),j=jbc_s,jbc_e), &
                          k=kbc_s,kbc_e)
                     write (fid(n),10) (((0.d0       ,i=ibc_s,ibc_e),j=jbc_s,jbc_e), &
                          k=kbc_s,kbc_e)
                  end if
                  
                  deallocate(cp)
                  if (viscous) then
                     deallocate(cf,yp)
                  end if
                  
               end if
               
!--------end loop (ibc1=1,nbcs)
            end do
            
!------loop over surfaces - END
         end do
         
!---- loop over harmonics - END
      end do
      
      
 14   format(3e22.14)
 10   format(3e22.14)
      
      return
      end      

!-----------------------------------------------------------------------
      subroutine wr_cgns_surf(nl,q,x,y,z,xdot,ydot,zdot,dist,si,sj,sk,xideri, &
        xiderj,xiderk,etaderi,etaderj,etaderk,zetaderi,zetaderj,zetaderk, &
        bctopo)
!-----------------------------------------------------------------------
!     
!-----------------------------------------------------------------------
      
      use parallelutils, only: surfblockstarts, surfblockends, &
           surfbcindex, surfblockindex
      use cgns
      use common_variables
      use cosa_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,ierr,bctopo(*),surfacenumber
      integer(kind=cosa_int) nl,iblock,iq,iblk,ixyz,iimt,ivmt,idist,n,fid(0:2*mharms)
      logical parallel,amcontrol
      real (kind=cosa_real) q(*),x(*),y(*),z(*),xdot(*),ydot(*),zdot(*),dist(*), &
        si(*),sj(*),sk(*),xideri(*),xiderj(*),xiderk(*),etaderi(*),etaderj(*), &
        etaderk(*),zetaderi(*),zetaderj(*),zetaderk(*)
      character*72 line,tag
      integer :: basenum
      double precision :: starttime, endtime, totaltime, maxtime, mintime

      call usingparallel(parallel)
      call amcontroller(amcontrol)

      call getmpitime(starttime)

      do n = 0,2*nharms

!------ build tecplot file(s) name

        if (.not.unsteady) then

          surftec='surf_steady.cgns'

        else if (harbal) then

          if (n.le.9) then
            write(tag,'(''00'',i1)') n
          else if (n.le.99) then
            write(tag,'(''0'',i2)') n
          else if (n.le.999) then
            write(tag,'(i3)') n
          end if
          surftec = 'surf_hb_'//tag
          surftec = trim(surftec)//'.cgns'

        else if (dualt) then

          if (itime.le.9) then
            write(tag,'(''000'',i1)') itime
          elseif (itime.le.99) then
            write(tag,'(''00'',i2)')  itime
          elseif (itime.le.999) then
            write(tag,'(''0'',i3)')   itime
          else
            write(tag,'(i4)')         itime
          end if
          tag = 'surf_t_'//tag
          surftec = trim(tag)//'.cgns'

        else if (rgkuns) then

          surftec='surf_tec_urk.cgns'

        end if

!------ open tecplot file(s) and write 2-line file header

        fid(n)    = 200 + n

        if(parallel) then

           call writeparallelsurfcgnsheader(fid(n),basenum,n,nl)

           call cgp_open_f(surftec,CG_MODE_MODIFY,fid(n),ierr)
           if(ierr .ne. CG_OK) then
              write(*,*) 'cgp_open_f error wr_cgns_surf'
              call cg_error_print_f()
           end if
  
        else

           surftec = trim(surftec)
           call cg_open_f(surftec,CG_MODE_WRITE,fid(n),ierr)
           if(ierr .ne. CG_OK) then
              write(*,*) 'cg_open_f error on file ',surftec
              call cg_error_print_f()
           end if
            
           call cg_base_write_f(fid(n),'gridbase',3,3,basenum,ierr)
           if(ierr .ne. CG_OK) then
              write(*,*) 'cg_base_write_f error'
              call cg_error_print_f()
           end if

        end if
        
      end do
      
      surfacenumber = 1

      do iblock = 1,mynblocks

        imax   = i_imax     (iblock,nl)
        jmax   = j_jmax     (iblock,nl)
        kmax   = k_kmax     (iblock,nl)
        iq     = 1 + off_p3 (iblock,nl) * npde * dim5
        iblk   = 1 + off_bct(iblock,nl)
        ixyz   = 1 + off_p2 (iblock,nl)        * dim5h
        iimt   = 1 + off_p2 (iblock,nl) * lmet * dim5h
        ivmt   = 1 + off_0  (iblock,nl) * 3    * dim5h
        idist  = 1 + off_p1 (iblock,nl)
        if (parallel) then

            if(1.lt.0) then

               call write_parallel_wr_cgns_surf_b_M1(fid,imax,jmax,kmax,npde, &
                    nharms,iblock,nl,lmet,q(iq),x(ixyz),y(ixyz),z(ixyz), &
                    xdot(ixyz),ydot(ixyz),zdot(ixyz),si(iimt),sj(iimt),sk(iimt), &
                    xideri(ivmt),etaderj(ivmt),zetaderk(ivmt),dist(idist), &
                    bctopo(iblk),nbcs(iblock),basenum)

            else

               call write_parallel_wr_cgns_surf_b_M2(fid,imax,jmax,kmax,npde, &
                    nharms,iblock,nl,lmet,q(iq),x(ixyz),y(ixyz),z(ixyz), &
                    xdot(ixyz),ydot(ixyz),zdot(ixyz),si(iimt),sj(iimt),sk(iimt), &
                    xideri(ivmt),xiderj(ivmt),xiderk(ivmt),etaderi(ivmt), &
                    etaderj(ivmt),etaderk(ivmt),zetaderi(ivmt), &
                    zetaderj(ivmt),zetaderk(ivmt),dist(idist), &
                    bctopo(iblk),nbcs(iblock),basenum)
               
            end if

      else
          if (1.lt.0) then
             call wr_cgns_surf_b_M1(fid,q(iq),x(ixyz),y(ixyz),z(ixyz), &
                  xdot(ixyz),ydot(ixyz),zdot(ixyz),si(iimt),sj(iimt),sk(iimt), &
                  xideri(ivmt),etaderj(ivmt),zetaderk(ivmt),dist(idist),bctopo(iblk), &
                  imax,jmax,kmax,lmet,nbcs(iblock),npde,nharms, &
                  surfacenumber,basenum)
          else
             call wr_cgns_surf_b_M2(fid,q(iq),x(ixyz),y(ixyz),z(ixyz), &
                  xdot(ixyz),ydot(ixyz),zdot(ixyz),si(iimt),sj(iimt),sk(iimt), &
                  xideri(ivmt),xiderj(ivmt),xiderk(ivmt),etaderi(ivmt), &
                  etaderj(ivmt),etaderk(ivmt),zetaderi(ivmt),zetaderj(ivmt), &
                  zetaderk(ivmt),dist(idist),bctopo(iblk), &
                  imax,jmax,kmax,lmet,nbcs(iblock),npde,nharms, &
                  surfacenumber,basenum)
          end if
       end if
      end do
      
      do n = 0,2*nharms
!------close cgns file(s)
         if(parallel) then
            if(allocated(surfblockstarts)) deallocate(surfblockstarts)
            if(allocated(surfblockends)) deallocate(surfblockends)
            if(allocated(surfblockindex)) deallocate(surfblockindex)            
            if(allocated(surfbcindex)) deallocate(surfbcindex)            
            call closeparallelcgnsfile(fid(n))
         else
            call cg_close_f(fid(n),ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) 'cg_close_f error'
               call cg_error_print_f()
            end if   
         end if
      end do

      call getmpitime(endtime)

      totaltime = endtime-starttime
      maxtime = totaltime
      mintime =  totaltime

      call realmaxreduce(maxtime,1)
      call realminreduce(mintime,1)

      if(amcontrol) then
         write(*,'(A,F8.2,A,2F8.2,A)') 'Writing surface file time: ',totaltime,' (',maxtime,mintime,') seconds'
      end if


!del   15   format(e21.14)

      return
      end

!-----------------------------------------------------------------------
      subroutine wr_cgns_surf_b_M1(fid,q,x,y,z,xdot,ydot,zdot,si,sj,sk,xideri, &
        etaderj,zetaderk,dist,bctopo, &
        imax,jmax,kmax,lmet,nbcs,npde,nharms,surfacenumber,basenum)
!-----------------------------------------------------------------------
       
      use cgns
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,npde,nharms,lmet,nbcs,surfacenumber
      integer(kind=cosa_int) i,imax1,j,jmax1,k,kmax1,ipde,n,nh,ibc_s,jbc_s,kbc_s, &
        ibc_e,jbc_e,kbc_e,basenum
      integer(kind=cosa_int) ibc1
      integer(kind=cosa_int) bctopo(10,nbcs)
      integer(kind=cosa_int) i1,i2,i3,ijkmax(3),idir,istrt(3),iend(3),bctyp, &
        inrout,ibcpt,ibcpt2,ibcn,ibcn2,ibcm,ibc,jbc,kbc,ibc2,jbc2, &
        kbc2,in,jn,kn,in2,jn2,kn2,im,jm,km,ic1,ic2,ic3,sysize,ioff1, &
        ioff2,joff1,joff2,koff1,koff2,ierr
      integer(cgsize_t) ::  sizes(3,3)
      integer fid(0:2*mharms),blocknum
      real (kind=cosa_real) sgnm,nx,ny,nz,kx,ky,kz,rhow,tw,muw,uw,u1,u2,uwr,u1r,u2r, &
           vw,v1,v2,vwr,v1r,v2r,ww,w1,w2,wwr,w1r,w2r,ueta,veta,weta, &
           etax,etay,etaz,dudx,dudy,dudz,dvdx,dvdy, &
           dvdz,dwdx,dwdy,dwdz,divv,txx,txy,txz,tyy,tyz,tzz,mach,tauwx, &
           tauwy,tauwz,tauwpx,tauwpy,tauwpz,dn1,tauw,utau
      real (kind=cosa_real) &
           q   (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           x   ( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           y   ( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           z   ( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           xdot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           ydot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           zdot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           si      (lmet, 0:imax+1, 0:jmax+1, 0:kmax+1,0:2*nharms*hbmove), &
           sj      (lmet, 0:imax+1, 0:jmax+1, 0:kmax+1,0:2*nharms*hbmove), &
           sk      (lmet, 0:imax+1, 0:jmax+1, 0:kmax+1,0:2*nharms*hbmove), &
           xideri  (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           etaderj (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           zetaderk(3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           dist( 0:imax  , 0:jmax  , 0:kmax)
      real (kind=cosa_real),allocatable :: cp(:,:,:),cf(:,:,:),yp(:,:,:)
      real (kind=cosa_real),allocatable :: localcoord(:,:,:)
      character*200 line1
      character*20 zonename
      integer :: xnum, ynum, znum, solnum, indexnum
      integer :: initialsurfacenumber

      ijkmax(1) = imax
      ijkmax(2) = jmax
      ijkmax(3) = kmax

      imax1 = imax+1
      jmax1 = jmax+1
      kmax1 = kmax+1

      initialsurfacenumber = surfacenumber

!---- loop over harmonics - STARTS
      do n = 0,2*nharms
        nh = n*hbmove

!     Required to deal with resetting everything when iterating harmonics.
        surfacenumber = initialsurfacenumber


!------ loop over surfaces - START
!------ MSC 12/03/23. In future nsurface loop can be removed HERE if we do
!                     not wish to distinguish separate surfaces/walls
!                     in surftec file. We will just check for BCs
!                     starting with 15 or 14.
!
!                     Keeping nsurface is instead essential to get
!                     forces and moments in routine force.
!                     Getting nsurface from input.f is trivial when
!                     all walls are either visc. (RANS) or inviscid
!                     (Euler). We just count them. 
!                     NOTE: we will need to use for BC 14 the same syntax as 
!                     for BC 1500s (14xx).
!
!                     Less trivial is case when
!                     a body has both visc. and inv. BCs (e.g. hub).
!                     For this, previous code structure with list of
!                     blocks with boundaries on given body was better.
        
        
        do ibc1=1,nbcs
           if (any(bctopo(1,ibc1).eq. &
                [1500,1501,1502,1503,1400,1401,1402,1403])) then
              
!     store boundary topology in mnemonic names
              bctyp    = bctopo(1,ibc1)
              idir     = bctopo(2,ibc1)
              inrout   = bctopo(3,ibc1)
              istrt(1) = bctopo(4,ibc1)
              iend(1)  = bctopo(5,ibc1)
              istrt(2) = bctopo(6,ibc1)
              iend(2)  = bctopo(7,ibc1)
              istrt(3) = bctopo(8,ibc1)
              iend(3)  = bctopo(9,ibc1)
              
!     set needed variables depending on whether the boundary is
!     the inner boundary (inrout = 1) or
!     the outer boundary (inrout > 1)
!     ibcpt : boundary condition location (first aux. cell)
!     ibcpt2: boundary condition location outside the block from
!     ibcpt (second aux. cell)
!     ibcn  : point to the inside of the block from ibcpt
!     ibcm  : location of the metrics
              
              if (inrout.eq.1) then
                 ibcpt  =  0
                 ibcpt2 = -1
                 ibcn   =  1
                 ibcn2  =  2
                 ibcm   =  1
                 sgnm   = - 1.d0
              else
                 ibcpt  = ijkmax(idir)
                 ibcpt2 = ijkmax(idir) + 1
                 ibcn   = ijkmax(idir) - 1
                 ibcn2  = ijkmax(idir) - 2
                 ibcm   = ijkmax(idir)
                 sgnm   = 1.d0
              end if
              
              if (idir.eq.1) then
                 ioff1 = (-1)**(1+1/inrout)
                 ioff2 = 2*ioff1
                 joff1 = 0
                 joff2 = 0
                 koff1 = 0
                 koff2 = 0
              else if (idir.eq.2) then
                 ioff1 = 0
                 ioff2 = 0
                 joff1 = (-1)**(1+1/inrout)
                 joff2 = 2*joff1
                 koff1 = 0
                 koff2 = 0
              else if (idir.eq.3) then
                 ioff1 = 0
                 ioff2 = 0
                 joff1 = 0
                 joff2 = 0
                 koff1 = (-1)**(1+1/inrout)
                 koff2 = 2*koff1
              end if
              
              ic1 = cyc (idir, 1)
              ic2 = cyc (idir, 2)
              ic3 = cyc (idir, 3)
              
              ibc_s  = ibcpt *krd(ic1,1) + istrt(ic2)*krd(ic2,1) + &
                   istrt(ic3)*krd(ic3,1)
              jbc_s  = ibcpt *krd(ic1,2) + istrt(ic2)*krd(ic2,2) + &
                   istrt(ic3)*krd(ic3,2)
              kbc_s  = ibcpt *krd(ic1,3) + istrt(ic2)*krd(ic2,3) + &
                   istrt(ic3)*krd(ic3,3)
              
              ibc_e  = ibcpt *krd(ic1,1) + iend(ic2)*krd(ic2,1) + &
                   iend(ic3)*krd(ic3,1)
              jbc_e  = ibcpt *krd(ic1,2) + iend(ic2)*krd(ic2,2) + &
                   iend(ic3)*krd(ic3,2)
              kbc_e  = ibcpt *krd(ic1,3) + iend(ic2)*krd(ic2,3) + &
                   iend(ic3)*krd(ic3,3)
              
              sizes(1,1) = (ibc_e-ibc_s)+1
              sizes(2,1) = (jbc_e-jbc_s)+1
              sizes(3,1) = (kbc_e-kbc_s)+1
              
              sizes(1,2) = sizes(1,1)-1
              sizes(2,2) = sizes(2,1)-1
              sizes(3,2) = sizes(3,1)-1
              
              sizes(1,3) = 0
              sizes(2,3) = 0
              sizes(3,3) = 0              
              
              write(zonename, "(A12,I6)") "surfaceblock",surfacenumber
              call cg_zone_write_f(fid(n),basenum,zonename,sizes, &
                   Structured,blocknum,ierr)
              if(ierr .ne. CG_OK) then
                 write(*,*) &
                      'cg_zone_write_f error in cgns surf M1'
                 call cg_error_print_f()
              end if
              
              allocate(cp(ibc_s:ibc_e,jbc_s:jbc_e,kbc_s:kbc_e))
              allocate(cf(ibc_s:ibc_e,jbc_s:jbc_e,kbc_s:kbc_e), &
                   yp(ibc_s:ibc_e,jbc_s:jbc_e,kbc_s:kbc_e))
              allocate(localcoord(ibc_s:ibc_e,jbc_s:jbc_e,kbc_s:kbc_e))

!------------ extract/calculate surface coefficients - START
              do i3 = istrt(ic3),iend(ic3)
                do i2 = istrt(ic2),iend(ic2)

                  ibc  = ibcpt *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
                  jbc  = ibcpt *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
                  kbc  = ibcpt *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

                  im   = ibcm  *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
                  jm   = ibcm  *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
                  km   = ibcm  *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

                  if (idir.eq.1) then
                    nx    = si(1,im,jm,km,nh)
                    ny    = si(2,im,jm,km,nh)
                    nz    = si(3,im,jm,km,nh)
                  else if (idir.eq.2) then
                    nx    = sj(1,im,jm,km,nh)
                    ny    = sj(2,im,jm,km,nh)
                    nz    = sj(3,im,jm,km,nh)
                  else if (idir.eq.3) then
                    nx    = sk(1,im,jm,km,nh)
                    ny    = sk(2,im,jm,km,nh)
                    nz    = sk(3,im,jm,km,nh)
                  end if

                  kx   = (-1)**(1+1/inrout) * nx
                  ky   = (-1)**(1+1/inrout) * ny
                  kz   = (-1)**(1+1/inrout) * nz

                  cp(ibc,jbc,kbc) = 2*(q(ibc,jbc,kbc,5,n)-1/gamma)/ &
                    machfs**2

                  if (viscous) then

                    ibc2 = ibcpt2*krd(ic1,1) +i2*krd(ic2,1) +i3*krd(ic3,1)
                    jbc2 = ibcpt2*krd(ic1,2) +i2*krd(ic2,2) +i3*krd(ic3,2)
                    kbc2 = ibcpt2*krd(Ic1,3) +i2*krd(ic2,3) +i3*krd(ic3,3)

                    in   = ibcn  *krd(ic1,1) +i2*krd(ic2,1) +i3*krd(ic3,1)
                    jn   = ibcn  *krd(ic1,2) +i2*krd(ic2,2) +i3*krd(ic3,2)
                    kn   = ibcn  *krd(ic1,3) +i2*krd(ic2,3) +i3*krd(ic3,3)

                    in2  = ibcn2 *krd(ic1,1) +i2*krd(ic2,1) +i3*krd(ic3,1)
                    jn2  = ibcn2 *krd(ic1,2) +i2*krd(ic2,2) +i3*krd(ic3,2)
                    kn2  = ibcn2 *krd(ic1,3) +i2*krd(ic2,3) +i3*krd(ic3,3)

                    if (idir.eq.1) then
                      etax  = xideri(1,im,jm,km,nh)
                      etay  = xideri(2,im,jm,km,nh)
                      etaz  = xideri(3,im,jm,km,nh)
                    else if (idir.eq.2) then
                      etax  = etaderj(1,im,jm,km,nh)
                      etay  = etaderj(2,im,jm,km,nh)
                      etaz  = etaderj(3,im,jm,km,nh)
                    else if (idir.eq.3) then
                      etax  = zetaderk(1,im,jm,km,nh)
                      etay  = zetaderk(2,im,jm,km,nh)
                      etaz  = zetaderk(3,im,jm,km,nh)
                    end if

                    rhow = q(ibc,jbc,kbc,1,n)
                    tw   = q(ibc,jbc,kbc,5,n) *gamma / q(ibc,jbc,kbc,1,n)
                    muw  = (stemp+1)/(tw+stemp) * (tw**1.5d0)

!------------------ MSC, 11/03/2023: dist not computed for laminar now
                    if (kom.or.kom_bsl.or.kom_sst) then
                      dn1  = dist(in ,jn ,kn )
                    else
                      dn1  = 0
                    end if

                    uw   = q(ibc,jbc,kbc,2,n)
                    u1   = q(ibc+ioff1,jbc+joff1,kbc+koff1,2,n)
                    u2   = q(ibc+ioff2,jbc+joff2,kbc+koff2,2,n)
                    if (moving) then
                      uwr  = xdot(ibc      ,jbc      ,kbc      ,nh)
                      u1r  = xdot(ibc+ioff1,jbc+joff1,kbc+koff1,nh)
                      u2r  = xdot(ibc+ioff2,jbc+joff2,kbc+koff2,nh)
                      uw   = uw - uwr
                      u1   = u1 - u1r
                      u2   = u2 - u2r
                    end if

                    vw   = q(ibc,jbc,kbc,3,n)
                    v1   = q(ibc+ioff1,jbc+joff1,kbc+koff1,3,n)
                    v2   = q(ibc+ioff2,jbc+joff2,kbc+koff2,3,n)
                    if (moving) then
                      vwr  = ydot(ibc      ,jbc      ,kbc      ,nh)
                      v1r  = ydot(ibc+ioff1,jbc+joff1,kbc+koff1,nh)
                      v2r  = ydot(ibc+ioff2,jbc+joff2,kbc+koff2,nh)
                      vw   = vw - vwr
                      v1   = v1 - v1r
                      v2   = v2 - v2r
                    end if

                    ww   = q(ibc,jbc,kbc,4,n)
                    w1   = q(ibc+ioff1,jbc+joff1,kbc+koff1,4,n)
                    w2   = q(ibc+ioff2,jbc+joff2,kbc+koff2,4,n)
                    if (moving) then
                      wwr  = zdot(ibc      ,jbc      ,kbc      ,nh)
                      w1r  = zdot(ibc+ioff1,jbc+joff1,kbc+koff1,nh)
                      w2r  = zdot(ibc+ioff2,jbc+joff2,kbc+koff2,nh)
                      ww   = ww - wwr
                      w1   = w1 - w1r
                      w2   = w2 - w2r
                    end if

                    ueta = (-1)**(1+1/inrout) * (-8*uw + 9*u1 - u2) / 3
                    veta = (-1)**(1+1/inrout) * (-8*vw + 9*v1 - v2) / 3
                    weta = (-1)**(1+1/inrout) * (-8*ww + 9*w1 - w2) / 3

                    dudx = ueta*etax
                    dudy = ueta*etay
                    dudz = ueta*etaz
                    dvdx = veta*etax
                    dvdy = veta*etay
                    dvdz = veta*etaz
                    dwdx = weta*etax
                    dwdy = weta*etay
                    dwdz = weta*etaz
                    divv = dudx + dvdy + dwdz

                    txx = 2*muw * (dudx - divv/3)
                    txy =   muw * (dudy + dvdx)
                    txz =   muw * (dudz + dwdx)
                    tyy = 2*muw * (dvdy - divv/3)
                    tyz =   muw * (dvdz + dwdy)
                    tzz = 2*muw * (dwdz - divv/3)

                    tauwx = txx*kx + txy*ky + txz*kz
                    tauwy = txy*kx + tyy*ky + tyz*kz
                    tauwz = txz*kx + tyz*ky + tzz*kz

                    tauwpx = tauwx - (tauwx*kx+tauwy*ky+tauwz*kz)*kx
                    tauwpy = tauwy - (tauwx*kx+tauwy*ky+tauwz*kz)*ky
                    tauwpz = tauwz - (tauwx*kx+tauwy*ky+tauwz*kz)*kz

                    tauw   = dsqrt(tauwpx**2+tauwpy**2+tauwpz**2)
                    utau   = dsqrt(tauw/rhow)

                    cf(ibc,jbc,kbc) = 2*tauw / (reyno*machfs)
                    yp(ibc,jbc,kbc) = dsqrt(reyno/machfs)* &
                      rhow*dn1*utau/muw

                 end if

                end do
              end do

              if(.not.viscous) then
                 
                 cf = 0.d0
                 yp = 0.d0
                 
              end if
              
!------------ extract/calculate surface coefficients - END
              
!------------write cgns surface patch
            call cg_sol_write_f(fid(n),basenum,surfacenumber,'SurfaceSolutions', &
                 Vertex,solnum,ierr)
!AJ     &           CellCenter,solnum,ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) 'cg_sol_write_f error'
               call cg_error_print_f()
            end if
              localcoord = x(ibc_s:ibc_e,jbc_s:jbc_e,kbc_s:kbc_e,n)
              call cg_coord_write_f(fid(n),basenum,surfacenumber,RealDouble, &
                   'CoordinateX',localcoord,xnum,ierr)
              if(ierr .ne. CG_OK) then
               write(*,*) 'X cg_coord_write_f error'
               call cg_error_print_f()
            end if
              localcoord = y(ibc_s:ibc_e,jbc_s:jbc_e,kbc_s:kbc_e,n)
            call cg_coord_write_f(fid(n),basenum,surfacenumber,RealDouble, &
                 'CoordinateY',localcoord,ynum,ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) 'Y cg_coord_write_f error'
               call cg_error_print_f()
            end if
              localcoord = z(ibc_s:ibc_e,jbc_s:jbc_e,kbc_s:kbc_e,n)
            call cg_coord_write_f(fid(n),basenum,surfacenumber,RealDouble, &
                 'CoordinateZ',localcoord,znum,ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) 'Z cg_coord_write_f error'
               call cg_error_print_f()
            end if
            call cg_field_write_f(fid(n),basenum,surfacenumber,solnum, &
                 RealDouble,'StaticPressureCoefficient', &
                 cp(ibc_s:ibc_e,jbc_s:jbc_e,kbc_s:kbc_e),indexnum,ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) 'cg_field_write_f StaticPressureCoefficient error'
               call cg_error_print_f()
            end if
            call cg_field_write_f(fid(n),basenum,surfacenumber,solnum, &
                 RealDouble,'SkinFrictionCoefficient', &
                 cf(ibc_s:ibc_e,jbc_s:jbc_e,kbc_s:kbc_e),indexnum,ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) 'cg_field_write_f SkinFrictionCoefficient error'
               call cg_error_print_f()
            end if
            call cg_field_write_f(fid(n),basenum,surfacenumber,solnum, &
                 RealDouble,'YPlus', &
                 yp(ibc_s:ibc_e,jbc_s:jbc_e,kbc_s:kbc_e),indexnum,ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) 'cg_field_write_f YPlus error'
               call cg_error_print_f()
            end if

            deallocate(cp,cf,yp,localcoord)
            
            surfacenumber = surfacenumber + 1

            end if

!-------- end loop (ibc1=1,nbcs)
          end do
  
!---- loop over harmonics - END
      end do


 14   format(3e22.14)
 10   format(3e22.14)

      return
      end

!-----------------------------------------------------------------------
      subroutine wr_cgns_surf_b_M2(fid,q,x,y,z,xdot,ydot,zdot,si,sj,sk,xideri, &
        xiderj,xiderk,etaderi,etaderj,etaderk,zetaderi,zetaderj,zetaderk,dist, &
        bctopo,imax,jmax,kmax,lmet,nbcs,npde,nharms,surfacenumber,basenum)
!-----------------------------------------------------------------------
       
      use cgns
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,npde,nharms,lmet,nbcs
      integer(kind=cosa_int) i,imax1,j,jmax1,k,kmax1,ipde,n,nh,ibc_s,jbc_s,kbc_s, &
        ibc_e,jbc_e,kbc_e
      integer(kind=cosa_int) ibc1,surfacenumber,blocknum,basenum
      integer(kind=cosa_int) bctopo(10,nbcs)
      integer(kind=cosa_int) i1,i2,i3,ijkmax(3),idir,istrt(3),iend(3),bctyp, &
        inrout,ibcpt,ibcpt2,ibcn,ibcn2,ibcm,ibc,jbc,kbc,ibc2,jbc2, &
        kbc2,in,jn,kn,in2,jn2,kn2,im,jm,km,ic1,ic2,ic3,sysize,ioff1, &
        ioff2,joff1,joff2,koff1,koff2,ierr
      integer(cgsize_t) ::  sizes(3,3)
      integer fid(0:2*mharms)
      real (kind=cosa_real) sgnm,nx,ny,nz,kx,ky,kz,rhow,tw,muw,uw,u1,u2,uwr,u1r,u2r, &
        vw,v1,v2,vwr,v1r,v2r,ww,w1,w2,wwr,w1r,w2r, &
        xix,xiy,xiz,etax,etay,etaz,zetax,zetay,zetaz, &
        dudxi,dudeta,dudzeta,dvdxi,dvdeta,dvdzeta,dwdxi,dwdeta,dwdzeta, &
        dudx,dudy,dudz,dvdx,dvdy, &
        dvdz,dwdx,dwdy,dwdz,divv,txx,txy,txz,tyy,tyz,tzz,mach,tauwx, &
        tauwy,tauwz,tauwpx,tauwpy,tauwpz,dn1,tauw,utau
      real (kind=cosa_real) &
           q   (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           x   ( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           y   ( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           z   ( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           xdot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           ydot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           zdot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           si      (lmet, 0:imax+1, 0:jmax+1, 0:kmax+1,0:2*nharms*hbmove), &
           sj      (lmet, 0:imax+1, 0:jmax+1, 0:kmax+1,0:2*nharms*hbmove), &
           sk      (lmet, 0:imax+1, 0:jmax+1, 0:kmax+1,0:2*nharms*hbmove), &
           xideri  (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           xiderj  (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           xiderk  (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           etaderi (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           etaderj (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           etaderk (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           zetaderi(3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           zetaderj(3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           zetaderk(3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           dist( 0:imax  , 0:jmax  , 0:kmax)
      real (kind=cosa_real),allocatable :: cp(:,:,:),cf(:,:,:),yp(:,:,:)
      real (kind=cosa_real),allocatable :: localcoord(:,:,:)
      character*200 line1
      character*20 zonename
      integer :: xnum, ynum, znum, solnum, indexnum
      integer :: initialsurfacenumber

      ijkmax(1) = imax
      ijkmax(2) = jmax
      ijkmax(3) = kmax

      imax1 = imax+1
      jmax1 = jmax+1
      kmax1 = kmax+1

      initialsurfacenumber = surfacenumber

!---- loop over harmonics - STARTS
      do n = 0,2*nharms
         nh = n*hbmove
         
!     Required to deal with resetting everything when iterating harmonics.
         surfacenumber = initialsurfacenumber

!------loop over surfaces - START
!------MSC 12/03/23. In future nbofy loop can be removed HERE if we do
!     not wish to distinguish separate surfaces/walls
!     in surftec file. We will just check for BCs
!     starting with 15 or 14.
!     
!     Keeping nbpdy is instead essential to get
!     forces and moments in routine force.
!     Getting nsurface from input.f is trivial when
!     all walls are either visc. (RANS) or inviscid
!     (Euler). We just count them. 
!     NOTE: we will need to use for BC 14 the same syntax as 
!     for BC 1500s (14xx).
!     
!     Less trivial is case when
!     a body has both visc. and inv. BCs (e.g. hub).
!     For this, previous code structure with list of
!     blocks with boundaries on given body was better.
         
         do ibc1=1,nbcs
            if (any(bctopo(1,ibc1).eq. &
                 [1500,1501,1502,1503,1400,1401,1402,1403])) then
                              
!     store boundary topology in mnemonic names
               bctyp    = bctopo(1,ibc1)
               idir     = bctopo(2,ibc1)
               inrout   = bctopo(3,ibc1)
               istrt(1) = bctopo(4,ibc1)
               iend(1)  = bctopo(5,ibc1)
               istrt(2) = bctopo(6,ibc1)
               iend(2)  = bctopo(7,ibc1)
               istrt(3) = bctopo(8,ibc1)
               iend(3)  = bctopo(9,ibc1)
               
!     set needed variables depending on whether the boundary is
!     the inner boundary (inrout = 1) or
!     the outer boundary (inrout > 1)
!     ibcpt : boundary condition location (first aux. cell)
!     ibcpt2: boundary condition location outside the block from
!     ibcpt (second aux. cell)
!     ibcn  : point to the inside of the block from ibcpt
!     ibcm  : location of the metrics
               
               if (inrout.eq.1) then
                  ibcpt  =  0
!     del            ibcpt2 = -1
                  ibcn   =  1
                  ibcn2  =  2
                  ibcm   =  1
!     del            sgnm   = - 1.d0
               else
                  ibcpt  = ijkmax(idir)
!     del            ibcpt2 = ijkmax(idir) + 1
                  ibcn   = ijkmax(idir) - 1
                  ibcn2  = ijkmax(idir) - 2
                  ibcm   = ijkmax(idir)
!     del            sgnm   = 1.d0
               end if
               
               if (idir.eq.1) then
                  ioff1 = (-1)**(1+1/inrout)
                  ioff2 = 2*ioff1
                  joff1 = 0
                  joff2 = 0
                  koff1 = 0
                  koff2 = 0
               else if (idir.eq.2) then
                  ioff1 = 0
                  ioff2 = 0
                  joff1 = (-1)**(1+1/inrout)
                  joff2 = 2*joff1
                  koff1 = 0
                  koff2 = 0
               else if (idir.eq.3) then
                  ioff1 = 0
                  ioff2 = 0
                  joff1 = 0
                  joff2 = 0
                  koff1 = (-1)**(1+1/inrout)
                  koff2 = 2*koff1
               end if
               
               ic1 = cyc (idir, 1)
               ic2 = cyc (idir, 2)
               ic3 = cyc (idir, 3)
               
               ibc_s  = ibcpt *krd(ic1,1) + istrt(ic2)*krd(ic2,1) + &
                    istrt(ic3)*krd(ic3,1)
               jbc_s  = ibcpt *krd(ic1,2) + istrt(ic2)*krd(ic2,2) + &
                    istrt(ic3)*krd(ic3,2)
               kbc_s  = ibcpt *krd(ic1,3) + istrt(ic2)*krd(ic2,3) + &
                    istrt(ic3)*krd(ic3,3)
               
               ibc_e  = ibcpt *krd(ic1,1) + iend(ic2)*krd(ic2,1) + &
                    iend(ic3)*krd(ic3,1)
               jbc_e  = ibcpt *krd(ic1,2) + iend(ic2)*krd(ic2,2) + &
                    iend(ic3)*krd(ic3,2)
               kbc_e  = ibcpt *krd(ic1,3) + iend(ic2)*krd(ic2,3) + &
                    iend(ic3)*krd(ic3,3)
               
               sizes(1,1) = (ibc_e-ibc_s)+1
               sizes(2,1) = (jbc_e-jbc_s)+1
               sizes(3,1) = (kbc_e-kbc_s)+1
               
               sizes(1,2) = sizes(1,1)-1
               sizes(2,2) = sizes(2,1)-1
               sizes(3,2) = sizes(3,1)-1
               
               sizes(1,3) = 0
               sizes(2,3) = 0
               sizes(3,3) = 0              
               
               write(zonename, "(A12,I6)") "surfaceblock",surfacenumber
               call cg_zone_write_f(fid(n),basenum,zonename,sizes, &
                    Structured,blocknum,ierr)
               if(blocknum .ne. surfacenumber) then
                  write(*,*) 'Zone number different than the expected: ', &
                       blocknum,surfacenumber
                  stop
               end if
               if(ierr .ne. CG_OK) then
                  write(*,*) &
                       'cg_zone_write_f error in cgns surf M2'
                  call cg_error_print_f()
               end if
               
               allocate(cp(ibc_s:ibc_e,jbc_s:jbc_e,kbc_s:kbc_e))
               allocate(cf(ibc_s:ibc_e,jbc_s:jbc_e,kbc_s:kbc_e), &
                    yp(ibc_s:ibc_e,jbc_s:jbc_e,kbc_s:kbc_e))
               allocate(localcoord(ibc_s:ibc_e,jbc_s:jbc_e,kbc_s:kbc_e))
!------------extract/calculate surface coefficients - START
               do i3 = istrt(ic3),iend(ic3)
                  do i2 = istrt(ic2),iend(ic2)
                     
                     ibc  = ibcpt *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
                     jbc  = ibcpt *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
                     kbc  = ibcpt *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)
                     
                     im   = ibcm  *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
                     jm   = ibcm  *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
                     km   = ibcm  *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)
                     
                     if (idir.eq.1) then
                        nx    = si(1,im,jm,km,nh)
                        ny    = si(2,im,jm,km,nh)
                        nz    = si(3,im,jm,km,nh)
                     else if (idir.eq.2) then
                        nx    = sj(1,im,jm,km,nh)
                        ny    = sj(2,im,jm,km,nh)
                        nz    = sj(3,im,jm,km,nh)
                     else if (idir.eq.3) then
                        nx    = sk(1,im,jm,km,nh)
                        ny    = sk(2,im,jm,km,nh)
                        nz    = sk(3,im,jm,km,nh)
                     end if
                     
                     kx   = (-1)**(1+1/inrout) * nx
                     ky   = (-1)**(1+1/inrout) * ny
                     kz   = (-1)**(1+1/inrout) * nz
                     
                     cp(ibc,jbc,kbc) = 2*(q(ibc,jbc,kbc,5,n)-1/gamma)/ &
                          machfs**2
                     
                     if (viscous) then
                        
!     del                ibc2 = ibcpt2*krd(ic1,1) +i2*krd(ic2,1) +i3*krd(ic3,1)
!     del                jbc2 = ibcpt2*krd(ic1,2) +i2*krd(ic2,2) +i3*krd(ic3,2)
!     del                kbc2 = ibcpt2*krd(Ic1,3) +i2*krd(ic2,3) +i3*krd(ic3,3)
                        
                        in   = ibcn  *krd(ic1,1) +i2*krd(ic2,1) +i3*krd(ic3,1)
                        jn   = ibcn  *krd(ic1,2) +i2*krd(ic2,2) +i3*krd(ic3,2)
                        kn   = ibcn  *krd(ic1,3) +i2*krd(ic2,3) +i3*krd(ic3,3)
                        
                        in2  = ibcn2 *krd(ic1,1) +i2*krd(ic2,1) +i3*krd(ic3,1)
                        jn2  = ibcn2 *krd(ic1,2) +i2*krd(ic2,2) +i3*krd(ic3,2)
                        kn2  = ibcn2 *krd(ic1,3) +i2*krd(ic2,3) +i3*krd(ic3,3)
                        
                        rhow = q(ibc,jbc,kbc,1,n)
                        tw   = q(ibc,jbc,kbc,5,n) *gamma / q(ibc,jbc,kbc,1,n)
                        muw  = (stemp+1)/(tw+stemp) * (tw**1.5d0)
                        
!------------------MSC, 11/03/2023: dist not computed for laminar now
                        if (kom.or.kom_bsl.or.kom_sst) then
                           dn1  = dist(in ,jn ,kn )
                        else
                           dn1  = 0
                        end if
                        
                        
                        if (idir.eq.1) then
                           
                           dudxi = (-1)**(1+1/inrout) * &
                                (-8*q(ibc,jbc,kbc,2,n) + &
                                9*q(ibc+ioff1,jbc,kbc,2,n) - &
                                q(ibc+ioff2,jbc,kbc,2,n) ) / 3
                           dvdxi = (-1)**(1+1/inrout) * &
                                (-8*q(ibc,jbc,kbc,3,n) + &
                                9*q(ibc+ioff1,jbc,kbc,3,n) - &
                                q(ibc+ioff2,jbc,kbc,3,n) ) / 3
                           dwdxi = (-1)**(1+1/inrout) * &
                                (-8*q(ibc,jbc,kbc,4,n) + &
                                9*q(ibc+ioff1,jbc,kbc,4,n) - &
                                q(ibc+ioff2,jbc,kbc,4,n) ) / 3
                           
                           if (i2.eq.istrt(ic2)) then
                              u2 = (  q(ibc,jbc,kbc,2,n) + q(ibc,jbc+1,kbc,2,n)) / 2
                              u1 = (3*q(ibc,jbc,kbc,2,n) - q(ibc,jbc+1,kbc,2,n)) / 2
                              v2 = (  q(ibc,jbc,kbc,3,n) + q(ibc,jbc+1,kbc,3,n)) / 2
                              v1 = (3*q(ibc,jbc,kbc,3,n) - q(ibc,jbc+1,kbc,3,n)) / 2
                              w2 = (  q(ibc,jbc,kbc,4,n) + q(ibc,jbc+1,kbc,4,n)) / 2
                              w1 = (3*q(ibc,jbc,kbc,4,n) - q(ibc,jbc+1,kbc,4,n)) / 2
                           else if (i2.eq.iend(ic2)) then
                              u2 = (3*q(ibc,jbc,kbc,2,n) - q(ibc,jbc-1,kbc,2,n)) / 2
                              u1 = (  q(ibc,jbc,kbc,2,n) + q(ibc,jbc-1,kbc,2,n)) / 2
                              v2 = (3*q(ibc,jbc,kbc,3,n) - q(ibc,jbc-1,kbc,3,n)) / 2
                              v1 = (  q(ibc,jbc,kbc,3,n) + q(ibc,jbc-1,kbc,3,n)) / 2
                              w2 = (3*q(ibc,jbc,kbc,4,n) - q(ibc,jbc-1,kbc,4,n)) / 2
                              w1 = (  q(ibc,jbc,kbc,4,n) + q(ibc,jbc-1,kbc,4,n)) / 2
                           else
                              u2 = (q(ibc,jbc  ,kbc,2,n) + q(ibc,jbc+1,kbc,2,n)) / 2
                              u1 = (q(ibc,jbc-1,kbc,2,n) + q(ibc,jbc  ,kbc,2,n)) / 2
                              v2 = (q(ibc,jbc  ,kbc,3,n) + q(ibc,jbc+1,kbc,3,n)) / 2
                              v1 = (q(ibc,jbc-1,kbc,3,n) + q(ibc,jbc  ,kbc,3,n)) / 2
                              w2 = (q(ibc,jbc  ,kbc,4,n) + q(ibc,jbc+1,kbc,4,n)) / 2
                              w1 = (q(ibc,jbc-1,kbc,4,n) + q(ibc,jbc  ,kbc,4,n)) / 2
                           end if
                           dudeta = u2 - u1
                           dvdeta = v2 - v1
                           dwdeta = w2 - w1
                           
                           if (i3.eq.istrt(ic3)) then
                              u2 = (  q(ibc,jbc,kbc,2,n) + q(ibc,jbc,kbc+1,2,n)) / 2
                              u1 = (3*q(ibc,jbc,kbc,2,n) - q(ibc,jbc,kbc+1,2,n)) / 2
                              v2 = (  q(ibc,jbc,kbc,3,n) + q(ibc,jbc,kbc+1,3,n)) / 2
                              v1 = (3*q(ibc,jbc,kbc,3,n) - q(ibc,jbc,kbc+1,3,n)) / 2
                              w2 = (  q(ibc,jbc,kbc,4,n) + q(ibc,jbc,kbc+1,4,n)) / 2
                              w1 = (3*q(ibc,jbc,kbc,4,n) - q(ibc,jbc,kbc+1,4,n)) / 2
                           else if (i3.eq.iend(ic3)) then
                              u2 = (3*q(ibc,jbc,kbc,2,n) - q(ibc,jbc,kbc-1,2,n)) / 2
                              u1 = (  q(ibc,jbc,kbc,2,n) + q(ibc,jbc,kbc-1,2,n)) / 2
                              v2 = (3*q(ibc,jbc,kbc,3,n) - q(ibc,jbc,kbc-1,3,n)) / 2
                              v1 = (  q(ibc,jbc,kbc,3,n) + q(ibc,jbc,kbc-1,3,n)) / 2
                              w2 = (3*q(ibc,jbc,kbc,4,n) - q(ibc,jbc,kbc-1,4,n)) / 2
                              w1 = (  q(ibc,jbc,kbc,4,n) + q(ibc,jbc,kbc-1,4,n)) / 2
                           else
                              u2 = (q(ibc,jbc,kbc  ,2,n) + q(ibc,jbc,kbc+1,2,n)) / 2
                              u1 = (q(ibc,jbc,kbc-1,2,n) + q(ibc,jbc,kbc  ,2,n)) / 2
                              v2 = (q(ibc,jbc,kbc  ,3,n) + q(ibc,jbc,kbc+1,3,n)) / 2
                              v1 = (q(ibc,jbc,kbc-1,3,n) + q(ibc,jbc,kbc  ,3,n)) / 2
                              w2 = (q(ibc,jbc,kbc  ,4,n) + q(ibc,jbc,kbc+1,4,n)) / 2
                              w1 = (q(ibc,jbc,kbc-1,4,n) + q(ibc,jbc,kbc  ,4,n)) / 2
                           end if
                           dudzeta = u2 - u1
                           dvdzeta = v2 - v1
                           dwdzeta = w2 - w1
                           
                           xix   = xideri  (1,im,jm,km,n)
                           xiy   = xideri  (2,im,jm,km,n)
                           xiz   = xideri  (3,im,jm,km,n)
                           etax  = etaderi (1,im,jm,km,n)
                           etay  = etaderi (2,im,jm,km,n)
                           etaz  = etaderi (3,im,jm,km,n)
                           zetax = zetaderi(1,im,jm,km,n)
                           zetay = zetaderi(2,im,jm,km,n)
                           zetaz = zetaderi(3,im,jm,km,n)
                           
                        else if (idir.eq.2) then
                           
                           if (i3.eq.istrt(ic3)) then
                              u2 = (  q(ibc+1,jbc,kbc,2,n) + q(ibc  ,jbc,kbc,2,n))/2
                              u1 = (3*q(ibc  ,jbc,kbc,2,n) - q(ibc+1,jbc,kbc,2,n))/2
                              v2 = (  q(ibc+1,jbc,kbc,3,n) + q(ibc  ,jbc,kbc,3,n))/2
                              v1 = (3*q(ibc  ,jbc,kbc,3,n) - q(ibc+1,jbc,kbc,3,n))/2
                              w2 = (  q(ibc+1,jbc,kbc,4,n) + q(ibc  ,jbc,kbc,4,n))/2
                              w1 = (3*q(ibc  ,jbc,kbc,4,n) - q(ibc+1,jbc,kbc,4,n))/2
                           else if (i3.eq.iend(ic3)) then
                              u2 = (3*q(ibc,jbc,kbc,2,n) - q(ibc-1,jbc,kbc,2,n)) / 2
                              u1 = (  q(ibc,jbc,kbc,2,n) + q(ibc-1,jbc,kbc,2,n)) / 2
                              v2 = (3*q(ibc,jbc,kbc,3,n) - q(ibc-1,jbc,kbc,3,n)) / 2
                              v1 = (  q(ibc,jbc,kbc,3,n) + q(ibc-1,jbc,kbc,3,n)) / 2
                              w2 = (3*q(ibc,jbc,kbc,4,n) - q(ibc-1,jbc,kbc,4,n)) / 2
                              w1 = (  q(ibc,jbc,kbc,4,n) + q(ibc-1,jbc,kbc,4,n)) / 2
                           else
                              u2 = (q(ibc+1,jbc,kbc,2,n) + q(ibc  ,jbc,kbc,2,n)) / 2
                              u1 = (q(ibc  ,jbc,kbc,2,n) + q(ibc-1,jbc,kbc,2,n)) / 2
                              v2 = (q(ibc+1,jbc,kbc,3,n) + q(ibc  ,jbc,kbc,3,n)) / 2
                              v1 = (q(ibc  ,jbc,kbc,3,n) + q(ibc-1,jbc,kbc,3,n)) / 2
                              w2 = (q(ibc+1,jbc,kbc,4,n) + q(ibc  ,jbc,kbc,4,n)) / 2
                              w1 = (q(ibc  ,jbc,kbc,4,n) + q(ibc-1,jbc,kbc,4,n)) / 2
                           end if
                           dudxi = u2 - u1
                           dvdxi = v2 - v1
                           dwdxi = w2 - w1
                           
                           dudeta = (-1)**(1+1/inrout) * &
                                (-8*q(ibc,jbc,kbc,2,n) + &
                                9*q(ibc,jbc+joff1,kbc,2,n) - &
                                q(ibc,jbc+joff2,kbc,2,n) ) / 3
                           dvdeta = (-1)**(1+1/inrout) * &
                                (-8*q(ibc,jbc,kbc,3,n) + &
                                9*q(ibc,jbc+joff1,kbc,3,n) - &
                                q(ibc,jbc+joff2,kbc,3,n) ) / 3
                           dwdeta = (-1)**(1+1/inrout) * &
                                (-8*q(ibc,jbc,kbc,4,n) + &
                                9*q(ibc,jbc+joff1,kbc,4,n) - &
                                q(ibc,jbc+joff2,kbc,4,n) ) / 3
                           
                           if (i2.eq.istrt(ic2)) then
                              u2 = (  q(ibc,jbc,kbc+1,2,n) + q(ibc,jbc,kbc  ,2,n))/2
                              u1 = (3*q(ibc,jbc,kbc  ,2,n) - q(ibc,jbc,kbc+1,2,n))/2
                              v2 = (  q(ibc,jbc,kbc+1,3,n) + q(ibc,jbc,kbc  ,3,n))/2
                              v1 = (3*q(ibc,jbc,kbc  ,3,n) - q(ibc,jbc,kbc+1,3,n))/2
                              w2 = (  q(ibc,jbc,kbc+1,4,n) + q(ibc,jbc,kbc  ,4,n))/2
                              w1 = (3*q(ibc,jbc,kbc  ,4,n) - q(ibc,jbc,kbc+1,4,n))/2
                           else if (i2.eq.iend(ic2)) then
                              u2 = (3*q(ibc,jbc,kbc,2,n) - q(ibc,jbc,kbc-1,2,n)) / 2
                              u1 = (  q(ibc,jbc,kbc,2,n) + q(ibc,jbc,kbc-1,2,n)) / 2
                              v2 = (3*q(ibc,jbc,kbc,3,n) - q(ibc,jbc,kbc-1,3,n)) / 2
                              v1 = (  q(ibc,jbc,kbc,3,n) + q(ibc,jbc,kbc-1,3,n)) / 2
                              w2 = (3*q(ibc,jbc,kbc,4,n) - q(ibc,jbc,kbc-1,4,n)) / 2
                              w1 = (  q(ibc,jbc,kbc,4,n) + q(ibc,jbc,kbc-1,4,n)) / 2
                           else
                              u2 = (q(ibc,jbc,kbc+1,2,n) + q(ibc,jbc,kbc  ,2,n)) / 2
                              u1 = (q(ibc,jbc,kbc  ,2,n) + q(ibc,jbc,kbc-1,2,n)) / 2
                              v2 = (q(ibc,jbc,kbc+1,3,n) + q(ibc,jbc,kbc  ,3,n)) / 2
                              v1 = (q(ibc,jbc,kbc  ,3,n) + q(ibc,jbc,kbc-1,3,n)) / 2
                              w2 = (q(ibc,jbc,kbc+1,4,n) + q(ibc,jbc,kbc  ,4,n)) / 2
                              w1 = (q(ibc,jbc,kbc  ,4,n) + q(ibc,jbc,kbc-1,4,n)) / 2
                           end if
                           dudzeta = u2 - u1
                           dvdzeta = v2 - v1
                           dwdzeta = w2 - w1
                           
                           xix   = xiderj  (1,im,jm,km,n)
                           xiy   = xiderj  (2,im,jm,km,n)
                           xiz   = xiderj  (3,im,jm,km,n)
                           etax  = etaderj (1,im,jm,km,n)
                           etay  = etaderj (2,im,jm,km,n)
                           etaz  = etaderj (3,im,jm,km,n)
                           zetax = zetaderj(1,im,jm,km,n)
                           zetay = zetaderj(2,im,jm,km,n)
                           zetaz = zetaderj(3,im,jm,km,n)
                           
                        else if (idir.eq.3) then
                           if (i2.eq.istrt(ic2)) then
                              u2 = (  q(ibc,jbc,kbc,2,n) + q(ibc+1,jbc,kbc,2,n)) / 2
                              u1 = (3*q(ibc,jbc,kbc,2,n) - q(ibc+1,jbc,kbc,2,n)) / 2
                              v2 = (  q(ibc,jbc,kbc,3,n) + q(ibc+1,jbc,kbc,3,n)) / 2
                              v1 = (3*q(ibc,jbc,kbc,3,n) - q(ibc+1,jbc,kbc,3,n)) / 2
                              w2 = (  q(ibc,jbc,kbc,4,n) + q(ibc+1,jbc,kbc,4,n)) / 2
                              w1 = (3*q(ibc,jbc,kbc,4,n) - q(ibc+1,jbc,kbc,4,n)) / 2
                           else if (i2.eq.iend(ic2)) then
                              u2 = (3*q(ibc,jbc,kbc,2,n) - q(ibc-1,jbc,kbc,2,n)) / 2
                              u1 = (  q(ibc,jbc,kbc,2,n) + q(ibc-1,jbc,kbc,2,n)) / 2
                              v2 = (3*q(ibc,jbc,kbc,3,n) - q(ibc-1,jbc,kbc,3,n)) / 2
                              v1 = (  q(ibc,jbc,kbc,3,n) + q(ibc-1,jbc,kbc,3,n)) / 2
                              w2 = (3*q(ibc,jbc,kbc,4,n) - q(ibc-1,jbc,kbc,4,n)) / 2
                              w1 = (  q(ibc,jbc,kbc,4,n) + q(ibc-1,jbc,kbc,4,n)) / 2
                           else
                              u2 = (q(ibc  ,jbc,kbc,2,n) + q(ibc+1,jbc,kbc,2,n)) / 2
                              u1 = (q(ibc-1,jbc,kbc,2,n) + q(ibc  ,jbc,kbc,2,n)) / 2
                              v2 = (q(ibc  ,jbc,kbc,3,n) + q(ibc+1,jbc,kbc,3,n)) / 2
                              v1 = (q(ibc-1,jbc,kbc,3,n) + q(ibc  ,jbc,kbc,3,n)) / 2
                              w2 = (q(ibc  ,jbc,kbc,4,n) + q(ibc+1,jbc,kbc,4,n)) / 2
                              w1 = (q(ibc-1,jbc,kbc,4,n) + q(ibc  ,jbc,kbc,4,n)) / 2
                           end if
                           dudxi = u2 - u1
                           dvdxi = v2 - v1
                           dwdxi = w2 - w1
                           
                           if (i3.eq.istrt(ic3)) then
                              u2 = (  q(ibc,jbc,kbc,2,n) + q(ibc,jbc+1,kbc,2,n)) / 2
                              u1 = (3*q(ibc,jbc,kbc,2,n) - q(ibc,jbc+1,kbc,2,n)) / 2
                              v2 = (  q(ibc,jbc,kbc,3,n) + q(ibc,jbc+1,kbc,3,n)) / 2
                              v1 = (3*q(ibc,jbc,kbc,3,n) - q(ibc,jbc+1,kbc,3,n)) / 2
                              w2 = (  q(ibc,jbc,kbc,4,n) + q(ibc,jbc+1,kbc,4,n)) / 2
                              w1 = (3*q(ibc,jbc,kbc,4,n) - q(ibc,jbc+1,kbc,4,n)) / 2
                           else if (i3.eq.iend(ic3)) then
                              u2 = (3*q(ibc,jbc,kbc,2,n) - q(ibc,jbc-1,kbc,2,n)) / 2
                              u1 = (  q(ibc,jbc,kbc,2,n) + q(ibc,jbc-1,kbc,2,n)) / 2
                              v2 = (3*q(ibc,jbc,kbc,3,n) - q(ibc,jbc-1,kbc,3,n)) / 2
                              v1 = (  q(ibc,jbc,kbc,3,n) + q(ibc,jbc-1,kbc,3,n)) / 2
                              w2 = (3*q(ibc,jbc,kbc,4,n) - q(ibc,jbc-1,kbc,4,n)) / 2
                              w1 = (  q(ibc,jbc,kbc,4,n) + q(ibc,jbc-1,kbc,4,n)) / 2
                           else
                              u2 = (q(ibc,jbc  ,kbc,2,n) + q(ibc,jbc+1,kbc,2,n)) / 2
                              u1 = (q(ibc,jbc-1,kbc,2,n) + q(ibc,jbc  ,kbc,2,n)) / 2
                              v2 = (q(ibc,jbc  ,kbc,3,n) + q(ibc,jbc+1,kbc,3,n)) / 2
                              v1 = (q(ibc,jbc-1,kbc,3,n) + q(ibc,jbc  ,kbc,3,n)) / 2
                              w2 = (q(ibc,jbc  ,kbc,4,n) + q(ibc,jbc+1,kbc,4,n)) / 2
                              w1 = (q(ibc,jbc-1,kbc,4,n) + q(ibc,jbc  ,kbc,4,n)) / 2
                           end if
                           dudeta = u2 - u1
                           dvdeta = v2 - v1
                           dwdeta = w2 - w1
                           
                           dudzeta = (-1)**(1+1/inrout) * &
                                (-8*q(ibc,jbc,kbc,2,n) + &
                                9*q(ibc,jbc,kbc+koff1,2,n) - &
                                q(ibc,jbc,kbc+koff2,2,n) ) / 3
                           dvdzeta = (-1)**(1+1/inrout) * &
                                (-8*q(ibc,jbc,kbc,3,n) + &
                                9*q(ibc,jbc,kbc+koff1,3,n) - &
                                q(ibc,jbc,kbc+koff2,3,n) ) / 3
                           dwdzeta = (-1)**(1+1/inrout) * &
                                (-8*q(ibc,jbc,kbc,4,n) + &
                                9*q(ibc,jbc,kbc+koff1,4,n) - &
                                q(ibc,jbc,kbc+koff2,4,n) ) / 3
                           
                           xix   = xiderk  (1,im,jm,km,n)
                           xiy   = xiderk  (2,im,jm,km,n)
                           xiz   = xiderk  (3,im,jm,km,n)
                           etax  = etaderk (1,im,jm,km,n)
                           etay  = etaderk (2,im,jm,km,n)
                           etaz  = etaderk (3,im,jm,km,n)
                           zetax = zetaderk(1,im,jm,km,n)
                           zetay = zetaderk(2,im,jm,km,n)
                           zetaz = zetaderk(3,im,jm,km,n)
                           
                        end if
                        
                        dudx = dudxi*xix + dudeta*etax + dudzeta*zetax
                        dudy = dudxi*xiy + dudeta*etay + dudzeta*zetay
                        dudz = dudxi*xiz + dudeta*etaz + dudzeta*zetaz
                        dvdx = dvdxi*xix + dvdeta*etax + dvdzeta*zetax
                        dvdy = dvdxi*xiy + dvdeta*etay + dvdzeta*zetay
                        dvdz = dvdxi*xiz + dvdeta*etaz + dvdzeta*zetaz
                        dwdx = dwdxi*xix + dwdeta*etax + dwdzeta*zetax
                        dwdy = dwdxi*xiy + dwdeta*etay + dwdzeta*zetay
                        dwdz = dwdxi*xiz + dwdeta*etaz + dwdzeta*zetaz
                        
                        divv = dudx + dvdy + dwdz
                        
                        txx = 2*muw * (dudx - divv/3)
                        txy =   muw * (dudy + dvdx)
                        txz =   muw * (dudz + dwdx)
                        tyy = 2*muw * (dvdy - divv/3)
                        tyz =   muw * (dvdz + dwdy)
                        tzz = 2*muw * (dwdz - divv/3)
                        
                        tauwx = txx*kx + txy*ky + txz*kz
                        tauwy = txy*kx + tyy*ky + tyz*kz
                        tauwz = txz*kx + tyz*ky + tzz*kz
                        
                        tauwpx = tauwx - (tauwx*kx+tauwy*ky+tauwz*kz)*kx
                        tauwpy = tauwy - (tauwx*kx+tauwy*ky+tauwz*kz)*ky
                        tauwpz = tauwz - (tauwx*kx+tauwy*ky+tauwz*kz)*kz
                        
                        tauw   = dsqrt(tauwpx**2+tauwpy**2+tauwpz**2)
                        utau   = dsqrt(tauw/rhow)
                        
                        cf(ibc,jbc,kbc) = 2*tauw / (reyno*machfs)
                        yp(ibc,jbc,kbc) = dsqrt(reyno/machfs)* &
                             rhow*dn1*utau/muw
                        
                     end if
                     
                  end do
               end do
               
               
               if(.not.viscous) then
                 
                  cf = 0.d0
                  yp = 0.d0
                  
               end if
!------------extract/calculate surface coefficients - END
               
!------------write cgns surface patch
               call cg_sol_write_f(fid(n),basenum,blocknum,'SurfaceSolutions', &
                    Vertex,solnum,ierr)
!     AJ     &            CellCenter,solnum,ierr)
               if(ierr .ne. CG_OK) then
                  write(*,*) 'cg_sol_write_f error'
                  call cg_error_print_f()
               end if
               localcoord = x(ibc_s:ibc_e,jbc_s:jbc_e,kbc_s:kbc_e,n)
               call cg_coord_write_f(fid(n),basenum,blocknum,RealDouble, &
                    'CoordinateX',localcoord,xnum,ierr)
               if(ierr .ne. CG_OK) then
                  write(*,*) 'X cg_coord_write_f error'
                  call cg_error_print_f()
               end if
               localcoord = y(ibc_s:ibc_e,jbc_s:jbc_e,kbc_s:kbc_e,n)      
               call cg_coord_write_f(fid(n),basenum,blocknum,RealDouble, &
                    'CoordinateY',localcoord,ynum,ierr)
               if(ierr .ne. CG_OK) then
                  write(*,*) 'Y cg_coord_write_f error'
                  call cg_error_print_f()
               end if
               localcoord = z(ibc_s:ibc_e,jbc_s:jbc_e,kbc_s:kbc_e,n)
               call cg_coord_write_f(fid(n),basenum,blocknum,RealDouble, &
                    'CoordinateZ',localcoord,znum,ierr)
               if(ierr .ne. CG_OK) then
                  write(*,*) 'Z cg_coord_write_f error'
                  call cg_error_print_f()
               end if
               call cg_field_write_f(fid(n),basenum,blocknum,solnum, &
                    RealDouble,'StaticPressureCoefficient', &
                    cp(ibc_s:ibc_e,jbc_s:jbc_e,kbc_s:kbc_e),indexnum,ierr)
               if(ierr .ne. CG_OK) then
                  write(*,*) 'cg_field_write_f StaticPressureCoefficient error'
                  call cg_error_print_f()
               end if
               call cg_field_write_f(fid(n),basenum,blocknum,solnum, &
                    RealDouble,'SkinFrictionCoefficient', &
                    cf(ibc_s:ibc_e,jbc_s:jbc_e,kbc_s:kbc_e),indexnum,ierr)
               if(ierr .ne. CG_OK) then
                  write(*,*) 'cg_field_write_f SkinFrictionCoefficient error'
                  call cg_error_print_f()
               end if
               call cg_field_write_f(fid(n),basenum,blocknum,solnum, &
                    RealDouble,'YPlus', &
                    yp(ibc_s:ibc_e,jbc_s:jbc_e,kbc_s:kbc_e),indexnum,ierr)
               if(ierr .ne. CG_OK) then
                  write(*,*) 'cg_field_write_f YPlus error'
                  call cg_error_print_f()
               end if              
               
               deallocate(cp,cf,yp,localcoord)
               
               surfacenumber = surfacenumber + 1
               
            end if
            
!--------end loop (ibc1=1,nbcs)
         end do
         
!---- loop over harmonics - END
      end do


 14   format(3e22.14)
 10   format(3e22.14)

      return
      end
