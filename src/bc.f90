!-----------------------------------------------------------------------
      subroutine bc(q,si,sj,sk,dist,nl,bctopo)
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) iblock,nl,imax,jmax,kmax,ijkmax,iblk,iq,iimt,idist
      integer(kind=cosa_int) bctopo(*)
      real(kind=cosa_real) q(*),si(*),sj(*),sk(*),dist(*)

      do iblock = 1,mynblocks
        imax   = i_imax     (iblock,nl)
        jmax   = j_jmax     (iblock,nl)
        kmax   = k_kmax     (iblock,nl)
        ijkmax = ijk_ijkmax (iblock,nl)
        iblk   = 1 + off_bct(iblock,nl)
        iq     = 1 + off_p3 (iblock,nl) * npde * dim5
        iimt   = 1 + off_p2 (iblock,nl) * lmet * dim5h
        idist  = 1 + off_p1 (iblock,nl)
        call bbc(q(iq),si(iimt),sj(iimt),sk(iimt),dist(idist), &
                 nl,bctopo(iblk),nbcs(iblock),imax,jmax,kmax, &
                 ijkmax,npde,nharms,lmet)
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine bbc(q,si,sj,sk,dist,nl,bctopo,nbcs,imax,jmax,kmax, &
                     ijkmax,npde,nharms,lmet)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,ijkmax,nbcs,npde,nharms,lmet
      integer(kind=cosa_int) nl,ibc
      integer(kind=cosa_int) bctopo(10,nbcs)
      real(kind=cosa_real) &
          q    (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
          dist ( 0:imax  , 0:jmax  , 0:kmax), &
          si   (lmet, 0:imax+1, 0:jmax+1, 0:kmax+1,0:2*nharms*hbmove), &
          sj   (lmet, 0:imax+1, 0:jmax+1, 0:kmax+1,0:2*nharms*hbmove), &
          sk   (lmet, 0:imax+1, 0:jmax+1, 0:kmax+1,0:2*nharms*hbmove)

!-----------------------------------------------------------------------
!     bc type                                                       mask
!
!     inviscid wall (all vars are extrapolated except normal vel.)    14
!     viscous wall                                                    15
!     subsonic inflow                                                 21
!     subsonic or supersonic free-stream                              30
!     extrapolation                                                    4
!     symmetry wrt xy-, or xz- or yz-plane                            80
!     subsonic outflow                                                91
!-----------------------------------------------------------------------
      do ibc = 1,nbcs

!------ planar symmetry (bctyp=80)
        if (bctopo(1,ibc).eq.80) then

          call bc_sym(q,bctopo(1,ibc),imax,jmax,kmax,npde,nharms)

!------ extrapolation boundary condition (bctyp=4)
        else if (bctopo(1,ibc).eq.4) then

          call bc_extrapol(q,bctopo(1,ibc),imax,jmax,kmax,npde,nharms)

!------ subsonic/supersonic free stream boundary condition (jbctyp=30)
!       From Jameson and Baker, AIAA 83-1929
        else if (bctopo(1,ibc).eq.30) then

          call bc_far0(q,si,sj,sk,bctopo(1,ibc),imax,jmax,kmax,npde, &
                       nharms,lmet)

        else if (bctopo(1,ibc).eq.31) then

          call bc_far1(nl,q,si,sj,sk,bctopo(1,ibc),imax,jmax,kmax,npde, &
                       nharms,lmet)

!------ inviscid wall boundary conditions (14)
!         It extrapolates density, pressure and tangential
!         flow velocity from interior nodes.

        else if (bctopo(1,ibc).eq.14) then

          call bc_walinv_14(q,si,sj,sk,bctopo(1,ibc),imax,jmax,kmax, &
                            npde,nharms,lmet)

!------ viscous wall (jbctyp=15)
        else if (bctopo(1,ibc)/100.eq.15) then

          call bc_walvis_15(q,si,sj,sk,dist,bctopo(1,ibc),imax,jmax, &
                            kmax,npde,nharms,lmet)

!------ subsonic inlet 
        else if (bctopo(1,ibc).eq.21) then

          call bc_insub(q,si,sj,sk,bctopo(1,ibc),imax,jmax,kmax, &
                        npde,nharms,lmet,nl,ijkmax)

!------ subsonic outlet
        else if (bctopo(1,ibc).eq.91) then

          call bc_outsub(q,si,sj,sk,bctopo(1,ibc),imax,jmax,kmax, &
                         npde,nharms,lmet,nl,ijkmax)

        end if

      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine bc_far0(q,si,sj,sk,bctopo,imax,jmax,kmax,npde,nharms, &
                         lmet)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,npde,nharms,lmet
      integer(kind=cosa_int) bctopo(10)
      integer(kind=cosa_int) n,nh,i1,i2,i3,ijkmax(3),idir,istrt(3),iend(3),bctyp, &
        inrout,ibcpt,ibc,jbc,kbc,ibcpt2,ibc2,jbc2,kbc2,ibcn,in,jn,kn, &
        ibcm,im,jm,km,ic1,ic2,ic3
      real (kind=cosa_real) &
          q    (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
          si   (lmet, 0:imax+1, 0:jmax+1, 0:kmax+1,0:2*nharms*hbmove), &
          sj   (lmet, 0:imax+1, 0:jmax+1, 0:kmax+1,0:2*nharms*hbmove), &
          sk   (lmet, 0:imax+1, 0:jmax+1, 0:kmax+1,0:2*nharms*hbmove)
      real (kind=cosa_real) qfst(5,0:2*mharms),angle,rotaz(2,2),ua,va,rhobc, &
           sinf,sext,sx,sy,sz,qext,qinf,rhoext,pext,a,aext,rext,rinf, &
           qdotn,mach,p,sgnm

      if (aircraft.or.vawt) then
        do n = 0,2*nharms
          qfst(1,n) = rhoinf
          qfst(2,n) = machfs * dcos(alpha) * dcos(beta)
          qfst(3,n) = machfs * dsin(alpha) * dcos(beta)
          qfst(4,n) = machfs * dsin(beta)
          qfst(5,n) = pinf
        end do
      else if (hawt) then
        if ((.not.unsteady).or.(.not.relframe)) then
          do n = 0,2*nharms
            qfst(1,n) = rhoinf
            qfst(2,n) = machfs * dsin(beta)
            qfst(3,n) = machfs * dsin(alpha) * dcos(beta)
            qfst(4,n) = machfs * dcos(alpha) * dcos(beta)
            qfst(5,n) = pinf
          end do
        else if (relframe) then
          ua =  machfs * dsin(beta)
          va =  machfs * dsin(alpha) * dcos(beta)
          do n = 0,2*nharms
            angle      = -omegas * zeit(n)
            rotaz(1,1) =  cos(angle)
            rotaz(1,2) = -sin(angle)
            rotaz(2,1) =  sin(angle)
            rotaz(2,2) =  cos(angle)
            qfst(1,n) = rhoinf
            qfst(2,n) = (ua*rotaz(1,1) + va*rotaz(1,2))
            qfst(3,n) = (ua*rotaz(2,1) + va*rotaz(2,2))
            qfst(4,n) = machfs * dcos(alpha) * dcos(beta)
            qfst(5,n) = pinf
          end do
        end if
      end if

      sinf    = pinf / rhoinf**gamma

!     store imax, jmax in ijkmax to enable boundary data location

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

!az         write(*,*) 'IDIR',idir
!az         write(*,*) bctopo(1),bctopo(2),bctopo(3),bctopo(4), &
!az                    bctopo(5),bctopo(6),bctopo(7)

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
        ibcm   =  1
        sgnm   = - 1.d0
      else
        ibcpt  = ijkmax(idir)
        ibcpt2 = ijkmax(idir) + 1
        ibcn   = ijkmax(idir) - 1
        ibcm   = ijkmax(idir)
        sgnm   = 1.d0
      end if

      ic1 = cyc (idir, 1)
      ic2 = cyc (idir, 2)
      ic3 = cyc (idir, 3)

      do n = 0,2*nharms
        nh = n*hbmove

        do i3 = istrt(ic3),iend(ic3)
          do i2 = istrt(ic2),iend(ic2)

            ibc  = ibcpt *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
            jbc  = ibcpt *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
            kbc  = ibcpt *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

            ibc2 = ibcpt2*krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
            jbc2 = ibcpt2*krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
            kbc2 = ibcpt2*krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

            in   = ibcn  *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
            jn   = ibcn  *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
            kn   = ibcn  *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

            im   = ibcm  *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
            jm   = ibcm  *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
            km   = ibcm  *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

            if (lomach) then
!---------- low speed treatment: freeze bc as preconditioned bc has yet to
!                                be implemented

              q(ibc ,jbc ,kbc ,1,n) = qfst(1,n)
              q(ibc ,jbc ,kbc ,2,n) = qfst(2,n)
              q(ibc ,jbc ,kbc ,3,n) = qfst(3,n)
              q(ibc ,jbc ,kbc ,4,n) = qfst(4,n)
              q(ibc ,jbc ,kbc ,5,n) = qfst(5,n)
              q(ibc2,jbc2,kbc2,1,n) = q(ibc,jbc,kbc,1,n)
              q(ibc2,jbc2,kbc2,2,n) = q(ibc,jbc,kbc,2,n)
              q(ibc2,jbc2,kbc2,3,n) = q(ibc,jbc,kbc,3,n)
              q(ibc2,jbc2,kbc2,4,n) = q(ibc,jbc,kbc,4,n)
              q(ibc2,jbc2,kbc2,5,n) = q(ibc,jbc,kbc,5,n)

            else
!---------- high speed treatment.

              if (idir.eq.1) then
                sx    = si(1,im,jm,km,nh) * sgnm
                sy    = si(2,im,jm,km,nh) * sgnm
                sz    = si(3,im,jm,km,nh) * sgnm
              else if (idir.eq.2) then
                sx    = sj(1,im,jm,km,nh) * sgnm
                sy    = sj(2,im,jm,km,nh) * sgnm
                sz    = sj(3,im,jm,km,nh) * sgnm
              else if (idir.eq.3) then
                sx    = sk(1,im,jm,km,nh) * sgnm
                sy    = sk(2,im,jm,km,nh) * sgnm
                sz    = sk(3,im,jm,km,nh) * sgnm
              end if

              qext  = q(in,jn,kn,2,n) * sx + q(in,jn,kn,3,n) * sy + &
                      q(in,jn,kn,4,n) * sz
              qinf  = qfst(2,n) * sx + qfst(3,n) * sy + qfst(4,n) * sz
              rhoext = q(in,jn,kn,1,n)
              pext = q(in,jn,kn,5,n)

!             One Dimensional Riemann invariants. The sign of U dot n 
!             handles the sign change in the Riemann invariants
!
              aext  = sqrt (gamma * pext / rhoext)
              rext  = qext + 2 * aext / (gamma-1)
              rinf  = qinf - 2 * ainf / (gamma-1)
              qdotn = (rext + rinf) / 2
              a     = (gamma-1) * (rext - rinf) / 4
              mach = qdotn / a
              sext  = pext / rhoext**gamma

              if (qdotn.le.0.d0) then

!               SUPERSONIC INFLOW: all data specified
                if (abs(mach).ge. 1.0d0) then
                  q(ibc ,jbc ,kbc ,1,n) = qfst(1,n)
                  q(ibc ,jbc ,kbc ,2,n) = qfst(2,n)
                  q(ibc ,jbc ,kbc ,3,n) = qfst(3,n)
                  q(ibc ,jbc ,kbc ,4,n) = qfst(4,n)
                  q(ibc ,jbc ,kbc ,5,n) = qfst(5,n)
                  q(ibc2,jbc2,kbc2,1,n) = q(ibc,jbc,kbc,1,n)
                  q(ibc2,jbc2,kbc2,2,n) = q(ibc,jbc,kbc,2,n)
                  q(ibc2,jbc2,kbc2,3,n) = q(ibc,jbc,kbc,3,n)
                  q(ibc2,jbc2,kbc2,4,n) = q(ibc,jbc,kbc,4,n)
                  q(ibc2,jbc2,kbc2,5,n) = q(ibc,jbc,kbc,5,n)
!               SUBSONIC INFLOW: one-dimensional Riemann Problem
                else if (abs(mach).lt.1.0d0) then
                  rhobc  = (a*a / (gamma*sinf))**(1.d0/(gamma-1))
                  p = rhobc * a * a / gamma
                  q(ibc ,jbc ,kbc ,1,n) = rhobc
                  q(ibc ,jbc ,kbc ,2,n) = qfst(2,n) + (qdotn-qinf)*sx
                  q(ibc ,jbc ,kbc ,3,n) = qfst(3,n) + (qdotn-qinf)*sy
                  q(ibc ,jbc ,kbc ,4,n) = qfst(4,n) + (qdotn-qinf)*sz
                  q(ibc ,jbc ,kbc ,5,n) = p
                  q(ibc2,jbc2,kbc2,1,n) = q(ibc,jbc,kbc ,1,n)
                  q(ibc2,jbc2,kbc2,2,n) = q(ibc,jbc,kbc ,2,n)
                  q(ibc2,jbc2,kbc2,3,n) = q(ibc,jbc,kbc ,3,n)
                  q(ibc2,jbc2,kbc2,4,n) = q(ibc,jbc,kbc ,4,n)
                  q(ibc2,jbc2,kbc2,5,n) = q(ibc,jbc,kbc ,5,n)
                end if

              else if (qdotn.gt.0.d0) then

!               SUPERSONIC OUTFLOW: all data extrapolated
                if (abs(mach).ge. 1.0d0) then
                  q(ibc ,jbc ,kbc ,1,n) = q(in ,jn ,kn ,1,n)
                  q(ibc ,jbc ,kbc ,2,n) = q(in ,jn ,kn ,2,n)
                  q(ibc ,jbc ,kbc ,3,n) = q(in ,jn ,kn ,3,n)
                  q(ibc ,jbc ,kbc ,4,n) = q(in ,jn ,kn ,4,n)
                  q(ibc ,jbc ,kbc ,5,n) = q(in ,jn ,kn ,5,n)
                  q(ibc2,jbc2,kbc2,1,n) = q(ibc,jbc,kbc,1,n)
                  q(ibc2,jbc2,kbc2,2,n) = q(ibc,jbc,kbc,2,n)
                  q(ibc2,jbc2,kbc2,3,n) = q(ibc,jbc,kbc,3,n)
                  q(ibc2,jbc2,kbc2,4,n) = q(ibc,jbc,kbc,4,n)
                  q(ibc2,jbc2,kbc2,5,n) = q(ibc,jbc,kbc,5,n)
!               SUBSONIC OUTFLOW: one-dimensional Riemann Problem
                else if (abs(mach).lt.1.0d0) then
                  rhobc  = (a*a / (gamma*sext))**(1.d0/(gamma-1))
                  p = rhobc * a * a / gamma
                  q(ibc ,jbc ,kbc ,1,n) = rhobc
                  q(ibc ,jbc ,kbc ,2,n) = &
                                       q(in,jn,kn,2,n) + (qdotn-qext)*sx
                  q(ibc ,jbc ,kbc ,3,n) = &
                                       q(in,jn,kn,3,n) + (qdotn-qext)*sy
                  q(ibc ,jbc ,kbc ,4,n) = &
                                       q(in,jn,kn,4,n) + (qdotn-qext)*sz
                  q(ibc ,jbc ,kbc ,5,n) = p
                  q(ibc2,jbc2,kbc2,1,n) = q(ibc,jbc,kbc ,1,n)
                  q(ibc2,jbc2,kbc2,2,n) = q(ibc,jbc,kbc ,2,n)
                  q(ibc2,jbc2,kbc2,3,n) = q(ibc,jbc,kbc ,3,n)
                  q(ibc2,jbc2,kbc2,4,n) = q(ibc,jbc,kbc ,4,n)
                  q(ibc2,jbc2,kbc2,5,n) = q(ibc,jbc,kbc ,5,n)
                end if

              end if

            end if

          end do
        end do

      end do

      if (kom.or.kom_bsl.or.kom_sst) then

        do n = 0,2*nharms
          nh = n*hbmove

          do i3 = istrt(ic3),iend(ic3)
            do i2 = istrt(ic2),iend(ic2)

              ibc  = ibcpt *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
              jbc  = ibcpt *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
              kbc  = ibcpt *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

              ibc2 = ibcpt2*krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
              jbc2 = ibcpt2*krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
              kbc2 = ibcpt2*krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

              in   = ibcn  *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
              jn   = ibcn  *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
              kn   = ibcn  *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

              im   = ibcm  *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
              jm   = ibcm  *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
              km   = ibcm  *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

              if (lomach) then
!------------ low speed treatment: freeze bc as preconditioned bc has yet
!                                  to be implemented

                q(ibc ,jbc ,kbc ,6,n) = tkefar
                q(ibc ,jbc ,kbc ,7,n) = omefar
                q(ibc2,jbc2,kbc2,6,n) = q(ibc,jbc,kbc,6,n)
                q(ibc2,jbc2,kbc2,7,n) = q(ibc,jbc,kbc,7,n)

              else
!------------ high speed treatment.

                if (idir.eq.1) then
                  sx    = si(1,im,jm,km,nh) * sgnm
                  sy    = si(2,im,jm,km,nh) * sgnm
                  sz    = si(3,im,jm,km,nh) * sgnm
                else if (idir.eq.2) then
                  sx    = sj(1,im,jm,km,nh) * sgnm
                  sy    = sj(2,im,jm,km,nh) * sgnm
                  sz    = sj(3,im,jm,km,nh) * sgnm
                else if (idir.eq.3) then
                  sx    = sk(1,im,jm,km,nh) * sgnm
                  sy    = sk(2,im,jm,km,nh) * sgnm
                  sz    = sk(3,im,jm,km,nh) * sgnm
                end if

                qext  = q(in,jn,kn,2,n) * sx + q(in,jn,kn,3,n) * sy + &
                        q(in,jn,kn,4,n) * sz

                qinf  = qfst(2,n) * sx + qfst(3,n) * sy + qfst(4,n) * sz
                rhoext = q(in,jn,kn,1,n)
                pext = q(in,jn,kn,5,n)

!               One Dimensional Riemann invariants. The sign of U dot n 
!               handles the sign change in the Riemann invariants
!
                aext  = sqrt (gamma * pext / rhoext)
                rext  = qext + 2 * aext / (gamma-1)
                rinf  = qinf - 2 * ainf / (gamma-1)
                qdotn = (rext + rinf) / 2
                a     = (gamma-1) * (rext - rinf) / 4
                mach = qdotn / a
                sext  = pext / rhoext**gamma

                if (qdotn.le.0.d0) then

!                 SUPERSONIC OR SUBSONIC INFLOW: all data specified
                  q(ibc ,jbc ,kbc ,6,n) = tkefar
                  q(ibc ,jbc ,kbc ,7,n) = omefar
                  q(ibc2,jbc2,kbc2,6,n) = q(ibc,jbc,kbc,6,n)
                  q(ibc2,jbc2,kbc2,7,n) = q(ibc,jbc,kbc,7,n)

                else if (qdotn.gt.0.d0) then

!                 SUPERSONIC OR SUBSONIC OUTFLOW: all data extrapolated
!msc              q(ibc ,jbc ,kbc ,6,n) = max(q(in,jn,kn,6,n),qmin(6))
!msc              q(ibc ,jbc ,kbc ,7,n) = max(q(in,jn,kn,7,n),qmin(7))
                  q(ibc ,jbc ,kbc ,6,n) = q(in ,jn ,kn ,6,n)
                  q(ibc ,jbc ,kbc ,7,n) = q(in ,jn ,kn ,7,n)
                  q(ibc2,jbc2,kbc2,6,n) = q(ibc,jbc,kbc,6,n)
                  q(ibc2,jbc2,kbc2,7,n) = q(ibc,jbc,kbc,7,n)

                end if

              end if

            end do
          end do

        end do

      end if

      return
      end

!-----------------------------------------------------------------------
      subroutine bc_far1(nl,q,si,sj,sk,bctopo,imax,jmax,kmax, &
                         npde,nharms,lmet)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,npde,nharms,lmet,nl
      integer(kind=cosa_int) bctopo(10)
      integer(kind=cosa_int) n,nh,i1,i2,i3,ijkmax(3),idir,istrt(3),iend(3),bctyp, &
        inrout,ibcpt,ibc,jbc,kbc,ibcpt2,ibc2,jbc2,kbc2,ibcn,in,jn,kn, &
        ibcm,im,jm,km,ic1,ic2,ic3
      real (kind=cosa_real) &
          q    (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
          si   (lmet, 0:imax+1, 0:jmax+1, 0:kmax+1,0:2*nharms*hbmove), &
          sj   (lmet, 0:imax+1, 0:jmax+1, 0:kmax+1,0:2*nharms*hbmove), &
          sk   (lmet, 0:imax+1, 0:jmax+1, 0:kmax+1,0:2*nharms*hbmove)
      real (kind=cosa_real) qfst(5,0:2*mharms),angle,rotaz(2,2),ua,va, &
           sx,sy,sz,mach,qext,qinf,rhobc,rhoext,pext,sgnm, &
           q2,a2,a,un,m2,ump2,mp2s,am,apr,unb,pb

      if (aircraft.or.vawt) then
        do n = 0,2*nharms
          qfst(1,n) = rhoinf
          qfst(2,n) = machfs * dcos(alpha) * dcos(beta)
          qfst(3,n) = machfs * dsin(alpha) * dcos(beta)
          qfst(4,n) = machfs * dsin(beta)
          qfst(5,n) = pinf
        end do
      else if (hawt) then
        if ((.not.unsteady).or.(.not.relframe)) then
          do n = 0,2*nharms
            qfst(1,n) = rhoinf
            qfst(2,n) = machfs * dsin(beta)
            qfst(3,n) = machfs * dsin(alpha) * dcos(beta)
            qfst(4,n) = machfs * dcos(alpha) * dcos(beta)
            qfst(5,n) = pinf
          end do
        else if (relframe) then
          ua =  machfs * dsin(beta)
          va =  machfs * dsin(alpha) * dcos(beta)
          do n = 0,2*nharms
            angle      = -omegas * zeit(n)
            rotaz(1,1) =  cos(angle)
            rotaz(1,2) = -sin(angle)
            rotaz(2,1) =  sin(angle)
            rotaz(2,2) =  cos(angle)
            qfst(1,n) = rhoinf
            qfst(2,n) = (ua*rotaz(1,1) + va*rotaz(1,2))
            qfst(3,n) = (ua*rotaz(2,1) + va*rotaz(2,2))
            qfst(4,n) = machfs * dcos(alpha) * dcos(beta)
            qfst(5,n) = pinf
          end do
        end if
      end if

!     store imax, jmax in ijkmax to enable boundary data location

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

!az         write(*,*) 'IDIR',idir
!az         write(*,*) bctopo(1),bctopo(2),bctopo(3),bctopo(4), &
!az                    bctopo(5),bctopo(6),bctopo(7)

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
        ibcm   =  1
        sgnm   = - 1.d0
      else
        ibcpt  = ijkmax(idir)
        ibcpt2 = ijkmax(idir) + 1
        ibcn   = ijkmax(idir) - 1
        ibcm   = ijkmax(idir)
        sgnm   = 1.d0
      end if

      ic1 = cyc (idir, 1)
      ic2 = cyc (idir, 2)
      ic3 = cyc (idir, 3)

      do n = 0,2*nharms
        nh = n*hbmove

        do i3 = istrt(ic3),iend(ic3)
          do i2 = istrt(ic2),iend(ic2)

            ibc  = ibcpt *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
            jbc  = ibcpt *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
            kbc  = ibcpt *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

            ibc2 = ibcpt2*krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
            jbc2 = ibcpt2*krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
            kbc2 = ibcpt2*krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

            in   = ibcn  *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
            jn   = ibcn  *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
            kn   = ibcn  *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

            im   = ibcm  *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
            jm   = ibcm  *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
            km   = ibcm  *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

            if (idir.eq.1) then
              sx    = si(1,im,jm,km,nh) * sgnm
              sy    = si(2,im,jm,km,nh) * sgnm
              sz    = si(3,im,jm,km,nh) * sgnm
            else if (idir.eq.2) then
              sx    = sj(1,im,jm,km,nh) * sgnm
              sy    = sj(2,im,jm,km,nh) * sgnm
              sz    = sj(3,im,jm,km,nh) * sgnm
            else if (idir.eq.3) then
              sx    = sk(1,im,jm,km,nh) * sgnm
              sy    = sk(2,im,jm,km,nh) * sgnm
              sz    = sk(3,im,jm,km,nh) * sgnm
            end if

            if (lomach) then

!------------ LSP treatment: 

              qext = q(in,jn,kn,2,n) * sx + q(in,jn,kn,3,n) * sy + &
                     q(in,jn,kn,4,n) * sz
              qinf = qfst(2,n) * sx + qfst(3,n) * sy + qfst(4,n) * sz
              pext = q(in,jn,kn,5,n)
              pinf = qfst(5,n)

              rhoext = q(in,jn,kn,1,n)
              q2   = q(in,jn,kn,2,n)**2 + q(in,jn,kn,3,n)**2 + &
                     q(in,jn,kn,4,n)**2
              a2   = gamma*pext/rhoext
              a    = sqrt(a2)
              un   = qext

!------------ LSP cut-off, and artificial sound speed apr
              m2   = q2/a2
              ump2 = (uprv / a)**2
!tmp          mp2s = dmin1(1.d0, &
!tmp                 dmax1(m2,cutoff_pgr(in,jn,n), &
!tmp                 cutoff_vis(in,jn,n),(1-mixpi)*ump2,epsp(nl)**2))
              mp2s = dmin1(1.d0, &
                     dmax1(m2,(1-mixpi)*ump2,epsp(nl)**2))
              am   = (1-mp2s) / 2
              apr  = sqrt(am**2 * un**2+a2*mp2s)

              unb  = (qext + qinf)/2 + (pext - pinf) / (2*rhoext*apr) + &
                     (qext - qinf)/2 * am*un/apr
              pb   = &
               (pext  + pinf )/2 + rhoext*mp2s*a2/apr * (qext - qinf)/2 - &
               (pext  - pinf )/2 * am*un/apr

            else
!---------- high speed treatment.

              qext  = q(in,jn,kn,2,n) * sx + q(in,jn,kn,3,n) * sy + &
                      q(in,jn,kn,4,n) * sz
              qinf  = qfst(2,n) * sx + qfst(3,n) * sy + qfst(4,n) * sz
              pext  = q(in,jn,kn,5,n)
              pinf  = qfst(5,n)

              rhoext = q(in,jn,kn,1,n)
              a      = sqrt (gamma * pext / rhoext)

              unb   = (qext + qinf)/2 + (pext - pinf) / (2*rhoext*a)
              pb    = (pext + pinf)/2 + rhoext*a * (qext - qinf)/2

            end if

            mach = unb / a

            if (unb.le.0.d0) then

!             SUPERSONIC INFLOW: all data specified
              if (abs(mach).ge. 1.0d0) then
                q(ibc ,jbc ,kbc ,1,n) = qfst(1,n)
                q(ibc ,jbc ,kbc ,2,n) = qfst(2,n)
                q(ibc ,jbc ,kbc ,3,n) = qfst(3,n)
                q(ibc ,jbc ,kbc ,4,n) = qfst(4,n)
                q(ibc ,jbc ,kbc ,5,n) = qfst(5,n)
                q(ibc2,jbc2,kbc2,1,n) = q(ibc,jbc,kbc,1,n)
                q(ibc2,jbc2,kbc2,2,n) = q(ibc,jbc,kbc,2,n)
                q(ibc2,jbc2,kbc2,3,n) = q(ibc,jbc,kbc,3,n)
                q(ibc2,jbc2,kbc2,4,n) = q(ibc,jbc,kbc,4,n)
                q(ibc2,jbc2,kbc2,5,n) = q(ibc,jbc,kbc,5,n)
!             SUBSONIC INFLOW: 
              else if (abs(mach).lt.1.0d0) then
                rhobc  = qfst(1,n) * (pb/pinf)**(1/gamma)
                q(ibc ,jbc ,kbc ,1,n) = rhobc
                q(ibc ,jbc ,kbc ,2,n) = qfst(2,n) + (unb-qinf)*sx
                q(ibc ,jbc ,kbc ,3,n) = qfst(3,n) + (unb-qinf)*sy
                q(ibc ,jbc ,kbc ,4,n) = qfst(4,n) + (unb-qinf)*sz
                q(ibc ,jbc ,kbc ,5,n) = pb
                q(ibc2,jbc2,kbc2,1,n) = q(ibc,jbc,kbc ,1,n)
                q(ibc2,jbc2,kbc2,2,n) = q(ibc,jbc,kbc ,2,n)
                q(ibc2,jbc2,kbc2,3,n) = q(ibc,jbc,kbc ,3,n)
                q(ibc2,jbc2,kbc2,4,n) = q(ibc,jbc,kbc ,4,n)
                q(ibc2,jbc2,kbc2,5,n) = q(ibc,jbc,kbc ,5,n)
              end if

            else if (unb.gt.0.d0) then

!             SUPERSONIC OUTFLOW: all data extrapolated
              if (abs(mach).ge. 1.0d0) then
                q(ibc ,jbc ,kbc ,1,n) = q(in ,jn ,kn ,1,n)
                q(ibc ,jbc ,kbc ,2,n) = q(in ,jn ,kn ,2,n)
                q(ibc ,jbc ,kbc ,3,n) = q(in ,jn ,kn ,3,n)
                q(ibc ,jbc ,kbc ,4,n) = q(in ,jn ,kn ,4,n)
                q(ibc ,jbc ,kbc ,5,n) = q(in ,jn ,kn ,5,n)
                q(ibc2,jbc2,kbc2,1,n) = q(ibc,jbc,kbc,1,n)
                q(ibc2,jbc2,kbc2,2,n) = q(ibc,jbc,kbc,2,n)
                q(ibc2,jbc2,kbc2,3,n) = q(ibc,jbc,kbc,3,n)
                q(ibc2,jbc2,kbc2,4,n) = q(ibc,jbc,kbc,4,n)
                q(ibc2,jbc2,kbc2,5,n) = q(ibc,jbc,kbc,5,n)
!             SUBSONIC OUTFLOW: 
              else if (abs(mach).lt.1.0d0) then
                rhobc  = q(in,jn,kn,1,n) * (pb/pext)**(1/gamma)
                q(ibc ,jbc ,kbc ,1,n) = rhobc
                q(ibc ,jbc ,kbc ,2,n) = &
                                     q(in,jn,kn,2,n) + (unb-qext)*sx
                q(ibc ,jbc ,kbc ,3,n) = &
                                     q(in,jn,kn,3,n) + (unb-qext)*sy
                q(ibc ,jbc ,kbc ,4,n) = &
                                     q(in,jn,kn,4,n) + (unb-qext)*sz
                q(ibc ,jbc ,kbc ,5,n) = pb
                q(ibc2,jbc2,kbc2,1,n) = q(ibc,jbc,kbc ,1,n)
                q(ibc2,jbc2,kbc2,2,n) = q(ibc,jbc,kbc ,2,n)
                q(ibc2,jbc2,kbc2,3,n) = q(ibc,jbc,kbc ,3,n)
                q(ibc2,jbc2,kbc2,4,n) = q(ibc,jbc,kbc ,4,n)
                q(ibc2,jbc2,kbc2,5,n) = q(ibc,jbc,kbc ,5,n)
              end if

            end if

          end do
        end do
      end do

      if (kom.or.kom_bsl.or.kom_sst) then

        do n = 0,2*nharms
          nh = n*hbmove

          do i3 = istrt(ic3),iend(ic3)
            do i2 = istrt(ic2),iend(ic2)

              ibc  = ibcpt *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
              jbc  = ibcpt *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
              kbc  = ibcpt *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

              ibc2 = ibcpt2*krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
              jbc2 = ibcpt2*krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
              kbc2 = ibcpt2*krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

              in   = ibcn  *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
              jn   = ibcn  *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
              kn   = ibcn  *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

              im   = ibcm  *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
              jm   = ibcm  *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
              km   = ibcm  *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

              if (idir.eq.1) then
                sx    = si(1,im,jm,km,nh) * sgnm
                sy    = si(2,im,jm,km,nh) * sgnm
                sz    = si(3,im,jm,km,nh) * sgnm
              else if (idir.eq.2) then
                sx    = sj(1,im,jm,km,nh) * sgnm
                sy    = sj(2,im,jm,km,nh) * sgnm
                sz    = sj(3,im,jm,km,nh) * sgnm
              else if (idir.eq.3) then
                sx    = sk(1,im,jm,km,nh) * sgnm
                sy    = sk(2,im,jm,km,nh) * sgnm
                sz    = sk(3,im,jm,km,nh) * sgnm
              end if

              if (lomach) then

!-------------- LSP treatment: 

                qext = q(in,jn,kn,2,n) * sx + q(in,jn,kn,3,n) * sy + &
                       q(in,jn,kn,4,n) * sz
                qinf  = qfst(2,n) * sx + qfst(3,n) * sy + qfst(4,n) * sz
                pext = q(in,jn,kn,5,n)
                pinf = qfst(5,n)

                rhoext = q(in,jn,kn,1,n)
                q2   = q(in,jn,kn,2,n)**2 + q(in,jn,kn,3,n)**2 + &
                       q(in,jn,kn,4,n)**2
                a2   = gamma*pext/rhoext
                a    = sqrt(a2)
                un   = qext

!-------------- LSP cut-off, and artificial sound speed apr
                m2   = q2/a2
                ump2 = (uprv / a)**2
!tmp            mp2s = dmin1(1.d0, &
!tmp                   dmax1(m2,cutoff_pgr(in,jn,n), &
!tmp                   cutoff_vis(in,jn,n),(1-mixpi)*ump2,epsp(nl)**2))
                mp2s = dmin1(1.d0, &
                       dmax1(m2,(1-mixpi)*ump2,epsp(nl)**2))
                am   = (1-mp2s) / 2
                apr  = sqrt(am**2 * un**2+a2*mp2s)

                unb  = (qext + qinf)/2 + (pext - pinf) / (2*rhoext*apr) + &
                       (qext - qinf)/2 * am*un/apr

              else

!-------------- high speed treatment.

                qext  = q(in,jn,kn,2,n) * sx + q(in,jn,kn,3,n) * sy + &
                        q(in,jn,kn,4,n) * sz
                qinf  = qfst(2,n) * sx + qfst(3,n) * sy + qfst(4,n) * sz
                pext = q(in,jn,kn,5,n)
                pinf = qfst(5,n)

                rhoext = q(in,jn,kn,1,n)
                a  = sqrt (gamma * pext / rhoext)

                unb   = (qext + qinf)/2 + (pext - pinf) / (2*rhoext*a)

              end if 

              if (unb.le.0.d0) then

!               SUPERSONIC OR SUBSONIC INFLOW: all data specified
                q(ibc ,jbc ,kbc ,6,n) = tkefar
                q(ibc ,jbc ,kbc ,7,n) = omefar
                q(ibc2,jbc2,kbc2,6,n) = q(ibc,jbc,kbc,6,n)
                q(ibc2,jbc2,kbc2,7,n) = q(ibc,jbc,kbc,7,n)

              else if (unb.gt.0.d0) then

!               SUPERSONIC OR SUBSONIC OUTFLOW: all data extrapolated
!msc            q(ibc ,jbc ,kbc ,6,n) = max(q(in,jn,kn,6,n),qmin(6))
!msc            q(ibc ,jbc ,kbc ,7,n) = max(q(in,jn,kn,7,n),qmin(7))
                q(ibc ,jbc ,kbc ,6,n) = q(in ,jn ,kn ,6,n)
                q(ibc ,jbc ,kbc ,7,n) = q(in ,jn ,kn ,7,n)
                q(ibc2,jbc2,kbc2,6,n) = q(ibc,jbc,kbc,6,n)
                q(ibc2,jbc2,kbc2,7,n) = q(ibc,jbc,kbc,7,n)

              end if

            end do
          end do

        end do

      end if

      return
      end

!-----------------------------------------------------------------------
      subroutine bc_insub(q,si,sj,sk,bctopo,imax,jmax,kmax, &
                          npde,nharms,lmet,nl,ijkmaxt)
!-----------------------------------------------------------------------
!---- subsonic inflow boundary condition
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,npde,nharms,lmet,ijkmaxt
      integer(kind=cosa_int) bctopo(10)
      integer(kind=cosa_int) n,nh,i1,i2,i3,ijkmax(3),idir,istrt(3),iend(3),bctyp, &
                inrout,ibcpt,ibcpt2,ibcn,ibcn2,ibcm,ibc,jbc,kbc,ibc2, &
                jbc2,kbc2,in,jn,kn,in2,jn2,kn2,im,jm,km,ic1,ic2,ic3, &
                icell,nl,count
      real (kind=cosa_real) &
           q(-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           si(lmet, 0:imax+1, 0:jmax+1, 0:kmax+1, 0:2*nharms*hbmove), &
           sj(lmet, 0:imax+1, 0:jmax+1, 0:kmax+1, 0:2*nharms*hbmove), &
           sk(lmet, 0:imax+1, 0:jmax+1, 0:kmax+1, 0:2*nharms*hbmove), &
           qfs(ijkmaxt-1,ijkmaxt-1,npde)
      real (kind=cosa_real) nx,ny,nz,sx,sy,sz,dot,con1,con2,con3,par,f,fp,rho, &
           rho1,rho2,rc,rc1,rc2,u,u1,u2,v,v1,v2,w,w1,w2,un,un1,un2,dun, &
           p,p1,p2,p0,t0,h0inl,q2,a2,a,mp2,ump2,lambdap3,esse, &
           rcun,rcun1,rcun2,modv,sgnm

!     store imax, jmax, kmax in ijkmax to enable boundary data location

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
        ibcm   =  1
        sgnm   =  1.d0
      else
        ibcpt  = ijkmax(idir)
        ibcpt2 = ijkmax(idir) + 1
        ibcn   = ijkmax(idir) - 1
        ibcn2  = ijkmax(idir) - 2
        ibcm   = ijkmax(idir)
        sgnm   = - 1.d0
      end if

      ic1 = cyc (idir, 1)
      ic2 = cyc (idir, 2)
      ic3 = cyc (idir, 3)

      if (bctyp.eq.21) then

!---- determine (different) state in 2 rows of auxiliary cells and 
!     DO NOT enforce flux. FIRST ORDER extrapolation of primitive or
!     characteristic variables to auxiliary cells.

        do icell = 1,2

          do n = 0,2*nharms
            nh = n*hbmove

            do i3 = istrt(ic3),iend(ic3)
              do i2 = istrt(ic2),iend(ic2)
  
                ibc  = ibcpt *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
                jbc  = ibcpt *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
                kbc  = ibcpt *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

                ibc2 = ibcpt2*krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
                jbc2 = ibcpt2*krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
                kbc2 = ibcpt2*krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

                in   = ibcn  *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
                jn   = ibcn  *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
                kn   = ibcn  *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

                in2  = ibcn2 *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
                jn2  = ibcn2 *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
                kn2  = ibcn2 *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

                im   = ibcm  *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
                jm   = ibcm  *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
                km   = ibcm  *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

!-----------------------------------------------------------------------
!               Hardwired values for qfs(:,1:5)

                qfs(jbc,kbc,1) = 1.008
                qfs(jbc,kbc,2) = 0.734486515085076
                qfs(jbc,kbc,3) = 1.0
                qfs(jbc,kbc,4) = 0.0
                qfs(jbc,kbc,5) = 0.0

!-----------------------------------------------------------------------

                if (idir.eq.1) then
                  nx    = sgnm * si(1,im,jm,km,nh)
                  ny    = sgnm * si(2,im,jm,km,nh)
                  nz    = sgnm * si(3,im,jm,km,nh)
                else if (idir.eq.2) then
                  nx    = sgnm * sj(1,im,jm,km,nh)
                  ny    = sgnm * sj(2,im,jm,km,nh)
                  nz    = sgnm * sj(3,im,jm,km,nh)
                else if (idir.eq.3) then
                  nx    = sgnm * sk(1,im,jm,km,nh)
                  ny    = sgnm * sk(2,im,jm,km,nh)
                  nz    = sgnm * sk(3,im,jm,km,nh)
                end if

                rho1 = q(in ,jn ,kn ,1,n)
                rho2 = q(in2,jn2,kn2,1,n)
                u1   = q(in ,jn ,kn ,2,n)
                u2   = q(in2,jn2,kn2,2,n)
                v1   = q(in ,jn ,kn ,3,n)
                v2   = q(in2,jn2,kn2,3,n)
                w1   = q(in ,jn ,kn ,4,n)
                w2   = q(in2,jn2,kn2,4,n)
                p1   = q(in ,jn ,kn ,5,n)
                p2   = q(in2,jn2,kn2,5,n)

                un1 = u1*nx + v1*ny + w1*nz
                un2 = u2*nx + v2*ny + w2*nz

                if (.not.lomach) then
                  rc1 = sqrt(gamma *p1*rho1)
                  rc2 = sqrt(gamma *p2*rho2)
                else
                  q2 = u1**2 + v1**2 + w1**2
                  a2 = gamma * p1 / rho1
                  a = sqrt(a2)
                  mp2 = q2/a2
                  ump2 = (uprv / a)**2
                  mp2 = dmin1(dmax1(mp2,ump2,epsp(nl)**2),1.d0)
                  lambdap3 = 0.5*(un1*(1.+mp2)+sqrt(un1**2*(1.-mp2)**2+ &
                         4.*a2*mp2))
                  esse = un1 - lambdap3 
                  rc1 = - rho1 * esse
                  q2 = u2**2 + v2**2 + w2**2
                  a2 = gamma * p2 / rho2
                  a = sqrt(a2)
                  mp2 = q2/a2
                  ump2 = (uprv / a)**2
                  mp2 = dmin1(dmax1(mp2,ump2,epsp(nl)**2),1.d0)
                  lambdap3 = 0.5*(un2*(1.+mp2)+sqrt(un2**2*(1.-mp2)**2+ &
                         4.*a2*mp2))
                  esse = un2 - lambdap3 
                  rc2 = - rho2 * esse
                end if

                rcun1 = rc1 * un1
                rcun2 = rc2 * un2

                if (icell.eq.1) then
                  p  = 2*p1 - p2
                  un = 2*un1 - un2
!---------------- extrapolate characteristic directly
                  rc = 2*rc1 - rc2
                  rcun = 2*rcun1 - rcun2
!---------------- extrapolate prim. vars. and compute characteristic
!---------------- msc 10/10/09: THIS OPTION IS CURRENTLY NOT USABLE WITH LSP
!cc               rho  = 2*rho1 - rho2
!cc               rc   = sqrt(gamma*p*rho)
!cc               rcun = rc * un
                else
                  p = 3*p1 - 2*p2
                  un = 3*un1 - 2*un2
!---------------- extrapolate characteristic directly
                  rc = 3*rc1 - 2*rc2
                  rcun = 3*rcun1 - 2*rcun2
!---------------- extrapolate prim. vars. and compute characteristic
!---------------- msc 10/10/09: THIS OPTION IS CURRENTLY NOT USABLE WITH LSP
!cc               rho  = 3*rho1 - 2*rho2
!cc               rc   = sqrt(gamma*p*rho)
!cc               rcun = rc * un
                end if

                t0 = qfs(jbc,kbc,1)
                p0 = qfs(jbc,kbc,2)
                h0inl = 1 / (gamma-1) * t0
                sx    = qfs(jbc,kbc,3)
                sy    = qfs(jbc,kbc,4)
                sz    = qfs(jbc,kbc,5)
                dot   = sx*nx + sy*ny +sz*nz 

                con1 = 0.5d0 / (h0inl*dot**2)
                con2 = 2 * con1*gamma / (gamma-1)
                con3 = p - rcun

                dun = 1.d0
                count = 0

                do while(abs(dun) .gt. 1.d-10)
                  count = count + 1
                  if(count.eq.10) then
                    write(*,*) 'exceeded count limit in inflow'
                    stop
                  endif
                  par = 1.d0 - con1*un**2
                  p   = p0 * par**(gamma / (gamma-1))
  
                  f  = p - rc*un - con3
                  fp = -con2*un*p / par - rc
                  dun = - f/fp
                  un  = un + dun
                enddo

                if (un .le. 0.) then
                  write(*,*) 'Attempted outflow at inflow !'
                endif
  
                p = +rc*un + con3
                modv = un/dot
                u = modv *sx
                v = modv *sy
                w = modv *sz
                rho = gamma *p / ((gamma-1) * (h0inl-0.5d0*modv**2))

                if (icell.eq.1) then
                  q(ibc ,jbc ,kbc ,1,n) = rho
                  q(ibc ,jbc ,kbc ,2,n) = u
                  q(ibc ,jbc ,kbc ,3,n) = v
                  q(ibc ,jbc ,kbc ,4,n) = w
                  q(ibc ,jbc ,kbc ,5,n) = p
                else
                  q(ibc2,jbc2,kbc2,1,n) = rho
                  q(ibc2,jbc2,kbc2,2,n) = u
                  q(ibc2,jbc2,kbc2,3,n) = v
                  q(ibc2,jbc2,kbc2,4,n) = w
                  q(ibc2,jbc2,kbc2,5,n) = p
                end if
              end do
            end do
          end do
        end do

        if (kom.or.kom_bsl.or.kom_sst) then

          do n = 0,2*nharms
            nh = n*hbmove

            do i3 = istrt(ic3),iend(ic3)
              do i2 = istrt(ic2),iend(ic2)

                ibc  = ibcpt *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
                jbc  = ibcpt *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
                kbc  = ibcpt *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

                ibc2 = ibcpt2*krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
                jbc2 = ibcpt2*krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
                kbc2 = ibcpt2*krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

                q(ibc ,jbc ,kbc ,6,n) = tkefar
                q(ibc ,jbc ,kbc ,7,n) = omefar
                q(ibc2,jbc2,kbc2,6,n) = q(ibc,jbc,kbc,6,n)
                q(ibc2,jbc2,kbc2,7,n) = q(ibc,jbc,kbc,7,n)

              end do
            end do

          end do

        end if

      end if

      return
      end

!-----------------------------------------------------------------------
      subroutine bc_outsub(q,si,sj,sk,bctopo,imax,jmax,kmax, &
                           npde,nharms,lmet,nl,ijkmaxt)
!-----------------------------------------------------------------------
!     subsonic outlet nc for internal flows
!     p_out(1) : static pressure in first (imax) auxiliary cell
!     p_out(2) : static pressure in second (imax1) auxiliary cell
!     p_out(3) : static pressure on outlet boundary
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,npde,nharms,lmet,ijkmaxt
      integer(kind=cosa_int) bctopo(10)
      integer(kind=cosa_int) n,nh,ipde,i1,i2,i3,ijkmax(3),idir,istrt(3),iend(3), &
                bctyp,inrout,ibcpt,ibcpt2,ibcn,ibcn2,ibcm,ibc,jbc,kbc, &
                ibc2,jbc2,kbc2,in,jn,kn,in2,jn2,kn2,im,jm,km,ic1,ic2, &
                ic3,nl,icell
      real (kind=cosa_real) &
           q(-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           si(lmet,0:imax+1,0:jmax+1,0:kmax+1,0:2*nharms*hbmove), &
           sj(lmet,0:imax+1,0:jmax+1,0:kmax+1,0:2*nharms*hbmove), &
           sk(lmet,0:imax+1,0:jmax+1,0:kmax+1,0:2*nharms*hbmove), &
           p_out(ijkmaxt-1,ijkmaxt-1,3)
      real (kind=cosa_real) nx,ny,nz,rhol,rhor,pl,pr,rhom1,rho0,rho1,rho2, &
           rhotkem1,rhotke0,rhotke1,rhotke2,tke1,tke(2),p(2),pm1,p0, &
           p1,p2,c0,c1,c2,w10,w11,rho,rhou,rhov,rhow,rhoe1,rhoe2, &
           w12,w20,w21,w22,um1,u0,u1,u2,vm1,v0,v1,v2,wm1,w0,w1,w2, &
           un0,un1,un2,ut0,ut1,ut2,sgnm,qcons(7,2)

!     store imax, jmax, kmax in ijkmax to enable boundary data location

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
        ibcm   =  1
        sgnm   = -1.d0
      else
        ibcpt  = ijkmax(idir)
        ibcpt2 = ijkmax(idir) + 1
        ibcn   = ijkmax(idir) - 1
        ibcn2  = ijkmax(idir) - 2
        ibcm   = ijkmax(idir)
        sgnm   = 1.d0
      end if

      ic1 = cyc (idir, 1)
      ic2 = cyc (idir, 2)
      ic3 = cyc (idir, 3)

      if (bctyp.eq.91) then

!---- first order extrapolation of rho, rho u and rho v
!     in first row of auxiliary cells; analytical
!     static pressure in first row of auxiliary cells; first 
!     order extrapolation of rho, rho u, and rho v in second row
!     of auxiliary cells; pressure in second row of auxiliary cells
!     such that Roe averaging gives analitical boundary pressure.
!     REMARK: using this bc, the subsonic nozzle flow with 
!     analitical solution (source flow) DEMONSTRATES SECOND ORDER 
!     of the whole code.

        do n = 0,2*nharms
          nh = n*hbmove

          do i3 = istrt(ic3),iend(ic3)
            do i2 = istrt(ic2),iend(ic2)
               
              ibc  = ibcpt *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
              jbc  = ibcpt *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
              kbc  = ibcpt *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

              ibc2 = ibcpt2*krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
              jbc2 = ibcpt2*krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
              kbc2 = ibcpt2*krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

              in   = ibcn  *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
              jn   = ibcn  *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
              kn   = ibcn  *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

              in2  = ibcn2 *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
              jn2  = ibcn2 *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
              kn2  = ibcn2 *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

              im   = ibcm  *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
              jm   = ibcm  *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
              km   = ibcm  *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

!-----------------------------------------------------------------------
!             Hardwired values for:
!             p_out(1) : static pressure in first (imax) auxiliary cell
!             p_out(2) : static pressure in second (imax1) auxiliary cell
!             p_out(3) : static pressure on outlet boundary
 
              p_out(jbc,kbc,1) = 0.714285714285714
              p_out(jbc,kbc,2) = 0.714285714285714
              p_out(jbc,kbc,3) = 0.714285714285714
!-----------------------------------------------------------------------

!             from prim 2 cons due to architecture of this particular bc
              if (kom.or.kom_bsl.or.kom_sst) then
                tke(1) = q(in ,jn ,kn ,6,n)
                tke(2) = q(in2,jn2,kn2,6,n)
              else
                tke(1) = 0.d0
                tke(2) = 0.d0
              end if
              do ipde = 1,npde
                if (ipde.eq.1) then
                  qcons(ipde,1) = q(in ,jn ,kn ,ipde,n)
                  qcons(ipde,2) = q(in2,jn2,kn2,ipde,n)
                else if (ipde.eq.5) then 
                  p(1) = q(in ,jn ,kn ,ipde,n)
                  p(2) = q(in2,jn2,kn2,ipde,n)
                else
                  qcons(ipde,1) = q(in ,jn ,kn ,ipde,n)* &
                                  q(in ,jn ,kn ,1   ,n)
                  qcons(ipde,2) = q(in2,jn2,kn2,ipde,n)* &
                                  q(in2,jn2,kn2,1   ,n)
                end if
              end do
              do icell=1,2
                rho  = qcons(1,icell)
                rhou = qcons(2,icell)
                rhov = qcons(3,icell)
                rhow = qcons(4,icell)
                qcons(5,icell) = p(icell)/(gamma-1) + &
                                 (rhou*rhou+rhov*rhov+rhow*rhow)/2/rho + &
                                 rho*tke(icell)
              end do

              do ipde = 1,npde
                if (ipde.ne.5)then
                  q(ibc,jbc,kbc,ipde,n) = 2*qcons(ipde,1)-qcons(ipde,2)
                end if
              end do

              if (kom.or.kom_bsl.or.kom_sst) then
                q(ibc ,jbc ,kbc ,6   ,n) = max(q(ibc ,jbc ,kbc ,6   ,n), &
                                               q(ibc ,jbc ,kbc ,1   ,n)* &
                                               qmin(6))
                q(ibc ,jbc ,kbc ,7   ,n) = max(q(ibc ,jbc ,kbc ,7   ,n), &
                                               q(ibc ,jbc ,kbc ,1   ,n)* &
                                               qmin(7))
                rhotke0  = q(ibc ,jbc ,kbc ,6,n)
              else
                rhotke0  = 0.d0
              end if

              q(ibc,jbc,kbc,5,n) = p_out(jbc,kbc,1)/(gamma-1) + &
                                   (q(ibc,jbc,kbc,2,n)**2 + &
                                    q(ibc,jbc,kbc,3,n)**2 + &
                                    q(ibc,jbc,kbc,4,n)**2) / 2 / &
                                   q(ibc,jbc,kbc,1,n) + &
                                   rhotke0

              do ipde = 1,npde
                if (ipde.ne.5) then
                  q(ibc2,jbc2,kbc2,ipde,n) = 3*qcons(ipde,1) - &
                                             2*qcons(ipde,2)
                end if
              end do

              if (kom.or.kom_bsl.or.kom_sst) then
                q(ibc2,jbc2,kbc2,6   ,n) = max(q(ibc2,jbc2,kbc2,6   ,n), &
                                               q(ibc2,jbc2,kbc2,1   ,n)* &
                                               qmin(6))
                q(ibc2,jbc2,kbc2,7   ,n) = max(q(ibc2,jbc2,kbc2,7   ,n), &
                                               q(ibc2,jbc2,kbc2,1   ,n)* &
                                               qmin(7))
                rhotke2  = qcons(6,2)
                rhotke1  = qcons(6,1)
                rhotkem1 = q(ibc2,jbc2,kbc2,6,n)
              else
                rhotke2  = 0.d0
                rhotke1  = 0.d0
                rhotkem1 = 0.d0
              end if

              rhol = ( 3*qcons(1,1) - qcons(1,2) ) / 2
              rhor = ( 3*q(ibc,jbc,kbc,1,n) - q(ibc2,jbc2,kbc2,1,n) )/ 2
              p1  = (gamma-1)* &
                    (qcons(5,1) - rhotke1 - &
                     (qcons(2,1)**2+qcons(3,1)**2+qcons(4,1)**2) / &
                     2 / qcons(1,1))
              p2 = (gamma-1)* &
                   (qcons(5,2) - rhotke2 - &
                   (qcons(2,2)**2+qcons(3,2)**2+qcons(4,2)**2) / &
                   2 / qcons(1,2))
              pl = ( 3*p1 - p2 ) / 2
              pr = (p_out(jbc,kbc,3) * (sqrt(rhol) + sqrt(rhor)) - &
                      pl * sqrt(rhol) ) / sqrt(rhor)
              pm1 = 3*p_out(jbc,kbc,1) - 2*pr
              q(ibc2,jbc2,kbc2,5,n) = pm1/(gamma-1) + &
                                      (q(ibc2,jbc2,kbc2,2,n)**2 + &
                                       q(ibc2,jbc2,kbc2,3,n)**2 + &
                                       q(ibc2,jbc2,kbc2,4,n)**2) / 2 / &
                                      q(ibc2,jbc2,kbc2,1,n) + &
                                      rhotkem1

!             transofrm auxiliary cell value from cons to prim
              do ipde=2,npde
                if (ipde.ne.5) then
                  q(ibc ,jbc ,kbc ,ipde,n) = q(ibc ,jbc ,kbc ,ipde,n) / &
                                             q(ibc ,jbc ,kbc ,1   ,n)
                  q(ibc2,jbc2,kbc2,ipde,n) = q(ibc2,jbc2,kbc2,ipde,n) / &
                                             q(ibc2,jbc2,kbc2,1   ,n)
                end if
              end do
              if (kom.or.kom_bsl.or.kom_sst) then
                tke(1) = q(ibc ,jbc ,kbc ,6,n)
                tke(2) = q(ibc2,jbc2,kbc2,6,n)
              else
                tke(1) = 0.d0
                tke(2) = 0.d0
              end if
              rho1  = q(ibc ,jbc ,kbc ,1,n)
              rho2  = q(ibc2,jbc2,kbc2,1,n)
              u1    = q(ibc ,jbc ,kbc ,2,n)
              u2    = q(ibc2,jbc2,kbc2,2,n)
              v1    = q(ibc ,jbc ,kbc ,3,n)
              v2    = q(ibc2,jbc2,kbc2,3,n)
              w1    = q(ibc ,jbc ,kbc ,4,n)
              w2    = q(ibc2,jbc2,kbc2,4,n)
              rhoe1 = q(ibc ,jbc ,kbc ,5,n)
              rhoe2 = q(ibc2,jbc2,kbc2,5,n)
              q(ibc ,jbc ,kbc ,5,n) = (gamma-1)* &
                (rhoe1 - rho1*(u1*u1+v1*v1+w1*w1)/2 - rho1*tke(1))
              q(ibc2,jbc2,kbc2,5,n) = (gamma-1)* &
                (rhoe2 - rho2*(u2*u2+v2*v2+w2*w2)/2 - rho2*tke(2))

            end do
          end do

        end do

      end if

      return
      end
!-----------------------------------------------------------------------
      subroutine bc_walinv_14(q,si,sj,sk,bctopo,imax,jmax,kmax,npde, &
                              nharms,lmet)
!-----------------------------------------------------------------------
!---- inviscid wall bc based on tangential velocity extrapolation
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,npde,nharms,lmet
      integer(kind=cosa_int) bctopo(10)
      integer(kind=cosa_int) n,nh,i1,i2,i3,ijkmax(3),idir,istrt(3),iend(3),bctyp, &
                inrout,ibcpt,ibcpt2,ibcn,ibcn2,ibcm,ibc,jbc,kbc,ibc2, &
                jbc2,kbc2,in,jn,kn,in2,jn2,kn2,im,jm,km,ic1,ic2,ic3
      real (kind=cosa_real) &
          q    (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
          si   (lmet, 0:imax+1, 0:jmax+1, 0:kmax+1,0:2*nharms*hbmove), &
          sj   (lmet, 0:imax+1, 0:jmax+1, 0:kmax+1,0:2*nharms*hbmove), &
          sk   (lmet, 0:imax+1, 0:jmax+1, 0:kmax+1,0:2*nharms*hbmove)
      real (kind=cosa_real) sgnm,nx,ny,nz,xwdot,ywdot,zwdot,rhow,rho1,rho2,uw, &
           u1,u2,vw,v1,v2,ww,w1,w2,pw,p1,p2,unw,un1,un2

!     store imax, jmax in ijkmax to enable boundary data location

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

      ic1 = cyc (idir, 1)
      ic2 = cyc (idir, 2)
      ic3 = cyc (idir, 3)

      do n = 0,2*nharms
        nh = n*hbmove

        do i3 = istrt(ic3),iend(ic3)
          do i2 = istrt(ic2),iend(ic2)

            ibc  = ibcpt *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
            jbc  = ibcpt *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
            kbc  = ibcpt *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

            ibc2 = ibcpt2*krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
            jbc2 = ibcpt2*krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
            kbc2 = ibcpt2*krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

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
              nx    = si(1,im,jm,km,nh)
              ny    = si(2,im,jm,km,nh)
              nz    = si(3,im,jm,km,nh)
              if (moving) then
                xwdot = si(5,im,jm,km,nh)
                ywdot = si(6,im,jm,km,nh)
                zwdot = si(7,im,jm,km,nh)
              end if
            else if (idir.eq.2) then
              nx    = sj(1,im,jm,km,nh)
              ny    = sj(2,im,jm,km,nh)
              nz    = sj(3,im,jm,km,nh)
              if (moving) then
                xwdot = sj(5,im,jm,km,nh)
                ywdot = sj(6,im,jm,km,nh)
                zwdot = sj(7,im,jm,km,nh)
              end if
            else if (idir.eq.3) then
              nx    = sk(1,im,jm,km,nh)
              ny    = sk(2,im,jm,km,nh)
              nz    = sk(3,im,jm,km,nh)
              if (moving) then
                xwdot = sk(5,im,jm,km,nh)
                ywdot = sk(6,im,jm,km,nh)
                zwdot = sk(7,im,jm,km,nh)
              end if
            end if

            rho1 = q(in ,jn ,kn ,1,n)
            u1   = q(in ,jn ,kn ,2,n)
            v1   = q(in ,jn ,kn ,3,n)
            w1   = q(in ,jn ,kn ,4,n)
            p1   = q(in ,jn ,kn ,5,n)

            rho2 = q(in2,jn2,kn2,1,n)
            u2   = q(in2,jn2,kn2,2,n)
            v2   = q(in2,jn2,kn2,3,n)
            w2   = q(in2,jn2,kn2,4,n)
            p2   = q(in2,jn2,kn2,5,n)

            rhow = (3*rho1 - rho2) / 2
            uw   = (3*u1   - u2  ) / 2
            vw   = (3*v1   - v2  ) / 2
            ww   = (3*w1   - w2  ) / 2
            pw   = (3*p1   - p2  ) / 2

            un1  = u1*nx + v1*ny + w1*nz
            un2  = u2*nx + v2*ny + w2*nz
            unw  = (3*un1-un2) / 2

            if (moving) then
              unw = unw - (xwdot*nx + ywdot*ny + zwdot*nz)
            endif 

            q(ibc ,jbc ,kbc ,1,n) = rhow
            q(ibc ,jbc ,kbc ,2,n) = uw - unw*nx
            q(ibc ,jbc ,kbc ,3,n) = vw - unw*ny
            q(ibc ,jbc ,kbc ,4,n) = ww - unw*nz
            q(ibc ,jbc ,kbc ,5,n) = pw
            q(ibc2,jbc2,kbc2,1,n) = q(ibc,jbc,kbc,1,n)
            q(ibc2,jbc2,kbc2,2,n) = q(ibc,jbc,kbc,2,n)
            q(ibc2,jbc2,kbc2,3,n) = q(ibc,jbc,kbc,3,n)
            q(ibc2,jbc2,kbc2,4,n) = q(ibc,jbc,kbc,4,n)
            q(ibc2,jbc2,kbc2,5,n) = q(ibc,jbc,kbc,5,n)

          end do
        end do

      end do

      if (kom.or.kom_bsl.or.kom_sst) then

        do n = 0,2*nharms
          nh = n*hbmove

          do i3 = istrt(ic3),iend(ic3)
            do i2 = istrt(ic2),iend(ic2)

              ibc  = ibcpt *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
              jbc  = ibcpt *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
              kbc  = ibcpt *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

              ibc2 = ibcpt2*krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
              jbc2 = ibcpt2*krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
              kbc2 = ibcpt2*krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

              in   = ibcn  *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
              jn   = ibcn  *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
              kn   = ibcn  *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

              in2  = ibcn2 *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
              jn2  = ibcn2 *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
              kn2  = ibcn2 *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

              im   = ibcm  *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
              jm   = ibcm  *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
              km   = ibcm  *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

              q(ibc ,jbc ,kbc ,6,n) = q(in ,jn ,kn ,6,n)
              q(ibc ,jbc ,kbc ,7,n) = q(in ,jn ,kn ,7,n)
              q(ibc2,jbc2,kbc2,6,n) = q(ibc,jbc,kbc,6,n)
              q(ibc2,jbc2,kbc2,7,n) = q(ibc,jbc,kbc,7,n)

            end do
          end do

        end do

      end if

      return
      end

!-----------------------------------------------------------------------
      subroutine cutman_init(cutopo)
!-----------------------------------------------------------------------
!     The subroutine initializes data exchange along cuts.
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
!     Consider that, in order to prepare data exchange, some constants
!     have to be adequately set
!     * maxproc ----> maximum number of MPI processes
!     * maxproc_cutman ---> maximum number of MPI processes one process exchange data
!     * mrequest_perproc ---> maximum number of exchanged points for each send/recv operation
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

      call allocate_cutman_init_memory()

      maxsendid_num = 0
      maxrecvid_num = 0
      if(.not. allocated(mysendid)) write(*,*) "ERROR"
      mysendid = -1 ! array syntax
      myrecvid = -1 ! array syntax 
      myrecvid_num = -1 ! array syntax
      mysendid_num = -1 ! array syntax

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
            call getownerid(iblk2,recvfromid)
            found = .false.
            do i_id=1,maxrecvid_num
               if(myrecvid(i_id) .eq. recvfromid) found = .true. 
            enddo
            if(.not.found) then
              maxrecvid_num = maxrecvid_num + 1
              myrecvid(maxrecvid_num) = recvfromid
              myrecvid_num(recvfromid) = maxrecvid_num
            endif
            maxrecv(recvfromid) = maxrecv(recvfromid) + &
                                  (length(ic3)+1)*(length(ic2)+1)*3
         else if(myblock2) then
            call getownerid(iblk1,sendtoid)
            found = .false.
            do i_id=1,maxsendid_num
               if(mysendid(i_id) .eq. sendtoid) found = .true. 
            enddo
            if(.not.found) then
              maxsendid_num = maxsendid_num + 1
              mysendid(maxsendid_num) = sendtoid
              mysendid_num(sendtoid) = maxsendid_num
            endif
            maxsend(sendtoid) = maxsend(sendtoid) + &
                                (length(ic3)+1)*(length(ic2)+1)*3
         end if
      end do

      mrequest_perproc = max(maxval(maxsend), maxval(maxrecv))
      maxproc_cutman = max(maxsendid_num,maxrecvid_num)

      call allocate_cutman_memory()

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
      subroutine cutman_q(nl,cutopo,q)
!-----------------------------------------------------------------------

      use parallelutils
      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,iblk1,icut,iblk2,imax1,jmax1,kmax1,imax2,jmax2,kmax2, &
                iq1,iq2
      integer(kind=cosa_int) cutopo(21,myncuts),a,b,c,i,reqnum
      logical myblock1, myblock2
      real(kind=cosa_real) q(*)
      integer(kind=cosa_int) sendrequestnum,receiverequestnum
      integer(kind=cosa_int) i_id,i_in_id
      integer(kind=cosa_int) sendid,recvid
      integer myid

      numrecv = 0  ! array syntax
      numsend = 0  ! array syntax
      sendrecv_datasize = npde*((2*nharms)+1)

      allocate(sendarray(sendrecv_datasize, mrequest_perproc, &
                         maxproc_cutman))
      allocate(receivearray(sendrecv_datasize, mrequest_perproc, &
                            maxproc_cutman))

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
            iq1    = 1 + off_p3 (iblk1,nl) * npde * dim5
         else
            iq1    = 1
         end if
         if(myblock2) then
            imax2  = i_imax     (iblk2,nl)
            jmax2  = j_jmax     (iblk2,nl)
            kmax2  = k_kmax     (iblk2,nl)
            iq2    = 1 + off_p3 (iblk2,nl) * npde * dim5
         else
            iq2    = 1
         end if
         call cut_q(imax1,jmax1,kmax1,q(iq1),imax2,jmax2,kmax2,q(iq2), &
              npde,nharms,cutopo(1,icut),receivearray,receiverequests, &
              numrecv,receiveindices,sendarray,sendrequests,numsend,iq1, &
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
           iq1   = receiveindices(i_in_id,4,reqnum)
           imax1 = receiveindices(i_in_id,5,reqnum)
           jmax1 = receiveindices(i_in_id,6,reqnum)
           kmax1 = receiveindices(i_in_id,7,reqnum)
           call copyqreceivedata(imax1,jmax1,kmax1,a,b,c,q(iq1), &
                receivearray,i_in_id,reqnum)
         enddo
      enddo

      call waitallmessages(sendrequests,maxsendid_num)

      deallocate(sendarray, receivearray)

      return

      end

!-----------------------------------------------------------------------
      subroutine cut_q(imax1,jmax1,kmax1,q1,imax2,jmax2,kmax2,q2,npde, &
        nharms,cutopo,receivearray,receiverequests,numrecv, &
        receiveindices,sendarray,sendrequests,numsend,iq1, &
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

      integer(kind=cosa_int) imax1,jmax1,kmax1,imax2,jmax2,kmax2,npde,nharms
      integer(kind=cosa_int) cutopo(21),mrequest,iblk1,iblk2
      integer(kind=cosa_int) ijkmax1(3),ijkmax2(3),istr1(3),iend1(3),istr2(3), &
           iend2(3),isgn1(3),isgn2(3),iord(3),len(3),idir1,idir2,inout1, &
           inout2,ibcpt,ibcinc,ibcpt2,inr,inrinc,inr2,l,ipde,n,ic1,ic2, &
           ic3,jc1,jc2,jc3,i2,i3,ii1,jj1,kk1,ii2,jj2,kk2,in1,jn1,kn1, &
           in2,jn2,kn2
      real(kind=cosa_real) &
           q1(-1:imax1+1,-1:jmax1+1,-1:kmax1+1,npde,0:2*nharms), &
           q2(-1:imax2+1,-1:jmax2+1,-1:kmax2+1,npde,0:2*nharms)
      integer(kind=cosa_int) tempindex, datasize
      integer(kind=cosa_int) receiverequestnum,iq1,sendrequestnum
      logical myblock1, myblock2
      integer(kind=cosa_int) mrequest_perproc,maxproc_cutman
      integer(kind=cosa_int) sendid,recvid,numrecvid,numsendid,maxproc
      integer(kind=cosa_int) sendrequests(maxproc_cutman),receiverequests(maxproc_cutman)
      integer(kind=cosa_int) receiveindices(mrequest_perproc,8,maxproc_cutman)
      real(kind=cosa_real) &
        sendarray   (npde*(2*nharms+1),mrequest_perproc,maxproc_cutman), &
        receivearray(npde*(2*nharms+1),mrequest_perproc,maxproc_cutman)
      integer(kind=cosa_int) numsend(maxproc_cutman),numrecv(maxproc_cutman)
      integer(kind=cosa_int) myrecvid_num(0:maxproc-1),mysendid_num(0:maxproc-1)
!--------------------------------------------------------------------------

!     Store imax, jmax, kmax in ijmax for looping
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
         ibcinc = - 1
      else
         ibcpt  = ijkmax1(idir1)
         ibcinc = 1
      end if
!     
      if (inout2 .eq. 1) then
         inr    =   1
         inrinc =   1
      else
         inr    =   ijkmax2(idir2) - 1
         inrinc = - 1
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

          ii2 = ii1 + ibcinc * krd (ic1, 1)
          jj2 = jj1 + ibcinc * krd (ic1, 2)
          kk2 = kk1 + ibcinc * krd (ic1, 3)

          in1 = inr                          * krd (jc1, 1) + &
                (istr2(jc2) + isgn2(jc2)*i2) * krd (jc2, 1) + &
                (istr2(jc3) + isgn2(jc3)*i3) * krd (jc3, 1)
          jn1 = inr                          * krd (jc1, 2) + &
                (istr2(jc2) + isgn2(jc2)*i2) * krd (jc2, 2) + &
                (istr2(jc3) + isgn2(jc3)*i3) * krd (jc3, 2)
          kn1 = inr                          * krd (jc1, 3) + &
                (istr2(jc2) + isgn2(jc2)*i2) * krd (jc2, 3) + &
                (istr2(jc3) + isgn2(jc3)*i3) * krd (jc3, 3)

          in2 = in1 + inrinc * krd (jc1, 1)
          jn2 = jn1 + inrinc * krd (jc1, 2)
          kn2 = kn1 + inrinc * krd (jc1, 3)

          if (myblock1 .and. myblock2) then

            do n = 0, 2*nharms
              do ipde = 1, npde
                 q1(ii1,jj1,kk1,ipde,n) = q2(in1,jn1,kn1,ipde,n)
                 q1(ii2,jj2,kk2,ipde,n) = q2(in2,jn2,kn2,ipde,n)
              end do
            end do

          else if(myblock1) then

            call getownerid(iblk2,recvid)
            numrecvid = myrecvid_num(recvid)
            numrecv(numrecvid) = numrecv(numrecvid)+1 
            receiveindices(numrecv(numrecvid),1,numrecvid) = ii1
            receiveindices(numrecv(numrecvid),2,numrecvid) = jj1
            receiveindices(numrecv(numrecvid),3,numrecvid) = kk1
            receiveindices(numrecv(numrecvid),4,numrecvid) = iq1
            receiveindices(numrecv(numrecvid),5,numrecvid) = imax1
            receiveindices(numrecv(numrecvid),6,numrecvid) = jmax1
            receiveindices(numrecv(numrecvid),7,numrecvid) = kmax1
            numrecv(numrecvid) = numrecv(numrecvid)+1
            receiveindices(numrecv(numrecvid),1,numrecvid) = ii2
            receiveindices(numrecv(numrecvid),2,numrecvid) = jj2
            receiveindices(numrecv(numrecvid),3,numrecvid) = kk2
            receiveindices(numrecv(numrecvid),4,numrecvid) = iq1
            receiveindices(numrecv(numrecvid),5,numrecvid) = imax1
            receiveindices(numrecv(numrecvid),6,numrecvid) = jmax1
            receiveindices(numrecv(numrecvid),7,numrecvid) = kmax1
            if (numrecv(numrecvid).gt. mrequest_perproc) then
               write(*,*) 'numrecv (cut_q)--> INCREASE MREQUEST_PERPROC'
               call abortmpi()
            end if

          else if(myblock2) then

            call getownerid(iblk1,sendid)
            numsendid = mysendid_num(sendid)
            numsend(numsendid) = numsend(numsendid)+1 
            tempindex = 1
            do n = 0, 2*nharms
              do ipde = 1, npde
                sendarray(tempindex,numsend(numsendid),numsendid) = &
                  q2(in1,jn1,kn1,ipde,n)
                tempindex = tempindex + 1
              end do
            end do
            numsend(numsendid) = numsend(numsendid)+1
            tempindex  = 1
            do n = 0, 2*nharms
              do ipde = 1, npde
                sendarray(tempindex,numsend(numsendid),numsendid) = &
                  q2(in2,jn2,kn2,ipde,n)
                tempindex = tempindex + 1
              end do
            end do
            if (numsend(numsendid).gt. mrequest_perproc) then
               write(*,*) 'numsend (cut_q)--> INCREASE MREQUEST_PERPROC'
               call abortmpi()
            end if

          end if

        end do
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine bc_walvis_15(q,si,sj,sk,dist,bctopo,imax,jmax,kmax, &
                              npde,nharms,lmet)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,npde,nharms,lmet
      integer(kind=cosa_int) bctopo(10)
      integer(kind=cosa_int) n,nh,i1,i2,i3,ijkmax(3),idir,istrt(3),iend(3),bctyp, &
                inrout,ibcpt,ibcpt2,ibcn,ibcn2,ibcm,ibc,jbc,kbc,ibc2, &
                jbc2,kbc2,in,jn,kn,in2,jn2,kn2,im,jm,km,ic1,ic2,ic3
      real (kind=cosa_real) &
          q    (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
          dist ( 0:imax  , 0:jmax  , 0:kmax), &
          si   ( lmet, 0:imax+1, 0:jmax+1, 0:kmax+1, 0:2*nharms*hbmove), &
          sj   ( lmet, 0:imax+1, 0:jmax+1, 0:kmax+1, 0:2*nharms*hbmove), &
          sk   ( lmet, 0:imax+1, 0:jmax+1, 0:kmax+1, 0:2*nharms*hbmove)
      real (kind=cosa_real) small,smasq,p1,p2,pw,rho1,rhow,tw,muw,uw,vw,ww,nx, &
           ny,nz,omew,uc,vc,wc,upar,uorth,dupardn,tauw,rkrpls,sr,tke1, &
           tke2,sgnm,mor

      small  = 1.d-12
      smasq  = 1.d-24
      uw     = 0.d0
      vw     = 0.d0
      ww     = 0.d0
      mor    = machfs / reyno

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

      ic1 = cyc (idir, 1)
      ic2 = cyc (idir, 2)
      ic3 = cyc (idir, 3)

      do n = 0,2*nharms
        nh = n*hbmove

        do i3 = istrt(ic3),iend(ic3)
          do i2 = istrt(ic2),iend(ic2)

            ibc  = ibcpt *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
            jbc  = ibcpt *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
            kbc  = ibcpt *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

            ibc2 = ibcpt2*krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
            jbc2 = ibcpt2*krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
            kbc2 = ibcpt2*krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

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
              nx    = si(1,im,jm,km,nh)
              ny    = si(2,im,jm,km,nh)
              nz    = si(3,im,jm,km,nh)
              if (moving) then
                uw = si(5,im,jm,km,nh)
                vw = si(6,im,jm,km,nh)
                ww = si(7,im,jm,km,nh)
              end if
            else if (idir.eq.2) then
              nx    = sj(1,im,jm,km,nh)
              ny    = sj(2,im,jm,km,nh)
              nz    = sj(3,im,jm,km,nh)
              if (moving) then
                uw =  sj(5,im,jm,km,nh)
                vw =  sj(6,im,jm,km,nh)
                ww =  sj(7,im,jm,km,nh)
              end if
            else if (idir.eq.3) then
              nx    = sk(1,im,jm,km,nh)
              ny    = sk(2,im,jm,km,nh)
              nz    = sk(3,im,jm,km,nh)
              if (moving) then
                uw =  sk(5,im,jm,km,nh)
                vw =  sk(6,im,jm,km,nh)
                ww =  sk(7,im,jm,km,nh)
              end if
            end if

            p1 = q(in ,jn ,kn ,5,n)
            p2 = q(in2,jn2,kn2,5,n)
!           linear extrapolation:
            pw = 0.5d0 * (3*p1-p2)
!           quadratic extrapolation based on zero pressure gradient
!           pw    = (9.d0 * p1 - p2) / 8.d0

            rho1 = q(in,jn,kn,1,n)
            rhow = (pw/p1)**(1/gamma) * q(in,jn,kn,1,n)

            q(ibc ,jbc ,kbc ,1,n) = rhow
            q(ibc ,jbc ,kbc ,2,n) = uw
            q(ibc ,jbc ,kbc ,3,n) = vw
            q(ibc ,jbc ,kbc ,4,n) = ww
            q(ibc ,jbc ,kbc ,5,n) = pw
            q(ibc2,jbc2,kbc2,1,n) = q(ibc,jbc,kbc,1,n)
            q(ibc2,jbc2,kbc2,2,n) = q(ibc,jbc,kbc,2,n)
            q(ibc2,jbc2,kbc2,3,n) = q(ibc,jbc,kbc,3,n)
            q(ibc2,jbc2,kbc2,4,n) = q(ibc,jbc,kbc,4,n)
            q(ibc2,jbc2,kbc2,5,n) = q(ibc,jbc,kbc,5,n)

          end do
        end do

      end do

      if (kom.or.kom_bsl.or.kom_sst) then

        if (wallwilc) then

          do n = 0,2*nharms
            nh = n*hbmove

            do i3 = istrt(ic3),iend(ic3)
              do i2 = istrt(ic2),iend(ic2)

                ibc  = ibcpt *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
                jbc  = ibcpt *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
                kbc  = ibcpt *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

                ibc2 = ibcpt2*krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
                jbc2 = ibcpt2*krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
                kbc2 = ibcpt2*krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

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
                  nx    = si(1,im,jm,km,nh)
                  ny    = si(2,im,jm,km,nh)
                  nz    = si(3,im,jm,km,nh)
                  if (moving) then
                    uw = si(5,im,jm,km,nh)
                    vw = si(6,im,jm,km,nh)
                    ww = si(7,im,jm,km,nh)
                  end if
                else if (idir.eq.2) then
                  nx    = sj(1,im,jm,km,nh)
                  ny    = sj(2,im,jm,km,nh)
                  nz    = sj(3,im,jm,km,nh)
                  if (moving) then
                    uw =  sj(5,im,jm,km,nh)
                    vw =  sj(6,im,jm,km,nh)
                    ww =  sj(7,im,jm,km,nh)
                  end if
                else if (idir.eq.3) then
                  nx    = sk(1,im,jm,km,nh)
                  ny    = sk(2,im,jm,km,nh)
                  nz    = sk(3,im,jm,km,nh)
                  if (moving) then
                    uw =  sk(5,im,jm,km,nh)
                    vw =  sk(6,im,jm,km,nh)
                    ww =  sk(7,im,jm,km,nh)
                  end if
                end if

                rhow = q(ibc,jbc,kbc,1,n)
                pw   = q(ibc,jbc,kbc,5,n)
                tw   = pw * gamma / rhow
                muw  = (stemp+1)/(tw+stemp) * (tw**1.5d0)

!-------------- wall BC based on Wilcox's rough wall approach
!               Calculate the velocity parallel to the wall at the cell 
!               centre off the wall.  
!
!                     U(parallel) = U - U(orthogonal) 
!
!               where U is the total velocity vector
!               and U(orthogonal) is the velocity normal to the wall.  
!               For fixed-body problems, velocity on the wall is zero.
!               For moving-body problems, velocity on the wall is 
!               nonzero.
!               uorth is the velocity component (with sign) normal to
!               the wall at the first cell centre off the wall.
!               upar is the modulus of the velocity component parallel
!               to the wall at the first cell centre off the wall.

                uc = q(in,jn,kn,2,n) - uw
                vc = q(in,jn,kn,3,n) - vw
                wc = q(in,jn,kn,4,n) - ww
                uorth = nx * uc + ny * vc + nz * wc
                upar  = sqrt( (uc-uorth*nx)**2 + (vc-uorth*ny)**2 + &
                              (wc-uorth*nz)**2 )

!               Tau(wall) = mu * du/dn

                dupardn = (upar      ) / (dist(in,jn,kn) + small)
                tauw    = muw * abs( dupardn )

!               Kr+         = utau *  rho(wall) * Kr / mu(wall)
!               Omega(wall) = utau**2 rho(wall) * Sr / mu(wall)
!               If Kr = 0, set Kr+ = 5 -> hydraulically smooth

                if (roughk.gt.0.d0) then
!msc              After Wilcox, AIAA J. 1988:
!msc              rkrpls = sqrt (rhow * tauw) * roughk / muw
!msc              After Ferrer and Munduate,, AIAA} paper 2009-269
                  rkrpls = &
                    max(1.d0,sqrt (rhow * tauw / mor) * roughk / muw)
                else
!msc              After Wilcox, AIAA J. 1988:
!msc              rkrpls = 5.d0
!msc              After Ferrer and Munduate,, AIAA} paper 2009-269
                  rkrpls = 1.d0
                end if
!
                if (rkrpls.ge.25.d0) then
                   sr =  100.d0 / rkrpls
                else
                   sr = 2500.d0 / (rkrpls * rkrpls)
                end if

                omew = abs(dupardn) * sr

                q(ibc ,jbc ,kbc ,6,n) = 0.d0
                q(ibc ,jbc ,kbc ,7,n) = omew
                q(ibc2,jbc2,kbc2,6,n) = q(ibc,jbc,kbc,6,n)
                q(ibc2,jbc2,kbc2,7,n) = q(ibc,jbc,kbc,7,n)

              end do
            end do

          end do

        else if (wallment) then

          do n = 0,2*nharms
            nh = n*hbmove

            do i3 = istrt(ic3),iend(ic3)
              do i2 = istrt(ic2),iend(ic2)

                ibc  = ibcpt *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
                jbc  = ibcpt *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
                kbc  = ibcpt *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

                ibc2 = ibcpt2*krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
                jbc2 = ibcpt2*krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
                kbc2 = ibcpt2*krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

                in   = ibcn  *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
                jn   = ibcn  *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
                kn   = ibcn  *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

                in2  = ibcn2 *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
                jn2  = ibcn2 *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
                kn2  = ibcn2 *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

                im   = ibcm  *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
                jm   = ibcm  *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
                km   = ibcm  *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

                rhow = q(ibc,jbc,kbc,1,n)
                pw   = q(ibc,jbc,kbc,5,n)
                tw   = pw * gamma / rhow
                muw  = (stemp+1)/(tw+stemp) * (tw**1.5d0)

!-------------- wall BC based on Menter (NASA TM 103975)
                omew = 60.d0*muw / (bkw*rhow*dist(in,jn,kn)**2+smasq) &
                        *machfs/reyno

                q(ibc ,jbc ,kbc ,6,n) = 0.d0
                q(ibc ,jbc ,kbc ,7,n) = omew
                q(ibc2,jbc2,kbc2,6,n) = q(ibc,jbc,kbc,6,n)
                q(ibc2,jbc2,kbc2,7,n) = q(ibc,jbc,kbc,7,n)

              end do
            end do

          end do

        end if

      end if

      return
      end

!-----------------------------------------------------------------------
      subroutine bc_extrapol(q,bctopo,imax,jmax,kmax,npde,nharms)
!-----------------------------------------------------------------------
!---- extrapolation from inside the domain
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,npde,nharms
      integer(kind=cosa_int) bctopo(10)
      integer(kind=cosa_int) n,nh,ipde,i1,i2,i3,ijkmax(3),idir,istrt(3),iend(3), &
                bctyp,inrout,ibcpt,ibcpt2,ibcn,ibcm,ibc,jbc,kbc,ibc2, &
                jbc2,kbc2,in,jn,kn,ic1,ic2,ic3
      real (kind=cosa_real) &
          q    (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms)

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
        ibcm   =  1
      else
        ibcpt  = ijkmax(idir)
        ibcpt2 = ijkmax(idir) + 1
        ibcn   = ijkmax(idir) - 1
        ibcm   = ijkmax(idir)
      end if

      ic1 = cyc (idir, 1)
      ic2 = cyc (idir, 2)
      ic3 = cyc (idir, 3)

      do n = 0,2*nharms
        do ipde = 1,npde

          do i3 = istrt(ic3),iend(ic3)
            do i2 = istrt(ic2),iend(ic2)

              ibc  = ibcpt *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
              jbc  = ibcpt *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
              kbc  = ibcpt *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

              ibc2 = ibcpt2*krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
              jbc2 = ibcpt2*krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
              kbc2 = ibcpt2*krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

              in   = ibcn  *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
              jn   = ibcn  *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
              kn   = ibcn  *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

              q(ibc ,jbc ,kbc ,ipde,n) = q(in,jn,kn,ipde,n)
              q(ibc2,jbc2,kbc2,ipde,n) = q(in,jn,kn,ipde,n)

            end do
          end do

        end do
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine bc_sym(q,bctopo,imax,jmax,kmax,npde,nharms)
!-----------------------------------------------------------------------
!     symmetry wrt xy-plane                                       idir 3
!     symmetry wrt xz-plane                                       idir 2
!     symmetry wrt yz-plane                                       idir 1
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,npde,nharms
      integer(kind=cosa_int) bctopo(10)
      integer(kind=cosa_int) ipde,n,i1,i2,i3,ijkmax(3),idir,istrt(3),iend(3),bctyp, &
                inrout,ibcpt,ibcpt2,ibcn,ibcn2,ibc,jbc,kbc,ibc2,jbc2, &
                kbc2,in,jn,kn,in2,jn2,kn2,ic1,ic2,ic3,isym(7)
      real (kind=cosa_real) &
          q    (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms)

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
      do ipde = 1,npde
         isym(ipde) = 1
      end do

      isym(idir+1) = -1

      do n = 0,2*nharms
        do ipde = 1,npde

          do i3 = istrt(ic3),iend(ic3)
            do i2 = istrt(ic2),iend(ic2)

              ibc  = ibcpt *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
              jbc  = ibcpt *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
              kbc  = ibcpt *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

              ibc2 = ibcpt2*krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
              jbc2 = ibcpt2*krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
              kbc2 = ibcpt2*krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

              in   = ibcn  *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
              jn   = ibcn  *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
              kn   = ibcn  *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

              in2  = ibcn2 *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
              jn2  = ibcn2 *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
              kn2  = ibcn2 *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

              q(ibc ,jbc ,kbc ,ipde,n) = isym(ipde) * &
                                         q(in ,jn ,kn ,ipde,n)
              q(ibc2,jbc2,kbc2,ipde,n) = isym(ipde) * &
                                         q(in2,jn2,kn2,ipde,n)

            end do
          end do

        end do
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine cutman_mut(nl,cutopo,mut)
!-----------------------------------------------------------------------

      use parallelutils
      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,iblk1,icut,iblk2,imax1,jmax1,kmax1,imax2,jmax2,kmax2, &
                imut1,imut2
      integer(kind=cosa_int) cutopo(21,myncuts),a,b,c,i,reqnum
      logical myblock1, myblock2
      real(kind=cosa_real) mut(*)
      integer(kind=cosa_int) sendrequestnum,receiverequestnum
      integer(kind=cosa_int) i_id,i_in_id
      integer(kind=cosa_int) sendid,recvid

      numrecv = 0  ! array syntax
      numsend = 0  ! array syntax
      sendrecv_datasize = 2*nharms+1

      allocate(sendarray(sendrecv_datasize, mrequest_perproc, &
                         maxproc_cutman))
      allocate(receivearray(sendrecv_datasize, mrequest_perproc, &
                            maxproc_cutman))

      do icut = 1,myncuts
        iblk1 = cutopo( 1,icut)
        iblk2 = cutopo(10,icut)
        call myblock(iblk1,myblock1,.true.)
        call myblock(iblk2,myblock2,.true.)
        if (myblock1) then
          imax1  = i_imax     (iblk1,nl)
          jmax1  = j_jmax     (iblk1,nl)
          kmax1  = k_kmax     (iblk1,nl)
          imut1    = 1 + off_p3 (iblk1,nl) * dim5
        else
          imut1    = 1
        end if
        if (myblock2) then
          imax2  = i_imax     (iblk2,nl)
          jmax2  = j_jmax     (iblk2,nl)
          kmax2  = k_kmax     (iblk2,nl)
          imut2    = 1 + off_p3 (iblk2,nl) * dim5
        else
          imut2    = 1
        end if
        call cut_q(imax1,jmax1,kmax1,mut(imut1),imax2,jmax2,kmax2, &
             mut(imut2),1,nharms,cutopo(1,icut),receivearray, &
             receiverequests,numrecv,receiveindices,sendarray, &
             sendrequests,numsend,imut1,mrequest_perproc,maxproc_cutman, &
             maxproc,myrecvid_num,mysendid_num)
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
           imut1 = receiveindices(i_in_id,4,reqnum)
           imax1 = receiveindices(i_in_id,5,reqnum)
           jmax1 = receiveindices(i_in_id,6,reqnum)
           kmax1 = receiveindices(i_in_id,7,reqnum)
           call copymutreceivedata(imax1,jmax1,kmax1,a,b,c,mut(imut1), &
                receivearray,i_in_id,reqnum)
         end do
      end do

      call waitallmessages(sendrequests,maxsendid_num)

      deallocate(sendarray, receivearray)

      return
      end

!-----------------------------------------------------------------------
      subroutine bc_res(nl,bctopo,q)
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) iblock,nl,imax,jmax,kmax,iblk,iq
      integer(kind=cosa_int) bctopo(*)
      real(kind=cosa_real) q(*)

      do iblock = 1,mynblocks
        imax   = i_imax     (iblock,nl)
        jmax   = j_jmax     (iblock,nl)
        kmax   = k_kmax     (iblock,nl)
        iblk   = 1 + off_bct(iblock,nl)
        iq     = 1 + off_p3 (iblock,nl) * npde * dim5
        call bbc_res(nl,bctopo(iblk),nbcs(iblock),q(iq), &
                     imax,jmax,kmax,npde,nharms)
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine bbc_res(nl,bctopo,nbcs,q,imax,jmax,kmax,npde,nharms)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,nbcs,npde,nharms
      integer(kind=cosa_int) nl,n,nh,i,j,ipde,ifictrow,ibc
      integer(kind=cosa_int) bctopo(10,nbcs)
      real(kind=cosa_real) &
          q    (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms)

!-----------------------------------------------------------------------
!     bc type                                                       mask
!
!     inviscid wall (all vars are extrapolated except normal vel.)    14
!     viscous wall                                                    15
!     subsonic or supersonic free-stream                              30
!     extrapolation                                                    4
!     symmetry wrt xy-plane                                           81
!     symmetry wrt xz-plane                                           82
!     symmetry wrt yz-plane                                           83
!-----------------------------------------------------------------------

      do ibc = 1,nbcs

        if (bctopo(1,ibc)/100.eq.15) then

          call bc_res_bnd(q,-1.d0,-1.d0,bctopo(1,ibc),imax,jmax,kmax, &
                          npde,nharms)

        else

          call bc_res_bnd(q, 1.d0, 1.d0,bctopo(1,ibc),imax,jmax,kmax, &
                          npde,nharms)

        end if

      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine bc_res_bnd(q,pm1,pm2,bctopo,imax,jmax,kmax,npde,nharms)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,npde,nharms
      integer(kind=cosa_int) bctopo(10)
      integer(kind=cosa_int) n,nh,ipde,i1,i2,i3,ijkmax(3),idir,istrt(3),iend(3), &
                bctyp,inrout,ibcpt,ibcpt2,ibcn,ibcn2,ibc,jbc,kbc,ibc2, &
                jbc2,kbc2,in,jn,kn,in2,jn2,kn2,ic1,ic2,ic3
      real (kind=cosa_real) &
          q    (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms)
      real (kind=cosa_real) pm1,pm2

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

      do n = 0,2*nharms
        do ipde=1,5

          do i3 = istrt(ic3),iend(ic3)
            do i2 = istrt(ic2),iend(ic2)

              ibc  = ibcpt *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
              jbc  = ibcpt *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
              kbc  = ibcpt *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

              ibc2 = ibcpt2*krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
              jbc2 = ibcpt2*krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
              kbc2 = ibcpt2*krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

              in   = ibcn  *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
              jn   = ibcn  *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
              kn   = ibcn  *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

              q(ibc,jbc,kbc,ipde,n) = pm1 * q(in ,jn ,kn ,ipde,n)

            end do
          end do

        end do
      end do

      if (kom.or.kom_bsl.or.kom_sst) then

        do n = 0,2*nharms
          do ipde=6,npde

            do i3 = istrt(ic3),iend(ic3)
              do i2 = istrt(ic2),iend(ic2)

                ibc  = ibcpt *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
                jbc  = ibcpt *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
                kbc  = ibcpt *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

                ibc2 = ibcpt2*krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
                jbc2 = ibcpt2*krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
                kbc2 = ibcpt2*krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

                in   = ibcn  *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
                jn   = ibcn  *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
                kn   = ibcn  *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)

                q(ibc,jbc,kbc,ipde,n) = pm2 * q(in ,jn ,kn ,ipde,n)

              end do
            end do

          end do
        end do

      end if

      return
      end

!-----------------------------------------------------------------------
      subroutine pcutman_init(percutopo)
!-----------------------------------------------------------------------
!     The subroutine initializes data exchange along periodic cuts.
!     Each rank of process I have to perform send/receive is associated
!     to a number ranging from 1 to maxpsendid_num/maxprecvid_num.
!     The inverse association is provided by another array.
!     More precisely, the results are:
!     * maxpsendid_num = number of processes which I have to send to 
!     * maxprecvid_num = number of processes which I have to receive from
!     * mypsendid(1:maxpsendid_num) ---> from number to rank
!     * myprecvid(1:maxprecvid_num) ---> from number to rank
!     * mypsendid_num(1:maxproc) ---> from rank to number
!     * myprecvid_num(1:maxproc) ---> from rank to number
!     Maybe send and receive arrays are always identical, I do not know...
!     Consider that, in order to prepare data exchange, some constants
!     have to be adequately set
!     * maxproc  ----> maximum number of MPI processes
!     * maxproc_pcutman ---> maximum number of MPI processes one process exchange data
!     * mrequest_pperproc ---> maximum number of exchanged points for each send/recv operation
!-----------------------------------------------------------------------

      use parallelutils
      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) iblk1,icut,iblk2
      integer(kind=cosa_int) percutopo(22,mynpcuts)
      logical myblock1, myblock2
      integer(kind=cosa_int) length(3),istr1(3),iend1(3),idir1,l,ic2,ic3
      integer(kind=cosa_int) recvfromid,sendtoid,i_id,rank,ierr,msize
      logical found
      integer(kind=cosa_int) maxsend(0:maxproc-1), maxrecv(0:maxproc-1)

      maxsend = 0
      maxrecv = 0

      call allocate_pcutman_init_memory()

      maxpsendid_num = 0
      maxprecvid_num = 0
      mypsendid = -1 ! array syntax
      myprecvid = -1 ! array syntax 
      myprecvid_num = -1 ! array syntax
      mypsendid_num = -1 ! array syntax

      do icut = 1,mynpcuts
         iblk1 = percutopo( 1,icut)
         iblk2 = percutopo(10,icut)
         idir1    = percutopo(2,icut)
         istr1(1) = percutopo(4,icut)
         iend1(1) = percutopo(5,icut)
         istr1(2) = percutopo(6,icut)
         iend1(2) = percutopo(7,icut)
         istr1(3) = percutopo(8,icut)
         iend1(3) = percutopo(9,icut)
         do l = 1, 3
            if (l .ne. idir1) then
               if (iend1(l) .gt. istr1(l)) then
                  iend1(l)  = iend1(l) + 1
               else
                  istr1(l)  = istr1(l) + 1
               end if
            end if            
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
            do i_id=1,maxprecvid_num
               if(myprecvid(i_id) .eq. recvfromid) found = .true. 
            enddo
            if(.not.found) then
              maxprecvid_num = maxprecvid_num + 1
              myprecvid(maxprecvid_num) = recvfromid
              myprecvid_num(recvfromid) = maxprecvid_num
            endif
            maxrecv(recvfromid) = maxrecv(recvfromid) + &
                                  (length(ic3)+1)*(length(ic2)+1)*3
         else if(myblock2) then
            call getownerid(iblk1,sendtoid)
!            print*,'iblk1: ',iblk1
            found = .false.
            do i_id=1,maxpsendid_num
               if(mypsendid(i_id) .eq. sendtoid) found = .true. 
            enddo
            if(.not.found) then
              maxpsendid_num = maxpsendid_num + 1
              mypsendid(maxpsendid_num) = sendtoid
              mypsendid_num(sendtoid) = maxpsendid_num
            endif
            maxsend(sendtoid) = maxsend(sendtoid) + &
                                (length(ic3)+1)*(length(ic2)+1)*3
         end if
      end do

      mrequest_pperproc = max(maxval(maxsend), maxval(maxrecv))
      maxproc_pcutman = max(maxpsendid_num,maxprecvid_num)
      call allocate_pcutman_memory()

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
      subroutine pcutman_q(nl,percutopo,q)
!-----------------------------------------------------------------------

      use parallelutils
      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,iblk1,icut,iblk2,imax1,jmax1,kmax1,imax2,jmax2,kmax2, &
                iq1,iq2
      integer(kind=cosa_int) percutopo(22,mynpcuts),a,b,c,i,reqnum
      logical myblock1, myblock2
      real(kind=cosa_real) q(*)
      integer(kind=cosa_int) sendrequestnum,receiverequestnum
      integer(kind=cosa_int) i_id,i_in_id
      integer(kind=cosa_int) sendid,recvid

      pnumrecv = 0  ! array syntax
      pnumsend = 0  ! array syntax
      sendrecv_datasize = npde*((2*nharms)+1)  

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
            iq1    = 1 + off_p3 (iblk1,nl) * npde * dim5
         else
            iq1    = 1
         end if
         if(myblock2) then
            imax2  = i_imax     (iblk2,nl)
            jmax2  = j_jmax     (iblk2,nl)
            kmax2  = k_kmax     (iblk2,nl)
            iq2    = 1 + off_p3 (iblk2,nl) * npde * dim5
         else
            iq2    = 1
         end if
         call pcut_q(imax1,jmax1,kmax1,q(iq1),imax2,jmax2,kmax2,q(iq2), &
              npde,nharms,percutopo(1,icut),receivearray,preceiverequests, &
              pnumrecv,preceiveindices,sendarray,psendrequests,pnumsend,iq1, &
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

      if (nharms.eq.0) then
         do i_id=1,maxprecvid_num
            call waitanymessages(preceiverequests,maxprecvid_num,reqnum)
            do i_in_id=1,pnumrecv(reqnum)
              a     = preceiveindices(i_in_id,1,reqnum)
              b     = preceiveindices(i_in_id,2,reqnum)
              c     = preceiveindices(i_in_id,3,reqnum)
              iq1   = preceiveindices(i_in_id,4,reqnum)
              imax1 = preceiveindices(i_in_id,5,reqnum)
              jmax1 = preceiveindices(i_in_id,6,reqnum)
              kmax1 = preceiveindices(i_in_id,7,reqnum)
              call copypqreceivedata(imax1,jmax1,kmax1,a,b,c,q(iq1), &
                   receivearray,i_in_id,reqnum)
            enddo
         enddo
      else
         do i_id=1,maxprecvid_num
            call waitanymessages(preceiverequests,maxprecvid_num,reqnum)
            do i_in_id=1,pnumrecv(reqnum)
              a     = preceiveindices(i_in_id,1,reqnum)
              b     = preceiveindices(i_in_id,2,reqnum)
              c     = preceiveindices(i_in_id,3,reqnum)
              iq1   = preceiveindices(i_in_id,4,reqnum)
              imax1 = preceiveindices(i_in_id,5,reqnum)
              jmax1 = preceiveindices(i_in_id,6,reqnum)
              kmax1 = preceiveindices(i_in_id,7,reqnum)
              call copypqhbreceivedata(imax1,jmax1,kmax1,a,b,c,q(iq1), &
                   receivearray,i_in_id,reqnum)
            enddo
         enddo
      end if

      call waitallmessages(psendrequests,maxpsendid_num)

      deallocate(sendarray, receivearray)

      return

      end

!-----------------------------------------------------------------------
      subroutine pcut_q(imax1,jmax1,kmax1,q1,imax2,jmax2,kmax2,q2,npde, &
        nharms,percutopo,receivearray,receiverequests,numrecv, &
        receiveindices,sendarray,sendrequests,numsend,iq1, &
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

      integer(kind=cosa_int) imax1,jmax1,kmax1,imax2,jmax2,kmax2,npde,nharms
      integer(kind=cosa_int) percutopo(22),mrequest,iblk1,iblk2
      integer(kind=cosa_int) ijkmax1(3),ijkmax2(3),istr1(3),iend1(3),istr2(3), &
           iend2(3),isgn1(3),isgn2(3),iord(3),irot,len(3),idir1,idir2, &
           inout1,inout2,ibcpt,ibcinc,ibcpt2,inr,inrinc,inr2,l,ipde,n, &
           n1,ic,ic1,ic2,ic3,jc1,jc2,jc3,i2,i3,ii1,jj1,kk1,ii2,jj2,kk2, &
           in1,jn1,kn1,in2,jn2,kn2
      real(kind=cosa_real) &
           q1(-1:imax1+1,-1:jmax1+1,-1:kmax1+1,npde,0:2*nharms), &
           q2(-1:imax2+1,-1:jmax2+1,-1:kmax2+1,npde,0:2*nharms)
      integer(kind=cosa_int) tempindex, datasize
      integer(kind=cosa_int) receiverequestnum,iq1,sendrequestnum
      logical myblock1, myblock2
      integer(kind=cosa_int) mrequest_perproc,maxproc_cutman
      integer(kind=cosa_int) sendid,recvid,numrecvid,numsendid,maxproc
      integer(kind=cosa_int) sendrequests(maxproc_cutman),receiverequests(maxproc_cutman)
      integer(kind=cosa_int) receiveindices(mrequest_perproc,8,maxproc_cutman)
      real(kind=cosa_real) &
        sendarray   (npde*(2*nharms+1),mrequest_perproc,maxproc_cutman), &
        receivearray(npde*(2*nharms+1),mrequest_perproc,maxproc_cutman)
      real(kind=cosa_real) tmp(7,0:2*mharms,4)
      integer(kind=cosa_int) numsend(maxproc_cutman),numrecv(maxproc_cutman)
      integer(kind=cosa_int) myrecvid_num(0:maxproc-1),mysendid_num(0:maxproc-1)
!--------------------------------------------------------------------------

!     Store imax, jmax, kmax in ijmax for looping
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
         ibcpt  =   0
         ibcinc = - 1
      else
         ibcpt  = ijkmax1(idir1)
         ibcinc = 1
      end if
!     
      if (inout2 .eq. 1) then
         inr    =   1
         inrinc =   1
      else
         inr    =   ijkmax2(idir2) - 1
         inrinc = - 1
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

          ii2 = ii1 + ibcinc * krd (ic1, 1)
          jj2 = jj1 + ibcinc * krd (ic1, 2)
          kk2 = kk1 + ibcinc * krd (ic1, 3)

          in1 = inr                          * krd (jc1, 1) + &
                (istr2(jc2) + isgn2(jc2)*i2) * krd (jc2, 1) + &
                (istr2(jc3) + isgn2(jc3)*i3) * krd (jc3, 1)
          jn1 = inr                          * krd (jc1, 2) + &
                (istr2(jc2) + isgn2(jc2)*i2) * krd (jc2, 2) + &
                (istr2(jc3) + isgn2(jc3)*i3) * krd (jc3, 2)
          kn1 = inr                          * krd (jc1, 3) + &
                (istr2(jc2) + isgn2(jc2)*i2) * krd (jc2, 3) + &
                (istr2(jc3) + isgn2(jc3)*i3) * krd (jc3, 3)

          in2 = in1 + inrinc * krd (jc1, 1)
          jn2 = jn1 + inrinc * krd (jc1, 2)
          kn2 = kn1 + inrinc * krd (jc1, 3)

          if (myblock1 .and. myblock2) then

            if (nharms.eq.0) then

!alt          n = 0
!alt          q1(ii1,jj1,kk1,   1,n) = q2(in1,jn1,kn1,   1,n)
!alt          q1(ii2,jj2,kk2,   1,n) = q2(in2,jn2,kn2,   1,n)
!alt          tmp(   2,n,1) = q2(in1,jn1,kn1,   2,n)
!alt          tmp(   3,n,1) = q2(in1,jn1,kn1,   3,n)
!alt          tmp(   2,n,2) = q2(in2,jn2,kn2,   2,n)
!alt          tmp(   3,n,2) = q2(in2,jn2,kn2,   3,n)
!alt          q1(ii1,jj1,kk1,2,n) =   rotmat(1,1) * tmp(   2,n,1) &
!alt                                - rotmat(1,2) * tmp(   3,n,1) * irot
!alt          q1(ii2,jj2,kk2,2,n) =   rotmat(1,1) * tmp(   2,n,2) &
!alt                                - rotmat(1,2) * tmp(   3,n,2) * irot
!alt          q1(ii1,jj1,kk1,3,n) = - rotmat(2,1) * tmp(   2,n,1) * irot &
!alt                                + rotmat(2,2) * tmp(   3,n,1)
!alt          q1(ii2,jj2,kk2,3,n) = - rotmat(2,1) * tmp(   2,n,2) * irot &
!alt                                + rotmat(2,2) * tmp(   3,n,2)
!alt          do ipde = 4, npde
!alt            q1(ii1,jj1,kk1,ipde,n) = q2(in1,jn1,kn1,ipde,n)
!alt            q1(ii2,jj2,kk2,ipde,n) = q2(in2,jn2,kn2,ipde,n)
!alt          end do

              n = 0
              do ipde = 1, npde
                q1(ii1,jj1,kk1,ipde,n) = q2(in1,jn1,kn1,ipde,n)
                q1(ii2,jj2,kk2,ipde,n) = q2(in2,jn2,kn2,ipde,n)
              end do
              tmp(   2,n,1) = q2(in1,jn1,kn1,   2,n)
              tmp(   3,n,1) = q2(in1,jn1,kn1,   3,n)
              tmp(   2,n,2) = q2(in2,jn2,kn2,   2,n)
              tmp(   3,n,2) = q2(in2,jn2,kn2,   3,n)
              q1(ii1,jj1,kk1,2,n) =   rotmat(1,1) * tmp(   2,n,1) &
                                    - rotmat(1,2) * tmp(   3,n,1) * irot
              q1(ii2,jj2,kk2,2,n) =   rotmat(1,1) * tmp(   2,n,2) &
                                    - rotmat(1,2) * tmp(   3,n,2) * irot
              q1(ii1,jj1,kk1,3,n) = - rotmat(2,1) * tmp(   2,n,1) * irot &
                                    + rotmat(2,2) * tmp(   3,n,1)
              q1(ii2,jj2,kk2,3,n) = - rotmat(2,1) * tmp(   2,n,2) * irot &
                                    + rotmat(2,2) * tmp(   3,n,2)

            else

!------------ from time-domain to frequency-domain
              do n = 0, 2*nharms
                do ipde=1,npde
                  tmp(ipde,n,1) = 0.d0
                  tmp(ipde,n,2) = 0.d0
                  do n1 = 0, 2*nharms
                    tmp(ipde,n,1) = tmp(ipde,n,1) + ehb(n,n1) * &
                                    q2(in1,jn1,kn1,ipde,n1)
                    tmp(ipde,n,2) = tmp(ipde,n,2) + ehb(n,n1) * &
                                    q2(in2,jn2,kn2,ipde,n1)
                  end do
                end do
              end do

              do ipde = 1, npde
                tmp(ipde,0,3) = tmp(ipde,0,1)
                tmp(ipde,0,4) = tmp(ipde,0,2)
              end do

              do n = 1, nharms
                do ipde=1,npde
                  do ic=1,2
                    tmp(ipde,2*n-1,ic+2) = &
                        tmp(ipde,2*n-1,ic) * cos(n*gtheta) &
                      - tmp(ipde,2*n  ,ic) * sin(n*gtheta*somega*irot)
                    tmp(ipde,2*n  ,ic+2) = &
                      + tmp(ipde,2*n-1,ic) * sin(n*gtheta*somega*irot) &
                      + tmp(ipde,2*n  ,ic) * cos(n*gtheta)
                  end do
                end do
              end do

              do n = 0, 2*nharms
                tmp(   2,n,1) = tmp(   2,n,3)
                tmp(   3,n,1) = tmp(   3,n,3)
                tmp(   2,n,2) = tmp(   2,n,4)
                tmp(   3,n,2) = tmp(   3,n,4)
                tmp(   2,n,3) = &
                  + rotmat(1,1) * tmp(   2,n,1) &
                  - rotmat(1,2) * tmp(   3,n,1) * irot
                tmp(   2,n,4) = &
                  + rotmat(1,1) * tmp(   2,n,2) &
                  - rotmat(1,2) * tmp(   3,n,2) * irot
                tmp(   3,n,3) = &
                  - rotmat(2,1) * tmp(   2,n,1) * irot &
                  + rotmat(2,2) * tmp(   3,n,1)
                tmp(   3,n,4) = &
                  - rotmat(2,1) * tmp(   2,n,2) * irot &
                  + rotmat(2,2) * tmp(   3,n,2)
              end do

!------------ from frequency-domain to time-domain
              do n = 0, 2*nharms
                do ipde=1,npde
                  q1(ii1,jj1,kk1,ipde,n) = 0.d0
                  q1(ii2,jj2,kk2,ipde,n) = 0.d0
                  do n1 = 0, 2*nharms
                    q1(ii1,jj1,kk1,ipde,n) = q1(ii1,jj1,kk1,ipde,n) + &
                                             eihb(n,n1) * tmp(ipde,n1,3)
                    q1(ii2,jj2,kk2,ipde,n) = q1(ii2,jj2,kk2,ipde,n) + &
                                             eihb(n,n1) * tmp(ipde,n1,4)
                  end do
                end do
              end do

            end if

          else if(myblock1) then

            call getownerid(iblk2,recvid)
            numrecvid = myrecvid_num(recvid)
            numrecv(numrecvid) = numrecv(numrecvid)+1 
            receiveindices(numrecv(numrecvid),1,numrecvid) = ii1
            receiveindices(numrecv(numrecvid),2,numrecvid) = jj1
            receiveindices(numrecv(numrecvid),3,numrecvid) = kk1
            receiveindices(numrecv(numrecvid),4,numrecvid) = iq1
            receiveindices(numrecv(numrecvid),5,numrecvid) = imax1
            receiveindices(numrecv(numrecvid),6,numrecvid) = jmax1
            receiveindices(numrecv(numrecvid),7,numrecvid) = kmax1
            numrecv(numrecvid) = numrecv(numrecvid)+1
            receiveindices(numrecv(numrecvid),1,numrecvid) = ii2
            receiveindices(numrecv(numrecvid),2,numrecvid) = jj2
            receiveindices(numrecv(numrecvid),3,numrecvid) = kk2
            receiveindices(numrecv(numrecvid),4,numrecvid) = iq1
            receiveindices(numrecv(numrecvid),5,numrecvid) = imax1
            receiveindices(numrecv(numrecvid),6,numrecvid) = jmax1
            receiveindices(numrecv(numrecvid),7,numrecvid) = kmax1
            if (numrecv(numrecvid).gt. mrequest_perproc) then
               write(*,*)'numrecv (pcut_q)--> INCREASE MREQUEST_PERPROC'
               call abortmpi()
            end if

          else if(myblock2) then

            call getownerid(iblk1,sendid)
            numsendid = mysendid_num(sendid)
            numsend(numsendid) = numsend(numsendid)+1 
            tempindex = 1

            if (nharms.eq.0) then

              n = 0
              sendarray(tempindex,numsend(numsendid),numsendid) = &
                q2(in1,jn1,kn1,   1,n)
              tempindex = tempindex + 1
              tmp(   2,n,1) = q2(in1,jn1,kn1,   2,n)
              tmp(   3,n,1) = q2(in1,jn1,kn1,   3,n)
              sendarray(tempindex,numsend(numsendid),numsendid) = &
                + rotmat(1,1) * tmp(   2,n,1) &
                - rotmat(1,2) * tmp(   3,n,1) * irot
              tempindex = tempindex + 1
              sendarray(tempindex,numsend(numsendid),numsendid) = &
                - rotmat(2,1) * tmp(   2,n,1) * irot &
                + rotmat(2,2) * tmp(   3,n,1)
              tempindex = tempindex + 1
              do ipde = 4, npde
                sendarray(tempindex,numsend(numsendid),numsendid) = &
                  q2(in1,jn1,kn1,ipde,n)
                tempindex = tempindex + 1
              end do

              numsend(numsendid) = numsend(numsendid)+1
              tempindex  = 1
 
              sendarray(tempindex,numsend(numsendid),numsendid) = &
                q2(in2,jn2,kn2,   1,n)
              tempindex = tempindex + 1
              tmp(   2,n,1) = q2(in2,jn2,kn2,   2,n)
              tmp(   3,n,1) = q2(in2,jn2,kn2,   3,n)
              sendarray(tempindex,numsend(numsendid),numsendid) = &
                + rotmat(1,1) * tmp(   2,n,1) &
                - rotmat(1,2) * tmp(   3,n,1) * irot
              tempindex = tempindex + 1
              sendarray(tempindex,numsend(numsendid),numsendid) = &
                - rotmat(2,1) * tmp(   2,n,1) * irot &
                + rotmat(2,2) * tmp(   3,n,1)
              tempindex = tempindex + 1
              do ipde = 4, npde
                sendarray(tempindex,numsend(numsendid),numsendid) = &
                  q2(in2,jn2,kn2,ipde,n)
                tempindex = tempindex + 1
              end do

            else
!------------ from time-domain to frequency-domain
              do n = 0, 2*nharms
                do ipde=1,npde
                  tmp(ipde,n,1) = 0.d0
                  tmp(ipde,n,2) = 0.d0
                  do n1 = 0, 2*nharms
                    tmp(ipde,n,1) = tmp(ipde,n,1) + ehb(n,n1) * &
                                    q2(in1,jn1,kn1,ipde,n1)
                    tmp(ipde,n,2) = tmp(ipde,n,2) + ehb(n,n1) * &
                                    q2(in2,jn2,kn2,ipde,n1)
                  end do
                end do
              end do

              do ipde = 1, npde
                tmp(ipde,0,3) = tmp(ipde,0,1)
                tmp(ipde,0,4) = tmp(ipde,0,2)
              end do

              do n = 1, nharms
                do ipde=1,npde
                  do ic=1,2
                    tmp(ipde,2*n-1,ic+2) = &
                        tmp(ipde,2*n-1,ic) * cos(n*gtheta) &
                      - tmp(ipde,2*n  ,ic) * sin(n*gtheta*somega*irot)
                    tmp(ipde,2*n  ,ic+2) = &
                      + tmp(ipde,2*n-1,ic) * sin(n*gtheta*somega*irot) &
                      + tmp(ipde,2*n  ,ic) * cos(n*gtheta)
                  end do
                end do
              end do

              do n = 0, 2*nharms
                tmp(   2,n,1) = tmp(   2,n,3)
                tmp(   3,n,1) = tmp(   3,n,3)
                tmp(   2,n,2) = tmp(   2,n,4)
                tmp(   3,n,2) = tmp(   3,n,4)
                tmp(   2,n,3) = &
                  + rotmat(1,1) * tmp(   2,n,1) &
                  - rotmat(1,2) * tmp(   3,n,1) * irot
                tmp(   2,n,4) = &
                  + rotmat(1,1) * tmp(   2,n,2) &
                  - rotmat(1,2) * tmp(   3,n,2) * irot
                tmp(   3,n,3) = &
                  - rotmat(2,1) * tmp(   2,n,1) * irot &
                  + rotmat(2,2) * tmp(   3,n,1)
                tmp(   3,n,4) = &
                  - rotmat(2,1) * tmp(   2,n,2) * irot &
                  + rotmat(2,2) * tmp(   3,n,2)
              end do

              do n = 0, 2*nharms
                do ipde = 1, npde
                  sendarray(tempindex,numsend(numsendid),numsendid) = &
                    tmp(ipde,n,3)
                  tempindex = tempindex + 1
                end do
              end do

              numsend(numsendid) = numsend(numsendid)+1
              tempindex  = 1
              do n = 0, 2*nharms
                do ipde = 1, npde
                  sendarray(tempindex,numsend(numsendid),numsendid) = &
                    tmp(ipde,n,4)
                  tempindex = tempindex + 1
                end do
              end do

            end if

            if (numsend(numsendid).gt. mrequest_perproc) then
              write(*,*)'numsend (pcut_q)--> INCREASE MREQUEST_PERPROC'
              call abortmpi()
            end if

          end if

        end do
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine pcutman_mut(nl,percutopo,mut)
!-----------------------------------------------------------------------

      use parallelutils
      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,iblk1,icut,iblk2,imax1,jmax1,kmax1,imax2,jmax2,kmax2, &
                imut1,imut2
      integer(kind=cosa_int) percutopo(22,mynpcuts),a,b,c,i,reqnum
      logical myblock1, myblock2
      real(kind=cosa_real) mut(*)
      integer(kind=cosa_int) sendrequestnum,receiverequestnum
      integer(kind=cosa_int) i_id,i_in_id
      integer(kind=cosa_int) sendid,recvid

      pnumrecv = 0  ! array syntax
      pnumsend = 0  ! array syntax
      sendrecv_datasize = 2*nharms+1

      allocate(sendarray(sendrecv_datasize, mrequest_pperproc, &
                         maxproc_pcutman))
      allocate(receivearray(sendrecv_datasize, mrequest_pperproc, &
                            maxproc_pcutman))

      do icut = 1,mynpcuts
        iblk1 = percutopo( 1,icut)
        iblk2 = percutopo(10,icut)
        call myblock(iblk1,myblock1,.true.)
        call myblock(iblk2,myblock2,.true.)
        if (myblock1) then
          imax1  = i_imax     (iblk1,nl)
          jmax1  = j_jmax     (iblk1,nl)
          kmax1  = k_kmax     (iblk1,nl)
          imut1    = 1 + off_p3 (iblk1,nl) * dim5
        else
          imut1    = 1
        end if
        if (myblock2) then
          imax2  = i_imax     (iblk2,nl)
          jmax2  = j_jmax     (iblk2,nl)
          kmax2  = k_kmax     (iblk2,nl)
          imut2    = 1 + off_p3 (iblk2,nl) * dim5
        else
          imut2    = 1
        end if
        call pcut_mut(imax1,jmax1,kmax1,mut(imut1),imax2,jmax2,kmax2, &
             mut(imut2),1,nharms,percutopo(1,icut),receivearray, &
             preceiverequests,pnumrecv,preceiveindices,sendarray, &
             psendrequests,pnumsend,imut1,mrequest_pperproc,maxproc_pcutman, &
             maxproc,myprecvid_num,mypsendid_num)
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
           imut1 = preceiveindices(i_in_id,4,reqnum)
           imax1 = preceiveindices(i_in_id,5,reqnum)
           jmax1 = preceiveindices(i_in_id,6,reqnum)
           kmax1 = preceiveindices(i_in_id,7,reqnum)
           call copypmutreceivedata(imax1,jmax1,kmax1,a,b,c,mut(imut1), &
                receivearray,i_in_id,reqnum)
         end do
      end do

      call waitallmessages(psendrequests,maxpsendid_num)

      deallocate(sendarray, receivearray)

      return
      end

!-----------------------------------------------------------------------
      subroutine pcut_mut(imax1,jmax1,kmax1,q1,imax2,jmax2,kmax2,q2,npde, &
        nharms,percutopo,receivearray,receiverequests,numrecv, &
        receiveindices,sendarray,sendrequests,numsend,iq1, &
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

      integer(kind=cosa_int) imax1,jmax1,kmax1,imax2,jmax2,kmax2,npde,nharms
      integer(kind=cosa_int) percutopo(22),mrequest,iblk1,iblk2
      integer(kind=cosa_int) ijkmax1(3),ijkmax2(3),istr1(3),iend1(3),istr2(3), &
           iend2(3),isgn1(3),isgn2(3),iord(3),irot,len(3),idir1,idir2, &
           inout1,inout2,ibcpt,ibcinc,ibcpt2,inr,inrinc,inr2,l,ipde,n, &
           ic1,ic2,ic3,jc1,jc2,jc3,i2,i3,ii1,jj1,kk1,ii2,jj2,kk2,in1, &
           jn1,kn1,in2,jn2,kn2
      real(kind=cosa_real) &
           q1(-1:imax1+1,-1:jmax1+1,-1:kmax1+1,npde,0:2*nharms), &
           q2(-1:imax2+1,-1:jmax2+1,-1:kmax2+1,npde,0:2*nharms)
      integer(kind=cosa_int) tempindex, datasize
      integer(kind=cosa_int) receiverequestnum,iq1,sendrequestnum
      logical myblock1, myblock2
      integer(kind=cosa_int) mrequest_perproc,maxproc_cutman
      integer(kind=cosa_int) sendid,recvid,numrecvid,numsendid,maxproc
      integer(kind=cosa_int) sendrequests(maxproc_cutman),receiverequests(maxproc_cutman)
      integer(kind=cosa_int) receiveindices(mrequest_perproc,8,maxproc_cutman)
      real(kind=cosa_real) &
        sendarray   (npde*(2*nharms+1),mrequest_perproc,maxproc_cutman), &
        receivearray(npde*(2*nharms+1),mrequest_perproc,maxproc_cutman)
      real(kind=cosa_real) tmp(2)
      integer(kind=cosa_int) numsend(maxproc_cutman),numrecv(maxproc_cutman)
      integer(kind=cosa_int) myrecvid_num(0:maxproc-1),mysendid_num(0:maxproc-1)
!--------------------------------------------------------------------------

!     Store imax, jmax, kmax in ijmax for looping
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
         ibcpt  =   0
         ibcinc = - 1
      else
         ibcpt  = ijkmax1(idir1)
         ibcinc = 1
      end if
!     
      if (inout2 .eq. 1) then
         inr    =   1
         inrinc =   1
      else
         inr    =   ijkmax2(idir2) - 1
         inrinc = - 1
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

          ii2 = ii1 + ibcinc * krd (ic1, 1)
          jj2 = jj1 + ibcinc * krd (ic1, 2)
          kk2 = kk1 + ibcinc * krd (ic1, 3)

          in1 = inr                          * krd (jc1, 1) + &
                (istr2(jc2) + isgn2(jc2)*i2) * krd (jc2, 1) + &
                (istr2(jc3) + isgn2(jc3)*i3) * krd (jc3, 1)
          jn1 = inr                          * krd (jc1, 2) + &
                (istr2(jc2) + isgn2(jc2)*i2) * krd (jc2, 2) + &
                (istr2(jc3) + isgn2(jc3)*i3) * krd (jc3, 2)
          kn1 = inr                          * krd (jc1, 3) + &
                (istr2(jc2) + isgn2(jc2)*i2) * krd (jc2, 3) + &
                (istr2(jc3) + isgn2(jc3)*i3) * krd (jc3, 3)

          in2 = in1 + inrinc * krd (jc1, 1)
          jn2 = jn1 + inrinc * krd (jc1, 2)
          kn2 = kn1 + inrinc * krd (jc1, 3)

          if (myblock1 .and. myblock2) then

            do n = 0, 2*nharms
              do ipde = 1, npde
                q1(ii1,jj1,kk1,ipde,n) = q2(in1,jn1,kn1,ipde,n)
                q1(ii2,jj2,kk2,ipde,n) = q2(in2,jn2,kn2,ipde,n)
              end do
            end do

          else if(myblock1) then

            call getownerid(iblk2,recvid)
            numrecvid = myrecvid_num(recvid)
            numrecv(numrecvid) = numrecv(numrecvid)+1 
            receiveindices(numrecv(numrecvid),1,numrecvid) = ii1
            receiveindices(numrecv(numrecvid),2,numrecvid) = jj1
            receiveindices(numrecv(numrecvid),3,numrecvid) = kk1
            receiveindices(numrecv(numrecvid),4,numrecvid) = iq1
            receiveindices(numrecv(numrecvid),5,numrecvid) = imax1
            receiveindices(numrecv(numrecvid),6,numrecvid) = jmax1
            receiveindices(numrecv(numrecvid),7,numrecvid) = kmax1
            numrecv(numrecvid) = numrecv(numrecvid)+1
            receiveindices(numrecv(numrecvid),1,numrecvid) = ii2
            receiveindices(numrecv(numrecvid),2,numrecvid) = jj2
            receiveindices(numrecv(numrecvid),3,numrecvid) = kk2
            receiveindices(numrecv(numrecvid),4,numrecvid) = iq1
            receiveindices(numrecv(numrecvid),5,numrecvid) = imax1
            receiveindices(numrecv(numrecvid),6,numrecvid) = jmax1
            receiveindices(numrecv(numrecvid),7,numrecvid) = kmax1
            if (numrecv(numrecvid).gt. mrequest_perproc) then
               write(*,*) &
               'numrecv (pcut_mut)--> INCREASE MREQUEST_PERPROC'
               call abortmpi()
            end if

          else if(myblock2) then

            call getownerid(iblk1,sendid)
            numsendid = mysendid_num(sendid)
            numsend(numsendid) = numsend(numsendid)+1 
            tempindex = 1
            do n = 0, 2*nharms
              do ipde = 1, npde
                sendarray(tempindex,numsend(numsendid),numsendid) = &
                  q2(in1,jn1,kn1,ipde,n)
                tempindex = tempindex + 1
              end do
            end do
            numsend(numsendid) = numsend(numsendid)+1
            tempindex  = 1
            do n = 0, 2*nharms
              do ipde = 1, npde
                sendarray(tempindex,numsend(numsendid),numsendid) = &
                  q2(in2,jn2,kn2,ipde,n)
                tempindex = tempindex + 1
              end do
            end do
            if (numsend(numsendid).gt. mrequest_perproc) then
               write(*,*) &
               'numsend (pcut_mut)--> INCREASE MREQUEST_PERPROC'
               call abortmpi()
            end if

          end if

        end do
      end do

      return
      end
