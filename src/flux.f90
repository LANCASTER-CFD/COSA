!-----------------------------------------------------------------------
      subroutine roflux(idir,nl,qr,ql,alphar,work,flux,sijk)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_variables
      use cosa_precision

      implicit none


      integer(kind=cosa_int) idir,nl,iblock,imax,jmax,kmax,ic,jc,kc,iqlr,iimt, &
                iflux
      real (kind=cosa_real) qr(*),ql(*),sijk(*),alphar(*),work(*),flux(*)

      do iblock = 1,mynblocks
        imax   = i_imax     (iblock,nl)
        jmax   = j_jmax     (iblock,nl)
        kmax   = k_kmax     (iblock,nl)
        if (idir.eq.1) then
          ic = imax
          jc = jmax - 1
          kc = kmax - 1
        else if (idir.eq.2) then
          ic = imax - 1
          jc = jmax
          kc = kmax - 1
        else if (idir.eq.3) then
          ic = imax - 1
          jc = jmax - 1
          kc = kmax
        end if
        iqlr   = 1 + off_p3 (iblock,nl) * npde * dim5
        iimt   = 1 + off_p2 (iblock,nl) * lmet * dim5h
        iflux  = 1 + off_0  (iblock,nl) * npde * dim5
        call roflux_b(qr(iqlr),ql(iqlr),alphar(iqlr),work(iqlr), &
          flux(iflux),sijk(iimt),imax,jmax,kmax,npde,nharms,lmet,ic,jc, &
          kc)
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine roflux_b(qr,ql,alphar,work,flux,sijk,imax,jmax,kmax, &
                          npde,nharms,lmet,ic,jc,kc)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none


      integer(kind=cosa_int) imax,jmax,kmax,npde,nharms,lmet
      integer(kind=cosa_int) ic,jc,kc,i,j,k,n,nh
      real (kind=cosa_real) kx,ky,kz,ds,ubn, &
           rhor,ur,vr,wr,q2r,pr,unr,h0r,rh0r,tker,omer, &
           rhol,ul,vl,wl,q2l,pl,unl,h0l,rh0l,tkel,omel, &
           drho,du,dv,dw,dp,dun,dtke,dome, &
           rho,u,v,w,un,h0,q2,a2,a,tke,ome,rm123,rm4,rm5, &
           rhols,rhors,denom,rdelta,fr(5),fl(5),df(5)
      real (kind=cosa_real) &
           qr    (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           ql    (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           alphar(-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           work  (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           flux  (   imax  ,   jmax  ,   kmax  ,npde,0:2*nharms), &
           sijk  (lmet, 0:imax+1, 0:jmax+1, 0:kmax+1,0:2*nharms*hbmove)
      
      ubn = 0.d0

      do n = 0,2*nharms
        nh = n*hbmove

!------ store right and left turbulent kinetic energy in
!       work(:,:1) and work(:,:,2) respectively

        if (kom.or.kom_bsl.or.kom_sst) then
          do k = 1,kc
            do j = 1,jc
              do i = 1,ic
                work(i,j,k,1,n) = qr(i,j,k,6,n)
                work(i,j,k,2,n) = ql(i,j,k,6,n)
              end do
            end do
          end do
        else
          do k = 1,kc
            do j = 1,jc
              do i = 1,ic
                work(i,j,k,1,n) = 0.d0
                work(i,j,k,2,n) = 0.d0
              end do
            end do
          end do
        end if

!------ calculate flux of mean flow variables

        do k = 1,kc
          do j = 1,jc
            do i = 1,ic

              kx = sijk(1,i,j,k,nh)
              ky = sijk(2,i,j,k,nh)
              kz = sijk(3,i,j,k,nh)
              ds = sijk(4,i,j,k,nh)
              if (hawt.and.rotating.and.(.not.tpitch)) then
                ubn  = sijk(8,i,j,k,nh)
              else if (moving.or.(hawt.and.rotating.and.tpitch)) then
                ubn  = sijk(5,i,j,k,nh) * kx + sijk(6,i,j,k,nh) * ky + &
                       sijk(7,i,j,k,nh) * kz
              end if

!------------ calculate fr using qr vector
              rhor   = qr(i,j,k,1,n)
              ur     = qr(i,j,k,2,n)
              vr     = qr(i,j,k,3,n)
              wr     = qr(i,j,k,4,n)
              pr     = qr(i,j,k,5,n)
              unr    = ur*kx + vr*ky + wr*kz - ubn
              q2r    = ur**2 + vr**2 + wr**2
              tker   = work(i,j,k,1,n)
              h0r    = gamma/(gamma-1)*pr/rhor + q2r/2 + tker
              rh0r   = rhor * h0r
              fr(1)  = rhor*unr
              fr(2)  = (rhor*ur*unr+pr*kx)
              fr(3)  = (rhor*vr*unr+pr*ky)
              fr(4)  = (rhor*wr*unr+pr*kz)
              fr(5)  = (rh0r * unr + pr * ubn)   

!------------ calculate fl using ql vector
              rhol   = ql(i,j,k,1,n)
              ul     = ql(i,j,k,2,n)
              vl     = ql(i,j,k,3,n)
              wl     = ql(i,j,k,4,n)
              pl     = ql(i,j,k,5,n)
              unl    = ul*kx + vl*ky + wl*kz - ubn
              q2l    = ul**2 + vl**2 + wl**2
              tkel   = work(i,j,k,2,n)
              h0l    = gamma/(gamma-1)*pl/rhol + q2l/2 + tkel
              rh0l   = rhol * h0l
              fl(1)  = rhol*unl
              fl(2)  = (rhol*ul*unl+pl*kx)
              fl(3)  = (rhol*vl*unl+pl*ky)
              fl(4)  = (rhol*wl*unl+pl*kz)
              fl(5)  = (rh0l * unl + pl * ubn)   

!------------ calculate delta vector
              drho   = rhor - rhol
              dp     = pr - pl
              du     = ur - ul
              dv     = vr - vl
              dw     = wr - wl
              dun    = unr-unl
              dtke   = tker - tkel

!------------ calculate Roe average values
              rhols  = sqrt(rhol)
              rhors  = sqrt(rhor)
              rho    = rhols*rhors
              denom  = rhols+rhors
              u      = (rhols*ul+rhors*ur)/denom
              v      = (rhols*vl+rhors*vr)/denom
              w      = (rhols*wl+rhors*wr)/denom
              un     = u*kx + v*ky + w*kz
              h0     = (rhols*h0l+rhors*h0r)/denom
              tke    = (rhols*tkel+rhors*tker)/denom
              q2     = u**2 + v**2 + w**2
              a2     = (gamma-1)*(h0-q2/2-tke)
              a      = sqrt(a2)

!------------ eigenvalues
              rm123  = abs(un     - ubn)
              rm4    = abs(un + a - ubn)
              rm5    = abs(un - a - ubn)

              if (efx_setup_f) then
!-------------- Entropy Fix based on Yee, Warming, and Harten (1985) and 
!               Yee (). Apply fix only to acoustic eigenvalues to preserve 
!               boundary layer calculations. Avoidance of the carbuncle 
!               phenomena requires fix applied to convective eigenvalue 
!               also. Details based on conversations with Jeff White.
!               Begin {N.B. Comment out this section to remove the entropy fix
                rdelta = cntrpy_f * &
                       ( entfx_f_w(1) * entfxctf_f + &
                         entfx_f_w(2) * (sqrt(u**2+v**2+w**2)+a) + &
                         entfx_f_w(3) * 0.2d0*(abs(un)+a+abs(ubn)) )
                rdelta = max(rdelta,1.d-12)
                if ((rm123.le.rdelta).and.efx_all_f) then
                  rm123 = 0.5d0*(rm123**2 + rdelta**2) / rdelta
                end if
                if (rm4.le.rdelta) then
                  rm4  = 0.5d0*(rm4**2  + rdelta**2) / rdelta
                end if
                if (rm5.le.rdelta) then
                  rm5  = 0.5d0*(rm5**2  + rdelta**2) / rdelta
                end if
!               End} Entropy Fix
              end if

              if (eig_limit_f) then
                if (rm123.le.eig_cutoff_f) rm123 = eig_cutoff_f
                if (eig_lim_all_f.and.(rm4.le.eig_cutoff_f)) &
                  rm4   = eig_cutoff_f
                if (eig_lim_all_f.and.(rm5.le.eig_cutoff_f)) &
                  rm5   = eig_cutoff_f
              end if

!------------ deltas of characteristic variables
              alphar(i,j,k,1,n) = rm123 * (drho - dp/a2)
              alphar(i,j,k,2,n) = rm123 * rho
              alphar(i,j,k,3,n) = rm4 * (dp/a2 + rho*dun/a) / 2
              alphar(i,j,k,4,n) = rm5 * (dp/a2 - rho*dun/a) / 2

!------------ calculate flux differences
              df(1) = alphar(i,j,k,1,n) + alphar(i,j,k,3,n) + &
                      alphar(i,j,k,4,n)
              df(2) = alphar(i,j,k,1,n) * u + &
                      alphar(i,j,k,2,n) * (du - dun*kx) + &
                      alphar(i,j,k,3,n) * ( u + a  *kx) + &
                      alphar(i,j,k,4,n) * ( u - a  *kx)
              df(3) = alphar(i,j,k,1,n) * v + &
                      alphar(i,j,k,2,n) * (dv - dun*ky) + &
                      alphar(i,j,k,3,n) * ( v + a  *ky) + &
                      alphar(i,j,k,4,n) * ( v - a  *ky)
              df(4) = alphar(i,j,k,1,n) * w + &
                      alphar(i,j,k,2,n) * (dw - dun*kz) + &
                      alphar(i,j,k,3,n) * ( w + a  *kz) + &
                      alphar(i,j,k,4,n) * ( w - a  *kz)
              df(5) = alphar(i,j,k,1,n) * (q2/2 + tke) + &
                      alphar(i,j,k,2,n) * &
                        (u*du + v*dv + w*dw - un*dun + dtke) + &
                      alphar(i,j,k,3,n) * ( h0 + a  *un) + &
                      alphar(i,j,k,4,n) * ( h0 - a  *un)

!------------ calculate total flux at interface
              flux(i,j,k,1,n) = ( fr(1) + fl(1) - df(1) ) / 2 * ds
              flux(i,j,k,2,n) = ( fr(2) + fl(2) - df(2) ) / 2 * ds
              flux(i,j,k,3,n) = ( fr(3) + fl(3) - df(3) ) / 2 * ds
              flux(i,j,k,4,n) = ( fr(4) + fl(4) - df(4) ) / 2 * ds
              flux(i,j,k,5,n) = ( fr(5) + fl(5) - df(5) ) / 2 * ds

            end do
          end do
        end do

!------ calculate flux of turbulent variables

        if (kom.or.kom_bsl.or.kom_sst) then

          do k = 1,kc
            do j = 1,jc
              do i = 1,ic

                kx = sijk(1,i,j,k,nh)
                ky = sijk(2,i,j,k,nh)
                kz = sijk(3,i,j,k,nh)
                ds = sijk(4,i,j,k,nh)
                if (hawt.and.rotating.and.(.not.tpitch)) then
                  ubn  = sijk(8,i,j,k,nh)
                else if (moving.or.(hawt.and.rotating.and.tpitch)) then
                  ubn  = sijk(5,i,j,k,nh) * kx + sijk(6,i,j,k,nh) * ky + &
                         sijk(7,i,j,k,nh) * kz
                end if

!-------------- calculate fr using qr vector
                rhor   = qr(i,j,k,1,n)
                ur     = qr(i,j,k,2,n)
                vr     = qr(i,j,k,3,n)
                wr     = qr(i,j,k,4,n)
                pr     = qr(i,j,k,5,n)
                omer   = qr(i,j,k,7,n)
                unr    = ur*kx + vr*ky + wr*kz - ubn
                tker   = work(i,j,k,1,n)
                fr(1)  = rhor * unr * tker
                fr(2)  = rhor * unr * omer

!-------------- calculate fl using ql vector
                rhol   = ql(i,j,k,1,n)
                ul     = ql(i,j,k,2,n)
                vl     = ql(i,j,k,3,n)
                wl     = ql(i,j,k,4,n)
                pl     = ql(i,j,k,5,n)
                omel   = ql(i,j,k,7,n)
                unl    = ul*kx + vl*ky + wl*kz - ubn
                tkel   = work(i,j,k,2,n)
                fl(1)  = rhol * unl * tkel
                fl(2)  = rhol * unl * omel

!azz start
!------------ calculate delta vector
              drho   = rhor - rhol
              dp     = pr - pl
              du     = ur - ul
              dv     = vr - vl
              dw     = wr - wl
              dun    = unr-unl
              dtke   = tker - tkel
              dome   = omer - omel

!------------ calculate Roe average values
              rhols  = sqrt(rhol)
              rhors  = sqrt(rhor)
              rho    = rhols*rhors
              denom  = rhols+rhors
              u      = (rhols*ul+rhors*ur)/denom
              v      = (rhols*vl+rhors*vr)/denom
              w      = (rhols*wl+rhors*wr)/denom
              un     = u*kx + v*ky + w*kz
              h0     = (rhols*h0l+rhors*h0r)/denom
              tke    = (rhols*tkel+rhors*tker)/denom
              ome    = (rhols*omel+rhors*omer)/denom
              q2     = u**2 + v**2 + w**2
              a2     = (gamma-1)*(h0-q2/2-tke)
              a      = sqrt(a2)

!------------ eigenvalues
              rm123  = abs(un     - ubn)
              rm4    = abs(un + a - ubn)
              rm5    = abs(un - a - ubn)

              if (efx_setup_t) then
!-------------- Entropy Fix based on Yee, Warming, and Harten (1985) and 
!               Yee (). Apply fix only to acoustic eigenvalues to preserve 
!               boundary layer calculations. Avoidance of the carbuncle 
!               phenomena requires fix applied to convective eigenvalue 
!               also. Details based on conversations with Jeff White.
!               Begin {N.B. Comment out this section to remove the entropy fix
                rdelta = cntrpy_t * &
                       ( entfx_t_w(1) * entfxctf_t + &
                         entfx_t_w(2) * (sqrt(u**2+v**2+w**2)+a) + &
                         entfx_t_w(3) * 0.2d0*(abs(un)+a+abs(ubn)) )
                rdelta = max(rdelta,1.d-12)
                if ((rm123.le.rdelta).and.efx_all_t) then
                  rm123 = 0.5d0*(rm123**2 + rdelta**2) / rdelta
                end if
                if (rm4.le.rdelta) then
                  rm4  = 0.5d0*(rm4**2  + rdelta**2) / rdelta
                end if
                if (rm5.le.rdelta) then
                  rm5  = 0.5d0*(rm5**2  + rdelta**2) / rdelta
                end if
!               End} Entropy Fix
              end if

              if (eig_limit_t) then
                if (rm123.le.eig_cutoff_t) rm123 = eig_cutoff_t
                if (eig_lim_all_t.and.(rm4.le.eig_cutoff_t)) &
                  rm4   = eig_cutoff_t
                if (eig_lim_all_t.and.(rm5.le.eig_cutoff_t)) &
                  rm5   = eig_cutoff_t
              end if

!------------ deltas of characteristic variables
              alphar(i,j,k,1,n) = rm123 * (drho - dp/a2)
              alphar(i,j,k,2,n) = rm123 * rho
              alphar(i,j,k,3,n) = rm4 * (dp/a2 + rho*dun/a) / 2
              alphar(i,j,k,4,n) = rm5 * (dp/a2 - rho*dun/a) / 2

!------------ calculate flux differences

!azz end

!-------------- calculate flux differences
                df(1)  = ( alphar(i,j,k,1,n) + alphar(i,j,k,3,n) + &
                           alphar(i,j,k,4,n) ) * tke + &
                         alphar(i,j,k,2,n) * dtke
                df(2)  = ( alphar(i,j,k,1,n) + alphar(i,j,k,3,n) + &
                           alphar(i,j,k,4,n) ) * ome + &
                         alphar(i,j,k,2,n) * dome

!-------------- calculate total flux at interface
                flux(i,j,k,6,n) = (fr(1)+fl(1)-df(1)) / 2 * ds
                flux(i,j,k,7,n) = (fr(2)+fl(2)-df(2)) / 2 * ds

              end do
            end do
          end do

        end if

      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine proflux(idir,nl,qr,ql,alphar,coff_tke,flux,sijk)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) idir,nl,iblock,imax,jmax,kmax,ic,jc,kc, &
                iqlr,iimt,iflux,icof
      real (kind=cosa_real) qr(*),ql(*),coff_tke(*),sijk(*),flux(*),alphar(*)

      do iblock = 1,mynblocks
        imax   = i_imax     (iblock,nl)
        jmax   = j_jmax     (iblock,nl)
        kmax   = k_kmax     (iblock,nl)
        if (idir.eq.1) then
          ic = imax
          jc = jmax - 1
          kc = kmax - 1
        else if (idir.eq.2) then
          ic = imax - 1
          jc = jmax
          kc = kmax - 1
        else if (idir.eq.3) then
          ic = imax - 1
          jc = jmax - 1
          kc = kmax
        end if
        iqlr   = 1 + off_p3 (iblock,nl) * npde * dim5
        iimt   = 1 + off_p2 (iblock,nl) * lmet * dim5h
        iflux  = 1 + off_0  (iblock,nl) * npde * dim5
        call proflux_b(nl,qr(iqlr),ql(iqlr),alphar(iqlr),coff_tke(iqlr), &
                       flux(iflux), &
                       sijk(iimt),imax,jmax,kmax, &
                       npde,nharms,lmet,ic,jc,kc)
      end do


      return
      end

!-----------------------------------------------------------------------
      subroutine proflux_b(nl,qr,ql,alphar,coff_tke,flux,sijk, &
                           imax,jmax,kmax,npde,nharms,lmet,ic,jc,kc)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,npde,nharms,lmet
      integer(kind=cosa_int) ic,jc,kc,i,j,k,n,nh,nl
      real (kind=cosa_real) kx,ky,kz,ds,ubn, &
           rhor,ur,vr,wr,q2r,pr,thethr,h0r,rh0r,tker,omer, &
           rhol,ul,vl,wl,q2l,pl,thethl,h0l,rh0l,tkel,omel, &
           drho,dp,du,dv,dw,dtheth,dtke,dome,rhols,rhors,denom, &
           rho,u,v,w,theth,h0,q2,a2,a,asos,tke,ome,rm123,rm4,rm5,qp2, &
           m2,mp2,mp2s,rdelta, &
           lambdap4, lambdap5, esse, erre, ump2, fr(5), fl(5),df(5)
      real (kind=cosa_real) &
           qr      (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           ql      (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           alphar  (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           coff_tke(-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           flux    (     imax,     jmax,     kmax,npde,0:2*nharms), &
           sijk    (lmet,0:imax+1,0:jmax+1,0:kmax+1,0:2*nharms*hbmove)

      ubn = 0.d0

      do n = 0,2*nharms
        nh = n*hbmove

!------ store right and left turbulent kinetic energy in
!       coff_tke(:,:3) and coff_tke(:,:,4) respectively

        if (kom.or.kom_bsl.or.kom_sst) then
          do k = 1,kc
            do j = 1,jc
              do i = 1,ic
                coff_tke(i,j,k,3,n) = qr(i,j,k,6,n)
                coff_tke(i,j,k,4,n) = ql(i,j,k,6,n)
              end do
            end do
          end do
        else
          do k = 1,kc
            do j = 1,jc
              do i = 1,ic
                coff_tke(i,j,k,3,n) = 0.d0
                coff_tke(i,j,k,4,n) = 0.d0
              end do
            end do
          end do
        end if

        do k = 1,kc
          do j = 1,jc
            do i = 1,ic

              kx = sijk(1,i,j,k,nh)
              ky = sijk(2,i,j,k,nh)
              kz = sijk(3,i,j,k,nh)
              ds = sijk(4,i,j,k,nh)
              if (hawt.and.rotating.and.(.not.tpitch)) then
                ubn  = sijk(8,i,j,k,nh)
              else if (moving.or.(hawt.and.rotating.and.tpitch)) then
                ubn  = sijk(5,i,j,k,nh)*kx + sijk(6,i,j,k,nh)*ky + &
                     sijk(7,i,j,k,nh) * kz
              end if

!------------ calculate fr using qr vector
              rhor   = qr(i,j,k,1,n)
              ur     = qr(i,j,k,2,n)
              vr     = qr(i,j,k,3,n)
              wr     = qr(i,j,k,4,n)
              pr     = qr(i,j,k,5,n)
              thethr = ur*kx + vr*ky + wr*kz - ubn
              q2r    = ur**2 + vr**2 + wr**2
              tker   = coff_tke(i,j,k,3,n)
              h0r    = gamma/(gamma-1)*pr/rhor + q2r/2 + tker
              rh0r   = rhor * h0r
              fr(1)  = rhor*thethr
              fr(2)  = (rhor*ur*thethr+pr*kx)
              fr(3)  = (rhor*vr*thethr+pr*ky)
              fr(4)  = (rhor*wr*thethr+pr*kz)
              fr(5)  = (rh0r * thethr + pr * ubn)

!------------ calculate fl using ql vector
              rhol   = ql(i,j,k,1,n)
              ul     = ql(i,j,k,2,n)
              vl     = ql(i,j,k,3,n)
              wl     = ql(i,j,k,4,n)
              pl     = ql(i,j,k,5,n)
              thethl = ul*kx + vl*ky + wl*kz - ubn
              q2l    = ul**2 + vl**2 + wl**2
              tkel   = coff_tke(i,j,k,4,n)
              h0l    = gamma/(gamma-1)*pl/rhol + q2l/2 + tkel
              rh0l   = rhol * h0l
              fl(1)  = rhol*thethl
              fl(2)  = (rhol*ul*thethl+pl*kx)
              fl(3)  = (rhol*vl*thethl+pl*ky)
              fl(4)  = (rhol*wl*thethl+pl*kz)
              fl(5)  = (rh0l * thethl + pl * ubn)

!------------ calculate delta vector
              drho   = rhor - rhol
              dp     = pr - pl
              du     = ur - ul
              dv     = vr - vl
              dw     = wr - wl
              dtheth = thethr - thethl
              dtke   = tker-tkel

!------------ calculate Roe average values
              rhols  = sqrt(rhol)
              rhors  = sqrt(rhor)
              rho    = rhols*rhors
              denom  = rhols+rhors
              u      = (rhols*ul+rhors*ur)/denom
              v      = (rhols*vl+rhors*vr)/denom
              w      = (rhols*wl+rhors*wr)/denom
              theth  = u*kx + v*ky + w*kz
              h0     = (rhols*h0l+rhors*h0r)/denom
              tke    = (rhols*tkel+rhors*tker)/denom
              q2     = u**2 + v**2 + w**2
              a2     = (gamma-1)*(h0-q2/2-tke)
              a      = sqrt(a2)

              if (moving.and.(lsp_vfr.eq.1)) then
                qp2 = (u-sijk(5,i,j,k,nh))**2 + (v-sijk(6,i,j,k,nh))**2 + &
                    (w-sijk(7,i,j,k,nh))**2
              else
                qp2 = q2
              end if

!------------ preconditioner bits and bops
              m2     = qp2/a2
              ump2   = (uprv / a)**2
              mp2    = dmin1(1.d0, &
                           dmax1(m2,ump2,coff_tke(i,j,k,1,n), &
                                 coff_tke(i,j,k,2,n),epsp(nl)**2))
              mp2s   = dmin1(1.d0, &
                dmax1(m2,coff_tke(i,j,k,1,n),coff_tke(i,j,k,2,n), &
                    (1-mixpi)*ump2,epsp(nl)**2))

!------------ preconditioned eigenvalues and other bits and bops
              lambdap4 = 0.5d0*((theth-ubn)*(1+mp2s)+ &
                            sqrt((theth-ubn)**2*(1-mp2s)**2+4*a2*mp2s))
              lambdap5 = 0.5d0*((theth-ubn)*(1+mp2s)- &
                            sqrt((theth-ubn)**2*(1-mp2s)**2+4*a2*mp2s))

              rm123  = abs(theth-ubn)
              rm4    = abs(lambdap4)
              rm5    = abs(lambdap5)

              if (efx_setup_f) then
!-------------- Entropy Fix based on Yee, Warming, and Harten (1985) and 
!               Yee (). Apply fix only to acoustic eigenvalues to preserve 
!               boundary layer calculations. Avoidance of the carbuncle 
!               phenomena requires fix applied to convective eigenvalue 
!               also. Details based on conversations with Jeff White.
!               Begin {N.B. Comment out this section to remove the entropy fix
                asos = 0.5 * sqrt((theth-ubn)**2*(1-mp2s)**2+4*a2*mp2s)
                rdelta = cntrpy_f * &
                     ( entfx_f_w(1) * entfxctf_f + &
                       entfx_f_w(2) * (sqrt(u**2+v**2+w**2)+asos) + &
                       entfx_f_w(3) * 0.2d0*(abs(theth)+asos+abs(ubn)) )
                rdelta = max(rdelta,1.d-12)
                if ((rm123.le.rdelta).and.efx_all_f) then
                  rm123 = 0.5d0*(rm123**2 + rdelta**2) / rdelta
                end if
                if (rm4.le.rdelta) then
                  rm4  = 0.5d0*(rm4**2  + rdelta**2) / rdelta
                end if
                if (rm5.le.rdelta) then
                  rm5  = 0.5d0*(rm5**2  + rdelta**2) / rdelta
                end if
!               End} Entropy Fix
              end if

              if (eig_limit_f) then
                if (rm123.le.eig_cutoff_f) rm123 = eig_cutoff_f
                if (eig_lim_all_f.and.(rm4.le.eig_cutoff_f)) &
                  rm4   = eig_cutoff_f
                if (eig_lim_all_f.and.(rm5.le.eig_cutoff_f)) &
                  rm5   = eig_cutoff_f
              end if

              esse   = theth - ubn - lambdap4
              erre   = theth - ubn - lambdap5

!------------ preconditioned deltas of characteristic variables
              alphar(i,j,k,1,n)   = (drho -dp/a2)*rm123
              alphar(i,j,k,2,n)   = rho*rm123
              alphar(i,j,k,4,n)   = (rho*erre/(a*(erre-esse))*dtheth + &
                      dp/(a*(erre-esse))) * rm4
              alphar(i,j,k,5,n)   = (rho*esse/(a*(erre-esse))*dtheth + &
                      dp/(a*(erre-esse))) * rm5

!------------ MIXED preconditioned delta fluxes
              df(1)    = alphar(i,j,k,1,n) + &
                (alphar(i,j,k,4,n)*(-esse)+alphar(i,j,k,5,n)*erre)/a/mp2
              df(2)    = alphar(i,j,k,1,n)*u &
                + alphar(i,j,k,2,n)*(du - dtheth*kx) &
                + a*kx*(alphar(i,j,k,4,n)-alphar(i,j,k,5,n)) &
               + u*(alphar(i,j,k,5,n)*erre-alphar(i,j,k,4,n)*esse)/a/mp2
              df(3)    = alphar(i,j,k,1,n)*v &
                + alphar(i,j,k,2,n)*(dv - dtheth*ky) &
                + a*ky*(alphar(i,j,k,4,n)-alphar(i,j,k,5,n)) &
               + v*(alphar(i,j,k,5,n)*erre-alphar(i,j,k,4,n)*esse)/a/mp2
              df(4)    = alphar(i,j,k,1,n)*w &
                + alphar(i,j,k,2,n)*(dw - dtheth*kz) &
                + a*kz*(alphar(i,j,k,4,n)-alphar(i,j,k,5,n)) &
               + w*(alphar(i,j,k,5,n)*erre-alphar(i,j,k,4,n)*esse)/a/mp2
              df(5)    = alphar(i,j,k,1,n)*(q2/2+tke) &
                 + alphar(i,j,k,2,n) * &
                    (u*du + v*dv + w*dw - theth*dtheth + dtke) &
                 + a*theth*(alphar(i,j,k,4,n)-alphar(i,j,k,5,n)) &
               + alphar(i,j,k,4,n)*(-esse*(a2/(gamma-1)+q2/2+tke)/mp2/a) &
               + alphar(i,j,k,5,n)*( erre*(a2/(gamma-1)+q2/2+tke)/mp2/a)

!------------ calculate total flux at interface
              flux(i,j,k,1,n) = (fr(1)+fl(1)-df(1)) / 2 * ds
              flux(i,j,k,2,n) = (fr(2)+fl(2)-df(2)) / 2 * ds
              flux(i,j,k,3,n) = (fr(3)+fl(3)-df(3)) / 2 * ds
              flux(i,j,k,4,n) = (fr(4)+fl(4)-df(4)) / 2 * ds
              flux(i,j,k,5,n) = (fr(5)+fl(5)-df(5)) / 2 * ds

            end do
          end do
        end do

!------ calculate flux of turbulent variables
        if (kom.or.kom_bsl.or.kom_sst) then

          do k = 1,kc
            do j = 1,jc
              do i = 1,ic

                kx = sijk(1,i,j,k,nh)
                ky = sijk(2,i,j,k,nh)
                kz = sijk(3,i,j,k,nh)
                ds = sijk(4,i,j,k,nh)
                if (hawt.and.rotating.and.(.not.tpitch)) then
                  ubn  = sijk(8,i,j,k,nh)
                else if (moving.or.(hawt.and.rotating.and.tpitch)) then
                  ubn  = sijk(5,i,j,k,nh) * kx + sijk(6,i,j,k,nh) * ky + &
                       sijk(7,i,j,k,nh) * kz
                end if

!-------------- calculate fr using qr vector
                rhor   = qr(i,j,k,1,n)
                ur     = qr(i,j,k,2,n)
                vr     = qr(i,j,k,3,n)
                wr     = qr(i,j,k,4,n)
                pr     = qr(i,j,k,5,n)
                thethr = ur*kx + vr*ky + wr*kz - ubn
                q2r    = ur**2 + vr**2 + wr**2
                tker   = coff_tke(i,j,k,3,n)
                omer   = qr(i,j,k,7,n)
                h0r    = gamma/(gamma-1)*pr/rhor + q2r/2 + tker
                fr(1)  = rhor*thethr*tker
                fr(2)  = rhor*thethr*omer

!-------------- calculate fl using ql vector
                rhol   = ql(i,j,k,1,n)
                ul     = ql(i,j,k,2,n)
                vl     = ql(i,j,k,3,n)
                wl     = ql(i,j,k,4,n)
                pl     = ql(i,j,k,5,n)
                thethl = ul*kx + vl*ky + wl*kz - ubn
                q2l    = ul**2 + vl**2 + wl**2
                tkel   = coff_tke(i,j,k,4,n)
                omel   = ql(i,j,k,7,n)
                h0l    = gamma/(gamma-1)*pl/rhol + q2l/2 + tkel
                fl(1)  = rhol*thethl*tkel
                fl(2)  = rhol*thethl*omel

!------------ calculate delta vector
              drho   = rhor - rhol
              dp     = pr - pl
              du     = ur - ul
              dv     = vr - vl
              dw     = wr - wl
              dtheth = thethr - thethl
              dtke   = tker-tkel
              dome   = omer - omel

!------------ calculate Roe average values
              rhols  = sqrt(rhol)
              rhors  = sqrt(rhor)
              rho    = rhols*rhors
              denom  = rhols+rhors
              u      = (rhols*ul+rhors*ur)/denom
              v      = (rhols*vl+rhors*vr)/denom
              w      = (rhols*wl+rhors*wr)/denom
              theth  = u*kx + v*ky + w*kz
              h0     = (rhols*h0l+rhors*h0r)/denom
              tke    = (rhols*tkel+rhors*tker)/denom
              ome    = (rhols*omel+rhors*omer)/denom
              q2     = u**2 + v**2 + w**2
              a2     = (gamma-1)*(h0-q2/2-tke)
              a      = sqrt(a2)

              if (moving.and.(lsp_vfr.eq.1)) then
                qp2 = (u-sijk(5,i,j,k,nh))**2 + (v-sijk(6,i,j,k,nh))**2 + &
                    (w-sijk(7,i,j,k,nh))**2
              else
                qp2 = q2
              end if

!------------ preconditioner bits and bops
              m2     = qp2/a2
              ump2   = (uprv / a)**2
              mp2    = dmin1(1.d0, &
                           dmax1(m2,ump2,coff_tke(i,j,k,1,n), &
                                 coff_tke(i,j,k,2,n),epsp(nl)**2))
              mp2s   = dmin1(1.d0, &
                dmax1(m2,coff_tke(i,j,k,1,n),coff_tke(i,j,k,2,n), &
                    (1-mixpi)*ump2,epsp(nl)**2))

!------------ preconditioned eigenvalues and other bits and bops
              lambdap4 = 0.5d0*((theth-ubn)*(1+mp2s)+ &
                            sqrt((theth-ubn)**2*(1-mp2s)**2+4*a2*mp2s))
              lambdap5 = 0.5d0*((theth-ubn)*(1+mp2s)- &
                            sqrt((theth-ubn)**2*(1-mp2s)**2+4*a2*mp2s))

              rm123  = abs(theth-ubn)
              rm4    = abs(lambdap4)
              rm5    = abs(lambdap5)

              if (efx_setup_t) then
!-------------- Entropy Fix based on Yee, Warming, and Harten (1985) and 
!               Yee (). Apply fix only to acoustic eigenvalues to preserve 
!               boundary layer calculations. Avoidance of the carbuncle 
!               phenomena requires fix applied to convective eigenvalue 
!               also. Details based on conversations with Jeff White.
!               Begin {N.B. Comment out this section to remove the entropy fix
                asos = 0.5 * sqrt((theth-ubn)**2*(1-mp2s)**2+4*a2*mp2s)
                rdelta = cntrpy_t * &
                     ( entfx_t_w(1) * entfxctf_t + &
                       entfx_t_w(2) * (sqrt(u**2+v**2+w**2)+asos) + &
                       entfx_t_w(3) * 0.2d0*(abs(theth)+asos+abs(ubn)) )
                rdelta = max(rdelta,1.d-12)
                if ((rm123.le.rdelta).and.efx_all_t) then
                  rm123 = 0.5d0*(rm123**2 + rdelta**2) / rdelta
                end if
                if (rm4.le.rdelta) then
                  rm4  = 0.5d0*(rm4**2  + rdelta**2) / rdelta
                end if
                if (rm5.le.rdelta) then
                  rm5  = 0.5d0*(rm5**2  + rdelta**2) / rdelta
                end if
!               End} Entropy Fix
              end if

              if (eig_limit_t) then
                if (rm123.le.eig_cutoff_t) rm123 = eig_cutoff_t
                if (eig_lim_all_t.and.(rm4.le.eig_cutoff_t)) &
                  rm4   = eig_cutoff_t
                if (eig_lim_all_t.and.(rm5.le.eig_cutoff_t)) &
                  rm5   = eig_cutoff_t
              end if

              esse   = theth - ubn - lambdap4
              erre   = theth - ubn - lambdap5

!------------ preconditioned deltas of characteristic variables
              alphar(i,j,k,1,n)   = (drho -dp/a2)*rm123
              alphar(i,j,k,2,n)   = rho*rm123
              alphar(i,j,k,4,n)   = (rho*erre/(a*(erre-esse))*dtheth + &
                      dp/(a*(erre-esse))) * rm4
              alphar(i,j,k,5,n)   = (rho*esse/(a*(erre-esse))*dtheth + &
                      dp/(a*(erre-esse))) * rm5

!-------------- MIXED preconditioned delta fluxes
                df(1)    = tke * (alphar(i,j,k,1,n) - &
                  (alphar(i,j,k,4,n)*esse - alphar(i,j,k,5,n)*erre)/a/mp2) &
                  + alphar(i,j,k,2,n) * dtke
                df(2)    = ome * (alphar(i,j,k,1,n) - &
                  (alphar(i,j,k,4,n)*esse - alphar(i,j,k,5,n)*erre)/a/mp2) &
                  + alphar(i,j,k,2,n) * dome

!-------------- calculate total flux at interface
                flux(i,j,k,6,n) = (fr(1)+fl(1)-df(1)) / 2 * ds
                flux(i,j,k,7,n) = (fr(2)+fl(2)-df(2)) / 2 * ds

              end do
            end do
          end do
        end if

      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine vflux(idir,nl,qf,dqxi,dqeta,dqzeta,tau1,flux,fdist, &
                       xider,etader,zetader,sijk)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) idir,nl,iblock,imax,jmax,kmax,ic,jc,kc,iq,iflux, &
                iimt,ivmt,ifdst
      real (kind=cosa_real) qf(*),dqxi(*),dqeta(*),dqzeta(*),flux(*),sijk(*), &
           xider(*),etader(*),zetader(*),tau1(*),fdist(*)

      do iblock = 1,mynblocks
        imax   = i_imax     (iblock,nl)
        jmax   = j_jmax     (iblock,nl)
        kmax   = k_kmax     (iblock,nl)
        if (idir.eq.1) then
          ic = imax
          jc = jmax - 1
          kc = kmax - 1
        else if (idir.eq.2) then
          ic = imax - 1
          jc = jmax
          kc = kmax - 1
        else if (idir.eq.3) then
          ic = imax - 1
          jc = jmax - 1
          kc = kmax
        end if
        iq     = 1 + off_p3 (iblock,nl) * npde * dim5
        iflux  = 1 + off_0  (iblock,nl) * npde * dim5
        iimt   = 1 + off_p2 (iblock,nl) * lmet * dim5h
        ivmt   = 1 + off_0  (iblock,nl) * 3    * dim5h
        ifdst  = 1 + off_0  (iblock,nl)
        call vflux_b(qf(iq),dqxi(iq),dqeta(iq),dqzeta(iq), &
             tau1(iq),flux(iflux), &
             fdist(ifdst),xider(ivmt),etader(ivmt),zetader(ivmt), &
             sijk(iimt),imax,jmax,kmax,npde,nharms,lmet,ic,jc,kc)
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine vflux_b(qf,dqxi,dqeta,dqzeta,tau1,flux,fdist,xider, &
        etader,zetader,sijk,imax,jmax,kmax,npde,nharms,lmet,ic,jc,kc)
!-----------------------------------------------------------------------
!     flux(:,:,:,1,:): on entry, turbulent viscosity mut
!     dqp(:,:,:,1:6,:): Reynolds Stress Tensor (tau1): 
!       txx in tau1(i,j,1,n)
!       txy in tau1(i,j,2,n)
!       txz in tau1(i,j,3,n)
!       tyy in tau1(i,j,4,n)
!       tyz in tau1(i,j,5,n)
!       tzz in tau1(i,j,6,n)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,npde,nharms,lmet
      integer(kind=cosa_int) ic,jc,kc,i,j,k,n,nh
      real (kind=cosa_real) &
        qf     (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
        dqxi   (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
        dqeta  (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
        dqzeta (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
        tau1   (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
        flux   (   imax  ,   jmax  ,   kmax  ,npde,0:2*nharms), &
        sijk   (lmet,0:imax+1,0:jmax+1,0:kmax+1,0:2*nharms*hbmove), &
        xider  (   3,  imax  ,  jmax  ,  kmax  ,0:2*nharms*hbmove), &
        etader (   3,  imax  ,  jmax  ,  kmax  ,0:2*nharms*hbmove), &
        zetader(   3,  imax  ,  jmax  ,  kmax  ,0:2*nharms*hbmove), &
        fdist  (       imax  ,  jmax  ,  kmax)

      real (kind=cosa_real) xix,xiy,xiz,etax,etay,etaz,zetax,zetay,zetaz,kx,ky, &
        kz,ds,rho,u,v,w,p,t,tke,ome,uxi,vxi,wxi,ueta,veta,weta,uzeta, &
        vzeta,wzeta,txi,teta,tzeta,tkexi,tkeeta,tkezeta,omexi,omeeta, &
        omezeta,divv,dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,dtdx, &
        dtdy,dtdz,dtkedx,dtkedy,dtkedz,domedx,domedy,domedz,mu,mut,ctke, &
        come,txx,txy,txz,tyy,tyz,tzz,qx,qy,qz,arg1,cdt,acdt,f1,sigk1, &
        sigkav,sigwav

      if (kom) then

        do n = 0,2*nharms
          nh = n*hbmove
          do k = 1,kc
            do j = 1,jc
              do i = 1,ic

                kx = sijk(1,i,j,k,nh)
                ky = sijk(2,i,j,k,nh)
                kz = sijk(3,i,j,k,nh)
                ds = sijk(4,i,j,k,nh)

                rho  = qf(i,j,k,1,n)
                p    = qf(i,j,k,5,n)
                t    = gamma*p/rho
                mu   = (stemp+1)/(t+stemp) * (t**1.5d0) 
                mut  = flux(i,j,k,1,n)

                xix     = xider  (1,i,j,k,nh)
                xiy     = xider  (2,i,j,k,nh)
                xiz     = xider  (3,i,j,k,nh)
                etax    = etader (1,i,j,k,nh)
                etay    = etader (2,i,j,k,nh)
                etaz    = etader (3,i,j,k,nh)
                zetax   = zetader(1,i,j,k,nh)
                zetay   = zetader(2,i,j,k,nh)
                zetaz   = zetader(3,i,j,k,nh)

                tkexi   = dqxi  (i,j,k,6,n)
                tkeeta  = dqeta (i,j,k,6,n)
                tkezeta = dqzeta(i,j,k,6,n)
                omexi   = dqxi  (i,j,k,7,n)
                omeeta  = dqeta (i,j,k,7,n)
                omezeta = dqzeta(i,j,k,7,n)

                dtkedx  = tkexi*xix + tkeeta*etax + tkezeta*zetax
                dtkedy  = tkexi*xiy + tkeeta*etay + tkezeta*zetay
                dtkedz  = tkexi*xiz + tkeeta*etaz + tkezeta*zetaz
                domedx  = omexi*xix + omeeta*etax + omezeta*zetax
                domedy  = omexi*xiy + omeeta*etay + omezeta*zetay
                domedz  = omexi*xiz + omeeta*etaz + omezeta*zetaz

                txx = tau1(i,j,k,1,n)
                txy = tau1(i,j,k,2,n)
                txz = tau1(i,j,k,3,n)
                tyy = tau1(i,j,k,4,n)
                tyz = tau1(i,j,k,5,n)
                tzz = tau1(i,j,k,6,n)

                flux(i,j,k,2,n) = -(txx*kx + txy *ky + txz *kz) * ds
                flux(i,j,k,3,n) = -(txy*kx + tyy *ky + tyz *kz) * ds
                flux(i,j,k,4,n) = -(txz*kx + tyz *ky + tzz *kz) * ds

                ctke = (mu + sigk * mut) * machfs/reyno
                come = (mu + sigw * mut) * machfs/reyno
                flux(i,j,k,5,n) = -ctke * &
                  (dtkedx*kx+dtkedy*ky+dtkedz*kz) * ds
                flux(i,j,k,6,n) = -ctke * &
                  (dtkedx*kx+dtkedy*ky+dtkedz*kz) * ds
                flux(i,j,k,7,n) = -come * &
                  (domedx*kx+domedy*ky+dtkedz*kz) * ds

             end do
            end do
          end do
        end do

      else if (kom_bsl.or.kom_sst) then

        if (kom_bsl) sigk1 = sigk1_bsl
        if (kom_sst) sigk1 = sigk1_sst

        do n = 0,2*nharms
          nh = n*hbmove
          do k = 1,kc
            do j = 1,jc
              do i = 1,ic

                kx = sijk(1,i,j,k,nh)
                ky = sijk(2,i,j,k,nh)
                kz = sijk(3,i,j,k,nh)
                ds = sijk(4,i,j,k,nh)

                rho  = qf(i,j,k,1,n)
                p    = qf(i,j,k,5,n)
                tke  = qf(i,j,k,6,n)
                ome  = qf(i,j,k,7,n)
                t    = gamma*p/rho
                mu   = (stemp+1)/(t+stemp) * (t**1.5d0)
                mut  = flux(i,j,k,1,n)

                xix     = xider  (1,i,j,k,nh)
                xiy     = xider  (2,i,j,k,nh)
                xiz     = xider  (3,i,j,k,nh)
                etax    = etader (1,i,j,k,nh)
                etay    = etader (2,i,j,k,nh)
                etaz    = etader (3,i,j,k,nh)
                zetax   = zetader(1,i,j,k,nh)
                zetay   = zetader(2,i,j,k,nh)
                zetaz   = zetader(3,i,j,k,nh)

                tkexi   = dqxi  (i,j,k,6,n)
                tkeeta  = dqeta (i,j,k,6,n)
                tkezeta = dqzeta(i,j,k,6,n)
                omexi   = dqxi  (i,j,k,7,n)
                omeeta  = dqeta (i,j,k,7,n)
                omezeta = dqzeta(i,j,k,7,n)

                dtkedx  = tkexi*xix + tkeeta*etax + tkezeta*zetax
                dtkedy  = tkexi*xiy + tkeeta*etay + tkezeta*zetay
                dtkedz  = tkexi*xiz + tkeeta*etaz + tkezeta*zetaz
                domedx  = omexi*xix + omeeta*etax + omezeta*zetax
                domedy  = omexi*xiy + omeeta*etay + omezeta*zetay
                domedz  = omexi*xiz + omeeta*etaz + omezeta*zetaz

                arg1 = dmax1(dsqrt(tke)/0.09d0/ome/(fdist(i,j,k)+epst), &
                             500*mu/rho/ome/(fdist(i,j,k)+epst)**2* &
                             machfs/reyno)
                cdt    = 2*rho*sigw2/ome* &
                         (dtkedx*domedx+dtkedy*domedy+dtkedz*domedz)
                acdt   = dmax1(cdt,1.d-20)
                arg1   = dmin1(arg1,4*rho*sigw2*tke/acdt/ &
                                    (fdist(i,j,k)+epst)**2)
                arg1   = dmin1(arg1**4,1.d2)
                f1     = (dexp(2*arg1)-1)/(dexp(2*arg1)+1)
                sigkav = f1*sigk1 + (1-f1)*sigk2
                sigwav = f1*sigw1 + (1-f1)*sigw2

                txx = tau1(i,j,k,1,n)
                txy = tau1(i,j,k,2,n)
                txz = tau1(i,j,k,3,n)
                tyy = tau1(i,j,k,4,n)
                tyz = tau1(i,j,k,5,n)
                tzz = tau1(i,j,k,6,n)

                flux(i,j,k,2,n) = -(txx*kx + txy *ky + txz *kz) * ds
                flux(i,j,k,3,n) = -(txy*kx + tyy *ky + tyz *kz) * ds
                flux(i,j,k,4,n) = -(txz*kx + tyz *ky + tzz *kz) * ds

                ctke = (mu + sigkav * mut) * machfs/reyno
                come = (mu + sigwav * mut) * machfs/reyno
                flux(i,j,k,5,n) = -ctke * &
                  (dtkedx*kx+dtkedy*ky+dtkedz*kz) * ds
                flux(i,j,k,6,n) = -ctke * &
                  (dtkedx*kx+dtkedy*ky+dtkedz*kz) * ds
                flux(i,j,k,7,n) = -come * &
                  (domedx*kx+domedy*ky+dtkedz*kz) * ds

             end do
           end do
          end do
        end do

      else

        do n = 0,2*nharms
          do k = 1,kc
            do j = 1,jc
              do i = 1,ic
                flux(i,j,k,1,n) = 0
                flux(i,j,k,2,n) = 0
                flux(i,j,k,3,n) = 0
                flux(i,j,k,4,n) = 0
                flux(i,j,k,5,n) = 0
              end do
            end do
          end do
        end do

      end if

      do n = 0,2*nharms
        nh = n*hbmove
        do k = 1,kc
          do j = 1,jc
            do i = 1,ic

              kx = sijk(1,i,j,k,nh)
              ky = sijk(2,i,j,k,nh)
              kz = sijk(3,i,j,k,nh)
              ds = sijk(4,i,j,k,nh)

              rho  = qf(i,j,k,1,n)
              u    = qf(i,j,k,2,n)
              v    = qf(i,j,k,3,n)
              w    = qf(i,j,k,4,n)
              p    = qf(i,j,k,5,n)
              t    = gamma*p/rho
              mu   = (stemp+1)/(t+stemp) * (t**1.5d0)
              mut  = flux(i,j,k,1,n)

              xix     = xider  (1,i,j,k,nh)
              xiy     = xider  (2,i,j,k,nh)
              xiz     = xider  (3,i,j,k,nh)
              etax    = etader (1,i,j,k,nh)
              etay    = etader (2,i,j,k,nh)
              etaz    = etader (3,i,j,k,nh)
              zetax   = zetader(1,i,j,k,nh)
              zetay   = zetader(2,i,j,k,nh)
              zetaz   = zetader(3,i,j,k,nh)

              uxi   = dqxi  (i,j,k,2,n)
              ueta  = dqeta (i,j,k,2,n)
              uzeta = dqzeta(i,j,k,2,n)
              vxi   = dqxi  (i,j,k,3,n)
              veta  = dqeta (i,j,k,3,n)
              vzeta = dqzeta(i,j,k,3,n)
              wxi   = dqxi  (i,j,k,4,n)
              weta  = dqeta (i,j,k,4,n)
              wzeta = dqzeta(i,j,k,4,n)
              txi   = t* ( dqxi  (i,j,k,5,n)/p - dqxi  (i,j,k,1,n)/rho )
              teta  = t* ( dqeta (i,j,k,5,n)/p - dqeta (i,j,k,1,n)/rho )
              tzeta = t* ( dqzeta(i,j,k,5,n)/p - dqzeta(i,j,k,1,n)/rho )

              dudx = uxi*xix + ueta*etax + uzeta*zetax
              dudy = uxi*xiy + ueta*etay + uzeta*zetay
              dudz = uxi*xiz + ueta*etaz + uzeta*zetaz
              dvdx = vxi*xix + veta*etax + vzeta*zetax
              dvdy = vxi*xiy + veta*etay + vzeta*zetay
              dvdz = vxi*xiz + veta*etaz + vzeta*zetaz
              dwdx = wxi*xix + weta*etax + wzeta*zetax
              dwdy = wxi*xiy + weta*etay + wzeta*zetay
              dwdz = wxi*xiz + weta*etaz + wzeta*zetaz
              dtdx = txi*xix + teta*etax + tzeta*zetax
              dtdy = txi*xiy + teta*etay + tzeta*zetay
              dtdz = txi*xiz + teta*etaz + tzeta*zetaz
              divv = dudx + dvdy + dwdz

              txx = machfs/reyno * 2*mu * (dudx - divv/3)
              txy = machfs/reyno * mu*(dudy + dvdx)
              txz = machfs/reyno * mu*(dudz + dwdx)
              tyy = machfs/reyno * 2*mu * (dvdy - divv/3)
              tyz = machfs/reyno * mu*(dvdz + dwdy)
              tzz = machfs/reyno * 2*mu * (dwdz - divv/3)
              qx  = -machfs/reyno * &
                     (mu/pranl+mut/prant) /(gamma-1) * dtdx
              qy  = -machfs/reyno * &
                     (mu/pranl+mut/prant) /(gamma-1) * dtdy
              qz  = -machfs/reyno * &
                     (mu/pranl+mut/prant) /(gamma-1) * dtdz

              flux(i,j,k,1,n) =  0
              flux(i,j,k,2,n) = flux(i,j,k,2,n) &
                 -(txx*kx + txy *ky + txz *kz) * ds
              flux(i,j,k,3,n) = flux(i,j,k,3,n) &
                 -(txy*kx + tyy *ky + tyz *kz) * ds
              flux(i,j,k,4,n) = flux(i,j,k,4,n) &
                 -(txz*kx + tyz *ky + tzz *kz) * ds
              flux(i,j,k,5,n) = flux(i,j,k,5,n) &
                 -( (u*txx + v*txy + w*txz)*kx + &
                    (u*txy + v*tyy + w*tyz)*ky + &
                    (u*txy + v*tyy + w*tyz)*kz ) * ds &
                 +( qx*kx + qy*ky + qz*kz) * ds

           end do
         end do
        end do
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine muscl(idir,nl,bctopo,q,qp,qm,dq,dqp,dqm)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) idir,nl,iblock,imax,jmax,kmax,iblk,iq
      integer(kind=cosa_int) bctopo(*)
      real (kind=cosa_real) q(*),qp(*),qm(*),dq(*),dqp(*),dqm(*)

      do iblock = 1,mynblocks
        imax   = i_imax     (iblock,nl)
        jmax   = j_jmax     (iblock,nl)
        kmax   = k_kmax     (iblock,nl)
        iblk   = 1 + off_bct(iblock,nl)
        iq     = 1 + off_p3 (iblock,nl) * npde * dim5
        if (idir.eq.1) then
          call muscl_bi(nl,bctopo(iblk),q(iq),qp(iq),qm(iq),dq(iq), &
                        dqp(iq),dqm(iq),nbcs(iblock),imax,jmax, &
                        kmax,npde,nharms)
        else if (idir.eq.2) then
          call muscl_bj(nl,bctopo(iblk),q(iq),qp(iq),qm(iq),dq(iq), &
                        dqp(iq),dqm(iq),nbcs(iblock),imax,jmax, &
                        kmax,npde,nharms)
        else if (idir.eq.3) then
          call muscl_bk(nl,bctopo(iblk),q(iq),qp(iq),qm(iq),dq(iq), &
                        dqp(iq),dqm(iq),nbcs(iblock),imax,jmax, &
                        kmax,npde,nharms)
        end if
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine muscl_bi(nl,bctopo,q,qp,qm,dq,dqp,dqm, &
                          nbcs,imax,jmax,kmax,npde,nharms)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,nbcs,npde,nharms
      integer(kind=cosa_int) i,j,k,ipde,ibc,n,nl,ipde_fo_beg, &
                ipde_fo_end,ipde_so_beg,ipde_so_end,bctyp,idir,inrout, &
                ibcpt,iqpm,jstrt,jend,kstrt,kend
      integer(kind=cosa_int) bctopo(10,nbcs)
      real (kind=cosa_real) &
           q  (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           qp (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           qm (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           dq (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           dqp(-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           dqm(-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms)
      real (kind=cosa_real) rkp,rkm
      logical dofirst,doscnd

!     nl_crs  : grid level on and below which first order is used
!     if nl_crs.gt.1 then, on levels nl.lt.nl_crs:
!       rkap = 1/3 : third order upwind biased
!       rkap = 0   : second order upwind biased
!       rkap = -1  : second order fully upwind
!     else
!       first order fully upwind on ALL levels
!     end if

      ipde_fo_beg =  0
      ipde_fo_end = -1
      ipde_so_beg =  1
      ipde_so_end =  npde
      dofirst     = .false.
      doscnd      = .true.

      if (nl.ge.nl_crs) then
        ipde_fo_beg =  1
        ipde_fo_end =  npde
        ipde_so_beg =  0
        ipde_so_end = -1
        dofirst     = .true.
        doscnd      = .false.
      else if (foturb) then
        ipde_fo_beg =  6
        ipde_fo_end =  npde
        ipde_so_beg =  1
        ipde_so_end =  5
        dofirst     = .true.
        doscnd      = .true.
      end if

      if (dofirst) then
!-----------------------------------------------------------------------
!     first order calculation of qp and qm
!-----------------------------------------------------------------------

        do n = 0,2*nharms
          do ipde = ipde_fo_beg,ipde_fo_end
            do k = 1,kmax-1
              do j = 1,jmax-1
                do i = 1,imax
                  qp(i,j,k,ipde,n) = q(i  ,j,k,ipde,n)
                  qm(i,j,k,ipde,n) = q(i-1,j,k,ipde,n)
                end do
              end do
            end do
          end do
        end do

      end if

      if (doscnd) then
!-----------------------------------------------------------------------
!     higher order calculation of qp and qm
!-----------------------------------------------------------------------

        do n=0,2*nharms
          do ipde = ipde_so_beg,ipde_so_end
            do k = 1,kmax-1
              do j = 1,jmax-1
                do i = 0,imax+1
                  dq(i,j,k,ipde,n) = q(i,j,k,ipde,n)-q(i-1,j,k,ipde,n)
                end do
              end do
            end do
          end do
        end do

!------ fix dq on i=1 and i=imax strong boundaries
        do ibc=1,nbcs

          bctyp  = bctopo(1,ibc)
          idir   = bctopo(2,ibc)
          if (idir.eq.1.and. &
              (bctyp.eq.13.or.bctyp.eq.14.or.bctyp/100.eq.15.or. &
               bctyp.eq.22.or.bctyp.eq.23.or.bctyp.eq.92.or. &
               bctyp.eq.93.or.bctyp.eq. 5)) &
                                                                    then
            inrout = bctopo(3,ibc)
            jstrt  = bctopo(6,ibc)
            jend   = bctopo(7,ibc)
            kstrt  = bctopo(8,ibc)
            kend   = bctopo(9,ibc)
            if (inrout.eq. 1) then
               ibcpt = 1
            else
               ibcpt = imax
            end if
            do n = 0,2*nharms
              do ipde = ipde_so_beg,ipde_so_end
                do k = kstrt,kend
                  do j = jstrt,jend
!old              do j = bctopo(6,ibc),bctopo(7,ibc)
                    dq(ibcpt,j,k,ipde,n) = (q(ibcpt  ,j,k,ipde,n)- &
                                            q(ibcpt-1,j,k,ipde,n)) * 2
                  end do
                end do
              end do
            end do

          end if

        end do

!------ limiters and calculation of qm & qp
        call limit_i(dq,dqp,dqm,epslimf,rkap,1,   5,ilimf, &
                     imax,jmax,kmax,npde,nharms)
        call limit_i(dq,dqp,dqm,epslimt,rkap,6,npde,ilimt, &
                     imax,jmax,kmax,npde,nharms)

!------ construction of qm & qp
        do n = 0,2*nharms
          do ipde = ipde_so_beg,ipde_so_end
            do k = 1,kmax-1
              do j = 1,jmax-1
                do i = 1,imax
                  qp(i,j,k,ipde,n) = &
                    q(i  ,j,k,ipde,n) + dqp(i,j,k,ipde,n)
                  qm(i,j,k,ipde,n) = &
                    q(i-1,j,k,ipde,n) + dqm(i,j,k,ipde,n)
                end do
              end do
            end do
          end do
        end do

      end if

!-----------------------------------------------------------------------
!     fix qp and qm on i=1 and i=imax strong boundaries
!-----------------------------------------------------------------------
      do ibc=1,nbcs

        bctyp  = bctopo(1,ibc)
        idir   = bctopo(2,ibc)
        if (idir.eq.1.and. &
            (bctyp.eq.13.or.bctyp.eq.14.or.bctyp/100.eq.15.or. &
             bctyp.eq.22.or.bctyp.eq.23.or.bctyp.eq.92.or. &
             bctyp.eq.93.or.bctyp.eq. 5)) &
                                                                    then
          inrout = bctopo(3,ibc)
          jstrt  = bctopo(6,ibc)
          jend   = bctopo(7,ibc)
          kstrt  = bctopo(8,ibc)
          kend   = bctopo(9,ibc)
          if (inrout.eq.1) then
             iqpm  = 1
             ibcpt = 0
          else
             iqpm  = imax
             ibcpt = imax
          end if

          do n=0,2*nharms
            do ipde = 1,npde
              do k = kstrt,kend
                do j = jstrt,jend
!old            do j = bctopo(6,ibc),bctopo(7,ibc)
                  qp(iqpm,j,k,ipde,n) = q(ibcpt,j,k,ipde,n)
                  qm(iqpm,j,k,ipde,n) = q(ibcpt,j,k,ipde,n)
                end do
              end do
            end do
          end do

        end if

      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine muscl_bj(nl,bctopo,q,qp,qm,dq,dqp,dqm, &
                          nbcs,imax,jmax,kmax,npde,nharms)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,nbcs,npde,nharms
      integer(kind=cosa_int) i,j,k,ipde,ibc,n,nl,ipde_fo_beg, &
                ipde_fo_end,ipde_so_beg,ipde_so_end,bctyp,idir,inrout, &
                jbcpt,jqpm,istrt,iend,kstrt,kend
      integer(kind=cosa_int) bctopo(10,nbcs)
      real (kind=cosa_real) &
           q  (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           qp (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           qm (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           dq (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           dqp(-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           dqm(-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms)
      real (kind=cosa_real) rkp,rkm
      logical dofirst,doscnd

!     nl_crs  : grid level on and below which first order is used
!     if nl_crs.gt.1 then, on levels nl.lt.nl_crs:
!       rkap = 1/3 : third order upwind biased
!       rkap = 0   : second order upwind biased
!       rkap = -1  : second order fully upwind
!     else
!       first order fully upwind on ALL levels
!     end if

      ipde_fo_beg =  0
      ipde_fo_end = -1
      ipde_so_beg =  1
      ipde_so_end =  npde
      dofirst     = .false.
      doscnd      = .true.

      if (nl.ge.nl_crs) then
        ipde_fo_beg =  1
        ipde_fo_end =  npde
        ipde_so_beg =  0
        ipde_so_end = -1
        dofirst     = .true.
        doscnd      = .false.
      else if (foturb) then
        ipde_fo_beg =  6
        ipde_fo_end =  npde
        ipde_so_beg =  1
        ipde_so_end =  5
        dofirst     = .true.
        doscnd      = .true.
      end if

      if (dofirst) then
!-----------------------------------------------------------------------
!     first order calculation of qp and qm
!-----------------------------------------------------------------------

        do n = 0,2*nharms
          do ipde = ipde_fo_beg,ipde_fo_end
            do k = 1,kmax-1
              do j = 1,jmax
                do i = 1,imax-1
                  qp(i,j,k,ipde,n) = q(i,j  ,k,ipde,n)
                  qm(i,j,k,ipde,n) = q(i,j-1,k,ipde,n)
                end do
              end do
            end do
          end do
        end do

      end if

      if (doscnd) then
!-----------------------------------------------------------------------
!     higher order calculation of qp and qm
!-----------------------------------------------------------------------

        do n=0,2*nharms
          do ipde = ipde_so_beg,ipde_so_end
            do k = 1,kmax-1
              do j = 0,jmax+1
                do i = 1,imax-1
                  dq(i,j,k,ipde,n) = q(i,j,k,ipde,n)-q(i,j-1,k,ipde,n)
                end do
              end do
            end do
          end do
        end do

!------ fix dq on j=1 and j=imax strong boundaries
        do ibc=1,nbcs

          bctyp  = bctopo(1,ibc)
          idir   = bctopo(2,ibc)
          if (idir.eq.2.and. &
              (bctyp.eq.13.or.bctyp.eq.14.or.bctyp/100.eq.15.or. &
               bctyp.eq.22.or.bctyp.eq.23.or.bctyp.eq.92.or. &
               bctyp.eq.93.or.bctyp.eq. 5)) &
                                                                    then
            inrout = bctopo(3,ibc)
            istrt  = bctopo(4,ibc)
            iend   = bctopo(5,ibc)
            kstrt  = bctopo(8,ibc)
            kend   = bctopo(9,ibc)
            if (inrout.eq.1) then
               jbcpt = 1
            else
               jbcpt = jmax
            end if
            do n = 0,2*nharms
              do ipde = ipde_so_beg,ipde_so_end
                do k = kstrt,kend
                  do i = istrt,iend
!old              do i = bctopo(4,ibc),bctopo(5,ibc)
                    dq(i,jbcpt,k,ipde,n) = (q(i,jbcpt  ,k,ipde,n)- &
                                            q(i,jbcpt-1,k,ipde,n)) * 2
                  end do
                end do
              end do
            end do

          end if

        end do

!------ limiters
        call limit_j(dq,dqp,dqm,epslimf,rkap,1,   5,ilimf, &
                     imax,jmax,kmax,npde,nharms)
        call limit_j(dq,dqp,dqm,epslimt,rkap,6,npde,ilimt, &
                     imax,jmax,kmax,npde,nharms)

!------ construction of qm & qp
        do n = 0,2*nharms
          do ipde = ipde_so_beg,ipde_so_end
            do k = 1,kmax-1
              do j = 1,jmax
                do i = 1,imax-1
                  qp(i,j,k,ipde,n) = &
                    q(i,j  ,k,ipde,n) + dqp(i,j,k,ipde,n)
                  qm(i,j,k,ipde,n) = &
                    q(i,j-1,k,ipde,n) + dqm(i,j,k,ipde,n)
                end do
              end do
            end do
          end do
        end do

      end if

!-----------------------------------------------------------------------
!     fix qp and qm on j=1 and j=imax strong boundaries
!-----------------------------------------------------------------------
      do ibc=1,nbcs

        bctyp  = bctopo(1,ibc)
        idir   = bctopo(2,ibc)
        if (idir.eq.2.and. &
            (bctyp.eq.13.or.bctyp.eq.14.or.bctyp/100.eq.15.or. &
             bctyp.eq.22.or.bctyp.eq.23.or.bctyp.eq.92.or. &
             bctyp.eq.93.or.bctyp.eq. 5)) &
                                                                    then
          inrout = bctopo(3,ibc)
          istrt  = bctopo(4,ibc)
          iend   = bctopo(5,ibc)
          kstrt  = bctopo(8,ibc)
          kend   = bctopo(9,ibc)
          if (inrout.eq.1) then
             jqpm  = 1
             jbcpt = 0
          else
             jqpm  = jmax
             jbcpt = jmax
          end if

          do n=0,2*nharms
            do ipde = 1,npde
              do k = kstrt,kend
                do i = istrt,iend
!old            do i = bctopo(4,ibc),bctopo(5,ibc)
                  qp(i,jqpm,k,ipde,n) = q(i,jbcpt,k,ipde,n)
                  qm(i,jqpm,k,ipde,n) = q(i,jbcpt,k,ipde,n)
                end do
              end do
            end do
          end do

        end if

      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine muscl_bk(nl,bctopo,q,qp,qm,dq,dqp,dqm, &
                          nbcs,imax,jmax,kmax,npde,nharms)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,nbcs,npde,nharms
      integer(kind=cosa_int) i,j,k,ipde,ibc,n,nl,ipde_fo_beg, &
                ipde_fo_end,ipde_so_beg,ipde_so_end,bctyp,idir,inrout, &
                kbcpt,kqpm,istrt,iend,jstrt,jend
      integer(kind=cosa_int) bctopo(10,nbcs)
      real (kind=cosa_real) &
           q  (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           qp (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           qm (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           dq (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           dqp(-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           dqm(-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms)
      real (kind=cosa_real) rkp,rkm
      logical dofirst,doscnd

!     nl_crs  : grid level on and below which first order is used
!     if nl_crs.gt.1 then, on levels nl.lt.nl_crs:
!       rkap = 1/3 : third order upwind biased
!       rkap = 0   : second order upwind biased
!       rkap = -1  : second order fully upwind
!     else
!       first order fully upwind on ALL levels
!     end if

      ipde_fo_beg =  0
      ipde_fo_end = -1
      ipde_so_beg =  1
      ipde_so_end =  npde
      dofirst     = .false.
      doscnd      = .true.

      if (nl.ge.nl_crs) then
        ipde_fo_beg =  1
        ipde_fo_end =  npde
        ipde_so_beg =  0
        ipde_so_end = -1
        dofirst     = .true.
        doscnd      = .false.
      else if (foturb) then
        ipde_fo_beg =  6
        ipde_fo_end =  npde
        ipde_so_beg =  1
        ipde_so_end =  5
        dofirst     = .true.
        doscnd      = .true.
      end if

      if (dofirst) then
!-----------------------------------------------------------------------
!     first order calculation of qp and qm
!-----------------------------------------------------------------------

        do n = 0,2*nharms
          do ipde = ipde_fo_beg,ipde_fo_end
            do k = 1,kmax
              do j = 1,jmax-1
                do i = 1,imax-1
                  qp(i,j,k,ipde,n) = q(i,j,k  ,ipde,n)
                  qm(i,j,k,ipde,n) = q(i,j,k-1,ipde,n)
                end do
              end do
            end do
          end do
        end do

      end if

      if (doscnd) then
!-----------------------------------------------------------------------
!     higher order calculation of qp and qm
!-----------------------------------------------------------------------

        do n=0,2*nharms
          do ipde = ipde_so_beg,ipde_so_end
            do k = 0,kmax+1
              do j = 1,jmax-1
                do i = 1,imax-1
                  dq(i,j,k,ipde,n) = q(i,j,k,ipde,n)-q(i,j,k-1,ipde,n)
                end do
              end do
            end do
          end do
        end do

!------ fix dq on k=1 and k=kmax strong boundaries
        do ibc=1,nbcs

          bctyp  = bctopo(1,ibc)
          idir   = bctopo(2,ibc)
          if (idir.eq.3.and. &
              (bctyp.eq.13.or.bctyp.eq.14.or.bctyp/100.eq.15.or. &
               bctyp.eq.22.or.bctyp.eq.23.or.bctyp.eq.92.or. &
               bctyp.eq.93.or.bctyp.eq. 5)) &
                                                                    then
            inrout = bctopo(3,ibc)
            istrt  = bctopo(4,ibc)
            iend   = bctopo(5,ibc)
            jstrt  = bctopo(6,ibc)
            jend   = bctopo(7,ibc)
            if (inrout.eq. 1) then
               kbcpt = 1
            else
               kbcpt = kmax
            end if
            do n = 0,2*nharms
              do ipde = ipde_so_beg,ipde_so_end
                do j = jstrt,jend
                  do i = istrt,iend
                    dq(i,j,kbcpt,ipde,n) = (q(i,j,kbcpt  ,ipde,n)- &
                                            q(i,j,kbcpt-1,ipde,n)) * 2
                  end do
                end do
              end do
            end do

          end if

        end do

!------ limiters and calculation of qm & qp
        call limit_k(dq,dqp,dqm,epslimf,rkap,1,   5,ilimf, &
                     imax,jmax,kmax,npde,nharms)
        call limit_k(dq,dqp,dqm,epslimt,rkap,6,npde,ilimt, &
                     imax,jmax,kmax,npde,nharms)

!------ construction of qm & qp
        do n = 0,2*nharms
          do ipde = ipde_so_beg,ipde_so_end
            do k = 1,kmax
              do j = 1,jmax-1
                do i = 1,imax-1
                  qp(i,j,k,ipde,n) = &
                    q(i,j,k  ,ipde,n) + dqp(i,j,k,ipde,n)
                  qm(i,j,k,ipde,n) = &
                    q(i,j,k-1,ipde,n) + dqm(i,j,k,ipde,n)
                end do
              end do
            end do
          end do
        end do

      end if

!-----------------------------------------------------------------------
!     fix qp and qm on i=1 and i=imax strong boundaries
!-----------------------------------------------------------------------
      do ibc=1,nbcs

        bctyp  = bctopo(1,ibc)
        idir   = bctopo(2,ibc)
        if (idir.eq.3.and. &
            (bctyp.eq.13.or.bctyp.eq.14.or.bctyp/100.eq.15.or. &
             bctyp.eq.22.or.bctyp.eq.23.or.bctyp.eq.92.or. &
             bctyp.eq.93.or.bctyp.eq. 5)) &
                                                                    then
          inrout = bctopo(3,ibc)
          istrt  = bctopo(4,ibc)
          iend   = bctopo(5,ibc)
          jstrt  = bctopo(6,ibc)
          jend   = bctopo(7,ibc)
          if (inrout.eq.1) then
             kqpm  = 1
             kbcpt = 0
          else
             kqpm  = kmax
             kbcpt = kmax
          end if

          do n=0,2*nharms
            do ipde = 1,npde
              do j = jstrt,jend
                do i = istrt,iend
                  qp(i,j,kqpm,ipde,n) = q(i,j,kbcpt,ipde,n)
                  qm(i,j,kqpm,ipde,n) = q(i,j,kbcpt,ipde,n)
                end do
              end do
            end do
          end do

        end if

      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine q_der(idir,nl,bctopo,q,dqxi,dqeta,dqzeta)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) idir,nl,iblock,imax,jmax,kmax,iblk,iq
      integer(kind=cosa_int) bctopo(*)
      real (kind=cosa_real) q(*),dqxi(*),dqeta(*),dqzeta(*)

      do iblock = 1,mynblocks
        imax   = i_imax     (iblock,nl)
        jmax   = j_jmax     (iblock,nl)
        kmax   = k_kmax     (iblock,nl)
        iblk   = 1 + off_bct(iblock,nl)
        iq     = 1 + off_p3 (iblock,nl) * npde * dim5
        if (idir.eq.1) then
          call bq_der_i(bctopo(iblk),nbcs(iblock),q(iq),dqxi(iq), &
            dqeta(iq),dqzeta(iq),imax,jmax,kmax,npde,nharms)
        else if (idir.eq.2) then
          call bq_der_j(bctopo(iblk),nbcs(iblock),q(iq),dqxi(iq), &
            dqeta(iq),dqzeta(iq),imax,jmax,kmax,npde,nharms)
        else if (idir.eq.3) then
          call bq_der_k(bctopo(iblk),nbcs(iblock),q(iq),dqxi(iq), &
            dqeta(iq),dqzeta(iq),imax,jmax,kmax,npde,nharms)
        end if
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine bq_der_i(bctopo,nbcs,q,dqxi,dqeta,dqzeta,imax,jmax, &
                          kmax,npde,nharms)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,nbcs,npde,nharms
      integer(kind=cosa_int) bctopo(10,nbcs)
      integer(kind=cosa_int) i,j,k,ipde,ibc,n,bctyp,idir,inrout,jstrt,jend,kstrt, &
                kend,iaux,in1,in2,iface
      real (kind=cosa_real) &
        q     (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
        dqxi  (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
        dqeta (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
        dqzeta(-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms)
      real (kind=cosa_real) q2,q1,sgnm

      do n = 0,2*nharms
        do ipde = 1,npde
          do k = 1,kmax-1
            do j = 1,jmax-1
              do i = 1,imax
                q2 = q(i  ,j,k,ipde,n)
                q1 = q(i-1,j,k,ipde,n)
                dqxi(i,j,k,ipde,n) = q2 - q1
                q2 = ( q(i-1,j,k,ipde,n)   + q(i,j,k,ipde,n) + &
                       q(i,j+1,k,ipde,n)   + q(i-1,j+1,k,ipde,n) ) / 4
                q1 = ( q(i-1,j-1,k,ipde,n) + q(i,j-1,k,ipde,n) + &
                       q(i,j,k,ipde,n)     + q(i-1,j,k,ipde,n)   ) / 4
                dqeta(i,j,k,ipde,n) = q2 - q1
                q2 = ( q(i-1,j,k,ipde,n)   + q(i,j,k,ipde,n) + &
                       q(i,j,k+1,ipde,n)   + q(i-1,j,k+1,ipde,n) ) / 4
                q1 = ( q(i-1,j,k-1,ipde,n) + q(i,j,k-1,ipde,n) + &
                       q(i,j,k,ipde,n)     + q(i-1,j,k,ipde,n)   ) / 4
                dqzeta(i,j,k,ipde,n) = q2 - q1
              end do
            end do
          end do
        end do
      end do

      do ibc = 1,nbcs

!------ correction on i=1 and i=imax boundaries
        bctyp  = bctopo(1,ibc)
        idir   = bctopo(2,ibc)
        if (idir.eq.1.and. &
            (bctyp.eq.13.or.bctyp.eq.14.or.bctyp/100.eq.15.or. &
             bctyp.eq.22.or.bctyp.eq.23.or.bctyp.eq.92.or. &
             bctyp.eq.93.or.bctyp.eq. 5)) &
                                                                    then
          inrout = bctopo(3,ibc)
          jstrt  = bctopo(6,ibc)
          jend   = bctopo(7,ibc)
          kstrt  = bctopo(8,ibc)
          kend   = bctopo(9,ibc)
          if (inrout.eq. 1) then
            iaux  = 0
            in1   = 1
            in2   = 2
            iface = 1
            sgnm  = 1.d0
          else
            iaux  = imax
            in1   = imax - 1
            in2   = imax - 2
            iface = imax
            sgnm  = - 1.d0
          end if
          do n = 0,2*nharms
            do ipde = 1,npde
              do k = kstrt,kend
                do j = jstrt,jend
                  dqxi(iface,j,k,ipde,n) = sgnm / 3 * &
                    (-q(iaux,j,k,ipde,n)*8 + &
                      q(in1,j,k,ipde,n)*9  - q(in2,j,k,ipde,n) )
                  if (j.eq.jstrt) then
                    q2 = ( q(iaux,j  ,k,ipde,n) + &
                           q(iaux,j+1,k,ipde,n) ) / 2
                    q1 = ( q(iaux,j  ,k,ipde,n)*3 - &
                           q(iaux,j+1,k,ipde,n) ) / 2
                  else if (j.eq.jend) then
                    q2 = ( q(iaux,j  ,k,ipde,n)*3 - &
                           q(iaux,j-1,k,ipde,n) ) / 2
                    q1 = ( q(iaux,j  ,k,ipde,n) + &
                           q(iaux,j-1,k,ipde,n) ) / 2
                  else
                    q2 = ( q(iaux,j  ,k,ipde,n) + &
                           q(iaux,j+1,k,ipde,n) ) / 2
                    q1 = ( q(iaux,j-1,k,ipde,n) + &
                           q(iaux,j  ,k,ipde,n) ) / 2
                  end if
                  dqeta(iface,j,k,ipde,n) = q2 - q1
                  if (k.eq.kstrt) then
                    q2 = ( q(iaux,j,k  ,ipde,n) + &
                           q(iaux,j,k+1,ipde,n) ) / 2
                    q1 = ( q(iaux,j,k  ,ipde,n)*3 - &
                           q(iaux,j,k+1,ipde,n) ) / 2
                  else if (k.eq.kend) then
                    q2 = ( q(iaux,j,k  ,ipde,n)*3 - &
                           q(iaux,j,k-1,ipde,n) ) / 2
                    q1 = ( q(iaux,j,k  ,ipde,n) + &
                           q(iaux,j,k-1,ipde,n) ) / 2
                  else
                    q2 = ( q(iaux,j,k  ,ipde,n) + &
                           q(iaux,j,k+1,ipde,n) ) / 2
                    q1 = ( q(iaux,j,k-1,ipde,n) + &
                           q(iaux,j,k  ,ipde,n) ) / 2
                  end if
                  dqzeta(iface,j,k,ipde,n) = q2 - q1
                end do
              end do
            end do
          end do
        else if (idir.eq.1.and.bctyp.eq.83) then
          inrout = bctopo(3,ibc)
          jstrt  = bctopo(6,ibc)
          jend   = bctopo(7,ibc)
          kstrt  = bctopo(8,ibc)
          kend   = bctopo(9,ibc)
          if (inrout.eq. 1) then
            iaux  = 0
            in1   = 1
            in2   = 2
            iface = 1
            sgnm  = 1.d0
          else
            iaux  = imax
            in1   = imax - 1
            in2   = imax - 2
            iface = imax
            sgnm  = - 1.d0
          end if
          do n = 0,2*nharms
            do ipde = 1,npde
              do k = kstrt,kend
                do j = jstrt,jend
                  dqxi(iface,j,k,ipde,n) = 0.d0
                  if (j.eq.jstrt) then
                    q2 = ( q(iaux,j  ,k,ipde,n) + &
                           q(iaux,j+1,k,ipde,n) ) / 2
                    q1 = ( q(iaux,j  ,k,ipde,n)*3 - &
                           q(iaux,j+1,k,ipde,n) ) / 2
                    dqeta(iface,j,k,ipde,n) = q2 - q1
                  else if (j.eq.jend) then
                    q2 = ( q(iaux,j  ,k,ipde,n)*3 - &
                           q(iaux,j-1,k,ipde,n) ) / 2
                    q1 = ( q(iaux,j  ,k,ipde,n) + &
                           q(iaux,j-1,k,ipde,n) ) / 2
                    dqeta(iface,j,k,ipde,n) = q2 - q1
                  end if
                  if (k.eq.kstrt) then
                    q2 = ( q(iaux,j,k  ,ipde,n) + &
                           q(iaux,j,k+1,ipde,n) ) / 2
                    q1 = ( q(iaux,j,k  ,ipde,n)*3 - &
                           q(iaux,j,k+1,ipde,n) ) / 2
                    dqzeta(iface,j,k,ipde,n) = q2 - q1
                  else if (k.eq.kend) then
                    q2 = ( q(iaux,j,k  ,ipde,n)*3 - &
                           q(iaux,j,k-1,ipde,n) ) / 2
                    q1 = ( q(iaux,j,k  ,ipde,n) + &
                           q(iaux,j,k-1,ipde,n) ) / 2
                    dqzeta(iface,j,k,ipde,n) = q2 - q1
                  end if
                end do
              end do
            end do
          end do
        end if

      end do

!dbg  call wrqder(dqxi,dqeta,idir,imax,jmax,npde,nharms)

      return
      end

!-----------------------------------------------------------------------
      subroutine bq_der_j(bctopo,nbcs,q,dqxi,dqeta,dqzeta,imax,jmax, &
                          kmax,npde,nharms)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,nbcs,npde,nharms
      integer(kind=cosa_int) bctopo(10,nbcs)
      integer(kind=cosa_int) i,j,k,ipde,ibc,n,bctyp,idir,inrout,istrt,iend,kstrt, &
                kend,jaux,jn1,jn2,jface
      real (kind=cosa_real) &
        q     (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
        dqxi  (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
        dqeta (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
        dqzeta(-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms)
      real (kind=cosa_real) q2,q1,sgnm

      do n = 0,2*nharms
        do ipde = 1,npde
          do k = 1,kmax-1
            do j = 1,jmax
              do i = 1,imax-1
                q2 = ( q(i,j-1,k,ipde,n)   + q(i,j,k,ipde,n) + &
                       q(i+1,j,k,ipde,n)   + q(i+1,j-1,k,ipde,n) ) / 4
                q1 = ( q(i-1,j-1,k,ipde,n) + q(i-1,j,k,ipde,n) + &
                       q(i,j,k,ipde,n)     + q(i,j-1,k,ipde,n)   ) / 4
                dqxi(i,j,k,ipde,n) = q2 - q1
                q2  = q(i,j  ,k,ipde,n)
                q1  = q(i,j-1,k,ipde,n)
                dqeta(i,j,k,ipde,n) = q2 -q1
                q2 = ( q(i,j-1,k,ipde,n)   + q(i,j,k,ipde,n) + &
                       q(i,j,k+1,ipde,n)   + q(i,j-1,k+1,ipde,n) ) / 4
                q1 = ( q(i,j-1,k-1,ipde,n) + q(i,j,k-1,ipde,n) + &
                       q(i,j,k,ipde,n)     + q(i,j-1,k,ipde,n)   ) / 4
                dqzeta(i,j,k,ipde,n) = q2 - q1
              end do
            end do
          end do
        end do
      end do

      do ibc = 1,nbcs

!------ correction on j=1 and j=jmax boundaries
        bctyp  = bctopo(1,ibc)
        idir   = bctopo(2,ibc)
        if (idir.eq.2.and. &
            (bctyp.eq.13.or.bctyp.eq.14.or.bctyp/100.eq.15.or. &
             bctyp.eq.22.or.bctyp.eq.23.or.bctyp.eq.92.or. &
             bctyp.eq.93.or.bctyp.eq. 5)) &
                                                                    then
          inrout = bctopo(3,ibc)
          istrt  = bctopo(4,ibc)
          iend   = bctopo(5,ibc)
          kstrt  = bctopo(8,ibc)
          kend   = bctopo(9,ibc)
          if (inrout.eq. 1) then
            jaux  = 0
            jn1   = 1
            jn2   = 2
            jface = 1
            sgnm  = 1.d0
          else
            jaux  = jmax
            jn1   = jmax - 1
            jn2   = jmax - 2
            jface = jmax
            sgnm  = - 1.d0
          end if
          do n = 0,2*nharms
            do ipde = 1,npde
              do k = kstrt,kend
                do i = istrt,iend
                  dqeta(i,jface,k,ipde,n) = sgnm / 3 * &
                    (-q(i,jaux,k,ipde,n)*8 + &
                      q(i,jn1,k,ipde,n)*9  - q(i,jn2,k,ipde,n) )
                  if (i.eq.istrt) then
                    q2 = ( q(i+1,jaux,k,ipde,n) + &
                           q(i  ,jaux,k,ipde,n) ) / 2
                    q1 = ( q(i  ,jaux,k,ipde,n)*3 - &
                           q(i+1,jaux,k,ipde,n) ) / 2
                  else if (i.eq.iend) then
                    q2 = ( q(i  ,jaux,k,ipde,n)*3 - &
                           q(i-1,jaux,k,ipde,n) ) / 2
                    q1 = ( q(i  ,jaux,k,ipde,n) + &
                           q(i-1,jaux,k,ipde,n) ) / 2
                  else
                    q2 = ( q(i+1,jaux,k,ipde,n) + &
                           q(i  ,jaux,k,ipde,n) ) / 2
                    q1 = ( q(i  ,jaux,k,ipde,n) + &
                           q(i-1,jaux,k,ipde,n) ) / 2
                  end if
                  dqxi(i,jface,k,ipde,n) = q2 - q1
                  if (k.eq.kstrt) then
                    q2 = ( q(i,jaux,k+1,ipde,n) + &
                           q(i,jaux,k  ,ipde,n) ) / 2
                    q1 = ( q(i,jaux,k  ,ipde,n)*3 - &
                           q(i,jaux,k+1,ipde,n) ) / 2
                  else if (k.eq.kend) then
                    q2 = ( q(i,jaux,k  ,ipde,n)*3 - &
                           q(i,jaux,k-1,ipde,n) ) / 2
                    q1 = ( q(i,jaux,k  ,ipde,n) + &
                           q(i,jaux,k-1,ipde,n) ) / 2
                  else
                    q2 = ( q(i,jaux,k+1,ipde,n) + &
                           q(i,jaux,k  ,ipde,n) ) / 2
                    q1 = ( q(i,jaux,k  ,ipde,n) + &
                           q(i,jaux,k-1,ipde,n) ) / 2
                  end if
                  dqzeta(i,jface,k,ipde,n) = q2 - q1
                end do
              end do
            end do
          end do
        else if (idir.eq.2.and.bctyp.eq.82) then
          inrout = bctopo(3,ibc)
          istrt  = bctopo(4,ibc)
          iend   = bctopo(5,ibc)
          kstrt  = bctopo(8,ibc)
          kend   = bctopo(9,ibc)
          if (inrout.eq. 1) then
            jaux  = 0
            jn1   = 1
            jn2   = 2
            jface = 1
            sgnm  = 1.d0
          else
            jaux  = jmax
            jn1   = jmax - 1
            jn2   = jmax - 2
            jface = jmax
            sgnm  = - 1.d0
          end if
          do n = 0,2*nharms
            do ipde = 1,npde
              do k = kstrt,kend
                do i = istrt,iend
                  dqeta(i,jface,k,ipde,n) = 0.d0
                  if (i.eq.istrt) then
                    q2 = ( q(i+1,jaux,k,ipde,n) + &
                           q(i  ,jaux,k,ipde,n) ) / 2
                    q1 = ( q(i  ,jaux,k,ipde,n)*3 - &
                           q(i+1,jaux,k,ipde,n) ) / 2
                    dqxi(i,jface,k,ipde,n) = q2 - q1
                  else if (i.eq.iend) then
                    q2 = ( q(i  ,jaux,k,ipde,n)*3 - &
                           q(i-1,jaux,k,ipde,n) ) / 2
                    q1 = ( q(i  ,jaux,k,ipde,n) + &
                           q(i-1,jaux,k,ipde,n) ) / 2
                    dqxi(i,jface,k,ipde,n) = q2 - q1
                  end if
                  if (k.eq.kstrt) then
                    q2 = ( q(i,jaux,k+1,ipde,n) + &
                           q(i,jaux,k  ,ipde,n) ) / 2
                    q1 = ( q(i,jaux,k  ,ipde,n)*3 - &
                           q(i,jaux,k+1,ipde,n) ) / 2
                    dqzeta(i,jface,k,ipde,n) = q2 - q1
                  else if (k.eq.kend) then
                    q2 = ( q(i,jaux,k  ,ipde,n)*3 - &
                           q(i,jaux,k-1,ipde,n) ) / 2
                    q1 = ( q(i,jaux,k  ,ipde,n) + &
                           q(i,jaux,k-1,ipde,n) ) / 2
                    dqzeta(i,jface,k,ipde,n) = q2 - q1
                  end if
                end do
              end do
            end do
          end do
        end if

      end do

!dbg  call wrqder(dqxi,dqeta,idir,imax,jmax,npde,nharms)

      return
      end

!-----------------------------------------------------------------------
      subroutine bq_der_k(bctopo,nbcs,q,dqxi,dqeta,dqzeta,imax,jmax, &
                          kmax,npde,nharms)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,nbcs,npde,nharms
      integer(kind=cosa_int) bctopo(10,nbcs)
      integer(kind=cosa_int) i,j,k,ipde,ibc,n,bctyp,idir,inrout,istrt,iend,jstrt, &
                jend,kaux,kn1,kn2,kface
      real (kind=cosa_real) &
        q     (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
        dqxi  (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
        dqeta (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
        dqzeta(-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms)
      real (kind=cosa_real) q2,q1,sgnm

      do n = 0,2*nharms
        do ipde = 1,npde
          do k = 1,kmax
            do j = 1,jmax-1
              do i = 1,imax-1
                q2 = ( q(i,j,k-1,ipde,n)   + q(i,j,k,ipde,n) + &
                       q(i+1,j,k,ipde,n)   + q(i+1,j,k-1,ipde,n) ) / 4
                q1 = ( q(i-1,j,k-1,ipde,n) + q(i-1,j,k,ipde,n) + &
                       q(i,j,k,ipde,n)     + q(i,j,k-1,ipde,n)   ) / 4
                dqxi(i,j,k,ipde,n) = q2 - q1
                q2 = ( q(i,j,k-1,ipde,n)   + q(i,j,k,ipde,n) + &
                       q(i,j+1,k,ipde,n)   + q(i,j+1,k-1,ipde,n) ) / 4
                q1 = ( q(i,j-1,k-1,ipde,n) + q(i,j-1,k,ipde,n) + &
                       q(i,j,k,ipde,n)     + q(i,j,k-1,ipde,n)   ) / 4
                dqeta(i,j,k,ipde,n) = q2 - q1
                q2 = q(i,j,k  ,ipde,n)
                q1 = q(i,j,k-1,ipde,n)
                dqzeta(i,j,k,ipde,n) = q2 - q1
              end do
            end do
          end do
        end do
      end do

      do ibc = 1,nbcs

!------ correction on k=1 and k=kmax boundaries
        bctyp  = bctopo(1,ibc)
        idir   = bctopo(2,ibc)
        if (idir.eq.3.and. &
            (bctyp.eq.13.or.bctyp.eq.14.or.bctyp/100.eq.15.or. &
             bctyp.eq.22.or.bctyp.eq.23.or.bctyp.eq.92.or. &
             bctyp.eq.93.or.bctyp.eq. 5)) &
                                                                    then
          inrout = bctopo(3,ibc)
          istrt  = bctopo(4,ibc)
          iend   = bctopo(5,ibc)
          jstrt  = bctopo(6,ibc)
          jend   = bctopo(7,ibc)
          if (inrout.eq. 1) then
            kaux  = 0
            kn1   = 1
            kn2   = 2
            kface = 1
            sgnm  = 1.d0
          else
            kaux  = kmax
            kn1   = kmax - 1
            kn2   = kmax - 2
            kface = kmax
            sgnm  = - 1.d0
          end if
          do n = 0,2*nharms
            do ipde = 1,npde
              do j = jstrt,jend
                do i = istrt,iend
                  dqzeta(i,j,kface,ipde,n) = sgnm / 3 * &
                    (-q(i,j,kaux,ipde,n)*8 + &
                      q(i,j,kn1,ipde,n)*9  - q(i,j,kn2,ipde,n) )
                  if (i.eq.istrt) then
                    q2 = ( q(i  ,j,kaux,ipde,n) + &
                           q(i+1,j,kaux,ipde,n) ) / 2
                    q1 = ( q(i  ,j,kaux,ipde,n)*3 - &
                           q(i+1,j,kaux,ipde,n) ) / 2
                  else if (i.eq.iend) then
                    q2 = ( q(i  ,j,kaux,ipde,n)*3 - &
                           q(i-1,j,kaux,ipde,n) ) / 2
                    q1 = ( q(i  ,j,kaux,ipde,n) + &
                           q(i-1,j,kaux,ipde,n) ) / 2
                  else
                    q2 = ( q(i  ,j,kaux,ipde,n) + &
                           q(i+1,j,kaux,ipde,n) ) / 2
                    q1 = ( q(i-1,j,kaux,ipde,n) + &
                           q(i  ,j,kaux,ipde,n) ) / 2
                  end if
                  dqxi(i,j,kface,ipde,n) = q2 - q1
                  if (j.eq.jstrt) then
                    q2 = ( q(i,j  ,kaux,ipde,n) + &
                           q(i,j+1,kaux,ipde,n) ) / 2
                    q1 = ( q(i,j  ,kaux,ipde,n)*3 - &
                           q(i,j+1,kaux,ipde,n) ) / 2
                  else if (j.eq.jend) then
                    q2 = ( q(i,j  ,kaux,ipde,n)*3 - &
                           q(i,j-1,kaux,ipde,n) ) / 2
                    q1 = ( q(i,j  ,kaux,ipde,n) + &
                           q(i,j-1,kaux,ipde,n) ) / 2
                  else
                    q2 = ( q(i,j  ,kaux,ipde,n) + &
                           q(i,j+1,kaux,ipde,n) ) / 2
                    q1 = ( q(i,j-1,kaux,ipde,n) + &
                           q(i,j  ,kaux,ipde,n) ) / 2
                  end if
                  dqeta(i,j,kface,ipde,n) = q2 - q1
                end do
              end do
            end do
          end do
        else if (idir.eq.1.and.bctyp.eq.81) then
          inrout = bctopo(3,ibc)
          istrt  = bctopo(4,ibc)
          iend   = bctopo(5,ibc)
          jstrt  = bctopo(6,ibc)
          jend   = bctopo(7,ibc)
          if (inrout.eq. 1) then
            kaux  = 0
            kn1   = 1
            kn2   = 2
            kface = 1
            sgnm  = 1.d0
          else
            kaux  = kmax
            kn1   = kmax - 1
            kn2   = kmax - 2
            kface = kmax
            sgnm  = - 1.d0
          end if
          do n = 0,2*nharms
            do ipde = 1,npde
              do j = jstrt,jend
                do i = istrt,iend
                  dqzeta(i,j,kface,ipde,n) = 0.d0
                  if (i.eq.istrt) then
                    q2 = ( q(i  ,j,kaux,ipde,n) + &
                           q(i+1,j,kaux,ipde,n) ) / 2
                    q1 = ( q(i  ,j,kaux,ipde,n)*3 - &
                           q(i+1,j,kaux,ipde,n) ) / 2
                    dqxi(i,j,kface,ipde,n) = q2 - q1
                  else if (i.eq.iend) then
                    q2 = ( q(i  ,j,kaux,ipde,n)*3 - &
                           q(i-1,j,kaux,ipde,n) ) / 2
                    q1 = ( q(i  ,j,kaux,ipde,n) + &
                           q(i-1,j,kaux,ipde,n) ) / 2
                    dqxi(i,j,kface,ipde,n) = q2 - q1
                  end if
                  if (j.eq.jstrt) then
                    q2 = ( q(i,j  ,kaux,ipde,n) + &
                           q(i,j+1,kaux,ipde,n) ) / 2
                    q1 = ( q(i,j  ,kaux,ipde,n)*3 - &
                           q(i,j+1,kaux,ipde,n) ) / 2
                    dqeta(i,j,kface,ipde,n) = q2 - q1
                  else if (j.eq.jend) then
                    q2 = ( q(i,j  ,kaux,ipde,n)*3 - &
                           q(i,j-1,kaux,ipde,n) ) / 2
                    q1 = ( q(i,j  ,kaux,ipde,n) + &
                           q(i,j-1,kaux,ipde,n) ) / 2
                    dqeta(i,j,kface,ipde,n) = q2 - q1
                  end if
                end do
              end do
            end do
          end do
        end if

      end do

!dbg  call wrqder(dqxi,dqeta,idir,imax,jmax,npde,nharms)

      return
      end

!-----------------------------------------------------------------------
      subroutine q_face(idir,nl,bctopo,q,qf,mut,mutf)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_variables
      use cosa_precision

      implicit none


      integer(kind=cosa_int) idir,nl,iblock,imax,jmax,kmax,iblk,iq,imut,imutf
      integer(kind=cosa_int) bctopo(*)
      real (kind=cosa_real) q(*),qf(*),mut(*),mutf(*)

      do iblock = 1,mynblocks
        imax   = i_imax     (iblock,nl)
        jmax   = j_jmax     (iblock,nl)
        kmax   = k_kmax     (iblock,nl)
        iblk   = 1 + off_bct(iblock,nl)
        iq     = 1 + off_p3 (iblock,nl) * npde * dim5
        imut   = 1 + off_p3 (iblock,nl) *        dim5
        imutf  = 1 + off_0  (iblock,nl) * npde * dim5
        if (idir.eq.1) then
          call bq_face_i(bctopo(iblk),nbcs(iblock),q(iq),qf(iq), &
                         mut(imut),mutf(imutf),imax,jmax,kmax,npde, &
                         nharms)
        else if (idir.eq.2) then
          call bq_face_j(bctopo(iblk),nbcs(iblock),q(iq),qf(iq), &
                         mut(imut),mutf(imutf),imax,jmax,kmax,npde, &
                         nharms)
        else if (idir.eq.3) then
          call bq_face_k(bctopo(iblk),nbcs(iblock),q(iq),qf(iq), &
                         mut(imut),mutf(imutf),imax,jmax,kmax,npde, &
                         nharms)
        end if
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine bq_face_i(bctopo,nbcs,q,qf,mut,mutf,imax,jmax,kmax, &
                           npde,nharms)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,nbcs,npde,nharms
      integer(kind=cosa_int) bctopo(10,nbcs)
      integer(kind=cosa_int) i,j,k,ipde,ibc,n,bctyp,idir,inrout,jstrt,jend,kstrt, &
                kend,ibcpt,iface
      real (kind=cosa_real) &
        q   (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
        qf  (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
        mut (-1:imax+1,-1:jmax+1,-1:kmax+1     ,0:2*nharms), &
        mutf(   imax  ,   jmax  ,   kmax  ,npde,0:2*nharms)

!---- store primitive-form of cell-face flow field in qf
      do n = 0,2*nharms
        do ipde = 1,npde
          do k = 1,kmax-1
            do j = 1,jmax-1
              do i = 1,imax
                qf(i,j,k,ipde,n) = &
                  (q(i-1,j,k,ipde,n) + q(i,j,k,ipde,n)) / 2
              end do
            end do
          end do
        end do
      end do

      do ibc = 1,nbcs

!------ correction on i=1 and i=imax boundaries
        bctyp  = bctopo(1,ibc)
        idir   = bctopo(2,ibc)
        if (idir.eq.1.and. &
            (bctyp.eq.13.or.bctyp.eq.14.or.bctyp/100.eq.15.or. &
             bctyp.eq.22.or.bctyp.eq.23.or.bctyp.eq.92.or. &
             bctyp.eq.93.or.bctyp.eq. 5)) &
                                                                    then
          inrout = bctopo(3,ibc)
          jstrt  = bctopo(6,ibc)
          jend   = bctopo(7,ibc)
          kstrt  = bctopo(8,ibc)
          kend   = bctopo(9,ibc)
          if (inrout.eq. 1) then
             ibcpt = 0
             iface = 1
          else
             ibcpt = imax
             iface = imax
          end if
          do n = 0,2*nharms
            do ipde = 1,npde
              do k = kstrt,kend
                do j = jstrt,jend
!old            do j = bctopo(6,ibc),bctopo(7,ibc)
                  qf(iface,j,k,ipde,n) = q(ibcpt,j,k,ipde,n)
                end do
              end do
            end do
          end do
        end if

      end do

      if (calmut) then

!------ store primitive-form of cell-face flow field in qf
        do n = 0,2*nharms
          do k = 1,kmax-1
            do j = 1,jmax-1
              do i = 1,imax
                mutf(i,j,k,1,n) = (mut(i-1,j,k,n) + mut(i,j,k,n)) / 2
              end do
            end do
          end do
        end do

        do ibc = 1,nbcs

!-------- correction on i=1 and i=imax boundaries
          bctyp  = bctopo(1,ibc)
          idir   = bctopo(2,ibc)
          if (idir.eq.1.and. &
              (bctyp.eq.13.or.bctyp.eq.14.or.bctyp/100.eq.15.or. &
               bctyp.eq.22.or.bctyp.eq.23.or.bctyp.eq.92.or. &
               bctyp.eq.93.or.bctyp.eq. 5)) &
                                                                    then
            inrout = bctopo(3,ibc)
            jstrt  = bctopo(6,ibc)
            jend   = bctopo(7,ibc)
            kstrt  = bctopo(8,ibc)
            kend   = bctopo(9,ibc)
            if (inrout.eq. 1) then
               ibcpt = 0
               iface = 1
            else
               ibcpt = imax
               iface = imax
            end if
            do n = 0,2*nharms
              do k = kstrt,kend
                do j = jstrt,jend
!old            do j = bctopo(6,ibc),bctopo(7,ibc)
                  mutf(iface,j,k,1,n) = mut(ibcpt,j,k,n)
                end do
              end do
            end do
          end if

        end do

      end if

!dbg  call wrqf(qf,idir,imax,jmax,npde,nharms)

      return
      end

!-----------------------------------------------------------------------
      subroutine bq_face_j(bctopo,nbcs,q,qf,mut,mutf,imax,jmax,kmax, &
                           npde,nharms)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,nbcs,npde,nharms
      integer(kind=cosa_int) bctopo(10,nbcs)
      integer(kind=cosa_int) i,j,k,ipde,ibc,n,bctyp,idir,inrout,istrt,iend,kstrt, &
                kend,jbcpt,jface
      real (kind=cosa_real) &
        q   (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
        qf  (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
        mut (-1:imax+1,-1:jmax+1,-1:kmax+1     ,0:2*nharms), &
        mutf(   imax  ,   jmax  ,   kmax  ,npde,0:2*nharms)

!---- store primitive-form of cell-face flow field in qf
      do n=0,2*nharms
        do ipde=1,npde
          do k = 1,kmax-1
            do j = 1,jmax
              do i = 1,imax-1
                qf(i,j,k,ipde,n) = &
                  (q(i,j-1,k,ipde,n) + q(i,j,k,ipde,n)) / 2
              end do
            end do
          end do
        end do
      end do

      do ibc = 1,nbcs

!------ correction on j=1 and j=imax boundaries
        bctyp  = bctopo(1,ibc)
        idir   = bctopo(2,ibc)
        if (idir.eq.2.and. &
            (bctyp.eq.13.or.bctyp.eq.14.or.bctyp/100.eq.15.or. &
             bctyp.eq.22.or.bctyp.eq.23.or.bctyp.eq.92.or. &
             bctyp.eq.93.or.bctyp.eq. 5)) &
                                                                    then
          inrout = bctopo(3,ibc)
          istrt  = bctopo(4,ibc)
          iend   = bctopo(5,ibc)
          kstrt  = bctopo(8,ibc)
          kend   = bctopo(9,ibc)
          if (inrout.eq. 1) then
             jbcpt = 0
             jface = 1
          else
             jbcpt = jmax
             jface = jmax
          end if
          do n = 0,2*nharms
            do ipde = 1,npde
              do k = kstrt,kend
                do i = istrt,iend
!old            do i = bctopo(4,ibc),bctopo(5,ibc)
                  qf(i,jface,k,ipde,n) = q(i,jbcpt,k,ipde,n)
                end do
              end do
            end do
          end do
        end if

      end do

      if (calmut) then

!------ store primitive-form of cell-face flow field in qf
        do n=0,2*nharms
          do k = 1,kmax-1
            do j = 1,jmax
              do i = 1,imax-1
                mutf(i,j,k,1,n) = (mut(i,j-1,k,n) + mut(i,j,k,n)) / 2
              end do
            end do
          end do
        end do

        do ibc = 1,nbcs

!-------- correction on j=1 and j=imax boundaries
          bctyp  = bctopo(1,ibc)
          idir   = bctopo(2,ibc)
          if (idir.eq.2.and. &
              (bctyp.eq.13.or.bctyp.eq.14.or.bctyp/100.eq.15.or. &
               bctyp.eq.22.or.bctyp.eq.23.or.bctyp.eq.92.or. &
               bctyp.eq.93.or.bctyp.eq. 5)) &
                                                                    then
            inrout = bctopo(3,ibc)
            istrt  = bctopo(4,ibc)
            iend   = bctopo(5,ibc)
            kstrt  = bctopo(8,ibc)
            kend   = bctopo(9,ibc)
            if (inrout.eq. 1) then
              jbcpt = 0
              jface = 1
            else
              jbcpt = jmax
              jface = jmax
            end if
            do n = 0,2*nharms
              do k = kstrt,kend
                do i = istrt,iend
!old            do i = bctopo(4,ibc),bctopo(5,ibc)
                  mutf(i,jface,k,1,n) = mut(i,jbcpt,k,n)
                end do
              end do
            end do
          end if

        end do

      end if

!dbg  call wrqf(qf,idir,imax,jmax,npde,nharms)

      return
      end

!-----------------------------------------------------------------------
      subroutine bq_face_k(bctopo,nbcs,q,qf,mut,mutf,imax,jmax,kmax, &
                           npde,nharms)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,nbcs,npde,nharms
      integer(kind=cosa_int) bctopo(10,nbcs)
      integer(kind=cosa_int) i,j,k,ipde,ibc,n,bctyp,idir,inrout,istrt,iend,jstrt, &
                jend,kbcpt,kface
      real (kind=cosa_real) &
        q   (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
        qf  (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
        mut (-1:imax+1,-1:jmax+1,-1:kmax+1     ,0:2*nharms), &
        mutf(   imax  ,   jmax  ,   kmax  ,npde,0:2*nharms)

!---- store primitive-form of cell-face flow field in qf
      do n = 0,2*nharms
        do ipde = 1,npde
          do k = 1,kmax
            do j = 1,jmax-1
              do i = 1,imax-1
                qf(i,j,k,ipde,n) = &
                  (q(i,j,k-1,ipde,n) + q(i,j,k,ipde,n)) / 2
              end do
            end do
          end do
        end do
      end do

      do ibc = 1,nbcs

!------ correction on k=1 and k=kmax boundaries
        bctyp  = bctopo(1,ibc)
        idir   = bctopo(2,ibc)
        if (idir.eq.3.and. &
            (bctyp.eq.13.or.bctyp.eq.14.or.bctyp/100.eq.15.or. &
             bctyp.eq.22.or.bctyp.eq.23.or.bctyp.eq.92.or. &
             bctyp.eq.93.or.bctyp.eq. 5)) &
                                                                    then
          inrout = bctopo(3,ibc)
          istrt  = bctopo(4,ibc)
          iend   = bctopo(5,ibc)
          jstrt  = bctopo(6,ibc)
          jend   = bctopo(7,ibc)
          if (inrout.eq. 1) then
             kbcpt = 0
             kface = 1
          else
             kbcpt = kmax
             kface = kmax
          end if
          do n = 0,2*nharms
            do ipde = 1,npde
              do j = jstrt,jend
                do i = istrt,iend
                  qf(i,j,kface,ipde,n) = q(i,j,kbcpt,ipde,n)
                end do
              end do
            end do
          end do
        end if

      end do

      if (calmut) then

!------ store primitive-form of cell-face flow field in qf
        do n = 0,2*nharms
          do k = 1,kmax
            do j = 1,jmax-1
              do i = 1,imax-1
                mutf(i,j,k,1,n) = (mut(i,j,k-1,n) + mut(i,j,k,n)) / 2
              end do
            end do
          end do
        end do

        do ibc = 1,nbcs

!-------- correction on k=1 and k=kmax boundaries
          bctyp  = bctopo(1,ibc)
          idir   = bctopo(2,ibc)
          if (idir.eq.3.and. &
              (bctyp.eq.13.or.bctyp.eq.14.or.bctyp/100.eq.15.or. &
               bctyp.eq.22.or.bctyp.eq.23.or.bctyp.eq.92.or. &
               bctyp.eq.93.or.bctyp.eq. 5)) &
                                                                    then
            inrout = bctopo(3,ibc)
            istrt  = bctopo(4,ibc)
            iend   = bctopo(5,ibc)
            jstrt  = bctopo(6,ibc)
            jend   = bctopo(7,ibc)
            if (inrout.eq. 1) then
               kbcpt = 0
               kface = 1
            else
               kbcpt = kmax
               kface = kmax
            end if
            do n = 0,2*nharms
              do j = jstrt,jend
                do i = istrt,iend
                  mutf(i,j,kface,1,n) = mut(i,j,kbcpt,n)
                end do
              end do
            end do
          end if

        end do

      end if

!dbg  call wrqf(qf,idir,imax,jmax,npde,nharms)

      return
      end

!-----------------------------------------------------------------------
      subroutine rtst(idir,nl,qf,dqxi,dqeta,dqzeta,xider,etader,zetader, &
                      tau1,mutf)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) idir,nl,imax,jmax,kmax
      integer(kind=cosa_int) iblock,ic,jc,kc,iq,ivmt,imutf
      real (kind=cosa_real) qf(*),dqxi(*),dqeta(*),dqzeta(*),tau1(*),mutf(*), &
           xider(*),etader(*),zetader(*)

      do iblock = 1,mynblocks
        imax   = i_imax     (iblock,nl)
        jmax   = j_jmax     (iblock,nl)
        kmax   = k_kmax     (iblock,nl)
        if (idir.eq.1) then
          ic = imax
          jc = jmax - 1
          kc = kmax - 1
        else if (idir.eq.2) then
          ic = imax - 1
          jc = jmax
          kc = kmax - 1
        else if (idir.eq.3) then
          ic = imax - 1
          jc = jmax - 1
          kc = kmax
        end if
        iq     = 1 + off_p3 (iblock,nl) * npde * dim5
        ivmt   = 1 + off_0  (iblock,nl) * 3    * dim5h
        imutf  = 1 + off_0  (iblock,nl) * npde * dim5
        call b_rtst(qf(iq),dqxi(iq),dqeta(iq),dqzeta(iq), &
             xider(ivmt),etader(ivmt),zetader(ivmt),tau1(iq), &
             mutf(imutf),imax,jmax,kmax,npde,nharms,ic,jc,kc)
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine b_rtst(qf,dqxi,dqeta,dqzeta,xider,etader,zetader,tau1, &
                       mutf,imax,jmax,kmax,npde,nharms,ic,jc,kc)
!-----------------------------------------------------------------------
!     flux(:,:,:,1,:): turbulent viscosity mut
!     dqp(:,:,:,1:6,:): Reynolds Stress Tensor (tau1): 
!       txx in tau1(i,j,1,n)
!       txy in tau1(i,j,2,n)
!       txz in tau1(i,j,3,n)
!       tyy in tau1(i,j,4,n)
!       tyz in tau1(i,j,5,n)
!       tzz in tau1(i,j,6,n)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,npde,nharms
      integer(kind=cosa_int) ic,jc,kc,i,j,k,n,nh,ipde
      real (kind=cosa_real) &
        qf     (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
        dqxi   (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
        dqeta  (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
        dqzeta (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
        tau1   (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
        mutf   (   imax  ,   jmax  ,   kmax  ,npde,0:2*nharms), &
        xider  (3,imax,jmax,kmax,0:2*nharms*hbmove), &
        etader (3,imax,jmax,kmax,0:2*nharms*hbmove), &
        zetader(3,imax,jmax,kmax,0:2*nharms*hbmove)
      real (kind=cosa_real) xix, xiy, xiz, etax, etay, etaz, zetax, zetay, &
        zetaz, kx, ky, kz, rho, p, t, tke, ome, mu, mut, &
        uxi, vxi, wxi, ueta, veta, weta, uzeta, vzeta, wzeta, txi, &
        teta, tzeta, divv, dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, &
        dwdy, dwdz, txx, tyy, tzz, txy, txz, tyz, fact

      fact = 2.d0/3.d0
      do n = 0,2*nharms
        nh = n*hbmove
        do k = 1,kc
          do j = 1,jc
            do i = 1,ic

              rho = qf(i,j,k,1,n)
              tke = qf(i,j,k,6,n)
              ome = qf(i,j,k,7,n)
              mut = mutf(i,j,k,1,n)

              xix   = xider  (1,i,j,k,nh)
              xiy   = xider  (2,i,j,k,nh)
              xiz   = xider  (3,i,j,k,nh)
              etax  = etader (1,i,j,k,nh)
              etay  = etader (2,i,j,k,nh)
              etaz  = etader (3,i,j,k,nh)
              zetax = zetader(1,i,j,k,nh)
              zetay = zetader(2,i,j,k,nh)
              zetaz = zetader(3,i,j,k,nh)

              uxi   = dqxi  (i,j,k,2,n)
              ueta  = dqeta (i,j,k,2,n)
              uzeta = dqzeta(i,j,k,2,n)
              vxi   = dqxi  (i,j,k,3,n)
              veta  = dqeta (i,j,k,3,n)
              vzeta = dqzeta(i,j,k,3,n)
              wxi   = dqxi  (i,j,k,4,n)
              weta  = dqeta (i,j,k,4,n)
              wzeta = dqzeta(i,j,k,4,n)

              dudx  = uxi*xix + ueta*etax + uzeta*zetax
              dudy  = uxi*xiy + ueta*etay + uzeta*zetay
              dudz  = uxi*xiz + ueta*etaz + uzeta*zetaz
              dvdx  = vxi*xix + veta*etax + vzeta*zetax
              dvdy  = vxi*xiy + veta*etay + vzeta*zetay
              dvdz  = vxi*xiz + veta*etaz + vzeta*zetaz
              dwdx  = wxi*xix + weta*etax + wzeta*zetax
              dwdy  = wxi*xiy + weta*etay + wzeta*zetay
              dwdz  = wxi*xiz + weta*etaz + wzeta*zetaz

              divv  = (dudx + dvdy + dwdz) * turcmp

              txx  = 2*mut * (dudx - divv/3)
              txy  = mut*(dudy + dvdx)
              txz  = mut*(dudz + dwdx)
              tyy  = 2*mut * (dvdy - divv/3)
              tyz  = mut*(dvdz + dwdy)
              tzz  = 2*mut * (dwdz - divv/3)

              tau1(i,j,k,1,n) = machfs/reyno*txx - fact*rho*tke
              tau1(i,j,k,2,n) = machfs/reyno*txy
              tau1(i,j,k,3,n) = machfs/reyno*txz
              tau1(i,j,k,4,n) = machfs/reyno*tyy - fact*rho*tke
              tau1(i,j,k,5,n) = machfs/reyno*tyz
              tau1(i,j,k,6,n) = machfs/reyno*tzz - fact*rho*tke

            end do
          end do
        end do
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine q_der_c(idir,nl,qf,dqc)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) idir,nl,iblock,imax,jmax,kmax,iq
      real (kind=cosa_real) qf(*),dqc(*)

      do iblock = 1,mynblocks
        imax   = i_imax     (iblock,nl)
        jmax   = j_jmax     (iblock,nl)
        kmax   = k_kmax     (iblock,nl)
        iq     = 1 + off_p3 (iblock,nl) * npde * dim5
        if (idir.eq.1) then
          call bq_der_ci(qf(iq),dqc(iq),imax,jmax,kmax,npde,nharms)
        else if (idir.eq.2) then
          call bq_der_cj(qf(iq),dqc(iq),imax,jmax,kmax,npde,nharms)
        else if (idir.eq.3) then
          call bq_der_ck(qf(iq),dqc(iq),imax,jmax,kmax,npde,nharms)
        end if
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine bq_der_ci(qf,dqxi,imax,jmax,kmax,npde,nharms)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,npde,nharms
      integer(kind=cosa_int) i,j,k,ipde,n
      real (kind=cosa_real) &
           qf  (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           dqxi(-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms)
      real (kind=cosa_real) q2, q1

      do n = 0,2*nharms
        do ipde = 2,4
          do k = 1,kmax-1
            do j = 1,jmax-1
              do i = 1,imax-1
                q2 = qf(i+1,j,k,ipde,n)
                q1 = qf(i  ,j,k,ipde,n)
                dqxi(i,j,k,ipde,n) = q2 - q1
              end do
            end do
          end do
        end do
      end do

      if (kom_bsl.or.kom_sst) then

        do n = 0,2*nharms
          do ipde = 6,7
            do k = 1,kmax-1
              do j = 1,jmax-1
                do i = 1,imax-1
                  q2 = qf(i+1,j,k,ipde,n)
                  q1 = qf(i,j,k,ipde,n)
                  dqxi(i,j,k,ipde,n) = q2 - q1
                end do
              end do
            end do
          end do
        end do

      end if

      return
      end

!-----------------------------------------------------------------------
      subroutine bq_der_cj(qf,dqeta,imax,jmax,kmax,npde,nharms)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,npde,nharms
      integer(kind=cosa_int) i,j,k,ipde,ibc,n
      real (kind=cosa_real) &
           qf   (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           dqeta(-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms)
      real (kind=cosa_real) q2, q1

      do n = 0,2*nharms
        do ipde = 2,4
          do k = 1,kmax-1
            do j = 1,jmax-1
              do i = 1,imax-1
                q2  = qf(i,j+1,k,ipde,n)
                q1  = qf(i,j  ,k,ipde,n)
                dqeta(i,j,k,ipde,n) = q2 -q1
              end do
            end do
          end do
        end do
      end do

      if (kom_bsl.or.kom_sst) then

        do n = 0,2*nharms
          do ipde = 6,7
            do k = 1,kmax-1
              do j = 1,jmax-1
                do i = 1,imax-1
                  q2  = qf(i,j+1,k,ipde,n)
                  q1  = qf(i,j  ,k,ipde,n)
                  dqeta(i,j,k,ipde,n) = q2 -q1
                end do
              end do
            end do
          end do
        end do

      end if

      return
      end

!-----------------------------------------------------------------------
      subroutine bq_der_ck(qf,dqzeta,imax,jmax,kmax,npde,nharms)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,npde,nharms
      integer(kind=cosa_int) i,j,k,ipde,ibc,n
      real (kind=cosa_real) &
           qf    (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           dqzeta(-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms)
      real (kind=cosa_real) q2, q1

      do n = 0,2*nharms
        do ipde = 2,4
          do k = 1,kmax-1
            do j = 1,jmax-1
              do i = 1,imax-1
                q2  = qf(i,j,k+1,ipde,n)
                q1  = qf(i,j,k  ,ipde,n)
                dqzeta(i,j,k,ipde,n) = q2 -q1
              end do
            end do
          end do
        end do
      end do

      if (kom_bsl.or.kom_sst) then

        do n = 0,2*nharms
          do ipde = 6,7
            do k = 1,kmax-1
              do j = 1,jmax-1
                do i = 1,imax-1
                  q2  = qf(i,j,k+1,ipde,n)
                  q1  = qf(i,j,k  ,ipde,n)
                  dqzeta(i,j,k,ipde,n) = q2 -q1
                end do
              end do
            end do
          end do
        end do

      end if

      return
      end

!-----------------------------------------------------------------------
      subroutine q_corns(nl,q,mut)
!-----------------------------------------------------------------------

      use cosa_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,iblock,imax,jmax,kmax,iq,imut
      real (kind=cosa_real) q(*),mut(*)

      do iblock = 1,mynblocks
        imax   = i_imax     (iblock,nl)
        jmax   = j_jmax     (iblock,nl)
        kmax   = k_kmax     (iblock,nl)
        iq     = 1 + off_p3 (iblock,nl) * npde * dim5
        imut   = 1 + off_p3 (iblock,nl) *        dim5
        call q_bcorns(q(iq),mut(imut),imax,jmax,kmax,npde,nharms)
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine q_bcorns(q,mut,imax,jmax,kmax,npde,nharms)
!-----------------------------------------------------------------------
!     This routine has to be called before q_der_i, q_der_j and q_der_l.
!     The values of q in the 8 corners of the domain need to be set
!     before calling q_der_i, q_der_j and q_der_k. This is to avoid
!     inaccuracy in the calculation of ... to be updated ...
!     dqeta(1,1),dqeta(imax,1),dqeta(1,jmax-1),dqeta(imax,jmax-1)
!     on the xi-faces, and inaccuracies in the calculation of
!     dqxi(1,1),dqxi(1,jmax),dqxi(imax-1,1),dqxi(imax-1,jmax)
!     on the eta-faces.
!-----------------------------------------------------------------------
!     ijkc: array of size 3 x 5 x 8 . First component is i,j,k; second
!           component is cell index of particulr cell of onsidered corner 
!           region; third component denotes corner region.
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,npde,nharms
      integer(kind=cosa_int) ipde,n,ic,ijkc(3,0:4,8)
      real (kind=cosa_real) &
           q  (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           mut(-1:imax+1,-1:jmax+1,-1:kmax+1,     0:2*nharms)
      real (kind=cosa_real) a,b,c,d

!---- region 1
      ijkc(1,0,1) = 0
      ijkc(2,0,1) = 0
      ijkc(3,0,1) = 0
      ijkc(1,1,1) = 0
      ijkc(2,1,1) = 1
      ijkc(3,1,1) = 1
      ijkc(1,2,1) = 1
      ijkc(2,2,1) = 1
      ijkc(3,2,1) = 0
      ijkc(1,3,1) = 1
      ijkc(2,3,1) = 0
      ijkc(3,3,1) = 1
      ijkc(1,4,1) = 1
      ijkc(2,4,1) = 1
      ijkc(3,4,1) = 1

!---- region 2
      ijkc(1,0,2) = imax
      ijkc(2,0,2) = 0
      ijkc(3,0,2) = 0
      ijkc(1,1,2) = imax
      ijkc(2,1,2) = 1
      ijkc(3,1,2) = 1
      ijkc(1,2,2) = imax-1
      ijkc(2,2,2) = 1
      ijkc(3,2,2) = 0
      ijkc(1,3,2) = imax-1
      ijkc(2,3,2) = 0
      ijkc(3,3,2) = 1
      ijkc(1,4,2) = imax-1
      ijkc(2,4,2) = 1
      ijkc(3,4,2) = 1

!---- region 3
      ijkc(1,0,3) = imax
      ijkc(2,0,3) = jmax
      ijkc(3,0,3) = 0
      ijkc(1,1,3) = imax
      ijkc(2,1,3) = jmax-1
      ijkc(3,1,3) = 1
      ijkc(1,2,3) = imax-1
      ijkc(2,2,3) = jmax-1
      ijkc(3,2,3) = 0
      ijkc(1,3,3) = imax-1
      ijkc(2,3,3) = jmax
      ijkc(3,3,3) = 1
      ijkc(1,4,3) = imax-1
      ijkc(2,4,3) = jmax-1
      ijkc(3,4,3) = 1

!---- region 4
      ijkc(1,0,4) = 0
      ijkc(2,0,4) = jmax
      ijkc(3,0,4) = 0
      ijkc(1,1,4) = 0
      ijkc(2,1,4) = jmax-1
      ijkc(3,1,4) = 1
      ijkc(1,2,4) = 1
      ijkc(2,2,4) = jmax-1
      ijkc(3,2,4) = 0
      ijkc(1,3,4) = 1
      ijkc(2,3,4) = jmax
      ijkc(3,3,4) = 1
      ijkc(1,4,4) = 1
      ijkc(2,4,4) = jmax-1
      ijkc(3,4,4) = 1

!---- region 5
      ijkc(1,0,5) = 0
      ijkc(2,0,5) = 0
      ijkc(3,0,5) = kmax
      ijkc(1,1,5) = 0
      ijkc(2,1,5) = 1
      ijkc(3,1,5) = kmax-1
      ijkc(1,2,5) = 1
      ijkc(2,2,5) = 1
      ijkc(3,2,5) = kmax
      ijkc(1,3,5) = 1
      ijkc(2,3,5) = 0
      ijkc(3,3,5) = kmax-1
      ijkc(1,4,5) = 1
      ijkc(2,4,5) = 1
      ijkc(3,4,5) = kmax-1

!---- region 6
      ijkc(1,0,6) = imax
      ijkc(2,0,6) = 0
      ijkc(3,0,6) = kmax
      ijkc(1,1,6) = imax
      ijkc(2,1,6) = 1
      ijkc(3,1,6) = kmax-1
      ijkc(1,2,6) = imax-1
      ijkc(2,2,6) = 1
      ijkc(3,2,6) = kmax
      ijkc(1,3,6) = imax-1
      ijkc(2,3,6) = 0
      ijkc(3,3,6) = kmax-1
      ijkc(1,4,6) = imax-1
      ijkc(2,4,6) = 1
      ijkc(3,4,6) = kmax-1

!---- region 7
      ijkc(1,0,7) = imax
      ijkc(2,0,7) = jmax
      ijkc(3,0,7) = kmax
      ijkc(1,1,7) = imax
      ijkc(2,1,7) = jmax-1
      ijkc(3,1,7) = kmax-1
      ijkc(1,2,7) = imax-1
      ijkc(2,2,7) = jmax-1
      ijkc(3,2,7) = kmax
      ijkc(1,3,7) = imax-1
      ijkc(2,3,7) = jmax
      ijkc(3,3,7) = kmax-1
      ijkc(1,4,7) = imax-1
      ijkc(2,4,7) = jmax-1
      ijkc(3,4,7) = kmax-1

!---- region 8
      ijkc(1,0,8) = 0
      ijkc(2,0,8) = jmax
      ijkc(3,0,8) = kmax
      ijkc(1,1,8) = 0
      ijkc(2,1,8) = jmax-1
      ijkc(3,1,8) = kmax-1
      ijkc(1,2,8) = 1
      ijkc(2,2,8) = jmax-1
      ijkc(3,2,8) = kmax
      ijkc(1,3,8) = 1
      ijkc(2,3,8) = jmax
      ijkc(3,3,8) = kmax-1
      ijkc(1,4,8) = 1
      ijkc(2,4,8) = jmax-1
      ijkc(3,4,8) = kmax-1

      do n = 0,2*nharms
        do ipde = 1,npde
          do ic = 1,8

            a = ( ( (ijkc(3,3,ic) - ijkc(3,1,ic)) * &
                    (ijkc(2,4,ic) - ijkc(2,1,ic)) - &
                    (ijkc(3,4,ic) - ijkc(3,1,ic)) * &
                    (ijkc(2,3,ic) - ijkc(2,1,ic)) ) / &
                  ( ijkc(3,2,ic) - ijkc(3,1,ic) ) * &
                  ( q(ijkc(1,2,ic),ijkc(2,2,ic),ijkc(3,2,ic),ipde,n) - &
                    q(ijkc(1,1,ic),ijkc(2,1,ic),ijkc(3,1,ic),ipde,n) ) + &
                  ( (ijkc(3,2,ic) - ijkc(3,1,ic)) * &
                    (ijkc(2,1,ic) - ijkc(2,4,ic)) + &
                    (ijkc(3,4,ic) - ijkc(3,1,ic)) * &
                    (ijkc(2,2,ic) - ijkc(2,1,ic)) ) / &
                  ( ijkc(3,2,ic) - ijkc(3,1,ic) ) * &
                  ( q(ijkc(1,3,ic),ijkc(2,3,ic),ijkc(3,3,ic),ipde,n) - &
                    q(ijkc(1,1,ic),ijkc(2,1,ic),ijkc(3,1,ic),ipde,n) ) + &
                  ( (ijkc(3,2,ic) - ijkc(3,1,ic)) * &
                    (ijkc(2,3,ic) - ijkc(2,1,ic)) + &
                    (ijkc(3,3,ic) - ijkc(3,1,ic)) * &
                    (ijkc(2,1,ic) - ijkc(2,2,ic)) ) / &
                  ( ijkc(3,2,ic) - ijkc(3,1,ic) ) * &
                  ( q(ijkc(1,4,ic),ijkc(2,4,ic),ijkc(3,4,ic),ipde,n) - &
                    q(ijkc(1,1,ic),ijkc(2,1,ic),ijkc(3,1,ic),ipde,n) ))/ &
                ( ( (ijkc(1,4,ic) - ijkc(1,1,ic)) * &
                    (ijkc(2,3,ic) - ijkc(2,1,ic)) + &
                    (ijkc(1,1,ic) - ijkc(1,3,ic)) * &
                    (ijkc(2,4,ic) - ijkc(2,1,ic)) ) + &
                  ( (ijkc(1,4,ic) - ijkc(1,1,ic)) * &
                    (ijkc(2,1,ic) - ijkc(2,2,ic)) + &
                    (ijkc(1,2,ic) - ijkc(1,1,ic)) * &
                    (ijkc(2,4,ic) - ijkc(2,1,ic)) ) * &
                  ( (ijkc(3,3,ic) - ijkc(3,1,ic)) ) / &
                  ( (ijkc(3,2,ic) - ijkc(3,1,ic)) ) + &
                  ( (ijkc(1,1,ic) - ijkc(1,2,ic)) * &
                    (ijkc(2,3,ic) - ijkc(2,1,ic)) + &
                    (ijkc(1,1,ic) - ijkc(1,3,ic)) * &
                    (ijkc(2,1,ic) - ijkc(2,2,ic)) ) * &
                  ( (ijkc(3,4,ic) - ijkc(3,1,ic)) ) / &
                  ( (ijkc(3,2,ic) - ijkc(3,1,ic)) ) )

            b = ( ( (ijkc(1,1,ic) - ijkc(1,3,ic)) + &
                    (ijkc(1,2,ic) - ijkc(1,1,ic)) * &
                    (ijkc(3,3,ic) - ijkc(3,1,ic)) / &
                    (ijkc(3,2,ic) - ijkc(3,1,ic)) ) * a + &
                  (ijkc(3,3,ic) - ijkc(3,1,ic)) / &
                  (ijkc(3,2,ic) - ijkc(3,1,ic)) * &
                  ( q(ijkc(1,1,ic),ijkc(2,1,ic),ijkc(3,1,ic),ipde,n) - &
                    q(ijkc(1,2,ic),ijkc(2,2,ic),ijkc(3,2,ic),ipde,n) ) + &
                  ( q(ijkc(1,3,ic),ijkc(2,3,ic),ijkc(3,3,ic),ipde,n) - &
                    q(ijkc(1,1,ic),ijkc(2,1,ic),ijkc(3,1,ic),ipde,n) ))/ &
                ( (ijkc(2,3,ic) - ijkc(2,1,ic)) + &
                  (ijkc(2,1,ic) - ijkc(2,2,ic)) * &
                  (ijkc(3,3,ic) - ijkc(3,1,ic)) / &
                  (ijkc(3,2,ic) - ijkc(3,1,ic)) )

            c = ( (ijkc(1,1,ic) - ijkc(1,2,ic)) * a + &
                  (ijkc(2,1,ic) - ijkc(2,2,ic)) * b + &
                  ( q(ijkc(1,2,ic),ijkc(2,2,ic),ijkc(3,2,ic),ipde,n) - &
                    q(ijkc(1,1,ic),ijkc(2,1,ic),ijkc(3,1,ic),ipde,n) ))/ &
                (ijkc(3,2,ic) - ijkc(3,1,ic))

            d =   q(ijkc(1,1,ic),ijkc(2,1,ic),ijkc(3,1,ic),ipde,n) - &
                  a*ijkc(1,1,ic) -b*ijkc(2,1,ic) -c*ijkc(3,1,ic)

            q(ijkc(1,0,ic),ijkc(2,0,ic),ijkc(3,0,ic),ipde,n) = &
              a*ijkc(1,0,ic) + b*ijkc(2,0,ic) + c*ijkc(3,0,ic) + d

          end do
        end do
      end do

      if (calmut) then

        do n = 0,2*nharms
          do ic = 1,8

            a = ( ( (ijkc(3,3,ic) - ijkc(3,1,ic)) * &
                    (ijkc(2,4,ic) - ijkc(2,1,ic)) - &
                    (ijkc(3,4,ic) - ijkc(3,1,ic)) * &
                    (ijkc(2,3,ic) - ijkc(2,1,ic)) ) / &
                  ( ijkc(3,2,ic) - ijkc(3,1,ic) ) * &
                  ( mut(ijkc(1,2,ic),ijkc(2,2,ic),ijkc(3,2,ic),n) - &
                    mut(ijkc(1,1,ic),ijkc(2,1,ic),ijkc(3,1,ic),n) ) + &
                  ( (ijkc(3,2,ic) - ijkc(3,1,ic)) * &
                    (ijkc(2,1,ic) - ijkc(2,4,ic)) + &
                    (ijkc(3,4,ic) - ijkc(3,1,ic)) * &
                    (ijkc(2,2,ic) - ijkc(2,1,ic)) ) / &
                  ( ijkc(3,2,ic) - ijkc(3,1,ic) ) * &
                  ( mut(ijkc(1,3,ic),ijkc(2,3,ic),ijkc(3,3,ic),n) - &
                    mut(ijkc(1,1,ic),ijkc(2,1,ic),ijkc(3,1,ic),n) ) + &
                  ( (ijkc(3,2,ic) - ijkc(3,1,ic)) * &
                    (ijkc(2,3,ic) - ijkc(2,1,ic)) + &
                    (ijkc(3,3,ic) - ijkc(3,1,ic)) * &
                    (ijkc(2,1,ic) - ijkc(2,2,ic)) ) / &
                  ( ijkc(3,2,ic) - ijkc(3,1,ic) ) * &
                  ( mut(ijkc(1,4,ic),ijkc(2,4,ic),ijkc(3,4,ic),n) - &
                    mut(ijkc(1,1,ic),ijkc(2,1,ic),ijkc(3,1,ic),n) ))/ &
                ( ( (ijkc(1,4,ic) - ijkc(1,1,ic)) * &
                    (ijkc(2,3,ic) - ijkc(2,1,ic)) + &
                    (ijkc(1,1,ic) - ijkc(1,3,ic)) * &
                    (ijkc(2,4,ic) - ijkc(2,1,ic)) ) + &
                  ( (ijkc(1,4,ic) - ijkc(1,1,ic)) * &
                    (ijkc(2,1,ic) - ijkc(2,2,ic)) + &
                    (ijkc(1,2,ic) - ijkc(1,1,ic)) * &
                    (ijkc(2,4,ic) - ijkc(2,1,ic)) ) * &
                  ( (ijkc(3,3,ic) - ijkc(3,1,ic)) ) / &
                  ( (ijkc(3,2,ic) - ijkc(3,1,ic)) ) + &
                  ( (ijkc(1,1,ic) - ijkc(1,2,ic)) * &
                    (ijkc(2,3,ic) - ijkc(2,1,ic)) + &
                    (ijkc(1,1,ic) - ijkc(1,3,ic)) * &
                    (ijkc(2,1,ic) - ijkc(2,2,ic)) ) * &
                  ( (ijkc(3,4,ic) - ijkc(3,1,ic)) ) / &
                  ( (ijkc(3,2,ic) - ijkc(3,1,ic)) ) )

            b = ( ( (ijkc(1,1,ic) - ijkc(1,3,ic)) + &
                    (ijkc(1,2,ic) - ijkc(1,1,ic)) * &
                    (ijkc(3,3,ic) - ijkc(3,1,ic)) / &
                    (ijkc(3,2,ic) - ijkc(3,1,ic)) ) * a + &
                  (ijkc(3,3,ic) - ijkc(3,1,ic)) / &
                  (ijkc(3,2,ic) - ijkc(3,1,ic)) * &
                  ( mut(ijkc(1,1,ic),ijkc(2,1,ic),ijkc(3,1,ic),n) - &
                    mut(ijkc(1,2,ic),ijkc(2,2,ic),ijkc(3,2,ic),n) ) + &
                  ( mut(ijkc(1,3,ic),ijkc(2,3,ic),ijkc(3,3,ic),n) - &
                    mut(ijkc(1,1,ic),ijkc(2,1,ic),ijkc(3,1,ic),n) ))/ &
                ( (ijkc(2,3,ic) - ijkc(2,1,ic)) + &
                  (ijkc(2,1,ic) - ijkc(2,2,ic)) * &
                  (ijkc(3,3,ic) - ijkc(3,1,ic)) / &
                  (ijkc(3,2,ic) - ijkc(3,1,ic)) )

            c = ( (ijkc(1,1,ic) - ijkc(1,2,ic)) * a + &
                  (ijkc(2,1,ic) - ijkc(2,2,ic)) * b + &
                  ( mut(ijkc(1,2,ic),ijkc(2,2,ic),ijkc(3,2,ic),n) - &
                    mut(ijkc(1,1,ic),ijkc(2,1,ic),ijkc(3,1,ic),n) ))/ &
                (ijkc(3,2,ic) - ijkc(3,1,ic))

            d =   mut(ijkc(1,1,ic),ijkc(2,1,ic),ijkc(3,1,ic),n) - &
                  a*ijkc(1,1,ic) -b*ijkc(2,1,ic) -c*ijkc(3,1,ic)

            mut(ijkc(1,0,ic),ijkc(2,0,ic),ijkc(3,0,ic),n) = &
              a*ijkc(1,0,ic) + b*ijkc(2,0,ic) + c*ijkc(3,0,ic) + d

          end do
        end do

      end if

      return
      end

!-----------------------------------------------------------------------
      subroutine limit_i(dq,dqp,dqm,epslim,rkap,npde_s,npde_e,ilim, &
                         imax,jmax,kmax,npde,nharms)
!-----------------------------------------------------------------------

      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,npde,nharms
      integer(kind=cosa_int) i,im,j,jm,k,km,ipde,n,npde_s,npde_e,ilim
      real (kind=cosa_real) &
           dq (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           dqp(-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           dqm(-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms)
      real (kind=cosa_real) epslim,rkap,rkp,rkm,s,epsq,dqmx,dqpx,sdm,sdp,arg,b, &
                    d1,d2

!     ilim = 1: smooth limiter. In ISAAC, this corresponds to
!               LIMIT=ILSMTH=1.
!     ilim = 2: min-mod limiter. In ISAAC, this corresponds to
!               LIMIY=ILMNMD=2.
!     ilim = 0: no limiter

      im = imax-1
      jm = jmax-1
      km = kmax-1

      if (ilim.eq.0) then
        rkp = 1+rkap
        rkm = 1-rkap
        do n = 0,2*nharms
          do ipde = npde_s,npde_e
            do k = 1,km
              do j = 1,jm
                do i = 1,imax
                  dqp(i,j,k,ipde,n) = &
                    - (rkp*dq(i  ,j,k,ipde,n)+rkm*dq(i+1,j,k,ipde,n)) /4
                  dqm(i,j,k,ipde,n) = &
                      (rkm*dq(i-1,j,k,ipde,n)+rkp*dq(i  ,j,k,ipde,n)) /4
                end do
              end do
            end do
          end do
        end do
      else if (ilim.eq.1) then
        epsq = epslim*epslim
        do n = 0,2*nharms
          do ipde = npde_s,npde_e
            do k = 1,km
              do j = 1,jm
                do i = 1,imax
                  s = (2*dq(i+1,j,k,ipde,n)*dq(i,j,k,ipde,n)+epslim) / &
                      (dq(i+1,j,k,ipde,n)**2+dq(i,j,k,ipde,n)**2+epslim)
                  rkp = 1+rkap*s
                  rkm = 1-rkap*s
                  dqp(i,j,k,ipde,n) = -s * &
                    (rkp*dq(i  ,j,k,ipde,n)+rkm*dq(i+1,j,k,ipde,n)) / 4
                  s = (2*dq(i,j,k,ipde,n)*dq(i-1,j,k,ipde,n)+epslim) / &
                      (dq(i,j,k,ipde,n)**2+dq(i-1,j,k,ipde,n)**2+epslim)
                  rkp = 1+rkap*s
                  rkm = 1-rkap*s
                  dqm(i,j,k,ipde,n) =  s * &
                    (rkm*dq(i-1,j,k,ipde,n)+rkp*dq(i  ,j,k,ipde,n)) / 4
                end do
              end do
            end do
          end do
        end do
      else if (ilim.eq.2) then
        rkp = 1+rkap
        rkm = 1-rkap
        b = (3.d0-rkap)/(1.d0-rkap)
        do n = 0,2*nharms
          do ipde = npde_s,npde_e
            do k = 1,km
              do j = 1,jm
                do i = 0,imax
                  dqmx = ( sign(1.d0,dq(i  ,j,k,ipde,n)) + &
                           sign(1.d0,dq(i+1,j,k,ipde,n)) ) / 2 * &
                         min( abs(dq(i  ,j,k,ipde,n)), &
                              abs(dq(i+1,j,k,ipde,n))*b)
                  dqpx = ( sign(1.d0,dq(i+1,j,k,ipde,n)) + &
                           sign(1.d0,dq(i  ,j,k,ipde,n)) ) / 2 * &
                         min( abs(dq(i+1,j,k,ipde,n)), &
                              abs(dq(i  ,j,k,ipde,n))*b)
                  dqp(i  ,j,k,ipde,n) = - (rkp*dqmx + rkm*dqpx) / 4
                  dqm(i+1,j,k,ipde,n) =   (rkm*dqmx + rkp*dqpx) / 4
                end do
              end do
            end do
          end do
        end do
      else if (ilim.eq.3) then
        epsq = epslim * epslim
        do n = 0,2*nharms
          do ipde = npde_s,npde_e
            do k = 1,km
              do j = 1,jm
                do i = 0,imax
                  d1   = dq(i+1,j,k,ipde,n)
                  d2   = dq(i  ,j,k,ipde,n)
                  dqmx = ( d1*(d2**2+2*epsq) + d2*(2*d1**2+epsq) ) / &
                         ( 2*d1**2 - d1*d2 + 2*d2**2 + 3*epsq )
                  dqpx = ( d2*(d1**2+2*epsq) + d1*(2*d2**2+epsq) ) / &
                         ( 2*d2**2 - d2*d1 + 2*d1**2 + 3*epsq )
                  dqp(i  ,j,k,ipde,n) = - dqpx / 2
                  dqm(i+1,j,k,ipde,n) =   dqmx / 2
                end do
              end do
            end do
          end do
        end do
      else if (ilim.eq.4) then
        epsq = epslim * epslim
        do n = 0,2*nharms
          do ipde = npde_s,npde_e
            do k = 1,km
              do j = 1,jm
                do i = 0,imax
                  d1   = dq(i+1,j,k,ipde,n)
                  d2   = dq(i  ,j,k,ipde,n)
                  dqmx = ( d1*(d2**2+epsq) + d2*(d1**2+epsq) ) / &
                         ( d1**2 + d2**2 + 2*epsq )
                  dqpx = dqmx
!                 dqpx = ( d2*(d1**2+epsq) + d1*(d2**2+epsq) ) / &
!                        ( d2**2 + d1**2 + 2*epsq )
                  dqp(i  ,j,k,ipde,n) = - dqpx / 2
                  dqm(i+1,j,k,ipde,n) =   dqmx / 2
                end do
              end do
            end do
          end do
        end do
      end if

      return
      end

!-----------------------------------------------------------------------
      subroutine limit_j(dq,dqp,dqm,epslim,rkap,npde_s,npde_e,ilim, &
                         imax,jmax,kmax,npde,nharms)
!-----------------------------------------------------------------------

      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,npde,nharms
      integer(kind=cosa_int) i,im,j,jm,k,km,ipde,n,npde_s,npde_e,ilim
      real (kind=cosa_real) &
           dq (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           dqp(-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           dqm(-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms)
      real (kind=cosa_real) epslim,rkap,rkp,rkm,s,epsq,dqmx,dqpx,sdm,sdp,arg,b, &
                    d1,d2

!     ilim = 1: smooth limiter. In ISAAC, this corresponds to
!               LIMIT=ILSMTH=1.
!     ilim = 2: min-mod limiter. In ISAAC, this corresponds to
!               LIMIY=ILMNMD=2.
!     ilim = 0: no limiter

      im = imax-1
      jm = jmax-1
      km = kmax-1

      if (ilim.eq.0) then
        rkp = 1+rkap
        rkm = 1-rkap
        do n = 0,2*nharms
          do ipde = npde_s,npde_e
            do k = 1,km
              do j = 1,jmax
                do i = 1,im
                  dqp(i,j,k,ipde,n) = &
                    - (rkp*dq(i,j  ,k,ipde,n)+rkm*dq(i,j+1,k,ipde,n)) /4
                  dqm(i,j,k,ipde,n) = &
                      (rkm*dq(i,j-1,k,ipde,n)+rkp*dq(i,j  ,k,ipde,n)) /4
                end do
              end do
            end do
          end do
        end do
      else if (ilim.eq.1) then
        epsq = epslim*epslim
        do n = 0,2*nharms
          do ipde = npde_s,npde_e
            do k = 1,km
              do j = 1,jmax
                do i = 1,im
                  s = (2*dq(i,j+1,k,ipde,n)*dq(i,j,k,ipde,n)+epslim) / &
                      (dq(i,j+1,k,ipde,n)**2+dq(i,j,k,ipde,n)**2+epslim)
                  rkp = 1+rkap*s
                  rkm = 1-rkap*s
                  dqp(i,j,k,ipde,n) = -s * &
                    (rkp*dq(i,j  ,k,ipde,n)+rkm*dq(i,j+1,k,ipde,n)) / 4
                  s = (2*dq(i,j,k,ipde,n)*dq(i,j-1,k,ipde,n)+epslim) / &
                      (dq(i,j,k,ipde,n)**2+dq(i,j-1,k,ipde,n)**2+epslim)
                  rkp = 1+rkap*s
                  rkm = 1-rkap*s
                  dqm(i,j,k,ipde,n) =  s * &
                    (rkm*dq(i,j-1,k,ipde,n)+rkp*dq(i,j  ,k,ipde,n)) / 4
                end do
              end do
            end do
          end do
        end do
      else if (ilim.eq.2) then
        rkp = 1+rkap
        rkm = 1-rkap
        b = (3.d0-rkap)/(1.d0-rkap)
        do n = 0,2*nharms
          do ipde = npde_s,npde_e
            do k = 1,km
              do j = 0,jmax
                do i = 1,im
                  dqmx = ( sign(1.d0,dq(i,j  ,k,ipde,n)) + &
                           sign(1.d0,dq(i,j+1,k,ipde,n)) ) / 2 * &
                         min( abs(dq(i,j  ,k,ipde,n)), &
                              abs(dq(i,j+1,k,ipde,n))*b)
                  dqpx = ( sign(1.d0,dq(i,j+1,k,ipde,n)) + &
                           sign(1.d0,dq(i,j  ,k,ipde,n)) ) / 2 * &
                         min( abs(dq(i,j+1,k,ipde,n)), &
                              abs(dq(i,j  ,k,ipde,n))*b)
                  dqp(i,j  ,k,ipde,n) = - (rkp*dqmx + rkm*dqpx) / 4
                  dqm(i,j+1,k,ipde,n) =   (rkm*dqmx + rkp*dqpx) / 4
                end do
              end do
            end do
          end do
        end do
      else if (ilim.eq.3) then
        epsq = epslim * epslim
        do n = 0,2*nharms
          do ipde = npde_s,npde_e
            do k = 1,km
              do j = 0,jmax
                do i = 1,im
                  d1   = dq(i,j+1,k,ipde,n)
                  d2   = dq(i,j  ,k,ipde,n)
                  dqmx = ( d1*(d2**2+2*epsq) + d2*(2*d1**2+epsq) ) / &
                         ( 2*d1**2 - d1*d2 + 2*d2**2 + 3*epsq )
                  dqpx = ( d2*(d1**2+2*epsq) + d1*(2*d2**2+epsq) ) / &
                         ( 2*d2**2 - d2*d1 + 2*d1**2 + 3*epsq )
                  dqp(i,j  ,k,ipde,n) = - dqpx / 2
                  dqm(i,j+1,k,ipde,n) =   dqmx / 2
                end do
              end do
            end do
          end do
        end do
      else if (ilim.eq.4) then
        epsq = epslim * epslim
        do n = 0,2*nharms
          do ipde = npde_s,npde_e
            do k = 1,km
              do j = 0,jmax
                do i = 1,im
                  d1   = dq(i,j+1,k,ipde,n)
                  d2   = dq(i,j  ,k,ipde,n)
                  dqmx = ( d1*(d2**2+epsq) + d2*(d1**2+epsq) ) / &
                         ( d1**2 + d2**2 + 2*epsq )
                  dqpx = dqmx
!                 dqpx = ( d2*(d1**2+epsq) + d1*(d2**2+epsq) ) / &
!                        ( d2**2 + d1**2 + 2*epsq )
                  dqp(i,j  ,k,ipde,n) = - dqpx / 2
                  dqm(i,j+1,k,ipde,n) =   dqmx / 2
                end do
              end do
            end do
          end do
        end do
      end if

      return
      end

!-----------------------------------------------------------------------
      subroutine limit_k(dq,dqp,dqm,epslim,rkap,npde_s,npde_e,ilim, &
                         imax,jmax,kmax,npde,nharms)
!-----------------------------------------------------------------------

      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,npde,nharms
      integer(kind=cosa_int) i,im,j,jm,k,km,ipde,n,npde_s,npde_e,ilim
      real (kind=cosa_real) &
           dq (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           dqp(-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           dqm(-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms)
      real (kind=cosa_real) epslim,rkap,rkp,rkm,s,epsq,dqmx,dqpx,sdm,sdp,arg,b, &
                    d1,d2

!     ilim = 1: smooth limiter. In ISAAC, this corresponds to
!               LIMIT=ILSMTH=1.
!     ilim = 2: min-mod limiter. In ISAAC, this corresponds to
!               LIMIY=ILMNMD=2.
!     ilim = 0: no limiter

      im = imax-1
      jm = jmax-1
      km = kmax-1

      if (ilim.eq.0) then
        rkp = 1+rkap
        rkm = 1-rkap
        do n = 0,2*nharms
          do ipde = npde_s,npde_e
            do k = 1,kmax
              do j = 1,jm
                do i = 1,im
                  dqp(i,j,k,ipde,n) = &
                    - (rkp*dq(i,j,k  ,ipde,n)+rkm*dq(i,j,k+1,ipde,n)) /4
                  dqm(i,j,k,ipde,n) = &
                      (rkm*dq(i,j,k-1,ipde,n)+rkp*dq(i,j,k  ,ipde,n)) /4
                end do
              end do
            end do
          end do
        end do
      else if (ilim.eq.1) then
        epsq = epslim*epslim
        do n = 0,2*nharms
          do ipde = npde_s,npde_e
            do k = 1,kmax
              do j = 1,jm
                do i = 1,im
                  s = (2*dq(i,j,k+1,ipde,n)*dq(i,j,k,ipde,n)+epslim) / &
                      (dq(i,j,k+1,ipde,n)**2+dq(i,j,k,ipde,n)**2+epslim)
                  rkp = 1+rkap*s
                  rkm = 1-rkap*s
                  dqp(i,j,k,ipde,n) = -s * &
                    (rkp*dq(i,j,k  ,ipde,n)+rkm*dq(i,j,k+1,ipde,n)) / 4
                  s = (2*dq(i,j,k,ipde,n)*dq(i,j,k-1,ipde,n)+epslim) / &
                      (dq(i,j,k,ipde,n)**2+dq(i,j,k-1,ipde,n)**2+epslim)
                  rkp = 1+rkap*s
                  rkm = 1-rkap*s
                  dqm(i,j,k,ipde,n) =  s * &
                    (rkm*dq(i,j,k-1,ipde,n)+rkp*dq(i,j,k  ,ipde,n)) / 4
                end do
              end do
            end do
          end do
        end do
      else if (ilim.eq.2) then
        rkp = 1+rkap
        rkm = 1-rkap
        b = (3.d0-rkap)/(1.d0-rkap)
        do n = 0,2*nharms
          do ipde = npde_s,npde_e
            do k = 0,kmax
              do j = 1,jm
                do i = 1,im
                  dqmx = ( sign(1.d0,dq(i,j,k  ,ipde,n)) + &
                           sign(1.d0,dq(i,j,k+1,ipde,n)) ) / 2 * &
                         min( abs(dq(i,j,k  ,ipde,n)), &
                              abs(dq(i,j,k+1,ipde,n))*b)
                  dqpx = ( sign(1.d0,dq(i,j,k+1,ipde,n)) + &
                           sign(1.d0,dq(i,j,k  ,ipde,n)) ) / 2 * &
                         min( abs(dq(i,j,k+1,ipde,n)), &
                              abs(dq(i,j,k  ,ipde,n))*b)
                  dqp(i,j,k  ,ipde,n) = - (rkp*dqmx + rkm*dqpx) / 4
                  dqm(i,j,k+1,ipde,n) =   (rkm*dqmx + rkp*dqpx) / 4
                end do
              end do
            end do
          end do
        end do
      else if (ilim.eq.3) then
        epsq = epslim * epslim
        do n = 0,2*nharms
          do ipde = npde_s,npde_e
            do k = 0,kmax
              do j = 1,jm
                do i = 1,im
                  d1   = dq(i,j,k+1,ipde,n)
                  d2   = dq(i,j,k  ,ipde,n)
                  dqmx = ( d1*(d2**2+2*epsq) + d2*(2*d1**2+epsq) ) / &
                         ( 2*d1**2 - d1*d2 + 2*d2**2 + 3*epsq )
                  dqpx = ( d2*(d1**2+2*epsq) + d1*(2*d2**2+epsq) ) / &
                         ( 2*d2**2 - d2*d1 + 2*d1**2 + 3*epsq )
                  dqp(i,j,k  ,ipde,n) = - dqpx / 2
                  dqm(i,j,k+1,ipde,n) =   dqmx / 2
                end do
              end do
            end do
          end do
        end do
      else if (ilim.eq.4) then
        epsq = epslim * epslim
        do n = 0,2*nharms
          do ipde = npde_s,npde_e
            do k = 0,kmax
              do j = 1,jm
                do i = 1,im
                  d1   = dq(i,j,k+1,ipde,n)
                  d2   = dq(i,j,k  ,ipde,n)
                  dqmx = ( d1*(d2**2+epsq) + d2*(d1**2+epsq) ) / &
                         ( d1**2 + d2**2 + 2*epsq )
                  dqpx = dqmx
!                 dqpx = ( d2*(d1**2+epsq) + d1*(d2**2+epsq) ) / &
!                        ( d2**2 + d1**2 + 2*epsq )
                  dqp(i,j,k  ,ipde,n) = - dqpx / 2
                  dqm(i,j,k+1,ipde,n) =   dqmx / 2
                end do
              end do
            end do
          end do
        end do
      end if

      return
      end

!-----------------------------------------------------------------------
      subroutine q_edges(nl,q,mut)
!-----------------------------------------------------------------------

      use cosa_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,iblock,imax,jmax,kmax,iq,imut
      real (kind=cosa_real) q(*),mut(*)

      do iblock = 1,mynblocks
        imax   = i_imax     (iblock,nl)
        jmax   = j_jmax     (iblock,nl)
        kmax   = k_kmax     (iblock,nl)
        iq     = 1 + off_p3 (iblock,nl) * npde * dim5
        imut   = 1 + off_p3 (iblock,nl) *        dim5
        call q_bedges(q(iq),mut(imut),imax,jmax,kmax,npde,nharms)
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine q_bedges(q,mut,imax,jmax,kmax,npde,nharms)
!-----------------------------------------------------------------------
!     This routine has to be called before q_der_i, q_der_j and q_der_k.
!     The values of q in the 12 edges of the domain need to be set
!     before calling q_der_i, q_der_j and q_der_j. In the xi-faces case,
!     this is to avoid inaccuracy in the calculation of
!     dqeta(1,1,:),dqeta(imax,1,:),dqeta(1,jmax-1,:),dqeta(imax,jmax-1,:)
!     and dqzeta(1,:,1),dqzeta(imax,:,1),dqzeta(1,:,kmax-1),
!     dqzeta(imax,:,kmax-1)
!-----------------------------------------------------------------------
!     ijkc: array of size 2 x 4 x 12 . First component is (i,j) or (i,k)
!           or (j,k); second
!           component is cell index of particulr cell of onsidered edge
!           region; third component denotes edge region.
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,npde,nharms
      integer(kind=cosa_int) ipde,n,ic,ijkc(2,0:3,12),i,j,k
      real (kind=cosa_real) &
            q  (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
            mut(-1:imax+1,-1:jmax+1,-1:kmax+1,     0:2*nharms)
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
      do n = 0,2*nharms
        do ipde = 1,npde
          do k=1,kmax-1
            do ic = 1,4

              b = ( (q(ijkc(1,3,ic),ijkc(2,3,ic),k,ipde,n) - &
                     q(ijkc(1,1,ic),ijkc(2,1,ic),k,ipde,n)) * &
                    (ijkc(1,2,ic) - ijkc(1,1,ic)) - &
                    (q(ijkc(1,2,ic),ijkc(2,2,ic),k,ipde,n) - &
                     q(ijkc(1,1,ic),ijkc(2,1,ic),k,ipde,n)) * &
                    (ijkc(1,3,ic) - ijkc(1,1,ic)) ) / &
                  ( (ijkc(2,3,ic) - ijkc(2,1,ic)) * &
                    (ijkc(1,2,ic) - ijkc(1,1,ic)) - &
                    (ijkc(2,2,ic) - ijkc(2,1,ic)) * &
                    (ijkc(1,3,ic) - ijkc(1,1,ic)) )

              a = ( q(ijkc(1,2,ic),ijkc(2,2,ic),k,ipde,n) - &
                    q(ijkc(1,1,ic),ijkc(2,1,ic),k,ipde,n) - &
                    b * (ijkc(2,2,ic)-ijkc(2,1,ic)) ) / &
                  ( ijkc(1,2,ic) - ijkc(1,1,ic) )

              c = q(ijkc(1,1,ic),ijkc(2,1,ic),k,ipde,n) - &
                  a*ijkc(1,1,ic) -b*ijkc(2,1,ic)

              q(ijkc(1,0,ic),ijkc(2,0,ic),k,ipde,n) = &
                a*ijkc(1,0,ic) + b*ijkc(2,0,ic) + c

            end do
          end do
        end do
      end do

      if (calmut) then

        do n = 0,2*nharms
          do k=1,kmax-1
            do ic = 1,4

              b = ( (mut(ijkc(1,3,ic),ijkc(2,3,ic),k,n) - &
                     mut(ijkc(1,1,ic),ijkc(2,1,ic),k,n)) * &
                    (ijkc(1,2,ic) - ijkc(1,1,ic)) - &
                    (mut(ijkc(1,2,ic),ijkc(2,2,ic),k,n) - &
                     mut(ijkc(1,1,ic),ijkc(2,1,ic),k,n)) * &
                    (ijkc(1,3,ic) - ijkc(1,1,ic)) ) / &
                  ( (ijkc(2,3,ic) - ijkc(2,1,ic)) * &
                    (ijkc(1,2,ic) - ijkc(1,1,ic)) - &
                    (ijkc(2,2,ic) - ijkc(2,1,ic)) * &
                    (ijkc(1,3,ic) - ijkc(1,1,ic)) )

              a = ( mut(ijkc(1,2,ic),ijkc(2,2,ic),k,n) - &
                    mut(ijkc(1,1,ic),ijkc(2,1,ic),k,n) - &
                    b * (ijkc(2,2,ic)-ijkc(2,1,ic)) ) / &
                  ( ijkc(1,2,ic) - ijkc(1,1,ic) )

              c = mut(ijkc(1,1,ic),ijkc(2,1,ic),k,n) - &
                  a*ijkc(1,1,ic) -b*ijkc(2,1,ic)

              mut(ijkc(1,0,ic),ijkc(2,0,ic),k,n) = &
                a*ijkc(1,0,ic) + b*ijkc(2,0,ic) + c

            end do
          end do
        end do

      end if

!---- regions 5 to 8
      do n = 0,2*nharms
        do ipde = 1,npde
          do j=1,jmax-1
            do ic = 5,8

              b = ( (q(ijkc(1,3,ic),j,ijkc(2,3,ic),ipde,n) - &
                     q(ijkc(1,1,ic),j,ijkc(2,1,ic),ipde,n)) * &
                    (ijkc(1,2,ic) - ijkc(1,1,ic)) - &
                    (q(ijkc(1,2,ic),j,ijkc(2,2,ic),ipde,n) - &
                     q(ijkc(1,1,ic),j,ijkc(2,1,ic),ipde,n)) * &
                    (ijkc(1,3,ic) - ijkc(1,1,ic)) ) / &
                  ( (ijkc(2,3,ic) - ijkc(2,1,ic)) * &
                    (ijkc(1,2,ic) - ijkc(1,1,ic)) - &
                    (ijkc(2,2,ic) - ijkc(2,1,ic)) * &
                    (ijkc(1,3,ic) - ijkc(1,1,ic)) )

              a = ( q(ijkc(1,2,ic),j,ijkc(2,2,ic),ipde,n) - &
                    q(ijkc(1,1,ic),j,ijkc(2,1,ic),ipde,n) - &
                    b * (ijkc(2,2,ic)-ijkc(2,1,ic)) ) / &
                  ( ijkc(1,2,ic) - ijkc(1,1,ic) )

              c = q(ijkc(1,1,ic),j,ijkc(2,1,ic),ipde,n) - &
                  a*ijkc(1,1,ic) -b*ijkc(2,1,ic)

              q(ijkc(1,0,ic),j,ijkc(2,0,ic),ipde,n) = &
                a*ijkc(1,0,ic) + b*ijkc(2,0,ic) + c

            end do
          end do
        end do
      end do

      if (calmut) then

        do n = 0,2*nharms
          do j=1,jmax-1
            do ic = 5,8

              b = ( (mut(ijkc(1,3,ic),j,ijkc(2,3,ic),n) - &
                     mut(ijkc(1,1,ic),j,ijkc(2,1,ic),n)) * &
                    (ijkc(1,2,ic) - ijkc(1,1,ic)) - &
                    (mut(ijkc(1,2,ic),j,ijkc(2,2,ic),n) - &
                     mut(ijkc(1,1,ic),j,ijkc(2,1,ic),n)) * &
                    (ijkc(1,3,ic) - ijkc(1,1,ic)) ) / &
                  ( (ijkc(2,3,ic) - ijkc(2,1,ic)) * &
                    (ijkc(1,2,ic) - ijkc(1,1,ic)) - &
                    (ijkc(2,2,ic) - ijkc(2,1,ic)) * &
                    (ijkc(1,3,ic) - ijkc(1,1,ic)) )

              a = ( mut(ijkc(1,2,ic),j,ijkc(2,2,ic),n) - &
                    mut(ijkc(1,1,ic),j,ijkc(2,1,ic),n) - &
                    b * (ijkc(2,2,ic)-ijkc(2,1,ic)) ) / &
                  ( ijkc(1,2,ic) - ijkc(1,1,ic) )

              c = mut(ijkc(1,1,ic),j,ijkc(2,1,ic),n) - &
                  a*ijkc(1,1,ic) -b*ijkc(2,1,ic)

              mut(ijkc(1,0,ic),j,ijkc(2,0,ic),n) = &
                a*ijkc(1,0,ic) + b*ijkc(2,0,ic) + c

            end do
          end do
        end do

      end if

!---- regions 9 to 12
      do n = 0,2*nharms
        do ipde = 1,npde
          do i=1,imax-1
            do ic = 9,12

              b = ( (q(i,ijkc(1,3,ic),ijkc(2,3,ic),ipde,n) - &
                     q(i,ijkc(1,1,ic),ijkc(2,1,ic),ipde,n)) * &
                    (ijkc(1,2,ic) - ijkc(1,1,ic)) - &
                    (q(i,ijkc(1,2,ic),ijkc(2,2,ic),ipde,n) - &
                     q(i,ijkc(1,1,ic),ijkc(2,1,ic),ipde,n)) * &
                    (ijkc(1,3,ic) - ijkc(1,1,ic)) ) / &
                  ( (ijkc(2,3,ic) - ijkc(2,1,ic)) * &
                    (ijkc(1,2,ic) - ijkc(1,1,ic)) - &
                    (ijkc(2,2,ic) - ijkc(2,1,ic)) * &
                    (ijkc(1,3,ic) - ijkc(1,1,ic)) )

              a = ( q(i,ijkc(1,2,ic),ijkc(2,2,ic),ipde,n) - &
                    q(i,ijkc(1,1,ic),ijkc(2,1,ic),ipde,n) - &
                    b * (ijkc(2,2,ic)-ijkc(2,1,ic)) ) / &
                  ( ijkc(1,2,ic) - ijkc(1,1,ic) )

              c = q(i,ijkc(1,1,ic),ijkc(2,1,ic),ipde,n) - &
                  a*ijkc(1,1,ic) -b*ijkc(2,1,ic)

              q(i,ijkc(1,0,ic),ijkc(2,0,ic),ipde,n) = &
                a*ijkc(1,0,ic) + b*ijkc(2,0,ic) + c

            end do
          end do
        end do
      end do

      if (calmut) then

        do n = 0,2*nharms
          do i=1,imax-1
            do ic = 9,12

              b = ( (mut(i,ijkc(1,3,ic),ijkc(2,3,ic),n) - &
                     mut(i,ijkc(1,1,ic),ijkc(2,1,ic),n)) * &
                    (ijkc(1,2,ic) - ijkc(1,1,ic)) - &
                    (mut(i,ijkc(1,2,ic),ijkc(2,2,ic),n) - &
                     mut(i,ijkc(1,1,ic),ijkc(2,1,ic),n)) * &
                    (ijkc(1,3,ic) - ijkc(1,1,ic)) ) / &
                  ( (ijkc(2,3,ic) - ijkc(2,1,ic)) * &
                    (ijkc(1,2,ic) - ijkc(1,1,ic)) - &
                    (ijkc(2,2,ic) - ijkc(2,1,ic)) * &
                    (ijkc(1,3,ic) - ijkc(1,1,ic)) )

              a = ( mut(i,ijkc(1,2,ic),ijkc(2,2,ic),n) - &
                    mut(i,ijkc(1,1,ic),ijkc(2,1,ic),n) - &
                    b * (ijkc(2,2,ic)-ijkc(2,1,ic)) ) / &
                  ( ijkc(1,2,ic) - ijkc(1,1,ic) )

              c = mut(i,ijkc(1,1,ic),ijkc(2,1,ic),n) - &
                  a*ijkc(1,1,ic) -b*ijkc(2,1,ic)

              mut(i,ijkc(1,0,ic),ijkc(2,0,ic),n) = &
                a*ijkc(1,0,ic) + b*ijkc(2,0,ic) + c

            end do
          end do
        end do

      end if

      return
      end
