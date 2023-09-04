!-----------------------------------------------------------------------
      subroutine turprod(nl,q,dqxi,dqeta,dqzeta,prd,divvel,delpls, &
        xideri,etaderi,zetaderi,xiderj,etaderj,zetaderj,xiderk,etaderk, &
        zetaderk)
!-----------------------------------------------------------------------
!     K Omega model.
!     Calculate main components of source terms of K ^ Omega eqautions.
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,iblock,imax,jmax,kmax,iq,iprd,ivmt,ivol
      real (kind=cosa_real) q(*),dqxi(*),dqeta(*),dqzeta(*),prd(*),divvel(*), &
           delpls(*),xideri(*),etaderi(*),zetaderi(*),xiderj(*), &
           etaderj(*),zetaderj(*),xiderk(*),etaderk(*),zetaderk(*)

      do iblock = 1,mynblocks
        imax   = i_imax     (iblock,nl)
        jmax   = j_jmax     (iblock,nl)
        kmax   = k_kmax     (iblock,nl)
        iq     = 1 + off_p3 (iblock,nl) * npde * dim5
        iprd   = 1 + off_p3 (iblock,nl) *        dim5
        ivmt   = 1 + off_0  (iblock,nl) * 3    * dim5h
        ivol   = 1 + off_p1 (iblock,nl)
        call turprod_b(q(iq),dqxi(iq),dqeta(iq),dqzeta(iq), &
          prd(iprd),divvel(iprd),delpls(iprd),xideri(ivmt), &
          etaderi(ivmt),zetaderi(ivmt),xiderj(ivmt),etaderj(ivmt), &
          zetaderj(ivmt),xiderk(ivmt),etaderk(ivmt),zetaderk(ivmt), &
          imax,jmax,kmax,npde,nharms)
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine turprod_b(q,dqxi,dqeta,dqzeta,prd,divvel,delpls,xideri, &
        etaderi,zetaderi,xiderj,etaderj,zetaderj,xiderk,etaderk, &
        zetaderk,imax,jmax,kmax,npde,nharms)
!-----------------------------------------------------------------------
!     K Omega model.
!     Calculate main components of source terms of K ^ Omega eqautions.
!     prod: term P_d as by eqs. (31) & (32) of Liu and Zheng, JCP, 1996
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,npde,nharms
      integer(kind=cosa_int) i,j,k,n,nh
      real (kind=cosa_real) &
        q      (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
        dqxi   (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
        dqeta  (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
        dqzeta (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
        prd    (-1:imax+1,-1:jmax+1,-1:kmax+1,     0:2*nharms), &
        divvel (-1:imax+1,-1:jmax+1,-1:kmax+1,     0:2*nharms), &
        delpls (-1:imax+1,-1:jmax+1,-1:kmax+1,     0:2*nharms), &
        xideri  (3   ,imax    ,jmax    ,kmax    ,0:2*nharms*hbmove), &
        xiderj  (3   ,imax    ,jmax    ,kmax    ,0:2*nharms*hbmove), &
        xiderk  (3   ,imax    ,jmax    ,kmax    ,0:2*nharms*hbmove), &
        etaderi (3   ,imax    ,jmax    ,kmax    ,0:2*nharms*hbmove), &
        etaderj (3   ,imax    ,jmax    ,kmax    ,0:2*nharms*hbmove), &
        etaderk (3   ,imax    ,jmax    ,kmax    ,0:2*nharms*hbmove), &
        zetaderi(3   ,imax    ,jmax    ,kmax    ,0:2*nharms*hbmove), &
        zetaderj(3   ,imax    ,jmax    ,kmax    ,0:2*nharms*hbmove), &
        zetaderk(3   ,imax    ,jmax    ,kmax    ,0:2*nharms*hbmove)

      real (kind=cosa_real) xix,xiy,xiz,etax,etay,etaz,zetax,zetay,zetaz,rho, &
           tke,ome,uxi,vxi,wxi,ueta,veta,weta,uzeta,vzeta,wzeta, &
           dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,divv,txx,txy, &
           txz,tyy,tyz,tzz,prod,fact

      fact = 2.d0/3.d0

      do n = 0,2*nharms
        nh = n*hbmove
        do k = 1,kmax-1
          do j = 1,jmax-1
            do i = 1,imax-1

              rho = q(i,j,k,1,n)
              tke = q(i,j,k,6,n)
              ome = q(i,j,k,7,n)

              xix   = (xideri  (1,i+1,j  ,k  ,nh) + &
                       xideri  (1,i  ,j  ,k  ,nh)) / 2
              xiy   = (xiderj  (2,i  ,j+1,k  ,nh) + &
                       xiderj  (2,i  ,j  ,k  ,nh)) / 2
              xiz   = (xiderk  (3,i  ,j  ,k+1,nh) + &
                       xiderk  (3,i  ,j  ,k  ,nh)) / 2
              etax  = (etaderi (1,i+1,j  ,k  ,nh) + &
                       etaderi (1,i  ,j  ,k  ,nh)) / 2
              etay  = (etaderj (2,i  ,j+1,k  ,nh) + &
                       etaderj (2,i  ,j  ,k  ,nh)) / 2
              etaz  = (etaderk (3,i  ,j  ,k+1,nh) + &
                       etaderk (3,i  ,j  ,k  ,nh)) / 2
              zetax = (zetaderi(1,i+1,j  ,k  ,nh) + &
                       zetaderi(1,i  ,j  ,k  ,nh)) / 2
              zetay = (zetaderj(2,i  ,j+1,k  ,nh) + &
                       zetaderj(2,i  ,j  ,k  ,nh)) / 2
              zetaz = (zetaderk(3,i  ,j  ,k+1,nh) + &
                       zetaderk(3,i  ,j  ,k  ,nh)) / 2

              uxi   = dqxi  (i,j,k,2,n)
              ueta  = dqeta (i,j,k,2,n)
              uzeta = dqzeta(i,j,k,2,n)
              vxi   = dqxi  (i,j,k,3,n)
              veta  = dqeta (i,j,k,3,n)
              vzeta = dqzeta(i,j,k,3,n)
              wxi   = dqxi  (i,j,k,4,n)
              weta  = dqeta (i,j,k,4,n)
              wzeta = dqzeta(i,j,k,4,n)

              dudx = uxi*xix + ueta*etax + uzeta*zetax
              dudy = uxi*xiy + ueta*etay + uzeta*zetay
              dudz = uxi*xiz + ueta*etaz + uzeta*zetaz
              dvdx = vxi*xix + veta*etax + vzeta*zetax
              dvdy = vxi*xiy + veta*etay + vzeta*zetay
              dvdz = vxi*xiz + veta*etaz + vzeta*zetaz
              dwdx = wxi*xix + weta*etax + wzeta*zetax
              dwdy = wxi*xiy + weta*etay + wzeta*zetay
              dwdz = wxi*xiz + weta*etaz + wzeta*zetaz

              divv = (dudx + dvdy + dwdz) * turcmp

              if (1.lt.0) then

!------------ prod is term P_d as by eqs. (31) and (32) of Liu and 
!               Zheng, JCP, 1996.
!               prod = tauxx*DUDX + tauyy*DVDY + tauzz*DWDZ +
!               tauxy*(dudy+dvdx) + tauxz*(dudz+dwdx) + tauyz*(dvdz+dwdy)
              txx = 2 * (dudx - divv/3)
              txy = (dudy + dvdx)
              txz = (dudz + dwdx)
              tyy = 2 * (dvdy - divv/3)
              tyz = (dvdz + dwdy)
              tzz = 2 * (dwdz - divv/3)
              prod = txx*dudx + txy*dudy + txz*dudz + &
                     txy*dvdx + tyy*dvdy + tyz*dvdz + &
                     txz*dwdx + tyz*dwdy + tzz*dwdz

              else if (1.gt.0) then

!------------ prod is term P_d as by eqs. (31) and (32) of Liu and 
!               Zheng, JCP, 1996.
!               This is alternative definition which may have less
!               round-off
              txx = dudx
              txy = (dudy + dvdx)/2
              txz = (dudz + dwdx)/2
              tyy = dvdy
              tyz = (dvdz + dwdy)/2
              tzz = dwdz
              prod =  2 * (txx**2+tyy**2+tzz**2 + &
                2*(txy**2+txz**2+tyz**2) - divv**2/3)

              end if

              prd   (i,j,k,n) = prod
              divvel(i,j,k,n) = divv
              delpls(i,j,k,n) = fact * max(divv,0.d0)

           end do
          end do
        end do
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine srckw(nl,q,mut,prd,divvel,res,vol)
!-----------------------------------------------------------------------
!     K Omega model.
!     Calculate source terms of K and Omega eqautions, and include them
!     in residual array
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,iblock,imax,jmax,kmax,iq,imut,ivol
      real (kind=cosa_real) q(*),mut(*),prd(*),divvel(*),res(*),vol(*)

      do iblock = 1,mynblocks
        imax   = i_imax     (iblock,nl)
        jmax   = j_jmax     (iblock,nl)
        kmax   = k_kmax     (iblock,nl)
        iq     = 1 + off_p3 (iblock,nl) * npde * dim5
        imut   = 1 + off_p3 (iblock,nl) *        dim5
        ivol   = 1 + off_p1 (iblock,nl)
        call srckw_b(q(iq),mut(imut),prd(imut),divvel(imut),res(iq), &
                     vol(ivol),nl,imax,jmax,kmax,npde,nharms)
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine srckw_b(q,mutg,prd,divvel,res,vol,nl,imax,jmax,kmax,npde, &
                         nharms)
!-----------------------------------------------------------------------
!     K Omega model.
!     Calculate source terms of K and Omega eqautions, and include them
!     in residual array
!     prod: term P_d as by eqs. (31) & (32) of Liu and Zheng, JCP, 1996
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,imax,jmax,kmax,npde,nharms
      integer(kind=cosa_int) i,j,k,n,nh
      real (kind=cosa_real) &
        q      (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
        res    (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
        mutg   (-1:imax+1,-1:jmax+1,-1:kmax+1     ,0:2*nharms), &
        prd    (-1:imax+1,-1:jmax+1,-1:kmax+1,     0:2*nharms), &
        divvel (-1:imax+1,-1:jmax+1,-1:kmax+1,     0:2*nharms), &
        vol    ( 0:imax  , 0:jmax  , 0:kmax                  )
      real (kind=cosa_real) rho,tke,ome,mut,divv,prod,prodw,tstdis,fact,mor

      fact = 2.d0/3.d0
      mor  = machfs/reyno

      do n = 0,2*nharms
        nh = n*hbmove
        do k = 1,kmax-1
          do j = 1,jmax-1
            do i = 1,imax-1

              rho = q(i,j,k,1,n)
              tke = q(i,j,k,6,n)
              ome = q(i,j,k,7,n)
              mut = mutg(i,j,k,n)

              divv  = divvel(i,j,k,n)
              prod  = prd(i,j,k,n)
              prodw = prd(i,j,k,n)

              if (lim_prodtke.and.((nl.eq.1).or.prd_all_lvls)) then
                tstdis = prdlim * bstrkw * rho * tke * ome
                if (mor*mut*prod.gt.tstdis) then
                  prod = tstdis / mor / mut
                  if (lim_prodomega) prodw = prod
                end if
              end if

              res(i,j,k,6,n) = res(i,j,k,6,n) - vol(i,j,k) * &
                (mor*mut*prod - fact*divv*rho*tke - &
                              bstrkw*rho*tke*ome)
              res(i,j,k,7,n) = res(i,j,k,7,n) - vol(i,j,k) * &
                (gkw*cmu*rho*prodw - gkw*fact*divv*rho*ome - &
                               bkw*rho*ome**2)
           end do
          end do
        end do
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine srckw_sst(nl,q,prd,divvel,res,mut,ttrm,vol)
!-----------------------------------------------------------------------
!     K Omega BSL/SST model.
!     Calculate source terms of K and Omega eqautions,
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,iblock,imax,jmax,kmax,iq,ivmt,ivol,imut,ittrm
      real (kind=cosa_real) q(*),res(*),ttrm(*),vol(*),mut(*),prd(*),divvel(*)

      do iblock = 1,mynblocks
        imax   = i_imax     (iblock,nl)
        jmax   = j_jmax     (iblock,nl)
        kmax   = k_kmax     (iblock,nl)
        iq     = 1 + off_p3 (iblock,nl) * npde * dim5
        imut   = 1 + off_p3 (iblock,nl) *        dim5
        ittrm  = 1 + off_p3 (iblock,nl) * 3    * dim5
        ivol   = 1 + off_p1 (iblock,nl)
        call srckw_sst_b(q(iq),prd(imut),divvel(imut),res(iq),mut(imut), &
          ttrm(ittrm),vol(ivol),nl,imax,jmax,kmax,npde,nharms)
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine srckw_sst_b(q,prd,divvel,res,mutg,ttrm,vol,nl,imax,jmax, &
                             kmax,npde,nharms)
!-----------------------------------------------------------------------
!     K Omega BSL/SST model.
!     Calculate source terms of K and Omega eqautions,
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,imax,jmax,kmax,npde,nharms
      integer(kind=cosa_int) i,j,k,n
      real (kind=cosa_real) &
        q     (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
        res   (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
        mutg  (-1:imax+1,-1:jmax+1,-1:kmax+1,     0:2*nharms), &
        prd   (-1:imax+1,-1:jmax+1,-1:kmax+1,     0:2*nharms), &
        divvel(-1:imax+1,-1:jmax+1,-1:kmax+1,     0:2*nharms), &
        ttrm  (-1:imax+1,-1:jmax+1,-1:kmax+1,3   ,0:2*nharms), &
        vol   ( 0:imax  , 0:jmax  , 0:kmax                  )
      real (kind=cosa_real) rho,tke,ome,mut,gav,bav,fact,mor,divv,prod,prodw, &
           tstdis,cdif

      fact = 2.d0/3.d0
      mor  = machfs/reyno

      do n = 0,2*nharms
        do k = 1,kmax-1
          do j = 1,jmax-1
            do i = 1,imax-1

              rho = q(i,j,k,1,n)
              tke = q(i,j,k,6,n)
              ome = q(i,j,k,7,n)
              mut = mutg(i,j,k,n)

              divv  = divvel(i,j,k,n)
              prod  = prd(i,j,k,n)
              prodw = prd(i,j,k,n)
              cdif  = ttrm(i,j,k,1,n)
              gav   = ttrm(i,j,k,2,n)
              bav   = ttrm(i,j,k,3,n)

              if (lim_prodtke.and.((nl.eq.1).or.prd_all_lvls)) then
                tstdis = prdlim * bstrkw * rho * tke * ome
                if (mor*mut*prod.gt.tstdis) then
                  prod = tstdis / mor / mut
                  if (lim_prodomega) prodw = prod
                end if
              end if

              res(i,j,k,6,n) = res(i,j,k,6,n) - vol(i,j,k)* &
                             (mor*mut*prod - fact*divv*rho*tke - &
                              bstrkw*rho*tke*ome)

!msc          res(i,j,k,7,n) = res(i,j,k,7,n) - vol(i,j,k)* &
!msc                         (gav*cmu*rho*prodw - &
!msc                          gav*fact*divv*(rho**2)*tke/mut/mor - &
!msc                          bav*rho*ome**2 + cdif)

              res(i,j,k,7,n) = res(i,j,k,7,n) - vol(i,j,k)* &
                             (gav*cmu*rho*prodw - gav*fact*divv*rho*ome - &
                              bav*rho*ome**2 + cdif)

           end do
          end do
        end do
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine src_uns(nl,q,vol,res)
!---------------------------------------------------------------------       
       
      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,iblock,imax,jmax,kmax,iq,ivol
      real (kind=cosa_real) q(*),res(*),vol(*)

      do iblock = 1,mynblocks
        imax   = i_imax(iblock,nl)
        jmax   = j_jmax(iblock,nl)
        kmax   = k_kmax(iblock,nl)
        iq     = 1 + off_p3 (iblock,nl) * npde
        ivol   = 1 + off_p1 (iblock,nl)
        call src_uns_b(q(iq),vol(ivol),res(iq),imax,jmax,kmax,npde)
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine src_uns_b(q,vol,res,imax,jmax,kmax,npde)
!---------------------------------------------------------------------       
       
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,npde
      integer(kind=cosa_int) i,j,k,ipde
      real (kind=cosa_real) &
           q  (-1:imax+1,-1:jmax+1,-1:kmax+1,npde), &
           res(-1:imax+1,-1:jmax+1,-1:kmax+1,npde), &
           vol( 0:imax  , 0:jmax  , 0:kmax       )

      if ((itime.eq.1).and.(.not.dual_prop_start)) then

        do ipde = 1,npde
          do k = 1,kmax-1
            do j = 1,jmax-1
              do i = 1,imax-1
                res(i,j,k,ipde) = res(i,j,k,ipde) + &
                                  q(i,j,k,ipde)*vol(i,j,k) / dt
              end do
            end do
          end do
        end do

      else

        do ipde = 1,npde
          do k = 1,kmax-1
            do j = 1,jmax-1
              do i = 1,imax-1
                res(i,j,k,ipde) = res(i,j,k,ipde) + &
                                  3*q(i,j,k,ipde)*vol(i,j,k) / 2 / dt
              end do
            end do
          end do
        end do

      end if

      return
      end

!-----------------------------------------------------------------------
      subroutine rhs_uns(nl,rhs,q1,q2,vol)
!---------------------------------------------------------------------       
       
      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,imax,jmax,kmax
      integer(kind=cosa_int) iblock,iq,ivol
      real (kind=cosa_real) rhs(*),q1(*),q2(*),vol(*)

      do iblock = 1,mynblocks
        imax   = i_imax(iblock,nl)
        jmax   = j_jmax(iblock,nl)
        kmax   = k_kmax(iblock,nl)
        iq     = 1 + off_p3 (iblock,nl) * npde
        ivol   = 1 + off_p1 (iblock,nl)
        call rhs_uns_b(rhs(iq),q1(iq),q2(iq),vol(ivol),imax,jmax,kmax, &
                           npde)
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine rhs_uns_b(rhs,q1,q2,vol,imax,jmax,kmax,npde)
!---------------------------------------------------------------------       
      
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,npde
      integer(kind=cosa_int) i,j,k,ipde
      real (kind=cosa_real) &
           rhs(-1:imax+1,-1:jmax+1,-1:kmax+1,npde), &
           q1 (-1:imax+1,-1:jmax+1,-1:kmax+1,npde), &
           q2 (-1:imax+1,-1:jmax+1,-1:kmax+1,npde), &
           vol( 0:imax  , 0:jmax  ,0:kmax        )

!     q1: q^n
!     q2: q^(n-1)

      if ((itime.eq.1).and.(.not.dual_prop_start)) then

        do ipde = 1,npde
          do k = 1,kmax-1
            do j = 1,jmax-1
              do i = 1,imax-1
                rhs(i,j,k,ipde) = q1(i,j,k,ipde) * vol(i,j,k) / dt
              end do
            end do
          end do
        end do

      else

        do ipde = 1,npde
          do k = 1,kmax-1
            do j = 1,jmax-1
              do i = 1,imax-1
                rhs(i,j,k,ipde) = &
                  (2*q1(i,j,k,ipde) - q2(i,j,k,ipde)/2) * vol(i,j,k) /dt
              end do
            end do
          end do
        end do

       end if

       return
       end

!-----------------------------------------------------------------------
      subroutine src_hb(nl,q,vol,res,qhb)
!---------------------------------------------------------------------       
       
      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,iblock,imax,jmax,kmax,iq,ivol
      real (kind=cosa_real) q(*),res(*),vol(*),qhb(*)

      do iblock = 1,mynblocks
        imax   = i_imax     (iblock,nl)
        jmax   = j_jmax     (iblock,nl)
        kmax   = k_kmax     (iblock,nl)
        iq     = 1 + off_p3 (iblock,nl) * npde * dim5
        ivol   = 1 + off_p1 (iblock,nl)
        call src_bhb(q(iq),vol(ivol),res(iq),qhb(iq),imax,jmax,kmax, &
                     npde,nharms)
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine src_bhb(q,vol,res,qhb,imax,jmax,kmax,npde,nharms)
!---------------------------------------------------------------------       
       
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,npde,nharms
      integer(kind=cosa_int) i,j,k,ipde,ni,nj
      real (kind=cosa_real) &
           res(-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           q  (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           qhb(-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           vol( 0:imax  , 0:jmax  , 0:kmax                  )

      do ni = 0,2*nharms
        do ipde = 1,npde
          do k = 1,kmax-1
            do j = 1,jmax-1
              do i = 1,imax-1
                qhb(i,j,k,ipde,ni) = 0.d0
              end do
            end do
          end do
        end do
      end do

      do ni = 0,2*nharms
        do nj = 0,2*nharms
          do ipde = 1,npde
            do k = 1,kmax-1
              do j = 1,jmax-1
                do i = 1,imax-1
                  qhb(i,j,k,ipde,ni) = qhb(i,j,k,ipde,ni) + &
                    omega*vol(i,j,k)*dhb(ni,nj)*q(i,j,k,ipde,nj)
                end do
              end do
            end do
          end do
        end do
      end do

      do ni = 0,2*nharms
        do ipde = 1,npde
          do k = 1,kmax-1
            do j = 1,jmax-1
              do i = 1,imax-1
                res(i,j,k,ipde,ni) = res(i,j,k,ipde,ni) + &
                  qhb(i,j,k,ipde,ni)
              end do
            end do
          end do
        end do
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine src_relf(nl,q,vol,res)
!---------------------------------------------------------------------       
       
      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,iblock,imax,jmax,kmax,iq,ivol
      real (kind=cosa_real) q(*),res(*),vol(*)

      do iblock = 1,mynblocks
        imax   = i_imax(iblock,nl)
        jmax   = j_jmax(iblock,nl)
        kmax   = k_kmax(iblock,nl)
        iq     = 1 + off_p3 (iblock,nl) * npde * dim5
        ivol   = 1 + off_p1 (iblock,nl)
        call src_relf_b(q(iq),vol(ivol),res(iq),omegas,imax,jmax,kmax, &
                        npde,nharms)
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine src_relf_b(q,vol,res,omegas,imax,jmax,kmax,npde,nharms)
!---------------------------------------------------------------------       
       
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,npde,nharms
      integer(kind=cosa_int) i,j,k,ipde,n
      real (kind=cosa_real) &
           q  (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           res(-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           vol( 0:imax  , 0:jmax  , 0:kmax       )
      real (kind=cosa_real) omegas

      do n = 0,2*nharms
        do k = 1,kmax-1
          do j = 1,jmax-1
            do i = 1,imax-1
              res(i,j,k,2,n) = res(i,j,k,2,n) - &
                q(i,j,k,1,n)*omegas*q(i,j,k,3,n)*vol(i,j,k)
              res(i,j,k,3,n) = res(i,j,k,3,n) + &
                q(i,j,k,1,n)*omegas*q(i,j,k,2,n)*vol(i,j,k)
            end do
          end do
        end do
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine lim_prod_mg(nl,q,prd,mut)
!-----------------------------------------------------------------------
!     K Omega models: limit prod term
!     prod: term P_d as by eqs. (31) & (32) of Liu and Zheng, JCP, 1996
!     05/08/2021: this ruting limits prod when one calculate muT on all
!                 grid terms and restrics prod from fine grid. This
!                 implies that k and omega production terms must be both
!                 limited, as prod appears both in k and onega eqs, i.e.
!                 there is no separate prodw (array) at present.
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,iblock,imax,jmax,kmax,iq,ivmt,imut
      real (kind=cosa_real) q(*),mut(*),prd(*)

      do iblock = 1,mynblocks
        imax   = i_imax     (iblock,nl)
        jmax   = j_jmax     (iblock,nl)
        kmax   = k_kmax     (iblock,nl)
        iq     = 1 + off_p3 (iblock,nl) * npde * dim5
        imut   = 1 + off_p3 (iblock,nl) *        dim5
        call lim_prod_mg_b(q(iq),prd(imut),mut(imut), &
          imax,jmax,kmax,npde,nharms)
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine lim_prod_mg_b(q,prd,mutg,imax,jmax,kmax,npde,nharms)
!-----------------------------------------------------------------------
!     K Omega models: limit prod term
!     prod: term P_d as by eqs. (31) & (32) of Liu and Zheng, JCP, 1996
!     05/08/2021: this ruting limits prod when one calculate muT on all
!                 grid terms and restrics prod from fine grid. This
!                 implies that k and omega production terms must be both
!                 limited, as prod appears both in k and onega eqs, i.e.
!                 there is no separate prodw (array) at present.
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,imax,jmax,kmax,npde,nharms
      integer(kind=cosa_int) i,j,k,n
      real (kind=cosa_real) &
        q     (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
        mutg  (-1:imax+1,-1:jmax+1,-1:kmax+1,     0:2*nharms), &
        prd   (-1:imax+1,-1:jmax+1,-1:kmax+1,     0:2*nharms)
      real (kind=cosa_real) rho,tke,ome,mut,mor,prod,tstdis

      mor  = machfs/reyno

      do n = 0,2*nharms
        do k = 1,kmax-1
          do j = 1,jmax-1
            do i = 1,imax-1

              rho = q(i,j,k,1,n)
              tke = q(i,j,k,6,n)
              ome = q(i,j,k,7,n)
              mut = mutg(i,j,k,n)

              prod  = prd(i,j,k,n)

              tstdis = prdlim * bstrkw * rho * tke * ome
              if (mor*mut*prod.gt.tstdis) then
                prd(i,j,k,n) = tstdis / mor / mut
              end if

           end do
          end do
        end do
      end do

      return
      end
