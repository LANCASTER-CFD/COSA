!-----------------------------------------------------------------------
      subroutine lsp(nl,q,lr,prec,preci,ipkpi,cutoff_pgr,cutoff_vis, &
                     dtvol,dtvolt,vol,ttrm,delpls,xdot,ydot,zdot,ar)
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none
!

      integer(kind=cosa_int) nl,iblock,imax,jmax,kmax,iq,ilr,iprec,icof,idtv,ivol, &
          ittrm,idlpl,ixyz
      real(kind=cosa_real) q(*),lr(*),prec(*),preci(*),ipkpi(*),cutoff_pgr(*), &
          cutoff_vis(*),vol(*),dtvol(*),dtvolt(*),ttrm(*),delpls(*), &
          xdot(*),ydot(*),zdot(*)
      real(kind=cosa_real) ar

      do iblock = 1,mynblocks
        imax   = i_imax     (iblock,nl)
        jmax   = j_jmax     (iblock,nl)
        kmax   = k_kmax     (iblock,nl)
        iq     = 1 + off_p3 (iblock,nl) * npde * dim5
        iprec  = 1 + off_m1 (iblock,nl) * npde * npde * dim5
        icof   = 1 + off_m1 (iblock,nl) *        dim5
        idtv   = 1 + off_p1 (iblock,nl) *        dim5
        ivol   = 1 + off_p1 (iblock,nl)
        ixyz   = 1 + off_p2 (iblock,nl)        * dim5h
        ittrm  = 1 + off_p3 (iblock,nl) * 3    * dim5
        idlpl  = 1 + off_p3 (iblock,nl) *        dim5
        call blsp(nl,q(iq),lr(iq),prec(iprec),preci(iprec), &
          ipkpi(iprec),cutoff_pgr(icof),cutoff_vis(icof),dtvol(idtv), &
          dtvolt(idtv),vol(ivol),ttrm(ittrm),delpls(idlpl),xdot(ixyz), &
          ydot(ixyz),zdot(ixyz),ar,imax,jmax,kmax,npde,nharms)
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine blsp(nl,q,lr,prec,preci,ipkpi,cutoff_pgr,cutoff_vis, &
                  dtvol,dtvolt,vol,ttrm,delpls,xdot,ydot,zdot,ar,imax, &
                  jmax,kmax,npde,nharms)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none
!
      integer(kind=cosa_int) imax,jmax,kmax,npde,nharms
      integer(kind=cosa_int) im,jm,km,i,j,k,l,ii,jj,kk,ii2,jj2,n,itest,nl,nh, &
                ipde,ipde1,ipde2
      real(kind=cosa_real) &
          q         (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
          lr        (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
          prec      (imax-1,jmax-1,kmax-1,npde,npde,0:2*nharms), &
          preci     (imax-1,jmax-1,kmax-1,npde,npde,0:2*nharms), &
          dtvol     (0:imax,0:jmax,0:kmax,0:2*nharms), &
          dtvolt    (0:imax,0:jmax,0:kmax,0:2*nharms), &
          ipkpi     (imax-1,jmax-1,kmax-1,npde,npde,0:2*nharms), &
          cutoff_pgr(imax-1,jmax-1,kmax-1,0:2*nharms), &
          vol       (0:imax,0:jmax,0:kmax), &
          cutoff_vis(imax-1,jmax-1,kmax-1,0:2*nharms), &
          ttrm      (-1:imax+1,-1:jmax+1,-1:kmax+1,3   ,0:2*nharms), &
          delpls    (-1:imax+1,-1:jmax+1,-1:kmax+1,     0:2*nharms), &
          xdot      (0:imax+1,0:jmax+1,0:kmax+1,0:2*nharms*hbmove), &
          ydot      (0:imax+1,0:jmax+1,0:kmax+1,0:2*nharms*hbmove), &
          zdot      (0:imax+1,0:jmax+1,0:kmax+1,0:2*nharms*hbmove)

      real(kind=cosa_real) rho,u,v,w,tke,ome,p,e0,q2,qp2,term1, &
          dtau,dtaut,h0,a2,a,mp2,ump2,ar,fact1,fact2,fact3,fact4, &
          a66,a77,entry,gam1,mp2m1,akom(2,2),tmp(7,7),prod(7,7)

      im = imax-1
      jm = jmax-1
      km = kmax-1

      gam1  = gamma-1

!---- 05/07/09: new expression of preconditioner that depends only on
!               mp2 and no longer on rhopp 
!               (Philipp Wessler development, to be verified)

      do n = 0,2*nharms
        nh = n*hbmove

        do k = 1,km
          do j = 1,jm
            do i = 1,im

            rho   = q(i,j,k,1,n)
            u     = q(i,j,k,2,n)/rho
            v     = q(i,j,k,3,n)/rho
            w     = q(i,j,k,4,n)/rho
            q2    = u**2+v**2+w**2
            e0    = q(i,j,k,5,n)/rho
            if (kom.or.kom_bsl.or.kom_sst) then
                tke   = q(i,j,k,6,n)/rho
            else
                tke   = 0.d0
            end if
            p     = (gamma-1)*(e0-0.5*q2-tke)*rho
            a2    = gamma*p/rho
            a     = sqrt(a2)

            if (moving.and.(lsp_vfr.eq.1)) then
              qp2 = (u-xdot(i,j,k,nh))**2 + (v-ydot(i,j,k,nh))**2 + &
                    (w-zdot(i,j,k,nh))**2 
            else
              qp2 = q2
            end if

            mp2   = qp2/a2
            ump2  = (uprv / a)**2
            mp2   = dmin1(1.d0, &
                          dmax1(mp2,ump2,cutoff_pgr(i,j,k,n), &
                                cutoff_vis(i,j,k,n),epsp(nl)**2))
            mp2m1 = mp2-1

            prec(i,j,k,1,1,n) = gam1*mp2m1*q2/2/a2+1
            prec(i,j,k,1,2,n) = -u*mp2m1*gam1/a2
            prec(i,j,k,1,3,n) = -v*mp2m1*gam1/a2
            prec(i,j,k,1,4,n) = -w*mp2m1*gam1/a2
            prec(i,j,k,1,5,n) = mp2m1*gam1/a2

            prec(i,j,k,2,1,n) = u*gam1*mp2m1*q2/2/a2
            prec(i,j,k,2,2,n) = 1-u**2*mp2m1*gam1/a2
            prec(i,j,k,2,3,n) = -u*v*mp2m1*gam1/a2
            prec(i,j,k,2,4,n) = -u*w*mp2m1*gam1/a2
            prec(i,j,k,2,5,n) = u*mp2m1*gam1/a2

            prec(i,j,k,3,1,n) = v*gam1*mp2m1*q2/2/a2
            prec(i,j,k,3,2,n) = -u*v*mp2m1*gam1/a2
            prec(i,j,k,3,3,n) = 1-v**2*mp2m1*gam1/a2
            prec(i,j,k,3,4,n) = -v*w*mp2m1*gam1/a2
            prec(i,j,k,3,5,n) = v*mp2m1*gam1/a2

            prec(i,j,k,4,1,n) = w*gam1*mp2m1*q2/2/a2
            prec(i,j,k,4,2,n) = -u*w*mp2m1*gam1/a2
            prec(i,j,k,4,3,n) = -v*w*mp2m1*gam1/a2
            prec(i,j,k,4,4,n) = 1-w**2*mp2m1*gam1/a2
            prec(i,j,k,4,5,n) = w*mp2m1*gam1/a2

            prec(i,j,k,5,1,n) = q2*mp2m1*(2*a2+gam1*(q2+2*tke))/4/a2
            prec(i,j,k,5,2,n) = -u*mp2m1*(2*a2+gam1*(q2+2*tke))/2/a2
            prec(i,j,k,5,3,n) = -v*mp2m1*(2*a2+gam1*(q2+2*tke))/2/a2
            prec(i,j,k,5,4,n) = -w*mp2m1*(2*a2+gam1*(q2+2*tke))/2/a2
            prec(i,j,k,5,5,n) = mp2+gam1*(q2+2*tke)*mp2m1/2/a2

            if (kom.or.kom_bsl.or.kom_sst) then

              ome   = q(i,j,k,7,n)/rho

              prec(i,j,k,1,6,n) = -mp2m1*gam1/a2
              prec(i,j,k,1,7,n) = 0

              prec(i,j,k,2,6,n) = -u*mp2m1*gam1/a2
              prec(i,j,k,2,7,n) = 0

              prec(i,j,k,3,6,n) = -v*mp2m1*gam1/a2
              prec(i,j,k,3,7,n) = 0

              prec(i,j,k,4,6,n) = -w*mp2m1*gam1/a2
              prec(i,j,k,4,7,n) = 0

              prec(i,j,k,5,6,n) = 1-mp2-gam1*(q2+2*tke)*mp2m1/2/a2
              prec(i,j,k,5,7,n) = 0

              prec(i,j,k,6,1,n) = tke*gam1*mp2m1*q2/2/a2
              prec(i,j,k,6,2,n) = -u*tke*mp2m1*gam1/a2
              prec(i,j,k,6,3,n) = -v*tke*mp2m1*gam1/a2
              prec(i,j,k,6,4,n) = -w*tke*mp2m1*gam1/a2
              prec(i,j,k,6,5,n) = tke*mp2m1*gam1/a2
              prec(i,j,k,6,6,n) = 1-tke*mp2m1*gam1/a2
              prec(i,j,k,6,7,n) = 0

              prec(i,j,k,7,1,n) = ome*gam1*mp2m1*q2/2/a2
              prec(i,j,k,7,2,n) = -u*ome*mp2m1*gam1/a2
              prec(i,j,k,7,3,n) = -v*ome*mp2m1*gam1/a2
              prec(i,j,k,7,4,n) = -w*ome*mp2m1*gam1/a2
              prec(i,j,k,7,5,n) = ome*mp2m1*gam1/a2
              prec(i,j,k,7,6,n) = -ome*mp2m1*gam1/a2
              prec(i,j,k,7,7,n) = 1

            end if
          
            end do
          end do
        end do
      end do

      if (debug) then

        do n = 0,2*nharms
          nh = n*hbmove

          do k = 1,km 
            do j = 1,jm
               do i = 1,im

               rho   = q(i,j,k,1,n)
               u     = q(i,j,k,2,n)/rho
               v     = q(i,j,k,3,n)/rho
               w     = q(i,j,k,4,n)/rho
               q2    = u**2+v**2+w**2
               e0    = q(i,j,k,5,n)/rho
               if (kom.or.kom_bsl.or.kom_sst) then
                 tke   = q(i,j,k,6,n)/rho
               else
                 tke   = 0.d0
               end if
               p     = gam1*(e0-0.5*q2-tke)*rho
               a2    = gamma*p/rho
               a     = sqrt(a2)

               if (moving.and.(lsp_vfr.eq.1)) then
                 qp2 = (u-xdot(i,j,k,nh))**2 + (v-ydot(i,j,k,nh))**2 + &
                       (w-zdot(i,j,k,nh))**2 
               else
                 qp2 = q2
               end if

               mp2   = qp2/a2
               ump2  = (uprv / a)**2
               mp2   = dmin1(1.d0, &
                            dmax1(mp2,ump2,cutoff_pgr(i,j,k,n), &
                                  cutoff_vis(i,j,k,n),epsp(nl)**2))
               mp2m1 = mp2-1

               preci(i,j,k,1,1,n) = -gam1*mp2m1*q2/2/a2/mp2+1
               preci(i,j,k,1,2,n) = u*mp2m1*gam1/a2/mp2
               preci(i,j,k,1,3,n) = v*mp2m1*gam1/a2/mp2
               preci(i,j,k,1,4,n) = w*mp2m1*gam1/a2/mp2
               preci(i,j,k,1,5,n) = -mp2m1*gam1/a2/mp2

               preci(i,j,k,2,1,n) = -u*gam1*mp2m1*q2/2/a2/mp2
               preci(i,j,k,2,2,n) = 1+u**2*mp2m1*gam1/a2/mp2
               preci(i,j,k,2,3,n) = u*v*mp2m1*gam1/a2/mp2
               preci(i,j,k,2,4,n) = u*w*mp2m1*gam1/a2/mp2
               preci(i,j,k,2,5,n) = -u*mp2m1*gam1/a2/mp2

               preci(i,j,k,3,1,n) = -v*gam1*mp2m1*q2/2/a2/mp2
               preci(i,j,k,3,2,n) = u*v*mp2m1*gam1/a2/mp2
               preci(i,j,k,3,3,n) = 1+v**2*mp2m1*gam1/a2/mp2
               preci(i,j,k,3,4,n) = v*w*mp2m1*gam1/a2/mp2
               preci(i,j,k,3,5,n) = -v*mp2m1*gam1/a2/mp2

               preci(i,j,k,4,1,n) = -w*gam1*mp2m1*q2/2/a2/mp2
               preci(i,j,k,4,2,n) = u*w*mp2m1*gam1/a2/mp2
               preci(i,j,k,4,3,n) = v*w*mp2m1*gam1/a2/mp2
               preci(i,j,k,4,4,n) = 1+w**2*mp2m1*gam1/a2/mp2
               preci(i,j,k,4,5,n) = -w*mp2m1*gam1/a2/mp2

               preci(i,j,k,5,1,n) = &
                -q2*mp2m1*(2*a2+gam1*(q2+2*tke))/4/a2/mp2
               preci(i,j,k,5,2,n) = u*mp2m1*(2*a2+gam1*(q2+2*tke))/2/a2/mp2
               preci(i,j,k,5,3,n) = v*mp2m1*(2*a2+gam1*(q2+2*tke))/2/a2/mp2
               preci(i,j,k,5,4,n) = w*mp2m1*(2*a2+gam1*(q2+2*tke))/2/a2/mp2
               preci(i,j,k,5,5,n) = -gam1*(q2/2+tke)*mp2m1/a2/mp2+1/mp2

               if (kom.or.kom_bsl.or.kom_sst) then

                 ome   = q(i,j,k,7,n)/rho

                 preci(i,j,k,1,6,n) = mp2m1*gam1/a2/mp2
                 preci(i,j,k,1,7,n) = 0
 
                 preci(i,j,k,2,6,n) = u*mp2m1*gam1/a2/mp2
                 preci(i,j,k,2,7,n) = 0

                 preci(i,j,k,3,6,n) = v*mp2m1*gam1/a2/mp2
                 preci(i,j,k,3,7,n) = 0

                 preci(i,j,k,4,6,n) = w*mp2m1*gam1/a2/mp2
                 preci(i,j,k,4,7,n) = 0

                 preci(i,j,k,5,6,n) = 1-1/mp2+(q2/2+tke)*mp2m1*gam1/a2/mp2
                 preci(i,j,k,5,7,n) = 0

                 preci(i,j,k,6,1,n) = -tke*gam1*mp2m1*q2/2/a2/mp2
                 preci(i,j,k,6,2,n) = u*tke*mp2m1*gam1/a2/mp2
                 preci(i,j,k,6,3,n) = v*tke*mp2m1*gam1/a2/mp2
                 preci(i,j,k,6,4,n) = w*tke*mp2m1*gam1/a2/mp2
                 preci(i,j,k,6,5,n) = -tke*mp2m1*gam1/a2/mp2
                 preci(i,j,k,6,6,n) = 1+tke*mp2m1*gam1/a2/mp2
                 preci(i,j,k,6,7,n) = 0

                 preci(i,j,k,7,1,n) = -ome*gam1*mp2m1*q2/2/a2/mp2
                 preci(i,j,k,7,2,n) = u*ome*mp2m1*gam1/a2/mp2
                 preci(i,j,k,7,3,n) = v*ome*mp2m1*gam1/a2/mp2
                 preci(i,j,k,7,4,n) = w*ome*mp2m1*gam1/a2/mp2
                 preci(i,j,k,7,5,n) = -ome*mp2m1*gam1/a2/mp2
                 preci(i,j,k,7,6,n) = ome*mp2m1*gam1/a2/mp2
                 preci(i,j,k,7,7,n) = 1

              end if

              end do
            end do
          end do
        end do

      end if

      if ((.not.rkimp).and.(kom.or.kom_bsl.or.kom_sst)) then

!------ TURBULENT STEADY OR FERK TD OR FERK HB CASE

!       ipkpi: INV (I + fact * P * A). This stands for
!       Inverse of (Identity Plus fact times Preconditioner times A),
!       where fact=\alpha_k times \dtau (Eqn. (17), ASME GT2014-25562).

        do n = 0,2*nharms
          nh = n*hbmove

          do k= 1,km
            do j = 1,jm
              do i = 1,im

              rho   = q(i,j,k,1,n)
              u     = q(i,j,k,2,n)/rho
              v     = q(i,j,k,3,n)/rho
              w     = q(i,j,k,4,n)/rho
              q2    = u**2+v**2+w**2
              e0    = q(i,j,k,5,n)/rho
              tke   = q(i,j,k,6,n)/rho
              ome   = q(i,j,k,7,n)/rho
              p     = gam1*(e0-0.5*q2-tke)*rho
              a2    = gamma*p/rho
              a     = sqrt(a2)

              if (moving.and.(lsp_vfr.eq.1)) then
                qp2 = (u-xdot(i,j,k,nh))**2 + (v-ydot(i,j,k,nh))**2 + &
                      (w-zdot(i,j,k,nh))**2 
              else
                qp2 = q2
              end if

              mp2   = qp2/a2
              ump2  = (uprv / a)**2
              mp2   = dmin1(1.d0, &
                            dmax1(mp2,ump2,cutoff_pgr(i,j,k,n), &
                                  cutoff_vis(i,j,k,n),epsp(nl)**2))
              mp2m1 = mp2-1
              dtau  = dtvol (i,j,k,n)*vol(i,j,k)
              dtaut = dtvolt(i,j,k,n)*vol(i,j,k)

              a66   = delpls(i,j,k,n)+bstrkw*ome  
              a77   = ttrm(i,j,k,2,n)*delpls(i,j,k,n)+2*ttrm(i,j,k,3,n)*ome
              fact1 = ar+ar**2*a77*dtaut

              term1= (1+ar*dtaut*a77)* &
                (a2-ar*dtaut*(-a2+tke*mp2m1*gam1)*a66) - &
                ar*bstrkw*dtaut*tke*ome*mp2m1*gam1       

              ipkpi(i,j,k,1,1,n) = 1
              ipkpi(i,j,k,1,2,n) = 0
              ipkpi(i,j,k,1,3,n) = 0
              ipkpi(i,j,k,1,4,n) = 0
              ipkpi(i,j,k,1,5,n) = 0
              ipkpi(i,j,k,1,6,n) = mp2m1*gam1*dtau*a66*fact1/term1
              ipkpi(i,j,k,1,7,n) = ar*bstrkw*tke*dtau*mp2m1*gam1/term1

              ipkpi(i,j,k,2,1,n) = 0
              ipkpi(i,j,k,2,2,n) = 1
              ipkpi(i,j,k,2,3,n) = 0
              ipkpi(i,j,k,2,4,n) = 0
              ipkpi(i,j,k,2,5,n) = 0
              ipkpi(i,j,k,2,6,n) = u*mp2m1*gam1*dtau*a66*fact1/term1
              ipkpi(i,j,k,2,7,n) = &
                u*ar*bstrkw*tke*dtau*mp2m1*gam1/term1

              ipkpi(i,j,k,3,1,n) = 0
              ipkpi(i,j,k,3,2,n) = 0
              ipkpi(i,j,k,3,3,n) = 1
              ipkpi(i,j,k,3,4,n) = 0
              ipkpi(i,j,k,3,5,n) = 0
              ipkpi(i,j,k,3,6,n) = v*mp2m1*gam1*dtau*a66*fact1/term1
              ipkpi(i,j,k,3,7,n) = &
                v*ar*bstrkw*tke*dtau*mp2m1*gam1/term1

              ipkpi(i,j,k,4,1,n) = 0
              ipkpi(i,j,k,4,2,n) = 0
              ipkpi(i,j,k,4,3,n) = 0
              ipkpi(i,j,k,4,4,n) = 1
              ipkpi(i,j,k,4,5,n) = 0
              ipkpi(i,j,k,4,6,n) = w*mp2m1*gam1*dtau*a66*fact1/term1
              ipkpi(i,j,k,4,7,n) = &
                w*ar*bstrkw*tke*dtau*mp2m1*gam1/term1

              ipkpi(i,j,k,5,1,n) = 0
              ipkpi(i,j,k,5,2,n) = 0
              ipkpi(i,j,k,5,3,n) = 0
              ipkpi(i,j,k,5,4,n) = 0
              ipkpi(i,j,k,5,5,n) = 1
              ipkpi(i,j,k,5,6,n) = &
                (a2/gam1+q2/2+tke)*mp2m1*gam1*dtau*a66*fact1/term1
              ipkpi(i,j,k,5,7,n) = &
                (a2/gam1+q2/2+tke)*ar*bstrkw*tke*dtau*mp2m1*gam1/term1

              ipkpi(i,j,k,6,1,n) = 0
              ipkpi(i,j,k,6,2,n) = 0
              ipkpi(i,j,k,6,3,n) = 0
              ipkpi(i,j,k,6,4,n) = 0
              ipkpi(i,j,k,6,5,n) = 0
              ipkpi(i,j,k,6,6,n) = &
                (a2+a2*ar*dtaut*a77-ar*bstrkw*dtaut*tke*ome*mp2m1*gam1)/ &
                term1
              ipkpi(i,j,k,6,7,n) = &
                ar*bstrkw*tke*dtaut*(-a2+tke*mp2m1*gam1)/term1

              ipkpi(i,j,k,7,1,n) = 0
              ipkpi(i,j,k,7,2,n) = 0
              ipkpi(i,j,k,7,3,n) = 0
              ipkpi(i,j,k,7,4,n) = 0
              ipkpi(i,j,k,7,5,n) = 0
              ipkpi(i,j,k,7,6,n) = ar*dtaut*ome*a66*mp2m1*gam1/term1
              ipkpi(i,j,k,7,7,n) = &
                a2*(1-ar*dtaut*a66*(-1+tke*mp2m1*gam1/a2))/term1

!------------ lr is the last row of ipkp [=(I+ar*dtau*Pc*A)]
              lr(i,j,k,1,n) = 0
              lr(i,j,k,2,n) = 0 
              lr(i,j,k,3,n) = 0
              lr(i,j,k,4,n) = 0
              lr(i,j,k,5,n) = 0
              lr(i,j,k,6,n) = -ar*dtaut*ome*mp2m1*gam1*a66/a2
              lr(i,j,k,7,n) = &
                ar*dtaut*(a77-bstrkw*tke*ome*mp2m1*gam1/a2)+1
 
              end do
            end do
          end do
        end do  

      else if (rkimp.and.(kom.or.kom_bsl.or.kom_sst)) then

!------ TURBULENT PIRK TD CASE

!       ipkpi: INV (I + fact * P * (1.5/dt*I + A)). This stands for
!       Inverse of (Identity Plus fact times Preconditioner times 
!       (1.5/dt*I + A)),
!       where fact=\alpha_k times \dtau (Eqn. (17), ASME GT2014-25562).

        fact1 = 1.5d0
        if ((itime.eq.1).and.(.not.dual_prop_start)) then
          fact1 = 1.0d0
        end if

        do n = 0,2*nharms
          nh = n*hbmove

          do k = 1,km
            do j = 1,jm
              do i = 1,im

              rho   = q(i,j,k,1,n)
              u     = q(i,j,k,2,n)/rho
              v     = q(i,j,k,3,n)/rho
              w     = q(i,j,k,4,n)/rho
              q2    = u**2+v**2+w**2
              e0    = q(i,j,k,5,n)/rho
              tke   = q(i,j,k,6,n)/rho
              ome   = q(i,j,k,7,n)/rho
              p     = gam1*(e0-0.5*q2-tke)*rho
              a2    = gamma*p/rho
              a     = sqrt(a2)

              if (moving.and.(lsp_vfr.eq.1)) then
                qp2 = (u-xdot(i,j,k,nh))**2 + (v-ydot(i,j,k,nh))**2 + &
                      (w-zdot(i,j,k,nh))**2 
              else
                qp2 = q2
              end if

              mp2   = qp2/a2
              ump2  = (uprv / a)**2
              mp2   = dmin1(1.d0, &
                            dmax1(mp2,ump2,cutoff_pgr(i,j,k,n), &
                                  cutoff_vis(i,j,k,n),epsp(nl)**2))
              mp2m1 = mp2-1
              dtau  = dtvol (i,j,k,n)*vol(i,j,k)
              dtaut = dtvolt(i,j,k,n)*vol(i,j,k)

              fact2 = delpls(i,j,k,n)+bstrkw*ome+fact1/dt  
              fact3 = ttrm(i,j,k,2,n)*delpls(i,j,k,n) &
                      + 2*ttrm(i,j,k,3,n)*ome + fact1/dt
              fact4 = 1+fact1*ar*dtau/dt

              term1 = fact4**4*((1+ar*dtaut*fact3)*((1+ &
                      fact1*ar*dtau*mp2/dt)*(1+ar*dtaut*fact2) &
                      +ar*mp2m1*gam1*tke*(fact1*dtau/dt &
                      -dtaut*fact2)/a2)-ar*bstrkw*tke*mp2m1 &
                      *gam1*dtaut*ome*fact4/a2)

              ipkpi(i,j,k,1,1,n) = fact4**3*((1+ar*dtaut*fact3)* &
               (a2*dt*fact4*(1+ar*fact2*dtaut)+fact1*ar*dtau* &
               (a2/gam1-q2/2+tke)*mp2m1*gam1* &
               (1+ar*fact2*dtaut)-ar*dt*fact2*fact4*tke*mp2m1* &
               gam1*dtaut)/(a2*dt)-(ar*bstrkw*tke*mp2m1 &
               *gam1*dtaut*ome*fact4)/a2)/term1

              ipkpi(i,j,k,1,2,n) = fact1*fact4**3*u*ar*mp2m1*gam1 &
              *dtau*(1+ar*dtaut*fact2)*(1+ar*fact3*dtaut)/(a2*dt)/term1

              ipkpi(i,j,k,1,3,n) = fact1*fact4**3*v*ar*mp2m1*gam1 &
              *dtau*(1+ar*dtaut*fact2)*(1+ar*fact3*dtaut)/(a2*dt)/term1

              ipkpi(i,j,k,1,4,n) = fact1*fact4**3*w*ar*mp2m1*gam1 &
              *dtau*(1+ar*dtaut*fact2)*(1+ar*fact3*dtaut)/(a2*dt)/term1

              ipkpi(i,j,k,1,5,n) = -fact1*fact4**3*ar*mp2m1*gam1 &
              *dtau*(1+ar*dtaut*fact2)*(1+ar*fact3*dtaut)/(a2*dt)/term1

              ipkpi(i,j,k,1,6,n) = (ar*fact2*fact4**3*mp2m1*gam1* &
                dtau*(1+ar*fact3*dtaut))/a2/term1

              ipkpi(i,j,k,1,7,n) = (ar*bstrkw*fact4**3*tke*mp2m1* &
                gam1*dtau)/a2/term1
 
              ipkpi(i,j,k,2,1,n) = -fact1*fact4**3*ar*q2*u*mp2m1*gam1* &
                dtau*(1+ar*dtaut*fact2)*(1+ar*fact3*dtaut)/2/a2/dt/term1

              ipkpi(i,j,k,2,2,n) = fact4**3*((ar*fact3*dtaut+1)* &
                (1/(2*a2*dt)* &
                (1+ar*fact2*dtaut)*(2*a2*dt*fact4+fact1*ar*dtau*mp2m1 &
                *gam1*(2*(a2/gam1+tke)+2*u**2))-fact4* &
                ((tke*mp2m1*gam1)/a2)*ar*fact2*dtaut) - &
                (ar*bstrkw*fact4*mp2m1*gam1*dtaut*tke*ome)/a2)/term1

              ipkpi(i,j,k,2,3,n) = fact1*fact4**3*ar*u*v*mp2m1*gam1* &
                dtau*(1+ar*dtaut*fact2)*(1+ar*fact3*dtaut)/(a2*dt)/term1

              ipkpi(i,j,k,2,4,n) = fact1*fact4**3*ar*u*w*mp2m1*gam1* &
                dtau*(1+ar*dtaut*fact2)*(1+ar*fact3*dtaut)/(a2*dt)/term1

              ipkpi(i,j,k,2,5,n) = -fact1*fact4**3*ar*u*mp2m1*gam1* &
                dtau*(1+ar*dtaut*fact2)*(1+ar*fact3*dtaut)/(a2*dt)/term1

              ipkpi(i,j,k,2,6,n) = ar*fact2*fact4**4*mp2m1*gam1* &
                dtau*(1+ar*fact3*dtaut)*u/a2/term1

              ipkpi(i,j,k,2,7,n) = ar*bstrkw*fact4**4*tke*mp2m1* &
                gam1*dtau*u/a2/term1

              ipkpi(i,j,k,3,1,n) = -fact1*fact4**3*ar*q2*v*mp2m1*gam1* &
                dtau*(1+ar*dtaut*fact2)*(1+ar*fact3*dtaut)/2/a2/dt/term1

              ipkpi(i,j,k,3,2,n) = fact1*fact4**3*ar*u*v*mp2m1*gam1*dtau* &
                (1+ar*dtaut*fact2)*(1+ar*fact3*dtaut)/(a2*dt)/term1

              ipkpi(i,j,k,3,3,n) = fact4**3*((ar*fact3*dtaut+1)* &
                (1/(2*a2*dt)* &
                (1+ar*fact2*dtaut)*(2*a2*dt*fact4+fact1*ar*dtau*mp2m1 &
                *gam1*(2*(a2/gam1+tke)+2*v**2))-fact4* &
                (tke*mp2m1*gam1/a2)*ar*fact2*dtaut) &
                -(ar*bstrkw*fact4*mp2m1*gam1*dtaut*tke*ome)/a2)/term1

              ipkpi(i,j,k,3,4,n) = fact1*fact4**3*ar*u*w*mp2m1*gam1*dtau* &
                (1+ar*dtaut*fact2)*(1+ar*fact3*dtaut)/(a2*dt)/term1

              ipkpi(i,j,k,3,5,n) = -fact1*fact4**3*ar*v*mp2m1*gam1*dtau* &
                (1+ar*dtaut*fact2)*(1+ar*fact3*dtaut)/(a2*dt)/term1

              ipkpi(i,j,k,3,6,n) = (ar*fact2*fact4**4*mp2m1*gam1* &
                dtau*(1+ar*fact3*dtaut)*v)/a2/term1

              ipkpi(i,j,k,3,7,n) = (ar*bstrkw*fact4**4*tke*mp2m1* &
                gam1*dtau*v)/a2/term1

              ipkpi(i,j,k,4,1,n) = -fact1*fact4**3*ar*q2*w*mp2m1*gam1* &
                dtau*(1+ar*dtaut*fact2)*(1+ar*fact3*dtaut)/2/a2/dt/term1

              ipkpi(i,j,k,4,2,n) = fact1*fact4**3*ar*u*w*mp2m1*gam1*dtau* &
                (1+ar*dtaut*fact2)*(1+ar*fact3*dtaut)/(a2*dt)/term1

              ipkpi(i,j,k,4,3,n) = fact1*fact4**3*ar*v*w*mp2m1*gam1*dtau* &
                (1+ar*dtaut*fact2)*(1+ar*fact3*dtaut)/(a2*dt)/term1

              ipkpi(i,j,k,4,4,n) = fact4**3*((ar*fact3*dtaut+1)* &
                (1/(2*a2*dt)* &
                (1+ar*fact2*dtaut)*(2*a2*dt*fact4+fact1*ar*dtau*mp2m1 &
                *gam1*(2*(a2/gam1+tke)+2*w**2))-fact4* &
                (tke*mp2m1*gam1/a2)*ar*fact2*dtaut) &
                -(ar*bstrkw*fact4*mp2m1*gam1*dtaut*tke*ome)/a2)/term1

              ipkpi(i,j,k,4,5,n) = -fact1*fact4**3*ar*w*mp2m1*gam1*dtau* &
                (1+ar*dtaut*fact2)*(1+ar*fact3*dtaut)/(a2*dt)/term1

              ipkpi(i,j,k,4,6,n) = (ar*fact2*fact4**4*mp2m1*gam1* &
                dtau*(1+ar*fact3*dtaut)*w)/a2/term1

              ipkpi(i,j,k,4,7,n) = (ar*bstrkw*fact4**4*tke*mp2m1* &
                gam1*dtau*w)/a2/term1

              ipkpi(i,j,k,5,1,n) = &
                -fact1*fact4**3*ar*q2*(a2/gam1+q2/2+tke)*mp2m1*gam1* &
                dtau*(1+ar*dtaut*fact2)*(1+ar*fact3*dtaut)/2/a2/dt/term1

              ipkpi(i,j,k,5,2,n) = &
                fact1*fact4**3*ar*u*(a2/gam1+q2/2+tke)*mp2m1*gam1* &
                dtau*(1+ar*dtaut*fact2)*(1+ar*fact3*dtaut)/(a2*dt)/term1

              ipkpi(i,j,k,5,3,n) = &
                fact1*fact4**3*ar*v*(a2/gam1+q2/2+tke)*mp2m1*gam1* &
                 dtau*(1+ar*dtaut*fact2)*(1+ar*fact3*dtaut)/a2/dt/term1

              ipkpi(i,j,k,5,4,n) = &
                fact1*fact4**3*ar*w*(a2/gam1+q2/2+tke)*mp2m1*gam1* &
                 dtau*(1+ar*dtaut*fact2)*(1+ar*fact3*dtaut)/a2/dt/term1

              ipkpi(i,j,k,5,5,n) = fact4**3*((1+ar*dtaut*fact3)*(fact4* &
                (1+ar*dtaut*fact2-ar*dtaut*fact2*(tke*mp2m1* &
                gam1)/a2)-(1+ar*dtaut*fact2)*fact1*ar*mp2m1* &
                gam1*dtau*q2/(2*a2*dt))-(ar*bstrkw*fact4*mp2m1* &
                gam1*dtaut*tke*ome)/a2)/term1

              ipkpi(i,j,k,5,6,n) = (ar*fact2*fact4**4*mp2m1*gam1* &
                dtau*(1+ar*fact3*dtaut)*(a2/gam1+q2/2+tke))/a2/term1

              ipkpi(i,j,k,5,7,n) = (ar*bstrkw*fact4**4*tke*mp2m1* &
                gam1*dtau*(a2/gam1+q2/2+tke))/a2/term1

              ipkpi(i,j,k,6,1,n) = -fact1/(2*a2*dt)*ar*fact4**4*tke &
                *mp2m1*q2*gam1*dtaut*(1+ar*fact3*dtaut &
                -ar*bstrkw*dtaut*ome)/term1

              ipkpi(i,j,k,6,2,n) = fact1/(a2*dt)*ar*fact4**4*tke*mp2m1 &
                *u*gam1*dtaut*(1+ar*fact3*dtaut &
                -ar*bstrkw*dtaut*ome)/term1

              ipkpi(i,j,k,6,3,n) = fact1/(a2*dt)*ar*fact4**4*tke*mp2m1 &
                *v*gam1*dtaut*(1+ar*fact3*dtaut &
                -ar*bstrkw*dtaut*ome)/term1

              ipkpi(i,j,k,6,4,n) = fact1/(a2*dt)*ar*fact4**4*tke*mp2m1 &
                *w*gam1*dtaut*(1+ar*fact3*dtaut &
                -ar*bstrkw*dtaut*ome)/term1

              ipkpi(i,j,k,6,5,n) = -fact1/(a2*dt)*ar*fact4**4*tke*mp2m1* &
               gam1*dtaut*(1+ar*fact3*dtaut-ar*bstrkw*dtaut*ome)/term1

              ipkpi(i,j,k,6,6,n) = fact4**4*(fact1*ar*(a2/gam1+tke)* &
                mp2m1*gam1*dtau/(a2*dt)*(1+ar*dtaut*fact3)+fact4* &
                (1+ar*dtaut*fact3-(ar*bstrkw*tke*mp2m1*gam1* &
                dtaut*ome)/a2))/term1

              ipkpi(i,j,k,6,7,n) = -1/(a2*dt)*ar*fact4**4*tke*dtaut* &
                (a2*bstrkw*dt*(1-tke*mp2m1*gam1/a2)+fact1*a2*ar &
                *bstrkw*dtau*mp2)/term1

              ipkpi(i,j,k,7,1,n) = -fact1/(2*a2*dt)*ome*ar*fact4**4* &
                mp2m1*q2*gam1*dtaut*(1+ar*fact2*dtaut)/term1

              ipkpi(i,j,k,7,2,n) = fact1/(a2*dt)*ome*ar*fact4**4*mp2m1* &
                u*gam1*dtaut*(1+ar*fact2*dtaut)/term1

              ipkpi(i,j,k,7,3,n) = fact1/(a2*dt)*ome*ar*fact4**4*mp2m1* &
                v*gam1*dtaut*(1+ar*fact2*dtaut)/term1

              ipkpi(i,j,k,7,4,n) = fact1/(a2*dt)*ome*ar*fact4**4*mp2m1* &
                w*gam1*dtaut*(1+ar*fact2*dtaut)/term1

              ipkpi(i,j,k,7,5,n) = -fact1/(a2*dt)*ome*ar*fact4**4*mp2m1 &
                *gam1*dtaut*(1+ar*fact2*dtaut)/term1

              ipkpi(i,j,k,7,6,n) = 1/a2*ome*ar*fact4**5*fact2*mp2m1 &
                *gam1*dtaut/term1 

              ipkpi(i,j,k,7,7,n) = fact4**4*(fact1*ar*(a2/gam1+tke)* &
                mp2m1 &
                *gam1*dtau/(a2*dt)*(1+ar*dtaut*fact2)+fact4* &
                (1+ar*dtaut*fact2-ar*dtaut*fact2*(tke*mp2m1* &
                gam1)/a2))/term1

!------------ lr is the last row of 
!             ipkp [=(I+ar*dtau*Pc*(1.5/(Delta t)I + A))]
              lr(i,j,k,1,n) = fact1*ar*dtaut*ome*q2*mp2m1*gam1/(2*a2*dt)
              lr(i,j,k,2,n) = -fact1*ar*dtaut*u*ome*mp2m1*gam1/(a2*dt)
              lr(i,j,k,3,n) = -fact1*ar*dtaut*v*ome*mp2m1*gam1/(a2*dt)
              lr(i,j,k,4,n) = -fact1*ar*dtaut*w*ome*mp2m1*gam1/(a2*dt)
              lr(i,j,k,5,n) = fact1*ar*dtaut*ome*mp2m1*gam1/(a2*dt)
              lr(i,j,k,6,n) = -ar*dtaut*ome*mp2m1*gam1*fact2/a2
              lr(i,j,k,7,n) = &
                1+ar*dtaut*(fact3-bstrkw*tke*ome*mp2m1*gam1/a2)

              end do
            end do
          end do
        end do  

      else if (rkimp.and.(.not.(kom.or.kom_bsl.or.kom_sst))) then

!------ INVISCID/LAMINAR PIRK TD CASE

        fact1 = 1.5d0
        if ((itime.eq.1).and.(.not.dual_prop_start)) then
          fact1 = 1.0d0
        end if

!       ipkpi: INV (I + fact2 * P) and stands for
!              (Identity Plus K Preconditioner) Inverse

        do n = 0,2*nharms
          nh = n*hbmove

          do k = 1,km
            do j = 1,jm
              do i = 1,im

              rho  = q(i,j,k,1,n)
              u    = q(i,j,k,2,n)/rho
              v    = q(i,j,k,3,n)/rho
              w    = q(i,j,k,4,n)/rho
              q2   = u**2+v**2+w**2
              e0   = q(i,j,k,5,n)/rho
              p    = gam1*(e0-q2/2)*rho
              a2   = gamma*p/rho
              a    = sqrt(a2)

              if (moving.and.(lsp_vfr.eq.1)) then
                qp2 = (u-xdot(i,j,k,nh))**2 + (v-ydot(i,j,k,nh))**2 + &
                      (w-zdot(i,j,k,nh))**2 
              else
                qp2 = q2
              end if

              mp2  = qp2/a2
              ump2 = (uprv / a)**2
              mp2   = dmin1(1.d0, &
                            dmax1(mp2,ump2,cutoff_pgr(i,j,k,n), &
                                  cutoff_vis(i,j,k,n),epsp(nl)**2))
              mp2m1 = mp2-1
 
              fact2 = fact1 * ar * dtvol(i,j,k,n) * vol(i,j,k) / dt
              term1 = a2*(1+fact2)*(1+fact2*mp2)

              ipkpi(i,j,k,1,1,n) =  1/(1+fact2) - &
                                 fact2*gam1*mp2m1*q2/2/term1
              ipkpi(i,j,k,1,2,n) = fact2*gam1*mp2m1*u/term1
              ipkpi(i,j,k,1,3,n) = fact2*gam1*mp2m1*v/term1
              ipkpi(i,j,k,1,4,n) = fact2*gam1*mp2m1*w/term1
              ipkpi(i,j,k,1,5,n) = -fact2*gam1*mp2m1/term1

              ipkpi(i,j,k,2,1,n) = -fact2*gam1*mp2m1*u*q2/2/term1
              ipkpi(i,j,k,2,2,n) = 1/(1+fact2) + &
                                 fact2*gam1*mp2m1*u**2 /term1
              ipkpi(i,j,k,2,3,n) = fact2*gam1*mp2m1*u*v/term1
              ipkpi(i,j,k,2,4,n) = fact2*gam1*mp2m1*u*w/term1
              ipkpi(i,j,k,2,5,n) = -fact2*gam1*mp2m1*u/term1

              ipkpi(i,j,k,3,1,n) = -fact2*gam1*mp2m1*v*q2/2/term1
              ipkpi(i,j,k,3,2,n) = fact2*gam1*mp2m1*u*v/term1
              ipkpi(i,j,k,3,3,n) = 1/(1+fact2) + &
                                 fact2*gam1*mp2m1*v**2 /term1
              ipkpi(i,j,k,3,4,n) = fact2*gam1*mp2m1*v*w/term1
              ipkpi(i,j,k,3,5,n) =  -fact2*gam1*mp2m1*v/term1

              ipkpi(i,j,k,4,1,n) = -fact2*gam1*mp2m1*w*q2/2/term1
              ipkpi(i,j,k,4,2,n) = fact2*gam1*mp2m1*u*w/term1
              ipkpi(i,j,k,4,3,n) = fact2*gam1*mp2m1*v*w/term1
              ipkpi(i,j,k,4,4,n) = 1/(1+fact2) + &
                                 fact2*gam1*mp2m1*w**2 /term1
              ipkpi(i,j,k,4,5,n) =  -fact2*gam1*mp2m1*w/term1

              ipkpi(i,j,k,5,1,n) = -fact2*mp2m1*q2/2* &
                                 (a2+gam1*q2/2)/term1
              ipkpi(i,j,k,5,2,n) = fact2*mp2m1*u* &
                                 (a2+gam1*q2/2)/term1
              ipkpi(i,j,k,5,3,n) = fact2*mp2m1*v* &
                                 (a2+gam1*q2/2)/term1
              ipkpi(i,j,k,5,4,n) = fact2*mp2m1*w* &
                                 (a2+gam1*q2/2)/term1
              ipkpi(i,j,k,5,5,n) = 1/(1+fact2*mp2) - &
                                 fact2*gam1*mp2m1*q2/2/term1

              end do
            end do
          end do
        end do  

      end if

      if (debug) then

        itest = 1

        if (itest.eq.1) then
!------ test 1: check product of prec and its inverse
        write(*,*) 'from lsp: itest = ', itest

        do n = 0,2*nharms
          do k = 1,km
            do j = 1,jm
              do i = 1,im
                do ipde = 1,npde
                  do ipde1 = 1,npde
                    prod(ipde,ipde1)= 0.d0
                    do ipde2 = 1,npde
                      prod(ipde,ipde1) = prod(ipde,ipde1) + &
                      prec(i,j,k,ipde,ipde2,n)*preci(i,j,k,ipde2,ipde1,n)
                    end do
                    if ((ipde.ne.ipde1).and. &
                        (abs(prod(ipde,ipde1)).gt.1.d-9)) then
                      write(*,*) 'lsp warning: nonzero prod at: ', i,j,k
                      write(*,*) 'prod = ', prod(ipde,ipde1), ipde, ipde1
!dbg                  do ipde2 = 1,npde
!dbg                    write(*,*) (prec(i,j,ipde2,ipde3,n),ipde3=1,4)
!dbg                  enddo
!dbg                  write(*,*)
!dbg                  do ipde2 = 1,npde
!dbg                    write(*,*) (preci(i,j,ipde2,ipde3,n),ipde3=1,4)
!dbg                  enddo
                    end if
                    if ((ipde.eq.ipde1).and. &
                        (abs(prod(ipde,ipde1)-1.).gt.1.d-14)) then
                      write(*,*) 'lsp warning: non-unitary prod at: ', &
                                 i,j,k
                      write(*,*) 'prod = ', prod(ipde,ipde1) , ipde
                    end if
                  end do
                end do
              end do
            end do  
          end do  
        end do

!------ test 2: check product of ipkpi and its inverse
        else if ((itest.eq.2).and. &
                 ((.not.rkimp).and.(kom.or.kom_bsl.or.kom_sst))) then

        write(*,*) 'from non-rkim turb lsp: itest = ', itest

        do n = 0,2*nharms
          do k = 1,km
            do j = 1,jm
              do i = 1,im
                
                akom(1,1) = delpls(i,j,k,n)+bstrkw*q(i,j,k,7,n)/q(i,j,k,1,n)
                akom(1,2) = bstrkw*q(i,j,k,6,n)/q(i,j,k,1,n)
                akom(2,1) = 0
                akom(2,2) = ttrm(i,j,k,2,n)*delpls(i,j,k,n) + &
                            2*ttrm(i,j,k,3,n)*q(i,j,k,7,n)/q(i,j,k,1,n)

                do ipde = 1,npde
                  tmp(ipde,1) = 0.d0
                  tmp(ipde,2) = 0.d0
                  tmp(ipde,3) = 0.d0
                  tmp(ipde,4) = 0.d0
                  tmp(ipde,5) = 0.d0
                  tmp(ipde,6) = prec(i,j,k,ipde,6,n) * akom(1,1)
                  tmp(ipde,7) = prec(i,j,k,ipde,6,n) * akom(1,2) + &
                                prec(i,j,k,ipde,7,n) * akom(2,2)
                end do

                do ipde = 1,5
                  do ipde1 = 1,npde
                    fact1 = ar * dtvol(i,j,k,n) * vol(i,j,k) 
                    prod(ipde,ipde1) = 0.d0
                    do ipde2 = 1,npde
                      if (ipde.eq.ipde2) then
                        entry = 1.d0
                      else
                        entry = 0.d0
                      end if
                      prod(ipde,ipde1) = prod(ipde,ipde1) + &
                        (entry + &
                         fact1*tmp(ipde,ipde2))*ipkpi(i,j,k,ipde2,ipde1,n)
                    end do

                    if ((ipde.ne.ipde1).and. &
                        (abs(prod(ipde,ipde1)).gt.1.d-9)) then
                      write(*,*) 'lsp warning: nonzero prod at: ', i,j,k
                      write(*,*) 'prod = ', prod(ipde,ipde1), ipde, ipde1
                    end if
                    if ((ipde.eq.ipde1).and. &
                        (abs(prod(ipde,ipde1)-1.d0).gt.1.d-14)) then
                      write(*,*) 'lsp warning: non-unitary prod at: ', &
                                 i,j,k
                      write(*,*) 'prod = ', prod(ipde,ipde1), ipde, ipde1
                    end if

                  end do
                end do

                do ipde = 6,npde
                  do ipde1 = 1,npde
                    fact1 = ar * dtvolt(i,j,k,n) * vol(i,j,k)
                    prod(ipde,ipde1) = 0.d0
                    do ipde2 = 1,npde
                      if (ipde.eq.ipde2) then
                        entry = 1.d0
                      else
                        entry = 0.d0
                      end if
                      prod(ipde,ipde1) = prod(ipde,ipde1) + &
                        (entry + &
                         fact1*tmp(ipde,ipde2))*ipkpi(i,j,k,ipde2,ipde1,n)
                    end do

                    if ((ipde.ne.ipde1).and. &
                        (abs(prod(ipde,ipde1)).gt.1.d-9)) then
                      write(*,*) 'lsp warning: nonzero prod at: ', i,j,k
                      write(*,*) 'prod = ', prod(ipde,ipde1), ipde, ipde1
                    end if
                    if ((ipde.eq.ipde1).and. &
                        (abs(prod(ipde,ipde1)-1.d0).gt.1.d-14)) then
                      write(*,*) 'lsp warning: non-unitary prod at: ', &
                                 i,j,k
                      write(*,*) 'prod = ', prod(ipde,ipde1),ipde,ipde1
                    end if

                  end do
                end do
              end do
            end do  
          end do 
        end do

        else if ((itest.eq.2).and. &
                 (rkimp.and.(kom.or.kom_bsl.or.kom_sst))) then

        write(*,*) 'from rkim turb lsp: itest = ', itest

        do n = 0,2*nharms
          do k = 1,km
            do j = 1,jm
              do i = 1,im

                akom(1,1) = delpls(i,j,k,n)+bstrkw*q(i,j,k,7,n)/q(i,j,k,1,n)
                akom(1,2) = bstrkw*q(i,j,k,6,n)/q(i,j,k,1,n)
                akom(2,1) = 0
                akom(2,2) = ttrm(i,j,k,2,n)*delpls(i,j,k,n) + &
                            2*ttrm(i,j,k,3,n)*q(i,j,k,7,n)/q(i,j,k,1,n)

                do ipde = 1,npde
                  do ipde1=1,npde
                    tmp(ipde,ipde1) = prec(i,j,k,ipde,ipde1,n)*fact1/dt
                  end do
                  tmp(ipde,6) = tmp(ipde,6) + &
                    prec(i,j,k,ipde,6,n)*akom(1,1)
                  tmp(ipde,7) = tmp(ipde,7) + &
                    prec(i,j,k,ipde,6,n) * akom(1,2) + &
                    prec(i,j,k,ipde,7,n) * akom(2,2)
                end do

                do ipde = 1,5
                  do ipde1 = 1,npde
                    dtau  = dtvol(i,j,k,n) * vol(i,j,k) 
                    prod (ipde,ipde1) = 0.d0
                    do ipde2 = 1,npde
                      if (ipde.eq.ipde2) then
                        entry = 1.d0
                      else
                        entry = 0.d0
                      end if
                      prod(ipde,ipde1) = prod(ipde,ipde1) + &
                        (entry + ar*dtau*tmp(ipde,ipde2))* &
                        ipkpi(i,j,k,ipde2,ipde1,n)
                    end do
                    if ((ipde.ne.ipde1).and. &
                         (abs(prod(ipde,ipde1)).gt.1.d-9)) then
                      write(*,*) 'lsp warning: nonzero prod at: ', i,j,k
                      write(*,*) 'prod = ', prod(ipde,ipde1), ipde, ipde1
                    end if
                    if ((ipde.eq.ipde1).and. &
                        (abs(prod(ipde,ipde1)-1.d0).gt.1.d-14)) then
                      write(*,*) 'lsp warning: non-unitary prod at: ', &
                                 i,j,k
                      write(*,*) 'prod = ', prod(ipde,ipde1),ipde,ipde1
                    end if
                  end do
                end do

                do ipde = 6,npde
                  do ipde1 = 1,npde
                    dtaut  = dtvolt(i,j,k,n) * vol(i,j,k) 
                    prod(ipde,ipde1) = 0.d0
                    do ipde2 = 1,npde
                      if (ipde.eq.ipde2) then
                        entry = 1.d0
                      else
                        entry = 0.d0
                      end if
                      prod(ipde,ipde1) = prod(ipde,ipde1) + &
                        (entry + ar*dtaut*tmp(ipde,ipde2))* &
                        ipkpi(i,j,k,ipde2,ipde1,n)
                    end do
                    if  ((ipde.ne.ipde1).and. &
                         (abs(prod(ipde,ipde1)).gt.1.d-9)) then
                      write(*,*) 'lsp warning: nonzero prod at: ', i,j,k
                      write(*,*) 'prod = ', prod(ipde,ipde1), ipde, ipde1
                    end if
                    if ((ipde.eq.ipde1).and. &
                        (abs(prod(ipde,ipde1)-1.d0).gt.1.d-14)) then
                      write(*,*) 'lsp warning: non-unitary prod at: ', &
                                 i,j,k
                      write(*,*) 'prod = ', prod(ipde,ipde1),ipde,ipde1
                    end if
                  end do
                end do
                ipde = npde
                do ipde1 = 1,npde
                  prod(ipde,ipde1) = 0.d0
                  do ipde2 = 1,npde
                     prod(ipde,ipde1) = prod(ipde,ipde1) + &
                       lr(i,j,k,ipde2,n)*ipkpi(i,j,k,ipde2,ipde1,n)  
                  end do
                  if  ((ipde.ne.ipde1).and. &
                       (abs(prod(ipde,ipde1)).gt.1.d-9)) then
                    write(*,*) 'lsp warning: nonzero lr-prod at: ',i,j,k
                    write(*,*) 'prod = ', prod(ipde,ipde1), ipde,ipde1
                  end if
                  if ((ipde.eq.ipde1).and. &
                      (abs(prod(ipde,ipde1)-1.d0).gt.1.d-14)) then
                    write(*,*) 'lsp warning: non-unitary lr-prod at: ', &
                               i,j,k
                    write(*,*) 'prod = ', prod(ipde,ipde1),ipde,ipde1
                  end if
                end do

              end do
            end do  
          end do  
        end do

        else if (itest.eq.2.and. &
                 (rkimp.and.(.not.(kom.or.kom_bsl.or.kom_sst)))) then

        write(*,*) 'from rkim non-turb lsp: itest = ', itest

        do n = 0,2*nharms
          do k = 1,km
            do j = 1,jm
              do i = 1,im
                do ipde = 1,npde
                  do ipde1 = 1,npde
                    fact2 = fact1 * ar * dtvol(i,j,k,n) * vol(i,j,k) / dt
                    prod(ipde,ipde1) = 0.d0
                    do ipde2 = 1,npde
                      if (ipde.eq.ipde2) then
                        entry = 1.d0
                      else
                        entry = 0.d0
                      end if
                      prod(ipde,ipde1) = prod(ipde,ipde1) + &
                                 (entry + fact2*prec(i,j,k,ipde,ipde2,n)) * &
                                 ipkpi(i,j,k,ipde2,ipde1,n)
                    end do
                    if ((ipde.ne.ipde1).and. &
                        (abs(prod(ipde,ipde1)).gt.1.d-9)) then
                      write(*,*) 'lsp warning: nonzero prod at: ', i,j,k
                      write(*,*) 'prod = ', prod(ipde,ipde1), ipde, ipde1
!dbg                  do ipde2 = 1,4
!dbg                    write(*,*) (prec(i,j,ipde2,ipde3,n),ipde3=1,4)
!dbg                  enddo
!dbg                  write(*,*)
!dbg                  do ipde2 = 1,4
!dbg                    write(*,*) (preci(i,j,ipde2,ipde3,n),ipde3=1,4)
!dbg                  enddo
                    end if
                    if ((ipde.eq.ipde1).and. &
                        (abs(prod(ipde,ipde1)-1.d0).gt.1.d-14)) then
                      write(*,*) 'lsp warning: non-unitary prod at: ', &
                                 i,j,k
                      write(*,*) 'prod = ', prod(ipde,ipde1)
                    end if
                  end do
                end do
              end do
            end do  
          end do  
        end do

        end if

      end if

      return      
      end 

!-----------------------------------------------------------------------
      subroutine lm_cutoff_pg(nl,q,cutoff_pgr)
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,iblock,imax,jmax,kmax,iq,icof
      real(kind=cosa_real) q(*),cutoff_pgr(*)
      real(kind=cosa_real) rho,t,u,v,w,q2,a2,a,mu,dp(4),dpmax,dx,dy,dz,ds,ccutoff
 
      do iblock = 1,mynblocks
        imax   = i_imax     (iblock,nl)
        jmax   = j_jmax     (iblock,nl)
        kmax   = k_kmax     (iblock,nl)
        iq     = 1 + off_p3 (iblock,nl) * npde * dim5
        icof   = 1 + off_m1 (iblock,nl) *        dim5
        call lm_bcutoff_pg(q(iq),cutoff_pgr(icof),imax,jmax,kmax, &
                             lmet,npde,nharms)
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine lm_bcutoff_pg(q,cutoff_pgr,imax,jmax,kmax,lmet,npde,nharms)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,lmet,npde,nharms
      integer(kind=cosa_int) im,i,jm,j,km,k,iface,n,nh
      real(kind=cosa_real) &
          q         (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
          cutoff_pgr(   imax-1,   jmax-1,   kmax-1,     0:2*nharms)
      real(kind=cosa_real) rho,t,u,v,w,q2,a2,a,p,dp(6),dpmax,ccutoff
 
      im = imax-1
      jm = jmax-1
      km = kmax-1

!---- pressure-gradient-dependent cutoff (Weiss and Darmofal)
      
        do n = 0,2*nharms
          do k = 1,km
            do j = 1,jm
              do i = 1,im

              rho = q(i,j,k,1,n)
              u   = q(i,j,k,2,n)
              v   = q(i,j,k,3,n)
              w   = q(i,j,k,4,n)
              p   = q(i,j,k,5,n)
              q2  = u**2 + v**2 + w**2
              a2  = gamma*p/rho
              a = sqrt(a2)
              dp(1) = dabs(q(i+1,j  ,k,  5,n) - q(i,j,k,5,n))
              dp(2) = dabs(q(i  ,j+1,k,  5,n) - q(i,j,k,5,n))
              dp(3) = dabs(q(i  ,j  ,k+1,5,n) - q(i,j,k,5,n))
              dp(4) = dabs(q(i-1,j  ,k,  5,n) - q(i,j,k,5,n))
              dp(5) = dabs(q(i  ,j-1,k,  5,n) - q(i,j,k,5,n))
              dp(6) = dabs(q(i  ,j  ,k-1,5,n) - q(i,j,k,5,n))
              dpmax = 0.d0
              do iface=1,6
                dpmax = dmax1(dpmax,dp(iface))
              end do
              cutoff_pgr(i,j,k,n) = dpmax/rho/a2

              end do
            end do
          end do
        end do

      return
      end

!-----------------------------------------------------------------------
      subroutine lm_cutoff_bi(coff_lr,cutoff_pgr,cutoff_vis,imax,jmax, &
                             kmax,npde,nharms)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,npde,nharms
      integer(kind=cosa_int) im,i,jm,j,km,k,iface,n
      real (kind=cosa_real) &
          coff_lr   (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
          cutoff_pgr(   imax-1,   jmax-1,   kmax-1,     0:2*nharms), &
          cutoff_vis(   imax-1,   jmax-1,   kmax-1,     0:2*nharms)
 
      im = imax-1
      jm = jmax-1
      km = kmax-1

      do n=0,2*nharms
        do k = 1,km
          do j = 1,jm
            coff_lr(1,j,k,1,n) = cutoff_pgr(1,j,k,n)
            coff_lr(1,j,k,2,n) = cutoff_vis(1,j,k,n)
            do i = 2,im
              coff_lr(i,j,k,1,n) = &
                           dmax1(cutoff_pgr(i,j,k,n),cutoff_pgr(i-1,j,k,n))
              coff_lr(i,j,k,2,n) = &
                           dmax1(cutoff_vis(i,j,k,n),cutoff_vis(i-1,j,k,n))
            end do
            coff_lr(imax,j,k,1,n) = cutoff_pgr(im,j,k,n)
            coff_lr(imax,j,k,2,n) = cutoff_vis(im,j,k,n)
          end do
        end do
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine lm_cutoff_bj(coff_lr,cutoff_pgr,cutoff_vis,imax,jmax, &
                             kmax,npde,nharms)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,npde,nharms
      integer(kind=cosa_int) im,i,jm,j,km,k,iface,n
      real (kind=cosa_real) &
          coff_lr   (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
          cutoff_pgr(   imax-1,   jmax-1,   kmax-1,     0:2*nharms), &
          cutoff_vis(   imax-1,   jmax-1,   kmax-1,     0:2*nharms)
 
      im = imax-1
      jm = jmax-1
      km = kmax-1

      do n=0,2*nharms
        do k = 1,km
          do i = 1,im
            coff_lr(i,1,k,1,n) = cutoff_pgr(i,1,k,n)
            coff_lr(i,1,k,2,n) = cutoff_vis(i,1,k,n)
            do j = 2,jm
              coff_lr(i,j,k,1,n) = &
                           dmax1(cutoff_pgr(i,j,k,n),cutoff_pgr(i,j-1,k,n))
              coff_lr(i,j,k,2,n) = &
                           dmax1(cutoff_vis(i,j,k,n),cutoff_vis(i,j-1,k,n))
            end do
            coff_lr(i,jmax,k,1,n) = cutoff_pgr(i,jm,k,n)
            coff_lr(i,jmax,k,2,n) = cutoff_vis(i,jm,k,n)
          end do
        end do
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine lm_cutoff_bk(coff_lr,cutoff_pgr,cutoff_vis,imax,jmax, &
                             kmax,npde,nharms)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,npde,nharms
      integer(kind=cosa_int) im,i,jm,j,km,k,iface,n
      real (kind=cosa_real) &
          coff_lr   (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
          cutoff_pgr(   imax-1,   jmax-1,   kmax-1,     0:2*nharms), &
          cutoff_vis(   imax-1,   jmax-1,   kmax-1,     0:2*nharms)
 
      im = imax-1
      jm = jmax-1
      km = kmax-1

      do n=0,2*nharms
        do j = 1,jm
          do i = 1,im
            coff_lr(i,j,1,1,n) = cutoff_pgr(i,j,1,n)
            coff_lr(i,j,1,2,n) = cutoff_vis(i,j,1,n)
            do k = 2,km
              coff_lr(i,j,k,1,n) = &
                           dmax1(cutoff_pgr(i,j,k,n),cutoff_pgr(i,j,k-1,n))
              coff_lr(i,j,k,2,n) = &
                           dmax1(cutoff_vis(i,j,k,n),cutoff_vis(i,j,k-1,n))
            end do
            coff_lr(i,j,kmax,1,n) = cutoff_pgr(i,j,km,n)
            coff_lr(i,j,kmax,2,n) = cutoff_vis(i,j,km,n)
          end do
        end do
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine lm_viscoff(nl,q,si,sj,sk,cutoff_vis)
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,iblock,imax,jmax,kmax,iq,icof,iimt
      real(kind=cosa_real) q(*),cutoff_vis(*),si(*),sj(*),sk(*)
 
      do iblock = 1,mynblocks
        imax   = i_imax     (iblock,nl)
        jmax   = j_jmax     (iblock,nl)
        kmax   = k_kmax     (iblock,nl)
        iq     = 1 + off_p3 (iblock,nl) * npde * dim5
        icof   = 1 + off_m1 (iblock,nl) *        dim5
        iimt   = 1 + off_p2 (iblock,nl) * lmet * dim5h
        call lm_bviscoff(q(iq),si(iimt),sj(iimt),sk(iimt), &
                         cutoff_vis(icof),imax,jmax,kmax,lmet,npde,nharms)
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine lm_bviscoff(q,si,sj,sk,cutoff_vis,imax,jmax,kmax, &
                             lmet,npde,nharms)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,lmet,npde,nharms
      integer(kind=cosa_int) im,i,jm,j,km,k,iface,n,nh,imet
      real(kind=cosa_real) &
          q         (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
          cutoff_vis(   imax-1,   jmax-1,   kmax-1,     0:2*nharms), &
          si        (lmet,0:imax+1,0:jmax+1,0:kmax+1,0:2*nharms*hbmove), &
          sj        (lmet,0:imax+1,0:jmax+1,0:kmax+1,0:2*nharms*hbmove), &
          sk        (lmet,0:imax+1,0:jmax+1,0:kmax+1,0:2*nharms*hbmove)
      real(kind=cosa_real) rho,t,u,v,w,a,p,mu,dx,dy,dz,ds, &
          ccutoff,cutoffx,cutoffy,cutoffz,rex,rey,rez,machx,machy,machz, &
          xix,xiy,xiz,etax,etay,etaz,zetax,zetay,zetaz,ux,uy,uz, &
          sim(4),sip(4),sjm(4),sjp(4),skm(4),skp(4),epsi
 
      im = imax-1
      jm = jmax-1
      km = kmax-1

      epsi = 1.d-7

      do n = 0,2*nharms
        nh = n*hbmove
        do k = 1,km
          do j = 1,jm
            do i = 1,im

            rho = q(i,j,k,1,n)
            u   = q(i,j,k,2,n)
            v   = q(i,j,k,3,n)
            w   = q(i,j,k,4,n)
            p   = q(i,j,k,5,n)
            t   = gamma*p/rho
            a   = sqrt(t)
            mu  = (stemp+1)/(t+stemp) * (t**1.5d0)

            do imet = 1,4
              sim(imet) = si(imet,i  ,j  ,k  ,nh)
              sip(imet) = si(imet,i+1,j  ,k  ,nh)
              sjm(imet) = sj(imet,i,  j  ,k  ,nh)
              sjp(imet) = sj(imet,i,  j+1,k  ,nh)
              skm(imet) = sk(imet,i,  j  ,k  ,nh)
              skp(imet) = sk(imet,i,  j  ,k+1,nh)
            end do

            xix   = 0.5d0 * (sim(1) + sip(1))
            xiy   = 0.5d0 * (sim(2) + sip(2))
            xiz   = 0.5d0 * (sim(3) + sip(3))

            etax  = 0.5d0 * (sjm(1) + sjp(1))
            etay  = 0.5d0 * (sjm(2) + sjp(2))
            etaz  = 0.5d0 * (sjm(3) + sjp(3))

            zetax = 0.5d0 * (skm(1) + skp(1))
            zetay = 0.5d0 * (skm(2) + skp(2))
            zetaz = 0.5d0 * (skm(3) + skp(3))

            ux = max(epsi,abs( u*xix    + v*xiy    + w*xiz  ))
            uy = max(epsi,abs( u*etax   + v*etay   + w*etaz ))
            uz = max(epsi,abs( u*zetax  + v*zetay  + w*zetaz))

            dy = (sim(4)+sip(4)) / 2
            dx = (sjm(4)+sjp(4)) / 2
            dz = (skm(4)+skp(4)) / 2

            rex = reyno/machfs * rho * ux * dx / mu
            machx = ux / a
            cutoffx = (machx**2 * (1-rex)) / &
                      (rex**2 + machx**2*rex*(1-rex))

            rey = reyno/machfs * rho * uy * dy / mu
            machy = uy / a
            cutoffy = (machy**2 * (1-rey)) / &
                      (rey**2 + machy**2*rey*(1-rey))

            rez = reyno/machfs * rho * uz * dz / mu
            machz = uz / a
            cutoffz = (machz**2 * (1-rez)) / &
                      (rez**2 + machz**2*rez*(1-rez))

            ccutoff = dmax1(cutoffx,cutoffy,cutoffz)
            cutoff_vis(i,j,k,n) = ccutoff

            end do
          end do
        end do
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine lm_cutoff(idir,nl,coff_lr,cutoff_pgr,cutoff_vis)
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) idir,nl,iblock,imax,jmax,kmax,iwrk,icof
      real (kind=cosa_real) coff_lr(*),cutoff_pgr(*),cutoff_vis(*)
 
      do iblock = 1,mynblocks
        imax   = i_imax     (iblock,nl)
        jmax   = j_jmax     (iblock,nl)
        kmax   = k_kmax     (iblock,nl)
        iwrk   = 1 + off_p3 (iblock,nl) * npde * dim5
        icof   = 1 + off_m1 (iblock,nl) *        dim5
        if (idir.eq.1) then
          call lm_cutoff_bi(coff_lr(iwrk),cutoff_pgr(icof), &
                            cutoff_vis(icof),imax,jmax,kmax,npde,nharms)
        else if (idir.eq.2) then
          call lm_cutoff_bj(coff_lr(iwrk),cutoff_pgr(icof), &
                            cutoff_vis(icof),imax,jmax,kmax,npde,nharms)
        else if (idir.eq.3) then
          call lm_cutoff_bk(coff_lr(iwrk),cutoff_pgr(icof), &
                            cutoff_vis(icof),imax,jmax,kmax,npde,nharms)
        end if
      end do

      return
      end
