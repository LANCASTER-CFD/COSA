!-----------------------------------------------------------------------
      subroutine update(nl,res,q,qold,vol,dtvol,dtvolt,prec,ipkpi,work1, &
                        qhb,lr,delpls,prd,ttrm,ar)
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none
!
      integer(kind=cosa_int) nl,iblock,imax,jmax,kmax,iq,idtv,ivol,iprec,idlpl, &
                ittrm
      real (kind=cosa_real) ar
      real (kind=cosa_real) q(*),qold(*),res(*),vol(*),dtvol(*),dtvolt(*), &
           prec(*),ipkpi(*),work1(*),qhb(*),delpls(*),prd(*),lr(*), &
           ttrm(*)
!
      if (harbal) then

        if (lomach) then

          do iblock = 1,mynblocks
            imax   = i_imax     (iblock,nl)
            jmax   = j_jmax     (iblock,nl)
            kmax   = k_kmax     (iblock,nl)
            iq     = 1 + off_p3 (iblock,nl) * npde * dim5
            idlpl  = 1 + off_p3 (iblock,nl) *        dim5
            idtv   = 1 + off_p1 (iblock,nl)        * dim5
            ivol   = 1 + off_p1 (iblock,nl)
            iprec  = 1 + off_m1 (iblock,nl) * npde * dim5 * npde
            ittrm  = 1 + off_p3 (iblock,nl) * 3    * dim5
            call update_lsp_bhb(res(iq),q(iq),qold(iq),vol(ivol), &
                   dtvol(idtv),dtvolt(idtv),prec(iprec),ipkpi(iprec), &
                   qhb(iq),lr(iq),delpls(idlpl),prd(idlpl), &
                   ttrm(ittrm),ar,imax,jmax,kmax,npde,nharms)
          end do

        else

          do iblock = 1,mynblocks
            imax   = i_imax     (iblock,nl)
            jmax   = j_jmax     (iblock,nl)
            kmax   = k_kmax     (iblock,nl)
            iq     = 1 + off_p3 (iblock,nl) * npde * dim5
            idlpl  = 1 + off_p3 (iblock,nl) *        dim5
            idtv   = 1 + off_p1 (iblock,nl)        * dim5
            ivol   = 1 + off_p1 (iblock,nl)
            ittrm  = 1 + off_p3 (iblock,nl) * 3    * dim5
            call update_bhb(res(iq),q(iq),qold(iq),vol(ivol), &
                   dtvol(idtv),dtvolt(idtv), &
                   work1(iq),qhb(iq),delpls(idlpl),prd(idlpl), &
                   ttrm(ittrm),ar,imax,jmax,kmax,npde,nharms)
          end do

        end if

      else

        if (lomach) then

          do iblock = 1,mynblocks
            imax   = i_imax     (iblock,nl)
            jmax   = j_jmax     (iblock,nl)
            kmax   = k_kmax     (iblock,nl)
            iq     = 1 + off_p3 (iblock,nl) * npde * dim5
            idlpl  = 1 + off_p3 (iblock,nl) *        dim5
            idtv   = 1 + off_p1 (iblock,nl)        * dim5
            ivol   = 1 + off_p1 (iblock,nl)
            iprec  = 1 + off_m1 (iblock,nl) * npde * dim5 * npde
            ittrm  = 1 + off_p3 (iblock,nl) * 3    * dim5
            call update_lsp_b(res(iq),q(iq),qold(iq),vol(ivol), &
                   dtvol(idtv),dtvolt(idtv),prec(iprec),ipkpi(iprec), &
                   work1(iq),lr(iq),delpls(idlpl),prd(idlpl), &
                   ttrm(ittrm),ar,imax,jmax,kmax,npde,nharms)
          end do

        else

          do iblock = 1,mynblocks
            imax   = i_imax     (iblock,nl)
            jmax   = j_jmax     (iblock,nl)
            kmax   = k_kmax     (iblock,nl)
            iq     = 1 + off_p3 (iblock,nl) * npde * dim5
            idlpl  = 1 + off_p3 (iblock,nl) *        dim5
            idtv   = 1 + off_p1 (iblock,nl)        * dim5
            ivol   = 1 + off_p1 (iblock,nl)
            iprec  = 1 + off_m1 (iblock,nl) * npde * dim5 * npde
            ittrm  = 1 + off_p3 (iblock,nl) * 3    * dim5
            call update_b(res(iq),q(iq),qold(iq),vol(ivol),dtvol(idtv), &
                   dtvolt(idtv),delpls(idlpl),prd(idlpl), &
                   ttrm(ittrm),ar,imax,jmax,kmax,npde,nharms)
          end do

        end if

      end if

      return
      end

!-----------------------------------------------------------------------
      subroutine update_lsp_b(res,q,qold,vol,dtvol,dtvolt,prec,ipkpi, &
                   work1,lr,delpls,prd,ttrm,ar,imax,jmax,kmax,npde, &
                   nharms)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none
!
      integer(kind=cosa_int) imax,jmax,kmax,npde,nharms
      integer(kind=cosa_int) im,jm,km,i,j,k,ipde,ipde1
      real (kind=cosa_real) ar, fact, fact1, term1, dtau, rho, ome, tke, &
           rhomega,drhomega,rhomegamin,a(2,2),tmp(7,7),somma,bav,gav
      real (kind=cosa_real) &
           q     (-1:imax+1,-1:jmax+1,-1:kmax+1,     npde), &
           qold  (-1:imax+1,-1:jmax+1,-1:kmax+1,     npde), &
           res   (-1:imax+1,-1:jmax+1,-1:kmax+1,     npde), &
           work1 (-1:imax+1,-1:jmax+1,-1:kmax+1,     npde), &
           lr    (-1:imax+1,-1:jmax+1,-1:kmax+1,     npde), &
           ttrm  (-1:imax+1,-1:jmax+1,-1:kmax+1,        3), &
           delpls(-1:imax+1,-1:jmax+1,-1:kmax+1          ), &
           prd   (-1:imax+1,-1:jmax+1,-1:kmax+1          ), &
           dtvol ( 0:imax  , 0:jmax  , 0:kmax            ), &
           dtvolt( 0:imax  , 0:jmax  , 0:kmax            ), &
           prec  (   imax-1,   jmax-1,   kmax-1,npde,npde), &
           ipkpi (   imax-1,   jmax-1,   kmax-1,npde,npde), &
           vol   ( 0:imax  , 0:jmax  , 0:kmax            )
!
      im = imax-1
      jm = jmax-1
      km = kmax-1

      do ipde = 1,npde
        do ipde1 = 1,npde
          tmp(ipde,ipde1) = 0.d0
        end do
      end do

      if (.not.rkimp) then

!-----------------------------------------------------------------------
!     STEADY OR TIME-DEPENDENT WITH FERK
!-----------------------------------------------------------------------

        if (kom.or.kom_bsl.or.kom_sst) then

          do k = 1,km
            do j = 1,jm
              do i = 1,im

                gav        = ttrm(i,j,k,2)
                bav        = ttrm(i,j,k,3)
          
                a(1,1) = delpls(i,j,k) + &
                         bstrkw * q(i,j,k,7) / q(i,j,k,1)
                a(1,2) = bstrkw * q(i,j,k,6) / q(i,j,k,1)
                a(2,1) = 0.d0
                a(2,2) = gav * delpls(i,j,k) + &
                         2 * bav * q(i,j,k,7) / q(i,j,k,1)

                do ipde = 1,npde
                  tmp(ipde,6) = prec(i,j,k,ipde,6) * a(1,1)
                  tmp(ipde,7) = prec(i,j,k,ipde,6) * a(1,2) + &
                                prec(i,j,k,ipde,7) * a(2,2)
                  work1(i,j,k,ipde) = &
                    tmp(ipde,6) * q(i,j,k,6) + tmp(ipde,7) * q(i,j,k,7)
                end do

              end do
            end do  
          end do  

!-------- update OMEGA and its residual
          ipde=7
          do k = 1,km
            do j = 1,jm
              do i = 1,im

                gav        = ttrm(i,j,k,2)

                term1 = ar * dtvolt(i,j,k)
                fact1 = term1 * vol(i,j,k)
                q(i,j,k,ipde) = 0
                do ipde1 = 1,npde
                  q(i,j,k,ipde) = q(i,j,k,ipde) + &
                    ipkpi(i,j,k,ipde,ipde1) * &
                    (qold(i,j,k,ipde1) + fact1 * work1(i,j,k,ipde1) - &
                                   term1 * res(i,j,k,ipde1))
                end do
                rhomegamin = gav*cmu*q(i,j,k,1)*sqrt(prd(i,j,k))
                if (q(i,j,k,ipde).lt.rhomegamin) then
                  q(i,j,k,ipde) = rhomegamin
                  somma = 0.d0
                  do ipde1 = 1,npde
                    somma = somma + lr(i,j,k,ipde1) * q(i,j,k,ipde1)
                  end do
                  res(i,j,k,ipde) = &
                    (-somma +qold(i,j,k,ipde) + &
                      fact1*work1(i,j,k,ipde) ) / term1
                end if

              end do
            end do
          end do

!-------- update TKE
          ipde=6
          do k = 1,km
            do j = 1,jm
              do i = 1,im
                term1 = ar * dtvolt(i,j,k)
                fact1 = term1 * vol(i,j,k)
                q(i,j,k,ipde) = 0.d0
                do ipde1 = 1,npde
                  q(i,j,k,ipde) = q(i,j,k,ipde) + &
                    ipkpi(i,j,k,ipde,ipde1) * &
                    (qold(i,j,k,ipde1) + fact1 * work1(i,j,k,ipde1) - &
                                   term1 * res(i,j,k,ipde1))
                end do
              end do
            end do
          end do

!-------- update NS variables
          do ipde = 1,5
            do k = 1,km
              do j = 1,jm
                do i = 1,im
                  term1 = ar * dtvol(i,j,k)
                  fact1 = term1 * vol(i,j,k)
                  q(i,j,k,ipde) = 0.d0
                  do ipde1 = 1,npde
                    q(i,j,k,ipde) = q(i,j,k,ipde) + &
                      ipkpi(i,j,k,ipde,ipde1) * &
                      (qold(i,j,k,ipde1) + fact1 * work1(i,j,k,ipde1) - &
                                     term1 * res(i,j,k,ipde1))
                  end do
                end do
              end do
            end do
          end do

        else

          do ipde = 1,npde
            do k = 1,km
              do j = 1,jm
                do i = 1,im
                  q(i,j,k,ipde) = qold(i,j,k,ipde) - &
                                ar*res(i,j,k,ipde)*dtvol(i,j,k)
                end do
              end do
            end do
          end do

        end if

      else if (rkimp) then

        fact = 1.5d0
        if ((itime.eq.1).and.(.not.dual_prop_start)) then
          fact = 1.0d0
        end if

        if ((irsop.eq.0).and.(nlevel.eq.1)) then

!-----------------------------------------------------------------------
!       TIME-DEPENDENT WITH SIRK and WITHOUT MG OR IRS
!-----------------------------------------------------------------------
!       this update section is FORMALLY identical to the FERK update

          if (kom.or.kom_bsl.or.kom_sst) then

            do k = 1,km
              do j = 1,jm
                do i = 1,im

                  gav        = ttrm(i,j,k,2)
                  bav        = ttrm(i,j,k,3)
          
                  a(1,1) = delpls(i,j,k) + &
                           bstrkw * q(i,j,k,7) / q(i,j,k,1)
                  a(1,2) = bstrkw * q(i,j,k,6) / q(i,j,k,1)
                  a(2,1) = 0.d0
                  a(2,2) = gav * delpls(i,j,k) + &
                           2 * bav * q(i,j,k,7) / q(i,j,k,1)

                  do ipde = 1,npde
                    tmp(ipde,6) = prec(i,j,k,ipde,6) * a(1,1)
                    tmp(ipde,7) = prec(i,j,k,ipde,6) * a(1,2) + &
                                  prec(i,j,k,ipde,7) * a(2,2)
                    work1(i,j,k,ipde) = &
                      tmp(ipde,6) * q(i,j,k,6) + &
                                    tmp(ipde,7) * q(i,j,k,7)
                  end do

                end do
              end do  
            end do  

!---------- update OMEGA and its residual
            ipde=7
            do k = 1,km
              do j = 1,jm
                do i = 1,im

                  gav        = ttrm(i,j,k,2)

                  term1 = ar * dtvolt(i,j,k)
                  fact1 = term1 * vol(i,j,k)
                  q(i,j,k,ipde) = 0
                  do ipde1 = 1,npde
                    q(i,j,k,ipde) = q(i,j,k,ipde) + &
                      ipkpi(i,j,k,ipde,ipde1) * &
                      (qold(i,j,k,ipde1) + fact1 * work1(i,j,k,ipde1) - &
                                     term1 * res(i,j,k,ipde1))
                  end do
                  rhomegamin = gav*cmu*q(i,j,k,1)*sqrt(prd(i,j,k))
                  if (q(i,j,k,ipde).lt.rhomegamin) then
                    q(i,j,k,ipde) = rhomegamin
                    somma = 0.d0
                    do ipde1 = 1,npde
                      somma = somma + lr(i,j,k,ipde1) * q(i,j,k,ipde1)
                    end do
                    res(i,j,k,ipde) = &
                      (-somma + qold(i,j,k,ipde) + &
                        fact1*work1(i,j,k,ipde) ) / term1
                  end if

                end do
              end do
            end do

!---------- update TKE
            ipde=6
            do k = 1,km
              do j = 1,jm
                do i = 1,im
                  term1 = ar * dtvolt(i,j,k)
                  fact1 = term1 * vol(i,j,k)
                  q(i,j,k,ipde) = 0.d0
                  do ipde1 = 1,npde
                    q(i,j,k,ipde) = q(i,j,k,ipde) + &
                      ipkpi(i,j,k,ipde,ipde1) * &
                      (qold(i,j,k,ipde1) + fact1 * work1(i,j,k,ipde1) - &
                                   term1 * res(i,j,k,ipde1))
                  end do
                end do
              end do
            end do

!---------- update NS variables
            do ipde = 1,5
              do k = 1,km
                do j = 1,jm
                  do i = 1,im
                    term1 = ar * dtvol(i,j,k)
                    fact1 = term1 * vol(i,j,k)
                    q(i,j,k,ipde) = 0.d0
                    do ipde1 = 1,npde
                      q(i,j,k,ipde) = q(i,j,k,ipde) + &
                        ipkpi(i,j,k,ipde,ipde1) * &
                        (qold(i,j,k,ipde1) + &
                         fact1 * work1(i,j,k,ipde1) - &
                         term1 * res(i,j,k,ipde1))
                    end do
                  end do
                end do
              end do
            end do

          else

            do ipde = 1,5
              do k = 1,km
                do j = 1,jm
                  do i = 1,im
                    term1 = ar * dtvol(i,j,k)
                    q(i,j,k,ipde) = ( &
                               ipkpi(i,j,k,ipde,1) * &
                               (-term1 * res(i,j,k,1) + qold(i,j,k,1)) + &
                               ipkpi(i,j,k,ipde,2) * &
                               (-term1 * res(i,j,k,2) + qold(i,j,k,2)) + &
                               ipkpi(i,j,k,ipde,3) * &
                               (-term1 * res(i,j,k,3) + qold(i,j,k,3)) + &
                               ipkpi(i,j,k,ipde,4) * &
                               (-term1 * res(i,j,k,4) + qold(i,j,k,4)) + &
                               ipkpi(i,j,k,ipde,5) * &
                               (-term1 * res(i,j,k,5) + qold(i,j,k,5)) )
                  end do
                end do
              end do
            end do

          end if

        else if ((irsop.gt.0).or.((irsop.eq.0).and.(nlevel.gt.1))) &
                                                                  then

!-----------------------------------------------------------------------
!       TIME-DEPENDENT WITH SIRK and WITH MG AND/OR IRS
!-----------------------------------------------------------------------

          if (kom.or.kom_bsl.or.kom_sst) then

            do k = 1,km
              do j = 1,jm
                do i = 1,im

                  gav        = ttrm(i,j,k,2)
                  bav        = ttrm(i,j,k,3)
                  dtau       = dtvolt(i,j,k) * vol(i,j,k)
          
                  a(1,1) = delpls(i,j,k) + &
                           bstrkw * q(i,j,k,7) / q(i,j,k,1)
                  a(1,2) = bstrkw * q(i,j,k,6) / q(i,j,k,1)
                  a(2,1) = 0.d0
                  a(2,2) = gav * delpls(i,j,k) + &
                           2 * bav * q(i,j,k,7) / q(i,j,k,1)

                  do ipde = 1,npde
                    do ipde1=1,npde
                      tmp(ipde,ipde1) = prec(i,j,k,ipde,ipde1) * fact/dt
                    end do
                    tmp(ipde,6) = tmp(ipde,6) + &
                                  prec(i,j,k,ipde,6) * a(1,1)
                    tmp(ipde,7) = tmp(ipde,7) + &
                                  prec(i,j,k,ipde,6) * a(1,2) + &
                                  prec(i,j,k,ipde,7) * a(2,2)
                    work1(i,j,k,ipde) = 0.d0
                    do ipde1=1,npde
                      work1(i,j,k,ipde) = work1(i,j,k,ipde) + &
                        tmp(ipde,ipde1) * q(i,j,k,ipde1)
                    end do
                  end do

                end do
              end do  
            end do  

!---------- update OMEGA and its residual
            ipde=7
            do k = 1,km
              do j = 1,jm
                do i = 1,im

                  gav        = ttrm(i,j,k,2)

                  term1 = ar * dtvolt(i,j,k)
                  fact1 = term1 * vol(i,j,k)
                  q(i,j,k,ipde) = 0
                  do ipde1 = 1,npde
                    q(i,j,k,ipde) = q(i,j,k,ipde) + &
                      ipkpi(i,j,k,ipde,ipde1) * &
                      (qold(i,j,k,ipde1) + fact1 * work1(i,j,k,ipde1) - &
                                     term1 * res(i,j,k,ipde1))
                  end do
                  rhomegamin = gav*cmu*q(i,j,k,1)*sqrt(prd(i,j,k))
                  if (q(i,j,k,ipde).lt.rhomegamin) then
                    q(i,j,k,ipde) = rhomegamin
                    somma = 0.d0
                    do ipde1 = 1,npde
                      somma = somma + lr(i,j,k,ipde1) * q(i,j,k,ipde1)
                    end do
                    res(i,j,k,ipde) = &
                      (-somma + qold(i,j,k,ipde) + &
                        fact1*work1(i,j,k,ipde) ) / term1
                  end if

                end do
              end do
            end do

!---------- update TKE
            ipde=6
            do k = 1,km
              do j = 1,jm
                do i = 1,im
                  term1 = ar * dtvolt(i,j,k)
                  fact1 = term1 * vol(i,j,k)
                  q(i,j,k,ipde) = 0.d0
                  do ipde1 = 1,npde
                    q(i,j,k,ipde) = q(i,j,k,ipde) + &
                      ipkpi(i,j,k,ipde,ipde1) * &
                      (qold(i,j,k,ipde1) + fact1 * work1(i,j,k,ipde1) - &
                                     term1 * res(i,j,k,ipde1))
                  end do
                end do
              end do
            end do

!---------- update NS variables
            do ipde = 1,5
              do k = 1,km
                do j = 1,jm
                  do i = 1,im
                    term1 = ar * dtvol(i,j,k)
                    fact1 = term1 * vol(i,j,k)
                    q(i,j,k,ipde) = 0.d0
                    do ipde1 = 1,npde
                      q(i,j,k,ipde) = q(i,j,k,ipde) + &
                        ipkpi(i,j,k,ipde,ipde1) * &
                        (qold(i,j,k,ipde1) + fact1 * work1(i,j,k,ipde1)- &
                                       term1 * res(i,j,k,ipde1))
                    end do
                  end do
                end do
              end do
            end do

          else

            do ipde = 1,5
              do k = 1,km
                do j = 1,jm
                  do i = 1,im
                    work1(i,j,k,ipde) = &
                         ( prec(i,j,k,ipde,1) * q(i,j,k,1) + &
                           prec(i,j,k,ipde,2) * q(i,j,k,2) + &
                           prec(i,j,k,ipde,3) * q(i,j,k,3) + &
                           prec(i,j,k,ipde,4) * q(i,j,k,4) + &
                           prec(i,j,k,ipde,5) * q(i,j,k,5) )
                  end do
                end do
              end do
            end do  

            do ipde = 1,5
              do k = 1,km
                do j = 1,jm
                  do i = 1,im
                    term1 = ar * dtvol(i,j,k)
                    fact1 = fact * vol(i,j,k) / dt
                    q(i,j,k,ipde) = ( &
                     ( ipkpi(i,j,k,ipde,1) * &
                       (qold(i,j,k,1) - term1 * (-fact1*work1(i,j,k,1) + &
                                               res(i,j,k,1)))        + &
                       ipkpi(i,j,k,ipde,2) * &
                       (qold(i,j,k,2) - term1 * (-fact1*work1(i,j,k,2) + &
                                               res(i,j,k,2)))        + &
                       ipkpi(i,j,k,ipde,3) * &
                       (qold(i,j,k,3) - term1 * (-fact1*work1(i,j,k,3) + &
                                               res(i,j,k,3)))        + &
                       ipkpi(i,j,k,ipde,4) * &
                       (qold(i,j,k,4) - term1 * (-fact1*work1(i,j,k,4) + &
                                               res(i,j,k,4)))        + &
                       ipkpi(i,j,k,ipde,5) * &
                       (qold(i,j,k,5) - term1 * (-fact1*work1(i,j,k,5) + &
                                               res(i,j,k,5)))) ) 
                  end do
                end do
              end do
            end do

          end if

        end if

      end if

      return
      end 

!-----------------------------------------------------------------------
      subroutine update_b(res,q,qold,vol,dtvol,dtvolt, &
                   delpls,prd,ttrm,ar,imax,jmax,kmax,npde,nharms)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none
!
      integer(kind=cosa_int) imax,jmax,kmax,npde,nharms
      integer(kind=cosa_int) im,jm,km,i,j,k,ipde,n
      real (kind=cosa_real) ar, fact, fact1, term1, dtau, rho, ome, tke, &
                    rhomega,drhomega,rhomegamin,gav,bav
      real (kind=cosa_real) &
           q     (-1:imax+1,-1:jmax+1,-1:kmax+1,npde), &
           delpls(-1:imax+1,-1:jmax+1,-1:kmax+1     ), &
           prd   (-1:imax+1,-1:jmax+1,-1:kmax+1     ), &
           qold  (-1:imax+1,-1:jmax+1,-1:kmax+1,npde), &
           res   (-1:imax+1,-1:jmax+1,-1:kmax+1,npde), &
           dtvol ( 0:imax  , 0:jmax  , 0:kmax       ), &
           dtvolt( 0:imax  , 0:jmax  , 0:kmax       ), &
           vol   ( 0:imax  , 0:jmax  , 0:kmax       ), &
           ttrm  (-1:imax+1,-1:jmax+1,-1:kmax+1,3   )
!
      im = imax-1
      jm = jmax-1
      km = kmax-1

!---- update solution

      if (.not.rkimp) then

        if (kom.or.kom_bsl.or.kom_sst) then
          do k = 1,km
            do j = 1,jm
              do i = 1,im
                rho        = q(i,j,k,1)
                ome        = q(i,j,k,7) / q(i,j,k,1)
                tke        = q(i,j,k,6) / q(i,j,k,1)
                dtau       = dtvolt(i,j,k) * vol(i,j,k)
                gav        = ttrm(i,j,k,2)
                bav        = ttrm(i,j,k,3)
                rhomega    = (qold(i,j,k,7) - &
                              ar*res(i,j,k,7)*dtvolt(i,j,k) + &
                              ar*dtau*rho*ome* &
                              (gav*delpls(i,j,k) + 2*bav*ome)) / &
                             (1 + ar*dtau* &
                              (gav*delpls(i,j,k) + 2*bav*ome))
                rhomegamin = gav*cmu*rho*sqrt(prd(i,j,k))
                if (rhomega.lt.rhomegamin) then
                  rhomega = rhomegamin
                end if
                q(i,j,k,7) = rhomega
                drhomega   = q(i,j,k,7) - rho*ome
                q(i,j,k,6) = (qold(i,j,k,6) - &
                              ar*res(i,j,k,6)*dtvolt(i,j,k) + &
                              ar*dtau*rho*tke* &
                              (delpls(i,j,k) + bstrkw*ome) - &
                              ar*dtau*bstrkw*tke*drhomega) / &
                             (1 + ar*dtau* &
                              (delpls(i,j,k) + bstrkw*ome))
              end do
            end do
          end do
        end if

        do ipde = 1,5
          do k = 1,km
            do j = 1,jm
              do i = 1,im
                q(i,j,k,ipde) = qold(i,j,k,ipde) - &
                                ar*res(i,j,k,ipde)*dtvol(i,j,k)
              end do
            end do
          end do
        end do

      else if (rkimp) then

        fact = 1.5d0
        if ((itime.eq.1).and.(.not.dual_prop_start)) then
          fact = 1.0d0
        end if

        if ((irsop.eq.0).and.(nlevel.eq.1)) then

          if (kom.or.kom_bsl.or.kom_sst) then
            do k = 1,km
              do j = 1,jm
                do i = 1,im
                  rho        = q(i,j,k,1)
                  ome        = q(i,j,k,7) / q(i,j,k,1)
                  tke        = q(i,j,k,6) / q(i,j,k,1)
                  dtau       = dtvolt(i,j,k) * vol(i,j,k)
                  gav        = ttrm(i,j,k,2)
                  bav        = ttrm(i,j,k,3)
                  rhomega    = (qold(i,j,k,7) - &
                                ar*res(i,j,k,7)*dtvolt(i,j,k) + &
                                ar*dtau*q(i,j,k,7)* &
                                (gav*delpls(i,j,k) + 2*bav*ome)) / &
                               (1 + ar*dtau* &
                                (gav*delpls(i,j,k) + 2*bav*ome) + &
                                fact*ar*dtau/dt)
                  rhomegamin = gav*cmu*rho*sqrt(prd(i,j,k))
                  if (rhomega.lt.rhomegamin) then
                    rhomega = rhomegamin
                  end if
                  q(i,j,k,7) = rhomega
                  drhomega = q(i,j,k,7) - rho*ome
                  q(i,j,k,6) = (qold(i,j,k,6) - &
                                ar*res(i,j,k,6)*dtvolt(i,j,k) + &
                                ar*dtau*q(i,j,k,6)* &
                                (delpls(i,j,k) + bstrkw*ome - &
                                bstrkw*drhomega/rho)) / &
                               (1 + ar*dtau* &
                                (delpls(i,j,k) + bstrkw*ome) + &
                                ar*fact*dtau/dt)
                end do
              end do
            end do
          end if

          do ipde = 1,5
            do k = 1,km
              do j = 1,jm
                do i = 1,im
                  q(i,j,k,ipde) = ( qold(i,j,k,ipde) - &
                                   ar*res(i,j,k,ipde)*dtvol(i,j,k) ) / &
                                 ( 1 + fact * ar * dtvol(i,j,k) * &
                                       vol(i,j,k) / dt )
                end do
              end do
            end do
          end do

        else if ((irsop.gt.0).or.((irsop.eq.0).and.nlevel.gt.1)) then

          if (kom.or.kom_bsl.or.kom_sst) then
            do k = 1,km
              do j = 1,jm
                do i = 1,im
                  rho        = q(i,j,k,1)
                  ome        = q(i,j,k,7) / q(i,j,k,1)
                  dtau       = dtvolt(i,j,k) * vol(i,j,k)
                  gav        = ttrm(i,j,k,2)
                  bav        = ttrm(i,j,k,3)
                  rhomega    = &
                    (qold(i,j,k,7) + ar*fact*dtau/dt*q(i,j,k,7) - &
                     ar*res(i,j,k,7)*dtvolt(i,j,k) + &
                     ar*dtau*q(i,j,k,7)* &
                     (gav*delpls(i,j,k) + 2*bav*ome)) / &
                    (1 + ar*dtau* &
                     (gav*delpls(i,j,k) + 2*bav*ome) + &
                     ar*fact*dtau/dt)
                  rhomegamin = gav*cmu*rho*sqrt(prd(i,j,k))
                  if (rhomega.lt.rhomegamin) then
                    rhomega = rhomegamin
                  end if
                  q(i,j,k,7)   = rhomega
                  drhomega   = q(i,j,k,7) - rho*ome
                  q(i,j,k,6)   = &
                    (qold(i,j,k,6) + ar*fact*dtau/dt*q(i,j,k,6) - &
                     ar*res(i,j,k,6)*dtvolt(i,j,k) + &
                     ar*dtau*q(i,j,k,6)* &
                     (delpls(i,j,k) + bstrkw*ome - &
                      bstrkw*drhomega/rho)) / &
                    (1 + ar*dtau* &
                     (delpls(i,j,k) + bstrkw*ome) + &
                     ar*fact*dtau/dt)
                end do
              end do
            end do
          end if

          do ipde = 1,5
            do k = 1,km
              do j = 1,jm
                do i = 1,im
                  q(i,j,k,ipde) = &
                    ( qold(i,j,k,ipde) - &
                      ar*dtvol(i,j,k) * &
                      (-fact*vol(i,j,k)/dt*q(i,j,k,ipde) + &
                       res(i,j,k,ipde) ) ) / &
                    ( 1 + fact * ar * dtvol(i,j,k) * vol(i,j,k) / dt )
                end do
              end do
            end do
          end do

        end if

      end if

      return
      end 

!-----------------------------------------------------------------------
      subroutine update_lsp_bhb(res,q,qold,vol,dtvol,dtvolt,prec,ipkpi, &
                  qhb0,lr,delpls,prd,ttrm,ar,imax,jmax,kmax,npde,nharms)
!-----------------------------------------------------------------------
!     mm    : numb. of rows    of matrix to be LU-factorized (DGETRF)
!     nn    : numb. of columns of matrix to be LU-factorized (DGETRF) OR
!             numb. of columns of matrix to be LU-factorized (DGETRS)
!     lda   : leading dimension of matrix to be LU-factorized (DGETRF),
!             .ge.max(1,nn) mm OR
!             leading dimension of matrix to be LU-factorized (DGETRS),
!             .ge.max(1,nn) nn
!     ldb   : leading dimension of RHS matrix (DGETRS) .ge.max(1,nn)
!     nrhs  : number of RHS's (DGETRS)
!     trans : DGETRS option to transpose matrix before mat.-vec. prod.
!             a) 'N': matrix is passed unaltered to DGETRS
!             b) 'T': transpose matrix is passed unaltered to DGETRS
!             c) 'C': conjugate transpose  is passed unaltered to DGETRS
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none
!
      integer(kind=cosa_int) imax,jmax,kmax,npde,nharms
      integer(kind=cosa_int) lda,info,mm,nn,ldb,nrhs
      character*1 trans
      integer(kind=cosa_int) im,jm,km,i,j,k,ipde,ipde1,ipde2,n,ni,nj,ni1, &
                nj1,sizeA,iend
      real (kind=cosa_real) ar, fact, fact1, term1
      real (kind=cosa_real) &
           q     (-1:imax+1,-1:jmax+1,-1:kmax+1,     npde,0:2*nharms), &
           qold  (-1:imax+1,-1:jmax+1,-1:kmax+1,     npde,0:2*nharms), &
           res   (-1:imax+1,-1:jmax+1,-1:kmax+1,     npde,0:2*nharms), &
           qhb0  (-1:imax+1,-1:jmax+1,-1:kmax+1,     npde,0:2*nharms), &
           dtvol ( 0:imax  , 0:jmax  , 0:kmax  ,          0:2*nharms), &
           dtvolt( 0:imax  , 0:jmax  , 0:kmax  ,          0:2*nharms), &
           lr    (-1:imax+1,-1:jmax+1,-1:kmax+1,     npde,0:2*nharms), &
           delpls(-1:imax+1,-1:jmax+1,-1:kmax+1,          0:2*nharms), &
           prd   (-1:imax+1,-1:jmax+1,-1:kmax+1,          0:2*nharms), &
           ttrm  (-1:imax+1,-1:jmax+1,-1:kmax+1,        3,0:2*nharms), &
           prec  (   imax-1,   jmax-1,   kmax-1,npde,npde,0:2*nharms), &
           ipkpi (   imax-1,   jmax-1,   kmax-1,npde,npde,0:2*nharms), &
           vol   ( 0:imax  , 0:jmax  , 0:kmax                       )
      integer(kind=cosa_int), allocatable:: ipiv(:)
      real (kind=cosa_real) tke,dtau,rhomega,rhomegamin,work(7), &
           tmp(7,7),somma,dps,rho,ome,drhomega,gav,bav
      real (kind=cosa_real), allocatable::hbstb(:,:),a(:,:,:), &
           id1(:,:),qhb(:),temp1(:,:),temp2(:,:),q_l(:,:),rho_l(:), &
           qold_l(:,:),prec_l(:,:,:),res_l(:,:)

      sizeA  = npde * (2*nharms + 1)
      iend   = sizeA - 1
      allocate( hbstb(0:iend,0:iend),id1(0:iend,0:iend),qhb(0:iend), &
                a(2,2,0:2*nharms), &
                temp1(0:2*nharms,0:2*nharms),temp2(0:iend,0:iend), &
                q_l(npde,0:2*nharms),rho_l(0:2*nharms), &
                qold_l(npde,0:2*nharms),res_l(npde,0:2*nharms), &
                prec_l(npde,npde,0:2*nharms),ipiv(0:iend) )

      im = imax-1
      jm = jmax-1
      km = kmax-1

!---- update solution

      if (.not.rkimp) then

!-----------------------------------------------------------------------
!     FERK INTEGRATION
!-----------------------------------------------------------------------

        if (kom.or.kom_bsl.or.kom_sst) then

          do n = 0,2*nharms
          
            do k = 1,km
              do j = 1,jm
                do i = 1,im

                gav        = ttrm(i,j,k,2,n)
                bav        = ttrm(i,j,k,3,n)

                a(1,1,n) = delpls(i,j,k,n) + &
                           bstrkw * q(i,j,k,7,n) / q(i,j,k,1,n)
                a(1,2,n) = bstrkw * q(i,j,k,6,n) / q(i,j,k,1,n)
                a(2,1,n) = 0.d0
                a(2,2,n) = gav * delpls(i,j,k,n) &
                                 + 2 * bav * q(i,j,k,7,n) / q(i,j,k,1,n)

                do ipde = 1,npde
                  tmp(ipde,6) = prec(i,j,k,ipde,6,n) * a(1,1,n)
                  tmp(ipde,7) = prec(i,j,k,ipde,6,n) * a(1,2,n) + &
                                      prec(i,j,k,ipde,7,n) * a(2,2,n)
!del              work1(i,j,ipde,n) =
                  work(ipde) = &
                    tmp(ipde,6) * q(i,j,k,6,n) + &
                    tmp(ipde,7) * q(i,j,k,7,n)
                end do

!-------------- update OMEGA and its residual
                ipde=7
                term1 = ar * dtvolt(i,j,k,n)
                fact1 = term1 * vol(i,j,k)
                q(i,j,k,ipde,n) = 0
                do ipde1 = 1,npde
                  q(i,j,k,ipde,n) = q(i,j,k,ipde,n) + &
                    ipkpi(i,j,k,ipde,ipde1,n) * &
                    (qold(i,j,k,ipde1,n) + fact1 * work(ipde1) - &
                                   term1 * res(i,j,k,ipde1,n))
                end do
                rhomegamin = gav*cmu*q(i,j,k,1,n)*sqrt(prd(i,j,k,n))
                if (q(i,j,k,ipde,n).lt.rhomegamin) then
                  q(i,j,k,ipde,n) = rhomegamin
                  somma = 0.d0
                  do ipde1 = 1,npde
                    somma = somma + lr(i,j,k,ipde1,n) * q(i,j,k,ipde1,n)
                  end do
                  res(i,j,k,ipde,n) = &
                    (-somma +qold(i,j,k,ipde,n) +fact1*work(ipde) ) / &
                    term1
                end if

!-------------- update TKE
                ipde=6
                q(i,j,k,ipde,n) = 0.d0
                do ipde1 = 1,npde
                  q(i,j,k,ipde,n) = q(i,j,k,ipde,n) + &
                    ipkpi(i,j,k,ipde,ipde1,n) * &
!del                (qold(i,j,ipde1,n) + fact1 * work1(i,j,ipde1,n) - &
                    (qold(i,j,k,ipde1,n) + fact1 * work(ipde1) - &
                                   term1 * res(i,j,k,ipde1,n))
                end do

!-------------- update NS variables
                do ipde = 1,5
                  term1 = ar * dtvol(i,j,k,n)
                  fact1 = term1 * vol(i,j,k)
                  q(i,j,k,ipde,n) = 0.d0
                  do ipde1 = 1,npde
                    q(i,j,k,ipde,n) = q(i,j,k,ipde,n) + &
                      ipkpi(i,j,k,ipde,ipde1,n) * &
!del                  (qold(i,j,ipde1,n) + fact1 * work1(i,j,ipde1,n) - &
                      (qold(i,j,k,ipde1,n) + fact1 * work(ipde1) - &
                                     term1 * res(i,j,k,ipde1,n))
                  end do
                end do

                end do
              end do
            end do

          end do

        else
!------ Euler/laminar FERK case

          do n = 0,2*nharms
            do ipde = 1,npde
              do k = 1,km
                do j = 1,jm
                  do i = 1,im
                    q(i,j,k,ipde,n) = qold(i,j,k,ipde,n) - &
                                    ar*res(i,j,k,ipde,n)*dtvol(i,j,k,n)
                  end do
                end do
              end do
            end do
          end do

        end if

      else if (rkimp) then

!-----------------------------------------------------------------------
!     PIRK INTEGRATION
!-----------------------------------------------------------------------

        if (kom.or.kom_bsl.or.kom_sst) then

          do nj = 0,iend
            do ni = 0,iend
              temp2(ni,nj) = 0.d0
              if (ni.eq.nj) then
                id1(ni,nj) = 1.d0
              else
                id1(ni,nj) = 0.d0
              end if
            end do
          end do

          if ((irsop.eq.0).and.(nlevel.eq.1)) then

!           TO BE COMPLETED

          else if ((irsop.gt.0).or.((irsop.eq.0).and.(nlevel.gt.1))) then

            trans = 'N'
            mm = sizeA
            nn = sizeA
            lda = mm
            ldb = nn
            nrhs = 1

            do k = 1,kmax-1
              do j = 1,jmax-1
                do i = 1,imax-1

                do ni = 0,2*nharms

                  gav        = ttrm(i,j,k,2,ni)
                  bav        = ttrm(i,j,k,3,ni)

                  a(1,1,ni) = delpls(i,j,k,ni) + &
                             bstrkw * q(i,j,k,7,ni) / q(i,j,k,1,ni)
                  a(1,2,ni) = bstrkw * q(i,j,k,6,ni) / q(i,j,k,1,ni)
                  a(2,1,ni) = 0.d0
                  a(2,2,ni) = gav * delpls(i,j,k,ni) &
                                   + 2 * bav * q(i,j,k,7,ni) / q(i,j,k,1,ni)

                  do nj = 0,2*nharms
                    temp1(ni,nj) = ar*omega*     dhb(ni,nj)
                  end do

!---------------- initialize local arrays
                  do ipde2 = 1,npde
                    do ipde1 = 1,npde
                      prec_l(ipde1,ipde2,ni) = prec(i,j,k,ipde1,ipde2,ni)
                    end do
                  end do
                  rho_l(ni) = q(i,j,k,1,ni)
                  do ipde1 = 1,npde
                    qold_l(ipde1,ni) = qold(i,j,k,ipde1,ni)
                    res_l (ipde1,ni) = res (i,j,k,ipde1,ni)
                    q_l(ipde1,ni) = q(i,j,k,ipde1,ni)
                  end do

                end do

                do ni = 0,2*nharms
                  do ipde = 1,npde
                    tmp(ipde,6) = prec_l(ipde,6,ni) * a(1,1,ni)
                    tmp(ipde,7) = prec_l(ipde,6,ni) * a(1,2,ni) + &
                                  prec_l(ipde,7,ni) * a(2,2,ni)
                  end do
                  do nj = 0,2*nharms
                    do ipde1 = 1,npde
                      if (ipde1.le.5) then
                        dtau = dtvol(i,j,k,ni) * vol(i,j,k)
                      else
                        dtau = dtvolt(i,j,k,ni) * vol(i,j,k)
                      end if
                      ni1 = (ipde1-1) + ni * npde
                      do ipde2 = 1,npde
                        nj1 = (ipde2-1) + nj * npde
                        temp2(ni1,nj1) = dtau*prec_l(ipde1,ipde2,ni) * &
                                         temp1(ni,nj)
!old                    if((ni.eq.nj).and.(ipde1.ge.5).and.(ipde2.ge.5)) then
!old                    if ((ni.eq.nj).and.(ipde2.ge.5)) then
                        if (nj.eq.ni) then
                          temp2(ni1,nj1) = temp2(ni1,nj1) + &
                                           ar * dtau*tmp(ipde1,ipde2)
                        end if
                        hbstb(ni1,nj1) = id1(ni1,nj1) + temp2(ni1,nj1)
                      end do
                    end do
                  end do
                end do

                call DGETRF(mm,nn,hbstb,lda,ipiv,info)
!dbg            write(*,*) 'info DGETRF is:',info

!-------------- complete RHS construction
!del            do ni = 0,iend
!del              qhb(ni) = 0.d0
!del            end do
                do ni = 0,2*nharms
                  do ipde1 = 1,npde
                    ni1 = (ipde1-1) + ni * npde
                    qhb(ni1) = 0.d0
                    do ipde2 = 1,npde
                      qhb(ni1) = qhb(ni1) + &
                        prec_l(ipde1,ipde2,ni) * qhb0(i,j,k,ipde2,ni)
                    end do

                    if (ipde1.le.5) then
                      qhb(ni1) = &
                        qold_l(ipde1,ni) + ar*dtvol(i,j,k,ni) * &
                        ( qhb(ni1) - res_l(ipde1,ni) )
                    else
                      qhb(ni1) = &
                        qold_l(ipde1,ni) + ar*dtvolt(i,j,k,ni) * &
                        ( qhb(ni1) - res_l(ipde1,ni) )
                    end if

                    tmp(ipde1,6) = prec_l(ipde1,6,ni) * a(1,1,ni)
                    tmp(ipde1,7) = prec_l(ipde1,6,ni) * a(1,2,ni) + &
                                   prec_l(ipde1,7,ni) * a(2,2,ni)

                    if (ipde1.le.5) then
                      dtau = dtvol(i,j,k,ni) * vol(i,j,k)
                    else
                      dtau = dtvolt(i,j,k,ni) * vol(i,j,k)
                    end if
                    qhb(ni1) = qhb(ni1) + &
                      ar*dtau * ( tmp(ipde1,6)*q_l(6,ni) + &
                                  tmp(ipde1,7)*q_l(7,ni) )

                  end do
                end do
 
!dbg            do k = 0,iend
!dbg              if (ipiv(k).ne.(k+1)) then
!dbg                write(3,*) 'using the pivoting'
!dbg              end if
!dbg            end do

                call DGETRS(trans,nn,nrhs,hbstb,lda,ipiv, &
                            qhb,ldb,info)
!dbg            write(*,*) 'info DGETRS is:',info                

                do ni = 0,2*nharms
                  do ipde = 1,npde
                    ni1 = (ipde-1) + ni * npde
                    if (ipde.eq.7) then
                      gav        = ttrm(i,j,k,2,ni)
                      rhomegamin = gav*cmu*rho_l(ni)*sqrt(prd(i,j,k,ni))
                      if (qhb(ni1).lt.rhomegamin) qhb(ni1) = rhomegamin
                    end if
                    q(i,j,k,ipde,ni) = qhb(ni1)
                  end do
                end do

                end do
              end do
            end do

          end if

        else

          do nj = 0,iend
            do ni = 0,iend
              temp2(ni,nj) = 0.d0
              if (ni.eq.nj) then
                id1(ni,nj) = 1.d0
              else
                id1(ni,nj) = 0.d0
              end if
            end do
          end do

          if ((irsop.eq.0).and.(nlevel.eq.1)) then

            trans = 'N'
            mm = sizeA
            nn = sizeA
            lda = mm
            ldb = nn
            nrhs = 1

            do k = 1,kmax-1
              do j = 1,jmax-1
                do i = 1,imax-1

                do ni = 0,2*nharms
                  dtau = dtvol(i,j,k,ni) * vol(i,j,k)
                  do nj = 0,2*nharms
                    temp1(ni,nj) = ar*omega*dtau*dhb(ni,nj)
                  end do

!---------------- initialize local arrays
                  do ipde2 = 1,npde
                    do ipde1 = 1,npde
                      prec_l(ipde1,ipde2,ni) = prec(i,j,k,ipde1,ipde2,ni)
                    end do
                    qold_l(ipde2,ni) = qold(i,j,k,ipde2,ni)
                    res_l (ipde2,ni) = res (i,j,k,ipde2,ni)
                  end do

                end do

                do ni = 0,2*nharms
                  do nj = 0,2*nharms
                    do ipde1 = 1,npde
                      ni1 = (ipde1-1) + ni * npde
                      do ipde2 = 1,npde
                        nj1 = (ipde2-1) + nj * npde
                        temp2(ni1,nj1) = prec_l(ipde1,ipde2,ni) * &
                                         temp1(ni,nj)
                        hbstb(ni1,nj1) = id1(ni1,nj1) + temp2(ni1,nj1)
                      end do
                    end do
                  end do
                end do

                call DGETRF(mm,nn,hbstb,lda,ipiv,info)
!dbg            write(*,*) 'info DGETRF is:',info

!-------------- calculate RHS of counterpsrt eq. (29) of gt-2011-45303
!               for case without IRS and without MG
                do ni = 0,2*nharms
                  do ipde1 = 1,npde
                    ni1 = (ipde1-1) + ni * npde
                    qhb(ni1) = qold_l(ipde1,ni) - &
                      ar*dtvol(i,j,k,ni) * res_l(ipde1,ni)
                  end do
                end do

                call DGETRS(trans,nn,nrhs,hbstb,lda,ipiv, &
                            qhb,ldb,info)
!dbg            write(*,*) 'info DGETRS is:',info

                do ni = 0,2*nharms
                  do ipde = 1,npde
                    ni1 = (ipde-1) + ni * npde
                    q(i,j,k,ipde,ni) = qhb(ni1)
                  end do
                end do

                end do
              end do
            end do

          else if ((irsop.gt.0).or.((irsop.eq.0).and.(nlevel.gt.1))) &
                   then

!-------- msc 14th August 2009: this update requires the solution at
!                               preceding RK stage, which is stored in
!                               q when entering this block. So there is
!                               no choice but to use work1. This could
!                               be avoided only for the case (irsop=0,
!                               nlevel=0), which is the previous block.
!                               Moreover, the updates of the pde's
!                               are now coupled and can no longer be
!                               performed separately like in the case
!                               with HB and withOUT LSP.

            trans = 'N'
            mm = sizeA
            nn = sizeA
            lda = mm
            ldb = nn
            nrhs = 1

            do k = 1,kmax-1
              do j = 1,jmax-1
                do i = 1,imax-1

                do ni = 0,2*nharms

                  dtau = dtvol(i,j,k,ni) * vol(i,j,k)
                  do nj = 0,2*nharms
                    temp1(ni,nj) = ar*omega*dtau*dhb(ni,nj)
                  end do

!---------------- initialize local arrays
                  do ipde2 = 1,npde
                    do ipde1 = 1,npde
                      prec_l(ipde1,ipde2,ni) = prec(i,j,k,ipde1,ipde2,ni)
                    end do
                    qold_l(ipde2,ni) = qold(i,j,k,ipde2,ni)
                    res_l (ipde2,ni) = res (i,j,k,ipde2,ni)
                    q_l(ipde2,ni) = q(i,j,k,ipde2,ni)
                  end do

                end do

                do ni = 0,2*nharms
                  do nj = 0,2*nharms
                    do ipde1 = 1,npde
                      ni1 = (ipde1-1) + ni * npde
                      do ipde2 = 1,npde
                        nj1 = (ipde2-1) + nj * npde
                        temp2(ni1,nj1) = prec_l(ipde1,ipde2,ni) * &
                                         temp1(ni,nj)
                        hbstb(ni1,nj1) = id1(ni1,nj1) + temp2(ni1,nj1)
                      end do
                    end do
                  end do
                end do

                call DGETRF(mm,nn,hbstb,lda,ipiv,info)
!dbg            write(*,*) 'info DGETRF is:',info

!-------------- calculate RHS of eq. (29) of gt-2011-45303
!del            do ni = 0,iend
!del              qhb(ni) = 0.d0
!del            end do
                do ni = 0,2*nharms
                  do ipde1 = 1,npde
                    ni1 = (ipde1-1) + ni * npde
                    qhb(ni1) = 0.d0
                    do ipde2 = 1,npde
                      qhb(ni1) = qhb(ni1) + &
                        prec_l(ipde1,ipde2,ni) * qhb0(i,j,k,ipde2,ni)
                    end do
                    qhb(ni1) = &
                      qold_l(ipde1,ni) + ar*dtvol(i,j,k,ni) * &
                      ( qhb(ni1) - res_l(ipde1,ni) )
                  end do
                end do

!dbg            do k = 0,iend
!dbg              if (ipiv(k).ne.(k+1)) then
!dbg                write(3,*) 'using the pivoting'
!dbg              end if
!dbg            end do

                call DGETRS(trans,nn,nrhs,hbstb,lda,ipiv, &
                            qhb,ldb,info)
!dbg            write(*,*) 'info DGETRS is:',info                

                do ni = 0,2*nharms
                  do ipde = 1,npde
                    ni1 = (ipde-1) + ni * npde
                    q(i,j,k,ipde,ni) = qhb(ni1)
                  end do
                end do

                end do
              end do
            end do

          end if

        end if

      end if

      deallocate( hbstb,id1,qhb,a,temp1,temp2, &
                  q_l,rho_l,qold_l,res_l,prec_l,ipiv)

      return
      end

!-----------------------------------------------------------------------
      subroutine update_bhb(res,q,qold,vol,dtvol,dtvolt,work1,qhb0, &
                   delpls,prd,ttrm,ar,imax,jmax,kmax,npde,nharms)
!-----------------------------------------------------------------------
!     mm    : numb. of rows    of matrix to be LU-factorized (DGETRF)
!     nn    : numb. of columns of matrix to be LU-factorized (DGETRF) OR
!             numb. of columns of matrix to be LU-factorized (DGETRS)
!     lda   : leading dimension of matrix to be LU-factorized (DGETRF),
!             .ge.max(1,nn) mm OR
!             leading dimension of matrix to be LU-factorized (DGETRS),
!             .ge.max(1,nn) nn
!     ldb   : leading dimension of RHS matrix (DGETRS) .ge.max(1,nn)
!     nrhs  : number of RHS's (DGETRS)
!     trans : DGETRS option to transpose matrix before mat.-vec. prod.
!             a) 'N': matrix is passed unaltered to DGETRS
!             b) 'T': transpose matrix is passed unaltered to DGETRS
!             c) 'C': conjugate transpose  is passed unaltered to DGETRS
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none
!
      integer(kind=cosa_int) imax,jmax,kmax,npde,nharms,off
      integer(kind=cosa_int) lda,info,mm,nn,ldb,nrhs
      character*1 trans
      integer(kind=cosa_int) im,jm,km,i,j,k,ipde,ipde1,ipde2,n,ni,nj,ni1,nj1
      real (kind=cosa_real) ar, fact, fact1, term1
      real (kind=cosa_real) &
           q     (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           qold  (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           res   (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           qhb0  (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           work1 (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           dtvol ( 0:imax  , 0:jmax  , 0:kmax  ,     0:2*nharms), &
           dtvolt( 0:imax  , 0:jmax  , 0:kmax  ,     0:2*nharms), &
           delpls(-1:imax+1,-1:jmax+1,-1:kmax+1,     0:2*nharms), &
           prd   (-1:imax+1,-1:jmax+1,-1:kmax+1,     0:2*nharms), &
           ttrm  (-1:imax+1,-1:jmax+1,-1:kmax+1,3   ,0:2*nharms), &
           vol   ( 0:imax  , 0:jmax  , 0:kmax                  )
      integer(kind=cosa_int), allocatable:: ipiv(:)
      real (kind=cosa_real) tke,dtau,rhomega,rhomegamin,dps,rho,ome,drhomega, &
                    gav,bav
      real (kind=cosa_real), allocatable::hbstb(:,:), &
           ID1(:,:),qhb(:),temp1(:,:),temp2(:,:),q_l(:,:), &
           qold_l(:,:),res_l(:,:),work2(:,:)

      allocate( hbstb(0:2*nharms,0:2*nharms), &
                ID1(0:2*nharms,0:2*nharms), &
                qhb(0:2*nharms),work2(0:2*nharms,1:npde), &
                ipiv(0:2*nharms))

      im = imax-1
      jm = jmax-1
      km = kmax-1

!---- update solution

      if (.not.rkimp) then

        do n = 0,2*nharms

          if (kom.or.kom_bsl.or.kom_sst) then
            do k = 1,km
              do j = 1,jm
                do i = 1,im
                  gav        = ttrm(i,j,k,2,n)
                  bav        = ttrm(i,j,k,3,n)
                  rho        = q(i,j,k,1,n)
                  ome        = q(i,j,k,7,n) / q(i,j,k,1,n)
                  tke        = q(i,j,k,6,n) / q(i,j,k,1,n)
                  dtau       = dtvolt(i,j,k,n) * vol(i,j,k)
                  rhomega    = (qold(i,j,k,7,n) - &
                                ar*res(i,j,k,7,n)*dtvolt(i,j,k,n) + &
                                ar*dtau*rho*ome* &
                                (gav*delpls(i,j,k,n) + 2*bav*ome)) / &
                               (1 + ar*dtau* &
                                (gav*delpls(i,j,k,n) + 2*bav*ome))
                  rhomegamin = gav*cmu*rho*sqrt(prd(i,j,k,n))
                  if (rhomega.lt.rhomegamin) then
                    rhomega = rhomegamin
                  end if
                  q(i,j,k,7,n) = rhomega
                  drhomega   = q(i,j,k,7,n) - rho*ome
                  q(i,j,k,6,n) = (qold(i,j,k,6,n) - &
                                ar*res(i,j,k,6,n)*dtvolt(i,j,k,n) + &
                                ar*dtau*rho*tke* &
                                (delpls(i,j,k,n) + bstrkw*ome) - &
                                ar*dtau*bstrkw*tke*drhomega) / &
                               (1 + ar*dtau* &
                                (delpls(i,j,k,n) + bstrkw*ome))
                end do
              end do
            end do
          end if

          do ipde = 1,5
            do k = 1,km
              do j = 1,jm
                do i = 1,im
                  q(i,j,k,ipde,n) = qold(i,j,k,ipde,n) - &
                                    ar*res(i,j,k,ipde,n)*dtvol(i,j,k,n)
                end do
              end do
            end do
          end do

        end do

      else if (rkimp) then

        do nj=0,2*nharms
          do ni=0,2*nharms
            if (ni.eq.nj) then
              ID1(ni,nj) = 1.d0
            else
              ID1(ni,nj) = 0.d0
            end if
          end do
        end do

        if ((irsop.eq.0).and.(nlevel.eq.1)) then
!         case without IRS and without MG

          if (kom.or.kom_bsl.or.kom_sst) then

!---------- UPDATE (rho omega)

!---------  1. calculate RHS of counterpsrt eq. (21) of gt-2011-45303
!              In present case of omega equation there are additional
!              terms.
!           2. store previous rho*omega in work1(:,:,7,:)
            do ni = 0,2*nharms
              do k = 1,kmax-1
                do j = 1,jmax-1
                  do i = 1,imax-1
                    gav            = ttrm(i,j,k,2,ni)
                    bav            = ttrm(i,j,k,3,ni)
                    ome            = q(i,j,k,7,ni) / q(i,j,k,1,ni)
                    dtau           = dtvolt(i,j,k,ni) * vol(i,j,k)
                    qhb0(i,j,k,7,ni) = qold(i,j,k,7,ni) - &
                      ar*res(i,j,k,7,ni)*dtvolt(i,j,k,ni) + &
                      ar*dtau*q(i,j,k,7,ni) * &
                      (gav*delpls(i,j,k,ni) + 2*bav*ome)
                    work1(i,j,k,7,ni) = q(i,j,k,7,ni)
                  end do
                end do
              end do
            end do

            trans = 'N'
            mm = 2*nharms + 1
            nn = 2*nharms + 1
            lda = mm
            ldb = nn
            nrhs = 1

            do k = 1,kmax-1
              do j = 1,jmax-1
                do i = 1,imax-1

                  do ni = 0,2*nharms
                    gav      = ttrm(i,j,k,2,ni)
                    bav      = ttrm(i,j,k,3,ni)
                    dtau     = dtvolt(i,j,k,ni) * vol(i,j,k)
                    ome      = q(i,j,k,7,ni) / q(i,j,k,1,ni)
                    dps      = delpls(i,j,k,ni)
                    do nj = 0,2*nharms
                      hbstb(ni,nj) = ID1(ni,nj) * &
                        (1 + ar*dtau * (gav*dps +  2*bav*ome) ) + &
                        ar*omega*dtau*dhb(ni,nj)
                    end do
                  end do

                  call DGETRF(mm,nn,hbstb,lda,ipiv,info)
!dbg              write(*,*) 'info DGETRF is:',info

!---------------- store RHS in work2.
                  do ni = 0,2*nharms
                    work2(ni,1) = qhb0(i,j,k,7,ni)
                  end do

                  call DGETRS(trans,nn,nrhs,hbstb,lda,ipiv, &
                              work2,ldb,info)
!dbg              write(*,*) 'info DGETRS is:',info

!---------------- 1. updated omega is in work2. Before copying it back to 
!                    q, apply limitation.
!                 2. store difference between new and previous values 
!                    of rho*omega in work1(:,:,6,:)
                  do ni = 0,2*nharms
                    gav        = ttrm(i,j,k,2,ni)
                    rho        = q(i,j,k,1,ni)
                    rhomegamin = gav*cmu*rho*sqrt(prd(i,j,k,ni))
                    if (work2(ni,1).lt.rhomegamin) then
                      work2(ni,1) = rhomegamin
                    end if
                    q(i,j,k,7,ni) = work2(ni,1)
                    work1(i,j,k,6,ni) = &
                      q(i,j,k,7,ni) - work1(i,j,k,7,ni)
                  end do

                end do
              end do
            end do

!---------- UPDATE (rho K)

!---------  calculate RHS of counterpsrt eq. (21) of gt-2011-45303
!           In present case of K equation there are additional terms.
            do ni = 0,2*nharms
              do k = 1,kmax-1
                do j = 1,jmax-1
                  do i = 1,imax-1
                    rho              = q(i,j,k,1,ni)
                    ome              = work1(i,j,k,7,ni) / q(i,j,k,1,ni)
                    dtau             = dtvolt(i,j,k,ni) * vol(i,j,k)
                    drhomega         = work1(i,j,k,6,ni)
                    qhb0(i,j,k,6,ni) = (qold(i,j,k,6,ni) - &
                      ar*res(i,j,k,6,ni)*dtvolt(i,j,k,ni) + &
                      ar*dtau*q(i,j,k,6,ni)* &
                      (delpls(i,j,k,ni) + bstrkw*ome - &
                      bstrkw*drhomega/rho))
                  end do
                end do
              end do
            end do

            do k = 1,kmax-1
              do j = 1,jmax-1
                do i = 1,imax-1

                  do ni = 0,2*nharms
                    dtau     = dtvolt(i,j,k,ni) * vol(i,j,k)
                    ome      = work1(i,j,k,7,ni) / q(i,j,k,1,ni)
                    dps      = delpls(i,j,k,ni)
                    do nj = 0,2*nharms
                      hbstb(ni,nj) = ID1(ni,nj) * &
                        (1 + ar*dtau * (    dps +    bstrkw*ome) ) + &
                        ar*omega*dtau*dhb(ni,nj)
                    end do
                  end do

                  call DGETRF(mm,nn,hbstb,lda,ipiv,info)
!dbg              write(*,*) 'info DGETRF is:',info

!---------------- store RHS in work2.
                  do ni = 0,2*nharms
                    work2(ni,1) = qhb0(i,j,k,6,ni)
                  end do

                  call DGETRS(trans,nn,nrhs,hbstb,lda,ipiv, &
                              work2,ldb,info)
!dbg              write(*,*) 'info DGETRS is:',info

!---------------- store updated rho*K in q(:,:,6,:).
                  do ni = 0,2*nharms
                    q(i,j,k,6,ni) = work2(ni,1)
                  end do

                end do
              end do
            end do

          end if

!-------- UPDATE Euler / NS variables

!-------- calculate RHS of counterpsrt eq. (21) of gt-2011-45303 for
          do ni = 0,2*nharms
            do ipde = 1,5
              do k = 1,kmax-1
                do j = 1,jmax-1
                  do i = 1,imax-1
                    qhb0(i,j,k,ipde,ni) = qold(i,j,k,ipde,ni) - &
                        ar*res(i,j,k,ipde,ni)*dtvol(i,j,k,ni)
                  end do
                end do
              end do
            end do
          end do

          trans = 'N'
          mm = 2*nharms + 1
          nn = 2*nharms + 1
          lda = mm
          ldb = nn
          nrhs = 5

          do k = 1,kmax-1
            do j = 1,jmax-1
              do i = 1,imax-1

                do ni = 0,2*nharms
                  dtau = dtvol(i,j,k,ni) * vol(i,j,k)
                  do nj = 0,2*nharms
                    hbstb(ni,nj) = ID1(ni,nj) + &
                                     ar*omega*dtau*dhb(ni,nj)
                  end do
                end do

                call DGETRF(mm,nn,hbstb,lda,ipiv,info)
!dbg            write(*,*) 'info DGETRF is:',info

!               transpose RHS matrix, to have pattern [(2nharms+1)xnpde]
                do ni = 0,2*nharms
                  do ipde = 1,5
                    work2(ni,ipde) = qhb0(i,j,k,ipde,ni)
                  end do
                end do

                call DGETRS(trans,nn,nrhs,hbstb,lda,ipiv, &
                            work2,ldb,info)
!dbg            write(*,*) 'info DGETRS is:',info

!-------------- copy new solution back to q
                do ni = 0,2*nharms
                  do ipde = 1,5
                    q(i,j,k,ipde,ni) = work2(ni,ipde)
                  end do
                end do

              end do
            end do
          end do

        else if ((irsop.gt.0).or.((irsop.eq.0).and.nlevel.gt.1)) then
!         case with IRS and/or MG

          if (kom.or.kom_bsl.or.kom_sst) then

!---------- UPDATE (rho omega)

!---------  1. calculate RHS of counterpsrt eq. (22) of gt-2011-45303
!              In present case of omega equation there are additional
!              terms.
!           2. store previous rho*omega in work1(:,:,7,:)
            do ni = 0,2*nharms
              do k = 1,kmax-1
                do j = 1,jmax-1
                  do i = 1,imax-1
                    gav            = ttrm(i,j,k,2,ni)
                    bav            = ttrm(i,j,k,3,ni)
                    ome            = q(i,j,k,7,ni) / q(i,j,k,1,ni)
                    dtau           = dtvolt(i,j,k,ni) * vol(i,j,k)
                    qhb0(i,j,k,7,ni) = &
                      qold(i,j,k,7,ni) + ar*dtvolt(i,j,k,ni) * &
                      (qhb0(i,j,k,7,ni) - res(i,j,k,7,ni)) + &
                      ar*dtau*q(i,j,k,7,ni)* &
                      (gav*delpls(i,j,k,ni) + 2*bav*ome)
                    work1(i,j,k,7,ni) = q(i,j,k,7,ni)
                  end do
                end do
              end do
            end do

            trans = 'N'
            mm = 2*nharms + 1
            nn = 2*nharms + 1
            lda = mm
            ldb = nn
            nrhs = 1

            do k = 1,kmax-1
              do j = 1,jmax-1
                do i = 1,imax-1

                  do ni = 0,2*nharms
                    gav      = ttrm(i,j,k,2,ni)
                    bav      = ttrm(i,j,k,3,ni)
                    dtau     = dtvolt(i,j,k,ni) * vol(i,j,k)
                    ome      = q(i,j,k,7,ni) / q(i,j,k,1,ni)
                    dps      = delpls(i,j,k,ni)
                    do nj = 0,2*nharms
                      hbstb(ni,nj) = ID1(ni,nj) * &
                        (1 + ar*dtau * (gav*dps +  2*bav*ome) ) + &
                        ar*omega*dtau*dhb(ni,nj)
                    end do
                  end do

                  call DGETRF(mm,nn,hbstb,lda,ipiv,info)
!dbg              write(*,*) 'info DGETRF is:',info

!---------------- store RHS in work2.
                  do ni = 0,2*nharms
                    work2(ni,1) = qhb0(i,j,k,7,ni)
                  end do

                  call DGETRS(trans,nn,nrhs,hbstb,lda,ipiv, &
                              work2,ldb,info)
!dbg              write(*,*) 'info DGETRS is:',info

!---------------- 1. updated omega is in work2. Before copying it back to 
!                    q, apply limitation.
!                 2. store difference between new and previous values 
!                    of rho*omega in work1(:,:,6,:)
                  do ni = 0,2*nharms
                    gav      = ttrm(i,j,k,2,ni)
                    rho        = q(i,j,k,1,ni)
                    rhomegamin = gav*cmu*rho*sqrt(prd(i,j,k,ni))
                    if (work2(ni,1).lt.rhomegamin) then
                      work2(ni,1) = rhomegamin
                    end if
                    q(i,j,k,7,ni) = work2(ni,1)
                    work1(i,j,k,6,ni) = q(i,j,k,7,ni) - &
                                        work1(i,j,k,7,ni)
                  end do

                end do
              end do
            end do

!---------- UPDATE (rho K)

!---------  calculate RHS of counterpsrt eq. (22) of gt-2011-45303
!           In present case of K equation there are additional terms.
            do ni = 0,2*nharms
              do k = 1,kmax-1
                do j = 1,jmax-1
                  do i = 1,imax-1
                    rho            = q(i,j,k,1,ni)
                    ome            = work1(i,j,k,7,ni) / q(i,j,k,1,ni)
                    dtau           = dtvolt(i,j,k,ni) * vol(i,j,k)
                    drhomega       = work1(i,j,k,6,ni)
                    qhb0(i,j,k,6,ni) = (qold(i,j,k,6,ni) + &
                      ar*dtvolt(i,j,k,ni) * &
                      (qhb0(i,j,k,6,ni) - res(i,j,k,6,ni)) + &
                      ar*dtau*q(i,j,k,6,ni) * &
                      (delpls(i,j,k,ni) + bstrkw*ome - &
                      bstrkw*drhomega/rho))
                  end do
                end do
              end do
            end do

            do k = 1,kmax-1
              do j = 1,jmax-1
                do i = 1,imax-1

                  do ni = 0,2*nharms
                    dtau     = dtvolt(i,j,k,ni) * vol(i,j,k)
                    ome      = work1(i,j,k,7,ni) / q(i,j,k,1,ni)
                    dps      = delpls(i,j,k,ni)
                    do nj = 0,2*nharms
                      hbstb(ni,nj) = ID1(ni,nj) * &
                        (1 + ar*dtau * (    dps +    bstrkw*ome) ) + &
                        ar*omega*dtau*dhb(ni,nj)
                    end do
                  end do

                  call DGETRF(mm,nn,hbstb,lda,ipiv,info)
!dbg              write(*,*) 'info DGETRF is:',info

!---------------- store RHS in work2.
                  do ni = 0,2*nharms
                    work2(ni,1) = qhb0(i,j,k,6,ni)
                  end do

                  call DGETRS(trans,nn,nrhs,hbstb,lda,ipiv, &
                              work2,ldb,info)
!dbg              write(*,*) 'info DGETRS is:',info

!---------------- store update rho*K in q(:,:,6,:).
                  do ni = 0,2*nharms
                    q(i,j,k,6,ni) = work2(ni,1)
                  end do

                end do
              end do
            end do

          end if

!-------- UPDATE Euler / NS variables

!-------- calculate RHS of eq. (22) of gt-2011-45303
          do ni = 0,2*nharms
            do ipde = 1,5
              do k = 1,kmax-1
                do j = 1,jmax-1
                  do i = 1,imax-1
                    qhb0(i,j,k,ipde,ni) = qold(i,j,k,ipde,ni) + &
                      (qhb0(i,j,k,ipde,ni) - res(i,j,k,ipde,ni)) * &
                      ar*dtvol(i,j,k,ni)
                  end do
                end do
              end do
            end do
          end do

          trans = 'N'
          mm = 2*nharms + 1
          nn = 2*nharms + 1
          lda = mm
          ldb = nn
          nrhs = 5

          do k = 1,kmax-1
            do j = 1,jmax-1
              do i = 1,imax-1

                do ni = 0,2*nharms
                  dtau = dtvol(i,j,k,ni) * vol(i,j,k)
                  do nj = 0,2*nharms
                    hbstb(ni,nj) = ID1(ni,nj) + &
                                     ar*omega*dtau*dhb(ni,nj)
                  end do
                end do

                call DGETRF(mm,nn,hbstb,lda,ipiv,info)
!dbg            write(*,*) 'info DGETRF is:',info

!               transpose RHS matrix, to have pattern [(2nharms+1)xnpde]
                do ni = 0,2*nharms
                  do ipde = 1,5
                    work2(ni,ipde) = qhb0(i,j,k,ipde,ni)
                  end do
                end do

!dbg            do k = 0,2*nharms
!dbg              if (ipiv(k).ne.(k+1)) then
!dbg                write(3,*) 'using the pivoting'
!dbg              end if
!dbg            end do

                call DGETRS(trans,nn,nrhs,hbstb,lda,ipiv, &
                            work2,ldb,info)
!dbg            write(*,*) 'info DGETRS is:',info

!-------------- copy new solution back to q
                do ni = 0,2*nharms
                  do ipde = 1,5
                    q(i,j,k,ipde,ni) = work2(ni,ipde)
                  end do   
                end do

              end do
            end do
          end do

        end if

      end if

      deallocate(hbstb,ID1,qhb,work2,ipiv)

      return
      end

!-----------------------------------------------------------------------
      subroutine prim2cons(nl,q,work)
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,iblock,imax,jmax,kmax,iq
      real (kind=cosa_real) q(*),work(*)

      do iblock = 1,mynblocks
        imax   = i_imax     (iblock,nl)
        jmax   = j_jmax     (iblock,nl)
        kmax   = k_kmax     (iblock,nl)
        iq     = 1 + off_p3 (iblock,nl) * npde * dim5
        call prim2cons_b(q(iq),work(iq),npde,imax,jmax,kmax,nharms)
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine prim2cons_b(q,work,npde,imax,jmax,kmax,nharms)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) npde,imax,jmax,kmax,nharms
      integer(kind=cosa_int) i,j,k,n,ipde
      real (kind=cosa_real) &
           q   (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           work(-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms)
      real (kind=cosa_real) rho,rhou,rhov,rhow,p,tke

      do n=0,2*nharms

        if (kom.or.kom_bsl.or.kom_sst) then

          do k = 1,kmax-1
            do j = 1,jmax-1
              do i = -1,imax+1
                work(i,j,k,1,n) = q(i,j,k,6,n)
              end do
            end do
          end do
          do k = 1,kmax-1
            do j = -1,0
              do i = 1,imax-1
                work(i,j,k,1,n) = q(i,j,k,6,n)
              end do
            end do
          end do
          do k = 1,kmax-1
            do j = jmax,jmax+1
              do i = 1,imax-1
                work(i,j,k,1,n) = q(i,j,k,6,n)
              end do
            end do
          end do
          do k = -1,0
            do j = 1,jmax-1
              do i = 1,imax-1
                work(i,j,k,1,n) = q(i,j,k,6,n)
              end do
            end do
          end do
          do k = kmax,kmax+1
            do j = 1,jmax-1
              do i = 1,imax-1
                work(i,j,k,1,n) = q(i,j,k,6,n)
              end do
            end do
          end do

        else

          do k = 1,kmax-1
            do j = 1,jmax-1
              do i = -1,imax+1
                work(i,j,k,1,n) = 0.d0
              end do
            end do
          end do
          do k = 1,kmax-1
            do j = -1,0
              do i = 1,imax-1
                work(i,j,k,1,n) = 0.d0
              end do
            end do
          end do
          do k = 1,kmax-1
            do j = jmax,jmax+1
              do i = 1,imax-1
                work(i,j,k,1,n) = 0.d0
              end do
            end do
          end do
          do k = -1,0
            do j = 1,jmax-1
              do i = 1,imax-1
                work(i,j,k,1,n) = 0.d0
              end do
            end do
          end do
          do k = kmax,kmax+1
            do j = 1,jmax-1
              do i = 1,imax-1
                work(i,j,k,1,n) = 0.d0
              end do
            end do
          end do

        end if

        do ipde=2,npde
          if (ipde.ne.5) then
            do k = 1,kmax-1
              do j = 1,jmax-1
                do i = -1,imax+1
                  q(i,j,k,ipde,n) = q(i,j,k,ipde,n)*q(i,j,k,1,n)
                end do
              end do
            end do
            do k = 1,kmax-1
              do j = -1,0
                do i = 1,imax-1
                  q(i,j,k,ipde,n) = q(i,j,k,ipde,n)*q(i,j,k,1,n)
                end do
              end do
            end do
            do k = 1,kmax-1
              do j = jmax,jmax+1
                do i = 1,imax-1
                  q(i,j,k,ipde,n) = q(i,j,k,ipde,n)*q(i,j,k,1,n)
                end do
              end do
            end do
            do k = -1,0
              do j = 1,jmax-1
                do i = 1,imax-1
                  q(i,j,k,ipde,n) = q(i,j,k,ipde,n)*q(i,j,k,1,n)
                end do
              end do
            end do
            do k = kmax,kmax+1
              do j = 1,jmax-1
                do i = 1,imax-1
                  q(i,j,k,ipde,n) = q(i,j,k,ipde,n)*q(i,j,k,1,n)
                end do
              end do
            end do
          end if
        end do

        do k = 1,kmax-1
          do j = 1,jmax-1
            do i = -1,imax+1
              rho  = q(i,j,k,1,n)
              rhou = q(i,j,k,2,n)
              rhov = q(i,j,k,3,n)
              rhow = q(i,j,k,4,n)
              p    = q(i,j,k,5,n)
              tke  = work(i,j,k,1,n)
              q(i,j,k,5,n) = p/(gamma-1) + &
                         (rhou*rhou+rhov*rhov+rhow*rhow)/2/rho + rho*tke
            end do
          end do
        end do
        do k = 1,kmax-1
          do j = -1,0
            do i = 1,imax-1
              rho  = q(i,j,k,1,n)
              rhou = q(i,j,k,2,n)
              rhov = q(i,j,k,3,n)
              rhow = q(i,j,k,4,n)
              p    = q(i,j,k,5,n)
              tke  = work(i,j,k,1,n)
              q(i,j,k,5,n) = p/(gamma-1) + &
                         (rhou*rhou+rhov*rhov+rhow*rhow)/2/rho + rho*tke
            end do
          end do
        end do
        do k = 1,kmax-1
          do j = jmax,jmax+1
            do i = 1,imax-1
              rho  = q(i,j,k,1,n)
              rhou = q(i,j,k,2,n)
              rhov = q(i,j,k,3,n)
              rhow = q(i,j,k,4,n)
              p    = q(i,j,k,5,n)
              tke  = work(i,j,k,1,n)
              q(i,j,k,5,n) = p/(gamma-1) + &
                         (rhou*rhou+rhov*rhov+rhow*rhow)/2/rho + rho*tke
            end do
          end do
        end do
        do k = -1,0
          do j = 1,jmax-1
            do i = 1,imax-1
              rho  = q(i,j,k,1,n)
              rhou = q(i,j,k,2,n)
              rhov = q(i,j,k,3,n)
              rhow = q(i,j,k,4,n)
              p    = q(i,j,k,5,n)
              tke  = work(i,j,k,1,n)
              q(i,j,k,5,n) = p/(gamma-1) + &
                         (rhou*rhou+rhov*rhov+rhow*rhow)/2/rho + rho*tke
            end do
          end do
        end do
        do k = kmax,kmax+1
          do j = 1,jmax-1
            do i = 1,imax-1
              rho  = q(i,j,k,1,n)
              rhou = q(i,j,k,2,n)
              rhov = q(i,j,k,3,n)
              rhow = q(i,j,k,4,n)
              p    = q(i,j,k,5,n)
              tke  = work(i,j,k,1,n)
              q(i,j,k,5,n) = p/(gamma-1) + &
                         (rhou*rhou+rhov*rhov+rhow*rhow)/2/rho + rho*tke
            end do
          end do
        end do

      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine cons2prim(nl,xtnd,q,work)
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,imax,jmax,kmax
      integer(kind=cosa_int) xtnd,iblock,iq
      real (kind=cosa_real) q(*),work(*)

      do iblock = 1,mynblocks
        imax   = i_imax     (iblock,nl)
        jmax   = j_jmax     (iblock,nl)
        kmax   = k_kmax     (iblock,nl)
        iq     = 1 + off_p3 (iblock,nl) * npde * dim5
        call cons2prim_b(q(iq),work(iq),xtnd,npde,imax,jmax,kmax,nharms)
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine cons2prim_b(q,work,xtnd,npde,imax,jmax,kmax,nharms)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) npde,imax,jmax,kmax,nharms
      integer(kind=cosa_int) i,istr,iend,j,jstr,jend,k,kstr,kend,n,ipde,xtnd
      real (kind=cosa_real) &
           q   (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           work(-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms)
      real (kind=cosa_real) rho,u,v,w,rhoe,tke

      istr = 1 - xtnd
      jstr = 1 - xtnd
      kstr = 1 - xtnd
      iend = imax - 1 + xtnd
      jend = jmax - 1 + xtnd
      kend = kmax - 1 + xtnd

      do n=0,2*nharms

        do ipde=2,npde
          if (ipde.ne.5) then
            do k = kstr,kend
              do j = jstr,jend
                do i = istr,iend
                  q(i,j,k,ipde,n) = q(i,j,k,ipde,n)/q(i,j,k,1,n)
                end do
              end do
            end do
          end if
        end do

        if (kom.or.kom_bsl.or.kom_sst) then

          do k = kstr,kend
            do j = jstr,jend
              do i = istr,iend
                work(i,j,k,1,n) = q(i,j,k,6,n)
              end do
            end do
          end do

        else

          do k = kstr,kend
            do j = jstr,jend
              do i = istr,iend
                work(i,j,k,1,n) = 0.d0
              end do
            end do
          end do

        end if

        do k = kstr,kend
          do j = jstr,jend
            do i = istr,iend
              rho  = q(i,j,k,1,n)
              u    = q(i,j,k,2,n)
              v    = q(i,j,k,3,n)
              w    = q(i,j,k,4,n)
              rhoe = q(i,j,k,5,n)
              tke  = work(i,j,k,1,n)
              q(i,j,k,5,n) = (gamma-1)* &
                (rhoe - rho*(u*u+v*v+w*w)/2 - rho*tke)
            end do
          end do
        end do

      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine fixq(nl,q)
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none
!
      integer(kind=cosa_int) nl,iblock,imax,jmax,kmax,iq
      real (kind=cosa_real) ar
      real (kind=cosa_real) q(*)
!
      do iblock = 1,mynblocks
        imax   = i_imax     (iblock,nl)
        jmax   = j_jmax     (iblock,nl)
        kmax   = k_kmax     (iblock,nl)
        iq     = 1 + off_p3 (iblock,nl) * npde * dim5
        call fixq_b(q(iq),imax,jmax,kmax,npde,nharms,mgit,nl)
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine fixq_b(q,imax,jmax,kmax,npde,nharms,ic,nl)
!-----------------------------------------------------------------------
!     positivity preservation for turbulent flows
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none
!
      integer(kind=cosa_int) imax,jmax,kmax,npde,nharms
      integer(kind=cosa_int) i,j,k,ipde,n,ic,nl
      real (kind=cosa_real) q(-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms)
      real (kind=cosa_real) rho,tke,ome

      do n = 0,2*nharms

        do k = 1,kmax-1
          do j = 1,jmax-1
            do i = 1,imax-1
              rho     = q(i,j,k,1,n)
              tke     = q(i,j,k,6,n) / rho
              ome     = q(i,j,k,7,n) / rho
              if (tke.lt.qmin(6)) then
                tke = qmin(6)
                q(i,j,k,6,n) = rho*tke
                if (debug) then
                  write(2,100) nl,i,j,k,ic
                  write(*,100) nl,i,j,k,ic
                end if
              end if
              if (ome.lt.qmin(7)) then
                ome = qmin(7)
                q(i,j,k,7,n) = rho*ome
                if (debug) then
                  write(2,101) nl,i,j,k,ic
                  write(*,101) nl,i,j,k,ic
                end if
              end if
            end do
          end do
        end do

      end do

 100  format(///,2x,'negative TKE on grid level',i2,' at i=',i3,', &
      j=',i3,' k=',i3,//,2x,'after',i6,'  MG iterations. Enforcing posit &
     &ivity.')
 101  format(///,2x,'negative OME on grid level',i2,' at i=',i3,', &
      j=',i3,' k=',i3,//,2x,'after',i6,'  MG iterations. Enforcing posit &
     &ivity.')

      return      
      end

!-----------------------------------------------------------------------
      subroutine tur_vis(q,mut,nl)
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables
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
        call tur_vis_b(q(iq),mut(imut),npde,imax,jmax,kmax,nharms)
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine tur_vis_b(q,mut,npde,imax,jmax,kmax,nharms)
!-----------------------------------------------------------------------
!     MSC, 09 August 2015: 
!     the loop ends below are such that:
!     a. first plane of regular auxiliary cells is computed correctly
!        since call bc always preceeds call tur_vis.
!     b. corner muT is zero as no value is ever assigned to corner q(6).
!     c. edge muT is zero if call tur_vis preceeds call q_edges. This
!        does not matter because edge muT is never used.
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) npde,imax,jmax,kmax,nharms
      integer(kind=cosa_int) i,j,k,n,ipde
      real (kind=cosa_real) &
           q  (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           mut(-1:imax+1,-1:jmax+1,-1:kmax+1     ,0:2*nharms)
      real (kind=cosa_real) rho,tke,ome

      do n=0,2*nharms

        do k = 0,kmax
          do j = 0,jmax
            do i = 0,imax
              rho = q(i,j,k,1,n)
              tke = q(i,j,k,6,n)
              ome = q(i,j,k,7,n)
              mut(i,j,k,n) = cmu * rho * tke / (ome+epst) * &
                             reyno / machfs
            end do
          end do
        end do

      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine tur_vis_sst(q,mut,dqxi,dqeta,dqzeta,xideri,etaderi, &
        zetaderi,xiderj,etaderj,zetaderj,xiderk,etaderk,zetaderk,dist, &
        nl)
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,iblock,imax,jmax,kmax,iq,imut,ivmt,idist
      real (kind=cosa_real) q(*),mut(*),dqxi(*),dqeta(*),dqzeta(*),xideri(*), &
           etaderi(*),zetaderi(*),xiderj(*),etaderj(*),zetaderj(*), &
           xiderk(*),etaderk(*),zetaderk(*),dist(*)

      do iblock = 1,mynblocks
        imax   = i_imax     (iblock,nl)
        jmax   = j_jmax     (iblock,nl)
        kmax   = k_kmax     (iblock,nl)
        iq     = 1 + off_p3 (iblock,nl) * npde * dim5
        imut   = 1 + off_p3 (iblock,nl) *        dim5
        ivmt   = 1 + off_0  (iblock,nl) * 3    * dim5h
        idist  = 1 + off_p1 (iblock,nl)
        call tur_vis_sst_b(q(iq),mut(imut),dqxi(iq),dqeta(iq), &
          dqzeta(iq),xideri(ivmt),etaderi(ivmt),zetaderi(ivmt), &
          xiderj(ivmt),etaderj(ivmt),zetaderj(ivmt),xiderk(ivmt), &
          etaderk(ivmt),zetaderk(ivmt),dist(idist),npde,imax,jmax,kmax, &
          nharms)
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine tur_vis_sst_b(q,mut,dqxi,dqeta,dqzeta,xideri,etaderi, &
        zetaderi,xiderj,etaderj,zetaderj,xiderk,etaderk,zetaderk,dist, &
        npde,imax,jmax,kmax,nharms)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) npde,imax,jmax,kmax,nharms
      integer(kind=cosa_int) i,j,k,n,ipde,nh,i1,j1,k1
      real (kind=cosa_real) &
        q       (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
        mut     (-1:imax+1,-1:jmax+1,-1:kmax+1     ,0:2*nharms), &
        dqxi    (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
        dqeta   (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
        dqzeta  (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
        xideri  (   3,  imax  ,  jmax  ,  kmax  ,0:2*nharms*hbmove), &
        etaderi (   3,  imax  ,  jmax  ,  kmax  ,0:2*nharms*hbmove), &
        zetaderi(   3,  imax  ,  jmax  ,  kmax  ,0:2*nharms*hbmove), &
        xiderj  (   3,  imax  ,  jmax  ,  kmax  ,0:2*nharms*hbmove), &
        etaderj (   3,  imax  ,  jmax  ,  kmax  ,0:2*nharms*hbmove), &
        zetaderj(   3,  imax  ,  jmax  ,  kmax  ,0:2*nharms*hbmove), &
        xiderk  (   3,  imax  ,  jmax  ,  kmax  ,0:2*nharms*hbmove), &
        etaderk (   3,  imax  ,  jmax  ,  kmax  ,0:2*nharms*hbmove), &
        zetaderk(   3,  imax  ,  jmax  ,  kmax  ,0:2*nharms*hbmove), &
        dist    (     0:imax  ,0:jmax  ,0:kmax)

      real (kind=cosa_real) rho,tke,ome,xix,xiy,xiz,etax,etay,etaz,zetax,zetay, &
        zetaz,uxi,ueta,uzeta,vxi,veta,vzeta,wxi,weta,wzeta,dudx,dudy,dudz, &
        dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,sxx,sxy,sxz,syy,syz,szz,s,arg2,f2, &
        mu,vort,p,t,fact

      do n=0,2*nharms

        nh = n*hbmove

        do k = 1,kmax-1
          do j = 1,jmax-1
            do i = 1,imax-1

              p   = (gamma-1) * &
                (q(i,j,k,5,n) - &
                 (q(i,j,k,2,n)**2 + q(i,j,k,3,n)**2 + q(i,j,k,4,n)**2)/ &
                 q(i,j,k,1,n)/2 - q(i,j,k,6,n))
              rho = q(i,j,k,1,n)
              tke = q(i,j,k,6,n) / rho
              ome = q(i,j,k,7,n) / rho
              t   = gamma*p/rho
              mu  = (stemp+1)/(t+stemp) * (t**1.5d0)
              
              xix   = (xideri  (1,i+1,j  ,k  ,nh) + &
                       xideri  (1,i  ,j  ,k  ,nh) ) / 2
              xiy   = (xiderj  (2,i  ,j+1,k  ,nh) + &
                       xiderj  (2,i  ,j  ,k  ,nh) ) / 2
              xiz   = (xiderk  (3,i  ,j  ,k+1,nh) + &
                       xiderk  (3,i  ,j  ,k  ,nh) ) / 2
              etax  = (etaderi (1,i+1,j  ,k  ,nh) + &
                       etaderi (1,i  ,j  ,k  ,nh) ) / 2
              etay  = (etaderj (2,i  ,j+1,k  ,nh) + &
                       etaderj (2,i  ,j  ,k  ,nh) ) / 2
              etaz  = (etaderk (3,i  ,j  ,k+1,nh) + &
                       etaderk (3,i  ,j  ,k  ,nh) ) / 2
              zetax = (zetaderi(1,i+1,j  ,k  ,nh) + &
                       zetaderi(1,i  ,j  ,k  ,nh) ) / 2
              zetay = (zetaderj(2,i  ,j+1,k  ,nh) + &
                       zetaderj(2,i  ,j  ,k  ,nh) ) / 2
              zetaz = (zetaderk(3,i  ,j  ,k+1,nh) + &
                       zetaderk(3,i  ,j  ,k  ,nh) ) / 2

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

              vort = sqrt( (dwdy-dvdz)**2 + (dudz-dwdx)**2 + &
                           (dvdx-dudy)**2 )

              sxx = dudx
              sxy = (dudy + dvdx)/2
              sxz = (dudz + dwdx)/2
              syy = dvdy
              syz = (dvdz + dwdy)/2
              szz = dwdz

              s   = sxx**2+syy**2+szz**2 + 2*(sxy**2+sxz**2+syz**2)

              arg2 = dmax1(2*dsqrt(tke)/0.09d0/ome/(dist(i,j,k)+epst), &
                           500*mu/rho/ome/(dist(i,j,k)+epst)**2* &
                           machfs/reyno)
              arg2   = arg2**2
              arg2   = dmin1(arg2,1.d2)
              f2     = (dexp(2*arg2)-1)/(dexp(2*arg2)+1)
              
!tmp          fact = vort
!tmp          fact = sqrt(s)
              fact = sqrt(2*s)

              mut(i,j,k,n) = a1 * rho * tke / dmax1(a1*ome,fact*f2) * &
                      reyno / machfs

            end do
          end do
        end do

        do k = 0,kmax,kmax
          do j = 1,jmax-1
            do i = 1,imax-1
              k1 = k + (-1)**(k/kmax)
              mut(i,j,k,n) =  mut(i,j,k1,n)
            end do
          end do
        end do
        do k = 1,kmax-1
          do j = 0,jmax,jmax
            do i = 1,imax-1
              j1 = j + (-1)**(j/jmax)
              mut(i,j,k,n) =  mut(i,j1,k,n)
            end do
          end do
        end do
        do k = 1,kmax-1
          do j = 1,jmax-1
            do i = 0,imax,imax
              i1 = i + (-1)**(i/imax)
              mut(i,j,k,n) =  mut(i1,j,k,n)
            end do
          end do
        end do

      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine lim_res_ome(nl,res,q,qold,vol,dtvolt,prec,ipkpi,work1, &
                             delpls,prd,ttrm,lr,ar)
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none
!
      integer(kind=cosa_int) nl,imax,jmax,kmax
      integer(kind=cosa_int) iblock,iq,idtv,ivol,iprec,idlpl,ittrm
      real (kind=cosa_real) ar
      real (kind=cosa_real) q(*),qold(*),res(*),vol(*),dtvolt(*),prec(*), &
           ipkpi(*),work1(*),delpls(*),prd(*),ttrm(*),lr(*)
!
      if (harbal) then

        do iblock = 1,mynblocks
          imax   = i_imax     (iblock,nl)
          jmax   = j_jmax     (iblock,nl)
          kmax   = k_kmax     (iblock,nl)
          iq     = 1 + off_p3 (iblock,nl) * npde * dim5
          idlpl  = 1 + off_p3 (iblock,nl) *        dim5
          idtv   = 1 + off_p1 (iblock,nl)        * dim5
          ivol   = 1 + off_p1 (iblock,nl)
          iprec  = 1 + off_m1 (iblock,nl) * npde * dim5 * npde
          ittrm  = 1 + off_p3 (iblock,nl) * 3    * dim5
          call lim_res_ome_bhb(res(iq),q(iq),qold(iq),vol(ivol), &
                 dtvolt(idtv),prec(iprec),ipkpi(iprec),work1(iq), &
                 delpls(idlpl),prd(idlpl),ttrm(ittrm),lr(iq), &
                 imax,jmax,kmax,npde,nharms)
        end do

      else

        do iblock = 1,mynblocks
          imax   = i_imax     (iblock,nl)
          jmax   = j_jmax     (iblock,nl)
          kmax   = k_kmax     (iblock,nl)
          iq     = 1 + off_p3 (iblock,nl) * npde * dim5
          idlpl  = 1 + off_p3 (iblock,nl) *        dim5
          idtv   = 1 + off_p1 (iblock,nl)        * dim5
          ivol   = 1 + off_p1 (iblock,nl)
          iprec  = 1 + off_m1 (iblock,nl) * npde * dim5 * npde
          ittrm  = 1 + off_p3 (iblock,nl) * 3    * dim5
          call lim_res_ome_b(res(iq),q(iq),  vol(ivol), &
                 dtvolt(idtv),prec(iprec),ipkpi(iprec),work1(iq), &
                 delpls(idlpl),prd(idlpl),ttrm(ittrm),lr(iq),ar, &
                 imax,jmax,kmax,npde)
        end do

      end if

      return
      end

!-----------------------------------------------------------------------
      subroutine lim_res_ome_b(res,q,vol,dtvolt,prec,ipkpi,work1, &
                          delpls,prd,ttrm,lr,ar,imax,jmax,kmax,npde)
!-----------------------------------------------------------------------
!     MSC, 30 April 2012: 
!     this routine is called ONLY if nlevel.gt.1. 
!     Hence, the residual to which the limitation of this routine is 
!     applied is always assumed to become zero at convergence,
!     regardless of whether IRS is used or not.
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none
!
      integer(kind=cosa_int) imax,jmax,kmax,npde
      integer(kind=cosa_int) im,jm,km,i,j,k,ipde,ipde1
      real (kind=cosa_real) fact,fact1,term1,dtau,rho,ome,den1,a(2,2),tmp(7,7), &
           qtmp(7),den2,rhomega,drhomega,rhomegamin,res_ome,gav,bav, &
           somma,ar
      real (kind=cosa_real) &
           q     (-1:imax+1,-1:jmax+1,-1:kmax+1,npde), &
           delpls(-1:imax+1,-1:jmax+1,-1:kmax+1     ), &
           prd   (-1:imax+1,-1:jmax+1,-1:kmax+1     ), &
           res   (-1:imax+1,-1:jmax+1,-1:kmax+1,npde), &
           work1 (-1:imax+1,-1:jmax+1,-1:kmax+1,npde), &
           ttrm  (-1:imax+1,-1:jmax+1,-1:kmax+1,3   ), &
           lr    (-1:imax+1,-1:jmax+1,-1:kmax+1,npde), &
           dtvolt( 0:imax  , 0:jmax  , 0:kmax       ), &
           vol   ( 0:imax  ,0:jmax   , 0:kmax       ), &
           prec  (   imax-1,   jmax-1,   kmax-1,npde,npde), &
           ipkpi (   imax-1,   jmax-1,   kmax-1,npde,npde)
!
      im = imax-1
      jm = jmax-1
      km = kmax-1

!---- update solution

      if (lomach) then

        if (.not.rkimp) then

          do k = 1,km
            do j = 1,jm
              do i = 1,im

                gav        = ttrm(i,j,k,2)
                bav        = ttrm(i,j,k,3)

                a(1,1) = delpls(i,j,k) + &
                         bstrkw * q(i,j,k,7) / q(i,j,k,1)
                a(1,2) = bstrkw * q(i,j,k,6) / q(i,j,k,1)
                a(2,1) = 0.d0
                a(2,2) = gav * delpls(i,j,k) + &
                         2 * bav * q(i,j,k,7) / q(i,j,k,1)

                do ipde = 1,npde
                  tmp(ipde,6) = prec(i,j,k,ipde,6) * a(1,1)
                  tmp(ipde,7) = prec(i,j,k,ipde,6) * a(1,2) + &
                                prec(i,j,k,ipde,7) * a(2,2)
                  work1(i,j,k,ipde) = &
                    tmp(ipde,6) * q(i,j,k,6) + tmp(ipde,7) * q(i,j,k,7)
                end do

!------------ update OMEGA residual
              ipde=7

              dtau       = dtvolt(i,j,k) * vol(i,j,k)

              rhomega = 0
              do ipde1 = 1,npde
                rhomega = rhomega + &
                  ipkpi(i,j,k,ipde,ipde1) * &
                  (q(i,j,k,ipde1) + ar * dtau * work1(i,j,k,ipde1) - &
                                ar * dtvolt(i,j,k) * res(i,j,k,ipde1))
              end do
              rhomegamin = gav*cmu*q(i,j,k,1)*sqrt(prd(i,j,k))
              if (rhomega.lt.rhomegamin) then
                qtmp(ipde) = rhomegamin
                do ipde1 = 1,6
                  qtmp(ipde1) = q(i,j,k,ipde1)
                end do
                somma = 0.d0
                do ipde1 = 1,npde
                  somma = somma + lr(i,j,k,ipde1) * qtmp(ipde1)
                end do
                res(i,j,k,ipde) = &
                  (-somma +q(i,j,k,ipde) + ar* dtau*work1(i,j,k,ipde) )/ &
                  dtvolt(i,j,k) / ar
              end if

              end do
            end do
          end do


        else if (rkimp) then

          fact = 1.5d0
          if ((itime.eq.1).and.(.not.dual_prop_start)) then
            fact = 1.0d0
          end if

          if ((irsop.gt.0).or.((irsop.eq.0).and.(nlevel.gt.1)))  then

!-----------------------------------------------------------------------
!           TIME-DEPENDENT WITH SIRK and WITH MG AND/OR IRS
!           (the conditioned above is the only possible with rkim)
!-----------------------------------------------------------------------

            do k = 1,km
              do j = 1,jm
                do i = 1,im

                gav        = ttrm(i,j,k,2)
                bav        = ttrm(i,j,k,3)
                dtau       = dtvolt(i,j,k) * vol(i,j,k)

                a(1,1) = delpls(i,j,k) + bstrkw * q(i,j,k,7) / q(i,j,k,1)
                a(1,2) = bstrkw * q(i,j,k,6) / q(i,j,k,1)
                a(2,1) = 0.d0
                a(2,2) = gav * delpls(i,j,k) &
                             + 2 * bav * q(i,j,k,7) / q(i,j,k,1)

                do ipde = 1,npde
                  do ipde1=1,npde
                    tmp(ipde,ipde1) = prec(i,j,k,ipde,ipde1) * fact/dt
                  end do
                  tmp(ipde,6) = tmp(ipde,6) + prec(i,j,k,ipde,6) * a(1,1)
                  tmp(ipde,7) = tmp(ipde,7) + &
                                prec(i,j,k,ipde,6) * a(1,2) + &
                                prec(i,j,k,ipde,7) * a(2,2)
                  work1(i,j,k,ipde) = 0.d0
                  do ipde1=1,npde
                    work1(i,j,k,ipde) = work1(i,j,k,ipde) + &
                      tmp(ipde,ipde1) * q(i,j,k,ipde1)
                  end do
                end do

!-------------- update OMEGA residual
                ipde=7

                rhomega = 0
                do ipde1 = 1,npde
                  rhomega = rhomega + &
                    ipkpi(i,j,k,ipde,ipde1) * &
                    (q(i,j,k,ipde1) + ar * dtau * work1(i,j,k,ipde1) - &
                                  ar * dtvolt(i,j,k) * res(i,j,k,ipde1))
                end do
                rhomegamin = gav*cmu*q(i,j,k,1)*sqrt(prd(i,j,k))
                if (rhomega.lt.rhomegamin) then
                  qtmp(ipde) = rhomegamin
                  do ipde1 = 1,6
                    qtmp(ipde1) = q(i,j,k,ipde1)
                  end do
                  somma = 0.d0
                  do ipde1 = 1,npde
                    somma = somma + lr(i,j,k,ipde1) * qtmp(ipde1)
                  end do
                  res(i,j,k,ipde) = &
                    (-somma +q(i,j,k,ipde) + ar * dtau*work1(i,j,k,ipde) ) &
                    / dtvolt(i,j,k) / ar
                end if

                end do
              end do
            end do

          end if

        end if

      else if (.not.lomach) then

        if (.not.rkimp) then

          do k = 1,km
            do j = 1,jm
              do i = 1,im
                gav        = ttrm(i,j,k,2)
                bav        = ttrm(i,j,k,3)
                dtau       = dtvolt(i,j,k) * vol(i,j,k)

                rho        = q(i,j,k,1)
                ome        = q(i,j,k,7) / q(i,j,k,1)
                den1       = 1 + ar * dtau*(gav*delpls(i,j,k) + 2*bav*ome)
                rhomega    = q(i,j,k,7) &
                            - ar * res(i,j,k,7)*dtvolt(i,j,k)/den1

                rhomegamin = gav*cmu*rho*sqrt(prd(i,j,k))
                if (rhomega.lt.rhomegamin) then
                  res(i,j,k,7) = den1*(q(i,j,k,7)-rhomegamin) / &
                             dtvolt(i,j,k) / ar
                end if

              end do
            end do
          end do

        else if (rkimp) then

          fact = 1.5d0
          if ((itime.eq.1).and.(.not.dual_prop_start)) then
            fact = 1.0d0
          end if

          if ((irsop.gt.0).or.((irsop.eq.0).and.(nlevel.gt.1))) then

!-----------------------------------------------------------------------
!           TIME-DEPENDENT WITH SIRK and WITH MG AND/OR IRS
!           (the conditioned above is the only possible with rkim)
!-----------------------------------------------------------------------

            do k = 1,km
              do j = 1,jm
                do i = 1,im
      
                gav        = ttrm(i,j,k,2)
                bav        = ttrm(i,j,k,3)
                dtau       = dtvolt(i,j,k) * vol(i,j,k)

                rho        = q(i,j,k,1)
                ome        = q(i,j,k,7) / q(i,j,k,1)
                den1       = 1 + dtau*(gav*delpls(i,j,k) + 2*bav*ome)
                den2       = ar*fact*dtau/dt

                rhomega    = q(i,j,k,7) - &
                             ar*res(i,j,k,7)*dtvolt(i,j,k)/(den1+den2)

                rhomegamin = gav*cmu*rho*sqrt(prd(i,j,k))
                if (rhomega.lt.rhomegamin) then
                  res(i,j,k,7) = &
                    ((den1+den2)*(q(i,j,k,7)-rhomegamin)) / &
                     dtvolt(i,j,k) / ar
                end if

                end do
              end do
            end do

          end if 

        end if

      end if

      return
      end 

!-----------------------------------------------------------------------
      subroutine lim_res_ome_bhb(res,q,qold,vol,dtvolt,prec,ipkpi,work1, &
                          delpls,prd,ttrm,lr,imax,jmax,kmax,npde,nharms)
!-----------------------------------------------------------------------
!     MSC, 30 April 2012: this routine is called ONLY if nlevel.gt.1. 
!     Hence, the residual to which the limitation of this routine is 
!     applied is always assumed to become zero at convergence,
!     regardless of whether IRS is used or not.
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none
!
      integer(kind=cosa_int) imax,jmax,kmax,npde,nharms
      integer(kind=cosa_int) im,jm,km,i,j,k,ipde,n,ni,nj,lda,info,mm,nn,ldb,nrhs
      real (kind=cosa_real) term1,dtau,rho,ome,den1,den2, rhomega,drhomega, &
           rhomegamin,rhomeganew,res_ome,dps,gav,bav
      real (kind=cosa_real) &
           q     (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           qold  (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           res   (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           work1 (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           ttrm  (-1:imax+1,-1:jmax+1,-1:kmax+1,3   ,0:2*nharms), &
           lr    (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           delpls(-1:imax+1,-1:jmax+1,-1:kmax+1     ,0:2*nharms), &
           prd   (-1:imax+1,-1:jmax+1,-1:kmax+1     ,0:2*nharms), &
           dtvolt( 0:imax  , 0:jmax  , 0:kmax       ,0:2*nharms), &
           vol   ( 0:imax  , 0:jmax  , 0:kmax                  ), &
           prec  (   imax-1,   jmax-1,   kmax-1,npde,npde,0:2*nharms), &
           ipkpi (   imax-1,   jmax-1,   kmax-1,npde,npde,0:2*nharms)
      integer(kind=cosa_int), allocatable:: ipiv(:)
      real (kind=cosa_real), allocatable::hbstb(:,:),work2(:),id1(:,:)
      character*1 trans
!
      im = imax-1
      jm = jmax-1
      km = kmax-1

!---- update solution

      if (lomach) then

        if (.not.rkimp) then


        else if (rkimp) then


        end if

      else if (.not.lomach) then

        allocate(  hbstb(0:2*nharms,0:2*nharms),work2(0:2*nharms), &
                   id1(0:2*nharms,0:2*nharms),ipiv(0:2*nharms) )

        if (.not.rkimp) then

          do n = 0,2*nharms
            do k = 1,km
              do j = 1,jm
                do i = 1,im
                  gav        = ttrm(i,j,k,2,n)
                  bav        = ttrm(i,j,k,3,n)
                  rho        = q(i,j,k,1,n)
                  ome        = q(i,j,k,7,n) / q(i,j,k,1,n)
                  dtau       = dtvolt(i,j,k,n) * vol(i,j,k)
                  rhomega    = qold(i,j,k,7,n) - &
                    res(i,j,k,7,n) * dtvolt(i,j,k,n) / &
                    (1 + dtau*(gav*delpls(i,j,k,n) + 2*bav*ome))

                  rhomegamin = gav*cmu*rho*sqrt(prd(i,j,k,n))

                  res_ome = res(i,j,k,7,n) - &
                    (1 + dtau*(gav*delpls(i,j,k,n) + 2*bav*ome)) * &
                    max(0.d0,(rhomegamin-rhomega)/dtvolt(i,j,k,n))
                  res(i,j,k,7,n) = res_ome
                end do
              end do
            end do
          end do

        else if (rkimp) then

          do nj=0,2*nharms
            do ni=0,2*nharms
              if (ni.eq.nj) then
                id1(ni,nj) = 1.d0
              else
                id1(ni,nj) = 0.d0
              end if
            end do
          end do

          trans = 'N'
          mm = 2*nharms + 1
          nn = 2*nharms + 1
          lda = mm
          ldb = nn
          nrhs = 1

          do k = 1,kmax-1
            do j = 1,jmax-1
              do i = 1,imax-1

                do ni = 0,2*nharms
                  gav      = ttrm(i,j,k,2,ni)
                  bav      = ttrm(i,j,k,3,ni)
                  dtau     = dtvolt(i,j,k,ni) * vol(i,j,k)
                  ome      = q(i,j,k,7,ni) / q(i,j,k,1,ni)
                  dps      = delpls(i,j,k,ni)
                  do nj = 0,2*nharms
                    hbstb(ni,nj) = id1(ni,nj) * &
                          (1 + dtau * (gav*dps +  2*bav*ome) ) + &
                          omega*dtau*dhb(ni,nj)
                  end do
                end do

                call DGETRF(mm,nn,hbstb,lda,ipiv,info)
!dbg            write(*,*) 'info DGETRF is:',info

!-------------- store RHS in work2.
                do n = 0,2*nharms
                  work2(n) = res(i,j,k,7,n) * dtvolt(i,j,k,n)
!old              work2(n) = res(i,j,k,7,n)
                end do

                call DGETRS(trans,nn,nrhs,hbstb,lda,ipiv, &
                            work2,ldb,info)
!dbg            write(*,*) 'info DGETRS is:',info

!-------------- 1. updated omega is in work2. Before copying it back to
!               q, apply limitation.
!               2. store difference between new and previous values
!               of rho*omega in work1(:,:,6,:)
                do n = 0,2*nharms
                  gav        = ttrm(i,j,k,2,n)
                  rho        = q(i,j,k,1,n)
                  rhomegamin = gav*cmu*rho*sqrt(prd(i,j,k,n))
                  rhomeganew = qold(i,j,k,7,n)-work2(n)
                  res_ome = res(i,j,k,7,n) - &
                            max(0.d0,(rhomegamin-rhomeganew))/ &
                            dtvolt(i,j,k,n)
                  res(i,j,k,7,n) = res_ome
                end do

              end do
            end do
          end do

        end if

      end if

      if (lomach) then

      else if (.not.lomach) then
        deallocate( hbstb,work2,id1,ipiv )
      end if

      return
      end 

!-----------------------------------------------------------------------
      subroutine update_ttrm(nl,q,dqxi,dqeta,dqzeta,xideri,etaderi, &
        zetaderi,xiderj,etaderj,zetaderj,xiderk,etaderk,zetaderk,ttrm, &
        dist)
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,iblock,imax,jmax,kmax,iq,ittrm,idist,ivmt
      real (kind=cosa_real) q(*),ttrm(*),dqxi(*),dqeta(*),dqzeta(*),dist(*), &
        xideri(*),etaderi(*),zetaderi(*),xiderj(*),etaderj(*), &
        zetaderj(*),xiderk(*),etaderk(*),zetaderk(*)

      do iblock = 1,mynblocks
        imax   = i_imax     (iblock,nl)
        jmax   = j_jmax     (iblock,nl)
        kmax   = k_kmax     (iblock,nl)
        iq     = 1 + off_p3 (iblock,nl) * npde * dim5
        ittrm  = 1 + off_p3 (iblock,nl) * 3    * dim5
        idist  = 1 + off_p1 (iblock,nl)
        ivmt   = 1 + off_0  (iblock,nl) * 3    * dim5h
        call update_ttrm_b(q(iq),dqxi(iq),dqeta(iq),dqzeta(iq), &
          xideri(ivmt),etaderi(ivmt),zetaderi(ivmt),xiderj(ivmt), &
          etaderj(ivmt),zetaderj(ivmt),xiderk(ivmt),etaderk(ivmt), &
          zetaderk(ivmt),ttrm(ittrm),dist(idist),npde,imax,jmax,kmax, &
          nharms)

      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine update_ttrm_b(q,dqxi,dqeta,dqzeta,xideri,etaderi, &
        zetaderi,xiderj,etaderj,zetaderj,xiderk,etaderk,zetaderk,ttrm, &
        dist,npde,imax,jmax,kmax,nharms)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) npde,imax,jmax,kmax,nharms
      integer(kind=cosa_int) i,j,k,n,nh,ipde
      real (kind=cosa_real) &
        q       (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
        dqxi    (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
        dqeta   (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
        dqzeta  (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
        ttrm    (-1:imax+1,-1:jmax+1,-1:kmax+1,3   ,0:2*nharms), &
        dist    (  0:imax,0:jmax,0:kmax                  ), &
        xideri  (3,  imax,  jmax,  kmax,0:2*nharms*hbmove), &
        etaderi (3,  imax,  jmax,  kmax,0:2*nharms*hbmove), &
        zetaderi(3,  imax,  jmax,  kmax,0:2*nharms*hbmove), &
        xiderj  (3,  imax,  jmax,  kmax,0:2*nharms*hbmove), &
        etaderj (3,  imax,  jmax,  kmax,0:2*nharms*hbmove), &
        zetaderj(3,  imax,  jmax,  kmax,0:2*nharms*hbmove), &
        xiderk  (3,  imax,  jmax,  kmax,0:2*nharms*hbmove), &
        etaderk (3,  imax,  jmax,  kmax,0:2*nharms*hbmove), &
        zetaderk(3,  imax,  jmax,  kmax,0:2*nharms*hbmove)

      real (kind=cosa_real) rho,tke,ome,p,t,mu,cdt,acdt,arg1,f1,gav,bav,xix,xiy, &
        xiz,etax,etay,etaz,zetax,zetay,zetaz,tkexi,tkeeta,tkezeta,omexi, &
        omeeta,omezeta,dtkedx,dtkedy,dtkedz,domedx,domedy,domedz

      if (kom) then

        do n=0,2*nharms

          do k = 1,kmax-1
            do j = 1,jmax-1
              do i = 1,imax-1
                ttrm(i,j,k,1,n) = 0.d0
                ttrm(i,j,k,2,n) = gkw
                ttrm(i,j,k,3,n) = bkw
              end do
            end do
          end do

        end do

      else if (kom_bsl.or.kom_sst) then

        do n=0,2*nharms
          nh = n*hbmove

          do k = 1,kmax-1
            do j = 1,jmax-1
              do i = 1,imax-1

                rho = q(i,j,k,1,n)
                tke = q(i,j,k,6,n)
                ome = q(i,j,k,7,n)
                p   = q(i,j,k,5,n)
                t   = gamma*p/rho
                mu  = (stemp+1)/(t+stemp) * (t**1.5d0)

                xix   = (xideri  (1,i+1,j  ,k  ,nh) + &
                         xideri  (1,i  ,j  ,k  ,nh) ) / 2
                xiy   = (xiderj  (2,i  ,j+1,k  ,nh) + &
                         xiderj  (2,i  ,j  ,k  ,nh) ) / 2
                xiz   = (xiderk  (3,i  ,j  ,k+1,nh) + &
                         xiderk  (3,i  ,j  ,k  ,nh) ) / 2
                etax  = (etaderi (1,i+1,j  ,k  ,nh) + &
                         etaderi (1,i  ,j  ,k  ,nh) ) / 2
                etay  = (etaderj (2,i  ,j+1,k  ,nh) + &
                         etaderj (2,i  ,j  ,k  ,nh) ) / 2
                etaz  = (etaderk (3,i  ,j  ,k+1,nh) + &
                         etaderk (3,i  ,j  ,k  ,nh) ) / 2
                zetax = (zetaderi(1,i+1,j  ,k  ,nh) + &
                         zetaderi(1,i  ,j  ,k  ,nh) ) / 2
                zetay = (zetaderj(2,i  ,j+1,k  ,nh) + &
                         zetaderj(2,i  ,j  ,k  ,nh) ) / 2
                zetaz = (zetaderk(3,i  ,j  ,k+1,nh) + &
                         zetaderk(3,i  ,j  ,k  ,nh) ) / 2

                tkexi   = dqxi  (i,j,k,6,n)
                tkeeta  = dqeta (i,j,k,6,n)
                tkezeta = dqzeta(i,j,k,6,n)
                omexi   = dqxi  (i,j,k,7,n)
                omeeta  = dqeta (i,j,k,7,n)
                omezeta = dqzeta(i,j,k,7,n)

                dtkedx = tkexi*xix + tkeeta*etax + tkezeta*zetax
                dtkedy = tkexi*xiy + tkeeta*etay + tkezeta*zetay
                dtkedz = tkexi*xiz + tkeeta*etaz + tkezeta*zetaz
                domedx = omexi*xix + omeeta*etax + omezeta*zetax
                domedy = omexi*xiy + omeeta*etay + omezeta*zetay
                domedz = omexi*xiz + omeeta*etaz + omezeta*zetaz

                arg1 = dmax1(dsqrt(tke)/0.09d0/ome/(dist(i,j,k)+epst), &
                             500*mu/rho/ome/(dist(i,j,k)+epst)**2* &
                             machfs/reyno)
                cdt    = 2*rho*sigw2/ome* &
                         (dtkedx*domedx+dtkedy*domedy+dtkedz*domedz)
                acdt   = dmax1(cdt,1.d-20)
                arg1   = dmin1(arg1,4*rho*sigw2*tke/acdt/ &
                                    (dist(i,j,k)+epst)**2)
                arg1   = arg1**4
                arg1   = dmin1(arg1,1.d2)
                f1     = (dexp(2*arg1)-1)/(dexp(2*arg1)+1)
                gav    = f1*g1    + (1-f1)*g2
                bav    = f1*b1    + (1-f1)*b2

                ttrm(i,j,k,1,n) = cdt*(1-f1)
                ttrm(i,j,k,2,n) = gav
                ttrm(i,j,k,3,n) = bav

              end do
            end do
          end do

        end do

      end if

      return
      end
