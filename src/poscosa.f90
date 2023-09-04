
      program poscosa

!-----------------------------------------------------------------------
!     Developed and maintained by: M. Sergio Campobasso
!
!     with contributions from:     Jernej      Drofelnik
!                                  Andreas     Piskopakis
!                                  Mohammad H. Baba-Ahmadi
!                                  Carlos      Perez-Arroyo
!-----------------------------------------------------------------------
!     Lancaster, March 2021
!-----------------------------------------------------------------------
!                  : 
!-----------------------------------------------------------------------
!     In this version force functionals can be computed only on a 
!     particular wall type (e.g. 14 or 15); wall class is 'wallmask'.
!     NOTE: viscous walls MUST always have wallmask 15, but we need
!           wallmask input for cases in which domain has both viscous
!           and inviscid walls.
!-----------------------------------------------------------------------

      use cosa_precision

      implicit none
      include 'common.block'
      include 'cosa.inc'

      integer(kind=cosa_int) msurf
      parameter (msurf = 2000)
      integer(kind=cosa_int) imax,jmax,kmax,ijkmax,ijkmaxg,iblock,nsurf, &
                surf_dim(msurf),ijkmax_surf(msurf,3)
      integer(kind=cosa_int) itime1,itime1s,itime1e,nl,ixyz,ixyzc,ixyzc1,iblk, &
                iblk_pos

      real (kind=cosa_real) qp1, qp2, si, sj, sk, x(*), y(*), z(*), xdot(*), &
           ydot(*), zdot(*), x0, y0, z0, xc(*), yc(*), zc(*), xyzc(*), &
           xyzcdot(*), dx, dy, dz,workx, worky, workz, &
           rad, dist, xwall, ywall, zwall, xgwall, ygwall, zgwall, &
           vol, xideri,  xiderj, xiderk, etaderi, etaderj, etaderk, &
           zetaderi, zetaderj, zetaderk, ntimes
      real (kind=cosa_real),allocatable::p(:,:,:),cf(:,:,:),yplus(:,:,:), &
           xsurf(:,:),ysurf(:,:),zsurf(:,:), &
           cl_p (:,:), cd_p (:,:),cz_p (:,:), cmx_p(:,:), &
           cmy_p(:,:), cm_p (:,:),cl_v (:,:), cd_v (:,:), &
           cz_v (:,:), cmx_v(:,:),cmy_v(:,:), cm_v (:,:)


      integer startcount, endcount, countrate

      integer(kind=cosa_int) bctopo(*),cutopo,percutopo
      logical amcontrol
      character*72 flowtec1

!     qqp1    : tecplot array of primitive variables (rho,u,v,p). It 
!               uses the memory of dqp, (0:imax+1,0:jmax+1,npde), but 
!               the entries actually used are (0:imax,0:jmax,npde).
!     qqp2    : tecplot array containing temperature and Mach number. It
!               uses the memory of dqm, (0:imax+1,0:jmax+1,npde), but 
!               the entries actually used are (0:imax,0:jmax,1:2).
      pointer (pqp1 , qp1), (pqp2 , qp2), &
              (psi  ,  si), (psj  ,  sj), (psk  ,  sk), &
              (px   ,   x), (py   ,   y), (pz   ,   z), &
              (px0  ,  x0), (py0  ,  y0), (pz0  ,  z0), &
              (pxc  ,  xc), (pyc  ,  yc), (pzc  ,  zc), &
              (pxyzc , xyzc), &
              (pdx  ,  dx), (pdy  ,  dy), (pdz  ,  dz), &
              (pxdot,xdot), (pydot,ydot), (pzdot,zdot), &
              (pxyzcdot,xyzcdot), &
              (pworkx , workx),(pworky , worky),(pworkz , workz), &
              (pvol , vol),   (prad , rad),   (pdist,dist), &
              (pxideri,xideri),  (pxiderj,xiderj),  (pxiderk,xiderk), &
              (petaderi,etaderi),(petaderj,etaderj),(petaderk,etaderk), &
              (pzetaderi,zetaderi),(pzetaderj,zetaderj), &
              (pzetaderk,zetaderk), &
              (pxwall,xwall),  (pywall,ywall),  (pzwall,zwall), &
              (pxgwall,xgwall),(pygwall,ygwall),(pzgwall,zgwall), &
              (pbctopo,bctopo),(pcutopo,cutopo),(ppercutopo,percutopo)

      code = 'poscosa'

      call gettime(startcount, countrate)

      call amcontroller(amcontrol)

      if (amcontrol) write(*,*) 'calling input'
      call input()

      call paralleldecomposition()

      call movedatafordecomposition()
 
      nlevel = 1

      if (dualt) then
!old    ntime = wrifl(2) - wrifl(1) + 1
        ntimes = wrifl(2) - wrifl(1) 
        ntime = int(ntimes/wrifl(3)) + 1
      else if (harbal) then
        ntime    = 2*nharms + 1
!tmp    nharms   = 0
        wrifl(1) = -1
      else if (.not.unsteady) then
        ntime = 1
      end if

      if (amcontrol) write(*,*) 'calling alloc'
      call alloc()

      nl = 1

      call setoffs (nl)

      px         = p_x(nl)
      py         = p_y(nl)
      pz         = p_z(nl)
      pxc        = p_xc(nl)
      pyc        = p_yc(nl)
      pzc        = p_zc(nl)
      pxyzc      = p_qp(nl)
      pqp1       = p_dqp(nl)
      pqp2       = p_dqm(nl)
      pbctopo    = p_bctopo(nl)
      pcutopo    = p_cutopo(nl)
      ppercutopo = p_percutopo(nl)
      if (moving) then
        px0      = p_x0(nl)
        py0      = p_y0(nl)
        pz0      = p_z0(nl)
        pdx      = p_dx(nl)
        pdy      = p_dy(nl)
        pdz      = p_dz(nl)
        pxdot    = p_xdot(nl)
        pydot    = p_ydot(nl)
        pzdot    = p_zdot(nl)
        pxyzcdot = p_qm(nl)
      end if

      call setdim_bc(bctopo,nl)
      call setdimcut(cutopo,percutopo,nl)

      ijkmaxg   = -1
      do iblock = 1,nblocks
        if (ijkmaxg.lt.ijk_ijkmax(iblock,nl)) then
         ijkmaxg   = ijk_ijkmax(iblock,nl)
        end if
      end do
      ijkmaxg = ijkmaxg**2
      allocate(p(ijkmaxg-1,msurf,1),xsurf(ijkmaxg-1,msurf), &
               ysurf(ijkmaxg-1,msurf),zsurf(ijkmaxg-1,msurf), &
               cl_p (ijkmaxg-1,msurf), cd_p (ijkmaxg-1,msurf), &
               cz_p (ijkmaxg-1,msurf), cmx_p(ijkmaxg-1,msurf), &
               cmy_p(ijkmaxg-1,msurf), cm_p (ijkmaxg-1,msurf), &
               cl_v (ijkmaxg-1,msurf), cd_v (ijkmaxg-1,msurf), &
               cz_v (ijkmaxg-1,msurf), cmx_v(ijkmaxg-1,msurf), &
               cmy_v(ijkmaxg-1,msurf), cm_v (ijkmaxg-1,msurf))

      if (wldata.and.viscous) then
        allocate(cf(ijkmaxg-1,msurf,1),yplus(ijkmaxg-1,msurf,1))
      end if

      if (.not.unsteady) then
        simtime = 0.d0
      end if

!---- read in vertex-based steady mesh (x,y,z)
      if (amcontrol) write(*,*) 'calling readgrid'
      call readgrid(nl,x,y,z,qp1,dx,dy,dz,qp2)
!---- calculate distances only if required
      if (viscous.and.wldata) then
!       determine xyzwall coordinates, private to each block
        pxwall   = p_xwall(nl)
        pywall   = p_ywall(nl)
        pzwall   = p_zwall(nl)
        if (amcontrol) write(*,*) 'calling extract_wall',nl
        call extract_wall(nl,x,y,z,xwall,ywall,zwall,bctopo)
!       merge local xyzwall into global xyzgwall
        pxgwall   = p_xgwall(nl)
        pygwall   = p_ygwall(nl)
        pzgwall   = p_zgwall(nl)
        if (amcontrol) write(*,*) 'calling merge_wall',nl
        call merge_wall(nl,xwall,ywall,zwall,xgwall,ygwall,zgwall)
        if (debug) then
          call print_local_wall(nl,xwall,ywall,zwall)
        end if
!       calculate distances from walls
        pdist    = p_dist(nl)
        pvol     = p_vol(nl)
        imax     = i_imax(1,nl)
        jmax     = j_jmax(1,nl)
        kmax     = k_kmax(1,nl)
        if (amcontrol) write(*,*) 'calling dist2wall',nl
        call dist2wall(nl,x,y,z,xgwall,ygwall,zgwall,dist,vol,bctopo)
      end if

      nsurf = 0
      do iblk_pos = 1,nblock_pos

        iblock   = iblock_pos(iblk_pos)
        imax     = i_imax(iblock,nl)
        jmax     = j_jmax(iblock,nl)
        kmax     = k_kmax(iblock,nl)
        ixyz     = 1 + off_p2 (iblock,nl)
        ixyzc1   = 1 + off_p3 (iblock,nl) * npde
        iblk     = 1 + off_bct(iblock,nl)

!------ build cell-centred mesh (xyzc)
        call vert2cell_gb_c(x(ixyz),y(ixyz),z(ixyz),xyzc(ixyzc1),imax, &
                            jmax,kmax,npde,0)

!------ extract coordinates of surface nodes in relative frame.
        call surf_coord(imax,jmax,kmax,npde,ijkmaxg,nsurf,msurf, &
                        nbcs(iblock),xyzc(ixyzc1),bctopo(iblk),xsurf, &
                        ysurf,zsurf,surf_dim,ijkmax_surf)

      end do

      if (unsteady) then
        itime1s = 1
        itime1e = ntime
      else
        itime1s = 1
        itime1e = 1
      end if

!cc   if (moving.and.harbal) then
!cc     call copy_array(1,1,1,1, 2,nl,x0,x, &
!cc                     1,1,1,1, 0)
!cc     call copy_array(1,1,1,1, 2,nl,y0,y, &
!cc                     1,1,1,1, 0)
!cc     call movegrid(nl,x,y,xdot,ydot,dx,dy,x0,y0,xdot2,ydot2)
!cc   end if

      psi      = p_si(nl)
      psj      = p_sj(nl)
      psk      = p_sk(nl)
      pxideri  = p_xideri(nl)
      pxiderj  = p_xiderj(nl)
      pxiderk  = p_xiderk(nl)
      petaderi = p_etaderi(nl)
      petaderj = p_etaderj(nl)
      petaderk = p_etaderk(nl)
      pzetaderi= p_zetaderi(nl)
      pzetaderj= p_zetaderj(nl)
      pzetaderk= p_zetaderk(nl)
      pvol     = p_vol(nl)
      prad     = p_rad(nl)

!old  if (.not.moving) then
      if ((.not.moving).or.(.not.unsteady.and.moving)) then
        if (moving) then
          pxdot    = p_xdot(nl)
          pydot    = p_ydot(nl)
          pzdot    = p_zdot(nl)
          px0      = p_x0(nl)
          py0      = p_y0(nl)
          pz0      = p_z0(nl)
          pdx      = p_dx(nl)
          pdy      = p_dy(nl)
          pdz      = p_dz(nl)
          call copy_array(1,1,1,1, 2,nl,x0,x, &
                          1,1,1,1, 0)
          call copy_array(1,1,1,1, 2,nl,y0,y, &
                          1,1,1,1, 0)
          call copy_array(1,1,1,1, 2,nl,z0,z, &
                          1,1,1,1, 0)
          call movegrid(nl,x,y,z,xdot,ydot,zdot,dx,dy,dz,x0,y0,z0)
          if (debug) then
            write(filename,'(''grid_unstd_itime.'',i1,''.lev.'',i1,''. &
     &dat'')') 0,nl
            call wrgrid(nl,x,y,z)
          end if
        end if

        call vols(nl,x,y,z,xdot,ydot,zdot,vol,rad)
        call imetrix(nl,x,y,z,xdot,ydot,zdot,si,sj,sk)
        if (viscous) then
          call volex(nl,vol)
          call cutman_v(nl,cutopo,vol)
          if (persec) call pcutman_v(nl,percutopo,vol)
          call coordex(nl,x,y,z)
          call cutman_x(nl,cutopo,x,y,z)
          if (persec) call pcutman_x(nl,percutopo,x,y,z)
          call vmetrix(nl,x,y,z,xideri,xiderj,xiderk,etaderi,etaderj, &
                       etaderk,zetaderi,zetaderj,zetaderk,vol)
        end if
      end if

      do itime1 = itime1s,itime1e

        if (dualt) then
          itime = wrifl(1) + ( itime1 - 1 ) * wrifl(3)
          if (itime.le.9) then
            write(flowtec1,'(''000'',i1)') itime
          elseif (itime.le.99) then
            write(flowtec1,'(''00'',i2)') itime
          elseif (itime.le.999) then
            write(flowtec1,'(''0'',i3)') itime
          else
            write(flowtec1,'(i4)') itime
          end if
          flowtec1 = 'flow_tec_'//flowtec1
          flowtec = trim(flowtec1)//'.dat'
        else if (harbal) then
          itime = wrifl(1) + itime1
          if (itime.le.9) then
            write(flowtec1,'(''00'',i1)') itime
          elseif (itime.le.99) then
            write(flowtec1,'(''0'',i2)') itime
          elseif (itime.le.999) then
            write(flowtec1,'(i3)') itime
          end if
          flowtec1 = 'flow_tec_hb_'//flowtec1
          flowtec = trim(flowtec1)//'.dat'
        else if (.not.unsteady) then
          flowtec = 'flow_tec_steady.dat'
        end if
        
        if (amcontrol) write(*,*) 'calling rdtec'
        call rdtec(nl,xc,yc,zc,qp1,qp2)
!       msc, 21/09/10: the time-dependent cell-centred grid coordinates
!                      xc and yc read in by rdtec are used for:
!                        a. writing grid coordinates of 
!                           new tecplot file with addition of cp and/or 
!                           vorticity.
!                        b. calculating required velocity derivatives
!                           in direction orthogonal to wall in 
!                           routine sfht.
        call cp_wall(nl,ijkmaxg,nsurf,msurf,p(1,1,1),qp1,bctopo, &
                     surf_dim)

!jd     if (moving.and.(dualt.or.(harbal.and.(itime1.eq.itime1s)))) then
        if (moving.and.(dualt.or.harbal)) then

          if (.not.harbal) zeit(0) = simtime

          pxdot    = p_xdot(nl)
          pydot    = p_ydot(nl)
          pzdot    = p_zdot(nl)
          px0      = p_x0(nl)
          py0      = p_y0(nl)
          pz0      = p_z0(nl)
          pdx      = p_dx(nl)
          pdy      = p_dy(nl)
          pdz      = p_dz(nl)
          pcutopo  = p_cutopo(nl)

          if (itime1.eq.1) then
            call copy_array(1,1,1,1, 2,nl,x0,x, &
                            1,1,1,1, 0)
            call copy_array(1,1,1,1, 2,nl,y0,y, &
                            1,1,1,1, 0)
            call copy_array(1,1,1,1, 2,nl,z0,z, &
                            1,1,1,1, 0)
          end if

          call bodypos()
          call movegrid(nl,x,y,z,xdot,ydot,zdot,dx,dy,dz,x0,y0,z0)
          call vols(nl,x,y,z,xdot,ydot,zdot,vol,rad)
          call imetrix(nl,x,y,z,xdot,ydot,zdot,si,sj,sk)
          if (viscous) then
            call volex(nl,vol)
            call cutman_v(nl,cutopo,vol)
            if (persec) call pcutman_v(nl,percutopo,vol)
            call coordex(nl,x,y,z)
            call cutman_x(nl,cutopo,x,y,z)
            if (persec) call pcutman_x(nl,percutopo,x,y,z)
            call vmetrix(nl,x,y,z,xideri,xiderj,xiderk,etaderi,etaderj, &
              etaderk,zetaderi,zetaderj,zetaderk,vol)
          end if
          do iblock =1,nblocks
            imax     = i_imax(iblock,nl)
            jmax     = j_jmax(iblock,nl)
            kmax     = k_kmax(iblock,nl)
            ixyz     = 1 + off_p2 (iblock,nl)        * dim5h
            ixyzc1   = 1 + off_p3 (iblock,nl) * npde * dim5

!---------- build cell-centred mesh velocities (xycdot)
            call vert2cell_gb_v(xdot(ixyz),ydot(ixyz),zdot(ixyz), &
                 xyzcdot(ixyzc1),imax,jmax,kmax,npde,nharms)
          end do

        end if

        if (wldata.and.viscous) then

!old      call sfht(nl,itime1,ijkmaxg,nsurf,msurf,surf_dim, &
!old                cf(1,1,itime1),yplus(1,1,itime1),xc,yc,zc,xyzcdot, &
!old                dist,si,sj,sk,xideri,xiderj,xiderk,etaderi,etaderj, &
!old                etaderk,zetaderi,zetaderj,zetaderk,qp1,qp2,bctopo)
          call sfht(nl,itime1,ijkmaxg,nsurf,msurf,surf_dim, &
                    cf(1,1,1),yplus(1,1,1),xc,yc,zc,xyzcdot, &
                    dist,si,sj,sk,xideri,xiderj,xiderk,etaderi,etaderj, &
                    etaderk,zetaderi,zetaderj,zetaderk,qp1,qp2,bctopo)

        end if

        if (viscous.and.blprof.and.(.not.unsteady)) then
          write(*,*) 'Aborting!!! Blprof is not yet implemented!'
          stop
!         call blprofil(nl,itime1,xc,yc,xycdot,si,sj,qp1,qp2,dist, &
!                       bctopo)
        end if

!------ msc August 2010. I temporarily comment the following since
!                        I do not have time to revise this right now.
!tmp    if (calvort) then
!tmp      call vorticity(qp1,qp2,imax,jmax,npde,xc,yc)
!tmp    end if
!tmp    if (calcp) then
!tmp      call cp_glob(qp1,qp2,imax,jmax,npde)
!tmp    end if
!tmp    if (calvort.or.calcp) then
!tmp      call writec(imax,jmax,npde,xc,yc,qp1,qp2,simutime(itime1))
!tmp    end if

!old  end do

!------ We require rotated mesh at cell centres in order to calculate
!       moment coefficients properly

        pworkx   = p_xc(nl)
        pworky   = p_yc(nl)
        pworkz   = p_zc(nl)
        call vert2cell_g(nl,x,y,z,workx,worky,workz)

        call forces(nl,ijkmaxg,nsurf,msurf,qp1,qp2,bctopo,surf_dim, &
                    si,sj,sk,xideri,xiderj,xiderk,etaderi,etaderj, &
                    etaderk,zetaderi,zetaderj,zetaderk,x,y,z,xyzcdot, &
                    cl_p,cl_v,cd_p,cd_v,cz_p,cz_v,cmx_p,cmx_v,cmy_p, &
                    cmy_v,cm_p,cm_v)


        call wr_cp(ijkmaxg,p,xsurf,ysurf,zsurf,nsurf,msurf,surf_dim, &
                   ijkmax_surf)
        if (wldata.and.viscous) then
          call wr_cfyp(ijkmaxg,cf,yplus,xsurf,ysurf,zsurf,nsurf,msurf, &
                       surf_dim,ijkmax_surf)
        end if

        call wr_forces(ijkmaxg,xsurf,ysurf,zsurf,nsurf,msurf,surf_dim, &
                       ijkmax_surf,cl_p,cl_v,cd_p,cd_v,cz_p,cz_v,cmx_p, &
                       cmx_v,cmy_p,cmy_v,cm_p,cm_v)

      end do

      call system('rm -f .lock')

      deallocate(p,xsurf,ysurf,zsurf)
      deallocate(cl_p,cd_p,cz_p,cmx_p,cmy_p,cm_p,cl_v,cd_v,cz_v,cmx_v, &
                 cmy_v,cm_v)

      if (wldata.and.viscous) then
        deallocate(cf,yplus)
      end if

      if(amcontrol) then
         call gettime(endcount, countrate)
         write(*,*) 'Elapsed time:',(endcount-startcount)*1.0/countrate, &
                    'seconds'
      end if

      stop
      end

!-----------------------------------------------------------------------
      subroutine sfht(nl,itime1,ijkmaxg,nsurf,msurf,surf_dim,cf,yplus, &
                      xc,yc,zc,xyzcdot,dist,si,sj,sk,xideri,xiderj, &
                      xiderk,etaderi,etaderj,etaderk,zetaderi,zetaderj, &
                      zetaderk,q1,q2,bctopo)
!-----------------------------------------------------------------------

      use cosa_precision

      implicit none

      include 'common.block'
      include 'cosa.inc'

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

      use cosa_precision

      implicit none
      include 'common.block'

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

      use cosa_precision

      implicit none
      include 'common.block'

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

      use cosa_precision

      implicit none

      include 'common.block'
      include 'cosa.inc'

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

      use cosa_precision

      implicit none
      include 'common.block'

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

      use cosa_precision

      implicit none
      include 'common.block'

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

      use cosa_precision

      implicit none
      include 'common.block'
      include 'cosa.inc'

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

      use cosa_precision

      implicit none
      include 'common.block'

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

      use cosa_precision

      implicit none

      include 'common.block'
      include 'cosa.inc'

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

      use cosa_precision

      implicit none
      include 'common.block'

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

      use cosa_precision

      implicit none
      include 'common.block'

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

      use cosa_precision

      implicit none
      include 'common.block'

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

      use cosa_precision

      implicit none
      include 'common.block'

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

      use cosa_precision

      implicit none
      include 'common.block'

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

      use cosa_precision

      implicit none
      include 'common.block'

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

      use cosa_precision

      implicit none

      include 'common.block'

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
       
      use cosa_precision

      implicit none

      include 'common.block'

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

      use cosa_precision

      implicit none
      include 'common.block'

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

      use cosa_precision

      implicit none

      include 'common.block'
      include 'cosa.inc'

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

      use cosa_precision

      implicit none
      include 'common.block'

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

      use cosa_precision

      implicit none
      include 'common.block'

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

      use cosa_precision

      implicit none
      include 'common.block'

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
       
      use cosa_precision

      implicit none

      include 'common.block'
      include 'cosa.inc'

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
       
      use cosa_precision

      implicit none

      include 'common.block'

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

