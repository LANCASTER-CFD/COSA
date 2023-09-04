
      program cosa

!-----------------------------------------------------------------------
!     F90 Implementation of COSA.
!     3D structured time-accurate Euler/RANS code featuring:
!        Roe's flux difference splitting
!        Central Implicit Residual Smoothing
!        Point-implicit Runge-Kutta update for unsteady problems
!        Low Speed Preconditioning
!        K-Omega turbulence model
!        Dual-time-stepping-integration for time-dependent problems
!        Harmonic Balance frequency-domain solver
!        Explicit Runge Kutta time-accurate integration
!-----------------------------------------------------------------------
!     Developed and maintained by: M. Sergio   Campobasso
!
!     with contributions from:     Adrian      Jackson
!                                  Jernej      Drofelnik
!                                  Minghan     Yan
!                                  Fabio       Gigante
!                                  Francesco   Salvadore
!                                  Michela     Botti
!                                  Giuseppe    Pace
!                                  Andreas     Piskopakis
!                                  Grant       McLelland
!                                  Carlos      Perez-AArroyo
!                                  Mohammad H. Baba-Ahmadi
!                                  Giuseppe    Onorato
!                                  Stefania S. Latorre
!                                  Francois    Fraysse
!-----------------------------------------------------------------------
!     Lancaster/Glasgow/Bari/Edinburgh  June 2023
!-----------------------------------------------------------------------

      use cosa_cgns_utilities
      use common_variables
      use cosa_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,niter
      integer(kind=cosa_int) istart,irk

      real(kind=cosa_real) q, q1, q2, qp(*), qm(*), dq, dqm, pr, res, rhs, mut, &
          ttrm, si, sj, sk,x(*), y(*), z(*), x2, y2, z2, &
          x0, y0, z0, workx, worky, workz, dx(*), dy(*), dz(*), &
          dx2, dy2, dz2, xdot, ydot, zdot, &
          xwall, ywall, zwall, xgwall, ygwall, zgwall, vol, dist, &
          fidist, fjdist, fkdist, rad, xideri, xiderj, xiderk, &
          etaderi, etaderj, etaderk, zetaderi, zetaderj, zetaderk

      real(kind=cosa_single) tstart(2),tend(2),ttotal,etime
      double precision starttime, endtime
      integer startcount, endcount, countrate
      logical amcontrol

      integer(kind=cosa_int) bctopo,cutopo,percutopo

      logical from_fine,to_coarse,done

      character*72 flowtec1

      pointer &
         (pq   ,   q),(pq1  ,  q1),(pq2  ,  q2),(pdq  ,  dq), &
         (pqp  ,  qp),(pqm  ,  qm),(ppr  ,  pr),(prhs , rhs), &
         (pmut , mut),(pttrm,ttrm), &
         (psi  ,  si),(psj  ,  sj),(psk  ,  sk), &
         (px0  ,  x0),(py0  ,  y0),(pz0  ,  z0), &
         (px   ,   x),(py   ,   y),(pz   ,   z), &
         (pdx  ,  dx),(pdy  ,  dy),(pdz  ,  dz), &
         (px2  ,  x2),(py2  ,  y2),(pz2  ,  z2), &
         (pdx2 , dx2),(pdy2 , dy2),(pdz2 , dz2), &
         (pxdot,xdot),(pydot,ydot),(pzdot,zdot), &
         (pvol , vol),(pdist,dist),(prad , rad), &
         (pworkx   ,   workx),(pworky   ,   worky),(pworkz   ,   workz), &
         (pxwall   ,   xwall),(pywall   ,   ywall),(pzwall   ,   zwall), &
         (pxgwall  ,  xgwall),(pygwall  ,  ygwall),(pzgwall  ,  zgwall), &
         (pfidist  ,  fidist),(pfjdist  ,  fjdist),(pfkdist  ,  fkdist), &
         (pxideri  ,  xideri),(pxiderj  ,  xiderj),(pxiderk  ,  xiderk), &
         (petaderi , etaderi),(petaderj , etaderj),(petaderk , etaderk), &
         (pzetaderi,zetaderi),(pzetaderj,zetaderj),(pzetaderk,zetaderk), &
         (pbctopo  ,  bctopo),(pcutopo  ,  cutopo), &
         (ppercutopo,percutopo)

      integer :: ierr

      call initialisempi()

      code = 'cosa'
      call gettime(startcount, countrate)
      call getmpitime(starttime)
      ttotal = etime(tstart)

      call amcontroller(amcontrol)

      if (amcontrol.and.debug) write(*,*) 'calling input'
      
      call input()
      
      call paralleldecomposition()

      call movedatafordecomposition()

      if(amcontrol) then
        open(37, file='hist.dat',       status='replace')
        open(2 , file='diagnostic.dat', status='replace')
        open(39, file='histg.dat',      status='replace')

        if (kom.or.kom_bsl.or.kom_sst) then
          open(41, file='histur.dat',     status='replace')
        end if

        if (debug) then
          open(8,file='debug.dat',status='replace')
        end if
      end if

      if (amcontrol.and.debug) write(*,*) 'calling alloc'
      call alloc()

      do nl=1,nlevel
        pbctopo    = p_bctopo(nl)
        pcutopo    = p_cutopo(nl)
        ppercutopo = p_percutopo(nl)
        call setoffs(nl)
        call setdim_bc(bctopo,nl)
        call setdimcut(cutopo,percutopo,nl)
        if (nl.eq.1) then
          px       = p_x(nl)
          py       = p_y(nl)
          pz       = p_z(nl)
          pdx      = p_dx(nl)
          pdy      = p_dy(nl)
          pdz      = p_dz(nl)
          pqp      = p_qp(nl)
          pqm      = p_qm(nl)
          if(using_cgns) then
             if (amcontrol.and.debug) write(*,*) 'calling readcgnsgrid'
             call readcgnsgrid(nl,x,y,z,qp,dx,dy,dz,qm)
          else
             if (amcontrol.and.debug) write(*,*) 'calling readgrid'
             call readgrid(nl,x,y,z,qp,dx,dy,dz,qm)
          end if
          if (debug) then
             write(filename,'(''grid_std.'',i1,''.dat'')') nl
             call wrgrid(nl,x,y,z)
          end if
        else if (nl.gt.1) then
          px       = p_x(nl-1)
          py       = p_y(nl-1)
          pz       = p_z(nl-1)
          px2      = p_x(nl)
          py2      = p_y(nl)
          pz2      = p_z(nl)
          pdx      = p_dx(nl-1)
          pdy      = p_dy(nl-1)
          pdz      = p_dz(nl-1)
          pdx2     = p_dx(nl)
          pdy2     = p_dy(nl)
          pdz2     = p_dz(nl)
          if (amcontrol.and.debug) write(*,*) 'calling coarsen',nl
          call coarsen(nl,x,y,z,dx,dy,dz,x2,y2,z2,dx2,dy2,dz2)
          if (debug) then
             write(filename,'(''grid_std.'',i1,''.dat'')') nl
             call wrgrid(nl,x2,y2,z2)
          end if
        end if
      end do

      if (kom.or.kom_bsl.or.kom_sst) then
        do nl=1,nlevel
          pbctopo  = p_bctopo(nl)
          px       = p_x(nl)
          py       = p_y(nl)
          pz       = p_z(nl)
          pxwall   = p_xwall(nl)
          pywall   = p_ywall(nl)
          pzwall   = p_zwall(nl)
!         Executed by all MPI procs., use local block-data
          if (amcontrol.and.debug) write(*,*) 'calling extract_wall',nl
          call extract_wall(nl,x,y,z,xwall,ywall,zwall,bctopo)
          pxgwall   = p_xgwall(nl)
          pygwall   = p_ygwall(nl)
          pzgwall   = p_zgwall(nl)
!         Executed by all MPI procs., merge local xywall into global xygwall
          if (amcontrol.and.debug) write(*,*) 'calling merge_wall',nl
          call merge_wall(nl,xwall,ywall,zwall,xgwall,ygwall,zgwall)
          if (debug) then
            call print_local_wall(nl,xgwall,ygwall,zgwall)
          end if
        end do
        do nl=1,nlevel
          pbctopo  = p_bctopo(nl)
          px       = p_x(nl)
          py       = p_y(nl)
          pz       = p_z(nl)
          pdist    = p_dist(nl)
          pfidist  = p_fidist(nl)
          pfjdist  = p_fjdist(nl)
          pfkdist  = p_fkdist(nl)
          pvol     = p_vol(nl)
          pxgwall  = p_xgwall(nl)
          pygwall  = p_ygwall(nl)
          pzgwall  = p_zgwall(nl)
          if (amcontrol.and.debug) write(*,*) 'calling dist2wall',nl
          call dist2wall(nl,x,y,z,xgwall,ygwall,zgwall,dist,vol,bctopo)
          if (kom_bsl.or.kom_sst) then
            call facedist(nl,dist,fidist,fjdist,fkdist)
          else if (kom) then
            pttrm     = p_ttrm(nl)
            call update_ttrm(nl,q,qp,qm,dqm,xideri,etaderi,zetaderi, &
              xiderj,etaderj,zetaderj,xiderk,etaderk,zetaderk,ttrm,dist)
          end if
        end do
      end if

      nl = 1

      pq       = p_q(nl)
      pq2      = p_q2(nl)
      px       = p_x(nl)
      py       = p_y(nl)
      pz       = p_z(nl)
      pmut      = p_mut(nl)

      if (irest.eq.0) then
        call init2fs(q,nl)
        if (kom.or.kom_bsl.or.kom_sst) then
          call tur_vis(q,mut,nl)
        end if
        if (amcontrol) then
          write(*,*) &
            '     *** CARRYING OUT FREESTREAM INITIALIZATION ***'
        end if
      else if (irest.eq.1) then
        call zeitdata_in
        if (amcontrol.and.debug) write(*,*) 'calling readrest'
        call readrest(q,q2,mut,nl)
      end if


      pcutopo    = p_cutopo(1)    ! FS_NEW
      ppercutopo = p_percutopo(1) ! FS_NEW
      call cutman_init(cutopo)    ! FS_NEW.msc

      if(persec) then
         call pcutman_init(percutopo) ! FS_NEW.msc
      end if

      if (.not.harbal) then
        call check_cuts(nl,cutopo,x,y,z)
      end if

!-----------------------------------------------------------------------
!---- steady multigrid integration scheme
!-----------------------------------------------------------------------
      if (.not.unsteady) then

        do nl = 1,nlevel
          px        = p_x(nl)
          py        = p_y(nl)
          pz        = p_z(nl)
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
          pvol      = p_vol(nl)
          prad      = p_rad(nl)
          pcutopo   = p_cutopo(nl)
          ppercutopo= p_percutopo(nl)
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
          if (viscous.or.(nlevel.gt.1)) then
            call volex(nl,vol)
            call cutman_v(nl,cutopo,vol)
            if (persec) call pcutman_v(nl,percutopo,vol)
          end if
          if ((nl.lt.nlevel).and.(.not.bln_prol)) then
            call vol_edges(nl,vol)
          end if
          if (viscous) then
            call coordex(nl,x,y,z)
            call cutman_x(nl,cutopo,x,y,z)
            if (persec) call pcutman_x(nl,percutopo,x,y,z)
!           call vmetrix_old(nl,x,y,z,si,sj,sk,xideri,xiderj,xiderk, &
!             etaderi,etaderj,etaderk,zetaderi,zetaderj,zetaderk,vol)
            call vmetrix(nl,x,y,z,xideri,xiderj,xiderk,etaderi,etaderj, &
              etaderk,zetaderi,zetaderj,zetaderk,vol)
          end if
          pworkx   = p_xc(nl)
          pworky   = p_yc(nl)
          pworkz   = p_zc(nl)
          call vert2cell_g(nl,x,y,z,workx,worky,workz)
          if (moving) then
            call vert2cell_g(nl,xdot,ydot,zdot,workx,worky,workz)
          end if
        end do

        nstage = nstagerk
        do irk=1,nstage
          alphas(irk) = ar4(irk)
        end do
        ncycle = lmax
        itime  = 1

        call mg()

        nl = 1
        if (wri_tec) call out_tec(nl)

!-----------------------------------------------------------------------
!---- Harmonic Balance (Duke) based on multigrid solver
!-----------------------------------------------------------------------
      else if (harbal) then

        do nl = 1,nlevel
          px        = p_x(nl)
          py        = p_y(nl)
          pz        = p_z(nl)
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
          pvol      = p_vol(nl)
          prad      = p_rad(nl)
          pcutopo   = p_cutopo(nl)
          ppercutopo= p_percutopo(nl)
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
            call bodypos()
            call movegrid(nl,x,y,z,xdot,ydot,zdot,dx,dy,dz,x0,y0,z0)
            if (debug) then
              write(filename,'(''grid_unstd_itime.'',i1,''.lev.'',i1,''. &
     &dat'')') 0,nl
              call wrgrid(nl,x,y,z)
            end if
          end if
          call check_cuts(nl,cutopo,x,y,z)
          call vols(nl,x,y,z,xdot,ydot,zdot,vol,rad)
          call imetrix(nl,x,y,z,xdot,ydot,zdot,si,sj,sk)
          if (viscous.or.(nlevel.gt.1)) then
            call volex(nl,vol)
            call cutman_v(nl,cutopo,vol)
            if (persec) call pcutman_v(nl,percutopo,vol)
          end if
          if ((nl.lt.nlevel).and.(.not.bln_prol)) then
            call vol_edges(nl,vol)
          end if
          if (viscous) then
            call coordex(nl,x,y,z)
            call cutman_x(nl,cutopo,x,y,z)
            if (persec) call pcutman_x(nl,percutopo,x,y,z)
!           call vmetrix_old(nl,x,y,z,si,sj,sk,xideri,xiderj,xiderk, &
!             etaderi,etaderj,etaderk,zetaderi,zetaderj,zetaderk,vol)
            call vmetrix(nl,x,y,z,xideri,xiderj,xiderk,etaderi,etaderj, &
              etaderk,zetaderi,zetaderj,zetaderk,vol)
          end if
          pworkx   = p_xc(nl)
          pworky   = p_yc(nl)
          pworkz   = p_zc(nl)
          call vert2cell_g(nl,x,y,z,workx,worky,workz)
          if (moving) then
            call vert2cell_g(nl,xdot,ydot,zdot,workx,worky,workz)
          end if
        end do

        nstage = nstagerk
        do irk=1,nstage
          alphas(irk) = ar4(irk)
        end do
        ncycle = lmax
        itime  = 1

        call mg()

        nl = 1
        if (wri_tec) call out_tec(nl)

!-----------------------------------------------------------------------
!---- DUAL TIME STEPPING + Runge-Kutta integration 
!-----------------------------------------------------------------------
      else if (dualt) then
        zeit(0) = simtime

        do nl = 1,nlevel
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
          pxiderj   = p_xiderj(nl)
          pxiderk   = p_xiderk(nl)
          petaderi  = p_etaderi(nl)
          petaderj  = p_etaderj(nl)
          petaderk  = p_etaderk(nl)
          pzetaderi = p_zetaderi(nl)
          pzetaderj = p_zetaderj(nl)
          pzetaderk = p_zetaderk(nl)
          pvol      = p_vol(nl)
          prad      = p_rad(nl)
          pcutopo   = p_cutopo(nl)
          ppercutopo= p_percutopo(nl)
          if (moving) then
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
!del        if (relframe) zeit(0) = 0.d0
            call bodypos()
            call movegrid(nl,x,y,z,xdot,ydot,zdot,dx,dy,dz,x0,y0,z0)
            if (debug) then
              write(filename,'(''grid_unstd_itime.'',i1,''.lev.'',i1,''. &
     &dat'')') 0,nl
              call wrgrid(nl,x,y,z)
            end if
          end if
          call vols(nl,x,y,z,xdot,ydot,zdot,vol,rad)
          call imetrix(nl,x,y,z,xdot,ydot,zdot,si,sj,sk)
          if (viscous.or.(nlevel.gt.1)) then
            call volex(nl,vol)
            call cutman_v(nl,cutopo,vol)
            if (persec) call pcutman_v(nl,percutopo,vol)
          end if
          if ((nl.lt.nlevel).and.(.not.bln_prol)) then
            call vol_edges(nl,vol)
          end if
          if (viscous) then
            call coordex(nl,x,y,z)
            call cutman_x(nl,cutopo,x,y,z)
            if (persec) call pcutman_x(nl,percutopo,x,y,z)
!           call vmetrix_old(nl,x,y,z,si,sj,sk,xideri,xiderj,xiderk, &
!             etaderi,etaderj,etaderk,zetaderi,zetaderj,zetaderk,vol)
            call vmetrix(nl,x,y,z,xideri,xiderj,xiderk,etaderi,etaderj, &
              etaderk,zetaderi,zetaderj,zetaderk,vol)
          end if
          pworkx   = p_xc(nl)
          pworky   = p_yc(nl)
          pworkz   = p_zc(nl)
          call vert2cell_g(nl,x,y,z,workx,worky,workz)
          if (moving) then
            call vert2cell_g(nl,xdot,ydot,zdot,workx,worky,workz)
          end if
        end do

        nl = 1

        if (calfor) then
          pbctopo  = p_bctopo(nl)
          pq       = p_q(nl)
          psi      = p_si(nl)
          psj      = p_sj(nl)
          psk      = p_sk(nl)
          px       = p_x(nl)
          py       = p_y(nl)
          pxdot    = p_xdot(nl)
          pydot    = p_ydot(nl)
          pzdot    = p_zdot(nl)
          pdist    = p_dist(nl)
          pdq      = p_dq(nl)
          call cons2prim(nl,0,q,dq)
          call bc(q,si,sj,sk,dist,nl,bctopo)
          call prim2cons(nl,q,dq)
          call output(nl,0,.false.)
        end if

        istart = 1
        if (dual_prop_start) then
          if(amcontrol) then
            write(*,*) 'Proper (2-state) DTS startup'
          end if
        else
          if(amcontrol) then
            write(*,*) 'First order (backward Euler) DTS startup'
          end if
        end if

        pq1      = p_q1(nl)
        prhs     = p_rhs(nl)

        do itime = istart,ntime

!-----------------------------------------------------------------------
!         WARNING: if irest=1 & (.not.dual_prop_start)
!                  there are 2 possibilities:
!                  a) simtime is found in restart. In this case simtime
!                     is incremented starting from itime=1, therefore 
!                     simtime = simtime0 + itime * dt
!                     where simtime0 is value found in restart file.
!                  b) simtime is NOT found in restart. In this case
!                     simtime is incremented starting from itime=2, so
!                     simtime = (itime-1) * dt
!-----------------------------------------------------------------------

!old      if (((irest.eq.1).and.dual_prop_start).or. &
!old          ((irest.eq.1).and.(.not.dual_prop_start).and. &
!old           (itime.gt.1).and.(.not.found_simtime)).or. &
!old          ((irest.eq.1).and.(.not.dual_prop_start).and. &
!old           (itime.ge.1).and.(found_simtime)).or. &
!old          ((irest.eq.0).and.(itime.gt.1))) then

          if (found_simtime.or.(itime.gt.1)) simtime = simtime + dt
          zeit(0) = simtime

          if (moving) call bodypos()
          if (moving.and.(.not.relframe)) then
            do nl = 1,nlevel
              px        = p_x(nl)
              py        = p_y(nl)
              pz        = p_z(nl)
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
              pvol      = p_vol(nl)
              prad      = p_rad(nl)
              pxdot     = p_xdot(nl)
              pydot     = p_ydot(nl)
              pzdot     = p_zdot(nl)
              px0       = p_x0(nl)
              py0       = p_y0(nl)
              pz0       = p_z0(nl)
              pdx       = p_dx(nl)
              pdy       = p_dy(nl)
              pdz       = p_dz(nl)
              pcutopo   = p_cutopo(nl)
              ppercutopo= p_percutopo(nl)
              call movegrid(nl,x,y,z,xdot,ydot,zdot,dx,dy,dz,x0,y0,z0)
              call vols(nl,x,y,z,xdot,ydot,zdot,vol,rad)
              call imetrix(nl,x,y,z,xdot,ydot,zdot,si,sj,sk)
              if (viscous.or.(nlevel.gt.1)) then
                call volex(nl,vol)
                call cutman_v(nl,cutopo,vol)
                if (persec) call pcutman_v(nl,percutopo,vol)
              end if
              if (viscous) then
                call coordex(nl,x,y,z)
                call cutman_x(nl,cutopo,x,y,z)
                if (persec) call pcutman_x(nl,percutopo,x,y,z)
!               call vmetrix_old(nl,x,y,z,si,sj,sk,xideri,xiderj, &
!                 xiderk,etaderi,etaderj,etaderk,zetaderi,zetaderj, &
!                 zetaderk,vol)
                call vmetrix(nl,x,y,z,xideri,xiderj,xiderk,etaderi, &
                  etaderj,etaderk,zetaderi,zetaderj,zetaderk,vol)
              end if
              if (debug) then
                write(filename,'(''grid_unstd_itime.'',i1,''.lev.'', &
     &i1,''.dat'')') itime,nl
                call wrgrid(nl,x,y,z)
              end if
              pworkx   = p_xc(nl)
              pworky   = p_yc(nl)
              pworkz   = p_zc(nl)
              call vert2cell_g(nl,x   ,y   ,z   ,workx,worky,workz)
              call vert2cell_g(nl,xdot,ydot,zdot,workx,worky,workz)
            end do
          end if

          if(amcontrol) then
            write(*,*)
            write(*,*) 'itime, simtime', itime, simtime
            write(*,*)
          end if

          nl             = 1
          call copy_array(1,npde,1,1, 3,nl,q1,q, &
                          1,npde,1,1,-1)

!         DUAL TIME-STEP
          ncycle = lmax
          nstage = nstagerk
          do irk=1,nstage
            alphas(irk) = ar4(irk)
          end do
          pvol     = p_vol(nl)
          call rhs_uns(nl,rhs,q1,q2,vol)
          if ((ramping1.or.ramping2).and.(itime.gt.ntime_ramp)) then
            ramping1 = .false.
            ramping2 = .false.
          end if
          call mg()

          call copy_array(1,npde,1,1, 3,nl,q2,q1, &
                          1,npde,1,1,-1)

!-------- write output file 
          if ((wri_tec).and.(itime.ge.wrifl(1)).and.(itime.le.wrifl(2)) &
              .and.(mod((itime-wrifl(1)),wrifl(3)).eq.0))            then
            call out_tec(nl)
          end if

!------ increment physical time-step
        end do

!-----------------------------------------------------------------------
!---- UNsteady 4-stage Runge-Kutta integration scheme
!-----------------------------------------------------------------------
      else if (rgkuns) then
        px        = p_x(nl)
        py        = p_y(nl)
        pz        = p_z(nl)
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
        pvol      = p_vol(nl)
        prad      = p_rad(nl)
        pcutopo   = p_cutopo(nl)
        ppercutopo= p_percutopo(nl)
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
            write(filename,'(''grid_unstd_itime.'',i1,''.lev.'',i1,''.da &
     &t'')') 0,nl
            call wrgrid(nl,x,y,z)
          end if
        end if
        call vols(nl,x,y,z,xdot,ydot,zdot,vol,rad)
        call imetrix(nl,x,y,z,xdot,ydot,zdot,si,sj,sk)
        if (viscous) then
          call volex(nl,vol)
          call cutman_v(nl,cutopo,vol)
          if (persec) call pcutman_v(nl,percutopo,vol)
        end if
        if (viscous) then
          call coordex(nl,x,y,z)
          call cutman_x(nl,cutopo,x,y,z)
          if (persec) call pcutman_x(nl,percutopo,x,y,z)
!         call vmetrix_old(nl,x,y,z,si,sj,sk,xideri,xiderj,xiderk, &
!           etaderi,etaderj,etaderk,zetaderi,zetaderj,zetaderk,vol)
          call vmetrix(nl,x,y,z,xideri,xiderj,xiderk,etaderi,etaderj, &
            etaderk,zetaderi,zetaderj,zetaderk,vol)
        end if
        pworkx   = p_xc(nl)
        pworky   = p_yc(nl)
        pworkz   = p_zc(nl)
        call vert2cell_g(nl,x,y,z,workx,worky,workz)
        if (moving) then
          call vert2cell_g(nl,xdot,ydot,zdot,workx,worky,workz)
        end if

        if (calfor) then
          pbctopo  = p_bctopo(nl)
          pq       = p_q(nl)
          psi      = p_si(nl)
          psj      = p_sj(nl)
          psk      = p_sk(nl)
          px       = p_x(nl)
          py       = p_y(nl)
          pxdot    = p_xdot(nl)
          pydot    = p_ydot(nl)
          pzdot    = p_zdot(nl)
          pdist    = p_dist(nl)
          pdq      = p_dq(nl)
          call cons2prim(nl,0,q,dq)
          call bc(q,si,sj,sk,dist,nl,bctopo)
          call prim2cons(nl,q,dq)
          call output(nl,0,.false.)
        end if

        nstage = nstagerk
        do irk=1,nstage
          alphas(irk) = ar4(irk)
          betas(irk)  = br4(irk)
        end do

        call rkuns()

        if (wri_tec) call out_tec(nl)

      end if

      call system('rm -f .lock')

      if(amcontrol) then
        close(2)
        close(37)
        close(39)

        if (kom.or.kom_bsl.or.kom_sst) then
          close(41)
        end if

        if (debug) then
          close(8)
        end if
      end if

      call barrier()

      if(amcontrol) then
         ttotal = etime(tend)
         write(*,15) tend(1)-tstart(1)
 15      format('User-state CPU-time = ',f14.2,' seconds')
         call gettime(endcount, countrate)
         write(*,*) 'Elapsed time: ',(endcount-startcount)*1.0/countrate
         call getmpitime(endtime)
         write(*,*) 'MPI Elapsed time: ',(endtime-starttime)
      end if


      call finalisempi()
      
      stop 
      end

!-----------------------------------------------------------------------
      subroutine mg()
!-----------------------------------------------------------------------

      use common_variables
      use cosa_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,niter
      logical from_fine,to_coarse,done

      if (debug) then
        open(8,  file='debug.dat',      status='old',position='append')
      end if

      do mgit=1,ncycle

        if (nlevel.eq.1) then
          niter=1
        else
          if (mgit.eq.1) then
            niter=nstart
          else
!del        niter=npre
            niter=npost
          endif
        endif

        nl=1
        from_fine = .false.
        if (nlevel.gt.1) then
          to_coarse = .true.
        else
          to_coarse = .false.
        end if
        call smooth(nl,from_fine,to_coarse,done,niter)
        if (((nlevel.gt.1).and. &
             (((mgit.eq.ncycle).and.done).or.(mgit.ne.ncycle))).or. &
            (nlevel.eq.1)) then
          call output(nl,mgit,done)
        end if
        if (done) return

        if (nlevel.gt.1) then

!-------- restrict solution and residual of fine onto coarse (full weighting)
          if (debug) then
            write(*,*) 'calling restrict on level', nl
            write(8,*) 'calling restrict on level', nl
          end if
          call restrict(nl)
          nl = nl+1
          if (nlevel.eq.2) then
            niter=ncrs
          else
            niter=npre
          endif
!-------- smooth on coarser grid (unknown: correction for finer grid)
          from_fine = .true.
          if (nlevel.gt.2) then
            to_coarse = .true.
          else
            to_coarse = .false.
          end if
          call smooth(nl,from_fine,to_coarse,done,niter)

          if (nlevel.eq.3) then

!---------- restrict solution and residual of fine onto coarse (full weighting
            if (debug) then
              write(*,*) 'calling restrict on level', nl
              write(8,*) 'calling restrict on level', nl
            end if
            call restrict(nl)
            nl = nl+1
            niter=ncrs
!---------- smooth on coarser grid (unknown: correction for finer grid)
            from_fine = .true.
            to_coarse = .false.
            call smooth(nl,from_fine,to_coarse,done,niter)

!---------- prolong correction of coarse onto fine (linear interpol)
            nl   = nl-1
            if (debug) then
              write(*,*) 'calling prolong on level', nl
              write(8,*) 'calling prolong on level', nl
            end if
            call prolong(nl)
!-----------smooth on medium grid (unknown: correction for finer grid)
            niter=npost
            from_fine = .false.
            to_coarse = .false.
            call smooth(nl,from_fine,to_coarse,done,niter)

          end if

!-------- prolong correction of coarse onto fine (linear interpol)
          nl   = nl-1
          if (debug) then
            write(*,*) 'calling prolong on level', nl
            write(8,*) 'calling prolong on level', nl
          end if
          call prolong(nl)

        end if

      end do

      if (nlevel.gt.1) then
        to_coarse = .false.
        from_fine = .false.
        niter =npost
        call smooth(nl,from_fine,to_coarse,done,niter)
        if (itime.gt.0) then
          call output(nl,ncycle,done)
        end if
      end if

      if (debug) then
        close(8)
      end if

      return
      end

!-----------------------------------------------------------------------
      subroutine prtest(nl,q,pr,ic)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,iblock,imax,jmax,kmax,iq,ipr,ic
      real (kind=cosa_real) q(*),pr(*)
 
      do iblock = 1,mynblocks
        imax   = i_imax     (iblock,nl)
        jmax   = j_jmax     (iblock,nl)
        kmax   = k_kmax     (iblock,nl)
        iq     = 1 + off_p3 (iblock,nl) * npde * dim5
        ipr    = 1 + off_p3 (iblock,nl) *        dim5
         call bprtest(q(iq),pr(ipr),ic,imax,jmax,kmax,npde,nharms, &
                     iblock,lowernblock,nl)
!old                 iblock,nl)
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine bprtest(q,pr,ic,imax,jmax,kmax,npde,nharms,iblock, &
                         lowernblock,nl)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,npde,nharms,iblock,lowernblock,nl
      integer(kind=cosa_int) i,j,k,ic,n,liblock
      real (kind=cosa_real) &
           q (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           pr(-1:imax+1,-1:jmax+1,-1:kmax+1,     0:2*nharms)
      real (kind=cosa_real) rho
 
      liblock = lowernblock+iblock-1

      do n = 0,2*nharms

        do k = 1,kmax-1
          do j = 1,jmax-1
            do i = -1,imax+1
              rho     = q(i,j,k,1,n)
              pr(i,j,k,n) = q(i,j,k,5,n)
              if (pr(i,j,k,n).le.qmin(5)) then
                write(2,100) nl,liblock,i,j,k,ic
                write(*,100) nl,liblock,i,j,k,ic
#if MPI
                call abortmpi()
#else
                stop
#endif
              else if (rho.le.qmin(1)) then
                write(2,101) nl,liblock,i,j,k,ic
                write(*,101) nl,liblock,i,j,k,ic
#if MPI
                call abortmpi()
#else
                stop
#endif
              end if
            end do
          end do
        end do

        do k = 1,kmax-1
          do j = -1,0
            do i = 1,imax-1
              rho     = q(i,j,k,1,n)
              pr(i,j,k,n) = q(i,j,k,5,n)
              if (pr(i,j,k,n).le.qmin(5)) then
                write(2,100) nl,liblock,i,j,k,ic
                write(*,100) nl,liblock,i,j,k,ic
#if MPI
                call abortmpi()
#else
                stop
#endif
              else if (rho.le.qmin(1)) then
                write(2,101) nl,liblock,i,j,k,ic
                write(*,101) nl,liblock,i,j,k,ic
#if MPI
                call abortmpi()
#else
                stop
#endif
              end if
            end do
          end do
        end do

        do k = 1,kmax-1
          do j = jmax,jmax+1
            do i = 1,imax-1
              rho     = q(i,j,k,1,n)
              pr(i,j,k,n) = q(i,j,k,5,n)
              if (pr(i,j,k,n).le.qmin(5)) then
                write(2,100) nl,liblock,i,j,k,ic
                write(*,100) nl,liblock,i,j,k,ic
#if MPI
                call abortmpi()
#else
                stop
#endif
              else if (rho.le.qmin(1)) then
                write(2,101) nl,liblock,i,j,k,ic
                write(*,101) nl,liblock,i,j,k,ic
#if MPI
                call abortmpi()
#else
                stop
#endif
              end if
            end do
          end do
        end do

        do k = -1,0
          do j = 1,jmax-1
            do i = 1,imax-1
              rho     = q(i,j,k,1,n)
              pr(i,j,k,n) = q(i,j,k,5,n)
              if (pr(i,j,k,n).le.qmin(5)) then
                write(2,100) nl,liblock,i,j,k,ic
                write(*,100) nl,liblock,i,j,k,ic
#if MPI
                call abortmpi()
#else
                stop
#endif
              else if (rho.le.qmin(1)) then
                write(2,101) nl,liblock,i,j,k,ic
                write(*,101) nl,liblock,i,j,k,ic
#if MPI
                call abortmpi()
#else
                stop
#endif
              end if
            end do
          end do
        end do

        do k = kmax,kmax+1
          do j = 1,jmax-1
            do i = 1,imax-1
              rho     = q(i,j,k,1,n)
              pr(i,j,k,n) = q(i,j,k,5,n)
              if (pr(i,j,k,n).le.qmin(5)) then
                write(2,100) nl,liblock,i,j,k,ic
                write(*,100) nl,liblock,i,j,k,ic
#if MPI
                call abortmpi()
#else
                stop
#endif
              else if (rho.le.qmin(1)) then
                write(2,101) nl,liblock,i,j,k,ic
                write(*,101) nl,liblock,i,j,k,ic
#if MPI
                call abortmpi()
#else
                stop
#endif
              end if
            end do
          end do
        end do

      end do

 100  format(///,2x,'negative pressure on grid level',i2,', block ',i5, &
                    ' at i=',i3,', j=',i3,', k=',i3,//,2x, &
                    'after',i6,' iterations')
 101  format(///,2x,'negative density on grid level',i2,', block ',i5, &
                    ' at i=',i3,', j=',i3,', k=',i3,//,2x, &
                    'after',i6,' iterations')

      return
      end

!---------------------------------------------------------------------
      subroutine resid(idir,nl,flux,res)
!---------------------------------------------------------------------       

      use common_variables
      use cosa_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) idir,nl,iblock,imax,jmax,kmax,iflux,ires
      real (kind=cosa_real) res(*),flux(*)

      do iblock = 1,mynblocks
        imax   = i_imax     (iblock,nl)
        jmax   = j_jmax     (iblock,nl)
        kmax   = k_kmax     (iblock,nl)
        ires   = 1 + off_p3 (iblock,nl) * npde * dim5
        iflux  = 1 + off_0  (iblock,nl) * npde * dim5
        call bresid(flux(iflux),res(ires),idir,imax,jmax,kmax,npde, &
                    nharms)
      end do

      return
      end

!---------------------------------------------------------------------
      subroutine bresid(flux,res,idir,imax,jmax,kmax,npde,nharms)
!---------------------------------------------------------------------       

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,npde,nharms
      integer(kind=cosa_int) i,j,k,idir,ipde,n
      real (kind=cosa_real) &
           res (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           flux(   imax  ,   jmax  ,   kmax  ,npde,0:2*nharms)

      if (idir.eq.1) then

        do n = 0,2*nharms
          do ipde = 1,npde
            do k = 1,kmax-1
              do j = 1,jmax-1
                do i = 1,imax-1
                  res(i,j,k,ipde,n) = res(i,j,k,ipde,n) + &
                    flux(i+1,j,k,ipde,n) - flux(i,j,k,ipde,n)
                end do
              end do
            end do
          end do
        end do

      else if (idir.eq.2) then

        do n = 0,2*nharms
          do ipde = 1,npde
            do k = 1,kmax-1
              do j = 1,jmax-1
                do i = 1,imax-1
                  res(i,j,k,ipde,n) = res(i,j,k,ipde,n) + &
                    flux(i,j+1,k,ipde,n) - flux(i,j,k,ipde,n)
                end do
              end do
            end do
          end do
        end do

      else if (idir.eq.3) then

        do n = 0,2*nharms
          do ipde = 1,npde
            do k = 1,kmax-1
              do j = 1,jmax-1
                do i = 1,imax-1
                  res(i,j,k,ipde,n) = res(i,j,k,ipde,n) + &
                    flux(i,j,k+1,ipde,n) - flux(i,j,k,ipde,n)
                end do
              end do
            end do
          end do
        end do

      end if

      return
      end

!-----------------------------------------------------------------------
      subroutine lambda(nl,ic,q,si,sj,sk,xdot,ydot,zdot,vol,mut,work, &
                        lambdas,cutoff_pgr,cutoff_vis)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,ic,iblock,imax,jmax,kmax,iq,iimt,ixyz,ivol,icof,imut
      real (kind=cosa_real) q(*),si(*),sj(*),sk(*),xdot(*),ydot(*),zdot(*), &
           vol(*),mut(*),work(*),lambdas(*),cutoff_pgr(*),cutoff_vis(*)

      do iblock = 1,mynblocks
        imax   = i_imax     (iblock,nl)
        jmax   = j_jmax     (iblock,nl)
        kmax   = k_kmax     (iblock,nl)
        iq     = 1 + off_p3 (iblock,nl) * npde * dim5
        iimt   = 1 + off_p2 (iblock,nl) * lmet * dim5h
        ivol   = 1 + off_p1 (iblock,nl)
        icof   = 1 + off_m1 (iblock,nl) *        dim5
        ixyz   = 1 + off_p2 (iblock,nl) *        dim5h
        imut   = 1 + off_p3 (iblock,nl) *        dim5
        call lambdab(nl,ic,q(iq),si(iimt),sj(iimt),sk(iimt),xdot(ixyz), &
          ydot(ixyz),zdot(ixyz),vol(ivol),mut(imut),work(iq), &
          lambdas(iq),cutoff_pgr(icof),cutoff_vis(icof), &
          imax,jmax,kmax,npde,nharms,lmet)
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine lambdab(nl,ic,q,si,sj,sk,xdot,ydot,zdot,vol,mut,work, &
                   lambdas,cutoff_pgr,cutoff_vis, &
                   imax,jmax,kmax,npde,nharms,lmet)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,npde,nharms,lmet
      integer(kind=cosa_int) i,j,k,ic,ipde,imet,n,nh,nl
      real (kind=cosa_real) &
       q         (     -1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
       work      (     -1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
       lambdas   (     -1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
       mut       (     -1:imax+1,-1:jmax+1,-1:kmax+1,     0:2*nharms), &
       si        (lmet,0:imax+1 ,0:jmax+1 ,0:kmax+1 ,0:2*nharms*hbmove), &
       sj        (lmet,0:imax+1 ,0:jmax+1 ,0:kmax+1 ,0:2*nharms*hbmove), &
       sk        (lmet,0:imax+1 ,0:jmax+1 ,0:kmax+1 ,0:2*nharms*hbmove), &
       xdot      (     0:imax+1 ,0:jmax+1 ,0:kmax+1 ,0:2*nharms*hbmove), &
       ydot      (     0:imax+1 ,0:jmax+1 ,0:kmax+1 ,0:2*nharms*hbmove), &
       zdot      (     0:imax+1 ,0:jmax+1 ,0:kmax+1 ,0:2*nharms*hbmove), &
       cutoff_pgr(     imax-1   ,jmax-1   ,kmax-1   ,     0:2*nharms), &
       cutoff_vis(     imax-1   ,jmax-1   ,kmax-1   ,     0:2*nharms), &
       vol       (     0:imax   ,0:jmax   ,0:kmax)
      real (kind=cosa_real) sim(4),sip(4),sjm(4),sjp(4),skm(4),skp(4), &
           p,pm,pp,rho,rhom,rhop,t,tm,tp,mu,mum,mup,mutc,mutm,mutp, &
           mui,muj,muk,muti,mutj,mutk,lamvi,lamvj,lamvk,lami,lamj,lamk, &
           qc(5),volc,cutoffc

!---- compute turbulent viscosity
      if (kom.or.kom_bsl.or.kom_sst) then

        do n = 0,2*nharms
          do k = 0,kmax
            do j = 0,jmax
              do i = 0,imax
                work(i,j,k,4,n) = mut(i,j,k,n)
              end do
            end do
          end do
        end do

      else

        do n = 0,2*nharms
          do k = 0,kmax
            do j = 0,jmax
              do i = 0,imax
                work(i,j,k,4,n) = 0.d0
              end do
            end do
          end do
        end do

      end if

!-----------------------------------------------------------------------
!     cases WITHOUT low-speed preconditioning
!-----------------------------------------------------------------------
      if (.not.lomach) then

!------ calculate lambdas

        if (.not.viscous) then

          do n = 0,2*nharms
            nh = n*hbmove

            do k = 1,kmax-1
              do j = 1,jmax-1
                do i = 1,imax-1

                  do ipde = 1,5
                    qc(ipde) = q(i,j,k,ipde,n)
                  end do

                  do imet = 1,4
                    sim(imet) = si(imet,i  ,j  ,k  ,nh)
                    sip(imet) = si(imet,i+1,j  ,k  ,nh)
                    sjm(imet) = sj(imet,i,  j  ,k  ,nh)
                    sjp(imet) = sj(imet,i,  j+1,k  ,nh)
                    skm(imet) = sk(imet,i,  j  ,k  ,nh)
                    skp(imet) = sk(imet,i,  j  ,k+1,nh)
                  end do

                  call eig_eul(qc,sim,sip,sjm,sjp,skm,skp, &
                    xdot(i,j,k,nh),ydot(i,j,k,nh),zdot(i,j,k,nh), &
                    lami,lamj,lamk)

                  lambdas(i,j,k,1,n) = lami
                  lambdas(i,j,k,2,n) = lamj
                  lambdas(i,j,k,3,n) = lamk

                end do
              end do
            end do

          end do

        else if (viscous) then

          do n = 0,2*nharms
            nh = n*hbmove

            do k = 1,kmax-1
              do j = 1,jmax-1
                do i = 1,imax-1

                  rho  = q(i  ,j  ,k  ,1,n)
                  p    = q(i  ,j  ,k  ,5,n)
                  t    = gamma*p/rho
                  mu   = (stemp+1.d0)/(t+stemp) * (t**1.5d0)
                  mutc = work(i  ,j  ,k  ,4,n)

                  rhom = q(i-1,j  ,k  ,1,n)
                  pm   = q(i-1,j  ,k  ,5,n)
                  tm   = gamma*pm/rhom
                  mum  = (stemp+1.d0)/(tm+stemp) * (tm**1.5d0)
                  mutm = work(i-1,j  ,k  ,4,n)

                  rhop = q(i+1,j  ,k  ,1,n)
                  pp   = q(i+1,j  ,k  ,5,n)
                  tp   = gamma*pp/rhom
                  mup  = (stemp+1.d0)/(tp+stemp) * (tp**1.5d0)
                  mutp = work(i+1,j  ,k  ,4,n)

                  mui  = (mum  + 2*mu   + mup ) / 4
                  muti = (mutm + 2*mutc + mutp) / 4

                  rhom = q(i  ,j-1,k  ,1,n)
                  pm   = q(i  ,j-1,k  ,5,n)
                  tm   = gamma*pm/rhom
                  mum  = (stemp+1.d0)/(tm+stemp) * (tm**1.5d0)
                  mutm = work(i  ,j-1,k  ,4,n)

                  rhop = q(i  ,j+1,k  ,1,n)
                  pp   = q(i  ,j+1,k  ,5,n)
                  tp   = gamma*pp/rhom
                  mup  = (stemp+1.d0)/(tp+stemp) * (tp**1.5d0)
                  mutp = work(i  ,j+1,k  ,4,n)

                  muj  = (mum  + 2*mu   + mup ) / 4
                  mutj = (mutm + 2*mutc + mutp) / 4

                  rhom = q(i  ,j  ,k-1,1,n)
                  pm   = q(i  ,j  ,k-1,5,n)
                  tm   = gamma*pm/rhom
                  mum  = (stemp+1.d0)/(tm+stemp) * (tm**1.5d0)
                  mutm = work(i  ,j  ,k-1,4,n)

                  rhop = q(i  ,j  ,k+1,1,n)
                  pp   = q(i  ,j  ,k+1,5,n)
                  tp   = gamma*pp/rhom
                  mup  = (stemp+1.d0)/(tp+stemp) * (tp**1.5d0)
                  mutp = work(i  ,j  ,k+1,4,n)

                  muk  = (mum  + 2*mu   + mup ) / 4
                  mutk = (mutm + 2*mutc + mutp) / 4

                  do ipde = 1,5
                    qc(ipde) = q(i,j,k,ipde,n)
                  end do

                  do imet = 1,4
                    sim(imet) = si(imet,i  ,j  ,k  ,nh)
                    sip(imet) = si(imet,i+1,j  ,k  ,nh)
                    sjm(imet) = sj(imet,i,  j  ,k  ,nh)
                    sjp(imet) = sj(imet,i,  j+1,k  ,nh)
                    skm(imet) = sk(imet,i,  j  ,k  ,nh)
                    skp(imet) = sk(imet,i,  j  ,k+1,nh)
                  end do

                  volc = vol(i,j,k)

                  call eig_ns(mui,muj,muk,muti,mutj,mutk,qc,sim,sip, &
                    sjm,sjp,skm,skp,xdot(i,j,k,nh),ydot(i,j,k,nh), &
                    zdot(i,j,k,nh),volc,lami,lamj,lamk,lamvi,lamvj, &
                    lamvk)

                  lambdas(i,j,k,1,n) = lami
                  lambdas(i,j,k,2,n) = lamj
                  lambdas(i,j,k,3,n) = lamk
                  work(i,j,k,1,n) = lamvi
                  work(i,j,k,2,n) = lamvj
                  work(i,j,k,3,n) = lamvk

                end do
              end do
            end do

          end do

        end if

!-----------------------------------------------------------------------
!     cases WITH low-speed preconditioning
!-----------------------------------------------------------------------
      else if (lomach) then

!------ calculate lambdas

        if (.not.viscous) then

          do n = 0,2*nharms
            nh = n*hbmove

            do k = 1,kmax-1
              do j = 1,jmax-1
                do i = 1,imax-1

                  cutoffc = dmax1(cutoff_pgr(i,j,k,n), &
                                  cutoff_vis(i,j,k,n))

                  do ipde = 1,5
                    qc(ipde) = q(i,j,k,ipde,n)
                  end do

                  do imet = 1,4
                    sim(imet) = si(imet,i  ,j  ,k  ,nh)
                    sip(imet) = si(imet,i+1,j  ,k  ,nh)
                    sjm(imet) = sj(imet,i,  j  ,k  ,nh)
                    sjp(imet) = sj(imet,i,  j+1,k  ,nh)
                    skm(imet) = sk(imet,i,  j  ,k  ,nh)
                    skp(imet) = sk(imet,i,  j  ,k+1,nh)
                  end do

                  call eig_eul_p(qc,sim,sip,sjm,sjp,skm,skp, &
                    xdot(i,j,k,nh),ydot(i,j,k,nh),zdot(i,j,k,nh), &
                    cutoffc,lami,lamj,lamk,nl)

                  lambdas(i,j,k,1,n) = lami
                  lambdas(i,j,k,2,n) = lamj
                  lambdas(i,j,k,3,n) = lamk

                end do
              end do
            end do

          end do

        else if (viscous) then

          do n = 0,2*nharms
            nh = n*hbmove

            do k = 1,kmax-1
              do j = 1,jmax-1
                do i = 1,imax-1

                  cutoffc = dmax1(cutoff_pgr(i,j,k,n), &
                                  cutoff_vis(i,j,k,n))

                  rho  = q(i  ,j  ,k  ,1,n)
                  p    = q(i  ,j  ,k  ,5,n)
                  t    = gamma*p/rho
                  mu   = (stemp+1.d0)/(t+stemp) * (t**1.5d0)
                  mutc = work(i  ,j  ,k  ,4,n)

                  rhom = q(i-1,j  ,k  ,1,n)
                  pm   = q(i-1,j  ,k  ,5,n)
                  tm   = gamma*pm/rhom
                  mum  = (stemp+1.d0)/(tm+stemp) * (tm**1.5d0)
                  mutm = work(i-1,j  ,k  ,4,n)

                  rhop = q(i+1,j  ,k  ,1,n)
                  pp   = q(i+1,j  ,k  ,5,n)
                  tp   = gamma*pp/rhom
                  mup  = (stemp+1.d0)/(tp+stemp) * (tp**1.5d0)
                  mutp = work(i+1,j  ,k  ,4,n)

                  mui  = (mum  + 2*mu   + mup ) / 4
                  muti = (mutm + 2*mutc + mutp) / 4

                  rhom = q(i  ,j-1,k  ,1,n)
                  pm   = q(i  ,j-1,k  ,5,n)
                  tm   = gamma*pm/rhom
                  mum  = (stemp+1.d0)/(tm+stemp) * (tm**1.5d0)
                  mutm = work(i  ,j-1,k  ,4,n)

                  rhop = q(i  ,j+1,k  ,1,n)
                  pp   = q(i  ,j+1,k  ,5,n)
                  tp   = gamma*pp/rhom
                  mup  = (stemp+1.d0)/(tp+stemp) * (tp**1.5d0)
                  mutp = work(i  ,j+1,k  ,4,n)

                  muj  = (mum  + 2*mu   + mup ) / 4
                  mutj = (mutm + 2*mutc + mutp) / 4

                  rhom = q(i  ,j  ,k-1,1,n)
                  pm   = q(i  ,j  ,k-1,5,n)
                  tm   = gamma*pm/rhom
                  mum  = (stemp+1.d0)/(tm+stemp) * (tm**1.5d0)
                  mutm = work(i  ,j  ,k-1,4,n)

                  rhop = q(i  ,j  ,k+1,1,n)
                  pp   = q(i  ,j  ,k+1,5,n)
                  tp   = gamma*pp/rhom
                  mup  = (stemp+1.d0)/(tp+stemp) * (tp**1.5d0)
                  mutp = work(i  ,j  ,k+1,4,n)

                  muk  = (mum  + 2*mu   + mup ) / 4
                  mutk = (mutm + 2*mutc + mutp) / 4

                  do ipde = 1,5
                    qc(ipde) = q(i,j,k,ipde,n)
                  end do

                  do imet = 1,4
                    sim(imet) = si(imet,i  ,j  ,k  ,nh)
                    sip(imet) = si(imet,i+1,j  ,k  ,nh)
                    sjm(imet) = sj(imet,i,  j  ,k  ,nh)
                    sjp(imet) = sj(imet,i,  j+1,k  ,nh)
                    skm(imet) = sk(imet,i,  j  ,k  ,nh)
                    skp(imet) = sk(imet,i,  j  ,k+1,nh)
                  end do

                  volc = vol(i,j,k)

                  call eig_ns_p(mui,muj,muk,muti,mutj,mutk,qc,sim,sip, &
                    sjm,sjp,skm,skp,xdot(i,j,k,nh),ydot(i,j,k,nh), &
                    zdot(i,j,k,nh),volc,cutoffc,lami,lamj,lamk,lamvi, &
                    lamvj,lamvk,nl)

                  lambdas(i,j,k,1,n) = lami
                  lambdas(i,j,k,2,n) = lamj
                  lambdas(i,j,k,3,n) = lamk
                  work(i,j,k,1,n) = lamvi
                  work(i,j,k,2,n) = lamvj
                  work(i,j,k,3,n) = lamvk

                end do
              end do
            end do

          end do
    
        end if

      end if

      return
      end

!-----------------------------------------------------------------------
      subroutine smooth(nl,from_fine,to_coarse,done,niter)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,niter
      integer(kind=cosa_int) niterp,iter,irk,irs_its,idir

      real(kind=cosa_real) q, qold, qp, qm, dq, dqp, dqm, pr, mut, ttrm, &
          prd, flux, res, rhs, prec, preci, ipkpi, dtvol, dtvolt, &
          si, sj, sk, x, y, z, xdot, ydot, zdot, vol, dist, fidist, &
          fjdist, fkdist, divvel, delpls, rad, betaxi, betaeta, betazeta, &
          tbetaxi, tbetaeta, tbetazeta, xideri, xiderj, xiderk, etaderi, &
          etaderj, etaderk, zetaderi, zetaderj, zetaderk, work1, work2, &
          cutoff_pgr, cutoff_vis

      integer(kind=cosa_int) bctopo,cutopo,percutopo

      real (kind=cosa_real) errmax

      logical from_fine,to_coarse,done

      pointer &
         (pq   ,   q), (pqold,qold), (pqp  ,  qp), (pqm  ,  qm), &
         (pdq  ,  dq), (pdqp , dqp), (pdqm , dqm), (pres , res), &
         (ppr  ,  pr), (prhs , rhs), (pflux,flux), &
         (psi  ,  si), (psj  ,  sj), (psk  ,  sk), &
         (px   ,   x), (py   ,   y), (pz   ,   z), &
         (pxdot,xdot), (pydot,ydot), (pzdot,zdot), &
         (pvol , vol), (prad , rad), &
         (pdtvol   ,   dtvol),(pdtvolt  ,  dtvolt), &
         (pmut     ,     mut),(pprd     ,     prd),(pttrm    ,    ttrm), &
         (pdelpls  ,  delpls),(pdivvel  ,  divvel),(pdist    ,    dist), &
         (pfidist  ,  fidist),(pfjdist  ,  fjdist),(pfkdist  ,  fkdist), &
         (pbctopo  ,  bctopo),(pcutopo  ,  cutopo), &
         (ppercutopo , percutopo), &
         (pxideri  ,  xideri),(pxiderj  ,  xiderj),(pxiderk  ,  xiderk), &
         (petaderi , etaderi),(petaderj , etaderj),(petaderk , etaderk), &
         (pzetaderi,zetaderi),(pzetaderj,zetaderj),(pzetaderk,zetaderk), &
         (pwork1   ,   work1),(pwork2   ,   work2), &
         (pprec    ,    prec),(ppreci   ,   preci),(pipkpi   ,   ipkpi), &
         (pbetaxi  ,  betaxi),(pbetaeta , betaeta),(pbetazeta,betazeta), &
         (ptbetaxi   ,   tbetaxi), (ptbetaeta  ,  tbetaeta), &
         (ptbetazeta , tbetazeta), &
         (pcutoff_pgr,cutoff_pgr), (pcutoff_vis,cutoff_vis)

      px         = p_x(nl)
      py         = p_y(nl)
      pz         = p_z(nl)
      pq         = p_q(nl)
      pqp        = p_qp(nl)
      pqm        = p_qm(nl)
      pflux      = p_flux(nl)
      pres       = p_res(nl)
      pdtvol     = p_dtvol(nl)
      pdtvolt    = p_dtvolt(nl)
      ppr        = p_pr(nl)
      pmut       = p_mut(nl)
      pttrm      = p_ttrm(nl)
      pprd       = p_prd(nl)
      psi        = p_si(nl)
      psj        = p_sj(nl)
      psk        = p_sk(nl)
      pxideri    = p_xideri(nl)
      pxiderj    = p_xiderj(nl)
      pxiderk    = p_xiderk(nl)
      petaderi   = p_etaderi(nl)
      petaderj   = p_etaderj(nl)
      petaderk   = p_etaderk(nl)
      pzetaderi  = p_zetaderi(nl)
      pzetaderj  = p_zetaderj(nl)
      pzetaderk  = p_zetaderk(nl)
      pxdot      = p_xdot(nl)
      pydot      = p_ydot(nl)
      pzdot      = p_zdot(nl)
      pvol       = p_vol(nl)
      pdist      = p_dist(nl)
      pfidist    = p_fidist(nl)
      pfjdist    = p_fjdist(nl)
      pfkdist    = p_fkdist(nl)
      pdelpls    = p_delpls(nl)
      pdivvel    = p_divvel(nl)
      prad       = p_rad(nl)
      pqold      = p_qold(nl)
      pdq        = p_dq(nl)
      pdqp       = p_dqp(nl)
      pdqm       = p_dqm(nl)
      prhs       = p_rhs(nl)
      pprec      = p_prec(nl)
      pbetaxi    = p_betaxi(nl)
      pbetaeta   = p_betaeta(nl)
      pbetazeta  = p_betazeta(nl)
      ptbetaxi   = p_tbetaxi(nl)
      ptbetaeta  = p_tbetaeta(nl)
      ptbetazeta = p_tbetazeta(nl)
      ppreci     = p_preci(nl)
      pwork1     = p_work1(nl)
      pwork2     = p_work2(nl)
      pipkpi     = p_ipkpi(nl)
      pbctopo    = p_bctopo(nl)
      pcutopo    = p_cutopo(nl)
      ppercutopo = p_percutopo(nl)
      pcutoff_pgr = p_cutoff_pgr(nl)
      pcutoff_vis = p_cutoff_vis(nl)

      done = .false.

      if(to_coarse) then
        niterp = niter+1
      else
        niterp = niter
      endif

      do iter = 1,niterp

!------ move initial solution to qold
        if (iter.le.niter) then
          call copy_array(1,npde,0,2*nharms, 3,nl,qold,q, &
                          1,npde,0,2*nharms,-1)
        end if

        do irk = 1,nstage

          call cons2prim(nl,0,q,dq)

          call bc(q,si,sj,sk,dist,nl,bctopo)
          call cutman_q(nl,cutopo,q)
          if (persec) call pcutman_q(nl,percutopo,q)

          if ((kom.or.kom_bsl.or.kom_sst).and. &
              (mgit.eq.1).and.(iter.eq.1).and.(irk.eq.1).and. &
!     &        ((nl.eq.1).or.(from_fine.and.(.not.prd_all_lvls)))) then &
              (from_fine.and.(.not.prd_all_lvls))) then
!---------- msc, 21/07/2021. If prd_all_lvls is false, it calculates
!             muT also on coarse grid levels (does not calculate turbul
!             prod term on coarser grids, uses instead fine grid
!             restriction.) &
!old          ((nl.eq.1).or.from_fine)) then
            if (kom.or.kom_bsl) then
                write(*,*) 'Recalculating mut (tur_vis)'
              call tur_vis(q,mut,nl)
            else if (kom_sst) then
!------------ comp./store (dqxi,dqeta,dqzeta) in (qp,qm,dq) respectivetly
              call copy_array(1,npde,0,2*nharms, 3,nl,res,q, &
                              1,npde,0,2*nharms, 3)
              calmut = .false.
              call q_face(1,nl,bctopo,res,dq,mut,flux)
              call q_der_c(1,nl,dq,qp)
              call q_face(2,nl,bctopo,res,dq,mut,flux)
              call q_der_c(2,nl,dq,qm)
              call q_face(3,nl,bctopo,res,dq,mut,flux)
              call q_der_c(3,nl,dq,dqm)
!------------ revert to conservative form needed for tur_vis_sst
              call prim2cons(nl,res,dq)
               write(*,*) 'Recalculating mut (tur_vis_sst)'
              call tur_vis_sst(res,mut,qp,qm,dqm,xideri,etaderi,zetaderi, &
                xiderj,etaderj,zetaderj,xiderk,etaderk,zetaderk,dist,nl)
              call zero(1,npde,0,2*nharms, 3,nl,res, &
                        1,npde,0,2*nharms, 3)
            end if
            call cutman_mut(nl,cutopo,mut)
            if (persec) call pcutman_mut(nl,percutopo,mut)
          end if

          call prtest(nl,q,pr,iter)
          if (lomach) then
!---------- Darmofal's pressuse-gradient-based cut-off:
!           cutoff_type=3: cut-off on coarser grids restricted from fine
!           cutoff_type=2: cut-off definition used on all grid levels
            if (((cutoff_type.eq.3).and. &
                 (((nlevel.gt.1).and.(nl.eq.1).and.(to_coarse)).or. &
                  (nlevel.eq.1))).or. &
                (cutoff_type.eq.2)) then
                call lm_cutoff_pg(nl,q,cutoff_pgr)
            end if
            if ((visprec).and. &
                (((nlevel.gt.1).and.(nl.eq.1).and.(to_coarse)).or. &
                 (nlevel.eq.1))) then
                call lm_viscoff(nl,q,si,sj,sk,cutoff_vis)
            end if
          end if

          if ((iter/iupdt*iupdt.eq.iter).and.(irk.eq.1)) then
!---------- local time stepping for all other cases
            call lambda(nl,iter,q,si,sj,sk,xdot,ydot,zdot,vol,mut,dq, &
                        qp,cutoff_pgr,cutoff_vis)
            call deltat(nl,iter,dtvol,dtvolt,qp,dq,betaxi,betaeta, &
                        betazeta,tbetaxi,tbetaeta,tbetazeta)
          end if
!msc: 03/06/2015: here we zero also the first auxiliary cell values of 
!                 res since zero values on the external bounaries are
!                 needed by routine cirs. This is to override action
!                 of bc_res called by restrict.
          call zero(1,npde,0,2*nharms, 3,nl,res, &
                    1,npde,0,2*nharms, 1)

          
          if (viscous) then
            calmut = .false.
            call q_edges(nl,q,mut)
          end if

          idir = 1
          call muscl(idir,nl,bctopo,q,qp,qm,dq,dqp,dqm)
          if (lomach) then
            call lm_cutoff(idir,nl,dqm,cutoff_pgr,cutoff_vis)
            call proflux(idir,nl,qp,qm,dqp,dqm,flux,si)
          else
            call roflux(idir,nl,qp,qm,dqp,dqm,flux,si)
          end if
          call resid(idir,nl,flux,res)
          if (viscous) then
            if (kom.or.kom_bsl.or.kom_sst) then
              calmut = .true.
            else
              calmut = .false.
            end if
            call q_face(idir,nl,bctopo,q,dq,mut,flux)
            call q_der(idir,nl,bctopo,q,qp,qm,dqm)
            if (kom.or.kom_bsl.or.kom_sst) then
              call rtst(idir,nl,dq,qp,qm,dqm,xideri,etaderi,zetaderi, &
                        dqp,flux)
            end if
            call vflux(idir,nl,dq,qp,qm,dqm,dqp,flux,fidist,xideri, &
                       etaderi,zetaderi,si)
            call resid(idir,nl,flux,res)
          end if

          idir = 2
          call muscl(idir,nl,bctopo,q,qp,qm,dq,dqp,dqm)
          if (lomach) then
            call lm_cutoff(idir,nl,dqm,cutoff_pgr,cutoff_vis)
            call proflux(idir,nl,qp,qm,dqp,dqm,flux,sj)
          else
            call roflux(idir,nl,qp,qm,dqp,dqm,flux,sj)
          end if
          call resid(idir,nl,flux,res)
          if (viscous) then
            if (kom.or.kom_bsl.or.kom_sst) then
              calmut = .true.
            else
              calmut = .false.
            end if
            call q_face(idir,nl,bctopo,q,dq,mut,flux)
            call q_der(idir,nl,bctopo,q,qp,qm,dqm)
            if (kom.or.kom_bsl.or.kom_sst) then
              call rtst(idir,nl,dq,qp,qm,dqm,xiderj,etaderj,zetaderj, &
                        dqp,flux)
            end if
            call vflux(idir,nl,dq,qp,qm,dqm,dqp,flux,fjdist,xiderj, &
                       etaderj,zetaderj,sj)
            call resid(idir,nl,flux,res)
          end if

          idir = 3
          call muscl(idir,nl,bctopo,q,qp,qm,dq,dqp,dqm)
          if (lomach) then
            call lm_cutoff(idir,nl,dqm,cutoff_pgr,cutoff_vis)
            call proflux(idir,nl,qp,qm,dqp,dqm,flux,sk)
          else
            call roflux(idir,nl,qp,qm,dqp,dqm,flux,sk)
          end if
          call resid(idir,nl,flux,res)
          if (viscous) then
            if (kom.or.kom_bsl.or.kom_sst) then
              calmut = .true.
            else
              calmut = .false.
            end if
            call q_face(idir,nl,bctopo,q,dq,mut,flux)
            call q_der(idir,nl,bctopo,q,qp,qm,dqm)
            if (kom.or.kom_bsl.or.kom_sst) then
              call rtst(idir,nl,dq,qp,qm,dqm,xiderk,etaderk,zetaderk, &
                        dqp,flux)
            end if
            call vflux(idir,nl,dq,qp,qm,dqm,dqp,flux,fkdist,xiderk, &
                       etaderk,zetaderk,sk)
            call resid(idir,nl,flux,res)
          end if

          if (kom.or.kom_bsl.or.kom_sst) then
            calmut = .false.
            call q_face(1,nl,bctopo,q,dq,mut,flux)
            call q_der_c(1,nl,dq,qp)
            call q_face(2,nl,bctopo,q,dq,mut,flux)
            call q_der_c(2,nl,dq,qm)
            call q_face(3,nl,bctopo,q,dq,mut,flux)
            call q_der_c(3,nl,dq,dqm)
!---------- msc, 21/07/2021. If prd_all_lvls is true, it calculates
!             production term on all grid levels and calculates muT on
!             coarser grids restricting muT from fine grid (nl=1).
            if ((nl.eq.1).or.prd_all_lvls) then
              call turprod(nl,q,qp,qm,dqm,prd,divvel,delpls,xideri, &
                etaderi,zetaderi,xiderj,etaderj,zetaderj,xiderk,etaderk, &
                zetaderk)
            end if
            if (kom) then
              call srckw(nl,q,mut,prd,divvel,res,vol)
            else if (kom_bsl.or.kom_sst) then
              call update_ttrm(nl,q,qp,qm,dqm,xideri,etaderi,zetaderi, &
                xiderj,etaderj,zetaderj,xiderk,etaderk,zetaderk,ttrm, &
                dist)
              call srckw_sst(nl,q,prd,divvel,res,mut,ttrm,vol)
            end if
          end if

          if (relframe) then
            call src_relf(nl,q,vol,res)
          end if

          call prim2cons(nl,q,dq)

          if (dualt) then
!-------- unsteady source term (absent in predict./correct.)
            if (.not.rkimp.or. &
                (rkimp.and.(irsop.gt.0.or.nlevel.gt.1))) then
              call src_uns(nl,q,vol,res)
            end if
          else if (harbal) then
!-------- unsteady source term
            if (.not.rkimp.or. &
                (rkimp.and.(irsop.gt.0.or.nlevel.gt.1))) then
              call src_hb(nl,q,vol,res,dq)
            end if
          end if

!-------- compute RMS of residual for hist.dat
          if ((nl.eq.1).and.(irk.eq.1).and.(iter.ne.niter+1).and. &
              ((mgit.eq.1).or.(mgit/iprint*iprint.eq.mgit).or. &
               (mgit.eq.ncycle))) then
            if (.not.harbal) then
              call resg_rms(nl,q,res,rhs,vol)
            else if (harbal) then
              call copy_array(1,npde,0,2*nharms, 3,nl,work1,res, &
                              1,npde,0,2*nharms,-1)
              call resg_rms_hb(nl,q,work1,vol)
            end if
          end if

!-------- compute OVERALL UNPRECONDITIONED residual on finest grid
          if ((nl.eq.1).and.(dualt)) then
            call add_array(1,npde,0,2*nharms, 3,nl,res,rhs, &
                           1,npde,0,2*nharms,-1,-1.d0)
          end if

!-------- apply low-speed preconditioning
          if (lomach) then
            call lsp(nl,q,dqp,prec,preci,ipkpi,cutoff_pgr,cutoff_vis, &
              dtvol,dtvolt,vol,ttrm,delpls,xdot,ydot,zdot,alphas(irk))
            call precres(nl,res,prec,work1)
          end if

!-------- correct multigrid forcing term on coarser grids (nl>1)
          if (from_fine.and.(iter.eq.1).and.(irk.eq.1)) then
            call add_array(1,npde,0,2*nharms, 3,nl,rhs,res, &
                           1,npde,0,2*nharms,-1,1.d0)
          end if

!-------- compute OVERALL PRECONDITIONED residual on coarser grids
          if (nl.gt.1) then
            call add_array(1,npde,0,2*nharms, 3,nl,res,rhs, &
                           1,npde,0,2*nharms,-1,-1.d0)
          end if

!-------- no update; restrict PRECONDITIONED and SMOOTHED OVERALL 
!         residual to next coarser mesh
          if (iter .eq. niter+1) then
            if (kom.or.kom_bsl.or.kom_sst) then
              if ((nlevel.gt.1).and.(nl.ne.nlevel).and. &
                  (.not.prd_all_lvls).and.lim_prodtke)  then
                call lim_prod_mg(nl,q,prd,mut)
              end if
              call lim_res_ome(nl,res,q,qold,vol,dtvolt,prec,ipkpi, &
                               work1,delpls,prd,ttrm,dqp,alphas(irk))
            end if
            return
          endif

!-------- smooth residual
          if ( ((irsop.gt.0).and.(irsop.le.3)).or. &
               ((irsop.gt.4).and.(irsop.le.6)) ) then
!msc: 04/06/2015: by calling cutman_q we apply cut condition based on
!                 initial residuals at first smoothing step. Commenting
!                 out following call results in using homogeneous
!                 Dirichlet BCs at first smoothing step
            if (cutcirs.eq.1) then
              call scalef(nl,res,dtvol,dtvolt)
              call copy_array(1,npde,0,2*nharms, 3,nl,work1,res, &
                              1,npde,0,2*nharms, 1)
              do irs_its = 1,2
                call cutman_q(nl,cutopo,res)
                call ccirs(nl,res,work1,betaxi,betaeta,betazeta, &
                           tbetaxi,tbetaeta,tbetazeta,work2)
              end do
              call scaleb(nl,res,dtvol,dtvolt)
            else if (cutcirs.eq.0) then
              call cirs(nl,res,dtvol,betaxi,betaeta,betazeta,dtvolt, &
                        tbetaxi,tbetaeta,tbetazeta,work2)
            end if
         
          end if

          call update(nl,res,q,qold,vol,dtvol,dtvolt,prec,ipkpi,work1, &
                      dq,dqp,delpls,prd,ttrm,alphas(irk))
          if (kom.or.kom_bsl.or.kom_sst) then
            if (ifixq.eq.1) call fixq(nl,q)
            if ((irk.eq.nstage).and. &
                ((nl.eq.1).or.(.not.prd_all_lvls))) then
              if (kom.or.kom_bsl) then
                call tur_vis(q,mut,nl)
              else if (kom_sst) then
                call tur_vis_sst(q,mut,qp,qm,dqm,xideri,etaderi, &
                  zetaderi,xiderj,etaderj,zetaderj,xiderk,etaderk, &
                  zetaderk,dist,nl)
              end if
              call cutman_mut(nl,cutopo,mut)
              if (persec) call pcutman_mut(nl,percutopo,mut)
            end if
          end if

!------ increment RK stage
        end do    

        errmax =dmax1(resrms(1),resrms(2),resrms(3),resrms(4),resrms(5))

        if (errmax.lt.toler) then
          done = .true.
          return
        end if

!---- increment RK cycle
      end do    

      return
      end

!---------------------------------------------------------------------
      subroutine ccirs(nl,res,ores,betaxi,betaeta,betazeta, &
                       tbetaxi,tbetaeta,tbetazeta,work)
!---------------------------------------------------------------------       

      use common_variables
      use cosa_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,iblock,imax,jmax,kmax,ijkmax,ires,idtv,iwrk
      real (kind=cosa_real) ores(*),res(*),betaxi(*),betaeta(*),betazeta(*), &
           tbetaxi(*),tbetaeta(*),tbetazeta(*),work(*)

      do iblock = 1,mynblocks
        imax   = i_imax     (iblock,nl)
        jmax   = j_jmax     (iblock,nl)
        kmax   = k_kmax     (iblock,nl)
        ijkmax = ijk_ijkmax (iblock,nl)
        ires   = 1 + off_p3 (iblock,nl) * npde * dim5
        idtv   = 1 + off_p1 (iblock,nl) *        dim5
        iwrk   = 1 + off_1d (iblock,nl) * npde
        call ccirs_b(res(ires),ores(ires),betaxi(idtv),betaeta(idtv), &
                betazeta(idtv),tbetaxi(idtv),tbetaeta(idtv), &
                tbetazeta(idtv),work(iwrk),imax,jmax,kmax,npde,nharms, &
                ijkmax)
      end do

      return
      end

!---------------------------------------------------------------------
      subroutine ccirs_b(res,ores,betaxi,betaeta,betazeta,tbetaxi, &
                    tbetaeta,tbetazeta,work,imax,jmax,kmax,npde,nharms, &
                    ijkmax)
!---------------------------------------------------------------------       

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,ijkmax,npde,nharms
      integer(kind=cosa_int) im,i,jm,j,km,k,ipde,n
      real (kind=cosa_real) &
           res      (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           ores     (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           betaxi   ( 0:imax  , 0:jmax  , 0:kmax  ,     0:2*nharms), &
           betaeta  ( 0:imax  , 0:jmax  , 0:kmax  ,     0:2*nharms), &
           betazeta ( 0:imax  , 0:jmax  , 0:kmax  ,     0:2*nharms), &
           tbetaxi  ( 0:imax  , 0:jmax  , 0:kmax  ,     0:2*nharms), &
           tbetaeta ( 0:imax  , 0:jmax  , 0:kmax  ,     0:2*nharms), &
           tbetazeta( 0:imax  , 0:jmax  , 0:kmax  ,     0:2*nharms), &
           work     ( ijkmax-1, npde)

      real (kind=cosa_real) epsirc

      im = imax-1
      jm = jmax-1
      km = kmax-1

      if (irsop.eq.1) then

        epsirc = 0.25d0 * (rat**2 - 1.d0)

!------ xi direction

        do n = 0,2*nharms
          do ipde = 1,npde
            do k = 1,km
              do j = 1,jm
                work(1,1) = 0
                work(1,2) = 1 + 2*epsirc
                work(1,3) = - epsirc
                work(1,4) = ores(1,j,k,ipde,n) + &
                  epsirc*res(0,j,k,ipde,n)
                do i = 2,im-1
                  work(i,1) = - epsirc
                  work(i,2) = 1 + 2*epsirc
                  work(i,3) = - epsirc
                  work(i,4) = ores(i,j,k,ipde,n)
                end do
                work(im,1) = - epsirc
                work(im,2) = 1 + 2*epsirc
                work(im,3) = 0
                work(im,4) = ores(im,j,k,ipde,n) + &
                  epsirc*res(imax,j,k,ipde,n)
                call tridi(work(1,1),work(1,2),work(1,3),work(1,4), &
                           ijkmax-1,im)
                do i=1,im
                  res(i,j,k,ipde,n) = work(i,4)
                end do
              end do
            end do
          end do
        end do

!------ eta direction

        do n = 0,2*nharms
          do ipde = 1,npde
            do k = 1,km
              do i = 1,im
                work(1,1) = 0
                work(1,2) = 1 + 2*epsirc
                work(1,3) = - epsirc
                work(1,4) = res(i,1,k,ipde,n) + &
                  epsirc*res(i,0,k,ipde,n)
                do j = 2,jm-1
                  work(j,1) = - epsirc
                  work(j,2) = 1 + 2*epsirc
                  work(j,3) = - epsirc
                  work(j,4) = res(i,j,k,ipde,n)
                end do
                work(jm,1) = - epsirc
                work(jm,2) = 1 + 2*epsirc
                work(jm,3) = 0
                work(jm,4) = res(i,jm,k,ipde,n) + &
                  epsirc*res(i,jmax,k,ipde,n)
                call tridi(work(1,1),work(1,2),work(1,3),work(1,4), &
                           ijkmax-1,jm)
                do j=1,jm
                  res(i,j,k,ipde,n) = work(j,4)
                end do
              end do
            end do
          end do
        end do

!------ zeta direction

        do n = 0,2*nharms
          do ipde = 1,npde
            do j = 1,jm
              do i = 1,im
                work(1,1) = 0
                work(1,2) = 1 + 2*epsirc
                work(1,3) = - epsirc
                work(1,4) = res(i,j,1,ipde,n) + &
                  epsirc*res(i,j,0,ipde,n)
                do k = 2,km-1
                  work(k,1) = - epsirc
                  work(k,2) = 1 + 2*epsirc
                  work(k,3) = - epsirc
                  work(k,4) = res(i,j,k,ipde,n)
                end do
                work(km,1) = - epsirc
                work(km,2) = 1 + 2*epsirc
                work(km,3) = 0
                work(km,4) = res(i,j,km,ipde,n) + &
                  epsirc*res(i,j,kmax,ipde,n)
                call tridi(work(1,1),work(1,2),work(1,3),work(1,4), &
                           ijkmax-1,km)
                do k=1,km
                  res(i,j,k,ipde,n) = work(k,4)
                end do
              end do
            end do
          end do
        end do

      else if (irsop.ge.2) then

!------ xi direction

        do n = 0,2*nharms
          do ipde = 1,5
            do k = 1,km
              do j = 1,jm
                epsirc = betaxi(1,j,k,n)
                work(1,1) = 0
                work(1,2) = 1 + 2*epsirc
                work(1,3) = - epsirc
                work(1,4) = ores(1,j,k,ipde,n) + &
                  epsirc*res(0,j,k,ipde,n)
                do i = 2,im-1
                  epsirc = betaxi(i,j,k,n)
                  work(i,1) = - epsirc
                  work(i,2) = 1 + 2*epsirc
                  work(i,3) = - epsirc
                  work(i,4) = ores(i,j,k,ipde,n)
                end do
                epsirc = betaxi(im,j,k,n)
                work(im,1) = - epsirc
                work(im,2) = 1 + 2*epsirc
                work(im,3) = 0
                work(im,4) = ores(im,j,k,ipde,n) + &
                  epsirc*res(imax,j,k,ipde,n)
                call tridi(work(1,1),work(1,2),work(1,3),work(1,4), &
                           ijkmax-1,im)
                do i=1,im
                  res(i,j,k,ipde,n) = work(i,4)
                end do
              end do
            end do
          end do
        end do

!------ eta direction

        do n = 0,2*nharms
          do ipde = 1,5
            do k = 1,km
              do i = 1,im
                epsirc = betaeta(i,1,k,n)
                work(1,1) = 0
                work(1,2) = 1 + 2*epsirc
                work(1,3) = - epsirc
                work(1,4) = res(i,1,k,ipde,n) + &
                  epsirc*res(i,0,k,ipde,n)
                do j = 2,jm-1
                  epsirc = betaeta(i,j,k,n)
                  work(j,1) = - epsirc
                  work(j,2) = 1 + 2*epsirc
                  work(j,3) = - epsirc
                  work(j,4) = res(i,j,k,ipde,n)
                end do
                epsirc = betaeta(i,jm,k,n)
                work(jm,1) = - epsirc
                work(jm,2) = 1 + 2*epsirc
                work(jm,3) = 0
                work(jm,4) = res(i,jm,k,ipde,n) + &
                  epsirc*res(i,jmax,k,ipde,n)
                call tridi(work(1,1),work(1,2),work(1,3),work(1,4), &
                           ijkmax-1,jm)
                do j=1,jm
                  res(i,j,k,ipde,n) = work(j,4)
                end do
              end do
            end do
          end do
        end do

!------ zeta direction

        do n = 0,2*nharms
          do ipde = 1,5
            do j = 1,jm
              do i = 1,im
                epsirc = betazeta(i,j,1,n)
                work(1,1) = 0
                work(1,2) = 1 + 2*epsirc
                work(1,3) = - epsirc
                work(1,4) = res(i,j,1,ipde,n) + &
                  epsirc*res(i,j,0,ipde,n)
                do k = 2,km-1
                  epsirc = betazeta(i,j,k,n)
                  work(k,1) = - epsirc
                  work(k,2) = 1 + 2*epsirc
                  work(k,3) = - epsirc
                  work(k,4) = res(i,j,k,ipde,n)
                end do
                epsirc = betazeta(i,j,km,n)
                work(km,1) = - epsirc
                work(km,2) = 1 + 2*epsirc
                work(km,3) = 0
                work(km,4) = res(i,j,km,ipde,n) + &
                  epsirc*res(i,j,kmax,ipde,n)
                call tridi(work(1,1),work(1,2),work(1,3),work(1,4), &
                           ijkmax-1,km)
                do k=1,km
                  res(i,j,k,ipde,n) = work(k,4)
                end do
              end do
            end do
          end do
        end do

        if (kom.or.kom_bsl.or.kom_sst) then

!-------- xi direction

          do n = 0,2*nharms
            do ipde = 6,npde
              do k = 1,km
                do j = 1,jm
                  epsirc = tbetaxi(1,j,k,n)
                  work(1,1) = 0
                  work(1,2) = 1 + 2*epsirc
                  work(1,3) = - epsirc
                  work(1,4) = ores(1,j,k,ipde,n) + &
                    epsirc*res(0,j,k,ipde,n)
                  do i = 2,im-1
                    epsirc = tbetaxi(i,j,k,n)
                    work(i,1) = - epsirc
                    work(i,2) = 1 + 2*epsirc
                    work(i,3) = - epsirc
                    work(i,4) = ores(i,j,k,ipde,n)
                  end do
                  epsirc = tbetaxi(im,j,k,n)
                  work(im,1) = - epsirc
                  work(im,2) = 1 + 2*epsirc
                  work(im,3) = 0
                  work(im,4) = ores(im,j,k,ipde,n) + &
                    epsirc*res(imax,j,k,ipde,n)
                  call tridi(work(1,1),work(1,2),work(1,3),work(1,4), &
                             ijkmax-1,im)
                  do i=1,im
                    res(i,j,k,ipde,n) = work(i,4)
                  end do
                end do
              end do
            end do
          end do

!-------- eta direction

          do n = 0,2*nharms
            do ipde = 6,npde
              do k = 1,km
                do i = 1,im
                  epsirc = tbetaeta(i,1,k,n)
                  work(1,1) = 0
                  work(1,2) = 1 + 2*epsirc
                  work(1,3) = - epsirc
                  work(1,4) = res(i,1,k,ipde,n) + &
                    epsirc*res(i,0,k,ipde,n)
                  do j = 2,jm-1
                    epsirc = tbetaeta(i,j,k,n)
                    work(j,1) = - epsirc
                    work(j,2) = 1 + 2*epsirc
                    work(j,3) = - epsirc
                    work(j,4) = res(i,j,k,ipde,n)
                  end do
                  epsirc = tbetaeta(i,jm,k,n)
                  work(jm,1) = - epsirc
                  work(jm,2) = 1 + 2*epsirc
                  work(jm,3) = 0
                  work(jm,4) = res(i,jm,k,ipde,n) + &
                    epsirc*res(i,jmax,k,ipde,n)
                  call tridi(work(1,1),work(1,2),work(1,3),work(1,4), &
                             ijkmax-1,jm)
                  do j=1,jm
                    res(i,j,k,ipde,n) = work(j,4)
                  end do
                end do
              end do
            end do
          end do

!-------- zeta direction

          do n = 0,2*nharms
            do ipde = 6,npde
              do j = 1,jm
                do i = 1,im
                  epsirc = tbetazeta(i,j,1,n)
                  work(1,1) = 0
                  work(1,2) = 1 + 2*epsirc
                  work(1,3) = - epsirc
                  work(1,4) = res(i,j,1,ipde,n) + &
                    epsirc*res(i,j,0,ipde,n)
                  do k = 2,km-1
                    epsirc = tbetazeta(i,j,k,n)
                    work(k,1) = - epsirc
                    work(k,2) = 1 + 2*epsirc
                    work(k,3) = - epsirc
                    work(k,4) = res(i,j,k,ipde,n)
                  end do
                  epsirc = tbetazeta(i,j,km,n)
                  work(km,1) = - epsirc
                  work(km,2) = 1 + 2*epsirc
                  work(km,3) = 0
                  work(km,4) = res(i,j,km,ipde,n) + &
                    epsirc*res(i,j,kmax,ipde,n)
                  call tridi(work(1,1),work(1,2),work(1,3),work(1,4), &
                             ijkmax-1,km)
                  do k=1,km
                    res(i,j,k,ipde,n) = work(k,4)
                  end do
                end do
              end do
            end do
          end do

        end if

      end if

      return
      end

!---------------------------------------------------------------------
      subroutine tridi(a,b,c,x,dimx,dim)
!---------------------------------------------------------------------
        
      use cosa_precision

      implicit none
      
      integer(kind=cosa_int) dimx,dim,j
      real(kind=cosa_real) a(dimx),b(dimx),c(dimx),d(dimx),x(dimx),r(dimx), &
           t(dimx)
      
      r(1) = c(1) / b(1)
      t(1) = x(1) / b(1)
      do j = 2,dim-1
         r(j) = c(j) / (b(j) - r(j-1)*a(j))
         t(j) = (x(j) - t(j-1)*a(j)) / (b(j) - r(j-1)*a(j)) 
      end do
      t(dim) = (x(dim) - t(dim-1)*a(dim)) / (b(dim) - r(dim-1)*a(dim)) 
      
      x(dim) = t(dim)
      do j=dim-1,1,-1
         x(j) = t(j) - x(j+1) * r(j)
      end do
      
      return
    end subroutine tridi

!-----------------------------------------------------------------------
      subroutine precres(nl,res,prec,work1)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_variables
      use cosa_precision

      implicit none
!
      integer(kind=cosa_int) nl,imax,jmax,kmax
      integer(kind=cosa_int) iblock,ires,iwrk,iprec
      real (kind=cosa_real) res(*),prec(*),work1(*)

      do iblock = 1,mynblocks
        imax   = i_imax     (iblock,nl)
        jmax   = j_jmax     (iblock,nl)
        kmax   = k_kmax     (iblock,nl)
        ires   = 1 + off_p3 (iblock,nl  ) * npde * dim5
        iwrk   = 1 + off_p3 (iblock,nl  ) * npde * dim5
        iprec  = 1 + off_m1 (iblock,nl) * npde * npde * dim5
        call prec_bres(res(ires),prec(iprec),work1(iwrk),imax,jmax,kmax, &
                       npde,nharms)
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine prec_bres(res,prec,work1,imax,jmax,kmax,npde,nharms)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none
!
      integer(kind=cosa_int) imax,jmax,kmax,npde,nharms
      integer(kind=cosa_int) im,jm,km,i,j,k,ipde,ipde1,n
      real (kind=cosa_real) &
           res  (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,     0:2*nharms), &
           prec (   imax-1,   jmax-1,   kmax-1,npde,npde,0:2*nharms), &
           work1(-1:imax+1,-1:jmax+1,-1:kmax+1,npde,     0:2*nharms)

      im = imax-1
      jm = jmax-1
      km = kmax-1

      do n = 0,2*nharms
        do ipde = 1,npde
          do k = 1,km
            do j = 1,jm
              do i = 1,im
                work1(i,j,k,ipde,n) = 0.d0
                do ipde1 = 1,npde
                  work1(i,j,k,ipde,n) = work1(i,j,k,ipde,n) + &
                    prec(i,j,k,ipde,ipde1,n) * res(i,j,k,ipde1,n)
                end do
              end do
            end do
          end do
        enddo
      end do  

      do n = 0,2*nharms
        do ipde = 1,npde
          do k = 1,km
            do j = 1,jm
              do i = 1,im
                res(i,j,k,ipde,n) = work1(i,j,k,ipde,n)
              end do
            end do
          end do
        end do  
      end do  

      return
      end

!-----------------------------------------------------------------------
      subroutine deltat_urk(nl,vol,dtvol)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,iblock,imax,jmax,kmax,ivol,idtv
      real (kind=cosa_real) dtvol(*), vol(*)

      do iblock = 1,mynblocks
        imax   = i_imax     (iblock,nl)
        jmax   = j_jmax     (iblock,nl)
        kmax   = k_kmax     (iblock,nl)
        ivol   = 1 + off_p1 (iblock,nl)
        idtv   = 1 + off_p1 (iblock,nl)
        call bdeltat_urk(vol(ivol),dtvol(idtv),imax,jmax,kmax)
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine bdeltat_urk(vol,dtvol,imax,jmax,kmax)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax
      integer(kind=cosa_int) i,j,k
      real (kind=cosa_real) dtvol(0:imax,0:jmax,0:kmax), &
                      vol(0:imax,0:jmax,0:kmax)

      do k = 1,kmax-1
        do j = 1,jmax-1
          do i = 1,imax-1
            dtvol(i,j,k) = dt / vol(i,j,k)
          end do
        end do
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine coarsen(nl,x,y,z,dx,dy,dz,x2,y2,z2,dx2,dy2,dz2)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,iblock,imax,jmax,kmax,imax2,jmax2,kmax2,ixyz,ixyz2, &
                idxyz,idxyz2
      real (kind=cosa_real) x(*),y(*),z(*),x2(*),y2(*),z2(*),dx(*),dy(*),dz(*), &
           dx2(*),dy2(*),dz2(*)

      do iblock = 1,mynblocks
        imax    = i_imax  (iblock,nl-1)         
        jmax    = j_jmax  (iblock,nl-1)
        kmax    = k_kmax  (iblock,nl-1)
        imax2   = i_imax  (iblock,nl  )         
        jmax2   = j_jmax  (iblock,nl  )
        kmax2   = k_kmax  (iblock,nl  )
        ixyz    = 1 + off_p2(iblock,nl-1)
        ixyz2   = 1 + off_p2(iblock,nl  )
        idxyz   = 1 + off_p1(iblock,nl-1)
        idxyz2  = 1 + off_p1(iblock,nl  )
        call coarsen_b(x(ixyz),y(ixyz),z(ixyz),dx(idxyz),dy(idxyz), &
             dz(idxyz),imax,jmax,kmax,x2(ixyz2),y2(ixyz2),z2(ixyz2), &
             dx2(idxyz2),dy2(idxyz2),dz2(idxyz2),imax2,jmax2,kmax2)
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine coarsen_b(x,y,z,dx,dy,dz,imax,jmax,kmax,x2,y2,z2,dx2, &
                 dy2,dz2,imax2,jmax2,kmax2)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,imax2,jmax2,kmax2
      integer(kind=cosa_int) i,j,k,i2,j2,k2
      real (kind=cosa_real) &
           x  (0:imax+1 ,0:jmax+1 ,0:kmax+1 ), &
           y  (0:imax+1 ,0:jmax+1 ,0:kmax+1 ), &
           z  (0:imax+1 ,0:jmax+1 ,0:kmax+1 ), &
           x2 (0:imax2+1,0:jmax2+1,0:kmax2+1), &
           y2 (0:imax2+1,0:jmax2+1,0:kmax2+1), &
           z2 (0:imax2+1,0:jmax2+1,0:kmax2+1), &
           dx (0:imax   ,0:jmax   ,0:kmax   ), &
           dy (0:imax   ,0:jmax   ,0:kmax   ), &
           dz (0:imax   ,0:jmax   ,0:kmax   ), &
           dx2(0:imax2  ,0:jmax2  ,0:kmax2  ), &
           dy2(0:imax2  ,0:jmax2  ,0:kmax2  ), &
           dz2(0:imax2  ,0:jmax2  ,0:kmax2  )

      do k=1,kmax,2
        k2 = 1+(k-1)/2
        do j=1,jmax,2
          j2 = 1+(j-1)/2
          do i=1,imax,2
            i2 = 1+(i-1)/2

            x2(i2,j2,k2) = x(i,j,k)
            y2(i2,j2,k2) = y(i,j,k)
            z2(i2,j2,k2) = z(i,j,k)

          end do
        end do
      end do

      if (moving.and.pitching) then

        do k=1,kmax,2
          k2 = 1+(k-1)/2
          do j=1,jmax,2
            j2 = 1+(j-1)/2
            do i=1,imax,2
              i2 = 1+(i-1)/2

              dx2(i2,j2,k2) = dx(i,j,k)
              dy2(i2,j2,k2) = dy(i,j,k)
              dz2(i2,j2,k2) = dz(i,j,k)

            end do
          end do
        end do

      end if

      return
      end

!-----------------------------------------------------------------------
      subroutine restrict(nl)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl, cutopo1, percutopo1
      real (kind=cosa_real) q1, q2, res1, rhs1, rhs2, vol1, vol2, cutoff_pgr1, &
        cutoff_pgr2, cutoff_vis1, cutoff_vis2, prd1, prd2, dpl1, dpl2, &
        div1, div2, mut1, mut2, bctopo1, mut
      pointer (pq1  ,  q1), (pq2  ,  q2), (pvol1,vol1), (pvol2,vol2), &
              (prhs1,rhs1), (prhs2,rhs2), (pres1,res1), (pprd1,prd1), &
              (pprd2,prd2), (pdpl1,dpl1), (pdpl2,dpl2), (pdiv1,div1), &
              (pdiv2,div2), (pmut1,mut1), (pmut2,mut2), &
              (pcutoff_pgr1,cutoff_pgr1),   (pcutoff_pgr2,cutoff_pgr2), &
              (pcutoff_vis1,cutoff_vis1),   (pcutoff_vis2,cutoff_vis2), &
              (pcutopo1,cutopo1),           (pbctopo1,bctopo1), &
              (ppercutopo1,percutopo1)

!     res1 : residuals of FINER grid
!     rhs2 : residuals of FINER grid restricted to COARSER grid. This is
!           the RHS of the residual equation to be solved on COARSER grid

      pq1   = p_q(nl)
      pres1 = p_res(nl)
      prhs1 = p_rhs(nl)
      pvol1 = p_vol(nl)
      pprd1 = p_prd(nl)
      pdpl1 = p_delpls(nl)
      pdiv1 = p_divvel(nl)
      pmut1 = p_mut(nl)
      pcutoff_pgr1 = p_cutoff_pgr(nl)
      pcutoff_vis1 = p_cutoff_vis(nl)
      pcutopo1     = p_cutopo(nl)
      ppercutopo1  = p_percutopo(nl)
      pbctopo1     = p_bctopo(nl)

      pq2   = p_q(nl+1)
      prhs2 = p_rhs(nl+1)
      pvol2 = p_vol(nl+1)
      pprd2 = p_prd(nl+1)
      pdpl2 = p_delpls(nl+1)
      pdiv2 = p_divvel(nl+1)
      pmut2 = p_mut(nl+1)
      pcutoff_pgr2 = p_cutoff_pgr(nl+1)
      pcutoff_vis2 = p_cutoff_vis(nl+1)

      call rest1(nl,q1,q2,vol1,vol2)

      if (lomach) then
        if (cutoff_type.eq.3) then
          call rest3(nl,cutoff_pgr1,cutoff_pgr2,vol1,vol2)
        end if
        if (visprec) then
          call rest3(nl,cutoff_vis1,cutoff_vis2,vol1,vol2)
        end if
      end if

      if (kom.or.kom_bsl.or.kom_sst) then
!------ msc, 21/07/2021. If prd_all_lvls is true, it calculates
        if (.not.prd_all_lvls) then
          call rest4(nl,prd1,prd2,vol1,vol2)
          call rest4(nl,dpl1,dpl2,vol1,vol2)
          call rest4(nl,div1,div2,vol1,vol2)
        else
          call rest4(nl,mut1,mut2,vol1,vol2)
        end if
        if (ho_rest) then
          call cutman_q(nl,cutopo1,res1)
          if (persec) call pcutman_q(nl,percutopo1,res1)
          call bc_res(nl,bctopo1,res1)
          calmut = .false.
          call q_edges(nl,res1,mut)
          call q_corns(nl,res1,mut)
        end if
      end if
      call rest2(nl,rhs1,rhs2,res1,ho_rest)

      return
      end

!-----------------------------------------------------------------------
      subroutine rest1(nl,q1,q2,vol1,vol2)
!-----------------------------------------------------------------------

      use cosa_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,iblock,imax1,jmax1,kmax1,imax2,jmax2,kmax2,iq1,iq2, &
                ivol1,ivol2
      real(kind=cosa_real) q1(*),vol1(*),q2(*),vol2(*)

      do iblock = 1,mynblocks
        imax1  = i_imax(iblock,nl)
        jmax1  = j_jmax(iblock,nl)
        kmax1  = k_kmax(iblock,nl)
        imax2  = i_imax(iblock,nl+1)
        jmax2  = j_jmax(iblock,nl+1)
        kmax2  = k_kmax(iblock,nl+1)
        iq1    = 1 + off_p3 (iblock,nl  ) * npde * dim5
        iq2    = 1 + off_p3 (iblock,nl+1) * npde * dim5
        ivol1  = 1 + off_p1 (iblock,nl  )
        ivol2  = 1 + off_p1 (iblock,nl+1)
        call rest1_b(imax1,jmax1,kmax1,imax2,jmax2,kmax2,npde,nharms, &
                     q1(iq1),q2(iq2),vol1(ivol1),vol2(ivol2))
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine rest1_b(imax1,jmax1,kmax1,imax2,jmax2,kmax2,npde, &
                       nharms,q1,q2,vol1,vol2)
!-----------------------------------------------------------------------

      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax1,jmax1,kmax1,imax2,jmax2,kmax2,npde,nharms
      integer(kind=cosa_int) i1,j1,k1,i2,j2,k2,ipde,n
      real(kind=cosa_real) &
          q1(-1:imax1+1,-1:jmax1+1,-1:kmax1+1,npde,0:2*nharms), &
          q2(-1:imax2+1,-1:jmax2+1,-1:kmax2+1,npde,0:2*nharms), &
          vol1(0:imax1,0:jmax1,0:kmax1), &
          vol2(0:imax2,0:jmax2,0:kmax2)

!     q1 : solution on FINER grid
!     q2 : solution on FINER grid restricted to COARSER grid.

      do n = 0,2*nharms
        do k1=1,kmax1-2,2
          do j1=1,jmax1-2,2
            do i1=1,imax1-2,2
              i2 = 1 + i1/2
              j2 = 1 + j1/2
              k2 = 1 + k1/2
              do ipde=1,npde
                q2(i2,j2,k2,ipde,n) = &
                  ( vol1(i1  ,j1  ,k1  ) * q1(i1  ,j1  ,k1  ,ipde,n) + &
                    vol1(i1  ,j1+1,k1  ) * q1(i1  ,j1+1,k1  ,ipde,n) + &
                    vol1(i1+1,j1  ,k1  ) * q1(i1+1,j1  ,k1  ,ipde,n) + &
                    vol1(i1+1,j1+1,k1  ) * q1(i1+1,j1+1,k1  ,ipde,n) + &
                    vol1(i1  ,j1  ,k1+1) * q1(i1  ,j1  ,k1+1,ipde,n) + &
                    vol1(i1  ,j1+1,k1+1) * q1(i1  ,j1+1,k1+1,ipde,n) + &
                    vol1(i1+1,j1  ,k1+1) * q1(i1+1,j1  ,k1+1,ipde,n) + &
                    vol1(i1+1,j1+1,k1+1) * q1(i1+1,j1+1,k1+1,ipde,n) ) / &
                  vol2(i2,j2,k2)
              end do
            end do
          end do
        end do
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine rest2(nl,rhs1,rhs2,res1,ho_rest)
!-----------------------------------------------------------------------

      use cosa_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,iblock,imax1,jmax1,kmax1,imax2,jmax2,kmax2,ires1, &
                irhs1,irhs2
      real(kind=cosa_real) rhs1(*),rhs2(*),res1(*)
      logical ho_rest

      do iblock = 1,mynblocks
        imax1  = i_imax(iblock,nl)
        jmax1  = j_jmax(iblock,nl)
        kmax1  = k_kmax(iblock,nl)
        imax2  = i_imax(iblock,nl+1)
        jmax2  = j_jmax(iblock,nl+1)
        kmax2  = k_kmax(iblock,nl+1)
        irhs1  = 1 + off_p3 (iblock,nl  ) * npde * dim5
        ires1  = 1 + off_p3 (iblock,nl  ) * npde * dim5
        irhs2  = 1 + off_p3 (iblock,nl+1) * npde * dim5
        if (ho_rest) then
          call rest2_ho_b(imax1,jmax1,kmax1,imax2,jmax2,kmax2,npde, &
                          nharms,rhs1(irhs1),rhs2(irhs2),res1(ires1))
        else
          call rest2_lo_b(imax1,jmax1,kmax1,imax2,jmax2,kmax2,npde, &
                          nharms,rhs1(irhs1),rhs2(irhs2),res1(ires1))
        end if
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine rest2_ho_b(imax1,jmax1,kmax1,imax2,jmax2,kmax2,npde, &
                            nharms,rhs1,rhs2,res1)
!-----------------------------------------------------------------------

      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax1,jmax1,kmax1,imax2,jmax2,kmax2,npde,nharms
      integer(kind=cosa_int) i1,j1,k1,i2,j2,k2,ipde,n
      real(kind=cosa_real) &
          rhs1(-1:imax1+1,-1:jmax1+1,-1:kmax1+1,npde,0:2*nharms), &
          rhs2(-1:imax2+1,-1:jmax2+1,-1:kmax2+1,npde,0:2*nharms), &
          res1(-1:imax1+1,-1:jmax1+1,-1:kmax1+1,npde,0:2*nharms)

!     rhs1 :
!     res1 :
!     rhs2 :

!cc   do n = 0,2*nharms
!cc     do ipde=1,npde
!cc       do k2 = 1,kmax2-1
!cc         do j2 = 1,jmax2-1
!cc           do i2 = 1,imax2-1
!cc             rhs2(i2,j2,k2,ipde,n) = 0
!cc           end do
!cc         end do
!cc       end do
!cc     end do
!cc   end do

      do n = 0,2*nharms
        do ipde=1,npde
          do k1=1,kmax1-2,2
            do j1=1,jmax1-2,2
              do i1=1,imax1-2,2
                i2 = 1 + i1/2
                j2 = 1 + j1/2
                k2 = 1 + k1/2
!cc   rhs2(i2,j2,k2,ipde,n) = rhs2(i2,j2,k2,ipde,n) -
      rhs2(i2,j2,k2,ipde,n) = - &
        (27*(res1(i1  ,j1  ,k1  ,ipde,n)+res1(i1  ,j1+1,k1  ,ipde,n) + &
             res1(i1+1,j1  ,k1  ,ipde,n)+res1(i1+1,j1+1,k1  ,ipde,n) + &
             res1(i1  ,j1  ,k1+1,ipde,n)+res1(i1  ,j1+1,k1+1,ipde,n) + &
             res1(i1+1,j1  ,k1+1,ipde,n)+res1(i1+1,j1+1,k1+1,ipde,n))+ &
          9*(res1(i1-1,j1  ,k1  ,ipde,n)+res1(i1-1,j1+1,k1  ,ipde,n) + &
             res1(i1+2,j1  ,k1  ,ipde,n)+res1(i1+2,j1+1,k1  ,ipde,n) + &
             res1(i1-1,j1  ,k1+1,ipde,n)+res1(i1-1,j1+1,k1+1,ipde,n) + &
             res1(i1+2,j1  ,k1+1,ipde,n)+res1(i1+2,j1+1,k1+1,ipde,n) + &
             res1(i1  ,j1-1,k1  ,ipde,n)+res1(i1+1,j1-1,k1  ,ipde,n) + &
             res1(i1  ,j1+2,k1  ,ipde,n)+res1(i1+1,j1+2,k1  ,ipde,n) + &
             res1(i1  ,j1-1,k1+1,ipde,n)+res1(i1+1,j1-1,k1+1,ipde,n) + &
             res1(i1  ,j1+2,k1+1,ipde,n)+res1(i1+1,j1+2,k1+1,ipde,n) + &
             res1(i1  ,j1  ,k1-1,ipde,n)+res1(i1  ,j1+1,k1-1,ipde,n) + &
             res1(i1  ,j1  ,k1+2,ipde,n)+res1(i1  ,j1+1,k1+2,ipde,n) + &
             res1(i1+1,j1  ,k1-1,ipde,n)+res1(i1+1,j1+1,k1-1,ipde,n) + &
             res1(i1+1,j1  ,k1+2,ipde,n)+res1(i1+1,j1+1,k1+2,ipde,n))+ &
          3*(res1(i1  ,j1-1,k1-1,ipde,n)+res1(i1+1,j1-1,k1-1,ipde,n) + &
             res1(i1  ,j1+2,k1-1,ipde,n)+res1(i1+1,j1+2,k1-1,ipde,n) + &
             res1(i1-1,j1  ,k1-1,ipde,n)+res1(i1-1,j1+1,k1-1,ipde,n) + &
             res1(i1+2,j1  ,k1-1,ipde,n)+res1(i1+2,j1+1,k1-1,ipde,n) + &
             res1(i1  ,j1-1,k1+2,ipde,n)+res1(i1+1,j1-1,k1+2,ipde,n) + &
             res1(i1  ,j1+2,k1+2,ipde,n)+res1(i1+1,j1+2,k1+2,ipde,n) + &
             res1(i1-1,j1  ,k1+2,ipde,n)+res1(i1-1,j1+1,k1+2,ipde,n) + &
             res1(i1+2,j1  ,k1+2,ipde,n)+res1(i1+2,j1+1,k1+2,ipde,n) + &
             res1(i1-1,j1-1,k1  ,ipde,n)+res1(i1-1,j1-1,k1+1,ipde,n) + &
             res1(i1-1,j1+2,k1  ,ipde,n)+res1(i1-1,j1+2,k1+1,ipde,n) + &
             res1(i1+2,j1-1,k1  ,ipde,n)+res1(i1+2,j1-1,k1+1,ipde,n) + &
             res1(i1+2,j1+2,k1  ,ipde,n)+res1(i1+2,j1+2,k1+1,ipde,n))+ &
            (res1(i1-1,j1-1,k1-1,ipde,n)+res1(i1+2,j1-1,k1-1,ipde,n) + &
             res1(i1-1,j1+2,k1-1,ipde,n)+res1(i1+2,j1+2,k1-1,ipde,n) + &
             res1(i1-1,j1-1,k1+2,ipde,n)+res1(i1+2,j1-1,k1+2,ipde,n) + &
             res1(i1-1,j1+2,k1+2,ipde,n)+res1(i1+2,j1+2,k1+2,ipde,n))) / &
        64
              end do
            end do
          end do
        end do
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine rest2_lo_b(imax1,jmax1,kmax1,imax2,jmax2,kmax2,npde, &
                            nharms,rhs1,rhs2,res1)
!-----------------------------------------------------------------------

      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax1,jmax1,kmax1,imax2,jmax2,kmax2,npde,nharms
      integer(kind=cosa_int) i1,j1,k1,i2,j2,k2,ipde,n
      real(kind=cosa_real) &
          rhs1(-1:imax1+1,-1:jmax1+1,-1:kmax1+1,npde,0:2*nharms), &
          rhs2(-1:imax2+1,-1:jmax2+1,-1:kmax2+1,npde,0:2*nharms), &
          res1(-1:imax1+1,-1:jmax1+1,-1:kmax1+1,npde,0:2*nharms)

!     rhs1 :
!     res1 :
!     rhs2 :

!cc   do n = 0,2*nharms
!cc     do ipde=1,npde
!cc       do k2 = 1,kmax2-1
!cc         do j2 = 1,jmax2-1
!cc           do i2 = 1,imax2-1
!cc             rhs2(i2,j2,k2,ipde,n) = 0
!cc           end do
!cc         end do
!cc       end do
!cc     end do
!cc   end do

      do n = 0,2*nharms
        do ipde=1,npde
          do k1=1,kmax1-2,2
            do j1=1,jmax1-2,2
              do i1=1,imax1-2,2
                i2 = 1 + i1/2
                j2 = 1 + j1/2
                k2 = 1 + k1/2
!cc             rhs2(i2,j2,k2,ipde,n) = rhs2(i2,j2,k2,ipde,n) -
                rhs2(i2,j2,k2,ipde,n) = - &
                  ( res1(i1  ,j1  ,k1  ,ipde,n) + &
                    res1(i1  ,j1+1,k1  ,ipde,n) + &
                    res1(i1+1,j1  ,k1  ,ipde,n) + &
                    res1(i1+1,j1+1,k1  ,ipde,n) + &
                    res1(i1  ,j1  ,k1+1,ipde,n) + &
                    res1(i1  ,j1+1,k1+1,ipde,n) + &
                    res1(i1+1,j1  ,k1+1,ipde,n) + &
                    res1(i1+1,j1+1,k1+1,ipde,n) )
              end do
            end do
          end do
        end do
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine rest3(nl,cutoff1,cutoff2,vol1,vol2)
!-----------------------------------------------------------------------

      use cosa_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,iblock,imax1,jmax1,kmax1,imax2,jmax2,kmax2,icof1, &
                icof2,ivol1,ivol2
      real(kind=cosa_real) cutoff1(*),cutoff2(*),vol1(*),vol2(*)

      do iblock = 1,mynblocks
        imax1  = i_imax(iblock,nl)
        jmax1  = j_jmax(iblock,nl)
        kmax1  = k_kmax(iblock,nl)
        imax2  = i_imax(iblock,nl+1)
        jmax2  = j_jmax(iblock,nl+1)
        kmax2  = k_kmax(iblock,nl+1)
        icof1  = 1 + off_m1 (iblock,nl  ) *        dim5
        icof2  = 1 + off_m1 (iblock,nl+1) *        dim5
        ivol1  = 1 + off_p1 (iblock,nl  )
        ivol2  = 1 + off_p1 (iblock,nl+1)
        call rest3_b(imax1,jmax1,kmax1,imax2,jmax2,kmax2,npde,nharms, &
                     cutoff1(icof1),cutoff2(icof2),vol1(ivol1), &
                     vol2(ivol2))
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine rest3_b(imax1,jmax1,kmax1,imax2,jmax2,kmax2,npde, &
                         nharms,cutoff1,cutoff2,vol1,vol2)
!-----------------------------------------------------------------------

      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax1,jmax1,kmax1,imax2,jmax2,kmax2,npde,nharms
      integer(kind=cosa_int) i1,j1,k1,i2,j2,k2,ipde,n
      real(kind=cosa_real) &
          cutoff1(imax1-1,jmax1-1,kmax1-1,0:2*nharms), &
          cutoff2(imax2-1,jmax2-1,kmax2-1,0:2*nharms), &
          vol1(0:imax1,0:jmax1,0:kmax1), &
          vol2(0:imax2,0:jmax2,0:kmax2)

!     cutoff1 : cutoff on FINER grid
!     cutoff2 : cutoff on FINER grid restricted to COARSER grid.

      do n = 0,2*nharms
        do k1=1,kmax1-2,2
          do j1=1,jmax1-2,2
            do i1=1,imax1-2,2
              i2 = 1 + i1/2
              j2 = 1 + j1/2
              k2 = 1 + k1/2
              cutoff2(i2,j2,k2,n) = &
                  ( vol1(i1  ,j1  ,k1  ) * cutoff1(i1  ,j1  ,k1  ,n) + &
                    vol1(i1  ,j1+1,k1  ) * cutoff1(i1  ,j1+1,k1  ,n) + &
                    vol1(i1+1,j1  ,k1  ) * cutoff1(i1+1,j1  ,k1  ,n) + &
                    vol1(i1+1,j1+1,k1  ) * cutoff1(i1+1,j1+1,k1  ,n) + &
                    vol1(i1  ,j1  ,k1+1) * cutoff1(i1  ,j1  ,k1+1,n) + &
                    vol1(i1  ,j1+1,k1+1) * cutoff1(i1  ,j1+1,k1+1,n) + &
                    vol1(i1+1,j1  ,k1+1) * cutoff1(i1+1,j1  ,k1+1,n) + &
                    vol1(i1+1,j1+1,k1+1) * cutoff1(i1+1,j1+1,k1+1,n) ) / &
                  vol2(i2,j2,k2)
            end do
          end do
        end do
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine rest4(nl,prd1,prd2,vol1,vol2)
!-----------------------------------------------------------------------

      use cosa_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,iblock,imax1,jmax1,kmax1,imax2,jmax2,kmax2,iprd1, &
                iprd2,ivol1,ivol2
      real(kind=cosa_real) prd1(*),prd2(*),vol1(*),vol2(*)

      do iblock = 1,mynblocks
        imax1  = i_imax(iblock,nl)
        jmax1  = j_jmax(iblock,nl)
        kmax1  = k_kmax(iblock,nl)
        imax2  = i_imax(iblock,nl+1)
        jmax2  = j_jmax(iblock,nl+1)
        kmax2  = k_kmax(iblock,nl+1)
        iprd1  = 1 + off_p3 (iblock,nl  ) *        dim5
        iprd2  = 1 + off_p3 (iblock,nl+1) *        dim5
        ivol1  = 1 + off_p1 (iblock,nl  )
        ivol2  = 1 + off_p1 (iblock,nl+1)
        call rest4_b(imax1,jmax1,kmax1,imax2,jmax2,kmax2,nharms, &
                     prd1(iprd1),prd2(iprd2),vol1(ivol1),vol2(ivol2))
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine rest4_b(imax1,jmax1,kmax1,imax2,jmax2,kmax2,nharms, &
                         prd1,prd2,vol1,vol2)
!-----------------------------------------------------------------------

      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax1,jmax1,kmax1,imax2,jmax2,kmax2,nharms
      integer(kind=cosa_int) i1,j1,k1,i2,j2,k2,n
      real(kind=cosa_real) &
          prd1(-1:imax1+1,-1:jmax1+1,-1:kmax1+1,0:2*nharms), &
          prd2(-1:imax2+1,-1:jmax2+1,-1:kmax2+1,0:2*nharms), &
          vol1(0:imax1   ,0:jmax1   ,0:kmax1), &
          vol2(0:imax2   ,0:jmax2   ,0:kmax2)

!     prd1    : production term on FINER grid
!     prd2    : production term restricted to COARSER grid.

      do n = 0,2*nharms
        do k1=1,kmax1-2,2
          do j1=1,jmax1-2,2
            do i1=1,imax1-2,2
              i2 = 1 + i1/2
              j2 = 1 + j1/2
              k2 = 1 + k1/2
              prd2(i2,j2,k2,n) = &
                  ( vol1(i1  ,j1  ,k1  ) * prd1(i1  ,j1  ,k1  ,n) + &
                    vol1(i1  ,j1+1,k1  ) * prd1(i1  ,j1+1,k1  ,n) + &
                    vol1(i1+1,j1  ,k1  ) * prd1(i1+1,j1  ,k1  ,n) + &
                    vol1(i1+1,j1+1,k1  ) * prd1(i1+1,j1+1,k1  ,n) + &
                    vol1(i1  ,j1  ,k1+1) * prd1(i1  ,j1  ,k1+1,n) + &
                    vol1(i1  ,j1+1,k1+1) * prd1(i1  ,j1+1,k1+1,n) + &
                    vol1(i1+1,j1  ,k1+1) * prd1(i1+1,j1  ,k1+1,n) + &
                    vol1(i1+1,j1+1,k1+1) * prd1(i1+1,j1+1,k1+1,n) ) / &
                  vol2(i2,j2,k2)
            end do
          end do
        end do
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine prolong(nl)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl
      real (kind=cosa_real) q1, q2,        qold2, work, vol1, vol2, dist, si2, &
        sj2, sk2, mut
      integer(kind=cosa_int) bctopo2,cutopo2,percutopo2
      pointer &
        (pq1   ,   q1),(pq2   ,   q2),               (pqold2,qold2), &
        (pvol1 , vol1),(pvol2 , vol2),(pdist , dist),(psi2  ,  si2), &
        (psj2  ,  sj2),(psk2  ,  sk2),(pwork , work), &
        (pbctopo2 , bctopo2),(pcutopo2 , cutopo2), &
        (ppercutopo2,percutopo2)

! --- xxxx  :
! --  xxxx  :

      pq1         = p_q      (nl)
      pvol1       = p_vol    (nl)

      pq2         = p_q      (nl+1)
      pqold2      = p_qold   (nl+1)
      pwork       = p_dq     (nl+1)
      pbctopo2    = p_bctopo (nl+1)
      pcutopo2    = p_cutopo (nl+1)
      ppercutopo2 = p_percutopo(nl+1)
      psi2        = p_si     (nl+1)
      psj2        = p_sj     (nl+1)
      psk2        = p_sk     (nl+1)
      pvol2       = p_vol    (nl+1)
      pdist       = p_dist   (nl+1)

      call rest1(nl,q1,qold2,vol1,vol2)

      call cons2prim(nl+1,0,qold2,work)
      call bc(qold2,si2,sj2,sk2,dist,nl+1,bctopo2)
      call cutman_q(nl+1,cutopo2,qold2)
      if (persec) call pcutman_q(nl+1,percutopo2,qold2)
      call prim2cons(nl+1,qold2,work)
      if (bln_prol) then
        calmut = .false.
        call q_edges(nl+1,qold2,mut)
        call q_corns(nl+1,qold2,mut)
      else
        calmut = .false.
        call q_edges(nl+1,qold2,mut)
        call q_corns(nl+1,qold2,mut)
      end if

      call cons2prim(nl+1,0,q2,work)
      call bc(q2,si2,sj2,sk2,dist,nl+1,bctopo2)
      call cutman_q(nl+1,cutopo2,q2)
      if (persec) call pcutman_q(nl+1,percutopo2,q2)
      call prim2cons(nl+1,q2,work)
      if (bln_prol) then
        calmut = .false.
        call q_edges(nl+1,q2,mut)
        call q_corns(nl+1,q2,mut)
      else
        calmut = .false.
        call q_edges(nl+1,q2,mut)
        call q_corns(nl+1,q2,mut)
      end if

      call add_array(1,npde,0,2*nharms, 3,nl+1,q2,qold2, &
                     1,npde,0,2*nharms, 1,-1.d0)

      if (.not.bln_prol) then
        call prol_lin(nl,q1,q2,vol1)
      else
        call prol_bln(nl,q1,q2,vol1,lim_corr)
      end if

      if (lim_corr(1).or.lim_corr(3)) call fixq(nl,q1)

      return
      end

!-----------------------------------------------------------------------
      subroutine prol_lin(nl,q1,q2,vol1)
!-----------------------------------------------------------------------

      use cosa_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,iblock,imax1,jmax1,kmax1,imax2,jmax2,kmax2,iq1,iq2, &
                ivol1
      real (kind=cosa_real) q1(*),q2(*),vol1(*)

      do iblock = 1,mynblocks
        imax1  = i_imax(iblock,nl)
        jmax1  = j_jmax(iblock,nl)
        kmax1  = k_kmax(iblock,nl)
        imax2  = i_imax(iblock,nl+1)
        jmax2  = j_jmax(iblock,nl+1)
        kmax2  = k_kmax(iblock,nl+1)
        iq1    = 1 + off_p3 (iblock,nl  ) * npde * dim5
        iq2    = 1 + off_p3 (iblock,nl+1) * npde * dim5
        ivol1  = 1 + off_p1 (iblock,nl  )

        call prol_lin_b(imax1,jmax1,kmax1,imax2,jmax2,kmax2,npde, &
          nharms,q1(iq1),q2(iq2),vol1(ivol1))
      end do

      return
      end
!-----------------------------------------------------------------------
      subroutine prol_bln(nl,q1,q2,vol1,lim_corr)
!-----------------------------------------------------------------------

      use cosa_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,iblock,imax1,jmax1,kmax1,imax2,jmax2,kmax2,iq1,iq2, &
                ivol1
      real (kind=cosa_real) q1(*),q2(*),vol1(*)
      logical lim_corr(3)

      do iblock = 1,mynblocks
        imax1  = i_imax(iblock,nl)
        jmax1  = j_jmax(iblock,nl)
        kmax1  = k_kmax(iblock,nl)
        imax2  = i_imax(iblock,nl+1)
        jmax2  = j_jmax(iblock,nl+1)
        kmax2  = k_kmax(iblock,nl+1)
        iq1    = 1 + off_p3 (iblock,nl  ) * npde * dim5
        iq2    = 1 + off_p3 (iblock,nl+1) * npde * dim5
        ivol1  = 1 + off_p1 (iblock,nl  )

        call prol_bln_b(imax1,jmax1,kmax1,imax2,jmax2,kmax2,npde, &
          nharms,q1(iq1),q2(iq2),vol1(ivol1),lim_corr)

      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine prol_lin_b(imax1,jmax1,kmax1,imax2,jmax2,kmax2,npde, &
        nharms,q1,q2,vol1)
!-----------------------------------------------------------------------
!     prolongation based on linear interpolations
!-----------------------------------------------------------------------

      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax1,jmax1,kmax1,imax2,jmax2,kmax2,npde,nharms
      integer(kind=cosa_int) i1,j1,k1,i2,j2,k2,ipde,n
      real (kind=cosa_real) &
        q1  (-1:imax1+1,-1:jmax1+1,-1:kmax1+1,npde,0:2*nharms), &
        q2  (-1:imax2+1,-1:jmax2+1,-1:kmax2+1,npde,0:2*nharms), &
        vol1( 0:imax1  , 0:jmax1  , 0:kmax1                  )
      real (kind=cosa_real) cnec,cnwc,csec,cswc,dcnc,dcsc,dcec,dcwc,cnef_l, &
           cnwf_l,csef_l,cswf_l,cnef_u,cnwf_u,csef_u,cswf_u,fact

!-----------------------------------------------------------------------
!     q2  : coarse grid correction to be extrapolated on finegrid
! --  q1  : fine grid solution 

!     cnec   : Correction North-East (NE) on Coarse grid
!     cnwc   : Correction North-West (NW) on Coarse grid
!     csec   : Correction South-East (SE) on Coarse grid
!     cswc   : Correction South-West (SW) on Coarse grid

!     dcnc   : Delta Correction North on Coarse grid
!     dcsc   : Delta Correction South on Coarse grid
!     dcec   : Delta Correction East  on Coarse grid
!     dcwc   : Delta Correction West  on Coarse grid

!     cswf_l : Correction SW at Fine grid cell centre (based on k2 plane)
!     csef_l : Correction SE at Fine grid cell centre (based on k2 plane)
!     cnwf_l : Correction NW at Fine grid cell centre (based on k2 plane)
!     cnef_l : Correction NE at Fine grid cell centre (based on k2 plane)

!     cswf_u : Correction SW at Fine grid cell centre (based on k2+1 plane)
!     csef_u : Correction SE at Fine grid cell centre (based on k2+1 plane)
!     cnwf_u : Correction NW at Fine grid cell centre (based on k2+1 plane)
!     cnef_u : Correction NE at Fine grid cell centre (based on k2+1 plane)

!-----------------------------------------------------------------------

      do n = 0,2*nharms
        do ipde=1,npde
          do k1=0,kmax1-1,2
            do j1=0,jmax1-1,2
              do i1=0,imax1-1,2
                i2 = i1/2
                j2 = j1/2
                k2 = k1/2
!-------------- k2 plane
!cc             deltas of corrections
                cnec = q2(i2+1,j2+1,k2  ,ipde,n)
                cnwc = q2(i2  ,j2+1,k2  ,ipde,n)
                csec = q2(i2+1,j2  ,k2  ,ipde,n)
                cswc = q2(i2  ,j2  ,k2  ,ipde,n)
                dcnc = cnec - cnwc
                dcsc = csec - cswc
                dcec = cnec - csec
                dcwc = cnwc - cswc
!cc             weighted interpolation
                cswf_l = cswc + &
                  dcsc/2 * vol1(i1  ,j1  ,k1  ) / &
                  (vol1(i1  ,j1  ,k1  )+vol1(i1+1,j1  ,k1  )) + &
                  dcwc/2 * vol1(i1  ,j1  ,k1  ) / &
                  (vol1(i1  ,j1  ,k1  )+vol1(i1  ,j1+1,k1  ))
                csef_l = csec - &
                  dcsc/2 * vol1(i1+1,j1  ,k1  ) / &
                  (vol1(i1  ,j1  ,k1  )+vol1(i1+1,j1  ,k1  )) + &
                  dcec/2 * vol1(i1+1,j1  ,k1  ) / &
                  (vol1(i1+1,j1  ,k1  )+vol1(i1+1,j1+1,k1  ))
                cnwf_l = cnwc + &
                  dcnc/2 * vol1(i1  ,j1+1,k1  ) / &
                  (vol1(i1  ,j1+1,k1  )+vol1(i1+1,j1+1,k1  )) - &
                  dcwc/2 * vol1(i1  ,j1+1,k1  ) / &
                  (vol1(i1  ,j1  ,k1  )+vol1(i1  ,j1+1,k1  ))
                cnef_l = cnec - &
                  dcnc/2 * vol1(i1+1,j1+1,k1  ) / &
                  (vol1(i1+1,j1+1,k1  )+vol1(i1  ,j1+1,k1  )) - &
                  dcec/2 * vol1(i1+1,j1+1,k1  ) / &
                  (vol1(i1+1,j1+1,k1  )+vol1(i1+1,j1  ,k1  ))
!-------------- UNweighted interpolation (uniform mesh)
!               cswf_l = cswc + 0.25 * ( dcsc +dcwc)
!               csef_l = csec + 0.25 * ( dcec -dcsc)
!               cnef_l = cnec + 0.25 * (-dcnc -dcec)
!               cnwf_l = cnwc + 0.25 * ( dcnc -dcwc)
!-------------- k2+1 plane
!cc             deltas of corrections
                cnec = q2(i2+1,j2+1,k2+1,ipde,n)
                cnwc = q2(i2  ,j2+1,k2+1,ipde,n)
                csec = q2(i2+1,j2  ,k2+1,ipde,n)
                cswc = q2(i2  ,j2  ,k2+1,ipde,n)
                dcnc = cnec - cnwc
                dcsc = csec - cswc
                dcec = cnec - csec
                dcwc = cnwc - cswc
!cc             weighted interpolation
                cswf_u = cswc + &
                  dcsc/2 * vol1(i1  ,j1  ,k1+1) / &
                  (vol1(i1  ,j1  ,k1+1)+vol1(i1+1,j1  ,k1+1)) + &
                  dcwc/2 * vol1(i1  ,j1  ,k1+1) / &
                  (vol1(i1  ,j1  ,k1+1)+vol1(i1  ,j1+1,k1+1))
                csef_u = csec - &
                  dcsc/2 * vol1(i1+1,j1  ,k1+1) / &
                  (vol1(i1  ,j1  ,k1+1)+vol1(i1+1,j1  ,k1+1)) + &
                  dcec/2 * vol1(i1+1,j1  ,k1+1) / &
                  (vol1(i1+1,j1  ,k1+1)+vol1(i1+1,j1+1,k1+1))
                cnwf_u = cnwc + &
                  dcnc/2 * vol1(i1  ,j1+1,k1+1) / &
                  (vol1(i1  ,j1+1,k1+1)+vol1(i1+1,j1+1,k1+1)) - &
                  dcwc/2 * vol1(i1  ,j1+1,k1+1) / &
                  (vol1(i1  ,j1  ,k1+1)+vol1(i1  ,j1+1,k1+1))
                cnef_u = cnec - &
                  dcnc/2 * vol1(i1+1,j1+1,k1+1) / &
                  (vol1(i1+1,j1+1,k1+1)+vol1(i1  ,j1+1,k1+1)) - &
                  dcec/2 * vol1(i1+1,j1+1,k1+1) / &
                  (vol1(i1+1,j1+1,k1+1)+vol1(i1+1,j1  ,k1+1))
!-------------- UNweighted interpolation (uniform mesh)
!               cswf_u = cswc + 0.25 * ( dcsc +dcwc)
!               csef_u = csec + 0.25 * ( dcec -dcsc)
!               cnef_u = cnec + 0.25 * (-dcnc -dcec)
!               cnwf_u = cnwc + 0.25 * ( dcnc -dcwc)
!-------------- correct fine grid solution
!msc, 24 July 2013: fact can probably be computed only once, because if
!                   grid is sufficiently smooth, fact will not change
!                   much within each step of this loop.
                fact =  vol1(i1  ,j1  ,k1  ) / 2 / &
                       (vol1(i1  ,j1  ,k1  ) + vol1(i1  ,j1  ,k1+1))
                q1(i1  ,j1  ,k1  ,ipde,n) = q1(i1  ,j1  ,k1  ,ipde,n) + &
                  cswf_l + fact/2 * (cswf_u-cswf_l)
                q1(i1  ,j1  ,k1+1,ipde,n) = q1(i1  ,j1  ,k1+1,ipde,n) + &
                  cswf_u - (1-fact)/2 * (cswf_u-cswf_l)
                fact =  vol1(i1+1,j1  ,k1  ) / 2 / &
                       (vol1(i1+1,j1  ,k1  ) + vol1(i1+1,j1  ,k1+1))
                q1(i1+1,j1  ,k1  ,ipde,n) = q1(i1+1,j1  ,k1  ,ipde,n) + &
                  csef_l + fact/2 * (csef_u-csef_l)
                q1(i1+1,j1  ,k1+1,ipde,n) = q1(i1+1,j1  ,k1+1,ipde,n) + &
                  csef_u - (1-fact)/2 * (csef_u-csef_l)
                fact =  vol1(i1+1,j1+1,k1  ) / 2 / &
                       (vol1(i1+1,j1+1,k1  ) + vol1(i1+1,j1+1,k1+1))
                q1(i1+1,j1+1,k1  ,ipde,n) = q1(i1+1,j1+1,k1  ,ipde,n) + &
                  cnef_l + fact/2 * (cnef_u-cnef_l)
                q1(i1+1,j1+1,k1+1,ipde,n) = q1(i1+1,j1+1,k1+1,ipde,n) + &
                  cnef_u - (1-fact)/2 * (cnef_u-cnef_l)
                fact =  vol1(i1  ,j1+1,k1  ) / 2 / &
                       (vol1(i1  ,j1+1,k1  ) + vol1(i1  ,j1+1,k1+1))
                q1(i1  ,j1+1,k1  ,ipde,n) = q1(i1  ,j1+1,k1  ,ipde,n) + &
                  cnwf_l + fact/2 * (cnwf_u-cnwf_l)
                q1(i1  ,j1+1,k1+1,ipde,n) = q1(i1  ,j1+1,k1+1,ipde,n) + &
                  cnwf_u - (1-fact)/2 * (cnwf_u-cnwf_l)
              end do
            end do
          end do
        end do
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine prol_bln_b(imax1,jmax1,kmax1,imax2,jmax2,kmax2,npde, &
        nharms,q1,q2,vol1,lim_corr)
!-----------------------------------------------------------------------
!     prolongation based on bilinear interpolations
!-----------------------------------------------------------------------

      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax1,jmax1,kmax1,imax2,jmax2,kmax2,npde,nharms
      integer(kind=cosa_int) i1,j1,k1,i2,j2,k2,ipde,n
      real (kind=cosa_real) &
        q1  (-1:imax1+1,-1:jmax1+1,-1:kmax1+1,npde,0:2*nharms), &
        q2  (-1:imax2+1,-1:jmax2+1,-1:kmax2+1,npde,0:2*nharms), &
        vol1( 0:imax1  , 0:jmax1  , 0:kmax1                  )
      real (kind=cosa_real) cnec_l,cnwc_l,csec_l,cswc_l,cnec_u,cnwc_u,csec_u, &
           cswc_u,cnef_l,cnwf_l,csef_l,cswf_l,cnef_u,cnwf_u,csef_u, &
           cswf_u,qsw_l,qse_l,qne_l,qnw_l,qsw_u,qse_u,qne_u,qnw_u,beta, &
           epsl
      logical lim_corr(3)

!-----------------------------------------------------------------------
!     q2  : coarse grid correction to be extrapolated on finegrid
! --  q1  : fine grid solution 

!     cswc_l : Correction SW at Fine grid cell centre on k2 plane
!     csec_l : Correction SE at Fine grid cell centre on k2 plane
!     cnwc_l : Correction NW at Fine grid cell centre on k2 plane
!     cnec_l : Correction NE at Fine grid cell centre on k2 plane

!     cswc_u : Correction SW at Fine grid cell centre on k2+1 plane
!     csec_u : Correction SE at Fine grid cell centre on k2+1 plane
!     cnwc_u : Correction NW at Fine grid cell centre on k2+1 plane
!     cnec_u : Correction NE at Fine grid cell centre on k2+1 plane

!     cswf_l : Correction SW at Fine grid cell centre (based on k2 plane)
!     csef_l : Correction SE at Fine grid cell centre (based on k2 plane)
!     cnwf_l : Correction NW at Fine grid cell centre (based on k2 plane)
!     cnef_l : Correction NE at Fine grid cell centre (based on k2 plane)

!     cswf_u : Correction SW at Fine grid cell centre (based on k2+1 plane)
!     csef_u : Correction SE at Fine grid cell centre (based on k2+1 plane)
!     cnwf_u : Correction NW at Fine grid cell centre (based on k2+1 plane)
!     cnef_u : Correction NE at Fine grid cell centre (based on k2+1 plane)

!-----------------------------------------------------------------------

      beta  = 100.d0
      epsl  = 1.d-12

      do n = 0,2*nharms
        do ipde=1,5
          do k1=0,kmax1-1,2
            do j1=0,jmax1-1,2
              do i1=0,imax1-1,2
                i2 = i1/2
                j2 = j1/2
                k2 = k1/2
!-------------- set coarse grid corrections
                cnec_l = q2(i2+1,j2+1,k2  ,ipde,n)
                cnwc_l = q2(i2  ,j2+1,k2  ,ipde,n)
                csec_l = q2(i2+1,j2  ,k2  ,ipde,n)
                cswc_l = q2(i2  ,j2  ,k2  ,ipde,n)
                cnec_u = q2(i2+1,j2+1,k2+1,ipde,n)
                cnwc_u = q2(i2  ,j2+1,k2+1,ipde,n)
                csec_u = q2(i2+1,j2  ,k2+1,ipde,n)
                cswc_u = q2(i2  ,j2  ,k2+1,ipde,n)
!-------------- fine grid corrections
                cswf_l = (27*cswc_l + 9*csec_l + 9*cnwc_l + 3*cnec_l + &
                           9*cswc_u + 3*csec_u + 3*cnwc_u +   cnec_u)/64
                csef_l = (27*csec_l + 9*cswc_l + 9*cnec_l + 3*cnwc_l + &
                           9*csec_u + 3*cswc_u + 3*cnec_u +   cnwc_u)/64
                cnef_l = (27*cnec_l + 9*cnwc_l + 9*csec_l + 3*cswc_l + &
                           9*cnec_u + 3*cnwc_u + 3*csec_u +   cswc_u)/64
                cnwf_l = (27*cnwc_l + 9*cnec_l + 9*cswc_l + 3*csec_l + &
                           9*cnwc_u + 3*cnec_u + 3*cswc_u +   csec_u)/64
                cswf_u = (27*cswc_u + 9*csec_u + 9*cnwc_u + 3*cnec_u + &
                           9*cswc_l + 3*csec_l + 3*cnwc_l +   cnec_l)/64
                csef_u = (27*csec_u + 9*cswc_u + 9*cnec_u + 3*cnwc_u + &
                           9*csec_l + 3*cswc_l + 3*cnec_l +   cnwc_l)/64
                cnef_u = (27*cnec_u + 9*cnwc_u + 9*csec_u + 3*cswc_u + &
                           9*cnec_l + 3*cnwc_l + 3*csec_l +   cswc_l)/64
                cnwf_u = (27*cnwc_u + 9*cnec_u + 9*cswc_u + 3*csec_u + &
                           9*cnwc_l + 3*cnec_l + 3*cswc_l +   csec_l)/64
!-------------- correct fine grid solution
                q1(i1  ,j1  ,k1  ,ipde,n) = &
                  q1(i1  ,j1  ,k1  ,ipde,n) + cswf_l
                q1(i1+1,j1  ,k1  ,ipde,n) = &
                  q1(i1+1,j1  ,k1  ,ipde,n) + csef_l
                q1(i1+1,j1+1,k1  ,ipde,n) = &
                  q1(i1+1,j1+1,k1  ,ipde,n) + cnef_l
                q1(i1  ,j1+1,k1  ,ipde,n) = &
                  q1(i1  ,j1+1,k1  ,ipde,n) + cnwf_l
                q1(i1  ,j1  ,k1+1,ipde,n) = &
                  q1(i1  ,j1  ,k1+1,ipde,n) + cswf_u
                q1(i1+1,j1  ,k1+1,ipde,n) = &
                  q1(i1+1,j1  ,k1+1,ipde,n) + csef_u
                q1(i1+1,j1+1,k1+1,ipde,n) = &
                  q1(i1+1,j1+1,k1+1,ipde,n) + cnef_u
                q1(i1  ,j1+1,k1+1,ipde,n) = &
                  q1(i1  ,j1+1,k1+1,ipde,n) + cnwf_u
              end do
            end do
          end do
        end do
      end do

      if (lim_corr(1)) then

        do n = 0,2*nharms
          do ipde=6,npde
            do k1=0,kmax1-1,2
              do j1=0,jmax1-1,2
                do i1=0,imax1-1,2
                  i2 = i1/2
                  j2 = j1/2
                  k2 = k1/2
!---------------- set coarse grid corrections
                  cnec_l = q2(i2+1,j2+1,k2  ,ipde,n)
                  cnwc_l = q2(i2  ,j2+1,k2  ,ipde,n)
                  csec_l = q2(i2+1,j2  ,k2  ,ipde,n)
                  cswc_l = q2(i2  ,j2  ,k2  ,ipde,n)
                  cnec_u = q2(i2+1,j2+1,k2+1,ipde,n)
                  cnwc_u = q2(i2  ,j2+1,k2+1,ipde,n)
                  csec_u = q2(i2+1,j2  ,k2+1,ipde,n)
                  cswc_u = q2(i2  ,j2  ,k2+1,ipde,n)
!---------------- fine grid corrections
                  cswf_l = (27*cswc_l + 9*csec_l + 9*cnwc_l + 3*cnec_l + &
                             9*cswc_u + 3*csec_u + 3*cnwc_u +   cnec_u)/ &
                           64
                  csef_l = (27*csec_l + 9*cswc_l + 9*cnec_l + 3*cnwc_l + &
                             9*csec_u + 3*cswc_u + 3*cnec_u +   cnwc_u)/ &
                           64
                  cnef_l = (27*cnec_l + 9*cnwc_l + 9*csec_l + 3*cswc_l + &
                             9*cnec_u + 3*cnwc_u + 3*csec_u +   cswc_u)/ &
                           64
                  cnwf_l = (27*cnwc_l + 9*cnec_l + 9*cswc_l + 3*csec_l + &
                             9*cnwc_u + 3*cnec_u + 3*cswc_u +   csec_u)/ &
                           64
                  cswf_u = (27*cswc_u + 9*csec_u + 9*cnwc_u + 3*cnec_u + &
                             9*cswc_l + 3*csec_l + 3*cnwc_l +   cnec_l)/ &
                           64
                  csef_u = (27*csec_u + 9*cswc_u + 9*cnec_u + 3*cnwc_u + &
                             9*csec_l + 3*cswc_l + 3*cnec_l +   cnwc_l)/ &
                           64
                  cnef_u = (27*cnec_u + 9*cnwc_u + 9*csec_u + 3*cswc_u + &
                             9*cnec_l + 3*cnwc_l + 3*csec_l +   cswc_l)/ &
                           64
                  cnwf_u = (27*cnwc_u + 9*cnec_u + 9*cswc_u + 3*csec_u + &
                             9*cnwc_l + 3*cnec_l + 3*cswc_l +   csec_l)/ &
                           64
!---------------- correct fine grid solution
                  q1(i1  ,j1  ,k1  ,ipde,n) = &
                    q1(i1  ,j1  ,k1  ,ipde,n) + cswf_l
                  q1(i1+1,j1  ,k1  ,ipde,n) = &
                    q1(i1+1,j1  ,k1  ,ipde,n) + csef_l
                  q1(i1+1,j1+1,k1  ,ipde,n) = &
                    q1(i1+1,j1+1,k1  ,ipde,n) + cnef_l
                  q1(i1  ,j1+1,k1  ,ipde,n) = &
                    q1(i1  ,j1+1,k1  ,ipde,n) + cnwf_l
                  q1(i1  ,j1  ,k1+1,ipde,n) = &
                    q1(i1  ,j1  ,k1+1,ipde,n) + cswf_u
                  q1(i1+1,j1  ,k1+1,ipde,n) = &
                    q1(i1+1,j1  ,k1+1,ipde,n) + csef_u
                  q1(i1+1,j1+1,k1+1,ipde,n) = &
                    q1(i1+1,j1+1,k1+1,ipde,n) + cnef_u
                  q1(i1  ,j1+1,k1+1,ipde,n) = &
                    q1(i1  ,j1+1,k1+1,ipde,n) + cnwf_u
                end do
              end do
            end do
          end do
        end do

      else if (lim_corr(2)) then

        do n = 0,2*nharms
          do ipde=6,npde
            do k1=0,kmax1-1,2
              do j1=0,jmax1-1,2
                do i1=0,imax1-1,2
                  i2 = i1/2
                  j2 = j1/2
                  k2 = k1/2
!---------------- set coarse grid corrections
                  cnec_l = q2(i2+1,j2+1,k2  ,ipde,n)
                  cnwc_l = q2(i2  ,j2+1,k2  ,ipde,n)
                  csec_l = q2(i2+1,j2  ,k2  ,ipde,n)
                  cswc_l = q2(i2  ,j2  ,k2  ,ipde,n)
                  cnec_u = q2(i2+1,j2+1,k2+1,ipde,n)
                  cnwc_u = q2(i2  ,j2+1,k2+1,ipde,n)
                  csec_u = q2(i2+1,j2  ,k2+1,ipde,n)
                  cswc_u = q2(i2  ,j2  ,k2+1,ipde,n)
!---------------- fine grid corrections
                  cswf_l = (27*cswc_l + 9*csec_l + 9*cnwc_l + 3*cnec_l + &
                             9*cswc_u + 3*csec_u + 3*cnwc_u +   cnec_u)/ &
                           64
                  csef_l = (27*csec_l + 9*cswc_l + 9*cnec_l + 3*cnwc_l + &
                             9*csec_u + 3*cswc_u + 3*cnec_u +   cnwc_u)/ &
                           64
                  cnef_l = (27*cnec_l + 9*cnwc_l + 9*csec_l + 3*cswc_l + &
                             9*cnec_u + 3*cnwc_u + 3*csec_u +   cswc_u)/ &
                           64
                  cnwf_l = (27*cnwc_l + 9*cnec_l + 9*cswc_l + 3*csec_l + &
                             9*cnwc_u + 3*cnec_u + 3*cswc_u +   csec_u)/ &
                           64
                  cswf_u = (27*cswc_u + 9*csec_u + 9*cnwc_u + 3*cnec_u + &
                             9*cswc_l + 3*csec_l + 3*cnwc_l +   cnec_l)/ &
                           64
                  csef_u = (27*csec_u + 9*cswc_u + 9*cnec_u + 3*cnwc_u + &
                             9*csec_l + 3*cswc_l + 3*cnec_l +   cnwc_l)/ &
                           64
                  cnef_u = (27*cnec_u + 9*cnwc_u + 9*csec_u + 3*cswc_u + &
                             9*cnec_l + 3*cnwc_l + 3*csec_l +   cswc_l)/ &
                           64
                  cnwf_u = (27*cnwc_u + 9*cnec_u + 9*cswc_u + 3*csec_u + &
                             9*cnwc_l + 3*cnec_l + 3*cswc_l +   csec_l)/ &
                           64
!---------------- limit cswf_l fine grid corrections
                  qsw_l = q1(i1  ,j1  ,k1  ,ipde,n) 
                  if ((qsw_l+cswf_l).le.0) then
                    if ((qsw_l+cswf_l/2).gt.0) then
                      cswf_l = cswf_l/2
                    else if ((qsw_l+cswf_l/4).gt.0) then
                      cswf_l = cswf_l/4
                    else
                      cswf_l = 0.d0   
                    end if
                  end if
!---------------- limit csef_l fine grid corrections
                  qse_l = q1(i1+1,j1  ,k1  ,ipde,n)
                  if ((qse_l+csef_l).le.0) then
                    if ((qse_l+csef_l/2).gt.0) then
                      csef_l = csef_l/2
                    else if ((qse_l+csef_l/4).gt.0) then
                      csef_l = csef_l/4
                    else
                      csef_l = 0.d0   
                    end if
                  end if
!---------------- limit cnef_l fine grid corrections
                  qne_l = q1(i1+1,j1+1,k1  ,ipde,n)
                  if ((qne_l+cnef_l).le.0) then
                    if ((qne_l+cnef_l/2).gt.0) then
                      cnef_l = cnef_l/2
                    else if ((qse_l+csef_l/4).gt.0) then
                      cnef_l = cnef_l/4
                    else
                      cnef_l = 0.d0   
                    end if
                  end if
!---------------- limit cnwf_l fine grid corrections
                  qnw_l = q1(i1  ,j1+1,k1  ,ipde,n)
                  if ((qnw_l+cnwf_l).le.0) then
                    if ((qnw_l+cnwf_l/2).gt.0) then
                      cnwf_l = cnwf_l/2
                    else if ((qnw_l+cnwf_l/4).gt.0) then
                      cnwf_l = cnwf_l/4
                    else
                      cnwf_l = 0.d0   
                    end if
                  end if
!---------------- limit cswf_u fine grid corrections
                  qsw_u = q1(i1  ,j1  ,k1+1,ipde,n) 
                  if ((qsw_u+cswf_u).le.0) then
                    if ((qsw_u+cswf_u/2).gt.0) then
                      cswf_u = cswf_u/2
                    else if ((qsw_u+cswf_u/4).gt.0) then
                      cswf_u = cswf_u/4
                    else
                      cswf_u = 0.d0   
                    end if
                  end if
!---------------- limit csef_u fine grid corrections
                  qse_u = q1(i1+1,j1  ,k1+1,ipde,n)
                  if ((qse_u+csef_u).le.0) then
                    if ((qse_u+csef_u/2).gt.0) then
                      csef_u = csef_u/2
                    else if ((qse_u+csef_u/4).gt.0) then
                      csef_u = csef_u/4
                    else
                      csef_u = 0.d0   
                    end if
                  end if
!---------------- limit cnef_u fine grid corrections
                  qne_u = q1(i1+1,j1+1,k1+1,ipde,n)
                  if ((qne_u+cnef_u).le.0) then
                    if ((qne_u+cnef_u/2).gt.0) then
                      cnef_u = cnef_u/2
                    else if ((qse_u+csef_u/4).gt.0) then
                      cnef_u = cnef_u/4
                    else
                      cnef_u = 0.d0   
                    end if
                  end if
!---------------- limit cnwf_u fine grid corrections
                  qnw_u = q1(i1  ,j1+1,k1+1,ipde,n)
                  if ((qnw_u+cnwf_u).le.0) then
                    if ((qnw_u+cnwf_u/2).gt.0) then
                      cnwf_u = cnwf_u/2
                    else if ((qnw_u+cnwf_u/4).gt.0) then
                      cnwf_u = cnwf_u/4
                    else
                      cnwf_u = 0.d0   
                    end if
                  end if
!---------------- correct fine grid solution
                  q1(i1  ,j1  ,k1  ,ipde,n) = &
                    q1(i1  ,j1  ,k1  ,ipde,n) + cswf_l
                  q1(i1+1,j1  ,k1  ,ipde,n) = &
                    q1(i1+1,j1  ,k1  ,ipde,n) + csef_l
                  q1(i1+1,j1+1,k1  ,ipde,n) = &
                    q1(i1+1,j1+1,k1  ,ipde,n) + cnef_l
                  q1(i1  ,j1+1,k1  ,ipde,n) = &
                    q1(i1  ,j1+1,k1  ,ipde,n) + cnwf_l
                  q1(i1  ,j1  ,k1+1,ipde,n) = &
                    q1(i1  ,j1  ,k1+1,ipde,n) + cswf_u
                  q1(i1+1,j1  ,k1+1,ipde,n) = &
                    q1(i1+1,j1  ,k1+1,ipde,n) + csef_u
                  q1(i1+1,j1+1,k1+1,ipde,n) = &
                    q1(i1+1,j1+1,k1+1,ipde,n) + cnef_u
                  q1(i1  ,j1+1,k1+1,ipde,n) = &
                    q1(i1  ,j1+1,k1+1,ipde,n) + cnwf_u
                end do
              end do
            end do
          end do
        end do

      else if (lim_corr(3)) then

        do n = 0,2*nharms
          do ipde=6,npde
            do k1=0,kmax1-1,2
              do j1=0,jmax1-1,2
                do i1=0,imax1-1,2
                  i2 = i1/2
                  j2 = j1/2
                  k2 = k1/2
!---------------- set coarse grid corrections
                  cnec_l = q2(i2+1,j2+1,k2  ,ipde,n)
                  cnwc_l = q2(i2  ,j2+1,k2  ,ipde,n)
                  csec_l = q2(i2+1,j2  ,k2  ,ipde,n)
                  cswc_l = q2(i2  ,j2  ,k2  ,ipde,n)
                  cnec_u = q2(i2+1,j2+1,k2+1,ipde,n)
                  cnwc_u = q2(i2  ,j2+1,k2+1,ipde,n)
                  csec_u = q2(i2+1,j2  ,k2+1,ipde,n)
                  cswc_u = q2(i2  ,j2  ,k2+1,ipde,n)
!---------------- fine grid corrections
                  cswf_l = (27*cswc_l + 9*csec_l + 9*cnwc_l + 3*cnec_l + &
                             9*cswc_u + 3*csec_u + 3*cnwc_u +   cnec_u)/ &
                           64
                  csef_l = (27*csec_l + 9*cswc_l + 9*cnec_l + 3*cnwc_l + &
                             9*csec_u + 3*cswc_u + 3*cnec_u +   cnwc_u)/ &
                           64
                  cnef_l = (27*cnec_l + 9*cnwc_l + 9*csec_l + 3*cswc_l + &
                             9*cnec_u + 3*cnwc_u + 3*csec_u +   cswc_u)/ &
                           64
                  cnwf_l = (27*cnwc_l + 9*cnec_l + 9*cswc_l + 3*csec_l + &
                             9*cnwc_u + 3*cnec_u + 3*cswc_u +   csec_u)/ &
                           64
                  cswf_u = (27*cswc_u + 9*csec_u + 9*cnwc_u + 3*cnec_u + &
                             9*cswc_l + 3*csec_l + 3*cnwc_l +   cnec_l)/ &
                           64
                  csef_u = (27*csec_u + 9*cswc_u + 9*cnec_u + 3*cnwc_u + &
                             9*csec_l + 3*cswc_l + 3*cnec_l +   cnwc_l)/ &
                           64
                  cnef_u = (27*cnec_u + 9*cnwc_u + 9*csec_u + 3*cswc_u + &
                             9*cnec_l + 3*cnwc_l + 3*csec_l +   cswc_l)/ &
                           64
                  cnwf_u = (27*cnwc_u + 9*cnec_u + 9*cswc_u + 3*csec_u + &
                             9*cnwc_l + 3*cnec_l + 3*cswc_l +   csec_l)/ &
                           64

!---------------- limit correction and correct fine grid solution
                  q1(i1  ,j1  ,k1  ,ipde,n) = q1(i1  ,j1  ,k1  ,ipde,n)+ &
            cswf_l / &
            (1 - beta*min(cswf_l/(q1(i1  ,j1  ,k1  ,ipde,n)+epsl),0.d0))
                  q1(i1+1,j1  ,k1  ,ipde,n) = q1(i1+1,j1  ,k1  ,ipde,n)+ &
            csef_l / &
            (1 - beta*min(csef_l/(q1(i1+1,j1  ,k1  ,ipde,n)+epsl),0.d0))
                  q1(i1+1,j1+1,k1  ,ipde,n) = q1(i1+1,j1+1,k1  ,ipde,n)+ &
            cnef_l / &
            (1 - beta*min(cnef_l/(q1(i1+1,j1+1,k1  ,ipde,n)+epsl),0.d0))
                  q1(i1  ,j1+1,k1  ,ipde,n) = q1(i1  ,j1+1,k1  ,ipde,n)+ &
            cnwf_l / &
            (1 - beta*min(cnwf_l/(q1(i1  ,j1+1,k1  ,ipde,n)+epsl),0.d0))
                  q1(i1  ,j1  ,k1+1,ipde,n) = q1(i1  ,j1  ,k1+1,ipde,n)+ &
            cswf_u / &
            (1 - beta*min(cswf_u/(q1(i1  ,j1  ,k1+1,ipde,n)+epsl),0.d0))
                  q1(i1+1,j1  ,k1+1,ipde,n) = q1(i1+1,j1  ,k1+1,ipde,n)+ &
            csef_u / &
            (1 - beta*min(csef_u/(q1(i1+1,j1  ,k1+1,ipde,n)+epsl),0.d0))
                  q1(i1+1,j1+1,k1+1,ipde,n) = q1(i1+1,j1+1,k1+1,ipde,n)+ &
            cnef_u / &
            (1 - beta*min(cnef_u/(q1(i1+1,j1+1,k1+1,ipde,n)+epsl),0.d0))
                  q1(i1  ,j1+1,k1+1,ipde,n) = q1(i1  ,j1+1,k1+1,ipde,n)+ &
            cnwf_u / &
            (1 - beta*min(cnwf_u/(q1(i1  ,j1+1,k1+1,ipde,n)+epsl),0.d0))
                end do
              end do
            end do
          end do
        end do

      end if

      return
      end

!-----------------------------------------------------------------------
      subroutine rkuns()
!-----------------------------------------------------------------------

      use common_variables
      use cosa_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) irk,idir,nl

      real(kind=cosa_real) q, q1, qold, qp, qm, dq, dqp, dqm, pr, flux, res, &
          rhs, work1, mut, prd, divvel, delpls, ttrm, &
          prec, preci, ipkpi, dtvol, dtvolt, &
          si, sj, sk, x, y, z, x0, y0, z0, dx, dy, dz, xdot, ydot, zdot, &
          workx, worky, workz, vol, &
          dist, fidist, fjdist, fkdist, rad, &
          xideri, xiderj, xiderk, etaderi, etaderj, etaderk, &
          zetaderi, zetaderj, zetaderk

      integer(kind=cosa_int) bctopo,cutopo,percutopo

      logical amcontrol

      real (kind=cosa_real) gtime

      pointer &
         (pq   ,   q), (pq1  ,  q1), (pqold,qold), (pdq  ,  dq), &
         (pdqp , dqp), (pdqm , dqm), (pqp  ,  qp), (pqm  ,  qm), &
         (pflux,flux), (pres , res), (ppr  ,  pr), &
         (psi  ,  si), (psj  ,  sj), (psk  ,  sk), &
         (px   ,   x), (py   ,   y), (pz   ,   z), &
         (pxdot,xdot), (pydot,ydot), (pzdot,zdot), &
         (px0  ,  x0), (py0  ,  y0), (pz0  ,  z0), &
         (pdx  ,  dx), (pdy  ,  dy), (pdz  ,  dz), &
         (pvol , vol), (prad , rad), (pdist,dist), &
         (pdtvol   ,   dtvol),(pdtvolt  ,  dtvolt), &
         (pmut     ,     mut),(pttrm    ,    ttrm),(pdelpls  ,  delpls), &
         (pdivvel  ,  divvel),(pprd     ,     prd), &
         (pfidist  ,  fidist),(pfjdist  ,  fjdist),(pfkdist  ,  fkdist), &
         (pworkx   ,   workx),(pworky   ,   worky),(pworkz   ,   workz), &
         (pxideri  ,  xideri),(pxiderj  ,  xiderj),(pxiderk  ,  xiderk), &
         (petaderi , etaderi),(petaderj , etaderj),(petaderk , etaderk), &
         (pzetaderi,zetaderi),(pzetaderj,zetaderj),(pzetaderk,zetaderk), &
         (pwork1   ,   work1), &
         (pprec    ,    prec),(ppreci   ,   preci),(pipkpi   ,   ipkpi), &
         (pbctopo  ,  bctopo),(pcutopo  ,  cutopo), &
         (ppercutopo , percutopo)

!---- Runge-Kutta integration scheme
      nl        = 1
      px        = p_x(nl)
      py        = p_y(nl)
      pz        = p_z(nl)
      pdx       = p_dx(nl)
      pdy       = p_dy(nl)
      pdz       = p_dz(nl)
      pq        = p_q(nl)
      pqp       = p_qp(nl)
      pqm       = p_qm(nl)
      pflux     = p_flux(nl)
      pmut      = p_mut(nl)
      pres      = p_res(nl)
      pdtvol    = p_dtvol(nl)
      pdtvolt   = p_dtvolt(nl)
      ppr       = p_pr(nl)
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
      pxdot     = p_xdot(nl)
      pydot     = p_ydot(nl)
      pzdot     = p_zdot(nl)
      px0       = p_x0(nl)
      py0       = p_y0(nl)
      pz0       = p_z0(nl)
      pvol      = p_vol(nl)
      pdist     = p_dist(nl)
      pfidist   = p_fidist(nl)
      pfjdist   = p_fjdist(nl)
      pfkdist   = p_fkdist(nl)
      pdelpls   = p_delpls(nl)
      pdivvel   = p_divvel(nl)
      pprd      = p_prd(nl)
      pttrm     = p_ttrm(nl)
      prad      = p_rad(nl)
      pqold     = p_qold(nl)
      pdq       = p_dq(nl)
      pdqp      = p_dqp(nl)
      pdqm      = p_dqm(nl)
      pq1       = p_q1(nl)
      pprec     = p_prec(nl)
      pworkx    = p_xc(nl)
      pworky    = p_yc(nl)
      pworkz    = p_zc(nl)
      ppreci    = p_preci(nl)
      pwork1    = p_work1(nl)
      pipkpi    = p_ipkpi(nl)
      pbctopo   = p_bctopo(nl)
      pcutopo   = p_cutopo(nl)
      ppercutopo= p_percutopo(nl)

      gtime = simtime

      do itime = 1,ntime

        call amcontroller(amcontrol)

!------ move initial solution to qold
        call copy_array(1,npde,0,0, 3,nl,qold,q, &
                        1,npde,0,0,-1)

        do irk = 1,nstage

          simtime = gtime + dt*br4(irk)

          if (moving) then
            call movegrid(nl,x,y,z,xdot,ydot,zdot,dx,dy,dz,x0,y0,z0)
            call vols(nl,x,y,z,xdot,ydot,zdot,vol,rad)
            call imetrix(nl,x,y,z,xdot,ydot,zdot,si,sj,sk)
            if (viscous) then
              call volex(nl,vol)
              call cutman_v(nl,cutopo,vol)
              if (persec) call pcutman_v(nl,percutopo,vol)
            end if
            if (viscous) then
              call coordex(nl,x,y,z)
              call cutman_x(nl,cutopo,x,y,z)
              if (persec) call pcutman_x(nl,percutopo,x,y,z)
!             call vmetrix_old(nl,x,y,z,si,sj,sk,xideri,xiderj,xiderk, &
!                          etaderi,etaderj,etaderk,zetaderi,zetaderj, &
!                          zetaderk,vol)
              call vmetrix(nl,x,y,z,xideri,xiderj,xiderk,etaderi, &
                etaderj,etaderk,zetaderi,zetaderj,zetaderk,vol)
            end if
            pworkx   = p_xc(nl)
            pworky   = p_yc(nl)
            pworkz   = p_zc(nl)
            call vert2cell_g(nl,x   ,y   ,z   ,workx,worky,workz)
            call vert2cell_g(nl,xdot,ydot,zdot,workx,worky,workz)
          end if
          if ((irk.eq.1).and.(itime.eq.1)) then
            call deltat_urk(nl,vol,dtvol)
          end if

          call cons2prim(nl,0,q,dq)

          call bc(q,si,sj,sk,dist,nl,bctopo)
          call cutman_q(nl,cutopo,q)
          if (persec) call pcutman_q(nl,percutopo,q)
          call prtest(nl,q,pr,itime)

          call zero(1,npde,0,0, 3,nl,res, &
                    1,npde,0,0,-1)

          idir = 1
          call muscl(idir,nl,bctopo,q,qp,qm,dq,dqp,dqm)
          call roflux(idir,nl,qp,qm,dqp,dqm,flux,si)
          call resid(idir,nl,flux,res)
          if (viscous) then
            if (kom.or.kom_bsl.or.kom_sst) then
              calmut = .true.
            else
              calmut = .false.
            end if
            call q_edges(nl,q,mut)
            call q_face(idir,nl,bctopo,q,dq,mut,flux)
            call q_der(idir,nl,bctopo,q,qp,qm,dqm)
            if (kom.or.kom_bsl.or.kom_sst) then
              call rtst(idir,nl,dq,qp,qm,dqm,xideri,etaderi,zetaderi, &
                        dqp,flux)
            end if
            call vflux(idir,nl,dq,qp,qm,dqm,dqp,flux,fidist,xideri, &
                       etaderi,zetaderi,si)
            call resid(idir,nl,flux,res)
          end if

          idir = 2
          call muscl(idir,nl,bctopo,q,qp,qm,dq,dqp,dqm)
          call roflux(idir,nl,qp,qm,dqp,dqm,flux,sj)
          call resid(idir,nl,flux,res)
          if (viscous) then
            if (kom.or.kom_bsl.or.kom_sst) then
              calmut = .true.
            else
              calmut = .false.
            end if
            call q_face(idir,nl,bctopo,q,dq,mut,flux)
            call q_der(idir,nl,bctopo,q,qp,qm,dqm)
            if (kom.or.kom_bsl.or.kom_sst) then
              call rtst(idir,nl,dq,qp,qm,dqm,xiderj,etaderj,zetaderj, &
                        dqp,flux)
            end if
            call vflux(idir,nl,dq,qp,qm,dqm,dqp,flux,fjdist,xiderj, &
                       etaderj,zetaderj,sj)
            call resid(idir,nl,flux,res)
          end if

          idir = 3
          call muscl(idir,nl,bctopo,q,qp,qm,dq,dqp,dqm)
          call roflux(idir,nl,qp,qm,dqp,dqm,flux,sk)
          call resid(idir,nl,flux,res)
          if (viscous) then
            if (kom.or.kom_bsl.or.kom_sst) then
              calmut = .true.
            else
              calmut = .false.
            end if
            call q_face(idir,nl,bctopo,q,dq,mut,flux)
            call q_der(idir,nl,bctopo,q,qp,qm,dqm)
            if (kom.or.kom_bsl.or.kom_sst) then
              call rtst(idir,nl,dq,qp,qm,dqm,xiderk,etaderk,zetaderk, &
                        dqp,flux)
            end if
            call vflux(idir,nl,dq,qp,qm,dqm,dqp,flux,fkdist,xiderk, &
                       etaderk,zetaderk,sk)
            call resid(idir,nl,flux,res)
          end if

          if (kom.or.kom_bsl.or.kom_sst) then
            calmut = .false.
            call q_face(1,nl,bctopo,q,dq,mut,flux)
            call q_der_c(1,nl,dq,qp)
            call q_face(2,nl,bctopo,q,dq,mut,flux)
            call q_der_c(2,nl,dq,qm)
            call q_face(3,nl,bctopo,q,dq,mut,flux)
            call q_der_c(3,nl,dq,dqm)
              call turprod(nl,q,qp,qm,dqm,prd,divvel,delpls,xideri, &
                etaderi,zetaderi,xiderj,etaderj,zetaderj,xiderk,etaderk, &
                zetaderk)
            if (kom) then
              call srckw(nl,q,mut,prd,divvel,res,vol)
            else if (kom_bsl.or.kom_sst) then
              call update_ttrm(nl,q,qp,qm,dqm,xideri,etaderi,zetaderi, &
                xiderj,etaderj,zetaderj,xiderk,etaderk,zetaderk,ttrm, &
                dist)
              call srckw_sst(nl,q,prd,divvel,res,mut,ttrm,vol)
            end if
          end if

          call prim2cons(nl,q,dq)

          if ((itime/iprint*iprint.eq.itime).or.(itime.eq.ntime)) then
            call resg_rms(nl,q,res,rhs,vol)
          end if
          call update(nl,res,q,qold,vol,dtvol,dtvolt,prec,ipkpi,work1, &
                      dq,dqp,delpls,prd,ttrm,alphas(irk))

!------ increment RK stage
        end do    

        if (calfor) then
          simtime = gtime + dt
          call force(nl,q,x,y,z,xdot,ydot,zdot,si,sj,sk,xideri,etaderj, &
                     zetaderk,bctopo,itime)
        end if

        gtime = gtime + dt

        if ((((itime/srest*srest.eq.itime).and.(itime.ne.(ntime))).or. &
             (itime.eq.(ntime))).and.(wri_res)) then
!-------- write restart file
          simtime = gtime
          if(amcontrol) then
            call zeitdata_out
          end if
          call writerest(nl,q,q1,mut)
        end if

        if ((itime/iprint*iprint.eq.itime).or.(itime.eq.ntime)) then
          if(amcontrol) then
            write(37,100) itime,resrms(1),resrms(2),resrms(3),resrms(5)
            write(*,100)  itime,resrms(1),resrms(2),resrms(3),resrms(5)
          end if
        end if

!---- increment RK cycle
      end do    

 100  format(i8,4(1x,e13.5))

      return
      end

!-----------------------------------------------------------------------
      subroutine eig_eul(q,sim,sip,sjm,sjp,skm,skp,xdot,ydot,zdot,lami, &
                         lamj,lamk)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      real (kind=cosa_real) q(5),sim(4),sip(4),sjm(4),sjp(4),skm(4),skp(4), &
           xdot,ydot,zdot
      real (kind=cosa_real) p,rho,u,v,w,a,xix,xiy,xiz,etax,etay,etaz,zetax, &
           zetay,zetaz,axi,aeta,azeta,lami,lamj,lamk


      rho = q(1)
      u   = q(2)
      v   = q(3)
      w   = q(4)
      p   = q(5)
      a   = sqrt(gamma*p/rho)

!---- msc, 14 June 2013: Below, the quantities 
!          (xix, xiy, xiz, etax, etay, etaz, zetax, zetay, zetaz)
!          are actually the true metrics multiplied by the cell volume.

      xix   = 0.5d0 * (sim(1)*sim(4) + sip(1)*sip(4))
      xiy   = 0.5d0 * (sim(2)*sim(4) + sip(2)*sip(4))
      xiz   = 0.5d0 * (sim(3)*sim(4) + sip(3)*sip(4))

      etax  = 0.5d0 * (sjm(1)*sjm(4) + sjp(1)*sjp(4))
      etay  = 0.5d0 * (sjm(2)*sjm(4) + sjp(2)*sjp(4))
      etaz  = 0.5d0 * (sjm(3)*sjm(4) + sjp(3)*sjp(4))

      zetax = 0.5d0 * (skm(1)*skm(4) + skp(1)*skp(4))
      zetay = 0.5d0 * (skm(2)*skm(4) + skp(2)*skp(4))
      zetaz = 0.5d0 * (skm(3)*skm(4) + skp(3)*skp(4))

      axi   = (sim(4)+sip(4)) / 2
      aeta  = (sjm(4)+sjp(4)) / 2
      azeta = (skm(4)+skp(4)) / 2

      if (moving) then
        u   = u - xdot
        v   = v - ydot
        w   = w - zdot
      end if

      lami = abs(u*xix    + v*xiy    + w*xiz)   + a*axi
      lamj = abs(u*etax   + v*etay   + w*etaz ) + a*aeta
      lamk = abs(u*zetax  + v*zetay  + w*zetaz) + a*azeta

      return
      end

!-----------------------------------------------------------------------
      subroutine eig_ns(mui,muj,muk,muti,mutj,mutk,q,sim,sip,sjm,sjp, &
                        skm,skp,xdot,ydot,zdot,vol,lami,lamj,lamk, &
                        lamvi,lamvj,lamvk)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      real (kind=cosa_real) q(5),sim(4),sip(4),sjm(4),sjp(4),skm(4),skp(4), &
           xdot,ydot,zdot
      real (kind=cosa_real) p,rho,t,a,u,v,w,mui,muj,muk,muti,mutj,mutk,axi,aeta, &
           azeta,xix,xiy,xiz,etax,etay,etaz,zetax,zetay,zetaz,vol, &
           lami,lamj,lamk,lamvi,lamvj,lamvk

      rho  = q(1)
      u    = q(2)
      v    = q(3)
      w    = q(4)
      p    = q(5)
      t    = gamma*p/rho
      a    = sqrt(t)

      xix   = 0.5d0 * (sim(1)*sim(4) + sip(1)*sip(4))
      xiy   = 0.5d0 * (sim(2)*sim(4) + sip(2)*sip(4))
      xiz   = 0.5d0 * (sim(3)*sim(4) + sip(3)*sip(4))

      etax  = 0.5d0 * (sjm(1)*sjm(4) + sjp(1)*sjp(4))
      etay  = 0.5d0 * (sjm(2)*sjm(4) + sjp(2)*sjp(4))
      etaz  = 0.5d0 * (sjm(3)*sjm(4) + sjp(3)*sjp(4))

      zetax = 0.5d0 * (skm(1)*skm(4) + skp(1)*skp(4))
      zetay = 0.5d0 * (skm(2)*skm(4) + skp(2)*skp(4))
      zetaz = 0.5d0 * (skm(3)*skm(4) + skp(3)*skp(4))

      axi   = (sim(4)+sip(4)) / 2
      aeta  = (sjm(4)+sjp(4)) / 2
      azeta = (skm(4)+skp(4)) / 2

      if (moving) then
        u   = u - xdot
        v   = v - ydot
        w   = w - zdot
      end if

      lami = abs(u*xix    + v*xiy    + w*xiz)   + a*axi
      lamj = abs(u*etax   + v*etay   + w*etaz ) + a*aeta
      lamk = abs(u*zetax  + v*zetay  + w*zetaz) + a*azeta

      lamvi = cdff * axi  *axi   * gamma*machfs/reyno * &
                         (mui/pranl + muti/prant) / (rho * vol)
      lamvj = cdff * aeta *aeta  * gamma*machfs/reyno * &
                         (muj/pranl + mutj/prant) / (rho * vol)
      lamvk = cdff * azeta*azeta * gamma*machfs/reyno * &
                         (muk/pranl + mutk/prant) / (rho * vol)

      return
      end

!-----------------------------------------------------------------------
      subroutine eig_eul_p(q,sim,sip,sjm,sjp,skm,skp,xdot,ydot,zdot, &
                           cutoff,lami,lamj,lamk,nl)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl
      real (kind=cosa_real) q(5),sim(4),sip(4),sjm(4),sjp(4),skm(4),skp(4), &
           xdot,ydot,zdot
      real (kind=cosa_real) p,rho,u,v,w,a,xix,xiy,xiz,etax,etay,etaz,zetax, &
           zetay,zetaz,axi,aeta,azeta,lami,lamj,lamk,cutoff,q2,a2,mp2, &
           ump2

      rho = q(1)
      u   = q(2)
      v   = q(3)
      w   = q(4)
      p   = q(5)
      a   = sqrt(gamma*p/rho)

      if (moving.and.(lsp_vfr.eq.1)) then
        q2    = (u-xdot)**2 + (v-ydot)**2 + (w-zdot)**2
      else
        q2    = u**2 + v**2 + w**2
      end if

      a2    = a*a
      mp2   = q2/a2
      ump2  = (uprv / a)**2
      mp2   = dmin1(dmax1(mp2,ump2,cutoff,epsp(nl)**2),1.d0)

!---- msc, 14 June 2013: Below, the quantities
!          (xix, xiy, xiz, etax, etay, etaz, zetax, zetay, zetaz)
!          are actually the true metrics multiplied by the cell volume.

      xix   = 0.5d0 * (sim(1)*sim(4) + sip(1)*sip(4))
      xiy   = 0.5d0 * (sim(2)*sim(4) + sip(2)*sip(4))
      xiz   = 0.5d0 * (sim(3)*sim(4) + sip(3)*sip(4))

      etax  = 0.5d0 * (sjm(1)*sjm(4) + sjp(1)*sjp(4))
      etay  = 0.5d0 * (sjm(2)*sjm(4) + sjp(2)*sjp(4))
      etaz  = 0.5d0 * (sjm(3)*sjm(4) + sjp(3)*sjp(4))

      zetax = 0.5d0 * (skm(1)*skm(4) + skp(1)*skp(4))
      zetay = 0.5d0 * (skm(2)*skm(4) + skp(2)*skp(4))
      zetaz = 0.5d0 * (skm(3)*skm(4) + skp(3)*skp(4))

      axi   = (sim(4)+sip(4)) / 2
      aeta  = (sjm(4)+sjp(4)) / 2
      azeta = (skm(4)+skp(4)) / 2

      if (moving) then
        u   = u - xdot
        v   = v - ydot
        w   = w - zdot
      end if

      lami = ( (1+mp2) * abs(u*  xix+v*  xiy+w*  xiz) + &
               sqrt((u*  xix+v*  xiy+w*  xiz)**2*(1-mp2)**2 + &
                    4*a2*mp2*  axi**2) ) / 2
      lamj = ( (1+mp2) * abs(u* etax+v* etay+w* etaz) + &
               sqrt((u* etax+v* etay+w* etaz)**2*(1-mp2)**2 + &
                    4*a2*mp2* aeta**2) ) / 2
      lamk = ( (1+mp2) * abs(u*zetax+v*zetay+w*zetaz) + &
               sqrt((u*zetax+v*zetay+w*zetaz)**2*(1-mp2)**2 + &
                    4*a2*mp2*azeta**2) ) / 2

      return
      end

!-----------------------------------------------------------------------
      subroutine eig_ns_p(mui,muj,muk,muti,mutj,mutk,q,sim,sip,sjm,sjp, &
                          skm,skp,xdot,ydot,zdot,vol,cutoff,lami,lamj, &
                          lamk,lamvi,lamvj,lamvk,nl)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl

      real (kind=cosa_real) q(5),sim(4),sip(4),sjm(4),sjp(4),skm(4),skp(4), &
           xdot,ydot,zdot
      real (kind=cosa_real) p,rho,t,a,u,v,w,mui,muj,muk,muti,mutj,mutk,axi,aeta, &
           azeta,xix,xiy,xiz,etax,etay,etaz,zetax,zetay,zetaz,vol, &
           lami,lamj,lamk,lamvi,lamvj,lamvk,cutoff,q2,a2,mp2,ump2

      rho  = q(1)
      u    = q(2)
      v    = q(3)
      w    = q(4)
      p    = q(5)
      t    = gamma*p/rho
      a    = sqrt(t)

      if (moving.and.(lsp_vfr.eq.1)) then
        q2    = (u-xdot)**2 + (v-ydot)**2 + (w-zdot)**2
      else
        q2    = u**2 + v**2 + w**2
      end if

      a2    = a*a
      mp2   = q2/a2
      ump2  = (uprv / a)**2
      mp2   = dmin1(dmax1(mp2,ump2,cutoff,epsp(nl)**2),1.d0)

!---- msc, 14 June 2013: Below, the quantities
!          (xix, xiy, xiz, etax, etay, etaz, zetax, zetay, zetaz)
!          are actually the true metrics multiplied by the cell volume.

      xix   = 0.5d0 * (sim(1)*sim(4) + sip(1)*sip(4))
      xiy   = 0.5d0 * (sim(2)*sim(4) + sip(2)*sip(4))
      xiz   = 0.5d0 * (sim(3)*sim(4) + sip(3)*sip(4))

      etax  = 0.5d0 * (sjm(1)*sjm(4) + sjp(1)*sjp(4))
      etay  = 0.5d0 * (sjm(2)*sjm(4) + sjp(2)*sjp(4))
      etaz  = 0.5d0 * (sjm(3)*sjm(4) + sjp(3)*sjp(4))

      zetax = 0.5d0 * (skm(1)*skm(4) + skp(1)*skp(4))
      zetay = 0.5d0 * (skm(2)*skm(4) + skp(2)*skp(4))
      zetaz = 0.5d0 * (skm(3)*skm(4) + skp(3)*skp(4))

      axi   = (sim(4)+sip(4)) / 2
      aeta  = (sjm(4)+sjp(4)) / 2
      azeta = (skm(4)+skp(4)) / 2

      if (moving) then
        u   = u - xdot
        v   = v - ydot
        w   = w - zdot
      end if

      lami = ( (1+mp2) * abs(u*  xix+v*  xiy+w*  xiz) + &
               sqrt((u*  xix+v*  xiy+w*  xiz)**2*(1-mp2)**2 + &
                    4*a2*mp2*  axi**2) ) / 2
      lamj = ( (1+mp2) * abs(u* etax+v* etay+w* etaz) + &
               sqrt((u* etax+v* etay+w* etaz)**2*(1-mp2)**2 + &
                    4*a2*mp2* aeta**2) ) / 2
      lamk = ( (1+mp2) * abs(u*zetax+v*zetay+w*zetaz) + &
               sqrt((u*zetax+v*zetay+w*zetaz)**2*(1-mp2)**2 + &
                    4*a2*mp2*azeta**2) ) / 2

      lamvi = cdff * axi  *axi   * gamma*machfs/reyno * &
                         (mui/pranl + muti/prant) / (rho * vol)
      lamvj = cdff * aeta *aeta  * gamma*machfs/reyno * &
                         (muj/pranl + mutj/prant) / (rho * vol)
      lamvk = cdff * azeta*azeta * gamma*machfs/reyno * &
                         (muk/pranl + mutk/prant) / (rho * vol)

      return
      end

!-----------------------------------------------------------------------
      subroutine resg_rms(nl,q,res,rhs,vol)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_variables
      use cosa_precision

      implicit none
!

      integer(kind=cosa_int) nl,iblock,imax,jmax,kmax,iq,ivol,ipde,ncell
      real(kind=cosa_real) q(*),res(*),rhs(*),vol(*)
      real(kind=cosa_real) resrms_b(8)

      ncell = 0

      do ipde=1,npde+1
        resrms(ipde) = 0
      end do

      do iblock = 1,mynblocks
        imax   = i_imax     (iblock,nl)
        jmax   = j_jmax     (iblock,nl)
        kmax   = k_kmax     (iblock,nl)
        iq     = 1 + off_p3 (iblock,nl) * npde
        ivol   = 1 + off_p1 (iblock,nl)
        call resg_brms(q(iq),res(iq),rhs(iq),vol(ivol),resrms_b,imax, &
                       jmax,kmax,npde)
        do ipde = 1,npde+1
          resrms(ipde) = resrms(ipde) + resrms_b(ipde)
        end do
        ncell = ncell + (imax-1) * (jmax-1) * (kmax-1)
      end do

      call combineresults(resrms,(npde+1),ncell)

      do ipde = 1,npde
        resrms(ipde) = sqrt( resrms(ipde)/ncell )
      end do
      resrms(npde+1) = sqrt( resrms(npde+1)/ncell/npde )

      return
      end

!-----------------------------------------------------------------------
      subroutine resg_brms(q,res,rhs,vol,resrms_b,imax,jmax,kmax,npde)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none
!
      integer(kind=cosa_int) imax,jmax,kmax,npde
      integer(kind=cosa_int) im,jm,km,i,j,k,ipde
      real(kind=cosa_real) fact
      real(kind=cosa_real) &
          q  (-1:imax+1,-1:jmax+1,-1:kmax+1,npde), &
          res(-1:imax+1,-1:jmax+1,-1:kmax+1,npde), &
          rhs(-1:imax+1,-1:jmax+1,-1:kmax+1,npde), &
          vol( 0:imax  , 0:jmax  , 0:kmax       )
      real(kind=cosa_real) resrms_b(npde+1)
!
      im = imax-1
      jm = jmax-1
      km = kmax-1

!---- calculate residual rms's

      do ipde=1,npde+1
        resrms_b(ipde) = 0
      end do

      if ((rkimp).and.(irsop.eq.0).and.(nlevel.eq.1)) then

        fact = 1.5d0
        if ((itime.eq.1).and.(.not.dual_prop_start)) then
          fact = 1.0d0
        end if

        do ipde = 1,npde
          do k = 1,km
            do j = 1,jm
              do i = 1,im
                resrms_b(ipde) = resrms_b(ipde) + &
                  (fact*q(i,j,k,ipde)*vol(i,j,k)/dt + &
                   res(i,j,k,ipde) - rhs(i,j,k,ipde))**2
              end do
            end do
          end do
        end do

      else if (.not.rgkuns) then

        do ipde = 1,npde
          do k = 1,km
            do j = 1,jm
              do i = 1,im
                resrms_b(ipde) = resrms_b(ipde) + &
                             ( res(i,j,k,ipde) - rhs(i,j,k,ipde))**2
              end do
            end do
          end do
        end do

      else

        do ipde = 1,npde
          do k = 1,km
            do j = 1,jm
              do i = 1,im
                resrms_b(ipde) = resrms_b(ipde) + res(i,j,k,ipde)**2
              end do
            end do
          end do
        end do

      end if

      do ipde = 1,npde
        resrms_b(npde+1) = resrms_b(npde+1) + resrms_b(ipde)
      end do

      return      
      end 

!-----------------------------------------------------------------------
      subroutine resg_rms_hb(nl,q,work1,vol)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_variables
      use cosa_precision

      implicit none
!
      integer(kind=cosa_int) nl,iblock,imax,jmax,kmax,iq,ivol,ipde,ncell
      real(kind=cosa_real) q(*),work1(*),vol(*)
      real(kind=cosa_real) resrms_b(8)

      ncell = 0

      do ipde=1,npde+1
        resrms(ipde) = 0
      end do

      do iblock = 1,mynblocks
        imax   = i_imax     (iblock,nl)
        jmax   = j_jmax     (iblock,nl)
        kmax   = k_kmax     (iblock,nl)
        iq     = 1 + off_p3 (iblock,nl) * npde * dim5
        ivol   = 1 + off_p1 (iblock,nl)
        call resg_brms_hb(q(iq),work1(iq),vol(ivol),resrms_b,imax,jmax, &
                          kmax,npde,nharms)
        do ipde = 1,npde+1
          resrms(ipde) = resrms(ipde) + resrms_b(ipde)
        end do
        ncell = ncell + (imax-1) * (jmax-1) * (kmax-1)
      end do

      call combineresults(resrms,(npde+1),ncell)

      do ipde = 1,npde
        resrms(ipde) = sqrt( resrms(ipde)/ncell/dim5 )
      end do
      resrms(npde+1) = sqrt( resrms(npde+1)/ncell/dim5/npde )

      return
      end

!-----------------------------------------------------------------------
      subroutine resg_brms_hb(q,work1,vol,resrms_b,imax,jmax,kmax,npde, &
                              nharms)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none
!
      integer(kind=cosa_int) imax,jmax,kmax,npde,nharms
      integer(kind=cosa_int) im,jm,km,i,j,k,ipde,n,ni,nj
      real(kind=cosa_real) ar, fact, fact1, term1
      real(kind=cosa_real) &
          q    (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
          work1(-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
          vol  ( 0:imax  , 0:jmax  , 0:kmax                  )
      real(kind=cosa_real) resrms_b(npde+1)
!
      im = imax-1
      jm = jmax-1
      km = kmax-1

!---- calculate residual rms's

      do ipde=1,npde+1
        resrms_b(ipde) = 0
      end do

      if ((rkimp).and.(irsop.eq.0).and.(nlevel.eq.1)) then

!---- this adds the hb coupled source terms. After defining res
!     with size (1:imax,1:jmax) the following bit can be 
!     done by src_HB

        do ni = 0,2*nharms
          do nj = 0,2*nharms
            do ipde = 1,npde
              do k = 1,kmax-1
                do j = 1,jmax-1
                  do i = 1,imax-1
                    work1(i,j,k,ipde,ni) = work1(i,j,k,ipde,ni) + &
                      omega * dhb(ni,nj)*q(i,j,k,ipde,nj)*vol(i,j,k)
                  end do
                end do
              end do
            end do
          end do
        end do

      end if

      do n = 0,2*nharms
        do ipde = 1,npde
          do k = 1,km
            do j = 1,jm
              do i = 1,im
                resrms_b(ipde) = resrms_b(ipde) + work1(i,j,k,ipde,n)**2
              end do
            end do
          end do
        end do
      end do

      do ipde = 1,npde
        resrms_b(npde+1) = resrms_b(npde+1) + resrms_b(ipde)
      end do

      return      
      end 


!-----------------------------------------------------------------------
      subroutine deltat(nl,ic,dtvol,dtvolt,lambdas,lambdasv,betaxi, &
                        betaeta,betazeta,tbetaxi,tbetaeta,tbetazeta)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,ic,iblock,imax,jmax,kmax,iq,ixyz,idtv
      real (kind=cosa_real) &
           dtvol(*),dtvolt(*),lambdas(*),lambdasv(*),betaxi(*), &
           betaeta(*),betazeta(*),tbetaxi(*),tbetaeta(*),tbetazeta(*)

      do iblock = 1,mynblocks
        imax   = i_imax     (iblock,nl)
        jmax   = j_jmax     (iblock,nl)
        kmax   = k_kmax     (iblock,nl)
        iq     = 1 + off_p3 (iblock,nl) * npde * dim5
        idtv   = 1 + off_p1 (iblock,nl) *        dim5
        ixyz   = 1 + off_p2 (iblock,nl) *        dim5h
        call bdeltat(nl,ic,dtvol(idtv),dtvolt(idtv),lambdas(iq), &
          lambdasv(iq),betaxi(idtv),betaeta(idtv),betazeta(idtv), &
          tbetaxi(idtv),tbetaeta(idtv),tbetazeta(idtv),imax,jmax,kmax, &
          npde,nharms,lmet)
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine bdeltat(nl,ic,dtvol,dtvolt,lambdas, &
                   lambdasv,betaxi,betaeta,betazeta, &
                   tbetaxi,tbetaeta,tbetazeta,imax,jmax,kmax, &
                   npde,nharms,lmet)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,npde,nharms,lmet
      integer(kind=cosa_int) i,j,k,ic,n,nh,nl
      real (kind=cosa_real) &
       lambdas  (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
       lambdasv (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
       betaxi   ( 0:imax   ,0:jmax   ,0:kmax   ,     0:2*nharms), &
       betaeta  ( 0:imax   ,0:jmax   ,0:kmax   ,     0:2*nharms), &
       betazeta ( 0:imax   ,0:jmax   ,0:kmax   ,     0:2*nharms), &
       tbetaxi  ( 0:imax   ,0:jmax   ,0:kmax   ,     0:2*nharms), &
       tbetaeta ( 0:imax   ,0:jmax   ,0:kmax   ,     0:2*nharms), &
       tbetazeta( 0:imax   ,0:jmax   ,0:kmax   ,     0:2*nharms), &
       dtvol    ( 0:imax   ,0:jmax   ,0:kmax   ,     0:2*nharms), &
       dtvolt   ( 0:imax   ,0:jmax   ,0:kmax   ,     0:2*nharms)
      real (kind=cosa_real) sqrt3,dcfli1,dcfli2,dcflit1,dcflit2, &
           lamvi,lamvj,lamvk,lami,lamj,lamk,c1,rji,rki,rij,rkj,rik,rjk, &
           bx,by,bz,bxd,byd,bzd,zeta,phixi,phieta,phizeta


      if (kom.or.kom_bsl.or.kom_sst) then
        if (ramping1.and.((mgit.eq.1).or.(mgit/100*100.eq.mgit))) then

          sqrt3   = sqrt(3.d0)

          dcfli1  = (vcfli(2) - vcfli(1)) * sqrt3**3
          dcfli2  =  vcfli(3) - vcfli(2)

          dcflit1 = (vcflit(2) - vcflit(1)) * sqrt3**3
          dcflit2 =  vcflit(3) - vcflit(2) 

          if (mgit.le.ncycr(1)) then
            cfl  = vcfli(1)
            cflt = vcflit(1)
          else if ((mgit.gt.ncycr(1)).and.(mgit.le.ncycr(2))) then
            cfl  = vcfli(1)  + dcfli1  * &
                    (float(mgit-ncycr(1))/(ncycr(2)-ncycr(1))/sqrt3)**3
            cflt = vcflit(1) + dcflit1 * &
                    (float(mgit-ncycr(1))/(ncycr(2)-ncycr(1))/sqrt3)**3
          else if ((mgit.gt.ncycr(2)).and.(mgit.le.ncycr(3))) then
            cfl  = vcfli(2)  + dcfli2  * &
                    (float(mgit-ncycr(2))/(ncycr(3)-ncycr(2)))
            cflt = vcflit(2) + dcflit2 * &
                    (float(mgit-ncycr(2))/(ncycr(3)-ncycr(2)))
          else if (mgit.gt.ncycr(3)) then
            cfl  = vcfli(3)
            cflt = vcflit(3)
          end if

        else if (ramping2) then

          dcfli1  = vcfli(2) - vcfli(1)
          dcfli2  = vcfli(3) - vcfli(2)

          dcflit1 = vcflit(2) - vcflit(1)
          dcflit2 = vcflit(3) - vcflit(2) 

          if (mgit.le.ncycr(0)) then
            cfl  = vcfli(1)
            cflt = vcflit(1)
          else if ((mgit.gt.ncycr(0)).and.(mgit.le.ncycr(1))) then
            cfl  = vcfli(1)  + dcfli1  * &
                    (float(mgit-ncycr(0))/(ncycr(1)-ncycr(0)))
            cflt = vcflit(1) + dcflit1 * &
                    (float(mgit-ncycr(0))/(ncycr(1)-ncycr(0)))
          else if ((mgit.gt.ncycr(1)).and.(mgit.le.ncycr(2))) then
            cfl  = vcfli(2)
            cflt = vcflit(2)
          else if ((mgit.gt.ncycr(2)).and.(mgit.le.ncycr(3))) then
            cfl  = vcfli(2)  + dcfli2  * &
                    (float(mgit-ncycr(2))/(ncycr(3)-ncycr(2)))
            cflt = vcflit(2) + dcflit2 * &
                    (float(mgit-ncycr(2))/(ncycr(3)-ncycr(2)))
          else if (mgit.gt.ncycr(3)) then
            cfl  = vcfli(3)
            cflt = vcflit(3)
          end if

        end if
      end if

!dbg  write(123,*) mgit, cfl, cflt,rat

!     zeta = 2.d0 / 3.d0
      zeta = 0.4d0
      c1   = 40.d0

!---- calculate local time step

      if ( (.not.viscous).or. &
           ((viscous).and.((irsop.eq.5).or.(irsop.eq.6))) )     then

        do n = 0,2*nharms
          nh = n*hbmove

          do k = 1,kmax-1
            do j = 1,jmax-1
              do i = 1,imax-1

                dtvol(i,j,k,n) = cfl / &
                  ( lambdas(i,j,k,1,n) + lambdas(i,j,k,2,n) + &
                    lambdas(i,j,k,3,n) )

              end do
            end do
          end do

          if (kom.or.kom_bsl.or.kom_sst) then
            do k = 1,kmax-1
              do j = 1,jmax-1
                do i = 1,imax-1
                  dtvolt(i,j,k,n) = cflt / cfl * dtvol(i,j,k,n)
                end do
              end do
            end do
          end if

        end do

      else if ( ((viscous).and.(irsop.ne.5)).or. &
                ((viscous).and.(irsop.ne.6)) )     then

        do n = 0,2*nharms
          nh = n*hbmove

          do k = 1,kmax-1
            do j = 1,jmax-1
              do i = 1,imax-1

                dtvol(i,j,k,n) = cfl / &
                  ( lambdas(i,j,k,1,n) + lambdas(i,j,k,2,n) + &
                    lambdas(i,j,k,3,n) + lambdasv(i,j,k,1,n) + &
                    lambdasv(i,j,k,2,n) + lambdasv(i,j,k,3,n) )

              end do
            end do
          end do

          if (kom.or.kom_bsl.or.kom_sst) then
            do k = 1,kmax-1
              do j = 1,jmax-1
                do i = 1,imax-1
                  dtvolt(i,j,k,n) = cflt / cfl * dtvol(i,j,k,n)
                end do
              end do
            end do
          end if

        end do

      end if

!------ calculate smoothing coefficients

      if (irsop.eq.6) then

        do n = 0,2*nharms
          nh = n*hbmove

          do k = 1,kmax-1
            do j = 1,jmax-1
              do i = 1,imax-1

                lami = lambdas(i,j,k,1,n)
                lamj = lambdas(i,j,k,2,n)
                lamk = lambdas(i,j,k,3,n)
                lamvi = lambdasv(i,j,k,1,n)
                lamvj = lambdasv(i,j,k,2,n)
                lamvk = lambdasv(i,j,k,3,n)

                rji = lamj / lami
                rki = lamk / lami
                rij = lami / lamj
                rkj = lamk / lamj
                rik = lami / lamk
                rjk = lamj / lamk

                bx   = 0.25d0*((rat / (1 + psirs*(rji+rki)))**2 - 1)
                bx   = dmax1(0.d0,bx)
                bxd  = 0.25d0 * c1 * lamvi / (lami+lamj+lamk)
                betaxi(i,j,k,n)   = dmax1(bx,bxd)

                by   = 0.25d0*((rat / (1 + psirs*(rij+rkj)))**2 - 1)
                by   = dmax1(0.d0,by)
                byd  = 0.25d0 * c1 * lamvj / (lami+lamj+lamk)
                betaeta(i,j,k,n)  = dmax1(by,byd)

                bz   = 0.25d0*((rat / (1 + psirs*(rik+rjk)))**2 - 1)
                bz   = dmax1(0.d0,bz)
                bzd  = 0.25d0 * c1 * lamvk / (lami+lamj+lamk)
                betazeta(i,j,k,n) = dmax1(bz,bzd)

              end do
            end do
          end do

          if (kom.or.kom_bsl.or.kom_sst) then

            do k = 1,kmax-1
              do j = 1,jmax-1
                do i = 1,imax-1

                  lami = lambdas(i,j,k,1,n)
                  lamj = lambdas(i,j,k,2,n)
                  lamk = lambdas(i,j,k,3,n)
                  lamvi = lambdasv(i,j,k,1,n)
                  lamvj = lambdasv(i,j,k,2,n)
                  lamvk = lambdasv(i,j,k,3,n)

                  rji = lamj / lami
                  rki = lamk / lami
                  rij = lami / lamj
                  rkj = lamk / lamj
                  rik = lami / lamk
                  rjk = lamj / lamk

                  bx   = 0.25d0*((ratt / (1 + psirs*(rji+rki)))**2 -1)
                  bx   = dmax1(0.d0,bx)
                  bxd  = 0.25d0 * c1 * lamvi / (lami+lamj+lamk)
                  tbetaxi(i,j,k,n)   = dmax1(bx,bxd)

                  by   = 0.25d0*((ratt / (1 + psirs*(rij+rkj)))**2 -1)
                  by   = dmax1(0.d0,by)
                  byd  = 0.25d0 * c1 * lamvj / (lami+lamj+lamk)
                  tbetaeta(i,j,k,n)  = dmax1(by,byd)

                  bz   = 0.25d0*((ratt / (1 + psirs*(rik+rjk)))**2 -1)
                  bz   = dmax1(0.d0,bz)
                  bzd  = 0.25d0 * c1 * lamvk / (lami+lamj+lamk)
                  tbetazeta(i,j,k,n) = dmax1(bz,bzd)

                end do
              end do
            end do

          end if

        end do

      else if (irsop.eq.5) then

        do n = 0,2*nharms
          nh = n*hbmove

          do k = 1,kmax-1
            do j = 1,jmax-1
              do i = 1,imax-1

                lami = lambdas(i,j,k,1,n)
                lamj = lambdas(i,j,k,2,n)
                lamk = lambdas(i,j,k,3,n)
                lamvi = lambdasv(i,j,k,1,n)
                lamvj = lambdasv(i,j,k,2,n)
                lamvk = lambdasv(i,j,k,3,n)

                phixi   = 1 + (lamj/lami)**zeta + (lamk/lami)**zeta
                phieta  = 1 + (lami/lamj)**zeta + (lamk/lamj)**zeta
                phizeta = 1 + (lami/lamk)**zeta + (lamj/lamk)**zeta

                bx   = ((rat * lami * phixi     / &
                         (lami+lamj+lamk))**2 -1) / 4
                bx   = dmax1(0.d0,bx)
                bxd  = c1 / 4 * lamvi / (lami+lamj+lamk)
                betaxi(i,j,k,n)   = dmax1(bx,bxd)

                by   = ((rat * lamj * phieta   / &
                         (lami+lamj+lamk))**2 -1) / 4
                by   = dmax1(0.d0,by)
                byd  = c1 / 4 * lamvj / (lami+lamj+lamk)
                betaeta(i,j,k,n)  = dmax1(by,byd)

                bz   = ((rat * lamk * phizeta  / &
                         (lami+lamj+lamk))**2 -1) / 4
                bz   = dmax1(0.d0,bz)
                bzd  = c1 / 4 * lamvk / (lami+lamj+lamk)
                betazeta(i,j,k,n) = dmax1(bz,bzd)

              end do
            end do
          end do

          if (kom.or.kom_bsl.or.kom_sst) then

            do k = 1,kmax-1
              do j = 1,jmax-1
                do i = 1,imax-1

                  lami = lambdas(i,j,k,1,n)
                  lamj = lambdas(i,j,k,2,n)
                  lamk = lambdas(i,j,k,3,n)
                  lamvi = lambdasv(i,j,k,1,n)
                  lamvj = lambdasv(i,j,k,2,n)
                  lamvk = lambdasv(i,j,k,3,n)

                  phixi   = 1 + (lamj/lami)**zeta + (lamk/lami)**zeta
                  phieta  = 1 + (lami/lamj)**zeta + (lamk/lamj)**zeta
                  phizeta = 1 + (lami/lamk)**zeta + (lamj/lamk)**zeta

                  bx   = ((ratt * lami * phixi     / &
                           (lami+lamj+lamk))**2 -1) / 4
                  bx   = dmax1(0.d0,bx)
                  bxd  = c1 / 4 * lamvi / (lami+lamj+lamk)
                  tbetaxi(i,j,k,n)   = dmax1(bx,bxd)

                  by   = ((ratt * lamj * phieta   / &
                           (lami+lamj+lamk))**2 -1) / 4
                  by   = dmax1(0.d0,by)
                  byd  = c1 / 4 * lamvj / (lami+lamj+lamk)
                  tbetaeta(i,j,k,n)  = dmax1(by,byd)

                  bz   = ((ratt * lamk * phizeta  / &
                         (lami+lamj+lamk))**2 -1) / 4
                  bz   = dmax1(0.d0,bz)
                  bzd  = c1 / 4 * lamvk / (lami+lamj+lamk)
                  tbetazeta(i,j,k,n) = dmax1(bz,bzd)

                end do
              end do
            end do

          end if

        end do

      else if (irsop.eq.2) then

        do n = 0,2*nharms
          nh = n*hbmove

          do k = 1,kmax-1
            do j = 1,jmax-1
              do i = 1,imax-1

                lami = lambdas(i,j,k,1,n)
                lamj = lambdas(i,j,k,2,n)
                lamk = lambdas(i,j,k,3,n)

                phixi   = 1 + (lamj/lami)**zeta + (lamk/lami)**zeta
                phieta  = 1 + (lami/lamj)**zeta + (lamk/lamj)**zeta
                phizeta = 1 + (lami/lamk)**zeta + (lamj/lamk)**zeta

                bx   = ((rat * lami * phixi     / &
                         (lami+lamj+lamk))**2 -1) / 4
                betaxi(i,j,k,n)   = dmax1(0.d0,bx)

                by   = ((rat * lamj * phieta   / &
                         (lami+lamj+lamk))**2 -1) / 4
                betaeta(i,j,k,n)  = dmax1(0.d0,by)

                bz   = ((rat * lamk * phizeta  / &
                         (lami+lamj+lamk))**2 -1) / 4
                betazeta(i,j,k,n) = dmax1(0.d0,bz)

              end do
            end do
          end do

          if (kom.or.kom_bsl.or.kom_sst) then

            do k = 1,kmax-1
              do j = 1,jmax-1
                do i = 1,imax-1

                  lami = lambdas(i,j,k,1,n)
                  lamj = lambdas(i,j,k,2,n)
                  lamk = lambdas(i,j,k,3,n)

                  phixi   = 1 + (lamj/lami)**zeta + (lamk/lami)**zeta
                  phieta  = 1 + (lami/lamj)**zeta + (lamk/lamj)**zeta
                  phizeta = 1 + (lami/lamk)**zeta + (lamj/lamk)**zeta

                  bx   = ((ratt * lami * phixi     / &
                           (lami+lamj+lamk))**2 -1) / 4
                  tbetaxi(i,j,k,n)   = dmax1(0.d0,bx)

                  by   = ((ratt * lamj * phieta   / &
                           (lami+lamj+lamk))**2 -1) / 4
                  tbetaeta(i,j,k,n)  = dmax1(0.d0,by)

                  bz   = ((ratt * lamk * phizeta  / &
                           (lami+lamj+lamk))**2 -1) / 4
                  tbetazeta(i,j,k,n) = dmax1(0.d0,bz)

                end do
              end do
            end do

          end if

        end do

      else if (irsop.eq.3) then

        do n = 0,2*nharms
          nh = n*hbmove

          do k = 1,kmax-1
            do j = 1,jmax-1
              do i = 1,imax-1

                lami = lambdas(i,j,k,1,n)
                lamj = lambdas(i,j,k,2,n)
                lamk = lambdas(i,j,k,3,n)

                rji = lamj / lami
                rki = lamk / lami
                rij = lami / lamj
                rkj = lamk / lamj
                rik = lami / lamk
                rjk = lamj / lamk

                bx   = 0.25d0*((rat / (1 + psirs*(rji+rki)))**2 - 1)
                betaxi(i,j,k,n)   = dmax1(0.d0,bx)

                by   = 0.25d0*((rat / (1 + psirs*(rij+rkj)))**2 - 1)
                betaeta(i,j,k,n)  = dmax1(0.d0,by)

                bz   = 0.25d0*((rat / (1 + psirs*(rik+rjk)))**2 - 1)
                betazeta(i,j,k,n) = dmax1(0.d0,bz)

              end do
            end do
          end do

          if (kom.or.kom_bsl.or.kom_sst) then

            do k = 1,kmax-1
              do j = 1,jmax-1
                do i = 1,imax-1

                  lami = lambdas(i,j,k,1,n)
                  lamj = lambdas(i,j,k,2,n)
                  lamk = lambdas(i,j,k,3,n)

                  rji = lamj / lami
                  rki = lamk / lami
                  rij = lami / lamj
                  rkj = lamk / lamj
                  rik = lami / lamk
                  rjk = lamj / lamk

                  bx  = 0.25d0*((ratt / (1 + psirs*(rji+rki)))**2 - 1)
                  tbetaxi(i,j,k,n)   = dmax1(0.d0,bx)

                  by  = 0.25d0*((ratt / (1 + psirs*(rij+rkj)))**2 - 1)
                  tbetaeta(i,j,k,n)  = dmax1(0.d0,by)

                  bz  = 0.25d0*((ratt / (1 + psirs*(rik+rjk)))**2 - 1)
                  tbetazeta(i,j,k,n) = dmax1(0.d0,bz)

                end do
              end do
            end do

          end if

        end do

      end if

      return
      end

!---------------------------------------------------------------------
      subroutine cirs(nl,res,dtvol,betaxi,betaeta,betazeta,dtvolt, &
                      tbetaxi,tbetaeta,tbetazeta,work)
!---------------------------------------------------------------------       

      use common_variables
      use cosa_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,iblock,imax,jmax,kmax,ijkmax,ires,idtv,iwrk
      real (kind=cosa_real) res(*),dtvol(*),betaxi(*),betaeta(*),betazeta(*), &
           dtvolt(*),tbetaxi(*),tbetaeta(*),tbetazeta(*),work(*)

      do iblock = 1,mynblocks
        imax   = i_imax     (iblock,nl)
        jmax   = j_jmax     (iblock,nl)
        kmax   = k_kmax     (iblock,nl)
        ijkmax = ijk_ijkmax (iblock,nl)
        ires   = 1 + off_p3 (iblock,nl) * npde * dim5
        idtv   = 1 + off_p1 (iblock,nl) *        dim5
        iwrk   = 1 + off_1d (iblock,nl) * npde
        call cirs_b(res(ires),dtvol(idtv),betaxi(idtv),betaeta(idtv), &
               betazeta(idtv),dtvolt(idtv),tbetaxi(idtv),tbetaeta(idtv), &
               tbetazeta(idtv),work(iwrk),imax,jmax,kmax,npde,nharms, &
               ijkmax)
      end do

      return
      end

!---------------------------------------------------------------------
      subroutine cirs_b(res,dtvol,betaxi,betaeta,betazeta,dtvolt, &
                   tbetaxi,tbetaeta,tbetazeta,work,imax,jmax,kmax,npde, &
                   nharms,ijkmax)
!---------------------------------------------------------------------       

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,ijkmax,npde,nharms
      integer(kind=cosa_int) im,i,jm,j,km,k,ipde,n
      real (kind=cosa_real) &
           res      (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           dtvol    ( 0:imax  , 0:jmax  , 0:kmax  ,     0:2*nharms), &
           betaxi   ( 0:imax  , 0:jmax  , 0:kmax  ,     0:2*nharms), &
           betaeta  ( 0:imax  , 0:jmax  , 0:kmax  ,     0:2*nharms), &
           betazeta ( 0:imax  , 0:jmax  , 0:kmax  ,     0:2*nharms), &
           dtvolt   ( 0:imax  , 0:jmax  , 0:kmax  ,     0:2*nharms), &
           tbetaxi  ( 0:imax  , 0:jmax  , 0:kmax  ,     0:2*nharms), &
           tbetaeta ( 0:imax  , 0:jmax  , 0:kmax  ,     0:2*nharms), &
           tbetazeta( 0:imax  , 0:jmax  , 0:kmax  ,     0:2*nharms), &
           work     ( ijkmax-1, npde)

      real (kind=cosa_real) epsirc

      im = imax-1
      jm = jmax-1
      km = kmax-1

      do n = 0,2*nharms
        do ipde = 1,5
          do k = 1,km
            do j = 1,jm
              do i = 1,im
                res(i,j,k,ipde,n) = res(i,j,k,ipde,n) * dtvol(i,j,k,n)
              end do
            end do
          end do
        end do
      end do

      if (kom.or.kom_bsl.or.kom_sst) then
        do n = 0,2*nharms
          do ipde = 6,npde
            do k = 1,km
              do j = 1,jm
                do i = 1,im
                  res(i,j,k,ipde,n) = res(i,j,k,ipde,n) * &
                                      dtvolt(i,j,k,n)
                end do
              end do
            end do
          end do
        end do
      end if

      if (irsop.eq.1) then

        epsirc = 0.25d0 * (rat**2 - 1.d0)

!------ xi direction

        do n = 0,2*nharms
          do ipde = 1,npde
            do k = 1,km
              do j = 1,jm
                work(1,1) = 0
                work(1,2) = 1 + 2*epsirc
                work(1,3) = - epsirc
                work(1,4) = res(1,j,k,ipde,n)
                do i = 2,im-1
                  work(i,1) = - epsirc
                  work(i,2) = 1 + 2*epsirc
                  work(i,3) = - epsirc
                  work(i,4) = res(i,j,k,ipde,n)
                end do
                work(im,1) = - epsirc
                work(im,2) = 1 + 2*epsirc
                work(im,3) = 0
                work(im,4) = res(im,j,k,ipde,n)
                call tridi(work(1,1),work(1,2),work(1,3),work(1,4), &
                           ijkmax-1,im)
                do i=1,im
                  res(i,j,k,ipde,n) = work(i,4)
                end do
              end do
            end do
          end do
        end do

!------ eta direction

        do n = 0,2*nharms
          do ipde = 1,npde
            do k = 1,km
              do i = 1,im
                work(1,1) = 0
                work(1,2) = 1 + 2*epsirc
                work(1,3) = - epsirc
                work(1,4) = res(i,1,k,ipde,n)
                do j = 2,jm-1
                  work(j,1) = - epsirc
                  work(j,2) = 1 + 2*epsirc
                  work(j,3) = - epsirc
                  work(j,4) = res(i,j,k,ipde,n)
                end do
                work(jm,1) = - epsirc
                work(jm,2) = 1 + 2*epsirc
                work(jm,3) = 0
                work(jm,4) = res(i,jm,k,ipde,n)
                call tridi(work(1,1),work(1,2),work(1,3),work(1,4), &
                           ijkmax-1,jm)
                do j=1,jm
                  res(i,j,k,ipde,n) = work(j,4)
                end do
              end do
            end do
          end do
        end do

!------ zeta direction

        do n = 0,2*nharms
          do ipde = 1,npde
            do j = 1,jm
              do i = 1,im
                work(1,1) = 0
                work(1,2) = 1 + 2*epsirc
                work(1,3) = - epsirc
                work(1,4) = res(i,j,1,ipde,n)
                do k = 2,km-1
                  work(k,1) = - epsirc
                  work(k,2) = 1 + 2*epsirc
                  work(k,3) = - epsirc
                  work(k,4) = res(i,j,k,ipde,n)
                end do
                work(km,1) = - epsirc
                work(km,2) = 1 + 2*epsirc
                work(km,3) = 0
                work(km,4) = res(i,j,km,ipde,n)
                call tridi(work(1,1),work(1,2),work(1,3),work(1,4), &
                           ijkmax-1,km)
                do k=1,km
                  res(i,j,k,ipde,n) = work(k,4)
                end do
              end do
            end do
          end do
        end do

      else if (irsop.ge.2) then

!------ xi direction

        do n = 0,2*nharms
          do ipde = 1,5
            do k = 1,km
              do j = 1,jm
                epsirc = betaxi(1,j,k,n)
                work(1,1) = 0
                work(1,2) = 1 + 2*epsirc
                work(1,3) = - epsirc
                work(1,4) = res(1,j,k,ipde,n)
                do i = 2,im-1
                  epsirc = betaxi(i,j,k,n)
                  work(i,1) = - epsirc
                  work(i,2) = 1 + 2*epsirc
                  work(i,3) = - epsirc
                  work(i,4) = res(i,j,k,ipde,n)
                end do
                epsirc = betaxi(im,j,k,n)
                work(im,1) = - epsirc
                work(im,2) = 1 + 2*epsirc
                work(im,3) = 0
                work(im,4) = res(im,j,k,ipde,n)
                call tridi(work(1,1),work(1,2),work(1,3),work(1,4), &
                           ijkmax-1,im)
                do i=1,im
                  res(i,j,k,ipde,n) = work(i,4)
                end do
              end do
            end do
          end do
        end do

!------ eta direction

        do n = 0,2*nharms
          do ipde = 1,5
            do k = 1,km
              do i = 1,im
                epsirc = betaeta(i,1,k,n)
                work(1,1) = 0
                work(1,2) = 1 + 2*epsirc
                work(1,3) = - epsirc
                work(1,4) = res(i,1,k,ipde,n)
                do j = 2,jm-1
                  epsirc = betaeta(i,j,k,n)
                  work(j,1) = - epsirc
                  work(j,2) = 1 + 2*epsirc
                  work(j,3) = - epsirc
                  work(j,4) = res(i,j,k,ipde,n)
                end do
                epsirc = betaeta(i,jm,k,n)
                work(jm,1) = - epsirc
                work(jm,2) = 1 + 2*epsirc
                work(jm,3) = 0
                work(jm,4) = res(i,jm,k,ipde,n)
                call tridi(work(1,1),work(1,2),work(1,3),work(1,4), &
                           ijkmax-1,jm)
                do j=1,jm
                  res(i,j,k,ipde,n) = work(j,4)
                end do
              end do
            end do
          end do
        end do

!------ zeta direction

        do n = 0,2*nharms
          do ipde = 1,5
            do j = 1,jm
              do i = 1,im
                epsirc = betazeta(i,j,1,n)
                work(1,1) = 0
                work(1,2) = 1 + 2*epsirc
                work(1,3) = - epsirc
                work(1,4) = res(i,j,1,ipde,n)
                do k = 2,km-1
                  epsirc = betazeta(i,j,k,n)
                  work(k,1) = - epsirc
                  work(k,2) = 1 + 2*epsirc
                  work(k,3) = - epsirc
                  work(k,4) = res(i,j,k,ipde,n)
                end do
                epsirc = betazeta(i,j,km,n)
                work(km,1) = - epsirc
                work(km,2) = 1 + 2*epsirc
                work(km,3) = 0
                work(km,4) = res(i,j,km,ipde,n)
                call tridi(work(1,1),work(1,2),work(1,3),work(1,4), &
                           ijkmax-1,km)
                do k=1,km
                  res(i,j,k,ipde,n) = work(k,4)
                end do
              end do
            end do
          end do
        end do

        if (kom.or.kom_bsl.or.kom_sst) then

!-------- xi direction

          do n = 0,2*nharms
            do ipde = 6,npde
              do k = 1,km
                do j = 1,jm
                  epsirc = tbetaxi(1,j,k,n)
                  work(1,1) = 0
                  work(1,2) = 1 + 2*epsirc
                  work(1,3) = - epsirc
                  work(1,4) = res(1,j,k,ipde,n)
                  do i = 2,im-1
                    epsirc = tbetaxi(i,j,k,n)
                    work(i,1) = - epsirc
                    work(i,2) = 1 + 2*epsirc
                    work(i,3) = - epsirc
                    work(i,4) = res(i,j,k,ipde,n)
                  end do
                  epsirc = tbetaxi(im,j,k,n)
                  work(im,1) = - epsirc
                  work(im,2) = 1 + 2*epsirc
                  work(im,3) = 0
                  work(im,4) = res(im,j,k,ipde,n)
                  call tridi(work(1,1),work(1,2),work(1,3),work(1,4), &
                             ijkmax-1,im)
                  do i=1,im
                    res(i,j,k,ipde,n) = work(i,4)
                  end do
                end do
              end do
            end do
          end do

!-------- eta direction

          do n = 0,2*nharms
            do ipde = 6,npde
              do k = 1,km
                do i = 1,im
                  epsirc = tbetaeta(i,1,k,n)
                  work(1,1) = 0
                  work(1,2) = 1 + 2*epsirc
                  work(1,3) = - epsirc
                  work(1,4) = res(i,1,k,ipde,n)
                  do j = 2,jm-1
                    epsirc = tbetaeta(i,j,k,n)
                    work(j,1) = - epsirc
                    work(j,2) = 1 + 2*epsirc
                    work(j,3) = - epsirc
                    work(j,4) = res(i,j,k,ipde,n)
                  end do
                  epsirc = tbetaeta(i,jm,k,n)
                  work(jm,1) = - epsirc
                  work(jm,2) = 1 + 2*epsirc
                  work(jm,3) = 0
                  work(jm,4) = res(i,jm,k,ipde,n)
                  call tridi(work(1,1),work(1,2),work(1,3),work(1,4), &
                             ijkmax-1,jm)
                  do j=1,jm
                    res(i,j,k,ipde,n) = work(j,4)
                  end do
                end do
              end do
            end do
          end do

!-------- zeta direction

          do n = 0,2*nharms
            do ipde = 6,npde
              do j = 1,jm
                do i = 1,im
                  epsirc = tbetazeta(i,j,1,n)
                  work(1,1) = 0
                  work(1,2) = 1 + 2*epsirc
                  work(1,3) = - epsirc
                  work(1,4) = res(i,j,1,ipde,n)
                  do k = 2,km-1
                    epsirc = tbetazeta(i,j,k,n)
                    work(k,1) = - epsirc
                    work(k,2) = 1 + 2*epsirc
                    work(k,3) = - epsirc
                    work(k,4) = res(i,j,k,ipde,n)
                  end do
                  epsirc = tbetazeta(i,j,km,n)
                  work(km,1) = - epsirc
                  work(km,2) = 1 + 2*epsirc
                  work(km,3) = 0
                  work(km,4) = res(i,j,km,ipde,n)
                  call tridi(work(1,1),work(1,2),work(1,3),work(1,4), &
                             ijkmax-1,km)
                  do k=1,km
                    res(i,j,k,ipde,n) = work(k,4)
                  end do
                end do
              end do
            end do
          end do

        end if

      end if

      do n = 0,2*nharms
        do ipde = 1,5
          do k = 1,km
            do j = 1,jm
              do i = 1,im
                res(i,j,k,ipde,n) = res(i,j,k,ipde,n) / dtvol(i,j,k,n)
              end do
            end do
          end do
        end do
      end do

      if (kom.or.kom_bsl.or.kom_sst) then
        do n = 0,2*nharms
          do ipde = 6,npde
            do k = 1,km
              do j = 1,jm
                do i = 1,im
                  res(i,j,k,ipde,n) = res(i,j,k,ipde,n) / &
                                      dtvolt(i,j,k,n)
                end do
              end do
            end do
          end do
        end do
      end if

      return
      end

!---------------------------------------------------------------------
      subroutine scalef(nl,res,dtvol,dtvolt)
!---------------------------------------------------------------------       

      use common_variables
      use cosa_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,iblock,imax,jmax,kmax,ires,idtv
      real (kind=cosa_real) res(*),dtvol(*),dtvolt(*)

      do iblock = 1,mynblocks
        imax   = i_imax     (iblock,nl)
        jmax   = j_jmax     (iblock,nl)
        kmax   = k_kmax     (iblock,nl)
        ires   = 1 + off_p3 (iblock,nl) * npde * dim5
        idtv   = 1 + off_p1 (iblock,nl) *        dim5
        call scalef_b(res(ires),dtvol(idtv),dtvolt(idtv), &
                      imax,jmax,kmax,npde,nharms)
      end do

      return
      end

!---------------------------------------------------------------------
      subroutine scalef_b(res,dtvol,dtvolt,imax,jmax,kmax,npde,nharms)
!---------------------------------------------------------------------       

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,npde,nharms
      integer(kind=cosa_int) im,i,jm,j,km,k,ipde,n
      real (kind=cosa_real) &
           res   (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           dtvol ( 0:imax  , 0:jmax  , 0:kmax  ,     0:2*nharms), &
           dtvolt( 0:imax  , 0:jmax  , 0:kmax  ,     0:2*nharms)

      im = imax-1
      jm = jmax-1
      km = kmax-1

      do n = 0,2*nharms
        do ipde = 1,5
          do k = 1,km
            do j = 1,jm
              do i = 1,im
                res(i,j,k,ipde,n) = res(i,j,k,ipde,n) * dtvol(i,j,k,n)
              end do
            end do
          end do
        end do
      end do

      if (kom.or.kom_bsl.or.kom_sst) then
        do n = 0,2*nharms
          do ipde = 6,npde
            do k = 1,km
              do j = 1,jm
                do i = 1,im
                  res(i,j,k,ipde,n) = res(i,j,k,ipde,n) * dtvolt(i,j,k,n)
                end do
              end do
            end do
          end do
        end do
      end if

      return
      end

!---------------------------------------------------------------------
      subroutine scaleb(nl,res,dtvol,dtvolt)
!---------------------------------------------------------------------       

      use common_variables
      use cosa_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,iblock,imax,jmax,kmax,ires,idtv
      real (kind=cosa_real) res(*),dtvol(*),dtvolt(*)

      do iblock = 1,mynblocks
        imax   = i_imax     (iblock,nl)
        jmax   = j_jmax     (iblock,nl)
        kmax   = k_kmax     (iblock,nl)
        ires   = 1 + off_p3 (iblock,nl) * npde * dim5
        idtv   = 1 + off_p1 (iblock,nl) *        dim5
        call scaleb_b(res(ires),dtvol(idtv),dtvolt(idtv), &
                      imax,jmax,kmax,npde,nharms)
      end do

      return
      end

!---------------------------------------------------------------------
      subroutine scaleb_b(res,dtvol,dtvolt,imax,jmax,kmax,npde,nharms)
!---------------------------------------------------------------------       

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,npde,nharms
      integer(kind=cosa_int) im,i,jm,j,km,k,ipde,n
      real (kind=cosa_real) &
           res   (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           dtvol ( 0:imax  , 0:jmax  , 0:kmax  ,     0:2*nharms), &
           dtvolt( 0:imax  , 0:jmax  , 0:kmax  ,     0:2*nharms)

      im = imax-1
      jm = jmax-1
      km = kmax-1

      do n = 0,2*nharms
        do ipde = 1,5
          do k = 1,km
            do j = 1,jm
              do i = 1,im
                res(i,j,k,ipde,n) = res(i,j,k,ipde,n) / dtvol(i,j,k,n)
              end do
            end do
          end do
        end do
      end do

      if (kom.or.kom_bsl.or.kom_sst) then
        do n = 0,2*nharms
          do ipde = 6,npde
            do k = 1,km
              do j = 1,jm
                do i = 1,im
                  res(i,j,k,ipde,n) = res(i,j,k,ipde,n) / dtvolt(i,j,k,n)
                end do
              end do
            end do
          end do
        end do
      end if

      return
      end
