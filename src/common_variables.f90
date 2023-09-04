!-----------------------------------------------------------------------
!
! Module to hold some of the common blocks from the F77 version of the 
! program. 
!
! Adrian Jackson, EPCC, The University of Edinburgh
! June 2023
!
!-----------------------------------------------------------------------

  module common_variables

    use cosa_precision

    implicit none

    public

    integer, parameter :: mharms=65
    integer, parameter :: msurface=10
    integer, parameter ::m_surblo=2000

    integer(kind=cosa_int) ::irest,srest,lmax,ncycle,iupdt, &
         iprint,ilimf,ilimt,irsop,cutcirs,nsurface,n_surblo(msurface), &
         i_surbl(msurface,m_surblo), & 
         profpar(4),nblock_pos,iblock_pos(msurface*m_surblo),wallmask,ncycr(0:3)
    integer(kind=cosa_int) :: krd(3,3),cyc(3,3),etpfxtyp_f,etpfxtyp_t,entfx_f_w(3),entfx_t_w(3)
    integer(kind=cosa_int) :: ntime,ntime_ramp,nsave,wrifl(3),itime,hbmove,nblade,somega
    integer(kind=cosa_int) :: nlevel,nl_crs,nstart,npre,npost,ncrs,nl_fmg,mgit
    integer(kind=cosa_int) :: nstagerk,nstage
    integer(kind=cosa_int) :: cutoff_type,mixpi,ifixq,lsp_vfr,turcmp
    
    
    real (kind=cosa_real):: xmom,ymom,zmom,gamma,machfs,alpha,beta,reyno,pranl, &
         lngth(3),tref,stemp,ainf,rhoinf,tinf,pinf,pi,cfl,cflt,vcfl(0:3),vcfli(0:3), &
         vcflt(0:3),vcflit(0:3),cdff,psirs,epslimf,epslimt,cntrpy_f,cntrpy_t, &
         entfxctf_f,entfxctf_t,rat,ratt,rkap,eig_cutoff_f, & 
         eig_cutoff_t,toler,resrms(8),qmin(7),omega,omegas,omegatp, &
         betaw,phitp,simtime,dt, &
         zeit(0:2*mharms),dhb(0:2*mharms,0:2*mharms), & 
         ahb(0:2*mharms,0:2*mharms),ehb(0:2*mharms,0:2*mharms), &
         eihb(0:2*mharms,0:2*mharms), &
         xh,yh,dtheta0,dh0x,dh0y,phipp,dtheta(0:2*mharms), &
         yhtp,zhtp,theta0tp,thetatp(0:2*mharms),xrotc,yrotc,zrotc,gtheta, &
         rotmat(2,2),epsp(12),lx,ly,uprv
    real (kind=cosa_real) alphas(5), betas(5), alphadum(5), ar4(4), br4(4), &
         ar5(5)
    real (kind=cosa_real) epst,prant,tkefar,mutfar,omefar,cmu,bstrkw,gkw,bkw, &
         sigk,sigw,roughk,prdlim,bstr,sigk1_bsl,sigk1_sst,sigw1,b1,g1, &
         sigk2,sigw2,b2,g2,a1
    
    character(len=10) :: code
    character(len=72) :: flowtec,surftec,filename
    character(len=200) :: cgns_filename
    
    logical :: aircraft,hawt,persec,vawt,relframe,calfor,exter,viscous, &
         foturb,kom,kom_bsl,kom_sst,ho_rest,prd_all_lvls,wallwilc,wallment, &
         lim_prodtke,lim_prodomega,lim_corr(3),unsteady,rgkuns, &
         dualt,dual_prop_start, &
         found_simtime,harbal,rkimp,moving,rotating,tpitch,pitching, &
         plunging,plupitching,lomach,visprec,debug,wri_res,wri_tec, &
         wldata,blprof,calvort,calcp,calmut,bln_prol,ramping1, &
         ramping2,tecbin,dgn_cell,efx_setup_f,efx_setup_t,eig_limit_f, &
         eig_limit_t,efx_all_f,efx_all_t,eig_lim_all_f,eig_lim_all_t, &
         using_cgns
    
!    common/comm1/gamma,machfs,alpha,beta,rat,ratt,cfl,cflt,vcfl(0:3), &
!         vcfli(0:3),vcflt(0:3),vcflit(0:3),cdff,psirs,epslimf,epslimt, &
!         cntrpy_f,cntrpy_t,etpfxtyp_f,etpfxtyp_t,entfx_f_w,entfx_t_w, &
!         entfxctf_f,entfxctf_t,rkap, &
!         eig_cutoff_f,eig_cutoff_t,xmom,ymom,zmom, &
!         lngth(3),reyno,pranl,tref,stemp,ainf,rhoinf,tinf,pinf,pi,toler, &
!         qmin,irest,srest,lmax,ncycle,ncycr,iupdt,iprint,ilimf,ilimt, &
!         irsop,cutcirs,nsurface,n_surblo(msurface),i_surbl(msurface,m_surblo), &
!         profpar,nblock_pos,iblock_pos(msurface*m_surblo),wallmask,flowtec, & 
!         surftec,filename,code,cgns_filename 
!    common/calind/krd(3,3),cyc(3,3)
!    common/resnorm/resrms
!    common/uns/omega,omegas,omegatp,betaw,phitp,simtime,dt, &
!         zeit(0:2*mharms),dhb(0:2*mharms,0:2*mharms), & 
!         ahb(0:2*mharms,0:2*mharms),ehb(0:2*mharms,0:2*mharms), &
!         eihb(0:2*mharms,0:2*mharms),ntime,ntime_ramp,nsave, & 
!         wrifl(3),itime,hbmove
!    common/movi/xh,yh,dtheta0,dh0x,dh0y,dtheta(0:2*mharms),phipp, &
!         yhtp,zhtp,theta0tp, &
!         thetatp(0:2*mharms),xrotc,yrotc,zrotc,gtheta,rotmat(2,2),nblade,somega 
!    common/losprec/epsp(12),lx,ly,uprv,cutoff_type,mixpi,lsp_vfr
!    common/logi/aircraft,hawt,persec,vawt,relframe,calfor,exter, &
!         viscous,foturb,kom,kom_bsl,kom_sst,ho_rest,prd_all_lvls,wallwilc, &
!         wallment,lim_prodtke,lim_prodomega,lim_corr, &
!         unsteady,dualt,rgkuns,harbal,found_simtime, &
!         rkimp,dual_prop_start,moving,rotating,tpitch,pitching,plunging, &
!         plupitching,lomach,visprec,debug,wri_res,wri_tec,wldata,blprof, &
!         calvort,calcp,calmut,bln_prol,ramping1,ramping2, &
!         tecbin,dgn_cell,efx_setup_f,efx_setup_t,eig_limit_f,eig_limit_t, &
!         efx_all_f,efx_all_t,eig_lim_all_f,eig_lim_all_t,using_cgns 
!    common/mgpars/nlevel,nl_crs,nstart,npre,npost,ncrs,nl_fmg,mgit
!    common/rkpars/alphas,betas,alphadum,nstagerk,nstage,ar4,br4,ar5
!    common/turpars/epst,prant,tkefar,mutfar,omefar,cmu,bstrkw,gkw, &
!         bkw,sigk,sigw,roughk,prdlim,bstr,sigk1_bsl,sigk1_sst,sigw1, &
!         b1,g1,sigk2,sigw2,b2,g2,a1,ifixq,turcmp
    
    save
    
    contains

   end module common_variables
