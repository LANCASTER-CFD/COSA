!-----------------------------------------------------------------------
!
! Module to hold some of the common blocks from the F77 version of the 
! program. 
!
! Adrian Jackson, EPCC, The University of Edinburgh
! June 2023
!
!-----------------------------------------------------------------------


module cosa_variables

  use cosa_precision

  implicit none
  
  public
  
  logical :: write_cgns = .false.
  logical :: write_plt = .false.
  logical :: write_szplt = .false.
  logical :: write_tec = .false.

  integer, parameter :: mblock=30000
  integer, parameter :: mlevel=3
  integer, parameter :: mcut=200000
  integer, parameter :: mbc=20
  integer, parameter :: mwls=4

  integer(kind=cosa_int) :: npde, irspde, nharms, dim5, dim5h, nblocks, &
       lmet, ncuts, cutdata(21,mcut), npcuts, pcutdata(22,mcut), &
       nbcs(mblock), bcdata(10,mbc,mblock), &
       i_imax(mblock,mlevel), j_jmax(mblock,mlevel), &
       k_kmax(mblock,mlevel), ijk_ijkmax(mblock,mlevel),&
       n_wall(mblock,mlevel), ng_wall(mlevel), &
       off_bct(mblock,mlevel), off_m1(mblock,mlevel), &
       off_0(mblock,mlevel),   off_1d(mblock,mlevel), &
       off_1dw(mblock,mlevel), &
       off_p1(mblock,mlevel),  off_p2(mblock,mlevel), & 
       off_p3(mblock,mlevel), off_1dwg(mblock,mlevel), &
       maxsendid_num, maxrecvid_num, &
       check_maxsendid_num, check_maxrecvid_num, &
       maxpsendid_num, maxprecvid_num,sendrecv_datasize 
  
!  integer :: npde, irspde, nharms, dim5, dim5h, nblocks, lmet, ncuts, & 
!!       cutdata, npcuts, pcutdata, nbcs, bcdata, i_imax, j_jmax, k_kmax, &
!       ijk_ijkmax, n_wall, ng_wall, off_bct, off_m1, off_0, off_p1, &
!       off_p2, off_p3, off_1d, off_1dw, off_1dwg, &
!       maxsendid_num , maxrecvid_num, &
!       sendrecv_datasize, &
!       check_maxrecvid_num, check_maxsendid_num, &
!       maxpsendid_num,  maxprecvid_num
  
  integer(kind=cosa_long) :: &
       p_q(mlevel),        p_q1(mlevel),       p_q2(mlevel), &
       p_qold(mlevel),     p_qp(mlevel),       p_qm(mlevel), & 
       p_dq(mlevel),       p_dqp(mlevel),      p_dqm(mlevel), &
       p_pr(mlevel),       p_flux(mlevel),     p_res(mlevel), &
       p_rhs(mlevel), &
       p_prec(mlevel),     p_preci(mlevel),    p_ipkpi(mlevel), &
       p_dtvol(mlevel),    p_dtvolt(mlevel),   p_rad(mlevel), &
       p_si(mlevel),       p_sj(mlevel),       p_sk(mlevel), &
       p_x(mlevel),        p_y(mlevel),        p_z(mlevel), &
       p_x0(mlevel),       p_y0(mlevel),       p_z0(mlevel), & 
       p_xc(mlevel),       p_yc(mlevel),       p_zc(mlevel), &
       p_dx(mlevel),       p_dy(mlevel),       p_dz(mlevel), &
       p_xdot(mlevel),     p_ydot(mlevel),     p_zdot(mlevel), &
       p_xwall(mlevel),    p_ywall(mlevel),    p_zwall(mlevel), &
       p_xgwall(mlevel),   p_ygwall(mlevel),   p_zgwall(mlevel), &
       p_delpls(mlevel),   p_divvel(mlevel),   p_prd(mlevel), &
       p_mut(mlevel),      p_ttrm(mlevel), & 
       p_vol(mlevel),      p_dist(mlevel), &
       p_betaxi(mlevel),   p_betaeta(mlevel),  p_betazeta(mlevel), &
       p_tbetaxi(mlevel),  p_tbetaeta(mlevel), p_tbetazeta(mlevel), &
       p_fidist(mlevel),   p_fjdist(mlevel),   p_fkdist(mlevel), & 
       p_xideri(mlevel),   p_xiderj(mlevel),   p_xiderk(mlevel), &
       p_etaderi(mlevel),  p_etaderj(mlevel),  p_etaderk(mlevel), &
       p_zetaderi(mlevel), p_zetaderj(mlevel), p_zetaderk(mlevel), &
       p_work1(mlevel),    p_work2(mlevel), &
       p_cutoff_pgr(mlevel),                   p_cutoff_vis(mlevel), &
       p_bctopo(mlevel),   p_cutopo(mlevel),   p_percutopo(mlevel)
  
!  integer(kind=8) :: &
!       p_q     ,p_q1    ,p_q2    ,p_qold  ,p_qp    ,p_qm    ,p_dq    , &
!       p_dqp   ,p_dqm   ,p_flux  ,p_res   ,p_rhs   ,p_pr    ,p_delpls, &
!       p_divvel,p_mut   ,p_prd   ,p_ttrm  ,p_prec  ,p_preci ,p_ipkpi , &
!       p_work1 ,p_work2 ,p_si    ,p_sj    ,p_sk    ,p_x     , &
!       p_y     ,p_z     ,p_x0    ,p_y0    ,p_z0    ,p_xc    ,p_yc    , &
!       p_zc    ,p_dx    ,p_dy    ,p_dz    ,p_xdot  ,p_ydot  ,p_zdot  , &
!       p_xwall ,p_ywall ,p_zwall ,p_xgwall,p_ygwall,p_zgwall,p_vol   , &
!       p_fidist    ,p_fjdist    ,p_fkdist    ,p_dist      ,p_betaxi    , &
!       p_betaeta   ,p_betazeta  ,p_tbetaxi   ,p_tbetaeta  ,p_tbetazeta , &
!       p_xideri    ,p_xiderj    ,p_xiderk    ,p_etaderi   ,p_etaderj   , &
!       p_etaderk   ,p_zetaderi  ,p_zetaderj  ,p_zetaderk  ,p_bctopo    , &
!       p_cutopo    ,p_percutopo ,p_cutoff_pgr,p_cutoff_vis,p_dtvol     , &
!       p_dtvolt    ,p_rad
  
  !mpi  Added for MPI parallelisation
!  integer :: lowernblock, uppernblock, mynblocks, maxblocks, myncuts, &
!       mynpcuts, mycuts, mypcuts, blockassignment, g_i_imax, g_j_jmax, &
!       g_k_kmax, g_ijk_ijkmax, typeint, typedouble, typechar, g_n_wall, &
!       g_nbcs, g_bcdata 
  
  !mpi Added for the MPI parallelisation
  integer(kind=cosa_int) :: lowernblock, uppernblock, mynblocks, maxblocks, &
       myncuts, mynpcuts, mycuts(mcut), mypcuts(mcut), &
       blockassignment(mblock), &
       g_i_imax(mblock,mlevel), & 
       g_j_jmax(mblock,mlevel), &
       g_k_kmax(mblock,mlevel), &
       g_ijk_ijkmax(mblock,mlevel), typeint, &
       typedouble, typechar, &
       g_n_wall(mblock,mlevel), &
       g_nbcs(mblock), &
       g_bcdata(10,mbc,mblock) 
   
  save
  
contains
  
end module cosa_variables
