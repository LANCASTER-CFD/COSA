      module cosa_cgns_utilities

      use cosa_precision

      implicit none
      
      private
      
      public :: get_cgns_data, readcgnsgrid
      
      save
      
      contains
      
      integer function find_family_name(family_names, family_name)
      
      use cosa_precision

      implicit none
      
      character(len=32) :: family_names(:)
      character(len=32) :: family_name
      integer :: number_of_families, i
      
      number_of_families = size(family_names)
      
      find_family_name = -1

      do i = 1, number_of_families
         if(trim(family_names(i)) .eq. trim(family_name))  then
            find_family_name = i
            continue
         end if
      end do

      end function



! Find the number of the surface referenced by this specific family names.
! The CGNS file can have many different family names, but not all of them 
! will link surfaces. Therefore, we have this routine to return a surface 
! id rather than an overall family id.
      integer function find_surface_name(family_names, is_surface, family_name)
      
      use cosa_precision

      implicit none
      
      character(len=32) :: family_names(:)
      logical :: is_surface(:)
      character(len=32) :: family_name
      integer :: number_of_families, i, surface_position
      
      number_of_families = size(family_names)
      
      surface_position = 1
      find_surface_name = -1

      do i = 1, number_of_families
         if(trim(family_names(i)) .eq. trim(family_name))  then
            find_surface_name = surface_position
            continue
         end if
         if(is_surface(i)) then
            surface_position = surface_position + 1
         end if
      end do

      end function

! Find the number of the surface referenced by this specific family names
! for a particular boundary condition type. The CGNS file can have many 
! different family names, but not all of them will link surfaces. Therefore, 
! we have this routine to return a surface id rather than an overall family id.
! We also split surfaces into two groups, 1400 and 1500 based groups (that 
! is the COSA BC types of the groups). For this routine we want a surface number
! that is just for one of the two groups, i.e. just for 1400 or just for 1500, so 
! we can get 1400,1401,1500,1501,1502,etc... This routine works out which type of 
! family this is, and then returns what position it should be.
      integer function find_surface_name_bc_type(family_names, family_types, &
           family_position, family_name)
      
      use cosa_precision

      implicit none
      
      character(len=32) :: family_names(:)
      integer(kind=cosa_int) :: family_types(:), family_position
      character(len=32) :: family_name
      integer :: number_of_families, i, bc_type
      logical :: surface_1500, surface_1400
      integer :: surface_1500_counter, surface_1400_counter      

      surface_1500 = .false.
      surface_1400 = .false.

      surface_1500_counter = 1
      surface_1400_counter = 1

      bc_type = family_types(family_position)

      if(bc_type .eq. 20) then
         write(*,*) 'Problem with getting the surface type. '
         write(*,*) 'Boundary type is ',bc_type
         write(*,*) 'Expecting it to be 21, or 22'
         write(*,*) '20 is too general, you should convert'
         write(*,*) 'to either 21 or 22.'
         stop
      end if
      if(bc_type .eq. 22) then
         surface_1500 = .true.
      else if(bc_type .eq. 21) then
         surface_1400 = .true.
      else
         write(*,*) 'Problem with getting the surface type. '
         write(*,*) 'Boundary type is ',bc_type
         write(*,*) 'Expecting it to be 21, or 22'
         stop
      end if


      number_of_families = size(family_names)
      
      find_surface_name_bc_type = -1

      do i = 1, number_of_families
         if(trim(family_names(i)) .eq. trim(family_name))  then
            if(family_types(i) .eq. 22) then
               find_surface_name_bc_type = surface_1500_counter
            else if(family_types(i) .eq. 21) then
               find_surface_name_bc_type = surface_1400_counter
            else
               write(*,*) 'Problem getting the surface type'
               write(*,*) 'The family name passed did not '
               write(*,*) 'match a suitable CGNS boundary condition, '
               write(*,*) 'i.e., 21 or 22.'
               stop
            end if
            continue
         end if
         if(family_types(i) .eq. 22) then
            surface_1500_counter = surface_1500_counter + 1
         else if(family_types(i) .eq. 21) then
            surface_1400_counter = surface_1400_counter + 1
         end if
      end do

      end function

      integer function get_character_position(string, character_to_find)
      
      use cosa_precision

      implicit none
      
      character(len=*) :: string
      character :: character_to_find
      integer :: string_length, i
      
      string_length = len(string)
      
      get_character_position = -1

      do i = 1, string_length
         if(string(i:i) .eq. character_to_find)  then
            get_character_position = i
            continue
         end if
      end do

      end function

      integer function get_first_number_position(string)

      use cosa_precision

      implicit none

      character(len=*) :: string
      integer :: string_length, i
      character :: temp_char

      string_length = len(string)
      
      get_first_number_position = -1

      do i = 1, string_length
         temp_char = string(i:i)
         if(temp_char .eq. '0' .or. temp_char .eq. '1' &
              .or. temp_char .eq. '2' &
              .or. temp_char .eq. '3' &
              .or. temp_char .eq. '4' &
              .or. temp_char .eq. '5' &
              .or. temp_char .eq. '6' &
              .or. temp_char .eq. '7' &
              .or. temp_char .eq. '8' &
              .or. temp_char .eq. '9')  then
            get_first_number_position = i
            exit
         end if
      end do


      end function

!-----------------------------------------------------------------------
      subroutine get_cgns_data()
!-----------------------------------------------------------------------

      use cgns
      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

#ifdef MPI
      include 'mpif.h'
      integer offset_kind
      parameter(offset_kind=MPI_OFFSET_KIND)
#else
      integer offset_kind
      parameter(offset_kind=8)
#endif

      integer nl,fid,ierr,i
      integer nbases,basenum
      integer(kind=offset_kind) offset
      integer :: maxbcsize,maxgcsize
      integer(kind=cosa_int) imax(mblock),jmax(mblock),kmax(mblock),ijkmax(mblock), &
           nwallv_1(mwls,mblock),nwallv_2(mwls,mblock), &
           nvwls(mblock)
      integer(kind=cosa_int), allocatable :: bcsize(:)
      integer(kind=cosa_int), allocatable :: bcrange(:,:,:),bctype(:,:),gcsize(:)
      integer(cgsize_t), allocatable :: gcrange(:,:,:,:)
      integer, allocatable :: gctransform(:,:,:)
      integer, allocatable :: donornum(:,:)
      character(len=32),allocatable :: blockname(:)
      character(len=32),allocatable :: family_names(:)
      integer(kind=cosa_int), allocatable :: family_types(:)
      logical parallel, amcontrol
      integer rank
      logical, allocatable :: pcut(:,:)
      real(kind=cosa_single), allocatable :: zrotangle(:,:)    
      real(kind=cosa_single), allocatable :: translength(:,:)     
      integer, dimension(25) :: cosabctypes
      integer(kind=cosa_int), allocatable :: blocknumbermapping(:)

      call amcontroller(amcontrol)
      call getmpiid(rank)      
      call setupcosabctypes(cosabctypes,lomach)

      if(amcontrol) then
         write(*,*) 'Now opening CGNS file ',trim(cgns_filename)
      end if            
      call cg_open_f(cgns_filename,CG_MODE_READ,fid,ierr)
      if(ierr.ne.0) then            
         if(amcontrol) then
            write(*,*) 'ERROR: cg_open_f',cgns_filename,' error ',ierr
         end if
         call cg_error_exit_f()
      end if
      call cg_nbases_f(fid,nbases,ierr)
      if(ierr.ne.0) then
         if(amcontrol) then
            write(*,*) 'ERROR: call cg_nbases_f'
         end if
         call cg_error_exit_f()
      end if
      if (nbases.ne.1 .and. amcontrol) then
         write(*, &
              '(''WARNING: Nbases='',i5,'' (and it should be 1)'')') nbases
      end if

      basenum = 1

      call readcgnsblockinformation(nblocks,basenum,fid,imax, &
           jmax,kmax,ijkmax,bcsize,bctype,bcrange,gcsize, &
           gcrange,gctransform,donornum, &
           maxbcsize,maxgcsize,blockname,pcut,zrotangle, &
           translength,mblock,blocknumbermapping)            

      if ((maxbcsize.eq.0).and.(maxgcsize.eq.0)) then
         write(*,*) 'WARNING: Problem with mesh connectivity'
         stop
      else
         
         call read_cgnsfamily(nsurface,msurface,family_names, &
              family_types,basenum,fid)

         call read_cgnsbc(bcsize,bctype,bcrange,maxbcsize,basenum, &
              family_names,family_types,n_surblo,msurface,m_surblo, &
              i_surbl,cosabctypes,nblocks,fid,blocknumbermapping, &
              amcontrol,rank)

         call read_cgnsgc(gcsize,gcrange,gctransform,maxgcsize, &
              blockname,donornum,basenum,nblocks,fid, &
              zrotangle,translength,pcut,blocknumbermapping)
      end if

      call convertcgnsboundaryinformationtocosa(nblocks,bcsize,nbcs)

      call convertcgnsblockinformationtocosa(nblocks,imax,jmax,kmax, &
           ijkmax,bcsize,bctype,bcrange,gcsize,gcrange,gctransform, &
           donornum,nwallv_1,nwallv_2,nvwls,nbcs,bcdata,cyc, &
           debug,kom,kom_bsl,kom_sst,viscous,code,mwls,mblock,mbc, &
           nsurface)

      call setdim(imax,jmax,kmax,ijkmax)
      if (((code.eq.'cosa').and.(kom.or.kom_bsl.or.kom_sst)).or. &
           ((code.eq.'poscosa').and.viscous)) then
         call setdim_vwall(nvwls,nwallv_1,nwallv_2)
         if(amcontrol.and.debug) then
            write(*,*) ' completed setdim_vwall'
         end if
      end if      
      
      call convertcgnsconnectivitydatatocosa(nblocks,gcsize, &
           gcrange,gctransform,maxgcsize,donornum,pcut,zrotangle, &
           translength,mcut,ncuts,npcuts,cutdata,pcutdata,cyc,persec, &
           debug)


      call cg_close_f(fid,ierr)
      if(ierr.ne.0) then
         write(*,*) 'ERROR: cg_close_f'
         call cg_error_exit_f()
      end if               
      
      if(allocated(family_names)) deallocate(family_names)
      if(allocated(family_types)) deallocate(family_types)

      return
      end

!-----------------------------------------------------------------------
      subroutine readcgnsgrid(nl,x,y,z,xyz,dx,dy,dz,dxyz)
!-----------------------------------------------------------------------

      use cgns
      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

#ifdef MPI
      include 'mpif.h'
      integer offset_kind
      parameter(offset_kind=MPI_OFFSET_KIND)
#else
      integer offset_kind
      parameter(offset_kind=8)
#endif

      integer nl,imax,jmax,kmax,iblock,ixyz,jxyz,jdxyz,fid
      integer ierr,basenum,blocknum
      real(kind=cosa_real) :: x(*),y(*),z(*),xyz(*),dx(*),dy(*),dz(*),dxyz(*)
      character*32 basename
      logical parallel, amcontrol

      call usingparallel(parallel)
      call amcontroller(amcontrol)

      if(parallel) then

         if(amcontrol.and.debug) then
            write(*,*) 'calling parallelreadcgnsgrid'
         end if

         call parallelreadcgnsgrid(x,y,z,xyz,dx,dy,dz,dxyz,nl)

      else

         call cg_open_f(cgns_filename,CG_MODE_READ,fid,ierr)
         if(ierr .ne. CG_OK) then
            write(*,*) 'cg_open_f error'
            call cg_error_print_f()
         end if

         basenum = 1

         do iblock = 1,mynblocks
            blocknum = iblock - 1 + lowernblock
            imax  = i_imax(iblock,nl)
            jmax  = j_jmax(iblock,nl)
            kmax  = k_kmax(iblock,nl)
            ixyz  = 1 + off_0 (iblock,nl)*3
            jxyz  = 1 + off_p2(iblock,nl)
            jdxyz = 1 + off_p1(iblock,nl)
            call read_bcgnsgrid(xyz(ixyz),dxyz(ixyz),imax,jmax,kmax, &
                 fid,basenum,blocknum)
            call xyz_offset(x(jxyz),y(jxyz),z(jxyz),xyz(ixyz),dx(jdxyz), &
                 dy(jdxyz),dz(jdxyz),dxyz(ixyz),imax,jmax,kmax)
         end do

         call cg_close_f(fid,ierr)
         if(ierr .ne. CG_OK) then
            write(*,*) 'cg_close_f error'
            call cg_error_print_f()
         end if

      end if

      return
      end

!-----------------------------------------------------------------------
      subroutine read_bcgnsgrid(xyz,dxyz,imax,jmax,kmax,fid,basenum, &
           blocknum)
!-----------------------------------------------------------------------

      use cgns
      use common_variables
      use cosa_precision

      implicit none

      integer imax,jmax,kmax
      integer i,j,k,fid,basenum,blocknum,zonenum,ierr
      integer(cgsize_t) minrange(3),maxrange(3)
      real(kind=cosa_real) &
          xyz (imax,jmax,kmax,3), &
          dxyz(imax,jmax,kmax,3)

      minrange(1) = 1
      minrange(2) = 1
      minrange(3) = 1

      maxrange(1) = imax
      maxrange(2) = jmax
      maxrange(3) = kmax

      zonenum = blocknum

      call cg_coord_read_f(fid,basenum,zonenum,'CoordinateX',RealDouble, &
           minrange,maxrange,xyz(:,:,:,1),ierr)
      if(ierr .ne. CG_OK) then
         write(*,*) 'X cg_coord_read_f error'
         call cg_error_print_f()
      end if
      call cg_coord_read_f(fid,basenum,zonenum,'CoordinateY',RealDouble, &
           minrange,maxrange,xyz(:,:,:,2),ierr)
      if(ierr .ne. CG_OK) then
         write(*,*) 'Y cg_coord_read_f error'
         call cg_error_print_f()
      end if
      call cg_coord_read_f(fid,basenum,zonenum,'CoordinateZ',RealDouble, &
           minrange,maxrange,xyz(:,:,:,3),ierr)
      if(ierr .ne. CG_OK) then
         write(*,*) 'Z cg_coord_read_f error'
         call cg_error_print_f()
      end if

      if (moving.and.pitching) then
         call cg_coord_read_f(fid,basenum,zonenum,'CoordinateDX',RealDouble, &
              minrange,maxrange,dxyz(:,:,:,1),ierr)
         if(ierr .ne. CG_OK) then
            write(*,*) 'moving and pitching X cg_coord_read_f error'
            call cg_error_print_f()
         end if
         call cg_coord_read_f(fid,basenum,zonenum,'CoordinateDY',RealDouble, &
              minrange,maxrange,dxyz(:,:,:,2),ierr)
         if(ierr .ne. CG_OK) then
            write(*,*) 'moving and pitching Y cg_coord_read_f error'
            call cg_error_print_f()
         end if
         call cg_coord_read_f(fid,basenum,zonenum,'CoordinateDZ',RealDouble, &
              minrange,maxrange,dxyz(:,:,:,3),ierr)
         if(ierr .ne. CG_OK) then
            write(*,*) 'moving and pitching Z cg_coord_read_f error'
            call cg_error_print_f()
         end if
      end if

      return
      end

!-----------------------------------------------------------------------
      subroutine setupcosabctypes(cosabctypes,lomach)
!-----------------------------------------------------------------------
! This setups the array cosabctypes to hardcode CGNS boundary condtion
! number mapping to COSA boundary condtion numbers. Each entry in the 
! array is one of the CGNS boundary condition types, and it is 
! initialised with the corresponding COSA boundary condition number.
! Entry 1 (i.e. cosabctypes(1)) is initialised to -1 as it is not used 
! in CGNS boundary condition numbering
!
!
!     CGNS BC Types
!     parameter (BCAxisymmetricWedge     = 2)
!     parameter (BCDegenerateLine        = 3)
!     parameter (BCDegeneratePoint       = 4)
!     parameter (BCDirichlet             = 5)
!     parameter (BCExtrapolate           = 6)
!     parameter (BCFarfield              = 7)
!     parameter (BCGeneral               = 8)
!     parameter (BCInflow                = 9)
!     parameter (BCInflowSubsonic        = 10)
!     parameter (BCInflowSupersonic      = 11)
!     parameter (BCNeumann               = 12)
!     parameter (BCOutflow               = 13)
!     parameter (BCOutflowSubsonic       = 14)
!     parameter (BCOutflowSupersonic     = 15)
!     parameter (BCSymmetryPlane         = 16)
!     parameter (BCSymmetryPolar         = 17)
!     parameter (BCTunnelInflow          = 18)
!     parameter (BCTunnelOutflow         = 19)
!     parameter (BCWall                  = 20)
!     parameter (BCWallInviscid          = 21)
!     parameter (BCWallViscous           = 22)
!     parameter (BCWallViscousHeatFlux   = 23)
!     parameter (BCWallViscousIsothermal = 24)
!     parameter (FamilySpecified         = 25)
!
!     COSA3D BC types:
!     inviscid wall (all vars are extrapolated except normal vel.)    1400 (23)
!     viscous wall                                                    1500 (22)
!     subsonic inflow                                                 21
!     subsonic or supersonic free-stream                              30
!     subsonic or supersonic free-stream,lsp                          31
!     extrapolation                                                    4
!     symmetry wrt xy-plane                                           81
!     symmetry wrt xz-plane                                           82
!     symmetry wrt yz-plane                                           83
!     subsonic outflow                                                91
!AJ 21 -> 14
!AJ 22 -> 15

      use cosa_precision

      implicit none

      integer, dimension(25), intent(inout) :: cosabctypes
      logical, intent(in) :: lomach

      cosabctypes(1) = -1 
      cosabctypes(2) = -1
      cosabctypes(3) = -1
      cosabctypes(4) = -1
      cosabctypes(5) = -1
      cosabctypes(6) = 4
      if(lomach) then
         cosabctypes(7) = 31
      else
         cosabctypes(7) = 30
      end if
      cosabctypes(8) = -1
      cosabctypes(9) = 21
      cosabctypes(10) = 20
      cosabctypes(11) = 21
      cosabctypes(12) = 1
      cosabctypes(13) = 4
      cosabctypes(14) = 90
      cosabctypes(15) = 90
      cosabctypes(16) = 80
      cosabctypes(17) = 80
      cosabctypes(18) = -1
      cosabctypes(19) = 1
      cosabctypes(20) = -1 !AJ 20 would be either 1400 or 1500 so make it an error
      cosabctypes(21) = 1400
      cosabctypes(22) = 1500
      cosabctypes(23) = -1
      cosabctypes(24) = -1
      cosabctypes(25) = -1 
      return
      end

!-----------------------------------------------------------------------
      subroutine readcgnsblockinformation(nblocks,basenum,fid, &
           imax,jmax,kmax,ijkmax, &
           bcsize,bctype,bcrange,gcsize,gcrange,gctransform, &
           donornum,maxbcsize,maxgcsize, &
           blockname,pcut,zrotangle,translength,mblock, &
           block_number_mapping)
!-----------------------------------------------------------------------
      
      use cgns
      use cosa_precision

      implicit none

      integer, intent(inout) :: basenum,fid
      integer(kind=cosa_int), intent(inout) :: nblocks
      integer iargc,iblock,ierr
      integer(kind=cosa_int), intent(inout) :: maxbcsize,maxgcsize
      integer(kind=cosa_int), intent(inout) :: imax(mblock), jmax(mblock)
      integer(kind=cosa_int), intent(inout) :: kmax(mblock), ijkmax(mblock)
      integer(kind=cosa_int), allocatable, intent(inout) :: bcsize(:),gcsize(:)
      integer(kind=cosa_int), allocatable, intent(inout) :: bcrange(:,:,:)
      integer(kind=cosa_int), allocatable, intent(inout) :: gctransform(:,:,:)
      integer(kind=cosa_int), allocatable, intent(inout) :: bctype(:,:),donornum(:,:)
      integer(cgsize_t), allocatable, intent(inout) :: gcrange(:,:,:,:)
      character(len=32),allocatable, intent(inout) :: blockname(:)
      logical, allocatable, intent(inout) :: pcut(:,:)
      real(kind=cosa_single), allocatable, intent(inout) :: zrotangle(:,:)
      real(kind=cosa_single), allocatable, intent(inout) :: translength(:,:)
      integer(kind=cosa_int), intent(in) :: mblock
      integer(kind=cosa_int), allocatable, intent(inout) :: block_number_mapping(:)
      integer :: maximax,maxjmax,maxkmax
      character*72 :: zonename
      integer*8, allocatable :: isize(:,:)
      logical amcontrol

      integer :: temp_nblocks
      integer :: isize_dimension
      integer :: character_position, actual_block_number
      integer :: rank

      call cg_index_dim_f(fid,basenum,1,isize_dimension,ierr)
      if(ierr.ne.0) then
         write(*,*) 'ERROR: call cg_index_dim_f'
         call cg_error_exit_f()
         call abortmpi()
      end if
      allocate(isize(isize_dimension,3))
      isize(:,:) = 1

      call amcontroller(amcontrol)    
      call getmpiid(rank)     

      maximax = 0
      maxjmax = 0
      maxkmax = 0
      maxbcsize = 0
      maxgcsize = 0
      
!     Get number of zones
      call cg_nzones_f(fid,basenum,temp_nblocks,ierr)
      nblocks = temp_nblocks
      if(ierr.ne.0) then
         write(*,*) 'ERROR: call cg_nzones_f'
         call cg_error_exit_f()
         call abortmpi()
      end if
      if(amcontrol) then
         write(*,'('' number of zones/blocks='',i6)') nblocks
      end if
      flush(6)
      if(nblocks .gt. mblock) then
         print*,'ERROR: nblock specified in CGNS is larger than mblock.'
         stop
      end if
      if (nblocks.lt.1) then
        write(*,*) 'Illegal value of nblocks. Aborting!'
        stop
      end if
      flush(6)
      allocate(block_number_mapping(nblocks),stat=ierr)
      if (ierr.ne.0) then
        print*,'ERROR: couldnt allocate memory for ', &
              'actual_block_mapping.'
        stop
      endif      
      allocate(bcsize(nblocks),stat=ierr)
      if (ierr.ne.0) then
         print*,'ERROR: couldnt allocate memory for bcsize.'
         stop
      endif
      allocate(gcsize(nblocks),stat=ierr)
      if (ierr.ne.0) then
         print*,'ERROR: couldnt allocate memory for gcsize.'
         stop
      endif
      allocate(blockname(nblocks),stat=ierr)
      if (ierr.ne.0) then
         print*,'ERROR: couldnt allocate memory for blockname.'
         stop
      endif
      
      
      do iblock=1,nblocks
!     Read general zone information
         call cg_zone_read_f(fid,basenum,iblock,zonename,isize,ierr)
         if(ierr.ne.0) then
            write(*,*) 'ERROR: call cg_zone_read_f'
            call cg_error_exit_f()
         end if
         blockname(iblock) = trim(zonename)
         
         character_position = get_first_number_position(zonename)  
         
         if(character_position .eq. -1) then
            write(*,*) 'Error finding _ in ',zonename
            stop
         end if
         
         read(zonename(character_position:),*,iostat=ierr)  actual_block_number
         
         if(actual_block_number .gt. nblocks) then
            write(*,*) 'Error, block number',actual_block_number, &
                 'exceeds the total number of blocks',nblocks
            stop
         end if

         block_number_mapping(iblock) = actual_block_number

!AJ Add to debug         if(amcontrol) then
!AJ            write(*,*) 'block mapping',iblock,actual_block_number,trim(blockname(iblock))
!AJ         end if

         actual_block_number = iblock

         imax(actual_block_number) = isize(1,1)
         jmax(actual_block_number) = isize(2,1)
         kmax(actual_block_number) = isize(3,1)
         maximax = max0(maximax, imax(actual_block_number))
         maxjmax = max0(maxjmax, jmax(actual_block_number))
         maxkmax = max0(maxkmax, kmax(actual_block_number))        

         ijkmax(actual_block_number) = max0(imax(actual_block_number),jmax(actual_block_number),kmax(actual_block_number))

         if ((imax(actual_block_number).le.1).or.(jmax(actual_block_number).le.1).or. &
              (kmax(actual_block_number).le.1)) then
            write(*,*) 'Inconsistency in problem dimensionality. Aborting! &
     &'
            stop
         end if
         
         call cg_nbocos_f(fid,basenum,actual_block_number,bcsize(actual_block_number),ierr)
         if(ierr.ne.0) then
            write(*,*) 'ERROR: call cg_nbocos_f'
            call cg_error_exit_f()
         end if
         maxbcsize = max0(maxbcsize,bcsize(actual_block_number))        

         call cg_n1to1_f(fid,basenum,actual_block_number,gcsize(actual_block_number),ierr)
         if(ierr.ne.0) then
            write(*,*) 'ERROR: cg_n1to1_f'
            call cg_error_exit_f()
         end if
         maxgcsize = max0(maxgcsize,gcsize(actual_block_number))
!       Initialising maxgcsize to 1 if there are no global connections
!       (i.e. this is single block simulation). Primarily this is to 
!        help with the logic below that allocates grid connectivity.        
         if(maxgcsize .eq. 0) then
            maxgcsize = 1
         end if
      enddo
            
      if ((maxbcsize.eq.0).and.(maxgcsize.eq.0)) then
         write(*,*) 'WARNING: Connectivity NOT found in CGNS file!'
      end if     
      
!     grid coordinates      
!     connectivity
      if ((maxbcsize.gt.0).and.(maxgcsize.gt.0)) then
         allocate(bctype(maxbcsize,nblocks),stat=ierr)
         if (ierr.ne.0) then
            print*,'ERROR: couldnt allocate memory for bctype.'
            stop
         endif
         allocate(bcrange(6,maxbcsize,nblocks),stat=ierr)
         if (ierr.ne.0) then
            print*,'ERROR: couldnt allocate memory for bcrange.'
            stop
         endif
         allocate(gcrange(6,2,maxgcsize,nblocks),stat=ierr)
         if (ierr.ne.0) then
            print*,'ERROR: couldnt allocate memory for gcrange.'
            stop
         endif
         allocate(gctransform(3,maxgcsize,nblocks),stat=ierr)
         if (ierr.ne.0) then
            print*,'ERROR: couldnt allocate memory for gctransform.'
            stop
         endif
         allocate(donornum(maxgcsize,nblocks),stat=ierr)
         if (ierr.ne.0) then
            print*,'ERROR: couldnt allocate memory for donornum.'
            stop
         endif
         allocate(zrotangle(maxgcsize,nblocks),stat=ierr)
         if (ierr.ne.0) then
            print*,'ERROR: couldnt allocate memory for zrotangle.'
            stop
         endif
         zrotangle = 0.d0
         allocate(translength(maxgcsize,nblocks),stat=ierr)
         if (ierr.ne.0) then
            print*,'ERROR: couldnt allocate memory for translength.'
            stop
         endif
         translength = 0.d0
         allocate(pcut(maxgcsize,nblocks),stat=ierr)
         if (ierr.ne.0) then
            print*,'ERROR: couldnt allocate memory for pcut.'
            stop
         endif
         pcut = .false.
      endif
            
      deallocate(isize)

      return 
      end


!-----------------------------------------------------------------------
      subroutine convertcgnsboundaryinformationtocosa(nblocks,bcsize, &
           nbcs)
!-----------------------------------------------------------------------

      use cosa_precision

      implicit none

      integer, intent(in) :: nblocks
      integer(kind=cosa_int), allocatable, intent(inout) :: bcsize(:)
      integer(kind=cosa_int), intent(inout) :: nbcs(nblocks)
      integer :: iblock, ibc

      do iblock = 1,nblocks
         nbcs(iblock) = bcsize(iblock)
      end do


      end
      
!-----------------------------------------------------------------------
      subroutine convertcgnsblockinformationtocosa(nblocks,imax,jmax, &
           kmax,ijkmax,bcsize,bctype,bcrange,gcsize,gcrange,gctransform, &
           donornum,nwallv_1,nwallv_2,nvwls,nbcs,bcdata,cyc, &
           debug,kom,kom_bsl,kom_sst,viscous,code,mwls,mblock,mbc, &
           nsurface)
!-----------------------------------------------------------------------
      
      use cgns
      use cosa_precision

      implicit none

      integer, intent(in) :: nblocks, mwls, mblock, mbc
      integer(kind=cosa_int), intent(in) :: cyc(3,3)
      integer iblock,ierr,ic,iwl,ibc,idir,i,ic2
      integer(kind=cosa_int), intent(inout) :: bcdata(10,mbc,mblock)
      integer(kind=cosa_int), intent(inout) :: imax(nblocks), jmax(nblocks)
      integer(kind=cosa_int), intent(inout) :: kmax(nblocks), ijkmax(nblocks)
      integer(kind=cosa_int), intent(inout) :: nbcs(nblocks)
      integer(kind=cosa_int), allocatable, intent(inout) :: bcsize(:),gcsize(:)
      integer(kind=cosa_int), allocatable, intent(inout) :: bcrange(:,:,:)
      integer(kind=cosa_int), allocatable, intent(inout) :: gctransform(:,:,:)
      integer(kind=cosa_int), allocatable, intent(inout) :: bctype(:,:),donornum(:,:)
      integer(cgsize_t), allocatable, intent(inout) :: gcrange(:,:,:,:) 
      integer(kind=cosa_int), intent(inout) :: nwallv_1(mwls,mblock)
      integer(kind=cosa_int), intent(inout) :: nwallv_2(mwls,mblock)
      integer(kind=cosa_int), intent(inout) :: nvwls(mblock)
      integer(kind=cosa_int), intent(inout) :: nsurface
      logical kom, kom_bsl, kom_sst, viscous
      character*10 code
      logical amcontrol, debug
      integer current_bc_total
      integer :: number_of_1500, number_of_1400

      call amcontroller(amcontrol)   

!---- set array dimensions on all grid levels
      do iblock = 1,mblock
         do iwl = 1,mwls
          nwallv_1(iwl,iblock) = 0
          nwallv_2(iwl,iblock) = 0
        end do
        nvwls(iblock) = 0
      end do

      do iblock = 1,nblocks

        if(amcontrol.and.debug) then
          write(*,*) 'block, imax, jmax, kmax, nbcs(iblock) : ', &
            iblock,imax(iblock),jmax(iblock),kmax(iblock),nbcs(iblock)
        end if

        current_bc_total = 0

        do ibc = 1,nbcs(iblock)
           if(bctype(ibc,iblock) .ne. -1 .and. bctype(ibc,iblock) .ne. 1) then
              current_bc_total = current_bc_total + 1
              bcdata(1,current_bc_total,iblock) = bctype(ibc,iblock)
              idir = 0
              do i = 1,3
                 if (bcrange(i,ibc,iblock).eq.bcrange(i+3,ibc,iblock)) then
                    idir = i
                 end if
              end do
              if (idir.lt.1.or.idir.gt.3) then
                 write(*,*) 'error in reconstructing bc data. Aborting!', &
                      idir,ibc,bcdata(1,ibc,iblock)
                 stop
              end if
              do i = 1,3
                 if (i.ne.idir) then
                    bcrange(i+3,ibc,iblock) = bcrange(i+3,ibc,iblock) - 1
                 end if
              end do
              bcdata(2,current_bc_total,iblock) = idir
!     Start of the istrt(idir) dimension
              bcdata(3,current_bc_total,iblock) = bcrange(idir,ibc,iblock)
!     Start of the first dimension (istrt(1))
              bcdata(4,current_bc_total,iblock) = bcrange(1,ibc,iblock)
!     End of the first dimension (iend(1))
              bcdata(5,current_bc_total,iblock) = bcrange(4,ibc,iblock)
!     Start of the second dimension (istrt(2))
              bcdata(6,current_bc_total,iblock) = bcrange(2,ibc,iblock)
!     End of the second dimension (iend(2))
              bcdata(7,current_bc_total,iblock) = bcrange(5,ibc,iblock)
!     Start of the third dimension (istrt(3))
              bcdata(8,current_bc_total,iblock) = bcrange(3,ibc,iblock)
!     End of the third dimension (iend(3))
              bcdata(9,current_bc_total,iblock) = bcrange(6,ibc,iblock)

              if(amcontrol.and.debug) then
                 write(*,'(A14,2I5,A5,9I5)') 'block, ibc : ',iblock,current_bc_total, &
                      'type ',bcdata(1,current_bc_total,iblock), &
                      bcdata(2,current_bc_total,iblock),bcdata(3,current_bc_total,iblock), &
                      bcdata(4,current_bc_total,iblock),bcdata(5,current_bc_total,iblock), &
                      bcdata(6,current_bc_total,iblock),bcdata(7,current_bc_total,iblock), &
                      bcdata(8,current_bc_total,iblock),bcdata(9,current_bc_total,iblock)
              end if
              
              if (bcdata(1,current_bc_total,iblock)/100.eq.15) then
                 if (((code.eq.'cosa').and.(kom.or.kom_bsl.or.kom_sst)).or. &
                      ((code.eq.'poscosa').and.viscous))             then
                 
                 nvwls(iblock) = nvwls(iblock) + 1
                 nwallv_1(nvwls(iblock),iblock) = 1
                 nwallv_2(nvwls(iblock),iblock) = 1              
                 
                 ic2 = cyc(idir,2)
                 nwallv_1(nvwls(iblock),iblock) = &
                      nwallv_1(nvwls(iblock),iblock) * &
                      (bcrange(ic2+3,ibc,iblock) - bcrange(ic2,ibc,iblock) + &
                      2)
                 
                 ic2 = cyc(idir,3)
                 nwallv_2(nvwls(iblock),iblock) = &
                      nwallv_2(nvwls(iblock),iblock) * &
                      (bcrange(ic2+3,ibc,iblock) - bcrange(ic2,ibc,iblock) + &
                      2)
                 
                 end if
              end if
           end if
        end do
        nbcs(iblock) = current_bc_total
      end do

!     These are to track the number of surfaces
!     If family names are being used this could be be tracked 
!     with nsurface
!     But if we are not, then we need to manually track
      number_of_1500 = 0
      number_of_1400 = 0

      if(nsurface .le. 0) then
        do iblock = 1,nblocks
          do ibc = 1,nbcs(iblock)
            if(bctype(ibc,iblock) .eq. 1500) then
              if(number_of_1500 .lt. 1) then
                 number_of_1500 = 1
              end if 
            else if(bctype(ibc,iblock) .eq. 1501) then
              if(number_of_1500 .lt. 2) then
                 number_of_1500 = 2
              end if
            else if(bctype(ibc,iblock) .eq. 1502) then
              if(number_of_1500 .lt. 3) then
                 number_of_1500 = 3
              end if
            else if(bctype(ibc,iblock) .eq. 1503) then
              if(number_of_1500 .lt. 4) then
                 number_of_1500 = 4
              end if
            else if(bctype(ibc,iblock) .eq. 1400) then
              if(number_of_1400 .lt. 1) then
                 number_of_1400 = 1
              end if
            else if(bctype(ibc,iblock) .eq. 1401) then
              if(number_of_1400 .lt. 2) then
                 number_of_1400 = 2
              end if
            else if(bctype(ibc,iblock) .eq. 1402) then
              if(number_of_1400 .lt. 3) then
                 number_of_1400 = 3
              end if
            else if(bctype(ibc,iblock) .eq. 1403) then
              if(number_of_1400 .lt. 4) then
                 number_of_1400 = 4
              end if
            end if
          end do
         end do
        nsurface = number_of_1500 + number_of_1400
      end if

      if(amcontrol.and.debug) then
         write(*,*) 'Final calculated number of surfaces is ', nsurface
      end if
                  
      return 
      end

!-----------------------------------------------------------------------
      subroutine convertcgnsconnectivitydatatocosa(nblocks,gcsize, &
           gcrange,gctransform,maxgcsize,donornum,pcut,zrotangle, &
           translength,mcut,ncuts,npcuts,cutdata,pcutdata,cyc,persec, &
           debug)
!-----------------------------------------------------------------------

      use cgns
      use cosa_precision

      implicit none

      integer(kind=cosa_int), intent(in) :: nblocks
      integer(kind=cosa_int), intent(in) :: maxgcsize
      integer(kind=cosa_int), intent(inout) :: ncuts, npcuts
      integer(kind=cosa_int), intent(in) :: mcut
      integer(kind=cosa_int), allocatable, intent(in) :: gcsize(:)
      integer(cgsize_t), allocatable, intent(in) :: gcrange(:,:,:,:)
      integer(kind=cosa_int), allocatable, intent(in) :: gctransform(:,:,:)
      integer(kind=cosa_int), allocatable, intent(in) :: donornum(:,:)
      integer(kind=cosa_int), intent(inout) :: cutdata(21,mcut), pcutdata(22,mcut)
      real(kind=cosa_single), allocatable, intent(in) :: zrotangle(:,:)
      real(kind=cosa_single), allocatable, intent(in) :: translength(:,:)
      real(kind=cosa_single) :: r4one, temp1
      logical, allocatable, intent(in) :: pcut(:,:)
      integer(kind=cosa_int), intent(in) :: cyc(3,3)
      integer :: iblock, igc, ic1, ic2, igc2, icut, nper, ierr
      integer :: idir, istrt_idir, current_cut
      integer :: donorigc(maxgcsize,nblocks)
      integer, allocatable :: temporange(:,:,:)
      logical amcontrol, period_msc, persec, debug

      call amcontroller(amcontrol)

      r4one = 1.0

      donorigc = -1

!     write connectivity
      allocate(temporange(6,maxgcsize,nblocks),stat=ierr)
      if (ierr.ne.0) then
        print*,'ERROR: couldnt allocate memory for temporange.'
        stop
      endif

      period_msc = .false.
      if (any(pcut.eqv..true.)) then
        period_msc = .true.
      end if

!     find interface (igc) index for donor block for the case if block 
!     and donor block are the same
      do iblock = 1,nblocks
        do igc = 1,gcsize(iblock)

          if (iblock.eq.donornum(igc,iblock)) then
!             write(*,*) 'gc_transform',iblock,igc,gctransform(:,igc,iblock)
            temporange(1,igc,iblock) = gcrange(1,2,igc,iblock)
            temporange(2,igc,iblock) = gcrange(2,2,igc,iblock)
            temporange(3,igc,iblock) = gcrange(3,2,igc,iblock)
            temporange(4,igc,iblock) = gcrange(4,2,igc,iblock)
            temporange(5,igc,iblock) = gcrange(5,2,igc,iblock)
            temporange(6,igc,iblock) = gcrange(6,2,igc,iblock)
            if (gctransform(1,igc,iblock).lt.0) then
              temporange(1,igc,iblock) = gcrange(4,2,igc,iblock)
              temporange(4,igc,iblock) = gcrange(1,2,igc,iblock)
            else if (gctransform(2,igc,iblock).lt.0) then
              temporange(2,igc,iblock) = gcrange(5,2,igc,iblock)
              temporange(5,igc,iblock) = gcrange(2,2,igc,iblock)
            else if (gctransform(3,igc,iblock).lt.0) then
              temporange(3,igc,iblock) = gcrange(6,2,igc,iblock)
              temporange(6,igc,iblock) = gcrange(3,2,igc,iblock)
            end if 

            do igc2 = 1,gcsize(iblock)
               if (all(temporange(:,igc,iblock).eq. &
                      gcrange(:,1,igc2,iblock))) then
                  donorigc(igc,iblock) = igc2
              end if
            end do

          end if

        end do
      end do

      deallocate(temporange)

      nper = 0

      if (.not.period_msc) then
         
         current_cut =  1
         
         do iblock = 1,nblocks
            do igc = 1,gcsize(iblock)
               if (iblock.gt.donornum(igc,iblock)) cycle

               if ((iblock.eq.donornum(igc,iblock)).and. &
                    (igc.gt.donorigc(igc,iblock))) cycle

               ic1 = 2 * current_cut - 1
               ic2 = 2 * current_cut
               
               cutdata( 1,ic1) = iblock
               cutdata( 4,ic1) = gcrange(1,1,igc,iblock)
               cutdata( 5,ic1) = gcrange(4,1,igc,iblock)
               cutdata( 6,ic1) = gcrange(2,1,igc,iblock)
               cutdata( 7,ic1) = gcrange(5,1,igc,iblock)
               cutdata( 8,ic1) = gcrange(3,1,igc,iblock)
               cutdata( 9,ic1) = gcrange(6,1,igc,iblock)
               
               idir = 0
               if(cutdata( 4,ic1) .eq. cutdata( 5,ic1)) then
                  idir = 1
                  istrt_idir = cutdata( 4,ic1)            
!     Modify start or end appropriately to account for cut input as
!     grid points but cut to be cell-centers. Do this only for the 
!     indices that the cut loops over, not the idir direction.
                  if(cutdata( 7,ic1) .gt. cutdata( 6,ic1)) then
                     cutdata( 7,ic1) = cutdata( 7,ic1) - 1
                  else
                     cutdata( 6,ic1) = cutdata( 6,ic1) - 1                  
                  end if
                  if(cutdata( 9,ic1) .gt. cutdata( 8,ic1)) then
                     cutdata( 9,ic1) = cutdata( 9,ic1) - 1
                  else
                     cutdata( 8,ic1) = cutdata( 8,ic1) - 1
                  end if               
               else if(cutdata( 6,ic1) .eq. cutdata( 7,ic1)) then
                  idir = 2
                  istrt_idir = cutdata( 6,ic1) 
!     Modify start or end appropriately to account for cut input as
!     grid points but cut to be cell-centers. Do this only for the 
!     indices that the cut loops over, not the idir direction.
                  if(cutdata( 5,ic1) .gt. cutdata( 4,ic1)) then
                     cutdata( 5,ic1) = cutdata( 5,ic1) - 1
                  else
                     cutdata( 4,ic1) = cutdata( 4,ic1) - 1                  
                  end if
                  if(cutdata( 9,ic1) .gt. cutdata( 8,ic1)) then
                     cutdata( 9,ic1) = cutdata( 9,ic1) - 1
                  else
                     cutdata( 8,ic1) = cutdata( 8,ic1) - 1
                  end if                              
               else if(cutdata( 8,ic1) .eq. cutdata( 9,ic1)) then
                  idir = 3
                  istrt_idir = cutdata( 8,ic1)
!     Modify start or end appropriately to account for cut input as
!     grid points but cut to be cell-centers. Do this only for the 
!     indices that the cut loops over, not the idir direction.
                  if(cutdata( 5,ic1) .gt. cutdata( 4,ic1)) then
                     cutdata( 5,ic1) = cutdata( 5,ic1) - 1
                  else
                     cutdata( 4,ic1) = cutdata( 4,ic1) - 1
                  end if                              
                  if(cutdata( 7,ic1) .gt. cutdata( 6,ic1)) then
                     cutdata( 7,ic1) = cutdata( 7,ic1) - 1
                  else
                     cutdata( 6,ic1) = cutdata( 6,ic1) - 1                  
                  end if
               else
                  write(*,*) 'error in reading cut data (idir1). &
     &Aborting!'
                  stop
               end if
               
               cutdata( 2,ic1) = idir
               cutdata( 3,ic1) = istrt_idir            
               
               cutdata(10, ic2) = cutdata( 1,ic1)
               cutdata(11, ic2) = cutdata( 2,ic1)
               cutdata(12, ic2) = cutdata( 3,ic1)
               cutdata(13, ic2) = cutdata( 4,ic1)
               cutdata(14, ic2) = cutdata( 5,ic1)
               cutdata(15, ic2) = cutdata( 6,ic1)
               cutdata(16, ic2) = cutdata( 7,ic1)
               cutdata(17, ic2) = cutdata( 8,ic1)
               cutdata(18, ic2) = cutdata( 9,ic1)
               
!     Move on to the second half of the cut data
               
               cutdata(13,ic1) = gcrange(1,2,igc,iblock)
               cutdata(14,ic1) = gcrange(4,2,igc,iblock)
               cutdata(15,ic1) = gcrange(2,2,igc,iblock)
               cutdata(16,ic1) = gcrange(5,2,igc,iblock)
               cutdata(17,ic1) = gcrange(3,2,igc,iblock)
               cutdata(18,ic1) = gcrange(6,2,igc,iblock)
               cutdata(19,ic1) = abs(gctransform(1,igc,iblock))
               cutdata(20,ic1) = abs(gctransform(2,igc,iblock))
               cutdata(21,ic1) = abs(gctransform(3,igc,iblock))
               
               if((cutdata(19,ic1) .lt. 1) .or. &
                    (cutdata(19,ic1) .gt. 3)) then
                  write(*,*) 'error in reading cut data (iord). &
     &Aborting!'
                  stop
               end if
               
               if(cutdata(20,ic1) .lt. 1 .or. &
                    cutdata(20,ic1) .gt. 3) then
                  write(*,*) 'error in reading cut data (iord). &
     &Aborting!'
                  stop
               end if
            
               if(cutdata(21,ic1) .lt. 1 .or. &
                    cutdata(21,ic1) .gt. 3) then
                  write(*,*) 'error in reading cut data (iord). &
     &Aborting!'
                  stop
               end if
               
               idir = 0
               if(cutdata(13,ic1) .eq. cutdata(14,ic1)) then
                  idir = 1
                  istrt_idir = cutdata(13,ic1)            
!     Modify start or end appropriately to account for cut input as
!     grid points but cut to be cell-centers. Do this only for the 
!     indices that the cut loops over, not the idir direction.
                  if(cutdata(16,ic1) .gt. cutdata(15,ic1)) then
                     cutdata(16,ic1) = cutdata(16,ic1) - 1
                  else
                     cutdata(15,ic1) = cutdata(15,ic1) - 1                  
                  end if
                  if(cutdata(18,ic1) .gt. cutdata(17,ic1)) then
                     cutdata(18,ic1) = cutdata(18,ic1) - 1
                  else
                     cutdata(17,ic1) = cutdata(17,ic1) - 1
                  end if               
               else if(cutdata(15,ic1) .eq. cutdata(16,ic1)) then
                  idir = 2
                  istrt_idir = cutdata(15,ic1) 
!     Modify start or end appropriately to account for cut input as
!     grid points but cut to be cell-centers. Do this only for the 
!     indices that the cut loops over, not the idir direction.
                  if(cutdata(14,ic1) .gt. cutdata(13,ic1)) then
                     cutdata(14,ic1) = cutdata(14,ic1) - 1
                  else
                     cutdata(13,ic1) = cutdata(13,ic1) - 1                  
                  end if
                  if(cutdata(18,ic1) .gt. cutdata(17,ic1)) then
                     cutdata(18,ic1) = cutdata(18,ic1) - 1
                  else
                     cutdata(17,ic1) = cutdata(17,ic1) - 1
                  end if
               else if(cutdata(17,ic1) .eq. cutdata(18,ic1)) then
                  idir = 3
                  istrt_idir = cutdata(17,ic1)
!     Modify start or end appropriately to account for cut input as
!     grid points but cut to be cell-centers. Do this only for the 
!     indices that the cut loops over, not the idir direction.
                  if(cutdata(14,ic1) .gt. cutdata(13,ic1)) then
                     cutdata(14,ic1) = cutdata(14,ic1) - 1
                  else
                     cutdata(13,ic1) = cutdata(13,ic1) - 1
                  end if                              
                  if(cutdata(16,ic1) .gt. cutdata(15,ic1)) then
                     cutdata(16,ic1) = cutdata(16,ic1) - 1
                  else
                     cutdata(15,ic1) = cutdata(15,ic1) - 1                  
                  end if
               else
                  write(*,*) 'error in reading cut data (idir2). &
     &Aborting!',cutdata(13,ic1),cutdata(14,ic1),cutdata(15,ic1), &
                       cutdata(16,ic1),cutdata(17,ic1),cutdata(18,ic1)
                  stop
               end if
               
               cutdata(10,ic1) = donornum(igc,iblock)
               cutdata(11,ic1) = idir
               cutdata(12,ic1) = istrt_idir
                              
               cutdata( 1,ic2) = cutdata(10,ic1)
               cutdata( 2,ic2) = cutdata(11,ic1)
               cutdata( 3,ic2) = cutdata(12,ic1)
               cutdata( 4,ic2) = cutdata(13,ic1)
               cutdata( 5,ic2) = cutdata(14,ic1)
               cutdata( 6,ic2) = cutdata(15,ic1)
               cutdata( 7,ic2) = cutdata(16,ic1)
               cutdata( 8,ic2) = cutdata(17,ic1)
               cutdata( 9,ic2) = cutdata(18,ic1)
               cutdata(19,ic2) = cyc(cutdata(19,ic1),1)
               cutdata(20,ic2) = cyc(cutdata(20,ic1),1)
               cutdata(21,ic2) = cyc(cutdata(21,ic1),1)
               
               if ((cutdata(19,ic1).eq.2) .and. &
                    (cutdata(20,ic1).eq.3)) then
                  cutdata(19,ic2) = cyc(cutdata(19,ic1),2)
                  cutdata(20,ic2) = cyc(cutdata(20,ic1),2)
                  cutdata(21,ic2) = cyc(cutdata(21,ic1),2)
               else if ((cutdata(19,ic1).eq.3) .and. &
                       (cutdata(20,ic1).eq.1)) then
                  cutdata(19,ic2) = cyc(cutdata(19,ic1),3)
                  cutdata(20,ic2) = cyc(cutdata(20,ic1),3)
                  cutdata(21,ic2) = cyc(cutdata(21,ic1),3)
               else
                  cutdata(19,ic2) = cyc(cutdata(19,ic1),1)
                  cutdata(20,ic2) = cyc(cutdata(20,ic1),1)
                  cutdata(21,ic2) = cyc(cutdata(21,ic1),1)
               end if

               current_cut = current_cut + 1
               
            end do
         end do
         ncuts = sum(gcsize)
         
      else if (period_msc) then

         nper = 0
         
         do iblock = 1,nblocks
            do igc = 1,gcsize(iblock)
               if(pcut(igc,iblock)) then
                  nper = nper + 1
               end if
            end do
         end do

         npcuts = nper/2
         if (npcuts.le.0 .and. nblocks .ne. 1) then
            write(*,*) 'Missing periodic boundary definition. Aborting!'
            stop
         end if

         current_cut = 1
         
         do iblock = 1,nblocks
            do igc = 1,gcsize(iblock)

              if(pcut(igc,iblock)) cycle

              if ((abs(zrotangle(igc,iblock)).gt.1d-7).or. &
                   (abs(translength(igc,iblock)).gt.1d-7)) cycle
                 
              if (iblock.gt.donornum(igc,iblock)) cycle
                 
              if ((iblock.eq.donornum(igc,iblock)).and. &
                   (igc.gt.donorigc(igc,iblock))) cycle

               ic1 = 2 * current_cut - 1
               ic2 = 2 * current_cut
               
               cutdata( 1,ic1) = iblock
               cutdata( 4,ic1) = gcrange(1,1,igc,iblock)
               cutdata( 5,ic1) = gcrange(4,1,igc,iblock)
               cutdata( 6,ic1) = gcrange(2,1,igc,iblock)
               cutdata( 7,ic1) = gcrange(5,1,igc,iblock)
               cutdata( 8,ic1) = gcrange(3,1,igc,iblock)
               cutdata( 9,ic1) = gcrange(6,1,igc,iblock)
               
               idir = 0
               if(cutdata( 4,ic1) .eq. cutdata( 5,ic1)) then
                  idir = 1
                  istrt_idir = cutdata( 4,ic1)            
!     Modify start or end appropriately to account for cut input as
!     grid points but cut to be cell-centers. Do this only for the 
!     indices that the cut loops over, not the idir direction.
                  if(cutdata( 7,ic1) .gt. cutdata( 6,ic1)) then
                     cutdata( 7,ic1) = cutdata( 7,ic1) - 1
                  else
                     cutdata( 6,ic1) = cutdata( 6,ic1) - 1                  
                  end if
                  if(cutdata( 9,ic1) .gt. cutdata( 8,ic1)) then
                     cutdata( 9,ic1) = cutdata( 9,ic1) - 1
                  else
                     cutdata( 8,ic1) = cutdata( 8,ic1) - 1
                  end if               
               else if(cutdata( 6,ic1) .eq. cutdata( 7,ic1)) then
                  idir = 2
                  istrt_idir = cutdata( 6,ic1) 
!     Modify start or end appropriately to account for cut input as
!     grid points but cut to be cell-centers. Do this only for the 
!     indices that the cut loops over, not the idir direction.
                  if(cutdata( 5,ic1) .gt. cutdata( 4,ic1)) then
                     cutdata( 5,ic1) = cutdata( 5,ic1) - 1
                  else
                     cutdata( 4,ic1) = cutdata( 4,ic1) - 1                  
                  end if
                  if(cutdata( 9,ic1) .gt. cutdata( 8,ic1)) then
                     cutdata( 9,ic1) = cutdata( 9,ic1) - 1
                  else
                     cutdata( 8,ic1) = cutdata( 8,ic1) - 1
                  end if                              
               else if(cutdata( 8,ic1) .eq. cutdata( 9,ic1)) then
                  idir = 3
                  istrt_idir = cutdata( 8,ic1)
!     Modify start or end appropriately to account for cut input as
!     grid points but cut to be cell-centers. Do this only for the 
!     indices that the cut loops over, not the idir direction.
                  if(cutdata( 5,ic1) .gt. cutdata( 4,ic1)) then
                     cutdata( 5,ic1) = cutdata( 5,ic1) - 1
                  else
                     cutdata( 4,ic1) = cutdata( 4,ic1) - 1
                  end if                              
                  if(cutdata( 7,ic1) .gt. cutdata( 6,ic1)) then
                     cutdata( 7,ic1) = cutdata( 7,ic1) - 1
                  else
                     cutdata( 6,ic1) = cutdata( 6,ic1) - 1                  
                  end if
               else
                  write(*,*) 'error in reading cut data (idir1). &
     &Aborting!'
                  stop
               end if
               
               cutdata( 2,ic1) = idir
               cutdata( 3,ic1) = istrt_idir            
               
               cutdata(10, ic2) = cutdata( 1,ic1)
               cutdata(11, ic2) = cutdata( 2,ic1)
               cutdata(12, ic2) = cutdata( 3,ic1)
               cutdata(13, ic2) = cutdata( 4,ic1)
               cutdata(14, ic2) = cutdata( 5,ic1)
               cutdata(15, ic2) = cutdata( 6,ic1)
               cutdata(16, ic2) = cutdata( 7,ic1)
               cutdata(17, ic2) = cutdata( 8,ic1)
               cutdata(18, ic2) = cutdata( 9,ic1)
               
!     Move on to the second half of the cut data
               
               cutdata(13,ic1) = gcrange(1,2,igc,iblock)
               cutdata(14,ic1) = gcrange(4,2,igc,iblock)
               cutdata(15,ic1) = gcrange(2,2,igc,iblock)
               cutdata(16,ic1) = gcrange(5,2,igc,iblock)
               cutdata(17,ic1) = gcrange(3,2,igc,iblock)
               cutdata(18,ic1) = gcrange(6,2,igc,iblock)
               cutdata(19,ic1) = abs(gctransform(1,igc,iblock))
               cutdata(20,ic1) = abs(gctransform(2,igc,iblock))
               cutdata(21,ic1) = abs(gctransform(3,igc,iblock))
               
               if((cutdata(19,ic1) .lt. 1) .or. &
                    (cutdata(19,ic1) .gt. 3)) then
                  write(*,*) 'error in reading cut data (iord). &
     &Aborting!'
                  stop
               end if
               
               if(cutdata(20,ic1) .lt. 1 .or. &
                    cutdata(20,ic1) .gt. 3) then
                  write(*,*) 'error in reading cut data (iord). &
     &Aborting!'
                  stop
               end if
            
               if(cutdata(21,ic1) .lt. 1 .or. &
                    cutdata(21,ic1) .gt. 3) then
                  write(*,*) 'error in reading cut data (iord). &
     &Aborting!'
                  stop
               end if
               
               idir = 0
               if(cutdata(13,ic1) .eq. cutdata(14,ic1)) then
                  idir = 1
                  istrt_idir = cutdata(13,ic1)            
!     Modify start or end appropriately to account for cut input as
!     grid points but cut to be cell-centers. Do this only for the 
!     indices that the cut loops over, not the idir direction.
                  if(cutdata(16,ic1) .gt. cutdata(15,ic1)) then
                     cutdata(16,ic1) = cutdata(16,ic1) - 1
                  else
                     cutdata(15,ic1) = cutdata(15,ic1) - 1                  
                  end if
                  if(cutdata(18,ic1) .gt. cutdata(17,ic1)) then
                     cutdata(18,ic1) = cutdata(18,ic1) - 1
                  else
                     cutdata(17,ic1) = cutdata(17,ic1) - 1
                  end if               
               else if(cutdata(15,ic1) .eq. cutdata(16,ic1)) then
                  idir = 2
                  istrt_idir = cutdata(15,ic1) 
!     Modify start or end appropriately to account for cut input as
!     grid points but cut to be cell-centers. Do this only for the 
!     indices that the cut loops over, not the idir direction.
                  if(cutdata(14,ic1) .gt. cutdata(13,ic1)) then
                     cutdata(14,ic1) = cutdata(14,ic1) - 1
                  else
                     cutdata(13,ic1) = cutdata(13,ic1) - 1                  
                  end if
                  if(cutdata(18,ic1) .gt. cutdata(17,ic1)) then
                     cutdata(18,ic1) = cutdata(18,ic1) - 1
                  else
                     cutdata(17,ic1) = cutdata(17,ic1) - 1
                  end if
               else if(cutdata(17,ic1) .eq. cutdata(18,ic1)) then
                  idir = 3
                  istrt_idir = cutdata(17,ic1)
!     Modify start or end appropriately to account for cut input as
!     grid points but cut to be cell-centers. Do this only for the 
!     indices that the cut loops over, not the idir direction.
                  if(cutdata(14,ic1) .gt. cutdata(13,ic1)) then
                     cutdata(14,ic1) = cutdata(14,ic1) - 1
                  else
                     cutdata(13,ic1) = cutdata(13,ic1) - 1
                  end if                              
                  if(cutdata(16,ic1) .gt. cutdata(15,ic1)) then
                     cutdata(16,ic1) = cutdata(16,ic1) - 1
                  else
                     cutdata(15,ic1) = cutdata(15,ic1) - 1                  
                  end if
               else
                  write(*,*) 'error in reading cut data (idir2). &
     &Aborting!'
                  stop
               end if
               
               cutdata(10,ic1) = donornum(igc,iblock)
               cutdata(11,ic1) = idir
               cutdata(12,ic1) = istrt_idir
               
               
               cutdata( 1,ic2) = cutdata(10,ic1)
               cutdata( 2,ic2) = cutdata(11,ic1)
               cutdata( 3,ic2) = cutdata(12,ic1)
               cutdata( 4,ic2) = cutdata(13,ic1)
               cutdata( 5,ic2) = cutdata(14,ic1)
               cutdata( 6,ic2) = cutdata(15,ic1)
               cutdata( 7,ic2) = cutdata(16,ic1)
               cutdata( 8,ic2) = cutdata(17,ic1)
               cutdata( 9,ic2) = cutdata(18,ic1)
               cutdata(19,ic2) = cyc(cutdata(19,ic1),1)
               cutdata(20,ic2) = cyc(cutdata(20,ic1),1)
               cutdata(21,ic2) = cyc(cutdata(21,ic1),1)
               
               if ((cutdata(19,ic1).eq.2) .and. &
                    (cutdata(20,ic1).eq.3)) then
                  cutdata(19,ic2) = cyc(cutdata(19,ic1),2)
                  cutdata(20,ic2) = cyc(cutdata(20,ic1),2)
                  cutdata(21,ic2) = cyc(cutdata(21,ic1),2)
               else if ((cutdata(19,ic1).eq.3) .and. &
                       (cutdata(20,ic1).eq.1)) then
                  cutdata(19,ic2) = cyc(cutdata(19,ic1),3)
                  cutdata(20,ic2) = cyc(cutdata(20,ic1),3)
                  cutdata(21,ic2) = cyc(cutdata(21,ic1),3)
               else
                  cutdata(19,ic2) = cyc(cutdata(19,ic1),1)
                  cutdata(20,ic2) = cyc(cutdata(20,ic1),1)
                  cutdata(21,ic2) = cyc(cutdata(21,ic1),1)
               end if

               current_cut = current_cut + 1
               
            end do
         end do

         ncuts =  sum(gcsize)/2 - nper/2
         ncuts = ncuts * 2         

         if(persec) then           
            
            current_cut = 1
            
            do iblock = 1,nblocks
               do igc = 1,gcsize(iblock)

                  if (pcut(igc,iblock)) then
                     
                     if (iblock.gt.donornum(igc,iblock)) cycle
                     
                     if ((iblock.eq.donornum(igc,iblock)).and. &
                          (igc.gt.donorigc(igc,iblock))) cycle
                     
                     ic1 = 2 * current_cut - 1
                     ic2 = 2 * current_cut                   
                     
                     pcutdata( 1,ic1) = iblock
                     pcutdata( 4,ic1) = gcrange(1,1,igc,iblock)
                     pcutdata( 5,ic1) = gcrange(4,1,igc,iblock)
                     pcutdata( 6,ic1) = gcrange(2,1,igc,iblock)
                     pcutdata( 7,ic1) = gcrange(5,1,igc,iblock)
                     pcutdata( 8,ic1) = gcrange(3,1,igc,iblock)
                     pcutdata( 9,ic1) = gcrange(6,1,igc,iblock)
                     
                     idir = 0
                     if(pcutdata( 4,ic1) .eq. pcutdata( 5,ic1)) then
                        idir = 1
                        istrt_idir = pcutdata( 4,ic1)            
!     Modify start or end appropriately to account for cut input as
!     grid points but cut to be cell-centers. Do this only for the 
!     indices that the cut loops over, not the idir direction.
                        if(pcutdata( 7,ic1) .gt. pcutdata( 6,ic1)) then
                           pcutdata( 7,ic1) = pcutdata( 7,ic1) - 1
                        else
                           pcutdata( 6,ic1) = pcutdata( 6,ic1) - 1                  
                        end if
                        if(pcutdata( 9,ic1) .gt. pcutdata( 8,ic1)) then
                           pcutdata( 9,ic1) = pcutdata( 9,ic1) - 1
                        else
                           pcutdata( 8,ic1) = pcutdata( 8,ic1) - 1
                        end if               
                     else if(pcutdata( 6,ic1) .eq. pcutdata( 7,ic1)) then
                        idir = 2
                        istrt_idir = pcutdata( 6,ic1) 
!     Modify start or end appropriately to account for cut input as
!     grid points but cut to be cell-centers. Do this only for the 
!     indices that the cut loops over, not the idir direction.
                        if(pcutdata( 5,ic1) .gt. pcutdata( 4,ic1)) then
                           pcutdata( 5,ic1) = pcutdata( 5,ic1) - 1
                        else
                           pcutdata( 4,ic1) = pcutdata( 4,ic1) - 1                  
                        end if
                        if(pcutdata( 9,ic1) .gt. pcutdata( 8,ic1)) then
                           pcutdata( 9,ic1) = pcutdata( 9,ic1) - 1
                        else
                           pcutdata( 8,ic1) = pcutdata( 8,ic1) - 1
                        end if                              
                     else if(pcutdata( 8,ic1) .eq. pcutdata( 9,ic1)) then
                        idir = 3
                        istrt_idir = pcutdata( 8,ic1)
!     Modify start or end appropriately to account for cut input as
!     grid points but cut to be cell-centers. Do this only for the 
!     indices that the cut loops over, not the idir direction.
                        if(pcutdata( 5,ic1) .gt. pcutdata( 4,ic1)) then
                           pcutdata( 5,ic1) = pcutdata( 5,ic1) - 1
                        else
                           pcutdata( 4,ic1) = pcutdata( 4,ic1) - 1
                        end if                              
                        if(pcutdata( 7,ic1) .gt. pcutdata( 6,ic1)) then
                           pcutdata( 7,ic1) = pcutdata( 7,ic1) - 1
                        else
                           pcutdata( 6,ic1) = pcutdata( 6,ic1) - 1                  
                        end if
                     else
                        write(*,*) 'error in reading pcut data (idir1). &
     &Aborting!'
                        stop
                     end if
                     
                     pcutdata( 2,ic1) = idir
                     pcutdata( 3,ic1) = istrt_idir            
                     
                     pcutdata(10, ic2) = pcutdata( 1,ic1)
                     pcutdata(11, ic2) = pcutdata( 2,ic1)
                     pcutdata(12, ic2) = pcutdata( 3,ic1)
                     pcutdata(13, ic2) = pcutdata( 4,ic1)
                     pcutdata(14, ic2) = pcutdata( 5,ic1)
                     pcutdata(15, ic2) = pcutdata( 6,ic1)
                     pcutdata(16, ic2) = pcutdata( 7,ic1)
                     pcutdata(17, ic2) = pcutdata( 8,ic1)
                     pcutdata(18, ic2) = pcutdata( 9,ic1)
                     
!     Move on to the second half of the pcut data
                     
                     pcutdata(13,ic1) = gcrange(1,2,igc,iblock)
                     pcutdata(14,ic1) = gcrange(4,2,igc,iblock)
                     pcutdata(15,ic1) = gcrange(2,2,igc,iblock)
                     pcutdata(16,ic1) = gcrange(5,2,igc,iblock)
                     pcutdata(17,ic1) = gcrange(3,2,igc,iblock)
                     pcutdata(18,ic1) = gcrange(6,2,igc,iblock)
                     pcutdata(19,ic1) = abs(gctransform(1,igc,iblock))
                     pcutdata(20,ic1) = abs(gctransform(2,igc,iblock))
                     pcutdata(21,ic1) = abs(gctransform(3,igc,iblock))                  
                     
                     if((pcutdata(19,ic1) .lt. 1) .or. &
                          (pcutdata(19,ic1) .gt. 3)) then
                        write(*,*) 'error in reading pcut data (iord). &
     &Aborting!'
                        stop
                     end if
                     
                     if(pcutdata(20,ic1) .lt. 1 .or. &
                          pcutdata(20,ic1) .gt. 3) then
                        write(*,*) 'error in reading pcut data (iord). &
     &Aborting!'
                        stop
                     end if
                     
                     if(pcutdata(21,ic1) .lt. 1 .or. &
                          pcutdata(21,ic1) .gt. 3) then
                        write(*,*) 'error in reading pcut data (iord). &
     &Aborting!'
                        stop
                     end if
                     
                     idir = 0
                     if(pcutdata(13,ic1) .eq. pcutdata(14,ic1)) then
                        idir = 1
                        istrt_idir = pcutdata(13,ic1)            
!     Modify start or end appropriately to account for cut input as
!     grid points but cut to be cell-centers. Do this only for the 
!     indices that the cut loops over, not the idir direction.
                        if(pcutdata(16,ic1) .gt. pcutdata(15,ic1)) then
                           pcutdata(16,ic1) = pcutdata(16,ic1) - 1
                        else
                           pcutdata(15,ic1) = pcutdata(15,ic1) - 1                  
                        end if
                        if(pcutdata(18,ic1) .gt. pcutdata(17,ic1)) then
                           pcutdata(18,ic1) = pcutdata(18,ic1) - 1
                        else
                           pcutdata(17,ic1) = pcutdata(17,ic1) - 1
                        end if               
                     else if(pcutdata(15,ic1) .eq. pcutdata(16,ic1)) then
                        idir = 2
                        istrt_idir = pcutdata(15,ic1) 
!     Modify start or end appropriately to account for cut input as
!     grid points but cut to be cell-centers. Do this only for the 
!     indices that the cut loops over, not the idir direction.
                        if(pcutdata(14,ic1) .gt. pcutdata(13,ic1)) then
                           pcutdata(14,ic1) = pcutdata(14,ic1) - 1
                        else
                           pcutdata(13,ic1) = pcutdata(13,ic1) - 1                  
                        end if
                        if(pcutdata(18,ic1) .gt. pcutdata(17,ic1)) then
                           pcutdata(18,ic1) = pcutdata(18,ic1) - 1
                        else
                           pcutdata(17,ic1) = pcutdata(17,ic1) - 1
                        end if
                     else if(pcutdata(17,ic1) .eq. pcutdata(18,ic1)) then
                        idir = 3
                        istrt_idir = pcutdata(17,ic1)
!     Modify start or end appropriately to account for cut input as
!     grid points but cut to be cell-centers. Do this only for the 
!     indices that the cut loops over, not the idir direction.
                        if(pcutdata(14,ic1) .gt. pcutdata(13,ic1)) then
                           pcutdata(14,ic1) = pcutdata(14,ic1) - 1
                        else
                           pcutdata(13,ic1) = pcutdata(13,ic1) - 1
                        end if                              
                        if(pcutdata(16,ic1) .gt. pcutdata(15,ic1)) then
                           pcutdata(16,ic1) = pcutdata(16,ic1) - 1
                        else
                           pcutdata(15,ic1) = pcutdata(15,ic1) - 1                  
                        end if
                     else
                        write(*,*) 'error in reading pcut data (idir2). &
     &Aborting!'
                        stop
                     end if
                     
                     pcutdata(10,ic1) = donornum(igc,iblock)
                     pcutdata(11,ic1) = idir
                     pcutdata(12,ic1) = istrt_idir
                     
                     pcutdata( 1,ic2) = pcutdata(10,ic1)
                     pcutdata( 2,ic2) = pcutdata(11,ic1)
                     pcutdata( 3,ic2) = pcutdata(12,ic1)
                     pcutdata( 4,ic2) = pcutdata(13,ic1)
                     pcutdata( 5,ic2) = pcutdata(14,ic1)
                     pcutdata( 6,ic2) = pcutdata(15,ic1)
                     pcutdata( 7,ic2) = pcutdata(16,ic1)
                     pcutdata( 8,ic2) = pcutdata(17,ic1)
                     pcutdata( 9,ic2) = pcutdata(18,ic1)
                     pcutdata(19,ic2) = cyc(pcutdata(19,ic1),1)
                     pcutdata(20,ic2) = cyc(pcutdata(20,ic1),1)
                     pcutdata(21,ic2) = cyc(pcutdata(21,ic1),1)
                     
                     if ((pcutdata(19,ic1).eq.2) .and. &
                          (pcutdata(20,ic1).eq.3)) then
                        pcutdata(19,ic2) = cyc(pcutdata(19,ic1),2)
                        pcutdata(20,ic2) = cyc(pcutdata(20,ic1),2)
                        pcutdata(21,ic2) = cyc(pcutdata(21,ic1),2)
                     else if ((pcutdata(19,ic1).eq.3) .and. &
                             (pcutdata(20,ic1).eq.1)) then
                        pcutdata(19,ic2) = cyc(pcutdata(19,ic1),3)
                        pcutdata(20,ic2) = cyc(pcutdata(20,ic1),3)
                        pcutdata(21,ic2) = cyc(pcutdata(21,ic1),3)
                     else
                        pcutdata(19,ic2) = cyc(pcutdata(19,ic1),1)
                        pcutdata(20,ic2) = cyc(pcutdata(20,ic1),1)
                        pcutdata(21,ic2) = cyc(pcutdata(21,ic1),1)
                     end if
                     
                     temp1 = 0.0
                     if (abs(zrotangle(igc,iblock)).gt.1d-7) then
                        temp1 = zrotangle(igc,iblock)
                     else if(abs(translength(igc,iblock)).gt.1d-7) then
                        temp1 = translength(igc,iblock)
                     else
                        write(*,*) 'Error writing the periodic cuts'
                        write(*,*) 'Neither rotational or translation'
                        stop
                     end if

                     pcutdata(22,ic2) = -int(sign(r4one,temp1))
                     pcutdata(22,ic1) = int(sign(r4one,temp1))
                     
                     current_cut = current_cut + 1
                     
                  end if
               end do
            end do
            
!AJ deal with the fact that some cuts are skipped in the code above
            npcuts  = (current_cut - 1) * 2
            
         end if            
      
      end if

      if(amcontrol.and.debug) then
         write(*,*) 'ncuts = ',ncuts/2         
         if(persec) then
            write(*,*) 'npcuts = ',npcuts/2
            write(*,*) 'current_cut = ',(current_cut-1)/2
         end if
      end if

      if(amcontrol.and.debug) then
         write(*,*) 'cuts'
         do icut = 1,ncuts
            write(*,'(A4,23I5)') 'icut',icut,cutdata(:,icut)
         end do
         write(*,*) ' '
      end if
      

      if(amcontrol.and.debug) then
         if(persec) then
            write(*,*) 'pcuts'
            do icut = 1,npcuts
               write(*,'(A4,23I5)') 'icut',icut,pcutdata(:,icut)
            end do
            write(*,*) ' '
         end if
      end if


      return
      end

!-----------------------------------------------------------------------
      subroutine read_cgnsgrid(imax,jmax,kmax,maximax,maxjmax,maxkmax, &
           x,y,z,basenum,nblocks,fid)
!-----------------------------------------------------------------------
      use cgns 
      use cosa_precision

      implicit none
      
      integer(kind=cosa_int) imax(nblocks),jmax(nblocks),kmax(nblocks)
      integer(kind=cosa_int) fid,nblocks,iblock
      integer(cgsize_t) minrange(3),maxrange(3)
      integer(kind=cosa_int) basenum,zonenum,ierr
      integer(kind=cosa_int) maximax,maxjmax,maxkmax
      
      real(kind=cosa_real) x(maximax,maxjmax,maxkmax,nblocks), &
           y(maximax,maxjmax,maxkmax,nblocks), &
           z(maximax,maxjmax,maxkmax,nblocks)
      real(kind=cosa_real),allocatable :: tempdata(:,:,:)
      
      
      do iblock = 1,nblocks
         
         minrange(1) = 1
         minrange(2) = 1
         minrange(3) = 1
         maxrange(1) = imax(iblock)
         maxrange(2) = jmax(iblock)
         maxrange(3) = kmax(iblock)
         
         zonenum = iblock
         
         allocate(tempdata(imax(iblock),jmax(iblock),kmax(iblock)), &
              stat=ierr)
         if (ierr.ne.0) then
            print*,'ERROR: couldnt allocate memory for tempdata.'
            stop
         endif
         
         call cg_coord_read_f(fid,basenum,zonenum,'CoordinateX', &
              RealDouble,minrange,maxrange,tempdata(:,:,:),ierr)
         if(ierr.ne.0) then
            write(*,*) 'ERROR: call cg_coord_read_f; CoordinateX'
            call cg_error_exit_f()
         end if
         
         x(1:imax(iblock),1:jmax(iblock),1:kmax(iblock),iblock) = &
              tempdata
         
         call cg_coord_read_f(fid,basenum,zonenum,'CoordinateY', &
              RealDouble,minrange,maxrange,tempdata(:,:,:),ierr)
         if(ierr.ne.0) then
            write(*,*) 'ERROR: call cg_coord_read_f; CoordinateY'
            call cg_error_exit_f()
         end if
         
         y(1:imax(iblock),1:jmax(iblock),1:kmax(iblock),iblock) = &
              tempdata
         
         call cg_coord_read_f(fid,basenum,zonenum,'CoordinateZ', &
              RealDouble,minrange,maxrange,tempdata(:,:,:),ierr)
         if(ierr.ne.0) then
            write(*,*) 'ERROR: call cg_coord_read_f; CoordinateZ'
            call cg_error_exit_f()
         end if
         
         z(1:imax(iblock),1:jmax(iblock),1:kmax(iblock),iblock) = &
              tempdata
         
         deallocate(tempdata)
         
      end do
      
      return
      
      end

!-----------------------------------------------------------------------
      subroutine read_cgnsfamily(nsurface,msurface,family_names, &
           family_types,basenum,fid)
!-----------------------------------------------------------------------
      
      use cgns
      use cosa_precision

      implicit none
      
!     bcrange is 6 as we are considering 3d data so there are 
!     lower and upper bounds for each of the 3 dimensions
      character(len=32),allocatable, intent(inout) :: family_names(:)
      integer(kind=cosa_int), allocatable, intent(inout) :: family_types(:)
      integer(kind=cosa_int) nsurface, msurface, fid, basenum, families, tempbcnum
      integer(kind=cosa_int) tempgcnum, ierr, tempsurf
      character*32 familyname 

      logical amcontrol

      call amcontroller(amcontrol)   

      call cg_nfamilies_f(fid, basenum, nsurface, ierr)
      if (ierr .ne. 0) then
         write(*,*) 'Problem reading cg_nfamilies_f',ierr
         call cg_error_exit_f()
      end if
               
      allocate(family_names(nsurface))
      allocate(family_types(nsurface))
   
      do families=1,nsurface
         call cg_family_read_f(fid, basenum, families, familyname, &
              tempbcnum, tempgcnum, ierr)
         if (ierr .ne. 0) then
            write(*,*) 'Problem reading cg_family_read_f',ierr,families
            stop                 
         end if                  
         
         family_names(families) = trim(familyname)

         call cg_fambc_read_f(fid, basenum, families, 1, familyname, &
              family_types(families), ierr)

      end do


      if(amcontrol) then
         write(*,*) 'Number of families ',nsurface
         do families=1,nsurface
            write(*,*) 'Family ',families,' is ',family_names(families), &
                 ' with type ',family_types(families)
         end do
      end if

      tempsurf = 0

      do families=1,nsurface
         if(any(family_types(families).eq.[21,22,25])) then
            tempsurf = tempsurf + 1
          end if
      end do

      nsurface = tempsurf
      
      if (nsurface.gt.msurface) then
         write(*,*) 'Exceeded max. number of bodies. Increase msurface and &
     &restart calculation. Aborting !'
         stop
      end if      


      return 

      end 

!-----------------------------------------------------------------------
      subroutine read_cgnsbc(bcsize,bctype,bcrange,maxbcsize,basenum, &
           family_names,family_types,n_surblo,msurface,m_surblo,i_surbl, &
           cosabctypes,nblocks,fid,block_number_mapping,amcontrol,rank)
!-----------------------------------------------------------------------
      
      use cgns
      use cosa_precision

      implicit none
      
      integer(kind=cosa_int), allocatable, intent(inout) :: bcsize(:)
!     bcrange is 6 as we are considering 3d data so there are 
!     lower and upper bounds for each of the 3 dimensions
      integer(kind=cosa_int), allocatable, intent(inout) :: bctype(:,:)
      integer(kind=cosa_int), allocatable, intent(inout) :: bcrange(:,:,:)
      character(len=32), allocatable, intent(in) :: family_names(:)
      integer(kind=cosa_int), allocatable, intent(inout) :: family_types(:)
      integer(kind=cosa_int), intent(inout) :: i_surbl(msurface,m_surblo)
      integer(kind=cosa_int), intent(inout) :: n_surblo(msurface)
      integer, dimension(25), intent(in) :: cosabctypes
      integer(kind=cosa_int), intent(in) :: block_number_mapping(nblocks)
      integer(kind=cosa_int) fid,nblocks,iblock,maxbcsize,ibc,ptsettype,basenum, &
           zonenum,datatype,numdataset,ierr,m_surblo, msurface
      logical, intent(in) :: amcontrol
!     NormalIndex is length 3 because we are considering 3d data
!     For 2d data this should be changed 
      integer(kind=cosa_int) NormalIndex(3)
      integer*8 npts
      integer(cgsize_t) :: NormalListFlag
      real :: NormalList(100)
      integer :: family_position, rank

      integer(cgsize_t), allocatable :: temp(:)

      logical :: expect_familyname      
      character*32 bcname 
      
      allocate(temp(6))

      do iblock = 1,nblocks
         
         zonenum = iblock
         
         do ibc = 1,bcsize(iblock)
            call cg_boco_info_f(fid,basenum,zonenum,ibc,bcname, &
                 bctype(ibc,iblock),ptsettype,npts, &
                 NormalIndex,NormalListFlag,datatype, &
                 numdataset,ierr)
            
            ! Convert from CGNS boundary condition types to COSA boundary condition types.
            ! We do this here so we can assign the surface boundary conditions numbers
            ! (i.e. 1500, 1501, 1502, etc...) based on the specific surface/body they are
            ! later on in this routine.
            if(bctype(ibc,iblock) .ne. 25 .and. bctype(ibc,iblock) .ne. 1 .and. cosabctypes(bctype(ibc,iblock)) .eq. -1) then
               write(*,*) 'Problem converting the boundary condition ', &
                    'type from cgns. CGNS boundary type is ', &
                    bctype(ibc,iblock),'  which is currently unmapped ', &
                    'in COSA'
               call abortmpi()
            else
            ! If the cgns boundary condition type is 25
            ! (FamilySpecified) then make sure we check 
            ! there is a family name specified.
               expect_familyname = .false.
               if(bctype(ibc,iblock) .eq. 25) then
                  expect_familyname = .true.
               end if
               bctype(ibc,iblock) = cosabctypes(bctype(ibc,iblock))
               if(ptsettype.eq.PointRange) then
                  call cg_boco_read_f(fid, basenum, zonenum, ibc, temp, &
                       NormalList, ierr)
                  if(ierr .ne. 0) then
                     write(*,*) 'Problem reading cg_ptset_read_f',ierr, &
                          iblock
                     stop
                  end if
                  bcrange(:,ibc,iblock) = temp(1:6)               
               else
                  write(*,*) 'ERROR: point set type not range, quiting!'
                  stop
               end if
            
               call cg_goto_f(fid, basenum, ierr, 'Zone_t', &
                    zonenum,'ZoneBC_t', 1, 'BC_t', &
                    ibc, 'end')
               if(ierr .ne. 0) then
                  write(*,*) 'Problem with cg_goto_f',ierr, &
                       iblock,ibc
                  stop
               end if
               call cg_famname_read_f(bcname, ierr)
               if(ierr .eq. CG_NODE_NOT_FOUND) then
               ! This is fine, a lot of boundary conditions won't have a family name (i.e. won't be a surface)
                  if(expect_familyname) then
                     write(*,*) 'Expecting a family name because the CGNS ', &
      'BC type is 25. Aborting.'
                     call abortmpi()
                  end if
               else if(ierr .ne. 0) then
               ! If we get to this point it's not just it can't be found, something else has gone wrong
                  write(*,*) 'Problem reading cg_famname_read_f',ierr, &
                       iblock,ibc
                  stop
               else
               ! If we find a family name, this means it is part of a surface so we need to add
               ! it to the specific surface data structure (i.e. i_surbl)            
                  family_position = find_family_name(family_names, bcname)
                  if(family_position .eq. -1) then
                     write(*,*) 'Found family name',bcname, &
                          'that does not match any family name we are ', &
                          'expecting'
                     call abortmpi()
                  end if
                  n_surblo(family_position) = n_surblo(family_position) + 1
                  if (n_surblo(family_position).gt.m_surblo) then
                     write(*,*) 'Exceeded max. number of wall blocks per body. In &
     &crease m_surblo and restart calculation. Aborting!'
                     stop
                  end if
                  i_surbl(family_position,n_surblo(family_position)) = zonenum
               ! If we have a surface boundary condition then renumber it to be 
               ! 1500, 1501, 1502, etc... depending on which surface it is.
               ! We want this to start from 1500 for the first surface/body
               ! and then increment by one from there. Provided family_position
               ! is correct, this should do that.
                  if(cosabctypes(family_types(family_position)) .eq. 1500) then
                     family_position = find_surface_name_bc_type(family_names, &
                       family_types, family_position, bcname)
                     if(family_position .eq. -1) then
                        write(*,*) 'Found family name',bcname, &
                             'that does not match any family name we are ', &
      'expecting'
                        stop
                     end if
                     bctype(ibc,iblock) = 1500 + (family_position - 1)
                  else if(cosabctypes(family_types(family_position)) .eq. 1400) then
                     write(*,*) '1400'
                     family_position = find_surface_name_bc_type(family_names, &
                          family_types, family_position, bcname)
                     if(family_position .eq. -1) then
                        write(*,*) 'Found family name',bcname, &
                             'that does not match any family name we are ', &
      'expecting'
                        stop
                     end if
                     bctype(ibc,iblock) = 1400 + (family_position - 1)
                  else
!     AJ We are using the CGNS type 1 of a CGNS family as a special value to
!     AJ mark up boundaries we can ignore because they are symmetric or periodic.
!     AJ TODO: We should revisit this, and potentially get rid of such specialisation
!     AJ TODO: by ensuring that this is fixed in CGNS initially.
                     if(family_types(family_position) .ne. 1 .and. &
                          cosabctypes(family_types(family_position)) .eq. -1) then
                        write(*,*) 'Problem converting the boundary condition ', &
      'type from cgns. CGNS boundary type is ', &
                          family_types(family_position),'  which is currently unmapped ', &
      'in COSA'
                        call abortmpi()
!     AJ We are using the CGNS type 1 of a CGNS family as a special value to
!     AJ mark up boundaries we can ignore because they are symmetric or periodic.
!     AJ TODO: We should revisit this, and potentially get rid of such specialisation
!     AJ TODO: by ensuring that this is fixed in CGNS initially.
                     else if(family_types(family_position) .eq. 1) then
                        bctype(ibc,iblock) = -1
                     else
!                        write(*,*) 'Not a family type ',family_types(family_position),
!     &                       cosabctypes(family_types(family_position))
                        bctype(ibc,iblock) = cosabctypes(family_types(family_position))
                     end if
                  end if
               end if
            end if
         end do
      end do
      
      deallocate(temp)

      return 

      end 

!-----------------------------------------------------------------------
      subroutine read_cgnsgc(gcsize,gcrange,gctransform,maxgcsize, &
           blockname,donornum,basenum,nblocks,fid, &
           zrotangle,translength,pcut,block_number_mapping)
!-----------------------------------------------------------------------
      
      use cgns
      use cosa_precision

      implicit none
      
!     gcrange is 6,2 as we are considering 3d data so there are 
!     lower and upper bounds for each of the 3 dimensions, 
!     and we get the ranges for this block and the donor block
      integer(kind=cosa_int) gcsize(nblocks)
      integer(cgsize_t) gcrange(6,2,maxgcsize,nblocks)
!               gctransform is 3 as we are doing 3d simulations
      integer(kind=cosa_int) gctransform(3,maxgcsize,nblocks)
      integer(kind=cosa_int) fid,nblocks,iblock,maxgcsize,igc,ptsettype
      integer(kind=cosa_int) basenum,zonenum,gclocation,gctype
      integer(kind=cosa_int), intent(in) :: block_number_mapping(nblocks)
      integer(kind=cosa_int) nzonetype,ndonorptstype,donordatatype,npntsdonor,ierr
      integer(kind=cosa_int) NormalIndex(3),donornum(maxgcsize,nblocks),inum
      integer*8 npts
      real*4    rotcenter(3),rotangle(3),translation(3) 
      real*4, allocatable, intent(inout) ::  zrotangle(:,:)
      real*4, allocatable, intent(inout) :: translength(:,:)
      character*32 gcname,donorname,blockname(nblocks)
      logical   pcut(maxgcsize,nblocks), found
      
      donornum = -1

      do iblock = 1,nblocks
         
         zonenum = iblock

         do igc = 1,gcsize(iblock)
            call cg_1to1_read_f(fid,basenum,zonenum,igc,gcname,donorname, &
                 gcrange(:,1,igc,iblock), &
                 gcrange(:,2,igc,iblock), &
                 gctransform(:,igc,iblock),ierr)

            gcname = trim(gcname)
            donorname = trim(donorname)
            
            call cg_1to1_periodic_read_f(fid,basenum,zonenum,igc, &
                 rotcenter,rotangle,translation, &
                 ierr) 
            if (ierr.eq.2) then
!     dbg        write(*,*) '1-to-1 Connectivity found. Cuts ONLY!'
            else
               pcut(igc,iblock) = .true.
!     dbg        write(*,*) '1-to-1 Connectivity found. Cuts and periodic cut
!     dbg &s!'
               if (any(abs(rotcenter).gt.1d-7)) write(*,*) 'WARNING, Rotation Ce &
     &nter NOT 0!'
               
               if (any(abs(rotangle).gt.1d-7)) then
                  if (abs(rotangle(1)).gt.1d-7) then
                     write(*,*) 'Rotation around X-axis currently not support &
     &ed by COSA. Please use Z-axis as rotational axis!'
                     write(*,*) rotangle(1)
                     stop
                  else if (abs(rotangle(2)).gt.1d-7) then
                     write(*,*) 'Rotation around Y-axis currently not support &
     &ed by COSA. Please use Z-axis as rotational axis!'
                     write(*,*) rotangle(2)
                     stop
                  else if (abs(rotangle(3)).gt.1d-7) then
                     zrotangle(igc,iblock) = rotangle(3)
                  end if 
               end if
               
               if (any(abs(translation).gt.1d-7)) then
                  if (abs(translation(1)).gt.1d-7) then
                     translength(igc,iblock) = translation(1)
                  else if (abs(translation(2)).gt.1d-7) then
                     translength(igc,iblock) = translation(2)
                  else if (abs(translation(3)).gt.1d-7) then
                     translength(igc,iblock) = translation(3)
                  end if 
               end if

            end if
            
            found =  .false.
            do inum = 1,nblocks
               if (blockname(inum).eq.donorname) then
                  donornum(igc,iblock) = inum
                  found = .true.
                  exit
               end if 
            end do
            if(.not. found) then
               write(*,*) donorname,'not found'
            end if


         end do
      end do
      
      return 
      
      end 

      end module cosa_cgns_utilities




!-----------------------------------------------------------------------
      subroutine input()
!-----------------------------------------------------------------------
 
! epslim        : limiter cutoff used within MUSCL extrapolations
! cntrpy        : entropy fix multiplicative constant:
!                 =0, NO entropy fix is applied
!                 >0,    entropy fix is applied
! entfxctf      : entropy fix cutoff used when etpfxtyp=0
! etpfxtyp      : entropy fix type:
!                 0, constant eigenvalue limiting (entfx_w=1)
!                 1, variable eigenvalue limiting (entfx_w=2)
!                 2, variable eigenvalue limiting (entfx_w=3)
! lmax          : number of MG iterations
! nsurface         : number of distinct bodies for which independent
!                 calculation of lift, drag, and possibly moment are
!                 requested
! n_surblo(msurface):
!                 array each element of which has number of blocks
!                 defining complete boundary of each of the nsurface
!                 bodies. The max. value of the actual number
!                 of bodies of the simulation (nsurface) is msurface, a
!                 parameter assigned in common.block.
! i_surbl(msurface,m_surblo):
!                 first matrix index refers to body (such index goes
!                 from 1 to nsurface). Second matrix index is pointer to
!                 block providing a boundary patch of the complete
!                 boundary of the current body.
! functional option:
!   can take two values: default or none. Informer case logical calfor is 
!   true, in latter case logical calfor is false. For default case, force 
!   and moment vectors are projected along directions depending on the logicals aircraft, hawt and vawt.
!   aircraft : freestream angles alpha and beta define freestream 
!              direction s; unit vector b is orthogonal to s and parallel
!              to xz (horizontal plane); unit vector a is orthogonal to 
!              both s and b. Lift is force component along a; drag is
!              force component along s; lateral force is force component
!              along b. Moments are computed about a given pole of 
!              coordinates xmom,ymom,zmom and projected along the 
!              axis x, y and z.
!   hawt     :
!   vawt     :
!   ROADWORKS
! cdff          : weight of viscous eigenvalues in calculation of
!                 local time step. Value of 1 may give faaster 
!                 convergence but lead to instability. Value of 4 is
!                 more conservative. Code expects 1 ge cdff le 4.
! ilim          : type of limiter:
!                 0, no limiter
!                 1, smooth limiter
!                 2, minmod limiter
!                 3, Venkatakrishnan smooth limiter
!                 4, Vanalbada smooth limiter
! ntime         : for dual-time-stepping, this is number of 
!                 time-intervals when:
!                 a. irest=1, d_prop_start=t.
!                 b. irest=1, d_prp_strt=f., fnd_simtime=t. 
!                    (e.g. time-dependent RK vortex TC analysis)
!               : for dual-time-stepping, this is ( number of 
!                 time-intervals - 1) when:
!                 d. irest=1, d_prp_strt=f., fnd_simtime=f.
!                    (Richardson extrapolation for cyl. shed vorticity)
!                 f. irest=0
! omega         : vibration frequency (Hertz)
! profpar(1)    : block from which bl profile has to be extracted
! profpar(2)    : i-value of profile
! profpar(3)    : j-start of progile
! profpar(4)    : j-end   of progile
! rkap          : order of MUSCL extrapolation when nl_crs.gt.1
! lim_corr1     : straightforward positivity enforcemnt of sol. after prol 
! lim_corr2     : correction of prolongation after Wasserman
! lim_corr3     : correction of prolongation after Eliasson
!
! irot          : a periodic cut is described by 2 lines of input.dat, 
!                 one defining the cut in block 1 coordinates
!                 the other in block 2 coordinates. In pcut_q
!                 The 'boundary point' is taken as the
!                 cut seen by block 1, so boundary cells of
!                 block 1 are set to values of interior cells of block 2.
!                 Suppose that periodic boundary 1 is brought onto 
!                 periodic boundary 2
!                 by performing a rotation in positive z direction.
!                 Then data of periodic boundary 2 are reported to
!                 periodic boundary 1 by performing a NEGATIVE rotation.
!                 In this circumstance, the convention is that irot in
!                 input.dat (on the second line) must be negative. COSA
!                 then sets irot of the first line equal to -irot=+1.

!-----------------------------------------------------------------------

!     Remark 1: theory on MUSCL extrapolation and limiters in chapter 21 
!               of Hirsch.

!-----------------------------------------------------------------------

      use cosa_cgns_utilities
      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax(mblock),jmax(mblock),kmax(mblock),ijkmax(mblock), &
                nwallv_1(mwls,mblock),nwallv_2(mwls,mblock), &
                nvwls(mblock)
      integer(kind=cosa_int) idir,istrt(3),iend(3),iord(3),irot,nbnd,ic1,ic2,iblk,i, &
                l,icut,ibc,n,nl,iblock,ibody,i_bodblo,iwl
      integer(kind=cosa_int) ierr
 
      real(kind=cosa_real) kappa

      character*10 argo,eqs_typ,cgns_type
      character*40 arg
      character*10000 line
      character*72 bcfile
      character*200 templine
      
!      data ar4      /0.25,0.333333,0.5,1.0/   ! Jameson's central scheme
!c     data ar4      /0.1084,0.2601,0.5051,1.0/! Van leer's upwind 
!c                                               (AIAA Journal 1995)
!c     data ar4      /0.125,0.306,0.587,1.0/   ! Eliasson's coefficiencts
!                                               (NUMECA upwind/central)
!      data br4      /0.0,0.25,0.333333,0.5/
!      data alphadum /1.0,0.0,0.0,0.0,0.0/
!      data ar5      /0.25,0.166666,0.375,0.5,1.0/
!      data cyc      /1,2,3,2,3,1,3,1,2/
!      data krd      /1,0,0,0,1,0,0,0,1/

      logical amcontrol

      nstagerk = 4

! Jameson's central scheme
      ar4 = (/0.25,0.333333,0.5,1.0/)
! Van leer's upwind (AIAA Journal 1995)
!  ar4 = (/0.1084,0.2601,0.5051,1.0/)
! Eliasson's coefficiencts (NUMECA upwind/central)
!  ar4 = (/0.125,0.306,0.587,1.0/)   
      br4 = (/0.0,0.25,0.333333,0.5/)
      alphadum = (/1.0,0.0,0.0,0.0,0.0/)
      ar5 = (/0.25,0.166666,0.375,0.5,1.0/)
      cyc = reshape((/1,2,3,2,3,1,3,1,2/),(/3,3/))
      krd = reshape((/1,0,0,0,1,0,0,0,1/),(/3,3/))

      call amcontroller(amcontrol)

      aircraft        = .false.
      hawt            = .false.
      vawt            = .false.
      persec          = .false.
      relframe        = .false.
      unsteady        = .false.
      moving          = .false.
      rotating        = .false.
      tpitch          = .false.
      pitching        = .false.
      plunging        = .false.
      plupitching     = .false.
      lomach          = .false.
      visprec         = .false.
      dualt           = .false.
      rgkuns          = .false.
      harbal          = .false.
      dual_prop_start = .false.
      debug           = .false.
      rkimp           = .false.
      viscous         = .false.
      wldata          = .false.
      blprof          = .false.
      calvort         = .false.
      calcp           = .false.
      foturb          = .false.
      calmut          = .false.
      bln_prol        = .true.
      found_simtime   = .false.
      ramping1        = .false.
      ramping2        = .false.
      lim_corr(1)     = .false.
      lim_corr(2)     = .false.
      lim_corr(3)     = .false.
      ho_rest         = .false.
      prd_all_lvls    = .false.
      tecbin          = .false.
      dgn_cell        = .false.

      efx_all_f       = .false.
      efx_all_t       = .false.
      eig_lim_all_f   = .false.
      eig_lim_all_t   = .false.

      using_cgns      = .false.

      ainf            = 1.d0
      rhoinf          = 1.d0
      tinf            = 1.d0

      pi      = dacos(-1.d0)
      nblade  = 0
      omega   = 0.d0
      omegas  = 0.d0
      omegatp = 0.d0
      phitp   = 0.d0
      gtheta  = 0.d0
      simtime = 0.d0
      xh      = 0.d0
      yh      = 0.d0
      dtheta0 = 0.d0
      dh0x    = 0.d0
      dh0y    = 0.d0
      machfs  = 0.d0
      alpha   = 0.d0
      beta    = 0.d0
!old  epsp    = 1.d-05
      lsp_vfr = 0
      mixpi   = 0
      uprv    = 0.d0
      iprint  = 10
      nharms  = 0
      nsave   = 1
      entfx_f_w = 0
      entfx_t_w = 0

      rotmat(1,1) = 1.d0
      rotmat(1,2) = 0.d0
      rotmat(2,1) = 0.d0
      rotmat(2,2) = 1.d0

      zeit(0)  = 0.d0
      wrifl(1) = 0
      wrifl(2) = 0

      prant   = 0.9d0

      qmin(1)  = 1.d-8
      qmin(2)  = 1.d0
      qmin(3)  = 1.d0
      qmin(4)  = 1.d0
      qmin(5)  = 1.d-8
      qmin(6)  = 1.d-12
      qmin(7)  = 1.d-8

      tref     = 0.d0

!---- K-Omega constants and parameters
      epst     = 1.d-12
      kappa    = 0.41d0
      cmu      = 1.d0
      sigk     = 0.5d0
      sigw     = 0.5d0
      bkw      = 3.d0 / 40.d0
      bstrkw   = 0.09d0
      gkw      = bkw/bstrkw - sigw*kappa**2/dsqrt(bstrkw)
!old  gkw      = 5.d0 / 9.d0

!---- K-Omega BSL/SST constants and parameters
      epst     = 1.d-12
      kappa    = 0.41d0
      bstr     = 0.09d0
      a1       = 0.31d0

      sigk1_bsl= 0.5d0
      sigk1_sst= 0.85d0
      sigk2    = 1.0d0

      sigw1    = 0.5d0
      b1       = 3.d0 / 40.d0
      g1       = b1/bstr - sigw1*kappa**2/dsqrt(bstr)

      sigw2    = 0.856d0
      b2       = 0.0828d0
      g2       = b2/bstr - sigw2*kappa**2/dsqrt(bstr)

      open(4,file='input.dat')
      rewind(4)

      read(4,100)
      read(4,'(a)') line
      call linarg(1,line,arg)
      read(arg,'(a)') cgns_type
      if (cgns_type .eq. 'cgns') then
         read(4,'(a)') templine
         cgns_filename = trim(templine)
         if(amcontrol) then
            write(*,*) 'cgns file name'
            write(*,*) '------------'      
            write(*,*) trim(cgns_filename)
            write(*,*) '------------'
         end if
         read(4,100)
         using_cgns = .true.
      else
         do iblock = 1,mblock
            do iwl = 1,mwls
               nwallv_1(iwl,iblock) = 0
               nwallv_2(iwl,iblock) = 0
            end do
            nvwls(iblock) = 0
         end do
      end if
      read(4,'(a)') line
      call linarg(1,line,arg)
      if (arg.eq.'y') then
        debug = .true.
      else if (arg.eq.'n') then
        debug = .false.
      else
        write(*,*) 'Unknown debugging option. Aborting!'
        stop
      end if
      call linarg(2,line,arg)
      read(arg,'(a)') eqs_typ
      if ((eqs_typ.eq.'euler')) then
        npde = 5
      else if (eqs_typ.eq.'laminar') then
        npde = 5
        viscous = .true.
        kom     = .false.
        kom_bsl = .false.
        kom_sst = .false.
      else if (eqs_typ.eq.'komega') then
        npde = 7
        viscous = .true.
        kom     = .true.
        kom_bsl = .false.
        kom_sst = .false.
      else if (eqs_typ.eq.'bsl') then
        npde = 7
        viscous = .true.
        kom     = .false.
        kom_bsl = .true.
        kom_sst = .false.
      else if (eqs_typ.eq.'sst') then
        npde = 7
        viscous = .true.
        kom     = .false.
        kom_bsl = .false.
        kom_sst = .true.
      else
        write(*,*) 'unknown flow model !'
        stop
      end if
      call linarg(3,line,arg)
      if (arg.eq.'external') then
        exter = .true.
      else if (arg.eq.'internal') then
        exter = .false.
      else
        write(*,*) 'Unknown flow type option. Aborting!'
        stop
      end if
      call linarg(4,line,arg)
      read(arg,'(a)') argo
      if (argo.eq.'aircraft') then
        aircraft=.true.
      else if (argo.eq.'saircraft') then
        aircraft=.true.
        persec  =.true.
      else if (argo.eq.'vawt') then
        vawt    =.true.
      else if ((argo.eq.'rhawt').or.(argo.eq.'shawt')) then
        hawt    =.true.
        if (argo.eq.'shawt') then
          persec = .true.
          call linarg(5,line,arg)
          read(arg,*) nblade
          if ((nblade.lt.0).or.(nblade.gt.100)) then
            write(*,*) 'Number of blades outside allowed range. Aborting &
     &!'
            stop
          end if
          gtheta = 2 * pi / nblade
          rotmat(1,1) =  cos(gtheta)
          rotmat(1,2) = -sin(gtheta)
          rotmat(2,1) =  sin(gtheta)
          rotmat(2,2) =  cos(gtheta)
        end if
      else
        write(*,*) 'Unknown configuration. Aborting!'
        stop
      end if
      if(amcontrol) then
        write(*,*) eqs_typ
      end if

      read(4,100)
      read(4,'(a)') line
      call linarg(1,line,arg)
      read(arg,*) gamma
      if (viscous) then
        call linarg(2,line,arg)
        read(arg,*) reyno
        call linarg(3,line,arg)
        read(arg,*) pranl
        call linarg(4,line,arg)
        read(arg,*) machfs
        if (exter) then
          call linarg(5,line,arg)
          read(arg,*) alpha
          call linarg(6,line,arg)
          read(arg,*) beta
        end if
        if(amcontrol) then
          write(*,*) 'gamma, machfs, alpha, beta, reyno', gamma, machfs, alpha, beta, reyno
        end if
      else if (.not.viscous) then
        if (.not.exter) then
          continue
        else 
          call linarg(2,line,arg)
          read(arg,*) machfs
          call linarg(3,line,arg)
          read(arg,*) alpha
          call linarg(4,line,arg)
          read(arg,*) beta
          if(amcontrol) then
            write(*,*) gamma, machfs, alpha, beta
          end if
        end if
      end if
      alpha = alpha / 180 * pi
      beta  = beta / 180 * pi
      pinf  = 1.d0/gamma

      if (kom.or.kom_bsl.or.kom_sst) then
        read(4,100)
!       prant  :  turbulent Pradtl number
!       tkefar : turbulence intensity squared
!       mutfar : (turbulent viscosity) / (laminar viscosity) at farfield
!                Note that (nondim. lam. visc. at farf.) = 1.d0
        read(4,'(a)') line
        call linarg(1,line,arg)
        read(arg,*) prant
        call linarg(2,line,arg)
        read(arg,*) tkefar
        call linarg(3,line,arg)
        read(arg,*) mutfar
        call linarg(4,line,arg)
        if (arg.eq.'wilcox') then
          wallwilc = .true.
          wallment = .false.
          call linarg(5,line,arg)
          read(arg,*) roughk
        else if (arg.eq.'menter') then
          wallwilc = .false.
          wallment = .true.
        end if
         
!       tkefar : convert to turbulence kinetic energy at farfield
!       omefar : omega at farfield
        tkefar = tkefar * machfs**2
        omefar = cmu * rhoinf * tkefar / mutfar * reyno / machfs

        read(4,100)
        read(4,'(a)') line
        call linarg(1,line,arg)
        if (arg.eq.'y') then
          turcmp = 1
        else if (arg.eq.'n') then
          turcmp = 0
        else
          write(*,*) 'Unknown type of positivity preserving type. Aborti &
     &ng!'
          stop
        end if
        call linarg(2,line,arg)
        if (arg.eq.'y') then
          lim_prodtke=.true.
        else if (arg.eq.'n') then
          lim_prodtke=.false.
        else
          write(*,*) 'Unknown type of tke limiting type. Aborting!'
          stop
        end if
        call linarg(3,line,arg)
        read(arg,*) prdlim
        call linarg(4,line,arg)
        if (arg.eq.'y') then
          lim_prodomega=.true.
        else if (arg.eq.'n') then
          lim_prodomega=.false.
        else
          write(*,*) 'Unknown type of omega limiting type. Aborting!'
          stop
        end if
        if (lim_prodomega.and.(.not.lim_prodtke)) then
          write(*,*) 'not advisable to limit omega prod term without lim &
     &iting k production term. Aborting!'
          stop
        end if
        call linarg(5,line,arg)
        if (arg.eq.'none') then
          ifixq = 0
        else if (arg.eq.'minimum') then
          ifixq = 1
        else
          write(*,*) 'Unknown type of solution positivity preserving typ &
     &e. Aborting!'
          stop
        end if
        call linarg(6,line,arg)
        if (arg.eq.'second') then
          foturb = .false.
        else if (arg.eq.'first') then
          foturb = .true.
        else
          write(*,*) 'Unknown type of turbulence order. Aborting!'
          stop
        end if

      end if

      read(4,100)
      read(4,'(a)') line
      call linarg(1,line,arg)
      if ((arg.ne.'steady').and.(arg.ne.'unsteady')) then
        write(*,*) 'Unknown type of flow mode. Aborting!'
        stop
      else if (arg.eq.'steady') then
        continue
      else if (arg.eq.'unsteady') then
        unsteady = .true.
        call linarg(2,line,arg)
        if (arg.eq.'rgkun') then
          rgkuns = .true.
          call linarg(3,line,arg)
          if (arg.eq.' ') then
            write(*,*) 'ERROR: missing physical time step. Aborting!'
            stop
          end if
          read(arg,*) dt
        else if (arg.eq.'dual') then
          dualt = .true.
          call linarg(3,line,arg)
          if (arg.eq.' ') then
            write(*,*) 'ERROR: missing physical time step. Aborting!'
            stop
          end if
          read(arg,*) dt
          call linarg(4,line,arg)
          if (arg.eq.' ') then
            write(*,*) 'Missing ntime for dual-time-stepping. Aborting!'
            stop
          else
            read(arg,*) ntime
          end if
          call linarg(5,line,arg)
          if (arg.eq.' ') then
            write(*,*) 'Missing saving frequency of unsteady restart fil &
     &e. Aborting!'
            stop
          else
            read(arg,*) nsave
          end if
          call linarg(6,line,arg)
          if (arg.eq.' ') then
            write(*,*) 'Missing restart order choice. Aborting!'
            stop
          else
            read(arg,'(a)') arg
            if (arg.eq.'second') then
!
            else if (arg.eq.'third') then
              write(*,*) '3rd time-accurate to be added!'
              stop
            else
              write(*,*) 'Unknown choice restart order. Aborting!'
              stop
            end if
          end if
          call linarg(7,line,arg)
          if (arg.eq.'rkex') then
            rkimp = .false.
          else if (arg.eq.'rkim') then
            rkimp = .true.
          else
            write(*,*) 'Unknown implicit/explicit RungeKutta. Aborting!'
            stop
          end if
        else if (arg.eq.'hb') then
          harbal = .true.
          call linarg(3,line,arg)
          if (arg.eq.'rkex') then
            rkimp = .false.
          else if (arg.eq.'rkim') then
            rkimp = .true.
          else
            write(*,*) 'Unknown implicit/explicit RungeKutta. Aborting!'
            stop
          end if
          call linarg(4,line,arg)
          if (arg.eq.' ') then
            write(*,*) 'Missing nharms for HB integration. Aborting!'
            stop
          else
            read(arg,*) nharms
            if (nharms.gt.mharms) then
              write(*,*) 'Exceeded max. number of harmonics. Aborting!'
              stop
            end if
          end if
        else
          write(*,*) 'unknown unsteady solver. Aborting!'
          stop
        end if
      end if

      if (hawt.and.(.not.unsteady)) then
        read(4,100)
        read(4,'(a)') line
        call linarg(1,line,arg)
        if (arg.eq.'fixed') then
          write(*,*) 'Simulation of motionless HAWT rotor.'
          write(*,*) 'If blades are NOT feathered, requested steady simu &
     &lation may be unstable.'
          if ((alpha.ne.0).or.(beta.ne.0)) then
            if (persec) then
              write(*,*) 'HAWT steady sector simulation is incompatible &
     &with misaligned freestream. Aborting!'
              stop
            end if
            write(*,*) 'WARNING: simulation outcome depends on azimuthal &
     & blade positions.'
          end if
        else if (arg.eq.'rotating') then
          moving   = .true.
          rotating = .true.
          call linarg(2,line,arg)
          if (arg.eq.' ') then
            write(*,*) 'ERROR: missing HAWT rotational frequency!'
            stop
          end if
          read(arg,*) omegas
          omegas = 2*pi*omegas
          if ((alpha.eq.0).and.(beta.eq.0)) then
            relframe=.true.
            call linarg(3,line,arg)
            if (arg.eq.' ') then
              write(*,*) 'Missing x-coordinate of rotor centre. Aborting &
     &!'
              stop
            end if
            read(arg,*) xrotc
            call linarg(4,line,arg)
            if (arg.eq.' ') then
              write(*,*) 'Missing y-coordinate of rotor centre. Aborting &
     &!'
              stop
            end if
            read(arg,*) yrotc
            call linarg(5,line,arg)
            if (arg.eq.' ') then
              write(*,*) 'Missing z-coordinate of rotor centre. Aborting &
     &!'
              stop
            end if
            read(arg,*) zrotc
          else
            write(*,*) 'HAWT simulation with nonzero alpha and/or beta c &
     &annot be steady. Aborting!'
            stop
          end if
        else
          write(*,*) 'Invalid HAWT motion parameter. Aborting!'
          stop
        end if
      end if

      if (unsteady) then
        read(4,100)
        read(4,'(a)') line
        call linarg(1,line,arg)
        if (arg.eq.'pitch') then
          moving = .true.
          pitching = .true.
          call linarg(2,line,arg)
          if (arg.eq.' ') then
            write(*,*) 'ERROR: missing frequency!'
            stop
          end if
          read(arg,*) omegas
          omegas = 2*pi*omegas
          call linarg(3,line,arg)
          if (arg.eq.' ') then
            write(*,*) 'ERROR: missing abscissa of hinge. Aborting!'
            stop
          end if
          read(arg,*) xh
          call linarg(4,line,arg)
          if (arg.eq.' ') then
            write(*,*) 'ERROR: missing ordinate of hinge. Aborting!'
            stop
          end if
          read(arg,*) yh
          call linarg(5,line,arg)
          if (arg.eq.' ') then
            write(*,*) 'ERROR: missing pitching amplitude!'
            stop
          end if
          read(arg,*) dtheta0
          dtheta0 = dtheta0 / 180.d0 * pi
        else if (arg.eq.'plunge') then
          moving = .true.
          plunging = .true.
          call linarg(2,line,arg)
          if (arg.eq.' ') then
            write(*,*) 'ERROR: missing frequency!'
            stop
          end if
          read(arg,*) omegas
          omegas = 2*pi*omegas
          call linarg(3,line,arg)
          if (arg.eq.' ') then
            write(*,*) 'ERROR: missing plunging x-amplitude!'
            stop
          end if
          read(arg,*) dh0x
          call linarg(4,line,arg)
          if (arg.eq.' ') then
            write(*,*) 'ERROR: missing plunging y-amplitude!'
            stop
          end if
          read(arg,*) dh0y
        else if (arg.eq.'plupit') then
          moving = .true.
          plupitching = .true.
          call linarg(2,line,arg)
          if (arg.eq.' ') then
            write(*,*) 'ERROR: missing frequency!'
            stop
          end if
          read(arg,*) omegas
          omegas = 2*pi*omegas
          call linarg(3,line,arg)
          if (arg.eq.' ') then
            write(*,*) 'ERROR: missing plunging x-amplitude!'
            stop
          end if
          read(arg,*) dh0x
          call linarg(4,line,arg)
          if (arg.eq.' ') then
            write(*,*) 'ERROR: missing plunging y-amplitude!'
            stop
          end if
          read(arg,*) dh0y
          call linarg(5,line,arg)
          if (arg.eq.' ') then
            write(*,*) 'ERROR: missing abscissa of hinge. Aborting!'
            stop
          end if
          read(arg,*) xh
          call linarg(6,line,arg)
          if (arg.eq.' ') then
            write(*,*) 'ERROR: missing ordinate of hinge. Aborting!'
            stop
          end if
          read(arg,*) yh
          call linarg(7,line,arg)
          if (arg.eq.' ') then
            write(*,*) 'ERROR: missing pitching amplitude!'
            stop
          end if
          read(arg,*) dtheta0
          dtheta0 = dtheta0 / 180.d0 * pi
          call linarg(8,line,arg)
          if (arg.eq.' ') then
            write(*,*) 'ERROR: missing pitch/plunge phase!'
            stop
          end if
          read(arg,*) phipp
          phipp = phipp / 180.d0 * pi
        else if (arg.eq.'rotating') then
          moving = .true.
          rotating = .true.
          call linarg(2,line,arg)
          if (arg.eq.' ') then
            write(*,*) 'ERROR: missing frequency!'
            stop
          end if
          read(arg,*) omegas
          omegas = 2*pi*omegas
          call linarg(3,line,arg)
          if (arg.eq.' ') then
            write(*,*) 'Missing x-coordinate of rotor centre. Aborting!'
            stop
          end if
          read(arg,*) xrotc
          call linarg(4,line,arg)
          if (arg.eq.' ') then
            write(*,*) 'Missing y-coordinate of rotor centre. Aborting!'
            stop
          end if
          read(arg,*) yrotc
          call linarg(5,line,arg)
          if (arg.eq.' ') then
            write(*,*) 'Missing z-coordinate of rotor centre. Aborting!'
            stop
          end if
          read(arg,*) zrotc
          if (hawt) then
            call linarg(6,line,arg)
            if (arg.eq.' ') then
              write(*,*) 'ERROR: missing rel/abs frame choice.'
              stop
            end if
            read(arg,*) arg
            if (arg.eq.'absolute') then
              relframe = .false.
              if (persec) then
                write(*,*) 'sector unsteady HAWT simulation cannot be pe &
     &rformed in absolute frame due to periodicity BCs. Aborting!'
                stop
              end if
            else if (arg.eq.'relative') then
              relframe = .true.
            else
              write(*,*) 'Unknown rel/abs frame choice. Aborting!'
              stop
            end if
            if ((alpha.eq.0.d0).and.(beta.eq.0.d0)) then
              write(*,*) 'Simulation of stalled (rotating) HAWT rotor.'
              if (persec) then
                write(*,*) 'WARNING: to avoid inaccuracies it may be mor &
     &e appropriate to use a whole-rotor grid or a single-blade grid wit &
     &h freestream boundaries only, i.e. without periodic boundaries.'
              end if
              if (harbal) then
                write(*,*) 'This simulation cannot be performed with the &
     & HB solver. Aborting!'
                stop
              end if
            else
!             Unsteady flow due to yaw or vertical wind component.
!             Stall may also exist.
              if (persec.and.(dualt.or.rgkuns)) then
                write(*,*) 'HAWT TD sector simulation is incompatible wi &
     &th misaligned freestream. Aborting!'
                stop
              end if
              continue
            end if
          end if
        else if (arg.eq.'fixed') then
          if (harbal) then
            write(*,*) 'Unsteady motionless grid problems cannot be perf &
     &ormed with HB solver. Aborting!'
            stop
          end if
          if (hawt) then
            write(*,*) 'TD simulation of possibly stalled flow past moti &
     &onless HAWT rotor'
            if ((alpha.ne.0).or.(beta.ne.0)) then
              write(*,*) 'WARNING: simulation outcome depends on azimuth &
     &al blade positions.'
              if (persec) then
                write(*,*) 'TD motionless HAWT rotor simulation cannot b &
      e performed with sector grid in presence of freestream misaligneme &
     &nt. Aborting!'
                stop
              end if
            end if
            if (persec) then
                write(*,*) 'WARNING: to avoid inaccuracies it may be mor &
      e appropriate to use a whole-rotor grid or a single-blade grid wit &
     &h freestream boundaries only, i.e. without periodic boundaries.'
            end if
          end if
        else
          write(*,*) 'Unknown type of grid motion. Aborting!'
          stop
        end if
      end if

      if ((dualt.or.harbal).and.hawt.and.rotating.and.(.not.persec)) then
        read(4,100)
        read(4,'(a)') line
        call linarg(1,line,arg)
        if (arg.eq.' ') then
          write(*,*) 'Missing tower motion type. Aborting!'
          stop
        else
          read(arg,'(a)') argo
          if (argo.eq.'fixed') then
            continue
          else if (argo.eq.'floating') then
            call linarg(2,line,arg)
            if (arg.eq.'pitch') then
              tpitch = .true.
              call linarg(3,line,arg)
              read(arg,*) theta0tp
              theta0tp = theta0tp /180.d0 * pi
              call linarg(4,line,arg)
              read(arg,*) omegatp
              omegatp  = 2*pi*omegatp
              call linarg(5,line,arg)
              read(arg,*) phitp
              phitp  = phitp /180.d0 * pi
              call linarg(6,line,arg)
              read(arg,*) zhtp
              call linarg(7,line,arg)
              read(arg,*) yhtp
              call linarg(8,line,arg)
              read(arg,*) betaw
              betaw = betaw / 180 * pi
            else
              write(*,*) 'Only pitch presently implemented. Aborting!'
              stop
            end if
          end if
        end if
      end if

      if (harbal) then
        if (nharms.gt.0) then
          if (tpitch) then
            omega = abs(omegatp)
          else
            omega = abs(omegas)
          end if
          dt = 2*pi/omega/(2*nharms+1)
          do n=0,2*nharms
            zeit(n) = n*dt
          end do
          call mathb(nharms)
        else if (nharms.eq.0) then
          unsteady = .false.
          moving   = .false.
        else
          write(*,*) 'Illegal number of HB harmonics. Aborting!'
          stop
        end if
      end if

      if (dualt) then
        read(4,100)
        read(4,'(a)') line
        call linarg(1,line,arg)
        if (arg.eq.' ') then
          write(*,*) 'Missing first time of solution write. Aborting!'
          stop
        else
          read(arg,*) wrifl(1)
          if (wrifl(1).le.0) then
            write(*,*) 'Illegal first time of solution write. Aborting!'
            stop
          end if
        end if
        call linarg(2,line,arg)
        if (arg.eq.' ') then
          write(*,*) 'Missing last time of solution write. Aborting!'
          stop
        else
          read(arg,*) wrifl(2)
          if (wrifl(2).lt.wrifl(1)) then
            write(*,*) 'Illegal last time of solution write. Aborting!'
            stop
          end if
        end if
        call linarg(3,line,arg)
        if (arg.eq.' ') then
          write(*,*) 'Missing frequency of solution write. Aborting!'
          stop
        else
          read(arg,*) wrifl(3)
          if ((wrifl(3).gt.(wrifl(2)-wrifl(1))).and. &
              (wrifl(2).ne.wrifl(1)))                               then
            write(*,*) 'Illegal frequency of solution write. Aborting!'
            stop
          end if
        end if
      end if

      read(4,100)
      if (.not.viscous) then
        read(4,*) irest,srest,vcfl(3),lmax,iupdt,toler
      else
        read(4,*) irest,srest,vcfl(3),cdff,lmax,iupdt,toler
        if ((cdff.lt.1.d0).or.(cdff.gt.4.d0)) then
          if(amcontrol) then
            write(*,*) 'Illegal value of constant cdff. Aborting!'
            stop
          end if
        end if
      end if
      if(amcontrol) then
        write(*,*) irest,srest,vcfl(3),lmax,iupdt,toler
      end if
      if ((unsteady).and.(.not.dualt)) then
        ntime = lmax
      end if
      if ((irest.ne.0).and.(irest.ne.1)) then
        write(*,*) 'Unknown restart option. Aborting'
        stop
      end if

      read(4,100)
      read(4,'(a)') line
      call linarg(1,line,arg)
      read(arg,*) rkap
      call linarg(2,line,arg)
      if (arg.eq.'cirs_v2v') then
        irsop = 6
      else if (arg.eq.'cirs_v1v') then
        irsop = 5
      else if (arg.eq.'cirs_v2') then
        irsop = 3
      else if (arg.eq.'cirs_v1') then
        irsop = 2
      else if (arg.eq.'cirs_c') then
        irsop = 1
      else if (arg.eq.'no_irs') then
        irsop = 0
      else
        write(*,*) 'Unknown Imp. Res. Smoothing option. Aborting!'
        stop
      end if
      if (irsop.gt.0) then
        call linarg(3,line,arg)
        read(arg,*) vcfli(3)
        call linarg(4,line,arg)
        read(arg,*) cutcirs
        if ((cutcirs.ne.0).and.(cutcirs.ne.1)) then
          write(*,*) 'Unknown IRS residual cut option. Aborting!'
          stop
        end if
        if ((irsop.eq.3).or.(irsop.eq.6)) then
          call linarg(5,line,arg)
          read(arg,*) psirs
        end if
      end if
      if(amcontrol) then
        write(*,*) 'rkap,irsop,vcfli(3),cutcirs,psirs', &
                    rkap,irsop,vcfli(3),cutcirs,psirs
      end if

      if (kom.or.kom_bsl.or.kom_sst) then
        read(4,100)
        read(4,'(a)') line
        call linarg(1,line,arg)
        read(arg,*) vcflt(3)
        call linarg(2,line,arg)
        read(arg,*) vcflit(3)
        call linarg(3,line,arg)
        if (arg(1:7).eq.'ramping') then
          if (arg(8:8).eq.'1') then
            ramping1 = .true.
          else if (arg(8:8).eq.'2') then
            ramping2 = .true.
          else
            write(*,*) 'Unknown ramping option. Aborting!'
            stop
          end if
          call linarg(4,line,arg) 
          read(arg,*) ncycr(3)
          call linarg(5,line,arg) 
          read(arg,*) ncycr(2)
          call linarg(6,line,arg) 
          read(arg,*) ncycr(1)
          if (ramping2) then
            call linarg(7,line,arg) 
            read(arg,*) ncycr(0)
          end if

          read(4,100)
          read(4,'(a)') line
          call linarg(1,line,arg) 
          read(arg,*) vcfli(2)
          call linarg(2,line,arg) 
          read(arg,*) vcflit(2)
          call linarg(3,line,arg) 
          read(arg,*) vcfli(1)
          call linarg(4,line,arg) 
          read(arg,*) vcflit(1)
          if (dualt) then
            call linarg(5,line,arg) 
            read(arg,*) ntime_ramp
          end if
        end if
      end if

      read(4,100)
      read(4,'(a)') line
      call linarg(1,line,arg)
      read(arg,*) ilimf
      call linarg(2,line,arg)
      read(arg,*) epslimf
      call linarg(3,line,arg)
      if (arg.eq.'y') then
         efx_setup_f = .true.
      else if (arg.eq.'n') then
         efx_setup_f = .false.
      else
         write(*,*) 'unknown option for E-fix-F set-up. Aborting!'
         stop
      end if
      if (efx_setup_f) then
        call linarg(4,line,arg)
        read(arg,*) cntrpy_f
        if (cntrpy_f.lt.0.d0.or.cntrpy_f.gt.1.d0) then
          write(*,*) 'Inadmittable coefficient for entropy fix F. Aborting!'
          write(*,*) 'cntrpy_F must be 0 <= cntrpy <= 1.'
          stop
        end if
        call linarg(5,line,arg)
        read(arg,*) etpfxtyp_f
        if (etpfxtyp_f.ne.0.and.etpfxtyp_f.ne.1.and.etpfxtyp_f.ne.2) then
          write(*,*) 'Inadmittable entropy fix F type. Aborting!'
          stop
        end if
        if (etpfxtyp_f.eq.0) entfx_f_w(1) = 1
        if (etpfxtyp_f.eq.1) entfx_f_w(2) = 1
        if (etpfxtyp_f.eq.2) entfx_f_w(3) = 1
        call linarg(6,line,arg)
        read(arg,*) entfxctf_f
        call linarg(7,line,arg)
        if (arg.eq.'y') then
          efx_all_f = .true.
        else if (arg.eq.'n') then
          efx_all_f = .false.
        else
          write(*,*) 'Invalid option for efx_all_f. Aborting!'
          stop
        end if
        if(amcontrol) then
          write(*,*) &
            'ilimf,epslimf,cntrpy_f,etpfxtyp_f,entfxctf_f', &
             ilimf,epslimf,cntrpy_f,etpfxtyp_f,entfxctf_f
        end if
      end if

      read(4,100)
      read(4,'(a)') line
      call linarg(1,line,arg)
      if (arg.eq.'y') then
        eig_limit_f = .true.
      else if (arg.eq.'n') then
        eig_limit_f = .false.
      else
        write(*,*) 'Invalid option for eig_limit_f. Aborting!'
        stop
      end if
      if (eig_limit_f) then
        call linarg(2,line,arg)
        read(arg,*) eig_cutoff_f
        if (eig_cutoff_f.lt.0.d0.or.eig_cutoff_f.gt.1.d0) then
          write(*,*) 'Invalid value for eig_cutoff_f. Aborting!'
          write(*,*) 'one must have 0 <= F eig_cutoff_f <= 1.'
          stop
        end if
        call linarg(3,line,arg)
        if (arg.eq.'y') then
          eig_lim_all_f = .true.
        else if (arg.eq.'n') then
          eig_lim_all_f = .false.
        else
          write(*,*) 'Invalid option for eig_lim_all_f. Aborting!'
          stop
        end if
      end if

      read(4,100)
      read(4,'(a)') line
      call linarg(1,line,arg)
      read(arg,*) ilimt
      call linarg(2,line,arg)
      read(arg,*) epslimt
      call linarg(3,line,arg)
      if (arg.eq.'y') then
        efx_setup_t = .true.
      else if (arg.eq.'n') then
        efx_setup_t = .false.
      else
        write(*,*) 'unknown option for E-fix-T set-up. Aborting!'
        stop
      end if
      if (efx_setup_t) then
        call linarg(4,line,arg)
        read(arg,*) cntrpy_t
        if (cntrpy_t.lt.0.d0.or.cntrpy_t.gt.1.d0) then
          write(*,*) 'Inadmittable coefficient for entropy fix T. Aborting!'
          write(*,*) 'cntrpy_T must be 0 <= cntrpy <= 1.'
          stop
        end if
        call linarg(5,line,arg)
        read(arg,*) etpfxtyp_t
        if (etpfxtyp_t.ne.0.and.etpfxtyp_t.ne.1.and.etpfxtyp_t.ne.2) then
          write(*,*) 'Inadmittable entropy fix T type. Aborting!'
          stop
        end if
        if (etpfxtyp_t.eq.0) entfx_t_w(1) = 1
        if (etpfxtyp_f.eq.1) entfx_t_w(2) = 1
        if (etpfxtyp_f.eq.2) entfx_t_w(3) = 1
        call linarg(6,line,arg)
        read(arg,*) entfxctf_t
        call linarg(7,line,arg)
        if (arg.eq.'y') then
          efx_all_t = .true.
        else if (arg.eq.'n') then
          efx_all_t = .false.
        else
          write(*,*) 'unknown option for efx_all_t set-up. Aborting!'
          stop
        end if
        if(amcontrol) then
          write(*,*) &
            'ilimt,epslimt,cntrpy_t,etpfxtyp_t,entfxctf_t', &
             ilimt,epslimt,cntrpy_t,etpfxtyp_t,entfxctf_t
        end if
      end if

      read(4,100)
      read(4,'(a)') line
      call linarg(1,line,arg)
      if (arg.eq.'y') then
        eig_limit_t = .true.
      else if (arg.eq.'n') then
        eig_limit_t = .false.
      else
        write(*,*) 'unknown option for T eigenval. lower bound. Aborting!'
        stop
      end if
      if (eig_limit_t) then
        call linarg(2,line,arg)
        read(arg,*) eig_cutoff_t
        if (eig_cutoff_t.lt.0.d0.or.eig_cutoff_t.gt.1.d0) then
          write(*,*) 'Invalid value for eig_cutoff_t. Aborting!'
          write(*,*) 'one must have 0 <= T eig_cutoff_t <= 1.'
          stop
        end if
        call linarg(3,line,arg)
        if (arg.eq.'y') then
          eig_lim_all_t = .true.
        else if (arg.eq.'n') then
          eig_lim_all_t = .false.
        else
          write(*,*) 'Invalid option for eig_lim_all_t. Aborting!'
          stop
        end if
      end if

!---- multigrid control parameters
      read(4,100)
      read(4,*) nlevel,nl_crs,nl_fmg,nstart,npre,npost,ncrs
      if ((nlevel.lt.1).or.(nlevel.gt.3)) then
        write(*,*) &
        '   Illegal value of multigrid levels. Aborting!'
        stop
      end if

!---- if MG, select interpolation type for prolongation
      if (nlevel.gt.1) then
        read(4,100)
        read(4,'(a)') line
        call linarg(1,line,arg)
        read(arg,'(a)') argo
        if (argo.eq.'linear') then
          bln_prol = .false.
        else if (argo.eq.'bilinear') then
          bln_prol = .true.
        else
          write(*,*) &
        '   Illegal interpolation type for prolongation. Aborting!'
          stop
        end if
!
        if (kom.or.kom_bsl.or.kom_sst) then
!-------- turb MG inp 1
          call linarg(2,line,arg)
          read(arg,'(a)') argo
          if (argo.eq.'ho_rest') then
            ho_rest = .true.
          else if (argo.eq.'lo_rest') then
            ho_rest = .false.
          else
            write(*,*) 'Illegal option for order of MG restriction. Abor &
     &ting!'
            stop
          end if
!-------- turb MG inp 2
          call linarg(3,line,arg)
          read(arg,'(a)') argo
          if (argo.eq.'lim_corr1') then
            lim_corr(1) = .true.
            lim_corr(2) = .false.
            lim_corr(3) = .false.
          else if (argo.eq.'lim_corr2') then
            lim_corr(1) = .false.
            lim_corr(2) = .true.
            lim_corr(3) = .false.
          else if (argo.eq.'lim_corr3') then
            lim_corr(1) = .false.
            lim_corr(2) = .false.
            lim_corr(3) = .true.
          else
            write(*,*) 'Illegal option for positivity preservation after &
     & FG correction. Aborting!'
            stop
          end if
!-------- turb MG inp 3
          call linarg(4,line,arg)
          read(arg,'(a)') argo
          if (argo.eq.'y') then
            prd_all_lvls = .true.
          else if (argo.eq.'n') then
            prd_all_lvls = .false.
          else
            write(*,*) 'Illegal option for input prd_all_lvls. Aborting! &
     &'
            stop
          end if
        end if
!
      end if

!---- low-speed-preconditioning control parameters
      read(4,100)
      read(4,'(a)') line
      call linarg(1,line,arg)
      if (arg.eq.'lomach') then
        lomach = .true.
        call linarg(2,line,arg)
        if (arg.eq.' ') then
          write(*,*) 'Missing cutoff option. Aborting!'
          stop
        else
          read(arg,*) cutoff_type
          if ((cutoff_type.lt.1).or.(cutoff_type.gt.3)) then
            write(*,*) 'Illegal cutoff type. Aborting!'
            stop
          end if
        end if
        call linarg(3,line,arg)
        if (arg.eq.' ') then
          write(*,*) 'Missing visc./inviscid prec. option. Aborting!'
          stop
        else if (arg.eq.'viscous') then
          if (.not.viscous) then
            write(*,*) 'Cannot use visc. preconditioning in inviscid run &
     &. Aborting!'
            stop
          else
            visprec = .true.
          end if
        else if (arg.eq.'inviscid') then
          visprec = .false.
        else
          write(*,*) 'Unknown type of visc./inviscid prec. mode. Abortin &
     &g!'
          stop
        end if

        if (moving.and.(.not.unsteady)) then
          call linarg(4,line,arg)
          if (arg.eq.' ') then
            write(*,*) 'Missing absolute/relative frame choice for pre &
     &conditioning parameter. Aborting!'
            stop
          else
            read(arg,*) lsp_vfr
            if ((lsp_vfr.lt.0).or.(lsp_vfr.gt.1)) then
              write(*,*) 'Illegal lsp_vfr choice. Aborting!'
              stop
            end if
          end if
        end if

        if (unsteady) then
          call linarg(4,line,arg)
          if (arg.eq.' ') then
            write(*,*) 'Missing preconditioner type (standard/mixed). Ab &
     &orting!'
            stop
          else
            if (arg.eq.'standard') then
              mixpi = 0
            else if (arg.eq.'mixed') then
              mixpi = 1
            else
              write(*,*) 'Illegal preconditioner type. Aborting!'
              stop
            end if
          end if
          call linarg(5,line,arg)
          if (arg.eq.' ') then
            write(*,*) 'Missing x length scale for unst. LSP. Aborting!'
            stop
          else
            read(arg,*) lx
          end if
          call linarg(6,line,arg)
          if (arg.eq.' ') then
            write(*,*) 'Missing y length scale for unst. LSP. Aborting!'
            stop
          else
            read(arg,*) ly
          end if
!         define Unsteady Preconditioned Reference Velocity
          uprv = dmax1(lx/(pi*dt),ly/(pi*dt))
          if (moving) then
            call linarg(7,line,arg)
            if (arg.eq.' ') then
              write(*,*) 'Missing absolute/relative frame choice for pre &
     &conditioning parameter. Aborting!'
              stop
            else
              read(arg,*) lsp_vfr
              if ((lsp_vfr.lt.0).or.(lsp_vfr.gt.1)) then
                write(*,*) 'Illegal lsp_vfr choice. Aborting!'
                stop
              end if
            end if
          end if
        end if
        if ((etpfxtyp_f.eq.1).or.(etpfxtyp_t.eq.1)) then
          write(*,*) 'Inadmittable entropy fix type in combination with &
     &LSP. Aborting!'
          stop
        end if
      end if

!---- constant cut-off on multigrid levels
      if (lomach) then
        read(4,100)
        read(4,'(a)') line
        do nl=1,nlevel
          call linarg(nl,line,arg)
          read(arg,*) epsp(nl)
        end do
      end if

!---- tref, a Sutherland's law constant.
      if (viscous) then
        read(4,100)
        read(4,'(a)') line
        call linarg(1,line,arg)
        read(arg,*) tref
!       define and nondimensionalize Sutherland's temperature
        stemp = 1.104d+2
        stemp = stemp/tref
        if(amcontrol) then
          write(*,*) 'tref', tref
        end if
      end if

!---- dgn_cell: allow (y) or prevent degenerate cells at rotor axes
      read(4,100)
      read(4,'(a)') line
      call linarg(1,line,arg)
      if (arg.eq.'y') then
        dgn_cell = .true.
      else if (arg.eq.'n') then
        dgn_cell = .false.
      else
        write(*,*) 'Unknown value for degenerate cell option. Aborting!'
        stop
      end if
      call linarg(2,line,arg)
      if (arg.eq.'y') then
        wri_res = .true.
      else if (arg.eq.'n') then
        wri_res = .false.
      else
        write(*,*) 'Unknown value for restart output. Aborting!'
        stop
      end if
      call linarg(3,line,arg)
      if (arg.eq.'y') then
        wri_tec = .true.
      else if (arg.eq.'n') then
        wri_tec = .false.
      else
        write(*,*) 'Unknown value for tecplot output. Aborting!'
        stop
      end if

!-----read flow tec file format preference      
        read(4,100) line
        read(4,'(a)') line
        call linarg(1,line,arg)
        if (arg.eq.' ') then
           write(*,*) 'Missing Flow Tec output format.' 
           write(*,*) 'Aborting!'
           stop
        else
           if(arg.eq.'cgns') then
              write_cgns = .true.
           else if(arg.eq.'plt') then
              write_plt = .true.
           else if(arg.eq.'szplt') then
              write_szplt = .true.
           else if(arg.eq.'bin') then
              write_tec = .true.
           else
              write(*,*) 'Unknown flow tec output format'
              write(*,*) 'Should be either cgns, plt, szplt'
              write(*,*) 'or bin. Aborting.'
              stop
           end if
        end if

!---- read desired functional
      read(4,100)
      read(4,'(a)') line
      call linarg(1,line,arg)
      read(arg,'(a)') argo
      if (argo.eq.'default') then
        calfor=.true.
        call linarg(2,line,arg)
        read(arg,*) xmom
        call linarg(3,line,arg)
        read(arg,*) ymom
        call linarg(4,line,arg)
        read(arg,*) zmom

!------ three reference lengths
        read(4,100)
        read(4,'(a)') line
        do i=1,3
          call linarg(i,line,arg)
          read(arg,*) lngth(i)
        end do
        if(amcontrol) then
          write(*,*) 'lref1, lref2, lref3', (lngth(i),i=1,3)
        end if
                                
        if(.not.using_cgns) then
!------read number of bodies and block indices defining each body
           read(4,'(a)') line
           call linarg(1,line,arg)
           read(arg,*) nsurface
           if (nsurface.gt.msurface) then
              write(*,*) 'Exceeded max. number of bodies. Increase msurface and &
     & restart calculation. Aborting!'
              stop
           end if
           do ibody = 1,nsurface
              read(4,'(a)') line
              call linarg(1,line,arg)
              read(arg,*) n_surblo(ibody)
              if (n_surblo(ibody).gt.m_surblo) then
                 write(*,*) 'Exceeded max. number of wall blocks per body. In &
     &crease m_surblo and restart calculation. Aborting!'
                 stop
              end if
              do i_bodblo = 1,n_surblo(ibody)
                 call linarg(1+i_bodblo,line,arg)
                 read(arg,*) i_surbl(ibody,i_bodblo)
              end do
           end do
        end if
      else if (argo.eq.'none') then
        continue
      else
        write(*,*) 'unknown functional. Aborting!'
        stop
      end if

      if(.not.using_cgns) then
!---- read overall number of blocks
         read(4,*) nblocks

!---- read block size and boundary data
         do iblock = 1,nblocks
            
            read(4,'(a)') line
            call linarg(1,line,arg)
            read(arg,*) imax(iblock)
            call linarg(2,line,arg)
            read(arg,*) jmax(iblock)
            call linarg(3,line,arg)
            read(arg,*) kmax(iblock)
            call linarg(4,line,arg)
            read(arg,*) nbcs(iblock)
            if(amcontrol.and.debug) then
               write(*,*) 'block, imax, jmax, kmax, nbcs(iblock) : ', &
                    iblock,imax(iblock),jmax(iblock),kmax(iblock),nbcs(iblock)
            end if
            
            ijkmax(iblock) = max0(imax(iblock),jmax(iblock),kmax(iblock))
            
            do ibc = 1,nbcs(iblock)
               read(4,*) istrt(1),istrt(2),istrt(3),iend(1),iend(2),iend(3), &
                    bcdata(1,ibc,iblock)
               idir = 0
               do i = 1,3
                  if (istrt(i).eq.iend(i)) idir = i
               end do
               if (idir.lt.1.or.idir.gt.3) then
                  write(*,*) 'error in reading bc data. Aborting!'
                  stop
               end if
               do i = 1,3
                  if (i.ne.idir) then
                     iend(i) = iend(i) - 1
                  end if
               end do
               bcdata(2,ibc,iblock) = idir
               bcdata(3,ibc,iblock) = istrt(idir)
               bcdata(4,ibc,iblock) = istrt(1)
               bcdata(5,ibc,iblock) = iend (1)
               bcdata(6,ibc,iblock) = istrt(2)
               bcdata(7,ibc,iblock) = iend (2)
               bcdata(8,ibc,iblock) = istrt(3)
               bcdata(9,ibc,iblock) = iend (3)
               if(amcontrol.and.debug) then
                  write(*,*) 'block, ibc : ',iblock,ibc
                  write(*,*) bcdata(1,ibc,iblock), &
                       bcdata(2,ibc,iblock),bcdata(3,ibc,iblock), &
                       bcdata(4,ibc,iblock),bcdata(5,ibc,iblock), &
                       bcdata(6,ibc,iblock),bcdata(7,ibc,iblock), &
                       bcdata(8,ibc,iblock),bcdata(9,ibc,iblock)
               end if
               
               if (bcdata(1,ibc,iblock)/100.eq.15) then
                  if (((code.eq.'cosa').and.(kom.or.kom_bsl.or.kom_sst)).or. &
                       ((code.eq.'poscosa').and.viscous))     then
                  
                     nvwls(iblock) = nvwls(iblock) + 1
                     nwallv_1(nvwls(iblock),iblock) = 1
                     nwallv_2(nvwls(iblock),iblock) = 1
                     
                     ic2 = cyc(idir,2)
                     nwallv_1(nvwls(iblock),iblock) = &
                          nwallv_1(nvwls(iblock),iblock) * &
                          (iend(ic2) - istrt(ic2) + 2)
                     
                     ic2 = cyc(idir,3)
                     nwallv_2(nvwls(iblock),iblock) = &
                          nwallv_2(nvwls(iblock),iblock) * &
                          (iend(ic2) - istrt(ic2) + 2)
                     
                  end if
               end if
            end do
         end do
         
!---- grid cuts
         read(4,*) ncuts
         if(amcontrol) then
            write(*,*) 'ncuts = ',ncuts
         end if
         if (ncuts.gt.0) then
            do icut = 1,ncuts
               
               ic1 = 2 * icut - 1
               ic2 = 2 * icut
               read(4,*) iblk, (istrt(l),l=1,3),(iend(i),i=1,3)
               
!     Determine CUT face. This is defined by constant index on that
!     face. Note that this will fail for the pathological case where
!     have two or more indices equal for a cut.
               
               idir = 0
               do l = 1, 3
                  if (istrt(l) .eq. iend(l)) idir = l
               end do
               if (idir.lt.1.or.idir.gt.3) then
                  write(*,*) 'error in reading cut data (idir1). Aborting!'
                  stop
               end if
               
!     Modify start or end appropriately to account for cut input as
!     grid points but cut to be cell-centers. Do this only for the 
!     indices that the cut loops over, not the idir direction.
               
               do l = 1, 3
                  if (l .ne. idir) then
                     if (iend(l) .gt. istrt(l)) then
                        iend(l)  = iend(l)  - 1
                     else
                        istrt(l) = istrt(l) - 1
                     end if
                  end if
               end do
               
               cutdata( 1,ic1) = iblk
               cutdata( 2,ic1) = idir
               cutdata( 3,ic1) = istrt(idir)
               cutdata( 4,ic1) = istrt(1)
               cutdata( 5,ic1) = iend (1)
               cutdata( 6,ic1) = istrt(2)
               cutdata( 7,ic1) = iend (2)
               cutdata( 8,ic1) = istrt(3)
               cutdata( 9,ic1) = iend (3)
               
               cutdata(10,ic2) = iblk
               cutdata(11,ic2) = idir
               cutdata(12,ic2) = istrt(idir)
               cutdata(13,ic2) = istrt(1)
               cutdata(14,ic2) = iend (1)
               cutdata(15,ic2) = istrt(2)
               cutdata(16,ic2) = iend (2)
               cutdata(17,ic2) = istrt(3)
               cutdata(18,ic2) = iend (3)
               
               read(4,*) iblk, (istrt(l),l=1,3),(iend(i),i=1,3), &
                    (iord(n),n=1,3)

               do l = 1, 3
                  if (iord(l) .lt. 1 .or. iord(l) .gt. 3) then
                     write(*,*) 'error in reading cut data (iord). Aborting!'
                     stop
                  end if
               end do
               
!     Determine CUT face. This is defined by constant index on that
!     face. Note that this will fail for the pathological case where
!     have two or more indices equal for a cut.
               
               idir = 0
               do l = 1, 3
                  if (istrt(l) .eq. iend(l)) idir = l
               end do
               if (idir.lt.1.or.idir.gt.3) then
                  write(*,*) 'error in reading cut data (idir2). Aborting!'
                  stop
               end if
               
               
!     Modify start or end appropriately to account for cut input as
!     grid points but cut to be cell-centers. Do this only for the 
!     indices that the cut loops over, not the idir direction.
               
               do l = 1, 3
                  if (l .ne. idir) then
                     if (iend(l) .gt. istrt(l)) then
                        iend(l)  = iend(l)  - 1
                     else
                        istrt(l) = istrt(l) - 1
                     end if
                  end if
               end do
               
               cutdata(10,ic1) = iblk
               cutdata(11,ic1) = idir
               cutdata(12,ic1) = istrt(idir)
               cutdata(13,ic1) = istrt(1)
               cutdata(14,ic1) = iend (1)
               cutdata(15,ic1) = istrt(2)
               cutdata(16,ic1) = iend (2)
               cutdata(17,ic1) = istrt(3)
               cutdata(18,ic1) = iend (3)
               cutdata(19,ic1) = iord (1)
               cutdata(20,ic1) = iord (2)
               cutdata(21,ic1) = iord (3)
               
               cutdata( 1,ic2) = iblk
               cutdata( 2,ic2) = idir
               cutdata( 3,ic2) = istrt(idir)
               cutdata( 4,ic2) = istrt(1)
               cutdata( 5,ic2) = iend (1)
               cutdata( 6,ic2) = istrt(2)
               cutdata( 7,ic2) = iend (2)
               cutdata( 8,ic2) = istrt(3)
               cutdata( 9,ic2) = iend (3)
               cutdata(19,ic2) = cyc(iord(1),1)
               cutdata(20,ic2) = cyc(iord(2),1)
               cutdata(21,ic2) = cyc(iord(3),1)
               
               if ((iord(1).eq.2).and.(iord(2).eq.3)) then
                  cutdata(19,ic2) = cyc(iord(1),2)
                  cutdata(20,ic2) = cyc(iord(2),2)
                  cutdata(21,ic2) = cyc(iord(3),2)
               else if ((iord(1).eq.3).and.(iord(2).eq.1)) then
                  cutdata(19,ic2) = cyc(iord(1),3)
                  cutdata(20,ic2) = cyc(iord(2),3)
                  cutdata(21,ic2) = cyc(iord(3),3)
               else
                  cutdata(19,ic2) = cyc(iord(1),1)
                  cutdata(20,ic2) = cyc(iord(2),1)
                  cutdata(21,ic2) = cyc(iord(3),1)
               end if
               
            end do
            ncuts = ncuts * 2
         end if
         
!---- periodic cuts
         if (persec) then
            
            read(4,*) npcuts
            if (npcuts.le.0) then
               write(*,*) 'Missing periodic boundary definition. Abortinh!'
               stop
            end if
            if(amcontrol) then
               write(*,*) 'npcuts = ',npcuts
            end if
!     del    if (npcuts.gt.0) then
            do icut = 1,npcuts
               
               ic1 = 2 * icut - 1
               ic2 = 2 * icut
               read(4,*) iblk, (istrt(l),l=1,3),(iend(i),i=1,3)
               
!     Determine CUT face. This is defined by constant index on that
!     face. Note that this will fail for the pathological case where
!     have two or more indices equal for a cut.
               
               idir = 0
               do l = 1, 3
                  if (istrt(l) .eq. iend(l)) idir = l
               end do
               if (idir.lt.1.or.idir.gt.3) then
                  write(*,*) 'error in reading pcut data (idir1). Aborting!'
                  stop
               end if
               
!     Modify start or end appropriately to account for cut input as
!     grid points but cut to be cell-centers. Do this only for the 
!     indices that the cut loops over, not the idir direction.
               
               do l = 1, 3
                  if (l .ne. idir) then
                     if (iend(l) .gt. istrt(l)) then
                        iend(l)  = iend(l)  - 1
                     else
                        istrt(l) = istrt(l) - 1
                     end if
                  end if
               end do
               
               pcutdata( 1,ic1) = iblk
               pcutdata( 2,ic1) = idir
               pcutdata( 3,ic1) = istrt(idir)
               pcutdata( 4,ic1) = istrt(1)
               pcutdata( 5,ic1) = iend (1)
               pcutdata( 6,ic1) = istrt(2)
               pcutdata( 7,ic1) = iend (2)
               pcutdata( 8,ic1) = istrt(3)
               pcutdata( 9,ic1) = iend (3)
               
               pcutdata(10,ic2) = iblk
               pcutdata(11,ic2) = idir
               pcutdata(12,ic2) = istrt(idir)
               pcutdata(13,ic2) = istrt(1)
               pcutdata(14,ic2) = iend (1)
               pcutdata(15,ic2) = istrt(2)
               pcutdata(16,ic2) = iend (2)
               pcutdata(17,ic2) = istrt(3)
               pcutdata(18,ic2) = iend (3)
               
               read(4,*) iblk, (istrt(l),l=1,3),(iend(i),i=1,3), &
                    (iord(n),n=1,3),irot
               
               do l = 1, 3
                  if (iord(l) .lt. 1 .or. iord(l) .gt. 3) then
                     write(*,*) 'error in reading pcut data (iord). Aborting! &
     &'
                     stop
                  end if
               end do
               
!     Determine CUT face. This is defined by constant index on that
!     face. Note that this will fail for the pathological case where
!     have two or more indices equal for a cut.
               
               idir = 0
               do l = 1, 3
                  if (istrt(l) .eq. iend(l)) idir = l
               end do
               if (idir.lt.1.or.idir.gt.3) then
                  write(*,*) 'error in reading pcut data (idir2). Aborting!'
                  stop
               end if
               
               
!     Modify start or end appropriately to account for cut input as
!     grid points but cut to be cell-centers. Do this only for the 
!     indices that the cut loops over, not the idir direction.
               
               do l = 1, 3
                  if (l .ne. idir) then
                     if (iend(l) .gt. istrt(l)) then
                        iend(l)  = iend(l)  - 1
                     else
                        istrt(l) = istrt(l) - 1
                     end if
                  end if
               end do
               
               pcutdata(10,ic1) = iblk
               pcutdata(11,ic1) = idir
               pcutdata(12,ic1) = istrt(idir)
               pcutdata(13,ic1) = istrt(1)
               pcutdata(14,ic1) = iend (1)
               pcutdata(15,ic1) = istrt(2)
               pcutdata(16,ic1) = iend (2)
               pcutdata(17,ic1) = istrt(3)
               pcutdata(18,ic1) = iend (3)
               pcutdata(19,ic1) = iord (1)
               pcutdata(20,ic1) = iord (2)
               pcutdata(21,ic1) = iord (3)
               
               pcutdata( 1,ic2) = iblk
               pcutdata( 2,ic2) = idir
               pcutdata( 3,ic2) = istrt(idir)
               pcutdata( 4,ic2) = istrt(1)
               pcutdata( 5,ic2) = iend (1)
               pcutdata( 6,ic2) = istrt(2)
               pcutdata( 7,ic2) = iend (2)
               pcutdata( 8,ic2) = istrt(3)
               pcutdata( 9,ic2) = iend (3)
               pcutdata(19,ic2) = cyc(iord(1),1)
               pcutdata(20,ic2) = cyc(iord(2),1)
               pcutdata(21,ic2) = cyc(iord(3),1)
               
               if ((iord(1).eq.2).and.(iord(2).eq.3)) then
                  pcutdata(19,ic2) = cyc(iord(1),2)
                  pcutdata(20,ic2) = cyc(iord(2),2)
                  pcutdata(21,ic2) = cyc(iord(3),2)
               else if ((iord(1).eq.3).and.(iord(2).eq.1)) then
                  pcutdata(19,ic2) = cyc(iord(1),3)
                  pcutdata(20,ic2) = cyc(iord(2),3)
                  pcutdata(21,ic2) = cyc(iord(3),3)
               else
                  pcutdata(19,ic2) = cyc(iord(1),1)
                  pcutdata(20,ic2) = cyc(iord(2),1)
                  pcutdata(21,ic2) = cyc(iord(3),1)
               end if
               
               pcutdata(22,ic2)   =  irot
               pcutdata(22,ic1)   = -irot
               
            end do
            npcuts = npcuts * 2
!     del    end if
            
         end if
         
      end if
      
!---- various initializations
      do n=0,nharms
        dtheta(n)  = 0.d0
      end do

      if (irsop.gt.0) then
        rat  = vcfli(3) / vcfl(3)
        if (kom.or.kom_bsl.or.kom_sst) ratt = vcflit(3) / vcflt(3)
        if (.not.(ramping1.or.ramping2)) then
          cfl  = vcfli(3)
          if (kom.or.kom_bsl.or.kom_sst) cflt = vcflit(3)
        end if
      else if (irsop.eq.0) then
        if (.not.(ramping1.or.ramping2)) then
          cfl  = vcfl(3)
          if (kom.or.kom_bsl.or.kom_sst) cflt = vcflt(3)
        else
          write(*,*) 'CFL ramping cannot be used without IRS. Aborting!'
          stop
        end if
      endif

      if (persec) then
        somega = dsign(1.d0,omegas)
      end if

!---- test on sizes, compatibilities and availabilities

      if ((machfs.lt.1.d-14).and.viscous) then
        write(*,*) 'Zero freestream velocity problems cannot be solved d &
      ue to present lack of suitable scaling of viscous stress terms. Ab &
     &orting!'
        stop
      end if

      if (rgkuns.and.(nlevel.gt.1)) then
        write(*,*) 'Unsteady Runge-Kutta solver not allowed with MG. Abo &
     &rting!'
        stop
      end if

      if (rgkuns.and.lomach) then
        write(*,*) 'Unsteady Runge-Kutta solver not allowed with precond &
     &itioning. Aborting!'
        stop
      end if

      if ((kom.or.kom_bsl.or.kom_sst).and.(nlevel.gt.1)) then
        if (.not.bln_prol) then
          write(*,*) 'Bilinear interpolation only can be used with turbu &
     &lent flow analyses. Aborting!'
          stop
        end if
        if ((.not.prd_all_lvls).and.lim_prodtke.and. &
            (.not.lim_prodomega)) then
          write(*,*) 'Presently not possible to limit only k production &
      terms on coarser grids when restricting production from finer &
     &finer grid. Aborting!'
          stop
        end if
      end if

      if (irsop.gt.0) then
        if ((irsop.eq.5).or.(irsop.eq.6)) then
          if (lomach) then
            write(*,*) 'IRS for removal of diffusion limit not available &
     & for low-speed flows. Aborting!'
            stop
          end if
        end if
      end if

      if(.not. using_cgns) then

!---- set array dimensions on all grid levels

         call setdim(imax,jmax,kmax,ijkmax)
         if (((code.eq.'cosa').and.(kom.or.kom_bsl.or.kom_sst)).or. &
              ((code.eq.'poscosa').and.viscous)) then
            call setdim_vwall(nvwls,nwallv_1,nwallv_2)
            if(amcontrol.and.debug) then
               write(*,*) ' completed setdim_vwall'
            end if
         end if         

      end if

!---- process POSCOSA input
      if (code.eq.'poscosa') then

        read(4,100)

        read(4,'(a)') line
        call linarg(1,line,arg)
        read(arg,*) nblock_pos
        if ((nblock_pos.lt.1).and.(nblock_pos.gt.nblocks)) then
          write(*,*) 'Illegal number of blocks to be postprocessed. Abor &
     &ting!'
          stop
        else
          do iblock = 1,nblock_pos
            call linarg(1+iblock,line,arg)
            read(arg,*) iblock_pos(iblock)
          end do
        end if

        read(4,'(a)') line
        call linarg(1,line,arg)
        if  (arg.eq.'y') then
          wldata = .true.
          call linarg(2,line,arg)
          read(arg,*) wallmask
        else if (arg.eq.'n') then
          wldata = .false.
        else
          write(*,*)'Unknown option of WLDATA status output. Aborting!'
          stop
        end if

        read(4,'(a)') line
        call linarg(1,line,arg)
        if  (arg.eq.'y') then
          blprof = .true.
        else if (arg.eq.'n') then
          blprof = .false.
        else
          write(*,*) 'Unknown option of BL profile output. Aborting!'
          stop
        end if
        if (blprof) then
          call linarg(2,line,arg)
          read(arg,*) profpar(1)
          call linarg(3,line,arg)
          read(arg,*) profpar(2)
          call linarg(4,line,arg)
          read(arg,*) profpar(3)
          call linarg(5,line,arg)
          read(arg,*) profpar(4)
        end if

        read(4,'(a)') line
        call linarg(1,line,arg)
        if  (arg.eq.'y') then
          calvort = .true.
        else if (arg.eq.'n') then
          calvort = .false.
        else
          write(*,*) 'Unknown option of vorticity output. Aborting!'
          stop
        end if
        call linarg(2,line,arg)
        if  (arg.eq.'y') then
          calcp = .true.
        else if (arg.eq.'n') then
          calcp = .false.
        else
          write(*,*) 'Unknown option of cp output. Aborting!'
          stop
        end if

        read(4,'(a)') line
        call linarg(1,line,arg)
        if  (arg.eq.'y') then
          tecbin = .true.
        else if (arg.eq.'n') then
          tecbin = .false.
        else
          write(*,*) 'Unknown option of TECPLOT file format. Aborting!'
          stop
        end if

      end if

      close(4)


!      Get the input mesh data, including connectivity from the CGNS file
      if(using_cgns) then
         call get_cgns_data()
      end if

      if (nblocks.lt.1.or.nblocks.gt.mblock) then
        write(*,*) 'Illegal value of nblocks. Aborting!'
        stop
      end if
      
 100  format(20a4)

      return 
      end

!-----------------------------------------------------------------------
      subroutine setdim(imax,jmax,kmax,ijkmax)
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax(mblock), jmax(mblock), kmax(mblock), &
                ijkmax(mblock)
      integer(kind=cosa_int) nl,ibc,iblock

      if (moving.and.harbal.and.(.not.relframe)) then
        hbmove = 1
      else
        hbmove = 0
      end if

      if (hawt.and.rotating) then
        lmet = 8
      else if (moving) then
        lmet = 7
      else
        lmet = 4
      end if

      if (irsop.ge.1) then
        irspde = npde
      end if

      do nl=1,nlevel
        do iblock = 1,nblocks
          i_imax(iblock,nl)     = (imax (iblock)-1) / 2**(nl-1) + 1
          j_jmax(iblock,nl)     = (jmax (iblock)-1) / 2**(nl-1) + 1
          k_kmax(iblock,nl)     = (kmax (iblock)-1) / 2**(nl-1) + 1
          ijk_ijkmax(iblock,nl) = (ijkmax(iblock)-1) / 2**(nl-1) + 1
        end do
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine setdim_vwall(nwls,nwall_1,nwall_2)
!-----------------------------------------------------------------------
!     n_wall is the array which should become private to each MPI
!     process, since it's used to allocate the memory of the walls of
!     individual processes (xwall, ywall).
!     ng_wall is the array which should stay global, since it's used
!     to allocate the memory of the whole array of wall coordinates
!     which will be seen by all MPI processes.
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nwall_1(mwls,mblock),nwall_2(mwls,mblock),nwls(mblock), &
                nl,iblock,iwl

      do nl=1,nlevel

        do iblock = 1,nblocks
          n_wall(iblock,nl) = 0
        end do

        do iblock = 1,nblocks
          do iwl = 1,nwls(iblock)
            n_wall(iblock,nl) = n_wall(iblock,nl) + &
                                ((nwall_1(iwl,iblock)-1)/2**(nl-1)+1)* &
                                ((nwall_2(iwl,iblock)-1)/2**(nl-1)+1)
          end do
        end do

        ng_wall(nl) = 0
        do iblock = 1,nblocks
          ng_wall(nl) = ng_wall(nl) + n_wall(iblock,nl)
        end do

      end do

      return
      end

      
!-----------------------------------------------------------------------
      subroutine readgrid(nl,x,y,z,xyz,dx,dy,dz,dxyz)
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

#ifdef MPI
      include 'mpif.h'
      integer offset_kind
      parameter(offset_kind=MPI_OFFSET_KIND)
#else
      integer offset_kind
      parameter(offset_kind=8)
#endif

      integer ierr
      integer(kind=cosa_int) nl,imax,jmax,kmax,iblock,ixyz,jxyz,jdxyz,fid
      integer(kind=offset_kind) offset
      real(kind=cosa_real) x(*),y(*),z(*),xyz(*),dx(*),dy(*),dz(*),dxyz(*)
      logical parallel

      call usingparallel(parallel)


      if(parallel) then
         
        call parallelreadgrid(x,y,z,xyz,dx,dy,dz,dxyz,nl)
        
      else

         fid = 1
         open (unit=fid,file ='mesh.dat',status='old',form='unformatted',iostat=ierr)

         if(ierr .ne. 0) then
            write(*,*) 'Problem opening mesh.dat'
            stop
         end if
         
         do iblock = 1,mynblocks
            imax  = i_imax(iblock,nl)
            jmax  = j_jmax(iblock,nl)
            kmax  = k_kmax(iblock,nl)
            !ixyz  = 1 + off_0 (iblock,nl)
            ixyz = 1
            jxyz  = 1 + off_p2(iblock,nl)
            jdxyz = 1 + off_p1(iblock,nl)
            call read_bgrid(xyz(ixyz),dxyz(ixyz),imax,jmax,kmax,fid)
            call xyz_offset(x(jxyz),y(jxyz),z(jxyz),xyz(ixyz),dx(jdxyz), &
                 dy(jdxyz),dz(jdxyz),dxyz(ixyz),imax,jmax,kmax)
         end do

         close(unit=fid)

      end if

      return
      end

!-----------------------------------------------------------------------
      subroutine read_bgrid(xyz,dxyz,imax,jmax,kmax,fid)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax
      integer(kind=cosa_int) i,j,k,fid
      real(kind=cosa_real) &
          xyz (imax,jmax,kmax,3), &
          dxyz(imax,jmax,kmax,3)

      read(fid) (((xyz(i,j,k,1),i=1,imax),j=1,jmax),k=1,kmax), &
                (((xyz(i,j,k,2),i=1,imax),j=1,jmax),k=1,kmax), &
                (((xyz(i,j,k,3),i=1,imax),j=1,jmax),k=1,kmax)

      if (moving.and.pitching) then
        read(fid) (((dxyz(i,j,k,1),i=1,imax),j=1,jmax),k=1,kmax), &
                  (((dxyz(i,j,k,2),i=1,imax),j=1,jmax),k=1,kmax), &
                  (((dxyz(i,j,k,3),i=1,imax),j=1,jmax),k=1,kmax)
      end if

      return
      end

!-----------------------------------------------------------------------
      subroutine xyz_offset(x,y,z,xyz,dx,dy,dz,dxyz,imax,jmax,kmax)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax
      integer(kind=cosa_int) i,j,k
      real(kind=cosa_real) &
          x   (0:imax+1,0:jmax+1,0:kmax+1), &
          y   (0:imax+1,0:jmax+1,0:kmax+1), &
          z   (0:imax+1,0:jmax+1,0:kmax+1), &
          dx  (0:imax,0:jmax,0:kmax), &
          dy  (0:imax,0:jmax,0:kmax), &
          dz  (0:imax,0:jmax,0:kmax), &
          xyz (imax,jmax,kmax,3), &
          dxyz(imax,jmax,kmax,3)

      do k=1,kmax
        do j=1,jmax
          do i=1,imax
            x(i,j,k) = xyz(i,j,k,1)
            y(i,j,k) = xyz(i,j,k,2)
            z(i,j,k) = xyz(i,j,k,3)
          end do
        end do
      end do

      if (moving.and.pitching) then
        do k=1,kmax
          do j=1,jmax
            do i=1,imax
              dx(i,j,k) = dxyz(i,j,k,1)
              dy(i,j,k) = dxyz(i,j,k,2)
              dz(i,j,k) = dxyz(i,j,k,3)
            end do
          end do
        end do
      end if

      return
      end



!-----------------------------------------------------------------------
      subroutine readrest(q,q2,mut,nl)
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,nl,iblock,iq,imut,fid
      real (kind=cosa_real) q(*),q2(*),mut(*)
      logical parallel
     
      call usingparallel(parallel)

      if(parallel) then

         call parallelreadrestart(q,q2,mut,nl)

      else
         
         fid=23
         open(fid,file='restart',form='unformatted',status='old')
         
         do iblock = 1,mynblocks
            iq   = 1 + off_p3  (iblock,nl) * npde * dim5
            imax = i_imax  (iblock,nl)
            jmax = j_jmax  (iblock,nl)
            kmax = k_kmax  (iblock,nl)
            if (kom.or.kom_bsl.or.kom_sst) then
               imut   = 1 + off_p3 (iblock,nl) *        dim5
            else
               imut = 1
            end if
            call read_brest(q(iq),q2(iq),mut(imut),imax,jmax,kmax,npde,nharms,fid)
         end do
         
         close(unit=fid)

      end if

      return
      end

!-----------------------------------------------------------------------
      subroutine read_brest(q,q2,mut,imax,jmax,kmax,npde,nharms,fid)
!-----------------------------------------------------------------------
!     This functionality is replicated in the parallelutils.f
!     file as well to provide the functionality in MPI-I/O when
!     that is used (i.e. when we're running the MPI code).
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,npde,nharms
      integer(kind=cosa_int) i,imax1,j,jmax1,k,kmax1,ipde,n,fid,ios1
      real (kind=cosa_real) &
           q (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           q2(-1:imax+1,-1:jmax+1,-1:kmax+1,npde), &
           mut(-1:imax+1,-1:jmax+1,-1:kmax+1 ,0:2*nharms)

      imax1 = imax+1
      jmax1 = jmax+1
      kmax1 = kmax+1

      read(fid) (((((q(i,j,k,ipde,n), &
        i=-1,imax1),j=-1,jmax1),k=-1,kmax1),ipde=1,npde),n=0,2*nharms)

      ios1 = 0

      if (dual_prop_start) then
        read(fid,iostat=ios1) ((((q2(i,j,k,ipde), &
          i=-1,imax1),j=-1,jmax1),k=-1,kmax1),ipde=1,npde)
        if (ios1.lt.0) then
          write(*,*) 'Error in reading restart (q1). Aborting'
          stop
        end if
      end if

      ios1 = 0

      if (kom.or.kom_bsl.or.kom_sst) then
        read(fid,iostat=ios1) ((((mut(i,j,k,n), &
          i=-1,imax1),j=-1,jmax1),k=-1,kmax1),n=0,2*nharms)
        if (ios1.lt.0) then
          write(*,*) 'Error in reading restart (mut). Aborting'
          stop
        end if
      end if

      return
      end

!-----------------------------------------------------------------------
      subroutine init2fs(q,nl)
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,nl,iblock,iq
      real (kind=cosa_real) q(*)

      do iblock = 1,mynblocks
        imax    = i_imax  (iblock,nl)
        jmax    = j_jmax  (iblock,nl)
        kmax    = k_kmax  (iblock,nl)
        iq      = 1 + off_p3(iblock,nl) * npde * dim5
        call init2fs_b(q(iq),imax,jmax,kmax,npde,nharms)
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine init2fs_b(q,imax,jmax,kmax,npde,nharms)
!-----------------------------------------------------------------------
!---- uniform flow initialization for external flow
!---- msc, 21/03/2010: in principle, the total energy should also
!                      include tkefar for turbulent flows. Here this
!                      fact is being neglected.
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax,npde,nharms
      integer(kind=cosa_int) i,j,k,n
      real (kind=cosa_real) q(-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms)
      real (kind=cosa_real) angle, rotaz(2,2), ua, va

      if (.not.hawt) then
        do n = 0,2*nharms
          do k = -1,kmax+1
            do j = -1,jmax+1
              do i = -1,imax+1
                q(i,j,k,1,n) = rhoinf
                q(i,j,k,2,n) = rhoinf * machfs * dcos(alpha) * dcos(beta)
                q(i,j,k,3,n) = rhoinf * machfs * dsin(alpha) * dcos(beta)
                q(i,j,k,4,n) = rhoinf * machfs * dsin(beta)
                q(i,j,k,5,n) = &
                  1/rhoinf/gamma/(gamma-1) + 0.5d0 * machfs**2
              end do
            end do
          end do
        end do

      else

        if (dualt.or.rgkuns.or.(.not.relframe)) then

!         In dualt and rgkuns, both rel and abs frames can be used.
!         In steady and HB cases, one can use absolute frame.
          do n = 0,2*nharms
            do k = -1,kmax+1
              do j = -1,jmax+1
                do i = -1,imax+1
                  q(i,j,k,1,n) = rhoinf
                  q(i,j,k,2,n) = rhoinf * machfs * dsin(beta)
                  q(i,j,k,3,n) = machfs * dsin(alpha) * dcos(beta)
                  q(i,j,k,4,n) = machfs * dcos(alpha) * dcos(beta)
                  q(i,j,k,5,n) = &
                    1/rhoinf/gamma/(gamma-1) + 0.5d0 * machfs**2
                end do
              end do
            end do
          end do

        else

!         This is for steady and HB calcs using rel frame.
          ua =  machfs * dsin(beta)
          va =  machfs * dsin(alpha) * dcos(beta)

          do n = 0,2*nharms
            angle      = -omegas * zeit(n)
            rotaz(1,1) =  cos(angle)
            rotaz(1,2) = -sin(angle)
            rotaz(2,1) =  sin(angle)
            rotaz(2,2) =  cos(angle)
            do k = -1,kmax+1
              do j = -1,jmax+1
                do i = -1,imax+1
                  q(i,j,k,1,n) = rhoinf
                  q(i,j,k,2,n) = rhoinf * &
                    (ua*rotaz(1,1) + va*rotaz(1,2))
                  q(i,j,k,3,n) = rhoinf * &
                    (ua*rotaz(2,1) + va*rotaz(2,2))
                  q(i,j,k,4,n) = rhoinf * machfs * dcos(alpha) * dcos(beta)
                  q(i,j,k,5,n) = &
                    1/rhoinf/gamma/(gamma-1) + 0.5d0 * machfs**2
                end do
              end do
            end do
          end do

        end if

      end if

      if (kom.or.kom_bsl.or.kom_sst) then
        do n = 0,2*nharms
          do k = 1,kmax-1
            do j = 1,jmax-1
              do i = 1,imax-1
                q(i,j,k,6,n) = rhoinf * tkefar
                q(i,j,k,7,n) = rhoinf * omefar
              end do
            end do
          end do
        end do
      end if

      return
      end

!-----------------------------------------------------------------------
      subroutine linarg(n,line,arg)
!-----------------------------------------------------------------------
!
      use cosa_precision

      implicit none
      character*40 arg
      character*10000 line
      character*10001 lin2
      integer(kind=cosa_int)    n, m, i1, i2
!
      lin2 = line//'!'
      i2   = 1
!
      do m = 1, n
        i1 = i2
!       detect beginning of m-th word
        do while(lin2(i1:i1).eq.' ' .and. lin2(i1+1:i1+1).ne.'!')
          i1 = i1+1
        enddo
        i2 = i1
!       detect end of m-th word
        do while(lin2(i2:i2).ne.' ' .and. lin2(i1+1:i1+1).ne.'!')
          i2 = i2+1
        enddo
      enddo

      arg = line(i1:i2)

      return
      end

!-----------------------------------------------------------------------
      subroutine alloc()
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

#ifdef MPI
      include 'mpif.h'
      integer offset_kind
      parameter(offset_kind=MPI_OFFSET_KIND)
#else
      integer offset_kind
      parameter(offset_kind=8)
#endif


      integer(kind=cosa_int) imax, jmax, kmax, ijkmax
      integer(kind=cosa_int) nl,dim4,offset,pop
      integer(kind=offset_kind) dimg

      if (harbal) then
        dim5  = 2*nharms+1
        dim5h = dim5**hbmove
      else
        dim5  = 1
        dim5h = 1
      end if

      pop = 0

      do nl=1,nlevel

        offset = 3
        call size_grid(dimg,nl,offset)
        call cosa_all(p_qp(nl)   ,'r8',dimg,npde,dim5 ,pop)
        call cosa_all(p_qm(nl)   ,'r8',dimg,npde,dim5 ,pop)
        call cosa_all(p_dqp(nl)  ,'r8',dimg,npde,dim5 ,pop)
        call cosa_all(p_dqm(nl)  ,'r8',dimg,npde,dim5 ,pop)
        if (code.eq.'cosa') then
          call cosa_all(p_pr(nl)   ,'r8',dimg,   1,dim5 ,pop)
          call cosa_all(p_q(nl)    ,'r8',dimg,npde,dim5 ,pop)
          call cosa_all(p_qold(nl) ,'r8',dimg,npde,dim5 ,pop)
          call cosa_all(p_res(nl)  ,'r8',dimg,npde,dim5 ,pop)
          call cosa_all(p_rhs(nl)  ,'r8',dimg,npde,dim5 ,pop)
          call cosa_all(p_dq(nl)   ,'r8',dimg,npde,dim5 ,pop)
          call cosa_all(p_work1(nl),'r8',dimg,npde,dim5,pop)
          if (dualt) then
            call cosa_all(p_q1(nl)   ,'r8',dimg,npde,   1 ,pop)
            call cosa_all(p_q2(nl)   ,'r8',dimg,npde,   1 ,pop)
          end if
          if (kom.or.kom_bsl.or.kom_sst) then
            call cosa_all(p_mut(nl)   ,'r8',dimg,1,dim5,pop)
            call cosa_all(p_divvel(nl),'r8',dimg,1,dim5,pop)
            call cosa_all(p_delpls(nl),'r8',dimg,1,dim5,pop)
            call cosa_all(p_prd(nl)   ,'r8',dimg,1,dim5,pop)
            call cosa_all(p_ttrm(nl)  ,'r8',dimg,3,dim5,pop)
          end if
        end if

        offset = 2
        call size_grid(dimg,nl,offset)
        call cosa_all(p_x(nl)  ,'r8',dimg,1   ,dim5h,pop)
        call cosa_all(p_y(nl)  ,'r8',dimg,1   ,dim5h,pop)
        call cosa_all(p_z(nl)  ,'r8',dimg,1   ,dim5h,pop)
        call cosa_all(p_xc(nl) ,'r8',dimg,1   ,dim5h,pop)
        call cosa_all(p_yc(nl) ,'r8',dimg,1   ,dim5h,pop)
        call cosa_all(p_zc(nl) ,'r8',dimg,1   ,dim5h,pop)
        call cosa_all(p_si(nl), 'r8',dimg,lmet,dim5h,pop)
        call cosa_all(p_sj(nl), 'r8',dimg,lmet,dim5h,pop)
        call cosa_all(p_sk(nl), 'r8',dimg,lmet,dim5h,pop)
        if (moving) then
          call cosa_all(p_x0(nl),   'r8',dimg,1,1    ,pop)
          call cosa_all(p_y0(nl),   'r8',dimg,1,1    ,pop)
          call cosa_all(p_z0(nl),   'r8',dimg,1,1    ,pop)
          call cosa_all(p_xdot(nl), 'r8',dimg,1,dim5h,pop)
          call cosa_all(p_ydot(nl), 'r8',dimg,1,dim5h,pop)
          call cosa_all(p_zdot(nl), 'r8',dimg,1,dim5h,pop)
        end if

        offset = 1
        call size_grid(dimg,nl,offset)
        call cosa_all(p_vol(nl),'r8',dimg,1,1    ,pop)
        if (moving.and.pitching) then
          call cosa_all(p_dx(nl),'r8',dimg,1,1,pop)
          call cosa_all(p_dy(nl),'r8',dimg,1,1,pop)
          call cosa_all(p_dz(nl),'r8',dimg,1,1,pop)
        end if
        if (code.eq.'cosa') then
          call cosa_all(p_dtvol(nl),'r8',dimg,1,dim5,pop)
          if (irsop.ge.2) then
            call cosa_all(p_betaxi(nl)  ,'r8',dimg,1,dim5,pop)
            call cosa_all(p_betaeta(nl) ,'r8',dimg,1,dim5,pop)
            call cosa_all(p_betazeta(nl),'r8',dimg,1,dim5,pop)
          end if
          if (kom.or.kom_bsl.or.kom_sst) then
            call cosa_all(p_dist(nl)  ,'r8',dimg,1,1   ,pop)
            call cosa_all(p_dtvolt(nl),'r8',dimg,1,dim5,pop)
            if (irsop.ge.2) then
              call cosa_all(p_tbetaxi(nl)  ,'r8',dimg,1,dim5,pop)
              call cosa_all(p_tbetaeta(nl) ,'r8',dimg,1,dim5,pop)
              call cosa_all(p_tbetazeta(nl),'r8',dimg,1,dim5,pop)
            end if
          end if
        end if
        if ((code.eq.'poscosa').and.viscous) then
          call cosa_all(p_dist(nl),'r8',dimg,1,1,pop)
        end if

        offset = 0
        call size_grid(dimg,nl,offset)
        if (viscous) then
          call cosa_all(p_xideri(nl)  ,'r8',dimg,3,dim5h,pop)
          call cosa_all(p_xiderj(nl)  ,'r8',dimg,3,dim5h,pop)
          call cosa_all(p_xiderk(nl)  ,'r8',dimg,3,dim5h,pop)
          call cosa_all(p_etaderi(nl) ,'r8',dimg,3,dim5h,pop)
          call cosa_all(p_etaderj(nl) ,'r8',dimg,3,dim5h,pop)
          call cosa_all(p_etaderk(nl) ,'r8',dimg,3,dim5h,pop)
          call cosa_all(p_zetaderi(nl),'r8',dimg,3,dim5h,pop)
          call cosa_all(p_zetaderj(nl),'r8',dimg,3,dim5h,pop)
          call cosa_all(p_zetaderk(nl),'r8',dimg,3,dim5h,pop)
          if (kom_bsl.or.kom_sst) then
            call cosa_all(p_fidist(nl),'r8',dimg,1,1,pop)
            call cosa_all(p_fjdist(nl),'r8',dimg,1,1,pop)
            call cosa_all(p_fkdist(nl),'r8',dimg,1,1,pop)
          end if
        end if

        if (code.eq.'cosa') then
          call cosa_all(p_flux(nl), 'r8',dimg,npde,dim5,pop)
        end if

        offset = -1
        call size_grid(dimg,nl,offset)
        if (1.lt.0) then
          call cosa_all(p_rad(nl),  'r8',dimg,1,dim5h,pop)
        end if

        if ((code.eq.'cosa').and.lomach) then
          dim4       = npde * npde
          call cosa_all(p_prec(nl)      ,'r8',dimg,dim4,dim5,pop)
          call cosa_all(p_cutoff_pgr(nl),'r8',dimg,1   ,dim5,pop)
          call cosa_all(p_cutoff_vis(nl),'r8',dimg,1   ,dim5,pop)
          if (rkimp.or.kom.or.kom_bsl.or.kom_sst) then
            call cosa_all(p_ipkpi(nl),'r8',dimg,dim4,dim5,pop)
          end if
          if (debug) then
            call cosa_all(p_preci(nl),'r8',dimg,dim4,dim5,pop)
          end if
        end if

        offset = -1
        call size_grid_1d(dimg,offset,mblock,ijk_ijkmax(1,nl),nblocks)
        if (irsop.ge.1) then
          dim4       = irspde
          call cosa_all(p_work2(nl),'r8',dimg,dim4,1,pop)
        end if

        call size_bnds(dimg,nl)
        call cosa_all(p_bctopo(nl),'i4',dimg,10,1,pop)

        dimg       = 21 * myncuts
        dimg       = dimg * mynblocks
        call cosa_all(p_cutopo(nl),'i4',dimg,1,1,pop)

        if (persec) then
          dimg       = 22 * mynpcuts
          dimg       = dimg * mynblocks
          call cosa_all(p_percutopo(nl),'i4',dimg,1,1,pop)
        end if

        if (((code.eq.'cosa').and.(kom.or.kom_bsl.or.kom_sst)).or. &
            ((code.eq.'poscosa').and.viscous)) then
          offset = 0
!         local memory allocation for each MPI process
          call size_grid_1d(dimg,offset,mblock,n_wall(1,nl),nblocks)
          call cosa_all(p_xwall(nl),'r8',dimg,   1,1,pop)
          call cosa_all(p_ywall(nl),'r8',dimg,   1,1,pop)
          call cosa_all(p_zwall(nl),'r8',dimg,   1,1,pop)
!         global memory allocation for complete array of wall nodes
          dimg = ng_wall(nl)
          call cosa_all(p_xgwall(nl),'r8',dimg,  1,1,pop)
          call cosa_all(p_ygwall(nl),'r8',dimg,  1,1,pop)
          call cosa_all(p_zgwall(nl),'r8',dimg,  1,1,pop)
        end if

      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine cosa_all(ptr,type,dimg,npde,nharms,pop)
!-----------------------------------------------------------------------

      use cosa_precision

      implicit none

#ifdef MPI
      include 'mpif.h'
      integer offset_kind
      parameter(offset_kind=MPI_OFFSET_KIND)
#else
      integer offset_kind
      parameter(offset_kind=8)
#endif

      integer(kind=cosa_int) npde,nharms,pop
      integer(kind=offset_kind) dimg,dim
      integer(kind=offset_kind) ptr,malloc
      character*2 type

      real (kind=cosa_real) r8
      pointer (pr8, r8)

      dim = dimg*npde*nharms

#ifdef XL
      if (type.eq.'r8') then
        ptr = malloc(%val(8*dim))
      else if (type.eq.'i4') then
        ptr = malloc(%val(4*dim))
      end if
#else
      if (type.eq.'r8') then
        ptr = malloc(8*dim)
      else if (type.eq.'i4') then
        ptr = malloc(4*dim)
      end if
#endif

      if(ptr .eq. 0 .and. dim .ne. 0) then
        write(*,*) 'Error, allocation has not worked, stopping',dim,ptr
#if MPI
        call abortmpi()
#else
        stop
#endif
      end if

      if (type.eq.'r8') then
        pr8 = ptr
        call zeroall(dim,r8)
      end if

!dbg            if (pop.eq.1) then
!dbg                  write(*,*) 'from all', dim1,dim2,dim3,dim4,dim,ptr
!dbg            end if

      return
      end

!-----------------------------------------------------------------------
      subroutine setdim_bc(bctopo,nl)
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,iblock,iblk,ibc
      integer(kind=cosa_int) bctopo(*)

      do iblock = 1,mynblocks
        iblk = 1 + off_bct(iblock,nl)
        call setdim_bbc(bctopo(iblk),iblock,nl)
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine setdim_bbc(bctopo,iblock,nl)
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,ibc,iblock,idbg
      integer(kind=cosa_int) bctopo(10,nbcs(iblock))

      logical amcontrol

      call amcontroller(amcontrol)

      do ibc=1,nbcs(iblock)
        bctopo(1,ibc) = bcdata(1,ibc,iblock)
        bctopo(2,ibc) = bcdata(2,ibc,iblock)
        if (nl.eq.1) then
          bctopo(3,ibc) = bcdata(3,ibc,iblock)
          bctopo(4,ibc) = bcdata(4,ibc,iblock)
          bctopo(5,ibc) = bcdata(5,ibc,iblock)
          bctopo(6,ibc) = bcdata(6,ibc,iblock)
          bctopo(7,ibc) = bcdata(7,ibc,iblock)
          bctopo(8,ibc) = bcdata(8,ibc,iblock)
          bctopo(9,ibc) = bcdata(9,ibc,iblock)
        else if (nl.gt.1) then
          bctopo(3,ibc) = bcdata(3,ibc,iblock) / 2**(nl-1) + 1
          bctopo(4,ibc) = (bcdata(4,ibc,iblock)-1) / 2**(nl-1) + 1
          bctopo(5,ibc) = (bcdata(5,ibc,iblock)-1) / 2**(nl-1) + 1
          bctopo(6,ibc) = (bcdata(6,ibc,iblock)-1) / 2**(nl-1) + 1
          bctopo(7,ibc) = (bcdata(7,ibc,iblock)-1) / 2**(nl-1) + 1
          bctopo(8,ibc) = (bcdata(8,ibc,iblock)-1) / 2**(nl-1) + 1
          bctopo(9,ibc) = (bcdata(9,ibc,iblock)-1) / 2**(nl-1) + 1
        end if
        if(amcontrol.and.debug) then
          write(*,*)
          write(*,*) 'From setdimbc for nl, iblock', nl, iblock
          write(*,*) 'ibc,nbcs(iblock)',ibc,nbcs(iblock)
          write(*,*) '(bctopo(idbg,ibc),idbg=1,7)', &
                      (bctopo(idbg,ibc),idbg=1,7)
          write(*,*)
        end if
      end do


      return
      end

!-----------------------------------------------------------------------
      subroutine mathb(nharms)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nharms
      integer(kind=cosa_int) i,j,k
      real(kind=cosa_real) fact

      do i=0,2*nharms
        do j=0,2*nharms
          ahb(i,j) = 0.d0
          dhb(i,j) = 0.d0
        end do
      end do

!---- matrix A, eq. 9, Liu et al., JCP 2006
      do i=1,nharms
        ahb( 2*(i-1)+1,2*(i-1)+2 ) =  i
        ahb( 2*(i-1)+2,2*(i-1)+1 ) = -i
      end do

!---- matrix E^(-1), eq. 13, Liu et al., JCP 2006
      do i=0,2*nharms
        eihb(i,0) = 1.d0
      end do
      do j=1,nharms
        do i=0,2*nharms
          eihb(i,2*(j-1)+1) = dcos(j*omega*zeit(i))
          eihb(i,2*(j-1)+2) = dsin(j*omega*zeit(i))
        end do
      end do

!---- matrix E, eq. 12, Liu et al., JCP 2006
      fact = 2.d0/(2*nharms+1)
      do j=0,2*nharms
        ehb(0,j) = fact/2
      end do
      do i=1,nharms
        do j=0,2*nharms
          ehb(2*(i-1)+1,j) = fact * dcos(i*omega*zeit(j))
          ehb(2*(i-1)+2,j) = fact * dsin(i*omega*zeit(j))
        end do
      end do

!---- compute dhb = [ eihb * ahb * ehb ] = [ E^(-1) * A * E ] 
      do i=0,2*nharms
        do j=0,2*nharms
          do k=0,2*nharms
            dhb(i,j) = dhb(i,j) + ahb(i,k) * ehb(k,j)
          end do
        end do
      end do
      do i=0,2*nharms
        do j=0,2*nharms
          ahb(i,j) = dhb(i,j)
          dhb(i,j) = 0.d0
        end do
      end do
      do i=0,2*nharms
        do j=0,2*nharms
          do k=0,2*nharms
            dhb(i,j) = dhb(i,j) + eihb(i,k) * ahb(k,j)
          end do
        end do
      end do

!dbg  write(*,*) 
!dbg  do i=0,2*nharms
!dbg    write(*,*) (dhb(i,j),j=0,2*nharms)
!dbg  end do

      return
      end

!-----------------------------------------------------------------------
      subroutine setdimcut(cutopo,percutopo,nl)
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) cutopo(21,ncuts),percutopo(22,npcuts)
      integer(kind=cosa_int) nl,icut,idbg
      logical amcontrol

      call amcontroller(amcontrol)

      do icut=1,myncuts
        cutopo( 1,icut) = cutdata( 1,icut)
        cutopo( 2,icut) = cutdata( 2,icut)
        cutopo(10,icut) = cutdata(10,icut)
        cutopo(11,icut) = cutdata(11,icut)
        cutopo(12,icut) = cutdata(12,icut)
        cutopo(19,icut) = cutdata(19,icut)
        cutopo(20,icut) = cutdata(20,icut)
        cutopo(21,icut) = cutdata(21,icut)
        if (nl.eq.1) then
          cutopo( 3,icut) = cutdata( 3,icut)
          cutopo( 4,icut) = cutdata( 4,icut)
          cutopo( 5,icut) = cutdata( 5,icut)
          cutopo( 6,icut) = cutdata( 6,icut)
          cutopo( 7,icut) = cutdata( 7,icut)
          cutopo( 8,icut) = cutdata( 8,icut)
          cutopo( 9,icut) = cutdata( 9,icut)
          cutopo(12,icut) = cutdata(12,icut)
          cutopo(13,icut) = cutdata(13,icut)
          cutopo(14,icut) = cutdata(14,icut)
          cutopo(15,icut) = cutdata(15,icut)
          cutopo(16,icut) = cutdata(16,icut)
          cutopo(17,icut) = cutdata(17,icut)
          cutopo(18,icut) = cutdata(18,icut)
        else if (nl.gt.1) then
          cutopo( 3,icut) =  cutdata( 3,icut)    / 2**(nl-1) + 1
          cutopo( 4,icut) = (cutdata( 4,icut)-1) / 2**(nl-1) + 1
          cutopo( 5,icut) = (cutdata( 5,icut)-1) / 2**(nl-1) + 1
          cutopo( 6,icut) = (cutdata( 6,icut)-1) / 2**(nl-1) + 1
          cutopo( 7,icut) = (cutdata( 7,icut)-1) / 2**(nl-1) + 1
          cutopo( 8,icut) = (cutdata( 8,icut)-1) / 2**(nl-1) + 1
          cutopo( 9,icut) = (cutdata( 9,icut)-1) / 2**(nl-1) + 1
          cutopo(12,icut) =  cutdata(12,icut)    / 2**(nl-1) + 1
          cutopo(13,icut) = (cutdata(13,icut)-1) / 2**(nl-1) + 1
          cutopo(14,icut) = (cutdata(14,icut)-1) / 2**(nl-1) + 1
          cutopo(15,icut) = (cutdata(15,icut)-1) / 2**(nl-1) + 1
          cutopo(16,icut) = (cutdata(16,icut)-1) / 2**(nl-1) + 1
          cutopo(17,icut) = (cutdata(17,icut)-1) / 2**(nl-1) + 1
          cutopo(18,icut) = (cutdata(18,icut)-1) / 2**(nl-1) + 1
        end if
!dbg    if (amcontrol) then
!dbg      write(*,*) 'nl,icut,(cutopo(idbg,icut),idbg=1,21)', &
!dbg                  nl,icut,(cutopo(idbg,icut),idbg=1,21)
!dbg    end if
      end do

!dbg  if(amcontrol) then
!dbg    write(*,*)
!dbg  end if

      if(persec) then
         do icut=1,mynpcuts
            percutopo( 1,icut) = pcutdata( 1,icut)
            percutopo( 2,icut) = pcutdata( 2,icut)
            percutopo(10,icut) = pcutdata(10,icut)
            percutopo(11,icut) = pcutdata(11,icut)
            percutopo(12,icut) = pcutdata(12,icut)
            percutopo(19,icut) = pcutdata(19,icut)
            percutopo(20,icut) = pcutdata(20,icut)
            percutopo(21,icut) = pcutdata(21,icut)
            percutopo(22,icut) = pcutdata(22,icut)
            if (nl.eq.1) then
               percutopo( 3,icut) = pcutdata( 3,icut)
               percutopo( 4,icut) = pcutdata( 4,icut)
               percutopo( 5,icut) = pcutdata( 5,icut)
               percutopo( 6,icut) = pcutdata( 6,icut)
               percutopo( 7,icut) = pcutdata( 7,icut)
               percutopo( 8,icut) = pcutdata( 8,icut)
               percutopo( 9,icut) = pcutdata( 9,icut)
               percutopo(12,icut) = pcutdata(12,icut)
               percutopo(13,icut) = pcutdata(13,icut)
               percutopo(14,icut) = pcutdata(14,icut)
               percutopo(15,icut) = pcutdata(15,icut)
               percutopo(16,icut) = pcutdata(16,icut)
               percutopo(17,icut) = pcutdata(17,icut)
               percutopo(18,icut) = pcutdata(18,icut)
            else if (nl.gt.1) then
               percutopo( 3,icut) =  pcutdata( 3,icut)    / 2**(nl-1) + 1
               percutopo( 4,icut) = (pcutdata( 4,icut)-1) / 2**(nl-1) + 1
               percutopo( 5,icut) = (pcutdata( 5,icut)-1) / 2**(nl-1) + 1
               percutopo( 6,icut) = (pcutdata( 6,icut)-1) / 2**(nl-1) + 1
               percutopo( 7,icut) = (pcutdata( 7,icut)-1) / 2**(nl-1) + 1
               percutopo( 8,icut) = (pcutdata( 8,icut)-1) / 2**(nl-1) + 1
               percutopo( 9,icut) = (pcutdata( 9,icut)-1) / 2**(nl-1) + 1
               percutopo(12,icut) =  pcutdata(12,icut)    / 2**(nl-1) + 1
               percutopo(13,icut) = (pcutdata(13,icut)-1) / 2**(nl-1) + 1
               percutopo(14,icut) = (pcutdata(14,icut)-1) / 2**(nl-1) + 1
               percutopo(15,icut) = (pcutdata(15,icut)-1) / 2**(nl-1) + 1
               percutopo(16,icut) = (pcutdata(16,icut)-1) / 2**(nl-1) + 1
               percutopo(17,icut) = (pcutdata(17,icut)-1) / 2**(nl-1) + 1
               percutopo(18,icut) = (pcutdata(18,icut)-1) / 2**(nl-1) + 1
            end if
!     dbg    if (amcontrol) then
!     dbg      write(*,*) 'nl,icut,(cutopo(idbg,icut),idbg=1,21)',
!     dbg &                nl,icut,(cutopo(idbg,icut),idbg=1,21)
!     dbg    end if
         end do
      end if
      

!dbg  if(amcontrol) then
!dbg    write(*,*)
!dbg  end if

      return
      end

!-----------------------------------------------------------------------
      subroutine size_grid(dimg,nl,offset)
!-----------------------------------------------------------------------

      use cosa_variables
      use cosa_precision

      implicit none

#ifdef MPI
      include 'mpif.h'
      integer offset_kind
      parameter(offset_kind=MPI_OFFSET_KIND)
#else
      integer offset_kind
      parameter(offset_kind=8)
#endif


      integer(kind=cosa_int) offset,iblock,nl
      integer(kind=offset_kind) dimg

      dimg = 0

      do iblock = 1,mynblocks
        dimg = dimg + (i_imax(iblock,nl)+offset) * &
                      (j_jmax(iblock,nl)+offset) * &
                      (k_kmax(iblock,nl)+offset)
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine size_grid_1d(dimg,offset,mblock,dim,nblocks)
!-----------------------------------------------------------------------

      use cosa_precision

      implicit none

#ifdef MPI
      include 'mpif.h'
      integer offset_kind
      parameter(offset_kind=MPI_OFFSET_KIND)
#else
      integer offset_kind
      parameter(offset_kind=8)
#endif

      integer(kind=cosa_int) mblock,dim(mblock)
      integer(kind=cosa_int) nblocks,offset,iblock
      integer(kind=offset_kind) dimg

      dimg = 0
      
      do iblock = 1,nblocks
        dimg = dimg + (dim(iblock)+offset)
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine size_bnds(dimg,nl)
!-----------------------------------------------------------------------

      use cosa_variables
      use cosa_precision

      implicit none

#ifdef MPI
      include 'mpif.h'
      integer offset_kind
      parameter(offset_kind=MPI_OFFSET_KIND)
#else
      integer offset_kind
      parameter(offset_kind=8)
#endif


      integer(kind=cosa_int) iblock,nl
      integer(kind=offset_kind) dimg

      dimg = 0

      do iblock = 1,nblocks
        dimg = dimg + nbcs(iblock)
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine setoffs(nl)
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,iblock

      off_1d    (1,nl) = 0
      off_1dw   (1,nl) = 0
      off_1dwg  (1,nl) = 0
      off_bct   (1,nl) = 0
      off_m1    (1,nl) = 0
      off_0     (1,nl) = 0
      off_p1    (1,nl) = 0
      off_p2    (1,nl) = 0
      off_p3    (1,nl) = 0

!mpi  This is purposefully left like this for the MPI code
!mpi  It is used to calculate the offsets into the xgwall, ygwall
!mpi  and zgwall arrays.
      do iblock=2,nblocks
        off_1dwg(iblock,nl) = off_1dwg(iblock-1,nl) + &
             g_n_wall(iblock-1,nl)

      end do

      if(mynblocks .gt. 1) then

        do iblock=2,mynblocks
          off_1d (iblock,nl) = off_1d (iblock-1,nl) + &
                               (ijk_ijkmax(iblock-1,nl)-1)
        end do

        do iblock=2,mynblocks
          off_1dw(iblock,nl) = off_1dw(iblock-1,nl) + n_wall(iblock-1,nl)
!dng      write(*,*) 'off_1dw(iblock,nl)', off_1dw(iblock,nl)
        end do

        do iblock=2,mynblocks
          off_bct(iblock,nl) = &
           off_bct(iblock-1,nl) + nbcs(iblock-1) * 10
        end do

        do iblock=2,mynblocks
          off_m1(iblock,nl) = off_m1(iblock-1,nl) + &
            (i_imax(iblock-1,nl)-1) * (j_jmax(iblock-1,nl)-1) * &
            (k_kmax(iblock-1,nl)-1)
        end do

        do iblock=2,mynblocks
          off_0(iblock,nl) = off_0(iblock-1,nl) + &
            (i_imax(iblock-1,nl)  ) * (j_jmax(iblock-1,nl)  ) * &
            (k_kmax(iblock-1,nl)  )
        end do

        do iblock=2,mynblocks
          off_p1(iblock,nl) = off_p1(iblock-1,nl) + &
            (i_imax(iblock-1,nl)+1) * (j_jmax(iblock-1,nl)+1) * &
            (k_kmax(iblock-1,nl)+1)
        end do

        do iblock=2,mynblocks
          off_p2(iblock,nl) = off_p2(iblock-1,nl) + &
            (i_imax(iblock-1,nl)+2) * (j_jmax(iblock-1,nl)+2) * &
            (k_kmax(iblock-1,nl)+2)
        end do

        do iblock=2,mynblocks
          off_p3(iblock,nl) = off_p3(iblock-1,nl) + &
            (i_imax(iblock-1,nl)+3) * (j_jmax(iblock-1,nl)+3) * &
            (k_kmax(iblock-1,nl)+3)
        end do

      end if

      return
      end

!-----------------------------------------------------------------------
      subroutine zeitdata_in
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) key1,ierr,fid

      fid = 15
      open(fid,file='zeit.dat',form='formatted',status='old', &
           iostat=ierr,err=100)

      read(fid,*,iostat=ierr,err=100) key1

      if (unsteady.and.(.not.harbal)) then

        if (key1.eq.9990999) then
          write(*,*) 'No simtime is read from zeit.dat.'
          write(*,*) 'Only one state at UNspecified time will be read in &
     &.'
          found_simtime = .false.
        else if (key1.eq.9991999) then
          read(fid,*,iostat=ierr,err=100) simtime
          write(*,*) 'Found simtime in zeit.dat.'
          write(*,*) 'Only one state at specified time will be read in.'
          found_simtime = .true.
        else if ((key1.eq.9992999).and.(dualt)) then
          read(fid,*,iostat=ierr,err=100) simtime
          write(*,*) 'Found simtime in zeit.dat.'
          write(*,*) 'Two states at specified times will be read in.'
          found_simtime = .true.
          dual_prop_start = .true.
        else
          write(*,*) key1
          write(*,*) 'Illegal key1 in zeit.dat. Aborting'
          stop
        end if

      end if

      close(fid)
      return

 100  if (ierr.eq.29) then 
        write(*,*) 'Error detected in zeit.dat. Iostat =', ierr 
        write(*,*) 'File not found!' 
        stop
      else 
        write(*,*) 'Error detected in zeit.dat. Iostat =', ierr 
        stop
      end if

      end
