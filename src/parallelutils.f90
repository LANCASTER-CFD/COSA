!-----------------------------------------------------------------------
!     File to contain all the specific routines used for the OpenMP 
!     and MPI parallelisations of the COSA code.
!     
!     Developed at EPCC, The University of Edinburgh
!     Adrian Jackson, August 2010 ... June 2023
!-----------------------------------------------------------------------

      module parallelutils
  
      use cosa_precision

      implicit none
  
      private
        
      public :: initialise_cut_mpi_memory, allocate_cutman_memory, &
                allocate_pcutman_memory, allocate_check_memory
      public :: deallocate_cutman_memory, deallocate_pcutman_memory, &
                deallocate_check_memory
      public :: allocate_cutman_init_memory, allocate_check_init_memory, &
                allocate_pcutman_init_memory
      public :: mysendid, myrecvid, myrecvid_num, mysendid_num
      public :: mypsendid, myprecvid, myprecvid_num, mypsendid_num
      public :: check_mysendid, check_myrecvid, check_myrecvid_num, &
                check_mysendid_num
      public :: numsend, numrecv, pnumsend, pnumrecv
      public :: check_numsend, check_numrecv
      public :: maxproc_cutman, maxproc, mrequest_perproc
      public :: maxproc_pcutman, mrequest_pperproc
      public :: maxproc_check_cutman, mrequest_check_perproc
      public :: sendarray, receivearray
      public :: sendrequests, receiverequests, receiveindices
      public :: psendrequests, preceiverequests, preceiveindices
      public :: check_sendrequests, check_receiverequests, &
                check_receiveindices
      public :: surfblockstarts, surfblockends, nsurfblocks
      public :: surfblockindex, surfbcindex
        
      integer(kind=cosa_int), allocatable :: mysendid(:), myrecvid(:), myrecvid_num(:)
      integer(kind=cosa_int), allocatable :: mysendid_num(:)
      integer(kind=cosa_int), allocatable :: check_mysendid(:), check_myrecvid(:)
      integer(kind=cosa_int), allocatable :: check_myrecvid_num(:), &
                              check_mysendid_num(:)
      integer(kind=cosa_int), allocatable :: myprecvid(:), mypsendid(:)
      integer(kind=cosa_int), allocatable :: mypsendid_num(:), myprecvid_num(:)
      integer(kind=cosa_int), allocatable :: numsend(:), numrecv(:)
      integer(kind=cosa_int), allocatable :: pnumsend(:), pnumrecv(:)
      integer(kind=cosa_int), allocatable :: check_numsend(:), check_numrecv(:)
      integer(kind=cosa_int) :: maxproc_cutman, maxproc, mrequest_perproc
      integer(kind=cosa_int) :: maxproc_pcutman, mrequest_pperproc
      integer(kind=cosa_int) :: maxproc_check_cutman, mrequest_check_perproc
      
      integer(kind=cosa_int), allocatable :: sendrequests(:),receiverequests(:)
      integer(kind=cosa_int), allocatable :: receiveindices(:,:,:)
      
      integer(kind=cosa_int), allocatable :: check_sendrequests(:), &
                              check_receiverequests(:)
      integer(kind=cosa_int), allocatable :: check_receiveindices(:,:,:)

      real(kind=cosa_real), allocatable :: sendarray(:,:,:),receivearray(:,:,:)
      
      integer(kind=cosa_int), allocatable :: psendrequests(:),preceiverequests(:)
      integer(kind=cosa_int), allocatable :: preceiveindices(:,:,:)
      
      integer(kind=cosa_int), allocatable :: surfblockstarts(:,:)
      integer(kind=cosa_int), allocatable :: surfblockends(:,:)
      integer(kind=cosa_int), allocatable :: surfblockindex(:)
      integer(kind=cosa_int), allocatable :: surfbcindex(:)
      integer(kind=cosa_int) :: nsurfblocks

      logical :: initialised = .false.
      
      save
      
      contains
  
      subroutine initialise_cut_mpi_memory(lmaxproc)
      
      use cosa_precision

      implicit none
      
      integer(kind=cosa_int) :: lmaxproc
      
      if(.not. initialised) then
         maxproc = lmaxproc
         initialised = .true.
      end if
      
      end subroutine initialise_cut_mpi_memory
      
      subroutine allocate_cutman_init_memory()
      
      use cosa_precision

      implicit none

      integer(kind=cosa_int) :: data_size

      if(initialised) then
         
         if(maxproc-1 .eq. 0) then 
            data_size = 1
         else 
            data_size = maxproc - 1
         end if
          
         allocate(myrecvid_num(0:data_size))
         allocate(mysendid_num(0:data_size))
         allocate(mysendid(0:data_size))
         allocate(myrecvid(0:data_size))
         
      end if    
      
      end subroutine allocate_cutman_init_memory

      subroutine allocate_cutman_memory()
      
      use cosa_precision

      implicit none
      
      if(initialised) then
         
         allocate(receiveindices(mrequest_perproc,8,maxproc_cutman))
         allocate(numsend(maxproc_cutman))
         allocate(numrecv(maxproc_cutman))
         allocate(sendrequests(maxproc_cutman))
         allocate(receiverequests(maxproc_cutman))
         
      end if
      
      end subroutine allocate_cutman_memory
      
      subroutine allocate_check_init_memory()

      use cosa_precision

      implicit none

      integer(kind=cosa_int) :: data_size

      if(initialised) then

         if(maxproc-1 .eq. 0) then 
            data_size = 1
         else 
            data_size = maxproc - 1
         end if
         
         allocate(check_myrecvid_num(0:data_size))
         allocate(check_mysendid_num(0:data_size))
         allocate(check_mysendid(0:data_size))
         allocate(check_myrecvid(0:data_size))
         
      end if    
      
      end subroutine allocate_check_init_memory

      subroutine allocate_check_memory()
      
      use cosa_precision

      implicit none
      
      if(initialised) then
         
         allocate(check_receiveindices(mrequest_check_perproc,16, &
                                       maxproc_check_cutman))
         allocate(check_numsend(maxproc_check_cutman))
         allocate(check_numrecv(maxproc_check_cutman))
         allocate(check_sendrequests(maxproc_check_cutman))
         allocate(check_receiverequests(maxproc_check_cutman))
  
      end if
      
      end subroutine allocate_check_memory
      
      subroutine allocate_pcutman_init_memory()

      use cosa_precision

      implicit none

      integer(kind=cosa_int) :: data_size

      if(initialised) then

         if(maxproc-1 .eq. 0) then 
            data_size = 1
         else 
            data_size = maxproc - 1
         end if
         
         allocate(myprecvid_num(0:data_size))
         allocate(mypsendid_num(0:data_size))
         allocate(mypsendid(0:data_size))
         allocate(myprecvid(0:data_size))
         
      end if    
      
      end subroutine allocate_pcutman_init_memory
      
      subroutine allocate_pcutman_memory()
      
      use cosa_precision

      implicit none
      
      if(initialised) then
         
         allocate(pnumrecv(mrequest_pperproc))
         allocate(pnumsend(mrequest_pperproc))
         allocate(psendrequests(maxproc_pcutman))
         allocate(preceiverequests(maxproc_pcutman))
         allocate(preceiveindices(mrequest_pperproc,8,maxproc_pcutman))
         
      end if
      
      end subroutine allocate_pcutman_memory
      
      subroutine deallocate_cutman_memory()
      
      use cosa_precision

      implicit none
      
      if(allocated(mysendid)) then
         deallocate(mysendid, myrecvid, mysendid_num, myrecvid_num)
         deallocate(numsend, numrecv)
         deallocate(sendrequests, receiverequests, receiveindices)
      end if
      
      end subroutine deallocate_cutman_memory
      
      subroutine deallocate_check_memory()
      
      use cosa_precision

      implicit none
      
      if(allocated(check_mysendid)) then
         deallocate(check_mysendid, check_myrecvid, check_mysendid_num, &
                    check_myrecvid_num)
         deallocate(check_numsend, check_numrecv)
         deallocate(check_sendrequests, check_receiverequests, &
                    check_receiveindices)
      end if
      
      end subroutine deallocate_check_memory
      
      subroutine deallocate_pcutman_memory()
      
      use cosa_precision

      implicit none
      
      if(allocated(mypsendid)) then
         deallocate(mypsendid, myprecvid, mypsendid_num, myprecvid_num)
         deallocate(pnumsend, pnumrecv)
         deallocate(psendrequests, preceiverequests, preceiveindices)
      end if
      
      end subroutine deallocate_pcutman_memory
      
      end module parallelutils


!-----------------------------------------------------------------------
      subroutine usingparallel(parallel)
!     This routine provides switching/logical functionality for any 
!     part of the COSA code which requires to chose functionality 
!     depending upon whether the parallel (MPI) code is being used
!     or not.  It enables us to keep the main source code clean of 
!     pre-processing flags and other parallel functionality which 
!     is restricted to this file by design.
!-----------------------------------------------------------------------

      
      logical parallel

#ifdef MPI
      parallel = .TRUE.
#else
      parallel = .FALSE.
#endif  

      end


!-----------------------------------------------------------------------
      subroutine gettime(count,countrate)
!     This routine is designed to provide a usable timer for the OpenMP 
!     parallelisation as well as a general timer for the program.  It 
!     calls gettimemax which in turn wraps the system_clock timer, which
!     returns wall clock time (i.e. how long the program has been 
!     running for) rather than cpu or user time.  It 
!     returns to the user the count (since a certain point) in countrate 
!     per seconds.  system_clock is a wrapping timer (i.e. it goes back 
!     zero when it reaches a certain point) which needs to be considered 
!     for long program runs.
!-----------------------------------------------------------------------

      use cosa_precision

      implicit none

      integer(kind=cosa_int) :: count, countrate, countmax

      call gettimemax(count, countrate, countmax)

      return
      end

!-----------------------------------------------------------------------
      subroutine gettimemax(count,countrate,countmax)
!     This wraps the FORTRAN library routine system_clock (described 
!     in the comments for the gettime routine).  It is generally 
!     designed to be called by the user through the gettime routine
!     but if you require access to the countmax variable you can 
!     call this one directly.
!-----------------------------------------------------------------------

      use cosa_precision

      implicit none


      integer(kind=cosa_int) :: count, countrate, countmax

      call system_clock(count, countrate, countmax)

      return
      end

!-----------------------------------------------------------------------
      subroutine getmpitime(ttime)
!     getmpitime wraps the MPI_Wtime function.  It is designed only to 
!     be used with MPI versions of the program so has conditional 
!     compilation parts so that the MPI libray functionality is not 
!     included if the main program is not being compiled with MPI 
!     enabled.
!     MPI_Wtime returns a count since a certain point.
!-----------------------------------------------------------------------

      use cosa_precision

      implicit none

#ifdef MPI
      include 'mpif.h'
#endif

      double precision ttime

#ifdef MPI
      ttime = MPI_Wtime()
#else
      ttime = 0.d0
#endif

      return
      end

!-----------------------------------------------------------------------
      subroutine initialisempi()
!     Wraps MPI_INIT routine.  Does not take or return any parameters
!     but it does print out the number of parallel processors being 
!     used and the id of this particular process calling this routine.
!     If this is compiled without the parallel functionality enabled 
!     this routine does nothing.
!-----------------------------------------------------------------------

      use parallelutils
      use cosa_variables
      use cosa_precision

      implicit none

#ifdef MPI
      include 'mpif.h'

      integer(kind=cosa_int) :: ierror, parallelsize, parallelid

      call MPI_INIT(ierror)

      call getmpisize(parallelsize)
      call getmpiid(parallelid)

      call initialise_cut_mpi_memory(parallelsize)

!mpi  These variables are used in some of the output routines to flag up
!mpi  the type of data being written to file so the data can be processed 
!mpi  later
      typeint = 0
      typedouble = 1
      typechar = 2

      if(parallelid .eq. 0) then
      write(*,*) 'Using ',parallelsize,' MPI processes'
      end if

#else

      call initialise_cut_mpi_memory(1)

#endif


      return
      end

!-----------------------------------------------------------------------
      subroutine finalisempi()
!     Wraps MPI_FINALIZE routine.  Does not take or return any parameters
!     If this is compiled without the parallel functionality enabled 
!     this routine does nothing.
!-----------------------------------------------------------------------

      use cosa_precision

      implicit none

#ifdef MPI
      include 'mpif.h'

      integer(kind=cosa_int) :: ierror

      call MPI_FINALIZE(ierror)
#endif      

      return
      end

!-----------------------------------------------------------------------
      subroutine abortmpi()
!     Wraps MPI_ABORT routine.  Does not take or return any parameters
!     If this is compiled without the parallel functionality enabled 
!     this routine does nothing.
!-----------------------------------------------------------------------

      use cosa_precision

      implicit none

#ifdef MPI
      include 'mpif.h'

      integer(kind=cosa_int) :: ierror

      call MPI_ABORT(MPI_COMM_WORLD,ierror)

#else
      stop

#endif      



      return
      end

!-----------------------------------------------------------------------
      subroutine barrier()
!     Wraps MPI_BARRIER routine.  Does not take or return any parameters
!     If this is compiled without the parallel functionality enabled 
!     this routine does nothing.
!-----------------------------------------------------------------------

      use cosa_precision

      implicit none

#ifdef MPI
      include 'mpif.h'

      integer(kind=cosa_int) :: ierror

      call MPI_BARRIER(MPI_COMM_WORLD,ierror)
#endif      

      return
      end



!-----------------------------------------------------------------------
      subroutine getmpisize(msize)
!     Wraps MPI_COMM_SIZE routine.  Returns an integer which is set to 
!     the number of processors/cores used if this is an MPI program
!     and to 1 if this isn't an MPI program.
!-----------------------------------------------------------------------

      use cosa_precision

      implicit none

#ifdef MPI
      include 'mpif.h'
#endif

      integer(kind=cosa_int) :: msize, ierror

#ifdef MPI
      call MPI_COMM_SIZE(MPI_COMM_WORLD,msize,ierror)
#else
      msize = 1
#endif
      

      return
      end

!-----------------------------------------------------------------------
      subroutine getmpiid(rank)
!     Wraps MPI_COMM_RANK routine.  Returns an integer which is set to 
!     the rank(id) of the call process if this is an MPI program
!     and to 0 if this isn't an MPI program.
!-----------------------------------------------------------------------

      use cosa_precision

      implicit none

#ifdef MPI
      include 'mpif.h'
#endif

      integer(kind=cosa_int) :: rank, ierror

#ifdef MPI
      call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)
#else
      rank = 0
#endif
      

      return
      end


!-----------------------------------------------------------------------
      subroutine amcontroller(trueorfalse)
!     A utility subroutine which returns true if the call process is 
!     of mpi rank 0 (or this is a serial/openmp program and false 
!     otherwise.  
!-----------------------------------------------------------------------

      use cosa_precision

      implicit none

      integer(kind=cosa_int) :: rank
      logical trueorfalse

      call getmpiid(rank)

      if(rank .eq. 0) then
         trueorfalse = .TRUE.
      else
         trueorfalse = .FALSE.
      end if

      return 
      end


!-----------------------------------------------------------------------
      subroutine paralleldecomposition()
!     The subroutine which calculates and constructs the block 
!     decomposition for the parallel program.  It can handle only
!     one process and also more processes than blocks.
!-----------------------------------------------------------------------

      use cosa_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) :: parallelsize, parallelid, i, blockoffset
      integer(kind=cosa_int) :: wholeblocks, remainderblocks, remblock
      integer(kind=cosa_int) :: iblk1, iblk2, icut

      call getmpisize(parallelsize)
      call getmpiid(parallelid)

      if(parallelsize .lt. 1) then
         write(*,*) 'error with paralleldecomposition parallel size ', &
              'less than 1'
         stop
      else if(parallelsize .eq. 1) then
         lowernblock = 1
         uppernblock = nblocks
         mynblocks   = nblocks
         myncuts     = ncuts
         mynpcuts    = npcuts
         maxblocks   = nblocks
         do i=1,nblocks
            blockassignment(i) = 0
         end do
         do i=1,ncuts
            mycuts(i) = i
         end do
         do i=1,npcuts
            mypcuts(i) = i
         end do
         return
      end if

      wholeblocks = nblocks/parallelsize

      if(wholeblocks .eq. 0) then
!     mpi     This is the case where there are less blocks than 
!     mpi     processors used so we assign one block to each 
!     mpi     processor until they are all assigned
         if((parallelid + 1) .le. nblocks) then
            lowernblock = parallelid + 1
            uppernblock = lowernblock
            mynblocks = 1
            maxblocks = 1
         else
            mynblocks = 0
            lowernblock = 0
            uppernblock = -1
            maxblocks = 1
         end if
         
!     mpi     Build a parallel id map for blocks to process ranks
         do i=1,nblocks
!     mpi     Ranks in MPI start at 0 so we must offset from the 
!     mpi     nblock number
            blockassignment(i) = i-1
         end do

      else
!     mpi     This is the case where there are more blocks than 
!     mpi     processors (which is the ideal case).  
!     mpi     By default we assume that the number of blocks 
!     mpi     exactly divides by the number of processors
         mynblocks = wholeblocks
!     mpi     Now we check whether the above assumption (i.e. the
!     mpi     number of blocks exactly divides by the number of 
!     mpi     processors) is correct.
         remainderblocks = mod(nblocks, parallelsize)
         if(remainderblocks .ne. 0) then
            maxblocks = mynblocks + 1
         else
            maxblocks = mynblocks
         end if
!     mpi     Check what blocks are remaining in the division 
!     mpi     and if there are blocks remaining then assign a 
!     mpi     block to each processor until they are all used up
         if((parallelid + 1) .le. remainderblocks) then
            mynblocks = mynblocks + 1
         end if

!     mpi     Build a parallel id map for blocks to process ranks
!     mpi     remblocks calculates point in the number of blocks when the 
!     mpi     over assignment of blocks stop.  We over assign blocks to 
!     mpi     fill up the first processes with more blocks than the rest 
!     mpi     to ensure all the blocks are assigned to someone.  This is 
!     mpi     because in this situation we have more blocks than processes 
!     mpi     and the number of blocks do not evenly divide into processes 
!     mpi     so we need to assign spare blocks (once we've done as even 
!     mpi     a division as possible) to processes to ensure all blocks are 
!     mpi     assigned.  The strategy we have taken is to increase the number 
!     mpi     of blocks the first processes get by one.  The remblock variable 
!     mpi     holds the point in the block list where we have assigned all the 
!     mpi     spare blocks and can go back to assigning the even number of blocks
!     mpi     (wholeblocks) to processes.
         remblock = (wholeblocks + 1) * remainderblocks
         do i=1,nblocks
            if(remainderblocks .ne. 0) then
               if(i .gt. remblock) then
                  blockassignment(i) = ((remblock/(wholeblocks+1)) + &
                       ((i-remblock-1)/wholeblocks))
               else
                  blockassignment(i) = ((i-1)/(wholeblocks+1))
               end if
            else
               blockassignment(i) = ((i-1)/wholeblocks)
            end if
         end do



!     mpi     Work out the blockoffset.  This is the amount of blocks 
!     mpi     up to the upper limit of the blocks this process owns.
!     mpi     This is slightly complicated because it needs to be able 
!     mpi     to deal with the case where processors have different 
!     mpi     numbers of blocks (i.e. some processors have one more 
!     mpi     block that others).

!     mpi     This is the situation where there are an uneven number 
!     mpi     of blocks but this processor has the lower number of blocks 
!     mpi     so we need to calculate the offset by working out the total
!     mpi     blocks used by processors with the higher block count then 
!     mpi     adding the number of blocks used by this processor and the 
!     mpi     preceding processors who have the lower number of blocks
         if(remainderblocks .ne. 0 .AND. &
              (parallelid + 1) .gt. remainderblocks) then
            blockoffset = ((mynblocks + 1) * remainderblocks) + &
                 ((parallelid + 1 - remainderblocks) * mynblocks)
!     mpi     This code covers the situtations of calculating the offset 
!     mpi     where either there are no remainder blocks (i.e. all 
!     mpi     processors have the same number of blocks) or there are 
!     mpi     remainder blocks but this processor has the higher number 
!     mpi     of blocks.
         else
            blockoffset = (parallelid + 1) * mynblocks
         end if
!     mpi     Calculate the number of the highest and lowest block this 
!     mpi     processors owns (i.e. the block range).  We are assuming 
!     mpi     that the blocks for each processors are contiguous (i.e. 
!     mpi     all the blocks are together).
         uppernblock = blockoffset
         lowernblock = uppernblock - mynblocks + 1
      end if
      
      if(ncuts .gt. mcut) then
         write(*,*) 'Error in paralleldecompositon, ncuts (',ncuts, &
              ') greater than ',mcut,' which is our current maximum ', &
              'size for the mycuts array, this needs increased.'
         stop
      end if
      
      if(npcuts .gt. mcut) then
         write(*,*) 'Error in paralleldecompositon, npcuts (',npcuts, &
              ') greater than ',mcut,' which is our current maximum ', &
              'size for the mypcuts array, this needs increased.'
         stop
      end if

      myncuts = 0

      do icut=1,ncuts
         iblk1 = cutdata(1,icut)
         iblk2 = cutdata(10,icut)
         if(((iblk1 .ge. lowernblock) .and. &
              (iblk1 .le. uppernblock)) &
              .or. ((iblk2 .ge. lowernblock) .and. &
              (iblk2 .le. uppernblock))) then
            myncuts = myncuts + 1
            mycuts(myncuts) = icut
         end if
      end do

      mynpcuts = 0

      do icut=1,npcuts
         iblk1 = pcutdata(1,icut)
         iblk2 = pcutdata(10,icut)
         if(((iblk1 .ge. lowernblock) .and. &
              (iblk1 .le. uppernblock)) &
              .or. ((iblk2 .ge. lowernblock) .and. &
              (iblk2 .le. uppernblock))) then
            mynpcuts = mynpcuts + 1
            mypcuts(mynpcuts) = icut
         end if
      end do

      return 
      end

!-----------------------------------------------------------------------
      subroutine movedatafordecomposition()
!     This routine is designed to move from a serial layout (as 
!     assumed in input.dat data) where a process has blocks 1 to nblocks
!     and cuts 1 to ncuts to a situation where each process only has a 
!     certain subset of cuts and blocks (as is required by the MPI 
!     parallel program).  As all the serial data is read in by each 
!     process in the parallel program what this routine does is move 
!     the blocks that this process has been assigned to the 
!     beginning of the arrays so they can be used in 
!     do iblock=1,mynblocks and icut=1,myncuts loops.
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) :: myblock, parallelsize, nl, ibc, iblock
      integer(kind=cosa_int) :: icut, currentcut, rank
      call getmpisize(parallelsize)
      call getmpiid(rank)

!     mpi     Store the whole i,j,k, and ijk max values for use later
!     mpi     (in routines like readgrid, or the surf tec file production)
      do iblock=1,nblocks
         do nl = 1,nlevel
            g_i_imax(iblock,nl) = i_imax(iblock,nl)
            g_j_jmax(iblock,nl) = j_jmax(iblock,nl)
            g_k_kmax(iblock,nl) = k_kmax(iblock,nl)
            g_ijk_ijkmax(iblock,nl) = ijk_ijkmax(iblock,nl)
            g_n_wall(iblock,nl) = n_wall(iblock,nl)
         end do
         g_nbcs(iblock) = nbcs(iblock)
      end do

!     mpi  Store the full boundary condition for later use in functionality such
!     mpi  as constructing surf tec files
      do iblock=1,nblocks
         do ibc = 1,nbcs(iblock)
            g_bcdata(1,ibc,iblock) = bcdata(1,ibc,iblock)
            g_bcdata(2,ibc,iblock) = bcdata(2,ibc,iblock)
            g_bcdata(3,ibc,iblock) = bcdata(3,ibc,iblock)
            g_bcdata(4,ibc,iblock) = bcdata(4,ibc,iblock)
            g_bcdata(5,ibc,iblock) = bcdata(5,ibc,iblock)
            g_bcdata(6,ibc,iblock) = bcdata(6,ibc,iblock)
            g_bcdata(7,ibc,iblock) = bcdata(7,ibc,iblock)
            g_bcdata(8,ibc,iblock) = bcdata(8,ibc,iblock)
            g_bcdata(9,ibc,iblock) = bcdata(9,ibc,iblock)      
         end do
      end do

      if(parallelsize .gt. 1) then

         myblock = 1

         do iblock = lowernblock,uppernblock
            nbcs(myblock) = nbcs(iblock)
            
            do nl = 1,nlevel
               i_imax(myblock,nl) = i_imax(iblock,nl) 
               j_jmax(myblock,nl) = j_jmax(iblock,nl)
               k_kmax(myblock,nl) = k_kmax(iblock,nl)
               ijk_ijkmax(myblock,nl) = ijk_ijkmax(iblock,nl)
               n_wall(myblock,nl) = n_wall(iblock,nl)
            end do
            
            do ibc = 1,nbcs(myblock)
               bcdata(1,ibc,myblock) = bcdata(1,ibc,iblock)
               bcdata(2,ibc,myblock) = bcdata(2,ibc,iblock)
               bcdata(3,ibc,myblock) = bcdata(3,ibc,iblock)
               bcdata(4,ibc,myblock) = bcdata(4,ibc,iblock)
               bcdata(5,ibc,myblock) = bcdata(5,ibc,iblock)
               bcdata(6,ibc,myblock) = bcdata(6,ibc,iblock)
               bcdata(7,ibc,myblock) = bcdata(7,ibc,iblock)
               bcdata(8,ibc,myblock) = bcdata(8,ibc,iblock)
               bcdata(9,ibc,myblock) = bcdata(9,ibc,iblock)
            end do
            
            myblock = myblock + 1

         end do
         
         do icut=1,myncuts

            currentcut = mycuts(icut)

            cutdata( 1,icut) = cutdata( 1,currentcut)
            cutdata( 2,icut) = cutdata( 2,currentcut)
            cutdata(10,icut) = cutdata(10,currentcut)
            cutdata(11,icut) = cutdata(11,currentcut)
            cutdata(12,icut) = cutdata(12,currentcut)
            cutdata(19,icut) = cutdata(19,currentcut)
            cutdata(20,icut) = cutdata(20,currentcut)
            cutdata(21,icut) = cutdata(21,currentcut)
            cutdata( 3,icut) = cutdata( 3,currentcut)
            cutdata( 4,icut) = cutdata( 4,currentcut)
            cutdata( 5,icut) = cutdata( 5,currentcut)
            cutdata( 6,icut) = cutdata( 6,currentcut)
            cutdata( 7,icut) = cutdata( 7,currentcut)
            cutdata( 8,icut) = cutdata( 8,currentcut)
            cutdata( 9,icut) = cutdata( 9,currentcut)
            cutdata(12,icut) = cutdata(12,currentcut)
            cutdata(13,icut) = cutdata(13,currentcut)
            cutdata(14,icut) = cutdata(14,currentcut)
            cutdata(15,icut) = cutdata(15,currentcut)
            cutdata(16,icut) = cutdata(16,currentcut)
            cutdata(17,icut) = cutdata(17,currentcut)
            cutdata(18,icut) = cutdata(18,currentcut)

         end do

         do icut=1,mynpcuts

            currentcut = mypcuts(icut)

            pcutdata( 1,icut) = pcutdata( 1,currentcut)
            pcutdata( 2,icut) = pcutdata( 2,currentcut)
            pcutdata(10,icut) = pcutdata(10,currentcut)
            pcutdata(11,icut) = pcutdata(11,currentcut)
            pcutdata(12,icut) = pcutdata(12,currentcut)
            pcutdata(19,icut) = pcutdata(19,currentcut)
            pcutdata(20,icut) = pcutdata(20,currentcut)
            pcutdata(21,icut) = pcutdata(21,currentcut)
            pcutdata(22,icut) = pcutdata(22,currentcut)
            pcutdata( 3,icut) = pcutdata( 3,currentcut)
            pcutdata( 4,icut) = pcutdata( 4,currentcut)
            pcutdata( 5,icut) = pcutdata( 5,currentcut)
            pcutdata( 6,icut) = pcutdata( 6,currentcut)
            pcutdata( 7,icut) = pcutdata( 7,currentcut)
            pcutdata( 8,icut) = pcutdata( 8,currentcut)
            pcutdata( 9,icut) = pcutdata( 9,currentcut)
            pcutdata(12,icut) = pcutdata(12,currentcut)
            pcutdata(13,icut) = pcutdata(13,currentcut)
            pcutdata(14,icut) = pcutdata(14,currentcut)
            pcutdata(15,icut) = pcutdata(15,currentcut)
            pcutdata(16,icut) = pcutdata(16,currentcut)
            pcutdata(17,icut) = pcutdata(17,currentcut)
            pcutdata(18,icut) = pcutdata(18,currentcut)

         end do

      end if

      return
      end

!-----------------------------------------------------------------------
      subroutine myblock(blockid, ismyblock, updateid)
!     myblock calculates whether the a given blockid is owned 
!     by the calling process.  If it is then the logical ismyblock 
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) :: blockid
      logical ismyblock
      logical updateid

      ismyblock = .false.

      if((blockid .ge. lowernblock) .and. (blockid .le. uppernblock)) &
           then
         ismyblock = .true.
         if(updateid) then
            blockid = blockid - (lowernblock - 1)
         end if
      end if

      return 
      end

!-----------------------------------------------------------------------
      subroutine myblockglobalid(blockid, updatedid)
!     myblockglobalid calculates the original global block id for 
!     the provided blockid. 
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) :: blockid, updatedid

      updatedid = blockid + (lowernblock - 1)

      return 
      end

!-----------------------------------------------------------------------
      subroutine gridbinaryoffset(offset, nl, integersize, doublesize)
!     gridoffset computes where the data for this process is in the 
!     grid data file (mesh.dat).  It is done this way as each block
!     read in separately from the grid data file and the blocks can 
!     vary in size so we have to go through each preceding block 
!     that this process does not own an work out how much of the grid 
!     data file they occupy to work out where this process data resides
!     in the file.
!-----------------------------------------------------------------------
      
      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none
      
#ifdef MPI
      include 'mpif.h'
      integer(kind=cosa_int) :: offset_kind
      parameter(offset_kind=MPI_OFFSET_KIND)
#else
      integer(kind=cosa_int) :: offset_kind
      parameter(offset_kind=8)
#endif

      integer(kind=offset_kind) :: offset
      integer(kind=cosa_int) :: nl
      integer(kind=cosa_int) :: i, j, k, n, imax, jmax, kmax
      integer(kind=cosa_int) :: integersize, doublesize, ierr
#ifdef MPI
      call mpi_type_size(MPI_DOUBLE_PRECISION, doublesize, ierr)
      call mpi_type_size(MPI_INTEGER, integersize, ierr)

      offset = 0

      do n = 1,(lowernblock-1)
!     This accounts for the integer before and after the data in the 
!     fortran binary format that records how many bytes are in the record.
         offset = offset + 2*integersize
         imax = g_i_imax(n,nl)
         jmax = g_j_jmax(n,nl)
         kmax = g_k_kmax(n,nl)
         do i = 1,imax
            do j = 1,jmax
               do k = 1,kmax
                  offset = offset + 3*doublesize
               end do
            end do
         end do
         if (moving.and.pitching) then
            offset = offset + 2*integersize
            do i = 1,imax
               do j = 1,jmax
                  do k = 1,kmax
                     offset = offset + 3*doublesize
                  end do
               end do
            end do
         end if
      end do
#endif
      return
      end

!-----------------------------------------------------------------------
      subroutine readrestartinitialoffset(totalblockamount, nl, &
                 integersize, doublesize)
!     readrestartinitialoffset computes where the data for this process is 
!     in the restart data file (restart).  It is done this way as each 
!     block read and write separately from the restart file and the 
!     blocks can vary in size so we have to go through each preceding 
!     block that this process does not own an work out how much of the 
!     restart data file they occupy to work out where this process data 
!     resides in the file.
!     This is further complicated as the data written to the restart 
!     file varies depending upon the specifics of the simulation being 
!     run.  
!-----------------------------------------------------------------------
      
      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none
      
#ifdef MPI
      include 'mpif.h'
#endif

      integer(kind=cosa_int) :: nl, k, imax, jmax, kmax
      integer*8 totalblockamount, blockamount
      integer(kind=cosa_int) :: doublesize, integersize
      integer(kind=cosa_int) :: ierr

      totalblockamount = 0
#ifdef MPI

!     mpi  Calculate the offset position in the file that this block 
!     mpi  will be written to.  It is based on the blocks in the file 
!     mpi  that will be written by processes that have lower blocks 
!     mpi  plus the blocks this process has already written.
      call mpi_type_size(MPI_DOUBLE_PRECISION, doublesize, ierr)
      call mpi_type_size(MPI_INTEGER, integersize, ierr)
     
      do k = 1,(lowernblock-1)
         imax = g_i_imax(k,nl)
         jmax = g_j_jmax(k,nl)
         kmax = g_k_kmax(k,nl)

         call restartblockoffset(blockamount, imax, jmax, kmax, &
              integersize, doublesize)
         
         totalblockamount = totalblockamount + blockamount

         call readrestartextraoffset(blockamount, imax, jmax, kmax, &
              integersize, doublesize)

         totalblockamount = totalblockamount + blockamount

      end do
#endif

      return
      end

!-----------------------------------------------------------------------
      subroutine writerestartinitialoffset(totalblockamount, nl, &
                 integersize, doublesize)
!     writerestartinitialoffset computes where the data for this process is 
!     in the restart data file (restart).  It is done this way as each 
!     block read and write separately from the restart file and the 
!     blocks can vary in size so we have to go through each preceding 
!     block that this process does not own an work out how much of the 
!     restart data file they occupy to work out where this process data 
!     resides in the file.
!     This is further complicated as the data written to the restart 
!     file varies depending upon the specifics of the simulation being 
!     run.  
!-----------------------------------------------------------------------
      
      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none
      
#ifdef MPI
      include 'mpif.h'
#endif

      integer(kind=cosa_int) :: nl, k, imax, jmax, kmax
      integer*8 totalblockamount, blockamount
      integer(kind=cosa_int) :: doublesize, integersize
      integer(kind=cosa_int) :: ierr

      totalblockamount = 0
#ifdef MPI

!     mpi  Calculate the offset position in the file that this block 
!     mpi  will be written to.  It is based on the blocks in the file 
!     mpi  that will be written by processes that have lower blocks 
!     mpi  plus the blocks this process has already written.
      call mpi_type_size(MPI_DOUBLE_PRECISION, doublesize, ierr)
      call mpi_type_size(MPI_INTEGER, integersize, ierr)
     
      do k = 1,(lowernblock-1)
         imax = g_i_imax(k,nl)
         jmax = g_j_jmax(k,nl)
         kmax = g_k_kmax(k,nl)

         call restartblockoffset(blockamount, imax, jmax, kmax, &
              integersize, doublesize)
         
         totalblockamount = totalblockamount + blockamount

         call writerestartextraoffset(blockamount, imax, jmax, kmax, &
              integersize, doublesize)

         totalblockamount = totalblockamount + blockamount

      end do
#endif

      return
      end

!-----------------------------------------------------------------------
      subroutine restartblockoffset(blockamount, imax, jmax, kmax, &
                 integersize, doublesize)
!     restartblockoffset calculates the amount of space in a file 
!     the basic block data occupies.
!-----------------------------------------------------------------------
      
      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none
      
      integer(kind=cosa_int) :: imax, jmax, kmax
      integer*8 blockamount 
      integer(kind=cosa_int) :: doublesize, integersize
      integer(kind=cosa_int) :: expr4, numberofints


      expr4 = npde
      numberofints = 2
      blockamount = ((2*nharms)+1) * &
                    (doublesize*expr4*(imax+3)*(jmax+3)*(kmax+3)) + &
                    numberofints*integersize


      return
      end

!-----------------------------------------------------------------------
      subroutine readrestartextraoffset(blockamount, imax, jmax, kmax, &
                 integersize, doublesize)
!     readrestartextraoffset calculates the amount of space in the restart 
!     file extra data requires for this simulation.  Some simulations will
!     not have any extra data, but some do so this routine calculates
!     what is required.
!-----------------------------------------------------------------------
      
      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none
      
      integer(kind=cosa_int) :: imax, jmax, kmax
      integer*8 blockamount
      integer(kind=cosa_int) :: expr4, doublesize, integersize, numberofints
      
      blockamount = 0

      if(dual_prop_start) then

         expr4 = npde
         numberofints = 2

         blockamount = doublesize*(imax+3)*(jmax+3)*(kmax+3)*expr4 + &
                       numberofints*integersize

      end if

      if (kom.or.kom_bsl.or.kom_sst) then

         expr4 = ((2*nharms)+1)
         numberofints = 2

         blockamount = blockamount + &
                       (doublesize*(imax+3)*(jmax+3)*(kmax+3)*expr4 + &
                       numberofints*integersize)
          
      end if

      return
      end


!-----------------------------------------------------------------------
      subroutine writerestartextraoffset(blockamount, imax, jmax, kmax, &
                 integersize, doublesize)
!     writerestartextraoffset calculates the amount of space in the restart 
!     file extra data requires for this simulation.  Some simulations will
!     not have any extra data, but some do so this routine calculates
!     what is required.
!-----------------------------------------------------------------------
      
      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none
      
      integer(kind=cosa_int) :: imax, jmax, kmax
      integer*8 blockamount
      integer(kind=cosa_int) :: expr4, doublesize, integersize, numberofints
      
      blockamount = 0

      if(dualt) then

         expr4 = npde
         numberofints = 2

         blockamount = doublesize*(imax+3)*(jmax+3)*(kmax+3)*expr4 + &
                       numberofints*integersize

      end if

      if (kom.or.kom_bsl.or.kom_sst) then

         expr4 = ((2*nharms)+1)
         numberofints = 2

         blockamount = blockamount + &
                       (doublesize*(imax+3)*(jmax+3)*(kmax+3)*expr4 + &
                       numberofints*integersize)

      end if

      return
      end 

!-----------------------------------------------------------------------
      subroutine newflowtecinitialoffset(totalblockamount, nl, iblock, &
                 integersize, doublesize, charactersize)
!-----------------------------------------------------------------------
!     newflowtecinitialoffset computes where the data for this process is 
!     in the flowtec data file.  It is done this way as each 
!     block is written separately to the flowtec files and the 
!     blocks can vary in size so we have to go through each preceding 
!     block that this process does not own and work out how much of the 
!     flowtec data files they occupy to work out where this process's data 
!     resides in the file.
!     This is further complicated as the data written to the flowtec 
!     files varies depending upon the specifics of the simulation being 
!     run.  
!-----------------------------------------------------------------------
      
      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none
      
#ifdef MPI
      include 'mpif.h'
#endif

      integer(kind=cosa_int) :: nl, k, imax, jmax, kmax, iblock, upperbound
      integer*8 totalblockamount
      integer(kind=cosa_int) :: blockamount, doublesize, integersize
      integer(kind=cosa_int) :: charactersize, ierr

      totalblockamount = 0
#ifdef MPI

      call mpi_type_size(MPI_DOUBLE_PRECISION, doublesize, ierr)
      call mpi_type_size(MPI_INTEGER, integersize, ierr)
      call mpi_type_size(MPI_CHARACTER, charactersize, ierr)
     
      upperbound = (lowernblock-1) + (iblock-1)

!     The 7 below is due to 7 integer control variables written by
!     routine writeparallelflowtecheader. This variables are read in
!     by the utilities that convert the MPI binary flowtec files to
!     other formats.
      totalblockamount = 7*integersize

!     If unsteady we have written the simtime into the file
      if(unsteady) then
         totalblockamount = totalblockamount + doublesize
      end if

      totalblockamount = totalblockamount + (integersize * 3 * nblocks)

      do k = 1,upperbound
         imax = g_i_imax(k,nl)
         jmax = g_j_jmax(k,nl)
         kmax = g_k_kmax(k,nl)

         call newflowtecblockoffset(blockamount, imax, jmax, kmax, &
              integersize, doublesize, charactersize)
         
         totalblockamount = totalblockamount + blockamount

      end do

#endif

      return
      end

!-----------------------------------------------------------------------
      subroutine newflowtecblockoffset(blockamount, imax, jmax, kmax, &
                 integersize, doublesize, charactersize)
!-----------------------------------------------------------------------
!     newflowtecblockoffset calculates the amount of space in a file 
!     a specific block occupies.
!-----------------------------------------------------------------------
      
      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none
      
      integer(kind=cosa_int) :: imax, jmax, kmax
      integer(kind=cosa_int) :: blockamount, doublesize, integersize, charactersize

      blockamount = 0
      
!---- msc/jd, 27 July 2015.
!     At present there are 14/12 variables in turbulent 
!     compressible/incompressible TECPLOT file and
!     10/8 variables in Euler/laminar compressible/incompressible 
!     TECPLOT file. If these numbers 
!     will change, the npde*2 factor in 2 if branches below will have to
!     be changed.
      if(kom.or.kom_bsl.or.kom_sst) then
         blockamount = blockamount + &
                       (doublesize*(imax+1)*(jmax+1)*(kmax+1)*npde*2)
      else
         blockamount = blockamount + &
                       (doublesize*(imax+1)*(jmax+1)*(kmax+1)*npde*2)
      end if

      return
      end

!-----------------------------------------------------------------------
      subroutine surftecinitialoffset(totalblockamount, nl, iblock, &
           nsurfblocks, surfblockstarts, surfblockends, surfblockindex, &
           surfbcindex, surfacenumber, integersize, doublesize, &
           charactersize)
!-----------------------------------------------------------------------
!     surftecinitialoffset computes where the data for this process is 
!     in the surftec data file.  It is done this way as each 
!     block is written separately to the surftec files and the 
!     blocks can vary in size so we have to go through each preceding 
!     block that this process does not own and work out how much of the 
!     surftec data files they occupy to work out where this process's data 
!     resides in the file.
!     This is further complicated as the data written to the surftec 
!     files varies depending upon the specifics of the simulation being 
!     run.  
!-----------------------------------------------------------------------
      
      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none
      
#ifdef MPI
      include 'mpif.h'
#endif

      integer(kind=cosa_int) :: nl, k, l, imax, jmax, kmax, iblock, upperbound
      integer*8 totalblockamount
      integer(kind=cosa_int) :: nsurfblocks, surfacenumber
      integer(kind=cosa_int) :: surfblockstarts(nsurfblocks,3)
      integer(kind=cosa_int) :: surfblockends(nsurfblocks,3)
      integer(kind=cosa_int) :: surfblockindex(nsurfblocks)
      integer(kind=cosa_int) :: surfbcindex(nsurfblocks)
      integer(kind=cosa_int) :: blockamount, doublesize, integersize
      integer(kind=cosa_int) :: charactersize, ierr

      totalblockamount = 0
#ifdef MPI

      call mpi_type_size(MPI_DOUBLE_PRECISION, doublesize, ierr)
      call mpi_type_size(MPI_INTEGER, integersize, ierr)
      call mpi_type_size(MPI_CHARACTER, charactersize, ierr)     

!     The 6 below is due to 6 integer control variables written by
!     routine writeparallelsurftecheader. This variables are read in
!     by the utilities that convert the MPI binary surftec files to
!     other formats.
      totalblockamount = 6*integersize

!     These are the surface block start end end index amounts
      totalblockamount = totalblockamount + (integersize * 6 * nsurfblocks)

      surfacenumber = 0

      upperbound = (lowernblock-1) + (iblock-1)

      do k = 1,upperbound
         do l=1,g_nbcs(k)

            call surftecblockoffset(k, l, blockamount, nsurfblocks, &
                 surfblockstarts, surfblockends, surfblockindex, surfbcindex, &
                 surfacenumber, integersize, doublesize, charactersize)
            
            totalblockamount = totalblockamount + blockamount
         end do         
      end do

!     AJ This is required because we are calculating the offset here, i.e. the 
!     AJ next block to be taken from the surfblockends etc... arrays. Up until
!     AJ this point we've calculated the last block to be seen.
      surfacenumber = surfacenumber + 1
      
#endif

      return
      end

!-----------------------------------------------------------------------
      subroutine surftecblockoffset(iblock, ibc, blockamount, &
                 nsurfblocks, surfblockstarts, surfblockends, &
                 surfblockindex, surfbcindex, surfacenumber, &
                 integersize, doublesize, charactersize)
!-----------------------------------------------------------------------
!     surftecblockoffset calculates the amount of space in a file 
!     a specific block occupies.
!-----------------------------------------------------------------------
      
      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none
      
      integer(kind=cosa_int) :: blockamount, doublesize, integersize, charactersize, ibc
      integer(kind=cosa_int) :: iblock, nsurfblocks, surfacenumber, isurface
      integer(kind=cosa_int) :: surfblockstarts(3, nsurfblocks)
      integer(kind=cosa_int) :: surfblockends(3, nsurfblocks)
      integer(kind=cosa_int) :: surfblockindex(nsurfblocks)
      integer(kind=cosa_int) :: surfbcindex(nsurfblocks)

      blockamount = 0
      do isurface = 1,nsurfblocks
         if(surfblockindex(isurface) .eq. iblock) then  
            if(surfbcindex(isurface) .eq. ibc) then
               blockamount = 6 * doublesize
               blockamount = &
                    (((surfblockends(1,isurface)-surfblockstarts(1,isurface))+1) * &
                    ((surfblockends(2,isurface)-surfblockstarts(2,isurface))+1) * &
                    ((surfblockends(3,isurface)-surfblockstarts(3,isurface))+1)) &
                    * blockamount
               surfacenumber = isurface
               exit
            end if
         end if
      end do
            
      return
      end

!-----------------------------------------------------------------------
      subroutine setuprestartiovalues1(imax, jmax, kmax, npde, nharms, &
           amount1)
!     This subroutine is used to calculate the amount of data written 
!     to the restart file in parts of the I/O code.  It is necessary 
!     to avoid the problems that occur when one component of the data 
!     to be written out is of zero size.  Without this code that would
!     cause no data to be written at all.  This code fixes that to 
!     ensure that having part of the data being zero extend does not 
!     restrict the other data being written out.
!     This particular routine is setup to calculate the amounts of data 
!     for the n = 0,2*nharms loop in the writerestartblock subroutine.
!-----------------------------------------------------------------------
      
      use common_variables
      use cosa_precision

      implicit none
      
      integer(kind=cosa_int) :: imax, jmax, kmax, npde, nharms
      integer(kind=cosa_int) :: amount1, expr1, expr2, expr3, expr4, expr5

      expr5 = (2*nharms) + 1

!     mpi  Workout the npde value separately as it is used in all these 
!     mpi  amounts in the same way.
      if(npde .eq. 0) then
         expr4 = 1
      else
         expr4 = npde
      end if
      
!     mpi  Workout the amount for the first write line i.e. 
!          (imax+2*jmax-1*npde)      
!     mpi  If all components are zero then it is correct to make the amount
!     mpi  written zero.
      if(((imax+3) .eq. 3) .AND. ((jmax+3) .eq. 3) .AND. &
         ((kmax+3) .eq. 3) .AND. (npde .eq. 0)) then
         amount1 = 0
      else
!     mpi  All components are not zero so go through them individually and 
!     mpi  fix them if individual ones are zero (i.e. set them to 1 so they 
!     mpi  don't affect the final multiplication)
         if((imax+3) .eq. 3) then
            expr1 = 1
         else
            expr1 = imax+3
         end if
         
         if((jmax+3) .eq. 3) then
            expr2 = 1
         else
            expr2 = jmax+3
         end if

         if((kmax+3) .eq. 3) then
            expr3 = 1
         else
            expr3 = kmax+3
         end if
         
!     mpi     Work out the total amount of data to be written
         amount1 = expr1*expr2*expr3*expr4*expr5

      end if

      return
      end 

!-----------------------------------------------------------------------
      subroutine setuprestartiovalues2(imax, jmax, kmax, npde, amount1)
!     This subroutine is used to calculate the amount of data written 
!     to the restart file in parts of the I/O code.  It is necessary 
!     to avoid the problems that occur when one component of the data 
!     to be written out is of zero size.  Without this code that would
!     cause no data to be written at all.  This code fixes that to 
!     ensure that having part of the data being zero extend does not 
!     restrict the other data being written out.
!     This particular routine is setup to calculate the amount of data 
!     in the if(dualt) part of the writerestartblock routine.
!-----------------------------------------------------------------------
      
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) :: imax, jmax, kmax, npde
      integer(kind=cosa_int) :: amount1, expr1, expr2, expr3, expr4

      if(((imax+3) .eq. 3) .AND. ((jmax+3) .eq. 3) .AND. &
         ((kmax+3) .eq. 3) .AND. (npde .eq. 0)) then
         amount1 = 0
      else

         if((imax+3) .eq. 3) then
            expr1 = 1
         else
            expr1 = imax+3
         end if
         
         if((jmax+3) .eq. 3) then
            expr2 = 1
         else
            expr2 = jmax+3
         end if
         
         if((kmax+3) .eq. 3) then
            expr3 = 1
         else
            expr3 = kmax+3
         end if

         if(npde .eq. 0) then
            expr4 = 1
         else
            expr4 = npde
         end if
         
         amount1 = expr1*expr2*expr3*expr4

      end if

      return
      end 

!-----------------------------------------------------------------------
      subroutine setuprestartiovalues3(imax, jmax, kmax, nharms, &
                 amount1)
!     This subroutine is used to calculate the amount of data written
!     to the restart file in parts of the I/O code.  It is necessary
!     to avoid the problems that occur when one component of the data
!     to be written out is of zero size.  Without this code that would
!     cause no data to be written at all.  This code fixes that to
!     ensure that having part of the data being zero extend does not
!     restrict the other data being written out.
!     This particular routine is setup to calculate the amounts of data
!     for the if(kom ....) loop in the writerestartblock subroutine.
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) :: imax, jmax, kmax, nharms
      integer(kind=cosa_int) :: amount1, expr1, expr2, expr3, expr4

      expr4 = (2*nharms) + 1

!     mpi  Workout the amount for the first write line i.e.
!          (imax+2*jmax-1*npde)
!     mpi  If all components are zero then it is correct to make the amount
!     mpi  written zero.
      if(((imax+3) .eq. 3) .AND. ((jmax+3) .eq. 3) .AND. &
         ((kmax+3) .eq. 3)) then
         amount1 = 0
      else
!     mpi  All components are not zero so go through them individually and
!     mpi  fix them if individual ones are zero (i.e. set them to 1 so they
!     mpi  don't affect the final multiplication)
         if((imax+3) .eq. 3) then
            expr1 = 1
         else
            expr1 = imax+3
         end if

         if((jmax+3) .eq. 3) then
            expr2 = 1
         else
            expr2 = jmax+3
         end if

         if((kmax+3) .eq. 3) then
            expr3 = 1
         else
            expr3 = kmax+3
         end if

!     mpi     Work out the total amount of data to be written
         amount1 = expr1*expr2*expr3*expr4

      end if

      return
      end

!-----------------------------------------------------------------------
      subroutine unformattedmovefilepointer(offset, fileid)
!     This routine moves offset lines through fileid 
!     It is used to find a position within a file for 
!     a process to read from.
!     It assumes that the file has been opened already using the fileid
!     provided above as the unit.  This particular routine also assumes 
!     that the file has been opened as 'unformatted'.  For formatted 
!     I/O use formattedmovefilepointer instead.
!-----------------------------------------------------------------------
      
      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none
      
      integer*8 offset
      integer(kind=cosa_int) :: fileid, i


      do i=1,offset
         read(fileid)
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine getownerid(blockid, ownerid)
!     Return the MPI rank of the process that owns a block.  This is
!     useful for working out which processes we should send data to.
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) :: blockid,ownerid

      ownerid = blockassignment(blockid)

      return 
      end

!-----------------------------------------------------------------------
      subroutine receiveblockdata(recvdata,receivefrom,length,request)
!     This routine is used to receive a vector of reals from another 
!     process. receivefrom is an integer specifying the block this 
!     data is coming from.  With this information we can look up which 
!     process owns than block and know who we will be receiving from.                     
!     We are using non-blocking communications so there is a request              
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) :: receivefrom,length,request
      integer(kind=cosa_int) :: receivefromid
      real(kind=cosa_real) :: recvdata(length)

      call getownerid(receivefrom,receivefromid)

      call mpireceive(receivefromid,recvdata,length,request)

      return 
      end

!-----------------------------------------------------------------------
      subroutine receiveblockdata_fromid(recvdata,receivefromid,length, &
                 request)
!     This routine is used to receive a vector of reals from another 
!     process. receivefrom is an integer specifying the block this 
!     data is coming from.  With this information we can look up which 
!     process owns than block and know who we will be receiving from.                     
!     We are using non-blocking communications so there is a request              
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) :: length,request
      integer(kind=cosa_int) :: receivefromid
      real(kind=cosa_real) :: recvdata(length)

      call mpireceive(receivefromid,recvdata,length,request)

      return 
      end

!-----------------------------------------------------------------------
      subroutine sendblockdata(senddata,sendto,length,request)
!     This routine is used to send vectors of reals to another 
!     process.  The input argument sendto is an integer specifying 
!     the block this data is from.  With this information we can look up 
!     which process owns than block and know who to send the data too. 
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) :: sendto,length
      integer(kind=cosa_int) :: sendtoid,request
      real(kind=cosa_real) :: senddata(length)
      
      call getownerid(sendto,sendtoid)

      call mpisend(sendtoid,senddata,length,request)
      
      return 
      end

!-----------------------------------------------------------------------
      subroutine sendblockdata_toid(senddata,sendtoid,length,request)
!     This routine is used to send vectors of reals to another 
!     process.  The input argument sendto is an integer specifying 
!     the block this data is from.  With this information we can look up 
!     which process owns than block and know who to send the data too. 
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) :: length
      integer(kind=cosa_int) :: sendtoid,request
      real(kind=cosa_real) :: senddata(length)
      
      call mpisend(sendtoid,senddata,length,request)
      
      return 
      end

!-----------------------------------------------------------------------
      subroutine mpisendwt(destrank,senddata,length,request,tag)
!     This wraps the MPI library function MPI_Isend.
!     Currently it is designed to send an array of reals using the 
!     MPI_COMM_WORLD communicator
!-----------------------------------------------------------------------

      use cosa_precision

      implicit none

#ifdef MPI
      include 'mpif.h'
#endif
      
      integer(kind=cosa_int) :: destrank,length, ierror, tag
      integer(kind=cosa_int) :: request
      real(kind=cosa_real) :: senddata(length)

#ifdef MPI
      call MPI_ISEND(senddata, length, MPI_DOUBLE_PRECISION, destrank, &
           tag, MPI_COMM_WORLD, request, ierror)
!      call MPI_SSEND(senddata, length, MPI_DOUBLE_PRECISION, destrank, 
!     &     tag, MPI_COMM_WORLD, ierror)

#endif


      return 
      end

!-----------------------------------------------------------------------
      subroutine mpisend(destrank,senddata,length,request)
!     This wraps the MPI library function MPI_Isend.
!     Currently it is designed to send an array of reals using the 
!     MPI_COMM_WORLD communicator
!-----------------------------------------------------------------------

      use cosa_precision

      implicit none

#ifdef MPI
      include 'mpif.h'
#endif
      
      integer(kind=cosa_int) :: destrank,length, ierror, tag
      integer(kind=cosa_int) :: request
      real(kind=cosa_real) :: senddata(length)


#ifdef MPI
      tag = 10

      call mpisendwt(destrank,senddata,length,request,tag)

#endif


      return 
      end

!-----------------------------------------------------------------------
      subroutine mpireceive(sourcerank,recvdata,length,request)
!     This wraps the MPI library function MPI_Irecv.
!     Currently it is designed to receive an array of reals using the 
!     MPI_COMM_WORLD communicator.
!-----------------------------------------------------------------------

      use cosa_precision

      implicit none

#ifdef MPI
      include 'mpif.h'
#endif
      
      integer(kind=cosa_int) :: sourcerank,length, ierror, tag
      integer(kind=cosa_int) :: request
      real(kind=cosa_real) :: recvdata(length)

#ifdef MPI
      tag = 10

      call mpireceivewt(sourcerank,recvdata,length,request,tag)
#endif


      return 
      end

!-----------------------------------------------------------------------
      subroutine mpireceivewt(sourcerank,recvdata,length,request,tag)
!     This wraps the MPI library function MPI_Irecv.
!     Currently it is designed to receive an array of reals using the 
!     MPI_COMM_WORLD communicator.
!-----------------------------------------------------------------------

      use cosa_precision

      implicit none

#ifdef MPI
      include 'mpif.h'
#endif
      
      integer(kind=cosa_int) :: sourcerank,length, ierror, tag
      integer(kind=cosa_int) :: request
      real(kind=cosa_real) :: recvdata(length)

#ifdef MPI

      call MPI_IRECV(recvdata, length, MPI_DOUBLE_PRECISION, sourcerank, &
           tag, MPI_COMM_WORLD, request, ierror)
      if(ierror /= MPI_SUCCESS) then
        write(0,*) 'MPI_RECV err: ',ierror,'; len: ',length,' req: ',request
      endif

!      call MPI_RECV(recvdata, length, MPI_DOUBLE_PRECISION, sourcerank, 
!     &     tag, MPI_COMM_WORLD,status, ierror)

#endif


      return 
      end

!-----------------------------------------------------------------------
      subroutine initialiserequests(requests,length)
!     Initialised an array of requests to null
!-----------------------------------------------------------------------

      use cosa_precision

      implicit none

#ifdef MPI
      include 'mpif.h'
#endif
      
      integer(kind=cosa_int) :: length
      integer(kind=cosa_int) :: requests(length)
#ifdef MPI
      requests = MPI_REQUEST_NULL
#endif

      return 
      end

!-----------------------------------------------------------------------
      subroutine waitonmessage(request)
!     Provides an external interface to the wait request for 
!     nonblock communications.
!-----------------------------------------------------------------------

      use cosa_precision

      implicit none

#ifdef MPI
      include 'mpif.h'
#endif
      
      integer(kind=cosa_int) :: request
#ifdef MPI
      integer(kind=cosa_int) :: ierror
      integer(kind=cosa_int) :: status(MPI_STATUS_SIZE)

      call mpiwait(request)
#endif

      return 
      end

!-----------------------------------------------------------------------
      subroutine waitallmessages(requests,requestnumbers)
!     Provides an external interface to the waitall request for 
!     nonblock communications.
!-----------------------------------------------------------------------

      use cosa_precision

      implicit none

#ifdef MPI
      include 'mpif.h'
#endif
      
      integer(kind=cosa_int) :: requestnumbers
      integer(kind=cosa_int) :: requests(requestnumbers)

#ifdef MPI
      call mpiwaitall(requests,requestnumbers)
#endif


      return 
      end

!-----------------------------------------------------------------------
      subroutine waitanymessages(requests,requestnumbers,mindex)
!     Provides an external interface to the waitall request for 
!     nonblock communications.
!-----------------------------------------------------------------------

      use cosa_precision

      implicit none

#ifdef MPI
      include 'mpif.h'
#endif
      
      integer(kind=cosa_int) :: requestnumbers
      integer(kind=cosa_int) :: requests(requestnumbers)
      integer(kind=cosa_int) :: mindex

#ifdef MPI
      call mpiwaitany(requests,requestnumbers,mindex)
#endif


      return 
      end

!-----------------------------------------------------------------------
      subroutine mpiwait(request)
!     This wraps the MPI library function MPI_Wait.
!-----------------------------------------------------------------------

      use cosa_precision

      implicit none

#ifdef MPI
      include 'mpif.h'
#endif
      
      integer(kind=cosa_int) :: request
#ifdef MPI
      integer(kind=cosa_int) :: ierror
      integer(kind=cosa_int) :: status(MPI_STATUS_SIZE)

      call MPI_WAIT(request, status, ierror)
#endif

      return 
      end

!-----------------------------------------------------------------------
      subroutine mpiwaitall(requests,requestnumber)
!     This wraps the MPI library function MPI_Waitall.
!-----------------------------------------------------------------------

      use cosa_precision

      implicit none

#ifdef MPI
      include 'mpif.h'
#endif

      integer(kind=cosa_int) :: requestnumber
      integer(kind=cosa_int) :: requests(requestnumber)
#ifdef MPI
      integer(kind=cosa_int) :: ierror
      integer(kind=cosa_int) :: statuses(MPI_STATUS_SIZE,requestnumber)
      call MPI_WAITALL(requestnumber,requests,statuses,ierror)
#endif


      return 
      end

!-----------------------------------------------------------------------
      subroutine mpiwaitany(requests,requestnumber,mindex)
!     This wraps the MPI library function MPI_Waitany.
!-----------------------------------------------------------------------

      use cosa_precision

      implicit none

#ifdef MPI
      include 'mpif.h'
#endif

      integer(kind=cosa_int) :: requestnumber
      integer(kind=cosa_int) :: requests(requestnumber)
      integer(kind=cosa_int) :: mindex
#ifdef MPI
      integer(kind=cosa_int) :: ierror
      integer(kind=cosa_int) ::statuses(MPI_STATUS_SIZE)
!      ierror = 0
!      statuses = 0
!      mindex = 0
      call MPI_WAITANY(requestnumber,requests,mindex,statuses,ierror)
#endif

      return 
      end

!-----------------------------------------------------------------------
      subroutine copymutreceivedata(imax1,jmax1,kmax1,ii1,jj1,kk1,dat, &
                 lreceivearray,i_req,i_proc)
!-----------------------------------------------------------------------
!     This routine is used in cutman_mut to place data received 
!     from other processes into the correct place in the dat array
!     on this process.  It mimics the serial copy functionality of
!     cut_mut.
!-----------------------------------------------------------------------

      use parallelutils
      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none
      
      integer(kind=cosa_int) :: imax1,jmax1,kmax1,ii1,jj1,kk1,receiverequestnum
      integer(kind=cosa_int) ::n,ipde,tempindex,i_req,i_proc
      real(kind=cosa_real) :: &
        dat(-1:imax1+1,-1:jmax1+1,-1:kmax1+1,1,0:2*nharms)
      real(kind=cosa_real) :: &
        lreceivearray(1*((2*nharms)+1),mrequest_perproc,maxproc_cutman)

      tempindex = 1
      do n = 0, 2*nharms
        do ipde = 1, 1
          dat(ii1,jj1,kk1,ipde,n) = lreceivearray(tempindex,i_req,i_proc)
          tempindex = tempindex + 1
        end do
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine copypmutreceivedata(imax1,jmax1,kmax1,ii1,jj1,kk1,dat, &
                 lreceivearray,i_req,i_proc)
!-----------------------------------------------------------------------
!     This routine is used in pcutman_mut to place data received 
!     from other processes into the correct place in the dat array
!     on this process.  It mimics the serial copy functionality of
!     pcut_mut.
!-----------------------------------------------------------------------

      use parallelutils
      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none
      
      integer(kind=cosa_int) :: imax1,jmax1,kmax1,ii1,jj1,kk1,receiverequestnum
      integer(kind=cosa_int) :: n,ipde,tempindex,i_req,i_proc
      real(kind=cosa_real) :: &
        dat(-1:imax1+1,-1:jmax1+1,-1:kmax1+1,1,0:2*nharms)
      real(kind=cosa_real) :: &
        lreceivearray(1*((2*nharms)+1),mrequest_pperproc,maxproc_pcutman)

      tempindex = 1
      do n = 0, 2*nharms
        do ipde = 1, 1
          dat(ii1,jj1,kk1,ipde,n) = lreceivearray(tempindex,i_req,i_proc)
          tempindex = tempindex + 1
        end do
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine copyqreceivedata(imax1,jmax1,kmax1,ii1,jj1,kk1,dat, &
                 lreceivearray,i_req,i_proc)
!-----------------------------------------------------------------------
!     This routine is used in cutman_q to place data received 
!     from other processes into the correct place in the dat array
!     on this process.  It mimics the serial copy functionality of
!     cut_q.
!-----------------------------------------------------------------------

      use parallelutils
      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none
      
      integer(kind=cosa_int) :: imax1,jmax1,kmax1,ii1,jj1,kk1
      integer(kind=cosa_int) :: n,ipde,tempindex,i_req,i_proc
      real(kind=cosa_real) :: dat(-1:imax1+1,-1:jmax1+1,-1:kmax1+1,npde,0:2*nharms)
      real(kind=cosa_real) :: lreceivearray(npde*((2*nharms)+1),mrequest_perproc, &
                   maxproc_cutman)

      tempindex = 1
      do n = 0, 2*nharms
        do ipde = 1, npde
          dat(ii1,jj1,kk1,ipde,n) = lreceivearray(tempindex,i_req,i_proc)
          tempindex = tempindex + 1
        end do
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine copypqreceivedata(imax1,jmax1,kmax1,ii1,jj1,kk1,dat, &
                 lreceivearray,i_req,i_proc)
!-----------------------------------------------------------------------
!     This routine is used in pcutman_q to place data received 
!     from other processes into the correct place in the dat array
!     on this process.  It mimics the serial copy functionality of
!     pcut_q.
!-----------------------------------------------------------------------

      use parallelutils
      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none
      
      integer(kind=cosa_int) :: imax1,jmax1,kmax1,ii1,jj1,kk1
      integer(kind=cosa_int) :: n,ipde,tempindex,i_req,i_proc
      real(kind=cosa_real) :: dat(-1:imax1+1,-1:jmax1+1,-1:kmax1+1,npde,0:2*nharms)
      real(kind=cosa_real) :: lreceivearray(npde*((2*nharms)+1),mrequest_pperproc, &
                   maxproc_pcutman)

      tempindex = 1
      do n = 0, 2*nharms
        do ipde = 1, npde
          dat(ii1,jj1,kk1,ipde,n) = lreceivearray(tempindex,i_req,i_proc)
          tempindex = tempindex + 1
        end do
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine copypqhbreceivedata(imax1,jmax1,kmax1,ii1,jj1,kk1,dat, &
                 lreceivearray,i_req,i_proc)
!-----------------------------------------------------------------------
!     This routine is used in pcutman_q to place data received 
!     from other processes into the correct place in the dat array
!     on this process.  It mimics the serial copy functionality of
!     hb branch of pcut_q.
!-----------------------------------------------------------------------

      use parallelutils
      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none
      
      integer(kind=cosa_int) :: imax1,jmax1,kmax1,ii1,jj1,kk1
      integer(kind=cosa_int) :: n,n1,ipde,tempindex,i_req,i_proc
      real(kind=cosa_real) :: dat(-1:imax1+1,-1:jmax1+1,-1:kmax1+1,npde,0:2*nharms)
      real(kind=cosa_real) :: lreceivearray(npde*((2*nharms)+1),mrequest_pperproc, &
                   maxproc_pcutman), tmp(7,0:2*mharms)

      tempindex = 1
      do n = 0, 2*nharms
        do ipde = 1, npde
          tmp(ipde,n) = lreceivearray(tempindex,i_req,i_proc)
          tempindex = tempindex + 1
        end do
      end do
!------------ from frequency-domain to time-domain
      do n = 0, 2*nharms
        do ipde=1,npde
          dat(ii1,jj1,kk1,ipde,n) = 0.d0
          do n1 = 0, 2*nharms
            dat(ii1,jj1,kk1,ipde,n) = dat(ii1,jj1,kk1,ipde,n) + &
                                      eihb(n,n1) * tmp(ipde,n1)
          end do
        end do
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine copyxreceivedata(imax1,jmax1,kmax1,ii1,jj1,kk1,dat, &
                 lreceivearray,i_req,i_proc)
!-----------------------------------------------------------------------
!     This routine is used in cutman_x to place data received 
!     from other processes into the correct place in the dat array
!     on this process.  It mimics the serial copy functionality of
!     cut_s.
!-----------------------------------------------------------------------

      use parallelutils
      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) :: imax1,jmax1,kmax1,ii1,jj1,kk1,i_req,i_proc
      integer(kind=cosa_int) :: n
      real(kind=cosa_real) :: dat(0:imax1+1,0:jmax1+1,0:kmax1+1,0:2*nharms*hbmove)
      real(kind=cosa_real) :: lreceivearray(2*nharms*hbmove+1,mrequest_perproc, &
                   maxproc_cutman)

      do n = 0, 2*nharms*hbmove
        dat(ii1,jj1,kk1,n) = lreceivearray(n+1,i_req,i_proc)
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine copypxreceivedata(imax1,jmax1,kmax1,ii1,jj1,kk1,dat, &
                 lreceivearray,i_req,i_proc)
!-----------------------------------------------------------------------
!     This routine is used in pcutman_x to place data received 
!     from other processes into the correct place in the dat array
!     on this process.  It mimics the serial copy functionality of
!     pcut_x.
!-----------------------------------------------------------------------

      use parallelutils
      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none
      
      integer(kind=cosa_int) :: imax1,jmax1,kmax1,ii1,jj1,kk1,i_req,i_proc
      integer(kind=cosa_int) :: n
      real(kind=cosa_real) :: dat(0:imax1+1,0:jmax1+1,0:kmax1+1,0:2*nharms*hbmove)
      real(kind=cosa_real) :: lreceivearray(2*nharms*hbmove+1,mrequest_pperproc, &
                   maxproc_pcutman)

      do n = 0, 2*nharms*hbmove
        dat(ii1,jj1,kk1,n) = lreceivearray(n+1,i_req,i_proc)
      end do

      return
      end

!-----------------------------------------------------------------------
      subroutine copyvreceivedata(imax1,jmax1,kmax1,ii1,jj1,kk1,dat, &
                 lreceivearray,i_req,i_proc)
!-----------------------------------------------------------------------
!     This routine is used in cutman_v to place data received 
!     from other processes into the correct place in the dat array
!     on this process.  It mimics the serial copy functionality of
!     cut_v.
!-----------------------------------------------------------------------

      use parallelutils
      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none
      
      integer(kind=cosa_int) :: imax1,jmax1,kmax1,ii1,jj1,kk1,i_req,i_proc
      real(kind=cosa_real) :: dat(0:imax1,0:jmax1,0:kmax1)
      real(kind=cosa_real) :: lreceivearray(mrequest_perproc,maxproc_cutman)

      dat(ii1,jj1,kk1) = lreceivearray(i_req,i_proc)

      return
      end

!-----------------------------------------------------------------------
      subroutine copypvreceivedata(imax1,jmax1,kmax1,ii1,jj1,kk1,dat, &
                 lreceivearray,i_req,i_proc)
!-----------------------------------------------------------------------
!     This routine is used in pcutman_v to place data received 
!     from other processes into the correct place in the dat array
!     on this process.  It mimics the serial copy functionality of
!     pcut_v.
!-----------------------------------------------------------------------

      use parallelutils
      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none
      
      integer(kind=cosa_int) :: imax1,jmax1,kmax1,ii1,jj1,kk1,i_req,i_proc
      real(kind=cosa_real) :: dat(0:imax1,0:jmax1,0:kmax1)
      real(kind=cosa_real) :: lreceivearray(mrequest_pperproc,maxproc_pcutman)

      dat(ii1,jj1,kk1) = lreceivearray(i_req,i_proc)

      return
      end

!-----------------------------------------------------------------------
      subroutine combineforces(cl,cd,cm,n)
!     combineforces is used to take the data in the cl,cd, and cm arrays
!     on each process and combine them so all processes have the full 
!     forces calculated across the whole simulation.  It uses an 
!     allreduce to do the communication and combining of the data.
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none
      
      integer(kind=cosa_int) :: n,i,j,datalength
      real(kind=cosa_real) :: cl(0:2*mharms,msurface),cd(0:2*mharms,msurface), &
                   cm(0:2*mharms,msurface),temparray(3)
      
!mpi  This uses more allreduce than is strictly necessary so is a place 
!mpi  that could be optimised (i.e. with some data packing the 
!mpi  do i=1,nsurface loop could be removed and a single all reduce 
!mpi  used to combine the data)
      do i=1,nsurface
         temparray(1) = cl(n,i)
         temparray(2) = cd(n,i)
         datalength = 2
!         if(functag.eq.3) then
            temparray(3) = cm(n,i)
            datalength = 3
!         end if
         call realsumallreduce(temparray,datalength)
         cl(n,i) = temparray(1)
         cd(n,i) = temparray(2)
!         if(functag.eq.3) then
            cm(n,i) = temparray(3)
!         end if
      end do

      return 
      end

!-----------------------------------------------------------------------
      subroutine newcombineforces(cl,cd,cm)
!     combineforces is used to take the data in the cl,cd, and cm arrays
!     on each process and combine them so all processes have the full 
!     forces calculated across the whole simulation.  It uses an 
!     allreduce to do the communication and combining of the data.
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none
      
      integer(kind=cosa_int)::  n,i,j,k,datalength
      real(kind=cosa_real) :: cl(0:2*mharms,msurface),cd(0:2*mharms,msurface), &
                   cm(0:2*mharms,msurface)
      real(kind=cosa_real) :: temparray(((2*nharms)+1)*nsurface*3)
      
      j = 1
      do k = 0,2*nharms
         do i=1,nsurface
            temparray(j) = cl(k,i)
            j = j +1
            temparray(j) = cd(k,i)
            j = j + 1
            temparray(j) = cm(k,i)
            j = j + 1
         end do
      end do
      call realsumallreduce(temparray,j-1)
      j =  1
      do k = 0,2*nharms
         do i=1,nsurface
            cl(k,i) = temparray(j)
            j = j +1
            cd(k,i) = temparray(j)
            j = j + 1
            cm(k,i) = temparray(j)
            j = j + 1
         end do
      end do

      return 
      end

!-----------------------------------------------------------------------
      subroutine combinemassflow(m_in,m_out,datalength)
!     combinemassflow combines the data in each local processes
!     m_in and m_out arrays to enable each process to know what 
!     the full massflow of the simulation is.  It packs the data into 
!     temporary arrays and then does 2 all reduces (one for each array)
!     to combine the results.
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none
      
      integer(kind=cosa_int) :: datalength,n
      real(kind=cosa_real) :: m_in(0:80),m_out(0:80)
      real(kind=cosa_real) :: tempm_in(81),tempm_out(81)

      do n=1,datalength
         tempm_in(n) = m_in(n-1)
         tempm_out(n) = m_out(n-1)
      end do

      call realsumallreduce(tempm_in,datalength)
      call realsumallreduce(tempm_out,datalength)

      do n=1,datalength
         m_in(n-1) = tempm_in(n)
         m_out(n-1) = tempm_out(n)
      end do

      return 
      end

!-----------------------------------------------------------------------
      subroutine combineresults(resrmslocal,lengthresrms,ncell)
!     combineresults is used to calculate the overall resrms data 
!     for the whole simulation using MPI all reduces to perform 
!     the combination of the data.
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none
      
      integer(kind=cosa_int) :: ncell, ncelldata(1), lengthresrms
      real(kind=cosa_real) :: resrmslocal(lengthresrms)

      ncelldata(1) = ncell
      call realsumallreduce(resrmslocal,lengthresrms)
      call integersumallreduce(ncelldata,1)
      ncell = ncelldata(1)

      return 
      end

!-----------------------------------------------------------------------
      subroutine parallelmergewall(nl,xwall,ywall,zwall,xgwall,ygwall, &
                                   zgwall)
!     parallelmergewall is used to take each local xwall, ywall and
!     for each block a process owns and then communicate them so
!     at the end each process has all the local xwall, ywall and zwall
!     data in xgwall, ygwall and zgwall. This is needed to enable each 
!     process to correctly calculate the dist2wall for it's portion of 
!     the simulations (this calculation relies on knowing the walls 
!     in the whole simulation).
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none
      
      integer(kind=cosa_int) :: nl,iblock,ownerid,mpiid,wallsize,myblockid,ixyzw
      integer(kind=cosa_int) :: ixyzwg,j
      real(kind=cosa_real) :: xwall(*),ywall(*),zwall(*),xgwall(*),ygwall(*), &
                   zgwall(*)
      real(kind=cosa_real), allocatable:: tempwall(:)

      call getmpiid(mpiid)

      do iblock = 1,nblocks
         call getownerid(iblock,ownerid)
         wallsize = g_n_wall(iblock,nl)
         ixyzwg = off_1dwg(iblock,nl)
         allocate(tempwall(wallsize))

         if(ownerid .eq. mpiid) then
            myblockid = iblock - (lowernblock - 1)
            ixyzw = 1 + off_1dw(myblockid,nl)
            do j = 1,wallsize
               tempwall(j) = xwall(ixyzw + j - 1)
            end do
         end if

         call realbroadcast(ownerid,wallsize,tempwall)

         if(ownerid .eq. mpiid) then
            do j = 1,wallsize
               xgwall(ixyzwg + j) = xwall(ixyzw - 1 + j)
            end do
         else
            do j = 1,wallsize
               xgwall(ixyzwg + j) = tempwall(j)
            end do
         end if

         if(ownerid .eq. mpiid) then
            do j = 1,wallsize
               tempwall(j) = ywall(ixyzw + j - 1)
            end do
         end if

         call realbroadcast(ownerid,wallsize,tempwall)

         if(ownerid .eq. mpiid) then
            do j = 1,wallsize
               ygwall(ixyzwg + j) = ywall(ixyzw - 1 + j)
            end do
         else
            do j = 1,wallsize
               ygwall(ixyzwg + j) = tempwall(j)
            end do
         end if

         if(ownerid .eq. mpiid) then
            do j = 1,wallsize
               tempwall(j) = zwall(ixyzw + j - 1)
            end do
         end if

         call realbroadcast(ownerid,wallsize,tempwall)

         if(ownerid .eq. mpiid) then
            do j = 1,wallsize
               zgwall(ixyzwg + j) = zwall(ixyzw - 1 + j)
            end do
         else
            do j = 1,wallsize
               zgwall(ixyzwg + j) = tempwall(j)
            end do
         end if

         deallocate(tempwall)
      end do

      return 
      end

!-----------------------------------------------------------------------
      subroutine integersumallreduce(reducedata,length)
!     This wraps the MPI library function MPI_Allreduce.
!     This particular wrapper is designed to allow arrays 
!     of integers to be summed.
!-----------------------------------------------------------------------

      use cosa_precision

      implicit none

#ifdef MPI
      include 'mpif.h'
#endif

      integer(kind=cosa_int) :: length, ierror, i      
      integer(kind=cosa_int) :: reducedata(length), tempdata(length)

#ifdef MPI
!      call MPI_ALLREDUCE(MPI_IN_PLACE, reducedata, length, MPI_INTEGER, 
!     &     MPI_SUM, MPI_COMM_WORLD, ierror)
      call MPI_ALLREDUCE(reducedata, tempdata, length, MPI_INTEGER, &
           MPI_SUM, MPI_COMM_WORLD, ierror)
      do i=1,length
         reducedata(i) = tempdata(i)
      end do
#endif


      return 
      end

!-----------------------------------------------------------------------
      subroutine realsumallreduce(reducedata,length)
!     This wraps the MPI library function MPI_Allreduce.
!     This particular wrapper is designed to allow arrays 
!     of reals to be summed.
!-----------------------------------------------------------------------

      use cosa_precision

      implicit none

#ifdef MPI
      include 'mpif.h'
#endif
      
      integer(kind=cosa_int) :: length, ierror, i
      real(kind=cosa_real) :: reducedata(length), tempdata(length)

#ifdef MPI
!      call MPI_ALLREDUCE(MPI_IN_PLACE, reducedata, length, 
!     &     MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror)
      call MPI_ALLREDUCE(reducedata, tempdata, length, &
           MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror)
      do i=1,length
         reducedata(i) = tempdata(i) 
      end do
#endif

      return 
      end

!-----------------------------------------------------------------------
      subroutine realsumreduce(reducedata,length)
!     This wraps the MPI library function MPI_Reduce.
!     This particular wrapper is designed to allow arrays 
!     of reals to be summed with the answer being constructed on rank 0.
!-----------------------------------------------------------------------

      use cosa_precision

      implicit none

#ifdef MPI
      include 'mpif.h'
#endif
      
      integer(kind=cosa_int) :: length, ierror, i
      real(kind=cosa_real) :: reducedata(length), tempdata(length)

#ifdef MPI
!      call MPI_REDUCE(MPI_IN_PLACE, reducedata, length, 
!     &     MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror)
      call MPI_REDUCE(reducedata, tempdata, length, &
           MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror)
      do i=1,length
         reducedata(i) = tempdata(i) 
      end do
#endif

      return 
      end

!-----------------------------------------------------------------------
      subroutine realmaxreduce(reducedata,length)
!     This wraps the MPI library function MPI_Reduce.
!     This particular wrapper is designed to allow the maximum of each
!     element in an array of reals to be selected with the answer being 
!     constructed on rank 0.
!-----------------------------------------------------------------------

      use cosa_precision

      implicit none

#ifdef MPI
      include 'mpif.h'
#endif
      
      integer(kind=cosa_int) :: length, ierror, i
      real(kind=cosa_real) :: reducedata(length), tempdata(length)

#ifdef MPI
      call MPI_REDUCE(reducedata, tempdata, length, &
           MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierror)
      do i=1,length
         reducedata(i) = tempdata(i)
      end do
#endif

      return 
      end

!-----------------------------------------------------------------------
      subroutine realminreduce(reducedata,length)
!     This wraps the MPI library function MPI_Reduce.
!     This particular wrapper is designed to allow the minimum of each 
!     element in an array of reals to be selected with the answer being 
!     constructed on rank 0.
!-----------------------------------------------------------------------

      use cosa_precision

      implicit none

#ifdef MPI
      include 'mpif.h'
#endif
      
      integer(kind=cosa_int) :: length, ierror, i
      real(kind=cosa_real) :: reducedata(length), tempdata(length)

#ifdef MPI
      call MPI_REDUCE(reducedata, tempdata, length, &
           MPI_DOUBLE_PRECISION, MPI_MIN, 0, MPI_COMM_WORLD, ierror)
      do i=1,length
         reducedata(i) = tempdata(i)
      end do
#endif

      return 
      end

!-----------------------------------------------------------------------
      subroutine realmaxlocallreduce(reducedata,length,rank,loc)
!     This wraps the MPI library function MPI_Allreduce.
!     This particular wrapper is designed to allow the maximum value in an
!     array of reals to be discovered an also returns the rank of the 
!     process that owns that value.  The return values are valid on
!     all ranks.
!-----------------------------------------------------------------------

      use cosa_precision

      implicit none

#ifdef MPI
      include 'mpif.h'
#endif
      
      integer(kind=cosa_int) :: length, ierror, i, rank
      integer(kind=cosa_int) :: loc(length)
      real(kind=cosa_real) :: reducedata(length), tempdata(2,length), tempdata_out(2,length)
#ifdef MPI
      do i=1,length
         tempdata(1,i)=reducedata(i)
         tempdata(2,i)=rank
      end do
!      call MPI_ALLREDUCE(MPI_IN_PLACE, tempdata, length, 
!     &     MPI_2DOUBLE_PRECISION, MPI_MAXLOC, MPI_COMM_WORLD, ierror)
      call MPI_ALLREDUCE(tempdata, tempdata_out, length, &
           MPI_2DOUBLE_PRECISION, MPI_MAXLOC, MPI_COMM_WORLD, ierror)
      do i=1,length
         reducedata(i) = tempdata_out(1,i)
         loc(i) = tempdata_out(2,i)
      end do
#endif

      return 
      end

!-----------------------------------------------------------------------
      subroutine realbroadcast(rootid,length,broadcastdata)
!     This wraps the MPI library function MPI_Broadcast.
!     This particular wrapper is designed to allow arrays 
!     of reals to be broadcast from a specified root process 
!     (specified by the the rootid parameter).
!-----------------------------------------------------------------------

      use cosa_precision

      implicit none

#ifdef MPI
      include 'mpif.h'
#endif
      
      integer(kind=cosa_int) :: length, ierror, rootid
      real(kind=cosa_real) :: broadcastdata(length)

#ifdef MPI
      call MPI_BCAST(broadcastdata, length, &
           MPI_DOUBLE_PRECISION, rootid, MPI_COMM_WORLD, ierror)
#endif

      return 
      end

!-----------------------------------------------------------------------
      subroutine integerbroadcast(rootid,length,broadcastdata)
!     This wraps the MPI library function MPI_Broadcast.
!     This particular wrapper is designed to allow arrays 
!     of reals to be broadcast from a specified root process 
!     (specified by the the rootid parameter).
!-----------------------------------------------------------------------

      use cosa_precision

      implicit none

#ifdef MPI
      include 'mpif.h'
#endif
      
      integer(kind=cosa_int) :: length, ierror, rootid
      integer(kind=cosa_int) :: broadcastdata(length)

#ifdef MPI
      call MPI_BCAST(broadcastdata, length, &
           MPI_INTEGER, rootid, MPI_COMM_WORLD, ierror)
#endif

      return 
      end

!-----------------------------------------------------------------------
      subroutine getreadonlyiomode(iomode)
!     Sets the I/O mode for the MPI I/O files to be read.  This will 
!     cause an error if the file does not exist.
!-----------------------------------------------------------------------

      use cosa_precision

      implicit none

#ifdef MPI
      include 'mpif.h'
#endif

      integer(kind=cosa_int) :: iomode
      
#ifdef MPI
      iomode = MPI_MODE_RDONLY
#endif

      return
      end

!-----------------------------------------------------------------------
      subroutine getcreatewriteiomode(iomode)
!     Sets the I/O mode for the MPI I/O files to be read/write and 
!     create (i.e. make a file if it does not exist and allow 
!     read-write access).
!-----------------------------------------------------------------------

      use cosa_precision

      implicit none

#ifdef MPI
      include 'mpif.h'
#endif

      integer(kind=cosa_int) :: iomode
      
#ifdef MPI
      iomode = ior(MPI_MODE_WRONLY,MPI_MODE_CREATE)
#endif

      return
      end

!-----------------------------------------------------------------------
      subroutine getappendiomode(iomode)
!     Sets the I/O mode for MPI I/O file to read/write.  This will 
!     append any write data to the end of the file and also may 
!     cause an error if a file does not exist in the filesystem and 
!     reads or writes are attempted on it.
!-----------------------------------------------------------------------

      use cosa_precision

      implicit none

#ifdef MPI
      include 'mpif.h'
#endif

      integer(kind=cosa_int) :: iomode
      
#ifdef MPI
      iomode = MPI_MODE_WRONLY
#endif

      return
      end

!-----------------------------------------------------------------------
      subroutine openfile(filehandle,filename,filemode)
!     Wraps the MPI-I/O MPI_FILE_OPEN subroutine 
!-----------------------------------------------------------------------

      use cosa_precision

      implicit none

#ifdef MPI
      include 'mpif.h'
#endif

      integer(kind=cosa_int) :: filehandle, filemode, ierror
      character*(*) filename

#ifdef MPI
      call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, filemode, &
           MPI_INFO_NULL, filehandle, ierror)
      if(ierror .ne. 0) then
         write(*,*) 'Problem opening the file ',filename
         call abortmpi()
      end if
#endif

      return
      end

!-----------------------------------------------------------------------
      subroutine deletefile(filename)
!     Wraps the MPI-I/O MPI_FILE_DELETE subroutine 
!-----------------------------------------------------------------------

      use cosa_precision

      implicit none

#ifdef MPI
      include 'mpif.h'
#endif

      integer(kind=cosa_int) :: ierror
      character*(*) filename
      logical :: amcontrol
      
#ifdef MPI

      call amcontroller(amcontrol)

!     mpi  MPI_FILE_DELETE is a non-collective operation so we only need 
!     mpi  to call it from one MPI process.  This uses the value returned 
!     mpi  from the amcontroller subroutine called above to determine who
!     mpi  runs the delete command.
      if(amcontrol) then
         call MPI_FILE_DELETE(filename,MPI_INFO_NULL,ierror)
      end if
#endif

      return
      end

!-----------------------------------------------------------------------
      subroutine closempiiofile(filehandle)
!     This routine provides file close functionality.  It wraps the 
!     MPI-I/O routine MPI_FILE_CLOSE.  It should be noted that if 
!     a file was opened with mode MPI_MODE_DELETE_ON_CLOSE then this 
!     call will also delete the file.
!-----------------------------------------------------------------------

      use cosa_precision

      implicit none

#ifdef MPI
      include 'mpif.h'
#endif

      integer(kind=cosa_int) :: filehandle, ierror

#ifdef MPI
      call MPI_FILE_CLOSE(filehandle,ierror)
#endif

      return
      end

!-----------------------------------------------------------------------
      subroutine setupfile(filehandle, displacement, type)
!     setupfile is designed to wrap the MPI-I/O function/subroutine
!     MPI_FILE_SET_VIEW which sets the position for this process
!     the the file specified.  This will then allow each process to 
!     write in different places in a file.
!-----------------------------------------------------------------------

      use cosa_precision

      implicit none

#ifdef MPI
      include 'mpif.h'
      integer(kind=cosa_int) :: offset_kind
      parameter(offset_kind=MPI_OFFSET_KIND)
#else
      integer(kind=cosa_int) :: offset_kind
      parameter(offset_kind=8)
#endif

      integer(kind=cosa_int) :: filehandle, ierror
      integer*8 displacement
      integer(kind=cosa_int) :: type
#ifdef MPI
      integer(kind=offset_kind) :: disp

      disp = int(displacement, kind=offset_kind)

      call MPI_FILE_SEEK(filehandle, disp, MPI_SEEK_SET, &
           ierror)

#endif

      return
      end

!-----------------------------------------------------------------------
      subroutine parallelreadgrid(x,y,z,xyz,dx,dy,dz,dxyz,nl)
!     This routine reads the mesh/grid file in using MPI-I/O.  It is 
!     bespoke for the mesh file and uses functionality that is 
!     also used in the serial I/O code.  This means if you change 
!     any I/O functionality here (which is not just the MPI-I/O code)
!     then it should also be changed in the same code in input.f
!     (i.e. subroutine readgrid)
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

#ifdef MPI
      include 'mpif.h'
      integer(kind=cosa_int) :: offset_kind
      parameter(offset_kind=MPI_OFFSET_KIND)
#else
      integer(kind=cosa_int) :: offset_kind
      parameter(offset_kind=8)
#endif

      integer(kind=cosa_int) :: nl,imax,jmax,kmax,iblock,ixyz,jxyz,jdxyz
      integer(kind=offset_kind) :: disp
      integer(kind=cosa_int) :: fid,integersize,doublesize,ierr,iomode
      real(kind=cosa_real) :: x(*),y(*),z(*),xyz(*),dx(*),dy(*),dz(*),dxyz(*)
      logical :: amcontrol
      double precision :: starttime, endtime, totaltime, maxtime, mintime

#ifdef MPI
      call getmpitime(starttime)

      call amcontroller(amcontrol)

!     mpi  We now open the file for reading the grid
      call getreadonlyiomode(iomode)
      call openfile(fid,'mesh.dat',iomode)

      call gridbinaryoffset(disp,nl,integersize,doublesize)

      do iblock = 1,mynblocks
         imax  = i_imax(iblock,nl)
         jmax  = j_jmax(iblock,nl)
         kmax  = k_kmax(iblock,nl)
         ixyz  = 1 + off_0 (iblock,nl)
         jxyz  = 1 + off_p2(iblock,nl)
         jdxyz = 1 + off_p1(iblock,nl)
         call parallelreadgridblock(x(jxyz),y(jxyz),z(jxyz),xyz(ixyz), &
              dx(jdxyz),dy(jdxyz),dz(jdxyz),dxyz(ixyz),imax,jmax,kmax, &
              fid,disp,integersize,doublesize)
      end do

      call closempiiofile(fid)
      
      call getmpitime(endtime)

      totaltime = endtime-starttime
      maxtime = totaltime
      mintime =  totaltime

      call realmaxreduce(maxtime,1)
      call realminreduce(mintime,1)

      if(amcontrol) then
         write(*,'(A,F8.2,A,2F8.2,A)') 'Read grid time: ',totaltime,' (',maxtime,mintime,') seconds'
      end if

#endif

      return
      end

!-----------------------------------------------------------------------
      subroutine parallelreadgridblock(x,y,z,xyz,dx,dy,dz,dxyz,imax, &
           jmax,kmax,fid,disp,integersize,doublesize)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

#ifdef MPI
      include 'mpif.h'
      integer(kind=cosa_int) :: offset_kind
      parameter(offset_kind=MPI_OFFSET_KIND)
#else
      integer(kind=cosa_int) :: offset_kind
      parameter(offset_kind=8)
#endif

      integer(kind=cosa_int) :: imax,jmax,kmax
      integer(kind=cosa_int) :: i,j,k,fid,doublesize,integersize,blocksize,ierror
      integer(kind=offset_kind) :: disp
      real(kind=cosa_real) :: &
          xyz (imax,jmax,kmax,3), &
          dxyz(imax,jmax,kmax,3), &
          x   (0:imax+1,0:jmax+1,0:kmax+1), &
          y   (0:imax+1,0:jmax+1,0:kmax+1), &
          z   (0:imax+1,0:jmax+1,0:kmax+1), &
          dx  (0:imax,0:jmax,0:kmax), &
          dy  (0:imax,0:jmax,0:kmax), &
          dz  (0:imax,0:jmax,0:kmax)

#ifdef MPI
      blocksize = imax*jmax*kmax*3

      disp = disp + integersize
      call setupfile(fid,disp,MPI_DOUBLE_PRECISION)
      call MPI_FILE_READ(fid, xyz(1,1,1,1), blocksize, &
           MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierror)
      disp = disp + blocksize*doublesize + integersize

      do k=1,kmax
        do j=1,jmax
          do i=1,imax
            x(i,j,k) = xyz(i,j,k,1)
          end do
        end do
      end do

      do k=1,kmax
        do j=1,jmax
          do i=1,imax
            y(i,j,k) = xyz(i,j,k,2)
          end do
        end do
      end do

      do k=1,kmax
        do j=1,jmax
          do i=1,imax
            z(i,j,k) = xyz(i,j,k,3)
          end do
        end do
      end do

      if (moving.and.pitching) then

         disp = disp + integersize
         call setupfile(fid,disp,MPI_DOUBLE_PRECISION)
         call MPI_FILE_READ(fid, dxyz(1,1,1,1), blocksize, &
              MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierror)
         disp = disp + blocksize*doublesize + integersize

         do k=1,kmax
            do j=1,jmax
               do i=1,imax
                  dx(i,j,k) = dxyz(i,j,k,1)
               end do
            end do
         end do
         
         do k=1,kmax
            do j=1,jmax
               do i=1,imax
                  dy(i,j,k) = dxyz(i,j,k,2)
               end do
            end do
         end do
         
         do k=1,kmax
            do j=1,jmax
               do i=1,imax
                  dz(i,j,k) = dxyz(i,j,k,3)
               end do
            end do
         end do
         
      end if
#endif

      return
      end

!-----------------------------------------------------------------------
      subroutine parallelreadcgnsgrid(x,y,z,xyz,dx,dy,dz,dxyz,nl)
!-----------------------------------------------------------------------
      
      use cgns
      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none
      
#ifdef MPI
      include 'mpif.h'
      integer(kind=cosa_int) :: offset_kind
      parameter(offset_kind=MPI_OFFSET_KIND)
#else
      integer(kind=cosa_int) :: offset_kind
      parameter(offset_kind=8)
#endif
      
      integer(kind=cosa_int) :: comm_info
      integer(kind=cosa_int) :: nl,imax,jmax,kmax,iblock,ixyz,jxyz,jdxyz,fid
      integer(kind=cosa_int) :: ierr,celldim,phydim,basenum,blocknum
      real(kind=cosa_real) :: x(*),y(*),z(*),xyz(*),dx(*),dy(*),dz(*),dxyz(*)
      double precision :: starttime, endtime, totaltime, maxtime, mintime
      logical :: amcontrol

      call getmpitime(starttime)

      call amcontroller(amcontrol)

      basenum = 1

#ifdef MPI
!      call cgp_queue_set_f(1, ierr)
!      if(ierr .ne. CG_OK) then
!         write(*,*) 'queue set error'
!         call cg_error_print_f()
!      end if

!      call mpi_info_create(comm_info, ierr)

!      call mpi_info_set(comm_info, "striping_unit", "1048576", ierr)

!      call cgp_mpi_info_f(comm_info, ierr)

      call cgp_pio_mode_f(CGP_INDEPENDENT, ierr)
      if(ierr .ne. CG_OK) then
         write(*,*) 'cgp_pio_mode_f error parallelread'
         call cg_error_print_f()
         call abortmpi()
      end if

      call cgp_open_f(cgns_filename,CG_MODE_READ,fid,ierr)
      if(ierr .ne. CG_OK) then
         write(*,*) 'cgp_open_f error parallelread'
         call cg_error_print_f()
      end if
      
      do iblock = 1,mynblocks

         blocknum = iblock - 1 + lowernblock
         imax  = i_imax(iblock,nl)
         jmax  = j_jmax(iblock,nl)
         kmax  = k_kmax(iblock,nl)
         ixyz  = 1 + off_0 (iblock,nl)
         jxyz  = 1 + off_p2(iblock,nl)
         jdxyz = 1 + off_p1(iblock,nl)
         call parallel_read_bcgnsgrid(xyz(ixyz),dxyz(ixyz),imax,jmax,kmax, &
              fid,basenum,blocknum)
         call xyz_offset(x(jxyz),y(jxyz),z(jxyz),xyz(ixyz),dx(jdxyz), &
              dy(jdxyz),dz(jdxyz),dxyz(ixyz),imax,jmax,kmax)         
         
      end do           

      call closeparallelcgnsfile(fid)

      call getmpitime(endtime)

      totaltime = endtime-starttime
      maxtime = totaltime
      mintime =  totaltime

      call realmaxreduce(maxtime,1)
      call realminreduce(mintime,1)

      if(amcontrol) then
         write(*,'(A,F8.2,A,2F8.2,A)') 'Read grid time: ',totaltime,' (',maxtime,mintime,') seconds'
      end if
!      call cgp_queue_flush_f(ierr)
!      if(ierr .ne. CG_OK) then
!         write(*,*) 'queue flush error'
!         call cg_error_print_f()
!      end if
#endif   

      return
      end


!-----------------------------------------------------------------------
      subroutine closeparallelcgnsfile(filehandle)
!     This routine provides file close functionality.  It wraps the 
!     parallel cgns file clode routine. 
!-----------------------------------------------------------------------

      use cgns
      use cosa_precision

      implicit none

      integer(kind=cosa_int) :: filehandle, ierror

      call cgp_close_f(filehandle,ierror)
      if(ierror .ne. CG_OK) then
         write(*,*) 'cgp_close_f error'
         call cg_error_print_f()
      end if

      return
      end


!-----------------------------------------------------------------------
      subroutine parallel_read_bcgnsgrid(xyz,dxyz,imax,jmax,kmax,fid, &
           basenum,blocknum)
!-----------------------------------------------------------------------

      use cgns
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) :: imax,jmax,kmax
      integer(kind=cosa_int) :: i,j,k,fid,basenum,blocknum,zonenum,ierr
      integer(kind=cosa_int) :: xnum,ynum,znum,dxnum,dynum,dznum
      integer(kind=cosa_int) :: dimvec(3)
      integer(cgsize_t) :: minrange(3),maxrange(3),sizes(3,3)
      real(kind=cosa_real) :: &
          xyz (imax,jmax,kmax,3), &
          dxyz(imax,jmax,kmax,3)
      real(kind=cosa_real) :: temp(imax,jmax,kmax)
      character*32 zonename
      integer(kind=cosa_int) :: rank

      minrange(1) = 1
      minrange(2) = 1
      minrange(3) = 1

      call getmpiid(rank)

      maxrange(1) = imax
      maxrange(2) = jmax
      maxrange(3) = kmax

      zonenum = blocknum

!      call cg_index_dim_f(fid,basenum,1,isize_dimension,ierr)
!      if(ierr.ne.0) then
!         write(*,*) 'ERROR: call cg_index_dim_f'
!         call cg_error_exit_f()
!         call abortmpi()
!      end if

!      call cg_zone_read_f(fid,basenum,zonenum,zonename,sizes,ierr)
!      write(*,'(A6,I4,A5,A15,A5,9I4,A9,3I4)') 'block',zonenum,'name',zonename,'sizes',sizes,
!     &     'maxrange',maxrange

#ifdef MPI
!      dimvec(1) = 1
!      dimvec(2) = 2
!      dimvec(3) = 3
!      call cgp_coord_multi_write_data_f(fid,basenum,zonenum,dimvec,
!     &     minrange,maxrange,xyz(:,:,:,1),xyz(:,:,:,2),xyz(:,:,:,3),
!     &     ierr)

      xnum = 1
      call cgp_coord_read_data_f(fid,basenum,zonenum,xnum, &
           minrange,maxrange,xyz(:,:,:,1),ierr)
      if(ierr .ne. CG_OK) then
         write(*,*) rank,'X cgp_coord_read_data_f error',maxrange
         call cg_error_print_f()
         call abortmpi()
      end if

      ynum = 2
      call cgp_coord_read_data_f(fid,basenum,zonenum,ynum, &
           minrange,maxrange,xyz(:,:,:,2),ierr)
      if(ierr .ne. CG_OK) then
         write(*,*) rank,'Y cgp_coord_read_data_f error',maxrange
         call cg_error_print_f()
         call abortmpi()
      end if

      znum = 3
      call cgp_coord_read_data_f(fid,basenum,zonenum,znum, &
           minrange,maxrange,xyz(:,:,:,3),ierr)
      if(ierr .ne. CG_OK) then
         write(*,*) 'Z cgp_coord_read_data_f error',maxrange
         call cg_error_print_f()
         call abortmpi()
      end if

      if (moving.and.pitching) then
         dxnum = 4
         call cgp_coord_read_data_f(fid,basenum,zonenum,dxnum, &
              minrange,maxrange,dxyz(:,:,:,1),ierr)
         if(ierr .ne. CG_OK) then
            write(*,*) 'moving and pitching X cgp_coord_read_data_f error'
            call cg_error_print_f()
            call abortmpi()
         end if

         dynum = 5
         call cgp_coord_read_data_f(fid,basenum,zonenum,dynum, &
              minrange,maxrange,dxyz(:,:,:,2),ierr)
         if(ierr .ne. CG_OK) then
            write(*,*) 'moving and pitching Y cgp_coord_read_data_f error'
            call cg_error_print_f()
            call abortmpi()
         end if

         dznum = 6
         call cgp_coord_read_data_f(fid,basenum,zonenum,dznum, &
              minrange,maxrange,dxyz(:,:,:,3),ierr)
         if(ierr .ne. CG_OK) then
            write(*,*) 'moving and pitching Z cgp_coord_read_data_f error'
            call cg_error_print_f()
            call abortmpi()
         end if
      end if      
#endif


      return
      end


!-----------------------------------------------------------------------
      subroutine parallelreadrestart(q,q2,mut,nl)
!     This routine reads the restart file in using MPI-I/O.  It is 
!     bespoke for the restart file and uses functionality that is 
!     also used in the serial I/O code.  This means if you change 
!     any I/O functionality here (which is not just the MPI-I/O code)
!     then it should also be changed in the same code in input.f
!     (i.e. subroutine readrest)
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

#ifdef MPI
      include 'mpif.h'
      integer(kind=cosa_int) :: offset_kind
      parameter(offset_kind=MPI_OFFSET_KIND)
#else
      integer(kind=cosa_int) :: offset_kind
      parameter(offset_kind=8)
#endif

      integer(kind=cosa_int) :: imax,jmax,kmax,nl,iblock,iq,imut
      real(kind=cosa_real) :: q(*),q2(*),mut(*)
      integer(kind=offset_kind) :: disp
      integer(kind=cosa_int) :: fid,integersize,doublesize,ierr,iomode

#ifdef MPI

!     mpi  We now open the file for writing the restart data
      call getreadonlyiomode(iomode)
      call openfile(fid,'restart',iomode)

      call readrestartinitialoffset(disp,nl,integersize,doublesize)

      do iblock = 1,mynblocks
        iq   = 1 + off_p3  (iblock,nl) * npde * dim5
        if (kom.or.kom_bsl.or.kom_sst) then
          imut   = 1 + off_p3 (iblock,nl) *        dim5
        else
          imut = 1
        end if
        imax = i_imax  (iblock,nl)
        jmax = j_jmax  (iblock,nl)
        kmax = k_kmax  (iblock,nl)
        call parallelreadrestartblock(q(iq),q2(iq),mut(imut),imax,jmax, &
             kmax,npde,nharms,fid,disp,integersize,doublesize)

      end do

      call closempiiofile(fid)
      
#endif

      return
      end

!-----------------------------------------------------------------------
      subroutine parallelreadrestartblock(q,q2,mut,imax,jmax,kmax,npde, &
           nharms,fid,disp,integersize,doublesize)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

#ifdef MPI
      include 'mpif.h'
      integer(kind=cosa_int) :: offset_kind
      parameter(offset_kind=MPI_OFFSET_KIND)
#else
      integer(kind=cosa_int) :: offset_kind
      parameter(offset_kind=8)
#endif

      integer(kind=cosa_int) :: imax,jmax,kmax,npde,nharms
      integer(kind=cosa_int) :: i,j,k,fid,doublesize,integersize,ierror
      integer(kind=cosa_int) :: blocksize,blocksize2,blocksize3
      integer(kind=offset_kind) :: disp
      real(kind=cosa_real) :: &
           q (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           q2(-1:imax+1,-1:jmax+1,-1:kmax+1,npde), &
           mut(-1:imax+1,-1:jmax+1,-1:kmax+1 ,0:2*nharms)

#ifdef MPI
      call setuprestartiovalues1(imax,jmax,kmax,npde,nharms,blocksize)
      call setuprestartiovalues2(imax,jmax,kmax,npde,blocksize2)
      call setuprestartiovalues3(imax,jmax,kmax,nharms,blocksize3)

      disp = disp + integersize
      call setupfile(fid,disp,MPI_DOUBLE_PRECISION)
      call MPI_FILE_READ(fid, q(-1,-1,-1,1,0), blocksize, &
           MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierror)
      disp = disp + blocksize*doublesize + integersize


      if (dual_prop_start) then

         disp = disp + integersize
         call setupfile(fid,disp,MPI_DOUBLE_PRECISION)
         call MPI_FILE_READ(fid, q2(-1,-1,-1,1), blocksize2, &
              MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierror)
         disp = disp + blocksize2*doublesize + integersize
         
      end if

      if (kom.or.kom_bsl.or.kom_sst) then

         disp = disp + integersize
         call setupfile(fid,disp,MPI_DOUBLE_PRECISION)
         call MPI_FILE_READ(fid, mut(-1,-1,-1,0), blocksize3, &
              MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierror)
         disp = disp + blocksize3*doublesize + integersize

      end if


#endif

      return
      end

!-----------------------------------------------------------------------
      subroutine parallelwriterestart(q,q1,mut,nl)
!     This routine write the restart file out using MPI-I/O.  It is 
!     bespoke for the restart file and uses functionality that is 
!     also used in the serial I/O code.  This means if you change 
!     any I/O functionality here (which is not just the MPI-I/O code)
!     then it should also be changed in the same code in output.f
!     (i.e. subroutine writerest)
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

#ifdef MPI
      include 'mpif.h'
#endif

      integer(kind=cosa_int) :: iq,imut,imax,jmax,kmax,nl
      integer*8 disp
      integer(kind=cosa_int) :: fid,integersize,doublesize,ierr,iblock,iomode
      real(kind=cosa_real) :: q(*),q1(*),mut(*)

#ifdef MPI
!     mpi  Delete any existing restart file
      call deletefile('restart')

!     mpi  We now open the file for writing the restart data
      call getcreatewriteiomode(iomode)
      call openfile(fid,'restart',iomode)

      call writerestartinitialoffset(disp,nl,integersize,doublesize)

      do iblock = 1,mynblocks
         imax = i_imax  (iblock,nl)
         jmax = j_jmax  (iblock,nl)
         kmax = k_kmax  (iblock,nl)
         iq   = 1 + off_p3  (iblock,nl) * npde * dim5
         if (kom.or.kom_bsl.or.kom_sst) then
           imut   = 1 + off_p3 (iblock,nl) *        dim5
         else
           imut = 1
         end if
         call writerestartblock(q(iq),q1(iq),mut(imut),imax,jmax,kmax, &
              npde,nharms,fid,iblock,disp,integersize,doublesize)
      end do

      call closempiiofile(fid)
      
#endif

      return
      end

!-----------------------------------------------------------------------
      subroutine writerestartblock (q,q1,mut,imax,jmax,kmax,npde, &
                 nharms,fid,iblock,disp,integersize,doublesize)
!     This routine writes parts of the restart file in using MPI-I/O.  
!     It is bespoke for the restart file and uses functionality that is 
!     also used in the serial I/O code.  This means if you change 
!     any I/O functionality here (which is not just the MPI-I/O code)
!     then it should also be changed in the same code in output.f
!     (i.e. subroutine write_brest).
!     We also assume that variables dualt and rgkuns do not change
!     between calls to this routine and on other processes as we have 
!     to calculate the amount of data written to the file by other 
!     processes to find the position that this call is writing data 
!     to.  There should be no reasons that these variables will change 
!     but it is worth noting that this is a constraint on this code.
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

#ifdef MPI 
      include 'mpif.h'
#endif

      integer(kind=cosa_int) :: imax,jmax,kmax,npde,nharms,ipde,iblock
      integer*8 disp
      integer(kind=cosa_int) :: i,j,k,n,np,fid,key1,doublesize,ierr,integersize
      integer(kind=cosa_int) :: linelength, im, jm, km, imax1, jmax1, kmax1
      integer(kind=cosa_int) :: blockamount, amount1
      real(kind=cosa_real) :: &
           q (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           q1(-1:imax+1,-1:jmax+1,-1:kmax+1,npde), &
           mut(-1:imax+1,-1:jmax+1,-1:kmax+1 ,0:2*nharms)

#ifdef MPI

      im = imax-1
      jm = jmax-1
      km = kmax-1
      imax1 = imax+1
      jmax1 = jmax+1
      kmax1 = kmax+1

      call setuprestartiovalues1(imax,jmax,kmax,npde,nharms,amount1)

      linelength = amount1*doublesize

      call setupfile(fid,disp,MPI_INTEGER)
      call mpi_file_write(fid, linelength, 1, &
           MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
      disp = disp + integersize
            
      call setupfile(fid,disp,MPI_DOUBLE_PRECISION)
      call mpi_file_write(fid, q(-1,-1,-1,1,0), linelength/doublesize, &
           MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)
      disp = disp + linelength
           
      call setupfile(fid,disp,MPI_INTEGER)
      call mpi_file_write(fid, linelength, 1, &
           MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
      disp = disp + integersize            

      if (dualt) then

         call setuprestartiovalues2(imax,jmax,kmax,npde,amount1)

         linelength = amount1*doublesize
            
         call setupfile(fid,disp,MPI_INTEGER)
         call mpi_file_write(fid, linelength, 1, &
              MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
         disp = disp + integersize
            
         call setupfile(fid,disp,MPI_DOUBLE_PRECISION)
         call mpi_file_write(fid, q1(-1,-1,-1,1), linelength/doublesize, &
              MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)
         disp = disp + linelength
            
         call setupfile(fid,disp,MPI_INTEGER)
         call mpi_file_write(fid, linelength, 1, &
              MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
         disp = disp + integersize            

      end if

      if (kom.or.kom_bsl.or.kom_sst) then

         call setuprestartiovalues3(imax,jmax,kmax,nharms,amount1)

         linelength = amount1*doublesize

         call setupfile(fid,disp,MPI_INTEGER)
         call mpi_file_write(fid, linelength, 1, &
              MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
         disp = disp + integersize

         call setupfile(fid,disp,MPI_DOUBLE_PRECISION)
         call mpi_file_write(fid, mut(-1,-1,-1,0),linelength/doublesize, &
              MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)
         disp = disp + linelength

         call setupfile(fid,disp,MPI_INTEGER)
         call mpi_file_write(fid, linelength, 1, &
              MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
         disp = disp + integersize

      end if


#endif

      return
      end

!-----------------------------------------------------------------------
      subroutine writeparallelflowtecheader(fid,flowtec,n,nl,harbal, &
           dualt,rgkuns,kom,kom_bsl,kom_sst,unsteady,simtime,zeit, &
           mharms)
!     This opens the flowtec file and writes the header text into it 
!     for the parallel flowtec write.
!-----------------------------------------------------------------------
      
      use cosa_variables
      use cosa_precision

      implicit none

#if MPI 
      include 'mpif.h'
#endif
      integer*8 disp
      character*72 flowtec
      integer(kind=cosa_int) :: fid,charactersize,integersize,doublesize,iomode
      integer(kind=cosa_int) :: harbalset,komset,unsteadyset,dualtset,rgkunsset
      integer(kind=cosa_int) :: ierr,n,nl,i,mharms
      integer(kind=cosa_int), allocatable :: blockindexes(:,:)
      character*100 line
      logical :: harbal,dualt,rgkuns,kom,kom_bsl,kom_sst,unsteady
      real(kind=cosa_real) :: simtime,zeit(0:2*mharms)

      logical :: amcontrol

#if MPI
      call mpi_type_size(MPI_CHARACTER, charactersize, ierr)
      call mpi_type_size(MPI_INTEGER, integersize, ierr)
      call mpi_type_size(MPI_DOUBLE_PRECISION, doublesize, ierr)

      call amcontroller(amcontrol)

!------open tecplot file(s)

      call deletefile(flowtec)
      call getcreatewriteiomode(iomode)
      call openfile(fid,flowtec,iomode)
      
      if(amcontrol) then

         disp = 0
         
!     mpi  For the new MPI I/O Flowtec write we don't add any of the 
!     mpi  text into the file header. This is added in the 
!     mpi  flowtec_bin2asc utility.
!     mpi  The following lines write whether harbal is set, whether 
!     mpi  unsteady is set and whether kom, kom_bsl, or kom_sst is set. 
!     mpi  This is used for the flowtec_bin2asc utility to covert the
!     mpi  binary MPI-I/O format flowtec files into ascii files.

         if (harbal) then
            harbalset = 1
         else
            harbalset = 0
         end if

         call setupfile(fid,disp,MPI_INTEGER)
         call mpi_file_write(fid,harbalset,1, &
              MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
         disp = disp + integersize

         if (dualt) then
            dualtset = 1
         else
           dualtset = 0
         end if

         call setupfile(fid,disp,MPI_INTEGER)
         call mpi_file_write(fid,dualtset,1, &
              MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
         disp = disp + integersize

         if (rgkuns) then
            rgkunsset = 1
         else
           rgkunsset = 0
         end if

         call setupfile(fid,disp,MPI_INTEGER)
         call mpi_file_write(fid,rgkunsset,1, &
              MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
         disp = disp + integersize


         if (kom.or.kom_bsl.or.kom_sst) then
            komset = 1
         else
            komset = 0
         end if

         call setupfile(fid,disp,MPI_INTEGER)
         call mpi_file_write(fid,komset,1, &
              MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
         disp = disp + integersize

         if (unsteady) then
            unsteadyset = 1
         else
           unsteadyset = 0
         end if

         call setupfile(fid,disp,MPI_INTEGER)
         call mpi_file_write(fid,unsteadyset,1, &
              MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
         disp = disp + integersize

         call setupfile(fid,disp,MPI_INTEGER)
         call mpi_file_write(fid,nblocks,1, &
              MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
         disp = disp + integersize

!     If unsteady write the simtime at the start of the file
         if(unsteady) then
            if(dualt) then
               call setupfile(fid,disp,MPI_DOUBLE_PRECISION)
               call mpi_file_write(fid,simtime,1, &
                    MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
               disp = disp + doublesize
            else if(harbal) then
               call setupfile(fid,disp,MPI_DOUBLE_PRECISION)
               call mpi_file_write(fid,zeit(n),1, &
                    MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
               disp = disp + doublesize
            end if
         end if

         call setupfile(fid,disp,MPI_INTEGER)
         call mpi_file_write(fid,n,1, &
              MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
         disp = disp + integersize

         allocate(blockindexes(3,nblocks))

         do i = 1,nblocks
            blockindexes(1,i) = g_i_imax(i,nl) + 1
            blockindexes(2,i) = g_j_jmax(i,nl) + 1
            blockindexes(3,i) = g_k_kmax(i,nl) + 1
         end do

         call setupfile(fid,disp,MPI_INTEGER)
         call mpi_file_write(fid,blockindexes,3*nblocks, &
              MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
         disp = disp + integersize

         deallocate(blockindexes)

      end if

#endif
      
      return
      end

!-----------------------------------------------------------------------
      subroutine writeparallelpltheader(titlename,varlist,flowtec,pwd,n, &
                 nl,fileformat,filetype,isdebug,isdouble)

!     This opens a plt or szplt file in parallel
!
!-----------------------------------------------------------------------
        
      use cosa_variables
      use cosa_precision

      implicit none

#if MPI 
      include 'mpif.h'
#endif

      character*72 titlename,flowtec,varlist
      character*500 pwd
      integer(kind=cosa_int) ierr,n,nl
      integer(kind=cosa_int) tecini142, tecmpiinit142
      integer(kind=cosa_int) fileformat, filetype, isdebug, isdouble


#if MPI

      if(tecini142(trim(titlename),trim(varlist),trim(flowtec),trim(pwd), &
           fileformat,filetype,isdebug,isdouble) .ne. 0) then
         write(*,*) 'error initialising tecini142'
         call abortmpi()
      end if

      if(tecmpiinit142(MPI_COMM_WORLD, 0) .ne. 0) then
         write(*,*) 'error initialising tecmpiinit142'
         call abortmpi()
      end if
      
#endif

      return
      end

!-----------------------------------------------------------------------
      subroutine writeparallelsurftecheader(fid,n,nl)
!     This opens the surftec file and writes the header information  
!     for the parallel surftec write.
!-----------------------------------------------------------------------
      
      use parallelutils, only: surfblockstarts, surfblockends, &
           surfblockindex, surfbcindex, nsurfblocks
      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

#if MPI 
      include 'mpif.h'
#endif

      integer*8 disp
      integer(kind=cosa_int) :: fid,charactersize,integersize,doublesize,iomode
      integer(kind=cosa_int) :: harbalset,dualtset,rgkunsset
      integer(kind=cosa_int) :: ierr,n,nl,iblk
      integer(kind=cosa_int) :: surfacenumber
      integer(kind=cosa_int) :: iblock, imax, jmax, kmax


      logical :: amcontrol

#if MPI
      call mpi_type_size(MPI_CHARACTER, charactersize, ierr)
      call mpi_type_size(MPI_INTEGER, integersize, ierr)
      call mpi_type_size(MPI_DOUBLE_PRECISION, doublesize, ierr)

      call amcontroller(amcontrol)

!------open surface tecplot file(s)

      call deletefile(surftec)
      call getcreatewriteiomode(iomode)
      call openfile(fid,surftec,iomode)
      
      if(amcontrol) then

         disp = 0
         
!     mpi  For the MPI I/O surftec write we don't add any of the 
!     mpi  text into the file header. This is added in the 
!     mpi  surftec_bin2asc utility.
!     mpi  The following lines write whether harbal is set, etc...
!     mpi  This is used for the surftec_bin2asc utility to covert the
!     mpi  binary MPI-I/O format surftec files into ascii files.

         if (harbal) then
            harbalset = 1
         else
            harbalset = 0
         end if

         call setupfile(fid,disp,MPI_INTEGER)
         call mpi_file_write(fid,harbalset,1, &
              MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
         disp = disp + integersize

         if (dualt) then
            dualtset = 1
         else
           dualtset = 0
         end if

         call setupfile(fid,disp,MPI_INTEGER)
         call mpi_file_write(fid,dualtset,1, &
              MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
         disp = disp + integersize

         if (rgkuns) then
            rgkunsset = 1
         else
           rgkunsset = 0
         end if

         call setupfile(fid,disp,MPI_INTEGER)
         call mpi_file_write(fid,rgkunsset,1, &
              MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
         disp = disp + integersize

         call setupfile(fid,disp,MPI_INTEGER)
         call mpi_file_write(fid,nblocks,1, &
              MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
         disp = disp + integersize

         call setupfile(fid,disp,MPI_INTEGER)
         call mpi_file_write(fid,n,1, &
              MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
         disp = disp + integersize

      end if
            
!     Only do this once for all harmonics. We check whether it's been done before by
!     checking if one of the arrays has already been allocated.
      if(.not.allocated(surfblockstarts)) then
       
         nsurfblocks = 0

         do iblock = 1,nblocks
            call calculatenumsurfblocks(g_bcdata(:,:,iblock),g_nbcs(iblock), &
                 nsurfblocks)
         end do
         
         allocate(surfblockstarts(3,nsurfblocks))
         allocate(surfblockends(3,nsurfblocks))
         allocate(surfblockindex(nsurfblocks))
         allocate(surfbcindex(nsurfblocks))
         
         surfacenumber = 1
         do iblock = 1,nblocks
            imax   = g_i_imax     (iblock,nl)
            jmax   = g_j_jmax     (iblock,nl)
            kmax   = g_k_kmax     (iblock,nl)            
            call calculatesurfblock(iblock,g_bcdata(:,:,iblock),g_nbcs(iblock), &
                 imax,jmax,kmax,nsurfblocks, &
                 surfblockstarts,surfblockends,surfblockindex,surfbcindex, &
                 surfacenumber)
         end do
      end if     

      if(amcontrol) then

         call setupfile(fid,disp,MPI_INTEGER)
         call mpi_file_write(fid,nsurfblocks,1, &
              MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
         disp = disp + integersize
         
         call setupfile(fid,disp,MPI_INTEGER)
         call mpi_file_write(fid,surfblockstarts,3*nsurfblocks, &
              MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
         disp = disp + integersize*3*nsurfblocks        

         call setupfile(fid,disp,MPI_INTEGER)
         call mpi_file_write(fid,surfblockends,3*nsurfblocks, &
              MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
         disp = disp + integersize*3*nsurfblocks

      end if      
      
#endif
      
      return
      end

!-----------------------------------------------------------------------
      subroutine writeparallelsurfpltheader(titlename,varlist,flowtec,pwd, &
                 n,nl,fileformat,filetype,isdebug,isdouble)

!     This opens the surftec file and writes the header information  
!     for the parallel surftec write.
!-----------------------------------------------------------------------
      
      use parallelutils, only: surfblockstarts, surfblockends, &
           surfblockindex, surfbcindex, nsurfblocks
      use cosa_variables
   !   use common_variables
      use cosa_precision

      implicit none

#if MPI 
      include 'mpif.h'
#endif

      character*72 titlename,flowtec,varlist
      character*500 pwd
      integer*8 disp
      integer(kind=cosa_int) :: fid,charactersize,integersize,doublesize,iomode
      integer(kind=cosa_int) :: harbalset,dualtset,rgkunsset
      integer(kind=cosa_int) :: ierr,n,nl,iblk
      integer(kind=cosa_int) :: tecini142, tecmpiinit142
      integer(kind=cosa_int) :: surfacenumber
      integer(kind=cosa_int) :: iblock, imax, jmax, kmax
      integer(kind=cosa_int) fileformat, filetype, isdebug, isdouble


      logical :: amcontrol

#if MPI            
!     Only do this once for all harmonics. We check whether it's been done before by
!     checking if one of the arrays has already been allocated.
      if(.not.allocated(surfblockstarts)) then
       
         nsurfblocks = 0

         do iblock = 1,nblocks
            call calculatenumsurfblocks(g_bcdata(:,:,iblock),g_nbcs(iblock), &
                 nsurfblocks)
         end do
         
         allocate(surfblockstarts(3,nsurfblocks))
         allocate(surfblockends(3,nsurfblocks))
         allocate(surfblockindex(nsurfblocks))
         allocate(surfbcindex(nsurfblocks))
         
         surfacenumber = 1
         do iblock = 1,nblocks
            imax   = g_i_imax     (iblock,nl)
            jmax   = g_j_jmax     (iblock,nl)
            kmax   = g_k_kmax     (iblock,nl)            
            call calculatesurfblock(iblock,g_bcdata(:,:,iblock),g_nbcs(iblock), &
                 imax,jmax,kmax,nsurfblocks, &
                 surfblockstarts,surfblockends,surfblockindex,surfbcindex, &
                 surfacenumber)
         end do
      end if

      if(tecini142(trim(titlename),trim(varlist),trim(flowtec),trim(pwd), &
           fileformat,filetype,isdebug,isdouble) .ne. 0) then
         write(*,*) 'error initialising tecini142'
         call abortmpi()
      end if

      if(tecmpiinit142(MPI_COMM_WORLD, 0) .ne. 0) then
         write(*,*) 'error initialising tecmpiinit142'
         call abortmpi()
      end if

      
#endif
      
      return
      end


!-----------------------------------------------------------------------
      subroutine calculatenumsurfblocks(bctopo,nbcs,nsurfblocks)
!-----------------------------------------------------------------------

      use cosa_precision

      implicit none


      integer(kind=cosa_int) :: isurface, ibc1
      integer(kind=cosa_int) :: nsurfblocks
      integer(kind=cosa_int) :: bctopo(10,nbcs),nbcs,nsurface

      logical :: amcontrol

      call amcontroller(amcontrol)

      do ibc1=1,nbcs
         if (any(bctopo(1,ibc1).eq. &
              [1500,1501,1502,1503,1400,1401,1402,1403])) then
            nsurfblocks = nsurfblocks + 1
         end if
      end do            

      return 
      end

!-----------------------------------------------------------------------
      subroutine calculatesurfblock(iblock,bctopo,nbcs,imax,jmax,kmax, &
           nsurfblocks,surfacestartindexes,surfaceendindexes, &
           surfaceblockindex,surfacebcindex,surfacenumber)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) :: bctopo(10,nbcs)
      integer(kind=cosa_int) :: ierr,nbcs,imax,jmax,kmax,iblock
      integer(kind=cosa_int) :: nsurfblocks
      integer(kind=cosa_int) :: surfacenumber
      integer(kind=cosa_int) :: ibc1, block
      integer(kind=cosa_int) :: surfacestartindexes(3,nsurfblocks)
      integer(kind=cosa_int) :: surfaceendindexes(3,nsurfblocks)
      integer(kind=cosa_int) :: surfaceblockindex(nsurfblocks)
      integer(kind=cosa_int) :: surfacebcindex(nsurfblocks)
      integer(kind=cosa_int) :: i1,i2,i3,ijkmax(3),idir,istrt(3),iend(3),bctyp, &
        inrout,ibcpt,ibcpt2,ibcn,ibcn2,ibcm,ibc,jbc,kbc,ibc2,jbc2, &
        kbc2,ic1,ic2,ic3,ioff1,ioff2,joff1,joff2,koff1,koff2,sgnm, &
        imax1,jmax1,kmax1
      logical :: amcontrol

      call amcontroller(amcontrol)

      ijkmax(1) = imax
      ijkmax(2) = jmax
      ijkmax(3) = kmax

      imax1 = imax+1
      jmax1 = jmax+1
      kmax1 = kmax+1

      if(surfacenumber .gt. nsurfblocks) then

         !AJ We've found all our surface blocks, so don't do anything
         !AJ But check there isn't a mistake
         do ibc1=1,nbcs
            
            if (any(bctopo(1,ibc1).eq. &
                 [1500,1501,1502,1503,1400,1401,1402,1403])) then            

               write(*,*) 'Found a surface block in calculatesurfblock '
               write(*,*) 'but we have already found all expected blocks.'
               write(*,*) 'This means something has gone wrong.'
               call abort()
               
            end if

         end do

      else

         do ibc1=1,nbcs
            
            if (any(bctopo(1,ibc1).eq. &
                 [1500,1501,1502,1503,1400,1401,1402,1403])) then            
               
               
!     store boundary topology in mnemonic names
               bctyp    = bctopo(1,ibc1)
               idir     = bctopo(2,ibc1)
               inrout   = bctopo(3,ibc1)
               istrt(1) = bctopo(4,ibc1)
               iend(1)  = bctopo(5,ibc1)
               istrt(2) = bctopo(6,ibc1)
               iend(2)  = bctopo(7,ibc1)
               istrt(3) = bctopo(8,ibc1)
               iend(3)  = bctopo(9,ibc1)               
               
!     set needed variables depending on whether the boundary is
!     the inner boundary (inrout = 1) or
!     the outer boundary (inrout > 1)
!     ibcpt : boundary condition location (first aux. cell)
!     ibcpt2: boundary condition location outside the block from
!     ibcpt (second aux. cell)
!     ibcn  : point to the inside of the block from ibcpt
!     ibcm  : location of the metrics
               
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
               
               surfacestartindexes(1,surfacenumber) = ibcpt &
                    *krd(ic1,1) + istrt(ic2)*krd(ic2,1) + &
                    istrt(ic3)*krd(ic3,1)
               surfacestartindexes(2,surfacenumber) = ibcpt &
                    *krd(ic1,2) + istrt(ic2)*krd(ic2,2) + &
                    istrt(ic3)*krd(ic3,2)
               surfacestartindexes(3,surfacenumber) = ibcpt &
                    *krd(ic1,3) + istrt(ic2)*krd(ic2,3) + &
                    istrt(ic3)*krd(ic3,3)
               
               surfaceendindexes(1,surfacenumber) = ibcpt &
                    *krd(ic1,1) + iend(ic2)*krd(ic2,1) + &
                    iend(ic3)*krd(ic3,1)
               surfaceendindexes(2,surfacenumber) = ibcpt &
                    *krd(ic1,2) + iend(ic2)*krd(ic2,2) + &
                    iend(ic3)*krd(ic3,2)
               surfaceendindexes(3,surfacenumber) = ibcpt &
                    *krd(ic1,3) + iend(ic2)*krd(ic2,3) + &
                    iend(ic3)*krd(ic3,3)
               
               surfaceblockindex(surfacenumber) = iblock
               surfacebcindex(surfacenumber) = ibc1
               
               surfacenumber = surfacenumber + 1
            end if
         end do            
         
      end if
      
      return
      end

!-----------------------------------------------------------------------
      subroutine writeparallelsurfcgnsheader(fid,basenum,n,nl)
!     This opens the surftec file and writes the header information  
!     for the parallel surftec write using cgns.
!-----------------------------------------------------------------------
      
      use parallelutils, only: surfblockstarts, surfblockends, &
           surfblockindex, surfbcindex, nsurfblocks

      use cgns
      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) :: fid(0:2*nharms), n, nl
      integer(kind=cosa_int) :: basenum, blocknum
      integer(kind=cosa_int) :: surfacenumber
      integer(kind=cosa_int) :: iblock, ierr, imax, jmax, kmax
      integer(kind=cosa_int) :: xnum, ynum, znum, flownum, solnum
      integer(cgsize_t) :: sizes(3,3)
      character*20 zonename

      logical :: amcontrol

      call amcontroller(amcontrol)

!     Only do this once for all harmonics. We check whether it's been done before by
!     checking if one of the arrays has already been allocated.
      if(.not.allocated(surfblockstarts)) then
       
         nsurfblocks = 0

         do iblock = 1,nblocks
            call calculatenumsurfblocks(g_bcdata(:,:,iblock),g_nbcs(iblock), &
                 nsurfblocks)
         end do
         
         allocate(surfblockstarts(3,nsurfblocks))
         allocate(surfblockends(3,nsurfblocks))
         allocate(surfblockindex(nsurfblocks))
         allocate(surfbcindex(nsurfblocks))
         
         surfacenumber = 1
         do iblock = 1,nblocks
            imax   = g_i_imax     (iblock,nl)
            jmax   = g_j_jmax     (iblock,nl)
            kmax   = g_k_kmax     (iblock,nl)            
            call calculatesurfblock(iblock,g_bcdata(:,:,iblock),g_nbcs(iblock), &
                 imax,jmax,kmax,nsurfblocks, &
                 surfblockstarts,surfblockends,surfblockindex,surfbcindex, &
                 surfacenumber)
         end do
      end if

      call cgp_pio_mode_f(CGP_INDEPENDENT, ierr)
      if(ierr .ne. CG_OK) then
         write(*,*) 'cgp_pio_mode_f error writeparallelsurfcgnsheader'
         call cg_error_print_f()
         call abortmpi()
      end if

      call cgp_open_f(surftec,CG_MODE_WRITE,fid(n),ierr)
      if(ierr .ne. CG_OK) then
         write(*,*) 'cg_open_f error'
         call cg_error_print_f()
      end if
      call cg_base_write_f(fid(n),'gridbase',3,3,basenum,ierr)
      if(ierr .ne. CG_OK) then
         write(*,*) 'cg_base_write_f error'
         call cg_error_print_f()
      end if

!     AJ Write all the surface blocks as zones.
      do iblock = 1,nsurfblocks
         
         imax = (surfblockends(1,iblock)-surfblockstarts(1,iblock))+1
         jmax = (surfblockends(2,iblock)-surfblockstarts(2,iblock))+1
         kmax = (surfblockends(3,iblock)-surfblockstarts(3,iblock))+1
         
         sizes(1,1) = imax
         sizes(2,1) = jmax
         sizes(3,1) = kmax
         
         sizes(1,2) = imax-1
         sizes(2,2) = jmax-1
         sizes(3,2) = kmax-1
         
         sizes(1,3) = 0
         sizes(2,3) = 0
         sizes(3,3) = 0

         write(zonename, "(A12,I6)") "surfaceblock",iblock
         call cg_zone_write_f(fid(n),basenum,zonename,sizes, &
              Structured,blocknum,ierr)
         if(ierr .ne. CG_OK) then
            write(*,*) &
                 'cg_zone_write_f error in writeparallelsurfcgnsheader'
            call cg_error_print_f()
         end if
         if(blocknum .ne. iblock) then
            write(*,*) 'CGNS block not the same number as our loop',blocknum,iblock
            call abortmpi()
         end if
         call cg_sol_write_f(fid(n),basenum,iblock,'SurfaceSolutions', &
              Vertex,solnum,ierr)
!     &        CellCenter,solnum,ierr)
         if(ierr .ne. CG_OK) then
            write(*,*) &
                 'cg_sol_write_f error in writeparallelcgnsheader'
            call cg_error_print_f()
         end if
         call cgp_coord_write_f(fid(n),basenum,iblock,RealDouble, &
              'CoordinateX',xnum,ierr)
         if(ierr .ne. CG_OK) then
            write(*,*) 'X cgp_coord_write_f error'
            call cg_error_print_f()
         end if
         call cgp_coord_write_f(fid(n),basenum,iblock,RealDouble, &
              'CoordinateY',ynum,ierr)
         if(ierr .ne. CG_OK) then
            write(*,*) 'Y cgp_coord_write_f error'
            call cg_error_print_f()
         end if
         call cgp_coord_write_f(fid(n),basenum,iblock,RealDouble, &
              'CoordinateZ',znum,ierr)
         if(ierr .ne. CG_OK) then
            write(*,*) 'Z cgp_coord_write_f error'
            call cg_error_print_f()
         end if
         call cgp_field_write_f(fid(n),basenum,iblock,solnum, &
              RealDouble,'StaticPressureCoefficient',flownum,ierr)
         if(ierr .ne. CG_OK) then
            write(*,*) 'cgp_field_write_f StaticPressureCoefficient error'
            call cg_error_print_f()
         end if
         call cgp_field_write_f(fid(n),basenum,iblock,solnum, &
              RealDouble,'SkinFrictionCoefficient',flownum,ierr)
         if(ierr .ne. CG_OK) then
            write(*,*) 'cgp_field_write_f SkinFrictionCoefficient error'
            call cg_error_print_f()
         end if
         call cgp_field_write_f(fid(n),basenum,iblock,solnum, &
              RealDouble,'YPlus',flownum,ierr)
         if(ierr .ne. CG_OK) then
            write(*,*) 'cgp_field_write_f YPlus error'
            call cg_error_print_f()
         end if
      end do
      
      
      call cgp_close_f(fid(n),ierr)
      if(ierr .ne. CG_OK) then
         write(*,*) 'cg_close_f error'
         call cg_error_print_f()
      end if
                        
      return
      end

!-----------------------------------------------------------------------
      subroutine write_parallel_wr_tec_b(fid,var1,var2,imax,jmax,kmax, &
                 npde,nharms,iblock,nl)
!     This routine writes parts of the flowtec files using MPI-I/O.
!     Mainly this routine wraps the bespoke code which writes the
!     data (in parallel_wr_tec_b) and also deals with setting up
!     the file and moving the file pointer or displacement to ensure the
!     block is written in the correct place.
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

#if MPI 
      include 'mpif.h'
#endif

      integer(kind=cosa_int) :: imax,jmax,kmax,npde,nharms,iblock,nl
      integer(kind=cosa_int) :: fid(0:2*nharms)
      real(kind=cosa_real) :: &
           var1( 0:imax  , 0:jmax  , 0:kmax  ,npde,0:2*nharms), &
           var2( 0:imax  , 0:jmax  , 0:kmax  ,npde,0:2*nharms)
      integer(kind=cosa_int) :: integersize,doublesize,charactersize
      integer*8 disp

      call newflowtecinitialoffset(disp,nl,iblock,integersize, &
           doublesize,charactersize)

      call parallel_wr_tec_b(fid,var1,var2,imax,jmax,kmax,npde,nharms, &
           integersize,doublesize,charactersize,disp)

      return 
      end

!-----------------------------------------------------------------------
      subroutine parallel_wr_tec_b(fid,var1,var2,imax,jmax,kmax,npde, &
                 nharms,integersize,doublesize,charactersize,disp)
!     This routine writes parts of the flowtec files using MPI-I/O.  
!     It is bespoke for the flowtec files and uses functionality that is 
!     also used in the serial I/O code.  This means if you change 
!     any I/O functionality here (which is not just the MPI-I/O code)
!     then it should also be changed in the same code in output.f
!     (i.e. subroutine wr_tec_b).
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

#if MPI 
      include 'mpif.h'
#endif

      integer(kind=cosa_int) :: imax,jmax,kmax,npde,nharms
      integer(kind=cosa_int) :: i,imax1,j,jmax1,k,kmax1,ipde,n,nh
      integer(kind=cosa_int) :: fid(0:2*nharms)
      real(kind=cosa_real) :: &
           var1( 0:imax  , 0:jmax  , 0:kmax  ,npde,0:2*nharms), &
           var2( 0:imax  , 0:jmax  , 0:kmax  ,npde,0:2*nharms), &
           tempdata( 2 * (imax+1) * (jmax+1) * (kmax+1) * npde)
      integer*8 disp,initialdisp
      integer(kind=cosa_int) :: ierr,integersize,doublesize,charactersize,datasize
      integer*8 tempindex
      character*150 line1
      
#ifdef MPI

      imax1 = imax+1
      jmax1 = jmax+1
      kmax1 = kmax+1

      datasize = imax1*jmax1*kmax1*npde*2

      initialdisp = disp

      do n = 0,2*nharms
        nh = n*hbmove

!     This line is included because the n loop iterates through files 
!     so we need to ensure that each time round the loop we start and 
!     the correct place in the associated file.
        disp = initialdisp

        tempindex = 1

        do ipde=1,npde
          do k=0,kmax
            do j=0,jmax
              do i=0,imax
                 tempdata(tempindex) = var1(i,j,k,ipde,n)
                 tempindex = tempindex + 1
              end do
            end do
          end do
        end do
        do ipde=1,npde
          do k=0,kmax
            do j=0,jmax
              do i=0,imax
                 tempdata(tempindex) = var2(i,j,k,ipde,n)
                 tempindex = tempindex + 1
              end do
            end do
          end do
        end do

        call setupfile(fid(n),disp,MPI_DOUBLE_PRECISION)
        call mpi_file_write(fid(n),tempdata(1),datasize, &
             MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
        disp = disp + datasize*doublesize

      end do

#endif

      return
      end


!-----------------------------------------------------------------------
      subroutine write_parallel_wr_plt_b(fid,var1,var2,imax,jmax,kmax, &
                 npde,nharms,iblock,nl)
!     This routine writes parts of the flowtec files using tecio.
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

#if MPI 
      include 'mpif.h'
#endif

      integer(kind=cosa_int) :: imax,jmax,kmax,npde,nharms,iblock,nl,n,nh
      integer(kind=cosa_int) :: imax1,jmax1,kmax1
      integer(kind=cosa_int) :: fid(0:2*nharms)
      real(kind=cosa_real) :: &
           var1( 0:imax  , 0:jmax  , 0:kmax  ,npde,0:2*nharms), &
           var2( 0:imax  , 0:jmax  , 0:kmax  ,npde,0:2*nharms)      
      character(len=5) :: blocknumname
      integer(kind=cosa_int) :: nfconns, fnmode, shrconn, isblock, valuelocation, isdouble
      integer(kind=cosa_int) :: tnfnodes, ncbfaces, tnbconns
      integer(kind=cosa_int) :: zonetype, strandid, parentzone
      integer(kind=cosa_int) :: teczne142, tecznemap142, tecdat142, tecfil142
      integer(kind=cosa_int) :: imaxmax, jmaxmax, kmaxmax
      integer(kind=cosa_int) :: mpiid
      integer(kind=cosa_int) :: Null(*)
      POINTER   (NullPtr,Null)

#ifdef MPI

      call getmpiid(mpiid)
      
      imax1 = imax+1
      jmax1 = jmax+1
      kmax1 = kmax+1

! Specify that we are using an Ordered zone type
      zonetype = 0

! Zones don't have parents
      parentzone = 0

! We are not part of a strand
      strandid = 0
      nfconns = 0
      fnmode = 0
      shrconn = 0
      imaxmax = 0
      jmaxmax = 0
      kmaxmax = 0

! This specifies if we are writing in block or point format
! 1 is block format.  Binary tecplot files must be block format, 
! so this must be 1.
      isblock = 1
      valuelocation = 0
      isdouble = 1
      tnfnodes = 0
      ncbfaces = 0
      tnbconns = 0

      do n = 0,2*nharms
        nh = n*hbmove

        if(tecfil142(n+1) .ne. 0) then
           write(*,*) 'error calling tecfil142'
           stop
        end if

        write (blocknumname, "(I5)") iblock

        if(teczne142('block'//trim(blocknumname)//char(0),zonetype, imax1, jmax1, kmax1, &
             imaxmax, jmaxmax, kmaxmax, simtime, strandid, parentzone, isblock, nfconns, & 
             fnmode, tnfnodes, ncbfaces, tnbconns, Null, Null, Null, shrconn) .ne. 0) then
           write(*,*) 'error setting up zone'
           call abortmpi()
        end if
        
        if(tecznemap142(1, mpiid) .ne. 0) then
           write(*,*) 'error setting up parallel mapping for zone'
           call abortmpi()
        end if

        if (kom.or.kom_bsl.or.kom_sst) then

           if(tecdat142(imax1*jmax1*kmax1*7,var1(0,0,0,1,n),isdouble) .ne. 0) then
              write(*,*) 'error writing block data'
              call abortmpi()
           end if
           if(tecdat142(imax1*jmax1*kmax1*7,var2(0,0,0,1,n),isdouble) .ne. 0) then
              write(*,*) 'error writing block data'
              call abortmpi()
           end if

        else

           if(tecdat142(imax1*jmax1*kmax1*5,var1(0,0,0,1,n),isdouble) .ne. 0) then
              write(*,*) 'error writing block data'
              call abortmpi()
           end if
           if(tecdat142(imax1*jmax1*kmax1*5,var2(0,0,0,1,n),isdouble) .ne. 0) then
              write(*,*) 'error writing block data'
              call abortmpi()
           end if
           
        end if

      end do

#endif

      return
      end



! AJ TODO Merge the tec and cgns versions of this to save code duplication
!-----------------------------------------------------------------------
      subroutine write_parallel_wr_tec_surf_b_M1(fid,imax, &
                 jmax,kmax,npde,nharms,iblock,nl,lmet,q,x,y,z,xdot,ydot, &
                  zdot,si,sj,sk,xideri,etaderj,zetaderk,dist,bctopo,nbcs)
!     This routine writes parts of the surftec files using MPI-I/O.
!     Mainly this routine wraps the bespoke code which writes the
!     data (in parallel_wr_surf_tec_b_M1) and also deals with setting up
!     the file and moving the file pointer or displacement to ensure the
!     block is written in the correct place.
!-----------------------------------------------------------------------

      use parallelutils, only: nsurfblocks, surfblockstarts, &
           surfblockends, surfblockindex, surfbcindex
      use common_variables
      use cosa_precision

      implicit none

#if MPI 
      include 'mpif.h'
#endif

      integer(kind=cosa_int) :: i,imax,jmax,kmax,npde,nharms,iblock,nl,lmet,nbcs
      integer(kind=cosa_int) :: fid(0:2*nharms)
      integer(kind=cosa_int) :: bctopo(10,nbcs)
      real(kind=cosa_real) :: &
           q   (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           x   ( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           y   ( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           z   ( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           xdot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           ydot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           zdot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           si      (lmet, 0:imax+1, 0:jmax+1, 0:kmax+1,0:2*nharms*hbmove), &
           sj      (lmet, 0:imax+1, 0:jmax+1, 0:kmax+1,0:2*nharms*hbmove), &
           sk      (lmet, 0:imax+1, 0:jmax+1, 0:kmax+1,0:2*nharms*hbmove), &
           xideri  (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           etaderj (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           zetaderk(3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           dist( 0:imax  , 0:jmax  , 0:kmax)
      integer(kind=cosa_int) :: integersize,doublesize,charactersize
      integer*8 disp
      integer(kind=cosa_int) :: surfacenumber

      disp = 0

      call surftecinitialoffset(disp,nl,iblock,nsurfblocks, &
           surfblockstarts,surfblockends,surfblockindex,surfbcindex, &
           surfacenumber,integersize,doublesize,charactersize)

      call parallel_wr_tec_surf_b_M1(fid,q,x,y,z,xdot,ydot,zdot,si,sj, &
           sk,xideri,etaderj,zetaderk,dist,bctopo,imax,jmax,kmax,npde, &
           nharms,nl,nbcs,lmet,integersize,doublesize, &
           charactersize,disp,surfblockstarts,surfblockends, &
           nsurfblocks,surfacenumber)

      return 
      end

!-----------------------------------------------------------------------
      subroutine parallel_wr_tec_surf_b_M1(fid,q,x,y,z,xdot,ydot,zdot, &
                 si,sj,sk,xideri,etaderj,zetaderk,dist,bctopo,imax,jmax, &
                 kmax,npde,nharms,nl,nbcs,lmet,integersize, &
                 doublesize,charactersize,disp,surfacestartindexes, &
                 surfaceendindexes,nsurfblocks,surfacenumber)
!     This routine writes parts of the flowtec files using MPI-I/O.  
!     It is bespoke for the flowtec files and uses functionality that is 
!     also used in the serial I/O code.  This means if you change 
!     any I/O functionality here (which is not just the MPI-I/O code)
!     then it should also be changed in the same code in output.f
!     (i.e. subroutine wr_tec_b).
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

#if MPI 
      include 'mpif.h'
#endif

      integer(kind=cosa_int) :: imax,jmax,kmax,npde,nharms,nl,lmet,nbcs
      integer(kind=cosa_int) :: i,imax1,j,jmax1,k,kmax1,ipde,n,nh,ibc1
      integer(kind=cosa_int) :: i1,i2,i3,ijkmax(3),idir,istrt(3),iend(3),bctyp, &
        inrout,ibcpt,ibcpt2,ibcn,ibcn2,ibcm,ibc,jbc,kbc,ibc2,jbc2, &
        kbc2,in,jn,kn,in2,jn2,kn2,im,jm,km,ic1,ic2,ic3,sysize,ioff1, &
        ioff2,joff1,joff2,koff1,koff2,ibc_s,jbc_s,kbc_s,ibc_e,jbc_e, &
        kbc_e
      integer(kind=cosa_int) :: fid(0:2*nharms)
      integer(kind=cosa_int) :: bctopo(10,nbcs)
      real(kind=cosa_real) :: sgnm,nx,ny,nz,kx,ky,kz,rhow,tw,muw,uw,u1,u2,uwr,u1r,u2r, &
          vw,v1,v2,vwr,v1r,v2r,ww,w1,w2,wwr,w1r,w2r,ueta,veta,weta, &
           etax,etay,etaz,dudx,dudy,dudz,dvdx,dvdy, &
           dvdz,dwdx,dwdy,dwdz,divv,txx,txy,txz,tyy,tyz,tzz,mach,tauwx, &
           tauwy,tauwz,tauwpx,tauwpy,tauwpz,dn1,tauw,utau
      real(kind=cosa_real) :: &
           q   (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           x   ( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           y   ( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           z   ( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           xdot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           ydot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           zdot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           si      (lmet, 0:imax+1, 0:jmax+1, 0:kmax+1,0:2*nharms*hbmove), &
           sj      (lmet, 0:imax+1, 0:jmax+1, 0:kmax+1,0:2*nharms*hbmove), &
           sk      (lmet, 0:imax+1, 0:jmax+1, 0:kmax+1,0:2*nharms*hbmove), &
           xideri  (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           etaderj (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           zetaderk(3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           dist( 0:imax  , 0:jmax  , 0:kmax)
      real(kind=cosa_real),allocatable :: cp(:,:,:),cf(:,:,:),yp(:,:,:)
      real(kind=cosa_real),allocatable :: tempdata(:)
      integer*8 disp,initialdisp
      integer(kind=cosa_int) :: ierr,integersize,doublesize,charactersize
      integer(kind=cosa_int) :: tempindex,datasize
      integer(kind=cosa_int) :: surfacestartindexes(3,nsurfblocks)
      integer(kind=cosa_int) :: surfaceendindexes(3,nsurfblocks)
      integer(kind=cosa_int) :: nsurfblocks, surfacenumber
      character*200 line1
      integer(kind=cosa_int) :: initialsurfacenumber
      
#ifdef MPI

      ijkmax(1) = imax
      ijkmax(2) = jmax
      ijkmax(3) = kmax

      imax1 = imax+1
      jmax1 = jmax+1
      kmax1 = kmax+1

      initialdisp = disp
      initialsurfacenumber = surfacenumber

      do n = 0,2*nharms
        nh = n*hbmove

!     This line is included because the n loop iterates through files 
!     so we need to ensure that each time round the loop we start and 
!     the correct place in the associated file.
        disp = initialdisp
        surfacenumber = initialsurfacenumber
           
        do ibc1=1,nbcs

           if (any(bctopo(1,ibc1).eq. &
                [1500,1501,1502,1503,1400,1401,1402,1403])) then
              
              ibc_s = surfacestartindexes(1,surfacenumber)
              jbc_s = surfacestartindexes(2,surfacenumber)
              kbc_s = surfacestartindexes(3,surfacenumber)
              
              ibc_e = surfaceendindexes(1,surfacenumber)
              jbc_e = surfaceendindexes(2,surfacenumber)
              kbc_e = surfaceendindexes(3,surfacenumber)
              
              datasize = 6
              if((ibc_e-ibc_s)+1 .ne. 0) datasize = datasize*((ibc_e-ibc_s)+1)
              if((jbc_e-jbc_s)+1 .ne. 0) datasize = datasize*((jbc_e-jbc_s)+1)
              if((kbc_e-kbc_s)+1 .ne. 0) datasize = datasize*((kbc_e-kbc_s)+1)

              allocate(tempdata(datasize))
              allocate(cp(ibc_s:ibc_e,jbc_s:jbc_e,kbc_s:kbc_e))
              if (viscous) then
                 allocate(cf(ibc_s:ibc_e,jbc_s:jbc_e,kbc_s:kbc_e), &
                      yp(ibc_s:ibc_e,jbc_s:jbc_e,kbc_s:kbc_e))
              end if
                            
!     store boundary topology in mnemonic names
              bctyp    = bctopo(1,ibc1)
              idir     = bctopo(2,ibc1)
              inrout   = bctopo(3,ibc1)
              istrt(1) = bctopo(4,ibc1)
              iend(1)  = bctopo(5,ibc1)
              istrt(2) = bctopo(6,ibc1)
              iend(2)  = bctopo(7,ibc1)
              istrt(3) = bctopo(8,ibc1)
              iend(3)  = bctopo(9,ibc1)
              
!             set needed variables depending on whether the boundary is
!             the inner boundary (inrout = 1) or
!             the outer boundary (inrout > 1)
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
                 sgnm   = - 1.d0
              else
                 ibcpt  = ijkmax(idir)
                 ibcpt2 = ijkmax(idir) + 1
                 ibcn   = ijkmax(idir) - 1
                 ibcn2  = ijkmax(idir) - 2
                 ibcm   = ijkmax(idir)
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
 
              do i3 = istrt(ic3),iend(ic3)
                 do i2 = istrt(ic2),iend(ic2)
                    
                    ibc  = ibcpt *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
                    jbc  = ibcpt *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
                    kbc  = ibcpt *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)
                    
                    im   = ibcm  *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
                    jm   = ibcm  *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
                    km   = ibcm  *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)                   


                    if (idir.eq.1) then
                       nx    = si(1,im,jm,km,nh)
                       ny    = si(2,im,jm,km,nh)
                       nz    = si(3,im,jm,km,nh)
                    else if (idir.eq.2) then
                       nx    = sj(1,im,jm,km,nh)
                       ny    = sj(2,im,jm,km,nh)
                       nz    = sj(3,im,jm,km,nh)
                    else if (idir.eq.3) then
                       nx    = sk(1,im,jm,km,nh)
                       ny    = sk(2,im,jm,km,nh)
                       nz    = sk(3,im,jm,km,nh)
                    end if
                    
                    kx   = (-1)**(1+1/inrout) * nx
                    ky   = (-1)**(1+1/inrout) * ny
                    kz   = (-1)**(1+1/inrout) * nz
                    
!     AJ Removed the machfs division here because if we are not viscous 
!     AJ machfs is zero.
                    cp(ibc,jbc,kbc) = 2*(q(ibc,jbc,kbc,5,n)-1/gamma)                                        
                    
                    if (viscous) then

                       cp(ibc,jbc,kbc) = cp(ibc,jbc,kbc)/ &
                            machfs**2

                       
                       ibc2 = ibcpt2*krd(ic1,1) +i2*krd(ic2,1) +i3*krd(ic3,1)
                       jbc2 = ibcpt2*krd(ic1,2) +i2*krd(ic2,2) +i3*krd(ic3,2)
                       kbc2 = ibcpt2*krd(Ic1,3) +i2*krd(ic2,3) +i3*krd(ic3,3)
                       
                       in   = ibcn  *krd(ic1,1) +i2*krd(ic2,1) +i3*krd(ic3,1)
                       jn   = ibcn  *krd(ic1,2) +i2*krd(ic2,2) +i3*krd(ic3,2)
                       kn   = ibcn  *krd(ic1,3) +i2*krd(ic2,3) +i3*krd(ic3,3)
                       
                       in2  = ibcn2 *krd(ic1,1) +i2*krd(ic2,1) +i3*krd(ic3,1)
                       jn2  = ibcn2 *krd(ic1,2) +i2*krd(ic2,2) +i3*krd(ic3,2)
                       kn2  = ibcn2 *krd(ic1,3) +i2*krd(ic2,3) +i3*krd(ic3,3)
                       
                       if (idir.eq.1) then
                          etax  = xideri(1,im,jm,km,nh)
                          etay  = xideri(2,im,jm,km,nh)
                          etaz  = xideri(3,im,jm,km,nh)
                       else if (idir.eq.2) then
                          etax  = etaderj(1,im,jm,km,nh)
                          etay  = etaderj(2,im,jm,km,nh)
                          etaz  = etaderj(3,im,jm,km,nh)
                       else if (idir.eq.3) then
                          etax  = zetaderk(1,im,jm,km,nh)
                          etay  = zetaderk(2,im,jm,km,nh)
                          etaz  = zetaderk(3,im,jm,km,nh)
                       end if
                       
                       rhow = q(ibc,jbc,kbc,1,n)
                       tw   = q(ibc,jbc,kbc,5,n) *gamma / q(ibc,jbc,kbc,1,n)
                       muw  = (stemp+1)/(tw+stemp) * (tw**1.5d0)
                       
!------------------MSC, 11/03/2023: dist not computed for laminar now
                       if (kom.or.kom_bsl.or.kom_sst) then
                          dn1  = dist(in ,jn ,kn )
                       else
                          dn1  = 0
                       end if
                       
                       uw   = q(ibc,jbc,kbc,2,n)
                       u1   = q(ibc+ioff1,jbc+joff1,kbc+koff1,2,n)
                       u2   = q(ibc+ioff2,jbc+joff2,kbc+koff2,2,n)
                       if (moving) then
                          uwr  = xdot(ibc      ,jbc      ,kbc      ,nh)
                          u1r  = xdot(ibc+ioff1,jbc+joff1,kbc+koff1,nh)
                          u2r  = xdot(ibc+ioff2,jbc+joff2,kbc+koff2,nh)
                          uw   = uw - uwr
                          u1   = u1 - u1r
                          u2   = u2 - u2r
                       end if
                       
                       vw   = q(ibc,jbc,kbc,3,n)
                       v1   = q(ibc+ioff1,jbc+joff1,kbc+koff1,3,n)
                       v2   = q(ibc+ioff2,jbc+joff2,kbc+koff2,3,n)
                       if (moving) then
                          vwr  = ydot(ibc      ,jbc      ,kbc      ,nh)
                          v1r  = ydot(ibc+ioff1,jbc+joff1,kbc+koff1,nh)
                          v2r  = ydot(ibc+ioff2,jbc+joff2,kbc+koff2,nh)
                          vw   = vw - vwr
                          v1   = v1 - v1r
                          v2   = v2 - v2r
                       end if
                       
                       ww   = q(ibc,jbc,kbc,4,n)
                       w1   = q(ibc+ioff1,jbc+joff1,kbc+koff1,4,n)
                       w2   = q(ibc+ioff2,jbc+joff2,kbc+koff2,4,n)
                       if (moving) then
                          wwr  = zdot(ibc      ,jbc      ,kbc      ,nh)
                          w1r  = zdot(ibc+ioff1,jbc+joff1,kbc+koff1,nh)
                          w2r  = zdot(ibc+ioff2,jbc+joff2,kbc+koff2,nh)
                          ww   = ww - wwr
                          w1   = w1 - w1r
                          w2   = w2 - w2r
                       end if
                       
                       ueta = (-1)**(1+1/inrout) * (-8*uw + 9*u1 - u2) / 3
                       veta = (-1)**(1+1/inrout) * (-8*vw + 9*v1 - v2) / 3
                       weta = (-1)**(1+1/inrout) * (-8*ww + 9*w1 - w2) / 3
                       
                       dudx = ueta*etax
                       dudy = ueta*etay
                       dudz = ueta*etaz
                       dvdx = veta*etax
                       dvdy = veta*etay
                       dvdz = veta*etaz
                       dwdx = weta*etax
                       dwdy = weta*etay
                       dwdz = weta*etaz
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
                       
                       cf(ibc,jbc,kbc) = 2*tauw / (reyno*machfs)
                       yp(ibc,jbc,kbc) = dsqrt(reyno/machfs)* &
                            rhow*dn1*utau/muw
                       
                    end if
                    
                 end do
              end do                 
              
              tempindex = 1
              
              do k=kbc_s,kbc_e
                 do j=jbc_s,jbc_e
                    do i=ibc_s,ibc_e
                       tempdata(tempindex) = x(i,j,k,n)
                       tempindex = tempindex + 1
                    end do
                 end do
              end do
              
              do k=kbc_s,kbc_e
                 do j=jbc_s,jbc_e
                    do i=ibc_s,ibc_e
                       tempdata(tempindex) = y(i,j,k,n)
                       tempindex = tempindex + 1
                    end do
                 end do
              end do
              
              do k=kbc_s,kbc_e
                 do j=jbc_s,jbc_e
                    do i=ibc_s,ibc_e
                       tempdata(tempindex) = z(i,j,k,n)
                       tempindex = tempindex + 1
                    end do
                 end do
              end do
              
              do k=kbc_s,kbc_e
                 do j=jbc_s,jbc_e
                    do i=ibc_s,ibc_e
                       tempdata(tempindex) = cp(i,j,k)
                       tempindex = tempindex + 1
                    end do
                 end do
              end do
              
              if(viscous) then
                 
                 do k=kbc_s,kbc_e
                    do j=jbc_s,jbc_e
                       do i=ibc_s,ibc_e
                          tempdata(tempindex) = cf(i,j,k)
                          tempindex = tempindex + 1
                       end do
                    end do
                 end do                    
                 
                 do k=kbc_s,kbc_e
                    do j=jbc_s,jbc_e
                       do i=ibc_s,ibc_e
                          tempdata(tempindex) = yp(i,j,k)
                          tempindex = tempindex + 1
                       end do
                    end do
                 end do                    

              else

                 do k=kbc_s,kbc_e
                    do j=jbc_s,jbc_e
                       do i=ibc_s,ibc_e
                          tempdata(tempindex) = 0.d0
                          tempindex = tempindex + 1
                       end do
                    end do
                 end do

                 do k=kbc_s,kbc_e
                    do j=jbc_s,jbc_e
                       do i=ibc_s,ibc_e
                          tempdata(tempindex) = 0.d0
                          tempindex = tempindex + 1
                       end do
                    end do
                 end do
                    
              end if                                       

              call setupfile(fid(n),disp,MPI_DOUBLE_PRECISION)
              call mpi_file_write(fid(n),tempdata(1),datasize, &
                   MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
              disp = disp + datasize*doublesize                
              
              if(allocated(tempdata)) deallocate(tempdata)
              if(allocated(cp)) deallocate(cp)
              if (viscous) then
                 if(allocated(cf)) deallocate(cf,yp)
              end if

              surfacenumber = surfacenumber + 1
              
           end if
        end do

      end do      

#endif

      return
      end

!-----------------------------------------------------------------------
      subroutine write_parallel_wr_plt_surf_b_M1(fid,imax, &
                 jmax,kmax,npde,nharms,iblock,nl,lmet,q,x,y,z,xdot,ydot, &
                  zdot,si,sj,sk,xideri,etaderj,zetaderk,dist,bctopo,nbcs)
!     This routine writes parts of the surftec files using tecplot.
!-----------------------------------------------------------------------

      use parallelutils, only: nsurfblocks, surfblockstarts, &
           surfblockends, surfblockindex, surfbcindex
      use common_variables
      use cosa_precision

      implicit none

#if MPI 
      include 'mpif.h'
#endif

      integer(kind=cosa_int) :: i,imax,jmax,kmax,npde,nharms,iblock,nl,lmet,nbcs
      integer(kind=cosa_int) :: fid(0:2*nharms)
      integer(kind=cosa_int) :: bctopo(10,nbcs)
      real(kind=cosa_real) :: &
           q   (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           x   ( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           y   ( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           z   ( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           xdot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           ydot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           zdot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           si      (lmet, 0:imax+1, 0:jmax+1, 0:kmax+1,0:2*nharms*hbmove), &
           sj      (lmet, 0:imax+1, 0:jmax+1, 0:kmax+1,0:2*nharms*hbmove), &
           sk      (lmet, 0:imax+1, 0:jmax+1, 0:kmax+1,0:2*nharms*hbmove), &
           xideri  (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           etaderj (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           zetaderk(3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           dist( 0:imax  , 0:jmax  , 0:kmax)
      integer(kind=cosa_int) :: integersize,doublesize,charactersize
      integer*8 disp
      integer(kind=cosa_int) :: surfacenumber

      disp = 0

      call surftecinitialoffset(disp,nl,iblock,nsurfblocks, &
           surfblockstarts,surfblockends,surfblockindex,surfbcindex, &
           surfacenumber,integersize,doublesize,charactersize)

      call parallel_wr_plt_surf_b_M1(fid,q,x,y,z,xdot,ydot,zdot,si,sj, &
           sk,xideri,etaderj,zetaderk,dist,bctopo,imax,jmax,kmax,npde, &
           nharms,nl,nbcs,lmet,surfblockstarts,surfblockends, &
           nsurfblocks,surfacenumber,iblock)

      return 
      end


!-----------------------------------------------------------------------
      subroutine parallel_wr_plt_surf_b_M1(fid,q,x,y,z,xdot,ydot,zdot, &
                 si,sj,sk,xideri,etaderj,zetaderk,dist,bctopo,imax,jmax, &
                 kmax,npde,nharms,nl,nbcs,lmet,surfacestartindexes, &
                 surfaceendindexes,nsurfblocks,surfacenumber,iblock)
!     This routine writes parts of the flowtec files using tecplot.  
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

#if MPI 
      include 'mpif.h'
#endif

      integer(kind=cosa_int) :: imax,jmax,kmax,npde,nharms,nl,lmet,nbcs
      integer(kind=cosa_int) :: i,imax1,j,jmax1,k,kmax1,ipde,n,nh,ibc1,iblock
      integer(kind=cosa_int) :: i1,i2,i3,ijkmax(3),idir,istrt(3),iend(3),bctyp, &
        inrout,ibcpt,ibcpt2,ibcn,ibcn2,ibcm,ibc,jbc,kbc,ibc2,jbc2, &
        kbc2,in,jn,kn,in2,jn2,kn2,im,jm,km,ic1,ic2,ic3,sysize,ioff1, &
        ioff2,joff1,joff2,koff1,koff2,ibc_s,jbc_s,kbc_s,ibc_e,jbc_e, &
        kbc_e
      integer(kind=cosa_int) :: fid(0:2*nharms)
      integer(kind=cosa_int) :: bctopo(10,nbcs)
      real(kind=cosa_real) :: sgnm,nx,ny,nz,kx,ky,kz,rhow,tw,muw,uw,u1,u2,uwr,u1r,u2r, &
           vw,v1,v2,vwr,v1r,v2r,ww,w1,w2,wwr,w1r,w2r,ueta,veta,weta, &
           etax,etay,etaz,dudx,dudy,dudz,dvdx,dvdy, &
           dvdz,dwdx,dwdy,dwdz,divv,txx,txy,txz,tyy,tyz,tzz,mach,tauwx, &
           tauwy,tauwz,tauwpx,tauwpy,tauwpz,dn1,tauw,utau
      real(kind=cosa_real) :: &
           q   (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           x   ( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           y   ( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           z   ( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           xdot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           ydot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           zdot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           si      (lmet, 0:imax+1, 0:jmax+1, 0:kmax+1,0:2*nharms*hbmove), &
           sj      (lmet, 0:imax+1, 0:jmax+1, 0:kmax+1,0:2*nharms*hbmove), &
           sk      (lmet, 0:imax+1, 0:jmax+1, 0:kmax+1,0:2*nharms*hbmove), &
           xideri  (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           etaderj (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           zetaderk(3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           dist( 0:imax  , 0:jmax  , 0:kmax)
      real(kind=cosa_real),allocatable :: cp(:,:,:),cf(:,:,:),yp(:,:,:),zeroarray(:)
      integer(kind=cosa_int) :: ierr
      integer(kind=cosa_int) :: tempindex,datasize
      integer(kind=cosa_int) :: surfacestartindexes(3,nsurfblocks)
      integer(kind=cosa_int) :: surfaceendindexes(3,nsurfblocks)
      integer(kind=cosa_int) :: initialsurfacenumber, nsurfblocks, surfacenumber
      character(len=5) :: blocknumname
      integer(kind=cosa_int) :: nfconns, fnmode, shrconn, isblock, valuelocation, isdouble
      integer(kind=cosa_int) :: tnfnodes, ncbfaces, tnbconns
      integer(kind=cosa_int) ::zonetype, strandid, parentzone
      integer(kind=cosa_int) :: teczne142, tecznemap142, tecdat142, tecfil142
      integer(kind=cosa_int) :: ibcmax, jbcmax, kbcmax
      integer(kind=cosa_int) :: imaxmax, jmaxmax, kmaxmax
      integer(kind=cosa_int) :: mpiid
      integer(kind=cosa_int) :: Null(*)
      POINTER   (NullPtr,Null)

      
#ifdef MPI

      call getmpiid(mpiid)    

! Specify that we are using an Ordered zone type
      zonetype = 0

! Zones don't have parents
      parentzone = 0

! We are not part of a strand
      strandid = 0
      nfconns = 0
      fnmode = 0
      shrconn = 0
      imaxmax = 0
      jmaxmax = 0
      kmaxmax = 0

! This specifies if we are writing in block or point format
! 1 is block format.  Binary tecplot files must be block format, 
! so this must be 1.
      isblock = 1
      valuelocation = 0
      isdouble = 1
      tnfnodes = 0
      ncbfaces = 0
      tnbconns = 0

      ijkmax(1) = imax
      ijkmax(2) = jmax
      ijkmax(3) = kmax

      imax1 = imax+1
      jmax1 = jmax+1
      kmax1 = kmax+1

      initialsurfacenumber = surfacenumber

      do n = 0,2*nharms
        nh = n*hbmove

!     This line is included because the n loop iterates through files 
!     so we need to ensure that each time round the loop we start and 
!     the correct place in the associated file.
        surfacenumber = initialsurfacenumber
           
        do ibc1=1,nbcs

           if (any(bctopo(1,ibc1).eq. &
                [1500,1501,1502,1503,1400,1401,1402,1403])) then
              
              ibc_s = surfacestartindexes(1,surfacenumber)
              jbc_s = surfacestartindexes(2,surfacenumber)
              kbc_s = surfacestartindexes(3,surfacenumber)
              
              ibc_e = surfaceendindexes(1,surfacenumber)
              jbc_e = surfaceendindexes(2,surfacenumber)
              kbc_e = surfaceendindexes(3,surfacenumber)
              
              datasize = 6
              if((ibc_e-ibc_s)+1 .ne. 0) datasize = datasize*((ibc_e-ibc_s)+1)
              if((jbc_e-jbc_s)+1 .ne. 0) datasize = datasize*((jbc_e-jbc_s)+1)
              if((kbc_e-kbc_s)+1 .ne. 0) datasize = datasize*((kbc_e-kbc_s)+1)

              allocate(cp(ibc_s:ibc_e,jbc_s:jbc_e,kbc_s:kbc_e))
              if (viscous) then
                 allocate(cf(ibc_s:ibc_e,jbc_s:jbc_e,kbc_s:kbc_e), &
                      yp(ibc_s:ibc_e,jbc_s:jbc_e,kbc_s:kbc_e))
              end if
                            
!     store boundary topology in mnemonic names
              bctyp    = bctopo(1,ibc1)
              idir     = bctopo(2,ibc1)
              inrout   = bctopo(3,ibc1)
              istrt(1) = bctopo(4,ibc1)
              iend(1)  = bctopo(5,ibc1)
              istrt(2) = bctopo(6,ibc1)
              iend(2)  = bctopo(7,ibc1)
              istrt(3) = bctopo(8,ibc1)
              iend(3)  = bctopo(9,ibc1)
              
!             set needed variables depending on whether the boundary is
!             the inner boundary (inrout = 1) or
!             the outer boundary (inrout > 1)
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
                 sgnm   = - 1.d0
              else
                 ibcpt  = ijkmax(idir)
                 ibcpt2 = ijkmax(idir) + 1
                 ibcn   = ijkmax(idir) - 1
                 ibcn2  = ijkmax(idir) - 2
                 ibcm   = ijkmax(idir)
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
 
              do i3 = istrt(ic3),iend(ic3)
                 do i2 = istrt(ic2),iend(ic2)
                    
                    ibc  = ibcpt *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
                    jbc  = ibcpt *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
                    kbc  = ibcpt *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)
                    
                    im   = ibcm  *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
                    jm   = ibcm  *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
                    km   = ibcm  *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)                   


                    if (idir.eq.1) then
                       nx    = si(1,im,jm,km,nh)
                       ny    = si(2,im,jm,km,nh)
                       nz    = si(3,im,jm,km,nh)
                    else if (idir.eq.2) then
                       nx    = sj(1,im,jm,km,nh)
                       ny    = sj(2,im,jm,km,nh)
                       nz    = sj(3,im,jm,km,nh)
                    else if (idir.eq.3) then
                       nx    = sk(1,im,jm,km,nh)
                       ny    = sk(2,im,jm,km,nh)
                       nz    = sk(3,im,jm,km,nh)
                    end if
                    
                    kx   = (-1)**(1+1/inrout) * nx
                    ky   = (-1)**(1+1/inrout) * ny
                    kz   = (-1)**(1+1/inrout) * nz
                    
!     AJ Removed the machfs division here because if we are not viscous 
!     AJ machfs is zero.
                    cp(ibc,jbc,kbc) = 2*(q(ibc,jbc,kbc,5,n)-1/gamma)                                        
                    
                    if (viscous) then

                       cp(ibc,jbc,kbc) = cp(ibc,jbc,kbc)/ &
                            machfs**2

                       
                       ibc2 = ibcpt2*krd(ic1,1) +i2*krd(ic2,1) +i3*krd(ic3,1)
                       jbc2 = ibcpt2*krd(ic1,2) +i2*krd(ic2,2) +i3*krd(ic3,2)
                       kbc2 = ibcpt2*krd(Ic1,3) +i2*krd(ic2,3) +i3*krd(ic3,3)
                       
                       in   = ibcn  *krd(ic1,1) +i2*krd(ic2,1) +i3*krd(ic3,1)
                       jn   = ibcn  *krd(ic1,2) +i2*krd(ic2,2) +i3*krd(ic3,2)
                       kn   = ibcn  *krd(ic1,3) +i2*krd(ic2,3) +i3*krd(ic3,3)
                       
                       in2  = ibcn2 *krd(ic1,1) +i2*krd(ic2,1) +i3*krd(ic3,1)
                       jn2  = ibcn2 *krd(ic1,2) +i2*krd(ic2,2) +i3*krd(ic3,2)
                       kn2  = ibcn2 *krd(ic1,3) +i2*krd(ic2,3) +i3*krd(ic3,3)
                       
                       if (idir.eq.1) then
                          etax  = xideri(1,im,jm,km,nh)
                          etay  = xideri(2,im,jm,km,nh)
                          etaz  = xideri(3,im,jm,km,nh)
                       else if (idir.eq.2) then
                          etax  = etaderj(1,im,jm,km,nh)
                          etay  = etaderj(2,im,jm,km,nh)
                          etaz  = etaderj(3,im,jm,km,nh)
                       else if (idir.eq.3) then
                          etax  = zetaderk(1,im,jm,km,nh)
                          etay  = zetaderk(2,im,jm,km,nh)
                          etaz  = zetaderk(3,im,jm,km,nh)
                       end if
                       
                       rhow = q(ibc,jbc,kbc,1,n)
                       tw   = q(ibc,jbc,kbc,5,n) *gamma / q(ibc,jbc,kbc,1,n)
                       muw  = (stemp+1)/(tw+stemp) * (tw**1.5d0)
                       
!------------------MSC, 11/03/2023: dist not computed for laminar now
                       if (kom.or.kom_bsl.or.kom_sst) then
                          dn1  = dist(in ,jn ,kn )
                       else
                          dn1  = 0
                       end if
                       
                       uw   = q(ibc,jbc,kbc,2,n)
                       u1   = q(ibc+ioff1,jbc+joff1,kbc+koff1,2,n)
                       u2   = q(ibc+ioff2,jbc+joff2,kbc+koff2,2,n)
                       if (moving) then
                          uwr  = xdot(ibc      ,jbc      ,kbc      ,nh)
                          u1r  = xdot(ibc+ioff1,jbc+joff1,kbc+koff1,nh)
                          u2r  = xdot(ibc+ioff2,jbc+joff2,kbc+koff2,nh)
                          uw   = uw - uwr
                          u1   = u1 - u1r
                          u2   = u2 - u2r
                       end if
                       
                       vw   = q(ibc,jbc,kbc,3,n)
                       v1   = q(ibc+ioff1,jbc+joff1,kbc+koff1,3,n)
                       v2   = q(ibc+ioff2,jbc+joff2,kbc+koff2,3,n)
                       if (moving) then
                          vwr  = ydot(ibc      ,jbc      ,kbc      ,nh)
                          v1r  = ydot(ibc+ioff1,jbc+joff1,kbc+koff1,nh)
                          v2r  = ydot(ibc+ioff2,jbc+joff2,kbc+koff2,nh)
                          vw   = vw - vwr
                          v1   = v1 - v1r
                          v2   = v2 - v2r
                       end if
                       
                       ww   = q(ibc,jbc,kbc,4,n)
                       w1   = q(ibc+ioff1,jbc+joff1,kbc+koff1,4,n)
                       w2   = q(ibc+ioff2,jbc+joff2,kbc+koff2,4,n)
                       if (moving) then
                          wwr  = zdot(ibc      ,jbc      ,kbc      ,nh)
                          w1r  = zdot(ibc+ioff1,jbc+joff1,kbc+koff1,nh)
                          w2r  = zdot(ibc+ioff2,jbc+joff2,kbc+koff2,nh)
                          ww   = ww - wwr
                          w1   = w1 - w1r
                          w2   = w2 - w2r
                       end if
                       
                       ueta = (-1)**(1+1/inrout) * (-8*uw + 9*u1 - u2) / 3
                       veta = (-1)**(1+1/inrout) * (-8*vw + 9*v1 - v2) / 3
                       weta = (-1)**(1+1/inrout) * (-8*ww + 9*w1 - w2) / 3
                       
                       dudx = ueta*etax
                       dudy = ueta*etay
                       dudz = ueta*etaz
                       dvdx = veta*etax
                       dvdy = veta*etay
                       dvdz = veta*etaz
                       dwdx = weta*etax
                       dwdy = weta*etay
                       dwdz = weta*etaz
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
                       
                       cf(ibc,jbc,kbc) = 2*tauw / (reyno*machfs)
                       yp(ibc,jbc,kbc) = dsqrt(reyno/machfs)* &
                            rhow*dn1*utau/muw
                       
                    end if
                    
                 end do
              end do                 
              
              if(tecfil142(n+1) .ne. 0) then
                 write(*,*) 'error calling tecfil142'
                 stop
              end if
              
              write (blocknumname, "(I5)") iblock
              
              ibcmax = ibc_e-ibc_s+1
              jbcmax = jbc_e-jbc_s+1
              kbcmax = kbc_e-kbc_s+1
              
              if(teczne142('block'//trim(blocknumname)//char(0),zonetype, ibcmax, jbcmax, kbcmax, &
                   imaxmax, jmaxmax, kmaxmax, simtime, strandid, parentzone, isblock, nfconns, & 
                   fnmode, tnfnodes, ncbfaces, tnbconns, Null, Null, Null, shrconn) .ne. 0) then
                 write(*,*) 'error setting up zone'
                 stop
              end if
              
              if(tecznemap142(1, mpiid) .ne. 0) then
                 write(*,*) 'error setting up parallel mapping for zone'
                 call abortmpi()
              end if              

              if(tecdat142(ibcmax*jbcmax*kbcmax,x(ibc_s,jbc_s,kbc_s,n),isdouble) .ne. 0) then
                 write(*,*) 'error writing block data'
              end if
              if(tecdat142(ibcmax*jbcmax*kbcmax,y(ibc_s,jbc_s,kbc_s,n),isdouble) .ne. 0) then
                 write(*,*) 'error writing block data'
              end if
              if(tecdat142(ibcmax*jbcmax*kbcmax,z(ibc_s,jbc_s,kbc_s,n),isdouble) .ne. 0) then
                 write(*,*) 'error writing block data'
              end if
              if(tecdat142(ibcmax*jbcmax*kbcmax,cp(ibc_s,jbc_s,kbc_s),isdouble) .ne. 0) then
                 write(*,*) 'error writing block data'
              end if
              
              if (viscous) then
                 if(tecdat142(ibcmax*jbcmax*kbcmax,cf(ibc_s,jbc_s,kbc_s),isdouble) .ne. 0) then
                    write(*,*) 'error writing block data'
                 end if
                 if(tecdat142(ibcmax*jbcmax*kbcmax,yp(ibc_s,jbc_s,kbc_s),isdouble) .ne. 0) then
                    write(*,*) 'error writing block data'
                 end if
              else
                 allocate(zeroarray(ibcmax*jbcmax*kbcmax))
                 zeroarray = 0.d0
                 if(tecdat142(ibcmax*jbcmax*kbcmax,zeroarray,isdouble) .ne. 0) then
                    write(*,*) 'error writing block data'
                 end if
                 if(tecdat142(ibcmax*jbcmax*kbcmax,zeroarray,isdouble) .ne. 0) then
                    write(*,*) 'error writing block data'
                 end if
                 deallocate(zeroarray)
              end if
              
              if(allocated(cp)) deallocate(cp)
              if (viscous) then
                 if(allocated(cf)) deallocate(cf,yp)
              end if

              surfacenumber = surfacenumber + 1
              
           end if
        end do

      end do      

#endif

      return
      end

! AJ TODO Merge the tec and cgns versions of this to save code duplication
!-----------------------------------------------------------------------
      subroutine write_parallel_wr_cgns_surf_b_M1(fid,imax, &
                 jmax,kmax,npde,nharms,iblock,nl,lmet,q,x,y,z,xdot,ydot, &
                 zdot,si,sj,sk,xideri,etaderj,zetaderk,dist,bctopo,nbcs, &
                 basenum)
!     This routine writes parts of the surface output files using MPI-I/O.
!     Mainly this routine wraps the bespoke code which writes the
!     data (in parallel_wr_surf_tec_b_M1) and also deals with setting up
!     the file and moving the file pointer or displacement to ensure the
!     block is written in the correct place.
!-----------------------------------------------------------------------

      use parallelutils, only: nsurfblocks, surfblockstarts, &
           surfblockends, surfblockindex, surfbcindex
      use common_variables
      use cosa_precision

      implicit none

#if MPI 
      include 'mpif.h'
#endif

      integer(kind=cosa_int) :: i,imax,jmax,kmax,npde,nharms,iblock,nl,lmet,nbcs
      integer(kind=cosa_int) :: fid(0:2*nharms)
      integer(kind=cosa_int) :: bctopo(10,nbcs)
      real(kind=cosa_real) :: &
           q   (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           x   ( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           y   ( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           z   ( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           xdot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           ydot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           zdot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           si      (lmet, 0:imax+1, 0:jmax+1, 0:kmax+1,0:2*nharms*hbmove), &
           sj      (lmet, 0:imax+1, 0:jmax+1, 0:kmax+1,0:2*nharms*hbmove), &
           sk      (lmet, 0:imax+1, 0:jmax+1, 0:kmax+1,0:2*nharms*hbmove), &
           xideri  (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           etaderj (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           zetaderk(3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           dist( 0:imax  , 0:jmax  , 0:kmax)
      integer(kind=cosa_int) :: integersize,doublesize,charactersize
      integer*8 disp
      integer(kind=cosa_int) :: surfacenumber
      integer(kind=cosa_int) :: basenum

!AJ We do not need the calculated offset here, but we do need the surfacenumber
!AJ variable. We could refactor and create a new subroutine that just produces that 
!AJ variable and not the displacement as well, but for the amount of compute time 
!AJ required it seems unnecessary
      call surftecinitialoffset(disp,nl,iblock,nsurfblocks, &
           surfblockstarts,surfblockends,surfblockindex,surfbcindex, &
           surfacenumber,integersize,doublesize,charactersize)

      call parallel_wr_cgns_surf_b_M1(fid,q,x,y,z,xdot,ydot,zdot,si,sj, &
           sk,xideri,etaderj,zetaderk,dist,bctopo,imax,jmax,kmax,npde, &
           nharms,nl,nbcs,lmet,integersize,doublesize, &
           charactersize,surfblockstarts,surfblockends, &
           nsurfblocks,surfacenumber,basenum)

      return 
      end

!-----------------------------------------------------------------------
      subroutine parallel_wr_cgns_surf_b_M1(fid,q,x,y,z,xdot,ydot,zdot, &
                 si,sj,sk,xideri,etaderj,zetaderk,dist,bctopo,imax,jmax, &
                 kmax,npde,nharms,nl,nbcs,lmet,integersize, &
                 doublesize,charactersize,surfacestartindexes, &
                 surfaceendindexes,nsurfblocks,surfacenumber,basenum)
!     This routine writes parts of the surface output files using MPI-I/O.  
!-----------------------------------------------------------------------

      use cgns
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) :: imax,jmax,kmax,npde,nharms,nl,lmet,nbcs
      integer(kind=cosa_int) :: i,imax1,j,jmax1,k,kmax1,ipde,n,nh,ibc1
      integer(kind=cosa_int) :: i1,i2,i3,ijkmax(3),idir,istrt(3),iend(3),bctyp, &
        inrout,ibcpt,ibcpt2,ibcn,ibcn2,ibcm,ibc,jbc,kbc,ibc2,jbc2, &
        kbc2,in,jn,kn,in2,jn2,kn2,im,jm,km,ic1,ic2,ic3,sysize,ioff1, &
        ioff2,joff1,joff2,koff1,koff2,ibc_s,jbc_s,kbc_s,ibc_e,jbc_e, &
        kbc_e
      integer(kind=cosa_int) :: fid(0:2*nharms)
      integer(kind=cosa_int) :: bctopo(10,nbcs)
      real(kind=cosa_real) :: sgnm,nx,ny,nz,kx,ky,kz,rhow,tw,muw,uw,u1,u2,uwr,u1r,u2r, &
           vw,v1,v2,vwr,v1r,v2r,ww,w1,w2,wwr,w1r,w2r,ueta,veta,weta, &
           etax,etay,etaz,dudx,dudy,dudz,dvdx,dvdy, &
           dvdz,dwdx,dwdy,dwdz,divv,txx,txy,txz,tyy,tyz,tzz,mach,tauwx, &
           tauwy,tauwz,tauwpx,tauwpy,tauwpz,dn1,tauw,utau
      real(kind=cosa_real) :: &
           q   (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           x   ( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           y   ( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           z   ( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           xdot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           ydot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           zdot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           si      (lmet, 0:imax+1, 0:jmax+1, 0:kmax+1,0:2*nharms*hbmove), &
           sj      (lmet, 0:imax+1, 0:jmax+1, 0:kmax+1,0:2*nharms*hbmove), &
           sk      (lmet, 0:imax+1, 0:jmax+1, 0:kmax+1,0:2*nharms*hbmove), &
           xideri  (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           etaderj (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           zetaderk(3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           dist( 0:imax  , 0:jmax  , 0:kmax)
      real(kind=cosa_real),allocatable :: cp(:,:,:),cf(:,:,:),yp(:,:,:)
      integer(kind=cosa_int) :: ierr,integersize,doublesize,charactersize
      integer(kind=cosa_int) :: tempindex,datasize
      integer(kind=cosa_int) :: surfacestartindexes(3,nsurfblocks)
      integer(kind=cosa_int) :: surfaceendindexes(3,nsurfblocks)
      integer(kind=cosa_int) :: nsurfblocks, surfacenumber
      integer(kind=cosa_int) :: xnum, ynum, znum, flownum, solnum
      integer(kind=cosa_int) :: basenum
      character*200 line1
      integer(cgsize_t) :: minrange(3), maxrange(3)     
      integer(kind=cosa_int) :: initialsurfacenumber

      ijkmax(1) = imax
      ijkmax(2) = jmax
      ijkmax(3) = kmax

      imax1 = imax+1
      jmax1 = jmax+1
      kmax1 = kmax+1

      initialsurfacenumber = surfacenumber

      do n = 0,2*nharms
        nh = n*hbmove
           
!     Required to deal with resetting everything when iterating harmonics.
        surfacenumber = initialsurfacenumber

        do ibc1=1,nbcs

           if (any(bctopo(1,ibc1).eq. &
                [1500,1501,1502,1503,1400,1401,1402,1403])) then
              
              ibc_s = surfacestartindexes(1,surfacenumber)
              jbc_s = surfacestartindexes(2,surfacenumber)
              kbc_s = surfacestartindexes(3,surfacenumber)
              
              ibc_e = surfaceendindexes(1,surfacenumber)
              jbc_e = surfaceendindexes(2,surfacenumber)
              kbc_e = surfaceendindexes(3,surfacenumber)
              
              minrange(1) = 1
              minrange(2) = 1
              minrange(3) = 1

              maxrange(1) = (ibc_e - ibc_s)+1
              maxrange(2) = (jbc_e - jbc_s)+1
              maxrange(3) = (kbc_e - kbc_s)+1

              allocate(cp(ibc_s:ibc_e,jbc_s:jbc_e,kbc_s:kbc_e))
              
              allocate(cf(ibc_s:ibc_e,jbc_s:jbc_e,kbc_s:kbc_e), &
                   yp(ibc_s:ibc_e,jbc_s:jbc_e,kbc_s:kbc_e))
                            
!     store boundary topology in mnemonic names
              bctyp    = bctopo(1,ibc1)
              idir     = bctopo(2,ibc1)
              inrout   = bctopo(3,ibc1)
              istrt(1) = bctopo(4,ibc1)
              iend(1)  = bctopo(5,ibc1)
              istrt(2) = bctopo(6,ibc1)
              iend(2)  = bctopo(7,ibc1)
              istrt(3) = bctopo(8,ibc1)
              iend(3)  = bctopo(9,ibc1)
              
!             set needed variables depending on whether the boundary is
!             the inner boundary (inrout = 1) or
!             the outer boundary (inrout > 1)
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
                 sgnm   = - 1.d0
              else
                 ibcpt  = ijkmax(idir)
                 ibcpt2 = ijkmax(idir) + 1
                 ibcn   = ijkmax(idir) - 1
                 ibcn2  = ijkmax(idir) - 2
                 ibcm   = ijkmax(idir)
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
 
              do i3 = istrt(ic3),iend(ic3)
                 do i2 = istrt(ic2),iend(ic2)
                    
                    ibc  = ibcpt *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
                    jbc  = ibcpt *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
                    kbc  = ibcpt *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)
                    
                    im   = ibcm  *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
                    jm   = ibcm  *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
                    km   = ibcm  *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)                   

                    if (idir.eq.1) then
                       nx    = si(1,im,jm,km,nh)
                       ny    = si(2,im,jm,km,nh)
                       nz    = si(3,im,jm,km,nh)
                    else if (idir.eq.2) then
                       nx    = sj(1,im,jm,km,nh)
                       ny    = sj(2,im,jm,km,nh)
                       nz    = sj(3,im,jm,km,nh)
                    else if (idir.eq.3) then
                       nx    = sk(1,im,jm,km,nh)
                       ny    = sk(2,im,jm,km,nh)
                       nz    = sk(3,im,jm,km,nh)
                    end if
                    
                    kx   = (-1)**(1+1/inrout) * nx
                    ky   = (-1)**(1+1/inrout) * ny
                    kz   = (-1)**(1+1/inrout) * nz
                    
!     AJ Removed the machfs division here because if we are not viscous 
!     AJ machfs is zero.
                    cp(ibc,jbc,kbc) = 2*(q(ibc,jbc,kbc,5,n)-1/gamma)                                        
                    
                    if (viscous) then

                       cp(ibc,jbc,kbc) = cp(ibc,jbc,kbc)/ &
                            machfs**2

                       
                       ibc2 = ibcpt2*krd(ic1,1) +i2*krd(ic2,1) +i3*krd(ic3,1)
                       jbc2 = ibcpt2*krd(ic1,2) +i2*krd(ic2,2) +i3*krd(ic3,2)
                       kbc2 = ibcpt2*krd(Ic1,3) +i2*krd(ic2,3) +i3*krd(ic3,3)
                       
                       in   = ibcn  *krd(ic1,1) +i2*krd(ic2,1) +i3*krd(ic3,1)
                       jn   = ibcn  *krd(ic1,2) +i2*krd(ic2,2) +i3*krd(ic3,2)
                       kn   = ibcn  *krd(ic1,3) +i2*krd(ic2,3) +i3*krd(ic3,3)
                       
                       in2  = ibcn2 *krd(ic1,1) +i2*krd(ic2,1) +i3*krd(ic3,1)
                       jn2  = ibcn2 *krd(ic1,2) +i2*krd(ic2,2) +i3*krd(ic3,2)
                       kn2  = ibcn2 *krd(ic1,3) +i2*krd(ic2,3) +i3*krd(ic3,3)
                       
                       if (idir.eq.1) then
                          etax  = xideri(1,im,jm,km,nh)
                          etay  = xideri(2,im,jm,km,nh)
                          etaz  = xideri(3,im,jm,km,nh)
                       else if (idir.eq.2) then
                          etax  = etaderj(1,im,jm,km,nh)
                          etay  = etaderj(2,im,jm,km,nh)
                          etaz  = etaderj(3,im,jm,km,nh)
                       else if (idir.eq.3) then
                          etax  = zetaderk(1,im,jm,km,nh)
                          etay  = zetaderk(2,im,jm,km,nh)
                          etaz  = zetaderk(3,im,jm,km,nh)
                       end if
                       
                       rhow = q(ibc,jbc,kbc,1,n)
                       tw   = q(ibc,jbc,kbc,5,n) *gamma / q(ibc,jbc,kbc,1,n)
                       muw  = (stemp+1)/(tw+stemp) * (tw**1.5d0)
                       
!------------------MSC, 11/03/2023: dist not computed for laminar now
                       if (kom.or.kom_bsl.or.kom_sst) then
                          dn1  = dist(in ,jn ,kn )
                       else
                          dn1  = 0
                       end if
                       
                       uw   = q(ibc,jbc,kbc,2,n)
                       u1   = q(ibc+ioff1,jbc+joff1,kbc+koff1,2,n)
                       u2   = q(ibc+ioff2,jbc+joff2,kbc+koff2,2,n)
                       if (moving) then
                          uwr  = xdot(ibc      ,jbc      ,kbc      ,nh)
                          u1r  = xdot(ibc+ioff1,jbc+joff1,kbc+koff1,nh)
                          u2r  = xdot(ibc+ioff2,jbc+joff2,kbc+koff2,nh)
                          uw   = uw - uwr
                          u1   = u1 - u1r
                          u2   = u2 - u2r
                       end if
                       
                       vw   = q(ibc,jbc,kbc,3,n)
                       v1   = q(ibc+ioff1,jbc+joff1,kbc+koff1,3,n)
                       v2   = q(ibc+ioff2,jbc+joff2,kbc+koff2,3,n)
                       if (moving) then
                          vwr  = ydot(ibc      ,jbc      ,kbc      ,nh)
                          v1r  = ydot(ibc+ioff1,jbc+joff1,kbc+koff1,nh)
                          v2r  = ydot(ibc+ioff2,jbc+joff2,kbc+koff2,nh)
                          vw   = vw - vwr
                          v1   = v1 - v1r
                          v2   = v2 - v2r
                       end if
                       
                       ww   = q(ibc,jbc,kbc,4,n)
                       w1   = q(ibc+ioff1,jbc+joff1,kbc+koff1,4,n)
                       w2   = q(ibc+ioff2,jbc+joff2,kbc+koff2,4,n)
                       if (moving) then
                          wwr  = zdot(ibc      ,jbc      ,kbc      ,nh)
                          w1r  = zdot(ibc+ioff1,jbc+joff1,kbc+koff1,nh)
                          w2r  = zdot(ibc+ioff2,jbc+joff2,kbc+koff2,nh)
                          ww   = ww - wwr
                          w1   = w1 - w1r
                          w2   = w2 - w2r
                       end if
                       
                       ueta = (-1)**(1+1/inrout) * (-8*uw + 9*u1 - u2) / 3
                       veta = (-1)**(1+1/inrout) * (-8*vw + 9*v1 - v2) / 3
                       weta = (-1)**(1+1/inrout) * (-8*ww + 9*w1 - w2) / 3
                       
                       dudx = ueta*etax
                       dudy = ueta*etay
                       dudz = ueta*etaz
                       dvdx = veta*etax
                       dvdy = veta*etay
                       dvdz = veta*etaz
                       dwdx = weta*etax
                       dwdy = weta*etay
                       dwdz = weta*etaz
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
                       
                       cf(ibc,jbc,kbc) = 2*tauw / (reyno*machfs)
                       yp(ibc,jbc,kbc) = dsqrt(reyno/machfs)* &
                            rhow*dn1*utau/muw
                       
                    end if
                    
                 end do
              end do                              
                            
              if(.not.viscous) then
                 
                 do k=kbc_s,kbc_e
                    do j=jbc_s,jbc_e
                       do i=ibc_s,ibc_e
                          cf(i,j,k) = 0.d0
                       end do
                    end do
                 end do                    
                 
                 do k=kbc_s,kbc_e
                    do j=jbc_s,jbc_e
                       do i=ibc_s,ibc_e
                          yp(i,j,k) = 0.d0
                       end do
                    end do
                 end do                    

              end if                                       

              xnum = 1
              call cgp_coord_write_data_f(fid(n),basenum,surfacenumber,xnum, &
                   minrange,maxrange,x(ibc_s:ibc_e,jbc_s:jbc_e,kbc_s:kbc_e,n),ierr)
              if(ierr .ne. CG_OK) then
                 write(*,*) 'X cgp_coord_write_data_f error'
                 call cg_error_print_f()
              end if
              ynum = 2
              call cgp_coord_write_data_f(fid(n),basenum,surfacenumber,ynum, &
                   minrange,maxrange,y(ibc_s:ibc_e,jbc_s:jbc_e,kbc_s:kbc_e,n),ierr)
              if(ierr .ne. CG_OK) then
                 write(*,*) 'Y cgp_coord_write_data_f error'
                 call cg_error_print_f()
              end if
              znum = 3
              call cgp_coord_write_data_f(fid(n),basenum,surfacenumber,znum, &
                   minrange,maxrange,z(ibc_s:ibc_e,jbc_s:jbc_e,kbc_s:kbc_e,n),ierr)
              if(ierr .ne. CG_OK) then
                 write(*,*) 'Z cgp_coord_write_data_f error'
                 call cg_error_print_f()
              end if
              solnum = 1
              flownum = 1
              call cgp_field_write_data_f(fid(n),basenum,surfacenumber,solnum, &
                   flownum,minrange,maxrange,cp(ibc_s:ibc_e,jbc_s:jbc_e,kbc_s:kbc_e),ierr)
              if(ierr .ne. CG_OK) then
                 write(*,*) 'cgp_field_write_data_f cp error'
                 call cg_error_print_f()
              end if
              flownum = flownum + 1
              call cgp_field_write_data_f(fid(n),basenum,surfacenumber,solnum, &
                   flownum,minrange,maxrange,cf(ibc_s:ibc_e,jbc_s:jbc_e,kbc_s:kbc_e),ierr)
              if(ierr .ne. CG_OK) then
                 write(*,*) 'cgp_field_write_data_f cf error'
                 call cg_error_print_f()
              end if
              flownum = flownum + 1
              call cgp_field_write_data_f(fid(n),basenum,surfacenumber,solnum, &
                   flownum,minrange,maxrange,yp(ibc_s:ibc_e,jbc_s:jbc_e,kbc_s:kbc_e),ierr)
              if(ierr .ne. CG_OK) then
                 write(*,*) 'cgp_field_write_data_f yp error'
                 call cg_error_print_f()
              end if
                            
              if(allocated(cp)) deallocate(cp)
              if(allocated(cf)) deallocate(cf,yp)

              surfacenumber = surfacenumber + 1
              
           end if
        end do

      end do      

      return
      end

! AJ TODO Merge the tec and cgns versions of this to save code duplication
!-----------------------------------------------------------------------
      subroutine write_parallel_wr_tec_surf_b_M2(fid,imax, &
                 jmax,kmax,npde,nharms,iblock,nl,lmet,q,x,y,z,xdot,ydot, &
                 zdot,si,sj,sk,xideri,xiderj,xiderk,etaderi,etaderj, &
                 etaderk,zetaderi,zetaderj,zetaderk,dist,bctopo,nbcs)
!     This routine writes parts of the surface output files using MPI-I/O.
!     Mainly this routine wraps the bespoke code which writes the
!     data (in parallel_wr_surf_tec_b_M2) and also deals with setting up
!     the file and moving the file pointer or displacement to ensure the
!     block is written in the correct place.
!-----------------------------------------------------------------------

      use parallelutils, only: nsurfblocks, surfblockstarts, &
           surfblockends, surfblockindex, surfbcindex
      use common_variables
      use cosa_precision

      implicit none

#if MPI 
      include 'mpif.h'
#endif

      integer(kind=cosa_int) :: i,imax,jmax,kmax,npde,nharms,iblock,nl,lmet,nbcs
      integer(kind=cosa_int) :: fid(0:2*nharms)
      integer(kind=cosa_int) :: bctopo(10,nbcs)
      real(kind=cosa_real) :: &
           q   (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           x   ( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           y   ( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           z   ( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           xdot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           ydot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           zdot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           si      (lmet, 0:imax+1, 0:jmax+1, 0:kmax+1,0:2*nharms*hbmove), &
           sj      (lmet, 0:imax+1, 0:jmax+1, 0:kmax+1,0:2*nharms*hbmove), &
           sk      (lmet, 0:imax+1, 0:jmax+1, 0:kmax+1,0:2*nharms*hbmove), &
           xideri  (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           xiderj  (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           xiderk  (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           etaderi (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           etaderj (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           etaderk (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           zetaderi(3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           zetaderj(3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           zetaderk(3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           dist( 0:imax  , 0:jmax  , 0:kmax)
      integer(kind=cosa_int) :: integersize,doublesize,charactersize
      integer*8 disp
      integer(kind=cosa_int) :: surfacenumber

      disp = 0

      call surftecinitialoffset(disp,nl,iblock,nsurfblocks, &
           surfblockstarts,surfblockends,surfblockindex,surfbcindex, &
           surfacenumber,integersize,doublesize,charactersize)

      call parallel_wr_tec_surf_b_M2(fid,q,x,y,z,xdot,ydot,zdot,si,sj, &
           sk,xideri,xiderj,xiderk,etaderi,etaderj,etaderk,zetaderi, &
           zetaderj,zetaderk,dist,bctopo,imax,jmax,kmax,npde, &
           nharms,nl,nbcs,lmet,integersize,doublesize, &
           charactersize,disp,surfblockstarts,surfblockends, &
           nsurfblocks,surfacenumber)

      return 
      end

!-----------------------------------------------------------------------
      subroutine parallel_wr_tec_surf_b_M2(fid,q,x,y,z,xdot,ydot,zdot, &
                 si,sj,sk,xideri,xiderj,xiderk,etaderi,etaderj,etaderk, &
                 zetaderi,zetaderj,zetaderk,dist,bctopo,imax,jmax, &
                 kmax,npde,nharms,nl,nbcs,lmet,integersize, &
                 doublesize,charactersize,disp,surfacestartindexes, &
                 surfaceendindexes,nsurfblocks,surfacenumber)
!     This routine writes parts of the surface output files using MPI-I/O
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

#if MPI 
      include 'mpif.h'
#endif

      integer(kind=cosa_int) :: imax,jmax,kmax,npde,nharms,nl,lmet,nbcs
      integer(kind=cosa_int) :: i,imax1,j,jmax1,k,kmax1,ipde,n,nh,ibc1
      integer(kind=cosa_int) :: i1,i2,i3,ijkmax(3),idir,istrt(3),iend(3),bctyp, &
        inrout,ibcpt,ibcpt2,ibcn,ibcn2,ibcm,ibc,jbc,kbc,ibc2,jbc2, &
        kbc2,in,jn,kn,in2,jn2,kn2,im,jm,km,ic1,ic2,ic3,sysize,ioff1, &
        ioff2,joff1,joff2,koff1,koff2,ibc_s,jbc_s,kbc_s,ibc_e,jbc_e, &
        kbc_e
      integer(kind=cosa_int) :: fid(0:2*nharms)
      integer(kind=cosa_int) :: bctopo(10,nbcs)
      real(kind=cosa_real) :: sgnm,nx,ny,nz,kx,ky,kz,rhow,tw,muw,uw,u1,u2,uwr,u1r,u2r, &
           vw,v1,v2,vwr,v1r,v2r,ww,w1,w2,wwr,w1r,w2r,ueta,veta,weta, &
           xix,xiy,xiz,etax,etay,etaz,zetax,zetay,zetaz,dudeta,dudzeta, &
           dudxi,dudx,dudy,dudz,dvdx,dvdy,dwdxi,dwdeta,dwdzeta,dvdeta, &
           dvdz,dwdx,dwdy,dwdz,divv,txx,txy,txz,tyy,tyz,tzz,mach,tauwx, &
           tauwy,tauwz,tauwpx,tauwpy,tauwpz,dn1,tauw,utau,dvdxi,dvdzeta
      real(kind=cosa_real) :: &
           q   (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           x   ( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           y   ( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           z   ( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           xdot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           ydot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           zdot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           si      (lmet, 0:imax+1, 0:jmax+1, 0:kmax+1,0:2*nharms*hbmove), &
           sj      (lmet, 0:imax+1, 0:jmax+1, 0:kmax+1,0:2*nharms*hbmove), &
           sk      (lmet, 0:imax+1, 0:jmax+1, 0:kmax+1,0:2*nharms*hbmove), &
           xideri  (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           xiderj  (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           xiderk  (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           etaderi (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           etaderj (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           etaderk (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           zetaderi(3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           zetaderj(3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           zetaderk(3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           dist( 0:imax  , 0:jmax  , 0:kmax)
      real(kind=cosa_real),allocatable :: cp(:,:,:),cf(:,:,:),yp(:,:,:)
      real(kind=cosa_real),allocatable :: tempdata(:)
      integer*8 disp,initialdisp
      integer(kind=cosa_int) :: ierr,integersize,doublesize,charactersize
      integer(kind=cosa_int) :: tempindex,datasize
      integer(kind=cosa_int) :: surfacestartindexes(3,nsurfblocks)
      integer(kind=cosa_int) :: surfaceendindexes(3,nsurfblocks)
      integer(kind=cosa_int) :: nsurfblocks, surfacenumber
      character*200 line1
      integer(kind=cosa_int) :: initialsurfacenumber
      
#ifdef MPI

      ijkmax(1) = imax
      ijkmax(2) = jmax
      ijkmax(3) = kmax

      imax1 = imax+1
      jmax1 = jmax+1
      kmax1 = kmax+1

      initialdisp = disp
      initialsurfacenumber = surfacenumber

      do n = 0,2*nharms
        nh = n*hbmove
           
!     This line is included because the n loop iterates through files 
!     so we need to ensure that each time round the loop we start and 
!     the correct place in the associated file.
        disp = initialdisp
        surfacenumber = initialsurfacenumber

        do ibc1=1,nbcs

           if (any(bctopo(1,ibc1).eq. &
                [1500,1501,1502,1503,1400,1401,1402,1403])) then
              
              ibc_s = surfacestartindexes(1,surfacenumber)
              jbc_s = surfacestartindexes(2,surfacenumber)
              kbc_s = surfacestartindexes(3,surfacenumber)
              
              ibc_e = surfaceendindexes(1,surfacenumber)
              jbc_e = surfaceendindexes(2,surfacenumber)
              kbc_e = surfaceendindexes(3,surfacenumber)
                                          
!     store boundary topology in mnemonic names
              bctyp    = bctopo(1,ibc1)
              idir     = bctopo(2,ibc1)
              inrout   = bctopo(3,ibc1)
              istrt(1) = bctopo(4,ibc1)
              iend(1)  = bctopo(5,ibc1)
              istrt(2) = bctopo(6,ibc1)
              iend(2)  = bctopo(7,ibc1)
              istrt(3) = bctopo(8,ibc1)
              iend(3)  = bctopo(9,ibc1)
              
!             set needed variables depending on whether the boundary is
!             the inner boundary (inrout = 1) or
!             the outer boundary (inrout > 1)
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

              datasize = 6
              if((ibc_e-ibc_s)+1 .ne. 0) datasize = datasize*((ibc_e-ibc_s)+1)
              if((jbc_e-jbc_s)+1 .ne. 0) datasize = datasize*((jbc_e-jbc_s)+1)
              if((kbc_e-kbc_s)+1 .ne. 0) datasize = datasize*((kbc_e-kbc_s)+1)

              allocate(tempdata(datasize))
              allocate(cp(ibc_s:ibc_e,jbc_s:jbc_e,kbc_s:kbc_e))
              if (viscous) then
                 allocate(cf(ibc_s:ibc_e,jbc_s:jbc_e,kbc_s:kbc_e), &
                      yp(ibc_s:ibc_e,jbc_s:jbc_e,kbc_s:kbc_e))
              end if      
 
              do i3 = istrt(ic3),iend(ic3)
                 do i2 = istrt(ic2),iend(ic2)
                    
                    ibc  = ibcpt *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
                    jbc  = ibcpt *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
                    kbc  = ibcpt *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)
                    
                    im   = ibcm  *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
                    jm   = ibcm  *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
                    km   = ibcm  *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)                   

                    if (idir.eq.1) then
                       nx    = si(1,im,jm,km,nh)
                       ny    = si(2,im,jm,km,nh)
                       nz    = si(3,im,jm,km,nh)
                    else if (idir.eq.2) then
                       nx    = sj(1,im,jm,km,nh)
                       ny    = sj(2,im,jm,km,nh)
                       nz    = sj(3,im,jm,km,nh)
                    else if (idir.eq.3) then
                       nx    = sk(1,im,jm,km,nh)
                       ny    = sk(2,im,jm,km,nh)
                       nz    = sk(3,im,jm,km,nh)
                    end if
                    
                    kx   = (-1)**(1+1/inrout) * nx
                    ky   = (-1)**(1+1/inrout) * ny
                    kz   = (-1)**(1+1/inrout) * nz
                    
!     AJ Removed the machfs division here because if we are not viscous 
!     AJ machfs is zero.
                    cp(ibc,jbc,kbc) = 2*(q(ibc,jbc,kbc,5,n)-1/gamma)                                        
                    
                    if (viscous) then

                       cp(ibc,jbc,kbc) = cp(ibc,jbc,kbc)/ &
                            machfs**2
                                              
                       in   = ibcn  *krd(ic1,1) +i2*krd(ic2,1) +i3*krd(ic3,1)
                       jn   = ibcn  *krd(ic1,2) +i2*krd(ic2,2) +i3*krd(ic3,2)
                       kn   = ibcn  *krd(ic1,3) +i2*krd(ic2,3) +i3*krd(ic3,3)
                       
                       in2  = ibcn2 *krd(ic1,1) +i2*krd(ic2,1) +i3*krd(ic3,1)
                       jn2  = ibcn2 *krd(ic1,2) +i2*krd(ic2,2) +i3*krd(ic3,2)
                       kn2  = ibcn2 *krd(ic1,3) +i2*krd(ic2,3) +i3*krd(ic3,3)
                       
                       rhow = q(ibc,jbc,kbc,1,n)
                       tw   = q(ibc,jbc,kbc,5,n) *gamma / q(ibc,jbc,kbc,1,n)
                       muw  = (stemp+1)/(tw+stemp) * (tw**1.5d0)
                       
!------------------MSC, 11/03/2023: dist not computed for laminar now
                       if (kom.or.kom_bsl.or.kom_sst) then
                          dn1  = dist(in ,jn ,kn )
                       else
                          dn1  = 0
                       end if

                       if (idir.eq.1) then
                          
                          dudxi = (-1)**(1+1/inrout) * &
                               (-8*q(ibc,jbc,kbc,2,n) + &
                               9*q(ibc+ioff1,jbc,kbc,2,n) - &
                               q(ibc+ioff2,jbc,kbc,2,n) ) / 3
                          dvdxi = (-1)**(1+1/inrout) * &
                               (-8*q(ibc,jbc,kbc,3,n) + &
                               9*q(ibc+ioff1,jbc,kbc,3,n) - &
                               q(ibc+ioff2,jbc,kbc,3,n) ) / 3
                          dwdxi = (-1)**(1+1/inrout) * &
                               (-8*q(ibc,jbc,kbc,4,n) + &
                               9*q(ibc+ioff1,jbc,kbc,4,n) - &
                               q(ibc+ioff2,jbc,kbc,4,n) ) / 3
                          
                          if (i2.eq.istrt(ic2)) then
                             u2 = (  q(ibc,jbc,kbc,2,n) + q(ibc,jbc+1,kbc,2,n)) / 2
                             u1 = (3*q(ibc,jbc,kbc,2,n) - q(ibc,jbc+1,kbc,2,n)) / 2
                             v2 = (  q(ibc,jbc,kbc,3,n) + q(ibc,jbc+1,kbc,3,n)) / 2
                             v1 = (3*q(ibc,jbc,kbc,3,n) - q(ibc,jbc+1,kbc,3,n)) / 2
                             w2 = (  q(ibc,jbc,kbc,4,n) + q(ibc,jbc+1,kbc,4,n)) / 2
                             w1 = (3*q(ibc,jbc,kbc,4,n) - q(ibc,jbc+1,kbc,4,n)) / 2
                          else if (i2.eq.iend(ic2)) then
                             u2 = (3*q(ibc,jbc,kbc,2,n) - q(ibc,jbc-1,kbc,2,n)) / 2
                             u1 = (  q(ibc,jbc,kbc,2,n) + q(ibc,jbc-1,kbc,2,n)) / 2
                             v2 = (3*q(ibc,jbc,kbc,3,n) - q(ibc,jbc-1,kbc,3,n)) / 2
                             v1 = (  q(ibc,jbc,kbc,3,n) + q(ibc,jbc-1,kbc,3,n)) / 2
                             w2 = (3*q(ibc,jbc,kbc,4,n) - q(ibc,jbc-1,kbc,4,n)) / 2
                             w1 = (  q(ibc,jbc,kbc,4,n) + q(ibc,jbc-1,kbc,4,n)) / 2
                          else
                             u2 = (q(ibc,jbc  ,kbc,2,n) + q(ibc,jbc+1,kbc,2,n)) / 2
                             u1 = (q(ibc,jbc-1,kbc,2,n) + q(ibc,jbc  ,kbc,2,n)) / 2
                             v2 = (q(ibc,jbc  ,kbc,3,n) + q(ibc,jbc+1,kbc,3,n)) / 2
                             v1 = (q(ibc,jbc-1,kbc,3,n) + q(ibc,jbc  ,kbc,3,n)) / 2
                             w2 = (q(ibc,jbc  ,kbc,4,n) + q(ibc,jbc+1,kbc,4,n)) / 2
                             w1 = (q(ibc,jbc-1,kbc,4,n) + q(ibc,jbc  ,kbc,4,n)) / 2
                          end if
                          dudeta = u2 - u1
                          dvdeta = v2 - v1
                          dwdeta = w2 - w1
                          
                          if (i3.eq.istrt(ic3)) then
                             u2 = (  q(ibc,jbc,kbc,2,n) + q(ibc,jbc,kbc+1,2,n)) / 2
                             u1 = (3*q(ibc,jbc,kbc,2,n) - q(ibc,jbc,kbc+1,2,n)) / 2
                             v2 = (  q(ibc,jbc,kbc,3,n) + q(ibc,jbc,kbc+1,3,n)) / 2
                             v1 = (3*q(ibc,jbc,kbc,3,n) - q(ibc,jbc,kbc+1,3,n)) / 2
                             w2 = (  q(ibc,jbc,kbc,4,n) + q(ibc,jbc,kbc+1,4,n)) / 2
                             w1 = (3*q(ibc,jbc,kbc,4,n) - q(ibc,jbc,kbc+1,4,n)) / 2
                          else if (i3.eq.iend(ic3)) then
                             u2 = (3*q(ibc,jbc,kbc,2,n) - q(ibc,jbc,kbc-1,2,n)) / 2
                             u1 = (  q(ibc,jbc,kbc,2,n) + q(ibc,jbc,kbc-1,2,n)) / 2
                             v2 = (3*q(ibc,jbc,kbc,3,n) - q(ibc,jbc,kbc-1,3,n)) / 2
                             v1 = (  q(ibc,jbc,kbc,3,n) + q(ibc,jbc,kbc-1,3,n)) / 2
                             w2 = (3*q(ibc,jbc,kbc,4,n) - q(ibc,jbc,kbc-1,4,n)) / 2
                             w1 = (  q(ibc,jbc,kbc,4,n) + q(ibc,jbc,kbc-1,4,n)) / 2
                          else
                             u2 = (q(ibc,jbc,kbc  ,2,n) + q(ibc,jbc,kbc+1,2,n)) / 2
                             u1 = (q(ibc,jbc,kbc-1,2,n) + q(ibc,jbc,kbc  ,2,n)) / 2
                             v2 = (q(ibc,jbc,kbc  ,3,n) + q(ibc,jbc,kbc+1,3,n)) / 2
                             v1 = (q(ibc,jbc,kbc-1,3,n) + q(ibc,jbc,kbc  ,3,n)) / 2
                             w2 = (q(ibc,jbc,kbc  ,4,n) + q(ibc,jbc,kbc+1,4,n)) / 2
                             w1 = (q(ibc,jbc,kbc-1,4,n) + q(ibc,jbc,kbc  ,4,n)) / 2
                          end if
                          dudzeta = u2 - u1
                          dvdzeta = v2 - v1
                          dwdzeta = w2 - w1
                          
                          xix   = xideri  (1,im,jm,km,n)
                          xiy   = xideri  (2,im,jm,km,n)
                          xiz   = xideri  (3,im,jm,km,n)
                          etax  = etaderi (1,im,jm,km,n)
                          etay  = etaderi (2,im,jm,km,n)
                          etaz  = etaderi (3,im,jm,km,n)
                          zetax = zetaderi(1,im,jm,km,n)
                          zetay = zetaderi(2,im,jm,km,n)
                          zetaz = zetaderi(3,im,jm,km,n)
                          
                       else if (idir.eq.2) then
                          
                          if (i3.eq.istrt(ic3)) then
                             u2 = (  q(ibc+1,jbc,kbc,2,n) + q(ibc  ,jbc,kbc,2,n))/2
                             u1 = (3*q(ibc  ,jbc,kbc,2,n) - q(ibc+1,jbc,kbc,2,n))/2
                             v2 = (  q(ibc+1,jbc,kbc,3,n) + q(ibc  ,jbc,kbc,3,n))/2
                             v1 = (3*q(ibc  ,jbc,kbc,3,n) - q(ibc+1,jbc,kbc,3,n))/2
                             w2 = (  q(ibc+1,jbc,kbc,4,n) + q(ibc  ,jbc,kbc,4,n))/2
                             w1 = (3*q(ibc  ,jbc,kbc,4,n) - q(ibc+1,jbc,kbc,4,n))/2
                          else if (i3.eq.iend(ic3)) then
                             u2 = (3*q(ibc,jbc,kbc,2,n) - q(ibc-1,jbc,kbc,2,n)) / 2
                             u1 = (  q(ibc,jbc,kbc,2,n) + q(ibc-1,jbc,kbc,2,n)) / 2
                             v2 = (3*q(ibc,jbc,kbc,3,n) - q(ibc-1,jbc,kbc,3,n)) / 2
                             v1 = (  q(ibc,jbc,kbc,3,n) + q(ibc-1,jbc,kbc,3,n)) / 2
                             w2 = (3*q(ibc,jbc,kbc,4,n) - q(ibc-1,jbc,kbc,4,n)) / 2
                             w1 = (  q(ibc,jbc,kbc,4,n) + q(ibc-1,jbc,kbc,4,n)) / 2
                          else
                             u2 = (q(ibc+1,jbc,kbc,2,n) + q(ibc  ,jbc,kbc,2,n)) / 2
                             u1 = (q(ibc  ,jbc,kbc,2,n) + q(ibc-1,jbc,kbc,2,n)) / 2
                             v2 = (q(ibc+1,jbc,kbc,3,n) + q(ibc  ,jbc,kbc,3,n)) / 2
                             v1 = (q(ibc  ,jbc,kbc,3,n) + q(ibc-1,jbc,kbc,3,n)) / 2
                             w2 = (q(ibc+1,jbc,kbc,4,n) + q(ibc  ,jbc,kbc,4,n)) / 2
                             w1 = (q(ibc  ,jbc,kbc,4,n) + q(ibc-1,jbc,kbc,4,n)) / 2
                          end if
                          dudxi = u2 - u1
                          dvdxi = v2 - v1
                          dwdxi = w2 - w1
                          
                          dudeta = (-1)**(1+1/inrout) * &
                               (-8*q(ibc,jbc,kbc,2,n) + &
                               9*q(ibc,jbc+joff1,kbc,2,n) - &
                               q(ibc,jbc+joff2,kbc,2,n) ) / 3
                          dvdeta = (-1)**(1+1/inrout) * &
                               (-8*q(ibc,jbc,kbc,3,n) + &
                               9*q(ibc,jbc+joff1,kbc,3,n) - &
                               q(ibc,jbc+joff2,kbc,3,n) ) / 3
                          dwdeta = (-1)**(1+1/inrout) * &
                               (-8*q(ibc,jbc,kbc,4,n) + &
                               9*q(ibc,jbc+joff1,kbc,4,n) - &
                               q(ibc,jbc+joff2,kbc,4,n) ) / 3
                          
                          if (i2.eq.istrt(ic2)) then
                             u2 = (  q(ibc,jbc,kbc+1,2,n) + q(ibc,jbc,kbc  ,2,n))/2
                             u1 = (3*q(ibc,jbc,kbc  ,2,n) - q(ibc,jbc,kbc+1,2,n))/2
                             v2 = (  q(ibc,jbc,kbc+1,3,n) + q(ibc,jbc,kbc  ,3,n))/2
                             v1 = (3*q(ibc,jbc,kbc  ,3,n) - q(ibc,jbc,kbc+1,3,n))/2
                             w2 = (  q(ibc,jbc,kbc+1,4,n) + q(ibc,jbc,kbc  ,4,n))/2
                             w1 = (3*q(ibc,jbc,kbc  ,4,n) - q(ibc,jbc,kbc+1,4,n))/2
                          else if (i2.eq.iend(ic2)) then
                             u2 = (3*q(ibc,jbc,kbc,2,n) - q(ibc,jbc,kbc-1,2,n)) / 2
                             u1 = (  q(ibc,jbc,kbc,2,n) + q(ibc,jbc,kbc-1,2,n)) / 2
                             v2 = (3*q(ibc,jbc,kbc,3,n) - q(ibc,jbc,kbc-1,3,n)) / 2
                             v1 = (  q(ibc,jbc,kbc,3,n) + q(ibc,jbc,kbc-1,3,n)) / 2
                             w2 = (3*q(ibc,jbc,kbc,4,n) - q(ibc,jbc,kbc-1,4,n)) / 2
                             w1 = (  q(ibc,jbc,kbc,4,n) + q(ibc,jbc,kbc-1,4,n)) / 2
                          else
                             u2 = (q(ibc,jbc,kbc+1,2,n) + q(ibc,jbc,kbc  ,2,n)) / 2
                             u1 = (q(ibc,jbc,kbc  ,2,n) + q(ibc,jbc,kbc-1,2,n)) / 2
                             v2 = (q(ibc,jbc,kbc+1,3,n) + q(ibc,jbc,kbc  ,3,n)) / 2
                             v1 = (q(ibc,jbc,kbc  ,3,n) + q(ibc,jbc,kbc-1,3,n)) / 2
                             w2 = (q(ibc,jbc,kbc+1,4,n) + q(ibc,jbc,kbc  ,4,n)) / 2
                             w1 = (q(ibc,jbc,kbc  ,4,n) + q(ibc,jbc,kbc-1,4,n)) / 2
                          end if
                          dudzeta = u2 - u1
                          dvdzeta = v2 - v1
                          dwdzeta = w2 - w1
                          
                          xix   = xiderj  (1,im,jm,km,n)
                          xiy   = xiderj  (2,im,jm,km,n)
                          xiz   = xiderj  (3,im,jm,km,n)
                          etax  = etaderj (1,im,jm,km,n)
                          etay  = etaderj (2,im,jm,km,n)
                          etaz  = etaderj (3,im,jm,km,n)
                          zetax = zetaderj(1,im,jm,km,n)
                          zetay = zetaderj(2,im,jm,km,n)
                          zetaz = zetaderj(3,im,jm,km,n)
                          
                       else if (idir.eq.3) then
                          if (i2.eq.istrt(ic2)) then
                             u2 = (  q(ibc,jbc,kbc,2,n) + q(ibc+1,jbc,kbc,2,n)) / 2
                             u1 = (3*q(ibc,jbc,kbc,2,n) - q(ibc+1,jbc,kbc,2,n)) / 2
                             v2 = (  q(ibc,jbc,kbc,3,n) + q(ibc+1,jbc,kbc,3,n)) / 2
                             v1 = (3*q(ibc,jbc,kbc,3,n) - q(ibc+1,jbc,kbc,3,n)) / 2
                             w2 = (  q(ibc,jbc,kbc,4,n) + q(ibc+1,jbc,kbc,4,n)) / 2
                             w1 = (3*q(ibc,jbc,kbc,4,n) - q(ibc+1,jbc,kbc,4,n)) / 2
                          else if (i2.eq.iend(ic2)) then
                             u2 = (3*q(ibc,jbc,kbc,2,n) - q(ibc-1,jbc,kbc,2,n)) / 2
                             u1 = (  q(ibc,jbc,kbc,2,n) + q(ibc-1,jbc,kbc,2,n)) / 2
                             v2 = (3*q(ibc,jbc,kbc,3,n) - q(ibc-1,jbc,kbc,3,n)) / 2
                             v1 = (  q(ibc,jbc,kbc,3,n) + q(ibc-1,jbc,kbc,3,n)) / 2
                             w2 = (3*q(ibc,jbc,kbc,4,n) - q(ibc-1,jbc,kbc,4,n)) / 2
                             w1 = (  q(ibc,jbc,kbc,4,n) + q(ibc-1,jbc,kbc,4,n)) / 2
                          else
                             u2 = (q(ibc  ,jbc,kbc,2,n) + q(ibc+1,jbc,kbc,2,n)) / 2
                             u1 = (q(ibc-1,jbc,kbc,2,n) + q(ibc  ,jbc,kbc,2,n)) / 2
                             v2 = (q(ibc  ,jbc,kbc,3,n) + q(ibc+1,jbc,kbc,3,n)) / 2
                             v1 = (q(ibc-1,jbc,kbc,3,n) + q(ibc  ,jbc,kbc,3,n)) / 2
                             w2 = (q(ibc  ,jbc,kbc,4,n) + q(ibc+1,jbc,kbc,4,n)) / 2
                             w1 = (q(ibc-1,jbc,kbc,4,n) + q(ibc  ,jbc,kbc,4,n)) / 2
                          end if
                          dudxi = u2 - u1
                          dvdxi = v2 - v1
                          dwdxi = w2 - w1
                          
                          if (i3.eq.istrt(ic3)) then
                             u2 = (  q(ibc,jbc,kbc,2,n) + q(ibc,jbc+1,kbc,2,n)) / 2
                             u1 = (3*q(ibc,jbc,kbc,2,n) - q(ibc,jbc+1,kbc,2,n)) / 2
                             v2 = (  q(ibc,jbc,kbc,3,n) + q(ibc,jbc+1,kbc,3,n)) / 2
                             v1 = (3*q(ibc,jbc,kbc,3,n) - q(ibc,jbc+1,kbc,3,n)) / 2
                             w2 = (  q(ibc,jbc,kbc,4,n) + q(ibc,jbc+1,kbc,4,n)) / 2
                             w1 = (3*q(ibc,jbc,kbc,4,n) - q(ibc,jbc+1,kbc,4,n)) / 2
                          else if (i3.eq.iend(ic3)) then
                             u2 = (3*q(ibc,jbc,kbc,2,n) - q(ibc,jbc-1,kbc,2,n)) / 2
                             u1 = (  q(ibc,jbc,kbc,2,n) + q(ibc,jbc-1,kbc,2,n)) / 2
                             v2 = (3*q(ibc,jbc,kbc,3,n) - q(ibc,jbc-1,kbc,3,n)) / 2
                             v1 = (  q(ibc,jbc,kbc,3,n) + q(ibc,jbc-1,kbc,3,n)) / 2
                             w2 = (3*q(ibc,jbc,kbc,4,n) - q(ibc,jbc-1,kbc,4,n)) / 2
                             w1 = (  q(ibc,jbc,kbc,4,n) + q(ibc,jbc-1,kbc,4,n)) / 2
                          else
                             u2 = (q(ibc,jbc  ,kbc,2,n) + q(ibc,jbc+1,kbc,2,n)) / 2
                             u1 = (q(ibc,jbc-1,kbc,2,n) + q(ibc,jbc  ,kbc,2,n)) / 2
                             v2 = (q(ibc,jbc  ,kbc,3,n) + q(ibc,jbc+1,kbc,3,n)) / 2
                             v1 = (q(ibc,jbc-1,kbc,3,n) + q(ibc,jbc  ,kbc,3,n)) / 2
                             w2 = (q(ibc,jbc  ,kbc,4,n) + q(ibc,jbc+1,kbc,4,n)) / 2
                             w1 = (q(ibc,jbc-1,kbc,4,n) + q(ibc,jbc  ,kbc,4,n)) / 2
                          end if
                          dudeta = u2 - u1
                          dvdeta = v2 - v1
                          dwdeta = w2 - w1
                          
                          dudzeta = (-1)**(1+1/inrout) * &
                               (-8*q(ibc,jbc,kbc,2,n) + &
                               9*q(ibc,jbc,kbc+koff1,2,n) - &
                               q(ibc,jbc,kbc+koff2,2,n) ) / 3
                          dvdzeta = (-1)**(1+1/inrout) * &
                               (-8*q(ibc,jbc,kbc,3,n) + &
                               9*q(ibc,jbc,kbc+koff1,3,n) - &
                               q(ibc,jbc,kbc+koff2,3,n) ) / 3
                          dwdzeta = (-1)**(1+1/inrout) * &
                               (-8*q(ibc,jbc,kbc,4,n) + &
                               9*q(ibc,jbc,kbc+koff1,4,n) - &
                               q(ibc,jbc,kbc+koff2,4,n) ) / 3
                          
                          xix   = xiderk  (1,im,jm,km,n)
                          xiy   = xiderk  (2,im,jm,km,n)
                          xiz   = xiderk  (3,im,jm,km,n)
                          etax  = etaderk (1,im,jm,km,n)
                          etay  = etaderk (2,im,jm,km,n)
                          etaz  = etaderk (3,im,jm,km,n)
                          zetax = zetaderk(1,im,jm,km,n)
                          zetay = zetaderk(2,im,jm,km,n)
                          zetaz = zetaderk(3,im,jm,km,n)
                          
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
                       
                       cf(ibc,jbc,kbc) = 2*tauw / (reyno*machfs)
                       yp(ibc,jbc,kbc) = dsqrt(reyno/machfs)* &
                            rhow*dn1*utau/muw
                       
                    end if
                    
                 end do
              end do

              
              tempindex = 1
              
              do k=kbc_s,kbc_e
                 do j=jbc_s,jbc_e
                    do i=ibc_s,ibc_e
                       tempdata(tempindex) = x(i,j,k,n)
                       tempindex = tempindex + 1
                    end do
                 end do
              end do
              
              do k=kbc_s,kbc_e
                 do j=jbc_s,jbc_e
                    do i=ibc_s,ibc_e
                       tempdata(tempindex) = y(i,j,k,n)
                       tempindex = tempindex + 1
                    end do
                 end do
              end do
              
              do k=kbc_s,kbc_e
                 do j=jbc_s,jbc_e
                    do i=ibc_s,ibc_e
                       tempdata(tempindex) = z(i,j,k,n)
                       tempindex = tempindex + 1
                    end do
                 end do
              end do
              
              do k=kbc_s,kbc_e
                 do j=jbc_s,jbc_e
                    do i=ibc_s,ibc_e
                       tempdata(tempindex) = cp(i,j,k)
                       tempindex = tempindex + 1
                    end do
                 end do
              end do
              
              if(viscous) then
                 
                 do k=kbc_s,kbc_e
                    do j=jbc_s,jbc_e
                       do i=ibc_s,ibc_e
                          tempdata(tempindex) = cf(i,j,k)
                          tempindex = tempindex + 1
                       end do
                    end do
                 end do                    
                 
                 do k=kbc_s,kbc_e
                    do j=jbc_s,jbc_e
                       do i=ibc_s,ibc_e
                          tempdata(tempindex) = yp(i,j,k)
                          tempindex = tempindex + 1
                       end do
                    end do
                 end do                    

              else

                 do k=kbc_s,kbc_e
                    do j=jbc_s,jbc_e
                       do i=ibc_s,ibc_e
                          tempdata(tempindex) = 0.d0
                          tempindex = tempindex + 1
                       end do
                    end do
                 end do

                 do k=kbc_s,kbc_e
                    do j=jbc_s,jbc_e
                       do i=ibc_s,ibc_e
                          tempdata(tempindex) = 0.d0
                          tempindex = tempindex + 1
                       end do
                    end do
                 end do
                    
              end if                          
 
              call setupfile(fid(n),disp,MPI_DOUBLE_PRECISION)
              call mpi_file_write(fid(n),tempdata(1),datasize, &
                   MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
              disp = disp + datasize*doublesize                
              
              if(allocated(tempdata)) deallocate(tempdata)
              if(allocated(cp)) deallocate(cp)
              if (viscous) then
                 if(allocated(cf)) deallocate(cf,yp)
              end if

              surfacenumber = surfacenumber + 1
              
           end if
        end do

      end do      

#endif

      return
      end

! AJ TODO Merge the tec and cgns versions of this to save code duplication
!-----------------------------------------------------------------------
      subroutine write_parallel_wr_plt_surf_b_M2(fid,imax, &
                 jmax,kmax,npde,nharms,iblock,nl,lmet,q,x,y,z,xdot,ydot, &
                 zdot,si,sj,sk,xideri,xiderj,xiderk,etaderi,etaderj, &
                 etaderk,zetaderi,zetaderj,zetaderk,dist,bctopo,nbcs)
!     This routine writes parts of the surface output files using tecplot.
!-----------------------------------------------------------------------

      use parallelutils, only: nsurfblocks, surfblockstarts, &
           surfblockends, surfblockindex, surfbcindex
      use common_variables
      use cosa_precision

      implicit none

#if MPI 
      include 'mpif.h'
#endif

      integer(kind=cosa_int) :: i,imax,jmax,kmax,npde,nharms,iblock,nl,lmet,nbcs
      integer(kind=cosa_int) :: fid(0:2*nharms)
      integer(kind=cosa_int) :: bctopo(10,nbcs)
      real(kind=cosa_real) :: &
           q   (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           x   ( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           y   ( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           z   ( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           xdot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           ydot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           zdot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           si      (lmet, 0:imax+1, 0:jmax+1, 0:kmax+1,0:2*nharms*hbmove), &
           sj      (lmet, 0:imax+1, 0:jmax+1, 0:kmax+1,0:2*nharms*hbmove), &
           sk      (lmet, 0:imax+1, 0:jmax+1, 0:kmax+1,0:2*nharms*hbmove), &
           xideri  (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           xiderj  (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           xiderk  (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           etaderi (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           etaderj (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           etaderk (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           zetaderi(3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           zetaderj(3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           zetaderk(3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           dist( 0:imax  , 0:jmax  , 0:kmax)
      integer(kind=cosa_int) :: integersize,doublesize,charactersize
      integer*8 disp
      integer(kind=cosa_int) :: surfacenumber

      disp = 0

      call surftecinitialoffset(disp,nl,iblock,nsurfblocks, &
           surfblockstarts,surfblockends,surfblockindex,surfbcindex, &
           surfacenumber,integersize,doublesize,charactersize)

      call parallel_wr_plt_surf_b_M2(fid,q,x,y,z,xdot,ydot,zdot,si,sj, &
           sk,xideri,xiderj,xiderk,etaderi,etaderj,etaderk,zetaderi, &
           zetaderj,zetaderk,dist,bctopo,imax,jmax,kmax,npde, &
           nharms,nl,nbcs,lmet,surfblockstarts,surfblockends, &
           nsurfblocks,surfacenumber,iblock)

      return 
      end

!-----------------------------------------------------------------------
      subroutine parallel_wr_plt_surf_b_M2(fid,q,x,y,z,xdot,ydot,zdot, &
                 si,sj,sk,xideri,xiderj,xiderk,etaderi,etaderj,etaderk, &
                 zetaderi,zetaderj,zetaderk,dist,bctopo,imax,jmax, &
                 kmax,npde,nharms,nl,nbcs,lmet,surfacestartindexes, &
                 surfaceendindexes,nsurfblocks,surfacenumber,iblock)
!     This routine writes parts of the surface output files using tecplot
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

#if MPI 
      include 'mpif.h'
#endif

      integer(kind=cosa_int) :: imax,jmax,kmax,npde,nharms,nl,lmet,nbcs
      integer(kind=cosa_int) :: i,imax1,j,jmax1,k,kmax1,ipde,n,nh,ibc1,iblock
      integer(kind=cosa_int) :: i1,i2,i3,ijkmax(3),idir,istrt(3),iend(3),bctyp, &
        inrout,ibcpt,ibcpt2,ibcn,ibcn2,ibcm,ibc,jbc,kbc,ibc2,jbc2, &
        kbc2,in,jn,kn,in2,jn2,kn2,im,jm,km,ic1,ic2,ic3,sysize,ioff1, &
        ioff2,joff1,joff2,koff1,koff2,ibc_s,jbc_s,kbc_s,ibc_e,jbc_e, &
        kbc_e
      integer(kind=cosa_int) :: fid(0:2*nharms)
      integer(kind=cosa_int) :: bctopo(10,nbcs)
      real(kind=cosa_real) :: sgnm,nx,ny,nz,kx,ky,kz,rhow,tw,muw,uw,u1,u2,uwr,u1r,u2r, &
           vw,v1,v2,vwr,v1r,v2r,ww,w1,w2,wwr,w1r,w2r,ueta,veta,weta, &
           xix,xiy,xiz,etax,etay,etaz,zetax,zetay,zetaz,dudeta,dudzeta, &
           dudxi,dudx,dudy,dudz,dvdx,dvdy,dwdxi,dwdeta,dwdzeta,dvdeta, &
           dvdz,dwdx,dwdy,dwdz,divv,txx,txy,txz,tyy,tyz,tzz,mach,tauwx, &
           tauwy,tauwz,tauwpx,tauwpy,tauwpz,dn1,tauw,utau,dvdxi,dvdzeta
      real(kind=cosa_real) :: &
           q   (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           x   ( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           y   ( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           z   ( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           xdot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           ydot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           zdot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           si      (lmet, 0:imax+1, 0:jmax+1, 0:kmax+1,0:2*nharms*hbmove), &
           sj      (lmet, 0:imax+1, 0:jmax+1, 0:kmax+1,0:2*nharms*hbmove), &
           sk      (lmet, 0:imax+1, 0:jmax+1, 0:kmax+1,0:2*nharms*hbmove), &
           xideri  (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           xiderj  (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           xiderk  (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           etaderi (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           etaderj (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           etaderk (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           zetaderi(3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           zetaderj(3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           zetaderk(3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           dist( 0:imax  , 0:jmax  , 0:kmax)
      real(kind=cosa_real),allocatable :: cp(:,:,:),cf(:,:,:),yp(:,:,:)
      real(kind=cosa_real),allocatable :: zeroarray(:)
      integer(kind=cosa_int) :: ierr
      integer(kind=cosa_int) :: tempindex,datasize
      integer(kind=cosa_int) :: surfacestartindexes(3,nsurfblocks)
      integer(kind=cosa_int) :: surfaceendindexes(3,nsurfblocks)
      integer(kind=cosa_int) :: nsurfblocks, surfacenumber
      integer(kind=cosa_int) :: initialsurfacenumber
      character(len=5) :: blocknumname
      integer(kind=cosa_int) :: nfconns, fnmode, shrconn, isblock, valuelocation, isdouble
      integer(kind=cosa_int) :: tnfnodes, ncbfaces, tnbconns
      integer(kind=cosa_int) :: zonetype, strandid, parentzone
      integer(kind=cosa_int) :: teczne142, tecznemap142, tecdat142, tecfil142
      integer(kind=cosa_int) :: ibcmax, jbcmax, kbcmax
      integer(kind=cosa_int) :: imaxmax, jmaxmax, kmaxmax
      integer(kind=cosa_int) :: mpiid
      integer(kind=cosa_int) :: Null(*)
      POINTER   (NullPtr,Null)

      
#ifdef MPI

      call getmpiid(mpiid)    

! Specify that we are using an Ordered zone type
      zonetype = 0

! Zones don't have parents
      parentzone = 0

! We are not part of a strand
      strandid = 0
      nfconns = 0
      fnmode = 0
      shrconn = 0
      imaxmax = 0
      jmaxmax = 0
      kmaxmax = 0

! This specifies if we are writing in block or point format
! 1 is block format.  Binary tecplot files must be block format, 
! so this must be 1.
      isblock = 1
      valuelocation = 0
      isdouble = 1
      tnfnodes = 0
      ncbfaces = 0
      tnbconns = 0

      ijkmax(1) = imax
      ijkmax(2) = jmax
      ijkmax(3) = kmax

      imax1 = imax+1
      jmax1 = jmax+1
      kmax1 = kmax+1

      initialsurfacenumber = surfacenumber

      do n = 0,2*nharms
        nh = n*hbmove
           
!     This line is included because the n loop iterates through files 
!     so we need to ensure that each time round the loop we start and 
!     the correct place in the associated file.
        surfacenumber = initialsurfacenumber

        do ibc1=1,nbcs

           if (any(bctopo(1,ibc1).eq. &
                [1500,1501,1502,1503,1400,1401,1402,1403])) then
              
              ibc_s = surfacestartindexes(1,surfacenumber)
              jbc_s = surfacestartindexes(2,surfacenumber)
              kbc_s = surfacestartindexes(3,surfacenumber)
              
              ibc_e = surfaceendindexes(1,surfacenumber)
              jbc_e = surfaceendindexes(2,surfacenumber)
              kbc_e = surfaceendindexes(3,surfacenumber)
                                          
!     store boundary topology in mnemonic names
              bctyp    = bctopo(1,ibc1)
              idir     = bctopo(2,ibc1)
              inrout   = bctopo(3,ibc1)
              istrt(1) = bctopo(4,ibc1)
              iend(1)  = bctopo(5,ibc1)
              istrt(2) = bctopo(6,ibc1)
              iend(2)  = bctopo(7,ibc1)
              istrt(3) = bctopo(8,ibc1)
              iend(3)  = bctopo(9,ibc1)
              
!             set needed variables depending on whether the boundary is
!             the inner boundary (inrout = 1) or
!             the outer boundary (inrout > 1)
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

              datasize = 6
              if((ibc_e-ibc_s)+1 .ne. 0) datasize = datasize*((ibc_e-ibc_s)+1)
              if((jbc_e-jbc_s)+1 .ne. 0) datasize = datasize*((jbc_e-jbc_s)+1)
              if((kbc_e-kbc_s)+1 .ne. 0) datasize = datasize*((kbc_e-kbc_s)+1)

              allocate(cp(ibc_s:ibc_e,jbc_s:jbc_e,kbc_s:kbc_e))
              if (viscous) then
                 allocate(cf(ibc_s:ibc_e,jbc_s:jbc_e,kbc_s:kbc_e), &
                      yp(ibc_s:ibc_e,jbc_s:jbc_e,kbc_s:kbc_e))
              end if      
 
              do i3 = istrt(ic3),iend(ic3)
                 do i2 = istrt(ic2),iend(ic2)
                    
                    ibc  = ibcpt *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
                    jbc  = ibcpt *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
                    kbc  = ibcpt *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)
                    
                    im   = ibcm  *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
                    jm   = ibcm  *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
                    km   = ibcm  *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)                   

                    if (idir.eq.1) then
                       nx    = si(1,im,jm,km,nh)
                       ny    = si(2,im,jm,km,nh)
                       nz    = si(3,im,jm,km,nh)
                    else if (idir.eq.2) then
                       nx    = sj(1,im,jm,km,nh)
                       ny    = sj(2,im,jm,km,nh)
                       nz    = sj(3,im,jm,km,nh)
                    else if (idir.eq.3) then
                       nx    = sk(1,im,jm,km,nh)
                       ny    = sk(2,im,jm,km,nh)
                       nz    = sk(3,im,jm,km,nh)
                    end if
                    
                    kx   = (-1)**(1+1/inrout) * nx
                    ky   = (-1)**(1+1/inrout) * ny
                    kz   = (-1)**(1+1/inrout) * nz
                    
!     AJ Removed the machfs division here because if we are not viscous 
!     AJ machfs is zero.
                    cp(ibc,jbc,kbc) = 2*(q(ibc,jbc,kbc,5,n)-1/gamma)                                        
                    
                    if (viscous) then

                       cp(ibc,jbc,kbc) = cp(ibc,jbc,kbc)/ &
                            machfs**2
                                              
                       in   = ibcn  *krd(ic1,1) +i2*krd(ic2,1) +i3*krd(ic3,1)
                       jn   = ibcn  *krd(ic1,2) +i2*krd(ic2,2) +i3*krd(ic3,2)
                       kn   = ibcn  *krd(ic1,3) +i2*krd(ic2,3) +i3*krd(ic3,3)
                       
                       in2  = ibcn2 *krd(ic1,1) +i2*krd(ic2,1) +i3*krd(ic3,1)
                       jn2  = ibcn2 *krd(ic1,2) +i2*krd(ic2,2) +i3*krd(ic3,2)
                       kn2  = ibcn2 *krd(ic1,3) +i2*krd(ic2,3) +i3*krd(ic3,3)
                       
                       rhow = q(ibc,jbc,kbc,1,n)
                       tw   = q(ibc,jbc,kbc,5,n) *gamma / q(ibc,jbc,kbc,1,n)
                       muw  = (stemp+1)/(tw+stemp) * (tw**1.5d0)
                       
!------------------MSC, 11/03/2023: dist not computed for laminar now
                       if (kom.or.kom_bsl.or.kom_sst) then
                          dn1  = dist(in ,jn ,kn )
                       else
                          dn1  = 0
                       end if

                       if (idir.eq.1) then
                          
                          dudxi = (-1)**(1+1/inrout) * &
                               (-8*q(ibc,jbc,kbc,2,n) + &
                               9*q(ibc+ioff1,jbc,kbc,2,n) - &
                               q(ibc+ioff2,jbc,kbc,2,n) ) / 3
                          dvdxi = (-1)**(1+1/inrout) * &
                               (-8*q(ibc,jbc,kbc,3,n) + &
                               9*q(ibc+ioff1,jbc,kbc,3,n) - &
                               q(ibc+ioff2,jbc,kbc,3,n) ) / 3
                          dwdxi = (-1)**(1+1/inrout) * &
                               (-8*q(ibc,jbc,kbc,4,n) + &
                               9*q(ibc+ioff1,jbc,kbc,4,n) - &
                               q(ibc+ioff2,jbc,kbc,4,n) ) / 3
                          
                          if (i2.eq.istrt(ic2)) then
                             u2 = (  q(ibc,jbc,kbc,2,n) + q(ibc,jbc+1,kbc,2,n)) / 2
                             u1 = (3*q(ibc,jbc,kbc,2,n) - q(ibc,jbc+1,kbc,2,n)) / 2
                             v2 = (  q(ibc,jbc,kbc,3,n) + q(ibc,jbc+1,kbc,3,n)) / 2
                             v1 = (3*q(ibc,jbc,kbc,3,n) - q(ibc,jbc+1,kbc,3,n)) / 2
                             w2 = (  q(ibc,jbc,kbc,4,n) + q(ibc,jbc+1,kbc,4,n)) / 2
                             w1 = (3*q(ibc,jbc,kbc,4,n) - q(ibc,jbc+1,kbc,4,n)) / 2
                          else if (i2.eq.iend(ic2)) then
                             u2 = (3*q(ibc,jbc,kbc,2,n) - q(ibc,jbc-1,kbc,2,n)) / 2
                             u1 = (  q(ibc,jbc,kbc,2,n) + q(ibc,jbc-1,kbc,2,n)) / 2
                             v2 = (3*q(ibc,jbc,kbc,3,n) - q(ibc,jbc-1,kbc,3,n)) / 2
                             v1 = (  q(ibc,jbc,kbc,3,n) + q(ibc,jbc-1,kbc,3,n)) / 2
                             w2 = (3*q(ibc,jbc,kbc,4,n) - q(ibc,jbc-1,kbc,4,n)) / 2
                             w1 = (  q(ibc,jbc,kbc,4,n) + q(ibc,jbc-1,kbc,4,n)) / 2
                          else
                             u2 = (q(ibc,jbc  ,kbc,2,n) + q(ibc,jbc+1,kbc,2,n)) / 2
                             u1 = (q(ibc,jbc-1,kbc,2,n) + q(ibc,jbc  ,kbc,2,n)) / 2
                             v2 = (q(ibc,jbc  ,kbc,3,n) + q(ibc,jbc+1,kbc,3,n)) / 2
                             v1 = (q(ibc,jbc-1,kbc,3,n) + q(ibc,jbc  ,kbc,3,n)) / 2
                             w2 = (q(ibc,jbc  ,kbc,4,n) + q(ibc,jbc+1,kbc,4,n)) / 2
                             w1 = (q(ibc,jbc-1,kbc,4,n) + q(ibc,jbc  ,kbc,4,n)) / 2
                          end if
                          dudeta = u2 - u1
                          dvdeta = v2 - v1
                          dwdeta = w2 - w1
                          
                          if (i3.eq.istrt(ic3)) then
                             u2 = (  q(ibc,jbc,kbc,2,n) + q(ibc,jbc,kbc+1,2,n)) / 2
                             u1 = (3*q(ibc,jbc,kbc,2,n) - q(ibc,jbc,kbc+1,2,n)) / 2
                             v2 = (  q(ibc,jbc,kbc,3,n) + q(ibc,jbc,kbc+1,3,n)) / 2
                             v1 = (3*q(ibc,jbc,kbc,3,n) - q(ibc,jbc,kbc+1,3,n)) / 2
                             w2 = (  q(ibc,jbc,kbc,4,n) + q(ibc,jbc,kbc+1,4,n)) / 2
                             w1 = (3*q(ibc,jbc,kbc,4,n) - q(ibc,jbc,kbc+1,4,n)) / 2
                          else if (i3.eq.iend(ic3)) then
                             u2 = (3*q(ibc,jbc,kbc,2,n) - q(ibc,jbc,kbc-1,2,n)) / 2
                             u1 = (  q(ibc,jbc,kbc,2,n) + q(ibc,jbc,kbc-1,2,n)) / 2
                             v2 = (3*q(ibc,jbc,kbc,3,n) - q(ibc,jbc,kbc-1,3,n)) / 2
                             v1 = (  q(ibc,jbc,kbc,3,n) + q(ibc,jbc,kbc-1,3,n)) / 2
                             w2 = (3*q(ibc,jbc,kbc,4,n) - q(ibc,jbc,kbc-1,4,n)) / 2
                             w1 = (  q(ibc,jbc,kbc,4,n) + q(ibc,jbc,kbc-1,4,n)) / 2
                          else
                             u2 = (q(ibc,jbc,kbc  ,2,n) + q(ibc,jbc,kbc+1,2,n)) / 2
                             u1 = (q(ibc,jbc,kbc-1,2,n) + q(ibc,jbc,kbc  ,2,n)) / 2
                             v2 = (q(ibc,jbc,kbc  ,3,n) + q(ibc,jbc,kbc+1,3,n)) / 2
                             v1 = (q(ibc,jbc,kbc-1,3,n) + q(ibc,jbc,kbc  ,3,n)) / 2
                             w2 = (q(ibc,jbc,kbc  ,4,n) + q(ibc,jbc,kbc+1,4,n)) / 2
                             w1 = (q(ibc,jbc,kbc-1,4,n) + q(ibc,jbc,kbc  ,4,n)) / 2
                          end if
                          dudzeta = u2 - u1
                          dvdzeta = v2 - v1
                          dwdzeta = w2 - w1
                          
                          xix   = xideri  (1,im,jm,km,n)
                          xiy   = xideri  (2,im,jm,km,n)
                          xiz   = xideri  (3,im,jm,km,n)
                          etax  = etaderi (1,im,jm,km,n)
                          etay  = etaderi (2,im,jm,km,n)
                          etaz  = etaderi (3,im,jm,km,n)
                          zetax = zetaderi(1,im,jm,km,n)
                          zetay = zetaderi(2,im,jm,km,n)
                          zetaz = zetaderi(3,im,jm,km,n)
                          
                       else if (idir.eq.2) then
                          
                          if (i3.eq.istrt(ic3)) then
                             u2 = (  q(ibc+1,jbc,kbc,2,n) + q(ibc  ,jbc,kbc,2,n))/2
                             u1 = (3*q(ibc  ,jbc,kbc,2,n) - q(ibc+1,jbc,kbc,2,n))/2
                             v2 = (  q(ibc+1,jbc,kbc,3,n) + q(ibc  ,jbc,kbc,3,n))/2
                             v1 = (3*q(ibc  ,jbc,kbc,3,n) - q(ibc+1,jbc,kbc,3,n))/2
                             w2 = (  q(ibc+1,jbc,kbc,4,n) + q(ibc  ,jbc,kbc,4,n))/2
                             w1 = (3*q(ibc  ,jbc,kbc,4,n) - q(ibc+1,jbc,kbc,4,n))/2
                          else if (i3.eq.iend(ic3)) then
                             u2 = (3*q(ibc,jbc,kbc,2,n) - q(ibc-1,jbc,kbc,2,n)) / 2
                             u1 = (  q(ibc,jbc,kbc,2,n) + q(ibc-1,jbc,kbc,2,n)) / 2
                             v2 = (3*q(ibc,jbc,kbc,3,n) - q(ibc-1,jbc,kbc,3,n)) / 2
                             v1 = (  q(ibc,jbc,kbc,3,n) + q(ibc-1,jbc,kbc,3,n)) / 2
                             w2 = (3*q(ibc,jbc,kbc,4,n) - q(ibc-1,jbc,kbc,4,n)) / 2
                             w1 = (  q(ibc,jbc,kbc,4,n) + q(ibc-1,jbc,kbc,4,n)) / 2
                          else
                             u2 = (q(ibc+1,jbc,kbc,2,n) + q(ibc  ,jbc,kbc,2,n)) / 2
                             u1 = (q(ibc  ,jbc,kbc,2,n) + q(ibc-1,jbc,kbc,2,n)) / 2
                             v2 = (q(ibc+1,jbc,kbc,3,n) + q(ibc  ,jbc,kbc,3,n)) / 2
                             v1 = (q(ibc  ,jbc,kbc,3,n) + q(ibc-1,jbc,kbc,3,n)) / 2
                             w2 = (q(ibc+1,jbc,kbc,4,n) + q(ibc  ,jbc,kbc,4,n)) / 2
                             w1 = (q(ibc  ,jbc,kbc,4,n) + q(ibc-1,jbc,kbc,4,n)) / 2
                          end if
                          dudxi = u2 - u1
                          dvdxi = v2 - v1
                          dwdxi = w2 - w1
                          
                          dudeta = (-1)**(1+1/inrout) * &
                               (-8*q(ibc,jbc,kbc,2,n) + &
                               9*q(ibc,jbc+joff1,kbc,2,n) - &
                               q(ibc,jbc+joff2,kbc,2,n) ) / 3
                          dvdeta = (-1)**(1+1/inrout) * &
                               (-8*q(ibc,jbc,kbc,3,n) + &
                               9*q(ibc,jbc+joff1,kbc,3,n) - &
                               q(ibc,jbc+joff2,kbc,3,n) ) / 3
                          dwdeta = (-1)**(1+1/inrout) * &
                               (-8*q(ibc,jbc,kbc,4,n) + &
                               9*q(ibc,jbc+joff1,kbc,4,n) - &
                               q(ibc,jbc+joff2,kbc,4,n) ) / 3
                          
                          if (i2.eq.istrt(ic2)) then
                             u2 = (  q(ibc,jbc,kbc+1,2,n) + q(ibc,jbc,kbc  ,2,n))/2
                             u1 = (3*q(ibc,jbc,kbc  ,2,n) - q(ibc,jbc,kbc+1,2,n))/2
                             v2 = (  q(ibc,jbc,kbc+1,3,n) + q(ibc,jbc,kbc  ,3,n))/2
                             v1 = (3*q(ibc,jbc,kbc  ,3,n) - q(ibc,jbc,kbc+1,3,n))/2
                             w2 = (  q(ibc,jbc,kbc+1,4,n) + q(ibc,jbc,kbc  ,4,n))/2
                             w1 = (3*q(ibc,jbc,kbc  ,4,n) - q(ibc,jbc,kbc+1,4,n))/2
                          else if (i2.eq.iend(ic2)) then
                             u2 = (3*q(ibc,jbc,kbc,2,n) - q(ibc,jbc,kbc-1,2,n)) / 2
                             u1 = (  q(ibc,jbc,kbc,2,n) + q(ibc,jbc,kbc-1,2,n)) / 2
                             v2 = (3*q(ibc,jbc,kbc,3,n) - q(ibc,jbc,kbc-1,3,n)) / 2
                             v1 = (  q(ibc,jbc,kbc,3,n) + q(ibc,jbc,kbc-1,3,n)) / 2
                             w2 = (3*q(ibc,jbc,kbc,4,n) - q(ibc,jbc,kbc-1,4,n)) / 2
                             w1 = (  q(ibc,jbc,kbc,4,n) + q(ibc,jbc,kbc-1,4,n)) / 2
                          else
                             u2 = (q(ibc,jbc,kbc+1,2,n) + q(ibc,jbc,kbc  ,2,n)) / 2
                             u1 = (q(ibc,jbc,kbc  ,2,n) + q(ibc,jbc,kbc-1,2,n)) / 2
                             v2 = (q(ibc,jbc,kbc+1,3,n) + q(ibc,jbc,kbc  ,3,n)) / 2
                             v1 = (q(ibc,jbc,kbc  ,3,n) + q(ibc,jbc,kbc-1,3,n)) / 2
                             w2 = (q(ibc,jbc,kbc+1,4,n) + q(ibc,jbc,kbc  ,4,n)) / 2
                             w1 = (q(ibc,jbc,kbc  ,4,n) + q(ibc,jbc,kbc-1,4,n)) / 2
                          end if
                          dudzeta = u2 - u1
                          dvdzeta = v2 - v1
                          dwdzeta = w2 - w1
                          
                          xix   = xiderj  (1,im,jm,km,n)
                          xiy   = xiderj  (2,im,jm,km,n)
                          xiz   = xiderj  (3,im,jm,km,n)
                          etax  = etaderj (1,im,jm,km,n)
                          etay  = etaderj (2,im,jm,km,n)
                          etaz  = etaderj (3,im,jm,km,n)
                          zetax = zetaderj(1,im,jm,km,n)
                          zetay = zetaderj(2,im,jm,km,n)
                          zetaz = zetaderj(3,im,jm,km,n)
                          
                       else if (idir.eq.3) then
                          if (i2.eq.istrt(ic2)) then
                             u2 = (  q(ibc,jbc,kbc,2,n) + q(ibc+1,jbc,kbc,2,n)) / 2
                             u1 = (3*q(ibc,jbc,kbc,2,n) - q(ibc+1,jbc,kbc,2,n)) / 2
                             v2 = (  q(ibc,jbc,kbc,3,n) + q(ibc+1,jbc,kbc,3,n)) / 2
                             v1 = (3*q(ibc,jbc,kbc,3,n) - q(ibc+1,jbc,kbc,3,n)) / 2
                             w2 = (  q(ibc,jbc,kbc,4,n) + q(ibc+1,jbc,kbc,4,n)) / 2
                             w1 = (3*q(ibc,jbc,kbc,4,n) - q(ibc+1,jbc,kbc,4,n)) / 2
                          else if (i2.eq.iend(ic2)) then
                             u2 = (3*q(ibc,jbc,kbc,2,n) - q(ibc-1,jbc,kbc,2,n)) / 2
                             u1 = (  q(ibc,jbc,kbc,2,n) + q(ibc-1,jbc,kbc,2,n)) / 2
                             v2 = (3*q(ibc,jbc,kbc,3,n) - q(ibc-1,jbc,kbc,3,n)) / 2
                             v1 = (  q(ibc,jbc,kbc,3,n) + q(ibc-1,jbc,kbc,3,n)) / 2
                             w2 = (3*q(ibc,jbc,kbc,4,n) - q(ibc-1,jbc,kbc,4,n)) / 2
                             w1 = (  q(ibc,jbc,kbc,4,n) + q(ibc-1,jbc,kbc,4,n)) / 2
                          else
                             u2 = (q(ibc  ,jbc,kbc,2,n) + q(ibc+1,jbc,kbc,2,n)) / 2
                             u1 = (q(ibc-1,jbc,kbc,2,n) + q(ibc  ,jbc,kbc,2,n)) / 2
                             v2 = (q(ibc  ,jbc,kbc,3,n) + q(ibc+1,jbc,kbc,3,n)) / 2
                             v1 = (q(ibc-1,jbc,kbc,3,n) + q(ibc  ,jbc,kbc,3,n)) / 2
                             w2 = (q(ibc  ,jbc,kbc,4,n) + q(ibc+1,jbc,kbc,4,n)) / 2
                             w1 = (q(ibc-1,jbc,kbc,4,n) + q(ibc  ,jbc,kbc,4,n)) / 2
                          end if
                          dudxi = u2 - u1
                          dvdxi = v2 - v1
                          dwdxi = w2 - w1
                          
                          if (i3.eq.istrt(ic3)) then
                             u2 = (  q(ibc,jbc,kbc,2,n) + q(ibc,jbc+1,kbc,2,n)) / 2
                             u1 = (3*q(ibc,jbc,kbc,2,n) - q(ibc,jbc+1,kbc,2,n)) / 2
                             v2 = (  q(ibc,jbc,kbc,3,n) + q(ibc,jbc+1,kbc,3,n)) / 2
                             v1 = (3*q(ibc,jbc,kbc,3,n) - q(ibc,jbc+1,kbc,3,n)) / 2
                             w2 = (  q(ibc,jbc,kbc,4,n) + q(ibc,jbc+1,kbc,4,n)) / 2
                             w1 = (3*q(ibc,jbc,kbc,4,n) - q(ibc,jbc+1,kbc,4,n)) / 2
                          else if (i3.eq.iend(ic3)) then
                             u2 = (3*q(ibc,jbc,kbc,2,n) - q(ibc,jbc-1,kbc,2,n)) / 2
                             u1 = (  q(ibc,jbc,kbc,2,n) + q(ibc,jbc-1,kbc,2,n)) / 2
                             v2 = (3*q(ibc,jbc,kbc,3,n) - q(ibc,jbc-1,kbc,3,n)) / 2
                             v1 = (  q(ibc,jbc,kbc,3,n) + q(ibc,jbc-1,kbc,3,n)) / 2
                             w2 = (3*q(ibc,jbc,kbc,4,n) - q(ibc,jbc-1,kbc,4,n)) / 2
                             w1 = (  q(ibc,jbc,kbc,4,n) + q(ibc,jbc-1,kbc,4,n)) / 2
                          else
                             u2 = (q(ibc,jbc  ,kbc,2,n) + q(ibc,jbc+1,kbc,2,n)) / 2
                             u1 = (q(ibc,jbc-1,kbc,2,n) + q(ibc,jbc  ,kbc,2,n)) / 2
                             v2 = (q(ibc,jbc  ,kbc,3,n) + q(ibc,jbc+1,kbc,3,n)) / 2
                             v1 = (q(ibc,jbc-1,kbc,3,n) + q(ibc,jbc  ,kbc,3,n)) / 2
                             w2 = (q(ibc,jbc  ,kbc,4,n) + q(ibc,jbc+1,kbc,4,n)) / 2
                             w1 = (q(ibc,jbc-1,kbc,4,n) + q(ibc,jbc  ,kbc,4,n)) / 2
                          end if
                          dudeta = u2 - u1
                          dvdeta = v2 - v1
                          dwdeta = w2 - w1
                          
                          dudzeta = (-1)**(1+1/inrout) * &
                               (-8*q(ibc,jbc,kbc,2,n) + &
                               9*q(ibc,jbc,kbc+koff1,2,n) - &
                               q(ibc,jbc,kbc+koff2,2,n) ) / 3
                          dvdzeta = (-1)**(1+1/inrout) * &
                               (-8*q(ibc,jbc,kbc,3,n) + &
                               9*q(ibc,jbc,kbc+koff1,3,n) - &
                               q(ibc,jbc,kbc+koff2,3,n) ) / 3
                          dwdzeta = (-1)**(1+1/inrout) * &
                               (-8*q(ibc,jbc,kbc,4,n) + &
                               9*q(ibc,jbc,kbc+koff1,4,n) - &
                               q(ibc,jbc,kbc+koff2,4,n) ) / 3
                          
                          xix   = xiderk  (1,im,jm,km,n)
                          xiy   = xiderk  (2,im,jm,km,n)
                          xiz   = xiderk  (3,im,jm,km,n)
                          etax  = etaderk (1,im,jm,km,n)
                          etay  = etaderk (2,im,jm,km,n)
                          etaz  = etaderk (3,im,jm,km,n)
                          zetax = zetaderk(1,im,jm,km,n)
                          zetay = zetaderk(2,im,jm,km,n)
                          zetaz = zetaderk(3,im,jm,km,n)
                          
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
                       
                       cf(ibc,jbc,kbc) = 2*tauw / (reyno*machfs)
                       yp(ibc,jbc,kbc) = dsqrt(reyno/machfs)* &
                            rhow*dn1*utau/muw
                       
                    end if
                    
                 end do
              end do

              if(tecfil142(n+1) .ne. 0) then
                 write(*,*) 'error calling tecfil142'
                 stop
              end if
              
              write (blocknumname, "(I5)") iblock
              
              ibcmax = ibc_e-ibc_s+1
              jbcmax = jbc_e-jbc_s+1
              kbcmax = kbc_e-kbc_s+1
              
              if(teczne142('block'//trim(blocknumname)//char(0),zonetype, ibcmax, jbcmax, kbcmax, &
                   imaxmax, jmaxmax, kmaxmax, simtime, strandid, parentzone, isblock, nfconns, & 
                   fnmode, tnfnodes, ncbfaces, tnbconns, Null, Null, Null, shrconn) .ne. 0) then
                 write(*,*) 'error setting up zone'
                 stop
              end if
              
              if(tecznemap142(1, mpiid) .ne. 0) then
                 write(*,*) 'error setting up parallel mapping for zone'
                 call abortmpi()
              end if              

              if(tecdat142(ibcmax*jbcmax*kbcmax,x(ibc_s,jbc_s,kbc_s,n),isdouble) .ne. 0) then
                 write(*,*) 'error writing block data'
              end if
              if(tecdat142(ibcmax*jbcmax*kbcmax,y(ibc_s,jbc_s,kbc_s,n),isdouble) .ne. 0) then
                 write(*,*) 'error writing block data'
              end if
              if(tecdat142(ibcmax*jbcmax*kbcmax,z(ibc_s,jbc_s,kbc_s,n),isdouble) .ne. 0) then
                 write(*,*) 'error writing block data'
              end if
              if(tecdat142(ibcmax*jbcmax*kbcmax,cp(ibc_s,jbc_s,kbc_s),isdouble) .ne. 0) then
                 write(*,*) 'error writing block data'
              end if
              
              if (viscous) then
                 if(tecdat142(ibcmax*jbcmax*kbcmax,cf(ibc_s,jbc_s,kbc_s),isdouble) .ne. 0) then
                    write(*,*) 'error writing block data'
                 end if
                 if(tecdat142(ibcmax*jbcmax*kbcmax,yp(ibc_s,jbc_s,kbc_s),isdouble) .ne. 0) then
                    write(*,*) 'error writing block data'
                 end if
              else
                 allocate(zeroarray(ibcmax*jbcmax*kbcmax))
                 zeroarray = 0.d0
                 if(tecdat142(ibcmax*jbcmax*kbcmax,zeroarray,isdouble) .ne. 0) then
                    write(*,*) 'error writing block data'
                 end if
                 if(tecdat142(ibcmax*jbcmax*kbcmax,zeroarray,isdouble) .ne. 0) then
                    write(*,*) 'error writing block data'
                 end if
                 deallocate(zeroarray)
              end if
              
              if(allocated(cp)) deallocate(cp)
              if (viscous) then
                 if(allocated(cf)) deallocate(cf,yp)
              end if

              surfacenumber = surfacenumber + 1
              
           end if
        end do

      end do      

#endif

      return
      end


! AJ TODO Merge the tec and cgns versions of this to save code duplication
!-----------------------------------------------------------------------
      subroutine write_parallel_wr_cgns_surf_b_M2(fid,imax, &
                 jmax,kmax,npde,nharms,iblock,nl,lmet,q,x,y,z,xdot,ydot, &
                 zdot,si,sj,sk,xideri,xiderj,xiderk,etaderi,etaderj, &
                 etaderk,zetaderi,zetaderj,zetaderk,dist,bctopo,nbcs, &
                 basenum)
!     This routine writes parts of the surface output files using CGNS
!     Mainly this routine wraps the bespoke code which writes the
!     data (in parallel_wr_surf_cgns_b_M2).
!-----------------------------------------------------------------------

      use parallelutils, only: nsurfblocks, surfblockstarts, &
           surfblockends, surfblockindex, surfbcindex
      use common_variables
      use cosa_precision

      implicit none

#if MPI 
      include 'mpif.h'
#endif

      integer(kind=cosa_int) :: i,imax,jmax,kmax,npde,nharms,iblock,nl,lmet,nbcs
      integer(kind=cosa_int) :: fid(0:2*nharms)
      integer(kind=cosa_int) :: bctopo(10,nbcs)
      real(kind=cosa_real) :: &
           q   (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           x   ( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           y   ( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           z   ( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           xdot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           ydot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           zdot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           si      (lmet, 0:imax+1, 0:jmax+1, 0:kmax+1,0:2*nharms*hbmove), &
           sj      (lmet, 0:imax+1, 0:jmax+1, 0:kmax+1,0:2*nharms*hbmove), &
           sk      (lmet, 0:imax+1, 0:jmax+1, 0:kmax+1,0:2*nharms*hbmove), &
           xideri  (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           xiderj  (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           xiderk  (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           etaderi (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           etaderj (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           etaderk (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           zetaderi(3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           zetaderj(3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           zetaderk(3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           dist( 0:imax  , 0:jmax  , 0:kmax)
      integer(kind=cosa_int) :: integersize,doublesize,charactersize
      integer(kind=cosa_int) :: surfacenumber, basenum
      integer*8 :: disp

!AJ We do not need the calculated offset here, but we do need the surfacenumber
!AJ variable. We could refactor and create a new subroutine that just produces that 
!AJ variable and not the displacement as well, but for the amount of compute time 
!AJ required it seems unnecessary
      call surftecinitialoffset(disp,nl,iblock,nsurfblocks, &
           surfblockstarts,surfblockends,surfblockindex,surfbcindex, &
           surfacenumber,integersize,doublesize,charactersize)

      call parallel_wr_cgns_surf_b_M2(fid,q,x,y,z,xdot,ydot,zdot,si,sj, &
           sk,xideri,xiderj,xiderk,etaderi,etaderj,etaderk,zetaderi, &
           zetaderj,zetaderk,dist,bctopo,imax,jmax,kmax,npde, &
           nharms,nl,nbcs,lmet,integersize,doublesize, &
           charactersize,surfblockstarts,surfblockends, &
           nsurfblocks,surfacenumber,basenum)

      return 
      end

!-----------------------------------------------------------------------
      subroutine parallel_wr_cgns_surf_b_M2(fid,q,x,y,z,xdot,ydot,zdot, &
                 si,sj,sk,xideri,xiderj,xiderk,etaderi,etaderj,etaderk, &
                 zetaderi,zetaderj,zetaderk,dist,bctopo,imax,jmax, &
                 kmax,npde,nharms,nl,nbcs,lmet,integersize, &
                 doublesize,charactersize,surfacestartindexes, &
                 surfaceendindexes,nsurfblocks,surfacenumber,basenum)
!     This routine writes parts of the surface output files using MPI-I/O.  
!-----------------------------------------------------------------------

      use cgns
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) :: imax,jmax,kmax,npde,nharms,nl,lmet,nbcs
      integer(kind=cosa_int) :: i,imax1,j,jmax1,k,kmax1,ipde,n,nh,ibc1
      integer(kind=cosa_int) :: i1,i2,i3,ijkmax(3),idir,istrt(3),iend(3),bctyp, &
        inrout,ibcpt,ibcpt2,ibcn,ibcn2,ibcm,ibc,jbc,kbc,ibc2,jbc2, &
        kbc2,in,jn,kn,in2,jn2,kn2,im,jm,km,ic1,ic2,ic3,sysize,ioff1, &
        ioff2,joff1,joff2,koff1,koff2,ibc_s,jbc_s,kbc_s,ibc_e,jbc_e, &
        kbc_e
      integer(kind=cosa_int) :: fid(0:2*nharms)
      integer(kind=cosa_int) :: bctopo(10,nbcs)
      real(kind=cosa_real) :: sgnm,nx,ny,nz,kx,ky,kz,rhow,tw,muw,uw,u1,u2,uwr,u1r,u2r, &
           vw,v1,v2,vwr,v1r,v2r,ww,w1,w2,wwr,w1r,w2r,ueta,veta,weta, &
           xix,xiy,xiz,etax,etay,etaz,zetax,zetay,zetaz,dudeta,dudzeta, &
           dudxi,dudx,dudy,dudz,dvdx,dvdy,dwdxi,dwdeta,dwdzeta,dvdeta, &
           dvdz,dwdx,dwdy,dwdz,divv,txx,txy,txz,tyy,tyz,tzz,mach,tauwx, &
           tauwy,tauwz,tauwpx,tauwpy,tauwpz,dn1,tauw,utau,dvdxi,dvdzeta
      real(kind=cosa_real) :: &
           q   (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           x   ( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           y   ( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           z   ( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           xdot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           ydot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           zdot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           si      (lmet, 0:imax+1, 0:jmax+1, 0:kmax+1,0:2*nharms*hbmove), &
           sj      (lmet, 0:imax+1, 0:jmax+1, 0:kmax+1,0:2*nharms*hbmove), &
           sk      (lmet, 0:imax+1, 0:jmax+1, 0:kmax+1,0:2*nharms*hbmove), &
           xideri  (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           xiderj  (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           xiderk  (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           etaderi (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           etaderj (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           etaderk (3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           zetaderi(3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           zetaderj(3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           zetaderk(3   ,   imax  ,   jmax  ,   kmax  ,0:2*nharms*hbmove), &
           dist( 0:imax  , 0:jmax  , 0:kmax)
      real(kind=cosa_real),allocatable :: cp(:,:,:),cf(:,:,:),yp(:,:,:)
      integer(kind=cosa_int) :: ierr,integersize,doublesize,charactersize
      integer(kind=cosa_int) :: surfacestartindexes(3,nsurfblocks)
      integer(kind=cosa_int) :: surfaceendindexes(3,nsurfblocks)
      integer(kind=cosa_int) :: nsurfblocks, surfacenumber
      integer(kind=cosa_int) :: xnum, ynum, znum, flownum, solnum, basenum
      character*200 line1
      integer(cgsize_t) :: minrange(3), maxrange(3)
      integer(kind=cosa_int) :: initialsurfacenumber

      ijkmax(1) = imax
      ijkmax(2) = jmax
      ijkmax(3) = kmax

      imax1 = imax+1
      jmax1 = jmax+1
      kmax1 = kmax+1

      initialsurfacenumber = surfacenumber

      do n = 0,2*nharms
        nh = n*hbmove
           
!     Required to deal with resetting everything when iterating harmonics.
        surfacenumber = initialsurfacenumber

        do ibc1=1,nbcs

           if (any(bctopo(1,ibc1).eq. &
                [1500,1501,1502,1503,1400,1401,1402,1403])) then
              
              ibc_s = surfacestartindexes(1,surfacenumber)
              jbc_s = surfacestartindexes(2,surfacenumber)
              kbc_s = surfacestartindexes(3,surfacenumber)
              
              ibc_e = surfaceendindexes(1,surfacenumber)
              jbc_e = surfaceendindexes(2,surfacenumber)
              kbc_e = surfaceendindexes(3,surfacenumber)

              minrange(1) = 1
              minrange(2) = 1
              minrange(3) = 1

              maxrange(1) = (ibc_e - ibc_s)+1
              maxrange(2) = (jbc_e - jbc_s)+1
              maxrange(3) = (kbc_e - kbc_s)+1

!              write(*,*) 'final surf size',ibc1,surfacenumber,maxrange
              
!     store boundary topology in mnemonic names
              bctyp    = bctopo(1,ibc1)
              idir     = bctopo(2,ibc1)
              inrout   = bctopo(3,ibc1)
              istrt(1) = bctopo(4,ibc1)
              iend(1)  = bctopo(5,ibc1)
              istrt(2) = bctopo(6,ibc1)
              iend(2)  = bctopo(7,ibc1)
              istrt(3) = bctopo(8,ibc1)
              iend(3)  = bctopo(9,ibc1)
              
!             set needed variables depending on whether the boundary is
!             the inner boundary (inrout = 1) or
!             the outer boundary (inrout > 1)
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

              allocate(cp(ibc_s:ibc_e,jbc_s:jbc_e,kbc_s:kbc_e))
              allocate(cf(ibc_s:ibc_e,jbc_s:jbc_e,kbc_s:kbc_e), &
                   yp(ibc_s:ibc_e,jbc_s:jbc_e,kbc_s:kbc_e))

 
              do i3 = istrt(ic3),iend(ic3)
                 do i2 = istrt(ic2),iend(ic2)
                    
                    ibc  = ibcpt *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
                    jbc  = ibcpt *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
                    kbc  = ibcpt *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)
                    
                    im   = ibcm  *krd(ic1,1) + i2*krd(ic2,1) + i3*krd(ic3,1)
                    jm   = ibcm  *krd(ic1,2) + i2*krd(ic2,2) + i3*krd(ic3,2)
                    km   = ibcm  *krd(ic1,3) + i2*krd(ic2,3) + i3*krd(ic3,3)                   

                    if (idir.eq.1) then
                       nx    = si(1,im,jm,km,nh)
                       ny    = si(2,im,jm,km,nh)
                       nz    = si(3,im,jm,km,nh)
                    else if (idir.eq.2) then
                       nx    = sj(1,im,jm,km,nh)
                       ny    = sj(2,im,jm,km,nh)
                       nz    = sj(3,im,jm,km,nh)
                    else if (idir.eq.3) then
                       nx    = sk(1,im,jm,km,nh)
                       ny    = sk(2,im,jm,km,nh)
                       nz    = sk(3,im,jm,km,nh)
                    end if
                    
                    kx   = (-1)**(1+1/inrout) * nx
                    ky   = (-1)**(1+1/inrout) * ny
                    kz   = (-1)**(1+1/inrout) * nz
                    
!     AJ Removed the machfs division here because if we are not viscous 
!     AJ machfs is zero.
                    cp(ibc,jbc,kbc) = 2*(q(ibc,jbc,kbc,5,n)-1/gamma)                                        
                    
                    if (viscous) then

                       cp(ibc,jbc,kbc) = cp(ibc,jbc,kbc)/ &
                            machfs**2
                                              
                       in   = ibcn  *krd(ic1,1) +i2*krd(ic2,1) +i3*krd(ic3,1)
                       jn   = ibcn  *krd(ic1,2) +i2*krd(ic2,2) +i3*krd(ic3,2)
                       kn   = ibcn  *krd(ic1,3) +i2*krd(ic2,3) +i3*krd(ic3,3)
                       
                       in2  = ibcn2 *krd(ic1,1) +i2*krd(ic2,1) +i3*krd(ic3,1)
                       jn2  = ibcn2 *krd(ic1,2) +i2*krd(ic2,2) +i3*krd(ic3,2)
                       kn2  = ibcn2 *krd(ic1,3) +i2*krd(ic2,3) +i3*krd(ic3,3)
                       
                       rhow = q(ibc,jbc,kbc,1,n)
                       tw   = q(ibc,jbc,kbc,5,n) *gamma / q(ibc,jbc,kbc,1,n)
                       muw  = (stemp+1)/(tw+stemp) * (tw**1.5d0)
                       
!------------------MSC, 11/03/2023: dist not computed for laminar now
                       if (kom.or.kom_bsl.or.kom_sst) then
                          dn1  = dist(in ,jn ,kn )
                       else
                          dn1  = 0
                       end if

                       if (idir.eq.1) then
                          
                          dudxi = (-1)**(1+1/inrout) * &
                               (-8*q(ibc,jbc,kbc,2,n) + &
                               9*q(ibc+ioff1,jbc,kbc,2,n) - &
                               q(ibc+ioff2,jbc,kbc,2,n) ) / 3
                          dvdxi = (-1)**(1+1/inrout) * &
                               (-8*q(ibc,jbc,kbc,3,n) + &
                               9*q(ibc+ioff1,jbc,kbc,3,n) - &
                               q(ibc+ioff2,jbc,kbc,3,n) ) / 3
                          dwdxi = (-1)**(1+1/inrout) * &
                               (-8*q(ibc,jbc,kbc,4,n) + &
                               9*q(ibc+ioff1,jbc,kbc,4,n) - &
                               q(ibc+ioff2,jbc,kbc,4,n) ) / 3
                          
                          if (i2.eq.istrt(ic2)) then
                             u2 = (  q(ibc,jbc,kbc,2,n) + q(ibc,jbc+1,kbc,2,n)) / 2
                             u1 = (3*q(ibc,jbc,kbc,2,n) - q(ibc,jbc+1,kbc,2,n)) / 2
                             v2 = (  q(ibc,jbc,kbc,3,n) + q(ibc,jbc+1,kbc,3,n)) / 2
                             v1 = (3*q(ibc,jbc,kbc,3,n) - q(ibc,jbc+1,kbc,3,n)) / 2
                             w2 = (  q(ibc,jbc,kbc,4,n) + q(ibc,jbc+1,kbc,4,n)) / 2
                             w1 = (3*q(ibc,jbc,kbc,4,n) - q(ibc,jbc+1,kbc,4,n)) / 2
                          else if (i2.eq.iend(ic2)) then
                             u2 = (3*q(ibc,jbc,kbc,2,n) - q(ibc,jbc-1,kbc,2,n)) / 2
                             u1 = (  q(ibc,jbc,kbc,2,n) + q(ibc,jbc-1,kbc,2,n)) / 2
                             v2 = (3*q(ibc,jbc,kbc,3,n) - q(ibc,jbc-1,kbc,3,n)) / 2
                             v1 = (  q(ibc,jbc,kbc,3,n) + q(ibc,jbc-1,kbc,3,n)) / 2
                             w2 = (3*q(ibc,jbc,kbc,4,n) - q(ibc,jbc-1,kbc,4,n)) / 2
                             w1 = (  q(ibc,jbc,kbc,4,n) + q(ibc,jbc-1,kbc,4,n)) / 2
                          else
                             u2 = (q(ibc,jbc  ,kbc,2,n) + q(ibc,jbc+1,kbc,2,n)) / 2
                             u1 = (q(ibc,jbc-1,kbc,2,n) + q(ibc,jbc  ,kbc,2,n)) / 2
                             v2 = (q(ibc,jbc  ,kbc,3,n) + q(ibc,jbc+1,kbc,3,n)) / 2
                             v1 = (q(ibc,jbc-1,kbc,3,n) + q(ibc,jbc  ,kbc,3,n)) / 2
                             w2 = (q(ibc,jbc  ,kbc,4,n) + q(ibc,jbc+1,kbc,4,n)) / 2
                             w1 = (q(ibc,jbc-1,kbc,4,n) + q(ibc,jbc  ,kbc,4,n)) / 2
                          end if
                          dudeta = u2 - u1
                          dvdeta = v2 - v1
                          dwdeta = w2 - w1
                          
                          if (i3.eq.istrt(ic3)) then
                             u2 = (  q(ibc,jbc,kbc,2,n) + q(ibc,jbc,kbc+1,2,n)) / 2
                             u1 = (3*q(ibc,jbc,kbc,2,n) - q(ibc,jbc,kbc+1,2,n)) / 2
                             v2 = (  q(ibc,jbc,kbc,3,n) + q(ibc,jbc,kbc+1,3,n)) / 2
                             v1 = (3*q(ibc,jbc,kbc,3,n) - q(ibc,jbc,kbc+1,3,n)) / 2
                             w2 = (  q(ibc,jbc,kbc,4,n) + q(ibc,jbc,kbc+1,4,n)) / 2
                             w1 = (3*q(ibc,jbc,kbc,4,n) - q(ibc,jbc,kbc+1,4,n)) / 2
                          else if (i3.eq.iend(ic3)) then
                             u2 = (3*q(ibc,jbc,kbc,2,n) - q(ibc,jbc,kbc-1,2,n)) / 2
                             u1 = (  q(ibc,jbc,kbc,2,n) + q(ibc,jbc,kbc-1,2,n)) / 2
                             v2 = (3*q(ibc,jbc,kbc,3,n) - q(ibc,jbc,kbc-1,3,n)) / 2
                             v1 = (  q(ibc,jbc,kbc,3,n) + q(ibc,jbc,kbc-1,3,n)) / 2
                             w2 = (3*q(ibc,jbc,kbc,4,n) - q(ibc,jbc,kbc-1,4,n)) / 2
                             w1 = (  q(ibc,jbc,kbc,4,n) + q(ibc,jbc,kbc-1,4,n)) / 2
                          else
                             u2 = (q(ibc,jbc,kbc  ,2,n) + q(ibc,jbc,kbc+1,2,n)) / 2
                             u1 = (q(ibc,jbc,kbc-1,2,n) + q(ibc,jbc,kbc  ,2,n)) / 2
                             v2 = (q(ibc,jbc,kbc  ,3,n) + q(ibc,jbc,kbc+1,3,n)) / 2
                             v1 = (q(ibc,jbc,kbc-1,3,n) + q(ibc,jbc,kbc  ,3,n)) / 2
                             w2 = (q(ibc,jbc,kbc  ,4,n) + q(ibc,jbc,kbc+1,4,n)) / 2
                             w1 = (q(ibc,jbc,kbc-1,4,n) + q(ibc,jbc,kbc  ,4,n)) / 2
                          end if
                          dudzeta = u2 - u1
                          dvdzeta = v2 - v1
                          dwdzeta = w2 - w1
                          
                          xix   = xideri  (1,im,jm,km,n)
                          xiy   = xideri  (2,im,jm,km,n)
                          xiz   = xideri  (3,im,jm,km,n)
                          etax  = etaderi (1,im,jm,km,n)
                          etay  = etaderi (2,im,jm,km,n)
                          etaz  = etaderi (3,im,jm,km,n)
                          zetax = zetaderi(1,im,jm,km,n)
                          zetay = zetaderi(2,im,jm,km,n)
                          zetaz = zetaderi(3,im,jm,km,n)
                          
                       else if (idir.eq.2) then
                          
                          if (i3.eq.istrt(ic3)) then
                             u2 = (  q(ibc+1,jbc,kbc,2,n) + q(ibc  ,jbc,kbc,2,n))/2
                             u1 = (3*q(ibc  ,jbc,kbc,2,n) - q(ibc+1,jbc,kbc,2,n))/2
                             v2 = (  q(ibc+1,jbc,kbc,3,n) + q(ibc  ,jbc,kbc,3,n))/2
                             v1 = (3*q(ibc  ,jbc,kbc,3,n) - q(ibc+1,jbc,kbc,3,n))/2
                             w2 = (  q(ibc+1,jbc,kbc,4,n) + q(ibc  ,jbc,kbc,4,n))/2
                             w1 = (3*q(ibc  ,jbc,kbc,4,n) - q(ibc+1,jbc,kbc,4,n))/2
                          else if (i3.eq.iend(ic3)) then
                             u2 = (3*q(ibc,jbc,kbc,2,n) - q(ibc-1,jbc,kbc,2,n)) / 2
                             u1 = (  q(ibc,jbc,kbc,2,n) + q(ibc-1,jbc,kbc,2,n)) / 2
                             v2 = (3*q(ibc,jbc,kbc,3,n) - q(ibc-1,jbc,kbc,3,n)) / 2
                             v1 = (  q(ibc,jbc,kbc,3,n) + q(ibc-1,jbc,kbc,3,n)) / 2
                             w2 = (3*q(ibc,jbc,kbc,4,n) - q(ibc-1,jbc,kbc,4,n)) / 2
                             w1 = (  q(ibc,jbc,kbc,4,n) + q(ibc-1,jbc,kbc,4,n)) / 2
                          else
                             u2 = (q(ibc+1,jbc,kbc,2,n) + q(ibc  ,jbc,kbc,2,n)) / 2
                             u1 = (q(ibc  ,jbc,kbc,2,n) + q(ibc-1,jbc,kbc,2,n)) / 2
                             v2 = (q(ibc+1,jbc,kbc,3,n) + q(ibc  ,jbc,kbc,3,n)) / 2
                             v1 = (q(ibc  ,jbc,kbc,3,n) + q(ibc-1,jbc,kbc,3,n)) / 2
                             w2 = (q(ibc+1,jbc,kbc,4,n) + q(ibc  ,jbc,kbc,4,n)) / 2
                             w1 = (q(ibc  ,jbc,kbc,4,n) + q(ibc-1,jbc,kbc,4,n)) / 2
                          end if
                          dudxi = u2 - u1
                          dvdxi = v2 - v1
                          dwdxi = w2 - w1
                          
                          dudeta = (-1)**(1+1/inrout) * &
                               (-8*q(ibc,jbc,kbc,2,n) + &
                               9*q(ibc,jbc+joff1,kbc,2,n) - &
                               q(ibc,jbc+joff2,kbc,2,n) ) / 3
                          dvdeta = (-1)**(1+1/inrout) * &
                               (-8*q(ibc,jbc,kbc,3,n) + &
                               9*q(ibc,jbc+joff1,kbc,3,n) - &
                               q(ibc,jbc+joff2,kbc,3,n) ) / 3
                          dwdeta = (-1)**(1+1/inrout) * &
                               (-8*q(ibc,jbc,kbc,4,n) + &
                               9*q(ibc,jbc+joff1,kbc,4,n) - &
                               q(ibc,jbc+joff2,kbc,4,n) ) / 3
                          
                          if (i2.eq.istrt(ic2)) then
                             u2 = (  q(ibc,jbc,kbc+1,2,n) + q(ibc,jbc,kbc  ,2,n))/2
                             u1 = (3*q(ibc,jbc,kbc  ,2,n) - q(ibc,jbc,kbc+1,2,n))/2
                             v2 = (  q(ibc,jbc,kbc+1,3,n) + q(ibc,jbc,kbc  ,3,n))/2
                             v1 = (3*q(ibc,jbc,kbc  ,3,n) - q(ibc,jbc,kbc+1,3,n))/2
                             w2 = (  q(ibc,jbc,kbc+1,4,n) + q(ibc,jbc,kbc  ,4,n))/2
                             w1 = (3*q(ibc,jbc,kbc  ,4,n) - q(ibc,jbc,kbc+1,4,n))/2
                          else if (i2.eq.iend(ic2)) then
                             u2 = (3*q(ibc,jbc,kbc,2,n) - q(ibc,jbc,kbc-1,2,n)) / 2
                             u1 = (  q(ibc,jbc,kbc,2,n) + q(ibc,jbc,kbc-1,2,n)) / 2
                             v2 = (3*q(ibc,jbc,kbc,3,n) - q(ibc,jbc,kbc-1,3,n)) / 2
                             v1 = (  q(ibc,jbc,kbc,3,n) + q(ibc,jbc,kbc-1,3,n)) / 2
                             w2 = (3*q(ibc,jbc,kbc,4,n) - q(ibc,jbc,kbc-1,4,n)) / 2
                             w1 = (  q(ibc,jbc,kbc,4,n) + q(ibc,jbc,kbc-1,4,n)) / 2
                          else
                             u2 = (q(ibc,jbc,kbc+1,2,n) + q(ibc,jbc,kbc  ,2,n)) / 2
                             u1 = (q(ibc,jbc,kbc  ,2,n) + q(ibc,jbc,kbc-1,2,n)) / 2
                             v2 = (q(ibc,jbc,kbc+1,3,n) + q(ibc,jbc,kbc  ,3,n)) / 2
                             v1 = (q(ibc,jbc,kbc  ,3,n) + q(ibc,jbc,kbc-1,3,n)) / 2
                             w2 = (q(ibc,jbc,kbc+1,4,n) + q(ibc,jbc,kbc  ,4,n)) / 2
                             w1 = (q(ibc,jbc,kbc  ,4,n) + q(ibc,jbc,kbc-1,4,n)) / 2
                          end if
                          dudzeta = u2 - u1
                          dvdzeta = v2 - v1
                          dwdzeta = w2 - w1
                          
                          xix   = xiderj  (1,im,jm,km,n)
                          xiy   = xiderj  (2,im,jm,km,n)
                          xiz   = xiderj  (3,im,jm,km,n)
                          etax  = etaderj (1,im,jm,km,n)
                          etay  = etaderj (2,im,jm,km,n)
                          etaz  = etaderj (3,im,jm,km,n)
                          zetax = zetaderj(1,im,jm,km,n)
                          zetay = zetaderj(2,im,jm,km,n)
                          zetaz = zetaderj(3,im,jm,km,n)
                          
                       else if (idir.eq.3) then
                          if (i2.eq.istrt(ic2)) then
                             u2 = (  q(ibc,jbc,kbc,2,n) + q(ibc+1,jbc,kbc,2,n)) / 2
                             u1 = (3*q(ibc,jbc,kbc,2,n) - q(ibc+1,jbc,kbc,2,n)) / 2
                             v2 = (  q(ibc,jbc,kbc,3,n) + q(ibc+1,jbc,kbc,3,n)) / 2
                             v1 = (3*q(ibc,jbc,kbc,3,n) - q(ibc+1,jbc,kbc,3,n)) / 2
                             w2 = (  q(ibc,jbc,kbc,4,n) + q(ibc+1,jbc,kbc,4,n)) / 2
                             w1 = (3*q(ibc,jbc,kbc,4,n) - q(ibc+1,jbc,kbc,4,n)) / 2
                          else if (i2.eq.iend(ic2)) then
                             u2 = (3*q(ibc,jbc,kbc,2,n) - q(ibc-1,jbc,kbc,2,n)) / 2
                             u1 = (  q(ibc,jbc,kbc,2,n) + q(ibc-1,jbc,kbc,2,n)) / 2
                             v2 = (3*q(ibc,jbc,kbc,3,n) - q(ibc-1,jbc,kbc,3,n)) / 2
                             v1 = (  q(ibc,jbc,kbc,3,n) + q(ibc-1,jbc,kbc,3,n)) / 2
                             w2 = (3*q(ibc,jbc,kbc,4,n) - q(ibc-1,jbc,kbc,4,n)) / 2
                             w1 = (  q(ibc,jbc,kbc,4,n) + q(ibc-1,jbc,kbc,4,n)) / 2
                          else
                             u2 = (q(ibc  ,jbc,kbc,2,n) + q(ibc+1,jbc,kbc,2,n)) / 2
                             u1 = (q(ibc-1,jbc,kbc,2,n) + q(ibc  ,jbc,kbc,2,n)) / 2
                             v2 = (q(ibc  ,jbc,kbc,3,n) + q(ibc+1,jbc,kbc,3,n)) / 2
                             v1 = (q(ibc-1,jbc,kbc,3,n) + q(ibc  ,jbc,kbc,3,n)) / 2
                             w2 = (q(ibc  ,jbc,kbc,4,n) + q(ibc+1,jbc,kbc,4,n)) / 2
                             w1 = (q(ibc-1,jbc,kbc,4,n) + q(ibc  ,jbc,kbc,4,n)) / 2
                          end if
                          dudxi = u2 - u1
                          dvdxi = v2 - v1
                          dwdxi = w2 - w1
                          
                          if (i3.eq.istrt(ic3)) then
                             u2 = (  q(ibc,jbc,kbc,2,n) + q(ibc,jbc+1,kbc,2,n)) / 2
                             u1 = (3*q(ibc,jbc,kbc,2,n) - q(ibc,jbc+1,kbc,2,n)) / 2
                             v2 = (  q(ibc,jbc,kbc,3,n) + q(ibc,jbc+1,kbc,3,n)) / 2
                             v1 = (3*q(ibc,jbc,kbc,3,n) - q(ibc,jbc+1,kbc,3,n)) / 2
                             w2 = (  q(ibc,jbc,kbc,4,n) + q(ibc,jbc+1,kbc,4,n)) / 2
                             w1 = (3*q(ibc,jbc,kbc,4,n) - q(ibc,jbc+1,kbc,4,n)) / 2
                          else if (i3.eq.iend(ic3)) then
                             u2 = (3*q(ibc,jbc,kbc,2,n) - q(ibc,jbc-1,kbc,2,n)) / 2
                             u1 = (  q(ibc,jbc,kbc,2,n) + q(ibc,jbc-1,kbc,2,n)) / 2
                             v2 = (3*q(ibc,jbc,kbc,3,n) - q(ibc,jbc-1,kbc,3,n)) / 2
                             v1 = (  q(ibc,jbc,kbc,3,n) + q(ibc,jbc-1,kbc,3,n)) / 2
                             w2 = (3*q(ibc,jbc,kbc,4,n) - q(ibc,jbc-1,kbc,4,n)) / 2
                             w1 = (  q(ibc,jbc,kbc,4,n) + q(ibc,jbc-1,kbc,4,n)) / 2
                          else
                             u2 = (q(ibc,jbc  ,kbc,2,n) + q(ibc,jbc+1,kbc,2,n)) / 2
                             u1 = (q(ibc,jbc-1,kbc,2,n) + q(ibc,jbc  ,kbc,2,n)) / 2
                             v2 = (q(ibc,jbc  ,kbc,3,n) + q(ibc,jbc+1,kbc,3,n)) / 2
                             v1 = (q(ibc,jbc-1,kbc,3,n) + q(ibc,jbc  ,kbc,3,n)) / 2
                             w2 = (q(ibc,jbc  ,kbc,4,n) + q(ibc,jbc+1,kbc,4,n)) / 2
                             w1 = (q(ibc,jbc-1,kbc,4,n) + q(ibc,jbc  ,kbc,4,n)) / 2
                          end if
                          dudeta = u2 - u1
                          dvdeta = v2 - v1
                          dwdeta = w2 - w1
                          
                          dudzeta = (-1)**(1+1/inrout) * &
                               (-8*q(ibc,jbc,kbc,2,n) + &
                               9*q(ibc,jbc,kbc+koff1,2,n) - &
                               q(ibc,jbc,kbc+koff2,2,n) ) / 3
                          dvdzeta = (-1)**(1+1/inrout) * &
                               (-8*q(ibc,jbc,kbc,3,n) + &
                               9*q(ibc,jbc,kbc+koff1,3,n) - &
                               q(ibc,jbc,kbc+koff2,3,n) ) / 3
                          dwdzeta = (-1)**(1+1/inrout) * &
                               (-8*q(ibc,jbc,kbc,4,n) + &
                               9*q(ibc,jbc,kbc+koff1,4,n) - &
                               q(ibc,jbc,kbc+koff2,4,n) ) / 3
                          
                          xix   = xiderk  (1,im,jm,km,n)
                          xiy   = xiderk  (2,im,jm,km,n)
                          xiz   = xiderk  (3,im,jm,km,n)
                          etax  = etaderk (1,im,jm,km,n)
                          etay  = etaderk (2,im,jm,km,n)
                          etaz  = etaderk (3,im,jm,km,n)
                          zetax = zetaderk(1,im,jm,km,n)
                          zetay = zetaderk(2,im,jm,km,n)
                          zetaz = zetaderk(3,im,jm,km,n)
                          
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
                       
                       cf(ibc,jbc,kbc) = 2*tauw / (reyno*machfs)
                       yp(ibc,jbc,kbc) = dsqrt(reyno/machfs)* &
                            rhow*dn1*utau/muw
                       
                    end if
                    
                 end do
              end do
              
              if(.not.viscous) then
                 
                 do k=kbc_s,kbc_e
                    do j=jbc_s,jbc_e
                       do i=ibc_s,ibc_e
                          cf(i,j,k) = 0.d0
                       end do
                    end do
                 end do                    
                 
                 do k=kbc_s,kbc_e
                    do j=jbc_s,jbc_e
                       do i=ibc_s,ibc_e
                          yp(i,j,k) = 0.d0
                       end do
                    end do
                 end do                    

              end if      

              xnum = 1
              call cgp_coord_write_data_f(fid(n),basenum,surfacenumber,xnum, &
                   minrange,maxrange,x(ibc_s:ibc_e,jbc_s:jbc_e,kbc_s:kbc_e,n),ierr)
              if(ierr .ne. CG_OK) then
                 write(*,*) 'X cgp_coord_write_data_f error'
                 call cg_error_print_f()
              end if
              ynum = 2
              call cgp_coord_write_data_f(fid(n),basenum,surfacenumber,ynum, &
                   minrange,maxrange,y(ibc_s:ibc_e,jbc_s:jbc_e,kbc_s:kbc_e,n),ierr)
              if(ierr .ne. CG_OK) then
                 write(*,*) 'Y cgp_coord_write_data_f error'
                 call cg_error_print_f()
              end if
              znum = 3
              call cgp_coord_write_data_f(fid(n),basenum,surfacenumber,znum, &
                   minrange,maxrange,z(ibc_s:ibc_e,jbc_s:jbc_e,kbc_s:kbc_e,n),ierr)
              if(ierr .ne. CG_OK) then
                 write(*,*) 'Z cgp_coord_write_data_f error'
                 call cg_error_print_f()
              end if
              solnum = 1
              flownum = 1
              call cgp_field_write_data_f(fid(n),basenum,surfacenumber,solnum, &
                   flownum,minrange,maxrange,cp(ibc_s:ibc_e,jbc_s:jbc_e,kbc_s:kbc_e),ierr)
              if(ierr .ne. CG_OK) then
                 write(*,*) 'cgp_field_write_data_f cp error'
                 call cg_error_print_f()
              end if
              flownum = flownum + 1
              call cgp_field_write_data_f(fid(n),basenum,surfacenumber,solnum, &
                   flownum,minrange,maxrange,cf(ibc_s:ibc_e,jbc_s:jbc_e,kbc_s:kbc_e),ierr)
              if(ierr .ne. CG_OK) then
                 write(*,*) 'cgp_field_write_data_f cf error'
                 call cg_error_print_f()
              end if
              flownum = flownum + 1
              call cgp_field_write_data_f(fid(n),basenum,surfacenumber,solnum, &
                   flownum,minrange,maxrange,yp(ibc_s:ibc_e,jbc_s:jbc_e,kbc_s:kbc_e),ierr)
              if(ierr .ne. CG_OK) then
                 write(*,*) 'cgp_field_write_data_f yp error'
                 call cg_error_print_f()
              end if
               
              if(allocated(cp)) deallocate(cp)              
              if(allocated(cf)) deallocate(cf,yp) 

              surfacenumber = surfacenumber + 1
              
           end if
        end do

      end do      

      return
      end

!-----------------------------------------------------------------------
      subroutine writeparallelflowteccgnsheader(fid,basenum,nl,n)
!     This opens the cgns file and writes the header text into it 
!     for the parallel flowtec write.
!-----------------------------------------------------------------------
      
      use cgns
      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

#if MPI 
      include 'mpif.h'
#endif
      integer(kind=cosa_int) :: fid(0:2*nharms)
      integer(kind=cosa_int) :: basenum, numblocks, blocknum
      integer(kind=cosa_int) :: imax, jmax, kmax, solnum, flownum, arrnum
      integer(kind=cosa_int) :: xnum, ynum, znum
      integer(kind=cosa_int) :: ierr,n,nl,i
      character*20 zonename
      integer(cgsize_t) :: sizes(3,3)

      logical :: amcontrol

      call amcontroller(amcontrol)

      call cgp_pio_mode_f(CGP_INDEPENDENT, ierr)
      if(ierr .ne. CG_OK) then
         write(*,*) 'cgp_pio_mode_f error writeparallelflowtecgnsheader'
         call cg_error_print_f()
         call abortmpi()
      end if

      call cgp_open_f(flowtec,CG_MODE_WRITE,fid(n),ierr)
      if(ierr .ne. CG_OK) then
         write(*,*) 'cg_open_f error'
         call cg_error_print_f()
      end if
      call cg_base_write_f(fid(n),'gridbase',3,3,basenum,ierr)
      if(ierr .ne. CG_OK) then
         write(*,*) 'cg_base_write_f error'
         call cg_error_print_f()
      end if

!     AJ Write all the blocks as zones.
      do i = 1,nblocks
         blocknum = i
         
         imax = g_i_imax(blocknum,nl)
         jmax = g_j_jmax(blocknum,nl)
         kmax = g_k_kmax(blocknum,nl)         
         
         sizes(1,1) = imax+1
         sizes(2,1) = jmax+1
         sizes(3,1) = kmax+1
         
         sizes(1,2) = imax
         sizes(2,2) = jmax
         sizes(3,2) = kmax
         
         sizes(1,3) = 0
         sizes(2,3) = 0
         sizes(3,3) = 0
         write(zonename, "(A5,I6)") "block",blocknum
         call cg_zone_write_f(fid(n),basenum,zonename,sizes, &
              Structured,blocknum,ierr)
         if(ierr .ne. CG_OK) then
            write(*,*) &
                 'cg_zone_write_f error in writeparallelflowtecgnsheader'
            call cg_error_print_f()
         end if
         call cg_sol_write_f(fid(n),basenum,blocknum,'FlowSolution1', &
              Vertex,solnum,ierr)
!     &        CellCenter,solnum,ierr)
         if(ierr .ne. CG_OK) then
            write(*,*) &
                 'cg_sol_write_f error in writeparallelflowtecgnsheader'
            call cg_error_print_f()
         end if
         call cgp_coord_write_f(fid(n),basenum,blocknum,RealDouble, &
              'CoordinateX',xnum,ierr)
         if(ierr .ne. CG_OK) then
            write(*,*) 'X cgp_coord_write_f error'
            call cg_error_print_f()
         end if
         call cgp_coord_write_f(fid(n),basenum,blocknum,RealDouble, &
              'CoordinateY',ynum,ierr)
         if(ierr .ne. CG_OK) then
            write(*,*) 'Y cgp_coord_write_f error'
            call cg_error_print_f()
         end if
         call cgp_coord_write_f(fid(n),basenum,blocknum,RealDouble, &
              'CoordinateZ',znum,ierr)
         if(ierr .ne. CG_OK) then
            write(*,*) 'Z cgp_coord_write_f error'
            call cg_error_print_f()
         end if
         call cgp_field_write_f(fid(n),basenum,blocknum,solnum, &
              RealDouble,'Density',flownum,ierr)
         if(ierr .ne. CG_OK) then
            write(*,*) 'cgp_field_write_f density error'
            call cg_error_print_f()
         end if
         call cgp_field_write_f(fid(n),basenum,blocknum,solnum, &
              RealDouble,'VelocityX',flownum,ierr)
         if(ierr .ne. CG_OK) then
            write(*,*) 'cgp_field_write_f velocityX error'
            call cg_error_print_f()
         end if
         call cgp_field_write_f(fid(n),basenum,blocknum,solnum, &
              RealDouble,'VelocityY',flownum,ierr)
         if(ierr .ne. CG_OK) then
            write(*,*) 'cgp_field_write_f velocityY error'
            call cg_error_print_f()
         end if
         call cgp_field_write_f(fid(n),basenum,blocknum,solnum, &
              RealDouble,'VelocityZ',flownum,ierr)
         if(ierr .ne. CG_OK) then
            write(*,*) 'cgp_field_write_f velocityZ error'
            call cg_error_print_f()
         end if
         call cgp_field_write_f(fid(n),basenum,blocknum,solnum, &
              RealDouble,'Pressure',flownum,ierr)
         if(ierr .ne. CG_OK) then
            write(*,*) 'cgp_field_write_f pressure error'
            call cg_error_print_f()
         end if
         call cgp_field_write_f(fid(n),basenum,blocknum,solnum, &
              RealDouble,'Temperature',flownum,ierr)
         if(ierr .ne. CG_OK) then
            write(*,*) 'cgp_field_write_f temperature error'
            call cg_error_print_f()
         end if
         call cgp_field_write_f(fid(n),basenum,blocknum,solnum, &
              RealDouble,'Mach',flownum,ierr)
         if(ierr .ne. CG_OK) then
            write(*,*) 'cgp_field_write_f mach error'
            call cg_error_print_f()
         end if
         if (kom.or.kom_bsl.or.kom_sst) then                     
            call cgp_field_write_f(fid(n),basenum,blocknum,solnum, &
                 RealDouble,'TurbulentEnergyKinetic',flownum, &
                 ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) 'cgp_field_write_f turbulent kinetic error'
               call cg_error_print_f()
            end if
            call cgp_field_write_f(fid(n),basenum,blocknum,solnum, &
                 RealDouble,'TurbulentDissipationRate',flownum,ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) 'cgp_field_write_f turbulent disrate error'
               call cg_error_print_f()
            end if
            call cgp_field_write_f(fid(n),basenum,blocknum,solnum, &
                 RealDouble,'TurbulentDistance',flownum,ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) 'cgp_field_write_f turbulent distance error'
               call cg_error_print_f()
            end if
            call cgp_field_write_f(fid(n),basenum,blocknum,solnum, &
                 RealDouble,'ViscosityEddy',flownum,ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) 'cgp_field_write_f viscosity eddy error'
               call cg_error_print_f()
            end if
         end if
      end do
      
      call cgp_close_f(fid(n),ierr)
      if(ierr .ne. CG_OK) then
         write(*,*) 'cg_close_f error'
         call cg_error_print_f()
      end if      

      return 
      end



!-----------------------------------------------------------------------
      subroutine write_parallel_wr_cgns_b(fid,q,mut,x,y,z,xdot,ydot, &
           zdot,dist,imax,jmax,kmax,npde,nharms,basenum,zonenum)
!     This routine writes parts of the flowtec files using parallel CGNS
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

#if MPI 
      include 'mpif.h'
#endif

      integer(kind=cosa_int) :: imax,jmax,kmax,npde,nharms,basenum,zonenum
      integer(kind=cosa_int) :: fid(0:2*nharms)
      real(kind=cosa_real) :: &
           q   (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           mut (-1:imax+1,-1:jmax+1,-1:kmax+1     ,0:2*nharms), &
           x   ( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           y   ( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           z   ( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           xdot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           ydot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           zdot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           dist( 0:imax  , 0:jmax  , 0:kmax)

      call parallel_wr_cgns_b(fid,q,mut,x,y,z,xdot,ydot,zdot,dist,imax, &
           jmax,kmax,npde,nharms,basenum,zonenum)

      return 
      end

!-----------------------------------------------------------------------
      subroutine parallel_wr_cgns_b(fid,q,mut,x,y,z,xdot,ydot, &
           zdot,dist,imax,jmax,kmax,npde,nharms,basenum,zonenum)
!     This routine writes parts of the flowtec files using MPI-I/O.  
!     It is bespoke for the flowtec files and uses functionality that is 
!     also used in the serial I/O code.  This means if you change 
!     any I/O functionality here (which is not just the MPI-I/O code)
!     then it should also be changed in the same code in output.f
!     (i.e. subroutine wr_tec_b).
!-----------------------------------------------------------------------

      use cgns
      use common_variables
      use cosa_precision

      implicit none

#if MPI 
      include 'mpif.h'
#endif

      integer(kind=cosa_int) :: imax,jmax,kmax,npde,nharms
      integer(kind=cosa_int) :: i,imax1,j,jmax1,k,kmax1,ipde,n,nh
      integer(kind=cosa_int) :: ierr,zonenum,basenum,solnum,flownum,arrnum
      integer(kind=cosa_int) :: fid(0:2*nharms)
      integer(kind=cosa_int) :: xnum,ynum,znum
      integer(cgsize_t) :: minrange(3), maxrange(3)
      real(kind=cosa_real) :: &
           q   (-1:imax+1,-1:jmax+1,-1:kmax+1,npde,0:2*nharms), &
           mut (-1:imax+1,-1:jmax+1,-1:kmax+1     ,0:2*nharms), &
           x   ( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           y   ( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           z   ( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           xdot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           ydot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           zdot( 0:imax+1, 0:jmax+1, 0:kmax+1     ,0:2*nharms*hbmove), &
           dist( 0:imax  , 0:jmax  , 0:kmax), &
           localx(0:imax  , 0:jmax  , 0:kmax), &
           localy(0:imax  , 0:jmax  , 0:kmax), &
           localz(0:imax  , 0:jmax  , 0:kmax), &
           velocityX(0:imax  , 0:jmax  , 0:kmax), &
           velocityY(0:imax  , 0:jmax  , 0:kmax), &
           velocityZ(0:imax  , 0:jmax  , 0:kmax), &
           temperature(0:imax  , 0:jmax  , 0:kmax), &
           density(0:imax  , 0:jmax  , 0:kmax), &
           pressure(0:imax  , 0:jmax  , 0:kmax), &
           turbulentkinetic(0:imax  , 0:jmax  , 0:kmax), &
           turbulentdissipation(0:imax  , 0:jmax  , 0:kmax), &
           mach(0:imax  , 0:jmax  , 0:kmax), &
           localmut(0:imax  , 0:jmax  , 0:kmax)

      imax1 = imax+1
      jmax1 = jmax+1
      kmax1 = kmax+1

      do n = 0,2*nharms
         nh = n*hbmove               
         
         minrange = 1
         maxrange(1) = imax+1
         maxrange(2) = jmax+1
         maxrange(3) = kmax+1
         
         if (kom.or.kom_bsl.or.kom_sst) then                     
            
            do k=0,kmax
               do j=0,jmax
                  do i=0,imax
                     localx(i,j,k) = x(i,j,k,nh)
                     localy(i,j,k) = y(i,j,k,nh)
                     localz(i,j,k) = z(i,j,k,nh)
                     if (moving) then
                        velocityX(i,j,k) = q(i,j,k,2,n) - xdot(i,j,k,nh)
                        velocityY(i,j,k) = q(i,j,k,3,n) - ydot(i,j,k,nh)
                        velocityZ(i,j,k) = q(i,j,k,4,n) - zdot(i,j,k,nh)
                     else
                        velocityX(i,j,k) = q(i,j,k,2,n)
                        velocityY(i,j,k) = q(i,j,k,3,n)
                        velocityZ(i,j,k) = q(i,j,k,4,n)
                     end if
! AJ TODO: Check whether these copies are necessary or can CGNS pull these ranges out automagically in the write calls
                     density(i,j,k) = q(i,j,k,1,n)
                     pressure(i,j,k) = q(i,j,k,5,n)
                     turbulentkinetic(i,j,k) = q(i,j,k,6,n)
                     turbulentdissipation(i,j,k) = q(i,j,k,7,n)
                     temperature(i,j,k) = gamma * q(i,j,k,5,n) / &
                          q(i,j,k,1,n)
                     mach(i,j,k) = &
                          sqrt(velocityX(i,j,k)**2+velocityY(i,j,k)**2+ &
                          velocityZ(i,j,k)**2)/ sqrt(temperature(i,j,k))
                     localmut(i,j,k) = mut(i,j,k,n)
                  end do
               end do
            end do
            xnum = 1
            call cgp_coord_write_data_f(fid(n),basenum,zonenum,xnum, &
                 minrange,maxrange,localx,ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) 'X cgp_coord_write_data_f error'
               call cg_error_print_f()
            end if
            ynum = 2
            call cgp_coord_write_data_f(fid(n),basenum,zonenum,ynum, &
                 minrange,maxrange,localy,ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) 'Y cgp_coord_write_data_f error'
               call cg_error_print_f()
            end if
            znum = 3
            call cgp_coord_write_data_f(fid(n),basenum,zonenum,znum, &
                 minrange,maxrange,localz,ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) 'Z cgp_coord_write_data_f error'
               call cg_error_print_f()
            end if
            solnum = 1
            flownum = 1
            call cgp_field_write_data_f(fid(n),basenum,zonenum,solnum, &
                 flownum,minrange,maxrange,density,ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) 'cgp_field_write_data_f density error'
               call cg_error_print_f()
            end if
            flownum = flownum + 1
            call cgp_field_write_data_f(fid(n),basenum,zonenum,solnum, &
                 flownum,minrange,maxrange,velocityX,ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) 'cgp_field_write_data_f velocityX error'
               call cg_error_print_f()
            end if
            flownum = flownum + 1
            call cgp_field_write_data_f(fid(n),basenum,zonenum,solnum, &
                 flownum,minrange,maxrange,velocityY,ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) 'cgp_field_write_data_f velocityY error'
               call cg_error_print_f()
            end if
            flownum = flownum + 1
            call cgp_field_write_data_f(fid(n),basenum,zonenum,solnum, &
                 flownum,minrange,maxrange,velocityZ,ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) 'cgp_field_write_data_f velocityZ error'
               call cg_error_print_f()
            end if
            flownum = flownum + 1
            call cgp_field_write_data_f(fid(n),basenum,zonenum,solnum, &
                 flownum,minrange,maxrange,pressure,ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) 'cgp_field_write_data_f pressure error'
               call cg_error_print_f()
            end if
            flownum = flownum + 1
            call cgp_field_write_data_f(fid(n),basenum,zonenum,solnum, &
                 flownum,minrange,maxrange,temperature,ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) 'cgp_field_write_data_f temperature error'
               call cg_error_print_f()
            end if
            flownum = flownum + 1
            call cgp_field_write_data_f(fid(n),basenum,zonenum,solnum, &
                 flownum,minrange,maxrange,mach,ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) 'cgp_field_write_data_f mach error'
               call cg_error_print_f()
            end if
            flownum = flownum + 1
             call cgp_field_write_data_f(fid(n),basenum,zonenum,solnum, &
                 flownum,minrange,maxrange,turbulentkinetic,ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) &
                    'cgp_field_write_data_f turbulent kinetic error'
               call cg_error_print_f()
            end if
            flownum = flownum + 1
            call cgp_field_write_data_f(fid(n),basenum,zonenum,solnum, &
                 flownum,minrange,maxrange,turbulentdissipation,ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) &
                    'cgp_field_write_data_f turbulent disrate error'
               call cg_error_print_f()
            end if
            flownum = flownum + 1
            call cgp_field_write_data_f(fid(n),basenum,zonenum,solnum, &
                 flownum,minrange,maxrange,dist,ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) &
                    'cgp_field_write_data_f turbulent distance error'
               call cg_error_print_f()
            end if
            flownum = flownum + 1
            call cgp_field_write_data_f(fid(n),basenum,zonenum,solnum, &
                 flownum,minrange,maxrange,localmut,ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) &
                    'cgp_field_write_data_f viscosity eddy error'
               call cg_error_print_f()
            end if
            
         else

            do k=0,kmax
               do j=0,jmax
                  do i=0,imax
                     localx(i,j,k) = x(i,j,k,nh)
                     localy(i,j,k) = y(i,j,k,nh)
                     localz(i,j,k) = z(i,j,k,nh)
                     if (moving) then
                        velocityX(i,j,k) = q(i,j,k,2,n) - xdot(i,j,k,nh)
                        velocityY(i,j,k) = q(i,j,k,3,n) - ydot(i,j,k,nh)
                        velocityZ(i,j,k) = q(i,j,k,4,n) - zdot(i,j,k,nh)
                     else
                        velocityX(i,j,k) = q(i,j,k,2,n)
                        velocityY(i,j,k) = q(i,j,k,3,n)
                        velocityZ(i,j,k) = q(i,j,k,4,n)
                     end if
                     density(i,j,k) = q(i,j,k,1,n)
                     pressure(i,j,k) = q(i,j,k,5,n)
                     temperature(i,j,k) = gamma * q(i,j,k,5,n) / &
                          q(i,j,k,1,n)
                     mach(i,j,k) = &
                          sqrt(velocityX(i,j,k)**2+velocityY(i,j,k)**2+ &
                          velocityZ(i,j,k)**2)/ sqrt(temperature(i,j,k))
                  end do
               end do
            end do

            xnum = 1
            call cgp_coord_write_data_f(fid(n),basenum,zonenum,xnum, &
                 minrange,maxrange,localx,ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) 'X cgp_coord_write_data_f error'
               call cg_error_print_f()
            end if
            ynum = 2
            call cgp_coord_write_data_f(fid(n),basenum,zonenum,ynum, &
                 minrange,maxrange,localy,ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) 'Y cgp_coord_write_data_f error'
               call cg_error_print_f()
            end if
            znum = 3
            call cgp_coord_write_data_f(fid(n),basenum,zonenum,znum, &
                 minrange,maxrange,localz,ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) 'Z cgp_coord_write_data_f error'
               call cg_error_print_f()
            end if
            solnum = 1
            flownum = 1
            call cgp_field_write_data_f(fid(n),basenum,zonenum,solnum, &
                 flownum,minrange,maxrange,density,ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) 'cgp_field_write_data_f density error'
               call cg_error_print_f()
            end if
            flownum = flownum + 1
            call cgp_field_write_data_f(fid(n),basenum,zonenum,solnum, &
                 flownum,minrange,maxrange,velocityX,ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) 'cgp_field_write_data_f velocityX error'
               call cg_error_print_f()
            end if
            flownum = flownum + 1
            call cgp_field_write_data_f(fid(n),basenum,zonenum,solnum, &
                 flownum,minrange,maxrange,velocityY,ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) 'cgp_field_write_data_f velocityY error'
               call cg_error_print_f()
            end if
            flownum = flownum + 1
            call cgp_field_write_data_f(fid(n),basenum,zonenum,solnum, &
                 flownum,minrange,maxrange,velocityZ,ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) 'cgp_field_write_data_f velocityZ error'
               call cg_error_print_f()
            end if
            flownum = flownum + 1
            call cgp_field_write_data_f(fid(n),basenum,zonenum,solnum, &
                 flownum,minrange,maxrange,pressure,ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) 'cgp_field_write_data_f pressure error'
               call cg_error_print_f()
            end if
            flownum = flownum + 1
            call cgp_field_write_data_f(fid(n),basenum,zonenum,solnum, &
                 flownum,minrange,maxrange,temperature,ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) 'cgp_field_write_data_f temperature error'
               call cg_error_print_f()
            end if
            call cgp_field_write_data_f(fid(n),basenum,zonenum,solnum, &
                 flownum,minrange,maxrange,mach,ierr)
            if(ierr .ne. CG_OK) then
               write(*,*) 'cgp_field_write_data_f mach error'
               call cg_error_print_f()
            end if

         end if
        
      end do

      return
      end


!-----------------------------------------------------------------------
      subroutine write_parallel_simtime_tec(fid,nl,nblocks)
!     Write simtime at the end of the flowtec file.
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

#if MPI 
      include 'mpif.h'
      integer(kind=cosa_int) :: offset_kind
      parameter(offset_kind=MPI_OFFSET_KIND)
#else
      integer(kind=cosa_int) :: offset_kind
      parameter(offset_kind=8)
#endif

      integer(kind=cosa_int) :: nl,nblocks
      integer(kind=cosa_int) :: fid
      integer(kind=cosa_int) :: ierr,integersize,doublesize,charactersize
      integer(kind=offset_kind) :: disp

#if MPI 
      call newflowtecinitialoffset(disp,nl,nblocks+1,integersize, &
           doublesize,charactersize)
      call setupfile(fid,disp,MPI_DOUBLE_PRECISION)
      call mpi_file_write(fid,simtime,1, &
           MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
#endif

      return
      end
