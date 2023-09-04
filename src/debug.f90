!-----------------------------------------------------------------------
       subroutine wrqf(qf,idir,imax,jmax,npde,nharms)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,npde,nharms
      integer(kind=cosa_int) i,im,imax1,j,jm,jmax1,ipde,idir,n
      real (kind=cosa_real) qf(0:imax+1,0:jmax+1,npde,0:2*nharms)
      character*72 fileqf
       
      if (idir.eq.1) then
        fileqf = 'qf_1.dat'
        open(335,file=fileqf,status='replace',form='formatted')
        do n = 0,2*nharms
          do ipde = 1,npde
            do j = 1,jmax-1
              do i = 1,imax
                write(335,*) i,j,ipde,n,qf(i,j,ipde,n)
              end do
            end do
          end do
        end do
        close(335)
      else if (idir.eq.2) then
        fileqf = 'qf_2.dat'
        open(335,file=fileqf,status='replace',form='formatted')
        do n = 0,2*nharms
          do ipde = 1,npde
            do j = 1,jmax
              do i = 1,imax-1
                write(335,*) i,j,ipde,n,qf(i,j,ipde,n)
              end do
            end do
          end do
        end do
        close(335)
      end if

      return
      end

!-----------------------------------------------------------------------
       subroutine wrqder(dqxi,dqeta,idir,imax,jmax,npde,nharms)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,npde,nharms
      integer(kind=cosa_int) i,im,imax1,j,jm,jmax1,ipde,idir,n
      real (kind=cosa_real) dqxi(imax,jmax,npde,0:2*nharms), &
           dqeta(imax,jmax,npde,0:2*nharms)
      character*72 filedqxi,filedqeta
       
      if (idir.eq.1) then
        filedqxi  = 'dqxi_face_1.dat'
        filedqeta = 'dqeta_face_1.dat'
        open(335,file=filedqxi,status='replace',form='formatted')
        open(337,file=filedqeta,status='replace',form='formatted')
        do n = 0,2*nharms
          do ipde = 1,npde
            do j = 1,jmax-1
              do i = 1,imax
                write(335,*) i,j,ipde,n,dqxi (i,j,ipde,n)
                write(337,*) i,j,ipde,n,dqeta(i,j,ipde,n)
              end do
            end do
          end do
        end do
        close(335)
      else if (idir.eq.2) then
        filedqxi  = 'dqxi_face_2.dat'
        filedqeta = 'dqeta_face_2.dat'
        open(335,file=filedqxi,status='replace',form='formatted')
        open(337,file=filedqeta,status='replace',form='formatted')
        do n = 0,2*nharms
          do ipde=1,npde
            do j = 1,jmax
              do i = 1,imax-1
                write(335,*) i,j,ipde,n,dqxi (i,j,ipde,n)
                write(337,*) i,j,ipde,n,dqeta(i,j,ipde,n)
              end do
            end do
          end do
        end do
        close(335)
      end if

      return
      end

!-----------------------------------------------------------------------
       subroutine wrgrid(nl,x,y,z)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,iblock,imax,jmax,kmax,ixyz,fid
      real(kind=cosa_real) x(*),y(*),z(*)
      character*72 line

      fid = 1
      open(fid,file=filename,status='replace')

      write(line,'(''TITLE = "MG level '',i1,''level"'')') nl
      write(fid,'(a)') line
      write(fid,*) 'Variables = "x","y","z"'

      do iblock = 1,nblocks
        imax    = i_imax  (iblock,nl)
        jmax    = j_jmax  (iblock,nl)
        kmax    = k_kmax  (iblock,nl)
        ixyz    = 1 + off_p2(iblock,nl)
        call wrgrid_b(iblock,imax,jmax,kmax,fid,x(ixyz),y(ixyz),z(ixyz))
      end do

      close(fid)

      return
      end

!-----------------------------------------------------------------------
       subroutine wrgrid_b(iblock,imax,jmax,kmax,fid,x,y,z)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) iblock,imax,jmax,kmax
      integer(kind=cosa_int) i,j,k,fid
      real(kind=cosa_real) x(0:imax+1,0:jmax+1,0:kmax+1), &
        y(0:imax+1,0:jmax+1,0:kmax+1),z(0:imax+1,0:jmax+1,0:kmax+1)
      character*72 line

      write(line,'(''Zone T = "block '',i1,''"'')') iblock
      write(1,'(a)') line
      write(1,11) imax,jmax,kmax
      do k=1,kmax
        do j=1,jmax
          do i=1,imax
            write(fid,10) x(i,j,k), y(i,j,k), z(i,j,k)
          end do
        end do
      end do
!
 11   format ('I = ',I4,',J = ',I4,'K = ',I4,',F = POINT')
 10   format (2e22.12)

      return
      end

!-----------------------------------------------------------------------
       subroutine wrdist(imax,jmax,kmax,x,y,z,dist,nl)
!-----------------------------------------------------------------------

      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) imax,jmax,kmax
      integer(kind=cosa_int) i,j,k,nl
      real(kind=cosa_real) &
          x   (0:imax+1,0:jmax+1,0:kmax+1), &
          y   (0:imax+1,0:jmax+1,0:kmax+1), &
          z   (0:imax+1,0:jmax+1,0:kmax+1), &
          dist(0:imax  ,0:jmax  ,0:kmax)
      character*72 line

      open(1,file=filename,status='replace')

      write(line,'(''TITLE = "MG level '',i1,''level"'')') nl
      write(1,'(a)') line
      write(1,*) 'Variables = "x","y","z","dist"'

      write(line,'(''Zone T = "block '',i1,''"'')') 1
      write(1,'(a)') line
      write(1,11) imax,jmax,kmax
      do k=1,kmax
        do j=1,jmax
          do i=1,imax
            write(1,12) x(i,j,k), y(i,j,k), z(i,j,k), dist(i,j,k)
          end do
        end do
      end do
!
 11   format ('I = ',I4,',J = ',I4,'K = ',I4,',F = POINT')
 12   format (3e22.12)

      close(unit=1)

      return
      end

!-----------------------------------------------------------------------
      subroutine print_local_wall(nl,xgwall,ygwall,zgwall)
!-----------------------------------------------------------------------

      use cosa_variables
      use common_variables
      use cosa_precision

      implicit none

      integer(kind=cosa_int) nl,iw
      real (kind=cosa_real) xgwall(ng_wall(nl)),ygwall(ng_wall(nl)), &
                    zgwall(ng_wall(nl))

      write(filename,'(''lwall_nl'',i1,''.dat'')') nl
      open(335,file=filename,status='replace',form='formatted')

      do iw = 1,ng_wall(nl)
        write(335,*) xgwall(iw),ygwall(iw),zgwall(iw)
      end do

      close(335)

      return
      end
