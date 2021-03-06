        subroutine ElmGMRElas (x,       disp,    shp,     shgl,
     &                         iBC,     BC,      shpb,    shglb,
     &                         shpif,   elasres, elasBDiag,
     &                         iper,    ilwork,  elaslhsK,
     &                         col,     row,     meshq, gcnormal)
c
c----------------------------------------------------------------------
c
c This routine computes the LHS mass matrix, the RHS residual
c vector, and the preconditioning matrix, for mesh-elastic solve
c
c----------------------------------------------------------------------
c
         use pointer_data
         use timedataC
         use readarrays ! read BLflt, BLgr, BLtnv, BLlist
c
        include "common.h"
        include "mpif.h"
c
        integer col(nshg+1), row(nnz*nshg)
        real*8  elaslhsK(nelas*nelas,nnz_tot),
     &          meshq(numel),
     &          meshV(numel)
c
        dimension gcnormal(nshg, nsd)
c
        dimension x(numnp,nsd),        disp(numnp,nsd),
     &            xtmp(numnp,nsd),     iBC(nshg),
     &            BC(nshg,ndofBC),
     &            elasres(nshg,nelas),
     &            elasBDiag(nshg,nelas,nelas),
     &            iper(nshg)
c
        dimension shp(MAXTOP,maxsh,MAXQPT),
     &            shgl(MAXTOP,nsd,maxsh,MAXQPT),
     &            shpb(MAXTOP,maxsh,MAXQPT),
     &            shglb(MAXTOP,nsd,maxsh,MAXQPT)
c
        real*8, dimension(maxtop,    maxsh,maxqpt) :: shpif
c
        dimension ilwork(nlwork)
c
        integer errorcount(2)
c
        integer listcounter, ngc, itnv, basevID, nv, vID, vID2
        real*8  iflt, igr
        dimension inormal(nsd)
c
        real*8, allocatable :: tmpshp(:,:), tmpshgl(:,:,:)
        real*8, allocatable :: Estiff(:,:,:)
c
c.... ------------------->   layer base elements   <-------------------
c
c.... calculate the normal of each growth curve based on new boundary positions
c
        xtmp = x + disp
        gcnormal = zero
c
c.... loop over the boundary elements
c
        boundary_blocks: do iblk = 1, nelblb
c
c.... set up the parameters
c
          iel    = lcblkb(1,iblk)
          lelCat = lcblkb(2,iblk)
          lcsyst = lcblkb(3,iblk)
          iorder = lcblkb(4,iblk)
          nenl   = lcblkb(5,iblk)  ! no. of vertices per element
          nenbl  = lcblkb(6,iblk)  ! no. of vertices per bdry. face
          nshl   = lcblkb(9,iblk)
          nshlb  = lcblkb(10,iblk)
          mattyp = lcblkb(7,iblk)
          ndofl  = lcblkb(8,iblk)
          npro   = lcblkb(1,iblk+1) - iel
c
          if(lcsyst.eq.itp_wedge_tri) lcsyst=nenbl ! may not be necessary
          ngaussb = nintb(lcsyst)
c
c.... collect wedge_tri or surfID option is on
c
          if((lcsyst.ne.itp_wedge_tri) .and. (useBLbaseSrfID.eq.0)) cycle
c
c.... compute and assemble non-unit normal
c
          call calc_gc_normal (xtmp,             shpb(lcsyst,1:nshl,:),
     &                         mienb(iblk)%p,    miBCB(iblk)%p,
     &                         gcnormal)
c
        enddo boundary_blocks ! end loop the boundary elements
c
c.... loop over the interface elements
c
        interface_blocks: do iblk = 1, nelblif
c
          iel     = lcblkif(1, iblk)
          npro    = lcblkif(1,iblk+1) - iel
          lcsyst0 = lcblkif(3, iblk)    ! element0 type
          lcsyst1 = lcblkif(4, iblk)    ! element1 type
          ipord   = lcblkif(5, iblk)    ! polynomial order
          nenl0   = lcblkif(6, iblk)    ! number of vertices per element0
          nenl1   = lcblkif(7, iblk)    ! number of vertices per element1
          mater0  = lcblkif(9, iblk)
          mater1  = lcblkif(10,iblk)
          nshl0   = lcblkif(iblkif_nshl0,iblk)
          nshl1   = lcblkif(iblkif_nshl1,iblk)
          itpid   = lcblkif(iblkif_topology,iblk)
          ngaussif = nintif(itpid)
c
c.... the 0 side
c
          lcsyst = lcsyst0
          nenl = nenl0
          nshl = nshl0

c          if(lcsyst.eq.itp_wedge_tri) lcsyst=nenbl ! may not be necessary
          ngaussb = nintb(lcsyst)
c
c.... collect wedge_tri or surfID option is on
c
          if((lcsyst.ne.itp_wedge_tri) .and. (useBLbaseSrfID.eq.0)) cycle
c
c.... compute and assemble non-unit normal
c
          call calc_gc_normal (xtmp,            shpif(lcsyst,1:nshl,:),
     &                         mienif0(iblk)%p, miBCB(iblk)%p,
     &                         gcnormal)
c
c.... end of the 0 side
c
c
c.... the 1 side
c
          lcsyst = lcsyst1
          nenl = nenl1
          nshl = nshl1

c          if(lcsyst.eq.itp_wedge_tri) lcsyst=nenbl ! may not be necessary
          ngaussb = nintb(lcsyst)
c
c.... collect wedge_tri or surfID option is on
c
          if((lcsyst.ne.itp_wedge_tri) .and. (useBLbaseSrfID.eq.0)) cycle
c
c.... compute and assemble non-unit normal
c
          call calc_gc_normal (xtmp,            shpif(lcsyst,1:nshl,:),
     &                         mienif1(iblk)%p, miBCB(iblk)%p,
     &                         gcnormal)
c
c.... end of the 1 side
c
        enddo interface_blocks ! end loop the interface elements
c
c.... communication
c
        if (numpe > 1) then
          call commu (gcnormal  , ilwork, nsd  , 'in ')
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
          call commu (gcnormal  , ilwork, nsd  , 'out ')
        endif
c
c.... end calculation of growth curve normal
c
c.... ---------------->   Re-position layered mesh   <-----------------
c
c.... loop over growth curves
c
        listconuter = 0
        do ngc = 1, numgc
          itnv = BLtnv(ngc) ! total number of vertices for this growth curve
c
c.... precaution
c
          if (itnv .lt. 2) then
            listconuter = listconuter + itnv
            cycle ! not loop over vertices
          endif
c
c.... prepare other paramteres
c
          iflt = BLflt(ngc) ! first layer thickness for this growth curve
          igr  = BLgr(ngc)  ! growth ratio for this growth curve
          basevID = BLlist(listconuter + 1)
          inormal = gcnormal(basevID,:)
c
c.... loop over vertices on this growth curve
c
          do nv = 2, itnv
            vID = BLlist(listconuter + nv)
            vID2= BLlist(listconuter + nv - 1) ! the previous one
            xtmp(vID,:) = xtmp(vID2,:) + iflt * inormal(:) * igr**(nv-2)
            disp(vID,:) = xtmp(vID,:) - x(vID,:)
c
c.... assign iBC and BC arrays
c
            call assign_bl_bc( disp, iBC, BC(:,ndof+2:ndof+5),
     &                         basevID, vID )
c
c.... end loop vertices on this growth curve
c
          enddo
c
          listconuter = listconuter + itnv ! update counter
c
c.... end loop growth curves
c
        enddo
c
c.... end re-position layered mesh
c
c.... -------------------->   interior elements   <--------------------
c
c.... loop over element blocks to compute element residuals
c
c
c.... initialize the arrays
c
        elasres = zero
        errorcount = zero
        if (lhs. eq. 1)    elaslhsK  = zero
        if (iprec .ne. 0)  elasBDiag = zero
c
c.... loop over the element-blocks
c
        do iblk = 1, nelblk
c
c.... set up the parameters
c
          iblkts = iblk            ! used in timeseries
          iel    = lcblk(1,iblk)
          lelCat = lcblk(2,iblk)
          lcsyst = lcblk(3,iblk)
          iorder = lcblk(4,iblk)
          nenl   = lcblk(5,iblk)   ! no. of vertices per element
          nshl   = lcblk(10,iblk)
          mattyp = lcblk(7,iblk)
          ndofl  = lcblk(8,iblk)
          nsymdl = lcblk(9,iblk)
          npro   = lcblk(1,iblk+1) - iel
          inum   = iel + npro - 1
          ngauss = nint(lcsyst)
c
c.... compute and assemble the residual and tangent matrix
          ndofelas = nshl * nelas
c
          allocate (Estiff(npro,ndofelas,ndofelas))
          allocate (tmpshp(nshl,MAXQPT))
          allocate (tmpshgl(nsd,nshl,MAXQPT))
c
          Estiff = zero
          tmpshp(1:nshl,:) = shp(lcsyst,1:nshl,:)
          tmpshgl(:,1:nshl,:) = shgl(lcsyst,:,1:nshl,:)
c
c.... Shape measure. Calculate the shape quality
c
          call shpMeasure(x, mien(iblk)%p, tmpshp, tmpshgl,
     &                    meshq(iel:iel+npro-1),
     &                    meshV(iel:iel+npro-1), errorcount )
c
          call AsIGMRElas (x,             disp,
     &                     tmpshp,        tmpshgl,
     &                     mien(iblk)%p,  elasres,
     &                     elasBDiag,     Estiff,
     &                     meshq(iel:iel+npro-1),
     &                     meshV(iel:iel+npro-1)   )
c
c.... satisfy the BCs on the implicit LHS
c
          call bc3LHSElas (iBC, BC(:, ndof+2:ndof+5),
     &                     mien(iblk)%p, Estiff)
c
c.... Fill-up the global sparse LHS mass matrix
c
          call fillsparseElas( mien(iblk)%p, Estiff,
     &                         elaslhsK, row, col)
c
          deallocate ( Estiff )
          deallocate ( tmpshp )
          deallocate ( tmpshgl )
c
c.... end of interior element loop
c
       enddo
c
      if (errorcount(1) .lt. 0 .or. errorcount(2) .lt. 0) then
        write(*,*) errorcount(1), " elements Meshq larger than one; ",
     &             errorcount(2), " elements Meshq smaller than zero."
      endif
c
      if(iabc==1) then               ! are there any axisym BCs
          call rotabc(elasres(1,2), iBC,  'in ')
      endif
c
c.... -------------------->   communications <-------------------------
c
      if (numpe > 1) then
          call commu (elasres  , ilwork, nelas  , 'in ')
c
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
c
          if (iprec .ne. 0) then
             call commu (elasBDiag,    ilwork,
     &                   nelas*nelas,  'in ')
          endif
      endif
c
c.... ---------------------->   post processing  <----------------------
c
c.... satisfy the BCs on the residual
c
      call bc3ResElas (iBC,     BC(:, ndof+2:ndof+5),
     &                 elasres, iper,    ilwork)
c
c.... satisfy the BCs on the block-diagonal preconditioner
c
      if (iprec .ne. 0) then
         call bc3BDgElas (iBC,       BC(:, ndof+2:ndof+5),
     &                    elasBDiag, iper,    ilwork)
      endif
c
c.... return
c
      return
      end

