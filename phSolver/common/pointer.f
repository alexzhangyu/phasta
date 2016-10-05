       module pointer_data
c
c.... maximum number of blocks
c
         parameter ( MAXBLK2 = 50000 ) ! Note compiler was complaining 
c                                       because MAXBLK in common.h be careful
c    					to chang both places
c
c.... data type definitions
c
         type r1d
           real*8, pointer :: p(:)
         end type
c
         type r2d
           real*8, pointer :: p(:,:)
         end type
c
         type r3d
           real*8, pointer :: p(:,:,:)
         end type
c
         type i1d
           integer, pointer :: p(:)
         end type
c
         type i2d
           integer, pointer :: p(:,:)
         end type
c
         type i2d64
           integer*8, pointer :: p(:,:)
         end type
c
         type i3d
           integer, pointer :: p(:,:,:)
         end type
c
c.... pointer declarations
c
         type (i1d), dimension(MAXBLK2) ::  mmat,  mmatb
         type (i2d), dimension(MAXBLK2) ::  mien
         type (i2d64), dimension(MAXBLK2) ::  mienG
         type (i2d), dimension(MAXBLK2) ::  mienb,  miBCB
         type (i2d), dimension(MAXBLK2) ::  mienif0, mienif1
         type (r2d), dimension(MAXBLK2) ::  mxmudmi
         type (r3d), dimension(MAXBLK2) ::  mBCB
         type (r3d), dimension(MAXBLK2) ::  b !for solid,left Cauchy_green tensor,added
         type (r3d), dimension(MAXBLK2) ::  b_dot!time derivative of b,added
         type (r3d), dimension(MAXBLK2) ::  b_af! b at time step n+af,added
         type (r3d), dimension(MAXBLK2) ::  bdy_b !left Cauchy_green tensor on the boundary,added
         type (r3d), dimension(MAXBLK2) ::  bdy_b_dot!time derivative of b on the boudary,added
         type (r3d), dimension(MAXBLK2) ::  bdy_b_af! b at time step n+af on the boudary,added
         
c          
         real*8, allocatable :: gmass(:)
       end module
c 
c
c
