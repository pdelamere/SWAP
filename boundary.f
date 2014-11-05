      module boundary

      USE global
     
      contains


c---------------------------------------------------------------------
      SUBROUTINE periodic(b)
c---------------------------------------------------------------------
CVD$F VECTOR
      !include 'incurv.h'

      real b(nx,ny,nz,3)

c x direction

      do 10 j=1,ny
         do 10 k=1,nz
            do 10 m=1,3
               b(1,j,k,m) = b(nx-1,j,k,m)
               b(nx,j,k,m) = b(2,j,k,m)
c               b(1,j,k,m) = b(2,j,k,m)
cc               b(nx,j,k,m) = b(2,j,k,m)
cc               b(nx-1,j,k,m) = b(nx-2,j,k,m)
c               b(nx,j,k,m) = b(nx-1,j,k,m)
 10            continue


c y direction

      do 20 i=1,nx
         do 20 k=1,nz
            do 20 m=1,3
               b(i,1,k,m) = b(i,ny-1,k,m)
               b(i,ny,k,m) = b(i,2,k,m)
 20            continue

c z direction

      do 30 i=1,nx
         do 30 j=1,ny
            do 30 m=1,3
               b(i,j,1,m) = b(i,j,nz-1,m)
               b(i,j,nz,m) = b(i,j,2,m)
 30            continue




      return
      end SUBROUTINE periodic
c---------------------------------------------------------------------


c---------------------------------------------------------------------
      SUBROUTINE tangent_B_zero(b) !normal derivative = 0
c The normal derivative of the tangential components is used to 
c determine the tangential boundary values.  ALSO, the normal
c derivative of the normal components temporarily sets the boundary
c values of the normal components.  These values are undated later
c by the requirement that divB=0.  This helps with the values on
c corners and edges.
c---------------------------------------------------------------------
      !include 'incurv.h'

      real b(nx,ny,nz,3)

c x surfaces
      do 10 j=2,ny-1
         do 10 k=2,nz-1

            b(1,j,k,1) = b(2,j,k,1)       !normal
c            b(1,j,k,1) = b(3,j,k,1)
c            b(nx,j,k,1) = b(nx-1,j,k,1)

            b(1,j,k,2) = b(2,j,k,2)       !tangential
            b(1,j,k,3) = b(2,j,k,3)
            b(nx,j,k,2) = b(nx-1,j,k,2)
            b(nx,j,k,3) = b(nx-1,j,k,3)


 10         continue

c y surfaces
       do 20 i=2,nx-1
          do 20 k=2,nz-1

             b(i,1,k,2) = b(i,2,k,2)      !normal
c             b(i,1,k,2) = b(i,3,k,2)
c             b(i,ny,k,2) = b(i,ny-1,k,2)

             b(i,1,k,1) = b(i,2,k,1)      !tangential
             b(i,1,k,3) = b(i,2,k,3)
             b(i,ny,k,1) = b(i,ny-1,k,1)
             b(i,ny,k,3) = b(i,ny-1,k,3)

 20          continue

c z surfaces
       do 30 i=2,nx-1
          do 30 j=2,ny-1

             b(i,j,1,3) = b(i,j,2,3)      !normal
c             b(i,j,1,3) = b(i,j,3,3)
c             b(i,j,nz,3) = b(i,j,nz-1,3)

             b(i,j,1,1) = b(i,j,2,1)      !tangential
             b(i,j,1,2) = b(i,j,2,2)
             b(i,j,nz,1) = b(i,j,nz-1,1)
             b(i,j,nz,2) = b(i,j,nz-1,2)
             
 30          continue


      return
      end SUBROUTINE tangent_B_zero
c---------------------------------------------------------------------


c---------------------------------------------------------------------
      SUBROUTINE copy_to_boundary(b)
c---------------------------------------------------------------------
      !include 'incurv.h'

      real b(nx,ny,nz,3)

c x surfaces      !periodic
      do 10 j=1,ny
         do 10 k=1,nz
c            do 10 m=1,3
c               b(2,j,k,1) = b(nx-1,j,k,1)       !normal
               b(1,j,k,1) = b(nx-1,j,k,1)       !normal
               b(nx,j,k,1) = b(2,j,k,1)

               b(1,j,k,2) = b(nx-1,j,k,2)       !tangential
c               b(2,j,k,2) = b(nx-1,j,k,2)       !tangential
               b(nx,j,k,2) = b(2,j,k,2)

               b(1,j,k,3) = b(nx-1,j,k,3)       !tangential
c               b(2,j,k,2) = b(nx-1,j,k,2)       !tangential
               b(nx,j,k,3) = b(2,j,k,3)
 10         continue

c y surfaces
       do 20 i=1,nx
          do 20 k=1,nz
             do 20 m=1,3
                b(i,1,k,m) = b(i,2,k,m)      !tangential
                b(i,ny,k,m) = b(i,ny-1,k,m)
 20          continue

c z surfaces
       do 30 i=1,nx
          do 30 j=1,ny
             do 30 m=1,3
                b(i,j,1,m) = b(i,j,2,m)      !tangential
                b(i,j,nz,m) = b(i,j,nz-1,m)
 30          continue

      return
      end SUBROUTINE copy_to_boundary
c---------------------------------------------------------------------


c---------------------------------------------------------------------
      SUBROUTINE periodic_scalar(b)
c---------------------------------------------------------------------
      !include 'incurv.h'

      real b(nx,ny,nz)


c x surfaces      !periodic
      do 10 j=1,ny
         do 10 k=1,nz
c            do 10 m=1,3
c               b(1,j,k) = b(2,j,k)       !tangential
c               b(nx,j,k) = b(nx-1,j,k)
               b(1,j,k) = b(nx-1,j,k)       !tangential
               b(nx,j,k) = b(2,j,k)
 10         continue



c y surfaces
       do 20 i=1,nx
          do 20 k=1,nz
c             do 20 m=1,3
                b(i,1,k) = b(i,ny-1,k)      !tangential
                b(i,ny,k) = b(i,2,k)
 20          continue

c z surfaces
       do 30 i=1,nx
          do 30 j=1,ny
c             do 30 m=1,3
                b(i,j,1) = b(i,j,nz-1)      !tangential
                b(i,j,nz) = b(i,j,2)
 30          continue




      return
      end SUBROUTINE periodic_scalar
c---------------------------------------------------------------------



c---------------------------------------------------------------------
      SUBROUTINE fix_normal_b(b)
c---------------------------------------------------------------------
      !include 'incurv.h'

      real b(nx,ny,nz,3)

c normal x components
      do 10 j=2,ny-1
         do 10 k=2,nz-1
            b(2,j,k,1) = dx*(b(2,j+1,k,2) - b(2,j,k,2))/dy +
     x                    dx*(b(2,j,k+1,3) - b(2,j,k,3))/dz_grid(k) +
     x                    b(3,j,k,1)

            b(nx-1,j,k,1) = b(nx-2,j,k,1) -
     x               dx*(b(nx-2,j+1,k,2) - b(nx-2,j,k,2))/dy -
     x               dx*(b(nx-2,j,k+1,3) - b(nx-2,j,k,3))/dz_grid(k)

 10         continue

c normal y components
c      do 20 i=2,nx-1
c         do 20 k=2,nz-1
c            b(i,2,k,2) = dy*(b(i+1,2,k,1) - b(i,2,k,1))/dx + 
c     x                    dy*(b(i,2,k+1,3) - b(i,2,k,3))/dz_grid(k) + 
c     x                    b(i,3,k,2)

c            b(i,ny,k,2) = b(i,ny-1,k,2) -
c     x                  dy*(b(i+1,ny-1,k,1) - b(i,ny-1,k,1))/dx -
c     x                  dy*(b(i,ny-1,k+1,3) - b(i,ny-1,k,3))/dz_grid(k)
c 20         continue

c normal z components
c      do 30 i=2,nx-1
c         do 30 j=2,ny-1
c            b(i,j,2,3) = dz_grid(2)*(b(i+1,j,2,1) - b(i,j,2,1))/dx + 
c     x                   dz_grid(2)*(b(i,j+1,2,2) - b(i,j,2,2))/dy +
c     x                   b(i,j,3,3)

c            b(i,j,nz,3) = b(i,j,nz-1,3) -
c     x               dz_grid(nz)*(b(i+1,j,nz-1,1) - b(i,j,nz-1,1))/dx +
c     x               dz_grid(nz)*(b(i,j+1,nz-1,2) - b(i,j,nz-1,2))/dy

c 30         continue


      return
      end SUBROUTINE fix_normal_b
c---------------------------------------------------------------------

c---------------------------------------------------------------------
      SUBROUTINE smooth_boundary(b)
c---------------------------------------------------------------------
      !include 'incurv.h'

      real b(nx,ny,nz,3)

c x surfaces
      do 10 j=2,ny-1
         do 10 k=2,nz-1
            do 10 m=1,3
               b(1,j,k,m) = b(2,j,k,m)       !tangential
               b(nx,j,k,m) = b(nx-1,j,k,m)
 10         continue

c y surfaces
       do 20 i=1,nx
          do 20 k=1,nz
             do 20 m=1,3
                b(i,1,k,m) = b(i,2,k,m)      !tangential
                b(i,ny,k,m) = b(i,ny-1,k,m)
 20          continue

c z surfaces
       do 30 i=1,nx
          do 30 j=1,ny
             do 30 m=1,3
                b(i,j,1,m) = b(i,j,2,m)      !tangential
                b(i,j,nz,m) = b(i,j,nz-1,m)
 30          continue

      return
      end SUBROUTINE smooth_boundary
c---------------------------------------------------------------------


c---------------------------------------------------------------------
      SUBROUTINE fix_tangential_E(E)
c---------------------------------------------------------------------
      !include 'incurv.h'

      real E(nx,ny,nz,3)

c     i = 2 & i = nx
      do 10 j=2,ny     !periodic boundary conditions
         do 10 k=2,nz
cc            E(2,j,k,1) = E(nx-1,j,k,1)  !normal component
cc            E(2,j,k,2) = E(nx-1,j,k,2)
cc            E(2,j,k,3) = E(nx-1,j,k,3)
c            E(nx,j,k,1) = E(3,j,k,1)  !normal component
c            E(nx,j,k,2) = E(3,j,k,2)
c            E(nx,j,k,3) = E(3,j,k,3)
c            E(nx,j,k,1) =   !normal component
            E(nx,j,k,2) = -2.3
            E(nx,j,k,3) = 0.0

 10         continue

c      write(*,*) 'E bnd...',E(23,8,14,1),E(23,8,14,2),E(23,8,14,3)

cc     j = 2 & j = ny
c      do 20 i=2,nx
c         do 20 k=2,nz
c            E(i,2,k,1) = E(i,3,k,1)
cc            E(i,2,k,2) = E(i,3,k,2)   !normal component
c            E(i,2,k,3) = E(i,3,k,3)
c            E(i,ny,k,1) = E(i,ny-1,k,1)
c            E(i,ny,k,2) = E(i,ny-1,k,2)  !normal component
c            E(i,ny,k,3) = E(i,ny-1,k,3)
c 20         continue

cc     k = 2 & k = nz
c      do 30 i=2,nx-1
c         do 30 j=2,ny-1
cc            E(i,j,2,1) = E(i,j,3,1)
cc            E(i,j,2,2) = E(i,j,3,2)
cc            E(i,j,nz,1) = E(i,j,nz-1,1)
cc            E(i,j,nz,2) = E(i,j,nz-1,2)
c            E(i,j,2,1) = 0.0
c            E(i,j,2,2) = 0.0
cc            E(i,j,2,3) = 0.0   !normal component
c            E(i,j,nz-1,1) = 0.0
c            E(i,j,nz-1,2) = 0.0
cc            E(i,j,nz,3) = 0.0  !normal component
c 30         continue

c      call periodic(E)

      return
      end SUBROUTINE fix_tangential_E
c---------------------------------------------------------------------

      end module boundary



















