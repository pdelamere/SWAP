      MODULE initial

      USE global
c      USE dimensions
c      USE inputs

      contains

c----------------------------------------------------------------------
      SUBROUTINE grd6_setup(b0,bt,b12,b1,b1p2,nu)
c----------------------------------------------------------------------

      real eoverm
c      parameter(eoverm = q/mO)
      real mO_q
c      parameter(mO_q = mO/q)
      real vol


      real b0r,a1,a2,omegar,kr,nfr,kdz
      real b0_1x, b0_2x, b0_1y, b0_2y
      real phi



      real b0(nx,ny,nz,3),
     x     bt(nx,ny,nz,3),
     x     b12(nx,ny,nz,3),
     x     b1(nx,ny,nz,3),
     x     b1p2(nx,ny,nz,3),
     x     nu(nx,ny,nz)


c      include 'incurv.h'

      eoverm = q/mO
      mO_q = mO/q

      phi = 2.0*PI/180.

      b0_1x = b0_top*eoverm*sin(phi)
      b0_2x = b0_bottom*eoverm*sin(phi)

      b0_1y = -b0_top*eoverm*cos(phi)
      b0_2y = -b0_bottom*eoverm*cos(phi)  

c      kr = 2.0/delz

c      a1 = kr**2*b0r/(alpha*nfr)
c      a2 = kr**2*b0r**2/(alpha*nfr)

c      omegar = 0.5*(a1 + sqrt(a1**2 + 4*a2))

      do i = 1,nx
         do j = 1,ny 
            do k = 1,nz
c               b0(i,j,k,1) = 0.5*(b0_1x+b0_2x) + 
c     x              0.5*(b0_1x-b0_2x)*tanh((qz(k)-qz(nz/2))/(Lo))
c               b0(i,j,k,2) = 0.5*(b0_1y+b0_2y) + 
c     x              0.5*(b0_1y-b0_2y)*tanh((qz(k)-qz(nz/2))/(Lo))
               b0(i,j,k,1) = 0.0
               b0(i,j,k,2) = 0.0
               b0(i,j,k,3) = b0_init*eoverm
            enddo
         enddo
      enddo
      
c      do 10 k=1,nz
c         kdz = 2.0/dz_grid(k)
cc         b0(k) = b0r
cc         b0(k) = 1.2*(0.5*(-omegar + (omegar/kdz)*
cc     x                 sqrt(kdz**2 + 4*alpha*nfr)))
c 10      continue

      do 20 i=1,nx
         do 20 j=1,ny
            do 20 k=1,nz
c               nf(i,j,k) = nfr  !/km^3
               nu(i,j,k) = nu_init
               do 20 m=1,3

                  bt(i,j,k,m) = b0(i,j,k,m)
                  b12(i,j,k,m) = b0(i,j,k,m)
                  b1(i,j,k,m) = b0(i,j,k,m)
                  b1p2(i,j,k,m) = b0(i,j,k,m)
                  vol = dx_cell(i)*dy_cell(j)*dz_cell(k)*km_to_m**3
                  input_Eb = input_Eb + 
     x                      (vol/(2.0*mu0))*(mO_q*b0(i,j,k,m))**2 

 20            continue

c      nu(1,:,:) = nu_init*50.0
c      nu(2,:,:) = nu_init*50.0
c      nu(3,:,:) = nu_init*10.0
c      nu(4,:,:) = nu_init*5.0
c      nu(5,:,:) = nu_init*2.0

c      nu(nx,:,:) = nu_init*50.0
c      nu(nx-1,:,:) = nu_init*50.0
c      nu(nx-2,:,:) = nu_init*10.0
c      nu(nx-3,:,:) = nu_init*5.0
c      nu(nx-4,:,:) = nu_init*2.0


c set resistive boundary 

c      do 50 i=1,nx
c         do 50 k=1,nz
c            nu(i,4,k) = 5.0
c            nu(i,3,k) = 10.0
c            nu(i,2,k) = 20.0
c            nu(i,1,k) = 20.0
c            nu(i,ny,k) = 20.0
c            nu(i,ny-1,k) = 10.0
c            nu(i,ny-2,k) = 5.0
c 50         continue

c      do 60 j=1,ny
c         do 60 k=1,nz
c            nu(4,j,k) = 5.0
c            nu(3,j,k) = 10.0
c            nu(2,j,k) = 20.0
c            nu(1,j,k) = 20.0
c            nu(nx,j,k) = 20.0
c            nu(nx-1,j,k) = 10.0
c            nu(nx-2,j,k) = 5.0
c 60         continue



c      open(30,file='nf.dat',status='unknown',form='unformatted')
c      write(30) nz
c      write(30) nf
c      close(30)

      open(40,file='b0.dat',status='unknown',form='unformatted')
      write(40) nz
      write(40) b0
      close(40)

      return
      end subroutine grd6_setup
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE grd7()
c----------------------------------------------------------------------
c      include 'incurv.h'

      parameter(nrgrd = 20)
      real zplus,zminus,xplus, xminus, yplus, yminus
      real zsf,xsf

      rk=(nz/2)! - 35
      rj=ny/2
      ri=nx/2



      do 20 j=1,ny
         qy(j) = j*dy
         dy_grid(j) = dy
 20            continue



c      do 10 i=1,nx
c         qx(i) = i*dx
c         dx_grid(i) = dx
c 10            continue

c==============stretch x direction=====================================
               
      xsf = 0.0
c up from center
      do 12 i = ri,ri+nrgrd
         dx_grid(i) = dx
 12   continue
      do 14 i = ri+nrgrd+1,nx
         dx_grid(i) = dx +
     x     xsf*dx*(i-(ri+nrgrd+1))/(nx-(ri+nrgrd+1)) 
 14   continue

c down from center
      do 16 i = ri-nrgrd,ri-1
         dx_grid(i) = dx
 16      continue
      do 17 i = 1,ri-nrgrd-1
         ind = ri-nrgrd-i
         dx_grid(ind) = dx + 
     x     xsf*dx*(ri-nrgrd-1-ind)/(ri-nrgrd-1)
 17   continue

      qx(1) = dx
      do 18 i=2,nx
         qx(i) = qx(i-1)+dx_grid(i)
 18   continue

      do 19 i = 1,nx-1
         dx_grid(i) = qx(i+1)-qx(i)
 19   continue
      dx_grid(nx) = dx_grid(nx-1)
c======================================================================

c      print*,'dx...',dx_grid


c==============stretch z direction=====================================

      zsf = 0.0  !z stretch factor
c up from center
      do 32 k = rk,rk+nrgrd
         dz_grid(k) = delz
 32   continue
      do 34 k = rk+nrgrd+1,nz
         dz_grid(k) = delz +
     x     zsf*delz*(k-(rk+nrgrd+1))/(nz-(rk+nrgrd+1)) 
c     x     2.0*sin((k-(rk+nrgrd+1))*0.5*pi/(nz-(rk+nrgrd+1)))**2 
c                                !dz_grid(k-1) + 0.01*delz 
 34   continue

c down from center
      do 36 k = rk-nrgrd,rk-1
         dz_grid(k) = delz
 36      continue
      do 37 k = 1,rk-nrgrd-1
         ind = rk-nrgrd-k
         dz_grid(ind) = delz + 
     x     zsf*delz*(rk-nrgrd-1-ind)/(rk-nrgrd-1)
c     x     2.0*sin((rk-nrgrd-1-ind)*(-0.5*pi)/(rk-nrgrd-1))**2 
c                                !dz_grid(ind+1) + 0.01*delz
 37   continue

      qz(1) = delz
      do 38 k=2,nz
c         write(*,*) 'dz_grid...',k,dz_grid(k)
         qz(k) = qz(k-1)+dz_grid(k)
 38   continue

      do 39 k = 1,nz-1
         dz_grid(k) = qz(k+1)-qz(k)
 39   continue
      dz_grid(nz) = dz_grid(nz-1)
c======================================================================


      dz_cell(1) = dz_grid(1)
      dz_cell(nz) = dz_grid(nz)
      zrat(1) = 0.5
      zrat(nz) = 0.5
      do 40 k=2,nz-1
         dz_cell(k) = ((qz(k+1) + qz(k))/2.0) -
     x                ((qz(k) + qz(k-1))/2.0)
               
         zplus = (qz(k+1) + qz(k))/2.0
         zminus = (qz(k) + qz(k-1))/2.0
         zrat(k) = (qz(k) - zminus)/(zplus - zminus)
 40   continue


      dx_cell(1) = dx_grid(1)
      dx_cell(nx) = dx_grid(nx)
      xrat(1) = 0.5
      xrat(nx) = 0.5
      do 50 i=2,nx-1
         dx_cell(i) = ((qx(i+1) + qx(i))/2.0) -
     x                ((qx(i) + qx(i-1))/2.0)
         xplus = (qx(i+1) + qx(i))/2.0
         xminus = (qx(i) + qx(i-1))/2.0
         xrat(i) = (qx(i) - xminus)/(xplus - xminus)
 50   continue

      dy_cell(1) = dy_grid(1)
      dy_cell(ny) = dy_grid(ny)
      yrat(1) = 0.5
      yrat(ny) = 0.5
      do 60 j=2,ny-1
         dy_cell(j) = ((qy(j+1) + qy(j))/2.0) -
     x                ((qy(j) + qy(j-1))/2.0)
         yplus = (qy(j+1) + qy(j))/2.0
         yminus = (qy(j) + qy(j-1))/2.0
         yrat(j) = (qy(j) - yminus)/(yplus - yminus)
 60   continue

c      call assign('assign -F system -N ultrix f:' //'c.coord.dat')
      open(40,file=trim(out_dir)//'c.coord.dat',status='unknown',
     x         form='unformatted')

c      open(40,file='coord.dat',status='unknown',form='unformatted')
      write(40) nx
      write(40) ny
      write(40) nz
      write(40) qx
      write(40) qy
      write(40) qz
      write(40) dz_grid
      write(40) dz_cell
      close(40)

      return
      end subroutine grd7
c----------------------------------------------------------------------



      end module initial


