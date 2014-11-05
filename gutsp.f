      MODULE gutsp

      USE global
      USE misc
      USE mpi
      USE boundary
      USE grid_interp

      contains

c----------------------------------------------------------------------
      SUBROUTINE remove_ion(xp,vp,vp1,ion_l)
c Removes particles from simulation that have gone out of bounds
c----------------------------------------------------------------------
CVD$R VECTOR


      !include 'incurv.h'

      real xp(Ni_max,3)
      real vp(Ni_max,3)
      real vp1(Ni_max,3)

      do 5 m=1,3   !remove ion energy from total input energy
         input_E = input_E-0.5*m_arr(l)*(vp(ion_l,m)*km_to_m)**2 /beta
 5    continue

      do 10 l=ion_l,Ni_tot-1
         do 10 m=1,3 
            xp(l,m) = xp(l+1,m)
            vp(l,m) = vp(l+1,m)
            vp1(l,m) = vp1(l+1,m)
            ijkp(l,m) = ijkp(l+1,m)
c            wquad(l,m) = wquad(l+1,m)
 10      continue

      do 20 m=1,8
         do 20 l=ion_l,Ni_tot-1
            wght(l,m) = wght(l+1,m)
 20         continue

      Ni_tot = Ni_tot - 1

      return
      end SUBROUTINE remove_ion
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE check_min_den(np,xp,vp,vp1,up,bt)
c----------------------------------------------------------------------
      !!include 'incurv.h'

      real np(nx,ny,nz),
     x     xp(Ni_max,3),
     x     vp(Ni_max,3),
     x     vp1(Ni_max,3),
     x     up(nx,ny,nz,3),
     x     bt(nx,ny,nz,3)

      real minden           !minimum wake density 
      real den_part ! density of 1 particle per cell
      real ak
      real btot,a1,a2,womega,phi,deltat
      integer npart,ipart
      

      den_part = 1/(beta*dx**3)

c      minden = nf_init/10.
      minden = 2.0*den_part
      do i = 2,nx-1
         do j = 2,ny-1
            do k = 2,nz-1
               ak = PI/dx
               btot = sqrt(bt(i,j,k,1)**2 + bt(i,j,k,2)**2 + 
     x              bt(i,j,k,3)**2)
               a1 = ak**2*Btot/(alpha*(np(i,j,k)))
               a2 = (ak*Btot)**2/(alpha*(np(i,j,k)))
               womega = 0.5*(a1 + sqrt(a1**2 + 4*a2))
               phi = womega/ak
               deltat = 0.1*dx/phi
c               if (deltat .le. dtsub) then
c                  write(*,*) 'Time stepping error...'
c               endif
               if (np(i,j,k).le.minden) then
                  npart = nint(minden/np(i,j,k))
                  do ipart = 1,npart 

                     l=Ni_tot + 1 !beginning array element for new borns    
                  
                     vp(l,1) = up(i,j,k,1)
                     vp(l,2) = up(i,j,k,2)
                     vp(l,3) = up(i,j,k,3)
                     xp(l,1) = qx(i) + (0.5-pad_ranf())*dx
                     xp(l,2) = qy(j) + (0.5-pad_ranf())*dy
                     xp(l,3) = qz(k) + (0.5-pad_ranf())*dz_grid(k)

                     ijkp(l,1) = nint(xp(l,1)/dx) !particle grid location index
                     ijkp(l,2) = nint(xp(l,2)/dy)

                     kk=1
                     do 100 while(xp(l,3) .gt. qz(kk)) !find k on non-uniform 
                        ijkp(l,3) = kk !grid
                        kk=kk+1
 100                 continue
                     kk=ijkp(l,3)
                     if (xp(l,3) .gt. (qz(kk)+(dz_grid(kk)/2))) then
                        ijkp(l,3) = kk+1
                     endif
                     mrat(l) = 1.0
                     m_arr(l) = mion
                     Ni_tot = Ni_tot + 1
                  enddo
               endif
            enddo
         enddo
      enddo

      return
      end SUBROUTINE check_min_den
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE extrapol_up(up,vp,vp1,np)
c This subroutine does the provisional extrapolation of the particle
c bulk flow velocity to time level n, and replaces up_n-3/2 
c with up_n-1/2
c----------------------------------------------------------------------
CVD$R VECTOR
      !!include 'incurv.h'

      real up(nx,ny,nz,3),
     x     vp(Ni_max,3),
     x     vp1(Ni_max,3),
     x     np(nx,ny,nz)
      
      real v_at_n(Ni_max,3)

      do 10 m=1,3
         do 10 l=1,Ni_tot
            v_at_n(l,m) = 1.5*vp(l,m) - 0.5*vp1(l,m)
 10      continue

      call update_up(v_at_n,np,up)

      return
      end SUBROUTINE extrapol_up
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE get_Ep(Ep,aj,np,up,btc,nu)
c----------------------------------------------------------------------
CVD$F VECTOR
      !!include 'incurv.h'

      real Ep(Ni_max,3),
     x     aj(nx,ny,nz,3),
     x     np(nx,ny,nz),
c     x     nf(nx,ny,nz),
     x     up(nx,ny,nz,3),
c     x     uf(nx,ny,nz,3),
     x     btc(nx,ny,nz,3),
     x     nu(nx,ny,nz)
c     x     gradP(nx,ny,nz,3)

      real ajc(nx,ny,nz,3),       !aj at cell center
     x     upc(nx,ny,nz,3),       !up at cell center
     x     ufc(nx,ny,nz,3)       !uf at cell center
c     x     gradPc(nx,ny,nz,3)     !gradP at cell center

      real aa(3), bb(3), cc(3)    !dummy vars for doing cross product
      real ntot                   !total plasma density
      real fnf,fnp                !fraction nf,np of total n
      real aj3(3),up3(3),uf3(3),  !vector field values at Ba position
     x     btc3(3),gradP3(3)
      real np_at_Ba               !particle density at particle pos
      real nf_at_Ba

      call face_to_center(aj,ajc)
      call face_to_center(up,upc)

c      call face_to_center(uf,ufc)
c      call face_to_center(gradP,gradPc)

      do 10 l=1,Ni_tot

         i=ijkp(l,1)!+wquad(l,1)
         j=ijkp(l,2)!+wquad(l,2)
         k=ijkp(l,3)!+wquad(l,3)

         ip = i+1
         jp = j+1
         kp = k+1

c         if (ip .ge. nx) ip = 2      !periodic boundary conditions
c         if (jp .ge. ny) jp = 2      !periodic boundary conditions
c         if (kp .ge. nz) kp = 2      !periodic boundary conditions

c         np_at_Ba = np(i,j,k)*wght(l,1) + np(ip,j,k)*wght(l,2) + 
c     x              np(i,j,kp)*wght(l,3) + np(ip,j,kp)*wght(l,4) + 
c     x              np(i,jp,k)*wght(l,5) + np(ip,jp,k)*wght(l,6) +
c     x              np(i,jp,kp)*wght(l,7) + np(ip,jp,kp)*wght(l,8)

c         nf_at_Ba = nf(i,j,k)*wght(l,1) + nf(ip,j,k)*wght(l,2) + 
c     x              nf(i,j,kp)*wght(l,3) + nf(ip,j,kp)*wght(l,4) + 
c     x              nf(i,jp,k)*wght(l,5) + nf(ip,jp,k)*wght(l,6) +
c     x              nf(i,jp,kp)*wght(l,7) + nf(ip,jp,kp)*wght(l,8)

c         ntot = nf_at_Ba + np_at_Ba
c         fnf = nf_at_Ba/ntot
c         fnp = np_at_Ba/ntot

         do 15 m=1,3 
            aj3(m) = ajc(i,j,k,m)*wght(l,1) + ajc(ip,j,k,m)*wght(l,2) 
     x          + ajc(i,j,kp,m)*wght(l,3) + ajc(ip,j,kp,m)*wght(l,4)
     x          + ajc(i,jp,k,m)*wght(l,5) + ajc(ip,jp,k,m)*wght(l,6)
     x          + ajc(i,jp,kp,m)*wght(l,7) + ajc(ip,jp,kp,m)*wght(l,8)

            up3(m) = upc(i,j,k,m)*wght(l,1) + upc(ip,j,k,m)*wght(l,2) 
     x          + upc(i,j,kp,m)*wght(l,3) + upc(ip,j,kp,m)*wght(l,4)
     x          + upc(i,jp,k,m)*wght(l,5) + upc(ip,jp,k,m)*wght(l,6)
     x          + upc(i,jp,kp,m)*wght(l,7) + upc(ip,jp,kp,m)*wght(l,8)

c            uf3(m) = ufc(i,j,k,m)*wght(l,1) + ufc(ip,j,k,m)*wght(l,2) 
c     x          + ufc(i,j,kp,m)*wght(l,3) + ufc(ip,j,kp,m)*wght(l,4)
c     x          + ufc(i,jp,k,m)*wght(l,5) + ufc(ip,jp,k,m)*wght(l,6)
c     x          + ufc(i,jp,kp,m)*wght(l,7) + ufc(ip,jp,kp,m)*wght(l,8)

            btc3(m) = btc(i,j,k,m)*wght(l,1) 
     x               + btc(ip,j,k,m)*wght(l,2) 
     x               + btc(i,j,kp,m)*wght(l,3) 
     x               + btc(ip,j,kp,m)*wght(l,4)
     x               + btc(i,jp,k,m)*wght(l,5) 
     x               + btc(ip,jp,k,m)*wght(l,6)
     x               + btc(i,jp,kp,m)*wght(l,7) 
     x               + btc(ip,jp,kp,m)*wght(l,8)

c            gradP3(m) = gradPc(i,j,k,m)*wght(l,1) 
c     x                + gradPc(ip,j,k,m)*wght(l,2) 
c     x                + gradPc(i,j,kp,m)*wght(l,3) 
c     x                + gradPc(ip,j,kp,m)*wght(l,4)
c     x                + gradPc(i,jp,k,m)*wght(l,5) 
c     x                + gradPc(ip,jp,k,m)*wght(l,6)
c     x                + gradPc(i,jp,kp,m)*wght(l,7) 
c     x                + gradPc(ip,jp,kp,m)*wght(l,8)

 15         continue

         do 20 m=1,3
c            aa(m) = aj3(m) - fnp*up3(m)  
c     x                     - fnf*uf3(m)
            aa(m) = aj3(m) - up3(m)
            bb(m) = btc3(m)                   
 20         continue

         cc(1) = aa(2)*bb(3) - aa(3)*bb(2)    !do cross product
         cc(2) = aa(3)*bb(1) - aa(1)*bb(3)
         cc(3) = aa(1)*bb(2) - aa(2)*bb(1)

         do 30 m=1,3
            Ep(l,m) = cc(m) 
c     x                + nu(i,j,k)*fnf*(uf3(m)-up3(m)) 
c     x                - gradP3(m) 
c     x                + nuei*aj3(m) 
c                     + etar(i,j,k,m)*aj3(m)
            Ep(l,m) = Ep(l,m)*mrat(l) 
 30      continue

 10      continue

      return
      end SUBROUTINE get_Ep
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE get_vplus_vminus(Ep,btc,vp,vplus,vminus)
c----------------------------------------------------------------------

      !!include 'incurv.h'

      real Ep(Ni_max,3),
     x     btc(nx,ny,nz,3),   !bt at cell center
     x     vp(Ni_max,3),      !particle velocities at t level n-1/2
     x     vplus(Ni_max,3),
     x     vminus(Ni_max,3)

      real a1, a2, a3,       !coefficients for calculating vplus
     x     a_d,              !denominator for a1,a2,a3
     x     B2,dt2,           !B*B,dt*dt
     x     Bx,By,Bz,         !B for cross product call
     x     vminus_x_B(3),    !v- x B
     x     vminus_dot_B     !v- . B

      real btc3(3)
      

      do 10 m=1,3
         do 10 l=1,Ni_tot 
            vminus(l,m) = vp(l,m) + 0.5*dt*Ep(l,m)
 10         continue

c      do 20 l=1,Ni_tot

c         i = ijkp(l,1)+wquad(l,1)
c         j = ijkp(l,2)+wquad(l,2)
c         k = ijkp(l,3)+wquad(l,3)
         
c         ip=i+1
c         jp=j+1
c         kp=k+1

c         if (ip .ge. nx) ip = 2      !periodic boundary conditions
c         if (jp .ge. ny) jp = 2      !periodic boundary conditions
c         if (kp .ge. nz) kp = 2      !periodic boundary conditions

c         do 25 m=1,3    !interpolate to particle position
c            btc3(l,m) = btc(i,j,k,m)*wght(l,1) 
c     x               + btc(ip,j,k,m)*wght(l,2) 
c     x               + btc(i,j,kp,m)*wght(l,3) 
c     x               + btc(ip,j,kp,m)*wght(l,4)
c     x               + btc(i,jp,k,m)*wght(l,5) 
c     x               + btc(ip,jp,k,m)*wght(l,6)
c     x               + btc(i,jp,kp,m)*wght(l,7) 
c     x               + btc(ip,jp,kp,m)*wght(l,8)
c 25         continue

c         vminus_x_B(l,1) = vminus(l,2)*btc3(l,3)*mrat(l) - 
c     x                     vminus(l,3)*btc3(l,2)*mrat(l)   
c         vminus_x_B(l,2) = vminus(l,3)*btc3(l,1)*mrat(l) - 
c     x                     vminus(l,1)*btc3(l,3)*mrat(l)   
c         vminus_x_B(l,3) = vminus(l,1)*btc3(l,2)*mrat(l) - 
c     x                     vminus(l,2)*btc3(l,1)*mrat(l)   

c         vminus_dot_B(l) = vminus(l,1)*btc3(l,1)*mrat(l) + 
c     x                     vminus(l,2)*btc3(l,2)*mrat(l) + 
c     x                     vminus(l,3)*btc3(l,3)*mrat(l)   

c 20   continue
   
      do 30 l=1,Ni_tot 

         i = ijkp(l,1)!+wquad(l,1)
         j = ijkp(l,2)!+wquad(l,2)
         k = ijkp(l,3)!+wquad(l,3)
   
         ip = i+1
         jp = j+1
         kp = k+1

c         if (ip .ge. nx) ip = 2 !periodic boundary conditions
c         if (jp .ge. ny) jp = 2 !periodic boundary conditions
c         if (kp .ge. nz) kp = 2 !periodic boundary conditions

         do 35 m=1,3
            btc3(m) = btc(i,j,k,m)*wght(l,1) 
     x           + btc(ip,j,k,m)*wght(l,2) 
     x           + btc(i,j,kp,m)*wght(l,3) 
     x           + btc(ip,j,kp,m)*wght(l,4)
     x           + btc(i,jp,k,m)*wght(l,5) 
     x           + btc(ip,jp,k,m)*wght(l,6)
     x           + btc(i,jp,kp,m)*wght(l,7) 
     x           + btc(ip,jp,kp,m)*wght(l,8)
            
 35      continue

         vminus_x_B(1) = vminus(l,2)*btc3(3)*mrat(l) - 
     x                     vminus(l,3)*btc3(2)*mrat(l)   
         vminus_x_B(2) = vminus(l,3)*btc3(1)*mrat(l) - 
     x                     vminus(l,1)*btc3(3)*mrat(l)   
         vminus_x_B(3) = vminus(l,1)*btc3(2)*mrat(l) - 
     x                     vminus(l,2)*btc3(1)*mrat(l)   

         vminus_dot_B = vminus(l,1)*btc3(1)*mrat(l) + 
     x                     vminus(l,2)*btc3(2)*mrat(l) + 
     x                     vminus(l,3)*btc3(3)*mrat(l)   

         Bx = btc3(1)*mrat(l) 
         By = btc3(2)*mrat(l) 
         Bz = btc3(3)*mrat(l) 
      
         B2 = Bx*Bx + By*By + Bz*Bz
         dt2 = dt*dt

         a_d = 1.0/(1.0 + (0.25*B2*dt2))
         a1 = (1 - (0.25*B2*dt2))*a_d
         a2 = dt*a_d
         a3 = 0.5*dt2*a_d

         do 40 m=1,3
            vplus(l,m) = a1*vminus(l,m) + a2*vminus_x_B(m) + 
     x           a3*vminus_dot_B*btc3(m)*mrat(l) 
 40      continue

 30   continue

      return
      end SUBROUTINE get_vplus_vminus
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE improve_up(vp1,vplus,vminus,up,np)
c The routine calculates v at time level n, and the associated bulk
c flow velocity up using the v+, v- technique.  The new up at
c time level n replaces the provisional extrapolation for up.
c----------------------------------------------------------------------
      !!include 'incurv.h'

      real vp1(Ni_max,3),    !particle velocities at t level n
     x     vplus(Ni_max,3),
     x     vminus(Ni_max,3),
     x     up(nx,ny,nz,3),
     x     np(nx,ny,nz)

      do 10 m=1,3
         do 10 l = 1,Ni_tot
            vp1(l,m) = 0.5*(vplus(l,m) + vminus(l,m))
 10         continue

      call update_up(vp1,np,up)

      return
      end SUBROUTINE improve_up
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE get_vp_final(Ep,vp,vp1,vplus)
c----------------------------------------------------------------------
      !!include 'incurv.h'

      real Ep(Ni_max,3),
     x     vp(Ni_max,3),    !particle velocities at t level n+1/2
     x     vp1(Ni_max,3),   !particle velocity at t level n-1/2
     x     vplus(Ni_max,3)

      do 10 m=1,3
         do 10 l = 1,Ni_tot
            vp1(l,m) = vp(l,m)  !to be used in extrapol_up for n-3/2
            vp(l,m) = vplus(l,m) + 0.5*dt*Ep(l,m)
 10         continue

      return
      end SUBROUTINE get_vp_final
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE move_ion_half(xp,vp,vp1,input_p)
c----------------------------------------------------------------------
      !!include 'incurv.h'

      real xp(Ni_max,3),
     x     vp(Ni_max,3),
     x     vp1(Ni_max,3),
     x     input_p(3)
      integer zindx(1)
      integer Ni_tot_in
      
      real dth                 !half time step

      real v2

      dth = dt/2.

      do 10 l=1,Ni_tot                   !make 1/2 time step advance

         xp(l,1) = xp(l,1) + dth*vp(l,1)
c         ijkp(l,1) = nint(xp(l,1)/dx) 

         xp(l,2) = xp(l,2) + dth*vp(l,2)
c         ijkp(l,2) = nint(xp(l,2)/dy) 

         xp(l,3) = xp(l,3) + dth*vp(l,3)
c         k=1
c         do 15 while(xp(l,3) .gt. qz(k))  !find k on non-uniform 
c            ijkp(l,3) = k                 !grid
c            k=k+1
c 15      continue
c         k=ijkp(l,3)
c         if (xp(l,3) .gt. (qz(k)+(dz_grid(k)/2))) then
c            ijkp(l,3) = k+1
c         endif

         if (xp(l,1) .gt. qx(nx-1)) then
           xp(l,1) = qx(1) + ( xp(l,1) - qx(nx-1) )
         else if (xp(l,1) .le. qx(1)) then
           xp(l,1) = qx(nx-1) - (qx(1) - xp(l,1))
         endif
c         ijkp(l,1) = nint(xp(l,1)/dx)

         if (xp(l,2) .gt. qy(ny-1)) then
           xp(l,2) = qy(1) + ( xp(l,2) - qy(ny-1) )
         else if (xp(l,2) .le. qy(1)) then
           xp(l,2) = qy(ny-1) - (qy(1) - xp(l,2))
         endif
c         ijkp(l,2) = nint(xp(l,2)/dy)

         if (xp(l,3) .gt. qz(nz-1)) then
           xp(l,3) = qz(1) + ( xp(l,3) - qz(nz-1) )
c           ijkp(l,3) = nint(xp(l,3)/delz)
         else if (xp(l,3) .le. qz(1)) then
           xp(l,3) = qz(nz-1) - (qz(1) - xp(l,3))
c           ijkp(l,3) = nint(xp(l,3)/delz)
         endif

 10      continue


c      where (xp(:,1) .gt. qx(nx-1))
cc         ijkp(:,1) = 1
cc         wquad(:,1) = 0
c         xp(:,1) = qx(1) + ( xp(:,1) - qx(nx-1) )
cc         ijkp(:,1) = nint(xp(:,1)/dx)
c      endwhere


c      where (xp(:,1) .le. qx(1)) 
cc         ijkp(:,1) = nx-1
cc         wquad(:,1) = -1
c         xp(:,1) = qx(nx-1) - (qx(1) - xp(:,1))
cc         ijkp(:,1) = nint(xp(:,1)/dx)
cc         xp(:,2) = pad_ranf()*qy(ny)
cc         ijkp(:,2) = ninit(xp(:,2)/dy
cc         vp(:,1) = -vsw
cc         vp(:,2) = 0.0
cc         vp(:,3) = 0.0
c      endwhere

c      where (xp(:,2) .ge. qy(ny-1))
cc         ijkp(:,2) = 1
cc         wquad(:,2) = 0
c         xp(:,2) = qy(1) + ( xp(:,2) - qy(ny-1) )
cc         ijkp(:,2) = nint(xp(:,2)/dy)
c      endwhere

c      where (xp(:,2) .le. qy(1)) 
cc         ijkp(:,2) = ny-1
cc         wquad(:,2) = -1
c         xp(:,2) = qy(ny-1) - (qy(1) - xp(:,2))
cc         ijkp(:,2) = nint(xp(:,2)/dy)
c      endwhere


cc periodic boundary conditions
c      where (xp(:,3) .le. qz(1)) 
cc         ijkp(:,3) = nz
cc         wquad(:,3) = -1.0
c         xp(:,3) = qz(nz-1) - (qz(1) - xp(:,3))
cc         ijkp(:,3) = nint(xp(:,3)/delz)
cc         vp(:,1) = -vp(:,1)
c      endwhere

cc periodic boundary conditions
c      where (xp(:,3) .ge. qz(nz-1))
cc         ijkp(:,3) = 1
cc         wquad(:,3) = 0.0
c         xp(:,3) = qz(1) + (xp(:,3) - qz(nz-1) )
cc         ijkp(:,3) = nint(xp(:,3)/delz)
cc         vp(:,1) = -vp(:,1)
c      endwhere


      return
      end SUBROUTINE move_ion_half
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE check_min_den_boundary(np,xp,vp,up)
c----------------------------------------------------------------------
      !!include 'incurv.h'

      real np(nx,ny,nz),
     x     xp(Ni_max,3),
     x     vp(Ni_max,3),
     x     up(nx,ny,nz,3)

      real minden           !minimum wake density 
      real den_part ! density of 1 particle per cell
      real ak
c      real btot,a1,a2,womega,phi,deltat
      integer npart,ipart
      integer flg
      real theta,phi,f,v,rnd
      

      den_part = 1/(beta*dx**3)

      k = nz-1   !top boundary
      minden = np_top-den_part

c      print *,'minden...',np_top,den_part
      do i = 2,nx-1
         do j = 2,ny-1

            if (np(i,j,k) .le. (np_top-den_part)) then
               npart = nint((np_top - np(i,j,k))/den_part)

               if (my_rank .eq. nint(pad_ranf()*procnum)) then
                  do ipart = 1,npart 
                     
                     l=Ni_tot + 1 !beginning array element for new borns    
                     
c                     phi = 2.0*pi*pad_ranf()
c                     flg = 0
c                     do 30 while (flg .eq. 0)
c                        theta = pi*pad_ranf()
c                        f = sin(theta)
c                        rnd = pad_ranf()
c                        if (f .ge. rnd) flg = 1
c 30                  continue
                     
c                     flg = 0
c                     do 40 while (flg .eq. 0)
c                        v = (100*pad_ranf())
c                        f = exp(-(v)**2 / vth**2)
c                        rnd = pad_ranf()
c                        if (f .ge. rnd) then 
c                           flg = 1
c                           
c                           vp(l,1) = vsw + v*cos(phi)*sin(theta)
c                           vp(l,2) = v*sin(phi)*sin(theta)
c                           vp(l,3) = v*cos(theta)
c                        endif
c                        
c 40                  continue
  
                     flg = 0
                     do 40 while (flg .eq. 0)
                        v = (2*vth_max*pad_ranf())-vth_max
                        f = exp(-(v)**2 / vth_top**2)
                        rnd = pad_ranf()
                        if (f .ge. rnd) then 
                           flg = 1
                           vx = v
                        endif
 40                  continue
                     flg = 0
                     do 42 while (flg .eq. 0)
                        v = (2*vth_max*pad_ranf())-vth_max
                        f = exp(-(v)**2 / vth_top**2)
                        rnd = pad_ranf()
                        if (f .ge. rnd) then 
                           flg = 1
                           vy = v
                        endif
 42                  continue
                     flg = 0
                     do 44 while (flg .eq. 0)
                        v = (2*vth_max*pad_ranf())-vth_max
                        f = exp(-(v)**2 / vth_top**2)
                        rnd = pad_ranf()
                        if (f .ge. rnd) then 
                           flg = 1
                           vz = v
                        endif
 44                  continue
                     
                     vp(l,1) = vsw*(tanh((qz(nz)-qz(nz/2))/(Lo)))+vx 
                     vp(l,2) = vy !*sin(phi)*sin(theta)
                     vp(l,3) = vz !*cos(theta)
                     
                     xp(l,1) = qx(i) + (0.5-pad_ranf())*dx
                     xp(l,2) = qy(j) + (0.5-pad_ranf())*dy
                     xp(l,3) = qz(k) + 2.0*(0.5-pad_ranf())*dz_grid(k)

                     
                     ijkp(l,1) = nint(xp(l,1)/dx) !particle grid location index
                     ijkp(l,2) = nint(xp(l,2)/dy)
                     
                     kk=1
                     do 100 while(xp(l,3) .gt. qz(kk)) !find k on non-uniform 
                        ijkp(l,3) = kk !grid
                        kk=kk+1
 100                 continue
                     kk=ijkp(l,3)
                     if (xp(l,3) .gt. (qz(kk)+(dz_grid(kk)/2))) then
                        ijkp(l,3) = kk+1
                     endif
                     mrat(l) = mion/m_top
                     m_arr(l) = m_top

c add energy
                     do m=1,3
                        input_E = input_E + 
     x                       0.5*m_arr(l)*(vp(l,m)*km_to_m)**2 /beta
c                        input_p(m) = input_p(m) - m_arr(l)*vp(l,m)/beta
                     enddo         

                     Ni_tot = Ni_tot + 1
                  enddo
               endif
            endif
         enddo
         
      enddo

      k = 2   !bottom boundary
      minden = np_bottom-den_part
      do i = 2,nx-1
         do j = 2,ny-1

            if (np(i,j,k) .le. (np_bottom-den_part)) then
               npart = nint((np_bottom - np(i,j,k))/den_part)
c               print *,'npart bottom...',npart
               if (my_rank .eq. nint(pad_ranf()*procnum)) then
                  do ipart = 1,npart 

                     l=Ni_tot + 1 !beginning array element for new borns    
                     
c                     phi = 2.0*pi*pad_ranf()
c                     flg = 0
c                     do 35 while (flg .eq. 0)
c                        theta = pi*pad_ranf()
c                        f = sin(theta)
c                        rnd = pad_ranf()
c                        if (f .ge. rnd) flg = 1
c 35                  continue
                     
c                     flg = 0
c                     do 45 while (flg .eq. 0)
c                        v = (100*pad_ranf())
c                        f = exp(-(v)**2 / vth**2)
c                        rnd = pad_ranf()
c                        if (f .ge. rnd) then 
c                           flg = 1
c                           
c                           vp(l,1) = vsw + v*cos(phi)*sin(theta)
c                           vp(l,2) = v*sin(phi)*sin(theta)
c                           vp(l,3) = v*cos(theta)
c                        endif
                        
c 45                  continue

  
                     flg = 0
                     do 46 while (flg .eq. 0)
                        v = (2*vth_max*pad_ranf())-vth_max
                        f = exp(-(v)**2 / vth_bottom**2)
                        rnd = pad_ranf()
                        if (f .ge. rnd) then 
                           flg = 1
                           vx = v
                        endif
 46                  continue
                     flg = 0
                     do 48 while (flg .eq. 0)
                        v = (2*vth_max*pad_ranf())-vth_max
                        f = exp(-(v)**2 / vth_bottom**2)
                        rnd = pad_ranf()
                        if (f .ge. rnd) then 
                           flg = 1
                           vy = v
                        endif
 48                  continue
                     flg = 0
                     do 50 while (flg .eq. 0)
                        v = (2*vth_max*pad_ranf())-vth_max
                        f = exp(-(v)**2 / vth_bottom**2)
                        rnd = pad_ranf()
                        if (f .ge. rnd) then 
                           flg = 1
                           vz = v
                        endif
 50                  continue
                     
                     vp(l,1) = vsw*(tanh((qz(1)-qz(nz/2))/(Lo)))+vx 
                     vp(l,2) = vy !*sin(phi)*sin(theta)
                     vp(l,3) = vz !*cos(theta)
                     
                     xp(l,1) = qx(i) - (0.5-pad_ranf())*dx
                     xp(l,2) = qy(j) - (0.5-pad_ranf())*dy
                     xp(l,3) = qz(k) - 2.0*(0.5-pad_ranf())*dz_grid(k)
                     
                     ijkp(l,1) = nint(xp(l,1)/dx) !particle grid location index
                     ijkp(l,2) = nint(xp(l,2)/dy)
                     
                     kk=1
                     do 110 while(xp(l,3) .gt. qz(kk)) !find k on non-uniform 
                        ijkp(l,3) = kk !grid
                        kk=kk+1
 110                 continue
                     kk=ijkp(l,3)
                     if (xp(l,3) .gt. (qz(kk)+(dz_grid(kk)/2))) then
                        ijkp(l,3) = kk+1
                     endif
                     mrat(l) = mion/m_bottom
                     m_arr(l) = m_bottom

c add energy
                     do m=1,3
                        input_E = input_E + 
     x                       0.5*m_arr(l)*(vp(l,m)*km_to_m)**2 /beta
c                        input_p(m) = input_p(m) - m_arr(l)*vp(l,m)/beta
                     enddo         


                     Ni_tot = Ni_tot + 1
                  enddo
               endif
            endif
         enddo
         
      enddo
      
      
      return
      end SUBROUTINE check_min_den_boundary
c----------------------------------------------------------------------



c----------------------------------------------------------------------
      SUBROUTINE check_min_den_boundary_1(np,xp,vp,up)
c----------------------------------------------------------------------
      !!include 'incurv.h'

      real np(nx,ny,nz),
     x     xp(Ni_max,3),
     x     vp(Ni_max,3),
     x     up(nx,ny,nz,3)

      real minden           !minimum wake density 
      real den_part ! density of 1 particle per cell
      real ak
c      real btot,a1,a2,womega,phi,deltat
      integer npart,ipart
      integer flg
      real theta,phi,f,v,rnd
      

      den_part = 1/(beta*dx**3)

      k = nz-1   !top boundary
      minden = np_top-den_part

      do i = 2,nx-1
         do j = 2,ny-1

            if (np(i,j,k) .le. (np_top-den_part)) then
               npart = nint((np_top - np(i,j,k))/den_part)

               if (my_rank .eq. nint(pad_ranf()*procnum)) then
                  do ipart = 1,npart 
                     
                     l=Ni_tot + 1 !beginning array element for new borns    
                     
c                     phi = 2.0*pi*pad_ranf()
c                     flg = 0
c                     do 30 while (flg .eq. 0)
c                        theta = pi*pad_ranf()
c                        f = sin(theta)
c                        rnd = pad_ranf()
c                        if (f .ge. rnd) flg = 1
c 30                  continue
                     
c                     flg = 0
c                     do 40 while (flg .eq. 0)
c                        v = (100*pad_ranf())
c                        f = exp(-(v)**2 / vth**2)
c                        rnd = pad_ranf()
c                        if (f .ge. rnd) then 
c                           flg = 1
c                           
c                           vp(l,1) = vsw + v*cos(phi)*sin(theta)
c                           vp(l,2) = v*sin(phi)*sin(theta)
c                           vp(l,3) = v*cos(theta)
c                        endif
c                        
c 40                  continue
  
                     flg = 0
                     do 40 while (flg .eq. 0)
                        
                        vx = (600*pad_ranf())-300
                        vy = (600*pad_ranf())-300
                        vz = (600*pad_ranf())-300
                        
                        v = sqrt(vx**2 + vy**2 + vz**2)
                        f = exp(-(v)**2 / vth**2)
                        rnd = pad_ranf()
                        if (f .ge. rnd) then 
                           flg = 1
                        endif
 40                  continue

c                     do 40 while (flg .eq. 0)
c                        v = (200*pad_ranf())-100
c                        f = exp(-(v)**2 / vth_top**2)
c                        rnd = pad_ranf()
c                        if (f .ge. rnd) then 
c                           flg = 1
c                           vx = v
c                        endif
c 40                  continue
c                     flg = 0
c                     do 42 while (flg .eq. 0)
c                        v = (200*pad_ranf())-100
c                        f = exp(-(v)**2 / vth_top**2)
c                        rnd = pad_ranf()
c                        if (f .ge. rnd) then 
c                           flg = 1
c                           vy = v
c                        endif
c 42                  continue
c                     flg = 0
c                     do 44 while (flg .eq. 0)
c                        v = (200*pad_ranf())-100
c                        f = exp(-(v)**2 / vth_top**2)
c                        rnd = pad_ranf()
c                        if (f .ge. rnd) then 
c                           flg = 1
c                           vz = v
c                        endif
c 44                  continue
                     
                     vp(l,1) = vsw*(tanh((qz(nz)-qz(nz/2))/(Lo)))+vx 
                     vp(l,2) = vy !*sin(phi)*sin(theta)
                     vp(l,3) = vz !*cos(theta)
                     
                     xp(l,1) = qx(i) + (0.5-pad_ranf())*dx
                     xp(l,2) = qy(j) + (0.5-pad_ranf())*dy
                     xp(l,3) = qz(k) + 2.0*(0.5-pad_ranf())*dz_grid(k)

                     
                     ijkp(l,1) = nint(xp(l,1)/dx) !particle grid location index
                     ijkp(l,2) = nint(xp(l,2)/dy)
                     
                     kk=1
                     do 100 while(xp(l,3) .gt. qz(kk)) !find k on non-uniform 
                        ijkp(l,3) = kk !grid
                        kk=kk+1
 100                 continue
                     kk=ijkp(l,3)
                     if (xp(l,3) .gt. (qz(kk)+(dz_grid(kk)/2))) then
                        ijkp(l,3) = kk+1
                     endif
                     mrat(l) = 1.0
                     m_arr(l) = mion

c add energy
                     do m=1,3
                        input_E = input_E + 
     x                       0.5*m_arr(l)*(vp(l,m)*km_to_m)**2 /beta
c                        input_p(m) = input_p(m) - m_arr(l)*vp(l,m)/beta
                     enddo         

                     Ni_tot = Ni_tot + 1
                  enddo
               endif
            endif
         enddo
         
      enddo

      k = 2   !bottom boundary
      minden = np_bottom-den_part
      do i = 2,nx-1
         do j = 2,ny-1

            if (np(i,j,k) .le. (np_bottom-den_part)) then
               npart = nint((np_bottom - np(i,j,k))/den_part)
c               print *,'npart bottom...',npart
               if (my_rank .eq. nint(pad_ranf()*procnum)) then
                  do ipart = 1,npart 

                     l=Ni_tot + 1 !beginning array element for new borns    
                     
c                     phi = 2.0*pi*pad_ranf()
c                     flg = 0
c                     do 35 while (flg .eq. 0)
c                        theta = pi*pad_ranf()
c                        f = sin(theta)
c                        rnd = pad_ranf()
c                        if (f .ge. rnd) flg = 1
c 35                  continue
                     
c                     flg = 0
c                     do 45 while (flg .eq. 0)
c                        v = (100*pad_ranf())
c                        f = exp(-(v)**2 / vth**2)
c                        rnd = pad_ranf()
c                        if (f .ge. rnd) then 
c                           flg = 1
c                           
c                           vp(l,1) = vsw + v*cos(phi)*sin(theta)
c                           vp(l,2) = v*sin(phi)*sin(theta)
c                           vp(l,3) = v*cos(theta)
c                        endif
                        
c 45                  continue

  
                     flg = 0
                     do 46 while (flg .eq. 0)
                        
                        vx = (600*pad_ranf())-300
                        vy = (600*pad_ranf())-300
                        vz = (600*pad_ranf())-300
                        
                        v = sqrt(vx**2 + vy**2 + vz**2)
                        f = exp(-(v)**2 / vth**2)
                        rnd = pad_ranf()
                        if (f .ge. rnd) then 
                           flg = 1
                        endif
 46                  continue
                     

c                     do 46 while (flg .eq. 0)
c                        v = (200*pad_ranf())-100
c                        f = exp(-(v)**2 / vth_bottom**2)
c                        rnd = pad_ranf()
c                        if (f .ge. rnd) then 
c                           flg = 1
c                           vx = v
c                        endif
c 46                  continue
c                     flg = 0
c                     do 48 while (flg .eq. 0)
c                        v = (200*pad_ranf())-100
c                        f = exp(-(v)**2 / vth_bottom**2)
c                        rnd = pad_ranf()
c                        if (f .ge. rnd) then 
c                           flg = 1
c                           vy = v
c                        endif
c 48                  continue
c                     flg = 0
c                     do 50 while (flg .eq. 0)
c                        v = (200*pad_ranf())-100
c                        f = exp(-(v)**2 / vth_bottom**2)
c                        rnd = pad_ranf()
c                        if (f .ge. rnd) then 
c                           flg = 1
c                           vz = v
c                        endif
c 50                  continue
                     
                     vp(l,1) = vsw*(tanh((qz(1)-qz(nz/2))/(Lo)))+vx 
                     vp(l,2) = vy !*sin(phi)*sin(theta)
                     vp(l,3) = vz !*cos(theta)
                     
                     xp(l,1) = qx(i) - (0.5-pad_ranf())*dx
                     xp(l,2) = qy(j) - (0.5-pad_ranf())*dy
                     xp(l,3) = qz(k) - 2.0*(0.5-pad_ranf())*dz_grid(k)
                     
                     ijkp(l,1) = nint(xp(l,1)/dx) !particle grid location index
                     ijkp(l,2) = nint(xp(l,2)/dy)
                     
                     kk=1
                     do 110 while(xp(l,3) .gt. qz(kk)) !find k on non-uniform 
                        ijkp(l,3) = kk !grid
                        kk=kk+1
 110                 continue
                     kk=ijkp(l,3)
                     if (xp(l,3) .gt. (qz(kk)+(dz_grid(kk)/2))) then
                        ijkp(l,3) = kk+1
                     endif
                     mrat(l) = 1.0
                     m_arr(l) = mion

c add energy
                     do m=1,3
                        input_E = input_E + 
     x                       0.5*m_arr(l)*(vp(l,m)*km_to_m)**2 /beta
c                        input_p(m) = input_p(m) - m_arr(l)*vp(l,m)/beta
                     enddo         


                     Ni_tot = Ni_tot + 1
                  enddo
               endif
            endif
         enddo
         
      enddo
      
      
      return
      end SUBROUTINE check_min_den_boundary_1
c----------------------------------------------------------------------





c----------------------------------------------------------------------
      SUBROUTINE get_interp_weights_2(xp)
c Weights are used for trilinear interpolation to/from main cell
c centers to particle positions.  For each particle there are 8
c grid points associated with the interpolation.  These 8 points
c are determined by the location of the particle within the main
c cell.  There are 8 sets of 8 grid points for each cell.
c----------------------------------------------------------------------
      !!include 'incurv.h'

      real xp(Ni_max,3)
      real x1,x2,y1,y2,z1,z2,vol


c      where((xp(:,1) .le. qx(ijkp(:,1))) .and. 
c     x    (xp(:,2) .le. qy(ijkp(:,2))) .and.
c     x    (xp(:,3) .le. qz(ijkp(:,3)))) 
c         wquad(:,1) = -1             
c         wquad(:,2) = -1
c         wquad(:,3) = -1  
c         wght(:,1) = abs(xp(:,1)-qx(ijkp(:,1)))*
c     x               abs(xp(:,2)-qy(ijkp(:,2)))*
c     x               abs(xp(:,3)-qz(ijkp(:,3)))/
c     x               (dx*dy*(qz(ijkp(:,3))-qz(ijkp(:,3)-1)))
c         wght(:,2) = abs(xp(:,1)-qx(ijkp(:,1)-1))*
c     x               abs(xp(:,2)-qy(ijkp(:,2)))*
c     x               abs(xp(:,3)-qz(ijkp(:,3)))/
c     x               (dx*dy*(qz(ijkp(:,3))-qz(ijkp(:,3)-1)))
c         wght(:,3) = abs(xp(:,1)-qx(ijkp(:,1)))*
c     x               abs(xp(:,2)-qy(ijkp(:,2)))*
c     x               abs(xp(:,3)-qz(ijkp(:,3)-1))/
c     x               (dx*dy*(qz(ijkp(:,3))-qz(ijkp(:,3)-1)))
c         wght(:,4) = abs(xp(:,1)-qx(ijkp(:,1)-1))*
c     x               abs(xp(:,2)-qy(ijkp(:,2)))*
c     x               abs(xp(:,3)-qz(ijkp(:,3)-1))/
c     x               (dx*dy*(qz(ijkp(:,3))-qz(ijkp(:,3)-1)))
c         wght(:,5) = abs(xp(:,1)-qx(ijkp(:,1)))*
c     x               abs(xp(:,2)-qy(ijkp(:,2)-1))*
c     x               abs(xp(:,3)-qz(ijkp(:,3)))/
c     x               (dx*dy*(qz(ijkp(:,3))-qz(ijkp(:,3)-1)))
c         wght(:,6) = abs(xp(:,1)-qx(ijkp(:,1)-1))*
c     x               abs(xp(:,2)-qy(ijkp(:,2)-1))*
c     x               abs(xp(:,3)-qz(ijkp(:,3)))/
c     x               (dx*dy*(qz(ijkp(:,3))-qz(ijkp(:,3)-1)))
c         wght(:,7) = abs(xp(:,1)-qx(ijkp(:,1)))*
c     x               abs(xp(:,2)-qy(ijkp(:,2)-1))*
c     x               abs(xp(:,3)-qz(ijkp(:,3)-1))/
c     x               (dx*dy*(qz(ijkp(:,3))-qz(ijkp(:,3)-1)))
c         wght(:,8) = abs(xp(:,1)-qx(ijkp(:,1)-1))*
c     x               abs(xp(:,2)-qy(ijkp(:,2)-1))*
c     x               abs(xp(:,3)-qz(ijkp(:,3)-1))/
c     x               (dx*dy*(qz(ijkp(:,3))-qz(ijkp(:,3)-1)))

c      endwhere


      do 10 l=1,Ni_tot

         i=ijkp(l,1)
         j=ijkp(l,2)
         k=ijkp(l,3)

c         if ((ijkp(l,1) .gt. nx) .or. (ijkp(l,2) .gt. ny) .or. 
c     x       (ijkp(l,3) .gt. nz) .or. (ijkp(l,1) .lt. 2) .or. 
c     x       (ijkp(l,2) .lt. 2) .or. (ijkp(l,3) .lt. 2)) then
c            write(*,*) 'Out of bounds...',i,j,k
cc            call remove_ion(xp,vp,vp1,l)
c            endif

c 111111111111111111111111111111111111111111111111111111111111111111111


      if ((xp(l,1) .le. qx(i)) .and. (xp(l,2) .le. qy(j)) .and. 
     x    (xp(l,3) .le. qz(k))) then
c         wquad(l,1) = -1
c         wquad(l,2) = -1
c         wquad(l,3) = -1
         vol = dx*dy*(qz(k)-qz(k-1))
         x1=abs(xp(l,1)-qx(i-1))
         x2=abs(xp(l,1)-qx(i))
         y1=abs(xp(l,2)-qy(j-1))
         y2=abs(xp(l,2)-qy(j))
         z1=abs(xp(l,3)-qz(k-1))
         z2=abs(xp(l,3)-qz(k))
         wght(l,1) = x2*y2*z2/vol
         wght(l,2) = x1*y2*z2/vol
         wght(l,3) = x2*y2*z1/vol
         wght(l,4) = x1*y2*z1/vol
         wght(l,5) = x2*y1*z2/vol
         wght(l,6) = x1*y1*z2/vol
         wght(l,7) = x2*y1*z1/vol
         wght(l,8) = x1*y1*z1/vol

      endif

c 222222222222222222222222222222222222222222222222222222222222222222222

      if ((xp(l,1) .gt. qx(i)) .and. (xp(l,2) .le. qy(j)) .and. 
     x    (xp(l,3) .le. qz(k))) then
c         wquad(l,1) = 0
c         wquad(l,2) = -1
c         wquad(l,3) = -1
         vol = dx*dy*(qz(k)-qz(k-1))
         x1=abs(xp(l,1)-qx(i))
         x2=abs(xp(l,1)-qx(i+1))
         y1=abs(xp(l,2)-qy(j-1))
         y2=abs(xp(l,2)-qy(j))
         z1=abs(xp(l,3)-qz(k-1))
         z2=abs(xp(l,3)-qz(k))
         wght(l,1) = x2*y2*z2/vol
         wght(l,2) = x1*y2*z2/vol
         wght(l,3) = x2*y2*z1/vol
         wght(l,4) = x1*y2*z1/vol
         wght(l,5) = x2*y1*z2/vol
         wght(l,6) = x1*y1*z2/vol
         wght(l,7) = x2*y1*z1/vol
         wght(l,8) = x1*y1*z1/vol

      endif

c 333333333333333333333333333333333333333333333333333333333333333333333

      if ((xp(l,1) .le. qx(i)) .and. (xp(l,2) .le. qy(j)) .and. 
     x    (xp(l,3) .gt. qz(k))) then
c         wquad(l,1) = -1
c         wquad(l,2) = -1
c         wquad(l,3) = 0
         vol = dx*dy*(qz(k+1)-qz(k))
         x1=abs(xp(l,1)-qx(i-1))
         x2=abs(xp(l,1)-qx(i))
         y1=abs(xp(l,2)-qy(j-1))
         y2=abs(xp(l,2)-qy(j))
         z1=abs(xp(l,3)-qz(k))
         z2=abs(xp(l,3)-qz(k+1))
         wght(l,1) = x2*y2*z2/vol
         wght(l,2) = x1*y2*z2/vol
         wght(l,3) = x2*y2*z1/vol
         wght(l,4) = x1*y2*z1/vol
         wght(l,5) = x2*y1*z2/vol
         wght(l,6) = x1*y1*z2/vol
         wght(l,7) = x2*y1*z1/vol
         wght(l,8) = x1*y1*z1/vol

      endif

c 444444444444444444444444444444444444444444444444444444444444444444444

      if ((xp(l,1) .gt. qx(i)) .and. (xp(l,2) .le. qy(j)) .and. 
     x    (xp(l,3) .gt. qz(k))) then
c         wquad(l,1) = 0
c         wquad(l,2) = -1
c         wquad(l,3) = 0
         vol = dx*dy*(qz(k+1)-qz(k))
         x1=abs(xp(l,1)-qx(i))
         x2=abs(xp(l,1)-qx(i+1))
         y1=abs(xp(l,2)-qy(j-1))
         y2=abs(xp(l,2)-qy(j))
         z1=abs(xp(l,3)-qz(k))
         z2=abs(xp(l,3)-qz(k+1))
         wght(l,1) = x2*y2*z2/vol
         wght(l,2) = x1*y2*z2/vol
         wght(l,3) = x2*y2*z1/vol
         wght(l,4) = x1*y2*z1/vol
         wght(l,5) = x2*y1*z2/vol
         wght(l,6) = x1*y1*z2/vol
         wght(l,7) = x2*y1*z1/vol
         wght(l,8) = x1*y1*z1/vol

      endif

c 555555555555555555555555555555555555555555555555555555555555555555555

      if ((xp(l,1) .le. qx(i)) .and. (xp(l,2) .gt. qy(j)) .and. 
     x    (xp(l,3) .le. qz(k))) then
c         wquad(l,1) = -1
c         wquad(l,2) = 0
c         wquad(l,3) = -1
         vol = dx*dy*(qz(k)-qz(k-1))
         x1=abs(xp(l,1)-qx(i-1))
         x2=abs(xp(l,1)-qx(i))
         y1=abs(xp(l,2)-qy(j))
         y2=abs(xp(l,2)-qy(j+1))
         z1=abs(xp(l,3)-qz(k-1))
         z2=abs(xp(l,3)-qz(k))
         wght(l,1) = x2*y2*z2/vol
         wght(l,2) = x1*y2*z2/vol
         wght(l,3) = x2*y2*z1/vol
         wght(l,4) = x1*y2*z1/vol
         wght(l,5) = x2*y1*z2/vol
         wght(l,6) = x1*y1*z2/vol
         wght(l,7) = x2*y1*z1/vol
         wght(l,8) = x1*y1*z1/vol

      endif

c 666666666666666666666666666666666666666666666666666666666666666666666

      if ((xp(l,1) .gt. qx(i)) .and. (xp(l,2) .gt. qy(j)) .and. 
     x    (xp(l,3) .le. qz(k))) then
c         wquad(l,1) = 0
c         wquad(l,2) = 0
c         wquad(l,3) = -1
         vol = dx*dy*(qz(k)-qz(k-1))
         x1=abs(xp(l,1)-qx(i))
         x2=abs(xp(l,1)-qx(i+1))
         y1=abs(xp(l,2)-qy(j))
         y2=abs(xp(l,2)-qy(j+1))
         z1=abs(xp(l,3)-qz(k-1))
         z2=abs(xp(l,3)-qz(k))
         wght(l,1) = x2*y2*z2/vol
         wght(l,2) = x1*y2*z2/vol
         wght(l,3) = x2*y2*z1/vol
         wght(l,4) = x1*y2*z1/vol
         wght(l,5) = x2*y1*z2/vol
         wght(l,6) = x1*y1*z2/vol
         wght(l,7) = x2*y1*z1/vol
         wght(l,8) = x1*y1*z1/vol

      endif

c 777777777777777777777777777777777777777777777777777777777777777777777

      if ((xp(l,1) .le. qx(i)) .and. (xp(l,2) .gt. qy(j)) .and. 
     x    (xp(l,3) .gt. qz(k))) then
c         wquad(l,1) = -1
c         wquad(l,2) = 0
c         wquad(l,3) = 0
         vol = dx*dy*(qz(k+1)-qz(k))
         x1=abs(xp(l,1)-qx(i-1))
         x2=abs(xp(l,1)-qx(i))
         y1=abs(xp(l,2)-qy(j))
         y2=abs(xp(l,2)-qy(j+1))
         z1=abs(xp(l,3)-qz(k))
         z2=abs(xp(l,3)-qz(k+1))
         wght(l,1) = x2*y2*z2/vol
         wght(l,2) = x1*y2*z2/vol
         wght(l,3) = x2*y2*z1/vol
         wght(l,4) = x1*y2*z1/vol
         wght(l,5) = x2*y1*z2/vol
         wght(l,6) = x1*y1*z2/vol
         wght(l,7) = x2*y1*z1/vol
         wght(l,8) = x1*y1*z1/vol

      endif

c 888888888888888888888888888888888888888888888888888888888888888888888

      if ((xp(l,1) .gt. qx(i)) .and. (xp(l,2) .gt. qy(j)) .and. 
     x    (xp(l,3) .gt. qz(k))) then
c         wquad(l,1) = 0
c         wquad(l,2) = 0
c         wquad(l,3) = 0
         vol = dx*dy*(qz(k+1)-qz(k))
         x1=abs(xp(l,1)-qx(i))
         x2=abs(xp(l,1)-qx(i+1))
         y1=abs(xp(l,2)-qy(j))
         y2=abs(xp(l,2)-qy(j+1))
         z1=abs(xp(l,3)-qz(k))
         z2=abs(xp(l,3)-qz(k+1))
         wght(l,1) = x2*y2*z2/vol
         wght(l,2) = x1*y2*z2/vol
         wght(l,3) = x2*y2*z1/vol
         wght(l,4) = x1*y2*z1/vol
         wght(l,5) = x2*y1*z2/vol
         wght(l,6) = x1*y1*z2/vol
         wght(l,7) = x2*y1*z1/vol
         wght(l,8) = x1*y1*z1/vol

      endif

 10   continue


      return
      end SUBROUTINE get_interp_weights_2
c----------------------------------------------------------------------



c----------------------------------------------------------------------
      SUBROUTINE get_interp_weights(xp)
c Weights are used for trilinear interpolation to/from main cell
c centers to particle positions.  For each particle there are 8
c grid points associated with the interpolation.  These 8 points
c are determined by the location of the particle within the main
c cell.  There are 8 sets of 8 grid points for each cell.
c----------------------------------------------------------------------
      !!include 'incurv.h'

      real xp(Ni_max,3)
      real x1,x2,y1,y2,z1,z2,vol


      do 10 l=1,Ni_tot

c         i = floor(xp(l,1)/dx) 
c         ijkp(l,1) = i

         i=0
 16      continue
         i = i + 1
         if (xp(l,1) .gt. qx(i)) go to 16 !find i on non-uniform 
         i = i-1
         ijkp(l,1)= i


         j = floor(xp(l,2)/dy) 
         ijkp(l,2) = j

         k=0
 15      continue
         k = k + 1
         if (xp(l,3) .gt. qz(k)) go to 15  !find k on non-uniform 
         k = k-1
         ijkp(l,3)= k

c         vol = 1.0/(dx*dy*(qz(k+1)-qz(k)))
         vol = 1.0/((qx(i+1)-qx(i))*(qy(j+1)-qy(j))*(qz(k+1)-qz(k)))
         x1=abs(xp(l,1)-qx(i))
         x2=abs(xp(l,1)-qx(i+1))
         y1=abs(xp(l,2)-qy(j))
         y2=abs(xp(l,2)-qy(j+1))
         z1=abs(xp(l,3)-qz(k))
         z2=abs(xp(l,3)-qz(k+1))
         wght(l,1) = x2*y2*z2*vol
         wght(l,2) = x1*y2*z2*vol
         wght(l,3) = x2*y2*z1*vol
         wght(l,4) = x1*y2*z1*vol
         wght(l,5) = x2*y1*z2*vol
         wght(l,6) = x1*y1*z2*vol
         wght(l,7) = x2*y1*z1*vol
         wght(l,8) = x1*y1*z1*vol


 10   continue

      return
      end SUBROUTINE get_interp_weights
c----------------------------------------------------------------------



c----------------------------------------------------------------------
      SUBROUTINE update_np(np)
c Weight density to eight nearest grid points.
c----------------------------------------------------------------------
      !!include 'incurv.h'

      real np(nx,ny,nz)

      real volb              !cell volume times beta
      real recvbuf(nx*ny*nz)
      integer count
      count = nx*ny*nz

c      real sumnp,vol

      do 10 i=1,nx
         do 10 j=1,ny
            do 10 k=1,nz
               np(i,j,k)=0.0
 10            continue

      do 20 l=1,Ni_tot

         i=ijkp(l,1)!+wquad(l,1)
         j=ijkp(l,2)!+wquad(l,2)
         k=ijkp(l,3)!+wquad(l,3)

c         if (i .lt. 1) i = nx-1
c         if (j .lt. 1) i = ny-1
c         if (k .lt. 1) i = nz-1

         ip = i+1
         jp = j+1
         kp = k+1

c         if (ip .ge. nx) ip = 2      !periodic boundary conditions
c         if (jp .ge. ny) jp = 2      !periodic boundary conditions
c         if (kp .ge. nz) kp = 2      !periodic boundary conditions

c         volb = 1.0/(dx*dy*(qz(k+1)-qz(k))*beta)
         volb = 1.0/(dx_grid(i)*dy_grid(j)*dz_grid(k)*beta)
c         volb = 1.0/((qx(ip)-qx(i))*(qy(jp)-qy(j))*(qz(kp)-qz(k))*beta)
c         volb = dx*dy*dz_cell(k)*beta

         np(i,j,k) = np(i,j,k) + wght(l,1)*volb
         np(ip,j,k) = np(ip,j,k) + wght(l,2)*volb
         np(i,j,kp) = np(i,j,kp) + wght(l,3)*volb
         np(ip,j,kp) = np(ip,j,kp) + wght(l,4)*volb
         np(i,jp,k) = np(i,jp,k) + wght(l,5)*volb
         np(ip,jp,k) = np(ip,jp,k) + wght(l,6)*volb
         np(i,jp,kp) = np(i,jp,kp) + wght(l,7)*volb
         np(ip,jp,kp) = np(ip,jp,kp) + wght(l,8)*volb

 20      continue

c use for periodic boundary conditions
         np(nx-1,:,:) = np(nx-1,:,:)+np(1,:,:)
         np(:,ny-1,:) = np(:,ny-1,:)+np(:,1,:)
         np(:,:,nz-1) = np(:,:,nz-1)+np(:,:,1)


         call MPI_Barrier(MPI_COMM_WORLD,ierr)

         call MPI_ALLREDUCE(np(:,:,:),recvbuf,count,
     x        MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)

         np(:,:,:) = reshape(recvbuf,(/nx,ny,nz/))


      call periodic_scalar(np)

      return
      end SUBROUTINE update_np
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE update_rho(mnp)
c Weight density to eight nearest grid points.
c----------------------------------------------------------------------
      !!include 'incurv.h'

      real mnp(nx,ny,nz)

      real volb              !cell volume times beta
      real recvbuf(nx*ny*nz)
      integer count
      count = nx*ny*nz

c      real sumnp,vol

      do 10 i=1,nx
         do 10 j=1,ny
            do 10 k=1,nz
               mnp(i,j,k)=0.0
 10            continue

      do 20 l=1,Ni_tot

         i=ijkp(l,1)!+wquad(l,1)
         j=ijkp(l,2)!+wquad(l,2)
         k=ijkp(l,3)!+wquad(l,3)

c         if (i .lt. 1) i = nx-1
c         if (j .lt. 1) i = ny-1
c         if (k .lt. 1) i = nz-1

         ip = i+1
         jp = j+1
         kp = k+1

c         if (ip .ge. nx) ip = 2      !periodic boundary conditions
c         if (jp .ge. ny) jp = 2      !periodic boundary conditions
c         if (kp .ge. nz) kp = 2      !periodic boundary conditions

c         volb = dx*dy*(qz(k+1)-qz(k))*beta
c         volb = dx*dy*dz_cell(k)*beta
         volb = (qx(ip)-qx(i))*(qy(jp)-qy(j))*(qz(kp)-qz(k))*beta

         mnp(i,j,k) = mnp(i,j,k) + (wght(l,1)/mrat(l))/volb
         mnp(ip,j,k) = mnp(ip,j,k) + (wght(l,2)/mrat(l))/volb
         mnp(i,j,kp) = mnp(i,j,kp) + (wght(l,3)/mrat(l))/volb
         mnp(ip,j,kp) = mnp(ip,j,kp) + (wght(l,4)/mrat(l))/volb
         mnp(i,jp,k) = mnp(i,jp,k) + (wght(l,5)/mrat(l))/volb
         mnp(ip,jp,k) = mnp(ip,jp,k) + (wght(l,6)/mrat(l))/volb
         mnp(i,jp,kp) = mnp(i,jp,kp) + (wght(l,7)/mrat(l))/volb
         mnp(ip,jp,kp) = mnp(ip,jp,kp) + (wght(l,8)/mrat(l))/volb

 20      continue

         mnp(:,:,:) = mion*mnp(:,:,:) !mass density

c use for periodic boundary conditions
         mnp(nx-1,:,:) = mnp(nx-1,:,:)+mnp(1,:,:)
         mnp(:,ny-1,:) = mnp(:,ny-1,:)+mnp(:,1,:)
         mnp(:,:,nz-1) = mnp(:,:,nz-1)+mnp(:,:,1)


         call MPI_Barrier(MPI_COMM_WORLD,ierr)

         call MPI_ALLREDUCE(mnp(:,:,:),recvbuf,count,
     x        MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)

         mnp(:,:,:) = reshape(recvbuf,(/nx,ny,nz/))


      call periodic_scalar(mnp)

      return
      end SUBROUTINE update_rho
c----------------------------------------------------------------------



c----------------------------------------------------------------------
      SUBROUTINE update_mixed()
c Weight density to eight nearest grid points.
c----------------------------------------------------------------------
      !!include 'incurv.h'

      real recvbuf(nx*ny*nz)
      integer count
      real mix_cnt(nx,ny,nz)
      count = nx*ny*nz


      do 10 i=1,nx
         do 10 j=1,ny
            do 10 k=1,nz
               mixed(i,j,k) = 0.0
               mix_cnt(i,j,k) = 0
 10            continue

      do 20 l=1,Ni_tot

         i=ijkp(l,1)
         j=ijkp(l,2)
         k=ijkp(l,3)

         mixed(i,j,k) = mixed(i,j,k) + mix_ind(l)
         mix_cnt(i,j,k) = mix_cnt(i,j,k) + 1

 20      continue

         call MPI_Barrier(MPI_COMM_WORLD,ierr)

         call MPI_ALLREDUCE(mixed(:,:,:),recvbuf,count,
     x        MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)

         mixed(:,:,:) = reshape(recvbuf,(/nx,ny,nz/))


         call MPI_Barrier(MPI_COMM_WORLD,ierr)

         call MPI_ALLREDUCE(mix_cnt(:,:,:),recvbuf,count,
     x        MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)

         mix_cnt(:,:,:) = reshape(recvbuf,(/nx,ny,nz/))

         mixed(:,:,:) = mixed(:,:,:)/mix_cnt(:,:,:)

      return
      end SUBROUTINE update_mixed
c----------------------------------------------------------------------



c----------------------------------------------------------------------
      SUBROUTINE update_up(vp,np,up)
c----------------------------------------------------------------------
      !!include 'incurv.h'

      real vp(Ni_max,3),
     x     np(nx,ny,nz),
     x     up(nx,ny,nz,3)

      real volb,nvolb      !np times vol times beta

      real recvbuf(nx*ny*nz)
      integer count
      count = nx*ny*nz


      do 10 m=1,3          !clear out temp variable ct
         do 10 i=1,nx
            do 10 j=1,ny
               do 10 k=1,nz
                  up(i,j,k,m)=0.0
                  ct(i,j,k,m)=0.0
 10               continue

c      where (ijkp(:,1) .eq. 1)
c         wquad(:,1) = 0.0
c      endwhere
         
      do 20 l=1,Ni_tot

         i=ijkp(l,1)!+wquad(l,1)
         j=ijkp(l,2)!+wquad(l,2)
         k=ijkp(l,3)!+wquad(l,3)

         ip = i+1
         jp = j+1
         kp = k+1

c         if (ip .ge. nx) ip = 2      !periodic boundary conditions
c         if (jp .ge. ny) jp = 2      !periodic boundary conditions
c         if (kp .ge. nz) kp = 2      !periodic boundary conditions


c         volb = dx*dy*(qz(kp)-qz(k))*beta
         volb = (qx(ip)-qx(i))*(qy(jp)-qy(j))*(qz(kp)-qz(k))*beta
c         volb = dx_grid(i)*dy_grid(j)*dz_grid(k)*beta
c         volb = dx*dy*dz_cell(k)*beta
      
c         np_at_Ba = np(i,j,k)*wght(l,1) + np(ip,j,k)*wght(l,2) + 
c     x              np(i,j,kp)*wght(l,3) + np(ip,j,kp)*wght(l,4) + 
c     x              np(i,jp,k)*wght(l,5) + np(ip,jp,k)*wght(l,6) +
c     x              np(i,jp,kp)*wght(l,7) + np(ip,jp,kp)*wght(l,8)


c         if (np(i,j,k) .gt. 0.0) then
         nvolb = 1.0/(np(i,j,k)*volb)
         ct(i,j,k,1) = ct(i,j,k,1) + vp(l,1)*wght(l,1)*nvolb
         ct(i,j,k,2) = ct(i,j,k,2) + vp(l,2)*wght(l,1)*nvolb
         ct(i,j,k,3) = ct(i,j,k,3) + vp(l,3)*wght(l,1)*nvolb
c         endif

c         if (np(ip,j,k) .gt. 0.0) then
         nvolb = 1.0/(np(ip,j,k)*volb)
         ct(ip,j,k,1) = ct(ip,j,k,1) + vp(l,1)*wght(l,2)*nvolb
         ct(ip,j,k,2) = ct(ip,j,k,2) + vp(l,2)*wght(l,2)*nvolb
         ct(ip,j,k,3) = ct(ip,j,k,3) + vp(l,3)*wght(l,2)*nvolb
c         endif

c         if (np(i,j,kp) .gt. 0.0) then
         nvolb = 1.0/(np(i,j,kp)*volb)
         ct(i,j,kp,1) = ct(i,j,kp,1) + vp(l,1)*wght(l,3)*nvolb
         ct(i,j,kp,2) = ct(i,j,kp,2) + vp(l,2)*wght(l,3)*nvolb
         ct(i,j,kp,3) = ct(i,j,kp,3) + vp(l,3)*wght(l,3)*nvolb
c         endif

c         if (np(ip,j,kp) .gt. 0.0) then
         nvolb = 1.0/(np(ip,j,kp)*volb)
         ct(ip,j,kp,1) = ct(ip,j,kp,1) + vp(l,1)*wght(l,4)*nvolb
         ct(ip,j,kp,2) = ct(ip,j,kp,2) + vp(l,2)*wght(l,4)*nvolb
         ct(ip,j,kp,3) = ct(ip,j,kp,3) + vp(l,3)*wght(l,4)*nvolb
c         endif

c         if (np(i,jp,k) .gt. 0.0) then
         nvolb = 1.0/(np(i,jp,k)*volb)
         ct(i,jp,k,1) = ct(i,jp,k,1) + vp(l,1)*wght(l,5)*nvolb
         ct(i,jp,k,2) = ct(i,jp,k,2) + vp(l,2)*wght(l,5)*nvolb
         ct(i,jp,k,3) = ct(i,jp,k,3) + vp(l,3)*wght(l,5)*nvolb
c         endif

c         if (np(ip,jp,k) .gt. 0.0) then
         nvolb = 1.0/(np(ip,jp,k)*volb)
         ct(ip,jp,k,1) = ct(ip,jp,k,1) + vp(l,1)*wght(l,6)*nvolb
         ct(ip,jp,k,2) = ct(ip,jp,k,2) + vp(l,2)*wght(l,6)*nvolb
         ct(ip,jp,k,3) = ct(ip,jp,k,3) + vp(l,3)*wght(l,6)*nvolb
c         endif

c         if (np(i,jp,kp) .gt. 0.0) then
         nvolb = 1.0/(np(i,jp,kp)*volb)
         ct(i,jp,kp,1) = ct(i,jp,kp,1) + vp(l,1)*wght(l,7)*nvolb
         ct(i,jp,kp,2) = ct(i,jp,kp,2) + vp(l,2)*wght(l,7)*nvolb
         ct(i,jp,kp,3) = ct(i,jp,kp,3) + vp(l,3)*wght(l,7)*nvolb
c         endif

c         if (np(ip,jp,kp) .gt. 0.0) then
         nvolb = 1.0/(np(ip,jp,kp)*volb)
         ct(ip,jp,kp,1) = ct(ip,jp,kp,1) + vp(l,1)*wght(l,8)*nvolb
         ct(ip,jp,kp,2) = ct(ip,jp,kp,2) + vp(l,2)*wght(l,8)*nvolb
         ct(ip,jp,kp,3) = ct(ip,jp,kp,3) + vp(l,3)*wght(l,8)*nvolb
c         endif


c         nvolb = np(i,j,k)*volb
c         ct(i,j,k,1) = ct(i,j,k,1) + vp(l,1)*wght(l,1)/nvolb
c         ct(i,j,k,2) = ct(i,j,k,2) + vp(l,2)*wght(l,1)/nvolb
c         ct(i,j,k,3) = ct(i,j,k,3) + vp(l,3)*wght(l,1)/nvolb

c         nvolb = np(ip,j,k)*volb
c         ct(ip,j,k,1) = ct(ip,j,k,1) + vp(l,1)*wght(l,2)/nvolb
c         ct(ip,j,k,2) = ct(ip,j,k,2) + vp(l,2)*wght(l,2)/nvolb
c         ct(ip,j,k,3) = ct(ip,j,k,3) + vp(l,3)*wght(l,2)/nvolb

c         nvolb = np(i,j,kp)*volb
c         ct(i,j,kp,1) = ct(i,j,kp,1) + vp(l,1)*wght(l,3)/nvolb
c         ct(i,j,kp,2) = ct(i,j,kp,2) + vp(l,2)*wght(l,3)/nvolb
c         ct(i,j,kp,3) = ct(i,j,kp,3) + vp(l,3)*wght(l,3)/nvolb

c         nvolb = np(ip,j,kp)*volb
c         ct(ip,j,kp,1) = ct(ip,j,kp,1) + vp(l,1)*wght(l,4)/nvolb
c         ct(ip,j,kp,2) = ct(ip,j,kp,2) + vp(l,2)*wght(l,4)/nvolb
c         ct(ip,j,kp,3) = ct(ip,j,kp,3) + vp(l,3)*wght(l,4)/nvolb

c         nvolb = np(i,jp,k)*volb
c         ct(i,jp,k,1) = ct(i,jp,k,1) + vp(l,1)*wght(l,5)/nvolb
c         ct(i,jp,k,2) = ct(i,jp,k,2) + vp(l,2)*wght(l,5)/nvolb
c         ct(i,jp,k,3) = ct(i,jp,k,3) + vp(l,3)*wght(l,5)/nvolb

c         nvolb = np(ip,jp,k)*volb
c         ct(ip,jp,k,1) = ct(ip,jp,k,1) + vp(l,1)*wght(l,6)/nvolb
c         ct(ip,jp,k,2) = ct(ip,jp,k,2) + vp(l,2)*wght(l,6)/nvolb
c         ct(ip,jp,k,3) = ct(ip,jp,k,3) + vp(l,3)*wght(l,6)/nvolb

c         nvolb = np(i,jp,kp)*volb
c         ct(i,jp,kp,1) = ct(i,jp,kp,1) + vp(l,1)*wght(l,7)/nvolb
c         ct(i,jp,kp,2) = ct(i,jp,kp,2) + vp(l,2)*wght(l,7)/nvolb
c         ct(i,jp,kp,3) = ct(i,jp,kp,3) + vp(l,3)*wght(l,7)/nvolb

c         nvolb = np(ip,jp,kp)*volb
c         ct(ip,jp,kp,1) = ct(ip,jp,kp,1) + vp(l,1)*wght(l,8)/nvolb
c         ct(ip,jp,kp,2) = ct(ip,jp,kp,2) + vp(l,2)*wght(l,8)/nvolb
c         ct(ip,jp,kp,3) = ct(ip,jp,kp,3) + vp(l,3)*wght(l,8)/nvolb

 20   continue

c use for periodic boundary conditions
      ct(nx-1,:,:,:) = ct(nx-1,:,:,:)+ct(1,:,:,:)
      ct(:,ny-1,:,:) = ct(:,ny-1,:,:)+ct(:,1,:,:)
      ct(:,:,nz-1,:) = ct(:,:,nz-1,:)+ct(:,:,1,:)


      call periodic(ct)

      call MPI_Barrier(MPI_COMM_WORLD,ierr)
         
      call MPI_ALLREDUCE(ct(:,:,:,1),recvbuf,count,
     x     MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
      ct(:,:,:,1) = reshape(recvbuf,(/nx,ny,nz/))
      
      call MPI_ALLREDUCE(ct(:,:,:,2),recvbuf,count,
     x     MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
      ct(:,:,:,2) = reshape(recvbuf,(/nx,ny,nz/))
      
      call MPI_ALLREDUCE(ct(:,:,:,3),recvbuf,count,
     x     MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
      ct(:,:,:,3) = reshape(recvbuf,(/nx,ny,nz/))
      

      do 30 i=1,nx-1      !interpolate back to contravarient positions
         do 30 j=1,ny-1
            do 30 k=1,nz-1
c               up(i,j,k,1) = ct(i,j,k,1)
c               up(i,j,k,2) = ct(i,j,k,2)
c               up(i,j,k,3) = ct(i,j,k,3)

               up(i,j,k,1) = 0.5*(ct(i,j,k,1)+ct(i+1,j,k,1))
               up(i,j,k,2) = 0.5*(ct(i,j,k,2)+ct(i,j+1,k,2))
               up(i,j,k,3) = 0.5*(ct(i,j,k,3)+ct(i,j,k+1,3))

 30            continue


      call periodic(up)

      return
      end SUBROUTINE update_up
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE update_up_2(vp,np,up)
c----------------------------------------------------------------------
c      include 'incurv.h'

      real vp(Ni_max,3),
     x     np(nx,ny,nz),
     x     up(nx,ny,nz,3)

c      real volb,nvolb      !np times vol times beta

      real cnt(nx,ny,nz)

      real recvbuf(nx*ny*nz)
      integer count
      count = nx*ny*nz


c      integer dest, source
c      real out_buf_z(nx,ny,3)
c      real in_buf_z(nx,ny,3)
c      integer cnt_buf_z
c      integer reqs(2)
c      integer stats(MPI_STATUS_SIZE,2)

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

c      cnt_buf_z = nx*ny*3

c      real recvbuf(nx*ny*nz)
c      integer count
c      count = nx*ny*nz
    

      do 10 m=1,3          !clear out temp variable ct
         do 10 i=1,nx
            do 10 j=1,ny
               do 10 k=1,nz
                  up(i,j,k,m)=0.0
                  ct(i,j,k,m)=0.0
 10            continue
               
               cnt(:,:,:) = 0.0
               
c      where (ijkp(:,1) .eq. 1)
c         wquad(:,1) = 0.0
c      endwhere
     

    
      do 20 l=1,Ni_tot

         i=ijkp(l,1)!+wquad(l,1)
         j=ijkp(l,2)!+wquad(l,2)
         k=ijkp(l,3)!+wquad(l,3)

         ip = i+1
         jp = j+1
         kp = k+1

c         nvolb = 1.0
c         if (np(i,j,k) .gt. 0.0) then
c         nvolb = np(i,j,k)*volb
         cnt(i,j,k) = cnt(i,j,k) + wght(l,1)
         ct(i,j,k,1) = ct(i,j,k,1) + vp(l,1)*wght(l,1) 
         ct(i,j,k,2) = ct(i,j,k,2) + vp(l,2)*wght(l,1) 
         ct(i,j,k,3) = ct(i,j,k,3) + vp(l,3)*wght(l,1) 
         
c         endif

c         if (np(ip,j,k) .gt. 0.0) then
c         nvolb = np(ip,j,k)*volb
         cnt(ip,j,k) = cnt(ip,j,k) + wght(l,2)
         ct(ip,j,k,1) = ct(ip,j,k,1) + vp(l,1)*wght(l,2) 
         ct(ip,j,k,2) = ct(ip,j,k,2) + vp(l,2)*wght(l,2) 
         ct(ip,j,k,3) = ct(ip,j,k,3) + vp(l,3)*wght(l,2) 
c         endif

c         if (np(i,j,kp) .gt. 0.0) then
c         nvolb = np(i,j,kp)*volb
         cnt(i,j,kp) = cnt(i,j,kp) + wght(l,3)
         ct(i,j,kp,1) = ct(i,j,kp,1) + vp(l,1)*wght(l,3) 
         ct(i,j,kp,2) = ct(i,j,kp,2) + vp(l,2)*wght(l,3) 
         ct(i,j,kp,3) = ct(i,j,kp,3) + vp(l,3)*wght(l,3) 
c         endif

c         if (np(ip,j,kp) .gt. 0.0) then
c         nvolb = np(ip,j,kp)*volb
         cnt(ip,j,kp) = cnt(ip,j,kp) + wght(l,4) 
         ct(ip,j,kp,1) = ct(ip,j,kp,1) + vp(l,1)*wght(l,4) 
         ct(ip,j,kp,2) = ct(ip,j,kp,2) + vp(l,2)*wght(l,4) 
         ct(ip,j,kp,3) = ct(ip,j,kp,3) + vp(l,3)*wght(l,4) 
c         endif

c         if (np(i,jp,k) .gt. 0.0) then
c         nvolb = np(i,jp,k)*volb
         cnt(i,jp,k) = cnt(i,jp,k) + wght(l,5)
         ct(i,jp,k,1) = ct(i,jp,k,1) + vp(l,1)*wght(l,5) 
         ct(i,jp,k,2) = ct(i,jp,k,2) + vp(l,2)*wght(l,5) 
         ct(i,jp,k,3) = ct(i,jp,k,3) + vp(l,3)*wght(l,5) 
c         endif

c         if (np(ip,jp,k) .gt. 0.0) then
c         nvolb = np(ip,jp,k)*volb
         cnt(ip,jp,k) = cnt(ip,jp,k) + wght(l,6)
         ct(ip,jp,k,1) = ct(ip,jp,k,1) + vp(l,1)*wght(l,6) 
         ct(ip,jp,k,2) = ct(ip,jp,k,2) + vp(l,2)*wght(l,6) 
         ct(ip,jp,k,3) = ct(ip,jp,k,3) + vp(l,3)*wght(l,6) 
c         endif

c         if (np(i,jp,kp) .gt. 0.0) then
c         nvolb = np(i,jp,kp)*volb
         cnt(i,jp,kp) = cnt(i,jp,kp) + wght(l,7)
         ct(i,jp,kp,1) = ct(i,jp,kp,1) + vp(l,1)*wght(l,7) 
         ct(i,jp,kp,2) = ct(i,jp,kp,2) + vp(l,2)*wght(l,7) 
         ct(i,jp,kp,3) = ct(i,jp,kp,3) + vp(l,3)*wght(l,7) 
c         endif

c         if (np(ip,jp,kp) .gt. 0.0) then
c         nvolb = np(ip,jp,kp)*volb
         cnt(ip,jp,kp) = cnt(ip,jp,kp) + wght(l,8)
         ct(ip,jp,kp,1) = ct(ip,jp,kp,1) + vp(l,1)*wght(l,8) 
         ct(ip,jp,kp,2) = ct(ip,jp,kp,2) + vp(l,2)*wght(l,8) 
         ct(ip,jp,kp,3) = ct(ip,jp,kp,3) + vp(l,3)*wght(l,8) 
c         endif

 20   continue


c use for periodic boundary conditions
      ct(nx-1,:,:,:) = ct(nx-1,:,:,:)+ct(1,:,:,:)
      cnt(nx-1,:,:) = cnt(nx-1,:,:)+cnt(1,:,:)

      ct(:,ny-1,:,:) = ct(:,ny-1,:,:)+ct(:,1,:,:)
      cnt(:,ny-1,:) = cnt(:,ny-1,:)+cnt(:,1,:)

      ct(:,:,nz-1,:) = ct(:,:,nz-1,:)+ct(:,:,1,:)
      cnt(:,:,nz-1) = cnt(:,:,nz-1)+cnt(:,:,1)

      call MPI_Barrier(MPI_COMM_WORLD,ierr)

c      call periodic(ct)
c      call periodic_scalar(cnt)

      where (cnt(:,:,:) .gt. 0.0)
         ct(:,:,:,1) = ct(:,:,:,1)/cnt(:,:,:)/procnum
         ct(:,:,:,2) = ct(:,:,:,2)/cnt(:,:,:)/procnum
         ct(:,:,:,3) = ct(:,:,:,3)/cnt(:,:,:)/procnum
      endwhere


c      call periodic(ct)


      call MPI_Barrier(MPI_COMM_WORLD,ierr)

         
      call MPI_ALLREDUCE(ct(:,:,:,1),recvbuf,count,
     x     MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
      ct(:,:,:,1) = reshape(recvbuf,(/nx,ny,nz/))
      
      call MPI_ALLREDUCE(ct(:,:,:,2),recvbuf,count,
     x     MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
      ct(:,:,:,2) = reshape(recvbuf,(/nx,ny,nz/))
      
      call MPI_ALLREDUCE(ct(:,:,:,3),recvbuf,count,
     x     MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
      ct(:,:,:,3) = reshape(recvbuf,(/nx,ny,nz/))

      call periodic(ct)

      do 30 i=1,nx-1      !interpolate back to contravarient positions
         do 30 j=1,ny-1
            do 30 k=1,nz-1
               up(i,j,k,1) = 0.5*(ct(i,j,k,1)+ct(i+1,j,k,1))
               up(i,j,k,2) = 0.5*(ct(i,j,k,2)+ct(i,j+1,k,2))
               up(i,j,k,3) = 0.5*(ct(i,j,k,3)+ct(i,j,k+1,3))

 30            continue

      call periodic(up)

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      return
      end SUBROUTINE update_up_2
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE separate_np(np,flg)
c Weight density to eight nearest grid points.
c----------------------------------------------------------------------
      !!include 'incurv.h'

      real np(nx,ny,nz)
      integer flg(Ni_max)

      real volb              !cell volume times beta
      real recvbuf(nx*ny*nz)
      integer count
      count = nx*ny*nz

c      real sumnp,vol

      do 10 i=1,nx
         do 10 j=1,ny
            do 10 k=1,nz
               np(i,j,k)=0.0
 10            continue

      do 20 l=1,Ni_tot

         i=ijkp(l,1)!+wquad(l,1)
         j=ijkp(l,2)!+wquad(l,2)
         k=ijkp(l,3)!+wquad(l,3)

c         if (i .lt. 1) i = nx-1
c         if (j .lt. 1) i = ny-1
c         if (k .lt. 1) i = nz-1

         ip = i+1
         jp = j+1
         kp = k+1

c         if (ip .ge. nx) ip = 2      !periodic boundary conditions
c         if (jp .ge. ny) jp = 2      !periodic boundary conditions
c         if (kp .ge. nz) kp = 2      !periodic boundary conditions

c         volb = dx*dy*(qz(k+1)-qz(k))*beta
c         volb = dx*dy*dz_cell(k)*beta
         volb = (qx(ip)-qx(i))*(qy(jp)-qy(j))*(qz(kp)-qz(k))*beta


         np(i,j,k) = np(i,j,k) + flg(l)*wght(l,1)/volb
         np(ip,j,k) = np(ip,j,k) + flg(l)*wght(l,2)/volb
         np(i,j,kp) = np(i,j,kp) + flg(l)*wght(l,3)/volb
         np(ip,j,kp) = np(ip,j,kp) + flg(l)*wght(l,4)/volb
         np(i,jp,k) = np(i,jp,k) + flg(l)*wght(l,5)/volb
         np(ip,jp,k) = np(ip,jp,k) + flg(l)*wght(l,6)/volb
         np(i,jp,kp) = np(i,jp,kp) + flg(l)*wght(l,7)/volb
         np(ip,jp,kp) = np(ip,jp,kp) + flg(l)*wght(l,8)/volb

 20      continue

         call MPI_Barrier(MPI_COMM_WORLD,ierr)

         call MPI_ALLREDUCE(np(:,:,:),recvbuf,count,
     x        MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)

         np(:,:,:) = reshape(recvbuf,(/nx,ny,nz/))

c use for periodic boundary conditions
      np(nx-1,:,:) = np(nx-1,:,:)+np(1,:,:)
      np(:,ny-1,:) = np(:,ny-1,:)+np(:,1,:)
      np(:,:,nz-1) = np(:,:,nz-1)+np(:,:,1)


      call periodic_scalar(np)

      return
      end SUBROUTINE separate_np
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE separate_up(vp,np,flg,up)
c----------------------------------------------------------------------
      !!include 'incurv.h'

      real vp(Ni_max,3),
     x     np(nx,ny,nz)
      integer flg(Ni_max)
      real up(nx,ny,nz,3)



      real volb,nvolb      !np times vol times beta

      real recvbuf(nx*ny*nz)
      integer count
      count = nx*ny*nz


      do 10 m=1,3          !clear out temp variable ct
         do 10 i=1,nx
            do 10 j=1,ny
               do 10 k=1,nz
                  up(i,j,k,m)=0.0
                  ct(i,j,k,m)=0.0
 10               continue

c      where (ijkp(:,1) .eq. 1)
c         wquad(:,1) = 0.0
c      endwhere
         
      do 20 l=1,Ni_tot

         i=ijkp(l,1)!+wquad(l,1)
         j=ijkp(l,2)!+wquad(l,2)
         k=ijkp(l,3)!+wquad(l,3)

         ip = i+1
         jp = j+1
         kp = k+1

c         volb = dx*dy*dz_cell(k)*beta
         volb = (qx(ip)-qx(i))*(qy(jp)-qy(j))*(qz(kp)-qz(k))*beta
      
         if (np(i,j,k) .gt. 0.0) then
         nvolb = np(i,j,k)*volb
         ct(i,j,k,1) = ct(i,j,k,1) + flg(l)*vp(l,1)*wght(l,1)/nvolb
         ct(i,j,k,2) = ct(i,j,k,2) + flg(l)*vp(l,2)*wght(l,1)/nvolb
         ct(i,j,k,3) = ct(i,j,k,3) + flg(l)*vp(l,3)*wght(l,1)/nvolb
         endif

         if (np(ip,j,k) .gt. 0.0) then
         nvolb = np(ip,j,k)*volb
         ct(ip,j,k,1) = ct(ip,j,k,1) + flg(l)*vp(l,1)*wght(l,2)/nvolb
         ct(ip,j,k,2) = ct(ip,j,k,2) + flg(l)*vp(l,2)*wght(l,2)/nvolb
         ct(ip,j,k,3) = ct(ip,j,k,3) + flg(l)*vp(l,3)*wght(l,2)/nvolb
         endif

         if (np(i,j,kp) .gt. 0.0) then
         nvolb = np(i,j,kp)*volb
         ct(i,j,kp,1) = ct(i,j,kp,1) + flg(l)*vp(l,1)*wght(l,3)/nvolb
         ct(i,j,kp,2) = ct(i,j,kp,2) + flg(l)*vp(l,2)*wght(l,3)/nvolb
         ct(i,j,kp,3) = ct(i,j,kp,3) + flg(l)*vp(l,3)*wght(l,3)/nvolb
         endif

         if (np(ip,j,kp) .gt. 0.0) then
         nvolb = np(ip,j,kp)*volb
         ct(ip,j,kp,1) = ct(ip,j,kp,1) + flg(l)*vp(l,1)*wght(l,4)/nvolb
         ct(ip,j,kp,2) = ct(ip,j,kp,2) + flg(l)*vp(l,2)*wght(l,4)/nvolb
         ct(ip,j,kp,3) = ct(ip,j,kp,3) + flg(l)*vp(l,3)*wght(l,4)/nvolb
         endif

         if (np(i,jp,k) .gt. 0.0) then
         nvolb = np(i,jp,k)*volb
         ct(i,jp,k,1) = ct(i,jp,k,1) + flg(l)*vp(l,1)*wght(l,5)/nvolb
         ct(i,jp,k,2) = ct(i,jp,k,2) + flg(l)*vp(l,2)*wght(l,5)/nvolb
         ct(i,jp,k,3) = ct(i,jp,k,3) + flg(l)*vp(l,3)*wght(l,5)/nvolb
         endif

         if (np(ip,jp,k) .gt. 0.0) then
         nvolb = np(ip,jp,k)*volb
         ct(ip,jp,k,1) = ct(ip,jp,k,1) + flg(l)*vp(l,1)*wght(l,6)/nvolb
         ct(ip,jp,k,2) = ct(ip,jp,k,2) + flg(l)*vp(l,2)*wght(l,6)/nvolb
         ct(ip,jp,k,3) = ct(ip,jp,k,3) + flg(l)*vp(l,3)*wght(l,6)/nvolb
         endif

         if (np(i,jp,kp) .gt. 0.0) then
         nvolb = np(i,jp,kp)*volb
         ct(i,jp,kp,1) = ct(i,jp,kp,1) + flg(l)*vp(l,1)*wght(l,7)/nvolb
         ct(i,jp,kp,2) = ct(i,jp,kp,2) + flg(l)*vp(l,2)*wght(l,7)/nvolb
         ct(i,jp,kp,3) = ct(i,jp,kp,3) + flg(l)*vp(l,3)*wght(l,7)/nvolb
         endif

         if (np(ip,jp,kp) .gt. 0.0) then
         nvolb = np(ip,jp,kp)*volb
        ct(ip,jp,kp,1) = ct(ip,jp,kp,1) + flg(l)*vp(l,1)*wght(l,8)/nvolb
        ct(ip,jp,kp,2) = ct(ip,jp,kp,2) + flg(l)*vp(l,2)*wght(l,8)/nvolb
        ct(ip,jp,kp,3) = ct(ip,jp,kp,3) + flg(l)*vp(l,3)*wght(l,8)/nvolb
         endif

 20   continue

c use for periodic boundary conditions
      ct(nx-1,:,:,:) = ct(nx-1,:,:,:)+ct(1,:,:,:)
      ct(:,ny-1,:,:) = ct(:,ny-1,:,:)+ct(:,1,:,:)
      ct(:,:,nz-1,:) = ct(:,:,nz-1,:)+ct(:,:,1,:)


      call periodic(ct)

      call MPI_Barrier(MPI_COMM_WORLD,ierr)
         
      call MPI_ALLREDUCE(ct(:,:,:,1),recvbuf,count,
     x     MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
      ct(:,:,:,1) = reshape(recvbuf,(/nx,ny,nz/))
      
      call MPI_ALLREDUCE(ct(:,:,:,2),recvbuf,count,
     x     MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
      ct(:,:,:,2) = reshape(recvbuf,(/nx,ny,nz/))
      
      call MPI_ALLREDUCE(ct(:,:,:,3),recvbuf,count,
     x     MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
      ct(:,:,:,3) = reshape(recvbuf,(/nx,ny,nz/))
      

      do 30 i=1,nx-1      !interpolate back to contravarient positions
         do 30 j=1,ny-1
            do 30 k=1,nz-1

               up(i,j,k,1) = 0.5*(ct(i,j,k,1)+ct(i+1,j,k,1))
               up(i,j,k,2) = 0.5*(ct(i,j,k,2)+ct(i,j+1,k,2))
               up(i,j,k,3) = 0.5*(ct(i,j,k,3)+ct(i,j,k+1,3))

 30            continue


      call periodic(up)

      return
      end SUBROUTINE separate_up
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE get_temperature(xp,vp,np,temp_p)
c----------------------------------------------------------------------
      !!include 'incurv.h'
      
      real xp(Ni_max,3),
     x     vp(Ni_max,3),
     x     np(nx,ny,nz),
     x     up_ave(nx,ny,nz,3),
     x     up2(nx,ny,nz,3),
     x     temp_p(nx,ny,nz)


      real mvp(Ni_max,3)

      real volb,nvolb      !np times vol times beta

      real recvbuf(nx*ny*nz)
      integer count
      count = nx*ny*nz


c      do 10 m=1,3          !clear out temp variable ct
c         do 10 i=1,nx
c            do 10 j=1,ny
c               do 10 k=1,nz
c                  up2(i,j,k,m)=0.0
c                  ct(i,j,k,m)=0.0
c 10               continue
      up2(:,:,:,:) = 0.0
      up_ave(:,:,:,:) = 0.0
      ct(:,:,:,:) = 0.0


      do m = 1,3 
         mvp(:,m) = vp(:,m)/sqrt(mrat(:))
      enddo


c      do m = 1,3 
c         mvp(:,m) = vp(:,m)/mrat(:)
c      enddo

      do 20 l=1,Ni_tot

         i=ijkp(l,1)!+wquad(l,1)
         j=ijkp(l,2)!+wquad(l,2)
         k=ijkp(l,3)!+wquad(l,3)

         ip = i+1
         jp = j+1
         kp = k+1

c         volb = dx*dy*dz_cell(k)*beta
         volb = (qx(ip)-qx(i))*(qy(jp)-qy(j))*(qz(kp)-qz(k))*beta
      
         if (np(i,j,k) .gt. 0.0) then
         nvolb = np(i,j,k)*volb
         ct(i,j,k,1) = ct(i,j,k,1) + mvp(l,1)**2*wght(l,1)/nvolb  
         ct(i,j,k,2) = ct(i,j,k,2) + mvp(l,2)**2*wght(l,1)/nvolb
         ct(i,j,k,3) = ct(i,j,k,3) + mvp(l,3)**2*wght(l,1)/nvolb
         endif

         if (np(ip,j,k) .gt. 0.0) then
         nvolb = np(ip,j,k)*volb
         ct(ip,j,k,1) = ct(ip,j,k,1) + mvp(l,1)**2*wght(l,2)/nvolb
         ct(ip,j,k,2) = ct(ip,j,k,2) + mvp(l,2)**2*wght(l,2)/nvolb
         ct(ip,j,k,3) = ct(ip,j,k,3) + mvp(l,3)**2*wght(l,2)/nvolb
         endif

         if (np(i,j,kp) .gt. 0.0) then
         nvolb = np(i,j,kp)*volb
         ct(i,j,kp,1) = ct(i,j,kp,1) + mvp(l,1)**2*wght(l,3)/nvolb
         ct(i,j,kp,2) = ct(i,j,kp,2) + mvp(l,2)**2*wght(l,3)/nvolb
         ct(i,j,kp,3) = ct(i,j,kp,3) + mvp(l,3)**2*wght(l,3)/nvolb
         endif

         if (np(ip,j,kp) .gt. 0.0) then
         nvolb = np(ip,j,kp)*volb
         ct(ip,j,kp,1) = ct(ip,j,kp,1) + mvp(l,1)**2*wght(l,4)/nvolb
         ct(ip,j,kp,2) = ct(ip,j,kp,2) + mvp(l,2)**2*wght(l,4)/nvolb
         ct(ip,j,kp,3) = ct(ip,j,kp,3) + mvp(l,3)**2*wght(l,4)/nvolb
         endif

         if (np(i,jp,k) .gt. 0.0) then
         nvolb = np(i,jp,k)*volb
         ct(i,jp,k,1) = ct(i,jp,k,1) + mvp(l,1)**2*wght(l,5)/nvolb
         ct(i,jp,k,2) = ct(i,jp,k,2) + mvp(l,2)**2*wght(l,5)/nvolb
         ct(i,jp,k,3) = ct(i,jp,k,3) + mvp(l,3)**2*wght(l,5)/nvolb
         endif

         if (np(ip,jp,k) .gt. 0.0) then
         nvolb = np(ip,jp,k)*volb
         ct(ip,jp,k,1) = ct(ip,jp,k,1) + mvp(l,1)**2*wght(l,6)/nvolb
         ct(ip,jp,k,2) = ct(ip,jp,k,2) + mvp(l,2)**2*wght(l,6)/nvolb
         ct(ip,jp,k,3) = ct(ip,jp,k,3) + mvp(l,3)**2*wght(l,6)/nvolb
         endif

         if (np(i,jp,kp) .gt. 0.0) then
         nvolb = np(i,jp,kp)*volb
         ct(i,jp,kp,1) = ct(i,jp,kp,1) + mvp(l,1)**2*wght(l,7)/nvolb
         ct(i,jp,kp,2) = ct(i,jp,kp,2) + mvp(l,2)**2*wght(l,7)/nvolb
         ct(i,jp,kp,3) = ct(i,jp,kp,3) + mvp(l,3)**2*wght(l,7)/nvolb
         endif

         if (np(ip,jp,kp) .gt. 0.0) then
         nvolb = np(ip,jp,kp)*volb
         ct(ip,jp,kp,1) = ct(ip,jp,kp,1) + mvp(l,1)**2*wght(l,8)/nvolb
         ct(ip,jp,kp,2) = ct(ip,jp,kp,2) + mvp(l,2)**2*wght(l,8)/nvolb
         ct(ip,jp,kp,3) = ct(ip,jp,kp,3) + mvp(l,3)**2*wght(l,8)/nvolb
         endif


 20   continue

c use for periodic boundary conditions
      ct(nx-1,:,:,:) = ct(nx-1,:,:,:)+ct(1,:,:,:)
      ct(:,ny-1,:,:) = ct(:,ny-1,:,:)+ct(:,1,:,:)
      ct(:,:,nz-1,:) = ct(:,:,nz-1,:)+ct(:,:,1,:)


      call periodic(ct)

      call MPI_Barrier(MPI_COMM_WORLD,ierr)
         
      call MPI_ALLREDUCE(ct(:,:,:,1),recvbuf,count,
     x     MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
      ct(:,:,:,1) = reshape(recvbuf,(/nx,ny,nz/))
      
      call MPI_ALLREDUCE(ct(:,:,:,2),recvbuf,count,
     x     MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
      ct(:,:,:,2) = reshape(recvbuf,(/nx,ny,nz/))
      
      call MPI_ALLREDUCE(ct(:,:,:,3),recvbuf,count,
     x     MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
      ct(:,:,:,3) = reshape(recvbuf,(/nx,ny,nz/))
      

      do 30 i=1,nx-1      !interpolate back to contravarient positions
         do 30 j=1,ny-1
            do 30 k=1,nz-1

               up2(i,j,k,1) = 0.5*(ct(i,j,k,1)+ct(i+1,j,k,1))
               up2(i,j,k,2) = 0.5*(ct(i,j,k,2)+ct(i,j+1,k,2))
               up2(i,j,k,3) = 0.5*(ct(i,j,k,3)+ct(i,j,k+1,3))

 30            continue


      call periodic(up2)


      ct(:,:,:,:) = 0.0

      do 40 l=1,Ni_tot

         i=ijkp(l,1)!+wquad(l,1)
         j=ijkp(l,2)!+wquad(l,2)
         k=ijkp(l,3)!+wquad(l,3)

         ip = i+1
         jp = j+1
         kp = k+1

c         volb = dx*dy*dz_cell(k)*beta
         volb = (qx(ip)-qx(i))*(qy(jp)-qy(j))*(qz(kp)-qz(k))*beta
      
         if (np(i,j,k) .gt. 0.0) then
         nvolb = np(i,j,k)*volb
         ct(i,j,k,1) = ct(i,j,k,1) + mvp(l,1)*wght(l,1)/nvolb  
         ct(i,j,k,2) = ct(i,j,k,2) + mvp(l,2)*wght(l,1)/nvolb
         ct(i,j,k,3) = ct(i,j,k,3) + mvp(l,3)*wght(l,1)/nvolb
         endif

         if (np(ip,j,k) .gt. 0.0) then
         nvolb = np(ip,j,k)*volb
         ct(ip,j,k,1) = ct(ip,j,k,1) + mvp(l,1)*wght(l,2)/nvolb
         ct(ip,j,k,2) = ct(ip,j,k,2) + mvp(l,2)*wght(l,2)/nvolb
         ct(ip,j,k,3) = ct(ip,j,k,3) + mvp(l,3)*wght(l,2)/nvolb
         endif

         if (np(i,j,kp) .gt. 0.0) then
         nvolb = np(i,j,kp)*volb
         ct(i,j,kp,1) = ct(i,j,kp,1) + mvp(l,1)*wght(l,3)/nvolb
         ct(i,j,kp,2) = ct(i,j,kp,2) + mvp(l,2)*wght(l,3)/nvolb
         ct(i,j,kp,3) = ct(i,j,kp,3) + mvp(l,3)*wght(l,3)/nvolb
         endif

         if (np(ip,j,kp) .gt. 0.0) then
         nvolb = np(ip,j,kp)*volb
         ct(ip,j,kp,1) = ct(ip,j,kp,1) + mvp(l,1)*wght(l,4)/nvolb
         ct(ip,j,kp,2) = ct(ip,j,kp,2) + mvp(l,2)*wght(l,4)/nvolb
         ct(ip,j,kp,3) = ct(ip,j,kp,3) + mvp(l,3)*wght(l,4)/nvolb
         endif

         if (np(i,jp,k) .gt. 0.0) then
         nvolb = np(i,jp,k)*volb
         ct(i,jp,k,1) = ct(i,jp,k,1) + mvp(l,1)*wght(l,5)/nvolb
         ct(i,jp,k,2) = ct(i,jp,k,2) + mvp(l,2)*wght(l,5)/nvolb
         ct(i,jp,k,3) = ct(i,jp,k,3) + mvp(l,3)*wght(l,5)/nvolb
         endif

         if (np(ip,jp,k) .gt. 0.0) then
         nvolb = np(ip,jp,k)*volb
         ct(ip,jp,k,1) = ct(ip,jp,k,1) + mvp(l,1)*wght(l,6)/nvolb
         ct(ip,jp,k,2) = ct(ip,jp,k,2) + mvp(l,2)*wght(l,6)/nvolb
         ct(ip,jp,k,3) = ct(ip,jp,k,3) + mvp(l,3)*wght(l,6)/nvolb
         endif

         if (np(i,jp,kp) .gt. 0.0) then
         nvolb = np(i,jp,kp)*volb
         ct(i,jp,kp,1) = ct(i,jp,kp,1) + mvp(l,1)*wght(l,7)/nvolb
         ct(i,jp,kp,2) = ct(i,jp,kp,2) + mvp(l,2)*wght(l,7)/nvolb
         ct(i,jp,kp,3) = ct(i,jp,kp,3) + mvp(l,3)*wght(l,7)/nvolb
         endif

         if (np(ip,jp,kp) .gt. 0.0) then
         nvolb = np(ip,jp,kp)*volb
         ct(ip,jp,kp,1) = ct(ip,jp,kp,1) + mvp(l,1)*wght(l,8)/nvolb
         ct(ip,jp,kp,2) = ct(ip,jp,kp,2) + mvp(l,2)*wght(l,8)/nvolb
         ct(ip,jp,kp,3) = ct(ip,jp,kp,3) + mvp(l,3)*wght(l,8)/nvolb
         endif


 40   continue

c use for periodic boundary conditions
      ct(nx-1,:,:,:) = ct(nx-1,:,:,:)+ct(1,:,:,:)
      ct(:,ny-1,:,:) = ct(:,ny-1,:,:)+ct(:,1,:,:)
      ct(:,:,nz-1,:) = ct(:,:,nz-1,:)+ct(:,:,1,:)


      call periodic(ct)

      call MPI_Barrier(MPI_COMM_WORLD,ierr)
         
      call MPI_ALLREDUCE(ct(:,:,:,1),recvbuf,count,
     x     MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
      ct(:,:,:,1) = reshape(recvbuf,(/nx,ny,nz/))
      
      call MPI_ALLREDUCE(ct(:,:,:,2),recvbuf,count,
     x     MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
      ct(:,:,:,2) = reshape(recvbuf,(/nx,ny,nz/))
      
      call MPI_ALLREDUCE(ct(:,:,:,3),recvbuf,count,
     x     MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
      ct(:,:,:,3) = reshape(recvbuf,(/nx,ny,nz/))
      

      do 50 i=1,nx-1      !interpolate back to contravarient positions
         do 50 j=1,ny-1
            do 50 k=1,nz-1

               up_ave(i,j,k,1) = 0.5*(ct(i,j,k,1)+ct(i+1,j,k,1))
               up_ave(i,j,k,2) = 0.5*(ct(i,j,k,2)+ct(i,j+1,k,2))
               up_ave(i,j,k,3) = 0.5*(ct(i,j,k,3)+ct(i,j,k+1,3))

 50            continue


      call periodic(up_ave)


      do i = 1,nx
         do j = 1,ny
            do k = 1,nz
               temp_p(i,j,k) = (1./3.)*1e6*mion*(
     x                  sqrt((up2(i,j,k,1) - up_ave(i,j,k,1)**2)**2 +
     x                       (up2(i,j,k,2) - up_ave(i,j,k,2)**2)**2 + 
     x                       (up2(i,j,k,3) - up_ave(i,j,k,3)**2)**2))  
c     x                             (up_ave(i,j,k,1)**2 +
c     x                              up_ave(i,j,k,2)**2 +
c     x                              up_ave(i,j,k,3)**2))
c               if ((temp_p(i,j,k)/1.6e-19 .gt. 1.0) .and. 
c     x            (temp_p(i,j,k)/1.6e-19 .lt. 500.0)) then
c               endif
            enddo
         enddo
      enddo
      


      return
      end SUBROUTINE get_temperature
c----------------------------------------------------------------------

c----------------------------------------------------------------------
      SUBROUTINE check_index()
c----------------------------------------------------------------------
      !!include 'incurv.h'

      do l = 1, Ni_tot
         if (ijkp(l,1) .gt. nx) then
            write(*,*) 'i OB...',ijkp(l,:),m_arr(l)
         endif

         if (ijkp(l,1) .lt. 1) then
            write(*,*) 'i OB...',ijkp(l,:),m_arr(l)
         endif


         if (ijkp(l,2) .gt. ny) then
            write(*,*) 'j OB...',ijkp(l,:),m_arr(l)
         endif

         if (ijkp(l,2) .lt. 1) then
            write(*,*) 'j OB...',ijkp(l,:),m_arr(l)
         endif


         if (ijkp(l,3) .gt. nz) then
            write(*,*) 'k OB...',ijkp(l,:),m_arr(l)
         endif

         if (ijkp(l,3) .lt. 1) then
            write(*,*) 'k OB...',ijkp(l,:),m_arr(l)
         endif



      enddo

      return
      end SUBROUTINE check_index
c----------------------------------------------------------------------

      end module gutsp















































