c----------------------------------------------------------------------
      SUBROUTINE crossf(aa,bbmf,cc)
c The cross product is formed at the main cell center.  a is assumed
c be main cell contravarient (cell face) and b is assumed to be
c main cell covarient (cell edge).  The result is main cell
c contravarient (cell face).

c Can only center vectors on main cell for loops 3 to n...and can
c only extrapolate back for loops 2 to n-1.  Must handle other cases
c separately.

c The magnetic field does not require extrapolation to cell centers
c on boundaries since dB/dx = 0 is the boundary condition.  That is
c just copy the interior values to the boundary.
c----------------------------------------------------------------------
CVD$R VECTOR

      !include 'incurv.h'

      real aa(nx,ny,nz,3)        !main cell contravarient vector 
      real bbmf(nx,ny,nz,3)      !main cell contravarient vector
      real cc(nx,ny,nz,3)        !cross product result, main cell
                                 !contravarient (cell face)

      real ax,ay,az,bx,by,bz    !dummy vars
      real temp                 !used to vectorize loop
c      real zfrc(nz)             !0.5*dz_grid(k)/dz_cell(k)
c      real ct(nx,ny,nz,3)       !temp main cell center cross product
      real aac(3),bbc(3)


c extrapolate(/interpolate) to main cell center and do cross product


      call periodic(aa)
      call periodic(bbmf)


c      do 5 k=1,nz
c         zfrc(k) = 0.5*dz_grid(k)/dz_cell(k)
c 5       continue



      do 10 k=2,nz-1      
         do 10 j=2,ny-1
            do 10 i=2,nx-1

               im = i-1         !assume daa/dxyz = 0 at boundary
               jm = j-1         !bbmf is given on boundary
               km = k-1

               ax = 0.5*(aa(i,j,k,1) + aa(im,j,k,1))
               bx = 0.5*(bbmf(i,j,k,1) + bbmf(im,j,k,1))

               ay = 0.5*(aa(i,j,k,2) + aa(i,jm,k,2))
               by = 0.5*(bbmf(i,j,k,2) + bbmf(i,jm,k,2))

               az = zrat(k)*(aa(i,j,k,3) - aa(i,j,km,3)) + aa(i,j,km,3)
               bz = zrat(k)*(bbmf(i,j,k,3) - bbmf(i,j,km,3))
     x                     + bbmf(i,j,km,3)

               ct(i,j,k,1) = ay*bz - az*by
               ct(i,j,k,2) = az*bx - ax*bz
               ct(i,j,k,3) = ax*by - ay*bx

 10            continue

       call periodic(ct)

c extrapolate back to main cell contravarient positions.
c ...just average across cells since cell edges are centered
c about the grid points.
      
      do 60 k=2,nz-1
         do 60 j=2,ny-1
            do 60 i=2,nx-1

               ip = i+1
               jp = j+1
               kp = k+1

               cc(i,j,k,1) = 0.5*(ct(i,j,k,1) + ct(ip,j,k,1))
               cc(i,j,k,2) = 0.5*(ct(i,j,k,2) + ct(i,jp,k,2))
               cc(i,j,k,3) = 0.5*(ct(i,j,k,3) + ct(i,j,kp,3))

 60            continue

      call periodic(cc)


      return
      end SUBROUTINE crossf
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE cov_to_contra(bt,btmf)
c Converts total magnetic field from main cell covarient positions
c to main cell contravarient positions.  This is then used in the
c fluid velocity update routines.  This routine assumes that cell 
c edges and cell centers are "topologically centered".  So the grid
c points do not reside at the cell centers...rather they are offset
c a little so that the cell edges are equidistant from the k and k-1
c grid points.  In extrapolating the coventient vectors to the 
c contravarient vector positions, this assymetry is accounted for
c using a linear interpolation of the k and k-1 values to the grid
c point location.
c----------------------------------------------------------------------
CVD$R VECTOR
      !include 'incurv.h'

      real bt(nx,ny,nz,3),   !main cell covarient
     x     btmf(nx,ny,nz,3)  !main cell contravarient

      real bx1, bx2, by1, by2, bz1, bz2  !main cell center fields
c      real zrat           !ratio for doing linear interpolation
                          !to grid point position.
      real zplus, zminus  !position of main cell edges up and down
      real b_j, b_jm, b_i, b_im !intermediate step in average process

      do 10 k=2,nz-1
         do 10 j=2,ny-1
            do 10 i=2,nx-1

               ip = i+1
               jp = j+1
               kp = k+1
               im = i-1
               jm = j-1
               km = k-1

c The x component of B resides at the k and k-1 edges, so this
c requires the non-uniform grid interpolation

c               zplus = (qz(k+1) + qz(k))/2.0
c               zminus = (qz(k) + qz(k-1))/2.0
c               zrat = (qz(k) - zminus)/(zplus - zminus)
   
               b_j = bt(i,j,km,1) 
     x               + zrat(k)*(bt(i,j,k,1) - bt(i,j,km,1)) 
               b_jm = bt(i,jm,km,1)
     x                + zrat(k)*(bt(i,jm,k,1) - bt(i,jm,km,1))
               bx1 = (b_j + b_jm)/2.0

               b_j = bt(ip,j,km,1) 
     x               + zrat(k)*(bt(ip,j,k,1) - bt(ip,j,km,1)) 
               b_jm = bt(ip,jm,km,1)
     x                + zrat(k)*(bt(ip,jm,k,1) - bt(ip,jm,km,1))
               bx2 = (b_j + b_jm)/2.0

               
               b_i = bt(i,j,km,2) 
     x               + zrat(k)*(bt(i,j,k,2) - bt(i,j,km,2)) 
               b_im = bt(im,j,km,2)
     x                + zrat(k)*(bt(im,j,k,2) - bt(im,j,km,2))           
               by1 = (b_i + b_im)/2.0

               b_i = bt(i,jp,km,2) 
     x               + zrat(k)*(bt(i,jp,k,2) - bt(i,jp,km,2)) 
               b_im = bt(im,jp,km,2)
     x                + zrat(k)*(bt(im,jp,k,2) - bt(im,jp,km,2))
               by2 = (b_i + b_im)/2.0


               bz1 = 0.25*(bt(i,j,k,3) + bt(i,jm,k,3) +
     x                     bt(im,jm,k,3) + bt(im,j,k,3))
               bz2 = 0.25*(bt(i,j,kp,3) + bt(i,jm,kp,3) +
     x                     bt(im,jm,kp,3) + bt(im,j,kp,3))

               btmf(i,j,k,1) = 0.5*(bx1+bx2)
               btmf(i,j,k,2) = 0.5*(by1+by2)
               btmf(i,j,k,3) = 0.5*(bz1+bz2)

 10            continue

c      call boundaries(btmf)
      call periodic(btmf)

      return
      end SUBROUTINE cov_to_contra
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE sw_part_setup_maxwl_1(np,vp,vp1,xp,input_p,up,np_t_flg,
     x                               np_b_flg)
c----------------------------------------------------------------------
      !include 'incurv.h'

      real np(nx,ny,nz)
      real vp(Ni_max,3)
      real vp1(Ni_max,3)
      real xp(Ni_max,3)
      real input_p(3)
      real up(nx,ny,nz,3)
      real phi,theta,rnd,f,v
      real rand
      real vx,vy,vz
      real dvx,dvz,v1
      integer np_t_flg(Ni_max)
      integer np_b_flg(Ni_max)

      integer flg
      real nprat

      v1 = 1.0

      np_t_flg(:) = 0
      np_b_flg(:) = 0


      nprat = np_bottom/np_top

c      Ni_tot = 120000

c      do n = 0,procnum-1
c      if (my_rank .eq. n) then
      Ni_tot = Ni_tot + Ni_tot*nint((1/nprat)-1)/2

      do 10 l = 1,Ni_tot
c         write(*,*) 'procnum, random number...',n,pad_ranf()

c         phi = 2.0*pi*pad_ranf()
c         flg = 0
c         do 30 while (flg .eq. 0)
c            theta = pi*pad_ranf()
c            f = sin(theta)
c            rnd = pad_ranf()
c            if (f .ge. rnd) flg = 1
c 30      continue


         xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))
         xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))

         
         flg = 0
         do 20 while (flg .eq. 0)
            xp(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz)-qz(1))
             rnd = (1-nprat)* 
     x             ((1.0+tanh((xp(l,3)-qz(nz/2))/(Lo)))/2.0) +
     x             nprat
c            write(*,*) 'rnd...',rnd
             if (pad_ranf() .le. rnd) flg = 1
             if (xp(l,3) .ge. qz(nz/2)) np_t_flg(l) = 1
             if (xp(l,3) .lt. qz(nz/2)) np_b_flg(l) = 1

 20      continue

         ijkp(l,1) = nint(xp(l,1)/dx) !particle grid location index
         ijkp(l,2) = nint(xp(l,2)/dy)
         
         k=1
         do 50 while(xp(l,3) .gt. qz(k)) !find k on non-uniform 
            ijkp(l,3) = k       !grid
            k=k+1
 50      continue
         k=ijkp(l,3)
         if (xp(l,3) .gt. (qz(k)+(dz_grid(k)/2))) then
            ijkp(l,3) = k+1
         endif

         vth = 0.5*(vth_top + vth_bottom) + 
     x     0.5*(vth_top - vth_bottom)*tanh((qz(ijkp(l,3))-qz(nz/2))/Lo)

         flg = 0
         do 40 while (flg .eq. 0)

            vy = (400*pad_ranf())-200
            vz = (400*pad_ranf())-200
            vx = (400*pad_ranf())-200
            
            v = sqrt(vx**2 + vy**2 + vz**2)
            f = exp(-(v)**2 / vth**2)
            rnd = pad_ranf()
            if (f .gt. rnd) then 
               flg = 1
            endif
 40      continue

         ii = ijkp(l,1)
         kk = ijkp(l,3)
         dvx = -2.0*v1*cos(PI*(qx(ii)-qx(nx/2))/qx(nx/2))*
     x    tanh((qz(kk)-qz(nz/2))/Lo)*(cosh((qz(kk)-qz(nz/2))/Lo))**(-2)
         dvz = v1*sin(PI*(qx(ii)-qx(nx/2))/qx(nx/2))*
     x       (cosh((qz(kk)-qz(nz/2))/Lo))**(-2)
         
         vp(l,1) = vsw*(tanh((qz(k)-qz(nz/2))/(Lo))) + vx !+dvx
         vp(l,2) = vy 
         vp(l,3) = vz !+dvz 

         if (xp(l,3) .gt. qz(nz/2)) mix_ind(l) = 1
         if (xp(l,3) .le. qz(nz/2)) mix_ind(l) = 0
            
         do 45 m=1,3
            vp1(l,m) = vp(l,m)
            input_E = input_E + 
     x           0.5*mion*(vp(l,m)*km_to_m)**2 /beta
            input_p(m) = input_p(m) + mion*vp(l,m) / beta
 45      continue
 10      continue
c      endif
c      enddo


      write(*,*) 'get interp weights...'
      call get_interp_weights(xp)
      write(*,*) 'update_np...'
      call update_np(np)
      write(*,*) 'update_up...'
      call update_up(vp,np,up)
      write(*,*) 'update_up complete...'

   

      return
      end SUBROUTINE sw_part_setup_maxwl_1
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE part_setup_maxwl_h(np,vp,vp1,xp,input_p,up,np_t_flg,
     x                               np_b_flg)
c----------------------------------------------------------------------
      !include 'incurv.h'

      real np(nx,ny,nz)
      real vp(Ni_max,3)
      real vp1(Ni_max,3)
      real xp(Ni_max,3)
      real input_p(3)
      real up(nx,ny,nz,3)
      real phi,theta,rnd,f,v
      real rand
      real vx,vy,vz
      real dvx,dvz,v1
      integer np_t_flg(Ni_max)
      integer np_b_flg(Ni_max)
      integer Ni_tot_1

      integer flg
      real nprat

      v1 = 1.0

c      np_t_flg(:) = 0
c      np_b_flg(:) = 0


      nprat = np_bottom/np_top

      write(*,*) 'Ni_tot...',Ni_tot

      dNi = Ni_tot/2
      Ni_tot_1 = Ni_tot
      Ni_tot = Ni_tot + dNi/m_heavy

      write(*,*) 'Ni_tot...',Ni_tot,my_rank,m_heavy


      do 10 l = Ni_tot_1,Ni_tot

         xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))
         xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))
         
         flg = 0
         do 20 while (flg .eq. 0)
            xp(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz)-qz(1))
            rnd = ((1.0-tanh((xp(l,3)-qz(nz/2))/(Lo)))/2.0)
            if (pad_ranf() .le. rnd) flg = 1
c            if (xp(l,3) .ge. qz(nz/2)) then
c               np_t_flg(l) = 1
c               m_arr(l) = m_top
c               mrat(l) = mion/m_top
c            endif
c            if (xp(l,3) .lt. qz(nz/2)) then
            np_b_flg(l) = 1
            m_arr(l) = m_heavy*mion
            mrat(l) = 1.0/m_heavy
c            endif
            
 20      continue

         ijkp(l,1) = nint(xp(l,1)/dx) !particle grid location index
         ijkp(l,2) = nint(xp(l,2)/dy)
         
         k=1
         do 50 while(xp(l,3) .gt. qz(k)) !find k on non-uniform 
            ijkp(l,3) = k       !grid
            k=k+1
 50      continue
         k=ijkp(l,3)
         if (xp(l,3) .gt. (qz(k)+(dz_grid(k)/2))) then
            ijkp(l,3) = k+1
         endif

         vth = 0.5*(vth_top + vth_bottom) + 
     x     0.5*(vth_top - vth_bottom)*tanh((qz(ijkp(l,3))-qz(nz/2))/Lo)

         flg = 0
         do 40 while (flg .eq. 0)
            v = (2*vth_max*pad_ranf())-vth_max
            f = exp(-(v)**2 / vth**2)
            rnd = pad_ranf()
            if (f .ge. rnd) then 
               flg = 1
               vx = v
            endif
 40      continue
         flg = 0
         do 42 while (flg .eq. 0)
            v = (2*vth_max*pad_ranf())-vth_max
            f = exp(-(v)**2 / vth**2)
            rnd = pad_ranf()
            if (f .ge. rnd) then 
               flg = 1
               vy = v
            endif
 42      continue
         flg = 0
         do 44 while (flg .eq. 0)
            v = (2*vth_max*pad_ranf())-vth_max
            f = exp(-(v)**2 / vth**2)
            rnd = pad_ranf()
            if (f .ge. rnd) then 
               flg = 1
               vz = v
            endif
 44      continue

         ii = ijkp(l,1)
         kk = ijkp(l,3)
         dvx = -2.0*v1*cos(PI*(qx(ii)-qx(nx/2))/qx(nx/2))*
     x    tanh((qz(kk)-qz(nz/2))/Lo)*(cosh((qz(kk)-qz(nz/2))/Lo))**(-2)
         dvz = v1*sin(PI*(qx(ii)-qx(nx/2))/qx(nx/2))*
     x       (cosh((qz(kk)-qz(nz/2))/Lo))**(-2)
         
         vp(l,1) = vsw*(tanh((qz(k)-qz(nz/2))/(Lo))) + vx !+dvx
         vp(l,2) = vy 
         vp(l,3) = vz !+dvz 

         if (xp(l,3) .gt. qz(nz/2)) mix_ind(l) = 1
         if (xp(l,3) .le. qz(nz/2)) mix_ind(l) = 0
            
         do 45 m=1,3
            vp1(l,m) = vp(l,m)
            input_E = input_E + 
     x           0.5*m_arr(l)*(vp(l,m)*km_to_m)**2 /beta
            input_p(m) = input_p(m) + m_arr(l)*vp(l,m) / beta
 45      continue
 10      continue
      
         write(*,*) 'Ni_tot...',Ni_tot

         Ni_tot_1 = Ni_tot+1
c         Ni_tot = Ni_tot + Ni_tot*np_top/np_bottom
         Ni_tot = Ni_tot + dNi/m_heavy

         write(*,*) 'Ni_tot...',Ni_tot,my_rank


         do 60 l = Ni_tot_1,Ni_tot

         xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))
         xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))

         
         flg = 0
         do 70 while (flg .eq. 0)
            xp(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz)-qz(1))
             rnd = (1.0+tanh((xp(l,3)-qz(nz/2))/(Lo)))/2.0
             if (pad_ranf() .le. rnd) flg = 1
c             if (xp(l,3) .ge. qz(nz/2)) then
             np_t_flg(l) = 1
             m_arr(l) = m_heavy*mion
             mrat(l) = 1.0/m_heavy
c             endif
c             if (xp(l,3) .lt. qz(nz/2)) then
c                np_b_flg(l) = 1
c                m_arr(l) = m_bottom
c                mrat(l) = mion/m_bottom
c             endif

 70      continue

         ijkp(l,1) = nint(xp(l,1)/dx) !particle grid location index
         ijkp(l,2) = nint(xp(l,2)/dy)
         
         k=1
         do 100 while(xp(l,3) .gt. qz(k)) !find k on non-uniform 
            ijkp(l,3) = k       !grid
            k=k+1
 100     continue
         k=ijkp(l,3)
         if (xp(l,3) .gt. (qz(k)+(dz_grid(k)/2))) then
            ijkp(l,3) = k+1
         endif

         vth = 0.5*(vth_top + vth_bottom) + 
     x     0.5*(vth_top - vth_bottom)*tanh((qz(ijkp(l,3))-qz(nz/2))/Lo)

         flg = 0
         do 90 while (flg .eq. 0)
            v = (2*vth_max*pad_ranf())-vth_max
            f = exp(-(v)**2 / vth**2)
            rnd = pad_ranf()
            if (f .ge. rnd) then 
               flg = 1
               vx = v
            endif
 90      continue
         flg = 0
         do 92 while (flg .eq. 0)
            v = (2*vth_max*pad_ranf())-vth_max
            f = exp(-(v)**2 / vth**2)
            rnd = pad_ranf()
            if (f .ge. rnd) then 
               flg = 1
               vy = v
            endif
 92      continue
         flg = 0
         do 94 while (flg .eq. 0)
            v = (2*vth_max*pad_ranf())-vth_max
            f = exp(-(v)**2 / vth**2)
            rnd = pad_ranf()
            if (f .ge. rnd) then 
               flg = 1
               vz = v
            endif
 94      continue

         ii = ijkp(l,1)
         kk = ijkp(l,3)
         dvx = -2.0*v1*cos(PI*(qx(ii)-qx(nx/2))/qx(nx/2))*
     x    tanh((qz(kk)-qz(nz/2))/Lo)*(cosh((qz(kk)-qz(nz/2))/Lo))**(-2)
         dvz = v1*sin(PI*(qx(ii)-qx(nx/2))/qx(nx/2))*
     x       (cosh((qz(kk)-qz(nz/2))/Lo))**(-2)
         
         vp(l,1) = vsw*(tanh((qz(k)-qz(nz/2))/(Lo))) + vx !+dvx
         vp(l,2) = vy 
         vp(l,3) = vz !+dvz 

         if (xp(l,3) .gt. qz(nz/2)) mix_ind(l) = 1
         if (xp(l,3) .le. qz(nz/2)) mix_ind(l) = 0
            
         do 95 m=1,3
            vp1(l,m) = vp(l,m)
            input_E = input_E + 
     x           0.5*m_arr(l)*(vp(l,m)*km_to_m)**2 /beta
            input_p(m) = input_p(m) + m_arr(l)*vp(l,m) / beta
 95      continue
 60   continue
c      endif
c      enddo



 
      write(*,*) 'get interp weights...'
      call get_interp_weights(xp)
      write(*,*) 'update_np...'
      call update_np(np)
      write(*,*) 'update_up...'
      call update_up(vp,np,up)
      write(*,*) 'update_up complete...'

   

      return
      end SUBROUTINE part_setup_maxwl_h
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE part_setup_maxwl_p(np,vp,vp1,xp,input_p,up,np_t_flg,
     x                               np_b_flg)
c----------------------------------------------------------------------
      !include 'incurv.h'

      real np(nx,ny,nz)
      real vp(Ni_max,3)
      real vp1(Ni_max,3)
      real xp(Ni_max,3)
      real input_p(3)
      real up(nx,ny,nz,3)
      real phi,theta,rnd,f,v
      real rand
      real vx,vy,vz
      real dvx,dvz,v1
      integer np_t_flg(Ni_max)
      integer np_b_flg(Ni_max)
      integer Ni_tot_1

      integer flg
      real nprat

      v1 = 1.0

      np_t_flg(:) = 0
      np_b_flg(:) = 0


      nprat = np_bottom/np_top


      do 10 l = 1,Ni_tot

         xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))
         xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))
         
         flg = 0
         do 20 while (flg .eq. 0)
            xp(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz)-qz(1))
            rnd = ((1.0-tanh((xp(l,3)-qz(nz/2))/(Lo)))/2.0)
            if (pad_ranf() .le. rnd) flg = 1
c            if (xp(l,3) .ge. qz(nz/2)) then
c               np_t_flg(l) = 1
c               m_arr(l) = m_top
c               mrat(l) = mion/m_top
c            endif
c            if (xp(l,3) .lt. qz(nz/2)) then
            np_b_flg(l) = 1
            m_arr(l) = m_bottom
            mrat(l) = mion/m_bottom
c            endif
            
 20      continue

         ijkp(l,1) = nint(xp(l,1)/dx) !particle grid location index
         ijkp(l,2) = nint(xp(l,2)/dy)
         
         k=1
         do 50 while(xp(l,3) .gt. qz(k)) !find k on non-uniform 
            ijkp(l,3) = k       !grid
            k=k+1
 50      continue
         k=ijkp(l,3)
         if (xp(l,3) .gt. (qz(k)+(dz_grid(k)/2))) then
            ijkp(l,3) = k+1
         endif

         vth = 0.5*(vth_top + vth_bottom) + 
     x     0.5*(vth_top - vth_bottom)*tanh((qz(ijkp(l,3))-qz(nz/2))/Lo)

         flg = 0
         do 40 while (flg .eq. 0)
            v = (2*vth_max*pad_ranf())-vth_max
            f = exp(-(v)**2 / vth**2)
            rnd = pad_ranf()
            if (f .ge. rnd) then 
               flg = 1
               vx = v
            endif
 40      continue
         flg = 0
         do 42 while (flg .eq. 0)
            v = (2*vth_max*pad_ranf())-vth_max
            f = exp(-(v)**2 / vth**2)
            rnd = pad_ranf()
            if (f .ge. rnd) then 
               flg = 1
               vy = v
            endif
 42      continue
         flg = 0
         do 44 while (flg .eq. 0)
            v = (2*vth_max*pad_ranf())-vth_max
            f = exp(-(v)**2 / vth**2)
            rnd = pad_ranf()
            if (f .ge. rnd) then 
               flg = 1
               vz = v
            endif
 44      continue

         ii = ijkp(l,1)
         kk = ijkp(l,3)
         dvx = -2.0*v1*cos(PI*(qx(ii)-qx(nx/2))/qx(nx/2))*
     x    tanh((qz(kk)-qz(nz/2))/Lo)*(cosh((qz(kk)-qz(nz/2))/Lo))**(-2)
         dvz = v1*sin(PI*(qx(ii)-qx(nx/2))/qx(nx/2))*
     x       (cosh((qz(kk)-qz(nz/2))/Lo))**(-2)
         
         vp(l,1) = vsw*(tanh((qz(k)-qz(nz/2))/(Lo))) + vx !+dvx
         vp(l,2) = vy 
         vp(l,3) = vz !+dvz 

         if (xp(l,3) .gt. qz(nz/2)) mix_ind(l) = 1
         if (xp(l,3) .le. qz(nz/2)) mix_ind(l) = 0
            
         do 45 m=1,3
            vp1(l,m) = vp(l,m)
            input_E = input_E + 
     x           0.5*m_arr(l)*(vp(l,m)*km_to_m)**2 /beta
            input_p(m) = input_p(m) + m_arr(l)*vp(l,m) / beta
 45      continue
 10      continue
      

         Ni_tot_1 = Ni_tot+1
         Ni_tot = Ni_tot + Ni_tot*np_top/np_bottom

         do 60 l = Ni_tot_1,Ni_tot

         xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))
         xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))

         
         flg = 0
         do 70 while (flg .eq. 0)
            xp(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz)-qz(1))
             rnd = (1.0+tanh((xp(l,3)-qz(nz/2))/(Lo)))/2.0
             if (pad_ranf() .le. rnd) flg = 1
c             if (xp(l,3) .ge. qz(nz/2)) then
             np_t_flg(l) = 1
             m_arr(l) = m_top
             mrat(l) = mion/m_top
c             endif
c             if (xp(l,3) .lt. qz(nz/2)) then
c                np_b_flg(l) = 1
c                m_arr(l) = m_bottom
c                mrat(l) = mion/m_bottom
c             endif

 70      continue

         ijkp(l,1) = nint(xp(l,1)/dx) !particle grid location index
         ijkp(l,2) = nint(xp(l,2)/dy)
         
         k=1
         do 100 while(xp(l,3) .gt. qz(k)) !find k on non-uniform 
            ijkp(l,3) = k       !grid
            k=k+1
 100     continue
         k=ijkp(l,3)
         if (xp(l,3) .gt. (qz(k)+(dz_grid(k)/2))) then
            ijkp(l,3) = k+1
         endif

         vth = 0.5*(vth_top + vth_bottom) + 
     x     0.5*(vth_top - vth_bottom)*tanh((qz(ijkp(l,3))-qz(nz/2))/Lo)

         flg = 0
         do 90 while (flg .eq. 0)
            v = (2*vth_max*pad_ranf())-vth_max
            f = exp(-(v)**2 / vth**2)
            rnd = pad_ranf()
            if (f .ge. rnd) then 
               flg = 1
               vx = v
            endif
 90      continue
         flg = 0
         do 92 while (flg .eq. 0)
            v = (2*vth_max*pad_ranf())-vth_max
            f = exp(-(v)**2 / vth**2)
            rnd = pad_ranf()
            if (f .ge. rnd) then 
               flg = 1
               vy = v
            endif
 92      continue
         flg = 0
         do 94 while (flg .eq. 0)
            v = (2*vth_max*pad_ranf())-vth_max
            f = exp(-(v)**2 / vth**2)
            rnd = pad_ranf()
            if (f .ge. rnd) then 
               flg = 1
               vz = v
            endif
 94      continue

         ii = ijkp(l,1)
         kk = ijkp(l,3)
         dvx = -2.0*v1*cos(PI*(qx(ii)-qx(nx/2))/qx(nx/2))*
     x    tanh((qz(kk)-qz(nz/2))/Lo)*(cosh((qz(kk)-qz(nz/2))/Lo))**(-2)
         dvz = v1*sin(PI*(qx(ii)-qx(nx/2))/qx(nx/2))*
     x       (cosh((qz(kk)-qz(nz/2))/Lo))**(-2)
         
         vp(l,1) = vsw*(tanh((qz(k)-qz(nz/2))/(Lo))) + vx !+dvx
         vp(l,2) = vy 
         vp(l,3) = vz !+dvz 

         if (xp(l,3) .gt. qz(nz/2)) mix_ind(l) = 1
         if (xp(l,3) .le. qz(nz/2)) mix_ind(l) = 0
            
         do 95 m=1,3
            vp1(l,m) = vp(l,m)
            input_E = input_E + 
     x           0.5*m_arr(l)*(vp(l,m)*km_to_m)**2 /beta
            input_p(m) = input_p(m) + m_arr(l)*vp(l,m) / beta
 95      continue
 60   continue
c      endif
c      enddo



 
      write(*,*) 'get interp weights...'
      call get_interp_weights(xp)
      write(*,*) 'update_np...'
      call update_np(np)
      write(*,*) 'update_up...'
      call update_up(vp,np,up)
      write(*,*) 'update_up complete...'

   

      return
      end SUBROUTINE part_setup_maxwl_p
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE sw_part_setup_maxwl_eq(np,vp,vp1,xp,input_p,up,
     x                               np_t_flg,np_b_flg)
c----------------------------------------------------------------------
      !include 'incurv.h'

      real np(nx,ny,nz)
      real vp(Ni_max,3)
      real vp1(Ni_max,3)
      real xp(Ni_max,3)
      real input_p(3)
      real up(nx,ny,nz,3)
      real phi,theta,rnd,f,v
      real rand
      real vx,vy,vz
      real dvx,dvz,v1
      integer np_t_flg(Ni_max)
      integer np_b_flg(Ni_max)
      integer Ni_tot_1
      real n1,n2,mr,nr,p0

      integer flg
      real nprat

      v1 = 1.0

      np_t_flg(:) = 0
      np_b_flg(:) = 0


      nprat = np_bottom/np_top

      do 10 l = 1,Ni_tot

         xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))
         xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))
         
         flg = 0
         do 20 while (flg .eq. 0)
            xp(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz)-qz(1))
            rnd = ((1.0-tanh((xp(l,3)-qz(nz/2))/(Lo)))/2.0)
            if (pad_ranf() .le. rnd) flg = 1
c            if (xp(l,3) .ge. qz(nz/2)) then
c               np_t_flg(l) = 1
c               m_arr(l) = m_top
c               mrat(l) = mion/m_top
c            endif
c            if (xp(l,3) .lt. qz(nz/2)) then
            np_b_flg(l) = 1
            m_arr(l) = m_bottom
            mrat(l) = mion/m_bottom
c            endif
            
 20      continue

         ijkp(l,1) = nint(xp(l,1)/dx) !particle grid location index
         ijkp(l,2) = nint(xp(l,2)/dy)
         
         k=1
         do 50 while(xp(l,3) .gt. qz(k)) !find k on non-uniform 
            ijkp(l,3) = k       !grid
            k=k+1
 50      continue
         k=ijkp(l,3)
         if (xp(l,3) .gt. (qz(k)+(dz_grid(k)/2))) then
            ijkp(l,3) = k+1
         endif

         vth = vth_bottom

c         vth = 0.5*(vth_top + vth_bottom) + 
c     x     0.5*(vth_top - vth_bottom)*tanh((qz(ijkp(l,3))-qz(nz/2))/Lo)

         flg = 0
         do 40 while (flg .eq. 0)
            v = (2*vth_max*pad_ranf())-vth_max
            f = exp(-(v)**2 / vth**2)
            rnd = pad_ranf()
            if (f .ge. rnd) then 
               flg = 1
               vx = v
            endif
 40      continue
         flg = 0
         do 42 while (flg .eq. 0)
            v = (2*vth_max*pad_ranf())-vth_max
            f = exp(-(v)**2 / vth**2)
            rnd = pad_ranf()
            if (f .ge. rnd) then 
               flg = 1
               vy = v
            endif
 42      continue
         flg = 0
         do 44 while (flg .eq. 0)
            v = (2*vth_max*pad_ranf())-vth_max
            f = exp(-(v)**2 / vth**2)
            rnd = pad_ranf()
            if (f .ge. rnd) then 
               flg = 1
               vz = v
            endif
 44      continue

         ii = ijkp(l,1)
         kk = ijkp(l,3)
         dvx = -2.0*v1*cos(PI*(qx(ii)-qx(nx/2))/qx(nx/2))*
     x    tanh((qz(kk)-qz(nz/2))/Lo)*(cosh((qz(kk)-qz(nz/2))/Lo))**(-2)
         dvz = v1*sin(PI*(qx(ii)-qx(nx/2))/qx(nx/2))*
     x       (cosh((qz(kk)-qz(nz/2))/Lo))**(-2)
         
         vp(l,1) = vsw*(tanh((qz(k)-qz(nz/2))/(Lo))) + vx !+dvx
         vp(l,2) = vy 
         vp(l,3) = vz !+dvz 

         if (xp(l,3) .gt. qz(nz/2)) mix_ind(l) = 1
         if (xp(l,3) .le. qz(nz/2)) mix_ind(l) = 0
            
         do 45 m=1,3
            vp1(l,m) = vp(l,m)
            input_E = input_E + 
     x           0.5*m_arr(l)*(vp(l,m)*km_to_m)**2 /beta
            input_p(m) = input_p(m) + m_arr(l)*vp(l,m) / beta
 45      continue
 10      continue
      

         Ni_tot_1 = Ni_tot+1
         Ni_tot = Ni_tot + Ni_tot*np_top/np_bottom



         do 60 l = Ni_tot_1,Ni_tot

         xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))
         xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))

         
         flg = 0
         do 70 while (flg .eq. 0)
            xp(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz)-qz(1))
             rnd = (1.0+tanh((xp(l,3)-qz(nz/2))/(Lo)))/2.0
             if (pad_ranf() .le. rnd) flg = 1
c             if (xp(l,3) .ge. qz(nz/2)) then
             np_t_flg(l) = 1
             m_arr(l) = m_top
             mrat(l) = mion/m_top
c             endif
c             if (xp(l,3) .lt. qz(nz/2)) then
c                np_b_flg(l) = 1
c                m_arr(l) = m_bottom
c                mrat(l) = mion/m_bottom
c             endif

 70      continue

         ijkp(l,1) = nint(xp(l,1)/dx) !particle grid location index
         ijkp(l,2) = nint(xp(l,2)/dy)
         
         k=1
         do 100 while(xp(l,3) .gt. qz(k)) !find k on non-uniform 
            ijkp(l,3) = k       !grid
            k=k+1
 100     continue
         k=ijkp(l,3)
         if (xp(l,3) .gt. (qz(k)+(dz_grid(k)/2))) then
            ijkp(l,3) = k+1
         endif


c update np and then determine eq vth for protons (vth_top)


c      write(*,*) 'get interp weights...'
c      call get_interp_weights(xp)
c      write(*,*) 'update_np...'
c      call update_np(np)
c      write(*,*) 'update_up...'


c 60   continue

c      do 62 l = Ni_tot_1,Ni_tot


c         ii = ijkp(l,1)
c         jj = ijkp(l,2)
         kk = ijkp(l,3)

         vth = vth_top

         n1 = np_bottom*(1.0-tanh((xp(l,3)-qz(nz/2))/(Lo)))/2.0
         n2 = np_top*(1.0+tanh((xp(l,3)-qz(nz/2))/(Lo)))/2.0
      
         p0 = m_bottom*np_bottom*vth_bottom**2

         if ((n1 .gt. 0) .and. (n2 .gt. 0)) then 
         vth = sqrt((p0 - n1*m_bottom*vth_bottom*vth_bottom)/(n2*m_top))
c            write(*,*) 'vth_top...',l,vth,xp(l,3)
         endif

c         if ((n2 .gt. 0) .and. (n1 .gt. 0)) then 
c            vth = (n1*m_bottom*vth_bottom**2/(n2*m_top))**(0.5)
c            write(*,*) 'vth_top...',l,vth,xp(l,3)
c         endif



c         vth = 0.5*(vth_top + vth_bottom) + 
c     x     0.5*(vth_top - vth_bottom)*tanh((qz(ijkp(l,3))-qz(nz/2))/Lo)

         flg = 0
         do 90 while (flg .eq. 0)
            v = (2*vth_max*pad_ranf())-vth_max
            f = exp(-(v)**2 / vth**2)
            rnd = pad_ranf()
            if (f .ge. rnd) then 
               flg = 1
               vx = v
            endif
 90      continue
         flg = 0
         do 92 while (flg .eq. 0)
            v = (2*vth_max*pad_ranf())-vth_max
            f = exp(-(v)**2 / vth**2)
            rnd = pad_ranf()
            if (f .ge. rnd) then 
               flg = 1
               vy = v
            endif
 92      continue
         flg = 0
         do 94 while (flg .eq. 0)
            v = (2*vth_max*pad_ranf())-vth_max
            f = exp(-(v)**2 / vth**2)
            rnd = pad_ranf()
            if (f .ge. rnd) then 
               flg = 1
               vz = v
            endif
 94      continue

         ii = ijkp(l,1)
         kk = ijkp(l,3)
         dvx = -2.0*v1*cos(PI*(qx(ii)-qx(nx/2))/qx(nx/2))*
     x    tanh((qz(kk)-qz(nz/2))/Lo)*(cosh((qz(kk)-qz(nz/2))/Lo))**(-2)
         dvz = v1*sin(PI*(qx(ii)-qx(nx/2))/qx(nx/2))*
     x       (cosh((qz(kk)-qz(nz/2))/Lo))**(-2)
         
         vp(l,1) = vsw*(tanh((qz(k)-qz(nz/2))/(Lo))) + vx !+dvx
         vp(l,2) = vy 
         vp(l,3) = vz !+dvz 

         if (xp(l,3) .gt. qz(nz/2)) mix_ind(l) = 1
         if (xp(l,3) .le. qz(nz/2)) mix_ind(l) = 0
            
         do 95 m=1,3
            vp1(l,m) = vp(l,m)
            input_E = input_E + 
     x           0.5*m_arr(l)*(vp(l,m)*km_to_m)**2 /beta
            input_p(m) = input_p(m) + m_arr(l)*vp(l,m) / beta
 95      continue
 60   continue
c      endif
c      enddo



 
      write(*,*) 'get interp weights...'
      call get_interp_weights(xp)
      write(*,*) 'update_np...'
      call update_np(np)
      write(*,*) 'update_up...'
      call update_up(vp,np,up)
      write(*,*) 'update_up complete...'

   

      return
      end SUBROUTINE sw_part_setup_maxwl_eq
c----------------------------------------------------------------------



c----------------------------------------------------------------------
      SUBROUTINE sw_part_setup_temp(np,vp,vp1,xp,input_p,up)
c----------------------------------------------------------------------
      !include 'incurv.h'

      real np(nx,ny,nz)
      real vp(Ni_max,3)
      real vp1(Ni_max,3)
      real xp(Ni_max,3)
      real input_p(3)
      real up(nx,ny,nz,3)
      real phi,theta,rnd,f,v
      real rand

c      Ni_tot = 120000

c      do n = 0,procnum-1
c      if (my_rank .eq. n) then
      do 10 l = 1,Ni_tot
c         write(*,*) 'procnum, random number...',n,pad_ranf()

         phi = 2.0*pi*pad_ranf()
         flg = 0
         do 30 while (flg .eq. 0)
            theta = pi*pad_ranf()
            f = sin(theta)
            rnd = pad_ranf()
            if (f .ge. rnd) flg = 1
 30      continue


         xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))
         xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))
         xp(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz)-qz(1))


         flg = 0
         do 40 while (flg .eq. 0)
            v = (100*pad_ranf())
c            f = (vth**2/exp(1.0))*v**2*exp(-(v)**2 / vth**2)
            f = exp(-(v)**2 / vth**2)
            rnd = pad_ranf()
            if (f .ge. rnd) then 
               flg = 1
               vp(l,1) = vsw + v*cos(phi)*sin(theta)
               vp(l,2) = v*sin(phi)*sin(theta)
               vp(l,3) = v*cos(theta)
            endif

c         vp(l,1) = -vsw
c         vp(l,2) = 0.0
c         vp(l,3) = 0.0

 40      continue


         ijkp(l,1) = nint(xp(l,1)/dx) !particle grid location index
         ijkp(l,2) = nint(xp(l,2)/dy)
         
         k=1
         do 50 while(xp(l,3) .gt. qz(k)) !find k on non-uniform 
            ijkp(l,3) = k       !grid
            k=k+1
 50      continue
         k=ijkp(l,3)
         if (xp(l,3) .gt. (qz(k)+(dz_grid(k)/2))) then
            ijkp(l,3) = k+1
         endif
         do 45 m=1,3
            vp1(l,m) = vp(l,m)
            input_E = input_E + 
     x           0.5*mion*(vp(l,m)*km_to_m)**2 /beta
            input_p(m) = input_p(m) + mion*vp(l,m) / beta
 45      continue
 10      continue
c      endif
c      enddo

 
      write(*,*) 'get interp weights...'
      call get_interp_weights(xp)
      write(*,*) 'update_np...'
      call update_np(np)
      write(*,*) 'update_up...'
      call update_up(vp,np,up)
      write(*,*) 'update_up complete...'

   

      return
      end SUBROUTINE sw_part_setup_temp
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE sw_part_setup_maxwl(np,vp,vp1,xp,input_p,up,np_t_flg,
     x                               np_b_flg)
c----------------------------------------------------------------------
      !include 'incurv.h'

      real np(nx,ny,nz)
      real vp(Ni_max,3)
      real vp1(Ni_max,3)
      real xp(Ni_max,3)
      real input_p(3)
      real up(nx,ny,nz,3)
      real phi,theta,rnd,f,v
      real rand
      real vx,vy,vz
      real dvx,dvz,v1
      integer np_t_flg(Ni_max)
      integer np_b_flg(Ni_max)
      integer Ni_tot_1

      integer flg
      real nprat

      v1 = 1.0

      np_t_flg(:) = 0
      np_b_flg(:) = 0


      nprat = np_bottom/np_top
      Ni_tot_heavy = Ni_tot

c add heavies to bottom first

      do 10 l = 1,Ni_tot

         xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))
         xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))
         
         flg = 0
         do 20 while (flg .eq. 0)
            xp(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz)-qz(1))
            rnd = ((1.0-tanh((xp(l,3)-qz(nz/2))/(Lo)))/2.0)
            if (pad_ranf() .le. rnd) flg = 1
c            if (xp(l,3) .ge. qz(nz/2)) then
c               np_t_flg(l) = 1
c               m_arr(l) = m_top
c               mrat(l) = mion/m_top
c            endif
c            if (xp(l,3) .lt. qz(nz/2)) then
            np_b_flg(l) = 1
            m_arr(l) = m_bottom
            mrat(l) = mion/m_bottom
c            endif
            
 20      continue

         ijkp(l,1) = nint(xp(l,1)/dx) !particle grid location index
         ijkp(l,2) = nint(xp(l,2)/dy)
         
         k=1
         do 30 while(xp(l,3) .gt. qz(k)) !find k on non-uniform 
            ijkp(l,3) = k       !grid
            k=k+1
 30      continue
         k=ijkp(l,3)
         if (xp(l,3) .gt. (qz(k)+(dz_grid(k)/2))) then
            ijkp(l,3) = k+1
         endif

         vth = vth_bottom

c         vth = 0.5*(vth_top + vth_bottom) + 
c     x     0.5*(vth_top - vth_bottom)*tanh((qz(ijkp(l,3))-qz(nz/2))/Lo)

         flg = 0
         do 40 while (flg .eq. 0)
            v = (2*vth_max*pad_ranf())-vth_max
            f = exp(-(v)**2 / vth**2)
            rnd = pad_ranf()
            if (f .ge. rnd) then 
               flg = 1
               vx = v
            endif
 40      continue
         flg = 0
         do 42 while (flg .eq. 0)
            v = (2*vth_max*pad_ranf())-vth_max
            f = exp(-(v)**2 / vth**2)
            rnd = pad_ranf()
            if (f .ge. rnd) then 
               flg = 1
               vy = v
            endif
 42      continue
         flg = 0
         do 44 while (flg .eq. 0)
            v = (2*vth_max*pad_ranf())-vth_max
            f = exp(-(v)**2 / vth**2)
            rnd = pad_ranf()
            if (f .ge. rnd) then 
               flg = 1
               vz = v
            endif
 44      continue

         ii = ijkp(l,1)
         kk = ijkp(l,3)
         dvx = -2.0*v1*cos(PI*(qx(ii)-qx(nx/2))/qx(nx/2))*
     x    tanh((qz(kk)-qz(nz/2))/Lo)*(cosh((qz(kk)-qz(nz/2))/Lo))**(-2)
         dvz = v1*sin(PI*(qx(ii)-qx(nx/2))/qx(nx/2))*
     x       (cosh((qz(kk)-qz(nz/2))/Lo))**(-2)
         
         vp(l,1) = vsw*(tanh((qz(k)-qz(nz/2))/(Lo))) + vx !+dvx
         vp(l,2) = vy 
         vp(l,3) = vz !+dvz 

         if (xp(l,3) .gt. qz(nz/2)) mix_ind(l) = 1
         if (xp(l,3) .le. qz(nz/2)) mix_ind(l) = 0
            
         do 45 m=1,3
            vp1(l,m) = vp(l,m)
            input_E = input_E + 
     x           0.5*m_arr(l)*(vp(l,m)*km_to_m)**2 /beta
            input_p(m) = input_p(m) + m_arr(l)*vp(l,m) / beta
 45      continue
 10      continue


c add 10% initial Ni_tot for uniform proton background

         Ni_tot_1 = Ni_tot + 1
         Ni_tot = Ni_tot + nint(Ni_tot*0.1)


         do 49 l = Ni_tot_1,Ni_tot
            xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))
            xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))
            xp(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz)-qz(1))         
            np_b_flg(l) = 0
            m_arr(l) = mion
            mrat(l) = 1.0

            ijkp(l,1) = nint(xp(l,1)/dx) !particle grid location index
            ijkp(l,2) = nint(xp(l,2)/dy)
            
            k=1
            do 50 while(xp(l,3) .gt. qz(k)) !find k on non-uniform 
               ijkp(l,3) = k    !grid
               k=k+1
 50         continue
            k=ijkp(l,3)
            if (xp(l,3) .gt. (qz(k)+(dz_grid(k)/2))) then
               ijkp(l,3) = k+1
            endif
            
            vth = vth_top
            
            flg = 0
            do 52 while (flg .eq. 0)
               v = (2*vth_max*pad_ranf())-vth_max
               f = exp(-(v)**2 / vth**2)
               rnd = pad_ranf()
               if (f .ge. rnd) then 
                  flg = 1
                  vx = v
               endif
 52         continue
            flg = 0
            do 54 while (flg .eq. 0)
               v = (2*vth_max*pad_ranf())-vth_max
               f = exp(-(v)**2 / vth**2)
               rnd = pad_ranf()
               if (f .ge. rnd) then 
                  flg = 1
                  vy = v
               endif
 54         continue
            flg = 0
            do 56 while (flg .eq. 0)
               v = (2*vth_max*pad_ranf())-vth_max
               f = exp(-(v)**2 / vth**2)
               rnd = pad_ranf()
               if (f .ge. rnd) then 
                  flg = 1
                  vz = v
               endif
 56         continue
            
            ii = ijkp(l,1)
            kk = ijkp(l,3)
            dvx = -2.0*v1*cos(PI*(qx(ii)-qx(nx/2))/qx(nx/2))*
     x           tanh((qz(kk)-qz(nz/2))/Lo)*(cosh((qz(kk)-qz(nz/2))/Lo))**(-2)
            dvz = v1*sin(PI*(qx(ii)-qx(nx/2))/qx(nx/2))*
     x           (cosh((qz(kk)-qz(nz/2))/Lo))**(-2)
            
            vp(l,1) = vsw*(tanh((qz(k)-qz(nz/2))/(Lo))) + vx !+dvx
            vp(l,2) = vy 
            vp(l,3) = vz        !+dvz 
            
            if (xp(l,3) .gt. qz(nz/2)) mix_ind(l) = 1
            if (xp(l,3) .le. qz(nz/2)) mix_ind(l) = 0
            
            do  m=1,3
               vp1(l,m) = vp(l,m)
               input_E = input_E + 
     x              0.5*m_arr(l)*(vp(l,m)*km_to_m)**2 /beta
               input_p(m) = input_p(m) + m_arr(l)*vp(l,m) / beta
            enddo
            

 49      enddo


c add protons to top half

         Ni_tot_1 = Ni_tot+1
         Ni_tot = Ni_tot + Ni_tot_heavy*np_top/np_bottom

         do 60 l = Ni_tot_1,Ni_tot

         xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))
         xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))

         
         flg = 0
         do 70 while (flg .eq. 0)
            xp(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz)-qz(1))
             rnd = (1.0+tanh((xp(l,3)-qz(nz/2))/(Lo)))/2.0
             if (pad_ranf() .le. rnd) flg = 1
c             if (xp(l,3) .ge. qz(nz/2)) then
             np_t_flg(l) = 1
             m_arr(l) = m_top
             mrat(l) = mion/m_top
c             endif
c             if (xp(l,3) .lt. qz(nz/2)) then
c                np_b_flg(l) = 1
c                m_arr(l) = m_bottom
c                mrat(l) = mion/m_bottom
c             endif

 70      continue

         ijkp(l,1) = nint(xp(l,1)/dx) !particle grid location index
         ijkp(l,2) = nint(xp(l,2)/dy)
         
         k=1
         do 100 while(xp(l,3) .gt. qz(k)) !find k on non-uniform 
            ijkp(l,3) = k       !grid
            k=k+1
 100     continue
         k=ijkp(l,3)
         if (xp(l,3) .gt. (qz(k)+(dz_grid(k)/2))) then
            ijkp(l,3) = k+1
         endif

         vth = vth_top

c         vth = 0.5*(vth_top + vth_bottom) + 
c     x     0.5*(vth_top - vth_bottom)*tanh((qz(ijkp(l,3))-qz(nz/2))/Lo)

         flg = 0
         do 90 while (flg .eq. 0)
            v = (2*vth_max*pad_ranf())-vth_max
            f = exp(-(v)**2 / vth**2)
            rnd = pad_ranf()
            if (f .ge. rnd) then 
               flg = 1
               vx = v
            endif
 90      continue
         flg = 0
         do 92 while (flg .eq. 0)
            v = (2*vth_max*pad_ranf())-vth_max
            f = exp(-(v)**2 / vth**2)
            rnd = pad_ranf()
            if (f .ge. rnd) then 
               flg = 1
               vy = v
            endif
 92      continue
         flg = 0
         do 94 while (flg .eq. 0)
            v = (2*vth_max*pad_ranf())-vth_max
            f = exp(-(v)**2 / vth**2)
            rnd = pad_ranf()
            if (f .ge. rnd) then 
               flg = 1
               vz = v
            endif
 94      continue

         ii = ijkp(l,1)
         kk = ijkp(l,3)
         dvx = -2.0*v1*cos(PI*(qx(ii)-qx(nx/2))/qx(nx/2))*
     x    tanh((qz(kk)-qz(nz/2))/Lo)*(cosh((qz(kk)-qz(nz/2))/Lo))**(-2)
         dvz = v1*sin(PI*(qx(ii)-qx(nx/2))/qx(nx/2))*
     x       (cosh((qz(kk)-qz(nz/2))/Lo))**(-2)
         
         vp(l,1) = vsw*(tanh((qz(k)-qz(nz/2))/(Lo))) + vx !+dvx
         vp(l,2) = vy 
         vp(l,3) = vz !+dvz 

         if (xp(l,3) .gt. qz(nz/2)) mix_ind(l) = 1
         if (xp(l,3) .le. qz(nz/2)) mix_ind(l) = 0
            
         do 95 m=1,3
            vp1(l,m) = vp(l,m)
            input_E = input_E + 
     x           0.5*m_arr(l)*(vp(l,m)*km_to_m)**2 /beta
            input_p(m) = input_p(m) + m_arr(l)*vp(l,m) / beta
 95      continue
 60   continue
c      endif
c      enddo




 
      write(*,*) 'get interp weights...'
      call get_interp_weights(xp)
      write(*,*) 'update_np...'
      call update_np(np)
      write(*,*) 'update_up...'
      call update_up(vp,np,up)
      write(*,*) 'update_up complete...'

   

      return
      end SUBROUTINE sw_part_setup_maxwl
c----------------------------------------------------------------------

c----------------------------------------------------------------------
      SUBROUTINE sw_part_setup(np,vp,vp1,xp,input_p,up)
c----------------------------------------------------------------------
      !include 'incurv.h'

      real np(nx,ny,nz)
      real vp(Ni_max,3)
      real vp1(Ni_max,3)
      real xp(Ni_max,3)
      real input_p(3)
      real up(nx,ny,nz,3)

c      Ni_tot = 120000

      do 10 l = 1,Ni_tot

         xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx)-qx(1))
         xp(l,2) = qy((ny/2)-8)+(1.0-pad_ranf())*(qy((ny/2)+8)-
     x                                        qy((ny/2)-8))
         xp(l,3) = qz((nz/2)-8)+(1.0-pad_ranf())*(qz((nz/2)+8)-
     x                                        qz((nz/2)-8))

c         i = nint(nx*pad_ranf())
c         j = nint((ny/2) + 16.*(0.5-pad_ranf()))
c         k = nint((nz/2) + 16.*(0.5-pad_ranf()))
cc         write(*,*) 'l...',l,i,j,k

c         xp(l,1) = qx(i)+dx*(0.5-pad_ranf())
c         xp(l,2) = qy(j)+dy*(0.5-pad_ranf())
c         xp(l,3) = qz(k)+dz_grid(k)*(0.5-pad_ranf())

         vp(l,1) = -vsw
         vp(l,2) = 0.0
         vp(l,3) = 0.0


         ijkp(l,1) = nint(xp(l,1)/dx) !particle grid location index
         ijkp(l,2) = nint(xp(l,2)/dy)
         
         k=1
         do 50 while(xp(l,3) .gt. qz(k)) !find k on non-uniform 
            ijkp(l,3) = k       !grid
            k=k+1
 50      continue
         k=ijkp(l,3)
         if (xp(l,3) .gt. (qz(k)+(dz_grid(k)/2))) then
            ijkp(l,3) = k+1
         endif
         do 45 m=1,3
            vp1(l,m) = vp(l,m)
            input_E = input_E + 
     x           0.5*mion*(vp(l,m)*km_to_m)**2 /beta
            input_p(m) = input_p(m) + mion*vp(l,m) / beta
 45      continue
 10      continue

        
 
c      write(*,*) 'get interp weights...'
      call get_interp_weights(xp)
c      write(*,*) 'update_np...'
      call update_np(np)
c      write(*,*) 'update_up...'
      call update_up(vp,np,up)
c      write(*,*) 'update_up complete...'

      return
      end SUBROUTINE sw_part_setup
c----------------------------------------------------------------------

