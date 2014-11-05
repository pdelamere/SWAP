      MODULE chem_rates

      USE global
      USE inputs
      USE gutsp
      

      contains


c----------------------------------------------------------------------
      SUBROUTINE Neut_Center(cx,cy,cz)
c Locates the cartisian coords of the center of the neutral cloud
c at a given time t.
c----------------------------------------------------------------------
      !include 'incurv.h'

c      t = m*dt + tstart          !0.2 reflects canister evacuation time
c      cx = qx(ri) + vsat*(t-0.2) !release point + cloud expansion
c      cx = qx(ri) + vsat*t       !release point 
c      cy = qy(rj) + dy/1e10      !second term to avoid division
c      cz = qz(rk)                !by zero.  That is to avoid

      x0 = dx/2
      y0 = dy/2
      z0 = dz_grid(nz/2)/2
 
      cx = qx(ri) + x0
      cy = qy(rj) + y0
      cz = qz(rk) + z0

                                 !centering the sat track on 
                                 !whole grid points, otherwise
                                 !danger of r = 0.
      return
      end SUBROUTINE Neut_Center
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE Ionize_pluto(np,vp,vp1,xp,m_tstep,input_p,up)
c Ionizes the neutral cloud with a 28 s time constant and fill particle
c arrays, np, vp, up (ion particle density, velocity, 
c and bulk velocity).   
c----------------------------------------------------------------------
CVD$R VECTOR
      !!include 'incurv.h'

      real np(nx,ny,nz),
     x     vp(Ni_max,3),
     x     vp1(Ni_max,3),
     x     xp(Ni_max,3),
     x     input_p(3),
     x     up(nx,ny,nz,3)

      real function pad_ranf      
      real uptot1(3),uptot2(3)
c      integer*4 dNi         !# of new born ions created in dt
c      real vol              !volume of grid cell

      real r                !dist of particle to neutral cloud center
      real t                !run time
      real v                !neutral cloud velocity, r/t
      real cx,cy,cz         !x,y,z coord of neutral cloud center
      real theta,phi        !spherical coords used to get velocities
      integer flg           !flag for while loop
      real zsum             !sum along z to find cell for new_borns
      real rnd              !random number
      real n_source,vol_source
      real r_xyz(200)       !radial distance
      real src(200)         !particle source distribution
      real Nofr(200)        !number of neutrals as func of r
      real nnofr(200)       !neutral density vs. r
      real npofr(200)       !plasma production rate vs. r
      real Np_total         !total number of ions /s
      real vol(200)         !volume of shell vs. r
      real intgl            !integral
      integer rijk
      real ddni

      real minden !minimum wake density of 1 particle per cell

      integer*4 ion_cnt(nx,ny,nz)  !keeps running count of ions 
                                   !in cell for calculating the bulk
                                   !flow velocity
      real vsum(nx,ny,nz,3) !total particle velocity in cell of the
                            !new born ions...for get bulk flow

c      parameter(tau=1.2e9)            !photoionization time constant
c                                     !3% per second
c      if (my_rank .eq. 0) then

      call Neut_Center(cx,cy,cz)

c get source density

      src(2) = Qo   !mol/s at Pluto
c      vrad = 0.4   !km/s
c      N_o = 5e34    !Number at source
c      Rp = 1200.0

      do i = 1,200 
         r_xyz(i) = i*dx 
      enddo

      vol(1) = (4./3.)*pi*r_xyz(1)**3
      do i = 2,199 
         vol(i) = (4./3.)*pi*(r_xyz(i+1)**3 - r_xyz(i)**3)
      enddo
      vol(200) = vol(199)


c      do i = 1,200
c         Nofr(i) =  (N_o - src(i)*tau_photo)*
c     x                exp(-r_xyz(i)/(tau_photo*vrad)) + 
c     x                src(i)*tau_photo
c         nnofr(i) = Nofr(i)/vol(i)
c         npofr(i) = nnofr(i)/tau_photo
c      enddo

      do i = 1,200
         nnofr(i) = Qo/(4*pi*r_xyz(i)**2*vrad)
c         write(*,*) 'nnofr1...',nnofr(i)
c         nnofr(i) = 6e13*(r_xyz(i)/1500.)**(-15)*1e15
c         write(*,*) 'nnofr2...',nnofr(i),r_xyz(i)
c         if (nnofr(i) .gt. 1e7*1e15) then 
c            nnofr(i) = 1e7*1e15
c         endif
         npofr(i) = nnofr(i)/tau_photo
      enddo

     


      l1=Ni_tot + 1  !beginning array element for new borns
      if (Ni_tot .eq. 0) then l1 = 1
      
      if ((Ni_tot+dNi) .gt. Ni_max) then Ni_tot = Ni_max
c      ddni=nint(dNi*(m_tstep/100.))
c      if (m_tstep .gt. 100) then ddni = dNi
c      write(*,*) 'ddni...',ddni,m_tstep,dNi*m_tstep/100.
      if ((Ni_tot+dNi) .le. Ni_max) then 
         Ni_tot = Ni_tot + dNi
c         Ni_tot = Ni_tot + ddni
      endif
c      write(*,*) 'Ni_tot in Ionize....',Ni_tot
      
      if (Ni_tot .le. Ni_max) then
         do 20 l=l1,Ni_tot !initialize new borns
            phi = 2.0*pi*pad_ranf()
            flg = 0
            do 30 while (flg .eq. 0)
               theta = pi*pad_ranf()
               f = sin(theta)
               rnd = pad_ranf()
               if (f .ge. rnd) flg = 1
 30         continue

c            flg = 0
c            do 40 while (flg .eq. 0)
            r = S_radius*dx*pad_ranf()
c               rijk = nint(r/dx)
            v = 0.0 !(10.0*pad_ranf())
c               f = exp(-(v-vo)**2 / vth**2)
c               f = npofr(rijk)/npofr(2)
c               rnd = pad_ranf()
c               if (f .ge. rnd) then 
c            flg = 1
c                  vp(l,1) = vsat + v*cos(phi)*sin(theta)
            vp(l,1) = v*cos(phi)*sin(theta)
            vp(l,2) = v*sin(phi)*sin(theta)
            vp(l,3) = v*cos(theta)
c                  r=v*t
c                  r = 15.0*dx*pad_ranf()
            xp(l,1) = cx + r*cos(phi)*sin(theta)
            xp(l,2) = cy + r*sin(phi)*sin(theta)
            xp(l,3) = cz + r*cos(theta)
c                  xp(l,1) = cx + 5.0*dx*(0.5-pad_ranf()) 
c                  xp(l,2) = cy + 5.0*dy*(0.5-pad_ranf()) 
c                  xp(l,3) = cz + 5.0*dz_cell(rk)*(0.5-pad_ranf()) 
            do 45 m=1,3
               vp1(l,m) = vp(l,m)
               input_E = input_E + 
     x              0.5*m_arr(l)*(vp(l,m)*km_to_m)**2 /beta
               input_p(m) = input_p(m) + m_arr(l)*vp(l,m) / beta
 45         continue
c                  endif

c 40         continue



            ijkp(l,1) = nint(xp(l,1)/dx) !particle grid location index
            ijkp(l,2) = nint(xp(l,2)/dy)

            k=1
            do 50 while(xp(l,3) .gt. qz(k))  !find k on non-uniform 
               ijkp(l,3) = k                 !grid
               k=k+1
 50         continue
            k=ijkp(l,3)
            if (xp(l,3) .gt. (qz(k)+(dz_grid(k)/2))) then
               ijkp(l,3) = k+1
            endif

c         write(*,*) 'xp,qz....',l,k,ijkp(l,3),xp(l,3),qz(k),
c     x                          qz(ijkp(l,3))

            mrat(l) = 1.0/m_pu
            m_arr(l) = mion*m_pu
 20      continue
         
         do 60 l = l1,Ni_tot
            if ((ijkp(l,1) .gt. nx) .or. (ijkp(l,2) .gt. ny) .or. 
     x           (ijkp(l,3) .gt. nz)) then
               call remove_ion(xp,vp,vp1,l)
               
            endif
 60      continue
         
         
c ionize neutral particles assuming constant ionization rate
c need to be careful ion ions formed near boundary...

      endif

c      endif   !for pickup on one processor
c maintain minimum wake density

c      minden = (2.0)/(beta*dx*dy*delz) 
      minden = nf_init/10.
c      print *,'minden...',minden,beta,dx*dy*delz
      do i = 2,nx-1
         do j = 2,ny-1
            do k = 2,nz-1
               if (np(i,j,k) .le. minden) then
c                  write(*,*) 'np...',np(i,j,k),minden
                                !choose a random processor
                  if (my_rank .eq. nint(pad_ranf()*procnum)) then

                     l=Ni_tot + 1 !beginning array element for new borns    
                  
c                     write(*,*) ' ',procnum
c           write(*,*) 'Wake depletion, adding particle',procnum,l1,i,j,k

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
                  endif

               endif
            enddo
         enddo
      enddo


c      write(*,*) 'Ni_tot after wake....',Ni_tot

      call get_interp_weights(xp)
      call update_np(np)
      call update_up(vp,np,up)


      
      return
      end SUBROUTINE Ionize_pluto
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE Ionize_pluto_mp(np,vp,vp1,xp,m_tstep,input_p,up)
c Ionizes the neutral cloud with a 28 s time constant and fill particle
c arrays, np, vp, up (ion particle density, velocity, 
c and bulk velocity).   
c----------------------------------------------------------------------
CVD$R VECTOR
      !!include 'incurv.h'

      real np(nx,ny,nz),
     x     vp(Ni_max,3),
     x     vp1(Ni_max,3),
     x     xp(Ni_max,3),
     x     input_p(3),
     x     up(nx,ny,nz,3)

      real function pad_ranf      
      real uptot1(3),uptot2(3)
c      integer*4 dNi         !# of new born ions created in dt
c      real vol              !volume of grid cell

      real r                !dist of particle to neutral cloud center
      real t                !run time
      real v                !neutral cloud velocity, r/t
      real cx,cy,cz         !x,y,z coord of neutral cloud center
      real theta,phi        !spherical coords used to get velocities
      integer flg           !flag for while loop
      real zsum             !sum along z to find cell for new_borns
      real rnd              !random number
      real n_source,vol_source
c      real r_xyz       !radial distance
c      real src(200)         !particle source distribution
c      real Nofr(200)        !number of neutrals as func of r
      real nnofr       !neutral density vs. r
      real npofr       !plasma production rate vs. r
      real Np_total         !total number of ions /s
      real vol         !volume of shell vs. r
      real vol_shell
      real vol_shell_min
      real intgl            !integral
      integer rijk
      real ddni
      integer cnt, l1

c      integer*4 ion_cnt(nx,ny,nz)  !keeps running count of ions 
                                   !in cell for calculating the bulk
                                   !flow velocity
c      real vsum(nx,ny,nz,3) !total particle velocity in cell of the
                            !new born ions...for get bulk flow

      call Neut_Center(cx,cy,cz)

c get source density

      vol = dx**3
      cnt = 0
      l1 = Ni_tot+1
      do i = 2,nx-1
         do j = 2,ny-1
            do k = 2,nz-1
               r = sqrt((qx(i)-cx)**2 + (qy(j)-cy)**2 + (qz(k)-cz)**2)
               if (r .le. dx*S_radius) then
                  nnofr = Qo/(4*pi*r**2*vrad)
                  npofr = vol*beta*nnofr*dt/tau_photo/procnum
c                  write(*,*) 'npofr...',npofr
                  if (npofr .ge. 1) then
                     
                     do ll = 1,nint(npofr)
                        l = Ni_tot + 1
                        vp(l,1) = 0.0
                        vp(l,2) = 0.0
                        vp(l,3) = 0.0                        

                       xp(l,1) = qx(i) + (pad_ranf()-0.5)*1.0*dx
                       xp(l,2) = qy(j) + (pad_ranf()-0.5)*1.0*dy
                       xp(l,3) = qz(k) + (pad_ranf()-0.5)*1.0*dz_grid(k)
                        
                        ijkp(l,1) = nint(xp(l,1)/dx) !particle grid location index
                        ijkp(l,2) = nint(xp(l,2)/dy)

                        kk=1
                       do 15 while((xp(l,3).gt.qz(kk)).and.(kk .le. nz)) !find k
                           ijkp(l,3) = kk !grid
                           kk=kk+1
 15                     continue
                        kk=ijkp(l,3)
                        if (xp(l,3) .gt. (qz(kk)+(dz_grid(kk)/2))) then
                           ijkp(l,3) = kk+1
                        endif
                        
                        mrat(l) = 1.0/m_pu
                        m_arr(l) = mion*m_pu
                        Ni_tot = l
                        cnt = cnt + 1
                        do m=1,3
                           vp1(l,m) = vp(l,m)
                           input_E = input_E + 
     x                          0.5*m_arr(l)*(vp(l,m)*km_to_m)**2 /beta
                         input_p(m) = input_p(m) + m_arr(l)*vp(l,m)/beta
                        enddo                     

                     enddo
                  endif
                  if (npofr .lt. 1) then
                  
                     if (npofr .gt. pad_ranf()) then
                        l = Ni_tot + 1
                        vp(l,1) = 0.0
                        vp(l,2) = 0.0
                        vp(l,3) = 0.0                        

                       xp(l,1) = qx(i) + (pad_ranf()-0.5)*1.0*dx
                       xp(l,2) = qy(j) + (pad_ranf()-0.5)*1.0*dy
                       xp(l,3) = qz(k) + (pad_ranf()-0.5)*1.0*dz_grid(k)
                        
                        ijkp(l,1) = nint(xp(l,1)/dx) !particle grid location index
                        ijkp(l,2) = nint(xp(l,2)/dy)
                        
                        kk=1
                       do 16 while((xp(l,3).gt.qz(kk)).and.(kk .le. nz)) !find k
                           ijkp(l,3) = kk !grid
                           kk=kk+1
 16                     continue
                        kk=ijkp(l,3)
                        if (xp(l,3) .gt. (qz(kk)+(dz_grid(kk)/2))) then
                           ijkp(l,3) = kk+1
                        endif

                        mrat(l) = 1.0/m_pu
                        m_arr(l) = mion*m_pu
                        Ni_tot = l
                        cnt = cnt + 1
                        do m=1,3
                           vp1(l,m) = vp(l,m)
                           input_E = input_E + 
     x                          0.5*m_arr(l)*(vp(l,m)*km_to_m)**2 /beta
                         input_p(m) = input_p(m) + m_arr(l)*vp(l,m)/beta
                        enddo                     
                     endif
                  endif
               endif
            enddo
         enddo
      enddo

      
      write(*,*) 'total new ions....',my_rank,cnt,dNi         
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

c      stop

      do 60 l = l1,Ni_tot
         if ((ijkp(l,1) .gt. nx) .or. (ijkp(l,2) .gt. ny) .or. 
     x        (ijkp(l,3) .gt. nz)) then
            write(*,*) 'removing ion...',my_rank,ijkp(l,:)
            call remove_ion(xp,vp,vp1,l)
            
         endif
 60   continue

c      stop

c      minden = nf_init/10.
c      do i = 2,nx-1
c         do j = 2,ny-1
c            do k = 2,nz-1
cc               if (np(i,j,k) .le. minden) then
c               do 62 while (np(i,j,k) .le. minden) 
c                  write(*,*) 'minden...',np(i,j,k),i,j,k
c                                !choose a random processor
c                  if (my_rank .eq. nint(pad_ranf()*procnum)) then

c                     l=Ni_tot + 1 !beginning array element for new borns    
                  
c                     vp(l,1) = up(i,j,k,1)
c                     vp(l,2) = up(i,j,k,2)
c                     vp(l,3) = up(i,j,k,3)
c                     xp(l,1) = qx(i) + (0.5-pad_ranf())*dx
c                     xp(l,2) = qy(j) + (0.5-pad_ranf())*dy
c                     xp(l,3) = qz(k) + (0.5-pad_ranf())*dz_grid(k)

c                     ijkp(l,1) = nint(xp(l,1)/dx) !particle grid location index
c                     ijkp(l,2) = nint(xp(l,2)/dy)

c                     kk=1
c                     do 100 while(xp(l,3) .gt. qz(kk)) !find k on non-uniform 
c                        ijkp(l,3) = kk !grid
c                        kk=kk+1
c 100                 continue
c                     kk=ijkp(l,3)
c                     if (xp(l,3) .gt. (qz(kk)+(dz_grid(kk)/2))) then
c                        ijkp(l,3) = kk+1
c                     endif
c                     mrat(l) = 1.0
c                     m_arr(l) = mion
c                     Ni_tot = Ni_tot + 1
c                  endif

c               endif
c 62            continue
c            enddo
c         enddo
c      enddo


c      write(*,*) 'Ni_tot after wake....',Ni_tot

      call get_interp_weights(xp)
      call update_np(np)
      call update_up(vp,np,up)


      
      return
      end SUBROUTINE Ionize_pluto_mp
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE Ionize_sw_mp(np,vp,vp1,xp,m_tstep,input_p,up)
c Ionizes the neutral cloud with a 28 s time constant and fill particle
c arrays, np, vp, up (ion particle density, velocity, 
c and bulk velocity).   
c----------------------------------------------------------------------
      !include 'incurv.h'

      real np(nx,ny,nz),
     x     vp(Ni_max,3),
     x     vp1(Ni_max,3),
     x     xp(Ni_max,3),
     x     input_p(3),
     x     up(nx,ny,nz,3)

      real function pad_ranf      
      real uptot1(3),uptot2(3)
c      integer*4 dNi         !# of new born ions created in dt
c      real vol              !volume of grid cell

      real r                !dist of particle to neutral cloud center
      real t                !run time
      real v                !neutral cloud velocity, r/t
      real cx,cy,cz         !x,y,z coord of neutral cloud center
      real theta,phi        !spherical coords used to get velocities
      integer flg           !flag for while loop
      real zsum             !sum along z to find cell for new_borns
      real rnd              !random number
      real n_source,vol_source
c      real r_xyz       !radial distance
c      real src(200)         !particle source distribution
c      real Nofr(200)        !number of neutrals as func of r
      real nnofr       !neutral density vs. r
      real npofr       !plasma production rate vs. r
      real Np_total         !total number of ions /s
      real vol         !volume of shell vs. r
      real vol_shell
      real vol_shell_min
      real intgl            !integral
      integer rijk
      real ddni
      integer cnt, l1
      
      real ndot
      ndot = 1e-6*nf_init*omega_p !normalized ionization rate

      dNi = ndot*dt*beta*beta_pu*dx*dy*delz*nz

      if (dNi .lt. 1.0) then 
         if (dNi .gt. pad_ranf()) then 
            dNi = 1.0
            write(*,*) 'new ions...',dNi
         endif
      endif

      dNi = nint(dNi)

c      write(*,*) 'dNi...',dNi,ndot*dt*beta*beta_pu*dx*dy*delz,Ni_tot
      l1 = Ni_tot+1

      do l = l1,l1+dNi-1

c         write(*,*) 'new ions....',l,dNi,Ni_tot
         theta = pad_ranf()*2*PI
         vp(l,1) = vsw+vsw*cos(theta) !+dvx
         vp(l,2) = vsw*sin(theta) !+dvz 
         vp(l,3) = 0.0

         xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))
         xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))
         xp(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz-1)-qz(1))

         i=0
 31      continue
         i = i + 1
         if (xp(l,1) .gt. qx(i)) go to 31 !find i on non-uniform 
         i = i-1
         ijkp(l,1)= i

         ijkp(l,2) = floor(xp(l,2)/dy) 
         
         k=0
 30      continue
         k = k + 1
         if (xp(l,3) .gt. qz(k)) go to 30 !find k on non-uniform 
         k = k-1
         ijkp(l,3)= k

         !add He++

c         if (pad_ranf() .lt. 0.1) then 
c            mrat(l) = 1.0/2.0
c            m_arr(l) = mion*2.0
c            beta_p(l) = beta_pu
c         endif

         !add protons

c         if (pad_ranf() .ge. 0.1) then 
            mrat(l) = 1.0
            m_arr(l) = mion
            beta_p(l) = beta_pu
c         endif

         

         do m=1,3
            vp1(l,m) = vp(l,m)
            input_E = input_E + 
     x           0.5*m_arr(l)*(vp(l,m)*km_to_m)**2 /
     x           beta*beta_p(l)
            input_p(m) = input_p(m) + m_arr(l)*vp(l,m)/
     x           beta*beta_p(l)
         enddo                     

      enddo
      
c      write(*,*) 'total new ions in ionize....',my_rank,cnt,dNi         
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)


      Ni_tot = Ni_tot + dNi


c      do 60 l = l1,Ni_tot
cc         write(*,*) 'l1...',l,l1,Ni_tot
c         if ((ijkp(l,1) .gt. nx) .or. (ijkp(l,2) .gt. ny) .or. 
c     x        (ijkp(l,3) .gt. nz)) then
c            write(*,*) 'removing ion in ionize...',my_rank,ijkp(l,:)
c            call remove_ion(xp,vp,vp1,l)
            
c         endif
c 60   continue

      call get_interp_weights(xp)
      call update_np(np)
      call update_up(vp,np,up)
      
      return
      end SUBROUTINE Ionize_sw_mp
c----------------------------------------------------------------------




c----------------------------------------------------------------------
      SUBROUTINE Ionize_sw_mp_1(np,vp,vp1,xp,m_tstep,input_p,up)
c Ionizes the neutral cloud with a 28 s time constant and fill particle
c arrays, np, vp, up (ion particle density, velocity, 
c and bulk velocity).   
c----------------------------------------------------------------------

      real np(nx,ny,nz),
     x     vp(Ni_max,3),
     x     vp1(Ni_max,3),
     x     xp(Ni_max,3),
     x     input_p(3),
     x     up(nx,ny,nz,3)

      real function pad_ranf      
      real uptot1(3),uptot2(3)
c      integer*4 dNi         !# of new born ions created in dt
c      real vol              !volume of grid cell

      real r                !dist of particle to neutral cloud center
      real t                !run time
      real v                !neutral cloud velocity, r/t
      real cx,cy,cz         !x,y,z coord of neutral cloud center
      real theta,phi        !spherical coords used to get velocities
      integer flg           !flag for while loop
      real zsum             !sum along z to find cell for new_borns
      real rnd              !random number
      real n_source,vol_source
c      real r_xyz       !radial distance
c      real src(200)         !particle source distribution
c      real Nofr(200)        !number of neutrals as func of r
      real nnofr       !neutral density vs. r
      real npofr       !plasma production rate vs. r
      real Np_total         !total number of ions /s
      real vol         !volume of shell vs. r
      real vol_shell
      real vol_shell_min
      real intgl            !integral
      integer rijk
      real ddni
      integer cnt, l1

      call Neut_Center(cx,cy,cz)

      dNi = 10.0
      l1 = Ni_tot+1

      do l = l1,l1+dNi
         theta = pad_ranf()*2*PI
         vp(l,1) = vsw+57.0*cos(theta) !+dvx
         vp(l,2) = 57.0*sin(theta) !+dvz 
         vp(l,3) = 0.0

c         vp(l,1) = 0.0
c         vp(l,2) = 0.0
c         vp(l,3) = 0.0                        

c         xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))


         flg = 0
         do 10 while (flg .eq. 0) 
         xp(l,1) = qx(nx/2-40)+(1.0-pad_ranf())*(80*dx)
         xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))
         xp(l,3) = qz(nz/2-40)+(1.0-pad_ranf())*(80*delz)

         r = sqrt((xp(l,1)-cx)**2 + (xp(l,2)-cy)**2 + (xp(l,3)-cz)**2)

         rnd = pad_ranf()
         if (exp(-r**2/(20*dx)**2) .gt. rnd) then
            flg = 1
         endif

         
 10      enddo

         i=0
 31      continue
         i = i + 1
         if (xp(l,1) .gt. qx(i)) go to 31 !find i on non-uniform 
         i = i-1
         ijkp(l,1)= i


         ijkp(l,2) = floor(xp(l,2)/dy) 
         
         k=0
 30      continue
         k = k + 1
         if (xp(l,3) .gt. qz(k)) go to 30 !find k on non-uniform 
         k = k-1
         ijkp(l,3)= k


c         ijkp(l,1) = nint(xp(l,1)/dx) !particle grid location index
c         ijkp(l,2) = nint(xp(l,2)/dy)

c         kk=1
c         do 15 while((xp(l,3).gt.qz(kk)).and.(kk .le. nz)) !find k
c            ijkp(l,3) = kk      !grid
c            kk=kk+1
c 15      continue
c         kk=ijkp(l,3)
c         if (xp(l,3) .gt. (qz(kk)+(dz_grid(kk)/2))) then
c            ijkp(l,3) = kk+1
c         endif

         if (pad_ranf() .lt. 0.5) then 
            mrat(l) = ion_amu/m_pu
            m_arr(l) = 1.67e-27*m_pu
         endif


         if (pad_ranf() .ge. 0.5) then 
            mrat(l) = ion_amu/48.0
            m_arr(l) = 1.67e-27*48.0
         endif

         

         do m=1,3
            vp1(l,m) = vp(l,m)
            input_E = input_E + 
     x           0.5*m_arr(l)*(vp(l,m)*km_to_m)**2 /beta
            input_p(m) = input_p(m) + m_arr(l)*vp(l,m)/beta
         enddo                     

      enddo
      
      write(*,*) 'total new ions....',my_rank,cnt,dNi         
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      Ni_tot = Ni_tot + dNi

      do 60 l = l1,Ni_tot
         if ((ijkp(l,1) .gt. nx) .or. (ijkp(l,2) .gt. ny) .or. 
     x        (ijkp(l,3) .gt. nz)) then
            write(*,*) 'removing ion...',my_rank,ijkp(l,:)
            call remove_ion(xp,vp,vp1,l)
            
         endif
 60   continue

      call get_interp_weights(xp)
      call update_np(np)
      call update_up(vp,np,up)
      
      return
      end SUBROUTINE Ionize_sw_mp_1
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE res_chex(xp,vp)
c----------------------------------------------------------------------
      !!include 'incurv.h'

      real xp(Ni_max,3)
      real vp(Ni_max,3)

      real cx,cy,cz,r,vrel
      real nn
      real chex_tau,chex_prob
      real sigma_chex

      PARAMETER (sigma_chex = 8e-26)  !10 A^2 in units of km^2
c      PARAMETER (pwl = 12)
c      PARAMETER (Ncol = 6e16*1e10)  !km^-2

      call Neut_Center(cx,cy,cz)
           
      do l = 1,Ni_tot 
c         if (mrat(l) .eq. 1.0) then 
c            r = sqrt((xp(l,1)-cx)**2 + (xp(l,2)-cy)**2 + 
c     x           (xp(l,3)-cz)**2)
            vrel = sqrt(vp(l,1)**2 + vp(l,2)**2 + vp(l,3)**2)
            nn = 10000e15 !Qo/(4*pi*r**2*vrad)
            chex_tau = 1./(nn*sigma_chex*vrel)
            chex_prob = dt/chex_tau
c            write(*,*) 'chex...',chex_prob,nn,vrel,dt
            
            if (pad_ranf() .lt. chex_prob) then
               write(*,*) 'chex...',l,chex_tau,chex_prob
               vp(l,:) = 0.0
               mrat(l) = 16.0/m_pu
               m_arr(l) = 1.67e-27*m_pu
            endif
c         endif
      enddo
      
      return
      end SUBROUTINE res_chex
c----------------------------------------------------------------------

      end MODULE chem_rates
