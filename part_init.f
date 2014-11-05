      MODULE part_init


      USE global
      USE dimensions
      USE misc
      USE gutsp
      USE mpi

      contains


c----------------------------------------------------------------------
      SUBROUTINE Energy_diag(vp,b0,b1,E,Evp,Euf,EB1,EB1x,EB1y,EB1z,
     x                       EE,EeP,nu,up,np)
c----------------------------------------------------------------------
      !include 'incurv.h'

      real vp(Ni_max,3),
c     x     uf(nx,ny,nz,3),
c     x     nf(nx,ny,nz),
     x     b0(nx,ny,nz,3),
     x     b1(nx,ny,nz,3),
     x     E(nx,ny,nz,3),
c     x     etemp(nx,ny,nz),
     x     nu(nx,ny,nz),
     x     up(nx,ny,nz,3),
     x     np(nx,ny,nz)

      real m_q
     

      real Evp                  !kinetic energy of particles
      real Euf                  !kinetic energy of fluid flow
      real EB1,EB1x,EB1y,EB1z   !Magnetic field energy 
      real EE                   !Electric field energy 
      real EeP                  !Electron pressure energy
      real total_E              !total energy
      real aveEvp               !average particle energy
      real norm_E               !normalized energy
      real vol                  !volume of cell
      real denf                 !fluid density

      real recvbuf
      integer count
      count = 1

      m_q = mion/q

      Euf = 0.0
      EB1 = 0.0
      EB1x = 0.0
      EB1y = 0.0
      EB1z = 0.0
      EE = 0.0
      EeP = 0.0                
      do 10 i=1,nx-1
         j = 2
c         do 10 j=1,ny-1
            do 10 k=1,nz-1
               vol = dx_cell(i)*dy_cell(j)*dz_cell(k)*km_to_m**3
               EB1x = EB1x + (vol/(2.0*mu0))*(m_q*b1(i,j,k,1))**2 
               EB1y = EB1y + (vol/(2.0*mu0))*(m_q*b1(i,j,k,2))**2 
               EB1z = EB1z + (vol/(2.0*mu0))*(m_q*b1(i,j,k,3))**2 
c               EeP = EeP + kboltz*etemp(i,j,k)
               do 10 m=1,3
                  denf = np(i,j,k)/(km_to_m**3)
                  Euf = Euf + 0.5*mO*denf*vol*(up(i,j,k,m)*km_to_m)**2
                  EB1 = EB1 + 
     x              (vol/(2.0*mu0))*(m_q*(b1(i,j,k,m)-b0(i,j,k,m)))**2
                  EE = EE + (epsilon*vol/2.0)*
     x                      (m_q*E(i,j,k,m)*km_to_m)**2
 10               continue

c      input_EeP = input_EeP + EeP

c      write(*,*) 'Energy diag...',Ni_tot,m_arr(2000000)
 
      Evp = 0.0
      do 15 l=1,Ni_tot
         do 15 m=1,3
            Evp = Evp + 0.5*m_arr(l)*(vp(l,m)*km_to_m)**2 /beta
 15   continue

c      write(*,*) 'Energy diag 2...',Ni_tot,m_arr(2000000)

      call MPI_Barrier(MPI_COMM_WORLD,ierr)

      call MPI_ALLREDUCE(Evp,recvbuf,count,
     x     MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)

      S_Evp = recvbuf

      call MPI_ALLREDUCE(input_E,recvbuf,count,
     x     MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)

      S_input_E = (recvbuf)


      total_E = S_Evp+EE+EB1
      aveEvp = S_Evp/S_input_E

c      write(*,*) 'Input energy (J).............',S_input_E
cc      write(*,*) 'Input EeP energy (J).........',input_EeP
c      write(*,*) 'Total vp energy (J)..........',S_Evp
c      write(*,*) 'Total up energy (J)..........',Euf
c      write(*,*) 'Total B energy (J)...........',EB1/S_input_E
c      write(*,*) 'Total E energy (J)...........',EE/S_input_E
cc      write(*,*) 'Total EeP energy (J).........',EeP
c      write(*,*) 'Total energy (J).............',total_E
cc      write(*,*) 'Total energy w/ eP (J).......',total_E+EeP
c      write(*,*) 'Energy thru boundaries.......',bndry_Eflux/S_input_E
c      write(*,*) 'Normalized particle energy...',aveEvp
      if (my_rank .eq. 0) then
      write(*,*) 'Normalized energy............',total_E/S_input_E,
     x   my_rank
      write(*,*) 'Normalized energy (bndry)....',
     x                (total_E)/(S_input_E+bndry_Eflux)
      endif
c      write(*,*) 'Normalized energy (no b1z)...',(S_Evp+Euf+EE+EB1x+
c     x                                            EB1y)/S_input_E
cc      write(*,*) 'Normalized energy (w/ eP)....',
cc     x                             (total_E+EeP)/(input_E + input_EeP)
c      write(*,*) ' '

      norm_E = total_E/S_input_E

c      if (prev_Etot .eq. 0.0) then prev_Etot = norm_E
c      do 20 i=1,nx 
c         do 20 j=1,ny
c            do 20 k=1,nz
c               nu(i,j,k) = nu(i,j,k) + 
c     x                 nu(i,j,k)*2.0*((norm_E - prev_Etot)/norm_E)
c 20            continue
      prev_Etot = norm_E

      return
      end SUBROUTINE Energy_diag
c----------------------------------------------------------------------



c----------------------------------------------------------------------
      SUBROUTINE sw_part_setup_maxwl_mix(np,vp,vp1,xp,input_p,up,
     x     np_t_flg,np_b_flg)
c----------------------------------------------------------------------
      !include 'incurv.h'
c      include mpif.h

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

c      nprat = np_bottom/np_top
c      Ni_tot_heavy = Ni_tot

c add cold populations first

      do 10 l = 1,Ni_tot

         xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))
         xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))
         
c         flg = 0
c         do 20 while (flg .eq. 0)
         xp(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz-1)-qz(1))
cc            rnd = ((1.0-tanh((xp(l,3)-qz(nz/2))/(Lo)))/2.0)  !for bottom
cc            rnd = (1.0+tanh((xp(l,3)-qz(nz/2))/(Lo)))/2.0 !for top
c            if (pad_ranf() .le. rnd) flg = 1
cc            if (xp(l,3) .ge. qz(nz/2)) then
cc               np_t_flg(l) = 1
cc               m_arr(l) = m_top
cc               mrat(l) = mion/m_top
cc            endif
cc            if (xp(l,3) .lt. qz(nz/2)) then
c            np_b_flg(l) = 1
         m_arr(l) = mion
         mrat(l) = 1.0

c            endif
            
c 20      continue


c         ijkp(l,1) = floor(xp(l,1)/dx) 

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
         
c         k=1
c         do 30 while(xp(l,3) .gt. qz(k)) !find k on non-uniform 
c            ijkp(l,3) = k       !grid
c            k=k+1
c 30      continue
c         k=ijkp(l,3)
c         if (xp(l,3) .gt. (qz(k)+(dz_grid(k)/2))) then
c            ijkp(l,3) = k+1
c         endif

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
c         dvx = -2.0*v1*cos(PI*(qx(ii)-qx(nx/2))/qx(nx/2))*
c     x    tanh((qz(kk)-qz(nz/2))/Lo)*(cosh((qz(kk)-qz(nz/2))/Lo))**(-2)
c         dvz = v1*sin(PI*(qx(ii)-qx(nx/2))/qx(nx/2))*
c     x       (cosh((qz(kk)-qz(nz/2))/Lo))**(-2)
         
c         vp(l,1) = vsw*(tanh((qz(k)-qz(nz/2))/(Lo))) + vx !+dvx
c         vp(l,1) =  50.0*exp(-(xp(l,3)-qz(nz/4*3))**2/(20.*delz)**2)+vx
c         vp(l,2) =  50.0*exp(-(xp(l,3)-qz(nz/4))**2/(20.*delz)**2)+vy 
         vp(l,1) = 57.0*exp(-(xp(l,3)-qz(nz/2))**2/(10.*delz)**2)+vx
         vp(l,2) = vy 
         vp(l,3) = vz 

         if (xp(l,3) .gt. qz(nz/2)) mix_ind(l) = 1
         if (xp(l,3) .le. qz(nz/2)) mix_ind(l) = 0
            
         do 45 m=1,3
            vp1(l,m) = vp(l,m)
            input_E = input_E + 
     x           0.5*m_arr(l)*(vp(l,m)*km_to_m)**2 /beta
            input_p(m) = input_p(m) + m_arr(l)*vp(l,m) / beta
 45      continue
 10      continue


c add hot and fast population

         Ni_tot_1 = Ni_tot + 1

         Ni_tot = Ni_tot + 0.0*Ni_tot


         do 49 l = Ni_tot_1,Ni_tot

            xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))
            xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))

c            flg = 0
c            do 22 while (flg .eq. 0)
            xp(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz-1)-qz(1))
c               rnd = ((1.0-tanh((xp(l,3)-qz(nz/2))/(Lo)))/2.0) !for bottom
c               rnd = (1.0+tanh((xp(l,3)-qz(nz/2))/(Lo)))/2.0 !for top
c               if (pad_ranf() .le. rnd) flg = 1
c               np_t_flg(l) = 1
            m_arr(l) = mion
            mrat(l) = 1.0
c 22         continue

c            ijkp(l,1) = floor(xp(l,1)/dx) 

            i=0
 51         continue
            i = i + 1
            if (xp(l,1) .gt. qx(i)) go to 51 !find i on non-uniform 
            i = i-1
            ijkp(l,1)= i


            ijkp(l,2) = floor(xp(l,2)/dy) 
         
            k=0
 50         continue
            k = k + 1
            if (xp(l,3) .gt. qz(k)) go to 50 !find k on non-uniform 
            k = k-1
            ijkp(l,3)= k

c            ijkp(l,1) = nint(xp(l,1)/dx) !particle grid location index
c            ijkp(l,2) = nint(xp(l,2)/dy)
            
c            k=1
c            do 50 while(xp(l,3) .gt. qz(k)) !find k on non-uniform 
c               ijkp(l,3) = k    !grid
c               k=k+1
c 50         continue
c            k=ijkp(l,3)
c            if (xp(l,3) .gt. (qz(k)+(dz_grid(k)/2))) then
c               ijkp(l,3) = k+1
c            endif
            
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
c            dvx = -2.0*v1*cos(PI*(qx(ii)-qx(nx/2))/qx(nx/2))*
c     x           tanh((qz(kk)-qz(nz/2))/Lo)*(cosh((qz(kk)-qz(nz/2))/Lo))**(-2)
c            dvz = v1*sin(PI*(qx(ii)-qx(nx/2))/qx(nx/2))*
c     x           (cosh((qz(kk)-qz(nz/2))/Lo))**(-2)
            
c            vp(l,1) = vsw*(tanh((qz(k)-qz(nz/2))/(Lo))) + vx !+dvx
            vp(l,1) = 1.0*vsw + vx !+dvx
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


c add pickup distribution

         Ni_tot_1 = Ni_tot + 1

         Ni_tot = Ni_tot + 0.0*Ni_tot


         do 69 l = Ni_tot_1,Ni_tot

            xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))
            xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))

c            flg = 0
c            do 22 while (flg .eq. 0)
            xp(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz-1)-qz(1))
c               rnd = ((1.0-tanh((xp(l,3)-qz(nz/2))/(Lo)))/2.0) !for bottom
c               rnd = (1.0+tanh((xp(l,3)-qz(nz/2))/(Lo)))/2.0 !for top
c               if (pad_ranf() .le. rnd) flg = 1
c               np_t_flg(l) = 1
            m_arr(l) = mion
            mrat(l) = 1.0
c 22         continue

c            ijkp(l,1) = floor(xp(l,1)/dx) 

            i=0
 71         continue
            i = i + 1
            if (xp(l,1) .gt. qx(i)) go to 71 !find i on non-uniform 
            i = i-1
            ijkp(l,1)= i


            ijkp(l,2) = floor(xp(l,2)/dy) 
            
            k=0
 70         continue
            k = k + 1
            if (xp(l,3) .gt. qz(k)) go to 70 !find k on non-uniform 
            k = k-1
            ijkp(l,3)= k

c            ijkp(l,1) = nint(xp(l,1)/dx) !particle grid location index
c            ijkp(l,2) = nint(xp(l,2)/dy)
            
c            k=1
c            do 70 while(xp(l,3) .gt. qz(k)) !find k on non-uniform 
c               ijkp(l,3) = k    !grid
c               k=k+1
c 70         continue
c            k=ijkp(l,3)
c            if (xp(l,3) .gt. (qz(k)+(dz_grid(k)/2))) then
c               ijkp(l,3) = k+1
c            endif
            
            vth = vth_top
            
            flg = 0
            do 72 while (flg .eq. 0)
               v = (2*vth_max*pad_ranf())-vth_max
               f = exp(-(v)**2 / vth**2)
               rnd = pad_ranf()
               if (f .ge. rnd) then 
                  flg = 1
                  vx = v
               endif
 72         continue
            flg = 0
            do 74 while (flg .eq. 0)
               v = (2*vth_max*pad_ranf())-vth_max
               f = exp(-(v)**2 / vth**2)
               rnd = pad_ranf()
               if (f .ge. rnd) then 
                  flg = 1
                  vy = v
               endif
 74         continue
            flg = 0
            do 76 while (flg .eq. 0)
               v = (2*vth_max*pad_ranf())-vth_max
               f = exp(-(v)**2 / vth**2)
               rnd = pad_ranf()
               if (f .ge. rnd) then 
                  flg = 1
                  vz = v
               endif
 76         continue
            
            ii = ijkp(l,1)
            kk = ijkp(l,3)

            theta = pad_ranf()*2*PI
            vp(l,1) = 1.0*vsw*cos(theta) + vx !+dvx
            vp(l,2) = vy 
            vp(l,3) = 1.0*vsw*sin(theta) + vz        !+dvz 
            
            if (xp(l,3) .gt. qz(nz/2)) mix_ind(l) = 1
            if (xp(l,3) .le. qz(nz/2)) mix_ind(l) = 0
            
            do  m=1,3
               vp1(l,m) = vp(l,m)
               input_E = input_E + 
     x              0.5*m_arr(l)*(vp(l,m)*km_to_m)**2 /beta
               input_p(m) = input_p(m) + m_arr(l)*vp(l,m) / beta
            enddo
            

 69      enddo


c      write(*,*) 'get interp weights...'
      call get_interp_weights(xp)
c      write(*,*) 'update_np...'
      call update_np(np)
c      write(*,*) 'update_up...'
      call update_up(vp,np,up)
c      write(*,*) 'update_up complete...'

   

      return
      end SUBROUTINE sw_part_setup_maxwl_mix
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE load_Maxwellian(np,vp,vp1,xp,input_p,up,vth,Ni_tot_1)
c----------------------------------------------------------------------

      real np(nx,ny,nz)
      real vp(Ni_max,3)
      real vp1(Ni_max,3)
      real xp(Ni_max,3)
      real input_p(3)
      real up(nx,ny,nz,3)
c      real phi,theta,rnd,f,v
c      real rand
c      real vx,vy,vz
c      real dvx,dvz,v1
c      integer np_t_flg(Ni_max)
c      integer np_b_flg(Ni_max)
      integer Ni_tot_1
      real vth

c      integer flg
c      real nprat

      v1 = 1.0

c      np_t_flg(:) = 0
c      np_b_flg(:) = 0

      do 10 l = 1,Ni_tot_1

         xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))
         xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))
         xp(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz-1)-qz(1))
         m_arr(l) = mion
         mrat(l) = 1.0

c         ijkp(l,1) = floor(xp(l,1)/dx) 

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
         
c         vth = vth_bottom

         vx = vsw+vth*sqrt(-alog(pad_ranf()))*cos(PI*pad_ranf())
         vy = vth*sqrt(-alog(pad_ranf()))*cos(PI*pad_ranf())
         vz = vth*sqrt(-alog(pad_ranf()))*cos(PI*pad_ranf())

         ii = ijkp(l,1)
         kk = ijkp(l,3)
         vp(l,1) = -0.0*(exp(-(xp(l,3)-qz(nz/2))**2/(10.*delz)**2)
     x        *exp(-(xp(l,1)-qx(nx/2))**2/(10.*dx)**2))+vx
         vp(l,2) = vy 
         vp(l,3) = vz 

c         if (xp(l,3) .gt. qz(nz/2)) mix_ind(l) = 1
c         if (xp(l,3) .le. qz(nz/2)) mix_ind(l) = 0
            
         do 45 m=1,3
            vp1(l,m) = vp(l,m)
            input_E = input_E + 
     x           0.5*m_arr(l)*(vp(l,m)*km_to_m)**2 /beta
            input_p(m) = input_p(m) + m_arr(l)*vp(l,m) / beta
 45      continue
 10      continue



      call get_interp_weights(xp)
      call update_np(np)
      call update_up(vp,np,up)

      return
      end SUBROUTINE load_Maxwellian
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE load_aniso_Maxwellian(np,vp,vp1,xp,input_p,up,vth, 
     x     Ni_tot_1)
c----------------------------------------------------------------------

      real np(nx,ny,nz)
      real vp(Ni_max,3)
      real vp1(Ni_max,3)
      real xp(Ni_max,3)
      real input_p(3)
      real up(nx,ny,nz,3)
c      real phi,theta,rnd,f,v
c      real rand
c      real vx,vy,vz
c      real dvx,dvz,v1
c      integer np_t_flg(Ni_max)
c      integer np_b_flg(Ni_max)
      integer Ni_tot_1
      real vth

      real aniso
      PARAMETER (aniso = 0.2)

c      integer flg
c      real nprat

      v1 = 1.0

c      np_t_flg(:) = 0
c      np_b_flg(:) = 0

      do 10 l = 1,Ni_tot_1

         xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))
         xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))
         xp(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz-1)-qz(1))
         m_arr(l) = mion
         mrat(l) = 1.0

c         ijkp(l,1) = floor(xp(l,1)/dx) 

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
         
c         vth = vth_bottom

         vx = vsw+vth*sqrt(-alog(pad_ranf()))*cos(PI*pad_ranf())
         vy = vth*sqrt(-alog(pad_ranf()))*cos(PI*pad_ranf())
         vz = vth*sqrt(-alog(pad_ranf()))*cos(PI*pad_ranf())

         ii = ijkp(l,1)
         kk = ijkp(l,3)
         vp(l,1) = -0.0*(exp(-(xp(l,3)-qz(nz/2))**2/(10.*delz)**2)
     x        *exp(-(xp(l,1)-qx(nx/2))**2/(10.*dx)**2))+vx
         vp(l,2) = vy 
         vp(l,3) = vz 

c         if (xp(l,3) .gt. qz(nz/2)) mix_ind(l) = 1
c         if (xp(l,3) .le. qz(nz/2)) mix_ind(l) = 0
            
         do 45 m=1,3
            vp1(l,m) = vp(l,m)
            input_E = input_E + 
     x           0.5*m_arr(l)*(vp(l,m)*km_to_m)**2 /beta
            input_p(m) = input_p(m) + m_arr(l)*vp(l,m) / beta
 45      continue
 10      continue

      do 50 l = Ni_tot_1+1,Ni_tot_1+0.06*Ni_tot_1

         xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))
         xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))
         xp(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz-1)-qz(1))
         m_arr(l) = mion
         mrat(l) = 1.0

         ijkp(l,1) = floor(xp(l,1)/dx) 

         i=0
 81      continue
         i = i + 1
         if (xp(l,1) .gt. qx(i)) go to 81 !find i on non-uniform 
         i = i-1
         ijkp(l,1)= i


         ijkp(l,2) = floor(xp(l,2)/dy) 
         
         k=0
 80      continue
         k = k + 1
         if (xp(l,3) .gt. qz(k)) go to 80 !find k on non-uniform 
         k = k-1
         ijkp(l,3)= k
         
c         vth = vth_bottom

         vx = vsw+1200.*sqrt(-alog(pad_ranf()))*cos(PI*pad_ranf())
         vy = 1200.*sqrt(-alog(pad_ranf()))*cos(PI*pad_ranf())
         vz = 500.*sqrt(-alog(pad_ranf()))*cos(PI*pad_ranf())

         ii = ijkp(l,1)
         kk = ijkp(l,3)
         vp(l,1) = -0.0*(exp(-(xp(l,3)-qz(nz/2))**2/(10.*delz)**2)
     x        *exp(-(xp(l,1)-qx(nx/2))**2/(10.*dx)**2))+vx
         vp(l,2) = vy 
         vp(l,3) = vz 

c         if (xp(l,3) .gt. qz(nz/2)) mix_ind(l) = 1
c         if (xp(l,3) .le. qz(nz/2)) mix_ind(l) = 0
            
         do 95 m=1,3
            vp1(l,m) = vp(l,m)
            input_E = input_E + 
     x           0.5*m_arr(l)*(vp(l,m)*km_to_m)**2 /beta
            input_p(m) = input_p(m) + m_arr(l)*vp(l,m) / beta
 95      continue
 50      continue

         Ni_tot = Ni_tot_1 + 0.06*Ni_tot_1


      do 20 l = Ni_tot+1,Ni_tot+0.06*Ni_tot

         xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))
         xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))
         xp(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz-1)-qz(1))
         m_arr(l) = 2*mion
         mrat(l) = 0.5

         ijkp(l,1) = floor(xp(l,1)/dx) 

         i=0
 28      continue
         i = i + 1
         if (xp(l,1) .gt. qx(i)) go to 28 !find i on non-uniform 
         i = i-1
         ijkp(l,1)= i


         ijkp(l,2) = floor(xp(l,2)/dy) 
         
         k=0
 25      continue
         k = k + 1
         if (xp(l,3) .gt. qz(k)) go to 25 !find k on non-uniform 
         k = k-1
         ijkp(l,3)= k
         
c         vth = vth_bottom

         vx = vsw+vth*sqrt(-alog(pad_ranf()))*cos(PI*pad_ranf())
         vy = vth*sqrt(-alog(pad_ranf()))*cos(PI*pad_ranf())
         vz = vth*sqrt(-alog(pad_ranf()))*cos(PI*pad_ranf())

         ii = ijkp(l,1)
         kk = ijkp(l,3)
         vp(l,1) = -0.0*(exp(-(xp(l,3)-qz(nz/2))**2/(10.*delz)**2)
     x        *exp(-(xp(l,1)-qx(nx/2))**2/(10.*dx)**2))+vx
         vp(l,2) = vy 
         vp(l,3) = vz 

c         if (xp(l,3) .gt. qz(nz/2)) mix_ind(l) = 1
c         if (xp(l,3) .le. qz(nz/2)) mix_ind(l) = 0
            
         do 35 m=1,3
            vp1(l,m) = vp(l,m)
            input_E = input_E + 
     x           0.5*m_arr(l)*(vp(l,m)*km_to_m)**2 /beta
            input_p(m) = input_p(m) + m_arr(l)*vp(l,m) / beta
 35      continue
 20      continue

         Ni_tot = Ni_tot + 0.06*Ni_tot
c         write(*,*) 'Ni_tot...',Ni_tot,Ni_tot_1
c         stop

      call get_interp_weights(xp)
      call update_np(np)
      call update_up(vp,np,up)

      return
      end SUBROUTINE load_aniso_Maxwellian
c----------------------------------------------------------------------





      end MODULE part_init









