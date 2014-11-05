      MODULE misc

      USE global
      USE gutsf
      USE boundary

      contains

c----------------------------------------------------------------------
      real FUNCTION pad_ranf()
c This is the random number generator that works on foo.
c----------------------------------------------------------------------
c      include 'incurv.h'

c      integer function irand
c      integer iflag, irnum
c      external irand

c      irnum = irand(iflag)

c      ranf = irnum/32767.0
c      ranf = irnum/2147483647.

      call random_number(pad_ranf)

      return
      end function pad_ranf
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      subroutine random_initialize ( seed_input )
c----------------------------------------------------------------------
!**********************************************************************
!
!! RANDOM_INITIALIZE initializes the FORTRAN90 random number seed.
!
!  Discussion:
!    If you don't initialize the FORTRAN90 random number generator
!    routine RANDOM_NUMBER, which is used by calls like
!
!      call random_number ( harvest = x )
!
!    then its behavior is not specified.  That may be OK for you.  But
!    you may want to be able to force the same sequence to be generated
!    each time, or to force a different sequence.
!
!    To control the sequence of random numbers, you need to set the seed.
!    In FORTRAN90, this is done by calling the RANDOM+SEED routine.
!    You can call it with no arguments, in fact.  But if you call
!    it with no arguments:
!
!      call random_seed ( )
!
!    then its behavior (or more particularly, the behavior of RANDOM_NUMBER)
!    is still not specified.  You might hope that the system will go out
!    and pick a nice random seed for you, but there's no guarantee.
!
!
!    For example, on the DEC ALPHA, if you compile a program that calls
!    RANDOM_NUMBER, then every time you run it, you get the same sequence
!    of "random" values.  If you compile a program that calls RANDOM_SEED
!    with no arguments, and then calls RANDOM_NUMBER, you still get the
!    same sequence each time you run the program.
!
!    In order to actually try to scramble up the random number generator
!    a bit, this routine goes through the tedious process of getting the
!    size of the random number seed, making up values based on the current
!    time, and setting the random number seed.
!
!    Unfortunately, the RANDOM_SEED routine has a very elastic definition.
!    It does not use a single scalar integer SEED.  Instead, it communicates
!    with the user through an integer vector whose size is not specified.
!    You actually have to "ask" the routine to tell you the size of this
!    vector.  Then you can fill up the vector with values to be used to
!    influence the seeding of the random number routine.  The details of
!    how the seed affects the sequence are also unspecified, but the only
!    thing we are reasonably confident about is that specifying the same
!    seed should result in the same sequence, and specifying different
!    seeds should result in different sequences!
!
!    I assume this is far more than you wanted to know.  (It's certainly
!    more than I wanted to know!)
!
!    The argument SEED is an input quantity only, so it is legal
!    to type
!
!      call random_initialize ( 0 )
!
!    or
!
!      call random_initialize ( 18867 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!     Input, integer ( kind = 4 ) SEED_INPUT, the user's "suggestion" for a seed
!     However, if the input value is 0, the routine will come up with
!     its own "suggestion", based on the system clock.
!     
c      implicit none
      
      integer count
      integer  count_max
      integer  count_rate
      logical, parameter :: debug = .false.
      integer  i
      integer  seed
      integer  seed_input
      integer, allocatable :: seed_vector(:)
      integer seed_size
      real    t
      integer, parameter :: warm_up = 100
      
      seed = seed_input
!     
!     Initialize the random seed routine.
!     
      call random_seed ( )
!     
!     Determine the size of the random number seed vector.
!     
      call random_seed ( size = seed_size )
!     
!     Allocate a vector of the right size to be used as a random seed.
!     
      allocate ( seed_vector(seed_size) )
!     
!     If the user supplied a SEED value, use that.
!     
!     Otherwise, use the system clock value to make up a value that is
!     likely to change based on when this routine is called.
!     
      if ( seed /= 0 ) then
         
         if ( debug ) then
            write (*,*) ' '
            write (*,*) 'RANDOM_INITIALIZE'
            write (*,*) 'Initialize RANDOM_NUMBER, 
     x                                 user SEED = ', seed
         end if
         
      else
         
         call system_clock ( count, count_rate, count_max )
         
      seed = count
      
      if ( debug ) then
         write (*,*) ' '
         write (*,* ) 'RANDOM_INITIALIZE'
         write (*,* ) '  Initialize RANDOM_NUMBER, arbitrary SEED = ',
     x                 seed
      end if
      
      end if
!     
!     Set the seed vector.  We don't know the significance of the
!     individual entries of the internal seed vector, so we'll just set
!     all entries to SEED.
!     
      seed_vector(1:seed_size) = seed
!     
!     Now call RANDOM_SEED, and tell it to use this seed vector.
!     
      call random_seed ( put = seed_vector(1:seed_size) )
!     
!     Free up the seed space.
!     
      deallocate ( seed_vector )
!     
!     Call the random number routine a bunch of times just to "warm it up".
!     
      do i = 1, warm_up
         call random_number ( harvest = t )
      enddo

      return
      end subroutine random_initialize
c----------------------------------------------------------------------





c----------------------------------------------------------------------
      SUBROUTINE Momentum_diag(up,uf,np,nf,E,b1,pup,puf,peb,input_p)
c----------------------------------------------------------------------
c      include 'incurv.h'

      real up(nx,ny,nz,3),
     x     uf(nx,ny,nz,3),
     x     np(nx,ny,nz),
     x     nf(nx,ny,nz),
     x     E(nx,ny,nz,3),
     x     b1(nx,ny,nz,3),
     x     pup(3),
     x     puf(3),
     x     peb(3),
     x     input_p(3)

      real vol
      real mom_flux
      real exb(nx,ny,nz,3)
      real npave(3),nfave(3)


      call crossf2(E,b1,exb)

      do 5 m=1,3
         pup(m) = 0
         puf(m) = 0
         peb(m) = 0
 5              continue

      do 10 i=2,nx-1
         do 10 j=2,ny-1
            do 10 k=2,nz-1
               ip = i+1
               jp = j+1
               kp = k+1
               if (ip .eq. nx) then ip = nx-1
               if (jp .eq. ny) then jp = ny-1
               if (kp .eq. nz) then kp = nz-1
               vol = dx_cell(i)*dy_cell(j)*dz_cell(k)
               npave(1) = 0.5*(np(i,j,k) + np(ip,j,k))
               npave(2) = 0.5*(np(i,j,k) + np(i,jp,k))
               npave(3) = 0.5*(np(i,j,k) + np(i,j,kp))
               nfave(1) = 0.5*(nf(i,j,k) + nf(ip,j,k))
               nfave(2) = 0.5*(nf(i,j,k) + nf(i,jp,k))
               nfave(3) = 0.5*(nf(i,j,k) + nf(i,j,kp))
               do 10 m=1,3
c                  pup(m) = pup(m) + npave(m)*vol*mBa*up(i,j,k,m)
                  pup(m) = pup(m) + np(i,j,k)*vol*mBa*up(i,j,k,m)
c                  puf(m) = puf(m) + nfave(m)*vol*mO*uf(i,j,k,m)
                  puf(m) = puf(m) + nf(i,j,k)*vol*mO*uf(i,j,k,m)
                  peb(m) = peb(m) + epsilon*1e3*exb(i,j,k,m)*vol*(mO/q)
 10               continue

c      write(*,*) 'Momentum conservation...'
c      write(*,*) '  Particles.............',pup(1),pup(2),pup(3)
c      write(*,*) '  Fluid.................',puf(1),puf(2),puf(3)
c      write(*,*) '  ExB...................',peb(1),peb(2),peb(3)
c      write(*,*) '  Normalized............',
c     x                     (pup(1)+puf(1)+peb(1))/input_p(1),
c     x                     (pup(2)+puf(2)+peb(2))/input_p(2),
c     x                     (pup(3)+puf(3)+peb(3))/input_p(3)

c Momentum flux through boundary faces

c i = 2 face

c      do 20 j=2,ny
c         do 20 k=2,nz
c            m=1
c            i=2
c            vol = uf(i,j,k,m)*dtsub*dy*dz_cell(k)
c            mom_flux = nf(i,j,k)*vol*mO*uf(i,j,k,m)
c            input_p(m) = input_p(m) + mom_flux  !+ sign since pos 
                                                !is flux into domain
c 20         continue

c i = nx face

c      do 30 j=2,ny
c         do 30 k=2,nz
c            m=1
c            i=nx
c            vol = uf(i,j,k,m)*dtsub*dy*dz_cell(k)
c            mom_flux = nf(i,j,k)*vol*mO*uf(i,j,k,m)
c            input_p(m) = input_p(m) - mom_flux  !- sign since pos 
                                                !is flux out domain
c 30         continue
c**********************
c j = 2 face

c      do 40 i=2,nx
c         do 40 k=2,nz
c            m=2
c            j=2
c            vol = uf(i,j,k,m)*dtsub*dy*dz_cell(k)
c            mom_flux = nf(i,j,k)*vol*mO*uf(i,j,k,m)
c            input_p(m) = input_p(m) + mom_flux  !+ sign since pos 
                                                !is flux into domain
c 40         continue

c j = ny face

c      do 50 i=2,nx
c         do 50 k=2,nz
c            m=2
c            j=ny
c            vol = uf(i,j,k,m)*dtsub*dy*dz_cell(k)
c            mom_flux = nf(i,j,k)*vol*mO*uf(i,j,k,m)
c            input_p(m) = input_p(m) - mom_flux  !- sign since pos 
                                                !is flux out domain
c 50         continue
c****************
c k = 2 face

c      do 60 i=2,nx
c         do 60 j=2,ny
c            m=3
c            k=2
c            vol = uf(i,j,k,m)*dtsub*dy*dz_cell(k)
c            mom_flux = nf(i,j,k)*vol*mO*uf(i,j,k,m)
c            input_p(m) = input_p(m) + mom_flux  !+ sign since pos 
                                                !is flux into domain
c 60         continue

c k = nz face

c      do 70 i=2,nx
c         do 70 j=2,ny
c            m=3
c            k=nz
c            vol = uf(i,j,k,m)*dtsub*dy*dz_cell(k)
c            mom_flux = nf(i,j,k)*vol*mO*uf(i,j,k,m)
c            input_p(m) = input_p(m) - mom_flux  !- sign since pos 
                                                !is flux out domain
c 70         continue

c      write(*,*) 'Normalized x momentum...',(pup(1)+puf(1))/input_p(1)
c      write(*,*) 'Normalized y momentum...',(pup(2)+puf(2))/input_p(2)
c      write(*,*) 'Normalized z momentum...',(pup(3)+puf(3))/input_p(3)

      return
      end SUBROUTINE Momentum_diag
c----------------------------------------------------------------------



c----------------------------------------------------------------------
      SUBROUTINE get_bndry_Eflux(b1,E)
c----------------------------------------------------------------------
c      include 'incurv.h'
     
      real b1(nx,ny,nz,3),
     x     E(nx,ny,nz,3)

      real vol
      real uf_flux
      real exb_flux

      real mO_q
      
      mO_q = mO/q


c Energy flux through boundary faces

cc i = 2 face 

c      do 20 j=2,ny
c         do 20 k=2,nz
c            m=1
c            i=2
c            vol = uf(i,j,k,m)*dtsub*dy*dz_cell(k)
c            uf_flux = 0.5*nf(i,j,k)*vol*mO*(uf(i,j,k,m)*km_to_m)**2
c            exb_flux = (mO_q)**2*(1.0/mu0)*dtsub*dy*dz_cell(k)*
c     x           (E(i,j,k,2)*b1(i,j,k,3) - E(i,j,k,3)*b1(i,j,k,2))*
c     x            km_to_m**2
c            bndry_Eflux = bndry_Eflux + uf_flux + exb_flux 
c                               !+ sign since pos is flux into domain
c 20         continue

cc i = nx face

c      do 30 j=2,ny
c         do 30 k=2,nz
c            m=1
c            i=nx
c            vol = uf(i,j,k,m)*dtsub*dy*dz_cell(k)
c            uf_flux = 0.5*nf(i,j,k)*vol*mO*(uf(i,j,k,m)*km_to_m)**2
c            exb_flux = (mO_q)**2*(1.0/mu0)*dtsub*dy*dz_cell(k)*
c     x           (E(i,j,k,2)*b1(i,j,k,3) - E(i,j,k,3)*b1(i,j,k,2))*
c     x           km_to_m**2
c            bndry_Eflux = bndry_Eflux - uf_flux - exb_flux 
c                               !- sign since neg is flux into domain
c 30         continue
cc**********************
cc j = 2 face

c      do 40 i=2,nx
c         do 40 k=2,nz
c            m=2
c            j=2
c            vol = uf(i,j,k,m)*dtsub*dx*dz_cell(k)
c            uf_flux = 0.5*nf(i,j,k)*vol*mO*(uf(i,j,k,m)*km_to_m)**2
c            exb_flux = (mO_q)**2*(1.0/mu0)*dtsub*dx*dz_cell(k)*
c     x            (-E(i,j,k,1)*b1(i,j,k,3)+E(i,j,k,3)*b1(i,j,k,1))*
c     x            km_to_m**2
c            bndry_Eflux = bndry_Eflux + uf_flux + exb_flux 
c                               !+ sign since neg is flux into domain
c 40         continue

cc j = ny face

c      do 50 i=2,nx
c         do 50 k=2,nz
c            m=2
c            j=ny
c            vol = uf(i,j,k,m)*dtsub*dx*dz_cell(k)
c            uf_flux = 0.5*nf(i,j,k)*vol*mO*(uf(i,j,k,m)*km_to_m)**2
c            exb_flux = (mO_q)**2*(1.0/mu0)*dtsub*dx*dz_cell(k)*
c     x            (-E(i,j,k,1)*b1(i,j,k,3)+E(i,j,k,3)*b1(i,j,k,1))*
c     x            km_to_m**2
c            bndry_Eflux = bndry_Eflux - uf_flux - exb_flux 
c                               !- sign since neg is flux into domain
c 50         continue
cc****************
c k = 2 face

      do 60 i=2,nx
         do 60 j=2,ny
            m=3
c            k=rk-20
            k=2
c            vol = uf(i,j,k,m)*dtsub*dx*dy
c            uf_flux = 0.5*nf(i,j,k)*vol*mO*(uf(i,j,k,m)*km_to_m)**2
            exb_flux = (mO_q)**2*(1.0/mu0)*dt*dx*dy*
     x             (E(i,j,k,1)*b1(i,j,k,2) - E(i,j,k,2)*b1(i,j,k,1))*
     x             km_to_m**3
            bndry_Eflux = bndry_Eflux + exb_flux 
                               !- sign since neg is flux into domain
 60         continue

c k = nz face

      do 70 i=2,nx
         do 70 j=2,ny
c            m=3
c            k=rk+20
            k=nz-1
c            vol = uf(i,j,k,m)*dtsub*dx*dy
c            uf_flux = 0.5*nf(i,j,k)*vol*mO*(uf(i,j,k,m)*km_to_m)**2
            exb_flux = (mO_q)**2*(1.0/mu0)*dt*dx*dy*
     x             (E(i,j,k,1)*b1(i,j,k,2) - E(i,j,k,2)*b1(i,j,k,1))*
     x             km_to_m**3
            bndry_Eflux = bndry_Eflux - exb_flux 
                               !- sign since neg is flux into domain
 70         continue

      return
      end SUBROUTINE get_bndry_Eflux
c----------------------------------------------------------------------



c----------------------------------------------------------------------
      SUBROUTINE get_beta()
c----------------------------------------------------------------------
c      include 'incurv.h'

      real vol

c divide particles up between procnum processors      
      vol = ((qx(nx-1)-qx(1))*(qy(ny-1)-qy(1))*(qz(nz)-qz(1)))/2
      beta = (Ni_tot_sys/vol)/np_top

      write(*,*) 'beta...',beta

      return
      end SUBROUTINE get_beta
cc----------------------------------------------------------------------


c----------------------------------------------------------------------
      subroutine get_np3(np,np3)
c----------------------------------------------------------------------
c      include 'incurv.h'

      real np(nx,ny,nz)
c      real nf(nx,ny,nz)
      real np3(nx,ny,nz,3)

      real nfp(nx,ny,nz)

      nfp = np

      do i = 2,nx-1
         do j = 2,ny-1
            do k = 2,nz-1
               np3(i,j,k,1) = 0.5*(nfp(i,j,k)+nfp(i+1,j,k))
               np3(i,j,k,2) = 0.5*(nfp(i,j,k)+nfp(i,j+1,k))
               np3(i,j,k,3) = 0.5*(nfp(i,j,k)+nfp(i,j,k+1))
            enddo
         enddo
      enddo

      call periodic(np3)

      return
      end subroutine get_np3
c----------------------------------------------------------------------


c----------------------------------------------------------------------
c      SUBROUTINE charge_exchange(np,xp,vp,vp1,m,chex_rate,
c     x                           input_p)
c----------------------------------------------------------------------
c      include 'incurv.h'

c      real np(nx,ny,nz)
c      real xp(Ni_max,3)
c      real vp(Ni_max,3)
c      real vp1(Ni_max,3)
c      real vn(3)
c      real input_p(3)

c      real cx,cy,cz          !neutral cloud center
c      real rx,ry,rz,r  
c      real t            !run time
c      real vr           !relative velocity between ions and neutrals
c      real sigma_chex
c      parameter (sigma_chex = 1.0e-24)   !km^2  check this
c      real vol
c      real dNcx
c      integer*4 nchex
c      real rnd
c      real nn           !neutral density
c      real nconst
c      real initial_E

c      nconst = vth*sqrt(pi)

c      call Neut_Center(cx,cy,cz)
      
c      initial_E = input_E
c      chex_rate = 0.0
c      nchex = 0      
c      do 10 l=1,Ni_tot

c         i=ijkp(l,1)
c         j=ijkp(l,2)
c         k=ijkp(l,3)

c         rx = qx(i) - cx
c         ry = qy(j) - cy
c         rz = qz(k) - cz
c         r = sqrt(rx**2 + ry**2 + rz**2)
c         vn(1) = vsat + rx/t
c         vn(2) = ry/t
c         vn(3) = rz/t
c         vr = sqrt((vp(l,1) - vn(1))**2 + 
c     x             (vp(l,2) - vn(2))**2 +
c     x             (vp(l,3) - vn(3))**2)

c         if (r .gt. 2.33*t) then    !2.33 km/s as distbn limit
                                    !otherwise float underflow
c            nn = 0.0
c         else
c            nn = (No/(4*pi*r*r*t*nconst)) *
c     x                exp(-(r-vo*t)**2 / (vth*t)**2)
c         endif

c            dNcx = dt*vr*sigma_chex*nn
cc            write(*,*) 'dNcx....',dNcx,dNcx/beta
c            rnd = pad_ranf()
c            if (rnd .lt. dNcx) then
c               nchex = nchex + 1
c               vol = dx*dy*dz_grid(k)
c               do 30 m=1,3             !remove neutral energy
c                  vp1(l,m) = vp(l,m)
c                  input_E = input_E -  
c     x                      0.5*mBa*(vp(l,m)*km_to_m)**2 /beta
c                  input_p(m) = input_p(m) - mBa*vp(l,m) / beta
c 30            continue
c               vp(l,1) = vn(1)
c               vp(l,2) = vn(2)
c               vp(l,3) = vn(3)
cc               xp(l,1) = xp(l,1)
cc               xp(l,2) = xp(l,2)
cc               xp(l,3) = xp(l,3)
c               do 40 m=1,3             !add ion energy
c                  vp1(l,m) = vp(l,m)
c                  input_E = input_E +  
c     x                      0.5*mBa*(vp(l,m)*km_to_m)**2 /beta
c                  input_p(m) = input_p(m) + mBa*vp(l,m) / beta
c 40               continue
c               endif
         
c 10      continue 


cc      write(*,*) 'nchex,chex_rate...',real(nchex)/beta,
cc     x            chex_rate/beta

c      input_chex = input_chex + (input_E - initial_E)
c      chex_rate = (real(nchex))/(dt*beta)
c      write(*,*) 'Normalized charge exchange energy gain...',
c     x            input_chex/(input_E - input_chex - input_bill),
c     x            input_E,input_chex,input_bill 
c      write(*,*) 'Charge exchange rate...',chex_rate


c      return
c      end SUBROUTINE charge_exchange
cc----------------------------------------------------------------------




      end module misc


