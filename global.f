      MODULE GLOBAL
      
      USE inputs
      USE dimensions
c      include 'mpif.h'

c contains simulation parameter list for a barium release

c simulation domain dimensions
c      PARAMETER (nx = 3, ny = 3, nz = 801)

c magnetic field and mass info for determining time step

c raw grid coordinate data
  
      real qx(nx), qy(ny), qz(nz)
                            !cartisian coordinates of grid points
      real lambda(nz)       !latitude of given grid point
      integer ri, rj, rk    !release coordinates
      real dz_grid(nz)      !dz for grid, variable in z direction
      real dz_cell(nz)      !dz for cell, variable in z direction 
      real dx_grid(nx),dx_cell(nx)
      real dy_grid(ny),dy_cell(ny)
      real xrat(nx), yrat(ny), zrat(nz)

      common /coords/ qx,qy,qz,lambda,ri,rj,rk,dz_grid,dz_cell,zrat,
     x                xrat,yrat

c Total number of ions produced at a given time during the simulation
      integer*4 Ni_tot, Nn_tot, Ni_tot_sys, Ni_tot_sw       

c Location (indices) of particles with the grid
      integer ijkp(Ni_max,3)
      logical in_bounds(Ni_max)
c      logical ionized(Ni_max)
      real mix_ind(Ni_max)
      common /ion_info/ Ni_tot, Nn_tot, ijkp, Ni_tot_sys, Ni_tot_sw,
     x       in_bounds,mix_ind

      integer*4 seed
      common /rndseed/ seed

c Indices for boundary cells
      integer ip,im,jp,jm,kp,km

c Total input energy (ions) 
      real input_E,input_Eb,input_EeP,prev_Etot, bndry_Eflux
      real input_chex, input_bill
      common /energy/ input_E, input_EeP, prev_Etot, bndry_Eflux, 
     x                input_chex, input_bill,input_Eb

c Momentum conservation
      real surf_tot(3),graduu_tot(3),ugradu_tot(3)
      common /momentum/ surf_tot,graduu_tot,ugradu_tot

c Dummy variables for doing vector operations and particle flux
c calculation....etc.

      real a(nx,ny,nz,3),c(nx,ny,nz,3),ct(nx,ny,nz,3)
      common /dummy/ a,c,ct

c Variable time step for fluid velocity update used in both
c particle update and fluid update

      real dtf,ntf,dtsub
      common /timestep/ ntf,dtsub

c Weight variables for trilinear interpolation

      real wght(Ni_max,8)
      common /weights/ wght

c      real wght(Ni_max,8), wquad(Ni_max,3)
c      common /weights/ wght, wquad


c variable for anomlous resistivity

c      real etar(nx,ny,nz,3)
      real nuei
      common /resist/ etar,nuei

c variable for particle scale

      real beta, dNi, dNi_sw
      common /scale/ beta, dNi, dNi_sw

c mass array for multi-ion species
      real mrat,m_arr
      common /mass/ mrat(Ni_max),m_arr(Ni_max)

c mixing array
      real mixed(nx,ny,nz)
      common /mix/ mixed


c parallel processor info
      integer procnum,my_rank
      common /procinfo/ procnum, my_rank

      
      END MODULE GLOBAL
