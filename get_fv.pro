;----------------------------------------------------------------
PRO read_part,file,nfrm,xp
;----------------------------------------------------------------

Ni_max=long(0)
nt=0l
ntout=0l
frm=0l

file = file+'.dat'
print,' reading...',file
openr,1,file,/f77_unformatted
readu,1,nt
readu,1,ntout
readu,1,Ni_max
print,nt,ntout,Ni_max

xp=fltarr(Ni_max,3,/nozero)

readu,1,frm
print,'  image #.....',frm
readu,1,xp
frmcnt = 1

while (frmcnt lt nfrm) do begin

   readu,1,frm
   print,'  image #.....',frm
   readu,1,xp
   frmcnt = frmcnt + 1

endwhile

close,1

return
end

;----------------------------------------------------------------


;----------------------------------------------------------------
PRO read_part_scalar,file,nfrm,xp
;----------------------------------------------------------------

Ni_max=long(0)
nt=0
ntout=0
frm=0

file = file+'.dat'
print,' reading...',file
openr,1,file,/f77_unformatted
readu,1,nt
readu,1,ntout
readu,1,Ni_max
print,nt,ntout,Ni_max

xp=fltarr(Ni_max,/nozero)

readu,1,frm
print,'  image #.....',frm
readu,1,xp
frmcnt = 1

while (frmcnt lt nfrm) do begin

   readu,1,frm
   print,'  image #.....',frm
   readu,1,xp
   frmcnt = frmcnt + 1

endwhile

close,1

return
end
;----------------------------------------------------------------


;----------------------------------------------------------------
pro get_dist,xyz_traj,lxyz,nlxyz,t
;----------------------------------------------------------------


mp=1.67e-27
nfile = '1'
nfrm = 1
;define look direction of instrument
;theta = 180.*!dtor
;phi = !pi/4
;dtheta= 20.*!dtor
;dphi = 20.*!dtor

rundir = '/Volumes/MacD97-2/hybrid/SWAP/run_test/'

part_dir = rundir+'part_1/'
fluid_fields_dir = rundir

print,part_dir

f_read_coord,fluid_fields_dir+'coord.dat',x,y,z,dzc,dzg,nx,ny,nz
f_read_3d_m,fluid_fields_dir+'npall_'+strtrim(string(nfile),2),nfrm,np
f_read_3d_vec_m_32,fluid_fields_dir+'b1all_'+strtrim(string(nfile),2),nfrm,b1
f_read_3d_vec_m_32,fluid_fields_dir+'upall_'+strtrim(string(nfile),2),nfrm,up

ztemp = fltarr(35)

dv = 2.0
nn = 4000.
alpha=1.9263418e-20
d3v = dv^3
vxyp = fltarr(nn/dv,nn/dv)
vxzp = fltarr(nn/dv,nn/dv)
vyzp = fltarr(nn/dv,nn/dv)
;fvp = fltarr(nn/dv,nn/dv,nn/dv)
;vxyi = fltarr(nn/dv,nn/dv)
;vxzi = fltarr(nn/dv,nn/dv)
;vyzi = fltarr(nn/dv,nn/dv)
;fvi = fltarr(nn/dv,nn/dv,nn/dv)
vx = -(nn/(2*dv))*dv + findgen(nn/dv)*dv
vy = -(nn/(2*dv))*dv + findgen(nn/dv)*dv
vz = -(nn/(2*dv))*dv + findgen(nn/dv)*dv

xpluto = nx/2
ypluto = ny/2 
zpluto = nz/2-1

;xpluto = xyz_traj(0)
;ypluto = xyz_traj(1)
;zpluto = xyz_traj(2)

dx = 300.0
ndx = 10.0
;Rp = 1200.
;svol = (2.0*ndx*dx*1e5)^3
b = b1(xpluto,ypluto,zpluto,*)
bx = b1(xpluto,ypluto,zpluto,0)
by = b1(xpluto,ypluto,zpluto,1)
bz = b1(xpluto,ypluto,zpluto,2)

;bmag = sqrt(bx^2 + by^2 + bz^2)
;bxhat = bx/bmag
;byhat = by/bmag
;bzhat = bz/bmag

bxy = b(1)/b(0)
bzx = b(2)/b(0)
bzy = b(2)/b(1)

;;look direction xyz
;lxhat = lxyz(0)
;lyhat = lxyz(1)
;lzhat = lxyz(2)

;lmag = sqrt(lxhat^2 + lyhat^2 + lzhat^2)
;lxhat = lxhat/lmag
;lyhat = lyhat/lmag
;lzhat = lzhat/lmag

;nlxhat = nlxyz(0)  ;axis perpedicular to the look direction
;nlyhat = nlxyz(1)
;nlzhat = nlxyz(2)

;look direction perp b
;nlxhat = bxhat
;nlyhat = byhat
;nlzhat = bzhat

;lxhat = 1.0
;lyhat = -lxhat*bxhat/byhat
;lzhat = 0.0

;lmag = sqrt(lxhat^2 + lyhat^2 + lzhat^2)
;lxhat = lxhat/lmag
;lyhat = lyhat/lmag
;lzhat = lzhat/lmag


for nfil = 0,5 do begin

    xfile = part_dir+'xp_'+strtrim(string(nfil),2)
    vfile = part_dir+'vp_'+strtrim(string(nfil),2)
    mratfile = part_dir+'mrat_'+strtrim(string(nfil),2)

    read_part,xfile,nfrm,xp
    read_part,vfile,nfrm,vp
    read_part_scalar,mratfile,nfrm,mrat

    i = xpluto
    j = ypluto
    k = zpluto

    wh = where((xp(*,0) ge x(i)-ndx*dx) and (xp(*,0) le x(i)+ndx*dx) and $
               (xp(*,1) ge y(j)-ndx*dx) and (xp(*,1) le y(j)+ndx*dx) and $ 
               (xp(*,2) ge z(k)-ndx*dx) and (xp(*,2) le z(k)+ndx*dx) and $
               (mrat(*) eq 1.0/1.0))

    if (wh(0) gt -1) then begin
    for l = 0ll,n_elements(wh)-1 do begin
       ii = fix(vp(wh(l),0)/dv) + (nn/(2*dv))
       jj = fix(vp(wh(l),1)/dv) + (nn/(2*dv))
       kk = fix(vp(wh(l),2)/dv) + (nn/(2*dv))
       vxyp(ii,jj) = 1.0 + vxyp(ii,jj)
       vxzp(ii,kk) = 1.0 + vxzp(ii,kk)
       vyzp(jj,kk) = 1.0 + vyzp(jj,kk)
;       fvp(ii,jj,kk) = 1.0 + fvp(ii,jj,kk)
    endfor
    endif

;    wh = where((xp(*,0) ge x(i)-ndx*dx) and (xp(*,0) le x(i)+ndx*dx) and $
;               (xp(*,1) ge y(j)-ndx*dx) and (xp(*,1) le y(j)+ndx*dx) and $ 
;               (xp(*,2) ge z(k)-ndx*dx) and (xp(*,2) le z(k)+ndx*dx) and $
;               (mrat(*) lt 1.0))
;    if (wh(0) gt -1) then begin
;    for l = 0ll,n_elements(wh)-1 do begin
;       ii = fix(vp(wh(l),0)/dv) + (nn/(2*dv))
;       jj = fix(vp(wh(l),1)/dv) + (nn/(2*dv))
;       kk = fix(vp(wh(l),2)/dv) + (nn/(2*dv))
;       vxyi(ii,jj) = 1.0 + vxyi(ii,jj)
;       vxzi(ii,kk) = 1.0 + vxzi(ii,kk)
;       vyzi(jj,kk) = 1.0 + vyzi(jj,kk)
;       fvi(ii,jj,kk) = 1.0 + fvi(ii,jj,kk)
;    endfor
;    endif
endfor

xx = reverse((x - x(55)))
xx = ((x - x(nx-55)))
yy = (y - y(ny/2))

;window,0,xsize=800,ysize=800
;contlevs = [0.02,0.05,0.1,0.2,0.5,1.0]*1e15
;contour,reform(np(*,*,zpluto)),levels=contlevs,xx/1e3,yy/1e3,$
;    /isotropic,xrange=[-140,40],yrange=[-70,70],xtitle='x 10!u3!n km',$
;    ytitle='y 10!u3!n km',/xsty,/ysty
;plots,[xx(xpluto)-ndx*dx,xx(xpluto)+ndx*dx]/1e3,([y(ypluto)-ndx*dx,y(ypluto)-ndx*dx]-y(ny/2))/1e3,thick=3
;plots,[xx(xpluto)-ndx*dx,xx(xpluto)+ndx*dx]/1e3,([y(ypluto)+ndx*dx,y(ypluto)+ndx*dx]-y(ny/2))/1e3,thick=3
;plots,[xx(xpluto)-ndx*dx,xx(xpluto)-ndx*dx]/1e3,([y(ypluto)-ndx*dx,y(ypluto)+ndx*dx]-y(ny/2))/1e3,thick=3
;plots,[xx(xpluto)+ndx*dx,xx(xpluto)+ndx*dx]/1e3,([y(ypluto)-ndx*dx,y(ypluto)+ndx*dx]-y(ny/2))/1e3,thick=3

;t1 = 25e3/cos(15*!dtor) - atan(15*!dtor)*x
t2 = -13.5e3/cos(15*!dtor) - atan(15*!dtor)*xx
;t3 = 7.5e3/cos(15*!dtor) - atan(15*!dtor)*x
 
;oplot,x/1e3,t1/1e3,linestyle=1,color=255,thick=4
oplot,xx/1e3,t2/1e3,linestyle=2,thick=1
;oplot,x/1e3,t3/1e3,linestyle=1,color=255,thick=4
;oplot,t.x/1e3,t.y/1e3,linestyle=3,thick=1

xrng = 900

device,decompose=0
loadct,39
;szf = size(fvp)
;fvp = smooth(fvp,2)
sz = size(vxyp)
lev=[1/exp(6),1/exp(5),1/exp(4),1/exp(3),1/exp(2),1/exp(1),1.0]*1.0*max(vxyp)

contour,smooth(vxyp,2),vx,vy,/isotropic,levels=lev,$
  xtitle='vx (km/s)',ytitle='vy (km/s)',$
  xrange=[-xrng,xrng],yrange=[-xrng,xrng],/xsty,/ysty

contour,smooth(vxzp,2),vx,vz,/isotropic,levels=lev,$
 xtitle='vx (km/s)',ytitle='vz (km/s)',$
  xrange=[-xrng,xrng],yrange=[-xrng,xrng],/xsty,/ysty

;xrng=400
vxzp_sm = vxzp 
vxzp_sm = smooth(vxzp,2)
vyzp_sm = median(vyzp,2)

;lev=[1/exp(4),1/exp(3),1/exp(2),1/exp(1),1.0]*1.1*max(vyzp_sm)
contour,smooth(vyzp,2),vy,vz,/isotropic,levels=lev,$
  xtitle='vy (km/s)',ytitle='vz (km/s)',$
  xrange=[-xrng,xrng],yrange=[-xrng,xrng],/xsty,/ysty

;vxy = reform(fvp(*,*,szf(3)/2))
;vxy_sm = smooth(vxy,10)

vxyp_sm = smooth(vxyp,20)
vxzp_sm = smooth(vxzp,20)

plot,vx,vxyp_sm(*,sz(2)/2),xrange=[-600,600],/xsty,xtitle='vx'
oplot,vx,1.0*max(vxyp_sm(*,sz(2)/2))*exp(-(vx-400)^2/1.^2),linestyle=1
plot,vz,vxzp_sm(sz(1)/2,*),xrange=[-600,600],/xsty,xtitle='vz'
oplot,vz,1.0*max(vxzp_sm(sz(1)/2,*))*exp(-vz^2/200.^2),linestyle=1
;contour,smooth(vyzi,2),vy,vz,/isotropic,nlev=10,xtitle='vy (km/s)',ytitle='vz (km/s)',$
;  xrange=[-xrng,xrng],yrange=[-xrng,xrng],/xsty,/ysty,/overplot,$
;  c_color=fsc_color('green')

;wh = where((abs(bzy*vy) le xrng) and (abs(vy) lt xrng))
;plots,vy(wh),bzy(wh)*vy(wh),linestyle=1
;arrow,0,0,lyhat*xrng/2,lzhat*xrng/2,color=fsc_color('red'),/data
;arrow,0,0,nlyhat*xrng/2,nlzhat*xrng/2,color=fsc_color('blue'),/data

;sz = size(fvp)
;scale3,xrange=[0,sz(1)],yrange=[0,sz(2)],zrange=[0,sz(3)]

;shade_volume, fvp,0.1*max(fvp),v,p  
;tv,polyshade(v,p,/t3d)           

;protons

wh = where(fvp gt 0.0)
ijkarr = array_indices(fvp,wh)
iarr = reform(ijkarr(0,*))
jarr = reform(ijkarr(1,*))
karr = reform(ijkarr(2,*))
vxarr = vx(iarr)
vyarr = vy(jarr)
vzarr = vz(karr)
varr = sqrt(vxarr^2 + vyarr^2 + vzarr^2)
vxhat = vxarr/varr
vyhat = vyarr/varr
vzhat = vzarr/varr

vhat = [vxhat,vyhat,vzhat]
lhat = [lxhat,lyhat,lzhat]
nlhat = [nlxhat,nlyhat,nlzhat]

lhat = [lxhat,lyhat,lzhat]
Earr = 0.0
cntarr = 0.0
flg = 'false'
for i = 0,n_elements(wh)-1 do begin
   vhat = [vxhat(i),vyhat(i),vzhat(i)]
   vdotl=transpose(vhat)#lhat
   vdotnl = transpose(vhat)#nlhat
;   print,vdotl
   if ((vdotl lt cos(45*!dtor)) and (abs(vdotnl) lt cos(80*!dtor))) then begin
      Earr = [Earr,0.5*mp*(vxarr(i)^2 + vyarr(i)^2 + vzarr(i)^2)*1e6/1.6e-19] 
;      print,Earr
      flg = 'true'
   endif
endfor

if (flg eq 'false') then begin
   print, 'No proton counts...'
   plot,[10,10000],[0,1],/nodata,title='Protons',/xlog
   goto, SKIP
endif
Earr = Earr(1:*)

dE = 10.
hmin = 1.0
hmax = 10000.

h = histogram(Earr,binsize=dE,min = hmin, max = hmax)
xE = dE*(findgen(n_elements(h)) + hmin)


;!p.multi=[0,1,1]
;window,1,xsize=400,ysize=400
plot,xE,h,psym=10,xtitle='Energy (eV)',ytitle='Counts',/xlog,$
   xrange=[10,10000],yrange=[0,max(h)*1.1],/ysty,/xsty,title='Protons'


SKIP:
;heavy ions

wh = where(fvi gt 0.0)
ijkarr = array_indices(fvi,wh)
iarr = reform(ijkarr(0,*))
jarr = reform(ijkarr(1,*))
karr = reform(ijkarr(2,*))
vxarr = vx(iarr)
vyarr = vy(jarr)
vzarr = vz(karr)
varr = sqrt(vxarr^2 + vyarr^2 + vzarr^2)
vxhat = vxarr/varr
vyhat = vyarr/varr
vzhat = vzarr/varr

vhat = [vxhat,vyhat,vzhat]
lhat = [lxhat,lyhat,lzhat]
nlhat = [nlxhat,nlyhat,nlzhat]

lhat = [lxhat,lyhat,lzhat]
Earr = 0.0
cntarr = 0.0
flg = 'false'
for i = 0,n_elements(wh)-1 do begin
   vhat = [vxhat(i),vyhat(i),vzhat(i)]
   vdotl=transpose(vhat)#lhat
   vdotnl = transpose(vhat)#nlhat
;   print,vdotl
   if ((vdotl lt 0.0) and (abs(vdotnl) lt cos(80*!dtor))) then begin
      Earr = [Earr,0.5*28.0*mp*(vxarr(i)^2 + vyarr(i)^2 + vzarr(i)^2)*1e6/1.6e-19] 
;      print,Earr
      flg = 'true'
   endif
endfor

if (flg eq 'false') then begin
   print, 'No pickup counts...'
   plot,[10,10000],[0,1],/nodata,title='Pickup ions',/xlog
   goto,SKIP_PICKUP
endif
Earr = Earr(1:*)

dE = 10
hmin = 1.0
hmax = 7500.

h = histogram(Earr,binsize=dE,min = hmin, max = hmax)

xE = dE*(findgen(n_elements(h))+hmin)

;!p.multi=[0,1,1]
;window,2,xsize=400,ysize=400
;wh1 = where(xE1 le 100)
;wh2 = where(xE2 gt 100)
plot,xE,h,psym=10,xtitle='Energy (eV)',ytitle='Counts',/xlog,$
   xrange=[10,10000],yrange=[0,max(h)*1.1],/ysty,/xsty,title='Pickup ions'


SKIP_PICKUP:

;device,/close

return
end
;----------------------------------------------------------------------


;----------------------------------------------------------------------
pro read_traj,t
;----------------------------------------------------------------------

t = {trajectory,Julian: 0.0, x: 0.0, y: 0.0, z: 0.0, vx: 0.0, vy: 0.0, $
     vz: 0.0,sx: 0.0, sy: 0.0, sz: 0.0, px: 0.0, py: 0.0, pz: 0.0} 

d = {trajectory,Julian: 0.0, x: 0.0, y: 0.0, z: 0.0, vx: 0.0, vy: 0.0, $
     vz: 0.0,sx: 0.0, sy: 0.0, sz: 0.0, px: 0.0, py: 0.0, pz: 0.0} 

close,1
openr, 1, 'trajectory1.csv'

junk = ' ' 
readf, 1, junk
print,junk


while not(eof(1)) do begin
   readf,1,d
   t = [t,d]
endwhile

t = t(1:*)

close,1

return
end
;----------------------------------------------------------------------


;----------------------------------------------------------------------
;main program
;----------------------------------------------------------------------

;set_plot,'ps'
;;device,/landscape
;!p.font=0
;device,filename='vdist.ps'
;!p.thick=2.0
;!x.thick=2.0
;!y.thick=2.0
!p.multi=[0,2,3]
;!p.charsize=1.5
;@x6x9
;device,/color


;read_traj,t
;rundir = '/Volumes/MacD97-2/hybrid/SWAP/run_test/'
rundir = './tmp1/'
f_read_coord,rundir+'c.coord.dat',x,y,z,dzc,dzg,nx,ny,nz

xyz_traj=[nx/2,ny/2,nz/2]
lxyz = [1.0,0.0,0.0]
nlxyz = [0.0,1.0,0.0]

get_dist,xyz_traj,lxyz,nlxyz

;stop

;t2 = -13.5e3/cos(15*!dtor) - atan(15*!dtor)*xx
;dx = 2200.

;xx = x - x(55)
;yy = y - y(ny/2)

;nstep = 10
;for i = 0,nstep-1 do begin 
;   isc = (nstep-i)*10
;   jsc = ny/2 - round((13.5e3/dx)/cos(15*!dtor)) - atan(15*!dtor)*(isc-(nx-55))
;   print,isc,jsc
;   xyz_traj=[isc,jsc,nz/2-35]
;   get_dist,xyz_traj,lxyz,nlxyz,t
;endfor

;evice,/close

stop
end
;----------------------------------------------------------------------
