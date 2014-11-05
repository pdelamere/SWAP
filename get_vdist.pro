;----------------------------------------------------------------
PRO read_part,file,nfrm,Ni_max,xp
;----------------------------------------------------------------

;Ni_max=long(0)
;nt=0l
;ntout=0l
frm=0l

file = file+'.dat'
print,' reading...',file
openr,1,file,/f77_unformatted
;readu,1,nt
;readu,1,ntout
;readu,1,Ni_max
;print,nt,ntout,Ni_max




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
PRO read_part_scalar,file,nfrm,Ni_max,xp
;----------------------------------------------------------------

;Ni_max=long(0)
;nt=0
;ntout=0
frm=0

file = file+'.dat'
print,' reading...',file
openr,1,file,/f77_unformatted
;readu,1,nt
;readu,1,ntout
;readu,1,Ni_max
;print,nt,ntout,Ni_max

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
pro get_dist,xcur,ycur,zcur,x,y,z,xp,vp,mrat,beta_p,np,ndx,w,oVid,vidStream,b1
;----------------------------------------------------------------

fnsz = 16

dv = 2.0
nn = 2000

vxyp = fltarr(nn/dv,nn/dv)
vxzp = fltarr(nn/dv,nn/dv)
vyzp = fltarr(nn/dv,nn/dv)
v_perp_par = fltarr(nn/dv,nn/dv)
;vpar = fltarr(nn/dv,nn/dv)

vx = -(nn/(2*dv))*dv + findgen(nn/dv)*dv
vy = -(nn/(2*dv))*dv + findgen(nn/dv)*dv
vz = -(nn/(2*dv))*dv + findgen(nn/dv)*dv
v_perp = -(nn/(2*dv))*dv + findgen(nn/dv)*dv
v_par = -(nn/(2*dv))*dv + findgen(nn/dv)*dv

dx = x(1)-x(0)

j = ycur
i = xcur
k = zcur
  
wh = where((xp(*,0) ge x(i)-ndx*dx) and (xp(*,0) le x(i)+ndx*dx) and $
           (xp(*,1) ge y(j)-ndx*dx) and (xp(*,1) le y(j)+ndx*dx) and $ 
           (xp(*,2) ge z(k)-ndx*dx) and (xp(*,2) le z(k)+ndx*dx))
print,n_elements(wh),z(k)

if (wh(0) gt -1) then begin
   e_arr = 0
   cnt_arr = 0
   for l = 0ll,n_elements(wh)-1 do begin
      ii = fix(vp(wh(l),0)/dv) + (nn/(2*dv))
      jj = fix(vp(wh(l),1)/dv) + (nn/(2*dv))
      kk = fix(vp(wh(l),2)/dv) + (nn/(2*dv))
;      vxyp(ii,jj) = 1/beta_p(wh(l)) + vxyp(ii,jj)
;      vxzp(ii,kk) = 1/beta_p(wh(l)) + vxzp(ii,kk)
;      vyzp(jj,kk) = 1/beta_p(wh(l)) + vyzp(jj,kk)
      vxyp(ii,jj) = 1.0 + vxyp(ii,jj)
      vxzp(ii,kk) = 1.0 + vxzp(ii,kk)
      vyzp(jj,kk) = 1.0 + vyzp(jj,kk)


      ii = fix(sqrt((vp(wh(l),0)-400.)^2 + vp(wh(l),1)^2)/dv) + (nn/(2*dv))
      ; do this in the plasma rest frame
;      v_perp_par(kk,ii) = 1/beta_p(wh(l)) + v_perp_par(kk,ii)
      v_perp_par(kk,ii) = 1.0 + v_perp_par(kk,ii)

;      vpp = reform(vp(wh(l),*))
;      vpp2 = sqrt(vpp(0)^2 + vpp(1)^2 + vpp(2)^2)
;      vpp1 = vpp/vpp2
      
;      vdotl = transpose(vpp1(*))#[-1,0,0]
      
;      if (vdotl gt cos(80*!dtor)) then begin
;         e_arr = [e_arr,(vpp(0)^2 + vpp(1)^2 + vpp(2)^2)/mrat(wh(l))]
;         cnt_arr = [cnt_arr,vdotl]
;;         cnt_arr = [cnt_arr,1.0]
;      endif
      
   endfor
endif


w.erase

sz = size(vxzp)
im1 = image(smooth(alog(vxzp>1),2),vx,vz,rgb_table=27,$
            axis_style=2,xtickdir=1,ytickdir=1,$
            xtitle='$v_x$',ytitle='$v_z$',layout=[2,2,1],/current,$
            font_size=fnsz,xrange=[-100,900],yrange=[-500,500],margin=0.15)

;   scl = 500./sz(1)
scl=1
im1.scale,scl,scl,1

im1 = image(smooth(alog(vxyp>1),2),vx,vz,rgb_table=27,$
            axis_style=2,xtickdir=1,ytickdir=1,$
            xtitle='$v_x$',ytitle='$v_y$',layout=[2,2,2],/current,$
            font_size=fnsz,xrange=[-100,900],yrange=[-500,500],margin=0.15)

im1.scale,scl,scl,1


;im1 = image(smooth(alog(vyzp>1),2),vx,vz,rgb_table=27,$
;            axis_style=2,xtickdir=1,ytickdir=1,$
;            xtitle='$v_y$',ytitle='$v_z$',layout=[2,2,3],/current,$
;            font_size=fnsz,xrange=[-500,500],yrange=[-500,500],margin=0.15)

;im1.scale,scl,scl,1

p = plot(b1(1,1,*,1)^2/b1(1,1,*,2)^2,layout=[2,2,3],margin=0.2,/current,$
        font_size=fnsz,ytitle='$B_y^2/B_o^2$',xtitle='z',yrange=[0,0.4])

im1 = image(smooth(alog(v_perp_par>1),2),v_par,v_perp,rgb_table=27,$
            axis_style=2,xtickdir=1,ytickdir=1,$
            ytitle='$v_{perp}$',xtitle='$v_{parallel}$',layout=[2,2,4],$
            /current,font_size=fnsz,yrange=[0,500],xrange=[-500,500],$
           aspect_ratio=1.0,margin=0.15)

im1.scale,scl,scl,1

for i = 0,10 do begin
   time = oVid.Put(vidStream, w.CopyWindow())
endfor

;im1 = image(smooth(np,2),rgb_table=33,layout=[2,2,4],/current,axis_style=2,$
;            xtickdir=1,ytickdir=1,font_size=fnsz)

;p = plot([0,i],[0,j],symbol="X",layout=[2,2,4],/current,overplot=1,$
;         linestyle=6,thick=2,sym_color='white')

;img = im1.CopyWindow()

;   tvscl,img,true=1
;   xinteranimate, frame = cnt, image = img

;   w.close


return
end
;----------------------------------------------------------------------



;----------------------------------------------------------------------
;main program
;----------------------------------------------------------------------

!p.multi=[0,1,1]

mp=1.67e-27
nfile = '1'
nfm = 16
procnum=1
ndx = 4.0

rio = 1800./40.

xsz=1100
ysz=1000

file = 'vdist.mp4'
width = xsz
height = ysz
frames = 180
fps = 30
speed = 2
samplerate = 22050L

        ; Create object and initialize video/audio streams
oVid = IDLffVideoWrite(file)
vidStream = oVid.AddVideoStream(width, height, fps)
        ; audStream = oVid.AddAudioStream(samplerate)


;restore,'/Volumes/Scratch/Dols_output/trajectory_0_24_27_31.sav'
;dir = '/Volumes/MacD97-2/hybrid/SWAP/run_test/'
dir = './tmp1/'
;dir = '/Volumes/MacD97-2/hybrid/3d_buf/run_test/'

read_para,dir

restore,filename=dir+'para.sav'

read_coords,dir,x,y,z

xx = (x - x(nx/2))/rio
yy = (y - y(ny/2))/rio

device,decomposed=0
loadct,39

w = window(dimensions=[xsz,ysz])

c_read_3d_m_32,dir,'c.np',nfm,np

plot,z,np(1,1,*)/1e15

for nfrm = 1,nfm do begin

nfil = 0
xfile = dir+'c.xp_'+strtrim(string(procnum),2)
vfile = dir+'c.vp_'+strtrim(string(procnum),2)
mratfile = dir+'c.mrat_'+strtrim(string(procnum),2)
beta_p_file = dir+'c.beta_p_'+strtrim(string(procnum),2)

read_part,xfile,nfrm,Ni_max,xp
read_part,vfile,nfrm,Ni_max,vp
read_part_scalar,mratfile,nfrm,Ni_max,mrat
read_part_scalar,beta_p_file,nfrm,Ni_max,beta_p


for nfil = 1,2 do begin
    xfile = dir+'c.xp_'+strtrim(string(procnum),2)
    vfile = dir+'c.vp_'+strtrim(string(procnum),2)
    mratfile = dir+'c.mrat_'+strtrim(string(procnum),2)
    beta_p_file = dir+'c.beta_p_'+strtrim(string(procnum),2)

    read_part,xfile,nfrm,Ni_max,xpp
    read_part,vfile,nfrm,Ni_max,vpp
    read_part_scalar,mratfile,nfrm,Ni_max,mratt
    read_part_scalar,beta_p_file,nfrm,Ni_max,beta_pp

    xp = [xp,xpp]
    vp = [vp,vpp]
    mrat = [mrat,mratt]
    beta_p = [beta_p,beta_pp]

endfor

c_read_3d_vec_m_32,dir,'c.b1',nfrm,b1

get_dist,1,1,nz/2,x,y,z,xp,vp,mrat,beta_p,np,ndx,w,oVid,vidStream,b1

endfor

oVid.Cleanup


end
;----------------------------------------------------------------------
