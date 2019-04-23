;NAME: straighten_loop

;PURPOSE: to staighten a coronal loop and save it to a datacube

;INPUTS: loop_datacube: a 3D image cube ([x,y,t]);

;KEYWORDS:  x0box & y0box: the position of the lower-left corner of the region you want to show (in pixel)
;           nxbox & nybox: the size of x-axis and y-axis of the region you choose to show
;           log: to set if the data should be in log scale

;OUTPUTS: datanew: the datacube containing straightened loop;
;         xfit & yfit: fitting result of the polynomial fitting;
;         xdata & ydata: the points you clicked

;NOTES: the procedure first need you to choose one frame you want to plot as the reference image,
;       then you should decide if the region should be rotated;
;       after that, click the left button to select points you want to trace the loop;
;       a polynomial fit (with degrees set by the user) will be done;
;       for each time, the pixels around the polyfit line will be converted to a straightened column.

;AUTHOR: Written by Zihao Yang (Peking University) on April 23rd, 2019;
;        Based on the code 'loop_width.pro' written by J. A. Klimchuk on April 22nd, 1992. The original 
;           code is implemented in SSW, under the branch of 'Yohkoh'.



pro straighten_loop, loop_datacube, x0box=x0box,y0box=y0box,nxbox=nxbox,nybox=nybox,log=log
restore,loop_datacube
help,data 
read,'which frame do you want to show? (0 for 1st image):',inum

int=alog10(data)
if not keyword_set(log) then int=data
window,0
plot_image,reform(int[*,*,inum])

if not keyword_set(x0box) then x0box=0
if not keyword_set(y0box) then y0box=0
if not keyword_set(nxbox) then nxbox=(size(data))[1]
if not keyword_set(nybox) then nybox=(size(data))[2]
nt=(size(data))[3]

sub=int[x0box:x0box+nxbox-1,y0box:y0box+nybox-1,inum]

read, 'Input magnification factor:', mag
sub_expand = rebin(sub, mag*nxbox, mag*nybox, /sample)

window,1, xpos=450,ypos=5,xsize=750,ysize=600
wset,1
tvscl,sub_expand

print, 'Loop must be vertically single valued;'
read, '   do you wish to rotate the image? (0-no; 1-yes)', irotate
if (irotate eq 1) then begin
   read, 'Input rotation angle (degrees clockwise):', angle
   sub_rotate = rot(sub, angle, /interp)
   sub = sub_rotate
   sub_comp = rebin(sub, mag*nxbox, mag*nybox, /sample)
   tvscl, sub_expand
endif

xtemp = intarr(100)
ytemp = xtemp
print, 'Select loop axis pixels with left button (starting with'
print, '  leftmost loop footpoint).'
print, '  Click right button when done.'
cursor, xx, yy, /device, /down
!err = 1               ; equal to 1 for left button clicks only
i = 0
while (!err eq 1) do begin
   xtemp(i) = xx/mag
   ytemp(i) = yy/mag
   print, i+1, ':   x, y  =', xtemp(i), ytemp(i)
   i = i + 1
   cursor, xx, yy, /device, /down
endwhile
xdata = xtemp(0: i-1)
ydata = ytemp(0: i-1)

; Compute polynomial fit to selected axis pixels and plot
;
read, 'Input degree of polynomial fit:', degree
coef = poly_fit(xdata, ydata, degree)
xlength = xdata(i-1) - xdata(0)
xfit = indgen(100)*xlength/100. + xdata(0)
yfit = coef(0)
for n = 1, degree do yfit = yfit + coef(n)*xfit^n
;
window,2,xpos=390,ypos=440,xsize=500,ysize=400
wset,2
plot, xfit, yfit, Title='Loop Axis', Xtitle='X', Ytitle='Y'
!psym=2
oplot,xdata,ydata
!psym=0

; Determine x range of data to be included in straighted loop image
; (based on an extension of the loop)
;
extend = 3.0                  ; amount of extension (in pixels)
dfdx = coef(1)                ; derivative of axis function
for n = 2, degree do dfdx = dfdx + n*coef(n)*xdata(0)^(n-1)
x0 = xdata(0) - extend/sqrt(1.0 + dfdx^2)
dfdx = coef(1)  
for n = 2, degree do dfdx = dfdx + n*coef(n)*xdata(i-1)^(n-1)
xmax = xdata(i-1) + extend/sqrt(1.0 + dfdx^2)
ymin = min(ydata) - extend
ymax = max(ydata) + extend

; Specify width of straightened loop image
;
jump1:
read, 'Input width of straightened loop image (odd):', nwide
noffset = (nwide - 1)/2
if (x0 lt noffset) or (xmax gt nxbox-noffset-1) or (ymin lt noffset) $
   or (ymax gt nybox-noffset-1) then begin
   print, 'Choose smaller width.'
   goto, jump1
endif else print, 'OK'

datatemp = fltarr(nwide,300,nt)
sub_temp=int[x0box:x0box+nxbox-1,y0box:y0box+nybox-1,0:nt-1]

x0_temp=x0
for tt=0,nt-1 do begin
j = 0
while (x0 le xmax) do begin       
;
; Determine direction of axis normal
   y0 = coef(0)
   for n = 1, degree do y0 = y0 + coef(n)*x0^n
   dfdx = coef(1)                       ; derivative of axis function
   for n = 2, degree do dfdx = dfdx + n*coef(n)*x0^(n-1)
   deltax = 1.0/sqrt(1.0 + 1.0/dfdx^2)
   deltay = -deltax/dfdx
   if (dfdx lt 0) then deltax = -deltax 
   if (dfdx lt 0) then deltay = -deltay
;
; Step along axis normal
   for i = 0, nwide-1 do begin     
      x = x0 + (i - noffset)*deltax   
      y = y0 + (i - noffset)*deltay
      xm = fix(x)
      ym = fix(y)
;
; Compute weighted mean intensity
      if (x eq xm) and (y eq ym) then datatemp(i,j,tt) = sub_temp(xm,ym,tt) $
      else begin             
         xp = xm + 1
         yp = ym + 1
         rmm = sqrt((x - xm)^2 + (y - ym)^2)
         rmp = sqrt((x - xm)^2 + (y - yp)^2)
         rpm = sqrt((x - xp)^2 + (y - ym)^2)
         rpp = sqrt((x - xp)^2 + (y - yp)^2)
         rnorm = 1/rmm + 1/rmp + 1/rpm + 1/rpp
         datatemp(i,j,tt) = (sub_temp(xm,ym,tt)/rmm + sub_temp(xm,yp,tt)/rmp + $
                       sub_temp(xp,ym,tt)/rpm + sub_temp(xp,yp,tt)/rpp)/rnorm
      endelse
   end
;save,filename='xy.sav',x,y,xm,ym
; Step along loop axis
   for i = 0, 9 do begin
      deltax = 0.1/sqrt(1.0 + dfdx^2)
      x0 = x0 + deltax
      dfdx = coef(1)
      for n = 2, degree do dfdx = dfdx + n*coef(n)*x0^(n-1)
   end
   j = j + 1
endwhile
x0=x0_temp
endfor

;
nlong = j
datanew = datatemp(0:nwide-1, 0:nlong-1,*)
;datanew_expand = rebin(datanew,mag*nwide,mag*nlong,nt,/sample)
;
window,3,xpos=0,ypos=390,xsize=200,ysize=450
wset,3
; tv, bytscl(datanew_expand, min=0.0, max=intmax)
tvscl,reform(datanew[*,*,inum])

save,filename='straightened_loop.sav',datanew,xfit,yfit,xdata,ydata

end
