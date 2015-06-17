pro explore_bmu

header=''

N=2920

bmu = fltarr(3,N)

dir_dat = '/mnt/meteo0/data/dargaville/rhuva/ACCESS/som_pak-3.1/'

;openr, lun, dir_dat+'ei_oper_an_sfc_15x15_6S105E49.5S165E_20100101_20111231.vis', /Get_LUN
openr, lun, dir_dat+'ei_oper_an_sfc_15x15_10.5S105E49.5S165E_20100101_20111231.vis', /Get_LUN

readf, lun, header 
readf, lun, bmu

free_lun, lun

errors=bmu(2,*)
stop
;min_error=min(errors)
;max_error=max(errors)

;LOOP THROUGH THE TIME STEPS IN THE ORIGINAL DATA AND SEARCH THE GRID
;FOR THE MAXIMUM ERROR WHEN COMPARED TO THE BMU.

header=0.
;mslp_aus=fltarr(41,30,N)
mslp_aus=fltarr(41,27,N)

;openr, lun_mslp, dir_dat+'ei_oper_an_sfc_15x15_6S105E49.5S165E_20100101_20111231.dat', /get_lun
openr, lun_mslp, dir_dat+'ei_oper_an_sfc_15x15_10.5S105E49.5S165E_20100101_20111231_zeros.dat', /get_lun

readf, lun_mslp, header
readf, lun_mslp, mslp_aus

free_lun, lun_mslp


;openr, lun_som, dir_dat+'MSLP_aus_region_6x5_alpha0.7_bubble.cod', /get_lun
openr, lun_som, dir_dat+'MSLP_aus_continent_6x5_alpha0.5_bubble.cod', /get_lun

xdim=0
ydim=0
nelements=0
junk=''
junk2=''
junk3=''
header=strarr(1)

readf, lun_som, header ;<-Read in the header
READS, header, junk, FORMAT='(A4)' ;<-Take the first four characters and send them to junk
nelements=fix(junk) ;Convert junk to int 
READS, header, junk2, xdim, FORMAT='(A9,x,I0)'  ;<-Take the 9th character and send it to xdim with type int
READS, header, junk3, ydim, FORMAT='(A11,x,I0)' ;<-Take the 11th character and send it to ydim with type int

nrow = xdim*ydim  ;<-The number of nodes in the SOM

som=fltarr(nelements,nrow)    ;<-Set size of data to be read in (if known)
readf, lun_som, som      ;<-If size of data is known, this will read in and automatically fill all elements from the
free_lun, lun_som        ;current LUN into the array nominated in the readf statement.
stop
;Reform to be size: Lat*Lon*Time or Lat*Lon*(no. of SOM nodes)
;som is currently 1230*12
som=reform(som,[41,27,nrow])  ;<-This appears to read the same way the .dat was written where lons where printed first.
som=reform(som,[41,27,xdim,ydim])
stop
som_comp=fltarr(41,27,N)
max_errors=fltarr(N)

for i=0,N-1 DO BEGIN

    x=bmu(0,i)
    y=bmu(1,i)
    som_comp(*,*,i)=abs(mslp_aus(*,*,i)-som(*,*,x,y))
    max_errors(i)=max(som_comp(*,*,i))

endfor

lt_20=where(max_errors lt 20, ngt)
lt_20_errors=max_errors(lt_20)

lt_20_bmu=bmu(*,lt_20)  ;max(lt_20_bmu(2,*))=

max_sub_bmu=where(max_errors eq max(lt_20_errors),ni) ; This gives , for a max of 

max_lt_20_bmu=where(bmu(2,*) eq max(lt_20_bmu(2,*))) ;This gives 
;print, bmu(*,22708)

;Order lt_20_bmu to judge how the worst few look (not just the worst):

lt_20_distances=lt_20_bmu(2,*)
lt_20_ordered=sort(lt_20_distances)  ;Sort the Euclidean distances into increasing order
lt_20_ordered_distances=lt_20_distances(lt_20_ordered)
lt_20_ordered_bmu=lt_20_bmu(*,lt_20_ordered)
stop

Ni=fix(n_elements(lt_20_distances)*0.95)

filter_bmu_test=fltarr(Ni)
points_gt_four=fltarr(Ni)

icount=0
for i=0,Ni-1 DO BEGIN

  filter_bmu=where(bmu(2,*) eq lt_20_ordered_bmu(2,icount)) ;Filter for the first (best) 95%
    if(n_elements(filter_bmu) gt 1) then BEGIN
        N=n_elements(filter_bmu)
        for j=0,N-1 DO BEGIN
            filter_bmu_test(icount)=filter_bmu(j)
            icount=icount+1
        endfor
        icount=icount-1
    endif ELSE BEGIN
        filter_bmu_test(icount)=filter_bmu
    endelse
;if(icount gt 10074) then stop
if(i ne icount) then i=icount

icount=icount+1
  

endfor

good_enough=som_comp(*,*,filter_bmu_test)   ;Find the corresponding array of differences
good_enough=reform(good_enough,[1107,Ni]) ;reform to 1-D of space to use where function

for j=0,Ni-1 DO BEGIN

    gt_four=where(good_enough(*,j) gt 4)  ;Determine how many grid points have a difference gt 1
    points_gt_four(j)=n_elements(gt_four) ;Sum the amount

endfor
stop
survivors=where(points_gt_four le points_gt_four(Ni-1)) ;Determine how many bmu have less than/equal to ... points gt 1
                                ;The above returns the array indeces for the sorted array                     
full_filtered_subs=filter_bmu_test(survivors) ;Filter for those that survive where statement
                                ;(where filter_bmu_test has the indeces in the orginal data that survived the first round of filtering)
full_filtered_bmu=bmu(*,full_filtered_subs)  ;Filter the original bmu array

stop


;Note that the first number for each row in bmu contains the x-dim
;number (column number in the SOM), second number is the y-dim (row number) in
;the SOM and the final number is the Euclidean distance to the best
;matching node for the corresponding time step in the data.


;FOR THE 6x5 ALPHA 0.7 BUBBLE RUN:-----------------------------------------------------------
;The best bmu (without filtering) occurs at time step 4621 (using qerrors)
;The worst bmu (without filtering) occurs at time step 2755 (using qerrors)

;The mean maximum discrepancy between each bmu and the corresponding
;map from the ERA-Interim data is 16.94hPa (median is 16.17). The
;maximum discrepancy, however, is 52.1hPa and the minimum 5.0hPa.

;There are 5287 time steps that have a maximum error less than 20hPa. Of
;these time steps the maximum (19.992) happens once: at time step
;5287 (when counting from zero in the original mslp
;data). The maximum Euclidean distance across the whole grid from those time steps that have
;a maximum error less than 20hPa happens at time step 6695 (in the
;original mslp jja data) and that distance is 283.181 (between node (4,0) and mslp_aus_jja(*,*,6695)).
;This is much more than the average Euclidean distance for all bmu, which is  163.355.








;FOR THE 6x4 ALPHA 0.001 BUBBLE RUN:-----------------------------------------------------------
;The mean maximum discrepancy between each bmu and the corresponding
;map from the ERA-Interim data is 17.54hPa (median is 16.9). The
;maximum discrepancy, however, is 58.29hPa and the minimum 3.82001hPa.

;There are 21917 time steps that have a maximum error less than 20hPa. Of
;these time steps the maximum (19.9999) happens once: at time step
;27286 (when counting from zero in the original mslp
;data). The maximum Euclidean distance across the whole grid from those time steps that have
;a maximum error less than 20hPa happens at time step 22708 (in the original mslp data) 
;and that distance is 325.484 (between node (4,0) and mslp_aus(*,*,22708)).
;This is much more than the average Euclidean distance for all bmu, which is 163.167.



end
