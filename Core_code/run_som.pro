pro run_som

;**************************************
;Use this to automate the training    *
;and creation of a SOM by calling the *
;SOM_PAK and iterating                *
;**************************************

;Alpha_type was tested and despite recommendations from the software
;the inverse_t tended to give higher quantised errors.


;****************************************************
;Hard-coded variables-change when mapping new fields:--------------------------------------------------------------------  
;****************************************************

data='MSLP_aus_continent_x'

;data='u80_aus_sample'


init_rlen='30680'    ;The number of samples (time-steps) in the data
final_rlen='306800'  ;10 times as many samples


;****************************************
;Invariant variables + SOM training code:--------------------------------------------------------------------------------
;****************************************


alpha_init_1=(findgen(9)+1)/1000.
alpha_init_2=(findgen(9)+1)/100.
alpha_init_3=[0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,1]

alpha_final_1=0.5*alpha_init_1
alpha_final_2=0.5*alpha_init_2
alpha_final_3=0.5*alpha_init_3

init_alpha=fltarr(n_elements(alpha_init_1)+n_elements(alpha_init_2)+n_elements(alpha_init_3))
init_alpha(0:n_elements(alpha_init_1)-1)=alpha_init_1
init_alpha(n_elements(alpha_init_1):n_elements(alpha_init_1)+n_elements(alpha_init_2)-1)=alpha_init_2
init_alpha(n_elements(alpha_init_1)+n_elements(alpha_init_2):n_elements(init_alpha)-1)=alpha_init_3

final_alpha=fltarr(n_elements(init_alpha))
final_alpha(0:n_elements(alpha_init_1)-1)=alpha_final_1
final_alpha(n_elements(alpha_init_1):n_elements(alpha_init_1)+n_elements(alpha_init_2)-1)=alpha_final_2
final_alpha(n_elements(alpha_init_1)+n_elements(alpha_init_2):n_elements(init_alpha)-1)=alpha_final_3

init_alpha=string(init_alpha)
final_alpha=string(final_alpha)



din=data+'.dat'

cout=data+'.cod'

;xdim=[3,4,5,6,7,8]
;ydim=xdim-2

xdim=[5,6,7]
ydim=xdim-1

xdim=string(xdim)
ydim=string(ydim)  

topol='rect'

;neigh='gaussian'

neigh='bubble'

;ntrials=20

ntrials=4

seed=indgen(ntrials)
seed=seed+1
seed=string(seed)

final_radius='1'

error=fltarr(1)

total_error=fltarr(ntrials,n_elements(xdim),n_elements(init_alpha))

for i=0,ntrials-1 DO BEGIN
print, i
    for j=0,n_elements(xdim)-1 DO BEGIN
print, j
        for k=0,n_elements(init_alpha)-1 DO BEGIN

            
spawn, './randinit '+'-din '+din+' '+'-cout '+cout+' '+'-xdim '+xdim(j)+' '+'-ydim '+ydim(j)+' '+'-topol '+topol+' '+'-neigh '+neigh+' '+'-rand '+seed(i)+' >&/dev/null'

spawn, './vsom '+'-din '+din+' '+'-cin '+cout+' '+'-cout '+cout+' '+'-rlen '+init_rlen+' '+'-alpha '+init_alpha(k)+' '+'-radius '+ydim(j)+' >&/dev/null'

spawn, './vsom '+'-din '+din+' '+'-cin '+cout+' '+'-cout '+cout+' '+'-rlen '+final_rlen+' '+'-alpha '+final_alpha(k)+' '+'-radius '+final_radius+' >&/dev/null'

spawn, './qerror '+'-din '+din+' '+'-cin '+cout+' > error.txt'

read_error,error

total_error(i,j,k)=error

       endfor
     
    endfor

endfor


min_error=min(total_error)

for i=0,ntrials-1 DO BEGIN

    for j=0,n_elements(xdim)-1 DO BEGIN

        for k=0,n_elements(init_alpha)-1 DO BEGIN

            if(total_error(i,j,k) eq min_error) then BEGIN
                print, 'THIS IS MIN SEED SUBSCRIPT: ',i
                seed_subscript=i
                print, 'THIS IS MIN XDIM SUBSCRIPT: ',j
                xdim_subscript=j
                print, 'THIS IS MIN ALPHA SUBSCRIPT: ',k
                alpha_subscript=k
            endif

        endfor
       
     endfor

endfor

;NOW RUN THE PROGRAM AGAIN, EXCEPT WITH THE BEST CONFIGURATION

            
spawn, './randinit '+'-din '+din+' '+'-cout '+cout+' '+'-xdim '+xdim(xdim_subscript)+' '+'-ydim '+ydim(xdim_subscript)+' '+'-topol '+topol+' '+'-neigh '+neigh+' '+'-rand '+seed(seed_subscript)+' >&/dev/null'

spawn, './vsom '+'-din '+din+' '+'-cin '+cout+' '+'-cout '+cout+' '+'-rlen '+init_rlen+' '+'-alpha '+init_alpha(alpha_subscript)+' '+'-radius '+ydim(xdim_subscript)+' >&/dev/null'

spawn, './vsom '+'-din '+din+' '+'-cin '+cout+' '+'-cout '+cout+' '+'-rlen '+final_rlen+' '+'-alpha '+final_alpha(alpha_subscript)+' '+'-radius '+final_radius+' >&/dev/null'


;TO BE CONSISTENT WITH ALEXANDER ET AL. 2010 RECORD THE 'BEST' SOM FOR EACH SOM SIZE:

get_lun,som_lun
openw,som_lun,'best_soms_'+data+'_bubble_II.txt'

for i=0,ntrials-1 DO BEGIN

    for j=0,n_elements(xdim)-1 DO BEGIN

        best_som=min(total_error(*,j,*))

        for k=0,n_elements(init_alpha)-1 DO BEGIN
            
            if(total_error(i,j,k) eq best_som) then BEGIN
                printf,som_lun,'For xdim: ',xdim(j),', The best SOM had seed and trial no.: ',i,' and initial alpha of: ',init_alpha(k)
            endif

        endfor

    endfor

endfor

free_lun, som_lun

stop
end

pro read_error,error

openr, error_lun, 'error.txt', /Get_Lun

all_error=''
error=''

readf,error_lun, all_error

free_lun, error_lun

error_pos=''

error_pos=strpos(all_error, 'is')
error_pos=error_pos+2

errors=strmid(all_error,error_pos,9)

error=float(errors)


end
