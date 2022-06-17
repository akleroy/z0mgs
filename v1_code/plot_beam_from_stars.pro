pro plot_beam_from_stars

  plot, findgen(10), xtitle='!6'

  for jj = 0, 1 do begin
     if jj eq 0 then begin
        res_str = 'gauss7p5'
     endif
     if jj eq 1 then begin
        res_str = 'gauss15'
     endif

     restore, '../measurements/star_stacks/star_stack_w1_'+res_str+'.idl', /v
     sz = size(stack)
     norm = stack
     ind = sort(val)
     for ii = 0, sz[3]-1 do begin
        norm[*,*,ind[ii]] = stack[*,*,ind[ii]]/val[ind[ii]]
     endfor  
     w1_stack = norm
     w1_beam = median(norm,dim=3)

     restore, '../measurements/star_stacks/star_stack_w2_'+res_str+'.idl', /v
     sz = size(stack)
     norm = stack
     ind = sort(val)
     for ii = 0, sz[3]-1 do begin
        norm[*,*,ind[ii]] = stack[*,*,ind[ii]]/val[ind[ii]]
     endfor  
     w2_stack = norm
     w2_beam = median(norm,dim=3)

     restore, '../measurements/star_stacks/star_stack_w3_'+res_str+'.idl', /v
     sz = size(stack)
     norm = stack
     ind = sort(val)
     for ii = 0, sz[3]-1 do begin
        norm[*,*,ind[ii]] = stack[*,*,ind[ii]]/val[ind[ii]]
     endfor  
     w3_stack = norm
     w3_beam = median(norm,dim=3)

     if res_str eq 'gauss15' then begin
        restore, '../measurements/star_stacks/star_stack_w4_'+res_str+'.idl', /v
        sz = size(stack)
        ind = sort(val)
        for ii = 0, sz[3]-1 do begin
           norm[*,*,ind[ii]] = stack[*,*,ind[ii]]/val[ind[ii]]
        endfor  
        w4_stack = norm
        w4_beam = median(norm,dim=3)
     endif else begin
        w4_beam = !values.f_nan*w1_beam
     endelse

     restore, '../measurements/star_stacks/star_stack_nuv_'+res_str+'.idl', /v
     sz = size(stack)
     norm = stack
     ind = sort(val)
     for ii = 0, sz[3]-1 do begin
        norm[*,*,ind[ii]] = stack[*,*,ind[ii]]/val[ind[ii]]
     endfor  
     nuv_stack = norm
     nuv_beam = median(norm,dim=3)

     restore, '../measurements/star_stacks/star_stack_fuv_'+res_str+'.idl', /v  
     sz = size(stack)
     norm = stack  
     ind = sort(val)
     for ii = 0, sz[3]-1 do begin
        norm[*,*,ind[ii]] = stack[*,*,ind[ii]]/val[ind[ii]]
     endfor  
     fuv_stack = norm
     fuv_beam = median(norm,dim=3)

     if jj eq 0 then begin
        w1_beam7p5 = w1_beam
        w2_beam7p5 = w2_beam
        w3_beam7p5 = w3_beam
        w4_beam7p5 = w4_beam
        fuv_beam7p5 = fuv_beam
        nuv_beam7p5 = nuv_beam
     endif else begin
        w1_beam15 = w1_beam
        w2_beam15 = w2_beam
        w3_beam15 = w3_beam
        w4_beam15 = w4_beam
        fuv_beam15 = fuv_beam
        nuv_beam15 = nuv_beam
     endelse

     x = findgen(sz[1])
     x -= mean(x)
     xgrid = x # (fltarr(sz[2])+1.)
     y = findgen(sz[2])
     y -= mean(y)
     ygrid =(fltarr(sz[1])+1.) # (y)
     r = sqrt(xgrid^2+ygrid^2)
     pix = 2.75

  endfor

  for jj = 0, 1 do begin

     if jj eq 0 then begin
        res_str = 'gauss15' 
        res_label = 'FWHM 15"'
        fwhm = 15.
     endif

     if jj eq 1 then begin
        res_str = 'gauss7p5' 
        res_label = 'FWHM 7.5"'
        fwhm = 7.5
     endif

     plot, findgen(10)
     loadct, 0
     psfile = '../plots/beam_from_stars_'+res_str+'.eps'
     pnfile = '../plots/beam_from_stars_'+res_str+'.png'
     ps, /def, /ps, xs=5, ys=5, /color, /encaps $
         , file=psfile
     
     plot $
        , [0], [0], /nodata $
        , ytitle='!6log!d10!n Fraction of Peak Intensity' $
        , xtitle='!6Radius ["]' $
        , xthick=5, ythick=5, charthick=3, charsize=1.25 $
        , xrange=[0., 100.], yrange=[1d-5, 1.], /ylo $
        , /xstyle
     
     for ii = 0, 20 do $
        oplot, ii*10*[1,1], [1d-10, 1d10], lines=1

     for ii = -10, 10 do $
        oplot, [0,1000], 10.^(ii)*[1, 1], lines=1

     
     for ii = 0, 5 do begin

        if ii eq 0 then begin
           name = 'wise1'
           if jj eq 0 then beam = w1_beam15
           if jj eq 1 then beam = w1_beam7p5
           color='firebrick'
        endif
        if ii eq 1 then begin
           name = 'wise2'
           if jj eq 0 then beam = w2_beam15
           if jj eq 1 then beam = w2_beam7p5
           color='salmon'
        endif
        if ii eq 2 then begin
           name = 'wise3'
           if jj eq 0 then beam = w3_beam15
           if jj eq 1 then beam = w3_beam7p5
           color='goldenrod'
        endif
        if ii eq 3 then begin
           name = 'wise4'
           if jj eq 0 then beam = w4_beam15
           if jj eq 1 then beam = w4_beam7p5
           color='lightseagreen'
        endif
        if ii eq 4 then begin
           name = 'nuv'
           if jj eq 0 then beam = nuv_beam15
           if jj eq 1 then beam = nuv_beam7p5
           color='dodgerblue'
        endif
        if ii eq 5 then begin
           name = 'fuv'
           if jj eq 0 then beam = fuv_beam15
           if jj eq 1 then beam = fuv_beam7p5
           color='orchid'
        endif
        
        ;oplot, r*pix, beam, color=cgcolor(color), ps=1, symsize=1
        bins = bin_data(r*pix, beam, /nan $
                        , xmin=-1.0, xmax=150., binsize=2.0)
        oplot, bins.xmid, (bins.ymean > 1d-6), color=cgcolor(color), thick=10

        fid = findgen(201)
        oplot, fid, (exp(-1.*(fid)^2/2./(fwhm/2.354)^2)), color=cgcolor('darkgray') $
               , thick=10, lines=2
                
     endfor

     al_legend, /top, /right $
                , box=1, clear=1, background='lightgray' $
                , lines=-99, [res_label] $
                , textcolor=[cgcolor('black')] $
                , charthick=3, charsize=1.25
    
     if res_str eq 'gauss15' then begin
        al_legend, /bottom, /left $
                   , box=0, clear=0, lines=-99 $
                   , ['WISE1','WISE2','WISE3','WISE4','NUV','FUV'] $
                   , textcolor=[cgcolor('firebrick'),cgcolor('salmon'),cgcolor('goldenrod'),cgcolor('lightseagreen'),cgcolor('dodgerblue'),cgcolor('purple')] $
                   , charthick=3, charsize=1.0
     endif else begin
        al_legend, pos=[70,1d-1] $
                   , box=0, clear=0, lines=-99 $
                   , ['WISE1','WISE2','WISE3','WISE4','NUV','FUV'] $
                   , textcolor=[cgcolor('firebrick'),cgcolor('salmon'),cgcolor('goldenrod'),cgcolor('lightseagreen'),cgcolor('dodgerblue'),cgcolor('purple')] $
                   , charthick=3, charsize=1.0
     endelse

     ps, /xw
     spawn, 'evince '+psfile+' &'       
     spawn, 'convert -density 300x300 '+psfile+' '+pnfile+' &'
     
  endfor

  stop

end
