function z0mgs_ellfunc, p, x=x, y=y, err=err
;+
;
; x - theta
; y - radius
; p - posang, major, minor
;
;-
  
  model = $
     p[1]*p[2] $
     / sqrt((p[1]*cos(x-p[0]))^2 + $
            (p[2]*sin(x-p[0]))^2)
  
  return, (y-model)/err

end
