local floor, sin, cos, pi2 = math.floor, math.sin, math.cos, math.pi * 2

local fft = {}

local function evalFFT (n, x, w, ox, dx, y, oy, dy, wp, dwp) 
  if n == 1 then 
    x[ox] = y[oy]
    x[ox + 1] = y[oy + 1]
  elseif n == 2 then 
    local ody = oy + dy
    x[ox] = y[oy] + y[ody]
    x[ox + 1] = y[oy + 1] + y[ody + 1]
    x[ox + dx] = y[oy] - y[ody]
    x[ox + dx + 1] = y[oy + 1] - y[ody + 1]
  else 
    local n05 = floor(0.5 * n)
    local nx = n05 * dx
    local onx = ox + nx
    local dy2 = dy + dy
    local dwp2 = dwp + dwp

    evalFFT(n05,w,x,ox,dx,y,oy,dy2,wp,dwp2)
    evalFFT(n05,w,x,onx,dx,y,oy + dy,dy2,wp,dwp2)

    x[ox] = w[ox] + w[onx]
    x[ox + 1] = w[ox + 1] + w[onx + 1]
    x[onx] = w[ox] - w[onx]
    x[onx + 1] = w[ox + 1] - w[onx + 1]

    local i = ox + dx
    local iw = dwp
    while i < onx do 
      local inx = i + nx
      local re = wp[iw] * w[inx] - wp[iw + 1] * w[inx + 1]
      local im = wp[iw] * w[inx + 1] + wp[iw + 1] * w[inx]
      x[i] = w[i] + re
      x[i + 1] = w[i + 1] + im
      x[inx] = w[i] - re
      x[inx + 1] = w[i + 1] - im
      iw = iw + dwp
      i = i + dx
    end
  end
end

fft.dft = function (x, y, skip, ox, oy) 
  skip, ox, oy = skip or 2, ox or 0, oy or 0

  local n = floor((#y - ox) / skip)
  local pi2dn = pi2 / n
  local wp = {}
  local re = cos(pi2dn)
  local im = sin(pi2dn)
  local cx = 1
  local cy = 1
  for i = 1, n do 
    local i2 = i + i
    wp[i2] = cx
    wp[i2 + 1] = cy
    local x2 = cx * re - cy * im
    cy = cx * im + cy * re
    cx = x2
  end

  evalFFT(n,x,{},ox,skip,y,oy,skip,wp,2)
end

fft.idft = function (x, y, skip, ox, oy) 
  skip, ox, oy = skip or 2, ox or 0, oy or 0

  local n = floor((#y - ox) / skip)
  local pi2dn = pi2 / n
  local wp = {}
  local re = cos(pi2dn)
  local im = -sin(pi2dn)
  local cx = 1
  local cy = 1
  for i = 1, n do 
    local i2 = i + i
    wp[i2] = cx
    wp[i2 + 1] = cy
    local x2 = cx * re - cy * im
    cy = cx * im + cy * re
    cx = x2
  end

  evalFFT(n,y,{},oy,skip,x,ox,skip,wp,2)

  local s = 1 / n
  for i = 1, n / 2 do 
    y[i] = y[i] * s
  end
end

fft._eval = evalFFT

return fft


