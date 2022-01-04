function odemodel!(dydt, initial, params, t)
    b, dm, kb, ku, f, thetar, k_cm, s0, gmax, cl, thetax, Kt, M, we, Km, vm, nx, Kq, Kp, vt, wr, wq, wp, hq, nr, ns = params
    cr, em, cp, cq, ct, et, cm, mt, mm, q, p, si, mq, mp, mr, r, a = initial

    dcr, dem, dcp, dcq, dct, det, dcm, dmt, dmm, dq, dp, dsi, dmq, dmp, dmr, dr, da = zeros(length(dydt))
    
    Kgamma = gmax/Kp
    gamma = gmax*a/(Kgamma+a)
    ttrate = (cq+cr+cp+ct+cm)gamma
    lam = ttrate/M
    vimp = (et*vt*s0/(Kt+s0))
    nucat = (em*vm*si/(Km+si))

    dydt[1] =  +r*mr*kb - cr*ku - cr*lam - gamma/nr*cr - f*cr
    dydt[2] = - lam*em + cm*gamma/nx
    dydt[3] =   +r*mp*kb - cp*ku - cp*lam - gamma/nx*cp - f*cp
    dydt[4] =  +r*mq*kb - cq*ku - cq*lam - gamma/nx*cq - f*cq
    dydt[5] =  +r*mt*kb - ct*ku - ct*lam - gamma/nx*ct - f*ct
    dydt[6] =  - lam*et + gamma/nx*ct
    dydt[7] =   +r*mm*kb - cm*ku - cm*lam - gamma/nx*cm - f*cm
    dydt[8] =   +(we*a/(thetax+a)) - mt*dm - mt*lam - r*mt*kb + ct*ku + gamma/nx*ct
    dydt[9] =   +(we*a/(thetax+a)) - mm*dm - mm*lam - r*mm*kb + cm*ku + gamma/nx*cm
    dydt[10] = - lam*q + gamma/nx*cq
    dydt[11] =  - lam*p + gamma/nx*cp
    dydt[12] = - lam*si + vimp - nucat
    dydt[13] =   +(wq*a/(thetax+a)/(1+(q/Kq)^hq)) - mq*dm - mq*lam - r*mq*kb + cq*ku + gamma/nx*cq
    dydt[14] = +(wp*a/(thetax+a)) - mp*dm - mp*lam - r*mp*kb + cp*ku + gamma/nx*cp
    dydt[15] =   +(wr*a/(thetar+a)) - mr*dm - mr*lam - r*mr*kb + cr*ku + gamma/nr*cr
    dydt[16] =  - lam*r - r*mr*kb - r*mt*kb - r*mm*kb - r*mq*kb - r*mp*kb + cr*ku + ct*ku + cm*ku + cq*ku + cp*ku + gamma/nr*cr + gamma/nr*cr + gamma/nx*ct + gamma/nx*cm + gamma/nx*cq + gamma/nx*cp
    dydt[17] =  +ns*nucat - ttrate - lam*a
end

function odemodelfull!(dydt, initial, params, t)
    b, dm, kb, ku, cl, k_cm, f, thetar, s0, gmax, thetax, Kt, M, we, Km, vm, nx, Kq, Kp, vt, wr, wq, wp, hq, nr, ns = params
    cr, em, cp, cq, ct, et, cm, zmr, zmp, zmq, zmt, zmm, mt, mm, q, p, si, mq, mp, mr, r, a = initial

    dcr, dem, dcp, dcq, dct, det, dcm, dzmr, dzmp, dzmq, dzmt, dzmm, dmt, dmm, dq, dp, dsi, dmq, dmp, dmr, dr, da = zeros(length(dydt))

    Kgamma = gmax/Kp
    gamma = gmax*a/(Kgamma+a)
    ttrate = (cq+cr+cp+ct+cm)gamma
    lam = ttrate/M
    vimp = (et*vt*s0/(Kt+s0))
    nucat = (em*vm*si/(Km+si))

    dydt[1] =  +r*mr*kb - cr*ku - cr*lam - gamma/nr*cr - f*cr
    dydt[2] = - lam*em + cm*gamma/nx
    dydt[3] =   +r*mp*kb - cp*ku - cp*lam - gamma/nx*cp - f*cp
    dydt[4] =  +r*mq*kb - cq*ku - cq*lam - gamma/nx*cq - f*cq
    dydt[5] =  +r*mt*kb - ct*ku - ct*lam - gamma/nx*ct - f*ct
    dydt[6] =  - lam*et + gamma/nx*ct
    dydt[7] =   +r*mm*kb - cm*ku - cm*lam - gamma/nx*cm - f*cm
    dydt[8] = +f*cr - b*zmr - lam*zmr
    dydt[9] = +f*cp - b*zmp - lam*zmp
    dydt[10] = +f*cq - b*zmq - lam*zmq
    dydt[11] = +f*ct - b*zmt - lam*zmt
    dydt[12] = +f*cm - b*zmm - lam*zmm
    dydt[13] =   +(we*a/(thetax+a)) - mt*dm - mt*lam - r*mt*kb + ct*ku + gamma/nx*ct
    dydt[14] =   +(we*a/(thetax+a)) - mm*dm - mm*lam - r*mm*kb + cm*ku + gamma/nx*cm
    dydt[15] = - lam*q + gamma/nx*cq
    dydt[16] =  - lam*p + gamma/nx*cp
    dydt[17] = - lam*si + vimp - nucat
    dydt[18] =   +(wq*a/(thetax+a)/(1+(q/Kq)^hq)) - mq*dm - mq*lam - r*mq*kb + cq*ku + gamma/nx*cq
    dydt[19] = +(wp*a/(thetax+a)) - mp*dm - mp*lam - r*mp*kb + cp*ku + gamma/nx*cp
    dydt[20] =   +(wr*a/(thetar+a)) - mr*dm - mr*lam - r*mr*kb + cr*ku + gamma/nr*cr
    dydt[21] =  - lam*r - r*mr*kb - r*mt*kb - r*mm*kb - r*mq*kb - r*mp*kb + cr*ku + ct*ku + cm*ku + cq*ku + cp*ku + gamma/nr*cr + gamma/nr*cr + gamma/nx*ct + gamma/nx*cm + gamma/nx*cq + gamma/nx*cp
    dydt[22] =  +ns*nucat - ttrate - lam*a

end

