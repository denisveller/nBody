
function [r_vec,v_vec] = elem_to_eci(a,e,i,RAAN,w,f,u)
    p = a*(1-e^2);
    [Rw,Vw] = Elements_to_Perifocal(p,e,f,u);
    
    [r_vec,v_vec] = Perifocal_to_ECI (Rw,Vw,i,RAAN,w);
end
function [Rw,Vw] = Elements_to_Perifocal (p,e, f,u)
    rmag = p/(1+e*cos(f));
    Rw = [rmag*cos(f),rmag*sin(f),0];
    Vw = sqrt(u/p)*[-1*sin(f),e+cos(f),0];
end
function [r_vec,v_vec] = Perifocal_to_ECI (Rw,Vw, i, RAAN, w)
    R11 = cos(RAAN)*cos(w) - sin(RAAN)*sin(w)*cos(i);
    R12 = -cos(RAAN)*sin(w) - sin(RAAN)*cos(w)*cos(i);
    R13 = sin(RAAN)*sin(i);
    R21 = sin(RAAN)*cos(w)+cos(RAAN)*sin(w)*cos(i);
    R22 = -sin(RAAN)*sin(w)+cos(RAAN)*cos(w)*cos(i);
    R23 = -cos(RAAN)*sin(i);
    R31 = sin(w)*sin(i);
    R32 = cos(w)*sin(i);
    R33 = cos(i);
    tMat = [R11,R12,R13;
            R21,R22,R23;
            R31,R32,R33];
    Rw = Rw.';
    Vw = Vw.';

    r_vec = tMat*Rw;
    v_vec = tMat*Vw;
    r_vec = r_vec.';
    v_vec = v_vec.';

end