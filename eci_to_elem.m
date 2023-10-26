function [a,e, i, RAAN, w, f] = eci_to_elem(r_vec,v_vec, u)
    rmag = norm(r_vec);
    vmag = norm(v_vec);

    h_vec = cross(r_vec,v_vec);
    e_vec = (1/u)*((vmag^2-(u/rmag))*r_vec - dot(r_vec,v_vec)*v_vec);
    n_vec = cross([0,0,1],h_vec);

    h = norm(h_vec);
    n = norm(n_vec);

    e = norm(e_vec);

    a = (h^2/u)/(1-e^2);

    i = acos(h_vec(3)/h);

    cos_RAAN = n_vec(1)/n;
    cos_argP = dot(n_vec,e_vec)/(n*e);   
    cos_v0 = dot(e_vec,r_vec)/(e*rmag);

    RAAN = acos(cos_RAAN);
    if(n_vec(2)<0)
        RAAN = 2*pi - RAAN;
    end

    w = acos(cos_argP);

    if(e_vec(3)<0)
        w = 2*pi - w;
    end

    f =acos(cos_v0);
    
    if(dot(r_vec,v_vec)<0)
        f = 2*pi - f;
    end
end