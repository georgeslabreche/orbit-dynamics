function [raan, incl, argp, a, e, anom_v] = r_v_to_param1(r,V)
% R_V_TO_PARAM Calculates the parameters of the orbit for the given values of r and V.

    %%%%%%%%%%%%%%%%%%%r_V_to_parameters file%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    global MUe
    global c1 c2 c3

    % Inclination.
    
    h = cross(r,V);
    mod_h = norm(h);
    vers_h = h/mod_h;
    incl = acos(dot(vers_h,c3));
    s_incl = sin(incl);

    if s_incl >= 1.e-10
        vect_N = cross(c3, vers_h);
        vers_N = vect_N / norm(vect_N);
        c_raan = dot(vers_N, c1);
        s_raan = dot(vers_N, c2);
        raan = atan2(s_raan, c_raan);

        if raan < 0
            raan = (2*pi) + raan;
        end

    else
        raan = 0;
        vers_N = c1;

    end

    % Eccentricity (0>e>1).

    vers_r = r / norm(r);
    vect_e = (cross(V,h)) / MUe - vers_r;
    e = norm(vect_e);

    if e>=1.e-10
        vers_e = vect_e / norm(vect_e);
    else
        vers_e = vers_N;

    end

    % Argument of perigee.

    vect_M = cross(vers_h, vers_N);
    vers_M = (cross(vers_h, vers_N)/norm(vect_M));
    c_argp = dot(vers_e,vers_N);
    s_argp = dot(vers_e,vers_M);
    argp = atan2(s_argp, c_argp);

    if argp < 0
        argp = (2*pi) + argp;
    end
    
    % Semi-major axis.

    mod_V = norm(V);
    mod_r = norm(r);
    En = (mod_V^2/2 - (MUe/mod_r)); % E is energy
    a = -MUe/(2*En);

    % True anomaly (0>anom_v>2*pi).

    mod_p = cross(vers_h, vers_e);
    vers_p = (cross(vers_h, vers_e)/norm(mod_p));
    c_anom_v = dot(vers_r,vers_e);
    s_anom_v = dot(vers_r,vers_p);
    anom_v = atan2(s_anom_v, c_anom_v);

    % Return calculated values.
    Params = [raan, incl, argp, a, e, anom_v];
end