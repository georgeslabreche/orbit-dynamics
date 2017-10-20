function [ position ] = eul_to_position(raan, incl, argp, xyz, Vxyz)
%EUL_TO_ROTMAT converts Euler Angles to a rotation matrix.

    global MUe 

    % First row of the rotation matrix.
    m1 = cos(raan) * cos(argp) - cos(incl) * sin(argp) * sin(raan);
    m2 = cos(argp) * sin(raan) + cos(incl) * cos(raan) * sin(argp);
    m3 = sin(incl) * sin(argp);

    % Second row of the rotation matrix.
    m4 = -sin(argp) * cos(raan) - cos(incl) * sin(raan) * cos(argp);
    m5 = -sin(raan) * sin(argp) + cos(incl) * cos(argp) * cos(raan);
    m6 = sin(incl) * cos(argp);

    % Third row of the rotation matrix.
    m7 = sin(raan) * sin(incl);
    m8 = -cos(raan) * sin(incl);
    m9 = cos(incl);

    RR = [m1, m2, m3;...
        m4, m5, m6;...
        m7, m8, m9];

    r = inv(RR) * xyz;
    V = inv(RR) * Vxyz;

    position = [r,V];
end


