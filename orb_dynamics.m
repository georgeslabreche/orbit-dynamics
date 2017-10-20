function dX_dt = orb_dynamics(t,X)
%ORB_DYNAMICS Calculates the derivatives of the position and veloctiy with respect to time.
    global MUe

    r = X(1:3,1);
    V = X(4:6,1);
    mod_r = norm(r);

    dV_dt = -MUe/mod_r^3 * r;

    dr_dt = V;


    dX_dt = [dr_dt;dV_dt];
end