%%-----------------------------------------------------------------------%%
%         R0008R Introduction to Space Mechanics and Electronics          %
%                    Lulea University of Technology                       %
%                      Teacher: Leonard Felicetti                         %
%------------------------------------------------------------------------%%
%                      Assignment 2: Orbit Dynamics                       %
%%-----------------------------------------------------------------------%%
%             Students: Georges Labreche and Natalie Lawton               %
%%-----------------------------------------------------------------------%%

% Clean up environment before starting.
clear all 
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               INIT              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The program automaticcaly exports all generated plots into PNG files.
% Define the plot export directories...
export_dir_root = 'img'; % Root export directory.
export_dir_orbits = strcat(export_dir_root, '/orbits'); % Orbit plots export directory.
export_dir_gt = strcat(export_dir_root, '/ground-tracks'); % Ground tracks export directory.
export_dir_ea = strcat(export_dir_root, '/elevation-azimuth'); % Elevation and azimuth plots export directory.
export_dir_ht = strcat(export_dir_root, '/hohmann-transfer'); % Hohmann Transfer plots export directory.

export_dirs = {export_dir_orbits, export_dir_gt, export_dir_ea, export_dir_ht};

% If an image export directory doesn't exist, create it.
for idx = 1:length(export_dirs)
    fn = fullfile(export_dirs{idx});
    if ~exist(fn, 'dir')
       mkdir(export_dirs{idx});
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               DIALOG              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct a question dialog asking the user if we should display the
% plots.
% Since there are a lot of plots we may just want to run the program
% without displaying all of the plots and just look at the exported images.
choice = questdlg(sprintf('Generated plots will be exported as PNG files in the project directory.\nWould you like to display them while they are being generated?'), ...
	'Options', ...
	'Yes','No', 'Cancel', 'Cancel');

% Handle dialog box response.
switch choice
    case 'Yes'
        plot_visible_orbit = 'on';
        plot_visible_gt = 'on';
        plot_visible_ea = 'on';
        plot_visible_op = 'on';
        plot_visible_ht = 'on';
        plot_visible_htop = 'on';
    case 'No'
        plot_visible_orbit = 'off';
        plot_visible_gt = 'off';
        plot_visible_ea = 'off';
        plot_visible_op = 'off';
        plot_visible_ht = 'off';
        plot_visible_htop = 'off';
    case 'Cancel'
        disp('Program canceled.');
        return;
    otherwise
        disp('Program canceled.');
        return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               INIT             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Global variables.
global MUe
global c1 c2 c3
global We Long_G0

% Choose how many orbits we want to simulate
NO = 10; % CHANGE NUMBER OF ORBITS HERE.

MUe = 398600.4;
c1 = [1 0 0];
c2 = [0 1 0];
c3 = [0 0 1];
We = (2*pi/86164);
Long_G0 = 0;
Rs = 6371; % Earth radius in KM.

gso_a = (MUe/We^2)^(1/3); % 42164 KM.
task4_a = 300+Rs;

% Loop through all orbits that need to be processed
for orbit = 1:7
    % Property matrices for LEO, GTO, GEO, Molniya, Tundra, MEO, and
    % task4LEO.
    % Putting their initalization inside the loop because unfortunately
    % some of these variables are overwritten further down the code so we 
    % need to re-init them at the beginning of every loop iteration.
    raan = [30 10 30 180 195 145.5 10]; %deg
    incl = [47.2 7 10 63.4 43 53.7 0.01]; %deg
    argp = [45 5 325 325 270 52.4 270]; %deg
    a    = [7978 24400 gso_a 26555 42164 26560 task4_a]; %km
    e    = [0.05 0.73 0.2 0.72 0.075 0.01 0.001];
    M    = [70 250 300 90 305 -37.3 0]; %deg

    % Determine orbit name and color based on loop index.
    orbit_name = '';
    color = '';

    switch orbit
       case 1 % LEO
          orbit_name = 'leo';
          color = [1 0 0]; % red
       case 2 % GTO
          orbit_name = 'gto';
          color = [0 1 0]; % green
       case 3 % GEO
          orbit_name = 'geo';
          color = [0 1 1]; % cyan
       case 4 % Molniya
          orbit_name = 'molniya';
          color = [1 0.5 0]; % orange
       case 5 % Tundra
          orbit_name = 'tundra';
          color = [1 0 1]; % magenta
       case 6 % MEO
          orbit_name = 'meo';
          color = [0 0 0]; % black
       case 7 % task4LEO
          orbit_name = 't4-meo';
          color = [0 0 1]; % blue
    end
    
    % Output feedback in Terminal.
    disp([' ']);
    disp(['Computing ', upper(orbit_name), ' orbit.']);
    
    raan = raan(orbit)*pi/180; %rad
    incl = incl(orbit)*pi/180; %rad
    argp = argp(orbit)*pi/180; %rad
    a = a(orbit);
    e = e(orbit);
    M = M(orbit)*pi/180; %rad
    
    time_0 = 0;
    time_F = NO*2*pi*sqrt(a^3/MUe);
    d_time = 10;
    tspan  = [time_0:d_time:time_F];

    % M = E-esinE
    % Iterative function to find the eccentric anomaly from the mean
    % anomaly.
    E0 = M;
    diff = 1;
    tol=(1.e-10);
    while abs(diff) >  tol
        E_new = (E0-e*sin(E0))-M;
        E_new1 = (1-e*cos(E0));
        diff = E_new/E_new1;
        if  abs(diff) >  tol
          E0 = E0 - diff;
        else
          E=E0;
        end
    end

    % Find the true anomaly.
    c_v = sqrt(1-e)*cos(E/2);
    s_v = sqrt(1+e)*sin(E/2);
    anom_v = atan2(s_v,c_v);

    % Find the cartesian coordinates using Eulers and the eccentric
    % anomaly.
    xyz = a*[cos(E)-e;sqrt(1-e^2)*sin(E);0];

    vers_xyz = xyz/norm(xyz);

    % Find the velocity vector.
    p = a*(1-e^2);
    vect_p = [0;1;0];
    h = sqrt(MUe*p);
    vect_anom_v = [-vers_xyz(2); vers_xyz(1); 0];
    Vxyz = (MUe/h)*(e*vect_p+vect_anom_v);
    mod_Vxyz = norm(Vxyz); % To check this is as expected.

    position = eul_to_position(raan, incl, argp, xyz, Vxyz);

    % Initial position vectors.
    r_0 = position(:,1);
    V_0 = position(:,2);

    % Integrate to find position and velocity.
    op = odeset('RelTol',1.e-10,'AbsTol',1.e-10);
    X_0 = [r_0' V_0'];

    [time, X]=ode45('orb_dynamics', tspan, X_0, op);

    pos = X(:,1:3);
    vel = X(:,4:6);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                 TASK 1: ORBITS                  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot the orbit and export it as PNG file. 
    export_filepath_orbit = strcat(export_dir_orbits, '/', 'orbit-', orbit_name);
    Orbits_3d('Orbit',time, pos, color, export_filepath_orbit, plot_visible_orbit)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %               TASK 2: GROUNDTRACKS              %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Output feedback in Terminal.
    disp(['Computing ', upper(orbit_name), ' ground track.']);

    % Set up zero matrices.
    Lat = zeros(1,length(time)) ; 
    Long = zeros(1,length(time));
    Long_G = zeros(1,length(time));


    for t_time = 1:length(time)
        Long_G(t_time) = Long_G0 + 180/pi*We*(time(t_time)-time_0); %rotation of earth
        Long(t_time) = 180/pi*atan2(pos(t_time, 2),pos(t_time, 1))-Long_G(t_time); % Satellite ground longitude
        Lat(t_time) = 180/pi*asin(pos(t_time,3)/norm(pos(t_time,:))); % Satellite ground latitude
        Long(t_time) = wrapTo180(Long(t_time));
    end

    gpos = [Lat', Long'];

    export_filepath_gt = strcat(export_dir_gt, '/', 'orbit-gt-', orbit_name);
    Earth_Map(gpos, time, color, export_filepath_gt, plot_visible_gt);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   TASK 3 VISIBILITY: ELEVATION AND AZIMUTH    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Coordinated for Kiruna, Malindi, Cape Canaveral, Tanegashima
    placet = [67 2 28 30]; 
    placeg = [20 40 -80 130];

    % Go through all locations to view the satellite from.
    for view = 1:length(placet)

        % Ground Station name.
        gs_name = '';

        % Get the name of the Ground Station from which we will plot visibility.
        switch view
            case 1 % Kiruna
                gs_name = 'kiruna';
            case 2 % Malindi
                gs_name = 'malindi';
            case 3 % Cape Carnaveral
                gs_name = 'cape-canaveral';
            case 4 % Tanegashima
                gs_name = 'tanegashima';
        end
        
        % Output feedback in Terminal.
        disp(['   - Computing elevation and azimuth view from ', upper(gs_name), '.']);

        RR_e = zeros(length(time),3);
        vect_R = zeros(length(time),3);
        vect_d = zeros(length(time),3);
        dir_d = zeros(length(time),3);
        RR_2G = zeros(3,3,length(time));
        
        elev = NaN(length(time),3);
        azi = NaN(length(time),3);

        Lat_G = placet(view)*pi/180; 
        Long_S = placeg(view)*pi/180;

        for t_time = 1:length(time)

            Long_G(t_time) = wrapTo180(Long_G0 + 180/pi*We*(time(t_time)-time_0)); 
            RR_e(t_time,:) = earth_rot(Lat_G, pi/180*Long_G(t_time), Long_S)';
            vect_R(t_time,:) = (Rs*RR_e(t_time,:))';
            vect_d(t_time,:) = pos(t_time,:) - vect_R(t_time,:);
            dir_d(t_time,:) = vect_d(t_time,:)/norm(vect_d(t_time,:));

            % Rotation matrices in three component parts North, Up and East
            N = [-sin(Lat_G)*cos(pi/180*Long_G(t_time)+Long_S);...
                -sin(Lat_G)*sin(pi/180*Long_G(t_time)+Long_S);...
                cos(Lat_G)];

            U = RR_e(t_time,:)';

            E = cross(N,U);

            % now find elevation and azimuth angles
            elev(t_time) = pi/2-acos(dot(dir_d(t_time,:),U));
            c_azi = dot(dir_d(t_time,:),N);
            s_azi = dot(dir_d(t_time,:),E);
            azi(t_time) = atan2(s_azi,c_azi);

            %Prevent graph from plotting when elevation is below horizon
            if elev(t_time) < 0
                elev(t_time)= NaN;
                azi(t_time) = NaN;
            end
        end

        % Plot elevation and azimuth.
        eaf = figure;
        set(eaf, 'Visible', plot_visible_ea);
        subplot(2,1,1);
        plot(time/60/60, elev*180/pi,...
            'LineWidth', 2,...
            'Color', color);
        title ('Elevation');
        ylabel('Angle (deg)');
        xlabel('Time (hours)');

        subplot(2,1,2);
        plot(time/60/60, azi*180/pi,...
            'LineWidth', 2,...
            'Color', color);
        title ('Azimuth');
        ylabel('Angle (deg)');
        xlabel('Time (hours)');

        export_filepath_ea = strcat(export_dir_ea, '/', orbit_name, '-ea-', gs_name);
        print(export_filepath_ea, '-dpng');
    end    

end % End loop that plots all orbits. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               TASK 4 PART A               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculating Hohmann transfer orbits.

h = 300; % Altitude of initial orbit.
e_in = 0.001; % Eccentricity of initial orbit.
gso_a = 42164; % Semi-major axis of final.
e_f = 0; % Eccentricity of final.
incl_f = pi/180*10; % Inclination of final.

r_in = (Rs + h)*(1-e_in); % Radius in first orbit.
V_in = sqrt(MUe/r_in); % Velcoity in first orbit.
En_in = -MUe/(2*r_in); % Energy of first orbit.

r_f = gso_a*(1-e_f);
V_f = sqrt(MUe/r_f);
En_f = -MUe/(2*r_f); % Energy of final orbit.

r_hohm = (r_in + r_f)/2; % Radius of hohmann.
En_hohm = -MUe/(2*r_hohm); % Energy of hohmann.
V_hohmP = sqrt(MUe*((2/r_in)-(1/r_hohm))); % Veloceity as perigee.
V_hohmA = sqrt(MUe*((2/r_f)-(1/r_hohm))); % Velocity at apogee.

hohm_am = r_in*V_hohmP; % Specific angular momentum trasfer orbit.
e_hohm = sqrt(1+((2*En_hohm*hohm_am^2)/MUe^2)); % ccentricity of hohmann transfer orbit.

del_V_in = V_hohmP - V_in;
del_V_f = V_f - V_hohmA;
del_V = norm(del_V_in)+norm(del_V_f);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               TASK 4 PART B               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m_in = 2500; %kg
Isp = 300; % s
g0 = 9.81; % gravity of earth

mp_1 = m_in*(1-exp((del_V_in*10^3)/(-Isp*g0))); % mass required at first change
mp_2 = (m_in-mp_1)*(1-exp((del_V_f*10^3)/(-Isp*g0))); % mass required at second change
m_T = mp_1 + mp_2; %total mass required

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               TASK 4 PART C               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
del_V_CA = sqrt(((V_hohmA^2)+(V_f^2))-(2*V_hohmA*V_f*cos(incl_f))); % change in velocity from GEO
del_V_CP = sqrt(((V_in^2)+(V_hohmP^2))-(2*V_in*V_hohmP*cos(incl_f))); % change in velocity from LEO

del_VTA = del_V_in + del_V_CA; % total change in velocity required combining orbital plane and height change at apogee
del_VTP = del_V_f + del_V_CP; % total change in velocity required combining orbital plane and height change at perigee

% Plotting these results
raan4 = [10 10 10]; %deg
argp4 = [270 0 0]; %deg
a4    = [r_in r_hohm r_f]; %km
e4    = [e_in e_hohm e_f];
M4    = [0 0 0]; %deg

% CAN CHANGE INCLINATION FOR PART A AND C HERE.
incl4a = [0 0 0]; %deg 
incl4c = [10 10 0]; %deg

% Can include new inclinations.
inclinations = {incl4a, incl4c};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               PLOT FOR TASK 4 PART A AND C                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For each defined inclination, compute the Hohmann Transfer and its
% orbital parameters
for idx = 1:length(inclinations)

    % Get inclination.
    incl4 = inclinations{idx};
    
    disp([' ']);
    disp('Computing Hohmann Transfer with inclination:');
    disp(incl4);

    % Empty orbit plot. Will be updated with Homann Transfer orbits...
    pos=NaN(1,3);
    Orbits_3d('Hohmann Transfer', 0, pos, [1 0 0], '', plot_visible_ht);

    time4 =[];
    pos5 =[];
    vel5 =[];

    for k=1:3
        raan = raan4(k)*pi/180; %rad
        incl = incl4(k)*pi/180; %rad
        argp = argp4(k)*pi/180; %rad
        a = a4(k);
        e = e4(k);
        M = M4(k)*pi/180; %rad

        % Choose how many need to be simulated.
        if k==1
            NO = 2;

        elseif k==3
            NO = 2;

        else
            NO = 0.5;
        end

        d_time = 50;

        if k ==1
            time_0 = 0;
            time_F = NO*2*pi*sqrt(a^3/MUe);
            tspan  = [time_0:d_time:time_F];
        elseif k==2
            time_0 = time4(end);
            time_F = time_0 + NO*2*pi*sqrt(a^3/MUe);
            tspan  = [time_0:d_time:time_F];
        elseif k==3
            time_0 = time4(end);
            time_F = time_0 + NO*2*pi*sqrt(a^3/MUe);
            tspan  = [time_0:d_time:time_F];
        end

        % M = E-esinE
        % Iterative function to find the eccentric anomaly from the mean
        % anomaly.
        E0 = M;
        diff = 1;
        tol=(1.e-10);
        while abs(diff) >  tol
            E_new = (E0-e*sin(E0))-M;
            E_new1 = (1-e*cos(E0));
            diff = E_new/E_new1;
            if  abs(diff) >  tol
              E0 = E0 - diff;
            else
              E=E0;
            end
        end

        % Find the true anomaly.
        c_v = sqrt(1-e)*cos(E/2);
        s_v = sqrt(1+e)*sin(E/2);
        anom_v = atan2(s_v,c_v);

        % Find the cartesian coordinates using Eulers and the eccentric
        % anomaly.
        xyz = a*[cos(E)-e;sqrt(1-e^2)*sin(E);0];

        vers_xyz = xyz/norm(xyz);

        % find the velocity vector.
        p = a*(1-e^2);
        vect_p = [0;1;0];
        h = sqrt(MUe*p);
        vect_anom_v = [-vers_xyz(2); vers_xyz(1); 0];
        Vxyz = (MUe/h)*(e*vect_p+vect_anom_v);
        mod_Vxyz = norm(Vxyz); % to check this is as expected.

        % Put into fucntion to make task 4 nicer.
        position = eul_to_position(raan, incl, argp, xyz, Vxyz);

        % Initial position vectors.
        r_0 = position(:,1);
        V_0 = position(:,2);

        % Integrate to find position and velocity.
        op = odeset('RelTol',1.e-10,'AbsTol',1.e-10);
        X_0 = [r_0' V_0'];

        [time, X]=ode45('orb_dynamics', tspan, X_0, op);

        time4 = [time4;time];

        pos4 = X(:,1:3);
        vel4 = X(:,4:6);

        pos5 = [pos5;X(:,1:3)];
        vel5 = [vel5;X(:,4:6)];


        hold on
        if k == 1
            scatter3(pos4(:,1),pos4(:,2),pos4(:,3), 'MarkerEdgeColor', [0 1 0], 'MarkerFaceColor', [0 1 0]);
        end
        if k ==2
            scatter3(pos4(:,1),pos4(:,2),pos4(:,3), 'MarkerEdgeColor', [1 0 0], 'MarkerFaceColor', [1 0 0]);
        end
        if k == 3
            scatter3(pos4(:,1),pos4(:,2),pos4(:,3), 'MarkerEdgeColor', [0 0 1], 'MarkerFaceColor', [0 0 1]);
        end
    end
    
    % Print resulting homann transfer pot as PNG export.
    export_filepath_ht = strcat(export_dir_ht, '/hohmann-transfer-with-inclination-', int2str(incl4(1)), '-', int2str(incl4(2)), '-', int2str(incl4(3)));
    print(export_filepath_ht, '-dpng')

    % TO FIND THE ORBITAL PARAMETERS FROM THE POSITION AND VELOCITY VECTORS
    disp('Computing Hohmann Transfer orbital parameters with inclination:');
    disp(incl4);

    OP = zeros(length(time4),6);
    for t_time = 1:length(time4)
        [OP(t_time,1),OP(t_time,2),OP(t_time,3),OP(t_time,4),OP(t_time,5),OP(t_time,6)] = r_v_to_param1(pos5(t_time,:), vel5(t_time,:));
    end
    
    time = time4;

    % Hohmann Transfer orbital parameter plots.
    fop = figure('Name', 'Orbital Parameters');
    set(fop, 'Visible', plot_visible_htop);

    subplot(2,3,1);
    plot(time/3600,OP(:,1)*180/pi);
    xlabel('time, hours');
    ylabel('raan, deg');
    grid;
    axis([0 time(end)/3600 0 360]);

    subplot(2,3,2);
    plot(time/3600,OP(:,2)*180/pi);
    xlabel('time, hours');
    ylabel('incl, deg');
    grid;
    axis([0 time(end)/3600 0 360]);

    subplot(2,3,3);
    plot(time/3600,OP(:,3)*180/pi);
    xlabel('time, hours');
    ylabel('argp, deg');
    grid;
    axis([0 time(end)/3600 0 360]);

    subplot(2,3,4);
    plot(time/3600,OP(:,4));
    xlabel('time, hours');
    ylabel('a, km');
    grid;

    subplot(2,3,5);
    plot(time/3600,OP(:,5));
    xlabel('time, hours');
    ylabel('e');
    grid;

    subplot(2,3,6);
    plot(time/3600,OP(:,6)*180/pi);
    xlabel('time, hours');
    ylabel('anom_v, deg');
    grid;
    axis([0 time(end)/3600 -181 181]);

    export_filepath_htop = strcat(export_dir_ht, '/ht-orbital-parameters-with-inclination-', int2str(incl4(1)), '-', int2str(incl4(2)), '-', int2str(incl4(3)));
    print(export_filepath_htop, '-dpng')

end

% Finished program. % 
disp('All done. Results have been exported to the image directory.');
