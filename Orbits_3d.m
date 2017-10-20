function Orbits_3d(name, time, pos, color, export_filepath, is_plot_visible)
%ORBITS_3D plots an orbit.
    load coast;
    load('topo.mat','topo','topomap1');

    f = figure('Name','Orbits');
    set(f, 'Visible', is_plot_visible);

    grid;

    xlabel('x-eci, km');
    ylabel('y-eci, km');
    zlabel('z-eci, km');
    axis equal;
    grid;

    npanels = 72;   % Number of globe panels around the equator deg/panel = 360/npanels
    alpha   = .5; % globe transparency level, 1 = opaque, through 0 = invisible
    image_file = 'land_ocean_ice_1024.jpg';
    erad    = 6371.0087714; % equatorial radius (km)
    prad    = 6371.0087714; % polar radius (km)
    [x, y, z] = ellipsoid(0, 0, 0, erad, erad, prad, npanels);
    globe = surf(x, y, -z);
    cdata = imread(image_file);
    set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', alpha,...
        'EdgeColor', 'none', ...
        'AmbientStrength', 0.5);

    hAnnotation = get(globe,'Annotation'); hLegendEntry = get(hAnnotation','LegendInformation'); set(hLegendEntry,'IconDisplayStyle','off');
    axis equal;
    grid;

    hold on
         % create plot  
         scatter3(pos(:,1),pos(:,2),pos(:,3), 5, 'MarkerEdgeColor', color, 'MarkerFaceColor', color);   
    hold off

    if ~strcmp(export_filepath, '') 
        print(export_filepath, '-dpng')
    end
end