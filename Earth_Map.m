function Earth_Map(gpos, time, color, export_filepath, is_plot_visible)
% EARTH_MAP Plots an orbit's ground tracks.

    global We Long_G0

    load('topo.mat');
    load coast;
    load('topo.mat','topo','topomap1');

    topo2=zeros(size(topo));
    for i_topo = 1:180
        topo2(:,i_topo) = topo(:,180+i_topo);
        topo2(:,180+i_topo) = topo(:,i_topo);
    end
    f = figure('Name','Ground Track');
    set(f, 'Visible', is_plot_visible);
    contour(-180:179,-89:90,topo2,[0 0],'b');
    ylabel('Latitude (deg)');
    xlabel('Longitude (deg)');
    %axis([-180 180 -90 90])

    grid
    hold on
        % Geolocate the groundstations.
        scatter(20, 67, 'k', 'filled'); % Kiruna
        text(20, 67, '  Kiruna' , 'Color', 'k');

        scatter(40, 2, 'k', 'filled'); % Malindi
        text(40, 2, '  Malindi' , 'Color', 'k');

        scatter(-80, 28, 'k', 'filled'); % Cape Carnaveral
        text(-80, 28, '  Cape Canaveral' , 'Color', 'k');

        scatter(130, 30, 'k', 'filled'); % Tanegashima
        text(130, 30, '  Tanegashima' , 'Color', 'k');

        % Draw the ground tracks.
        scatter(gpos(:,2), gpos(:,1), 1, 'MarkerEdgeColor', color, 'MarkerFaceColor', color);
    hold off

    % Export the plot as a PNG image.
    print(export_filepath, '-dpng')
end

