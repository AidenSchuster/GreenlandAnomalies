% Load Data
cd('C:\Users\ajs82292\Desktop\Research\Matlab\Source\Greenland_Melt') ;
addpath('seawater','C:\Users\ajs82292\Desktop\Research\Matlab\Source\seawater') ;
load("02cleanNODC_updated2.mat")
load('OMG_data.mat')
load("x_coast.mat")
load("y_coast.mat")
s = settings;
s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = "light";
% Commonly Changed Variables for box size and day interval (set run to one
% below if changing rectangle size)
recwidth = 10*2 ; % km doubled as it is 10 "radius" Change value for different sized cast boxes
reclength = 20*2 ; % km doubled as it is 20 "radius" Change value for different sized cast boxes
day_range = 15 ; % adjust this as needed, actual day range is 2*number listed
radius = 20 ; % km circle for opean ocean
min_count = 4 ; % the minimum # of data points needed in order to plot the statistics of a cast (std dev, anomalies ect.)
aspect_ratio = cosd(65) ; % Aspect Ratio at 65 N
target  = 100 ; % km from coast
%combine OMG_data into NODC data
length_NODC = length(lat) ;
lat = [lat,OMG_lat] ;
lon = [lon,OMG_lon] ;
dep = [dep,OMG_depth] ;
mon = [mon,OMG_mon] ;
yea = [yea,OMG_yea] ;
sal = [sal,OMG_sal] ;
temp = [temp,OMG_temp] ;
day = [day,OMG_day] ;
clear OMG_temp OMG_sal OMG_day OMG_mon OMG_yea OMG_lat OMG_lon OMG_depth
% extend coastline to -30 using ETOPO depth data and 0-5 countour (ice surface 15 arcseconds)
[Z_eto,R] = readgeoraster('extendedcoast.tiff') ;
[rows,cols] = size(Z_eto) ;
[rowGrid, colGrid] = ndgrid(1:rows, 1:cols);
[YY_eto, XX_eto] = intrinsicToGeographic(R, colGrid, rowGrid);
Z_eto = double(Z_eto) ;
coast_range = (Z_eto>=0) & (Z_eto<=5) ; % 0-5 m coastline
Z_masked = Z_eto ;
depth_min = 0 ;
depth_max = 5 ;
run = 2 ;
for i = 1:1:1
    if run == 1
figure;
[C,h] = contour(XX_eto, YY_eto, Z_masked,[depth_min,depth_max]) ;
xlabel('Longitude');
ylabel('Latitude');
%C(wholeNumberIndex) = NaN ;
zeroidx = (C == 0) ;
for i = 1:1:length(C)
    if C(1,i) == 0
C(:,i) = NaN ;
    end
end
indices = find(C > 0 & C < 10,1);
[row, col] = ind2sub(size(C), indices);
C = C(:,1:col-1) ;
cx = [cx;C(1,:)'] ;
cy = [cy;C(2,:)'] ;
save 'cx-cy.mat' cx cy
    end
end
load 'cx-cy.mat' cx cy
clear rowGrid colGrid cols R te rows depth_min depth_max C Z_masked coast_range h indices zeroidx
% Extend even furth east (0)
[Z_eto_E,R] = readgeoraster('extendedEcoast.tiff') ;
[rows,cols] = size(Z_eto_E) ;
[rowGrid, colGrid] = ndgrid(1:rows, 1:cols);
[YY_eto_E, XX_eto_E] = intrinsicToGeographic(R, colGrid, rowGrid);
Z_eto_E = double(Z_eto_E) ;
coast_range = (Z_eto_E>=0) & (Z_eto_E<=5) ; % 0-5 m coastline
Z_masked = Z_eto_E ;
depth_min = 0 ;
depth_max = 5 ;
run = 2 ;
for i = 1:1:1
    if run == 1
figure;
[C,h] = contour(XX_eto_E, YY_eto_E, Z_masked,[depth_min,depth_max]) ;
xlabel('Longitude');
ylabel('Latitude');
%C(wholeNumberIndex) = NaN ;
zeroidx = (C == 0) ;
for i = 1:1:length(C)
    if C(1,i) == 0
C(:,i) = NaN ;
    end
end
indices = find(C > 0 & C < 10,1);
[row, col] = ind2sub(size(C), indices);
C = C(:,1:col-1) ;
cx = [cx;C(1,:)'] ;
cy = [cy;C(2,:)'] ;
save 'cx-cy.mat' cx cy
    end
end
load 'cx-cy.mat' cx cy
clear rowGrid colGrid cols R te rows depth_min depth_max C Z_masked coast_range h indices zeroidx
% Extend it further north using same ETOPO dataset and 0-5 contour
[Z_eto_N,R] = readgeoraster('extendedNcoast.tiff') ;
[rows,cols] = size(Z_eto_N) ;
[rowGrid, colGrid] = ndgrid(1:rows, 1:cols);
[YY_eto_N, XX_eto_N] = intrinsicToGeographic(R, colGrid, rowGrid);
Z_eto_N = double(Z_eto_N) ;
coast_range = (Z_eto_N>=0) & (Z_eto_N<=5) ; % 0-5 m coastline
Z_masked = Z_eto_N ;
depth_min = 0 ;
depth_max = 5 ;
run = 2 ;
for i = 1:1:1
    if run == 1
figure;
[C,h] = contour(XX_eto_N, YY_eto_N, Z_masked,[depth_min,depth_max]) ;
xlabel('Longitude');
ylabel('Latitude');
%C(wholeNumberIndex) = NaN ;
zeroidx = (C == 0) ;
for i = 1:1:length(C)
    if C(1,i) == 0
C(:,i) = NaN ;
    end
end
indices = find(C > 0 & C < 10,1);
[row, col] = ind2sub(size(C), indices);
C = C(:,1:col-1) ;
cx = [cx;C(1,:)'] ;
cy = [cy;C(2,:)'] ;
save 'cx-cy.mat' cx cy
    end
end
load 'cx-cy.mat' cx cy
%Extend even furth north
[Z_eto_NN,R] = readgeoraster('extendedNNcoast.tiff') ;
[rows,cols] = size(Z_eto_NN) ;
[rowGrid, colGrid] = ndgrid(1:rows, 1:cols);
[YY_eto_NN, XX_eto_NN] = intrinsicToGeographic(R, colGrid, rowGrid);
Z_eto_NN = double(Z_eto_NN) ;
coast_range = (Z_eto_NN>=0) & (Z_eto_NN<=5) ; % 0-5 m coastline
Z_masked = Z_eto_NN ;
depth_min = 0 ;
depth_max = 5 ;
run = 2 ;
for i = 1:1:1
    if run == 1
figure;
[C,h] = contour(XX_eto_NN, YY_eto_NN, Z_masked,[depth_min,depth_max]) ;
xlabel('Longitude');
ylabel('Latitude');
%C(wholeNumberIndex) = NaN ;
zeroidx = (C == 0) ;
for i = 1:1:length(C)
    if C(1,i) == 0
C(:,i) = NaN ;
    end
end
indices = find(C > 0 & C < 10,1);
[row, col] = ind2sub(size(C), indices);
C = C(:,1:col-1) ;
cx = [cx;C(1,:)'] ;
cy = [cy;C(2,:)'] ;
save 'cx-cy.mat' cx cy
    end
end
load 'cx-cy.mat' cx cy
clear rowGrid colGrid cols R te rows depth_min depth_max C Z_masked coast_range h indices zeroidx YY_eto_NN XX_eto_NN Z_eto_NN c_lon c_lat
%%
% Defining Coastline
run = 2 ;
for i = 1:1:1
    if run == 1
clf
hold on
plot(cx,cy,'k')
xlim([-75,-37])
ylim([59,77])
daspect([1 aspect_ratio 1])
x_coast = [];
y_coast = [];
run = true;
while run
    [x, y, button] = ginput(1);
    if isempty(button) || button == 13
        run = false;
    else
        plot(x, y, 'ro'); % Plot the point
        xlim([x-3, x+3])
        ylim([y-3, y+3])
        x_coast = [x_coast; x];
        y_coast = [y_coast; y];
    end
end
save("x_coast.mat","x_coast")
save("y_coast.mat","y_coast")
load("x_coast.mat")
load("y_coast.mat")
    end
end
% refference coastline (for angles rectangles)
run = 2 ;
for i = 1:1:1
    if run == 1
clf
hold on
plot(cx,cy,'k')
xlim([-75,-30])
ylim([59,77])
daspect([1 aspect_ratio 1])
x_reff = [];
y_reff = [];
run = true;
while run
    [x, y, button] = ginput(1);
    if isempty(button) || button == 13
        run = false;
    else
        plot(x, y, 'ro'); % Plot the point
        xlim([x-10, x+10])
        ylim([y-10, y+10])
        x_reff = [x_reff; x];
        y_reff = [y_reff; y];
    end
end
save("x_reff.mat","x_reff")
save("y_reff.mat","y_reff")
    end
load("x_reff.mat")
load("y_reff.mat")    
end
% added coastline (DON'T THINK THIS IS NEEDED?)
%run = 2 ;
%for i = 1:1:1
%    if run == 1
%clf
%hold on
%plot(cx,cy,'k')
%plot(x_coast,y_coast,'r')
%xlim([-45,-30])
%ylim([59,77])
%daspect([1 aspect_ratio 1])
%extra_x = [] ;
%extra_y = [] ;
%run = true;
%while run
%    [x, y, button] = ginput(1);
%    if isempty(button) || button == 13
%        run = false;
%    else
%        plot(x, y, 'ro'); % Plot the point
%        xlim([x-3, x+3])
%        ylim([y-3, y+3])
%        extra_x = [extra_x; x];
%        extra_y = [extra_y; y];
%    end
%end
%save 'extra_coast.mat' extra_x extra_y
%    end
%end
%load 'extra_coast.mat' extra_x extra_y
%reff_cord = [] ;
%combine coasts
%x_coast = [x_coast;extra_x] ;
%y_coast = [y_coast;extra_y] ;
for i = 1:1:length(x_reff)-1
    x = [(x_reff(i)+x_reff(i+1))/2,(y_reff(i)+y_reff(i+1))/2] ;   %middle point of refference lines
    reff_cord(i,:) = x ;
end
for i = 1:1:length(x_reff)-1
    x = (y_reff(i+1,1)-y_reff(i,1))/(x_reff(i+1,1)-x_reff(i,1)) ;
    reff_slope(i,:) = x ;
end
num_points = 10 ;
%%
for i=1:1:length(x_reff)-1
x_reff_new(i,:) = linspace(x_reff(i),x_reff(i+1),num_points) ; % each row is a different line segment
slope_new(i,:) = repmat(reff_slope(i,:),1,num_points) ;
y_reff_new(i,:) = reff_slope(i,:)*(x_reff_new(i,:)-x_reff(i))+y_reff(i) ;
end
x_reff_new = reshape(x_reff_new,1,[]) ;
y_reff_new = reshape(y_reff_new,1,[]) ;
slope_new = reshape(slope_new,1,[]) ;
clear('x','y','button')
%NOTHING CHANGED PAST HERE
run = 2 ;
if run == 1 
% extended coastline code 
combined_y = zeros(size(y_coast));
combined_x = zeros(size(x_coast));
% Initial plot to allow zooming and panning
figure(1);
daspect([1 aspect_ratio 1]);
plot(x_coast, y_coast, '-b', 'LineWidth', 2); % Original coastline
xlim([-80, -30]);
ylim([55, 80]);
title('Zoom and Pan to Your Desired View, Then Press Any Key to Start');
xlabel('Longitude');
ylabel('Latitude');
grid on;
hold on;
% Pause to allow zoom and pan
disp('Zoom and pan as needed, then press any key to start...');
pause; % Wait for key press

% Loop through each point on the coastline
for i = 1:(length(y_coast)-1)
    % Calculate the bearing (azimuth) from the current point to the next
    [azm] = azimuth(y_coast(i), x_coast(i), y_coast(i+1), x_coast(i+1), 'degrees');
    
    % Calculate both +90 and -90 perpendicular bearings
    bearing_perpendicular1 = mod(azm + 90, 360); % Clockwise perpendicular
    bearing_perpendicular2 = mod(azm - 90, 360); % Counterclockwise perpendicular
    
    % Use the reckon function to compute the new points 100 km away
    [lat_new1, lon_new1] = reckon(y_coast(i), x_coast(i), km2deg(100), bearing_perpendicular1, 'degrees');
    [lat_new2, lon_new2] = reckon(y_coast(i), x_coast(i), km2deg(100), bearing_perpendicular2, 'degrees');
    
    % Display options on a plot for visual inspection
    figure(1);
    daspect([1 aspect_ratio 1])
    plot(cx,cy,'k')
    xlim([-80,-30])
    ylim([55,80])
    plot(x_coast, y_coast, '-b', 'LineWidth', 2); % Original coastline
    hold on;
    p1 = plot(lon_new1, lat_new1, 'or', 'MarkerSize', 8); % +90 option (Red)
    p2 = plot(lon_new2, lat_new2, 'og', 'MarkerSize', 8); % -90 option (Green)
    current_point = plot(x_coast(i), y_coast(i), 'ok', 'MarkerSize', 10, 'MarkerFaceColor', 'k'); % Current point (Black)
    legend([p1, p2, current_point], {'+90 Option (Red)', '-90 Option (Green)', 'Current Point (Black)'}, 'Location', 'best');
    title(sprintf('Point %d of %d: Select Option', i, length(y_coast)-1));
    hold off;
      choice = questdlg(sprintf('Which option is correct for point %d?', i), ...
                      'Select Direction', ...
                      '+90', '-90', 'Skip', '+90');
      switch choice
        case '+90'
            combined_y(i) = lat_new1;
            combined_x(i) = lon_new1;
        case '-90'
            combined_y(i) = lat_new2;
            combined_y(i) = lon_new2;
        case 'Skip'
            % Do nothing (keep NaN for this point)
            continue;
       end
end
% Final plot of the selected coastal zone
figure(2);
geoplot(y_coast,x_coast, '-b', 'LineWidth', 2); % Original coastline
hold on;
geoplot(combined_y,combined_x, '--r', 'LineWidth', 2); % Selected coastal zone
legend('Original Coastline', '100 km Coastal Line');
title('Final Coastal Zone');
% remove intersecting points by hand (needs to be run after you've run the previous portion
combined_x(end) = [] ;
combined_y(end) = [] ;
%hand select points to delete (doesn't need to be run more than once)
%cursorMode = datacursormode(gcf); % run first
%dataTips = cursorMode.getCursorInfo;
xValues = zeros(length(dataTips), 1);
yValues = zeros(length(dataTips), 1);
for i = 1:length(dataTips)
    xValues(i) = dataTips(i).Position(1);
    yValues(i) = dataTips(i).Position(2);
end
for i = 1:length(xValues)
    delete_idx(i,1) = dataTips(i).DataIndex ; 
end
combined_x(delete_idx) = [] ;
combined_y(delete_idx) = [] ;
save combined_x.mat combined_x
save combined_y.mat combined_y
end
load combined_x.mat
load combined_y.mat

% NOTHING CHANGED PAST HERE FOR REWORKING COAST
polygon_x = [x_coast] ; % may need to add more points but I think this is fine?
polygon_y = [y_coast] ;
combined_x = [x_coast',flip(combined_x')] ; % may need to add more points but I think this is fine?
combined_y = [y_coast',flip(combined_y')] ;
in = inpolygon(lon,lat,combined_x,combined_y) ; % All coasts within coastal section 
% Determine which refference point to base slope of rectangle boxes
reff_x = x_reff_new' ;
reff_y = y_reff_new' ;
coastal_lon = lon(in) ;
coastal_lat = lat(in) ;
for i= 1:length(lon(in))
    temp_x = [repmat(coastal_lon(1,i),length(reff_x),1),reff_x] ;
    temp_y = [repmat(coastal_lat(1,i),length(reff_y),1),reff_y] ;
for j = 1:length(temp_y)
    temp_distance(j,:) = sw_dist(temp_y(j,:),temp_x(j,:),'km') ;
    [min_distance(:,i),idx(:,i)] = min(temp_distance) ; %idx tells which refference point each coastal cast is closest too and should have their slope based upon 
end
end
idx = idx(1,:) ;
clear temp_x temp_y temp_distance min_distance y_reff_new x_reff_new
%% Defining parallelogram verticies for casts (length and width need work)
% 40x20 km boxes for coastal profiles
target_width = recwidth/2 ;
target_length = reclength/2 ;
run = 2 ;
if run == 1 % 10 km towards and away from refference line
width_nm = km2nm(target_width*2);
width_deg = nm2deg(width_nm);
length_nm = km2nm(target_length*2);
length_deg = nm2deg(length_nm);
offset_dist = km2deg(10) ;
    for i = 1:length(x_reff)-1
    azm(i) = azimuth(y_reff(i), x_reff(i), y_reff(i+1), x_reff(i+1), 'degrees');
    % Calculate both +90 and -90 perpendicular bearings
    bearing_perpendicular1 = mod(azm + 90, 360); % Clockwise perpendicular
    bearing_perpendicular2 = mod(azm - 90, 360); % Counterclockwise perpendicular
    [lat1,lon1] = reckon(y_reff(i),x_reff(i), offset_dist, azm(i),'degrees') ;
    if i > 1
        [lat2, lon2] = reckon(y_reff(i), x_reff(i), -offset_dist, azm(i-1), 'degrees'); % Backward point
    else
        lat2 = NaN; % Placeholder for the first iteration
        lon2 = NaN;
    end
    x_reff_mid1(:,i) = lon1; % Clockwise perpendicular
    y_reff_mid1(:,i) = lat1;
    x_reff_mid2(:,i) = lon2; % Counterclockwise perpendicular
    y_reff_mid2(:,i) = lat2;
    end
    x_reff_mid2(1) = [] ;
    x_reff_mid2(end+1) = NaN ;
    y_reff_mid2(1) = [] ;
    y_reff_mid2(end+1) = NaN ;
    x_reff_mid = [x_reff_mid1;x_reff_mid2] ;
    y_reff_mid = [y_reff_mid1;y_reff_mid2] ;
    %find index of nearest refference line
    for i = 1:length(coastal_lon)
         temp_dist = zeros(2, length(x_reff) - 1); % Two rows for distances (forward and backward points)
        for j = 1:length(x_reff)-1
        % Calculate distances for both rows (forward and backward points)
        for k = 1:2 % Loop through each row of combined midpoints
            temp_coords = [coastal_lon(i), x_reff_mid(k, j); coastal_lat(i), y_reff_mid(k, j)];
            temp_dist(k, j) = sw_dist(temp_coords(1, :), temp_coords(2, :), 'km');
        end
        end
    % Find the minimum distance and its corresponding column
    [~, min_idx] = min(temp_dist(:)); 
    [row_idx, col_idx] = ind2sub(size(temp_dist), min_idx); % Convert linear index back to row/column
    min_indices(i) = col_idx;
    end
    for i = 1:length(coastal_lon)
    [lat_side,lon_side] = reckon(coastal_lat(i),coastal_lon(i),width_deg/2,bearing_perpendicular1(min_indices(i))) ; % gets coordinate a point on side of the rectangle
    side_coords = [lat_side;lon_side] ; 
    vert1 = reckon(side_coords(1,1),side_coords(2,1),length_deg/2,bearing_perpendicular1(min_indices(i))+90) ; % 1/2 on same side, use these to find 3/4 on opposite side
    vert2 = reckon(side_coords(1,1),side_coords(2,1),length_deg/2,bearing_perpendicular1(min_indices(i))-90) ;
    vert3 = reckon(vert1(1),vert1(2),width_deg,mod(bearing_perpendicular1(min_indices(i))+180,360));
    vert4 = reckon(vert2(1),vert2(2),width_deg,mod(bearing_perpendicular1(min_indices(i))+180,360));
    vertices{i} = [vert1;vert2;vert3;vert4] ;
        for j = 1:4 % order vertices for inpolygon use
        angles{i}(j) = atan2(vertices{i}(j,2) - coastal_lon(i), vertices{i}(j,1) - coastal_lat(i));    
        end
    [~, sorted_indices] = sort(angles{i});
    vert{i} = vertices{i}(sorted_indices, :);
    end    
save("vert.mat","vert")
end
load("vert.mat","vert")
clear width_deg width_nm length_nm length_deg  bearing_perpendicular1 bearing_perpendicular2 azm temp_coords temp_coords temp_dist min_idx angles vertices offset_dist
clear sorted indices vert1 vert2 vert3 vert4 side_coords lat_side lon_side min_indices x_reff_mid1 x_reff_mid2 x_reff_mid y_reff_mid1 y_reff_mid2 y_reff_mid lon1 lon2 lat1 lat2
%% get rid of compartment
% Interpolate every 1 m for the region (will eventually need to exclude
% fjord stuff via inpolygon indexing)
dep = cellfun(@double, dep, 'UniformOutput', false);
max_length = max(cellfun(@max,dep)) ; %might be source of problem as its based off length not max value
DepInterval = (0:1:max_length); % 1 m intervals
run = 2 ; % will rerun interpolated, temp and salinity as well as cast box indicies, when changing box size you need to set to 1
for i = 1:1:length(dep)
    if run == 1 
        interp_temp{1,i} = interp1(dep{1,i},temp{1,i}, DepInterval, 'linear') ;
    end
end
for i = 1:1
    if run == 1
save ("interp_temp.mat", "interp_temp", '-v7.3');
    end
end
    load('interp_temp.mat')
%Salinity 
for i = 1:1:length(dep)
   if run == 1
    interp_sal{1,i} = interp1(dep{1,i},sal{1,i}, DepInterval, 'linear') ;
   end
end
for i = 1:1
    if run == 1
 save ("interp_sal.mat", "interp_sal", '-v7.3');
    end
end
load('interp_sal.mat')
%% Test Section
% remove fjord casts and create variables to use 
in_idx = inpolygon(lon,lat,polygon_x,polygon_y) ; % need this for the interp section
lon_fj = lon(in_idx) ;
lat_fj = lat(in_idx) ;
lon_a = lon(~in_idx) ;
lat_a = lat(~in_idx) ;
temp_a = interp_temp(~in_idx) ;
sal_a = interp_sal(~in_idx) ;
%ranges_a = ranges(:,~in_idx) ;
mon_a = mon(~in_idx) ;
day_a = day(~in_idx) ;
yea_a = yea(~in_idx) ;
mon_fj = mon(in_idx) ;
day_fj = day(in_idx) ;
yea_fj = yea(in_idx) ;
% format into matrices (hopefully can get rid of this subsection)
interp_temp_a = interp_temp(~in_idx) ; % gets rid of fjord casts
interp_sal_a = interp_sal(~in_idx) ;
numCells = length(interp_temp_a) ;
rows = length(DepInterval) ;
interp_sal_mat =  zeros(rows, numCells) ;
interp_temp_mat =  zeros(rows, numCells) ;
run = 2 ;
    if run == 1 
for i = 1:length(interp_temp_a)
insidearray1 = interp_temp_a{i}' ;
insidearray2 = interp_sal_a{i}' ;
interp_temp_mat(:,i) = insidearray1 ;
interp_sal_mat(:,i) = insidearray2 ;
end
save("interp_sal_mat.mat",'interp_sal_mat','-v7.3')
save("interp_temp_mat.mat",'interp_temp_mat','-v7.3')
    end
load("interp_sal_mat.mat")
load("interp_temp_mat.mat")
% Clean Sal_mat data
%create coastal variables to clean vert along with everything else
in_temp = inpolygon(lon_a,lat_a,combined_x,combined_y) ;
coast_sal_temporary = interp_sal_mat(:,in_temp) ;
%% Selecting Fjord coordinates/area % defining areas of individual fjords
run = 2 ;
if run == 1
    hold on
    plot(cx,cy,'k')
    scatter(lon_fj,lat_fj,0.7,'r')
    daspect([1 aspect_ratio 1])
        xlim([-80,-35])
        ylim([55,80])
        fjord_vert = {} ;
        while true
            zoom on;
            pan on;
            disp('Adjust view as needed, then press Enter to start drawing a polygon.');
            pause;
            zoom off;
            pan off;
            h = drawpolygon('Color', 'yellow'); % Create an interactive polygon
            wait(h); % Wait until the polygon is finished
            vertices = h.Position;
            fjord_vert{end + 1} = vertices;
            delete(h); % Remove the polygon object from the figure
            plot(vertices(:,1), vertices(:,2), '-o', 'Color', 'blue');
                    choice = questdlg('Do you want to draw another polygon?', ...
                          'Continue Drawing', 'Y', 'N','Y');
        % Check for exit condition
        if strcmp(choice, 'No') % Exit if 'No' is selected
            break; % Exit the loop
        end
        end
        % Replot all previously stored vertices
        for i = 1:length(fjord_vert)
            plot(fjord_vert{i}(:, 1), fjord_vert{i}(:, 2), 'o', 'Color', 'yellow', 'MarkerFaceColor', 'yellow');
        end
    hold off
    save('fjord_vert.mat', 'fjord_vert');
end
load fjord_vert.mat fjord_vert
%fjord_vert = fjord_vert(1:29) ; % cutoff empty cell
%% Selecting center points for fjord boxes for anomaly calculations
%select center points
run = 2 ;
if run == 1
    hold on
    scatter(lon_fj,lat_fj,0.7,'r')
    scatter(lon_a,lat_a,0.7,'r')
    isobath_interval = -200:-100:-3000;  % Define isobath intervals down to the minimum depth
    iso_interval = 200:100:3000 ;
    [~, ContourMatrix] = contour(XX_eto, YY_eto, Z_eto, isobath_interval);
    [~,CountourM] = contour(XX,YY,ZZ, iso_interval) ;
    [~,ContourE] = contour(XX_eto_E,YY_eto_E,Z_eto_E,isobath_interval) ;
    daspect([1 aspect_ratio 1])
    xlim([-80,-30])
    ylim([55,80])
    center_points = [];
        for i = 1:length(fjord_vert)
            fill(fjord_vert{i}(:,1),fjord_vert{i}(:,2),'b') ;
        end
    plot(cx,cy,'k') % plotting coast on top for better visibility
    % zoom and pan until ready, then hit enter to plot point
    disp('Click on the map to select points. Press Enter when done.');
    while true
            zoom on;
            pan on;
            disp('Adjust view as needed, then press Enter to select a point.');
            pause;
            zoom off;
            pan off;
        [lon_cen, lat_cen, button] = ginput(1); % Select a point
        scatter(lon_cen,lat_cen,2,'g')
        if isempty(button) % Stop if Enter is pressed
        break;
        end
    center_points = [center_points; lon_cen, lat_cen]; % Store selected point
    end
    hold off
save('center_points.mat', "center_points")
end
load center_points.mat 
clear isobath_interval iso_interval
%% Make rectangles 20x10 km using reckon around a center point, can adjust length in special cases
run = 2;
if run == 1
    % Convert width and length from km to degrees
    width_nm = km2nm(target_width*2);
    width_deg = nm2deg(width_nm);
    length_nm = km2nm(target_length*2);
    length_deg = nm2deg(length_nm);
    hold on
    scatter(lon_fj, lat_fj, 0.7, 'r')
    scatter(lon_a, lat_a, 0.7, 'r')
    isobath_interval = -200:-100:-3000;  % Define isobath intervals down to the minimum depth
    iso_interval = 200:100:3000 ;
    [~, ContourMatrix] = contour(XX_eto, YY_eto, Z_eto, isobath_interval);
    [~, ContourMatrix_N] = contour(XX_eto_N, YY_eto_N, Z_eto_N, isobath_interval);
    [~,CountourM] = contour(XX,YY,ZZ, iso_interval) ;
    [~,ContourE] = countour(XX_eto_E,YY_eto_E,Z_eto_E,isobath_interval) ;
    for i = 1:length(fjord_vert)
        fill(fjord_vert{i}(:, 1), fjord_vert{i}(:, 2), 'b');
    end
    plot(cx,cy,'k') % plotting coast on top for better visibility
    daspect([1, aspect_ratio, 1])
    xlim([-80, -30])
    ylim([55, 80])
    for i = 1:size(center_points, 1)
        lon_cen = center_points(i, 1);
        lat_cen = center_points(i, 2);
        angle_deg = 0;  % Start angle for each rectangle
        scatter(center_points(:, 1), center_points(:, 2), 15, 'g', 'filled')
        while true
            h = findobj(gca, 'Tag', sprintf('rectangle%d', i));
            delete(h);
            [lat1, lon1] = reckon(lat_cen, lon_cen, width_deg / 2, angle_deg);
            [lat2, lon2] = reckon(lat_cen, lon_cen, width_deg / 2, angle_deg + 180);
            [lat3, lon3] = reckon(lat1, lon1, length_deg / 2, angle_deg + 90);
            [lat4, lon4] = reckon(lat1, lon1, length_deg / 2, angle_deg - 90);
            [lat5, lon5] = reckon(lat2, lon2, length_deg / 2, angle_deg + 90);
            [lat6, lon6] = reckon(lat2, lon2, length_deg / 2, angle_deg - 90);
            rectangle_lat{i} = [lat3, lat4, lat6, lat5, lat3];
            rectangle_lon{i} = [lon3, lon4, lon6, lon5, lon3];
            plot(rectangle_lon{i}, rectangle_lat{i}, 'g-', 'LineWidth', 2, 'Tag', sprintf('rectangle%d', i));
            drawnow;
            angle_input = input('Enter new angle (or press Enter to confirm): ', 's');
            if isempty(angle_input)
                disp(['Angle confirmed for rectangle ', num2str(i)]);
                break;  % Exit the loop for this rectangle
            else
                % Update the angle for next iteration
                angle_deg = str2double(angle_input);
                if isnan(angle_deg)
                    warning('Invalid angle. Try again.');
                    angle_deg = 0;
                end
            end
        end
        % Ask if it's a special case
        special_case = questdlg('Is this a special box case?', 'Special Case', 'Yes', 'No', 'No');

        if strcmp(special_case, 'Yes')
            while true
                % Prompt to enter new length for the special case
                new_length_input = input('Enter new length (or press Enter to confirm): ', 's');
                if isempty(new_length_input)
                    disp(['Length confirmed for special box ', num2str(i)]);
                    break;  % Exit the loop for this special case
                else
                    % Update length if valid
                    new_length_km = str2double(new_length_input);
                    if isnan(new_length_km)
                        warning('Invalid length. Try again.');
                    else
                        % Convert new length from km to degrees
                        length_deg = nm2deg(km2nm(new_length_km));
                        
                        % rectangle with the new length
                        [lat3, lon3] = reckon(lat1, lon1, length_deg / 2, angle_deg + 90);
                        [lat4, lon4] = reckon(lat1, lon1, length_deg / 2, angle_deg - 90);
                        [lat5, lon5] = reckon(lat2, lon2, length_deg / 2, angle_deg + 90);
                        [lat6, lon6] = reckon(lat2, lon2, length_deg / 2, angle_deg - 90);
                        rectangle_lat{i} = [lat3, lat4, lat6, lat5, lat3];
                        rectangle_lon{i} = [lon3, lon4, lon6, lon5, lon3];
                        
                        % Update the plot
                        h = findobj(gca, 'Tag', sprintf('rectangle%d', i));
                        delete(h);
                        plot(rectangle_lon{i}, rectangle_lat{i}, 'g-', 'LineWidth', 2, 'Tag', sprintf('rectangle%d', i));
                        drawnow;
                    end
                end
            end
        end
    end
    save rectangle_lat.mat rectangle_lat
    save rectangle_lon.mat rectangle_lon
end
load rectangle_lon.mat
load rectangle_lat.mat
for i = 1:length(rectangle_lon)
fjord_box_cords{i} = [rectangle_lon{i}',rectangle_lat{i}'] ;
end
clear isobath_interval iso_interval rectangle_lon rectangle_lat
%% Extend fjord area to meet anomaly box (if special case i.e. close enough to justify it) manually done since small #
run = 2;
if run == 1 
    hold on
    scatter(lon_fj, lat_fj, 0.7, 'r')
    scatter(lon_a, lat_a, 0.7, 'r')
    for i = 1:length(fjord_vert)
        fill(fjord_vert{i}(:, 1), fjord_vert{i}(:, 2), 'b');
        plot(fjord_box_cords{i}(:, 1), fjord_box_cords{i}(:, 2), 'g')
        scatter(fjord_vert{i}(1,1),fjord_vert{i}(1,2),'k','filled') % lets me know which point i need to add first
    end
    plot(cx, cy, 'k')
    daspect([1, aspect_ratio, 1])
    xlim([-80, -30])
    ylim([55, 80])
    selected_points = {};
    continue_selecting = true;
    while continue_selecting
        zoom on
        pan on
        %msgbox('Use zoom and pan to adjust the view, then close this box to select points.');
        pause
        waitfor(gcf, 'CurrentCharacter', char(13)); % Wait for Enter key
        zoom off
        pan off
        disp('Click on the points you want to add. Press Enter when done with each set.');
        [lon_new, lat_new] = ginput;
        selected_points{end+1} = [lon_new, lat_new];
        plot(lon_new, lat_new, 'ro-', 'LineWidth', 1.5, 'MarkerSize', 6);
        more_points = questdlg('Would you like to adjust the view and add another set of points?', 'Continue', 'Yes', 'No', 'No');
        if strcmp(more_points, 'No')
            continue_selecting = false;
        end
    end
   hold off
    save selected_points.mat selected_points
end
clear lat_fj lon_fj
load selected_points.mat
%% edit fjord_vert points based on the points selected (don't really need this, faster to edit by hand)
run = 2; % only run if you've reselected boxes or # of fjords
if run == 1
fjord_affected = [1,2,3,4,5,6,7,9,10,11,12,13,15,17,18,19,20,21,22,23,24,26,27,28,30] ; % idx of fjords that will have points altered.(manual)
for i = 1:length(fjord_affected) 
fjord_vert{fjord_affected(i)} = [selected_points{i}(1,:);fjord_vert{fjord_affected(i)}] ; % adds first selected point in front of fjord vert
fjord_vert{fjord_affected(i)} = [fjord_vert{fjord_affected(i)};selected_points{i}(2,:)] ; % adds second selected point to end of fjord vert
end
save fjord_vert.mat fjord_vert
end
load fjord_vert.mat fjord_vert
clear fjord_affected selected_points
% group fjord data by individual fjords
fj_find_combined = [];
for i = 1:length(fjord_vert)
    in_fj{i} = inpolygon(lon, lat, fjord_vert{i}(:,1), fjord_vert{i}(:,2));
    fjord_find = find(in_fj{i});
    fj_find_combined = [fj_find_combined, fjord_find]; % Fj profile indices
end
OMG_fj_profiles = fj_find_combined >= length_NODC ; % OMG profiles, (should be any profile after NODC length)
clear center_points lat1 lon1 lat2 lon2 lat3 lon3 lat4 lon4 width_nm width_deg length_nm length_deg clear center_points special_case iso_interval isobath_interval fjord_find in_fj selected_points
%% clean fjord data loop through collumns and calculate derivatives (delta y /delta x)
% fjord sal
    sal_mat_fj = cellfun(@(c) c(:), interp_sal, 'UniformOutput', false) ; % fjord data 
    sal_mat_fj = [sal_mat_fj{:}] ; % all data, idx fjord profiles, clean and replace with nan for removed profiles
    fjord_sal_mat_fj = sal_mat_fj(:,fj_find_combined) ; % fjord data
%fjord temp
    temp_mat_fj = cellfun(@(c) c(:), interp_temp, 'UniformOutput', false) ; % fjord temp data 
    temp_mat_fj = [temp_mat_fj{:}] ; % all data, idx fjord profiles, clean and replace with nan for removed profiles
    fjord_temp_mat_fj = temp_mat_fj(:,fj_find_combined) ; % fjord data
    lon_fj =  lon(fj_find_combined) ;
    lat_fj = lat(fj_find_combined) ;
    mon_fj = mon(fj_find_combined) ;
    day_fj = day(fj_find_combined) ;
    yea_fj = yea(fj_find_combined) ;
    diff_result = NaN(size(fjord_sal_mat_fj)) ;
    % remove profiles with certain % of NaN in top 100 meters
    top_50 = fjord_sal_mat_fj(1:50,:) ;
    nan_count = sum(isnan(top_50), 1);  % Count NaNs along rows for each column
    threshold = 25 ; % max number of NaN values permitted in top 50 m
    too_many_nans_fj = nan_count >= threshold ;
    fjord_sal_mat_fj = fjord_sal_mat_fj(:,~too_many_nans_fj) ;
    fjord_temp_mat_fj = fjord_temp_mat_fj(:,~too_many_nans_fj) ;
    lon_fj = lon_fj(:,~too_many_nans_fj) ; % for tracking
    lat_fj = lat_fj(:,~too_many_nans_fj) ;
    mon_fj = mon_fj(:,~too_many_nans_fj) ;
    day_fj = day_fj(:,~too_many_nans_fj) ;
    yea_fj = yea_fj(:,~too_many_nans_fj) ;
    run = 2 ;
if run == 1
    for i  = 1:length(fjord_sal_mat_fj)
    non_nan = find(~isnan(fjord_sal_mat_fj(:,i))) ;
    non_nan_store{i} = non_nan ; 
    temporary = fjord_sal_mat_fj(:,i) ;
    diff_sal = diff(temporary(non_nan),1,1) ;
    diff_result(non_nan(1:end-1),i) = abs(diff_sal ./ diff(non_nan,1,1)) ;
    end
   threshold_25 = 12 ; % 0-25 m
   threshold_50 = 2 ; % 26-50 m
   threshold = 0.3 ; % 51-end m
     top_25 = diff_result(1:25,:) ;
     top_50 = diff_result(26:50,:) ;
     bottom = diff_result(51:end,:) ;
     top_25_remove = top_25 >= threshold_25 ; % index of top 25 m that has derivatives greater than threshold
     top_50_remove = top_50 >= threshold_50 ; % idx of 26-50m that has derivatives greater than threshold
     bottom_remove = bottom >= threshold ; % idx of bottom that has derivatives greater than threshold
     remove_fj = [top_25_remove;top_50_remove;bottom_remove] ; % need to account for both points that create a true remove value
     for i = 1:length(remove_fj) 
     idx = find(remove_fj(:,i) == 1) ;
        for j = 1:length(idx)
            if size(idx) >= 1
        value_below_idx = find(non_nan_store{i} == idx(j,:)) + 1 ;
        value_below = non_nan_store{i}(value_below_idx) ;
        remove_fj(value_below, i) = 1; % Mark the value from non_nan_store for removal
            end
        end
     end
save 'remove_fj.mat' 'remove_fj'
clear top_50 bottom bottom_remove diff_sal temporary non_nan value_below value_below_idx non_nan_store top_50_remove top_25_remove top_25 
% compare removed points to OMG list, ideally you are removing little to no OMG data since they've been cleaned already 
remove_fj_any = any(remove_fj, 1);
num_total = sum(remove_fj_any) ;
compare_OMG = remove_fj_any & OMG_fj_profiles ;
compare_OMG_nan = OMG_fj_profiles & too_many_nans_fj ;
num_OMG_nan = sum(compare_OMG_nan) ;
num_OMG = sum(compare_OMG);
disp(['Number of total Fjord profiles impacted: ', num2str(num_total)]) ;
disp(['Number of OMG profiles impacted by derivative cleaning: ', num2str(num_OMG)]);
disp(['Number of OMG profiles impacted by NaN cleaning: ', num2str(num_OMG_nan)]) ; % if you change derivative cleaning the 66 # needs to change too, run without removing any nans and set that value
end
load remove_fj.mat remove_fj
clear compare_OMG remove_fj_any num_total num_OMG too_many_nans num_OMG_nan compare_OMG_nan diff_result threshold_50 threshold threshold_25
% remove these points
fjord_sal_mat_fj(remove_fj) = NaN ;
fjord_temp_mat_fj(remove_fj) = NaN ;
%interp_temp_mat_fj(remove_fj) = NaN ; (not a thing yet)
%% Clean fjord box anomaly profiles (by same standards as open/coastal) (needs to be light/conservative)
 box_find_combined = [] ;
 for i = 1:length(fjord_box_cords)
 temp_idx{i} = inpolygon(lon,lat,fjord_box_cords{i}(:,1),fjord_box_cords{i}(:,2)) ;
 box_find = find(temp_idx{i});
 box_find_combined = [box_find_combined, box_find]; % Fj profile indices (contains duplicates)
 end
   box_find_combined = unique(box_find_combined) ; 
   OMG_box_profiles = box_find_combined >= length_NODC ; %idx of all box profiles that are also OMG
    box_idx = any(cell2mat(temp_idx'), 1);
    box_sal_mat = sal_mat_fj(:,box_idx) ; % only box salinities 
    box_temp_mat = temp_mat_fj(:,box_idx) ; % only box temps
top_100 = box_sal_mat(1:100,:) ;
nan_count = sum(isnan(top_100), 1);  % Count NaNs along rows for each column
    lon_box =  lon(box_find_combined) ;
    lat_box = lat(box_find_combined) ;
    mon_box = mon(box_find_combined) ;
    day_box = day(box_find_combined) ;
    yea_box = yea(box_find_combined) ;
run = 2 ;
if run == 1
% derivatives 
diff_result = NaN(size(box_sal_mat)) ;
for i  = 1:size(box_sal_mat,2)
    non_nan = find(~isnan(box_sal_mat(:,i))) ;
    non_nan_store{i} = non_nan ;
    temporary = box_sal_mat(:,i) ;
    diff_sal = diff(temporary(non_nan),1,1) ;
    diff_result(non_nan(1:end-1),i) = abs(diff_sal ./ diff(non_nan,1,1)) ;
    end
   threshold_50 = 2 ; % 0-50 m
   threshold = 0.3 ; % 50-end m
     top_50 = diff_result(1:50,:) ;
     bottom = diff_result(51:end,:) ;
     top_50_remove = top_50 >= threshold_50 ; % idx of bottom that has derivatives greater than threshold
     bottom_remove = bottom >= threshold ; % idx of bottom that has derivatives greater than threshold
     remove_box = [top_50_remove;bottom_remove] ; % need to account for both points that create a true remove value
     for i = 1:size(remove_box,2) 
     idx = find(remove_box(:,i) == 1) ;
        for j = 1:length(idx)
            if size(idx) >= 1
        value_below_idx = find(non_nan_store{i} == idx(j,:)) + 1 ;
        value_below = non_nan_store{i}(value_below_idx) ;
        remove_box(value_below, i) = 1; % Mark the value from non_nan_store for removal
            end
        end
     end
% Test
remove_box_any = any(remove_box, 1);
%num_total = sum(remove_box_any) ;
%compare_OMG = remove_box_any & OMG_box_profiles ;
%compare_OMG_nan = OMG_box_profiles & too_many_nans ;
%num_OMG_nan = sum(compare_OMG_nan) ;
%num_OMG = sum(compare_OMG);
%disp(['Number of total box profiles impacted: ', num2str(num_total)]) ;
%disp(['Number of box OMG profiles impacted by derivative cleaning: ', num2str(num_OMG)]);
%disp(['Number of box OMG profiles impacted by NaN cleaning: ', num2str(num_OMG_nan)]) ; % if you change derivative cleaning the 66 # needs to change too, run without removing any nans and set that value
save remove_box.mat remove_box
end
load remove_box.mat remove_box
box_sal_mat(remove_box) = NaN ; % eliminates points 
box_temp_mat(remove_box) = NaN ;
% NaN filter
threshold = 50 ; % max number of NaN values permitted in top 100 m
too_many_nans_box = nan_count > threshold ;
box_sal_mat = box_sal_mat(:,~too_many_nans_box) ;
box_temp_mat = box_temp_mat(:,~too_many_nans_box) ;
    lon_box = lon_box(:,~too_many_nans_box) ; % for tracking
    lat_box = lat_box(:,~too_many_nans_box) ;
    mon_box = mon_box(:,~too_many_nans_box) ;
    day_box = day_box(:,~too_many_nans_box) ;
    yea_box = yea_box(:,~too_many_nans_box) ;
clear temp_idx box_idx threshold_50 threshold top_50 bottom threshold bottom_remove top_50 threshold_50 diff_sal temporary non_nan value_below value_below_idx diff_result non_nan_store remove top_50_remove length_NODC num_OMG
clear num_OMG_nan 
%% loop through collumns and calculate derivatives (delta y /delta x)
run = 2 ;
if run == 1
    diff_result = NaN(size(interp_sal_mat)) ;
    for i  = 1:length(interp_sal_mat)
    non_nan = find(~isnan(interp_sal_mat(:,i))) ;
    non_nan_store{i} = non_nan ;
    temporary = interp_sal_mat(:,i) ;
    diff_sal = diff(temporary(non_nan),1,1) ;
    diff_result(non_nan(1:end-1),i) = abs(diff_sal ./ diff(non_nan,1,1)) ;
    end
   threshold_50 = 2 ; % 0-50 m
   threshold = 0.3 ; % 50-end m
     top_50 = diff_result(1:50,:) ;
     bottom = diff_result(51:end,:) ;
     top_50_remove = top_50 >= threshold_50 ; % idx of bottom that has derivatives greater than threshold
     bottom_remove = bottom >= threshold ; % idx of bottom that has derivatives greater than threshold
     remove = [top_50_remove;bottom_remove] ; % need to account for both points that create a true remove value
     for i = 1:length(remove) 
     idx = find(remove(:,i) == 1) ;
        for j = 1:length(idx)
            if size(idx) >= 1
        value_below_idx = find(non_nan_store{i} == idx(j,:)) + 1 ;
        value_below = non_nan_store{i}(value_below_idx) ;
        remove(value_below, i) = 1; % Mark the value from non_nan_store for removal
            end
        end
     end
save 'remove.mat' 'remove'
clear top_50 bottom threshold bottom_remove top_50 threshold_50 diff_sal temporary non_nan value_below value_below_idx diff_result non_nan_store remove top_50_remove
%same thing but for coast (for vert idx)
    diff_result = NaN(size(coast_sal_temporary)) ;
    for i  = 1:length(coast_sal_temporary)
    non_nan = find(~isnan(coast_sal_temporary(:,i))) ;
    non_nan_store{i} = non_nan ;
    temporary = coast_sal_temporary(:,i) ;
    diff_sal = diff(coast_sal_temporary(non_nan),1,1) ;
    diff_result(non_nan(1:end-1),i) = abs(diff_sal ./ diff(non_nan,1,1)) ;
    end
   threshold_50 = 2 ; % 0-50 m
   threshold = 0.3 ; % 50-end m
     top_50 = diff_result(1:50,:) ;
     bottom = diff_result(51:end,:) ;
     top_50_remove = top_50 >= threshold_50 ; % idx of bottom that has derivatives greater than threshold
     bottom_remove = bottom >= threshold ; % idx of bottom that has derivatives greater than threshold
     remove_coast = [top_50_remove;bottom_remove] ; % need to account for both points that create a true remove value
     for i = 1:length(remove_coast) 
     idx = find(remove_coast(:,i) == 1) ;
        for j = 1:length(idx)
            if size(idx) >= 1
        value_below_idx = find(non_nan_store{i} == idx(j,:)) + 1 ;
        value_below = non_nan_store{i}(value_below_idx) ;
        remove_coast(value_below, i) = 1; % Mark the value from non_nan_store for removal
            end
        end
     end
save 'remove_coast.mat' 'remove_coast'
clear remove_coast top_50 bottom threshold bottom_remove top_50 threshold_50 diff_sal temporary non_nan value_below value_below_idx diff_result non_nan_store top_50_remove
end
load 'remove.mat'
load 'remove_coast.mat'
interp_sal_mat(remove) = NaN ;
interp_temp_mat(remove) = NaN ;
coast_sal_temporary(remove_coast) = NaN ;
% remove anomalous profiles from interp_sal and interp_temp
temp_sal = interp_sal_mat(50:end, :) ;
coast_temp_sal = coast_sal_temporary(50:end,:) ;
bad_profiles = temp_sal <= 30; % finds casts with sals >30 at depths greater than 50 m
bad_coast_profiles = coast_temp_sal <= 30; % finds casts with sals >30 at depths greater than 50 m
remove_profiles = any(bad_profiles,1) ;
remove_coast_profiles = any(bad_coast_profiles) ;
% certain % of non-NaN values in 0-100 m
top_100 = interp_sal_mat(1:100,:) ;
top_100_coast = coast_sal_temporary(1:100,:) ;
nan_count = sum(isnan(top_100), 1);  % Count NaNs along rows for each column
coast_nan_count = sum(isnan(top_100_coast), 1);  % Count NaNs along rows for each column
threshold = 50 ; % max number of NaN values permitted in top 100 m
too_many_nans = nan_count > threshold ;
coast_many_nans = coast_nan_count > threshold ;
% combined idx and remove from dataset
vert_idx = bad_coast_profiles | coast_many_nans ;
vert_idx = any(vert_idx,1) ;
unique_col = remove_profiles | too_many_nans ;
vert(:,vert_idx) = [] ;
coastal_lat(:,vert_idx) = [] ;
coastal_lon(:,vert_idx) = [] ;
interp_sal_mat(:, unique_col) = [];
interp_temp_mat(:,unique_col) = [];
clear salinity_diff threshold abs_diff anomalous_profiles too_many_nans nan_count top_100 bad_profiles remove_profiles abs_diff_expanded in_temp coast_abs_diff coast_sal_temporary
clear insidearray1 insidearray2 rows numCells sal_a temp_a col coast_diff coast_abs_diff_expanded bad_coast_profiles %(need interp_sal_a/temp for open casts)
clear coast_temp_sal remove_coast_profiles coast_many_nans coast_nan_count vert_idx
% edit the variables in preperation for indexing
lon_temp = lon_a ; % for use in coastal indexing section
lat_temp = lat_a ;
lon_a(:,unique_col) = [] ;
lat_a(:,unique_col) = [] ;
mon_a(:,unique_col) = [] ;
yea_a(:,unique_col) = [] ;
day_a(:,unique_col) = [] ;
%%
% Create 15 day indicies regardless of year 
in_a = inpolygon(lon_a,lat_a,combined_x,combined_y) ; % All coasts within coastal section 
date = ([yea; mon; day]) ;
date(4,:) = datenum(0,date(2,:),date(3,:));% everything set in the year 0 to make for easy sorting by day and month
datenum = datenum(0,date(2,:),date(3,:)) ;% everything set in the year 0 to make for easy sorting by day and month
datenum_a = datenum(~in_idx) ;
datenum_a(:,unique_col) = [] ; 
date_a = date(:,~in_idx) ;
date_a(:,unique_col) = [] ; 
datenum_coast = datenum_a(in_a) ;
day_range_2 = day_range*2 ; % For dynamic naming
middle_idx = date_a(4,in_a) >= day_range & date_a(4,in_a) <= 365-day_range; % idx for coastal ranges, will need to go back for open ocean stuff
earlyjan_idx = (date_a(4,in_a) < day_range) ; % only includes earliest part of Jan that overlaps with interval 
latedec_idx = (date_a(4,in_a) > 365-day_range); % only includes latest part of Dec that overlaps with interval 
Jan_range = 365-day_range + (datenum_coast(earlyjan_idx) ); % Earliest Day to include
Jan_range(2,:) = datenum_coast(earlyjan_idx)+day_range ; % Latest Day to include
Middle_range = datenum_coast(middle_idx)- day_range ;
Middle_range(2,:) = datenum_coast(middle_idx)+ day_range ;
Dec_range = datenum_coast(latedec_idx)- day_range ;
Dec_range(2,:) = datenum_coast(latedec_idx)-365+day_range;
range = zeros(2,length(coastal_lat)) ;
range(:,earlyjan_idx) = Jan_range ;
range(:,middle_idx) = Middle_range ;
range(:,latedec_idx) = Dec_range ;
% Open Ocean Casts
open = inpolygon(lon_a,lat_a,combined_x,combined_y) ; % All coasts within coastal section 
datenum_open = datenum_a(~open) ;
lat_open = lat_a(~open) ;
lon_open = lon_a(~open) ;
middle_idx = date_a(4,~open) >= day_range & date_a(4,~open) <= 365-day_range; % idx for coastal ranges, will need to go back for open ocean stuff
earlyjan_idx = (date_a(4,~open) < day_range) ; % only includes earliest part of Jan that overlaps with interval 
latedec_idx = (date_a(4,~open) > 365-day_range); % only includes latest part of Dec that overlaps with interval 
Jan_range_open(1,:) = 365-day_range + (datenum_open(earlyjan_idx) ); % Earliest Day to include
Jan_range_open(2,:) = datenum_open(earlyjan_idx)+day_range ; % Latest Day to include
Middle_range_open = datenum_open(middle_idx)- day_range ;
Middle_range_open(2,:) = datenum_open(middle_idx)+ day_range ;
Dec_range_open = datenum_open(latedec_idx)- day_range ;
Dec_range_open(2,:) = datenum_open(latedec_idx)-365+day_range;
range_open = zeros(2,length(lat_open)) ;
range_open(:,earlyjan_idx) = Jan_range_open ;
range_open(:,middle_idx) = Middle_range_open ;
range_open(:,latedec_idx) = Dec_range_open ;
clear polygon_x polygon_y
clear middle_idx latedec_idx earlyjan_idx day_a
%%
% Get rid of garbage casts
%temp_sal = temp_sal(:,in) ;
%[~, col] = find(temp_sal <= 30); % finds casts with sals >30 at depths greater than 50 m
%unique_col = unique(col, 'stable');
%vert(:,unique_col) = [] ;
%coastal_lon(:,unique_col) = [] ;
%coastal_lat(:,unique_col) = [] ;
%range(:,unique_col) = [] ;
%clear temp_sal
% Create indicies for coastal casts within box and date range 
run = 2 ;
if run == 1
        for i = 1:length(vert)
            poly_idx = inpolygon(lon_a,lat_a,vert{i}(:,2),vert{i}(:,1)) ; % each collumn corresponds to a different cast
                if datenum_coast(i) > 14 && datenum_coast(i) < 350 % Middle ranges
                   coastal_idx = datenum_a >= range(1,i) & datenum_a <= range(2,i) ; 
                elseif datenum_coast(i) <= 14 || datenum_coast(i) >= 351 % Early Jan and Late Dec 
                   coastal_idx = datenum_a >= range(1,i) | datenum_a <= range(2,i) ;
                end
                % combine indices
                combined_idx = poly_idx & coastal_idx ;
                coast_find{i} = find(combined_idx == 1) ;
        end
save("coast_find.mat",'coast_find')
end
load('coast_find.mat') ;
clear combined_idx coastal_idx poly_idx  
% Create find index of values that are found within each box to simplify
% calculations later (I think  I can get rid of this little section?)
% Preallocate cell array to store indices for each column
%indices_cell = cell(1, size(poly_idx, 2));
% Iterate over each column of poly_idx
%for i = 1:size(poly_idx, 2)
    % Find the row indices where the value is 1 in the current column
 %   column_indices = find(poly_idx(:, i) == 1);
    % Store the row indices in the cell array
  %  indices_cell{i} = column_indices;
%end
clear index_minus index_minus_vert index_plus index_plus_vert intersect_x intersect_y idx_remove offset min_minus min_plus num_points reff_cord reff_slope reff_x reff_y 
clear row run slope_coast slope_new target target_length target_width tolerance vert W column_indices exten indices_cell
dummy_lon = lon_a ; % just for purposes of rerunning, delete after
dummy_lat = lat_a ;
%% Use sw_dist to create index for open casts within box and date range 
%refference lat and lon
start_pattern = [true,false] ;
pattern = [true,true,false] ;  %pattern for opendist/sw_dist data sorting
num_repeat = ceil(256629/length(pattern)) ; % # should be equal to length of open_dist -2 have to change if you change length!
pattern = repmat(pattern,1,num_repeat) ;
pattern = [start_pattern,pattern,true] ; % concats start and end values
clear num_repeat start_pattern lon_temp
lon_a = dummy_lon ;
lat_a = dummy_lat ;% for rerunning
if mod(length(lon_a),2) == 1
lon_a = [lon_a,0] ; % add additional garbage point to make even (remove after done)
lat_a = [lat_a,0] ;
end
run = 2 ;
if run == 1
    for i = 1:length(lon_open)
reff_lon = [repmat(lon_open(i),1,length(lon_a)/2)] ; % half the length of lat/lon_a to interleave
reff_lat = [repmat(lat_open(i),1,length(lat_a)/2)] ;
lon_reshaped = reshape(lon_a,2,[]) ;
lat_reshaped = reshape(lat_a,2,[]) ;
lon_inter = reshape([reff_lon;lon_reshaped],1,[]) ;
lon_inter = [lon_inter,lon_open(i)] ;
lat_inter = reshape([reff_lat;lat_reshaped],1,[]) ;
lat_inter = [lat_inter,lat_open(i)] ;
        open_dist = sw_dist(lon_inter,lat_inter,'km') ; % calculates distance for every point to the target cast (i) (needs to be edited, only storing one value)
        open_dist_final = open_dist(pattern) ;%remove trash distances
        circle_idx = open_dist_final <= recwidth ; % all casts less than 20km away
                if datenum_open(i) > 14 && datenum_open(i) < 350 % Middle ranges
                   open_idx = datenum_a >= range_open(1,i) & datenum_a <= range_open(2,i) ; 
                elseif datenum_open(i) <= 14 || datenum_open(i) >= 351 % Early Jan and Late Dec 
                   open_idx = datenum_a >= range_open(1,i) | datenum_a <= range_open(2,i) ;
                end
                        if length(circle_idx) < length(open_idx)
                        circle_idx = [circle_idx,0] ; % to match idx length
                        end
                % combine indices
                combined_idx = circle_idx & open_idx ;
                open_find{i} = find(combined_idx == 1) ;
    end
    save open_find.mat open_find 
end
load open_find.mat
clear lon_inter lat_inter lon_reshaped lat_reshaped reff_lat reff_lon pattern open_idx open_dist_final open_dist num_repeat start_pattern
clear idx temp temp_lon poly_idx date_idx_cell range Dec_range Jan_range circle_idx_open dec_idx date_idx_cell_open Dec_open individual_array jan_idx Jan_open Middle_idx Middle_open non_zero_elements
clear Middle_range Jan_indices Jan_range_open Dec_range_open Middle_range_open ranges_open reff_lon reff_lat off_coast circle_idx_mat day lengthdeg mon numCells sal 
clear i j max_length poly_off radius reclength_coast dist dep
%% remove dummy lon/lat
lon_a = dummy_lon ;
lat_a = dummy_lat ;
clear dummy_lon dummy_lat
% remove from indicies, profiles that contain less than 4 total profiles (3 not including the source profile)
cutoff = 4 ; % number of profiles needed to make a good average
numElements_open = cellfun(@numel, open_find);
numElements_coast = cellfun(@numel, coast_find);
numElements = [numElements_open,numElements_coast] ; % for plotting
open_idx = numElements_open >= cutoff ;
coast_idx = numElements_coast >= cutoff ;
open_find = open_find(open_idx) ;
lat_open = lat_open(open_idx) ; % only will include lat and lon for coast and open that we'll have anomalies for, will leave lat_a lon_a for all together
lon_open = lon_open(open_idx) ;
coast_find = coast_find(coast_idx) ;
coastal_lat = coastal_lat(coast_idx) ;
coastal_lon = coastal_lon(coast_idx) ;
clear numElements_open numElements_coast numElements cutoff
%% coast anomaly setup
in_a = inpolygon(lon_a,lat_a,combined_x,combined_y) ; % All coasts within coastal section 
coastal_sal = interp_sal_mat(:,in_a) ;
coastal_temp = interp_temp_mat(:,in_a) ;
% Coast means/std
run = 2 ;
if run == 1
for i = 1:1
    coastal_temp_mat = [] ;
    coastal_sal_mat = [] ;
    for j = 1:length(coast_find)
    coastal_temp_mat = interp_temp_mat(:,coast_find{j}) ;
    coastal_sal_mat = interp_sal_mat(:,coast_find{j}) ;
    coastal_sal_mean{j} = mean(coastal_sal_mat,2,'omitnan') ; 
    coastal_temp_mean{j} = mean(coastal_temp_mat,2,'omitnan') ;
    coastal_temp_std{j} = std(coastal_temp_mat,0,2,'omitnan') ;
    coastal_sal_std{j} = std(coastal_sal_mat,0,2,'omitnan') ;
    end
end
% Expand
numCells = length(coastal_temp_mean) ;
rows = length(DepInterval) ;
coast_temp_avg =  zeros(rows, numCells) ;
coast_sal_avg =  zeros(rows, numCells) ;
for i = 1:length(coastal_sal_mean) 
insidearray1 = coastal_temp_mean{i}' ;
insidearray2 = coastal_sal_mean{i}' ;
coast_temp_avg(:,i) = insidearray1 ;
coast_sal_avg(:,i) = insidearray2 ;
end
clear coastal_sal_mat coastal_temp_mat clear coastal_sal_mean coastal_temp_mean insidearray2 insidearray1
for i = 1:length(coast_sal_avg(1,:)) 
coast_sal_anom(:,i) = coastal_sal(:,i) - coast_sal_avg(:,i) ;
coast_temp_anom(:,i) = coastal_temp(:,i) - coast_temp_avg(:,i) ;
end
save coast_sal_anom.mat coast_sal_anom
save coast_temp_anom.mat coast_temp_anom
end
load coast_sal_anom.mat % can be used later, will clear now for space
load coast_temp_anom.mat
%combined idx (may or may not need this)
%[~, col] = find(coastal_sal <= 30 | coastal_sal >= 40); % this is in the code twice because it needs to be able to run if run = 2
%unique_col = unique(col, 'stable');
%unique_col_idx_coast = false(1,length(coastal_sal)) ;
%unique_col_idx_coast(unique_col) = true ;
%combined_idx_coast = coast_idx | unique_col_idx_coast ;
%coastal_lat(:, combined_idx_coast) = [];
%coastal_lon(:, combined_idx_coast) = [];
clear coast_sal_avg coast_temp_avg range_open range coastal_sal_std coastal_temp_std date_a watdep twd
%% Open casts
% Coast means/std
in_a = inpolygon(lon_a,lat_a,combined_x,combined_y) ; % All coasts within coastal section 
open_sal = interp_sal_mat(:,~in_a) ;
open_temp = interp_temp_mat(:,~in_a) ;
[~, col] = find(open_sal <= 30| open_sal >= 40);
unique_col = unique(col, 'stable');
unique_col_idx_open = false(1,length(open_sal)) ;
unique_col_idx_open(unique_col) = true ;

run = 2 ;
if run == 1
    % clear for open_find run (everything thats not needed)
    clearvars -except open_find interp_sal_mat interp_temp_mat open_sal open_temp
        for j = 1:length(open_find)
    open_temp_mat = [] ;
    open_sal_mat = [] ;
    open_temp_mat = interp_temp_mat(:,open_find{j}) ;
    open_sal_mat = interp_sal_mat(:,open_find{j}) ;
    open_sal_mean = mean(open_sal_mat,2,'omitnan') ; 
    open_temp_mean = mean(open_temp_mat,2,'omitnan') ;
    open_sal_anom(:,j) = open_sal(:,j) - open_sal_mean ;
    open_temp_anom(:,j) = open_temp(:,j) - open_temp_mean ;
    %open_temp_std{j} = std(open_temp_mat,0,2,'omitnan') ;
    %open_sal_std{j} = std(open_sal_mat,0,2,'omitnan') ; %only run either
    %anomaly or std at a time
        end
   save("open_temp_anom.mat",'open_temp_anom','-v7.3')
   save("open_sal_anom.mat",'open_sal_anom','-v7.3')
  % save("open_temp_std.mat",'open_temp_std')
   %save("open_sal_std",'open_sal_std')
end
load("open_temp_anom.mat") 
load("open_sal_anom.mat")

clear col unique_col
%load("open_temp_std.mat")
%load("open_sal_std.mat")
clear coast_find open_find open_temp_mat open_sal_mat % will need these back eventually
%% Fjord Anomaly Calculations (include profiles + or - 15 days)
run = 2 ;
if run == 1
for i = 1:length(fjord_box_cords)
fj_box_profiles{i} = inpolygon(lon_box,lat_box,fjord_box_cords{i}(:,1),fjord_box_cords{i}(:,2)) ; % gives index of all profiles within the fjord anomaly boxes.
fj_profiles{i} = inpolygon(lon_fj,lat_fj,fjord_vert{i}(:,1),fjord_vert{i}(:,2)) ; % idx of profiles inside defined fjords
end
fj_combined = [] ;
fj_box_combined = [] ;
for i = 1:length(fj_profiles)
fj_box_idx = find(fj_box_profiles{i} == 1) ; % create index instead of logical index
fj_profile_idx = find(fj_profiles{i} == 1) ;
fj_combined = [fj_combined, fj_profile_idx] ;
fj_box_combined = [fj_box_combined,fj_box_idx] ;
fj_box_idx_back{i} = fj_box_idx ;
fj_profile_back{i} = find(fj_profiles{i} == 1) ;
end
clear fj_box_idx 
% get box dates
box_datenum = datenum(fj_box_combined) ; 
box_datenum = box_datenum(:,~too_many_nans_box) ;
datenum_fj_temp = datenum(fj_combined) ;
for i = 1:length(fj_profile_back)
fj_datenum{i} = datenum_fj_temp(fj_profile_back{i}) ;
fj_coords{i} = [lon(fj_profile_back{i});lat(fj_profile_back{i})] ; % matching coordinates for anomalies
end
for i = 1:length(fj_datenum)
    for j= 1:length(fj_datenum{i})
    if fj_datenum{i}(j) > 15 && fj_datenum{i}(j) < 351  % normal cases
        date_idx{i}{j} = abs(fj_datenum{i}(j) - box_datenum) <= 15;  % within ±15 days
    elseif fj_datenum{i}(j) <= 15  % start in January, wrap to December
        date_idx{i}{j} = abs(fj_datenum{i}(j) - box_datenum) <= 15 | abs((fj_datenum{i}(j) + 365) - box_datenum) <= 15;
    elseif fj_datenum{i}(j) >= 351  % start in December, wrap to January
        date_idx{i}{j} = abs(fj_datenum{i}(j) - box_datenum) <= 15 | abs((fj_datenum{i}(j) - 365) - box_datenum) <= 15;
    end
    end
end
% replace logical index in date_idx cells with index from fj_box_combined
for i = 1:length(date_idx)
    for j = 1:length(date_idx{i})
        fj_anom_date_idx = find(date_idx{i}{j} == 1) ;
        fj_anom_idx{i}{j} = intersect(fj_anom_date_idx,fj_box_idx_back{i}) ;
    end
end
% compare date casts and box casts, only keep if a value is present in both (fj_profile_back contains idx for fj profiles, fj_anom_idx for box profiles that fit the correct box and date)_
clear date_idx fj_anom_date_idx box_datenum fj_box_idx_back fj_box_combined fj_combined fj_profiles fj_box_profiles
% calculate means, and subtract fjord values for anomaly
%interp_sal_mat_fj = cellfun(@(c) c(:), interp_sal, 'UniformOutput', false) ;
%interp_sal_mat_fj = [interp_sal_mat_fj{:}] ;
for i = 1:length(fj_anom_idx) 
    for j = 1:length(fj_anom_idx{i})
        % Extract relevant data for the current fj_anom_idx
        relevant_data = box_sal_mat(:, fj_anom_idx{i}{j});
        relevant_temp_data = box_temp_mat(:, fj_anom_idx{i}{j});
        % Initialize the anomaly column with NaN
        fj_anoms{i}(:, j) = NaN(size(fjord_sal_mat_fj, 1), 1);
        fj_temp_anoms{i}(:, j) = NaN(size(fjord_temp_mat_fj, 1), 1);
        % Check validity of rows 1:99 (at least 2 non-NaN values)
        valid_rows_1_149 = (sum(~isnan(relevant_data(1:149, :)), 2) >= 2);
        valid_rows_1_149_temp = (sum(~isnan(relevant_temp_data(1:149, :)), 2) >= 2);
        % Check validity of rows 100:end (at least 3 non-NaN values)
        valid_rows_150_end = (sum(~isnan(relevant_data(150:end, :)), 2) >= 3);
        valid_rows_150_end_temp = (sum(~isnan(relevant_temp_data(150:end, :)), 2) >= 3);
        % Combine row validity
        valid_rows = [valid_rows_1_149; valid_rows_150_end];       
        valid_temp_rows = [valid_rows_1_149_temp; valid_rows_150_end_temp];     
        % Calculate salinity anomalies only for valid depths
        if any(valid_rows)
            % Mean and std of valid depths (row-wise), ignoring NaNs
            mean_values = mean(relevant_data(valid_rows, :), 2, 'omitnan');    
            std_values = std(relevant_data(valid_rows,:), 2, 'omitnan') ;
            % Subtract the mean from fjord_sal_mat_fj for valid depths
            fj_anoms{i}(valid_rows, j) = fjord_sal_mat_fj(valid_rows, fj_profile_back{i}(j)) - mean_values;
            fj_z_anoms{i}(valid_rows,j) =fj_anoms{i}(valid_rows,j) ./ std_values ;
        end
        % same for temp anoms
        if any(valid_rows)
            % Mean of valid depths (row-wise), ignoring NaNs
            mean_temp_values = mean(relevant_temp_data(valid_temp_rows, :), 2, 'omitnan');            
            % Subtract the mean from fjord_sal_mat_fj for valid depths
            fj_temp_anoms{i}(valid_temp_rows, j) = fjord_temp_mat_fj(valid_temp_rows, fj_profile_back{i}(j)) - mean_temp_values;
        end
    end
end
%make into a mat
valid_anoms = fj_anoms(~cellfun('isempty', fj_anoms));
fj_anoms_mat = horzcat(valid_anoms{:});
valid_temp_anoms = fj_temp_anoms(~cellfun('isempty', fj_temp_anoms));
fj_temp_anoms_mat = horzcat(valid_temp_anoms{:});
% derivative cleaning for anomalies (only replace failed points with NaN)
diff_result = NaN(size(fj_anoms_mat)) ;
    for i  = 1:length(fj_anoms_mat)
    non_nan = find(~isnan(fj_anoms_mat(:,i))) ;
    non_nan_store{i} = non_nan ;
    temporary = fj_anoms_mat(:,i) ;
    diff_sal = diff(temporary(non_nan),1,1) ;
    diff_result(non_nan(1:end-1),i) = abs(diff_sal ./ diff(non_nan,1,1)) ;
    end
   threshold_25 = 12 ; % 0-25 m
   threshold_50 = 2 ; % 0-50 m
   threshold = 0.1 ; % 50-end m
     top_25 = diff_result(1:25,:) ;
     top_50 = diff_result(26:90,:) ;
     bottom = diff_result(91:end,:) ;
     top_25_remove = top_25 >= threshold_25 ; % index of top 25 m that has derivatives greater than threshold
     top_50_remove = top_50 >= threshold_50 ; % idx of bottom that has derivatives greater than threshold
     bottom_remove = bottom >= threshold ; % idx of bottom that has derivatives greater than threshold
     remove_anom = [top_25_remove;top_50_remove;bottom_remove] ; % need to account for both points that create a true remove value
     for i = 1:length(remove_anom) 
     idx = find(remove_anom(:,i) == 1) ;
        for j = 1:length(idx)
            if size(idx) >= 1
        value_below_idx = find(non_nan_store{i} == idx(j,:)) + 1 ;
        value_below = non_nan_store{i}(value_below_idx) ;
        remove_anom(value_below, i) = 1; % Mark the value from non_nan_store for removal
            end
        end
     end
%remove those anomaly points that fail the threshold
%remove_anom = logical(remove_anom) ;
fj_anoms_mat(remove_anom) = NaN ;
fj_temp_anoms_mat(remove_anom) = NaN ;
save fj_anoms_mat.mat fj_anoms_mat
save fj_temp_anoms_mat.mat fj_temp_anoms_mat
save fj_coords.mat fj_coords
end
load fj_anoms_mat.mat 
load fj_coords.mat
load fj_temp_anoms_mat.mat
%valid_anoms = fj_anoms(~cellfun('isempty', fj_anoms));
%fj_anoms_mat = horzcat(valid_anoms{:});
clear interp_sal interp_temp fj_profile_back  valid_rows valid_rows_1_149 valid_rows_150_end mean_values non_nan diff_result threshold_50 bottom top_50_remove %fj_anom_idx
clear bottom_remove remove_anom idx value_below_idx value_below top_25 top_25_remove valid_temp_anoms valid_rows_150_end_temp valid_rows_1_149_temp
clear mean_temp_values threshold threshold_25 too_many_nans_box too_many_nans_fj top_100_coast top_50 OMG_fj_profiles OMG_box_profiles
clear box_find_combined box_sal_mat date datenum day_box fjord_box_cords in lat_box lon_box lat_temp remove remove_box remove_coast % may need these
clear unique_col_idx_open yea_box yea box_temp_mat box_find exten2 sal_mat_fj
clear XX_eto YY_eto Z_eto XX_eto_E YY_eto_E Z_eto_E XX_eto_N YY_eto_N Z_eto_N
%% Potential Density and Spicity (needs to be rerun with higher optimization)
run = 2 ;
if run == 1
% Pressure (db) from Depth and Density (can clear after getting potential temp)
numpoints = length(DepInterval) ;
for i = 1:length(lat_fj)
k = sw_pres(DepInterval,repmat(lat_fj(i),1,numpoints)) ;
press_fj(:,i) = k ;
end
for i = 1:length(lat_a)
k = sw_pres(DepInterval,repmat(lat_a(i),1,numpoints)) ;
press_a(:,i) = k ;
end    
dens_fj = sw_dens(fjord_sal_mat_fj,fjord_temp_mat_fj,press_fj) ; % kg/m^3
dens_a = sw_dens(interp_sal_mat,interp_temp_mat,press_a) ; % kg/m^3
% Vectorize and combine
clear  k coastal_temp open_temp interp_temp_a interp_sal_a dens_a dens_fj %interp_sal_mat interp_temp_mat
% Split data_a into two chunks to proccess seperatly for more efficient memory usage
half = round(size(interp_sal_mat,2)/2) ; % halfway point
sal_1 = interp_sal_mat(:,1:half) ;
temp_1 = interp_temp_mat(:,1:half) ;
press_1 = press_a(:,1:half) ;
sal_2 = interp_sal_mat(:,half+1:end) ;
temp_2 = interp_temp_mat(:,half+1:end) ;
press_2 = press_a(:,half+1:end) ;
%temporarily clear the non-split variables for space, will recombined when combining spice
clear interp_sal_mat interp_temp_mat press_a half
% fjord spicity (get profiles of spicity for each fjord profile (using sw_pspi())
spice_fj = sw_pspi(fjord_sal_mat_fj,fjord_temp_mat_fj,press_fj,0) ; % last # is pressure refference (kg/m^3)
spice_fj = double(spice_fj) ;
run = 1 ;
spice_1 = sw_pspi(sal_1,temp_1,press_1,0) ; % last # is pressure refference (kg/m^3)
spice_2 = sw_pspi(sal_2,temp_2,press_2,0) ; % last # is pressure refference (kg/m^3)
spice_a = [spice_1,spice_2] ;
save('spice_a.mat', 'spice_a', '-v7.3');
interp_sal_mat = [sal_1,sal_2] ;
interp_temp_mat = [temp_1,temp_2] ;
end
%load spice_a.mat spice_a
clear press_a press_fj spice_1 spice_2 sal_1 sal_2 press_1 press_1
%%
%top 300 m of both open and coastal
sal_combined = [open_sal,coastal_sal] ;
coastal_sal = coastal_sal(1:300,:) ;
coastal_temp = coastal_temp(1:300,:) ;
open_sal = open_sal(1:300,:) ;
open_temp = open_temp(1:300,:) ;
length_open = length(open_sal) ;
sal_combined = sal_combined(1:300,:) ;
fj_anom_combined = fj_anoms_mat(1:300,:) ;
fj_temp_anom_combined = fj_temp_anoms_mat(1:300,:) ;
fj_combined = fjord_sal_mat_fj(1:300,:) ;
fj_temp_combined = fjord_temp_mat_fj(1:300,:) ;
%May_a = mon_a == 5 ;
%May_fj = mon_fj == 5 ;
%Jun_a = mon_a == 6 ;
%Jun_fj = mon_fj == 6 ;
%Jul_a = mon_a == 7 ;
%Jul_fj = mon_fj == 7 ;
%Aug_a = mon_a == 8 ;
%5Aug_fj = mon_fj == 8 ;
%Sep_a = mon_a == 9 ;
%Sep_fj = mon_fj == 9 ;
%Oct_a = mon_a == 10 ;
%Oct_fj = mon_fj == 10 ;

sal_anom_combined = [open_sal_anom,coast_sal_anom] ;
sal_anom_combined = sal_anom_combined(1:300, :);
sal_anom_combined = [sal_anom_combined,fj_anom_combined] ;
lat_combined = [lat_open,coastal_lat,lat_fj] ; 
lon_combined = [lon_open,coastal_lon,lon_fj ] ;

temp_anom_combined = [open_temp_anom,coast_temp_anom] ;
temp_anom_combined = temp_anom_combined(1:300, :);
temp_anom_combined = [temp_anom_combined,fj_temp_anom_combined] ;
% extropolate and replace NaN from 10 m to surface (avg of top three values
% or lacking that top value gets repeated to surface)
for i = 1:size(sal_anom_combined, 2)
    top10 = sal_anom_combined(1:10, i);  % Extract the top 10 rows of the column
    nonNaN_values = top10(~isnan(top10));  % Extract non-NaN values
    %if length(nonNaN_values) >= 3
        % Average the top three non-NaN values
        %avg_value = mean(nonNaN_values(1:3));
        %sal_anom_combined(1:10, i) = fillmissing(top10, 'constant', avg_value);
    if length(nonNaN_values) >= 1
        % If fewer than 3 non-NaN values, use the top-most non-NaN value
        top_most_value = nonNaN_values(1);
        sal_anom_combined(1:10, i) = fillmissing(top10, 'constant', top_most_value);
    end
end
% same for fj anomalies
for i = 1:size(fj_anom_combined, 2)
    top10 = fj_anom_combined(1:10, i);  % Extract the top 10 rows of the column
    top10_temp = fj_temp_anom_combined(1:10, i);
    nonNaN_values = top10(~isnan(top10));  % Extract non-NaN values
    nonNaN_temp_values = top10_temp(~isnan(top10_temp)) ;
    %if length(nonNaN_values) >= 3
        % Average the top three non-NaN values
        %avg_value = mean(nonNaN_values(1:3));
        %sal_anom_combined(1:10, i) = fillmissing(top10, 'constant', avg_value);
    if length(nonNaN_values) >= 1
        % use the top-most non-NaN value
        top_most_value = nonNaN_values(1);
        fj_anom_combined(1:10, i) = fillmissing(top10, 'constant', top_most_value);
    end
    if length(nonNaN_temp_values) >=1
    top_most_value_temp = nonNaN_temp_values(1) ;
    fj_temp_anom_combined(1:10,i) = fillmissing(top10_temp, 'constant', top_most_value_temp);
    end
end
% same for open salinity 
for i = 1:size(open_sal, 2)
    top10 = open_sal(1:10, i);  % Extract the top 10 rows of the column
    top10_temp = open_temp(1:10,i) ;
    nonNaN_values = top10(~isnan(top10));  % Extract non-NaN values
    nonNaN_temp = top10_temp(~isnan(top10_temp)) ;
    %if length(nonNaN_values) >= 3
        % Average the top three non-NaN values
        %avg_value = mean(nonNaN_values(1:3));
        %open_sal(1:10, i) = fillmissing(top10, 'constant', avg_value);
    if length(nonNaN_values) >= 1
        % If fewer than 3 non-NaN values, use the top-most non-NaN value
        top_most_value = nonNaN_values(1);
        open_sal(1:10, i) = fillmissing(top10, 'constant', top_most_value);
    end
    if length(nonNaN_temp) >= 1
        top_most_value = nonNaN_temp(1) ;
        open_temp(1:10,i) = fillmissing(top10_temp,'constant',top_most_value) ;
    end
end
% same for coastal_salinity
for i = 1:size(coastal_sal, 2)
    top10 = coastal_sal(1:10, i);  % Extract the top 10 rows of the column
    top10_temp = coastal_temp(1:10,i) ;
    nonNaN_values = top10(~isnan(top10));  % Extract non-NaN values
    nonNaN_temp = top10_temp(~isnan(top10_temp)) ;
    %if length(nonNaN_values) >= 3
        % Average the top three non-NaN values
        %avg_value = mean(nonNaN_values(1:3));
        %coastal_sal(1:10, i) = fillmissing(top10, 'constant', avg_value);
    if length(nonNaN_values) >= 1
        % If fewer than 3 non-NaN values, use the top-most non-NaN value
        top_most_value = nonNaN_values(1);
        coastal_sal(1:10, i) = fillmissing(top10, 'constant', top_most_value);
    end
    if length(nonNaN_values) >= 1
        top_most_value = nonNaN_values(1) ;
        coastal_temp(1:10,i) = fillmissing(top10_temp,'constant',top_most_value) ;
    end
end
% Same for fjord salinity
for i = 1:size(fj_combined, 2)
    top10 = fj_combined(1:10, i);  % Extract the top 10 rows of the column
    top10_temp = fj_temp_combined(1:10,i) ;
    nonNaN_values = top10(~isnan(top10));  % Extract non-NaN values
    nonNaN_temp_values = top10_temp(~isnan(top10_temp)) ;
    if length(nonNaN_values) >= 1
        % If fewer than 3 non-NaN values, use the top-most non-NaN value
        top_most_value = nonNaN_values(1);
        fj_combined(1:10, i) = fillmissing(top10, 'constant', top_most_value);
    end
    if length(nonNaN_temp_values) >=1
    top_most_value_temp = nonNaN_temp_values(1);
    fj_temp_combined(1:10, i) = fillmissing(top10_temp, 'constant', top_most_value_temp);
    end
end

% combined everything for another PCA  
open_lat = lat_a(~in_a) ; %unaltered length
open_lon = lon_a(~in_a) ;
lat_coastal = lat_a(in_a) ;
lon_coastal =  lon_a(in_a) ;

all_sal = [open_sal,coastal_sal,fj_combined] ;
all_temp = [open_temp,coastal_temp,fj_temp_combined] ;
all_lon = [open_lon,lon_coastal,lon_fj] ;
all_lat = [open_lat,lat_coastal,lat_fj] ;

clear top_most_value top10 nonNaN_values avg_value open_sal_anom coast_sal_anom coast_temp_anom open_temp_anom
clear top_most_value_temp top10_temp nonNaN_temp_values nonNaN_temp
%% Canadian area selection (not based off anything? maybe countour bathy to verfiy)
% Defining Coastline
run = 2 ;
for i = 1:1:1
    if run == 1
clf
hold on
plot(cx,cy,'k')
scatter(lon_open(:,year_mon_open), lat_open(:,year_mon_open));
xlim([-80,-50])
ylim([50,77])
daspect([1 aspect_ratio 1])
x_canada = [];
y_canada = [];
run = true;
while run
    [x, y, button] = ginput(1);
    if isempty(button) || button == 13
        run = false;
    else
        plot(x, y, 'ro'); % Plot the point
        xlim([x-10, x+10])
        ylim([y-10, y+10])
        x_canada = [x_canada; x];
        y_canada = [y_canada; y];
    end
end
save("x_canada.mat","x_canada")
save("y_canada.mat","y_canada")
    end
end
load("x_canada.mat")
load("y_canada.mat")
canada = inpolygon(lon_open,lat_open,x_canada,y_canada) ;
copy = open_sal ;
%% Invert
month_selected = 9 ; % for sprintf
year_selected = 2012 ; % for sprintf, replace with desired year
five_year_range = [year_selected - 0, year_selected + 0] ; % for scarcer fjord data
% Invert for PCA
if size(sal_anom_combined,2) > 301
sal_anom_combined = sal_anom_combined' ;
can_invert = canada' ;
end
if size(coastal_sal,2) > 301
coastal_sal = coastal_sal' ;
open_sal = open_sal' ;
all_temp = all_temp' ;
all_sal = all_sal' ;
temp_anom_combined = temp_anom_combined' ;
sal_anom_combined = sal_anom_combined' ;
end

% remove all nan columns for fj_combined
%valid_cols_idx = any(~isnan(fj_combined), 1); % use this to create lon/lat for anomalies
%fj_combined = fj_combined(:, valid_cols_idx);

if size(fj_combined,2) > 301
fj_combined = fj_combined' ;
fj_anom_combined = fj_anom_combined';
fj_temp_combined = fj_temp_combined' ;
fj_temp_anom_combined = fj_temp_anom_combined' ;
%spice_fj = spice_fj' ;
end
% Fj_combined single to double
fj_combined = double(fj_combined) ;
fj_temp_combined = double(fj_temp_combined) ;
fj_temp_anom_combined = double(fj_temp_combined) ;
all_sal = double(all_sal) ;
all_temp = double(all_temp) ;
% create variables for non-anomaly use
lat_coast_n = lat_a(in_a) ;
lon_coast_n = lon_a(in_a) ;
lat_open_n = lat_a(~in_a) ;
lon_open_n = lon_a(~in_a) ;
% month and year idx's
coastal_mon = mon_a(in_a) ;
coastal_mon = coastal_mon(coast_idx) ;
open_mon = mon_a(~in_a) ;
open_mon = open_mon(open_idx);
%can_mon = open_mon(canada) ;
%open_mon = open_mon(~canada) ;
coastal_yea = yea_a(in_a) ;
coastal_yea = coastal_yea(coast_idx) ;
open_yea = yea_a(~in_a) ;
open_yea = open_yea(open_idx) ;
%can_yea = open_yea(canada) ;
%open_yea = open_yea(~canada) ;
all_mon = [open_mon,coastal_mon,mon_fj] ;
all_yea = [open_yea,coastal_yea,yea_fj] ;
anom_mon = [];
anom_yea = [];
% month
month_string = {'January', 'February', 'March', 'April', 'May', 'June','July', 'August', 'September', 'October', 'November', 'December'} ;
mon_comb = [open_mon,coastal_mon] ;
month_s = mon_comb == month_selected ; % # of desired month;
month_s_fj = mon_fj == month_selected ; % don't need seperate one for anomaly since i didn't remove them, just set to NAN
month_open = open_mon == month_selected ;
month_coast = coastal_mon == month_selected ;
month_all = all_mon == month_selected ;
% Salinity Month
month_coast_n_test = mon_a(in_a) ;
month_open_n_test = mon_a(~in_a) ;
month_coast_n = month_coast_n_test == month_selected ;
month_open_n = month_open_n_test == month_selected ;
% year
yea_combined = [open_yea, coastal_yea] ;
year = yea_combined == year_selected ;
year_open = open_yea == year_selected ;
year_coast = coastal_yea == year_selected ;
five_yea_s_fj = yea_fj >= five_year_range(1,1) & yea_fj <= five_year_range(1,2) ; % don't need seperate
year_all = all_yea == year_selected ;
% Salinity Year
year_coast_n_test = yea_a(in_a) ;
year_open_n_test = yea_a(~in_a) ;
year_coast_n = year_coast_n_test == year_selected ;
year_open_n = year_open_n_test == year_selected ;
%combine
year_mon = year & month_s ;
year_mon_open = year_open & month_open ;
year_mon_coast = year_coast & month_coast ;
year_mon_fj = five_yea_s_fj & month_s_fj ;
all_yea_mon = month_all & year_all ;
% Salinity combine
yeamon_n_coast = year_coast_n & month_coast_n ;
yeamon_n_open = year_open_n & month_open_n ;
%clear open_mon coastal_mon %combined_idx_open combined_idx_coast
%copy = sal_anom_combined ; % for copying
%in_combined = inpolygon(lon_combined,lat_combined,combined_x,combined_y) ;
%open_sal_anom = sal_anom_combined(~in_combined,:) ;
%coast_sal_anom = sal_anom_combined(in_combined,:) ;
%lon_coast = lon_combined(in_combined) ;
%lon_open = lon_combined(~in_combined) ;
%lat_coast = lat_combined(in_combined) ;
%lat_open = lat_combined(~in_combined) ;
% Canada
%canada_n = inpolygon(lon_open_n,lat_open_n,x_canada,y_canada) ;
%can_lon_n = lon_open_n(canada_n) ;
%can_lat_n = lat_open_n(canada_n) ;
%lat_open_n = lat_open_n(~canada_n) ;
%lon_open_n = lon_open_n(~canada_n) ;
%can_lon = lon_open(canada) ;
%can_lat = lat_open(canada) ;
%lon_open = lon_open(~canada) ;
%lat_open = lat_open(~canada) ;
%yeamon_n_open = yeamon_n_open(~canada_n) ;
% sal/anom/temp/anom
    %if size(canada,1) == 1
    %can_invert = canada' ;
    %can_n_invert = canada_n' ;
    %end
    %can_sal = open_sal(can_n_invert,:) ; % copy so we can rerun this section as needed
    %can_sal_anom = open_sal_anom(can_invert,:) ;
    %open_sal_anom = open_sal_anom(~can_invert,:) ;
    %open_sal = open_sal(~can_n_invert,:) ;
clear length_open yea year_coast month_coast month_open year_open year_coast_n year_open_n month_coast_n month_open_n month_coast_n_test year_open_n_test year_coast_n_test
clear month_open_n_test canada canada_n can_invert can_n_invert yea_s_fj month_s_fj
clear can_sal can_sal_anom can_invert can_n_invert can_lon can_lat canada_n can_lat_n can_lon_n can_mon can_yea % clearing all canada variables (will get rid of them entirely eventually)
%% Helheim (or other restricted region anom PCA) PCA and GMM
% sal
% location index (helheim is fjord_vert{32}
hel_idx = inpolygon(lon_fj,lat_fj,fjord_vert{32}(:,1),fjord_vert{32}(:,2))' ;
fj_combined_hel = fj_combined(hel_idx,:) ;
fj_combined_test = fj_combined_hel(:,1:100) ; % reduced depths
fj_temp_hel = fj_temp_combined(hel_idx,:) ;
fj_temp_combined_test = fj_temp_hel(:,1:100) ;
run = 1 ;
if run == 1
valid_rows = all(~isnan(fj_combined_test), 2) & all(~isnan(fj_temp_combined_test), 2);
sal = fj_combined_test(valid_rows,:) ; % should just be able to change this
% Find the first column where there are less then 3 non-nan values and
% truncate
last_nan_col = find(sum(~isnan(sal), 1) < 3, 1);
if ~isempty(last_nan_col)
    % Cut off the columns from the first NaN column onwards
    sal = sal(:, 1:last_nan_col-1);
end
mu_sal = mean(sal,1) ;
std_sal = std(sal) ;
sal_norm = (sal-mu_sal)./std_sal ;
 % standard normalization of data
[coeff, score, latent , tsquared] = eof225(sal_norm,NaN,50); % Renato's Function (50 is number he gave) (very slow so reduce NaN's as much as possible)
first_PC = score(:,1) ; % first principal component
first_coeff = coeff(:,1); % first pc coeff
second_PC = score(:,2) ; % second
third_PC = score(:,3) ;
explained = 100 * latent / sum(latent);
clear last_nan_col

coeff_fj_sal = coeff ;
score_fj_sal = score ;
latent_fj_sal = latent ;
tsqaured_fj_sal = tsquared ;
explained_fj_sal = explained ;
end
% Fjord Temp PCA
run = 1 ;
if run == 1
temp = fj_temp_combined_test(valid_rows,:) ; 
% Find the first column where there are less then 3 non-nan values and
% truncate
last_nan_col = find(sum(~isnan(temp), 1) < 3, 1);
if ~isempty(last_nan_col)
    % Cut off the columns from the first NaN column onwards
    temp = temp(:, 1:last_nan_col-1);
end
mu_temp = mean(temp,1) ;
std_temp = std(temp) ;
temp_norm = (temp-mu_temp)./std_temp ;
% standard normalization of data
[coeff, score, latent , tsquared] = eof225(temp_norm,NaN,50); % Renato's Function (50 is number he gave) (very slow so reduce NaN's as much as possible)
first_PC = score(:,1) ; % first principal component
first_coeff = coeff(:,1); % first pc coeff
second_PC = score(:,2) ; % second
third_PC = score(:,3) ;
explained = 100 * latent / sum(latent);
clear last_nan_col

coeff_fj_temp = coeff ;
score_fj_temp = score ;
latent_fj_temp = latent ;
tsqaured_fj_temp = tsquared ;
explained_fj_temp = explained ;
end
% Train GMM on PCA Data (Depth Indepdent)
%depth independent
% select random training data (wasn't working right, do later) 
k = 5 ; % number of clusters
feature_matrix = [score_fj_temp(:,1:3),score_fj_sal(:,1:3)] ;
lat_fj_test = lat_fj(hel_idx) ;
lat_fj_test = lat_fj_test(valid_rows') ;
lon_fj_test = lon_fj(hel_idx) ;
lon_fj_test = lon_fj_test(valid_rows') ;
mon_fj_test = mon_fj(hel_idx) ;
mon_fj_test = mon_fj_test(valid_rows') ;
yea_fj_test = yea_fj(hel_idx) ;
yea_fj_test = yea_fj_test(valid_rows) ;
%Model
run = 1 ;
if run == 1
options = statset('MaxIter', 500, 'Display', 'final');  % Increase iterations to 500
indepen_model = fitgmdist(feature_matrix, k, 'Options', options);
%save indepen_model.mat indepen_model
end
%load indepen_model
cluster_labels = cluster(indepen_model, feature_matrix);
cluster_probs = posterior(indepen_model, feature_matrix);
clear index
%% Helheim Depth Dependent
hel_idx = inpolygon(lon_fj,lat_fj,fjord_vert{32}(:,1),fjord_vert{32}(:,2))' ;
fj_combined_hel = fj_combined(hel_idx,:) ;
fj_combined_test = fj_combined_hel(:,1:100) ; % reduced depths
fj_temp_hel = fj_temp_combined(hel_idx,:) ;
fj_temp_combined_test = fj_temp_hel(:,1:100) ;
num_splits = 5 ; % should be evenly divisible by the total length
disp([num2str(size(fj_combined_test, 2)/2) ' meter segments']); % length of depth segments (should be whole number)

% PCA 
% sal
run = 1 ;
if run == 1
valid_rows = all(~isnan(fj_combined_test), 2) & all(~isnan(fj_temp_combined_test), 2);
sal_init = fj_combined_test(valid_rows,:) ; % should just be able to change this
sal = interleave_matrix(sal_init,num_splits) ; % splits and interleaves collumns for depth dependent PCA

% Find the first column where there are less then 3 non-nan values and
% truncate
last_nan_col = find(sum(~isnan(sal), 1) < 3, 1);
if ~isempty(last_nan_col)
    % Cut off the columns from the first NaN column onwards
    sal = sal(:, 1:last_nan_col-1);
end
mu_sal = mean(sal,1) ;
std_sal = std(sal) ;
sal_norm = (sal-mu_sal)./std_sal ;
[coeff, score, latent , tsquared] = eof225(sal,NaN,50); % Renato's Function (50 is number he gave) (very slow so reduce NaN's as much as possible)
first_PC = score(:,1) ; % first principal component
first_coeff = coeff(:,1); % first pc coeff
second_PC = score(:,2) ; % second
third_PC = score(:,3) ;
explained = 100 * latent / sum(latent);
clear last_nan_col

coeff_fj_sal = coeff ;
score_fj_sal = score ;
latent_fj_sal = latent ;
tsqaured_fj_sal = tsquared ;
explained_fj_sal = explained ;
end
% Fjord Temp PCA
run = 1 ;
if run == 1
temp_init = fj_temp_combined_test(valid_rows,:) ; 
temp = interleave_matrix(temp_init,num_splits) ; % splits and interleaves collumns for depth dependent PCA
% Find the first column where there are less then 3 non-nan values and
% truncate
last_nan_col = find(sum(~isnan(temp), 1) < 3, 1);
if ~isempty(last_nan_col)
    % Cut off the columns from the first NaN column onwards
    temp = temp(:, 1:last_nan_col-1);
end
mu_temp = mean(temp,1) ;
std_temp = std(temp) ;
temp_norm = (temp-mu_temp)./std_temp ;
[coeff, score, latent , tsquared] = eof225(temp,NaN,50); % Renato's Function (50 is number he gave) (very slow so reduce NaN's as much as possible)
first_PC = score(:,1) ; % first principal component
first_coeff = coeff(:,1); % first pc coeff
second_PC = score(:,2) ; % second
third_PC = score(:,3) ;
explained = 100 * latent / sum(latent);
clear last_nan_col

coeff_fj_temp = coeff ;
score_fj_temp = score ;
latent_fj_temp = latent ;
tsqaured_fj_temp = tsquared ;
explained_fj_temp = explained ;
end
% GMM Training (Depth Dependent)
run = 1 ;
if run == 1
%depth independent
% select random training data (wasn't working right, do later) 
k = 3 ; % number of clusters
feature_matrix = [score_fj_temp(:,1:3),score_fj_sal(:,1:3)] ;
lat_fj_test = lat_fj(hel_idx) ;
lat_fj_test = lat_fj_test(valid_rows') ;
lon_fj_test = lon_fj(hel_idx) ;
lon_fj_test = lon_fj_test(valid_rows') ;
mon_fj_test = mon_fj(hel_idx) ;
mon_fj_test = mon_fj_test(valid_rows') ;
yea_fj_test = yea_fj(hel_idx) ;
yea_fj_test = yea_fj_test(valid_rows) ;
%Model
options = statset('MaxIter', 1000, 'Display', 'final');  % Increase iterations to 500
depth_model = fitgmdist(feature_matrix, k, 'Options', options);
clear index
save depth_model.mat depth_model
end
load depth_model.mat depth_model
lat_fj_test = lat_fj(hel_idx) ;
lat_fj_test = lat_fj_test(valid_rows') ;
lon_fj_test = lon_fj(hel_idx) ;
lon_fj_test = lon_fj_test(valid_rows') ;
mon_fj_test = mon_fj(hel_idx) ;
mon_fj_test = mon_fj_test(valid_rows') ;
yea_fj_test = yea_fj(hel_idx) ;
yea_fj_test = yea_fj_test(valid_rows) ;
feature_matrix = [score_fj_temp(:,1:3),score_fj_sal(:,1:3)] ;
cluster_labels = cluster(depth_model, feature_matrix);
cluster_probs = posterior(depth_model, feature_matrix);

%% Fjord sal anom PCA
run = 2 ;
if run == 1
sal_anom = fj_anom_combined ; % should just be able to change this
% Find the first column where there are less then 3 non-nan values and
% truncate
last_nan_col = find(sum(~isnan(sal_anom), 1) < 3, 1);
if ~isempty(last_nan_col)
    % Cut off the columns from the first NaN column onwards
    sal_anom = sal_anom(:, 1:last_nan_col-1);
end
[coeff, score, latent , tsquared] = eof225(sal_anom,NaN,50); % Renato's Function (50 is number he gave) (very slow so reduce NaN's as much as possible)
first_PC_anom = score(:,1) ; % first principal component
first_coeff_anom = coeff(:,1); % first pc coeff
second_PC_anom = score(:,2) ; % second
third_PC_anom = score(:,3) ;
explained_anom = 100 * latent / sum(latent);
clear last_nan_col

coeff_fj_anom = coeff ;
score_fj_anom = score ;
latent_fj_anom = latent ;
tsqaured_fj_anom = tsquared ;
explained_fj_anom = explained_anom ;
save fj_anom_PCA.mat coeff_fj_anom score_fj_anom latent_fj_anom tsqaured_fj_anom explained_fj_anom
end
load fj_anom_PCA.mat coeff_fj_anom score_fj_anom latent_fj_anom tsqaured_fj_anom
% Fjord Temp anom PCA
run = 2 ;
if run == 1
temp_anom = fj_temp_anom_combined ; 
% Find the first column where there are less then 3 non-nan values and
% truncate
last_nan_col = find(sum(~isnan(temp_anom), 1) < 3, 1);
if ~isempty(last_nan_col)
    % Cut off the columns from the first NaN column onwards
    temp_anom = temp_anom(:, 1:last_nan_col-1);
end
[coeff, score, latent , tsquared] = eof225(temp_anom,NaN,50); % Renato's Function (50 is number he gave) (very slow so reduce NaN's as much as possible)
first_PC_anom = score(:,1) ; % first principal component
first_coeff_anom = coeff(:,1); % first pc coeff
second_PC_anom = score(:,2) ; % second
third_PC_anom = score(:,3) ;
explained_anom = 100 * latent / sum(latent);
clear last_nan_col

coeff_fj_anom_temp = coeff ;
score_fj_anom_temp = score ;
latent_fj_anom_temp = latent ;
tsqaured_fj_anom_temp = tsquared ;
explained_fj_anom_temp = explained_anom ;
save fj_anom_temp_PCA.mat coeff_fj_anom_temp score_fj_anom_temp latent_fj_anom_temp tsqaured_fj_anom_temp explained_fj_anom_temp
end
load fj_anom_temp_PCA.mat coeff_fj_anom_temp score_fj_anom_temp latent_fj_anom_temp tsqaured_fj_anom_temp explained_fj_anom_temp
%% Test PCA's only looking at stuff with measurements to 100 m) for both temp and sal
% sal
fj_combined_test = fj_combined(:,1:100) ; % reduced depths
fj_temp_combined_test = fj_temp_combined(:,1:100) ;
run = 1 ;
if run == 1
valid_rows = all(~isnan(fj_combined_test), 2) & all(~isnan(fj_temp_combined_test), 2);
sal = fj_combined_test(valid_rows,:) ; % should just be able to change this
% Find the first column where there are less then 3 non-nan values and
% truncate
last_nan_col = find(sum(~isnan(sal), 1) < 3, 1);
if ~isempty(last_nan_col)
    % Cut off the columns from the first NaN column onwards
    sal = sal(:, 1:last_nan_col-1);
end
mu_sal = mean(sal,1) ;
std_sal = std(sal) ;
sal_norm = (sal-mu_sal)./std_sal ;
 % standard normalization of data
[coeff, score, latent , tsquared] = eof225(sal_norm,NaN,50); % Renato's Function (50 is number he gave) (very slow so reduce NaN's as much as possible)
first_PC = score(:,1) ; % first principal component
first_coeff = coeff(:,1); % first pc coeff
second_PC = score(:,2) ; % second
third_PC = score(:,3) ;
explained = 100 * latent / sum(latent);
clear last_nan_col

coeff_fj_sal = coeff ;
score_fj_sal = score ;
latent_fj_sal = latent ;
tsqaured_fj_sal = tsquared ;
explained_fj_sal = explained ;
end
% Fjord Temp PCA
run = 1 ;
if run == 1
temp = fj_temp_combined_test(valid_rows,:) ; 
% Find the first column where there are less then 3 non-nan values and
% truncate
last_nan_col = find(sum(~isnan(temp), 1) < 3, 1);
if ~isempty(last_nan_col)
    % Cut off the columns from the first NaN column onwards
    temp = temp(:, 1:last_nan_col-1);
end
mu_temp = mean(temp,1) ;
std_temp = std(temp) ;
temp_norm = (temp-mu_temp)./std_temp ;
% standard normalization of data
[coeff, score, latent , tsquared] = eof225(temp_norm,NaN,50); % Renato's Function (50 is number he gave) (very slow so reduce NaN's as much as possible)
first_PC = score(:,1) ; % first principal component
first_coeff = coeff(:,1); % first pc coeff
second_PC = score(:,2) ; % second
third_PC = score(:,3) ;
explained = 100 * latent / sum(latent);
clear last_nan_col

coeff_fj_temp = coeff ;
score_fj_temp = score ;
latent_fj_temp = latent ;
tsqaured_fj_temp = tsquared ;
explained_fj_temp = explained ;
end
% Train GMM on PCA Data (Depth Indepdent)
%depth independent
% select random training data (wasn't working right, do later) 
k = 8 ; % number of clusters
feature_matrix = [score_fj_temp(:,1:4),score_fj_sal(:,1:4)] ;
lat_fj_test = lat_fj(valid_rows) ;
lon_fj_test = lon_fj(valid_rows) ;
mon_fj_test = mon_fj(valid_rows) ;
yea_fj_test = yea_fj(valid_rows) ;
%Model
run = 1 ;
if run == 1
options = statset('MaxIter', 500, 'Display', 'final');  % Increase iterations to 500
indepen_model = fitgmdist(feature_matrix, k, 'Options', options);
save indepen_model.mat indepen_model
end
load indepen_model
cluster_labels = cluster(indepen_model, feature_matrix);
cluster_probs = posterior(indepen_model, feature_matrix);
clear index
%% Depth Dependent Setup, PCA, and GMM training
%setup (split and concatenate profiles by a selected depth)
fj_combined_test = fj_combined(:,1:100) ; % reduced depths
fj_temp_combined_test = fj_temp_combined(:,1:100) ;
num_splits = 5 ; % should be evenly divisible by the total length
disp([num2str(size(fj_combined_test, 2)/2) ' meter segments']); % length of depth segments (should be whole number)

% PCA 
% sal
run = 1 ;
if run == 1
valid_rows = all(~isnan(fj_combined_test), 2) & all(~isnan(fj_temp_combined_test), 2);
sal_init = fj_combined_test(valid_rows,:) ; % should just be able to change this
sal = interleave_matrix(sal_init,num_splits) ; % splits and interleaves collumns for depth dependent PCA

% Find the first column where there are less then 3 non-nan values and
% truncate
last_nan_col = find(sum(~isnan(sal), 1) < 3, 1);
if ~isempty(last_nan_col)
    % Cut off the columns from the first NaN column onwards
    sal = sal(:, 1:last_nan_col-1);
end
mu_sal = mean(sal,1) ;
std_sal = std(sal) ;
sal_norm = (sal-mu_sal)./std_sal ;
[coeff, score, latent , tsquared] = eof225(sal,NaN,50); % Renato's Function (50 is number he gave) (very slow so reduce NaN's as much as possible)
first_PC = score(:,1) ; % first principal component
first_coeff = coeff(:,1); % first pc coeff
second_PC = score(:,2) ; % second
third_PC = score(:,3) ;
explained = 100 * latent / sum(latent);
clear last_nan_col

coeff_fj_sal = coeff ;
score_fj_sal = score ;
latent_fj_sal = latent ;
tsqaured_fj_sal = tsquared ;
explained_fj_sal = explained ;
end
% Fjord Temp PCA
run = 1 ;
if run == 1
temp_init = fj_temp_combined_test(valid_rows,:) ; 
temp = interleave_matrix(temp_init,num_splits) ; % splits and interleaves collumns for depth dependent PCA
% Find the first column where there are less then 3 non-nan values and
% truncate
last_nan_col = find(sum(~isnan(temp), 1) < 3, 1);
if ~isempty(last_nan_col)
    % Cut off the columns from the first NaN column onwards
    temp = temp(:, 1:last_nan_col-1);
end
mu_temp = mean(temp,1) ;
std_temp = std(temp) ;
temp_norm = (temp-mu_temp)./std_temp ;
[coeff, score, latent , tsquared] = eof225(temp,NaN,50); % Renato's Function (50 is number he gave) (very slow so reduce NaN's as much as possible)
first_PC = score(:,1) ; % first principal component
first_coeff = coeff(:,1); % first pc coeff
second_PC = score(:,2) ; % second
third_PC = score(:,3) ;
explained = 100 * latent / sum(latent);
clear last_nan_col

coeff_fj_temp = coeff ;
score_fj_temp = score ;
latent_fj_temp = latent ;
tsqaured_fj_temp = tsquared ;
explained_fj_temp = explained ;
end
% GMM Training (Depth Dependent)
run = 1 ;
if run == 1
%depth independent
% select random training data (wasn't working right, do later) 
k = 5 ; % number of clusters
feature_matrix = [score_fj_temp(:,1:3),score_fj_sal(:,1:3)] ;
lat_fj_test = lat_fj(valid_rows) ;
lon_fj_test = lon_fj(valid_rows) ;
mon_fj_test = mon_fj(valid_rows) ;
yea_fj_test = yea_fj(valid_rows) ;
%Model
options = statset('MaxIter', 1000, 'Display', 'final');  % Increase iterations to 500
depth_model = fitgmdist(feature_matrix, k, 'Options', options);
clear index
save depth_model.mat depth_model
end
load depth_model.mat depth_model
lat_fj_test = lat_fj(valid_rows) ;
lon_fj_test = lon_fj(valid_rows) ;
mon_fj_test = mon_fj(valid_rows) ;
yea_fj_test = yea_fj(valid_rows) ;
feature_matrix = [score_fj_temp(:,1:3),score_fj_sal(:,1:3)] ;
cluster_labels = cluster(depth_model, feature_matrix);
cluster_probs = posterior(depth_model, feature_matrix);
%% PCA change month/index as desired
sal_anom = open_sal_anom(year_mon_open, :); % should just be able to change this
% Find the first column where there are less then 3 non-nan values and
% truncate
last_nan_col = find(sum(~isnan(sal_anom), 1) < 3, 1);
if ~isempty(last_nan_col)
    % Cut off the columns from the first NaN column onwards
    sal_anom = sal_anom(:, 1:last_nan_col-1);
end
[coeff, score, latent , tsquared] = eof225(sal_anom,NaN,50); % Renato's Function (50 is number he gave) (very slow so reduce NaN's as much as possible)
first_PC_anom = score(:,1) ; % first principal component
first_coeff_anom = coeff(:,1); % first pc coeff
second_PC_anom = score(:,2) ; % second
third_PC_anom = score(:,3) ;
explained_anom = 100 * latent / sum(latent);
clear last_nan_col
%% coastal raw salinity data 
sal = coastal_sal(yeamon_n_coast,:) ; % change this to open/coast
mean_sal = nanmean(sal, 1); % mean ignoring NaNs
sal_minus = (sal - mean_sal) ;
% Find the first column where there are less then 3 non-nan values and
% truncate
last_nan_col = find(sum(~isnan(sal_minus), 1) < 3, 1);
if ~isempty(last_nan_col)
    % Cut off the columns from the first NaN column onwards
    sal_minus = sal_minus(:, 1:last_nan_col-1);
end
[coeff_coast_sal, score_coast_sal, latent_coast_sal , ~] = eof225(sal_minus,NaN,50); % Renato's Function 
first_PC_coast_sal = score_coast_sal(:,1) ; % first principal component
first_coeff_coast_sal = coeff_coast_sal(:,1); % first pc coeff
second_PC_coast_sal = score_coast_sal(:,2) ; % second
third_PC_coast_sal = score_coast_sal(:,3) ;
explained_coast_sal = 100 * latent_coast_sal / sum(latent_coast_sal);
PC_coast_lon = lon_coast_n(yeamon_n_coast) ;
PC_coast_lat = lat_coast_n(yeamon_n_coast) ;
PC_coast_sal = coastal_sal(yeamon_n_coast,:) ;
clear last_nan_col
%
%% open raw salinity data 
sal = open_sal(yeamon_n_open,:) ; 
mean_sal = nanmean(sal, 1); % mean ignoring NaNs
sal_minus = (sal - mean_sal) ;
% Find the first column where there are less then 3 non-nan values and
% truncate
last_nan_col = find(sum(~isnan(sal_minus), 1) < 3, 1);
if ~isempty(last_nan_col)
    % Cut off the columns from the first NaN column onwards
    sal_minus = sal_minus(:, 1:last_nan_col-1);
end
[coeff_open_sal, score_open_sal, latent_open_sal , ~] = eof225(sal_minus,NaN,50); % Renato's Function 
first_PC_open_sal = score_open_sal(:,1) ; % first principal component
first_coeff_open_sal = coeff_open_sal(:,1); % first pc coeff
second_PC_open_sal = score_open_sal(:,2) ; % second
third_PC_open_sal = score_open_sal(:,3) ;
explained_open_sal = 100 * latent_open_sal / sum(latent_open_sal);
% get coords for plotting
PC_open_lon = lon_open_n(yeamon_n_open) ;
PC_open_lat = lat_open_n(yeamon_n_open) ;
PC_open_sal = open_sal(yeamon_n_open,:) ;
clear last_nan_col
%% Fjord salinity PCA
sal = fj_combined(year_mon_fj,:) ; % five year span centered around year_selected
mean_sal = nanmean(sal, 1); % mean ignoring NaNs
sal_minus = (sal - mean_sal) ;
% Find the first column where there are less then 3 non-nan values and truncate
last_nan_col = find(sum(~isnan(sal_minus), 1) < 3, 1);
if ~isempty(last_nan_col)
    % Cut off the columns from the first NaN column onwards
    sal_minus = sal_minus(:, 1:last_nan_col-1);
end
[coeff_fj, score_fj, latent_fj , ~] = eof225(sal_minus,NaN,50); % Renato's Function
first_PC_fj = score_fj(:,1) ; % first principal component
first_coeff_fj = coeff_fj(:,1); % first pc coeff
second_PC_fj = score_fj(:,2) ;
third_PC_fjk = score_fj(:,3) ;
explained_fj = 100 * latent_fj / sum(latent_fj);
clear last_nan_col
%% Fjord Tempeature PCA
temp = fj_combined(year_mon_fj,:) ; % five year span centered around year_selected
mean_temp = nanmean(temp, 1); % mean ignoring NaNs
temp_minus = (temp - mean_temp) ;
% Find the first column where there are less then 3 non-nan values and truncate
last_nan_col = find(sum(~isnan(temp_minus), 1) < 3, 1);
if ~isempty(last_nan_col)
    % Cut off the columns from the first NaN column onwards
    temp_minus = temp_minus(:, 1:last_nan_col-1);
end
[coeff_fj_t, score_fj_t, latent_fj_t , ~] = eof225(temp_minus,NaN,50); % Renato's Function
first_PC_fj_t = score_fj_t(:,1) ; % first principal component
first_coeff_fj_t = coeff_fj_t(:,1); % first pc coeff
second_PC_fj_t = score_fj_t(:,2) ;
third_PC_fj_t = score_fj_t(:,3) ;
explained_fj_t = 100 * latent_fj_t/ sum(latent_fj_t);
clear last_nan_col temp mean_temp temp_minus
%% Fjord Spiciness PCA
temp = spice_fj(year_mon_fj,:) ; % five year span centered around year_selected
mean_temp = nanmean(temp, 1); % mean ignoring NaNs
temp_minus = (temp - mean_temp) ;
% Find the first column where there are less then 3 non-nan values and truncate
last_nan_col = find(sum(~isnan(temp_minus), 1) < 3, 1);
if ~isempty(last_nan_col)
    % Cut off the columns from the first NaN column onwards
    temp_minus = temp_minus(:, 1:last_nan_col-1);
end
[coeff_fj_s, score_fj_s, latent_fj_s , ~] = eof225(temp_minus,NaN,50); % Renato's Function
first_PC_fj_s = score_fj_s(:,1) ; % first principal component
first_coeff_fj_s = coeff_fj_s(:,1); % first pc coeff
second_PC_fj_s = score_fj_s(:,2) ;
third_PC_fj_s = score_fj_s(:,3) ;
explained_fj_s = 100 * latent_fj_s/ sum(latent_fj_s);
clear last_nan_col temp mean_temp temp_minus
%%
% Potential Temp = [] ;
press_reff = 10 ; % not sure what a good refference would be
ptmp_a = sw_ptmp(interp_sal_mat,interp_temp_mat,press_a,press_reff) ;
clear interp_sal interp_temp numCells rows open_sal_mean open_temp_mean insidearray1 insidearray2 open_sal open_temp k num_points press_reff copy
%%  Statistics
%number of observations for each cast box
for i = 1:length(cast_idx_coast)
temp = lat_a(cast_idx_coast{i}) ;
count_coast(:,i) = numel(temp) ;
end
for i = 1:length(cast_idx_open)
temp = lat_a(cast_idx_open{i}) ;
count_open(:,i) = numel(temp) ;
end
% concatenate coastal_lat/lon with open_lat/lon (refference
% Plots_Greenland.mat to plot this information)
run = 1 ; % do not need to run this most times
for i = 1:1:1
    if run == 1
coastal_lat = lat_a(in_a) ;
coastal_lon =  lon_a(in_a) ;
coastal_mon = mon_a(in_a) ;
open_lat = lat_a(~in_a) ;
open_lon = lon_a(~in_a) ;
open_mon = mon_a(~in_a) ;
lat_comb = [coastal_lat,open_lat] ; % lat with coastal then open
lon_comb = [coastal_lon,open_lon] ; % coastal then open
count_comb = [count_coast,count_open] ; % coastal then open
mon_comb = [coastal_mon,open_mon] ;
    end
end
% month idx's
May = find(mon_comb == 5) ;
Jun = find(mon_comb == 6) ;
Jul = find(mon_comb == 7) ;
Aug = find(mon_comb == 8) ;
Sep = find(mon_comb == 9) ;
Oct = find(mon_comb == 10) ;
clear temp coast_lat coast_lon count_coast count_open coastal_mon open_mon 
clear open_sal_anom open_sal_std open_sal_avg open_temp_std open_temp_avg open_temp_anom coast_sal_anom coastal_sal_std coast_sal_avg coast_temp_anom coast_temp_avg coastal_temp_std % for clearing memory, will need these eventually
%% Boxes for avg counts
increment = 0.1 ; % degrees
run = 1 ; % no need to run every time
for i = 1:1:1
    if run == 1
lat_min = 55 ; % change for different view
lat_max = 80 ;
lon_min = -85 ;
lon_max = -35 ;
x = lon_min:increment:lon_max ;
y = lat_min:increment:lat_max ;
[XX,YY] = meshgrid(x,y) ;
clear lat_min lat_max lon_min lon_max
    end
end
numLat = length(y) - 1;
numLon = length(x) - 1;
for i = 1:numLat
    for j = 1:numLon
         latCorners = [y(i), y(i+1), y(i+1), y(i)];
         lonCorners = [x(j), x(j), x(j+1), x(j+1)];
         boxes{i,j} = [latCorners; lonCorners] ;
    end
end
cast_boxes = boxes(:)' ;
month = Oct ; % change to change month for following code
for i = 1:length(cast_boxes)
in_box = inpolygon(lon_comb(month),lat_comb(month),cast_boxes{i}(2,:),cast_boxes{i}(1,:)) ;
avg = mean(count_comb(in_box)) ;
box_avg(:,i) = avg ; 
end
clear x y latCorners LonCorners XX YY ZZ increment boxes % in_box avg
%% Number of Obersvations at each depth (DELETE)
for i= 1:length(DepInterval)
    for j= 1:length(SW_poly_temp_clean)
        if ~isempty(SW_poly_temp_clean{j}) && size(SW_poly_temp_clean{j}, 1) >= i ;
Depth_obs = nnz(~isnan(SW_poly_temp_clean{j}(i,:))) ;
Depth_obs_cell{j}(i,:) = Depth_obs ;
    end
    end
end
%% Statistics
% Frequency of Counts
for i = 1:numel(SW_poly_temp_cell)
    % Get the current cell array
    subCellArray = SW_poly_temp_cell{i};
    % Count the number of elements in the current cell array
    counts(i) = numel(subCellArray);
end
%Plot Histogram
Mo = mode(counts) ;
avg = mean(counts) ;
Std_counts = std(counts) ;
clf
hold on
histogram(counts, 'BinMethod', 'integers');
xlabel('Number of Casts per box');
ylabel('Frequency');
title(sprintf('Frequency of Casts with %d day intervals and %d x %d km box', day_range_2, width_2,length_2))
hold off
% Plot Example Temperature Anomaly
clf
hold on
plot(SW_temp_Anmly_cell{1,1280},DepInterval)
axis ij
xlabel('Temperature Anomaly')
ylabel('Depth (M)')
title(sprintf('Temperature Anomaly vs Depth %d day interval %d x %d box',day_range_2,width_2,length_2))
hold off
% Plot Example Salinity Anomaly Only with depths with Std devs > min_count
clf
hold on
plot(SW_sal_Anmly_cell{1,1280},DepInterval)
axis ij
xlabel('Salinity Anomaly')
ylabel('Depth (M)')
title(sprintf('Salinity Anomaly vs Depth %d day interval %d x %d box',day_range_2,width_2,length_2))
hold off
% Plot Std Dev Temp
hold on
for i = 1:length(DepInterval)
    if Depth_obs_cell{1,1280}(i,1)>=min_count
plot(Std_temp{1,1280}(i,:),DepInterval)
    end
end
axis ij
xlabel('Temperature Std dev')
ylabel('Depth (M)')
title(sprintf('Temperature Std dev vs Depth %d day interval %d x %d box',day_range_2,width_2,length_2))
hold off
% Std dev Sal Example Station
hold on
for i = 1:length(DepInterval)
    if Depth_obs_cell{1,1280}(i,1)>=min_count
plot(Std_sal{1280}(i,1),DepInterval)
    end
end
axis ij
xlabel('Salinity Std dev')
ylabel('Depth (M)')
title(sprintf('Salinity Std dev vs Depth %d day interval %d x %d box',day_range_2,width_2,length_2))
hold off
%Plot Example Station
clf
hold on
scatter(SW_poly_lon_cell{1280},SW_poly_lat_cell{1280},4,'MarkerFaceColor','b','MarkerEdgeAlpha','0')
scatter(W_lon(1,1280),W_lat(1,1280),4,'MarkerFaceColor','r','MarkerEdgeAlpha','0')
plot(rotvertcell{1280}(1,:),rotvertcell{1280}(2,:), 'b' )
xlim([-80,-35])
ylim([55,80])
plot(cx,cy, 'k')
plot(x_Wcoast,y_SWcoast, 'r')
plot(x5_Wcoast,y_SWcoast,'r')
hold off

% Open Ocean Casts (will need to make to loop, this is just an example)
clf
hold on
plot(cx,cy,'k')
xlim([-80,-35])
ylim([55,80])
run = 2 ;
for i = 1:1:1
    if run == 1 
        Open_Ocean = ginput(4);
        save ("OpeanOcean.mat", "Open_Ocean");
    end
end
load('OpeanOcean.mat','Open_Ocean')
plot(Open_Ocean(:,1),Open_Ocean(:,2))
Open_Casts_idx = inpolygon(lon,lat,Open_Ocean(:,1),Open_Ocean(:,2)) ; % needs to be replaced with all non-coastal casts
scatter(lon(Open_Casts_idx),lat(Open_Casts_idx),'MarkerFaceColor','b','MarkerEdgeAlpha','0')
lon_open = lon(Open_Casts_idx) ;
lat_open = lat(Open_Casts_idx) ;
%Find Distances from central points to all other points in limited box
%(needs work)
for i = 1:length(lat_open)
open_distances = distance(lat_open(342),lon_open(342),lat_open(i),lon_open(i)) ; % will need to be made to loop
open_distance_cell{i} = open_distances ;
end
for i = 1:length(lat_open)
open_distance_km{i} = deg2km(open_distance_cell{i}) ;
end
open_distance = cell2mat(open_distance_km) ;
open_idx = find(open_distance<radius) ; % can then run stats like histogram, means, anomalies, std dev
%Plot to check
aspect_ratio = cosd(62) ; % rough approx
clf
hold on
plot(cx,cy,'k')
xlim([-80,-35])
ylim([55,80])
daspect([1 aspect_ratio 1])
scatter(lon(Open_Casts_idx),lat(Open_Casts_idx),'MarkerFaceColor','b','MarkerEdgeAlpha','0')
scatter(lon_open(open_idx),lat_open(open_idx), 'MarkerFaceColor','r','MarkerEdgeAlpha','0')
scatter(lon_open(342),lat_open(342),'MarkerFaceColor','g','MarkerEdgeAlpha','0')
hold off
