% Load Data
cd('C:\Users\ajs82292\Desktop\Research\Matlab\Source') ;
addpath('seawater','C:\Users\ajs82292\Desktop\Research\Matlab\Source\seawater') ;
load("02cleanNODC.mat")
load("x_coast.mat")
load("y_coast.mat")
%% Commonly Changed Variables for box size and day interval (set run to one
% below if changing rectangle size)
recwidth = 10*2 ; % km doubled as it is 10 "radius" Change value for different sized cast boxes
reclength = 20*2 ; % km doubled as it is 20 "radius" Change value for different sized cast boxes
day_range = 15 ; % adjust this as needed, actual day range is 2*number listed
radius = 20 ; % km circle for opean ocean
min_count = 4 ; % the minimum # of data points needed in order to plot the statistics of a cast (std dev, anomalies ect.)
aspect_ratio = cosd(65) ; % Aspect Ratio at 65 N
%% Defining Coastlines (need to deliniate between west and east facing later on
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
clear('x','y','button')
%%
% West and East facing coasts bases on slope
for i = 1:length(x_coast)-1
W = (y_coast(i+1,1)-y_coast(i,1))/(x_coast(i+1,1)-x_coast(i,1)) ;
slope(i,1) = W ;
inverse(i,1) = -1/W ;
W = y_coast(i,1)-slope(i)*x_coast(i,1) ; 
intercept(i,1 ) = W ;
end
clear('W')
% Split into West/East Facing coastlines
Ecoast_idx = slope >= 0 ;
Wcoast_idx = slope <= 0 ;
x_Ecoast = x_coast(Ecoast_idx) ;
y_Ecoast = y_coast(Ecoast_idx) ;
x_Wcoast = x_coast(Wcoast_idx) ;
y_Wcoast = y_coast(Wcoast_idx) ;
slope_w = slope(Wcoast_idx) ;
inverse_W = inverse(Wcoast_idx) ;
slope_E = slope(Ecoast_idx) ;
inverse_E = inverse(Ecoast_idx) ;
% Creation of Parallelograms on West facing Coast
width_2 = recwidth ; % for dynamic plot titles
length_2 = reclength ; % for dynamic plot titles
recwidth_coast = 10*2 ; % do not change, for defining larger polygon
reclength_coast = 20*2 ; % do not change, for defining larger polygon
lengthdeg = km2deg(reclength) ;
widthdeg_lon = 111.120 .*cosd(lat) ;
widthdeg_lon = recwidth./widthdeg_lon ;
widthdeg_coast = recwidth_coast./widthdeg_lon ; % For larger polygon
lengthdeg_coast = km2deg(reclength_coast) ;% For larger polygon
Westwidthdeg = 111.120* cosd(y_Wcoast) ;
Westwidthdeg = recwidth./Westwidthdeg ;
Westwidthdeg_coast = 111.320* cosd(x_Wcoast) ;
Westwidthdeg_coast = recwidth_coast./Westwidthdeg_coast ;
% East
E_widthdeg = 111.120* cosd(y_Ecoast) ;
E_widthdeg = recwidth./E_widthdeg ;
E_widthdeg_coast = 111.120* cosd(y_Ecoast) ;
E_widthdeg_coast = recwidth_coast./E_widthdeg_coast ;
%%
% create points along lines
dist = 2.5 ; 
num_points = 1000 ; 
% West
for i = 1:length(x_Wcoast)-1
W = linspace(x_Wcoast(i,1), x_Wcoast(i,1) - dist, num_points); % minus longitude for western facing coasts
x_perp_W{i,1} = W' ; 
W =  inverse_W(i) * (x_perp_W{i} - x_Wcoast(i,1)) + y_Wcoast(i,1);
y_perp_W{i,1} = W ;
end
for i = 1:length(x_Wcoast)-1
    for j = 2:length(x_Wcoast)
W = linspace(x_Wcoast(i+1,1), x_Wcoast(i+1,1) - dist, num_points); % minus longitude for western facing coasts
x_perp_W2{i} = W' ; 
W =  inverse_W(i) * (x_perp_W2{i} - x_Wcoast(i+1,1)) + y_Wcoast(i+1,1);
y_perp_W2{i,1} = W ;
    end
end
%East
for i = 1:length(x_Ecoast)-1
E = linspace(x_Ecoast(i,1), x_Ecoast(i,1) + dist, num_points); % minus longitude for western facing coasts
x_perp_E{i} = E' ; 
E =  inverse_E(i) * (x_perp_E{i} - x_Ecoast(i,1)) + y_Ecoast(i,1);
y_perp_E{i,1} = E ;
end
for i = 1:length(x_Ecoast)-1
E = linspace(x_Ecoast(i+1,1), x_Ecoast(i+1,1) + dist, num_points); % minus longitude for western facing coasts
x_perp_E2{i,1} = E' ; 
E =  inverse_E(i) * (x_perp_E2{i} - x_Ecoast(i+1,1)) + y_Ecoast(i+1,1);
y_perp_E2{i,1} = E ;
end
clear('W','E')
%% 
% SW dist to limit coastal area (not done, only have one point, need
% distances for every single coastline point)
for i = 1:length(x_perp_W)
ref_Wx = [x_perp_W{i},repmat(x_Wcoast(i,1),1000,1)];
ref_Wy = [y_perp_W{i},repmat(y_Wcoast(i,1),1000,1)];
ref_W2x = [x_perp_W2{i},repmat(x_Wcoast(i,1),1000,1)];
ref_W2y = [y_perp_W2{i},repmat(y_Wcoast(i,1),1000,1)];
end
for i = 1:length(x_perp_E)
ref_Ex = [x_perp_E{i},repmat(x_Ecoast(i,1),1000,1)];
ref_Ey = [y_perp_E{i},repmat(y_Ecoast(i,1),1000,1)];
ref_E2x = [x_perp_E2{i},repmat(x_Ecoast(i,1),1000,1)];
ref_E2y = [y_perp_E2{i},repmat(y_Ecoast(i,1),1000,1)];
end
% need to do W2,E,E2 as well
for i = 1:1:length(ref_Wy)
    W_dist(i,1) = sw_dist(ref_Wx(i,:),ref_Wy(i,:),'km') ;
    W2_dist(i,1) = sw_dist(ref_W2x(i,:),ref_W2y(i,:),'km') ;
     E_dist(i,1) = sw_dist(ref_Ex(i,:),ref_Ey(i,:),'km') ;
    E2_dist(i,1) = sw_dist(ref_E2x(i,:),ref_E2y(i,:),'km') ;
end
clear('ref_Wx','ref_W2x','ref_Wy','ref_W2y','ref_Ex','ref_E2x','ref_Ey','ref_E2y','x_perp_E2','y_perp_E','x_perp_W','x_perp_W2','y_perp_W','y_perp_E2','y_perp_W2')
%% Find Casts within defined area
% Western Facing Coasts
x5_Wcoast = x_Wcoast - (5*Westwidthdeg_coast)/2 ; % This should not change
x5_Ecoast = Eastgcoast(:,1) + (5*E_widthdeg_coast)/2 ;
y5_Scoast = Southgcoast(:,2) - (5*Southlengthdeg_coast)/2 ;
extendedcoast_W = x5_Wcoast - (Westwidthdeg)/2; % This will change depending on size of cast boxes
extendedcoast_E = x5_Ecoast + (E_widthdeg)/2;
extendedcoast_S = y5_Scoast - (Southlengthdeg_coast)/2 ;
W_xcombined = [x5_Wcoast,x_Wcoast] ; % counterclockwise verticies order
W_ycombined = [y_Wcoast,y_Wcoast] ; % counterclockwise verticies order
E_xcombined = [x5_Ecoast', x_Ecoast'];
E_x ={E_xcombined} ;
E_ycombined = [y_Ecoast',y_Ecoast'] ;
E_y = {E_ycombined} ;
W_xcombined_coast = [extendedcoast_W,x_Wcoast] ;
for i = 1:length(W_xcombined)-1
    W = [W_xcombined(i+1,1),W_xcombined(i,1),W_xcombined(i,2),W_xcombined(i+1,2)] ;
    W_x{i} = W ;
    W = [W_ycombined(i+1,1),W_ycombined(i,1),W_ycombined(i,2),W_ycombined(i+1,2)] ;
    W_y{i} = W ;
end
for i = 1:length(W_xcombined_coast)-1
    W = [W_xcombined_coast(i+1,1),W_xcombined_coast(i,1),W_xcombined_coast(i,2),W_xcombined_coast(i+1,2)] ;
    W_x_coast{i} = W ;
    W = [W_ycombined(i+1,1),W_ycombined(i,1),W_ycombined(i,2),W_ycombined(i+1,2)] ;
    W_y_coast{i} = W ;
end
for i = 1:length(W_xcombined)-1
    W = inpolygon(lon,lat, W_x{i},W_y{i}) ; %idx of all core coastal areas compartmentalized for individual verticies rotation
    W_idx{i} = W ;
end
 for i = 1:length(W_xcombined)-1
    W = inpolygon(lon,lat, W_x_coast{i},W_y_coast{i}) ; %idx of all extended coastal areas (not sure if they need to be compartmentalized)
    W_idx_coast{i} = W ;
 end
W_idx_combined = any(cat(3, W_idx{:}), 3); % combined indicies of each compartment (will be used for a unified coastal idx)
% Eastern Facing Coasts
for i = 1:length(E_x)-1
E_idx = inpolygon(lon,lat, E_x{i},E_y{i}) ; %idx of all core coastal areas compartmentalized for individual verticies rotation
end
%Gap Coasts (intebetween methods)
Gap_plus = ([x5_Wcoast(1,1),Westgcoast(1,2);Southgcoast(1,1),y5_Scoast(1,1)]) ;
Gap_minus = ([Southgcoast(2,1),y5_Scoast(2,1);x5_Ecoast(1,1),Eastgcoast(1,2)]) ;
clear('W','lengthdegcoast','W_x','W_y','width','widthdeg_coast')
%% Isolate variables Lat and Lon for each compartment (Finish this later)
for i = 1:length(W_idx)
    W = lon(W_idx{i}) ;
W_lon{i} = W ;
    W = lat(W_idx{i}) ;
W_lat{i} = W ;
W = widthdeg_lon(W_idx{i}) ;
W_widthdeg{i} = W ;
end
for i = 1:length(E_idx)
    E = lon(E_idx{i}) ;
E_lon{i} = E ;
    E = lat(E_idx{i}) ;
E_lat{i} = E ;
E = widthdeg_lon(E_idx{i}) ;
E_widthdeg{i} = E ;
end
%Extended Coast Variables (might be easier to just make a variable with all
%data contained within?
for i = 1:length(W_idx_coast)
    W = lon(W_idx_coast{i}) ;
W_lon_coast{i} = W ;
    W = lat(W_idx_coast{i}) ;
W_lat_coast{i} = W ;
    W = sal(W_idx_coast{i}) ;
W_sal_coast{i} = W ;
    W = temp(W_idx_coast{i}) ;
W_temp_coast{i} = W ;
    W = dep(W_idx_coast{i}) ;
W_dep_coast{i} = W ;
    W = yea(W_idx_coast{i});
W_yea_coast{i} = W ;
    W = day(W_idx_coast{i});
W_day_coast{i} = W ;
    W = mon(W_idx_coast{i}) ;
W_mon_coast{i} = W ;
    W =  watdep(W_idx_coast{i});
W_watdep_coast{i} = W ; 
end
clear('W')
%% Defining parallelogram verticies for casts (length and width need work) (working on now)
lengthrad = deg2rad(lengthdeg) ;
lengthdeg = ones(1,length(W_widthdeg)).*lengthdeg ;
half_length = lengthdeg./1.6 ; % mess with the divisibles to get the km right (this is within 1 km after limited testing)
% West half widths
for i = 1:length(W_widthdeg)
W = W_widthdeg{i}./3; %  mess with the divisibles to get the km right (this is within 1 km after limited testing)
W_half_width{i} = W ;
end
% East Half Widths
for i = 1:length(E_widthdeg)
E = E_widthdeg./3; %  mess with the divisibles to get the km right (this is within 1 km after limited testing)
E_half_width = E ;
end
% South Half Width

% West Compartments (-longitude)
for i = 1:numel(W_lon)
    % Calculate the coordinates of the vertices
    W_vertex1 = [W_lon{i} - W_half_width{i}; W_lat{i} + half_length(i)]; %careful, longitude is negative so watch signs
    W_vertex2 = [W_lon{i} - W_half_width{i}; W_lat{i} - half_length(i)];
    W_vertex3 = [W_lon{i} + W_half_width{i}; W_lat{i} - half_length(i)];
    W_vertex4 = [W_lon{i} + W_half_width{i}; W_lat{i} + half_length(i)];
end
%East Compartment
for i = 1:numel(E_lon)
    % Calculate the coordinates of the vertices
    W_vertex1 = [W_lon{i} - W_half_width{i}; W_lat{i} + half_length(i)]; %careful, longitude is negative so watch signs
    W_vertex2 = [W_lon{i} - W_half_width{i}; W_lat{i} - half_length(i)];
    W_vertex3 = [W_lon{i} + W_half_width{i}; W_lat{i} - half_length(i)];
    W_vertex4 = [W_lon{i} + W_half_width{i}; W_lat{i} + half_length(i)];
end
for i= 1:1:1 % needs to be changed ---------------------------------------------------------------------------------------------------------------------------------------
    % Arrange vertices in counterclockwise order
    vertices = [vertex1(:,i), vertex2(:,i), vertex3(:,i), vertex4(:,i)];
    % Store in cell array
    vertices_cell{i} = vertices;
end
% Rotate Verticies to allign with angle to coast
SW_angle = atan2d(dx,dy) ; % degrees
bearing = SW_angle - 90 ; % find angle in relation to y-axis
R = [cosd(bearing), -sind(bearing); sind(bearing), cosd(bearing)];
for i = 1:length(vertices_cell)
center_x = vertices_cell{i}(1,:) - W_lon(i);
center_y = vertices_cell{i}(2,:) - W_lat(i);
SW_rotated = R * [center_x;center_y] ;
SW_rotated_cell{i} = SW_rotated ;
SW_xrotated = SW_rotated_cell{i}(1,:) + W_lon(i) ;
SW_yrotated = SW_rotated_cell{i}(2,:) + W_lat(i) ;
SW_xrotatedcell{i} = SW_xrotated ;
SW_yrotatedcell{i} = SW_yrotated ;
SW_rotated_verticies = [SW_xrotatedcell{i};SW_yrotatedcell{i}] ;
SW_rotvertcell{i} = SW_rotated_verticies ;
end
clear('SW_rotated_verticies','SW_yrotated','SW_xrotated','SW_rotated_cell','SW_yrotated','SW_xrotated','SW_yrotatedcell','SW_xrotatedcell','SW_rotated_cell','SW_rotated')
% Clear Unneeded variables
clear('vertex1')
clear('vertex2')
clear('vertex3')
clear('vertex4')
clear('R')
% Rectangle Calculations
% Interpolate every 1 m for the SW coastal region
max_length = max(cellfun(@numel,W_dep_coast)) ;
DepInterval = (0:1:max_length)';
run = 2 ; % will rerun interpolated, temp and salinity as well as cast box indicies, when changing box size you need to set to 1
for i = 1:1:length(W_dep_coast)
    if run == 1 
        SW_int_temp{1,i} = interp1(W_dep_coast{1,i},W_temp_coast{1,i}, DepInterval, 'linear') ;
    end
end
for i = 1:1
    if run == 1
save ("SW_int_temp.mat", "SW_int_temp");
    end
end
    load('SW_int_temp.mat')
%Salinity 
for i = 1:1:length(W_dep_coast)
   if run == 1
    SW_int_sal{1,i} = interp1(W_dep_coast{1,i},W_sal_coast{1,i}, DepInterval, 'linear') ;
   end
end
for i = 1:1
    if run == 1
 save ("SW_int_sal.mat", "SW_int_sal");
    end
end
load('SW_int_sal.mat')
% Create 15 day indicies regardless of year 
SW_date = ([W_yea_coast; W_mon_coast; W_day_coast]) ;
SW_date(4,:) = datenum(0,SW_date(2,:),SW_date(3,:));% everything set in the year 0 to make for easy sorting by day and month
SW_datenum = datenum(0,SW_date(2,:),SW_date(3,:)) ;% everything set in the year 0 to make for easy sorting by day and month
day_range_2 = day_range*2 ; % For dynamic naming
middle_idx = SW_date(4,:) >= day_range & SW_date(4,:) <= 365-day_range;
earlyjan_idx = (SW_date(4,:) < day_range) ; % only includes earliest part of Jan that overlaps with interval 
latedec_idx = (SW_date(4,:) > 365-day_range); % only includes latest part of Dec that overlaps with interval 
Jan_range = 365-day_range + (SW_datenum(earlyjan_idx) ); % Earliest Day to include
Jan_range(2,:) = SW_datenum(earlyjan_idx)+day_range ; % Latest Day to include
Middle_range = SW_datenum(middle_idx)- day_range ;
Middle_range(2,:) = SW_datenum(middle_idx)+ day_range ;
Dec_range = SW_datenum(latedec_idx)- day_range ;
Dec_range(2,:) = SW_datenum(latedec_idx)-365+day_range;
SW_ranges = [Jan_range,Middle_range,Dec_range] ;
clear('middle_idx')
clear('earlyjan_idx')
clear('latedec_idx')
% Create indicies for individual rectangles 
SW_poly_idxcell = cell(length(W_lon_coast), 1); % preallocate
for i = 1:length(SW_rotvertcell)
    for j = 1:length(W_lon_coast)
        if run == 1 ;
            SW_poly_idx = inpolygon(W_lon_coast,W_lat_coast,SW_rotvertcell{i}(1,:),SW_rotvertcell{i}(2,:)) ; % extended area +core casts that fall into core cast boxes
SW_poly_idxcell{j} = SW_poly_idx ;
    end
    end
end
for i = 1:1
    if run == 1
 save("SW_poly_idx.mat",'SW_poly_idxcell')
    end
end
    load('SW_poly_idx.mat')
clear('SW_poly_idx')
% Create indicies for dates
jan_idx = false(1,length(SW_ranges)) ;
Middle_idx = false(1,length(SW_ranges)) ;
dec_idx = false(1,length(SW_ranges)) ;
% Overlapping January Dates
for i = 1:length(SW_ranges)
    for j = 1:length(Jan_range)
    if SW_datenum(i) >= Jan_range(1,j) | SW_datenum(i) <= Jan_range(2,j)
        jan_idx(i) = 1 ;
        Jan_cell{j} = jan_idx ;
    end
    end
end
% Middle Dates
for i = 1:length(SW_ranges)
    for j = 1:length(Middle_range)
    if SW_datenum(i) >= Middle_range(1,j) & SW_datenum(i) <= Middle_range(2,j)
Middle_idx(i) = 1;
Middle{j} = Middle_idx ;
    end
    end
end
% Overlapping Dec Dates
for i = 1:length(SW_ranges)
    for j = 1:length(Dec_range)
         if SW_datenum(i) >= Dec_range(1,j) | SW_datenum(i) <= Dec_range(2,j)
        dec_idx(i) = 1 ;
       Dec_cell{j} = dec_idx ;
    end
    end
end
% Concatenate Overlapping [Jan,Middle,Dec] cell arrays (Looks Correct but need to test
SW_date_idx_cell =  ([Jan_cell,Middle,Dec_cell]) ;
% Combine Date and Polygon indicies for temp and sal
for i = 1:length(SW_poly_idxcell)
    SW_poly_temp = SW_int_temp(SW_poly_idxcell{i} & SW_date_idx_cell{i}) ;
    SW_poly_temp_cell{i} = SW_poly_temp ;
end
for i = 1:length(SW_poly_idxcell)
    SW_poly_sal = SW_int_sal(SW_poly_idxcell{i} & SW_date_idx_cell{i}) ;
    SW_poly_sal_cell{i} = SW_poly_sal ;
end
for i = 1:length(SW_poly_idxcell)
    SW= W_lon_coast(SW_poly_idxcell{i} & SW_date_idx_cell{i}) ;
    SW_poly_lon_cell{i} = SW ;
end
for i = 1:length(SW_poly_idxcell)
    SW= W_lat_coast(SW_poly_idxcell{i} & SW_date_idx_cell{i}) ;
    SW_poly_lat_cell{i} = SW ;
end
% Calculate Means
counts = zeros(1,numel(SW_poly_temp_cell) );
max_length = max(counts) ;
for i = 1:length(SW_poly_sal_cell)
    insidearray = cell2mat(SW_poly_sal_cell{i}) ;
    SW_sal_expanded{1,i} = insidearray ;
end
for i = 1:length(SW_poly_temp_cell)
    insidearray = cell2mat(SW_poly_temp_cell{i}) ;
    SW_temp_expanded{1,i} = insidearray ;
end
% Means 
for i = 1:length(SW_sal_expanded)
    SW_means = mean(SW_sal_expanded{i},2,'omitnan') ;
    SW_sal_means_cell{i} = SW_means ;
end
for i = 1:length(SW_temp_expanded)
    SW_means = mean(SW_temp_expanded{i},2,'omitnan') ;
    SW_temp_means_cell{i} = SW_means ;
end
% Calculate Std Dev
for i = 1:length(W_lon)
    SW = cell2mat(SW_poly_temp_cell{i}) ;
    SW_poly_temp_clean{i} = SW ;
end
for i = 1:length(W_lon)
    SW = cell2mat(SW_poly_sal_cell{i}) ;
    SW_poly_sal_clean{i} = SW ; 
end
for i = 1:length(W_lon)
Std = std(SW_poly_sal_clean{i},0,2,'omitnan') ;
Std_sal{i} = Std ;
end
for i = 1:length(W_lon)
Std = std(SW_poly_temp_clean{i},0,2,'omitnan') ;
Std_temp{i} = Std ;
end
% Remove unused casts from int_sal and int_temp
not_empty = ~cellfun(@isempty,SW_sal_means_cell) ;
SW_int_sal = SW_int_sal(not_empty) ;
SW_int_temp = SW_int_temp(not_empty) ;
SW_sal_means_cell = SW_sal_means_cell(not_empty) ;
SW_temp_means_cell = SW_temp_means_cell(not_empty) ;
% Calculate Anomalies
for i = 1:length(SW_int_temp)
SW_Anmly = SW_sal_means_cell{i} - SW_int_sal{i} ;
SW_sal_Anmly_cell{i} = SW_Anmly ;
end
for i = 1:length(SW_int_temp)
SW_Anmly = SW_temp_means_cell{i} - SW_int_temp{i} ;
SW_temp_Anmly_cell{i} = SW_Anmly ;
end
% Number of Obersvations at each depth
for i= 1:length(DepInterval)
    for j= 1:length(SW_poly_temp_clean)
        if ~isempty(SW_poly_temp_clean{j}) && size(SW_poly_temp_clean{j}, 1) >= i ;
Depth_obs = nnz(~isnan(SW_poly_temp_clean{j}(i,:))) ;
Depth_obs_cell{j}(i,:) = Depth_obs ;
    end
    end
end
clear('SW')
clear('Std')
clear('Jan_cell')
clear('Jan_idx')
clear('Jan_range')
clear('Dec_cell')
clear('Dec_idx')
clear('Dec_range')
clear('Middle')
clear('Middle_idx')
clear('Middle_range')
clear('SW_poly_temp')
clear('SW_poly_sal')
clear('SW_Anmly')
clear('dec_idx')
clear('')
clear('Depth_obs')
clear('dx')
clear('dy')
clear('half_length')
clear('W_half_width')
clear('insidearray')
clear('lengthdeg')
clear('lengthdeg_coast')
clear('mon')
clear('day')
clear('yea')
clear('notempty')
clear('lengthrad')
clear('SW_poly_idxcell') % might need this one
clear('W_widthdeg')

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
plot(SW_rotvertcell{1280}(1,:),SW_rotvertcell{1280}(2,:), 'b' )
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
