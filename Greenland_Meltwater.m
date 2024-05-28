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
target  = 100 ; % km from coast
%% Defining Coastline
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
run = 2 ;
for i = 1:1:1
    if run == 1
clf
hold on
plot(cx,cy,'k')
xlim([-75,-37])
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
reff_cord = [] ;
for i = 1:1:length(x_reff)-1
    x = [(x_reff(i)+x_reff(i+1))/2,(y_reff(i)+y_reff(i+1))/2] ;   %middle point of refference lines
    reff_cord(i,:) = x ;
end
for i = 1:1:length(x_reff)-1
    x = (y_reff(i+1,1)-y_reff(i,1))/(x_reff(i+1,1)-x_reff(i,1)) ;
    reff_slope(i,:) = x ;
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
width_2 = recwidth ; % for dynamic plot titles
length_2 = reclength ; % for dynamic plot titles
recwidth_coast = 10*2 ; % do not change, for defining larger polygon
reclength_coast = 20*2 ; % do not change, for defining larger polygon
lengthdeg = km2deg(reclength) ;
widthdeg_lon = 111.120 .*cosd(lat) ;
widthdeg_lon = recwidth./widthdeg_lon ;
widthdeg_coast = recwidth_coast./widthdeg_lon ; % For larger polygon
lengthdeg_coast = km2deg(reclength_coast) ;% For larger polygon
%%
% create points along lines
dist = 3 ; 
num_points = 5000 ; 
% All Coastline + and - for 1/2
for i = 1:length(x_coast)-1
W = linspace(x_coast(i,1), x_coast(i,1) - dist, num_points); % minus longitude for western facing coasts
x_perp_minus{i,1} = W' ;
W =  inverse(i) * (x_perp_minus{i} - x_coast(i,1)) + y_coast(i,1);
y_perp_minus{i,1} = W ;
end
for i = 1:length(x_coast)-1
W = linspace(x_coast(i,1), x_coast(i,1) + dist, num_points); % minus longitude for western facing coasts
x_perp_plus{i,1} = W' ;
W =  inverse(i) * (x_perp_plus{i} - x_coast(i,1)) + y_coast(i,1);
y_perp_plus{i,1} = W ;
end
for i = 1:length(x_coast)-1
W = linspace(x_coast(i+1,1), x_coast(i+1,1) + dist, num_points); % minus longitude for western facing coasts
x_perp_plus2{i,1} = W' ; 
W =  inverse(i) * (x_perp_plus2{i} - x_coast(i+1,1)) + y_coast(i+1,1);
y_perp_plus2{i,1} = W ;
end
for i = 1:length(x_coast)-1
W = linspace(x_coast(i+1,1), x_coast(i+1,1) - dist, num_points); % minus longitude for western facing coasts
x_perp_minus2{i,1} = W' ; 
W =  inverse(i) * (x_perp_minus2{i} - x_coast(i+1,1)) + y_coast(i+1,1);
y_perp_minus2{i,1} = W ;
end
clear W
%% 
% setting up sw_dist
for i = 1:length(x_perp_plus)
%plus
ref_plusx{i,1} = [x_perp_plus{i},repmat(x_coast(i,1),num_points,1)];
ref_plusy{i,1} = [y_perp_plus{i},repmat(y_coast(i,1),num_points,1)];
ref_plus2x{i,1} = [x_perp_plus2{i},repmat(x_coast(i+1,1),num_points,1)];
ref_plus2y{i,1}= [y_perp_plus2{i},repmat(y_coast(i+1,1),num_points,1)];
%minus
ref_minusx{i,1} = [x_perp_minus{i},repmat(x_coast(i,1),num_points,1)];
ref_minusy{i,1} = [y_perp_minus{i},repmat(y_coast(i,1),num_points,1)];
ref_minus2x{i,1} = [x_perp_minus2{i},repmat(x_coast(i+1,1),num_points,1)];
ref_minus2y{i,1}= [y_perp_minus2{i},repmat(y_coast(i+1,1),num_points,1)];
end
% sw_dist (only do for index plus minus when appropriate(not currently done))
for i = 1:1:length(ref_plusy)
    for j = 1:1:length(ref_plusy{1})
    dist_plus{i,1}(j,1) = sw_dist(ref_plusy{i}(j,:),ref_plusx{i}(j,:),'km') ;
    dist_plus2{i,1}(j,1) = sw_dist(ref_plus2y{i}(j,:),ref_plus2x{i}(j,:),'km') ;
    dist_minus{i,1}(j,1) = sw_dist(ref_minusy{i}(j,:),ref_minusx{i}(j,:),'km') ;
    dist_minus2{i,1}(j,1) = sw_dist(ref_minus2y{i}(j,:),ref_minus2x{i}(j,:),'km') ;
    end
    j = 1 ;
end
%% 
% Find point closest to x km within y km
for i = 1:1:length(dist_minus)
    for j =1:1:length(ref_minusy{1})
differences_minus{i,1}(j,1) = abs(dist_minus{i,1}(j,1)-target)  ;
differences_minus2{i,1}(j,1) = abs(dist_minus2{i,1}(j,1)-target) ;
    end
    j= 1 ;
[min_minus(i,1),index_minus(i,1)] = min(differences_minus{i}) ;
[min_minus2(i,1),index_minus2(i,1)] = min(differences_minus2{i}) ;
end
for i = 1:1:length(dist_plus)
    for j =1:1:length(ref_plusy{1})
differences_plus{i,1}(j,1) = abs(dist_plus{i,1}(j,1)-target)  ;
differences_plus2{i,1}(j,1) = abs(dist_plus2{i,1}(j,1)-target) ;
    end
    j= 1 ;
[min_plus(i,1),index_plus(i,1)] = min(differences_plus{i}) ;
[min_plus2(i,1),index_plus2(i,1)] = min(differences_plus2{i}) ;
end
% Index back in to be left with cordinates for plus and minus 1/2
for i = 1:1:length(dist_minus)
index = index_minus(i) ;
   exten_minus(i,:) = [x_perp_minus{i}(index),y_perp_minus{i}(index)] ;
index = index_minus2(i) ;
   exten_minus2(i,:) = [x_perp_minus2{i}(index),y_perp_minus2{i}(index)] ;
end
for i = 1:1:length(dist_plus)
index = index_plus(i) ;
   exten_plus(i,:) = [x_perp_plus{i}(index),y_perp_plus{i}(index)] ;
index = index_plus2(i) ;
   exten_plus2(i,:) = [x_perp_plus2{i}(index),y_perp_plus2{i}(index)] ;
end
% Test inpolygon to remove stuff inside coastline (need to keep 1/2
% seperate for avging later
polygon_x = [-67;x_coast;-40] ; % needs to change if we redo coastline
polygon_y = [95;y_coast;95] ;
in_plus = inpolygon(exten_plus(:,1),exten_plus(:,2),polygon_x,polygon_y) ;
in_plus2 = inpolygon(exten_plus2(:,1),exten_plus2(:,2),polygon_x,polygon_y) ;
in_minus = inpolygon(exten_minus(:,1),exten_minus(:,2),polygon_x,polygon_y) ;
in_minus2 = inpolygon(exten_minus2(:,1),exten_minus2(:,2),polygon_x,polygon_y) ;
% removing the inside coast points so I can stich together and average
% valid points
exten_plus_x = nan(size(exten_plus(:,1))) ;
exten_plus_x(~in_plus) = exten_plus(~in_plus,1) ;
exten_plus_y = nan(size(exten_plus(:,1))) ;
exten_plus_y(~in_plus) = exten_plus(~in_plus,2) ;
exten_plus_ = ([exten_plus_x,exten_plus_y]) ;
exten_plus2_x = nan(size(exten_plus2(:,1))) ;
exten_plus2_x(~in_plus2) = exten_plus2(~in_plus2,1) ;
exten_plus2_y = nan(size(exten_plus2(:,1))) ;
exten_plus2_y(~in_plus2) = exten_plus2(~in_plus2,2) ;
exten_plus2_ = ([exten_plus2_x,exten_plus2_y]) ;
exten_minus_x = nan(size(exten_minus(:,1))) ;
exten_minus_x(~in_minus) = exten_minus(~in_minus,1) ;
exten_minus_y = nan(size(exten_minus(:,1))) ;
exten_minus_y(~in_minus) = exten_minus(~in_minus,2) ;
exten_minus_ = ([exten_minus_x,exten_minus_y]) ;
exten_minus2_x = nan(size(exten_minus2(:,1))) ;
exten_minus2_x(~in_minus2) = exten_minus(~in_minus2,1) ;
exten_minus2_y = nan(size(exten_minus2(:,1))) ;
exten_minus2_y(~in_minus2) = exten_minus2(~in_minus2,2) ;
exten_minus2_ = ([exten_minus2_x,exten_minus2_y]) ;
% Stich plus and minus together 
nan_idx = isnan(exten_plus_) ;
exten_plus_(nan_idx) = exten_minus_(nan_idx) ;
exten = exten_plus_ ;
nan_idx2 = isnan(exten_plus2_) ;
exten_plus2_(nan_idx2) = exten_minus2_(nan_idx2) ;
exten2 = exten_plus2_  ;
% concatenate lon,lat from both 1/2 for ordering STILL NEED TO ORDER
offset = [NaN,NaN] ;
exten = [exten;offset] ;
exten2 = [offset;exten2] ;
off_x = [exten(:,1),exten2(:,1)] ;
off_y = [exten(:,2),exten2(:,2)] ;
% average offset off coast points to result in one point per point
off_coast = [mean(off_x,2,'omitnan'),mean(off_y,2,'omitnan')] ;
clear exten_plus_x exten_plus_y exten_plus2_x exten_plus2_y  exten_minus_x exten_minus_y exten_minus2_x exten_minus2_y 
clear in_minus2 in_minus in_plus in_plus2 
clear exten_minus exten_minus2 exten_plus exten_plus2 exten_plus_ exten_plus2_ exten_minus_ dist_plus dist_minus dist_minus2 dist_plus2 exten_minus2 index i index_minus index_plus index_plus2
clear  j inverse intercept min_plus2 min_plus min_minus2 min_minus ref_minus2x ref_minusx ref_minus2y ref_plus2x ref_plusx run slope x_perp_minus x_perp_minus2 x_perp_plus
clear x_perp_plus2 y_reff y_perp_minus y_perp_minus2 y_perp_plus y_perp_plus2 x_reff differences_plus2 differences_plus differences_minus differences_minus2 exten_minus2_ index_minus2 length_2 
clear ref_minusy ref_plusy ref_plus2y dist lengthdeg_coast off_x off_y
%% Find Casts within defined area (concatenate off coast and coast to make one polygon
%hand removing verticies that cause self-intersection
intersect_x = [-55.2664; -55.1943; -51.6715; -51.2208; -50.9012; -40.8574;-39; -39.1132; -38.8448; -39.3609;-39.1381]; % Hand picked for 100 km, needds to be edited if coastline or size changes
intersect_y = [66.734; 66.1034; 61.8243; 61.7822; 61.3371;62.2753; 63.5719; 63.7534; 64.2143; 64.1142;64.2335] ; % Hand picked for 100 km, needds to be edited if coastline or size changes
tolerance = 1e-4 ;
idx_remove = [] ;
for i = 1:length(intersect_y)
   idx = find(abs(off_coast(:,1) - intersect_x(i)) < tolerance & abs(off_coast(:,2) - intersect_y(i)) < tolerance);
    % Append the found indices to idx_remove
    idx_remove = [idx_remove;idx] ;
end
off_coast(idx_remove,:) = [] ;
%combine into polygon
poly_off = [-66.6413,75.7] ;
combined_x = [off_coast(:,1);flip(x_coast)] ;
combined_y = [off_coast(:,2);flip(y_coast)] ;
in = inpolygon(lon,lat,combined_x,combined_y) ;
%% Determine which refference point to base slope of rectangles on WORKING ON NOW
reff_x = reff_cord(:,1) ;
reff_y = reff_cord(:,2) ;
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
%%
clear temp_x temp_y temp_distance min_distance
%for i = 1:length(W_idx)
%    W = lon(W_idx{i}) ;
%W_lon{i} = W ;
%    W = lat(W_idx{i}) ;
%W_lat{i} = W ;
%W = widthdeg_lon(W_idx{i}) ;
%W_widthdeg{i} = W ;
%end
%for i = 1:length(E_idx)
%    E = lon(E_idx{i}) ;
%E_lon{i} = E ;
%    E = lat(E_idx{i}) ;
%E_lat{i} = E ;
%E = widthdeg_lon(E_idx{i}) ;
%E_widthdeg{i} = E ;
%end
%Extended Coast Variables (might be easier to just make a variable with all
%data contained within?
%for i = 1:length(W_idx_coast)
%    W = lon(W_idx_coast{i}) ;
%W_lon_coast{i} = W ;
%    W = lat(W_idx_coast{i}) ;
%W_lat_coast{i} = W ;
%    W = sal(W_idx_coast{i}) ;
%W_sal_coast{i} = W ;
%    W = temp(W_idx_coast{i}) ;
%W_temp_coast{i} = W ;
%    W = dep(W_idx_coast{i}) ;
%W_dep_coast{i} = W ;
%    W = yea(W_idx_coast{i});
%W_yea_coast{i} = W ;
%    W = day(W_idx_coast{i});
%W_day_coast{i} = W ;
%    W = mon(W_idx_coast{i}) ;
%W_mon_coast{i} = W ;
%    W =  watdep(W_idx_coast{i});
%W_watdep_coast{i} = W ; 
%end
%clear('W')
%% Defining parallelogram verticies for casts (length and width need work) (working on now)
lengthrad = deg2rad(lengthdeg) ;
lengthdeg = ones(1,length(W_widthdeg)).*lengthdeg ;
half_length = lengthdeg./1.6 ; % mess with the divisibles to get the km right (this is within 1 km after very limited testing)
% West half widths
for i = 1:length(W_widthdeg)
W = W_widthdeg{i}./3; %  mess with the divisibles to get the km right (this is within 1 km after very limited testing)
W_half_width{i} = W ;
end
% East Half Widths
for i = 1:length(E_widthdeg)
E = E_widthdeg./3; %  mess with the divisibles to get the km right (this is within 1 km after very limited testing)
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
