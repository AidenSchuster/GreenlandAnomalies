% Load Data
cd('C:\Users\ajs82292\Desktop\Research\Matlab\Source\Greenland_Melt') ;
addpath('seawater','C:\Users\ajs82292\Desktop\Research\Matlab\Source\seawater') ;
load("02cleanNODC_updated.mat")
load('OMG_data.mat')
load("x_coast.mat")
load("y_coast.mat")
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
% Extend it further north using same ETOPO dataset and 0-5 contour
[Z_eto,R] = readgeoraster('extendedNcoast.tiff') ;
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

%
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
% added coastline
run = 2 ;
for i = 1:1:1
    if run == 1
clf
hold on
plot(cx,cy,'k')
plot(x_coast,y_coast,'r')
xlim([-45,-30])
ylim([59,77])
daspect([1 aspect_ratio 1])
extra_x = [] ;
extra_y = [] ;
run = true;
while run
    [x, y, button] = ginput(1);
    if isempty(button) || button == 13
        run = false;
    else
        plot(x, y, 'ro'); % Plot the point
        xlim([x-3, x+3])
        ylim([y-3, y+3])
        extra_x = [extra_x; x];
        extra_y = [extra_y; y];
    end
end
save 'extra_coast.mat' extra_x extra_y
    end
end
load 'extra_coast.mat' extra_x extra_y
reff_cord = [] ;
%combine coasts
x_coast = [x_coast;extra_x] ;
y_coast = [y_coast;extra_y] ;
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
% Slope and inverse
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
%widthdeg_coast = recwidth_coast./widthdeg_lon ; % Think I can get rid of
%this
%lengthdeg_coast = km2deg(reclength_coast) ;% For larger polygon
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
clear widthdeg_lon extra_x extra_y
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
polygon_x = [-45;x_coast;-30] ; % needs to change if we redo coastline
polygon_y = [78.1036;y_coast;78.1036] ;
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
clear x_perp_plus2 y_perp_minus y_perp_minus2 y_perp_plus y_perp_plus2 differences_plus2 differences_plus differences_minus differences_minus2 exten_minus2_ index_minus2 length_2 
clear ref_minusy ref_plusy ref_plus2y dist lengthdeg_coast off_x off_y nan_idx2 nan_idx num_points reff_cords 
% Find Casts within defined area (concatenate off coast and coast to make one polygon
%hand removing verticies that cause self-intersection
intersect_x = [-55.2664; -55.1943; -51.6715; -51.2208; -50.9012; -40.8574;-39; -39.1132; -38.8448; -39.3609;-39.1381;-40.181;-39.1862;-38.2976]; % Hand picked for 100 km, needds to be edited if coastline or size changes
intersect_y = [66.734; 66.1034; 61.8243; 61.7822; 61.3371;62.2753; 63.5719; 63.7534; 64.2143; 64.1142;64.2335;61.8788;64.23;64.7891] ; % Hand picked for 100 km, needds to be edited if coastline or size changes
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
inverse_new = -1./slope_new ;
inverse_new = inverse_new(idx) ;
target_width = recwidth/2 ;
target_length = reclength/2 ;
dist = .3 ;
num_points = 1000 ;
%10 km towards and away from refference line (max distance off by 2.5 km)
% minus
run = 2 ;
if run == 1
for i = 1:length(coastal_lon)
W = linspace(coastal_lon(1,i), coastal_lon(1,i) - dist, num_points); % minus longitude for western facing coasts
x_perp_minus{1,i} = W' ;
W =  inverse_new(i) * (x_perp_minus{i} - coastal_lon(1,i)) + coastal_lat(1,i);
y_perp_minus{1,i} = W ;
ref_minusx{1,i} = [x_perp_minus{i},repmat(coastal_lon(1,i),num_points,1)];
ref_minusy{1,i} = [y_perp_minus{i},repmat(coastal_lat(1,i),num_points,1)];
end
for i = 1:1:length(ref_minusy)
for j = 1:1:length(ref_minusy{1})
dist_minus{1,i}(j,1) = sw_dist(ref_minusy{i}(j,:),ref_minusx{i}(j,:),'km') ;
differences_minus{1,i}(j,1) = abs(dist_minus{1,i}(j,1)-target_width)  ;
[min_minus(1,i),index_minus(1,i)] = min(differences_minus{i}) ;
end
j = 1 ;
end
for i = 1:1:length(dist_minus)
index = index_minus(i) ;
   exten_minus(:,i) = [x_perp_minus{i}(index),y_perp_minus{i}(index)] ; %cords for side poins ~10km away
end
clear x_perp_minus y_perp_minus ref_minusx ref_minusy index
% plus
for i = 1:length(coastal_lon)
W = linspace(coastal_lon(1,i), coastal_lon(1,i) + dist, num_points); % minus longitude for western facing coasts
x_perp_plus{1,i} = W' ;
W =  inverse_new(i) * (x_perp_plus{i} - coastal_lon(1,i)) + coastal_lat(1,i);
y_perp_plus{1,i} = W ;
ref_plusx{1,i} = [x_perp_plus{i},repmat(coastal_lon(1,i),num_points,1)];
ref_plusy{1,i} = [y_perp_plus{i},repmat(coastal_lat(1,i),num_points,1)];
end
for i = 1:1:length(ref_plusy)
for j = 1:1:length(ref_plusy{1})
dist_plus{1,i}(j,1) = sw_dist(ref_plusy{i}(j,:),ref_plusx{i}(j,:),'km') ;
differences_plus{1,i}(j,1) = abs(dist_plus{1,i}(j,1)-target_width)  ;
[min_plus(1,i),index_plus(1,i)] = min(differences_plus{i}) ;
end
j = 1 ;
end
for i = 1:1:length(dist_plus)
index = index_plus(i) ;
   exten_plus(:,i) = [x_perp_plus{i}(index),y_perp_plus{i}(index)] ; %cords for side poins ~10km away
end
clear x_perp_plus y_perp_plus ref_plusx ref_plusy index
% make verticies points target reclength away
%vert1 (uses exten_minus)
slope_coast = slope_new(idx)';
for i = 1:length(exten_minus)
W = linspace(exten_minus(1,i), exten_minus(1,i) - dist, num_points); % minus longitude for western facing coasts
bodge_xminus{1,i} = W' ;
W =  slope_coast(i) * (bodge_xminus{i} - exten_minus(1,i)) + exten_minus(2,i);
bodge_yminus{1,i} = W ;
ref_minusx{1,i} = [bodge_xminus{i},repmat(exten_minus(1,i),num_points,1)];
ref_minusy{1,i} = [bodge_yminus{i},repmat(exten_minus(2,i),num_points,1)];
end
for i = 1:1:length(ref_minusy)
for j = 1:1:length(ref_minusy{1})
dist_minus{1,i}(j,1) = sw_dist(ref_minusy{i}(j,:),ref_minusx{i}(j,:),'km') ;
differences_minus{1,i}(j,1) = abs(dist_minus{1,i}(j,1)-target_length)  ;
[min_minus(1,i),index_minus_vert(1,i)] = min(differences_minus{i}) ;
end
j = 1 ;
end
for i = 1:1:length(dist_minus)
index = index_minus_vert(i) ;
   vert_1(:,i) = [bodge_xminus{i}(index),bodge_yminus{i}(index)] ; %cords for side poins ~10km away
end
clear index ref_minusx ref_minusy bodge_yminus bodge_xminus
% Vert2 (uses exten_minus)
for i = 1:length(exten_minus)
W = linspace(exten_minus(1,i), exten_minus(1,i) + dist, num_points); % minus longitude for western facing coasts
bodge_xplus{1,i} = W' ;
W =  slope_coast(i) * (bodge_xplus{i} - exten_minus(1,i)) + exten_minus(2,i);
bodge_yplus{1,i} = W ;
ref_plusx{1,i} = [bodge_xplus{i},repmat(exten_minus(1,i),num_points,1)];
ref_plusy{1,i} = [bodge_yplus{i},repmat(exten_minus(2,i),num_points,1)];
end
for i = 1:1:length(ref_plusy)
for j = 1:1:length(ref_plusy{1})
dist_plus{1,i}(j,1) = sw_dist(ref_plusy{i}(j,:),ref_plusx{i}(j,:),'km') ;
differences_plus{1,i}(j,1) = abs(dist_plus{1,i}(j,1)-target_length)  ;
[min_plus(1,i),index_plus_vert(1,i)] = min(differences_plus{i}) ;
end
j = 1 ;
end
for i = 1:1:length(dist_plus)
index = index_plus_vert(i) ;
   vert_2(:,i) = [bodge_xplus{i}(index),bodge_yplus{i}(index)] ; %cords for side poins ~10km away
end
clear index ref_plusx ref_plusy bodge_yplus bodge_xplus
%Vert 3 (uses exten_plus) WORKING ON NOW
for i = 1:length(exten_minus)
W = linspace(exten_plus(1,i), exten_plus(1,i) - dist, num_points); % minus longitude for western facing coasts
bodge_xminus{1,i} = W' ;
W =  slope_coast(i) * (bodge_xminus{i} - exten_plus(1,i)) + exten_plus(2,i);
bodge_yminus{1,i} = W ;
ref_minusx{1,i} = [bodge_xminus{i},repmat(exten_plus(1,i),num_points,1)];
ref_minusy{1,i} = [bodge_yminus{i},repmat(exten_plus(2,i),num_points,1)];
end
for i = 1:1:length(ref_minusy)
for j = 1:1:length(ref_minusy{1})
dist_minus{1,i}(j,1) = sw_dist(ref_minusy{i}(j,:),ref_minusx{i}(j,:),'km') ;
differences_minus{1,i}(j,1) = abs(dist_minus{1,i}(j,1)-target_length)  ;
[min_minus(1,i),index_minus_vert(1,i)] = min(differences_minus{i}) ;
end
j = 1 ;
end
for i = 1:1:length(dist_minus)
index = index_minus_vert(i) ;
   vert_3(:,i) = [bodge_xminus{i}(index),bodge_yminus{i}(index)] ; %cords for side poins ~10km away
end
clear index ref_minusx ref_minusy bodge_yminus bodge_xminus
% Vert 4 (uses exten_plus)
for i = 1:length(exten_plus)
W = linspace(exten_plus(1,i), exten_plus(1,i) + dist, num_points); % minus longitude for western facing coasts
bodge_xplus{1,i} = W' ;
W =  slope_coast(i) * (bodge_xplus{i} - exten_plus(1,i)) + exten_plus(2,i);
bodge_yplus{1,i} = W ;
ref_plusx{1,i} = [bodge_xplus{i},repmat(exten_plus(1,i),num_points,1)];
ref_plusy{1,i} = [bodge_yplus{i},repmat(exten_plus(2,i),num_points,1)];
end
for i = 1:1:length(ref_plusy)
for j = 1:1:length(ref_plusy{1})
dist_plus{1,i}(j,1) = sw_dist(ref_plusy{i}(j,:),ref_plusx{i}(j,:),'km') ;
differences_plus{1,i}(j,1) = abs(dist_plus{1,i}(j,1)-target_length)  ;
[min_plus(1,i),index_plus_vert(1,i)] = min(differences_plus{i}) ;
end
j = 1 ;
end
for i = 1:1:length(dist_plus)
index = index_plus_vert(i) ;
   vert_4(:,i) = [bodge_xplus{i}(index),bodge_yplus{i}(index)] ; %cords for side poins ~10km away
end
clear exten_minus exten2 exten_plus i 
%
%combine verticies (not ordered)
for i= 1:1:length(vert_1)
vertices = [vert_1(:,i), vert_2(:,i), vert_3(:,i), vert_4(:,i)];
vert_combined{i} = vertices ;
end
clear verticies
% order the verticies based on angle relative to y-axis from central point
for  i = 1:length(vert_combined)
    for j = 1:length(vert_combined{i})
angles(j) = atan2(vert_combined{i}(2,j)-coastal_lat(i),vert_combined{i}(1,j)-coastal_lon(i)) ;
angle_vert{i} = angles ;
    end
    j = 1 ;
end
for  i = 1:length(vert_combined)
    for j = 1:length(vert_combined{i})
    [~, sortOrder] = sort(angle_vert{i});
order{i} = sortOrder ;
    end
end
for i = 1:length(order)
    v = vert_combined{i} ;
    o = order{i} ;
    s = v(:,o) ;
    vert{i} = s ;
end
save("vert.mat","vert")  
end
load("vert.mat","vert")
clear v o s vert_combined vert_1 vert_2 vert_3 vert_4 vertices angles angles angle_vert bodge_yplus bodge_yminus bodge_xminus bodge_xplus order sortOrder inverse_new ref_plusx ref_minusy ref_plusy ref_minusx
clear differences_plus differences_minus dist_plus dist_minus % need to tweak some stuff here eventually
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
save ("interp_temp.mat", "interp_temp");
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
 save ("interp_sal.mat", "interp_sal");
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
    [~,CountourM] = contour(XX,YY,ZZ, iso_interval) ;
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
%% group fjord data by individual fjords (move until after cleaning
% eventually)
for i = 1:length(fjord_vert)
in_fj{i} = inpolygon(lon,lat,fjord_vert{i}(:,1),fjord_vert{i}(:,2)) ;
end
clear center_points lat1 lon1 lat2 lon2 lat3 lon3 lat4 lon4 width_nm width_deg length_nm length_deg clear center_points special_case iso_interval isobath_interval
%% clean fjord data loop through collumns and calculate derivatives (delta y /delta x)
run = 2 ;
fjord_sal_mat_fj = cellfun(@(c) c(:), interp_sal, 'UniformOutput', false) ; % fjord data 
fjord_sal_mat_fj = [fjord_sal_mat_fj{:}] ; % all data, but cleaned specifically for fjord standards
if run == 1
    diff_result = NaN(size(interp_sal_mat)) ;
    for i  = 1:length(interp_sal_mat)
    non_nan = find(~isnan(interp_sal_mat(:,i))) ;
    non_nan_store{i} = non_nan ;
    temporary = interp_sal_mat(:,i) ;
    diff_sal = diff(temporary(non_nan),1,1) ;
    diff_result(non_nan(1:end-1),i) = abs(diff_sal ./ diff(non_nan,1,1)) ;
    end
   threshold_15 = 2 ; % 0-50 m
   threshold = 0.3 ; % 50-end m
     top_15 = diff_result(1:50,:) ;
     bottom = diff_result(51:end,:) ;
     top_15_remove = top_15 >= threshold_15 ; % idx of bottom that has derivatives greater than threshold
     bottom_remove = bottom >= threshold ; % idx of bottom that has derivatives greater than threshold
     remove = [top_15_remove;bottom_remove] ; % need to account for both points that create a true remove value
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
clear top_15 bottom threshold bottom_remove top_15 threshold_15 diff_sal temporary non_nan value_below value_below_idx diff_result non_nan_store remove top_15_remove

end




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
   threshold_15 = 2 ; % 0-50 m
   threshold = 0.3 ; % 50-end m
     top_15 = diff_result(1:50,:) ;
     bottom = diff_result(51:end,:) ;
     top_15_remove = top_15 >= threshold_15 ; % idx of bottom that has derivatives greater than threshold
     bottom_remove = bottom >= threshold ; % idx of bottom that has derivatives greater than threshold
     remove = [top_15_remove;bottom_remove] ; % need to account for both points that create a true remove value
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
clear top_15 bottom threshold bottom_remove top_15 threshold_15 diff_sal temporary non_nan value_below value_below_idx diff_result non_nan_store remove top_15_remove
%same thing but for coast (for vert idx)
    diff_result = NaN(size(coast_sal_temporary)) ;
    for i  = 1:length(coast_sal_temporary)
    non_nan = find(~isnan(coast_sal_temporary(:,i))) ;
    non_nan_store{i} = non_nan ;
    temporary = coast_sal_temporary(:,i) ;
    diff_sal = diff(coast_sal_temporary(non_nan),1,1) ;
    diff_result(non_nan(1:end-1),i) = abs(diff_sal ./ diff(non_nan,1,1)) ;
    end
   threshold_15 = 2 ; % 0-50 m
   threshold = 0.3 ; % 50-end m
     top_15 = diff_result(1:50,:) ;
     bottom = diff_result(51:end,:) ;
     top_15_remove = top_15 >= threshold_15 ; % idx of bottom that has derivatives greater than threshold
     bottom_remove = bottom >= threshold ; % idx of bottom that has derivatives greater than threshold
     remove_coast = [top_15_remove;bottom_remove] ; % need to account for both points that create a true remove value
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
clear remove_coast top_15 bottom threshold bottom_remove top_15 threshold_15 diff_sal temporary non_nan value_below value_below_idx diff_result non_nan_store top_15_remove
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
for i = 1:length(vert)
        if run == 1  
            poly_idx = inpolygon(lon_a,lat_a,vert{i}(1,:),vert{i}(2,:)) ; % each collumn corresponds to a different cast
                if datenum_coast(i) > 14 && datenum_coast(i) < 350 % Middle ranges
                   coastal_idx = datenum_a >= range(1,i) & datenum_a <= range(2,i) ; 
                elseif datenum_coast(i) <= 14 || datenum_coast(i) >= 351 % Early Jan and Late Dec 
                   coastal_idx = datenum_a >= range(1,i) | datenum_a <= range(2,i) ;
                end
                % combine indices
                combined_idx = poly_idx & coastal_idx ;
                coast_find{i} = find(combined_idx == 1) ;
                save("coast_find.mat",'coast_find')
        end
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
num_repeat = ceil(129752/length(pattern)) ; % # should be equal to length of open_dist -2 have to change if you change length!
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
for i = 1:1
    if run == 1
    open_temp_mat = [] ;
    open_sal_mat = [] ;
        for j = 1:length(open_find)
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
end
load("open_temp_anom.mat") 
load("open_sal_anom.mat")

clear col unique_col
%load("open_temp_std.mat")
%load("open_sal_std.mat")
clear coast_find open_find open_temp_mat open_sal_mat % will need these back eventually
%% Fjord Anomaly Calculations (include profiles + or - 15 days)
run = 1 ;
if run == 1
for i = 1:length(fjord_box_cords)
fj_box_profiles{i} = inpolygon(lon,lat,fjord_box_cords{i}(:,1),fjord_box_cords{i}(:,2)) ; % gives index of all profiles within the fjord anomaly boxes.
fj_profiles{i} = inpolygon(lon,lat,fjord_vert{i}(:,1),fjord_vert{i}(:,2)) ; % idx of profiles inside defined fjords
end
fj_combined = [] ;
fj_box_combined = [] ;
for i = 1:length(fj_profiles)
fj_box_idx = find(fj_box_profiles{i} == 1) ; % create index instead of logical index
fj_box_combined = [fj_box_combined,fj_box_idx] ;
fj_box_idx_back{i} = fj_box_idx ;
fj_profile_back{i} = find(fj_profiles{i} == 1) ;
end
clear fj_box_idx 
% get box dates
box_datenum = datenum(fj_box_combined) ; 
for i = 1:length(fj_profile_back)
fj_datenum{i} = datenum(fj_profile_back{i}) ;
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
        fj_anom_date_idx = fj_box_combined(date_idx{i}{j}) ;
        fj_anom_idx{i}{j} = intersect(fj_anom_date_idx,fj_box_idx_back{i}) ;
    end
end
% compare date casts and box casts, only keep if a value is present in both (fj_profile_back contains idx for fj profiles, fj_anom_idx for box profiles that fit the correct box and date)_
clear date_idx fj_anom_date_idx fj_datenum box_datenum fj_box_idx_back fj_box_combined fj_combined fj_profiles fj_box_profiles
% calculate means, and subtract fjord values for anomaly
%interp_sal_mat_fj = cellfun(@(c) c(:), interp_sal, 'UniformOutput', false) ;
%interp_sal_mat_fj = [interp_sal_mat_fj{:}] ;
for i = 1:length(fj_anom_idx) 
    for j = 1:length(fj_anom_idx{i})
        if length(fj_anom_idx{i}{j}) >= 2 % reduced limit to calculate anomalies
            fj_anoms{i}(:,j) = interp_sal_mat_fj(:,fj_profile_back{i}(j)) - mean(interp_sal_mat_fj(:,fj_anom_idx{i}{j}),2,'omitnan') ;
        else
            fj_anoms{i}(:,j) = NaN(size(interp_sal_mat_fj, 1), 1); 
        end
    end
end
save fj_anoms.mat fj_anoms
save fj_coords.mat fj_coords
end
load fj_anoms.mat 
load fj_coords.mat
clear interp_sal interp_temp fj_profile_back %fj_anom_idx
%%
% Pressure (db) from Depth and Density (can clear after getting potential temp)
%numpoints = length(DepInterval) ;
%for i = 1:length(lat_a)
%k = sw_pres(DepInterval,repmat(lat_a(i),1,numpoints)) ;
%press_a(:,i) = k ;
%end
%coastal_sal = interp_sal_mat(:,in_a) ; (don't think i need this any
%longer)
%coastal_temp = interp_temp_mat(:,in_a) ;
%coast_press = press_a(:,in_a) ;
%open_press = press_a(:,~in_a) ;
%open_dens = sw_dens(open_sal,open_temp,open_press) ; % kg/m^3
%coastal_dens = sw_dens(coastal_sal,coastal_temp,coast_press) ;
clear press_a
% Vectorize and combine
clear  coastal_temp open_temp interp_temp_a interp_sal_a %interp_sal_mat interp_temp_mat
%%
%top 300 m of both open and coastal
sal_combined = [open_sal,coastal_sal] ;
coastal_sal = coastal_sal(1:300,:) ;
open_sal = open_sal(1:300,:) ;
length_open = length(open_sal) ;
sal_combined = sal_combined(1:300,:) ;
May_a = mon_a == 5 ;
Jun_a = mon_a == 6 ;
Jul_a = mon_a == 7 ;
Aug_a = mon_a == 8 ;
Sep_a = mon_a == 9 ;
Oct_a = mon_a == 10 ;
sal_anom_combined = [open_sal_anom,coast_sal_anom] ;
sal_anom_combined = sal_anom_combined(1:300, :);
lat_combined = [lat_open,coastal_lat] ; % this coastal_lat/lon is NOT correct
lon_combined = [lon_open,coastal_lon ] ;
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
% same for open salinity 
for i = 1:size(open_sal, 2)
    top10 = open_sal(1:10, i);  % Extract the top 10 rows of the column
    nonNaN_values = top10(~isnan(top10));  % Extract non-NaN values
    %if length(nonNaN_values) >= 3
        % Average the top three non-NaN values
        %avg_value = mean(nonNaN_values(1:3));
        %open_sal(1:10, i) = fillmissing(top10, 'constant', avg_value);
    if length(nonNaN_values) >= 1
        % If fewer than 3 non-NaN values, use the top-most non-NaN value
        top_most_value = nonNaN_values(1);
        open_sal(1:10, i) = fillmissing(top10, 'constant', top_most_value);
    end
end
% same for coastal_salinity
for i = 1:size(coastal_sal, 2)
    top10 = coastal_sal(1:10, i);  % Extract the top 10 rows of the column
    nonNaN_values = top10(~isnan(top10));  % Extract non-NaN values
    %if length(nonNaN_values) >= 3
        % Average the top three non-NaN values
        %avg_value = mean(nonNaN_values(1:3));
        %coastal_sal(1:10, i) = fillmissing(top10, 'constant', avg_value);
    if length(nonNaN_values) >= 1
        % If fewer than 3 non-NaN values, use the top-most non-NaN value
        top_most_value = nonNaN_values(1);
        coastal_sal(1:10, i) = fillmissing(top10, 'constant', top_most_value);
    end
end
clear top_most_value top10 nonNaN_values avg_value Sep_a Aug_a Oct_a Jul_a Jun_a May_a open_sal_anom coast_sal_anom coast_temp_anom open_temp_anom
%% Canadian area selection
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
month_selected = 7 ; % for sprintf
year_selected = 2012 ; % for sprintf, replace with desired year
if size(sal_anom_combined,2) > 301
sal_anom_combined = sal_anom_combined' ;
can_invert = canada' ;
end
if size(coastal_sal,2) > 301
coastal_sal = coastal_sal' ;
open_sal = open_sal' ;
end
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
can_mon = open_mon(canada) ;
open_mon = open_mon(~canada) ;
coastal_yea = yea_a(in_a) ;
coastal_yea = coastal_yea(coast_idx) ;
open_yea = yea_a(~in_a) ;
open_yea = open_yea(open_idx) ;
can_yea = open_yea(canada) ;
open_yea = open_yea(~canada) ;
% month
month_string = {'January', 'February', 'March', 'April', 'May', 'June','July', 'August', 'September', 'October', 'November', 'December'} ;
mon_comb = [open_mon,coastal_mon] ;
month_s = mon_comb == month_selected ; % # of desired month;
month_open = open_mon == month_selected ;
month_coast = coastal_mon == month_selected ;
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
% Salinity Year
year_coast_n_test = yea_a(in_a) ;
year_open_n_test = yea_a(~in_a) ;
year_coast_n = year_coast_n_test == year_selected ;
year_open_n = year_open_n_test == year_selected ;
%combine
year_mon = year & month_s ;
year_mon_open = year_open & month_open ;
year_mon_coast = year_coast & month_coast ;
% Salinity combine
yeamon_n_coast = year_coast_n & month_coast_n ;
yeamon_n_open = year_open_n & month_open_n ;
%clear open_mon coastal_mon %combined_idx_open combined_idx_coast
%copy = sal_anom_combined ; % for copying
in_combined = inpolygon(lon_combined,lat_combined,combined_x,combined_y) ;
open_sal_anom = sal_anom_combined(~in_combined,:) ;
coast_sal_anom = sal_anom_combined(in_combined,:) ;
lon_coast = lon_combined(in_combined) ;
lon_open = lon_combined(~in_combined) ;
lat_coast = lat_combined(in_combined) ;
lat_open = lat_combined(~in_combined) ;
% Canada
canada_n = inpolygon(lon_open_n,lat_open_n,x_canada,y_canada) ;
can_lon_n = lon_open_n(canada_n) ;
can_lat_n = lat_open_n(canada_n) ;
lat_open_n = lat_open_n(~canada_n) ;
lon_open_n = lon_open_n(~canada_n) ;
can_lon = lon_open(canada) ;
can_lat = lat_open(canada) ;
lon_open = lon_open(~canada) ;
lat_open = lat_open(~canada) ;
yeamon_n_open = yeamon_n_open(~canada_n) ;
% sal/anom/temp/anom
run = 1 ;
if run == 1
    run = 2 ;
    if size(canada,1) == 1
    can_invert = canada' ;
    can_n_invert = canada_n' ;
    end
    can_sal = open_sal(can_n_invert,:) ; % copy so we can rerun this section as needed
    can_sal_anom = open_sal_anom(can_invert,:) ;
    open_sal_anom = open_sal_anom(~can_invert,:) ;
    open_sal = open_sal(~can_n_invert,:) ;
end
clear length_open yea year_coast month_coast month_open year_open year_coast_n year_open_n month_coast_n month_open_n month_coast_n_test year_open_n_test year_coast_n_test
clear month_open_n_test canada canada_n can_invert can_n_invert
%% Subsections and removals (maybe get rid of this?)
% Index casts that have NO nan's through 300 m
%sal_anom_combined = copy ;
%sal_anom_combined = sal_anom_combined(:,1:10) ;
%NaN_idx = isnan(sal_anom_combined) ;
%columnSums = sum(NaN_idx,2); % Sum along rows to get column sums
%NoNaN = columnSums == 0; % Find columns with no true values
%sal_anom_combined = sal_anom_combined(NoNaN,:) ;
%clear columnSums NaN_idx NoNaN
%% PCA change month/index as desired
sal_anom = coast_sal_anom(year_mon_coast, :); % should just be able to change this
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
%% corresponding raw salinity data (includes all salinity profiles)
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
[coeff_n, score_n, latent_n , ~] = eof225(sal_minus,NaN,50); % Renato's Function (50 is number he gave) (very slow so reduce NaN's as much as possible)
first_PC_n = score_n(:,1) ; % first principal component
first_coeff_n = coeff_n(:,1); % first pc coeff
second_PC_n = score_n(:,2) ; % second
third_PC_n = score_n(:,3) ;
explained_n = 100 * latent_n / sum(latent_n);
clear last_nan_col
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
