%recreating FINESST OMG NODC plot
% Set figure size
figure('Units', 'inches', 'Position', [1, 1, 8, 6]); % Adjust size as needed

% Existing plotting code
omg_lat = lat(length_NODC+1:end);
omg_lon = lon(length_NODC+1:end);
nodc_lat = lat(1:length_NODC);
nodc_lon = lon(1:length_NODC);
clf
hold on
daspect([1 aspect_ratio 1])
%scatter(nodc_lon, nodc_lat, 1.5, 'b','filled')
%scatter(omg_lon, omg_lat, 7, 'red','filled')
plot(cx, cy, 'k')
xlim([-100, -10])
ylim([55, 85])
% Plot Mooring Lat Lon on top (optional and not used in FINNEST Figure)
if iscell(lat_moor)
lat_moor = cell2mat(lat_moor) ;
lon_moor = cell2mat(lon_moor) ;
year_moor = cell2mat(year_moor) ;
month_moor = cell2mat(month_moor) ;
end
run = 1 ;
if run == 1
scatter(lon_moor,lat_moor,7,'green','filled')
end
% Export the plot with high resolution
set(gca,'box','on')
%exportgraphics(gcf, 'finesst_plot_mooring.png', 'Resolution', 300); % Save as PNG with 300 DPI
%% TS plot figure for FINESST
%work to setup
load glacier_FINESST.mat glacier
glacier_coords = [66.350, -38.2000] ;
glacier_idx = inpolygon(lon_fj, lat_fj, glacier(:, 1), glacier(:, 2));

hel_fj = mon_fj(glacier_idx) ;
Aug = hel_fj == 8 ;
Sep = hel_fj == 9 ;
%fjord_sal_mat_fj = double(fjord_sal_mat_fj) ;
%fjord_temp_mat_fj = double(fjord_temp_mat_fj) ;

temp_glacier = temp(in_idx) ; % not interpolated vertically
temp_glacier = temp_glacier(:,glacier_idx) ;
sal_glacier = sal(in_idx) ;
sal_glacier = sal_glacier(:,glacier_idx) ;
dep_glacier = dep(in_idx) ;
dep_glacier = dep_glacier(:,glacier_idx) ;

%Spline Interpolate
maxValues = cellfun(@max, dep_glacier) ;
for i = 1:length(dep_glacier)
    DepInterval_temp{i} = (1:maxValues(i))' ; % array until the last depth value
    spline_temp{1,i} = spline(dep_glacier{1,i},temp_glacier{1,i},DepInterval_temp{i}) ; 
    spline_sal{1,i} = spline(dep_glacier{1,i},sal_glacier{1,i},DepInterval_temp{i}) ;
end
% cell to matrix (don't think this works with my approach without adding Nans
%for i = 1:length(spline_temp)
%insidearray1 = spline_temp{i}' ;
%insidearray2 = spline_temp{i}' ;
%spline_temp_mat(:,i) = insidearray1 ;
%spline_sal_mat(:,i) = insidearray2 ;
%end

spline_temp_Sep = spline_temp(:,Sep) ;
spline_sal_Sep = spline_sal(:,Sep) ;
spline_dep_Sep = DepInterval_temp(:,Sep) ;
all_hel_sal = fjord_sal_mat_fj(:,glacier_idx) ;
all_hel_temp = fjord_temp_mat_fj(:,glacier_idx) ;
all_hel_sal_Sep = all_hel_sal(:,Sep) ;
all_hel_sal_Aug = all_hel_sal(:,Aug) ;
all_hel_temp_Sep = all_hel_temp(:,Sep) ;
all_hel_temp_Aug = all_hel_temp(:,Aug) ;

hel_lon_Sep = lon_fj(:,glacier_idx) ;
hel_lat_Sep = lat_fj(:,glacier_idx) ;
hel_lon_Sep = hel_lon_Sep(:,Sep) ;
hel_lat_Sep = hel_lat_Sep(:,Sep) ;

for i = 1:length(hel_lat_Sep)
    temp_gl = [glacier_coords;hel_lat_Sep(i),hel_lon_Sep(i)] ;
    temp_x = temp_gl(:,2) ;
    temp_y = temp_gl(:,1) ;
dist_glacier(1,i) = sw_dist(temp_y,temp_x,'km') ;
end
%five_sal = fjord_sal_mat_fj(51:70,hel_idx) ;
%five_temp = fjord_temp_mat_fj(51:70,hel_idx) ;
%five_filter = all(~isnan(five_sal), 1) ;
%five_sal = five_sal(:,five_filter) ;
%five_temp = five_temp(:,five_filter) ;
%five_dep = (51:70)' ;
%five_sal_sep = five_sal(:,Sep(five_filter)) ;
%five_temp_sep = five_temp(:,Sep(five_filter)) ;
%seven_sal = fjord_sal_mat_fj(71:90,hel_idx) ;
%seven_temp = fjord_temp_mat_fj(71:90,hel_idx) ;
%nine_sal = fjord_sal_mat_fj(91:110,hel_idx) ;
%nine_temp = fjord_temp_mat_fj(91:110,hel_idx) ;
%% Rough plot with too much data
clf
hold on
colorbar;
%caxis([min(dist_glacier), max(dist_glacier)]);  % Set color range
%colormap();

% Get the color for each data point based on dist_glacier
for i = 1:3:size(spline_temp_Sep,2)
    %color_value = dist_glacier(i);  % Scalar value for color
    %color = colormap('jet');  % Get the colormap
    %color_index = round((color_value - min(dist_glacier)) / (max(dist_glacier) - min(dist_glacier)) * (size(color, 1) - 1)) + 1;  % Normalize and index the colormap
    plot(spline_sal_Sep{i}(:,1), spline_temp_Sep{i}(:,1));
end

title('September Near Glacier Helheim TS')
hold off
%plot(seven_sal(:,hel_idx),seven_temp(:,hel_idx))
% clear hel_idx Aug Sep
%% Polished Plot (have to select a plot from the above grouping)
addpath 'C:\Users\ajs82292\Desktop\Research\Matlab\Source\Greenland_Melt\gsw_matlab_v3_06_16'
addpath 'C:\Users\ajs82292\Desktop\Research\Matlab\Source\Greenland_Melt\gsw_matlab_v3_06_16\html'
addpath 'C:\Users\ajs82292\Desktop\Research\Matlab\Source\Greenland_Melt\gsw_matlab_v3_06_16\library'
addpath 'C:\Users\ajs82292\Desktop\Research\Matlab\Source\Greenland_Melt\gsw_matlab_v3_06_16\pdf'
addpath 'C:\Users\ajs82292\Desktop\Research\Matlab\Source\Greenland_Melt\gsw_matlab_v3_06_16\thermodynamics_from_t'
load profile5.mat
load profile6.mat

%find index of profiles selected
five_x = profile5.Target.XData(:) ;
match_profile5 = cellfun(@(c) isequal(c, five_x), spline_sal_Sep);
five_idx = find(match_profile5 == 1) ;

six_x = profile6.Target.XData(:) ;
match_profile6 = cellfun(@(c) isequal(c, six_x), spline_sal_Sep);
six_idx = find(match_profile6 == 1) ;
clear six_x five_x

five_dep = spline_dep_Sep(five_idx) ;
six_dep =  spline_dep_Sep(six_idx) ;
figure
hold on
isopycnals = [17,18,19,20,21,22,23,24,25,26,27,28] ;
gsw_SA_CT_plot_Aiden(profile5.Target.XData,profile5.Target.YData,0,isopycnals,'Helheim Fjord TS plot')
scatter(profile5.Target.XData,profile5.Target.YData,50,five_dep{1},'filled')
scatter(profile6.Target.XData,profile6.Target.YData,50,six_dep{1},'filled')
colormap('jet');
a = colorbar;
a.Label.String = 'Depth (m)';
caxis([1,150])
ylim([-1.5,2])
xlim([20,36])
set(gca,'box','on')
exportgraphics(gcf, 'test.png', 'Resolution', 300); % Save as PNG with 300 DPI
%% bathmetry setup
hold on
isobath_interval = -200:-100:-3000;  % Define isobath intervals down to the minimum depth
iso_interval = 200:100:3000 ;
[~, ContourMatrix] = contour(XX_eto, YY_eto, Z_eto, isobath_interval);
[~,CountourM] = contour(XX,YY,ZZ, iso_interval) ;
daspect([1 aspect_ratio 1])
xlim([-80,-30])
ylim([55,80])
plot(cx,cy,'k')
hold off
clear isobath_interval iso_interval
%% plot fjords and anomaly boxes 
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
for i = 1:length(fjord_box_cords)
    plot(fjord_box_cords{i}(:, 1), fjord_box_cords{i}(:, 2), 'k')
end
plot(cx,cy,'k')
daspect([1 aspect_ratio 1])
xlim([-80,-30])
ylim([55,80])
hold off
clear CountourM ContourMatrix
%% Plot Fjord Anomalies
clf
hold on
for i = 1:length(fj_anoms)
    for j = 1:size(fj_anoms{i},2)
    plot(fj_anoms{i}(:,j),DepInterval)
    end
end
axis ij
ylim([0,300]) 
xlabel('Salinity')
ylabel('Depth M')
title('Salinity Anomalies for Fjord Profiles')
hold off
%%
%Plot Angled Rectangles
clf
hold on
plot(SW_rotvertcell{1}(1,:),SW_rotvertcell{1}(2,:) )
scatter(W_lon(1,1),W_lat(1,1))
plot(cx,cy,' k')
xlim([-60,-45])
ylim([60,70])
plot(x_Wcoast,y_SWcoast, 'r')
plot(x5_Wcoast,y_SWcoast,'r')
hold off
%%
% Test of non-overlapping lat-lons should line up with blue having some
% values with no red but not the alternative
clf
hold on
plot(cx,cy,'k')
daspect([1 aspect_ratio 1])
xlim([-80,-35])
ylim([55,80])
scatter(coastal_lon(:,year_mon_coast),coastal_lat(:,year_mon_coast),'MarkerFaceColor','r')
%scatter(lon_open(:,year_mon_open),lat_open(:,year_mon_open),'MarkerFaceColor','r')
scatter(lon_coast_n(:,yeamon_n_coast),lat_coast_n(:,yeamon_n_coast),'b')
%scatter(lon_open_n(:,yeamon_n_open),lat_open_n(:,yeamon_n_open),'b')
hold off
%% Final Extended Coast plot
clf
hold on
plot(cx,cy,'k')
plot(x_coast,y_coast,'b')
daspect([1 aspect_ratio 1])
xlim([-80,-35])
ylim([55,80])
plot(off_coast(:,1),off_coast(:,2),'r')
scatter(reff_x,reff_y,'g')
%hold off
%% Plot restricted Coastal Casts
clf
hold on
plot(cx,cy,'k')
plot(x_coast,y_coast,'b')
daspect([1 aspect_ratio 1])
xlim([-80,-35])
ylim([55,80])
plot(off_coast(:,1),off_coast(:,2),'r')
scatter(lon(in),lat(in),3,'MarkerFaceColor','g','MarkerEdgeAlpha','.05')
plot(combined_x,combined_y,'r')
hold off
%% Plot Rotated Boxes
clf
hold on
plot(cx,cy,'k')
plot(x_coast,y_coast,'b')
daspect([1 aspect_ratio 1])
xlim([-80,-35])
ylim([55,80])
plot(off_coast(:,1),off_coast(:,2),'b')
scatter(lon(in),lat(in),3,'MarkerFaceColor','g','MarkerEdgeAlpha','.05')
plot(rotvertcell{1}(1,:),rotvertcell{1}(2,:), 'k')
plot(rotvertcell{7}(1,:),rotvertcell{7}(2,:), 'k')
plot(rotvertcell{10}(1,:),rotvertcell{10}(2,:), 'k')
scatter(reff_x,reff_y,'k')
hold off
%% box
box_number = 10000 ;
hold on
clf
hold on
plot(cx,cy,'k')
plot(x_coast,y_coast,'b')
xlim([-80,-35])
ylim([55,80])
plot(vert{box_number}(1,:),vert{box_number}(2,:),'r')
scatter(reff_x,reff_y,'k')
scatter(coastal_lon(1,box_number),coastal_lat(1,box_number),'b.')
scatter(exten_minus(1,1),exten_minus(2,1),'k')
scatter(exten_plus(1,1),exten_plus(2,1),'k')
daspect([1 aspect_ratio 1])
hold off
%%
test = inpolygon(coastal_lon,coastal_lat,vert{1}(1,:),vert{1}(2,:)) ; % test points
hold on
clf
hold on
plot(cx,cy,'k')
plot(x_coast,y_coast,'b')
xlim([-80,-35])
ylim([55,80])
scatter(reff_x,reff_y,'k')
plot(vert{1}(1,:),vert{1}(2,:),'r')
scatter(coastal_lon(test),coastal_lat(test),3,'MarkerFaceColor','b','MarkerEdgeAlpha','.05')
hold off
%% Open Ocean Casts
cast = 74806 ;
clf
hold on
plot(cx,cy,'k')
daspect([1 aspect_ratio 1])
xlim([-80,-35])
ylim([55,80])
scatter(lon_a(open_find{cast}),lat_a(open_find{cast}),10,'MarkerFaceColor','b','MarkerEdgeAlpha','.05')
scatter(lon_open(cast),lat_open(cast),10,'MarkerFaceColor','r','MarkerEdgeAlpha','.05')
hold off
%% TS Diagram (need to select a subsect or something)
clf
hold on
for i = 1:10:length(interp_sal_mat(1,:))
plot(interp_sal_mat(:,i),ptmp_a(:,i),'x','Linestyle','none')
xlabel('Salinity')
ylabel('Potential Temperature')
end
hold off
%% 0.1 x 0.1 boxes
clf
hold on
daspect([1 aspect_ratio 1])
xlim([-80,-35])
ylim([55,80])
cmap = hot ;
for i = 1:length(box_avg)
        % Fill the box with the value
        fill(cast_boxes{i}(2,:), cast_boxes{i}(1,:), box_avg(i), 'EdgeColor', 'none');
end
c = colorbar ;
colormap('cool')
caxis([0 50]);
ylabel(c,'Avg # of profiles')
title('October') % check value of month
plot(cx,cy,'k','MarkerSize',200)
hold off
%% Plot Fjord Casts
clf
hold on
daspect([1 aspect_ratio 1])
xlim([-80,-35])
ylim([55,80])
plot(cx,cy,'k')
scatter(lon_fj,lat_fj,'r') ;
%% Plot Fjord Salinites colored by groupings
colors = lines(40);
clf
hold on
for i = 1:length(in_fj)
    salinityData = fjord_sal(in_fj{i});
    for j = 1:length(salinityData)
plot(salinityData{j}, DepInterval, 'Color', colors(i, :))
    end
end
axis ij
hold off
%% PCA correct (includes coeff and sal_anom/salinity) Coastal!
folderPath = 'C:\Users\ajs82292\Desktop\Research\Weekly Meeting\Images\10-24-24' ; % change depending on folderlocation
DepInterval_custom = DepInterval(1:length(sal_anom)) ;
hold on
daspect([1 aspect_ratio 1])
xlim([-80,-30])
ylim([55,80])
scatter(lon_coast(:,year_mon_coast), lat_coast(:,year_mon_coast), 50, first_PC_anom, 'filled');
colorbar;
% caxis([-3 3]);
colormap('cool')
plot(cx,cy,'k') ;
xlabel('Longitude');
ylabel('Latitude');
title(sprintf('First Principal Component of Salinity Anomaly for %s %d', month_string{month_selected}, year_selected));
hold off
filename = sprintf('coast_anom_PC_%s_%d.jpg',month_string{month_selected},year_selected) ;
print(fullfile(folderPath, filename), '-djpeg') ;
% Plot Coeff
figure
hold on
plot(first_coeff_anom,DepInterval_custom)
title('1st Coefficients vs Depth ')
axis ij
hold off
filename = sprintf('coast_anom_coeff_%s_%d.jpg',month_string{month_selected},year_selected) ;
print(fullfile(folderPath, filename), '-djpeg') ;
% Plot sal_anom
figure
hold on
plot(sal_anom,DepInterval_custom) ;
ylim([0,300])
axis ij
hold off
title('0-300m Salinity Anomaly')
xlabel('Salinity Anomaly')
ylabel('Depth')
hold off
filename = sprintf('coast_sal_anom_%s_%d.jpg',month_string{month_selected},year_selected) ;
print(fullfile(folderPath, filename), '-djpeg') ;
% Create Scree Plot
figure
hold on
plot(1:length(explained_anom), explained_anom, 'o-', 'LineWidth', 2);
xlabel('Principal Component');
xlim([0,5])
xticks(1:5); 
ylabel('Variance Explained (%)');
title('Scree Plot');
grid on;
hold off
filename = sprintf('coast_anom_scree_%s_%d.jpg',month_string{month_selected},year_selected) ;
print(fullfile(folderPath, filename), '-djpeg') ;
clear filename
%% Sal Anom Open PCAS
folderPath = 'C:\Users\ajs82292\Desktop\Research\Weekly Meeting\Images\10-24-24' ; % change depending on folderlocation
DepInterval_custom = DepInterval(1:length(sal_anom)) ;
hold on
daspect([1 aspect_ratio 1])
xlim([-80,-30])
ylim([55,80])
scatter(lon_open(:,year_mon_open), lat_open(:,year_mon_open), 50, first_PC_anom, 'filled');
colorbar;
% caxis([-3 3]);
colormap('cool')
plot(cx,cy,'k') ;
xlabel('Longitude');
ylabel('Latitude');
title(sprintf('First Principal Component of Salinity Anomaly for %s %d', month_string{month_selected}, year_selected));
hold off
filename = sprintf('open_anom_PC_%s_%d.jpg',month_string{month_selected},year_selected) ;
print(fullfile(folderPath, filename), '-djpeg') ;
% Plot Coeff
figure
hold on
plot(first_coeff_anom,DepInterval_custom)
title('1st Coefficients vs Depth ')
axis ij
hold off
filename = sprintf('open_anom_coeff_%s_%d.jpg',month_string{month_selected},year_selected) ;
print(fullfile(folderPath, filename), '-djpeg') ;
% Plot sal_anom
figure
hold on
plot(sal_anom,DepInterval_custom) ;
ylim([0,300])
axis ij
hold off
title('0-300m Salinity Anomaly')
xlabel('Salinity Anomaly')
ylabel('Depth')
hold off
filename = sprintf('open_sal_anom_%s_%d.jpg',month_string{month_selected},year_selected) ;
print(fullfile(folderPath, filename), '-djpeg') ;
% Create Scree Plot
figure
hold on
plot(1:length(explained_anom), explained_anom, 'o-', 'LineWidth', 2);
xlabel('Principal Component');
xlim([0,5])
xticks(1:5); 
ylabel('Variance Explained (%)');
title('Scree Plot');
grid on;
hold off
filename = sprintf('open_anom_scree_%s_%d.jpg',month_string{month_selected},year_selected) ;
print(fullfile(folderPath, filename), '-djpeg') ;
clear filename
%% Salinity PCA's (coastal only!)
folderPath = 'C:\Users\ajs82292\Desktop\Research\Weekly Meeting\Images\10-24-24' ; % change depending on folderlocation
DepInterval_custom = DepInterval(1:300) ;
hold on
daspect([1 aspect_ratio 1])
xlim([-80,-30])
ylim([55,80])
scatter(lon_coast_n(:,yeamon_n_coast), lat_coast_n(:,yeamon_n_coast), 50, first_PC_n, 'filled');
colorbar;
% caxis([-3 3]);
colormap('cool')
plot(cx,cy,'k') ;
xlabel('Longitude');
ylabel('Latitude');
title(sprintf('First Principal Component of Salinity for coastal casts %s %d', month_string{month_selected}, year_selected));
hold off
filename = sprintf('sal_coast_PC_%s_%d.jpg',month_string{month_selected},year_selected) ;
print(fullfile(folderPath, filename), '-djpeg') ;
% Plot Coeff
figure
hold on
plot(first_coeff_n,DepInterval_custom)
title('1st Coefficients vs Depth ')
axis ij
hold off
filename = sprintf('sal_coast_coeff_%s_%d.jpg',month_string{month_selected},year_selected) ;
print(fullfile(folderPath, filename), '-djpeg') ;
% Create Scree Plot
figure
hold on
plot(1:length(explained_n), explained_n, 'o-', 'LineWidth', 2);
xlabel('Principal Component');
xlim([0,5])
xticks(1:5); 
ylabel('Variance Explained (%)');
title('Scree Plot');
grid on;
hold off
filename = sprintf('sal_coast_scree_%s_%d.jpg',month_string{month_selected},year_selected) ;
print(fullfile(folderPath, filename), '-djpeg') ;
clear filename
% Plot Salinity
if size(coastal_sal,2) < 301
coastal_sal = coastal_sal' ;
open_sal = open_sal' ;
end
figure
hold on
plot(coastal_sal(:,yeamon_n_coast),DepInterval_custom) ;
title(sprintf('Salinity vs Depth for coastal casts %s %d', month_string{month_selected}, year_selected))
axis ij
xlabel('Salinity')
ylabel('Depth')
hold off
filename = sprintf('sal_coast_%s_%d.jpg',month_string{month_selected},year_selected) ;
print(fullfile(folderPath, filename), '-djpeg') ;
clear step_s sal_combined_reduced
if size(coastal_sal,2) > 301
coastal_sal = coastal_sal' ;
open_sal = open_sal' ;
end
%% Salinity PCA's (open only!)
folderPath = 'C:\Users\ajs82292\Desktop\Research\Weekly Meeting\Images\10-24-24' ; % change depending on folderlocation
DepInterval_custom = DepInterval(1:300) ;
hold on
daspect([1 aspect_ratio 1])
xlim([-80,-30])
ylim([55,80])
scatter(lon_open_n(:,yeamon_n_open), lat_open_n(:,yeamon_n_open), 50, first_PC_n, 'filled');
colorbar;
% caxis([-3 3]);
colormap('cool')
plot(cx,cy,'k') ;
xlabel('Longitude');
ylabel('Latitude');
title(sprintf('First Principal Component of open ocean Salinity for %s %d', month_string{month_selected}, year_selected));
hold off
filename = sprintf('sal_open_PC_%s_%d.jpg',month_string{month_selected},year_selected) ;
print(fullfile(folderPath, filename), '-djpeg') ;
% Plot Coeff
figure
hold on
plot(first_coeff_n,DepInterval_custom)
title('1st Coefficients vs Depth ')
axis ij
hold off
filename = sprintf('sal_open_coeff_%s_%d.jpg',month_string{month_selected},year_selected) ;
print(fullfile(folderPath, filename), '-djpeg') ;
% Create Scree Plot
figure
hold on
plot(1:length(explained_n), explained_n, 'o-', 'LineWidth', 2);
xlabel('Principal Component');
xlim([0,5])
xticks(1:5); 
ylabel('Variance Explained (%)');
title('Scree Plot');
grid on;
hold off
filename = sprintf('sal_open_scree_%s_%d.jpg',month_string{month_selected},year_selected) ;
print(fullfile(folderPath, filename), '-djpeg') ;
clear filename
% Plot Salinity
if size(open_sal,2) < 301
open_sal = open_sal' ;
end
figure
hold on
plot(open_sal(:,yeamon_n_open),DepInterval_custom) ;
title(sprintf('Salinity vs Depth for open ocean %s %d', month_string{month_selected}, year_selected))
axis ij
xlabel('Salinity')
ylabel('Depth')
hold off
filename = sprintf('sal_open_%s_%d.jpg',month_string{month_selected},year_selected) ;
print(fullfile(folderPath, filename), '-djpeg') ;
clear step_s sal_combined_reduced
if size(open_sal,2) > 301
open_sal = open_sal' ;
end
%% Salinity PCA (fjords)
clf
folderPath = 'C:\Users\ajs82292\Desktop\Research\Weekly Meeting\Images\12-12-24' ; % change depending on folderlocation
DepInterval_custom_change = DepInterval(1:length(first_coeff_fj)) ;
DepInterval_custom = DepInterval(1:300) ;
hold on
daspect([1 aspect_ratio 1])
xlim([-80,-30])
ylim([55,80])
scatter(lon_fj(:,year_mon_fj), lat_fj(:,year_mon_fj), 50, first_PC_fj, 'filled');
colorbar;
% caxis([-3 3]);
colormap('hot')
plot(cx,cy,'k') ;
xlabel('Longitude');
ylabel('Latitude');
title(sprintf('First Principal Component of Fjord Salinity for %s %d - %d', month_string{month_selected}, five_year_range(1),five_year_range(2)));
hold off
filename = sprintf('sal_fjord_PC_%s_%d.jpg',month_string{month_selected},year_selected) ;
print(fullfile(folderPath, filename), '-djpeg') ;
% Plot Coeff
figure
hold on
plot(first_coeff_fj,DepInterval_custom_change)
title('1st Coefficients vs Depth ')
axis ij
hold off
filename = sprintf('sal_fj_coeff_%s_%d.jpg',month_string{month_selected},year_selected) ;
print(fullfile(folderPath, filename), '-djpeg') ;
% Create Scree Plot
figure
hold on
plot(1:length(explained_fj), explained_fj, 'o-', 'LineWidth', 2);
xlabel('Principal Component');
xlim([0,5])
xticks(1:5); 
ylabel('Variance Explained (%)');
title('Scree Plot');
grid on;
hold off
filename = sprintf('sal_fj_scree_%s_%d.jpg',month_string{month_selected},year_selected) ;
print(fullfile(folderPath, filename), '-djpeg') ;
clear filename
% Plot Salinity
if size(fj_combined,2) < 301
fj_combined = fj_combined' ;
end
figure
hold on
plot(fj_combined(:,year_mon_fj),DepInterval_custom) ;
title(sprintf('Salinity vs Depth for open ocean %s %d - %d', month_string{month_selected}, five_year_range(1),five_year_range(2))) ;
axis ij
xlabel('Salinity')
ylabel('Depth')
hold off
filename = sprintf('sal_open_%s_%d.jpg',month_string{month_selected},year_selected) ;
print(fullfile(folderPath, filename), '-djpeg') ;
clear step_s sal_combined_reduced
if size(fj_combined,2) > 301
fj_combined = fj_combined' ;
end
%% Temperature PCA (fjords)
clf
folderPath = 'C:\Users\ajs82292\Desktop\Research\Weekly Meeting\Images\12-12-24' ; % change depending on folderlocation
DepInterval_custom_change = DepInterval(1:length(first_coeff_fj_t)) ;
DepInterval_custom = DepInterval(1:300) ;
hold on
daspect([1 aspect_ratio 1])
xlim([-80,-30])
ylim([55,80])
scatter(lon_fj(:,year_mon_fj), lat_fj(:,year_mon_fj), 50, first_PC_fj_t, 'filled');
colorbar;
% caxis([-3 3]);
colormap('hot')
plot(cx,cy,'k') ;
xlabel('Longitude');
ylabel('Latitude');
title(sprintf('First Principal Component of Fjord Temperature for %s %d - %d', month_string{month_selected}, five_year_range(1),five_year_range(2)));
hold off
filename = sprintf('temp_fjord_PC_%s_%d.jpg',month_string{month_selected},year_selected) ;
print(fullfile(folderPath, filename), '-djpeg') ;
% Plot Coeff
figure
hold on
plot(first_coeff_fj_t,DepInterval_custom_change)
title('1st Coefficients vs Depth ')
axis ij
hold off
filename = sprintf('temp_fj_coeff_%s_%d.jpg',month_string{month_selected},year_selected) ;
print(fullfile(folderPath, filename), '-djpeg') ;
% Create Scree Plot
figure
hold on
plot(1:length(explained_fj_t), explained_fj_t, 'o-', 'LineWidth', 2);
xlabel('Principal Component');
xlim([0,5])
xticks(1:5); 
ylabel('Variance Explained (%)');
title('Scree Plot');
grid on;
hold off
filename = sprintf('temp_fj_scree_%s_%d.jpg',month_string{month_selected},year_selected) ;
print(fullfile(folderPath, filename), '-djpeg') ;
clear filename
% Plot Temperature
if size(fj_temp_combined,2) < 301
fj_temp_combined = fj_temp_combined' ;
end
figure
hold on
plot(fj_temp_combined(:,year_mon_fj),DepInterval_custom) ;
title(sprintf('Temperature vs Depth for fjord %s %d - %d', month_string{month_selected}, five_year_range(1),five_year_range(2))) ;
axis ij
xlabel('Temperature (C)')
ylabel('Depth')
hold off
filename = sprintf('temp_fj_%s_%d.jpg',month_string{month_selected},year_selected) ;
print(fullfile(folderPath, filename), '-djpeg') ;
clear step_s sal_combined_reduced
if size(fj_temp_combined,2) > 301
fj_temp_combined = fj_combined' ;
end
%% EOF 1 vs 2 for coast, open, and fjord
%setup of GMM
k = 1 ; % # number of clusters
c_sal_combined = [first_PC_coast_sal,second_PC_coast_sal] ;
c_model = fitgmdist(c_sal_combined,k) ;
c_idx = cluster(c_model,c_sal_combined) ;

o_sal_combined = [first_PC_open_sal,second_PC_open_sal] ;
o_model = fitgmdist(o_sal_combined,k+1) ;
o_idx = cluster(o_model,o_sal_combined) ;

% start plot
clf
sgtitle(sprintf('Salinity %s %d',month_string{month_selected},year_selected))
%coast
subplot(2,2,1)
hold on
gscatter(c_sal_combined(:,1), c_sal_combined(:,2), c_idx, 'br', 'o','filled');
%scatter(first_PC_coast_sal,second_PC_coast_sal,'b','filled')
%scatter selected points using idx
%scatter(first_PC_coast_sal(c_idx),second_PC_coast_sal(c_idx),'r','filled')
xlabel('First PC')
ylabel('Second PC') 
title('Coastal Profiles')

% Coastal Salinities
DepInterval_custom = DepInterval(1:size(PC_coast_sal,2))' ;
subplot(2,2,3)
hold on
for i = 1:size(PC_coast_sal,1)
    plot(PC_coast_sal(i,:),DepInterval_custom,'b')
end
%plot selected profiles using idx
c2_idx = c_idx == 2 ;
temp_sal = PC_coast_sal(c2_idx,:) ;
for i = 1:size(temp_sal,1)
    plot(temp_sal(i,:),DepInterval_custom,'r')
end
axis ij
xlabel('Salinity (psu)')
ylabel('Depth (m)')

%open
subplot(2,2,2)
hold on
gscatter(o_sal_combined(:,1), o_sal_combined(:,2), o_idx, 'br', 'o','filled');
%scatter selected points using idx
xlabel('First PC')
ylabel('Second PC')
title('Open Profiles')

%open salinties
DepInterval_custom = DepInterval(1:size(PC_open_sal,2))' ;
subplot(2,2,4)
hold on
for i = 1:size(PC_open_sal,1)
    plot(PC_open_sal(i,:),DepInterval_custom,'b')
end
%plot selected profiles using idx
o2_idx = o_idx == 2 ;
temp_sal = PC_open_sal(o2_idx,:) ;
for i = 1:size(temp_sal,1)
    plot(temp_sal(i,:),DepInterval_custom,'r')
end
axis ij
xlabel('Salinity (psu)')
ylabel('Depth (m)')
clear temp_sal c_idx o_idx
%% Reconstruct a salinity profile using the 1st principal component
DepInterval_custom = DepInterval(1:300) ;
number = 1 ; % which profile
number_profile = coastal_sal(yeamon_n_coast,:) ;
run = 2 ;
if run == 1
coeff_n = coeff_n';
end
z_reconstructed = score_n(number,1) * coeff_n(:,1)' ;  % Multiply the score by the transpose of coeff
x_reconstructed = z_reconstructed .* std_sal(number) + mean_sal(number) ;
z_reconstructed_2 = score_n(number, 1:2) * coeff_n(:, 1:2)';  % for PC's 1/2
x_reconstructed_2 = z_reconstructed_2 .* std_sal(number,:) + mean_sal(number,:) ;
z_3 = score_n(number, 1:3) * coeff_n(:, 1:3)';  % for PC's 1/2/3
x_3 = z_3 .* std_sal(number,:) + mean_sal(number,:) ;
hold on
plot(number_profile(1,:),DepInterval_custom,'b')
plot(x_reconstructed,DepInterval_custom,'r')
plot(x_reconstructed_2,DepInterval_custom,'g')
plot(x_3,DepInterval_custom,'k')
axis ij
legend({'Raw Salinity', '1st PC Only', 'Top 2 PCs','Top 3'});
hold off
clear number_profile number z_reconstructed z_reconstructed_2 x_reconstructed x_reconstructed_2 x_3 z_3
%% Plot 1/10 salinities
DepInterval_custom = DepInterval(1:300) ;
step_s = 10 ;
reduced = sal_combined(:,1:step_s:end) ;
figure
plot(reduced,DepInterval_custom) ;
xlim([0,40])
axis ij
xlabel('Salinity');
ylabel('Depth');
title('Salinity profiles with new cleaning changes');
%% PLot individual fjords salinity and salinity anomaly
fjord = 32 ;
DepInterval_custom = DepInterval(1:300) ;
temp_idx = inpolygon(lon_fj,lat_fj,fjord_vert{fjord}(:,1),fjord_vert{fjord}(:,2)) ;
temp_box_idx = inpolygon(lon,lat,fjord_box_cords{fjord}(:,1),fjord_box_cords{fjord}(:,2)) ;
temp_lon = lon_fj(temp_idx) ;
temp_lat = lat_fj(temp_idx) ;
if size(fj_combined,2) < 301
fj_combined = fj_combined' ;
end
temp_fj = fj_combined(:,temp_idx) ;
hold on
plot(temp_fj,DepInterval_custom)
axis ij
title('Helheim Fjord Salinities')
if size(fj_combined,2) > 301
fj_combined = fj_combined' ;
end
hold off
%anomalies
figure
if size(fj_anom_combined,2) < 301
fj_anom_combined = fj_anom_combined' ;
end
temp_anom = fj_anom_combined(:,temp_idx) ;
hold on
plot(temp_anom,DepInterval_custom)
axis ij
title('Helheim Fjord Salinity Anomalies')
%map
figure
hold on
daspect([1 aspect_ratio 1])
fill(fjord_vert{fjord}(:,1),fjord_vert{fjord}(:,2),'b')
plot(fjord_box_cords{fjord}(:,1),fjord_box_cords{fjord}(:,2))
plot(cx,cy,'k')
scatter(temp_lon,temp_lat,10,'red','filled')
scatter(lon(temp_box_idx),lat(temp_box_idx),5,'b','filled')
clear temp_idx temp_lon temp_lat
%% Plot abs derivatives (line 622)
hold on
step_s = 5 ;
xline_half = 0.5 * ones(1,length(DepInterval)) ;
xline_half_2 = 0.3 * ones(1,length(DepInterval)) ;
xline_two = 2 * ones(1,length(DepInterval)) ;
reduced = diff_result(:,1:step_s:end) ;
plot(reduced,DepInterval) ;
plot(xline_half,DepInterval,'--r')
plot(xline_half_2,DepInterval,'--g')
plot(xline_two,DepInterval,'--b')
xlim([0,40])
axis ij
xlabel('Salinity Derivatives');
ylabel('Depth');
title('Salinity Derivatives vs Depth');
xlim([0,5])
ylim([0,300]) 
hold off
clear step_s reduced xline_half xline_two
%% Plot fjord derivatives (line 935)
figure
set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]); % Full-screen
folderPath = 'C:\Users\ajs82292\Desktop\Research\Weekly Meeting\Images\11-21-24' ; % change depending on folderlocation
hold on
xline_25 = threshold_25 * ones(1,25) ;
xline_half_2 = 0.3 * ones(1,length(DepInterval)-50) ;
xline_two = 2 * ones(1, numel(DepInterval(26:50)));
h1 = plot(diff_result,DepInterval,'b') ;
h2= plot(diff_result(:,OMG_fj_profiles),DepInterval,'r') ;
plot(xline_half_2,DepInterval(51:end),'--k')
plot(xline_two,DepInterval(26:50),'--k')
plot(xline_25,DepInterval(1:25),'--k')
xlim([0,40])
axis ij
xlabel('Salinity Derivatives');
ylabel('Depth');
title('Salinity Derivatives vs Depth');
legend([h1(1), h2(1)], 'NODC Profiles', 'OMG Profiles'); % Explicit legend assignment
xlim([0,12])
ylim([0,300]) 
hold off
filename = sprintf('deriv t25 %d 25-50 %f b50 %c', threshold_25,threshold_50, threshold);
filename = strcat(filename, '.jpeg');
print(fullfile(folderPath, filename), '-djpeg') ;
clear xline_half xline_two h1 h2 filename folderpath
% Plot profiles w removed data from Fjord profiles (935)
figure
set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]); % Full-screen
folderPath = 'C:\Users\ajs82292\Desktop\Research\Weekly Meeting\Images\11-21-24' ; % change depending on folderlocation
hold on
h1 = plot(fjord_sal_mat_fj(:,remove_fj_any), DepInterval, 'b'); % Plot NODC Profiles
h2 = plot(fjord_sal_mat_fj(:,compare_OMG), DepInterval, 'r');   % Plot OMG Profiles
xlim([0,40])
axis ij
xlabel('Salinity');
ylabel('Depth');
title('Profiles with Data Points Cleaned vs Depth');
legend([h1(1), h2(1)], 'NODC Profiles', 'OMG Profiles'); % Explicit legend assignment
%xlim([0,5])
%ylim([0,300]) 
hold off
filename = sprintf('salinity t25 %d 25-50 %f b50 %c',threshold_25, threshold_50, threshold);
filename = strcat(filename, '.jpeg');
print(fullfile(folderPath, filename), '-djpeg') ;
clear h1 h2 filename folderPath
%% Plot Eliminated Profiles due to derivative threshold (line 669 stop) (with dots at the point on the line that will get elimed
figure
hold on
step_s = 10 ;
col = any(remove, 1);
reduced = interp_sal_mat(:,col) ;
remove_reduced = remove(:,col) ;
remove_reduced = remove(:, 1:step_s:end);  
[row_idx, col_idx] = find(remove_reduced == 1) ;
reduced_mat = interp_sal_mat(:,1:step_s:end) ;
plot(reduced_mat,DepInterval) ;
for i = 1:length(row_idx)
    scatter(reduced_mat(row_idx(i), col_idx(i)), DepInterval(row_idx(i)), 'g', 'filled');
end
xlim([0,40])
ylim([0,50])
axis ij
xlabel('Salinity');
ylabel('Depth');
title('Salinity profiles with points eliminated');
hold off
%clear reduced step_s col row_idx col_idx col valid_idx
%% Fjord Anomaly Derivatives
DepInterval_custom = DepInterval(1:300) ;
%result = diff_result(:,temp_idx) ;
threshold_point25 = repmat(threshold_25,length(DepInterval)) ;
threshold_point = repmat(threshold_50,length(DepInterval)) ;
threshold_point2 = repmat(threshold,length(DepInterval)) ;
figure
hold on
plot(diff_result,DepInterval)
plot(threshold_point25,DepInterval,'--k')
plot(threshold_point,DepInterval,'--k')
plot(threshold_point2,DepInterval,'--k')
axis ij
ylim([0,300])
xlabel('Salinity Anomaly');
ylabel('Depth');
title('Cleaned Fjord Sal Anom Derivatives');
hold off
clear temp_idx sal fjord threshold_point2 threshold_point
%% Fjord Spiciness Profiles
fjord = 32 ; 
temp_idx = inpolygon(lon_fj,lat_fj,fjord_vert{fjord}(:,1),fjord_vert{fjord}(:,2)) ;
% seperate by month to see if we can see a noticeable difference
clf
hold on
plot(spice_fj(:,temp_idx),DepInterval)
axis ij
xlabel('Spiciness (kg/m^3)')
ylabel('Depth (m)')
title('Helheim Spiciness')
hold off
clear temp_idx
%% Plot profiles with eliminated points in top 50 m's geogrpahic location
top_50 = remove(1:50,:) ;
col = any(top_50, 1);
figure
hold on 
daspect([1 aspect_ratio 1])
scatter(lon_a(col),lat_a(col),'r')
plot(cx,cy,'k') ;
xlim([-80,-30])
ylim([55,80])
xlabel('Longitude');
ylabel('Latitude');
title('Location of profiles with points removed in top 50m')
hold off
%% Find and plot salinites with 0-1 m measurements >15 (should only be very near coast measurements) (eventually add back in if it produces good results
sal_combined = sal_combined(1:2,:) ;
[~,col] = find(sal_combined <= 15) ;
unique_col = unique(col, 'stable');
clear col
%plot
hold on
scatter(lon_a(unique_col),lat_a(unique_col),'r') 
plot(cx,cy,'k') ;
xlim([-80,-30])
ylim([55,80])
xlabel('Longitude');
ylabel('Latitude');
daspect([1 aspect_ratio 1])
hold off

