% Plot coastline
clf
hold on
plot(cx,cy,'k')
scatter(lon,lat,3,'MarkerFaceColor','b','MarkerEdgeAlpha','.05')
scatter(OMG_lon,OMG_lat,7,'MarkerFaceColor','r','MarkerEdgeAlpha','.05')
xlim([-75,-30])
ylim([55,80])
daspect([1 aspect_ratio 1])
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
hold off
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
%% find sal anom issue
test = sal_anom_combined >= 5 | sal_anom_combined <= -5 ;
column_keep = any(test,1) ;
sal_anom_combined(:,~column_keep) = [] ;
lat_combined(:,~column_keep) = [] ;
lon_combined(:,~column_keep) = [] ;
%plot 
hold on
daspect([1 aspect_ratio 1])
xlim([-80,-30])
ylim([55,80])
scatter(lon_combined,lat_combined,50,'k')
plot(cx,cy,'k') ;
xlabel('Longitude');
ylabel('Latitude');
hold off
%% example bad data with profile 80,064
hold on 
daspect([1 aspect_ratio 1])
xlim([-80,-30])
ylim([55,80])
scatter(lon_combined(80064),lat_combined(80064),50,'k')
plot(cx,cy,'k') ;
xlabel('Longitude');
ylabel('Latitude');
hold off

hold on 
plot(open_sal(:,80064),DepInterval) ;
xlabel('Salinity');
ylabel('Depth');
axis ij
hold off