% Plot coastline
clf
hold on
plot(cx,cy,'k')
scatter(lon,lat,3,'MarkerFaceColor','g','MarkerEdgeAlpha','.05')
xlim([-80,-35])
ylim([55,80])
daspect([1 aspect_ratio 1])
plot(x_coast,y_coast,'r')
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
% Test of Inpolygon
clf
hold on
plot(cx,cy,'k')
daspect([1 aspect_ratio 1])
xlim([-80,-35])
ylim([55,80])
plot(polygon_x,polygon_y,'b')
scatter(exten_plus(~in_plus,1),exten_plus(~in_plus,2),'r')
scatter(exten_plus2(~in_plus2,1),exten_plus2(~in_plus2,2),'r')
scatter(exten_minus(~in_minus,1),exten_minus(~in_minus,2),'r')
scatter(exten_minus(~in_minus2,1),exten_minus2(~in_minus2,2),'r')
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
hold on
clf
hold on
plot(cx,cy,'k')
plot(x_coast,y_coast,'b')
xlim([-80,-35])
ylim([55,80])
scatter(reff_x,reff_y,'k')
scatter(coastal_lon(1,1),coastal_lat(1,1),'c')
scatter(vert{1}(1,1),vert{1}(2,1),'r')
scatter(vert{1}(1,2),vert{1}(2,2),'b')
scatter(vert{1}(1,3),vert{1}(2,3),'y')
scatter(vert{1}(1,4),vert{1}(2,4),'g')
scatter(exten_minus(1,1),exten_minus(2,1),'k')
scatter(exten_plus(1,1),exten_plus(2,1),'k')
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
%%
clf
hold on
plot(cx,cy,'k')
daspect([1 aspect_ratio 1])
xlim([-80,-35])
ylim([55,80])
plot(polygon_x,polygon_y,'r')
scatter(lon(~in_idx),lat(~in_idx),3,'MarkerFaceColor','b','MarkerEdgeAlpha','.05')
hold off
%% Open Ocean Casts
clf
hold on
plot(cx,cy,'k')
daspect([1 aspect_ratio 1])
xlim([-80,-35])
ylim([55,80])
scatter(lon_a(cast_idx_open{2}),lat_a(cast_idx_open{2}),3,'MarkerFaceColor','b','MarkerEdgeAlpha','.05')
scatter(lon_a(2),lat_a(2),3,'MarkerFaceColor','r','MarkerEdgeAlpha','.05')
hold off
%% Plot Anomaly example w/ std dev
number = 1006 ;
clf
hold on
daspect([1 aspect_ratio 1])
xlim([-80,-35])
ylim([55,80])
plot(cx,cy,'k')
scatter(lon_a(cast_idx_coast{number}),lat_a(cast_idx_coast{number}),10,'MarkerFaceColor','b','MarkerEdgeAlpha','0.5')
scatter(coastal_lon(number),coastal_lat(number),10,'MarkerFaceColor','r','MarkerEdgeAlpha','0.5')
hold off
%% continued (Temp anom)
clf
hold on
plot(coast_temp_anom(:,number),DepInterval,'k')
axis ij
xlabel('Temperature Anomaly')
ylabel('Depth (M)')
title('Temperature Anomaly vs Depth')
hold off
%% contin (sal anom)
clf
hold on
plot(coast_sal_anom(:,number),DepInterval,'k')
axis ij
xlabel('Salinity Anomaly')
ylabel('Depth (M)')
title('Salinity Anomaly vs Depth')
hold off
%% contin (std dev sal)
clf
hold on
plot(coastal_sal_std{:,number},DepInterval,'k')
axis ij
xlabel('Salinity Std')
ylabel('Depth (M)')
title('Salinity Std vs Depth')
hold off
%% contin (std dev temp)
clf
hold on
plot(coastal_temp_std{:,number},DepInterval,'k')
axis ij
xlabel('Temperature Std')
ylabel('Depth (M)')
title('Temperature Std vs Depth')
hold off
%% Meshgrid for std (do one where every cast is included, and then seperate by month)
idx = Aug ; % for easy changing of month idx
tri = delaunay(lon_comb(idx), lat_comb(idx)) ;
clf 
hold on
trisurf(tri,lon_comb(idx),lat_comb(idx),count_comb(idx),'FaceColor', 'interp', 'EdgeColor', 'none') ;
daspect([1 aspect_ratio 1])
xlim([-80,-35])
ylim([55,80])
cb= colorbar()
caxis([0 20]);
plot(cx,cy,'k','MarkerSize',200)
xlabel('Lon');
ylabel('Lat');
ylabel(cb, '# of casts contained'); % Label the colorbar
zlabel('# of casts contained');
title('Aug casts')
hold off
