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