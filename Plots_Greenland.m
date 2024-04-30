% Plot Script 
%Plot Coast Line With extended area
clf
hold on
plot(cx,cy,'k')
scatter(lon,lat,3,'MarkerFaceColor','g','MarkerEdgeAlpha','.05')
xlim([-80,-35])
ylim([55,80])
daspect([1 aspect_ratio 1])
plot(Southgcoast(:,1),Southgcoast(:,2),'r')
plot(Southgcoast(:,1),y5_Scoast,'r')
plot(Southgcoast(:,1),extendedcoast_S,'--b')
plot(Westgcoast(:,1),Westgcoast(:,2),'r')
plot(x5_Wcoast,Westgcoast(:,2),'r')
plot(extendedcoast_W,Westgcoast(:,2),'--b')
plot(Eastgcoast(:,1),Eastgcoast(:,2),'r')
plot(x5_Ecoast,Eastgcoast(:,2),'r')
plot(extendedcoast_E,Eastgcoast(:,2),'--b')
plot(Gap_plus(:,1),Gap_plus(:,2), 'r')
plot(Gap_minus(:,1),Gap_minus(:,2), 'r')
hold off
%% Plot Section on Line 131
clf
hold on 
plot(x_Wcoast,y_SWcoast, 'r')
plot(x5_Wcoast,y_SWcoast,'r')
plot(cx,cy, 'k')
xlim([-60,-45])
ylim([60,70])
scatter(W_lon_coast,W_lat_coast,1,'MarkerFaceColor','b','MarkerEdgeAlpha','0') ;
scatter(W_lon,W_lat,1,'MarkerFaceColor','r','MarkerEdgeAlpha','0') ;
hold off
clf
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