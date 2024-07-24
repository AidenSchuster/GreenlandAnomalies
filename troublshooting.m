load( "bedmachinev5_decimated.mat")
load topo_east_coast.mat

nanElements = isnan(depth) ;
num_bed = numel(depth) - sum(nanElements(:)) ;


% topo
te=find(YY(1,:)<53);
YY(:,te)=[];
XX(:,te)=[];
ZZ(:,te)=[];

te=find(XX(:,1)<-84);
YY(te,:)=[];
XX(te,:)=[];
ZZ(te,:)=[];
nanElements = isnan(ZZ) ;
num_topo = numel(ZZ) - sum(nanElements(:)) ;

clf
hold on
plot(XX,YY,'b.')
plot(XX_eto,YY_eto,'g.')
plot(XX_bed,YY_bed,'r.')
xlim([-80,-35])
ylim([55,80])
plot(cx,cy,'k')
aspect_ratio = cosd(65) ; % Aspect Ratio at 65 N
daspect([1 aspect_ratio 1])
hold off

