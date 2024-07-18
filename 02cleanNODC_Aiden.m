clear all
close all


cd C:\UGA\IMCS\NODC\greenland

load 01GREENLAND_temp_sal
whos yea mon day lon lat dep watdep temp sal

watdep=double(watdep);
lat=double(lat);
lon=double(lon);
    
% First, find true water depth
load topo_east_coast


te=find(YY(1,:)<53);
YY(:,te)=[];
XX(:,te)=[];
ZZ(:,te)=[];

te=find(XX(:,1)<-84);
YY(te,:)=[];
XX(te,:)=[];
ZZ(te,:)=[];

    
TUTU=1
if TUTU==2
twd=griddata(XX,YY,ZZ,lon,lat);
save 02cleanNODC_temporary twd
end
load 02cleanNODC_temporary

plot(lon,lat,'.')
hold on
    
% Start cleaning data. First, automated methods used by
% Lentz (2003), Climatology of salty intrusions

% 1) Eliminating profiles where deepest sample exceeded the water depth by
% more than 10 m

% NOTE TO AIDEN: YOU MAY HAVE TO UPDATE THE TOPOGRAPHY HERE WITH
% BEDMACHINE. THIS WORKS FINE OFFSHORE, BUT IN THE FJORDS I AM AFRAID THIS
% WILL DELETE A LOT OF DATA BECAUSE THE TOPOGRAPHY DATA INSIDE THE FJORDS
% IN THIS DATA SET IS COMPLETELY WRONG. SO YOU CAN EITHER NOT APPLY THIS
% CRITERIUM INSIDE THE FJORDS, OR UPDATE THE TOPOGRAPHY WITH BEDMACHINE
% WHEN FINDING twd (total water depth) INSIDE THE FJORDS

dif=watdep-twd;
ind=find(dif<10);

whos yea mon day lon lat dep watdep temp sal

yea=yea(ind);
mon=mon(ind);
day=day(ind);
lon=lon(ind);
lat=lat(ind);
dep=dep(ind);
watdep=watdep(ind);
temp=temp(ind);
sal=sal(ind);
twd=twd(ind);

whos yea mon day lon lat dep watdep temp sal twd

plot(lon,lat,'.r')
hold on

% 2) 
% Now, clear profiles with the following characteristics:
% 0  and 25 isobath with 1 or less points
% 25 and 50 isobath with 2 or less points
% 50 and 75 isobath with 3 or less points
% offshore of 75 isobath with 4 or less points

clear ind cont
cont=1;
te=find(twd>=0 & twd<25);
for j=1:length(te)
   dd=dep{te(j)}; 
   if length(dd)>1
       ind(cont)=te(j);
       cont=cont+1;
   end
end
te=find(twd>=25 & twd<50);
for j=1:length(te)
   dd=dep{te(j)}; 
   if length(dd)>2
       ind(cont)=te(j);
       cont=cont+1;
   end
end
te=find(twd>=50 & twd<75);
for j=1:length(te)
   dd=dep{te(j)}; 
   if length(dd)>3
       ind(cont)=te(j);
       cont=cont+1;
   end
end
te=find(twd>=75);
for j=1:length(te)
   dd=dep{te(j)}; 
   if length(dd)>4
       ind(cont)=te(j);
       cont=cont+1;
   end
end

ind1=sort(ind);
clear ind
ind=ind1;
clear ind1

whos yea mon day lon lat dep watdep temp sal twd

yea=yea(ind);
mon=mon(ind);
day=day(ind);
lon=lon(ind);
lat=lat(ind);
dep=dep(ind);
watdep=watdep(ind);
temp=temp(ind);
sal=sal(ind);
twd=twd(ind);

whos yea mon day lon lat dep watdep temp sal twd

plot(lon,lat,'.g')

%3)
% Now, want to impose that minimum separation between samples over the
% shelf is less than 25 m in the vertical
% shelf is defined as water less than 200 m

clear ind cont
ind=find(twd>200);
cont=length(ind)+1;

te=find(twd<200);
for j=1:length(te)
    dd=dep{te(j)};
    if max(diff(dd))<=25
        ind(cont)=te(j);
        cont=cont+1;
    end
end

ind1=sort(ind);
clear ind
ind=ind1;
clear ind1

whos yea mon day lon lat dep watdep temp sal twd

whos ind

yea=yea(ind);
mon=mon(ind);
day=day(ind);
lon=lon(ind);
lat=lat(ind);
dep=dep(ind);
watdep=watdep(ind);
temp=temp(ind);
sal=sal(ind);
twd=twd(ind);

whos yea mon day lon lat dep watdep temp sal twd

plot(lon,lat,'.k')

%4)
% Now, want to garantee that density inversions between adjacent samples
% are less than 0.05 km m^-3. Apply this criterion to all profiles (shelf,
% slope and offshore)
clear ind cont

cont=1;
for j=1:length(temp)
    if mod(j,1000)==0;display([j length(temp)])
    end
    
    dd=dep{j};
    tt=temp{j};
    ss=sal{j};
    rho=sw_pden(ss,tt,dd,0);
    if min(diff(rho))>-0.05
        ind(cont)=j;
        cont=cont+1;
    end
end

ind1=sort(ind);
clear ind
ind=ind1;
clear ind1

whos ind
whos yea mon day lon lat dep watdep temp sal twd

yea=yea(ind);
mon=mon(ind);
day=day(ind);
lon=lon(ind);
lat=lat(ind);
dep=dep(ind);
watdep=watdep(ind);
temp=temp(ind);
sal=sal(ind);
twd=twd(ind);

whos yea mon day lon lat dep watdep temp sal twd

plot(lon,lat,'.m')

% Now, clean values that are obviously absurd

% First, salinity

for j=1:length(sal)
    ss=sal{j};
    te=find(abs(ss)>10000);
    ss(te)=nan;
    sal{j}=ss;
end


clear ind cont

cont=1;
for j=1:length(sal)
    ss=sal{j};
    te=find(ss<0 | ss > 40);
    te1=find(ss<=37);
    if isempty(te) & ~isempty(te1)
        ind(cont)=j;
        cont=cont+1;
    end
end

whos ind
whos yea mon day lon lat dep watdep temp sal twd

yea=yea(ind);
mon=mon(ind);
day=day(ind);
lon=lon(ind);
lat=lat(ind);
dep=dep(ind);
watdep=watdep(ind);
temp=temp(ind);
sal=sal(ind);
twd=twd(ind);

whos yea mon day lon lat dep watdep temp sal twd

plot(lon,lat,'.k')

% Now, temperature

for j=1:length(temp)
    tt=temp{j};
    te=find(abs(tt)>10000);
    tt(te)=nan;
    temp{j}=tt;
end


clear ind cont
cont=1;
for j=1:length(temp)
    tt=temp{j};
    te=find(tt<-3 | tt > 40);
    if isempty(te) 
        ind(cont)=j;
        cont=cont+1;
    end
end

   
whos yea mon day lon lat dep watdep temp sal twd
whos ind


yea=yea(ind);
mon=mon(ind);
day=day(ind);
lon=lon(ind);
lat=lat(ind);
dep=dep(ind);
watdep=watdep(ind);
temp=temp(ind);
sal=sal(ind);
twd=twd(ind);

whos yea mon day lon lat dep watdep temp sal twd

plot(lon,lat,'.c')


% Now compute statistics of profiles

%shelf, slope and off
te=find(twd<200);
per(1)=length(te)/length(twd)*100;
num(1)=length(te);
te=find(twd>=200 & twd < 1000);
per(2)=length(te)/length(twd)*100;
num(2)=length(te);
te=find(twd>=1000);
per(3)=length(te)/length(twd)*100;
num(3)=length(te);

subplot(2,2,1)
shelfx=[0 1 1 0 0];
shelfy=[0 0 per(1) per(1) 0];
fill(shelfx,shelfy,'k')
hold on
slopex=[1 2 2 1 1];
slopey=[0 0 per(2) per(2) 0];
fill(slopex,slopey,'r')
offx=[2 3 3 2 2];
offy=[0 0 per(3) per(3) 0];
fill(offx,offy,'b')
ylabel('% of valid profiles')

print -djpeg -r300 02cleanNODC_shelf_slope_off

close all

% month of the year separetely for shelf and slope
clear shelf slope mopershelf moperslope
shelf=find(twd<200);
slope=find(twd>=200 & twd < 1000);
for j=1:12
    te=find(mon(shelf)==j);
    mopershelf(j)=length(te)/length(shelf)*100;
    monumshelf(j)=length(te);
    te1=find(mon(slope)==j);
    moperslope(j)=length(te1)/length(slope)*100;
    monumslope(j)=length(te1);
end

clear tex tey
for j=1:12
    tex=[j-0.45 j+0.45 j+0.45 j-0.45 j-0.45];
    teyshelf=[0 0 mopershelf(j) mopershelf(j) 0];
    teyslope=[0 0 moperslope(j) moperslope(j) 0];
    
    subplot(2,1,1)
    fill(tex,teyshelf,'k')
    hold on
    axis([0 13 0 30])
    if j==1;
        text(1,15,'shelf');
        ylabel('% of profiles');
    end
    
    subplot(2,1,2)
    fill(tex,teyslope,'k')
    hold on
    axis([0 13 0 30])
    if j==1;
        text(1,15,'slope');
        ylabel('% of profiles');
        xlabel('Month');
    end
    
end
print -djpeg -r300 02cleanNODC_shelf_slope_by_month_per

close all
clear tex tey
for j=1:12
    tex=[j-0.45 j+0.45 j+0.45 j-0.45 j-0.45];
    teyshelf=[0 0 monumshelf(j) monumshelf(j) 0];
    teyslope=[0 0 monumslope(j) monumslope(j) 0];
    
    subplot(2,1,1)
    fill(tex,teyshelf,'k')
    hold on
    axis([0 13 0 800])
    if j==1;
        text(1,350,'shelf');
        ylabel('number of profiles');
    end
    
    subplot(2,1,2)
    fill(tex,teyslope,'k')
    hold on
    axis([0 13 0 2500])
    if j==1;
        text(1,350,'slope');
        ylabel('number of profiles');
        xlabel('Month');
    end
    
end
print -djpeg -r300 02cleanNODC_shelf_slope_by_month_num

close all



% sorting by date
%month
clear ind
[mon2,ind]=sort(mon);
yea=yea(ind);
mon=mon(ind);
day=day(ind);
lon=lon(ind);
lat=lat(ind);
dep=dep(ind);
watdep=watdep(ind);
temp=temp(ind);
sal=sal(ind);
twd=twd(ind);
whos yea mon day lon lat dep watdep temp sal twd

%day
    clear ind
    te=find(mon==1);
    [day2,ind]=sort(day(te));
        yea1=yea(te(ind));
        mon1=mon(te(ind));
        day1=day(te(ind));
        lon1=lon(te(ind));
        lat1=lat(te(ind));
        dep1=dep(te(ind));
        watdep1=watdep(te(ind));
        temp1=temp(te(ind));
        sal1=sal(te(ind));
        twd1=twd(te(ind));

for j=2:12
    clear ind
    te=find(mon==j);
    [day2,ind]=sort(day(te));
        tu=length(yea1);
        yea1(tu+1:tu+length(te))=yea(te(ind));
        mon1(tu+1:tu+length(te))=mon(te(ind));
        day1(tu+1:tu+length(te))=day(te(ind));
        lon1(tu+1:tu+length(te))=lon(te(ind));
        lat1(tu+1:tu+length(te))=lat(te(ind));
        dep1(tu+1:tu+length(te))=dep(te(ind));
        watdep1(tu+1:tu+length(te))=watdep(te(ind));
        temp1(tu+1:tu+length(te))=temp(te(ind));
        sal1(tu+1:tu+length(te))=sal(te(ind));
        twd1(tu+1:tu+length(te))=twd(te(ind));
end

clear yea mon day lon lat dep watdep temp sal twd
yea=yea1;
mon=mon1;
day=day1;
lon=lon1;
lat=lat1;
dep=dep1;
watdep=watdep1;
temp=temp1;
sal=sal1;
twd=twd1;
clear yea1 mon1 day1 lon1 lat1 dep1 watdep1 temp1 sal1 twd1

shelf=find(twd<200);
slope=find(twd>=200 & twd < 1000);

close all
plot(lon,lat,'.')
hold on
plot(lon(shelf),lat(shelf),'.r')
plot(lon(slope),lat(slope),'.k')
io=axis;
plot(cx,cy,'c')
%contour(XX,YY,ZZ,[ 500 1000 3000],'g')
axis(io)

print -djpeg -r300 02cleanNODC_map


figure
subplot(2,1,1)
hist(yea,1910:10:2020)
axis tight
ylabel('year')

subplot(2,1,2)
hist(mon,1:12)
axis tight
xlabel('month')

print -djpeg -r300 02cleanNODC_histogram_year


whos yea mon day lon lat dep watdep temp sal twd XX YY ZZ cx cy

save 02cleanNODC yea mon day lon lat dep watdep temp sal twd XX YY ZZ cx cy




for j=1:length(temp)
    t0(j)=temp{j}(1);
    s0(j)=sal{j}(1);
end

te=find(isfinite(t0+s0));

TUTU=1
if TUTU==2
xO=-65:0.2:-40;
yO=55:0.2:67;
[XO,YO]=meshgrid(xO,yO);
TOPO_0=griddata(XX,YY,ZZ,XO,YO);
save obj_interp_topo xO yO XO YO TOPO_0
end
load obj_interp_topo

loc=find(isnan(TOPO_0));


%ZOt = barnes(lon(te),lat(te),t0(te),XO,YO,1,1);
%ZOt(loc)=nan;

%ZOs = barnes(lon(te),lat(te),s0(te),XO,YO,1,1);
%ZOs(loc)=nan;


te1=find(isfinite(t0+s0));
ii=[30 35];
io=[-65.00        -40.00         55.00         67.00];

figure
te=find(mon(te1)>=7 & mon(te1)<=8 & yea(te1) >=1950 & yea(te1) <=1960);
subplot(2,3,1)
scatter(lon(te1(te)),lat(te1(te)),5,s0(te1(te)),'filled')
caxis(ii)
hold on
plot(cx,cy,'k')
axis(io)

te=find(mon(te1)>=7 & mon(te1)<=8 & yea(te1) >=1960 & yea(te1) <=1970);
subplot(2,3,2)
scatter(lon(te1(te)),lat(te1(te)),5,s0(te1(te)),'filled')
caxis(ii)
hold on
plot(cx,cy,'k')
axis(io)

te=find(mon(te1)>=7 & mon(te1)<=8 & yea(te1) >=1970 & yea(te1) <=1980);
subplot(2,3,3)
scatter(lon(te1(te)),lat(te1(te)),5,s0(te1(te)),'filled')
caxis(ii)
hold on
plot(cx,cy,'k')
axis(io)

te=find(mon(te1)>=7 & mon(te1)<=8 & yea(te1) >=1980 & yea(te1) <=1990);
subplot(2,3,4)
scatter(lon(te1(te)),lat(te1(te)),5,s0(te1(te)),'filled')
caxis(ii)
hold on
plot(cx,cy,'k')
axis(io)

te=find(mon(te1)>=7 & mon(te1)<=8 & yea(te1) >=1990 & yea(te1) <=2000);
subplot(2,3,5)
scatter(lon(te1(te)),lat(te1(te)),5,s0(te1(te)),'filled')
caxis(ii)
hold on
plot(cx,cy,'k')
axis(io)

te=find(mon(te1)>=7 & mon(te1)<=8 & yea(te1) >=2000 & yea(te1) <=2013);
subplot(2,3,6)
scatter(lon(te1(te)),lat(te1(te)),5,s0(te1(te)),'filled')
caxis(ii)
hold on
plot(cx,cy,'k')
axis(io)




%ZOs = barnes(lon(te),lat(te),s0(te),XO,YO,1,1);
%ZOs(loc)=nan;

te1=find(isfinite(t0+s0));
ii=[30 35];
io=[-65.00        -40.00         55.00         67.00];
inc=[2007:2013];
%inc=[2001:2007];

close all
for j=1:6
    
te=find(mon(te1)>=7 & mon(te1)<=9 & yea(te1) >=inc(j) & yea(te1) <inc(j+1));
ZOs = barnes(lon(te1(te)),lat(te1(te)),s0(te1(te)),XO,YO,1,1);
ZOs(loc)=nan;
subplot(2,3,j)
pcolor(XO,YO,ZOs)
shading flat
caxis(ii)
hold on
plot(cx,cy,'k')
axis(io)
%scatter(lon(te1(te)),lat(te1(te)),5,s0(te1(te)),'filled')
plot(lon(te1(te)),lat(te1(te)),'.','color',[1 1 1]*0.7,'markersize',5)
title([num2str(inc(j)) ', 7-9'])
end
po=get(gca,'position');
cb=colorbar;
po1=get(gca,'position');
set(gca,'position',po)


print -djpeg -r300 Z:\IMCS\NODC\greenland\surf_salinity_per_year2
%print -djpeg -r300 Z:\IMCS\NODC\greenland\surf_salinity_per_year1



te1=find(isfinite(t0+s0));
ii=[30 35];
io=[-65.00        -40.00         55.00         67.00];
inc=[1950:10:2000 2013];

close all
for j=1:6
    
te=find(mon(te1)>=7 & mon(te1)<=9 & yea(te1) >=inc(j) & yea(te1) <=inc(j+1));
ZOs = barnes(lon(te1(te)),lat(te1(te)),s0(te1(te)),XO,YO,1,1);
ZOs(loc)=nan;
subplot(2,3,j)
pcolor(XO,YO,ZOs)
shading flat
caxis(ii)
hold on
plot(cx,cy,'k')
axis(io)
%scatter(lon(te1(te)),lat(te1(te)),5,s0(te1(te)),'filled')
plot(lon(te1(te)),lat(te1(te)),'.','color',[1 1 1]*0.7,'markersize',5)
title([num2str(inc(j)) ' to ' num2str(inc(j+1)) ', 7-9'])
end
po=get(gca,'position');
cb=colorbar;
po1=get(gca,'position');
set(gca,'position',po)


print -djpeg -r300 Z:\IMCS\NODC\greenland\surf_salinity_per_decade

