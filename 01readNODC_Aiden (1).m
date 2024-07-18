clear all
close all



cd C:\UGA\IMCS\NODC\greenland\profiles

direc=ls;
direc(1:2,:)=[];

BB=1;
for DD=1:size(direc,1)
        
    ld=['cd ' direc(DD,:)];
    eval(ld)


name=ls
name(1:2,:)=[];

for j=1:size(name,1)
    if mod(j,50)==0
        display([num2str(DD) '   ' num2str(size(direc,1)) '    '  num2str(j) '   ' num2str(size(name,1))])
    end
    
 nh=[name(j,:)];

 lat(BB)=ncread(nh,'lat');
 lon(BB)=ncread(nh,'lon');
 date=double(ncread(nh,'date'));
 temp{BB}=ncread(nh,'Temperature');
 sal{BB}=ncread(nh,'Salinity');
 depth{BB}=ncread(nh,'z');
 uba=num2str(date);
 year(BB)=str2num(uba(1:4));
 month(BB)=str2num(uba(5:6));
 day(BB)=str2num(uba(7:8));
 

BB=BB+1;
end

cd ..

end
cd ..

whos lat lon temp sal depth year month day

save 01Greenland_temporary lat lon temp sal depth year month day



% checking if there are profiles with depths 'inverted'. Also,
% eliminate profiles with less than two points
cont=1;
for j=1:size(depth,2)
    dt=depth{j};
    iit=min(find(diff(dt)<0));
    if isempty(iit)  & length(dt)>1 
        ind(cont)=j;
        cont=cont+1;
    end
end

whos year month day lon lat depth temp sal ind

year=year(ind);
month=month(ind);
day=day(ind);
lon=lon(ind);
lat=lat(ind);
depth=depth(ind);
temp=temp(ind);
sal=sal(ind);

dept=depth;
sals=sal;
temt=temp;

clear depth sal temp

whos year month day lon lat dept temt sals ind
clear ind cont ds dt iis iit j ss tt


% now, need to select only profiles where both T and S exist
clear ind
cont=1;
for j=1:size(sals,2)
    if mod(j,1)==1000;display([j size(sals,2)])
    end
t=temt{j};
s=sals{j};

    if sum(isfinite(t))==sum(isfinite(s))
        ind(cont)=j;
        cont=cont+1;
    end
end

whos year month day lon lat dept temt sals ind

year=year(ind);
month=month(ind);
day=day(ind);
lon=lon(ind);
lat=lat(ind);
dept=dept(ind);
temt=temt(ind);
sals=sals(ind);


whos year month day lon lat dept temt sals ind
clear ind cont ds dt iis iit j ss tt



for j=1:size(dept,2)
    watdept(j)=max(dept{j});
end





ind=find(isfinite(watdept));
year=year(ind);
month=month(ind);
day=day(ind);
lon=lon(ind);
lat=lat(ind);
dept=dept(ind);
temt=temt(ind);
sals=sals(ind);
watdept=watdept(ind);



whos  year month day lon lat dept temt sals ind watdept 
clear ind cont



watdep=watdept;
dep=dept;
temp=temt;
sal=sals;

yea=year;
mon=month;

whos yea mon day lon lat dep watdep temp sal
save 01GREENLAND_temp_sal yea mon day lon lat dep watdep temp sal

