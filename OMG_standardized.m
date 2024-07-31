load('all_OMG_CTD_from_csv_ver_1.3.mat')
% normalize OMG Data
for i = 1:1:1908
Data{1,i}= (struct2cell(cast{1,i})) ; 
end
for j = 1:1:7
for i = 1:1:1908
OMG_data{j,i} = Data{1,i}{j,1} ;
end
end
%adding human time
format = ('mmmm,dd,yyyy') ;
for i= 1:1:1908
DataTime{1,i} = datevec(OMG_data{2,i},format) ;
end
for i= 1:1:length(Data)
    OMG_lon(1,i) = OMG_data{3,i}(1,1) ;
    OMG_lat(1,i) = OMG_data{4,i}(1,1) ;
    OMG_yea(1,i) = DataTime{1,i}(1,1) ;
    OMG_mon(1,i) = DataTime{1,i}(1,2) ;
    OMG_day(1,i) = DataTime{1,i}(1,3) ;
    OMG_temp{i} = OMG_data{5,i} ;
    OMG_sal{i} = OMG_data{6,i} ;
    OMG_depth{i} = OMG_data{7,i} ;
end
% cut off at -30
idx = (OMG_lon <= -30) ;
OMG_lon = OMG_lon(idx) ;
OMG_lat = OMG_lat(idx) ;
OMG_yea = OMG_yea(idx) ;
OMG_mon = OMG_mon(idx) ;
OMG_day = OMG_day(idx) ;
OMG_temp = OMG_temp(idx) ;
OMG_sal = OMG_sal(idx) ;
OMG_depth = OMG_depth(idx) ;
save OMG_data.mat OMG_lon OMG_depth OMG_sal OMG_temp OMG_day OMG_mon OMG_yea OMG_lat 