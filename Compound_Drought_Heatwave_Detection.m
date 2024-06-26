%%Analyzing Compound Drought and Heatwaves


%%Deriving the weekly drought indices
lon = ncread("/Users/smonda15/Ongoing Work/Publication_Compound_Drought_Heatwave_Over_Phoneix/Datasets/CPC Global Precipitation/precip.1979.nc",'lon');
lat = ncread("/Users/smonda15/Ongoing Work/Publication_Compound_Drought_Heatwave_Over_Phoneix/Datasets/CPC Global Precipitation/precip.1979.nc",'lat');

folderPath = '/Users/smonda15/Ongoing Work/Publication_Compound_Drought_Heatwave_Over_Phoneix/Datasets/CPC Global Precipitation';
ncFiles = dir(fullfile(folderPath, '*.nc'));


b=zeros(259200,1);
d=zeros(259200,1);

c=1979;

%for i = 1:length(ncFiles)
for i=1:45
    
    filePath = fullfile(folderPath, ncFiles(i).name);
    
    
    precipData = ncread(filePath, 'precip');
    a       = size(precipData);
    
  
    precipData = reshape(precipData,[a(1)*a(2),a(3)]);
    isLeapYear = (rem(c, 4) == 0 && rem(c, 100) ~= 0) || (rem(c, 400) == 0);
    if isLeapYear
        precipData1 = precipData(:,92:274);
        precipData2 = precipData(:,1:364);
    else 
        precipData1 = precipData(:,91:273);
        precipData2 = precipData(:,1:364);
    end
    b  =cat(2,b,precipData1);
    d  =cat(2,d,precipData2);
    

    c=c+1;
    i
end

b(:,1)=[];
d(:,1)=[];


SWNA             =  readmatrix('/Users/smonda15/Ongoing Work/Publication_Compound_Drought_Heatwave_Over_Phoneix/Datasets/Shapefiles/Southwest_north_America_Grid_Locations.csv');
SWNA_latlon      =  SWNA(:,1:2);
SWNA_index       =  SWNA(:,3);


SW_summer_pcp  =  b(SWNA_index,:);
SW_precip      =  d(SWNA_index,:);




%Converting Daily Precipitation data to weekly scale (considering the whole
%year)
numLocations = 2256;
numDays = 16380;
numYears = 2023 - 1979 + 1;

% Reshape the data to have one year per layer
numDaysPerYear = 364;
dataMatrixReshaped = reshape(SW_precip, numLocations, numDaysPerYear, numYears);

% Determine the number of weeks (assuming 12 weeks per summer, adjust as necessary)
numWeeksPerYear = 52;

% Initialize the weekly data matrix
weeklyDataMatrix = zeros(numLocations, numWeeksPerYear, numYears);

% Loop over each year and each week to calculate weekly averages
for year = 1:numYears
    for week = 1:numWeeksPerYear
        % Define the range of days for this week
        startDay = (week - 1) * 7 + 1;
        endDay = week * 7;
        if endDay > numDaysPerYear
            endDay = numDaysPerYear;
        end

        % Calculate the weekly average for each location
        weeklyDataMatrix(:, week, year) = sum(dataMatrixReshaped(:, startDay:endDay, year), 2);
    end
end

% Reshape the weekly data matrix back to two dimensions if needed
weeklyDataMatrix_r = reshape(weeklyDataMatrix, numLocations, numWeeksPerYear * numYears);



weeklyDataMatrix_r(find(weeklyDataMatrix_r<0))=0;


SW_SPI12 = zeros(2288,2256);
%SW_SPI2 = zeros(2288,2256);
%SW_SPI2 = zeros(585,439);
for i=1:2256
    SW_SPI12(:,i)= SPI(weeklyDataMatrix_r(i,:)',48,52);
    %SW_SPI2(:,i)= SPI(weeklyDataMatrix(i,:)',12,52);
end

SW_SPI3 = zeros(2288,2256);
%SW_SPI2 = zeros(2288,2256);
%SW_SPI2 = zeros(585,439);
for i=1:2256
    SW_SPI3(:,i)= SPI(weeklyDataMatrix_r(i,:)',12,52);
    %SW_SPI2(:,i)= SPI(weeklyDataMatrix(i,:)',12,52);
end

%Conisdering Moderate Drought (Threshold = -1.0 for moderate, -0.5 for
%mild, -1.5 for severe and -2.0 for extreme)

SW_drought3 = zeros(2288,2256);

for i=1:2288
    for j=1:2256
        if SW_SPI3(i,j)<-0.5
            SW_drought3(i,j)=1;
        end
    end
end

SW_drought12 = zeros(2288,2256);

for i=1:2288
    for j=1:2256
        if SW_SPI12(i,j)<-0.5
            SW_drought12(i,j)=1;
        end
    end
end



%Deriving The Heatwave Characteristics

lon = ncread("/Users/smonda15/Ongoing Work/Publication_Compound_Drought_Heatwave_Over_Phoneix/Datasets/CPC Global Temperature/tmax.1979.nc",'lon');
lat = ncread("/Users/smonda15/Ongoing Work/Publication_Compound_Drought_Heatwave_Over_Phoneix/Datasets/CPC Global Temperature/tmax.1979.nc",'lat');

folderPath = '/Users/smonda15/Ongoing Work/Publication_Compound_Drought_Heatwave_Over_Phoneix/Datasets/CPC Global Temperature';
ncFiles = dir(fullfile(folderPath, '*.nc'));


b=zeros(259200,1);
d=zeros(259200,1);
c=1979;

%for i = 1:length(ncFiles)
  
for i = 1:45
    filePath = fullfile(folderPath, ncFiles(i).name);
    
    
    maxtempData = ncread(filePath, 'tmax');
    a       = size(maxtempData);
    
  
    maxtempData = reshape(maxtempData,[a(1)*a(2),a(3)]);
    isLeapYear = (rem(c, 4) == 0 && rem(c, 100) ~= 0) || (rem(c, 400) == 0);
    if isLeapYear
        maxtempData1 = maxtempData(:,92:274);
        maxtempData2 = maxtempData(:,90:276);
    else 
        maxtempData1 = maxtempData(:,91:273);
        maxtempData2 = maxtempData(:,89:275);
    end
    b  =cat(2,b,maxtempData1);
    d  =cat(2,d,maxtempData2);

    c=c+1;
    i
end

b(:,1)  =   [];
d(:,1)  =   [];

V=zeros(259200,2);
for i=1:360
    V(((720*(i-1)+1):(720*i)),1)=lon(1:720,1);
    V(((720*(i-1)+1):(720*i)),2)=lat(i,1);
end


%Changing the coordinate system so that the southwest shapefile and CPC
%coordinates can be projected together
V_new = zeros(259200,2);
for i=1:259200
    if V(i,1)>180
        V_new(i,1)=-(360-V(i,1));
        V_new(i,2)= V(i,2);
    else
        V_new(i,1)= V(i,1);
        V_new(i,2)= V(i,2);
    end
end



V_new = zeros(259200,2);
for i=1:259200
    if V(i,1)>180
        V_new(i,1)=-(360-V(i,1));
        V_new(i,2)= V(i,2);
    else
        V_new(i,1)= V(i,1);
        V_new(i,2)= V(i,2);
    end
end


%Extracting Grid-Locations Corresponding to Southwest USA

SWNA             =  readmatrix('/Users/smonda15/Ongoing Work/Publication_Compound_Drought_Heatwave_Over_Phoneix/Datasets/Shapefiles/Southwest_north_America_Grid_Locations.csv');
SWNA_latlon      =  SWNA(:,1:2);
SWNA_index       =  SWNA(:,3);



%Extracting Temperature Data Corresponding to Southwest USA
SWNAheat         =  b(SWNA_index ,:);
SWNAheat1        =  d(SWNA_index ,:);

%Treating the missing values
%Missing data present at time step 599, 623 and 1257. Taking average of
%598 and 600 th timestep for filling in 599th time step and doing the same
%for others

SWNAheat(:,758)   =  SWNAheat(:,757);
SWNAheat(:,762)   =  SWNAheat(:,763);
SWNAheat(:,760)   =  0.5*(SWNAheat(:,758)+SWNAheat(:,762));
SWNAheat(:,759)   =  0.5*(SWNAheat(:,758)+SWNAheat(:,760));
SWNAheat(:,761)   =  0.5*(SWNAheat(:,760)+SWNAheat(:,762));
SWNAheat(:,1206)  =  0.5*(SWNAheat(:,1205)+SWNAheat(:,1207));
SWNAheat(:,1230)  =  0.5*(SWNAheat(:,1229)+SWNAheat(:,1231));
SWNAheat(:,1454)  =  0.5*(SWNAheat(:,1453)+SWNAheat(:,1455));
SWNAheat(:,2501)  =  0.5*(SWNAheat(:,2500)+SWNAheat(:,2502));

SWNAheat1(:,776)   =  SWNAheat1(:,775);
SWNAheat1(:,780)   =  SWNAheat1(:,781);
SWNAheat1(:,778)   =  0.5*(SWNAheat1(:,780)+SWNAheat1(:,776));
SWNAheat1(:,777)   =  0.5*(SWNAheat1(:,776)+SWNAheat1(:,778));
SWNAheat1(:,779)   =  0.5*(SWNAheat1(:,778)+SWNAheat1(:,780));
SWNAheat1(:,1232)  =  0.5*(SWNAheat1(:,1231)+SWNAheat1(:,1233));
SWNAheat1(:,1256)  =  0.5*(SWNAheat1(:,1255)+SWNAheat1(:,1257));
SWNAheat1(:,1484)  =  0.5*(SWNAheat1(:,1483)+SWNAheat1(:,1485));
SWNAheat1(:,2555)  =  0.5*(SWNAheat1(:,2554)+SWNAheat1(:,2556));


%Calculating spatially averaged and temporally averaged anomaly

SWNAheat_r   = reshape(SWNAheat,[2256,183,45]);

%Detection of Calendar Percentile and Heatwave Detection

SWNAheat_pool      =   zeros(2256,183,225);
for i=1:2256
    for j=1:183
        t=zeros(2,1);

        for k=1:45
            s  = (SWNAheat1(i,(183*(k-1)+j):(183*(k-1)+j+4)))';
            t  = (cat(1,t,s));
        end
      t(1:2)=[];  
     SWNAheat_pool(i,j,:)=t';
    end
    i
end

%Calculating the 90th and 95th percentile and other percentiles

SWheat_95 = zeros(2256,183);
SWheat_90 = zeros(2256,183);
SWheat_80 = zeros(2256,183);
SWheat_75 = zeros(2256,183);
SWheat_25 = zeros(2256,183);


for i=1:2256
    for j=1:183
        SWheat_95(i,j)=prctile(SWNAheat_pool(i,j,:),95);
        SWheat_90(i,j)=prctile(SWNAheat_pool(i,j,:),90);
        SWheat_80(i,j)=prctile(SWNAheat_pool(i,j,:),80);
        SWheat_25(i,j)=prctile(SWNAheat_pool(i,j,:),25);
        SWheat_75(i,j)=prctile(SWNAheat_pool(i,j,:),75);
    end
end

%Calculating Spatiotemporal Data with heatwave occurrence as 1 and not
%occurrence as 0

SWheat_r    =  reshape(SWNAheat,[2256,183,45]);
SWheat_b    =  zeros(2256,183,45);
SWheat_d    =  zeros(2256,183,45);
SWheat_p    =  zeros(2256,183,45);

for i=1:2256
    for j=1:183
        for k=1:45
            if SWheat_r(i,j,k)>SWheat_95(i,j)
                SWheat_b(i,j,k)=1;
            else 
                SWheat_b(i,j,k)=0;
            end
            if SWheat_r(i,j,k)>SWheat_90(i,j)
                SWheat_d(i,j,k)=1;
            else 
                SWheat_d(i,j,k)=0;
            end
            if SWheat_r(i,j,k)>SWheat_80(i,j)
                SWheat_p(i,j,k)=1;
            else 
                SWheat_p(i,j,k)=0;
            end
        end
    end
end

SWheat_d1 = (reshape(SWheat_d,[2256,183*45]))';
SWheat_b1 = (reshape(SWheat_b,[2256,183*45]))';
SWheat_p1 = (reshape(SWheat_p,[2256,183*45]))';
%Calculating Occurrence of Heatwaves as Binary Time series        
heatwaveEvents_95 = findHeatwaves(SWheat_b);


heatwaveEvents_90 = findHeatwaves(SWheat_d);

heatwaveEvents_80 = findHeatwaves(SWheat_p);


heatwaveEvents_95 = (reshape(heatwaveEvents_95,[2256,183*45]))';
heatwaveEvents_90 = (reshape(heatwaveEvents_90,[2256,183*45]))';
heatwaveEvents_80 = (reshape(heatwaveEvents_80,[2256,183*45]))';

hw_temp_95        = (heatwaveEvents_95.*SWNAheat'); 
hw_temp_90        = (heatwaveEvents_90.*SWNAheat');
hw_temp_80        = (heatwaveEvents_80.*SWNAheat');

%hw_temp_95 & hw_temp_90 can be used to calculate the heatwave severity as
%well as the 

%Calculating Duration of each Heatwaves

heatwaveDuration_95       = DurationHeatwaves(SWheat_b);

heatwaveDuration_90       = DurationHeatwaves(SWheat_d);
heatwaveDuration_80       = DurationHeatwaves(SWheat_p);

heatwaveDuration_95       = (reshape(heatwaveDuration_95,[2256,183*45]))';
heatwaveDuration_90       = (reshape(heatwaveDuration_90,[2256,183*45]))';
heatwaveDuration_80       = (reshape(heatwaveDuration_80,[2256,183*45]))';

%Compound Drought and Heatwave
heatwaveEvents_80_r = zeros(2256,183,45);

for i=1:45
    m = 183*(i-1)+1;
    n = 183*i;

    heatwaveEvents_80_r(:,:,i)= (heatwaveEvents_80(m:n,:))';
end

heatwaveEvents_80_r = heatwaveEvents_80_r(:,1:182,:);

heatwaveEvents_80_r = reshape(heatwaveEvents_80_r,[2256,182*45]);

heatwaveEvents_80_rw = zeros(2256,1170);

for i=1:1170
    k = 7*(i-1)+1;
    l = 7*i;
    for j=1:2256
    heatwaveEvents_80_rw(j,i) = sum(heatwaveEvents_80_r(j,k:l),2);
    end
end

%heatwaveEvents_90_rwsum  =  sum(heatwaveEvents_90_rw,3);


%Extracting data from 1980 to match with drought

heatwaveEvents_80_rw_1980      = heatwaveEvents_80_rw(:,27:1170);



%Calculating the heatwave severity

SWheat_90_r = SWheat_90';
SWheat_95_r = SWheat_95';
SWheat_80_r = SWheat_80';

SWheat_severity_90 = zeros(8235,2256);
SWheat_severity_95 = zeros(8235,2256);
SWheat_severity_80 = zeros(8235,2256);

for i=1:45
    for j=1:2256
        k = 183*(i-1)+1;
        l = 183*i;
        SWheat_severity_90(k:l,j)=SWNAheat(j,k:l)'-(SWheat_90_r(:,j));
        SWheat_severity_95(k:l,j)=SWNAheat(j,k:l)'-(SWheat_95_r(:,j));
        SWheat_severity_80(k:l,j)=SWNAheat(j,k:l)'-(SWheat_80_r(:,j));
    end
end

SWheat_hwseverity_90 = zeros(8235,2256);

SWheat_hwseverity_95 = zeros(8235,2256);
SWheat_hwseverity_80 = zeros(8235,2256);

for i=1:2256
    a = find(heatwaveDuration_90(:,i)>0);
    if a~=0
        b = heatwaveDuration_90(a,i);
        for j=1:size(a,1)

            
        SWheat_hwseverity_90(a(j,1),i) = sum(SWheat_severity_90(a(j,1):(a(j,1)+b(j,1)-1),i));
        end
    end
end


for i=1:2256
    a = find(heatwaveDuration_95(:,i)>0);
    if a~=0
        b = heatwaveDuration_95(a,i);
        for j=1:size(a,1)
        SWheat_hwseverity_95(a(j,1),i) = sum(SWheat_severity_95(a(j,1):(a(j,1)+b(j,1)-1),i));
        end
    end
end


for i=1:2256
    a = find(heatwaveDuration_80(:,i)>0);
    if a~=0
        b = heatwaveDuration_80(a,i);
        for j=1:size(a,1)
        SWheat_hwseverity_80(a(j,1),i) = sum(SWheat_severity_80(a(j,1):(a(j,1)+b(j,1)-1),i));
        end
    end
end

%Calculating Average Heatwave Severity for  whole summer


SWNAheat_avg_sev_80 = zeros(2256,1);
SWNAheat_avg_sev_90 = zeros(2256,1);
SWNAheat_avg_sev_95 = zeros(2256,1);

for i=1:2256
    a = find(heatwaveDuration_80(:,i)>0);
    if a~=0
        b = SWheat_hwseverity_80(a,i);
        SWNAheat_avg_sev_80(i,1) = mean(b);
    else 
        SWNAheat_avg_sev_80(i,1) = 0;
    end
end

for i=1:2256
    a = find(heatwaveDuration_90(:,i)>0);
    if a~=0
        b = SWheat_hwseverity_90(a,i);
        SWNAheat_avg_sev_90(i,1) = mean(b);
    else 
        SWNAheat_avg_sev_90(i,1) = 0;
    end
end

for i=1:2256
    a = find(heatwaveDuration_95(:,i)>0);
    if a~=0
        b = SWheat_hwseverity_95(a,i);
        SWNAheat_avg_sev_95(i,1) = mean(b);
    else 
        SWNAheat_avg_sev_95(i,1) = 0;
    end
end

%Calculating average severity for 2023

SWNAheat_avg_sev_80_2023 = zeros(2256,1);
SWNAheat_avg_sev_90_2023 = zeros(2256,1);
SWNAheat_avg_sev_95_2023 = zeros(2256,1);

for i=1:2256
    a = find(heatwaveDuration_80(8053:8235,i)>0);
    if a~=0
        b = SWheat_hwseverity_80(a+8052,i);
        SWNAheat_avg_sev_80_2023(i,1) = mean(b);
    else 
        SWNAheat_avg_sev_80_2023(i,1) = 0;
    end
end

for i=1:2256
    a = find(heatwaveDuration_90(8053:8235,i)>0);
    if a~=0
        b = SWheat_hwseverity_90(a+8052,i);
        SWNAheat_avg_sev_90_2023(i,1) = mean(b);
    else 
        SWNAheat_avg_sev_90_2023(i,1) = 0;
    end
end

for i=1:2256
    a = find(heatwaveDuration_95(8053:8235,i)>0);
    if a~=0
        b = SWheat_hwseverity_95(a+8052,i);
        SWNAheat_avg_sev_95_2023(i,1) = mean(b);
    else 
        SWNAheat_avg_sev_95_2023(i,1) = 0;
    end
end


%Calculating the decadal change characteristics of heatwaves

%First Decade
%Calculating Frequency of Heatwaves
heatwave_Frequency_951 = zeros(2256,1);
heatwave_Frequency_901 = zeros(2256,1);
heatwave_Frequency_801 = zeros(2256,1);

for i=1:2256
    a=find(heatwaveDuration_90(1:2013,i)>0);
    a=size(a);
    b=find(heatwaveDuration_95(1:2013,i)>0);
    b=size(b);
    c=find(heatwaveDuration_80(1:2013,i)>0);
    c=size(c);
    heatwave_Frequency_901(i,1)=a(1,1);
    heatwave_Frequency_951(i,1)=b(1,1);
    heatwave_Frequency_801(i,1)=c(1,1);
end

heatwave_Frequency_952 = zeros(2256,1);
heatwave_Frequency_902 = zeros(2256,1);
heatwave_Frequency_802 = zeros(2256,1);
for i=1:2256
    a=find(heatwaveDuration_90(2014:4026,i)>0);
    a=size(a);
    b=find(heatwaveDuration_95(2014:4026,i)>0);
    b=size(b);
    c=find(heatwaveDuration_80(2014:4026,i)>0);
    c=size(c);
    heatwave_Frequency_902(i,1)=a(1,1);
    heatwave_Frequency_952(i,1)=b(1,1);
    heatwave_Frequency_802(i,1)=c(1,1);
end

heatwave_Frequency_953 = zeros(2256,1);
heatwave_Frequency_903 = zeros(2256,1);
heatwave_Frequency_803 = zeros(2256,1);

for i=1:2256
    a=find(heatwaveDuration_90(4027:6039,i)>0);
    a=size(a);
    b=find(heatwaveDuration_95(4027:6039,i)>0);
    b=size(b);
    c=find(heatwaveDuration_80(4027:6039,i)>0);
    c=size(c);

    heatwave_Frequency_903(i,1)=a(1,1);
    heatwave_Frequency_953(i,1)=b(1,1);
    heatwave_Frequency_803(i,1)=c(1,1);
end

heatwave_Frequency_954 = zeros(2256,1);
heatwave_Frequency_904 = zeros(2256,1);
heatwave_Frequency_804 = zeros(2256,1);

for i=1:2256
    a=find(heatwaveDuration_90(6040:8235,i)>0);
    a=size(a);
    b=find(heatwaveDuration_95(6040:8235,i)>0);
    b=size(b);
    c=find(heatwaveDuration_80(6040:8235,i)>0);
    c=size(c);
    heatwave_Frequency_904(i,1)=a(1,1);
    heatwave_Frequency_954(i,1)=b(1,1);
    heatwave_Frequency_804(i,1)=c(1,1);
end

%isnan(heatwave_Frequency_90);



%Calculating Decadal Change in Duration of Heatwaves

HW_avg_duration_901 =zeros(2256,1);
HW_avg_duration_951 =zeros(2256,1);
HW_avg_duration_801 =zeros(2256,1);

for i=1:2256
    a=find(heatwaveDuration_90(1:2013,i)>0);
    a=heatwaveDuration_90(a,i);
    a=mean(a,1);
    b=find(heatwaveDuration_95(1:2013,i)>0);
    b=heatwaveDuration_95(b,i);
    b=mean(b,1);
    c=find(heatwaveDuration_80(1:2013,i)>0);
    c=heatwaveDuration_80(c,i);
    c=mean(c,1);
    HW_avg_duration_901(i,1)=a;
    HW_avg_duration_951(i,1)=b;
    HW_avg_duration_801(i,1)=c;
end

HW_avg_duration_902 =zeros(2256,1);
HW_avg_duration_952 =zeros(2256,1);
HW_avg_duration_802 =zeros(2256,1);

for i=1:2256
    a=find(heatwaveDuration_90(2014:4026,i)>0);
    a=heatwaveDuration_90(a+2013,i);
    a=mean(a,1);
    b=find(heatwaveDuration_95(2014:4026,i)>0);
    b=heatwaveDuration_95(b+2013,i);
    b=mean(b,1);
    c=find(heatwaveDuration_80(2014:4026,i)>0);
    c=heatwaveDuration_80(c+2013,i);
    c=mean(c,1);

    HW_avg_duration_902(i,1)=a;
    HW_avg_duration_952(i,1)=b;
    HW_avg_duration_802(i,1)=c;
end

HW_avg_duration_903 =zeros(2256,1);
HW_avg_duration_953 =zeros(2256,1);
HW_avg_duration_803 =zeros(2256,1);

for i=1:2256
    a=find(heatwaveDuration_90(4027:6039,i)>0);
    a=heatwaveDuration_90(a+4026,i);
    a=mean(a,1);
    b=find(heatwaveDuration_95(4027:6039,i)>0);
    b=heatwaveDuration_95(b+4026,i);
    b=mean(b,1);
    c=find(heatwaveDuration_80(4027:6039,i)>0);
    c=heatwaveDuration_80(c+4026,i);
    c=mean(c,1);
    HW_avg_duration_903(i,1)=a;
    HW_avg_duration_953(i,1)=b;
    HW_avg_duration_803(i,1)=c;
end


HW_avg_duration_904 =zeros(2256,1);
HW_avg_duration_954 =zeros(2256,1);
HW_avg_duration_804 =zeros(2256,1);
for i=1:2256
    a=find(heatwaveDuration_90(6040:8235,i)>0);
    a=heatwaveDuration_90(a+6039,i);
    a=mean(a,1);
    b=find(heatwaveDuration_95(6040:8235,i)>0);
    b=heatwaveDuration_95(b+6039,i);
    b=mean(b,1);
    c=find(heatwaveDuration_80(6040:8235,i)>0);
    c=heatwaveDuration_80(c+6039,i);
    c=mean(c,1);
    HW_avg_duration_904(i,1)=a;
    HW_avg_duration_954(i,1)=b;
    HW_avg_duration_804(i,1)=c;
end



%Calculating decadal change in the severity of heatwaves

HW_avg_severity_901 =zeros(2256,1);
HW_avg_severity_951 =zeros(2256,1);
HW_avg_severity_801 =zeros(2256,1);

for i=1:2256
    a=find(heatwaveDuration_90(1:2013,i)>0);
    a=SWheat_hwseverity_90(a,i);
    a=mean(a,1);
    b=find(heatwaveDuration_95(1:2013,i)>0);
    b=SWheat_hwseverity_95(b,i);
    b=mean(b,1);
    c=find(heatwaveDuration_80(1:2013,i)>0);
    c=SWheat_hwseverity_80(c,i);
    c=mean(c,1);
    HW_avg_severity_901(i,1)=a;
    HW_avg_severity_951(i,1)=b;
    HW_avg_severity_801(i,1)=c;
end

HW_avg_severity_902 =zeros(2256,1);
HW_avg_severity_952 =zeros(2256,1);
HW_avg_severity_802 =zeros(2256,1);

for i=1:2256
    a=find(heatwaveDuration_90(2014:4026,i)>0);
    a=SWheat_hwseverity_90(a+2013,i);
    a=mean(a,1);
    b=find(heatwaveDuration_95(2014:4026,i)>0);
    b=SWheat_hwseverity_95(b+2013,i);
    b=mean(b,1);
    c=find(heatwaveDuration_80(2014:4026,i)>0);
    c=SWheat_hwseverity_80(c+2013,i);
    c=mean(c,1);

    HW_avg_severity_902(i,1)=a;
    HW_avg_severity_952(i,1)=b;
    HW_avg_severity_802(i,1)=c;
end

HW_avg_severity_903 =zeros(2256,1);
HW_avg_severity_953 =zeros(2256,1);
HW_avg_severity_803 =zeros(2256,1);

for i=1:2256
    a=find(heatwaveDuration_90(4027:6039,i)>0);
    a=SWheat_hwseverity_90(a+4026,i);
    a=mean(a,1);
    b=find(heatwaveDuration_95(4027:6039,i)>0);
    b=SWheat_hwseverity_95(b+4026,i);
    b=mean(b,1);
    c=find(heatwaveDuration_80(4027:6039,i)>0);
    c=SWheat_hwseverity_80(c+4026,i);
    c=mean(c,1);
    HW_avg_severity_903(i,1)=a;
    HW_avg_severity_953(i,1)=b;
    HW_avg_severity_803(i,1)=c;
end


HW_avg_severity_904 =zeros(2256,1);
HW_avg_severity_954 =zeros(2256,1);
HW_avg_severity_804 =zeros(2256,1);
for i=1:2256
    a=find(heatwaveDuration_90(6040:8235,i)>0);
    a=SWheat_hwseverity_90(a+6039,i);
    a=mean(a,1);
    b=find(heatwaveDuration_95(6040:8235,i)>0);
    b=SWheat_hwseverity_95(b+6039,i);
    b=mean(b,1);
    c=find(heatwaveDuration_80(6040:8235,i)>0);
    c=SWheat_hwseverity_80(c+6039,i);
    c=mean(c,1);
    HW_avg_severity_904(i,1)=a;
    HW_avg_severity_954(i,1)=b;
    HW_avg_severity_804(i,1)=c;
end

%Changing Drought data to match the same dimension as the heatwaves
%SW_drought3_r       = SW_drought3';
%SW_drought3_r       = reshape(SW_drought3_r,[2256,52,44]);
%SW_drought3_r       = SW_drought3_r(:,14:39,:);
%SW_drought3_r       = reshape(SW_drought3_r,[2256,26*44]);

%If we are considering long-term drought in Compound Drought & Heatwave
SW_drought12_r        = SW_drought12';

SW_drought12_summer   = zeros(2256,1144);

for i=1:44
    k = 52*(i-1)+1;
    l = 52*i;
    m = 26*(i-1)+1;
    n = 26*i;

    p = SW_drought12_r(:,k:l);

    SW_drought12_summer(:,m:n) = p(:,14:39);
end

%If we are considering short-term drought in compound drought and heatwave
SW_drought3_r        = SW_drought3';

SW_drought3_summer   = zeros(2256,1144);

for i=1:44
    k = 52*(i-1)+1;
    l = 52*i;
    m = 26*(i-1)+1;
    n = 26*i;

    p = SW_drought3_r(:,k:l);

    SW_drought3_summer(:,m:n) = p(:,14:39);
end



%Analyzing Dependency among compound Drought and Heatwaves

heatwaveEvents80_w = zeros(2256,1144);

for i=1:2256
    for j=1:1144
        if heatwaveEvents_80_rw_1980(i,j)>=2
            heatwaveEvents80_w(i,j)=1;
        end
    end
end

%heatwaveEvents90_wsum = (sum(heatwaveEvents90_w,1))';

%SW_drought3_r_sum     = (sum(SW_drought3_r,1))';

Compound_dr_hw  = zeros(2256,1144);

for i=1:2256
    for j=1:1144
        if heatwaveEvents80_w(i,j)==1 && SW_drought3_summer(i,j)==1
            Compound_dr_hw(i,j)=1;
        end
    end
end

Compound_dr_hw_sum    = (sum(Compound_dr_hw,1))';

save('Nighttime_CDHW_only_occurrence_matrix','Compound_dr_hw');


%%Calculate the heatwave component of Compound Drought & Heatwave severity

SWheat_normf  =  zeros(2256,183);

for i=1:2256
    for j=1:183
        SWheat_normf(i,j)=SWheat_75(i,j)-SWheat_25(i,j);
    end
end

SWheat_numerator = zeros(2256,8235);

for i=1:45
    for j=1:2256
        k = 183*(i-1)+1;
        l = 183*i;

        SWheat_numerator(j,k:l) = SWNAheat(j,k:l)-SWheat_25(j,:);
    end
end




CDHW_hw_sev = zeros(2256,8235);

for i=1:45
    k = 183*(i-1)+1;
    l = 183*i;
    for j=1:2256

    CDHW_hw_sev(j,k:l) =  SWheat_numerator(j,k:l)./SWheat_normf(j,:);
    end

end

%Multiplying the severity with heatwave occurrence to extract the severity
%corresponding to only heatwave occurrence
heatwaveEvents_80_T =  heatwaveEvents_80';

CDHW_hw_sev_occ = zeros(2256,8235);

for i=1:8235
  
    for j=1:2256

    CDHW_hw_sev_occ(j,i) =  CDHW_hw_sev(j,i)*heatwaveEvents_80_T(j,i);
    end

end




%Extracting the severity from 1980

CDHW_hw_sev_1980 = CDHW_hw_sev_occ(:,184:8235);



CDHW_hw_sev_r    = zeros(2256,183,44);

for i=1:44
    k  =  183*(i-1)+1;
    l  =  183*i;

    CDHW_hw_sev_r(:,:,i) =  CDHW_hw_sev_1980(:,k:l);
end


CDHW_hw_sev_r   = CDHW_hw_sev_r(:,1:182,:);

CDHW_hw_sev_r   = reshape(CDHW_hw_sev_r,[2256,182*44]);

CDHW_hw_sev_rw  = zeros(2256,1144);

for i=1:1144
    for j=1:2256
        k = 7*(i-1)+1;
        l = 7*i;

        CDHW_hw_sev_rw(j,i)= sum(CDHW_hw_sev_r(j,k:l));
    end
end

%Extracting summer drought characteristics (lets keep it in comments for
%the time)
SW_SPI3_summer = zeros(1144,2256);

for i=1:44
    m   = 52*(i-1)+1;
    n   = 52*i;

    p   = 26*(i-1)+1;
    q   = 26*i;
    l   = SW_SPI3(m:n,:);
    SW_SPI3_summer(p:q,:) = l(14:39,:);
end


SW_SPI3_summer = SW_SPI3_summer';


CDHW_hw_dr3_sev  = zeros(2256,1144);

for i=1:2256
    for j=1:1144

        CDHW_hw_dr3_sev(i,j) = CDHW_hw_sev_rw(i,j)*SW_SPI3_summer(i,j)*SW_drought3_summer(i,j)*(-1);
    end
end



save('nighttime_CDHW_severity_matrix_05092024.mat','CDHW_hw_dr3_sev');



%Calculating Spatiotemporal Characteristics of Compound Drought & Heatwave

Compound_dr_hw    = Compound_dr_hw';
CDHW_hw_dr3_sev   = CDHW_hw_dr3_sev';




numLocations = 2256;
numWeeks     = 26;

numDroughtEvents = zeros(1, numLocations);
cDroughtSeverity = zeros(numWeeks, numLocations);
cDroughtDuration = zeros(numWeeks, numLocations);
cDroughtIntensity = zeros(numWeeks, numLocations);




for loc = 1:numLocations
    droughtPeriods = Compound_dr_hw(1119:1144, loc);      %for 2023, 1119:1144
    droughtspi     = CDHW_hw_dr3_sev(1119:1144,loc);
    % ... [any additional code for specific periods] ...

    % Find the start and end indices of each drought event
    startIndices = find(diff([0; droughtPeriods]) == 1);
    endIndices = find(diff([droughtPeriods; 0]) == -1);

    % Initialize the number of events for this location
    numEvents = 0;
    
    for i = 1:length(startIndices)
        eventDuration = endIndices(i) - startIndices(i) + 1;
        
        % Check if the event duration is at least two weeks
        if eventDuration >= 1
            numEvents = numEvents + 1; % Counting the valid event
            cDroughtDuration(startIndices(i), loc) = eventDuration;
            
            eventSeverity = sum(droughtspi(startIndices(i):endIndices(i)));
            eventIntensity = eventSeverity / eventDuration;

            cDroughtSeverity(startIndices(i), loc) = eventSeverity;
            cDroughtIntensity(startIndices(i), loc) = eventIntensity;
        end
    end

    % Store the number of valid drought events for this location
    numDroughtEvents(loc) = numEvents;
end



%Calculating Average Drought Duration, Severity and Intensity
numDroughtEvents         =  numDroughtEvents';
SW_dr_avg_duration       =  zeros(2256,1);
SW_dr_avg_severity       =  zeros(2256,1);
SW_dr_avg_intensity      =  zeros(2256,1);

for i=1:2256
    a  =  find(cDroughtDuration(:,i)>0);
    b  =  cDroughtDuration(a,i);
    b  =  mean(b,1);
    SW_dr_avg_duration(i,1)=b;
    c  =  cDroughtSeverity(a,i);
    c  =  mean(c,1);
    SW_dr_avg_severity(i,1)=c;
    d  =  cDroughtIntensity(a,i);
    d  =  mean(d,1);
    SW_dr_avg_intensity(i,1)=d;
end

    
a =find(isnan(SW_dr_avg_duration));

SW_dr_avg_duration(a,:) = 0;

SW_dr_avg_intensity(a,:) = 0;

SW_dr_avg_severity(a,:) = 0;




%Calculating the decadal changes in Compound Drought and Heatwave
%Characteristics

%1980-1990
numLocations = 2256;
numWeeks     = 286;

numDroughtEvents1 = zeros(1, numLocations);
cDroughtSeverity1 = zeros(numWeeks, numLocations);
cDroughtDuration1 = zeros(numWeeks, numLocations);
cDroughtIntensity1 = zeros(numWeeks, numLocations);




for loc = 1:numLocations
    droughtPeriods = Compound_dr_hw(1:286, loc);      %for 2020, 2081:2132 & %2237:2288
    droughtspi     = CDHW_hw_dr3_sev(1:286,loc);
    % ... [any additional code for specific periods] ...

    % Find the start and end indices of each drought event
    startIndices = find(diff([0; droughtPeriods]) == 1);
    endIndices = find(diff([droughtPeriods; 0]) == -1);

    % Initialize the number of events for this location
    numEvents = 0;
    
    for i = 1:length(startIndices)
        eventDuration = endIndices(i) - startIndices(i) + 1;
        
        % Check if the event duration is at least two weeks
        if eventDuration >= 1
            numEvents = numEvents + 1; % Counting the valid event
            cDroughtDuration1(startIndices(i), loc) = eventDuration;
            
            eventSeverity = sum(droughtspi(startIndices(i):endIndices(i)));
            eventIntensity = eventSeverity / eventDuration;

            cDroughtSeverity1(startIndices(i), loc) = eventSeverity;
            cDroughtIntensity1(startIndices(i), loc) = eventIntensity;
        end
    end

    % Store the number of valid drought events for this location
    numDroughtEvents1(loc) = numEvents;
end



%Calculating Average Drought Duration, Severity and Intensity
numDroughtEvents1         =  numDroughtEvents1';
SW_dr_avg_duration1       =  zeros(2256,1);
SW_dr_avg_severity1       =  zeros(2256,1);
SW_dr_avg_intensity1      =  zeros(2256,1);

for i=1:2256
    a  =  find(cDroughtDuration1(:,i)>0);
    b  =  cDroughtDuration1(a,i);
    b  =  mean(b,1);
    SW_dr_avg_duration1(i,1)=b;
    c  =  cDroughtSeverity1(a,i);
    c  =  mean(c,1);
    SW_dr_avg_severity1(i,1)=c;
    d  =  cDroughtIntensity1(a,i);
    d  =  mean(d,1);
    SW_dr_avg_intensity1(i,1)=d;
end

    
a =find(isnan(SW_dr_avg_duration1));

SW_dr_avg_duration1(a,:) = 0;

SW_dr_avg_intensity1(a,:) = 0;

SW_dr_avg_severity1(a,:) = 0;


%1991-2001
numLocations = 2256;
numWeeks     = 286;

numDroughtEvents2 = zeros(1, numLocations);
cDroughtSeverity2 = zeros(numWeeks, numLocations);
cDroughtDuration2 = zeros(numWeeks, numLocations);
cDroughtIntensity2 = zeros(numWeeks, numLocations);




for loc = 1:numLocations
    droughtPeriods = Compound_dr_hw(287:572, loc);      %for 2020, 2081:2132 & %2237:2288
    droughtspi     = CDHW_hw_dr3_sev(287:572,loc);
    % ... [any additional code for specific periods] ...

    % Find the start and end indices of each drought event
    startIndices = find(diff([0; droughtPeriods]) == 1);
    endIndices = find(diff([droughtPeriods; 0]) == -1);

    % Initialize the number of events for this location
    numEvents = 0;
    
    for i = 1:length(startIndices)
        eventDuration = endIndices(i) - startIndices(i) + 1;
        
        % Check if the event duration is at least two weeks
        if eventDuration >= 1
            numEvents = numEvents + 1; % Counting the valid event
            cDroughtDuration2(startIndices(i), loc) = eventDuration;
            
            eventSeverity = sum(droughtspi(startIndices(i):endIndices(i)));
            eventIntensity = eventSeverity / eventDuration;

            cDroughtSeverity2(startIndices(i), loc) = eventSeverity;
            cDroughtIntensity2(startIndices(i), loc) = eventIntensity;
        end
    end

    % Store the number of valid drought events for this location
    numDroughtEvents2(loc) = numEvents;
end



%Calculating Average Drought Duration, Severity and Intensity
numDroughtEvents2         =  numDroughtEvents2';
SW_dr_avg_duration2       =  zeros(2256,1);
SW_dr_avg_severity2       =  zeros(2256,1);
SW_dr_avg_intensity2      =  zeros(2256,1);

for i=1:2256
    a  =  find(cDroughtDuration2(:,i)>0);
    b  =  cDroughtDuration2(a,i);
    b  =  mean(b,1);
    SW_dr_avg_duration2(i,1)=b;
    c  =  cDroughtSeverity2(a,i);
    c  =  mean(c,1);
    SW_dr_avg_severity2(i,1)=c;
    d  =  cDroughtIntensity2(a,i);
    d  =  mean(d,1);
    SW_dr_avg_intensity2(i,1)=d;
end

    
a =find(isnan(SW_dr_avg_duration2));

SW_dr_avg_duration2(a,:) = 0;

SW_dr_avg_intensity2(a,:) = 0;

SW_dr_avg_severity2(a,:) = 0;


%2002-2012
numLocations = 2256;
numWeeks     = 286;

numDroughtEvents3 = zeros(1, numLocations);
cDroughtSeverity3 = zeros(numWeeks, numLocations);
cDroughtDuration3 = zeros(numWeeks, numLocations);
cDroughtIntensity3 = zeros(numWeeks, numLocations);




for loc = 1:numLocations
    droughtPeriods = Compound_dr_hw(573:858, loc);      %for 2020, 2081:2132 & %2237:2288
    droughtspi     = CDHW_hw_dr3_sev(573:858,loc);
    % ... [any additional code for specific periods] ...

    % Find the start and end indices of each drought event
    startIndices = find(diff([0; droughtPeriods]) == 1);
    endIndices = find(diff([droughtPeriods; 0]) == -1);

    % Initialize the number of events for this location
    numEvents = 0;
    
    for i = 1:length(startIndices)
        eventDuration = endIndices(i) - startIndices(i) + 1;
        
        % Check if the event duration is at least two weeks
        if eventDuration >= 1
            numEvents = numEvents + 1; % Counting the valid event
            cDroughtDuration3(startIndices(i), loc) = eventDuration;
            
            eventSeverity = sum(droughtspi(startIndices(i):endIndices(i)));
            eventIntensity = eventSeverity / eventDuration;

            cDroughtSeverity3(startIndices(i), loc) = eventSeverity;
            cDroughtIntensity3(startIndices(i), loc) = eventIntensity;
        end
    end

    % Store the number of valid drought events for this location
    numDroughtEvents3(loc) = numEvents;
end



%Calculating Average Drought Duration, Severity and Intensity
numDroughtEvents3         =  numDroughtEvents3';
SW_dr_avg_duration3       =  zeros(2256,1);
SW_dr_avg_severity3       =  zeros(2256,1);
SW_dr_avg_intensity3      =  zeros(2256,1);

for i=1:2256
    a  =  find(cDroughtDuration3(:,i)>0);
    b  =  cDroughtDuration3(a,i);
    b  =  mean(b,1);
    SW_dr_avg_duration3(i,1)=b;
    c  =  cDroughtSeverity3(a,i);
    c  =  mean(c,1);
    SW_dr_avg_severity3(i,1)=c;
    d  =  cDroughtIntensity3(a,i);
    d  =  mean(d,1);
    SW_dr_avg_intensity3(i,1)=d;
end

    
a =find(isnan(SW_dr_avg_duration3));

SW_dr_avg_duration3(a,:) = 0;

SW_dr_avg_intensity3(a,:) = 0;

SW_dr_avg_severity3(a,:) = 0;


%2013-2023
numLocations = 2256;
numWeeks     = 286;

numDroughtEvents4 = zeros(1, numLocations);
cDroughtSeverity4 = zeros(numWeeks, numLocations);
cDroughtDuration4 = zeros(numWeeks, numLocations);
cDroughtIntensity4 = zeros(numWeeks, numLocations);




for loc = 1:numLocations
    droughtPeriods = Compound_dr_hw(859:1144, loc);      %for 2020, 2081:2132 & %2237:2288
    droughtspi     = CDHW_hw_dr3_sev(859:1144,loc);
    % ... [any additional code for specific periods] ...

    % Find the start and end indices of each drought event
    startIndices = find(diff([0; droughtPeriods]) == 1);
    endIndices = find(diff([droughtPeriods; 0]) == -1);

    % Initialize the number of events for this location
    numEvents = 0;
    
    for i = 1:length(startIndices)
        eventDuration = endIndices(i) - startIndices(i) + 1;
        
        % Check if the event duration is at least two weeks
        if eventDuration >= 1
            numEvents = numEvents + 1; % Counting the valid event
            cDroughtDuration4(startIndices(i), loc) = eventDuration;
            
            eventSeverity = sum(droughtspi(startIndices(i):endIndices(i)));
            eventIntensity = eventSeverity / eventDuration;

            cDroughtSeverity4(startIndices(i), loc) = eventSeverity;
            cDroughtIntensity4(startIndices(i), loc) = eventIntensity;
        end
    end

    % Store the number of valid drought events for this location
    numDroughtEvents4(loc) = numEvents;
end



%Calculating Average Drought Duration, Severity and Intensity
numDroughtEvents4         =  numDroughtEvents4';
SW_dr_avg_duration4       =  zeros(2256,1);
SW_dr_avg_severity4       =  zeros(2256,1);
SW_dr_avg_intensity4      =  zeros(2256,1);

for i=1:2256
    a  =  find(cDroughtDuration4(:,i)>0);
    b  =  cDroughtDuration4(a,i);
    b  =  mean(b,1);
    SW_dr_avg_duration4(i,1)=b;
    c  =  cDroughtSeverity4(a,i);
    c  =  mean(c,1);
    SW_dr_avg_severity4(i,1)=c;
    d  =  cDroughtIntensity4(a,i);
    d  =  mean(d,1);
    SW_dr_avg_intensity4(i,1)=d;
end

    
a =find(isnan(SW_dr_avg_duration4));

SW_dr_avg_duration4(a,:) = 0;

SW_dr_avg_intensity4(a,:) = 0;

SW_dr_avg_severity4(a,:) = 0;




%%Analyzing 2023 Compound Drought and Heatwave Severity time series of 2023
%%with respect to other decades


%Getting Grid-locations corresponding to 4-corner states and Mexico

SWNA             =  readtable('/Users/smonda15/Ongoing Work/Publication_Compound_Drought_Heatwave_Over_Phoneix/Datasets/Shapefiles/Southwest_north_America_Grid_Locations.csv');
%SWNA_latlon      =  SWNA(:,1:2);
%SWNA_index       =  SWNA(:,3);
States           =  SWNA(:,19);

index_AZ = find(strcmp(States.NAME, 'Arizona'));
index_UT = find(strcmp(States.NAME, 'Utah'));
index_NM = find(strcmp(States.NAME, 'New Mexico'));
index_CO = find(strcmp(States.NAME, 'Colorado'));
index_MX = find(strcmp(States.NAME, ''));


%Day & Night severity 

day_cdhw_az_2023 = mean(day_cdhw(1119:1144,index_AZ),2);
day_cdhw_mx_2023 = mean(day_cdhw(1119:1144,index_MX),2);
day_cdhw_nm_2023 = mean(day_cdhw(1119:1144,index_NM),2);


night_cdhw_az_2023 = mean(night_cdhw(1119:1144,index_AZ),2);
night_cdhw_mx_2023 = mean(night_cdhw(1119:1144,index_MX),2);
night_cdhw_nm_2023 = mean(night_cdhw(1119:1144,index_NM),2);


AZ = [day_cdhw_az_2023,night_cdhw_az_2023];
NM = [day_cdhw_nm_2023,night_cdhw_nm_2023];
MX = [day_cdhw_mx_2023,night_cdhw_mx_2023];

%
CDHW_hw_dr3_sev_AZ  =  CDHW_hw_dr3_sev(:,index_AZ);
CDHW_hw_dr3_sev_UT  =  CDHW_hw_dr3_sev(:,index_UT);
CDHW_hw_dr3_sev_NM  =  CDHW_hw_dr3_sev(:,index_NM);
CDHW_hw_dr3_sev_CO  =  CDHW_hw_dr3_sev(:,index_CO);
CDHW_hw_dr3_sev_MX  =  CDHW_hw_dr3_sev(:,index_MX);

%Converting to Areally averaged values

CDHW_hw_dr3_sev_AZ_avg = mean(CDHW_hw_dr3_sev_AZ,2);
CDHW_hw_dr3_sev_UT_avg = mean(CDHW_hw_dr3_sev_UT,2);
CDHW_hw_dr3_sev_NM_avg = mean(CDHW_hw_dr3_sev_NM,2);
CDHW_hw_dr3_sev_CO_avg = mean(CDHW_hw_dr3_sev_CO,2);
CDHW_hw_dr3_sev_MX_avg = mean(CDHW_hw_dr3_sev_MX,2);

% Arizona
CDHW_hw_dr3_sev_AZ_avg_1990 = CDHW_hw_dr3_sev_AZ_avg(27:286);
CDHW_hw_dr3_sev_AZ_avg_2000 = CDHW_hw_dr3_sev_AZ_avg(287:546);
CDHW_hw_dr3_sev_AZ_avg_2010 = CDHW_hw_dr3_sev_AZ_avg(547:806);
CDHW_hw_dr3_sev_AZ_avg_2020 = CDHW_hw_dr3_sev_AZ_avg(807:1066);

CDHW_hw_dr3_sev_AZ_avg_1990 = reshape(CDHW_hw_dr3_sev_AZ_avg_1990,[26,10]);
CDHW_hw_dr3_sev_AZ_avg_2000 = reshape(CDHW_hw_dr3_sev_AZ_avg_2000,[26,10]);
CDHW_hw_dr3_sev_AZ_avg_2010 = reshape(CDHW_hw_dr3_sev_AZ_avg_2010,[26,10]);
CDHW_hw_dr3_sev_AZ_avg_2020 = reshape(CDHW_hw_dr3_sev_AZ_avg_2020,[26,10]);

CDHW_hw_dr3_sev_AZ_avg_1990 = mean(CDHW_hw_dr3_sev_AZ_avg_1990,2);
CDHW_hw_dr3_sev_AZ_avg_2000 = mean(CDHW_hw_dr3_sev_AZ_avg_2000,2);
CDHW_hw_dr3_sev_AZ_avg_2010 = mean(CDHW_hw_dr3_sev_AZ_avg_2010,2);
CDHW_hw_dr3_sev_AZ_avg_2020 = mean(CDHW_hw_dr3_sev_AZ_avg_2020,2);

CDHW_hw_dr3_sev_AZ_avg_2023 = CDHW_hw_dr3_sev_AZ_avg(1119:1144);

CDHW_AZ  = [CDHW_hw_dr3_sev_AZ_avg_1990,CDHW_hw_dr3_sev_AZ_avg_2000,CDHW_hw_dr3_sev_AZ_avg_2010,CDHW_hw_dr3_sev_AZ_avg_2020,CDHW_hw_dr3_sev_AZ_avg_2023];



%Utah

CDHW_hw_dr3_sev_UT_avg_1990 = CDHW_hw_dr3_sev_UT_avg(27:286);
CDHW_hw_dr3_sev_UT_avg_2000 = CDHW_hw_dr3_sev_UT_avg(287:546);
CDHW_hw_dr3_sev_UT_avg_2010 = CDHW_hw_dr3_sev_UT_avg(547:806);
CDHW_hw_dr3_sev_UT_avg_2020 = CDHW_hw_dr3_sev_UT_avg(807:1066);

CDHW_hw_dr3_sev_UT_avg_1990 = reshape(CDHW_hw_dr3_sev_UT_avg_1990,[26,10]);
CDHW_hw_dr3_sev_UT_avg_2000 = reshape(CDHW_hw_dr3_sev_UT_avg_2000,[26,10]);
CDHW_hw_dr3_sev_UT_avg_2010 = reshape(CDHW_hw_dr3_sev_UT_avg_2010,[26,10]);
CDHW_hw_dr3_sev_UT_avg_2020 = reshape(CDHW_hw_dr3_sev_UT_avg_2020,[26,10]);

CDHW_hw_dr3_sev_UT_avg_1990 = mean(CDHW_hw_dr3_sev_UT_avg_1990,2);
CDHW_hw_dr3_sev_UT_avg_2000 = mean(CDHW_hw_dr3_sev_UT_avg_2000,2);
CDHW_hw_dr3_sev_UT_avg_2010 = mean(CDHW_hw_dr3_sev_UT_avg_2010,2);
CDHW_hw_dr3_sev_UT_avg_2020 = mean(CDHW_hw_dr3_sev_UT_avg_2020,2);

CDHW_hw_dr3_sev_UT_avg_2023 = CDHW_hw_dr3_sev_UT_avg(1119:1144);

CDHW_UT  = [CDHW_hw_dr3_sev_UT_avg_1990,CDHW_hw_dr3_sev_UT_avg_2000,CDHW_hw_dr3_sev_UT_avg_2010,CDHW_hw_dr3_sev_UT_avg_2020,CDHW_hw_dr3_sev_UT_avg_2023];


%New Mexico

CDHW_hw_dr3_sev_NM_avg_1990 = CDHW_hw_dr3_sev_NM_avg(27:286);
CDHW_hw_dr3_sev_NM_avg_2000 = CDHW_hw_dr3_sev_NM_avg(287:546);
CDHW_hw_dr3_sev_NM_avg_2010 = CDHW_hw_dr3_sev_NM_avg(547:806);
CDHW_hw_dr3_sev_NM_avg_2020 = CDHW_hw_dr3_sev_NM_avg(807:1066);

CDHW_hw_dr3_sev_NM_avg_1990 = reshape(CDHW_hw_dr3_sev_NM_avg_1990,[26,10]);
CDHW_hw_dr3_sev_NM_avg_2000 = reshape(CDHW_hw_dr3_sev_NM_avg_2000,[26,10]);
CDHW_hw_dr3_sev_NM_avg_2010 = reshape(CDHW_hw_dr3_sev_NM_avg_2010,[26,10]);
CDHW_hw_dr3_sev_NM_avg_2020 = reshape(CDHW_hw_dr3_sev_NM_avg_2020,[26,10]);

CDHW_hw_dr3_sev_NM_avg_1990 = mean(CDHW_hw_dr3_sev_NM_avg_1990,2);
CDHW_hw_dr3_sev_NM_avg_2000 = mean(CDHW_hw_dr3_sev_NM_avg_2000,2);
CDHW_hw_dr3_sev_NM_avg_2010 = mean(CDHW_hw_dr3_sev_NM_avg_2010,2);
CDHW_hw_dr3_sev_NM_avg_2020 = mean(CDHW_hw_dr3_sev_NM_avg_2020,2);

CDHW_hw_dr3_sev_NM_avg_2023 = CDHW_hw_dr3_sev_NM_avg(1119:1144);

CDHW_NM  = [CDHW_hw_dr3_sev_NM_avg_1990,CDHW_hw_dr3_sev_NM_avg_2000,CDHW_hw_dr3_sev_NM_avg_2010,CDHW_hw_dr3_sev_NM_avg_2020,CDHW_hw_dr3_sev_NM_avg_2023];

%Colorado

CDHW_hw_dr3_sev_CO_avg_1990 = CDHW_hw_dr3_sev_CO_avg(27:286);
CDHW_hw_dr3_sev_CO_avg_2000 = CDHW_hw_dr3_sev_CO_avg(287:546);
CDHW_hw_dr3_sev_CO_avg_2010 = CDHW_hw_dr3_sev_CO_avg(547:806);
CDHW_hw_dr3_sev_CO_avg_2020 = CDHW_hw_dr3_sev_CO_avg(807:1066);

CDHW_hw_dr3_sev_CO_avg_1990 = reshape(CDHW_hw_dr3_sev_CO_avg_1990,[26,10]);
CDHW_hw_dr3_sev_CO_avg_2000 = reshape(CDHW_hw_dr3_sev_CO_avg_2000,[26,10]);
CDHW_hw_dr3_sev_CO_avg_2010 = reshape(CDHW_hw_dr3_sev_CO_avg_2010,[26,10]);
CDHW_hw_dr3_sev_CO_avg_2020 = reshape(CDHW_hw_dr3_sev_CO_avg_2020,[26,10]);

CDHW_hw_dr3_sev_CO_avg_1990 = mean(CDHW_hw_dr3_sev_CO_avg_1990,2);
CDHW_hw_dr3_sev_CO_avg_2000 = mean(CDHW_hw_dr3_sev_CO_avg_2000,2);
CDHW_hw_dr3_sev_CO_avg_2010 = mean(CDHW_hw_dr3_sev_CO_avg_2010,2);
CDHW_hw_dr3_sev_CO_avg_2020 = mean(CDHW_hw_dr3_sev_CO_avg_2020,2);

CDHW_hw_dr3_sev_CO_avg_2023 = CDHW_hw_dr3_sev_CO_avg(1119:1144);

CDHW_CO  = [CDHW_hw_dr3_sev_CO_avg_1990,CDHW_hw_dr3_sev_CO_avg_2000,CDHW_hw_dr3_sev_CO_avg_2010,CDHW_hw_dr3_sev_CO_avg_2020,CDHW_hw_dr3_sev_CO_avg_2023];


%Mexico

CDHW_hw_dr3_sev_MX_avg_1990 = CDHW_hw_dr3_sev_MX_avg(27:286);
CDHW_hw_dr3_sev_MX_avg_2000 = CDHW_hw_dr3_sev_MX_avg(287:546);
CDHW_hw_dr3_sev_MX_avg_2010 = CDHW_hw_dr3_sev_MX_avg(547:806);
CDHW_hw_dr3_sev_MX_avg_2020 = CDHW_hw_dr3_sev_MX_avg(807:1066);

CDHW_hw_dr3_sev_MX_avg_1990 = reshape(CDHW_hw_dr3_sev_MX_avg_1990,[26,10]);
CDHW_hw_dr3_sev_MX_avg_2000 = reshape(CDHW_hw_dr3_sev_MX_avg_2000,[26,10]);
CDHW_hw_dr3_sev_MX_avg_2010 = reshape(CDHW_hw_dr3_sev_MX_avg_2010,[26,10]);
CDHW_hw_dr3_sev_MX_avg_2020 = reshape(CDHW_hw_dr3_sev_MX_avg_2020,[26,10]);

CDHW_hw_dr3_sev_MX_avg_1990 = mean(CDHW_hw_dr3_sev_MX_avg_1990,2);
CDHW_hw_dr3_sev_MX_avg_2000 = mean(CDHW_hw_dr3_sev_MX_avg_2000,2);
CDHW_hw_dr3_sev_MX_avg_2010 = mean(CDHW_hw_dr3_sev_MX_avg_2010,2);
CDHW_hw_dr3_sev_MX_avg_2020 = mean(CDHW_hw_dr3_sev_MX_avg_2020,2);

CDHW_hw_dr3_sev_MX_avg_2023 = CDHW_hw_dr3_sev_MX_avg(1119:1144);

CDHW_MX  = [CDHW_hw_dr3_sev_MX_avg_1990,CDHW_hw_dr3_sev_MX_avg_2000,CDHW_hw_dr3_sev_MX_avg_2010,CDHW_hw_dr3_sev_MX_avg_2020,CDHW_hw_dr3_sev_MX_avg_2023];



%Creating Exceedence probability plot considering both spatial extent and
%area averaged (averaged over the area of occurrence)of CDHW

SWNA                =  readtable('/Users/smonda15/Ongoing Work/Publication_Compound_Drought_Heatwave_Over_Phoneix/Datasets/Shapefiles/Southwest_north_America_Grid_Locations.csv');
%SWNA_latlon      =  SWNA(:,1:2);
%SWNA_index       =  SWNA(:,3);
States              =  SWNA(:,19);

index_AZ            = find(strcmp(States.NAME, 'Arizona'));


Compound_dr_hw_AZ   = Compound_dr_hw(:,index_AZ);
CDHW_hw_dr3_sev_AZ  =  CDHW_hw_dr3_sev(:,index_AZ);


Compound_dr_hw_AZ_ar_avg   = sum(Compound_dr_hw_AZ,2);
CDHW_hw_dr3_sev_AZ_ar_avg  = sum(CDHW_hw_dr3_sev_AZ,2);

CDHW_hw_dr3_sev_AZ_ar_avg1 = zeros(1144,1);

for i=1:1144
    if Compound_dr_hw_AZ_ar_avg(i,1)~=0
        CDHW_hw_dr3_sev_AZ_ar_avg1(i,1) = CDHW_hw_dr3_sev_AZ_ar_avg(i,1)/Compound_dr_hw_AZ_ar_avg(i,1);
    else
        CDHW_hw_dr3_sev_AZ_ar_avg1(i,1) = 0;
    end
end


CDHW = [CDHW_hw_dr3_sev_AZ_ar_avg1,Compound_dr_hw_AZ_ar_avg];

CDHW_AZ_max_sev = zeros(44,1);

for i=1:44
    j      =  26*(i-1)+1;
    k      =  26*i;

    p      =  CDHW_hw_dr3_sev_AZ_ar_avg1(j:k,1);

    indexp =  find(p>3);

    q      =  Compound_dr_hw_AZ_ar_avg(j:k,1);

    q_p    = q(indexp,1);

    if size(indexp,1)~=0
        CDHW_AZ_max_sev(i,1) = max(q_p);
    else
        CDHW_AZ_max_sev(i,1)=0;
    end
end

CDHW_AZ_max_sev1 = zeros(44,1);
for i=1:44
    j   =  26*(i-1)+1;
    k   =  26*i;

    CDHW_AZ_max_sev1(i,1)=max(CDHW_hw_dr3_sev_AZ_ar_avg(j:k,1));
end






Parameters_CDHW_weibul_AZ=wblfit(CDHW_AZ_max_sev);
Return_Period_wb_AZ  = zeros(44,1);
for i=1:44
p = wblcdf(CDHW_AZ_max_sev(i,1),Parameters_CDHW_weibul_AZ(1,1),Parameters_CDHW_weibul_AZ(1,2));

Return_Period_wb_AZ(i,1)=1/(1-p);

end




Parameters_CDHW_weibul_AZ1=wblfit(CDHW_AZ_max_sev1);
Return_Period_wb_AZ1  = zeros(44,1);
for i=1:44
p = wblcdf(CDHW_AZ_max_sev1(i,1),Parameters_CDHW_weibul_AZ1(1,1),Parameters_CDHW_weibul_AZ1(1,2));

Return_Period_wb_AZ1(i,1)=1/(1-p);

end











% 9Arizona
CDHW_hw_dr3_sev_AZ_avg = mean(CDHW_hw_dr3_sev_AZ,2);
CDHW_hw_dr3_sev_AZ_avg_1990 = CDHW_hw_dr3_sev_AZ_avg(27:286);
CDHW_hw_dr3_sev_AZ_avg_2000 = CDHW_hw_dr3_sev_AZ_avg(287:546);
CDHW_hw_dr3_sev_AZ_avg_2010 = CDHW_hw_dr3_sev_AZ_avg(547:806);
CDHW_hw_dr3_sev_AZ_avg_2020 = CDHW_hw_dr3_sev_AZ_avg(807:1066);

CDHW_hw_dr3_sev_AZ_avg_1990 = reshape(CDHW_hw_dr3_sev_AZ_avg_1990,[26,10]);
CDHW_hw_dr3_sev_AZ_avg_2000 = reshape(CDHW_hw_dr3_sev_AZ_avg_2000,[26,10]);
CDHW_hw_dr3_sev_AZ_avg_2010 = reshape(CDHW_hw_dr3_sev_AZ_avg_2010,[26,10]);
CDHW_hw_dr3_sev_AZ_avg_2020 = reshape(CDHW_hw_dr3_sev_AZ_avg_2020,[26,10]);

CDHW_hw_dr3_sev_AZ_avg_1990 = mean(CDHW_hw_dr3_sev_AZ_avg_1990,2);
CDHW_hw_dr3_sev_AZ_avg_2000 = mean(CDHW_hw_dr3_sev_AZ_avg_2000,2);
CDHW_hw_dr3_sev_AZ_avg_2010 = mean(CDHW_hw_dr3_sev_AZ_avg_2010,2);
CDHW_hw_dr3_sev_AZ_avg_2020 = mean(CDHW_hw_dr3_sev_AZ_avg_2020,2);

CDHW_hw_dr3_sev_AZ_avg_2023 = CDHW_hw_dr3_sev_AZ_avg(1119:1144);

CDHW_AZ  = [CDHW_hw_dr3_sev_AZ_avg_1990,CDHW_hw_dr3_sev_AZ_avg_2000,CDHW_hw_dr3_sev_AZ_avg_2010,CDHW_hw_dr3_sev_AZ_avg_2020,CDHW_hw_dr3_sev_AZ_avg_2023];


%Plotting severity vs return period

CDHW_max_sev_AZ                =    CDHW_max_sev(index_AZ,:);

Return_Period_wb_all_years_AZ  =    Return_Period_wb_all_years(index_AZ,:);


CDHW_max_sev_AZ_mean                = (mean(CDHW_max_sev_AZ,1))';
Return_Period_wb_all_years_AZ_mean  = (mean(Return_Period_wb_all_years_AZ,1))';

loglog(Return_Period_wb_all_years_AZ_mean, CDHW_max_sev_AZ_mean , 'o', 'MarkerSize', 8, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');


b = [CDHW_max_sev_AZ_mean,Return_Period_wb_all_years_AZ_mean];

CDHW_max_sev_MX                =    CDHW_max_sev(index_MX,:);

Return_Period_wb_all_years_MX  =    Return_Period_wb_all_years(index_MX,:);

CDHW_max_sev_MX_mean                = (mean(CDHW_max_sev_MX,1))';
Return_Period_wb_all_years_MX_mean  = (mean(Return_Period_wb_all_years_MX,1))';

% Your data vectors
x = 1980:2023;  % This is the x-axis common to both y-axes
y1 = CDHW_max_sev_MX_mean;  % First set of data
y2 = Return_Period_wb_all_years_MX_mean;  % Second set of data, which we'll plot against a second y-axis

% Create the figure
figure;

% Left y-axis
yyaxis left;
plot(x, y1, '-b', 'LineWidth', 2);  % '-b' makes the line blue
ylabel('Left Y-Axis Data');

% Right y-axis
yyaxis right;
plot(x, y2, '-r', 'LineWidth', 2);  % '-r' makes the line red
ylabel('Right Y-Axis Data');

% X-axis and title
xlabel('Common X-Axis');
title('Dual Y-Axis Plot');

% Customize axes
xlim([min(x) max(x)]);  % Set the x-axis limits to the range of x
yyaxis left;
ylim([0 max(y1)]);  % Set the y-axis limits for left y-axis to the range of y1
yyaxis right;
ylim([0 max(y2)]);  % Set the y-axis limits for right y-axis to the range of y2

% Adding grid and legend
grid on;
legend('Data Y1', 'Data Y2');


%Preparing the severity values for extreme value analysis

%Considering Compound_dr_hw and CDHW_hw_dr3_sev is already computed

numLocations = 2256;
numWeeks     = 1144;

numDroughtEvents4 = zeros(1, numLocations);
cDroughtSeverity = zeros(numWeeks, numLocations);
cDroughtDuration = zeros(numWeeks, numLocations);
cDroughtIntensity = zeros(numWeeks, numLocations);




for loc = 1:numLocations
    droughtPeriods = Compound_dr_hw(:, loc);      %for 2020, 2081:2132 & %2237:2288
    droughtspi     = CDHW_hw_dr3_sev(:,loc);
    % ... [any additional code for specific periods] ...

    % Find the start and end indices of each drought event
    startIndices = find(diff([0; droughtPeriods]) == 1);
    endIndices = find(diff([droughtPeriods; 0]) == -1);

    % Initialize the number of events for this location
    numEvents = 0;
    
    for i = 1:length(startIndices)
        eventDuration = endIndices(i) - startIndices(i) + 1;
        
        % Check if the event duration is at least two weeks
        if eventDuration >= 1
            numEvents = numEvents + 1; % Counting the valid event
            cDroughtDuration(startIndices(i), loc) = eventDuration;
            
            eventSeverity = sum(droughtspi(startIndices(i):endIndices(i)));
            eventIntensity = eventSeverity / eventDuration;

            cDroughtSeverity(startIndices(i), loc) = eventSeverity;
            cDroughtIntensity(startIndices(i), loc) = eventIntensity;
        end
    end

    % Store the number of valid drought events for this location
    numDroughtEvents4(loc) = numEvents;
end



%Calculating Average Drought Duration, Severity and Intensity
numDroughtEvents4         =  numDroughtEvents4';
SW_dr_avg_duration4       =  zeros(2256,1);
SW_dr_avg_severity4       =  zeros(2256,1);
SW_dr_avg_intensity4      =  zeros(2256,1);

for i=1:2256
    a  =  find(cDroughtDuration(:,i)>0);
    b  =  cDroughtDuration(a,i);
    b  =  mean(b,1);
    SW_dr_avg_duration4(i,1)=b;
    c  =  cDroughtSeverity(a,i);
    c  =  mean(c,1);
    SW_dr_avg_severity4(i,1)=c;
    d  =  cDroughtIntensity(a,i);
    d  =  mean(d,1);
    SW_dr_avg_intensity4(i,1)=d;
end

    
a =find(isnan(SW_dr_avg_duration4));

SW_dr_avg_duration4(a,:) = 0;

SW_dr_avg_intensity4(a,:) = 0;

SW_dr_avg_severity4(a,:) = 0;


%Extracting maximum severity corresponding to each grid-location each year

CDHW_max_sev = zeros(2256,44);

for i=1:2256
    for j=1:44
        k = 26*(j-1)+1;
        l = 26*j;
        a = cDroughtSeverity(k:l,i);
        b = find(a>0);
        if size(b,1)~=0
           CDHW_max_sev(i,j)= max(a(b,1));
        else
            CDHW_max_sev(i,j)=0;
        end
    end
end

% Calculate grid-locations with no-compound drought -heatwave for each year
% Initialize an array to hold the count of zeros for each year
zero_count_per_year = zeros(1, 44);

% Loop through each year (column) and count the zeros
for year = 1:44
    zero_count_per_year(year) = sum(CDHW_max_sev(:, year) == 0);
end

% Display the count of zeros for each year
disp(zero_count_per_year);

zero_count_per_year    =  zero_count_per_year';
zero_count_per_year1   =  (zero_count_per_year/2256)*100;



%Extreme Value Analysis

% Load your data: assuming it's in a matrix called 'data'
% where rows represent grid-locations and columns represent years
% data = load('your_data_file.mat');

% Initialize a matrix to store the annual maxima for each location
annual_maxima = max(data, [], 2);  % Assuming 'data' has zeros for no events


%Calculating occurrence of Compound Drought & Severity for each
%grid-location
zero_count_per_grid = zeros(2256,1);
for i=1:2256
    a = find(CDHW_max_sev(i,:)==0);
    zero_count_per_grid(i,1) = size(a,2);
end

%Doing Extreme Value Analysis

a = CDHW_max_sev(1535,:);

a_n = a(find(a>0));

a_n = a_n';

PARMHAT = gevfit(a_n);



% Assuming PARMHAT is obtained from your previous code:
K = PARMHAT(1);
sigma = PARMHAT(2);
mu = PARMHAT(3);

% Value x for which you want to calculate the return period
x = a_n(23,1); % Replace 'your_desired_value' with the actual value

% Calculate the exceedance probability using the GEV CDF
if xi ~= 0
    p = 1 - exp(-((1 + xi * ((x - mu) / sigma)).^(-1/xi)));
else
    % Limit as xi approaches 0 (Gumbel distribution)
    p = 1 - exp(-exp(-(x - mu) / sigma));
end

% Calculate the return period
T = 1 / p;



P = gevcdf(x,K,sigma,mu);
T1 = 1/(1-P);

%Normalizing Maximum Severity Values

CDHW_max_min = zeros(2256,2);

for i=1:2256
    a = (CDHW_max_sev(i,:))';
    a = a(find(a>0));
    if size(a,1)~=0
        CDHW_max_min(i,1) =  max(a);
        CDHW_max_min(i,2) =  min(a);
    else 
        CDHW_max_min(i,1)=0;
        CDHW_max_min(i,2)=0;
    end
end

CDHW_max_sev_norm = zeros(2256,44);

for i=1:2256
    for j=1:44
        if  CDHW_max_sev(i,j)~=0

        CDHW_max_sev_norm(i,j)=(CDHW_max_sev(i,j)-CDHW_max_min(i,2))/(CDHW_max_min(i,1)-CDHW_max_min(i,2));
        else
          CDHW_max_sev_norm(i,j) = CDHW_max_sev(i,j);  
    end
end
end






Parameters_CDHW = zeros(2256,3);

for i=1:2256

    a = (CDHW_max_sev(i,:))';
    a = a(find(a>0));
    Parameters_CDHW(i,:)=gevfit(a);
end


Return_Period = zeros(2256,1);

for i=1:2256
    x = CDHW_max_sev(i,44);
    p = gevcdf(x,Parameters_CDHW(i,1),Parameters_CDHW(i,2),Parameters_CDHW(i,3));
    Return_Period(i,1)=1/(1-p);
end


Return_Period1 = Return_Period(find(Return_Period<1000));




%Checking the distribution of the data% Your severity data

CDHW_model_selection_AIC = zeros(2256,4);

CDHW_model_selection_BIC = zeros(2256,4);

for i=1:2256
data = CDHW_max_sev_norm(i,:);
data = data';
data = data(data>0);

if size(data,1)>2

% GEV Fit
[paramGEV, ciGEV] = gevfit(data);
llGEV = -gevlike(paramGEV, data);

% Gumbel (as special case of GEV or using evfit)
[paramGumbel, ciGumbel] = evfit(data);
llGumbel = -evlike(paramGumbel, data);

% Weibull Fit
[paramWeibull, ciWeibull] = wblfit(data);
llWeibull = -wbllike(paramWeibull, data);

% Exponential Fit
[paramExp, ciExp] = expfit(data);
llExp = -explike(paramExp, data);

% Compare using AIC and BIC
n = numel(data); % Number of data points
k = 2; % Number of parameters for each distribution except GEV (which has 3)
aic = @(ll, k) 2*k - 2*ll;
bic = @(ll, k, n) k*log(n) - 2*ll;

AIC = [aic(llGEV, 3), aic(llGumbel, k), aic(llWeibull, k), aic(llExp, k)];
BIC = [bic(llGEV, 3, n), bic(llGumbel, k, n), bic(llWeibull, k, n), bic(llExp, k, n)];

CDHW_model_selection_AIC(i,:)= AIC;
CDHW_model_selection_BIC(i,:)= BIC;
else
CDHW_model_selection_AIC(i,:)= 0;
CDHW_model_selection_BIC(i,:)= 0;
end
end

CDHW_model_selection_AIC_mean = mean(CDHW_model_selection_AIC,1); 
CDHW_model_selection_BIC_mean = mean(CDHW_model_selection_BIC,1); 

% Print or plot the results for comparison
%disp('AIC values:');
%disp(AIC);
%disp('BIC values:');
%disp(BIC);



Parameters_CDHW_weibul = zeros(2256,2);

for i=1:2256

    a = (CDHW_max_sev(i,:))';
    a = a(find(a>0));
    Parameters_CDHW_weibul(i,:)=wblfit(a);
end


Return_Period_wb = zeros(2256,1);

for i=1:2256
    x = CDHW_max_sev(i,44);
    p = wblcdf(x,Parameters_CDHW_weibul(i,1),Parameters_CDHW_weibul(i,2));
    Return_Period_wb(i,1)=1/(1-p);
end


Return_Period1 = Return_Period(find(Return_Period<1000));


%Plotting the spatial Map

% Sample Data Preparation
%latlon = rand(2256, 2) * 180; % Replace with your actual lat-lon data
%latlon(:,1) = latlon(:,1) - 90; % Latitude adjustment
%latlon(:,2) = latlon(:,2) - 180; % Longitude adjustment
%return_periods = rand(2256, 1) * 100; % Replace with your actual return period data

% Begin Visualization
figure;

% Focusing Map on Southwestern North America
% These lat and lon boundaries can be adjusted to better fit your specific data extent
latlim = [20 45];  % Approximate latitude limits for Southwestern North America
lonlim = [-125 -100];  % Approximate longitude limits for Southwestern North America
worldmap(latlim, lonlim); % Sets the map extent to Southwestern North America
load coastlines % Loading coastlines for background

% Plotting adjusted coastlines within the specified limits
plotm(coastlat, coastlon, 'k') % Plotting coastlines

% Data Scatter Plot
% Ensure SWNA_latlon and Return_Period_wb variables are correctly loaded and prepared
scatterm(SWNA_latlon(:,2), SWNA_latlon(:,1), 40, Return_Period1, 'filled'); % Lat, Lon, Size, Data, Marker
%caxis([10 160]);
% Colorbar and its customization
cb = colorbar;
cb.Label.String = 'Return Period (years)';
cb.FontSize = 12; % Adjust font size as needed

% Enhancements for Publication Quality
title('Return Period Visualization for Southwestern North America', 'FontSize', 14);
% Note: For the worldmap function, labels for longitude and latitude are automatically managed.
set(gca, 'FontSize', 12); % Sets the default text sizes

% Colormap
colormap(jet); % You can choose other colormaps like 'hot', 'jet', etc.

% Save the Figure with High Resolution
print('SWNAReturnPeriodMap','-dpng','-r300'); % Saves the figure as a PNG with 300 DPI




% Sample data setup - replace these with your actual matrices
% Assuming size(CDHW_model_selection_AIC) == [2256, 4] for example

% Prepare figure
figure;

% Model names for x-ticks
modelNames = {'GEV', 'Gumbel', 'Weibull', 'Exponential'};

% Colors for AIC and BIC
colors = {'b', 'r'}; % Blue for AIC, Red for BIC

% Loop through each model to plot AIC and BIC side by side
for iModel = 1:length(modelNames)
    % AIC values
    aicValues = CDHW_model_selection_AIC(:, iModel);
    bplot1 = boxplot(aicValues, 'positions', [iModel-0.2], 'Widths', 0.35, 'Colors', colors{1}, 'Symbol', '');
    set(bplot1, {'linew'}, {2}); % Set line width
    hold on;
    
    % BIC values
    bicValues = CDHW_model_selection_BIC(:, iModel);
    bplot2 = boxplot(bicValues, 'positions', [iModel+0.2], 'Widths', 0.35, 'Colors', colors{2}, 'Symbol', '');
    set(bplot2, {'linew'}, {2}); % Set line width
end

% Customizations
set(gca, 'XTick', 1:length(modelNames), 'XTickLabel', modelNames, 'FontSize', 12);
ylabel('Criterion Value', 'FontSize', 12);
title('Comparison of Model Selection Criteria', 'FontSize', 14);

% Legend
h = zeros(2, 1);
h(1) = plot(NaN,NaN,'b-');
h(2) = plot(NaN,NaN,'r-');
legend(h, 'AIC','BIC', 'Location', 'best', 'FontSize', 12);

grid on;
hold off;

% Save the figure with high resolution
print('ModelSelectionComparison','-dpng','-r300');




Return_Period_wb_all_years = zeros(2256,44);

for i=1:2256
    for j=1:44
    x = CDHW_max_sev(i,j);
    p = wblcdf(x,Parameters_CDHW_weibul(i,1),Parameters_CDHW_weibul(i,2));
    Return_Period_wb_all_years(i,j)=1/(1-p);
    end
end



save('nighttime_CDHW_Unusuality_28032024.mat','Return_Period_wb_all_years','CDHW_max_sev');
save('Daytime_CDHW_Unusuality_28032024.mat','Return_Period_wb_all_years','CDHW_max_sev');


%Comparing time series of daytime and nighttime unusuality over different
%states of Southwestern North America


SWNA             =  readtable('/Users/smonda15/Ongoing Work/Publication_Compound_Drought_Heatwave_Over_Phoneix/Datasets/Shapefiles/Southwest_north_America_Grid_Locations.csv');
%SWNA_latlon      =  SWNA(:,1:2);
%SWNA_index       =  SWNA(:,3);
States           =  SWNA(:,19);

index_AZ = find(strcmp(States.NAME, 'Arizona'));
index_UT = find(strcmp(States.NAME, 'Utah'));
index_NM = find(strcmp(States.NAME, 'New Mexico'));
index_CO = find(strcmp(States.NAME, 'Colorado'));
index_MX = find(strcmp(States.NAME, ''));


%Return_Period_d          =  load('Daytime_CDHW_Unusuality.mat');
Return_Period_d          =  Return_Period_wb_all_years;

Return_Period_n          =  load('nighttime_CDHW_Unusuality.mat');

%Return_Period_d          =  Return_Period_d.Return_Period_wb_all_years;
Return_Period_n          =  Return_Period_n.Return_Period_wb_all_years;

%Mexico
Return_Period_d_MX       =  Return_Period_d(index_MX,:);
Return_Period_n_MX       =  Return_Period_n(index_MX,:);

Return_Period_d_MX_med = zeros(44,1);
Return_Period_n_MX_med = zeros(44,1);

for i=1:44
    a = (Return_Period_d_MX(:,i));
    a = a(a>5);

    Return_Period_d_MX_med(i,1) = median(a);
end

Return_Period_d_MX_med1  = movmean(Return_Period_d_MX_med,5);

severity_MX              = CDHW_max_sev(index_MX,:);
severity_MX_med          = (median(severity_MX,1))';

severity_rp_MX           = [movmean(Return_Period_d_MX_med,7),movmean(severity_MX_med,7)];

% Create the plot
figure; % Opens a new figure window
semilogx(movmean(Return_Period_d_MX_med,5), movmean(severity_MX_med,5), 'o-'); % Creates a logarithmic plot with markers and lines
xlabel('Vector 1'); % Label for the x-axis
ylabel('Vector 2'); % Label for the y-axis
title('Logarithmic Plot of Vector 1 vs Vector 2'); % Title for the plot
grid on; % Adds a grid for easier data visualization


for i=1:44
    a = (Return_Period_n_MX(:,i));
    a = a(a>5);

    Return_Period_n_MX_med(i,1) = mean(a);
end


%Arizona
Return_Period_d_AZ       =  Return_Period_d(index_AZ,:);
Return_Period_n_AZ       =  Return_Period_n(index_AZ,:);

Return_Period_d_AZ_med = zeros(44,1);
Return_Period_n_AZ_med = zeros(44,1);

for i=1:44
    a = (Return_Period_d_AZ(:,i));
    a = a(a>2);

    Return_Period_d_AZ_med(i,1) = mode(a);
end

severity_AZ              = CDHW_max_sev(index_AZ,:);
severity_AZ_med          = (median(severity_AZ,1))';

%New Mexico 
Return_Period_d_NM       =  Return_Period_d(index_NM,:);
Return_Period_d_NM_med   = zeros(44,1);

for i=1:44
    a = (Return_Period_d_NM(:,i));
    a = a(a>1.5);

    Return_Period_d_NM_med(i,1) = median(a);
end

severity_NM              = CDHW_max_sev(index_NM,:);
severity_NM_med          = (median(severity_NM,1))';
 
%Using Generalized Extreme Value Distribution 






%Extracting Land-use Land -class corresponding to different locations of
%different states






%Different States

latlon_AZ  = SWNA_latlon(index_AZ,:);
latlon_CO  = SWNA_latlon(index_CO,:);
latlon_MX  = SWNA_latlon(index_MX,:);
latlon_UT  = SWNA_latlon(index_UT,:);
latlon_NM  = SWNA_latlon(index_NM,:);



% Read the raster and its spatial reference
[LULC, R] = readgeoraster('/Users/smonda15/Ongoing Work/Publication_Compound_Drought_Heatwave_Over_Phoneix/Datasets/Shapefiles/NLCD_Reprojected.tif');

% Print the raster reference object info
disp(R);

% Assuming your points are in latitude and longitude,
% check if they fall within the bounding box of the raster.
bbox = [min(R.LongitudeLimits), min(R.LatitudeLimits); 
        max(R.LongitudeLimits), max(R.LatitudeLimits)];
disp('Bounding box of the raster:');
disp(bbox);

% Check the first few points
disp('First few points:');
disp(latlon_AZ(1:min(end,5), :));


% Initialize an array to hold the land use values for each point
landUseModes_AZ = zeros(size(latlon_AZ,1),1);

% Loop through each point in latlon_AZ
for k = 1:size(latlon_AZ,1)
    % Convert each point's longitude and latitude to raster indices
    [col, row] = map2pix(R, latlon_AZ(k,1), latlon_AZ(k,2)); % Ensure latlon_AZ(k,1) is longitude and latlon_AZ(k,2) is latitude

    % Round the indices since they must be integers to access the raster matrix
    col = round(col);
    row = round(row);

    % Check if the computed indices are within the bounds of the raster
    if (row >= 1 && row <= size(LULC, 1) && col >= 1 && col <= size(LULC, 2))
        % Extract the land use value from the raster
        landUseModes_AZ(k) = LULC(row, col);
    else
        % Assign NaN if the point falls outside the raster
        landUseModes_AZ(k) = NaN;
    end
end

% Now, landUseModes_AZ contains the land use code for each point in latlon_AZ













% Sample data setup - replace these with your actual data
years = (1980:2023)'; % Years from 1980 to 2023
% Note: The following random data generation is for illustration purposes.
% You should use your actual data for Return_Period_d_MX_med and Return_Period_n_MX_med.
%Return_Period_d_MX_med = rand(44, 1) * 100; % Example daytime return period data
%Return_Period_n_MX_med = rand(44, 1) * 100; % Example nighttime return period data

% Calculate 5-year running means
windowSize = 5; % Define the window size for the running mean
d_MX_5yr_mean = movmean(Return_Period_d_MX_med, windowSize);
n_MX_5yr_mean = movmean(Return_Period_n_MX_med, windowSize);

% Begin Visualization
figure;

% Plotting both original time series and their 5-year running means
plot(years, Return_Period_d_MX_med, '--b', 'LineWidth', 1); % Original daytime data in blue
hold on;
plot(years, Return_Period_n_MX_med, '--r', 'LineWidth', 1); % Original nighttime data in red

% Plotting 5-year running means
plot(years, d_MX_5yr_mean, '-b', 'LineWidth', 2); % 5-year running mean for daytime in dashed blue
plot(years, n_MX_5yr_mean, '-r', 'LineWidth', 2); % 5-year running mean for nighttime in dashed red

% Enhancements for publication quality
title('Comparative Time Series of Return Periods', 'FontSize', 16);
xlabel('Year', 'FontSize', 14);
ylabel('Return Period (years)', 'FontSize', 14);
legend({'Daytime', 'Nighttime', '5-year Mean Daytime', '5-year Mean Nighttime'}, 'FontSize', 12, 'Location', 'best');
set(gca, 'FontSize', 12); % Sets the default text sizes
xlim([1980 2023]); % Set x-axis limits to cover the full range of years
grid on; % Add grid for better readability

% Save the Figure with High Resolution
print('TimeSeriesComparisonWithRunningMean','-dpng','-r300'); % Saves the figure as a PNG file with 300 DPI


RP = readmatrix('/Users/smonda15/Ongoing Work/Publication_Compound_Drought_Heatwave_Over_Phoneix/Results/Results Southwest North America/SouthwestNA_2023_CDHW_longdrReturn_Period.csv');

RPday = RP(:,3);
RPnight = RP(:,4);

a = find(isnan(RPnight));


%% Plotting cumuliative anomaly of Max, Min Temperature and Precipitation

%Making an attempt to generate the time series image with the same aspect
%ratio

img = imread('/Users/smonda15/Ongoing Work/Publication_Compound_Drought_Heatwave_Over_Phoneix/Results/Results Southwest North America/Plots/SouthwestNA_Tmin_JUL_AUG_2023 copy 3.png');
[height, width, ~] = size(img);
aspectRatio = width / height;
disp(['Aspect Ratio: ', num2str(aspectRatio)]);



% Sample data (replace with your actual data)
vector1 = cumMean_2023; % Replace with your first vector
vector2 = cumMean_2023; % Replace with your second vector
vector3 = SW_pcp_anm_2023_cum; % Replace with your third vector

% Load the image to get its dimensions
imgInfo = imfinfo('/Users/smonda15/Ongoing Work/Publication_Compound_Drought_Heatwave_Over_Phoneix/Results/Results Southwest North America/Plots/SouthwestNA_Tmax_JUL_AUG_2023 (2).png'); % Replace 'example.png' with your image file name
imageWidth = imgInfo.Width;
imageHeight = imgInfo.Height;

% Set the figure size
figure('Position', [100, 100, imageWidth, imageHeight]);
hold on; % Keep all plots on the same figure

% Time or X-axis values
time = 1:183; % Adjust if you have specific time or x-axis values

% Plotting the vectors
% Replace 'vector1', 'vector2', 'vector3' with your actual data vectors
plot(time, vector1, 'LineWidth', 2, 'DisplayName', 'Series 1');
plot(time, vector2, 'LineWidth', 2, 'DisplayName', 'Series 2');
plot(time, vector3, 'LineWidth', 2, 'DisplayName', 'Series 3');

% Enhancing the plot
xlabel('Time', 'FontSize', 14); % Replace 'Time' with your actual x-axis label
ylabel('Cumulative Anomaly', 'FontSize', 14); % Adjust label accordingly
legend('show', 'Location', 'best'); % Show legend
grid on; % Add grid for better readability
box on; % Enclose plot in a box

% Adjusting axes for better visibility
ax = gca; % Current axes
ax.FontSize = 12; % Adjust font size
ax.LineWidth = 1.5; % Adjust line width for the axes

% Save the figure with high resolution
saveas(gcf, 'HighQualityPlot.png', 'png'); % Saves the figure as a PNG file



%Dual Axis Plots

% Time or X-axis values
time = 1:92; % Adjust if you have specific time or x-axis values

% Create figure
figure;
hold on; % Keep all plots on the same figure

% Plotting the first two vectors on the left y-axis
yyaxis left; % Activate left y-axis
plot(time, vector1, 'b-', 'LineWidth', 2, 'DisplayName', 'Series 1'); % Blue line
plot(time, vector2, 'g-', 'LineWidth', 2, 'DisplayName', 'Series 2'); % Green line
ylabel('Left Axis Label', 'FontSize', 14); % Adjust left y-axis label

% Plotting the third vector on the right y-axis
yyaxis right; % Activate right y-axis
plot(time, vector3, 'r-', 'LineWidth', 2, 'DisplayName', 'Series 3'); % Red line
ylabel('Right Axis Label', 'FontSize', 14); % Adjust right y-axis label

% Enhancing the plot
xlabel('Time', 'FontSize', 14); % Replace 'Time' with your actual x-axis label
title('Dual Y-Axis Time Series Visualization', 'FontSize', 16); % Adjust title
legend('show', 'Location', 'best'); % Show legend
grid on; % Add grid for better readability
box on; % Enclose plot in a box

% Adjusting axes for better visibility
ax = gca; % Current axes
ax.FontSize = 12; % Adjust font size
ax.LineWidth = 1.5; % Adjust line width for the axes

% Save the figure with high resolution
saveas(gcf, 'DualAxisHighQualityPlot.png', 'png'); % Saves the figure as a PNG file




%Using matlab to place images properly:

% Read the images
img1 = imread('/Users/smonda15/Ongoing Work/Publication_Compound_Drought_Heatwave_Over_Phoneix/Results/Results Southwest North America/Plots/SouthwestNA_PCP_JUL_AUG_2023 copy.png');
img2 = imread('/Users/smonda15/Ongoing Work/Publication_Compound_Drought_Heatwave_Over_Phoneix/Results/Results Southwest North America/Plots/SouthwestNA_Tmax_JUL_AUG_2023 copy 3.png');
img3 = imread('/Users/smonda15/Ongoing Work/Publication_Compound_Drought_Heatwave_Over_Phoneix/Results/Results Southwest North America/Plots/SouthwestNA_Tmin_JUL_AUG_2023 copy 2.png');
img4 = imread('/Users/smonda15/Ongoing Work/Publication_Compound_Drought_Heatwave_Over_Phoneix/Results/Results Southwest North America/Plots/Cum_Anm_July_Sep_2023.png');




% Create figure
figure;

% Subplot 1
subplot(2,2,1);
imshow(img1);
title('a', 'FontSize', 14, 'FontWeight', 'bold'); % Label for the first image
% Adding a text box
annotation('textbox', [0.15, 0.83, 0.1, 0.1], 'String', 'T_m_a_x', ...
    'FitBoxToText', 'on', 'LineStyle', 'none', 'FontSize', 12);

% Subplot 2
subplot(2,2,2);
imshow(img2);
title('b', 'FontSize', 14, 'FontWeight', 'bold'); % Label for the second image
% Adding a text box
annotation('textbox', [0.55, 0.83, 0.1, 0.1], 'String', 'T_m_i_n', ...
    'FitBoxToText', 'on', 'LineStyle', 'none', 'FontSize', 12);

% Subplot 3
subplot(2,2,3);
imshow(img3);
title('c', 'FontSize', 14, 'FontWeight', 'bold'); % Label for the third image
% Adding a text box
annotation('textbox', [0.15, 0.38, 0.1, 0.1], 'String', 'P', ...
    'FitBoxToText', 'on', 'LineStyle', 'none', 'FontSize', 12);

% Subplot 4
subplot(2,2,4);
imshow(img4);
title('d', 'FontSize', 14, 'FontWeight', 'bold'); % Label for the fourth image

% Adjust subplot spacing
set(gcf, 'Position', [100, 100, 800, 600]); % Adjust figure size


%Trying the eps plots in matlab



% Prepare the figure
figure;
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 10, 8]); % Set the figure size

% Subplot positions (left, bottom, width, height) in normalized units
% Adjust these values to reduce space between subplots and ensure equal size
positions = [
    0.1, 0.6, 0.43, 0.39;  % Top-left
    0.54, 0.6, 0.43, 0.39; % Top-right
    0.1, 0.11, 0.43, 0.39; % Bottom-left
    0.54, 0.11, 0.43, 0.39; % Bottom-right
];


% File names (assuming you have converted your images from EPS to PNG)
files = {
    '/Users/smonda15/Ongoing Work/Publication_Compound_Drought_Heatwave_Over_Phoneix/Results/Results Southwest North America/Plots/SouthwestNA_Tmax_JUL_AUG_2023 copy 4.png',
    '/Users/smonda15/Ongoing Work/Publication_Compound_Drought_Heatwave_Over_Phoneix/Results/Results Southwest North America/Plots/SouthwestNA_Tmin_JUL_AUG_2023 copy 3.png',
    '/Users/smonda15/Ongoing Work/Publication_Compound_Drought_Heatwave_Over_Phoneix/Results/Results Southwest North America/Plots/SouthwestNA_PCP_JUL_AUG_2023 copy 2.png',
    '/Users/smonda15/Ongoing Work/Publication_Compound_Drought_Heatwave_Over_Phoneix/Results/Results Southwest North America/Plots/Cum_Anm_Apr_Sep1_2023.png'

};

% Labels for each subplot
labels = {'a', 'b', 'c', 'd'};

% Loop through each subplot
for i = 1:4
    % Create subplot
    subplot('Position', positions(i, :));
    % Read and show the image
    img = imread(files{i}); % Corrected: use curly braces for cell array indexing
% Read image file
    imshow(img);
    set(gca, 'Position', positions(i, :)); % Adjust position again after imshow
    
    % Remove axis
    axis off;
    
    % Add labels directly to images to avoid overlap with tight layout
    text(20, 20, labels{i}, 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'white'); % Adjust text position as needed
    
    % Add text boxes to the first three plots
    if i < 4
        dim = [positions(i,1), positions(i,2) - 0.1, 0.1, 0.05]; % Position for annotation
        annotation('textbox', dim, 'String', 'Your Text Here', ...
            'FitBoxToText', 'on', 'LineStyle', 'none', 'FontSize', 12, ...
            'BackgroundColor', 'white', 'Color', 'black');
    end
end

% Adjust figure and print settings for saving
set(gcf, 'PaperPositionMode', 'auto');
print('-depsc2', 'OutputFigure.eps'); % Save the figure as an EPS file






%Extracting Grid-cells having negative precipitation anomaly

PCP_anm             =    readmatrix('/Users/smonda15/Ongoing Work/Publication_Compound_Drought_Heatwave_Over_Phoneix/Results/Results Southwest North America/Precipitation_Anomaly_Temporally_Averaged.csv');

PCP_anm_July_Sep    =    PCP_anm(:,6);

PCP_latlon          =    PCP_anm(:,1:2);

index_PCP           =    find(PCP_anm_July_Sep<0);

index_PCP_latlon    =    PCP_latlon(index_PCP,:);

index_PCP_anm       =    PCP_anm_July_Sep(index_PCP,1);

T_anm               =    readmatrix('/Users/smonda15/Ongoing Work/Publication_Compound_Drought_Heatwave_Over_Phoneix/Results/Results Southwest North America/Temperature_Anomaly_Temporally_Averaged.csv');
Tmax_anm_Jul_Sep    =    T_anm(:,6);

index_Tmax          =    find(Tmax_anm_Jul_Sep>1.5);
index_Tmax_latlon   =    PCP_latlon(index_Tmax,:);
index_Tmax_anm      =    Tmax_anm_Jul_Sep(index_Tmax,1);




%Plotting Population against Precipitation anomalies

PCP_pop             =    readmatrix('/Users/smonda15/Ongoing Work/Publication_Compound_Drought_Heatwave_Over_Phoneix/Results/Results Southwest North America/Precipitation_negative_Anm_Population.csv');
PCP_population      =    PCP_pop(:,4)  









%Plotting of Figure 3:

dataSheet1    = readmatrix('/Users/smonda15/Ongoing Work/Publication_Compound_Drought_Heatwave_Over_Phoneix/Results/Results Southwest North America/Statewide_Severity_CDHW.xlsx', 'Sheet', 'Spatial_Extent_CDHW'); % By sheet name
dataSheet2    = readmatrix('/Users/smonda15/Ongoing Work/Publication_Compound_Drought_Heatwave_Over_Phoneix/Results/Results Southwest North America/Statewide_Severity_CDHW_night.xlsx', 'Sheet', 'Spatial_Extent_CDHW'); % By sheet index


datasheet     = [dataSheet1,dataSheet2];

AZ_severity   = [dataSheet1(:,2), dataSheet2(:,2)];

% Assuming vector1 and vector2 are your data vectors
vector1 = AZ_severity(:,1); % Replace with your actual data
vector2 = AZ_severity(:,2); % Replace with your actual data

%Getting the dimensions of the image
% Read the image
img = imread('/Users/smonda15/Ongoing Work/Publication_Compound_Drought_Heatwave_Over_Phoneix/Results/Results Southwest North America/Plots/Figure 3_Incomplete.png'); % Replace 'yourimage.jpg' with your image file

% Get the size of the image
[height, width, depth] = size(img);


% Create a new figure. You can adjust the size to match another figure
% For example, if the other figure is 560x420 pixels:
figure('Position', [100, 100, 2475, 1024]); 

% Plot both vectors
hold on; % Allows multiple plots in the same figure
plot(vector1, 'LineWidth', 2, 'DisplayName', 'Vector 1'); % Adjust line properties as needed
plot(vector2, 'LineWidth', 2, 'DisplayName', 'Vector 2'); % Adjust line properties as needed
hold off;

% Enhancing the plot
xlabel('X-axis Label', 'FontSize', 14); % Replace X-axis Label with your label
ylabel('Y-axis Label', 'FontSize', 14); % Replace Y-axis Label with your label
title('Your Title Here', 'FontSize', 16); % Replace Your Title Here with your title
legend('show', 'Location', 'best'); % Shows the legend
grid on; % Adds a grid
set(gca, 'FontSize', 12); % Adjust font size as needed

% Setting axis limits if necessary
xlim([1 44]); % Assuming your vectors have indices from 1 to 26
% ylim([min_value max_value]); % Uncomment and replace min_value and max_value with actual values if needed

% Make the plot look better by setting the box on and making the line thicker
set(gca, 'Box', 'on', 'LineWidth', 2);

% If you want to match this plot's size to another figure's size exactly,
% you can set the 'Position' property of the figure to the same values.
% You can find the position of an existing figure by getting its 'Position' attribute:
% pos = get(other_figure_handle, 'Position');
% and then setting this figure's position to 'pos'.



% Assuming vector1 and vector2 are your data vectors
vector1 = AZ_severity(:,1); % Replace with your actual data
vector2 = AZ_severity(:,2); % Replace with your actual data

% Calculate 5-year running means
windowSize = 5; % 5 years for the running mean
movAvg1 = movmean(vector1, windowSize);
movAvg2 = movmean(vector2, windowSize);

% Calculate trend lines (linear regression)
x = (1:length(vector1))'; % Independent variable (time)
p1 = polyfit(x, vector1, 1); % p1 contains the slope and intercept for vector1
p2 = polyfit(x, vector2, 1); % p2 contains the slope and intercept for vector2
trendLine1 = polyval(p1, x);
trendLine2 = polyval(p2, x);

% Read the image for dimensions
img = imread('/Users/smonda15/Ongoing Work/Publication_Compound_Drought_Heatwave_Over_Phoneix/Results/Results Southwest North America/Plots/Figure 3_Incomplete.png');

% Get the size of the image
[height, width, depth] = size(img);

% Create a new figure
figure('Position', [100, 100, 2475, 1024]);

% Plot both vectors and their 5-year running means and trend lines
hold on;
plot(vector1, 'LineWidth', 2, 'DisplayName', 'Vector 1');
plot(vector2, 'LineWidth', 2, 'DisplayName', 'Vector 2');
plot(movAvg1, 'LineWidth', 2, 'DisplayName', '5-Year Running Mean (Vector 1)', 'Color', 'cyan');
plot(movAvg2, 'LineWidth', 2, 'DisplayName', '5-Year Running Mean (Vector 2)', 'Color', 'magenta');
plot(x, trendLine1, '--', 'LineWidth', 2, 'DisplayName', 'Trend Line (Vector 1)', 'Color', 'blue');
plot(x, trendLine2, '--', 'LineWidth', 2, 'DisplayName', 'Trend Line (Vector 2)', 'Color', 'red');
hold off;

% Enhancing the plot
xlabel('X-axis Label', 'FontSize', 14); % Replace X-axis Label with your label
ylabel('Y-axis Label', 'FontSize', 14); % Replace Y-axis Label with your label
title('Your Title Here', 'FontSize', 16); % Replace Your Title Here with your title
legend('show', 'Location', 'best'); % Shows the legend
grid on; % Adds a grid
set(gca, 'FontSize', 12);

% Setting axis limits
xlim([1 length(vector1)]); % Assuming your vectors match the full range

% Make the plot look better
set(gca, 'Box', 'on', 'LineWidth', 2);



%Changing the dimension of the Figure File

openfig('/Users/smonda15/Ongoing Work/Publication_Compound_Drought_Heatwave_Over_Phoneix/Results/Results Southwest North America/Plots/Daytime_CDHW_model_Selection.fig'); % Replace 'your_figure.fig' with the path to your file
fig = gcf; % Get the handle to the current figure
fig.Position = [fig.Position(1) fig.Position(2) 1024 1024]; % Change the size while keeping the current position






%Analyzing relationship between severity and return period

SWNA                =  readtable('/Users/smonda15/Ongoing Work/Publication_Compound_Drought_Heatwave_Over_Phoneix/Datasets/Shapefiles/Southwest_north_America_Grid_Locations.csv');
%SWNA_latlon      =  SWNA(:,1:2);
%SWNA_index       =  SWNA(:,3);
States              =  SWNA(:,19);


index_AZ = find(strcmp(States.NAME, 'Arizona'));
index_UT = find(strcmp(States.NAME, 'Utah'));
index_NM = find(strcmp(States.NAME, 'New Mexico'));
index_CO = find(strcmp(States.NAME, 'Colorado'));
index_MX = find(strcmp(States.NAME, ''));



LULC                  = readmatrix('/Users/smonda15/Ongoing Work/Publication_Compound_Drought_Heatwave_Over_Phoneix/Results/Results Southwest North America/LULC_output_values1.csv');

nonnan_indices_LULC   = find(~isnan(LULC(:,3)));
LULC_nonnan           = LULC(nonnan_indices_LULC,3);

nonzero_indices_LULC  = find(LULC_nonnan~=0);
LULC_nozero           = LULC_nonnan(nonzero_indices_LULC,1);

latlon_nonnan         = LULC(nonnan_indices_LULC,1:2);
latlon_nonzero        = latlon_nonnan(nonzero_indices_LULC,:);



LULC_latlon          = readtable('/Users/smonda15/Ongoing Work/Publication_Compound_Drought_Heatwave_Over_Phoneix/Datasets/Shapefiles/4-corner_states_LULC_latlon.csv');
indices_4corner      = LULC_latlon(:,3);
indices_4corner1      = indices_4corner.Values ;


LULC_4_corner        = LULC_nonzero(indices_4corner1);

LULC_unique          = unique(LULC_4_corner);

CDHW_max_sev_day     = CDHW_max_sev(nonnan_indices_LULC,:);
CDHW_max_sev_day     = CDHW_max_sev_day(nonzero_indices_LULC,:);
CDHW_max_sev_day     = CDHW_max_sev_day(indices_4corner1,:);

Return_Period_day    = Return_Period_wb_all_years(nonnan_indices_LULC,:);
Return_Period_day    = Return_Period_day(nonzero_indices_LULC,:);
Return_Period_day    = Return_Period_day(indices_4corner1,:);


CDHW_max_sev_night     = CDHW_max_sev(nonnan_indices_LULC,:);
CDHW_max_sev_night     = CDHW_max_sev_night(nonzero_indices_LULC,:);
CDHW_max_sev_night     = CDHW_max_sev_night(indices_4corner1,:);

Return_Period_night   = Return_Period_wb_all_years(nonnan_indices_LULC,:);
Return_Period_night    = Return_Period_night(nonzero_indices_LULC,:);
Return_Period_night    = Return_Period_night(indices_4corner1,:);


%Extracting indices corresponding to different LULC

indices_Forest                 = find(LULC_4_corner == 41 | LULC_4_corner == 42 | LULC_4_corner == 43);
indices_developed              = find(LULC_4_corner == 21 | LULC_4_corner == 22 | LULC_4_corner == 23| LULC_4_corner==24);
indices_shurub                 = find(LULC_4_corner ==52);
indices_grassland              = find(LULC_4_corner == 71);
indices_agriculture            = find(LULC_4_corner ==81| LULC_4_corner==82);
indices_wetlands               = find(LULC_4_corner == 90);


CDHW_max_sev_day_forest        = CDHW_max_sev_day(indices_Forest,:);
CDHW_max_sev_night_forest      = CDHW_max_sev_night(indices_Forest,:);
CDHW_max_sev_day_Forest_av     = (mean(CDHW_max_sev_day_forest,1))';
CDHW_max_sev_night_Forest_av   = (mean(CDHW_max_sev_night_forest,1))';


Return_Period_day_forest       = Return_Period_day(indices_Forest,:);
Return_Period_night_forest     = Return_Period_night(indices_Forest,:);
Return_Period_day_forest_av    = (mean(Return_Period_day_forest,1))';
Return_Period_night_forest_av  = (mean(Return_Period_night_forest,1))';



CDHW_max_sev_day_developed        = CDHW_max_sev_day(indices_developed,:);
CDHW_max_sev_night_developed      = CDHW_max_sev_night(indices_developed,:);
CDHW_max_sev_day_developed_av     = (mean(CDHW_max_sev_day_developed,1))';
CDHW_max_sev_night_developed_av   = (mean(CDHW_max_sev_night_developed,1))';


Return_Period_day_developed       = Return_Period_day(indices_developed,:);
Return_Period_night_developed     = Return_Period_night(indices_developed,:);
Return_Period_day_developed_av    = (mean(Return_Period_day_developed,1))';
Return_Period_night_developed_av  = (mean(Return_Period_night_developed,1))';


CDHW_max_sev_day_shurub        = CDHW_max_sev_day(indices_shurub,:);
CDHW_max_sev_night_shurub      = CDHW_max_sev_night(indices_shurub,:);
CDHW_max_sev_day_shurub_av     = (mean(CDHW_max_sev_day_shurub,1))';
CDHW_max_sev_night_shurub_av   = (mean(CDHW_max_sev_night_shurub,1))';


Return_Period_day_shurub       = Return_Period_day(indices_shurub,:);
Return_Period_night_shurub     = Return_Period_night(indices_shurub,:);
Return_Period_day_shurub_av    = (mean(Return_Period_day_shurub,1))';
Return_Period_night_shurub_av  = (mean(Return_Period_night_shurub,1))';


CDHW_max_sev_day_grassland        = CDHW_max_sev_day(indices_grassland,:);
CDHW_max_sev_night_grassland      = CDHW_max_sev_night(indices_grassland,:);
CDHW_max_sev_day_grassland_av     = (mean(CDHW_max_sev_day_grassland,1))';
CDHW_max_sev_night_grassland_av   = (mean(CDHW_max_sev_night_grassland,1))';


Return_Period_day_grassland       = Return_Period_day(indices_grassland,:);
Return_Period_night_grassland     = Return_Period_night(indices_grassland,:);
Return_Period_day_grassland_av    = (mean(Return_Period_day_grassland,1))';
Return_Period_night_grassland_av  = (mean(Return_Period_night_grassland,1))';


CDHW_max_sev_day_agriculture        = CDHW_max_sev_day(indices_agriculture,:);
CDHW_max_sev_night_agriculture      = CDHW_max_sev_night(indices_agriculture,:);
CDHW_max_sev_day_agriculture_av     = (mean(CDHW_max_sev_day_agriculture,1))';
CDHW_max_sev_night_agriculture_av   = (mean(CDHW_max_sev_night_agriculture,1))';


Return_Period_day_agriculture       = Return_Period_day(indices_agriculture,:);
Return_Period_night_agriculture     = Return_Period_night(indices_agriculture,:);
Return_Period_day_agriculture_av    = (mean(Return_Period_day_agriculture,1))';
Return_Period_night_agriculture_av  = (mean(Return_Period_night_agriculture,1))';


CDHW_max_sev_day_wetlands        = CDHW_max_sev_day(indices_wetlands,:);
CDHW_max_sev_night_wetlands      = CDHW_max_sev_night(indices_wetlands,:);
CDHW_max_sev_day_wetlands_av     = (mean(CDHW_max_sev_day_wetlands,1))';
CDHW_max_sev_night_wetlands_av   = (mean(CDHW_max_sev_night_wetlands,1))';


Return_Period_day_wetlands       = Return_Period_day(indices_wetlands,:);
Return_Period_night_wetlands     = Return_Period_night(indices_wetlands,:);
Return_Period_day_wetlands_av    = (mean(Return_Period_day_wetlands,1))';
Return_Period_night_wetlands_av  = (mean(Return_Period_night_wetlands,1))';



SWNA                =  readtable('/Users/smonda15/Ongoing Work/Publication_Compound_Drought_Heatwave_Over_Phoneix/Datasets/Shapefiles/Southwest_north_America_Grid_Locations.csv');
%SWNA_latlon      =  SWNA(:,1:2);
%SWNA_index       =  SWNA(:,3);
States              =  SWNA(:,19);


index_AZ = find(strcmp(States.NAME, 'Arizona'));
index_UT = find(strcmp(States.NAME, 'Utah'));
index_NM = find(strcmp(States.NAME, 'New Mexico'));
index_CO = find(strcmp(States.NAME, 'Colorado'));
index_MX = find(strcmp(States.NAME, ''));


load('Daytime_CDHW_Unusuality_28032024.mat');
CHDW_max_sev_day               =    CDHW_max_sev;
Return_Period_day              =    Return_Period_wb_all_years;

load('nighttime_CDHW_Unusuality_28032024.mat');
CHDW_max_sev_night              =    CDHW_max_sev;
Return_Period_night              =    Return_Period_wb_all_years;


CDHW_max_sev_day_AZ                =    CHDW_max_sev_day(index_AZ,:);
CDHW_max_sev_night_AZ              =    CHDW_max_sev_night(index_AZ,:);
Return_Period_day_AZ               =    Return_Period_day(index_AZ,:);
Return_Period_night_AZ             =    Return_Period_night(index_AZ,:);

CDHW_max_sev_day_AZ_mean           = (mean(CDHW_max_sev_day_AZ,1))';
Return_Period_day_AZ_mean          = (mean(Return_Period_day_AZ ,1))';
CDHW_max_sev_night_AZ_mean         = (mean(CDHW_max_sev_night_AZ,1))';
Return_Period_night_AZ_mean        = (mean(Return_Period_night_AZ ,1))';

CDHW_max_sev_day_night_AZ          = [CDHW_max_sev_day_AZ_mean,CDHW_max_sev_night_AZ_mean];
Return_Period_day_night_AZ         = [Return_Period_day_AZ_mean,Return_Period_night_AZ_mean];

createfigure(CDHW_max_sev_day_night_AZ);
createfigure(Return_Period_day_night_AZ);

CDHW_max_sev_day_MX                =    CHDW_max_sev_day(index_MX,:);
CDHW_max_sev_night_MX              =    CHDW_max_sev_night(index_MX,:);
Return_Period_day_MX               =    Return_Period_day(index_MX,:);
Return_Period_night_MX             =    Return_Period_night(index_MX,:);

CDHW_max_sev_day_MX_mean           = (mean(CDHW_max_sev_day_MX,1))';
Return_Period_day_MX_mean          = (mean(Return_Period_day_MX ,1))';
CDHW_max_sev_night_MX_mean         = (mean(CDHW_max_sev_night_MX,1))';
Return_Period_night_MX_mean        = (mean(Return_Period_night_MX ,1))';

CDHW_max_sev_day_night_MX          = [CDHW_max_sev_day_MX_mean,CDHW_max_sev_night_MX_mean];
Return_Period_day_night_MX         = [Return_Period_day_MX_mean,Return_Period_night_MX_mean];

createfigure(CDHW_max_sev_day_night_MX);
createfigure(Return_Period_day_night_MX);


CDHW_max_sev_day_NM                =    CHDW_max_sev_day(index_NM,:);
CDHW_max_sev_night_NM              =    CHDW_max_sev_night(index_NM,:);
Return_Period_day_NM               =    Return_Period_day(index_NM,:);
Return_Period_night_NM             =    Return_Period_night(index_NM,:);

CDHW_max_sev_day_NM_mean           = (mean(CDHW_max_sev_day_NM,1))';
Return_Period_day_NM_mean          = (mean(Return_Period_day_NM ,1))';
CDHW_max_sev_night_NM_mean         = (mean(CDHW_max_sev_night_NM,1))';
Return_Period_night_NM_mean        = (mean(Return_Period_night_NM ,1))';

CDHW_max_sev_day_night_NM          = [CDHW_max_sev_day_NM_mean,CDHW_max_sev_night_NM_mean];
Return_Period_day_night_NM         = [Return_Period_day_NM_mean,Return_Period_night_NM_mean];

createfigure(CDHW_max_sev_day_night_NM);
createfigure(Return_Period_day_night_NM);

loglog(Return_Period_wb_all_years_AZ_mean, CDHW_max_sev_AZ_mean , 'o', 'MarkerSize', 8, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');


b = [CDHW_max_sev_AZ_mean,Return_Period_wb_all_years_AZ_mean];

CDHW_max_sev_MX                =    CDHW_max_sev(index_MX,:);

Return_Period_wb_all_years_MX  =    Return_Period_wb_all_years(index_MX,:);

CDHW_max_sev_MX_mean                = (mean(CDHW_max_sev_MX,1))';
Return_Period_wb_all_years_MX_mean  = (mean(Return_Period_wb_all_years_MX,1))';

% Your data vectors
x = 1980:2023;  % This is the x-axis common to both y-axes
y1 = CDHW_max_sev_MX_mean;  % First set of data
y2 = Return_Period_wb_all_years_MX_mean;  % Second set of data, which we'll plot against a second y-axis

% Create the figure
figure;

% Left y-axis
yyaxis left;
plot(x, y1, '-b', 'LineWidth', 2);  % '-b' makes the line blue
ylabel('Left Y-Axis Data');

% Right y-axis
yyaxis right;
plot(x, y2, '-r', 'LineWidth', 2);  % '-r' makes the line red
ylabel('Right Y-Axis Data');

% X-axis and title
xlabel('Common X-Axis');
title('Dual Y-Axis Plot');

% Customize axes
xlim([min(x) max(x)]);  % Set the x-axis limits to the range of x
yyaxis left;
ylim([0 max(y1)]);  % Set the y-axis limits for left y-axis to the range of y1
yyaxis right;
ylim([0 max(y2)]);  % Set the y-axis limits for right y-axis to the range of y2

% Adding grid and legend
grid on;
legend('Data Y1', 'Data Y2');




%Plotting both day and night severity and return period
% Assuming 'day_severity', 'night_severity', 'day_return_period', 'night_return_period' are your vectors
% x-axis is the common axis representing years from 1980 to 2023

x = 1980:2023;
day_severity = CDHW_max_sev_day_MX_mean ; % your day severity data
night_severity = CDHW_max_sev_night_MX_mean; % your night severity data
day_return_period = Return_Period_day_MX_mean; % your day return period data
night_return_period = Return_Period_night_MX_mean; % your night return period data

% Create the figure
figure;

% Left y-axis for severity
yyaxis left;
% Plotting day severity with a solid line and night severity with a dashed line
plot(x, day_severity, '-b', 'LineWidth', 2);
hold on; % Hold on to plot on the same figure
plot(x, night_severity, '--b', 'LineWidth', 2);

ylabel('Severity');
% Set the y-axis limits for left y-axis to the range of severity data
ylim([min([day_severity, night_severity]) max([day_severity, night_severity])]);
yticks(linspace(min([day_severity, night_severity]), max([day_severity, night_severity]), 5)); % 5 y-ticks for clear reading

% Right y-axis for return period
yyaxis right;
% Plotting day return period with a solid line and night return period with a dashed line
plot(x, day_return_period, '-r', 'LineWidth', 2);
plot(x, night_return_period, '--r', 'LineWidth', 2);

ylabel('Return Period (years)');
% Set the y-axis limits for right y-axis to the range of return period data
ylim([min([day_return_period, night_return_period]) max([day_return_period, night_return_period])]);
yticks(linspace(min([day_return_period, night_return_period]), max([day_return_period, night_return_period]), 5)); % 5 y-ticks

% X-axis and title
xlabel('Year');
title('Severity and Return Period Over Years');

% Customizing the tick marks for better readability
xticks(min(x):5:max(x)); % Tick marks every 5 years

% Legend
legend({'Day Severity', 'Night Severity', 'Day Return Period', 'Night Return Period'}, 'Location', 'northwest');

% Customize grid
grid on;
set(gca, 'GridLineStyle', '--'); % Dashed grid lines

% Set the figure size for better aspect ratio
set(gcf, 'Position', [100, 100, 800, 600]);

% Enhance line visibility
ax = gca;
ax.LineWidth = 1.5;

% Save the figure with high resolution
saveas(gcf, 'Severity_Return_Period.png', 'png');



x = 1980:2023;
day_severity = CDHW_max_sev_day_MX_mean; % your day severity data
night_severity = CDHW_max_sev_night_MX_mean; % your night severity data
day_return_period = Return_Period_day_MX_mean; % your day return period data
night_return_period = Return_Period_night_MX_mean; % your night return period data

% Create the figure
figure;

% Left y-axis for severity
yyaxis left;
plot(x, day_severity, '-b', 'LineWidth', 2);
hold on;
plot(x, night_severity, '--b', 'LineWidth', 2);

ylabel('Severity');

% Define limits for y-axis (ensure they are not equal and are valid)
left_y_lim = [min([day_severity; night_severity]), max([day_severity; night_severity])];
if left_y_lim(1) == left_y_lim(2)
    left_y_lim = left_y_lim + [-1, 1]; % Add a small buffer if limits are equal
end
ylim(left_y_lim);

% Right y-axis for return period
yyaxis right;
plot(x, day_return_period, '-r', 'LineWidth', 2);
plot(x, night_return_period, '--r', 'LineWidth', 2);

ylabel('Return Period (years)');

% Define limits for y-axis (ensure they are not equal and are valid)
right_y_lim = [min([day_return_period; night_return_period]), max([day_return_period; night_return_period])];
if right_y_lim(1) == right_y_lim(2)
    right_y_lim = right_y_lim + [-1, 1]; % Add a small buffer if limits are equal
end
ylim(right_y_lim);

% X-axis and title
xlabel('Year');
title('Severity and Return Period Over Years');

% Customizing the tick marks for better readability
xticks(min(x):5:max(x)); % Tick marks every 5 years

% Legend
legend({'Day Severity', 'Night Severity', 'Day Return Period', 'Night Return Period'}, 'Location', 'northwest');

% Customize grid
grid on;
set(gca, 'GridLineStyle', '--'); % Dashed grid lines

% Set the figure size for better aspect ratio
set(gcf, 'Position', [100, 100, 800, 600]);

% Enhance line visibility
ax = gca;
ax.LineWidth = 1.5;

% Save the figure with high resolution
saveas(gcf, 'Severity_Return_Period.png', 'png');



%Trying Stacked plot, top one for severity and bottom one of return period


x = 1980:2023;
day_severity = CDHW_max_sev_day_MX_mean; % your day severity data
night_severity = CDHW_max_sev_night_MX_mean; % your night severity data
day_return_period = Return_Period_day_MX_mean; % your day return period data
night_return_period = Return_Period_night_MX_mean; % your night return period data

% Create the figure
figure;

% Subplot for severity
subplot(2, 1, 1); % This means 2 rows, 1 column, first subplot
hold on;
plot(x, day_severity, '-b', 'LineWidth', 2, 'DisplayName', 'Day Severity');
plot(x, night_severity, '--b', 'LineWidth', 2, 'DisplayName', 'Night Severity');
ylabel('Severity');
xlabel('Year');
title('Severity Over Years');
legend('Location', 'northwest');
grid on;
hold off;

% Subplot for return period
subplot(2, 1, 2); % Second subplot
hold on;
plot(x, day_return_period, '-r', 'LineWidth', 2, 'DisplayName', 'Day Return Period');
plot(x, night_return_period, '--r', 'LineWidth', 2, 'DisplayName', 'Night Return Period');
ylabel('Return Period (years)');
xlabel('Year');
title('Return Period Over Years');
legend('Location', 'northwest');
grid on;
hold off;

% Adjust the layout to make sure there's no overlap and titles/labels are clearly visible
set(gcf, 'Position', [100, 100, 800, 600]);




%Exploring Relationship between HW amplification and Elevation

night_ampl    =    (night_ampl-1)*100;
day_ampl      =    (day_ampl-1)*100;

index_600     =    find(elev<600);
index_1200    =    find(elev<1600 & elev>800);
index_2400    =    find(elev<2400 & elev>1600);
index_3500    =    find(elev<3500 & elev>2400);

nigh_600      =    night_ampl(index_600,:);
nigh_1200     =    night_ampl(index_1200,:);
nigh_2400     =    night_ampl(index_2400,:);
nigh_3500     =    night_ampl(index_3500,:);


day_600       =    day_ampl(index_600,:);
day_1200      =    day_ampl(index_1200,:);
day_2400      =    day_ampl(index_2400,:);
day_3500      =    day_ampl(index_3500,:);








