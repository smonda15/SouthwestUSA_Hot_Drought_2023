
%Considering Utah, Colorado, Idaho, Nevada, California, Oregon, Texas,
%Arizona, Oklahoma, Nebraska, Wyoming, kansas and Mexico as Southwest North America


lon = ncread("/Users/smonda15/Ongoing Work/Publication_Compound_Drought_Heatwave_Over_Phoneix/Datasets/CPC Global Temperature/tmax.1979.nc",'lon');
lat = ncread("/Users/smonda15/Ongoing Work/Publication_Compound_Drought_Heatwave_Over_Phoneix/Datasets/CPC Global Temperature/tmax.1979.nc",'lat');

folderPath = '/Users/smonda15/Ongoing Work/Publication_Compound_Drought_Heatwave_Over_Phoneix/Datasets/CPC Global Temperature';
ncFiles = dir(fullfile(folderPath, '*.nc'));


b=zeros(259200,1);
d=zeros(259200,1);
c=1979;

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

%Extracting Long term mean of temperature to create a deseasonalized
%anomaly
heat_ltm          = ncread('/Users/smonda15/Ongoing Work/Publication_Compound_Drought_Heatwave_Over_Phoneix/Datasets/CPC Global Long Term Mean/tmax.day.ltm.1991-2020.nc','tmax');

heat_ltm          = reshape(heat_ltm,[720*360,365]);
SWNAheat_ltm      = heat_ltm(SWNA_index,91:273);

SWNAheat_anm      = zeros(2256,183,45);

for i=1:45
    SWNAheat_anm(:,:,i)=SWNAheat_r(:,:,i)-SWNAheat_ltm;
end

SWNAheat_anm1     = reshape(SWNAheat_anm,[2256,183*45]);

SWNAheat_anm_2020 = SWNAheat_anm(:,:,42);
SWNAheat_anm_2023 = SWNAheat_anm(:,:,45);

SWNAheat_anm_2020(1610,:)=0.5*SWNAheat_anm_2020(1609,:)+0.5*SWNAheat_anm_2020(1611,:);
SWNAheat_anm_2023(1610,:)=0.5*SWNAheat_anm_2023(1609,:)+0.5*SWNAheat_anm_2023(1611,:);

%Getting the temperature anomaly for 1981-1990, 1991-2000, 2001-2010
SWNAheat_anm_1981_1990 = SWNAheat_anm(:,:,3:12);
SWNAheat_anm_1991_2000 = SWNAheat_anm(:,:,13:22);
SWNAheat_anm_2001_2010 = SWNAheat_anm(:,:,23:32);
SWNAheat_anm_2011_2020 = SWNAheat_anm(:,:,33:42);

SWNAheat_anm_1981_1990 = mean(SWNAheat_anm_1981_1990,3);
SWNAheat_anm_1991_2000 = mean(SWNAheat_anm_1991_2000,3);
SWNAheat_anm_2001_2010 = mean(SWNAheat_anm_2001_2010,3);
SWNAheat_anm_2011_2020 = mean(SWNAheat_anm_2011_2020,3);


SWNAheat_anm_2001_2010(1610,:)=0.5*SWNAheat_anm_2001_2010(1609,:)+0.5*SWNAheat_anm_2001_2010(1611,:);
SWNAheat_anm_2011_2020(1610,:)=0.5*SWNAheat_anm_2011_2020(1609,:)+0.5*SWNAheat_anm_2011_2020(1611,:);

SWNAheat_anm_1981_1990 = (mean(SWNAheat_anm_1981_1990,1))';
SWNAheat_anm_1991_2000 = (mean(SWNAheat_anm_1991_2000,1))';
SWNAheat_anm_2001_2010 = (mean(SWNAheat_anm_2001_2010,1))';
SWNAheat_anm_2011_2020 = (mean(SWNAheat_anm_2011_2020,1))';

SWNAheat_anm_2023      = (mean(SWNAheat_anm_2023,1))';
SWNAheat_anm_2020      = (mean(SWNAheat_anm_2020,1))';

runningMean1981_1990 = movmean(SWNAheat_anm_1981_1990 , [2 0]); % Includes the current and the previous 2 days
runningMean1991_2000 = movmean(SWNAheat_anm_1991_2000, [2 0]);
runningMean2001_2010 = movmean(SWNAheat_anm_2001_2010,[2 0]);
runningMean2011_2020 = movmean(SWNAheat_anm_2011_2020,[2 0]);

runningMean_2023     = movmean(SWNAheat_anm_2023,[2 0]);

cumMean1981_1990 = cumsum(SWNAheat_anm_1981_1990); % Includes the current and the previous 2 days
cumMean1991_2000 = cumsum(SWNAheat_anm_1991_2000);
cumMean2001_2010 = cumsum(SWNAheat_anm_2001_2010);
cumMean2011_2020 = cumsum(SWNAheat_anm_2011_2020);

cumMean_2023     = cumsum(SWNAheat_anm_2023);
cumMean_2020     = cumsum(SWNAheat_anm_2020);



%Plotting spatially averaged Temperature Anomaly of 2023 with respect to
%other decades
t =1:183;
% Plotting the time series
figure; % Creates a new figure window
plot(t, cumMean_2023 , 'LineWidth', 2); % Plot ts1 with a line width of 2
hold on; % Holds the current plot so that new plots are added to the same figure
%plot(t, cumMean_2020, '--', 'LineWidth', 2); % Plot ts2 with dashed lines
plot(t, cumMean1981_1990, '--', 'LineWidth', 2); % Plot ts2 with dashed lines
plot(t, cumMean1991_2000, ':', 'LineWidth', 2); % Plot ts3 with dotted lines
plot(t, cumMean2001_2010, '-.', 'LineWidth', 2); % Plot ts4 with dash-dot lines
plot(t, cumMean2011_2020, 'LineWidth', 2); % Plot ts5 with a different line style or color

% Enhancing the plot
title('Comparative Time Series Plot'); % Adds a title
xlabel('Time (s)'); % Adds an x-axis label
ylabel('Data Value'); % Adds a y-axis label
legend('2023','1981-1990', '1991-2000', '2001-2010', '2011-2020', 'Location', 'best'); % Adds a legend
grid on; % Adds a grid for better readability

% Optional: Set the axes properties for more control
ax = gca; % Gets the current axes
ax.FontSize = 12; % Sets the font size
ax.LineWidth = 1.5; % Sets the axes line width
ax.Box = 'on'; % Turns the box border on

% Adjust axis limits if necessary
% xlim([min_time max_time]);
% ylim([min_data_value max_data_value]);

hold off; % Releases the hold on the current figure





%Spatially averaged Temperature Anomaly
SWNAheat_anm_SATA   =    mean(SWNAheat_anm,1);
SWNAheat_anm_SATA   =    reshape(SWNAheat_anm_SATA,[1*183,45]);


%Creating Dual Line plot to provide the comparative characteristics of 2020
%and 2023 Temperature anomaly

% Assuming SWNAheat_anm_SATA is your matrix with dimensions 183x45
% And the columns correspond to years 1979-2023 in order

% Extracting data for 2020 and 2023
data2020 = (mean(SWNAheat_anm_2020,1))';
data2023 = (mean(SWNAheat_anm_2023,1))';
dates = 1:183;

% Calculate three-day running mean
runningMean2020 = movmean(data2020, [2 0]); % Includes the current and the previous 2 days
runningMean2023 = movmean(data2023, [2 0]);

figure;
% Plot original data in dashed lines
plot(dates, data2020, 'b--', 'LineWidth', 1); % 2020 data in blue dashed
hold on;
plot(dates, data2023, 'r--', 'LineWidth', 1); % 2023 data in red dashed

% Plot running mean in solid bold lines
plot(dates, runningMean2020, 'b-', 'LineWidth', 2); % 2020 running mean in blue solid
plot(dates, runningMean2023, 'r-', 'LineWidth', 2); % 2023 running mean in red solid
hold off;

% Enhancements for clarity and publication quality
xlabel('Date', 'FontSize', 12);
ylabel('Temperature Anomaly (°C)', 'FontSize', 12);
title('Daily Maximum Temperature Anomalies: 2020 vs 2023', 'FontSize', 14);
legend({'2020', '2023', '2020 3-day Mean', '2023 3-day Mean'}, 'Location', 'northeastoutside');
grid on;
datetick('x', 'mmm', 'keeplimits');

% Adjusting aesthetics
set(gca, 'FontSize', 10);
set(gcf, 'Color', 'w'); % Set background color to white


%Temporally averaged Spatial Anomalies
SWNAheat_anm_2020_APR_JUN_TMIN        = SWNAheat_anm_2020(:,1:91);
SWNAheat_anm_2023_APR_JUN_TMIN        = SWNAheat_anm_2023(:,1:91);
SWNAheat_anm_2020_APR_JUN_TMIN_mean   = mean(SWNAheat_anm_2020_APR_JUN_TMIN,2); 
SWNAheat_anm_2023_APR_JUN_TMIN_mean   = mean(SWNAheat_anm_2023_APR_JUN_TMIN,2);

SWNAheat_anm_2020_JUL_SEP_TMIN        = SWNAheat_anm_2020(:,92:183);
SWNAheat_anm_2023_JUL_SEP_TMIN        = SWNAheat_anm_2023(:,92:183);
SWNAheat_anm_2020_JUL_SEP_TMIN_mean   = mean(SWNAheat_anm_2020_JUL_SEP_TMIN,2);
SWNAheat_anm_2023_JUL_SEP_TMIN_mean   = mean(SWNAheat_anm_2023_JUL_SEP_TMIN,2);







%Heatwave Characteristics calculation start


%Calculating Calender Day 90th and 95th percentile

%Creating sample population for each calender day (5 days each year =
%45*5days=225 days

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

%Calculating the 90th and 95th percentile

SWheat_95 = zeros(2256,183);
SWheat_90 = zeros(2256,183);
SWheat_80 = zeros(2256,183);


for i=1:2256
    for j=1:183
        SWheat_95(i,j)=prctile(SWNAheat_pool(i,j,:),95);
        SWheat_90(i,j)=prctile(SWNAheat_pool(i,j,:),90);
    end
end

%Calculating Spatiotemporal Data with heatwave occurrence as 1 and not
%occurrence as 0

SWheat_r    =  reshape(SWNAheat,[2256,183,45]);
SWheat_b    =  zeros(2256,183,45);
SWheat_d    =  zeros(2256,183,45);

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
        end
    end
end

SWheat_d1 = (reshape(SWheat_d,[2256,183*45]))';
SWheat_b1 = (reshape(SWheat_b,[2256,183*45]))';

%Calculating Occurrence of Heatwaves as Binary Time series        
heatwaveEvents_95 = findHeatwaves(SWheat_b);


heatwaveEvents_90 = findHeatwaves(SWheat_d);


heatwaveEvents_95 = (reshape(heatwaveEvents_95,[2256,183*45]))';
heatwaveEvents_90 = (reshape(heatwaveEvents_90,[2256,183*45]))';


hw_temp_95        = (heatwaveEvents_95.*SWNAheat'); 
hw_temp_90        = (heatwaveEvents_90.*SWNAheat');


%Calculating Duration of each Heatwaves

heatwaveDuration_95       = DurationHeatwaves(SWheat_b);

heatwaveDuration_90       = DurationHeatwaves(SWheat_d);

heatwaveDuration_95       = (reshape(heatwaveDuration_95,[2256,183*45]))';
heatwaveDuration_90       = (reshape(heatwaveDuration_90,[2256,183*45]))';

%Calculating Interarrival time of heatwaves

HW_duration_95            = reshape(heatwaveDuration_95',[2256,183,45]);

HW_duration_90            = reshape(heatwaveDuration_90',[2256,183,45]);

HW_intrarrival_time_90    = zeros(2256,45);
HW_intrarrival_time_95    = zeros(2256,45);

for i=1:2256
    for j=1:45
        a = find(HW_duration_90(i,:,j)>0);
        if size(a,1)==0
            HW_intrarrival_time_90(i,j)=0;
        else
        b = HW_duration_90(i,a,j);
        sum=0;
        for k=2:size(a,2)
            sum=sum+a(1,k)-a(1,k-1)-b(1,k-1);
        end
        sum = sum/size(a,1);
        HW_intrarrival_time_90(i,j)=sum;
        end
    end
end









%Calculating Frequency of Heatwaves
heatwave_Frequency_95 = zeros(2256,1);
heatwave_Frequency_90 = zeros(2256,1);

for i=1:2256
    a=find(heatwaveDuration_90(:,i)>0);
    a=size(a);
    b=find(heatwaveDuration_95(:,i)>0);
    b=size(b);
    heatwave_Frequency_90(i,1)=a(1,1);
    heatwave_Frequency_95(i,1)=b(1,1);
end

isnan(heatwave_Frequency_90);


%Calculating the inter-arrival time of the heatwaves

%heatwaveDuration95 can be used to calculate the interarrival time for each
%of the heatwaves





%Calculating Frequency of Heatwaves for 2020 and 2023
heatwave_Frequency_95 = zeros(2256,1);
heatwave_Frequency_90 = zeros(2256,1);

for i=1:2256
    a=find(heatwaveDuration_90(7504:7686,i)>0);
    a=size(a);
    b=find(heatwaveDuration_95(7504:7686,i)>0);
    b=size(b);
    heatwave_Frequency_90(i,1)=a(1,1);
    heatwave_Frequency_95(i,1)=b(1,1);
end
%Calculating Average Duration of Heatwaves

HW_avg_duration_90 =zeros(2256,1);
HW_avg_duration_95 =zeros(2256,1);

for i=1:2256
    a=find(heatwaveDuration_90(:,i)>0);
    a=heatwaveDuration_90(a,i);
    a=mean(a,1);
    b=find(heatwaveDuration_95(:,i)>0);
    b=heatwaveDuration_95(b,i);
    b=mean(b,1);
    HW_avg_duration_90(i,1)=a;
    HW_avg_duration_95(i,1)=b;
end

%Calculating Duration of Heatwaves during 2020 and 2023
heatwaveDuration_2023_90=heatwaveDuration_90(7504:7686,:);
heatwaveDuration_2023_95=heatwaveDuration_95(7504:7686,:);
HW_avg_duration_90 =zeros(2256,1);
HW_avg_duration_95 =zeros(2256,1);

for i=1:2256
    a=find(heatwaveDuration_2023_90(:,i)>0);
    a=heatwaveDuration_2023_90(a,i);
    a=mean(a,1);
    b=find(heatwaveDuration_2023_95(:,i)>0);
    b=heatwaveDuration_2023_95(b,i);
    b=mean(b,1);
    HW_avg_duration_90(i,1)=a;
    HW_avg_duration_95(i,1)=b;
end

a=find(isnan(HW_avg_duration_95));
HW_avg_duration_95(a,:)=0;

a=find(isnan(HW_avg_duration_90));
HW_avg_duration_90(a,:)=0;

%Calculating Severity and Intensity of Heatwaves
HW_severity_90 = zeros(8235,2256);
HW_severity_95 = zeros(8235,2256);

HW_intensity_90 = zeros(8235,2256);
HW_intensity_95 = zeros(8235,2256);

for i=1:2256
    for j=1:8235
        if heatwaveDuration_90(j,i)>0
            HW_severity_90(j,i)   = max(hw_temp_90(j:(j+heatwaveDuration_90(j,i)-1),i));
            HW_intensity_90(j,i)  = HW_severity_90(j,i);
        end
        if heatwaveDuration_95(j,i)>0
            HW_severity_95(j,i)   = max(hw_temp_95(j:(j+heatwaveDuration_95(j,i)-1),i));
            HW_intensity_95(j,i)  = HW_severity_95(j,i);

          
        end
    end
end


%Calculating Average Severity and Intensity of heatwaves over Southwest USA
HW_avg_amp_95 = zeros(2256,1);
HW_avg_amp_90 = zeros(2256,1);

HW_max_amp_95 = zeros(2256,1);
HW_max_amp_90 = zeros(2256,1);



for i=1:2256
   a=find(HW_intensity_90(:,i)>0);
   a=HW_intensity_90(a,i);
    a=max(a);
    HW_max_amp_90(i,1)=a;

    c=find(HW_intensity_90(:,i)>0);
    c=HW_intensity_90(c,i);
    c=mean(c,1);
    HW_avg_amp_90(i,1)=c;

    b=find(HW_intensity_95(:,i)>0);
    b=HW_intensity_95(b,i);
    b=max(b);
    HW_max_amp_95(i,1)=b;


   

    d=find(HW_intensity_95(:,i)>0);
    d=HW_intensity_95(d,i);
    d=mean(d,1);
    HW_avg_amp_95(i,1)=d;
    
end



%%Calculating Average Average and Maximum Amplitude of heatwaves over
%%Southwest USA for 2020 and 2023
SWNAheat = SWNAheat';

HW_intensity_2020_95  = HW_intensity_95(7504:7686,:);
HW_intensity_2020_90  = HW_intensity_90(7504:7686,:);

HW_avg_amp_95 = zeros(2256,1);
HW_avg_amp_90 = zeros(2256,1);

HW_max_amp_95 = zeros(2256,1);
HW_max_amp_90 = zeros(2256,1);



for i=1:2256
   a=find(HW_intensity_2020_90(:,i)>0);
   if size(a,1)==0
       HW_max_amp_90(i,1)=max(SWNAheat(7504:7686,i));
   else

   a=HW_intensity_2020_90(a,i);
    a=max(a);
    HW_max_amp_90(i,1)=a;
   end

    c=find(HW_intensity_2020_90(:,i)>0);
    if size(c,1)==0
        HW_avg_amp_90(i,1)=mean(SWNAheat(7504:7686,i));
    else
    c=HW_intensity_2020_90(c,i);
    c=mean(c,1);
    HW_avg_amp_90(i,1)=c;
    end

    b=find(HW_intensity_2020_95(:,i)>0);
    if size(b,1)==0
        HW_max_amp_95(i,1)=max(SWNAheat(7504:7686,i));
    else
    b=HW_intensity_2020_95(b,i);
    b=max(b);
    HW_max_amp_95(i,1)=b;
    end


   

    d=find(HW_intensity_2020_95(:,i)>0);
    if size(b,1)==0
        HW_avg_amp_95(i,1)=mean(SWNAheat(7504:7686,i));
    else
    d=HW_intensity_2020_95(d,i);
    d=mean(d,1);
    HW_avg_amp_95(i,1)=d;
    end
    
end











%%Measuring Decadal Change in Heatwave Characteristics

lon = ncread("/Users/smonda15/Ongoing Work/Publication_Compound_Drought_Heatwave_Over_Phoneix/Datasets/CPC Global Temperature/tmax.1979.nc",'lon');
lat = ncread("/Users/smonda15/Ongoing Work/Publication_Compound_Drought_Heatwave_Over_Phoneix/Datasets/CPC Global Temperature/tmax.1979.nc",'lat');

folderPath = '/Users/smonda15/Ongoing Work/Publication_Compound_Drought_Heatwave_Over_Phoneix/Datasets/CPC Global Temperature';
ncFiles = dir(fullfile(folderPath, '*.nc'));


b=zeros(259200,1);
d=zeros(259200,1);
c=1979;

for i = 1:length(ncFiles)
    
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

%Calculating the 90th and 95th percentile

SWheat_95 = zeros(2256,183);
SWheat_90 = zeros(2256,183);


for i=1:2256
    for j=1:183
        SWheat_95(i,j)=prctile(SWNAheat_pool(i,j,:),95);
        SWheat_90(i,j)=prctile(SWNAheat_pool(i,j,:),90);
    end
end

%Calculating Spatiotemporal Data with heatwave occurrence as 1 and not
%occurrence as 0

SWheat_r    =  reshape(SWNAheat,[2256,183,45]);
SWheat_b    =  zeros(2256,183,45);
SWheat_d    =  zeros(2256,183,45);

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
        end
    end
end

SWheat_d1 = (reshape(SWheat_d,[2256,183*45]))';
SWheat_b1 = (reshape(SWheat_b,[2256,183*45]))';

%Calculating Occurrence of Heatwaves as Binary Time series        
heatwaveEvents_95 = findHeatwaves(SWheat_b);


heatwaveEvents_90 = findHeatwaves(SWheat_d);


heatwaveEvents_95 = (reshape(heatwaveEvents_95,[2256,183*45]))';
heatwaveEvents_90 = (reshape(heatwaveEvents_90,[2256,183*45]))';


hw_temp_95        = (heatwaveEvents_95.*SWNAheat'); 
hw_temp_90        = (heatwaveEvents_90.*SWNAheat');


%Calculating Duration of each Heatwaves

heatwaveDuration_95       = DurationHeatwaves(SWheat_b);

heatwaveDuration_90       = DurationHeatwaves(SWheat_d);

heatwaveDuration_95       = (reshape(heatwaveDuration_95,[2256,183*45]))';
heatwaveDuration_90       = (reshape(heatwaveDuration_90,[2256,183*45]))';






%Dividing Decades in 11 (1:2013),11 (2014:4026),11 (4027:6039), 12 (6040:8235) Years

%First Decade
%Calculating Frequency of Heatwaves
heatwave_Frequency_951 = zeros(2256,1);
heatwave_Frequency_901 = zeros(2256,1);

for i=1:2256
    a=find(heatwaveDuration_90(1:2013,i)>0);
    a=size(a);
    b=find(heatwaveDuration_95(1:2013,i)>0);
    b=size(b);
    heatwave_Frequency_901(i,1)=a(1,1);
    heatwave_Frequency_951(i,1)=b(1,1);
end

heatwave_Frequency_952 = zeros(2256,1);
heatwave_Frequency_902 = zeros(2256,1);

for i=1:2256
    a=find(heatwaveDuration_90(2014:4026,i)>0);
    a=size(a);
    b=find(heatwaveDuration_95(2014:4026,i)>0);
    b=size(b);
    heatwave_Frequency_902(i,1)=a(1,1);
    heatwave_Frequency_952(i,1)=b(1,1);
end

heatwave_Frequency_953 = zeros(2256,1);
heatwave_Frequency_903 = zeros(2256,1);

for i=1:2256
    a=find(heatwaveDuration_90(4027:6039,i)>0);
    a=size(a);
    b=find(heatwaveDuration_95(4027:6039,i)>0);
    b=size(b);
    heatwave_Frequency_903(i,1)=a(1,1);
    heatwave_Frequency_953(i,1)=b(1,1);
end

heatwave_Frequency_954 = zeros(2256,1);
heatwave_Frequency_904 = zeros(2256,1);

for i=1:2256
    a=find(heatwaveDuration_90(6040:8235,i)>0);
    a=size(a);
    b=find(heatwaveDuration_95(6040:8235,i)>0);
    b=size(b);
    heatwave_Frequency_904(i,1)=a(1,1);
    heatwave_Frequency_954(i,1)=b(1,1);
end

%isnan(heatwave_Frequency_90);



%Calculating Decadal Change in Duration of Heatwaves

HW_avg_duration_901 =zeros(2256,1);
HW_avg_duration_951 =zeros(2256,1);

for i=1:2256
    a=find(heatwaveDuration_90(1:2013,i)>0);
    a=heatwaveDuration_90(a,i);
    a=mean(a,1);
    b=find(heatwaveDuration_95(1:2013,i)>0);
    b=heatwaveDuration_95(b,i);
    b=mean(b,1);
    HW_avg_duration_901(i,1)=a;
    HW_avg_duration_951(i,1)=b;
end

HW_avg_duration_902 =zeros(2256,1);
HW_avg_duration_952 =zeros(2256,1);

for i=1:2256
    a=find(heatwaveDuration_90(2014:4026,i)>0);
    a=heatwaveDuration_90(a+2013,i);
    a=mean(a,1);
    b=find(heatwaveDuration_95(2014:4026,i)>0);
    b=heatwaveDuration_95(b+2013,i);
    b=mean(b,1);
    HW_avg_duration_902(i,1)=a;
    HW_avg_duration_952(i,1)=b;
end

HW_avg_duration_903 =zeros(2256,1);
HW_avg_duration_953 =zeros(2256,1);

for i=1:2256
    a=find(heatwaveDuration_90(4027:6039,i)>0);
    a=heatwaveDuration_90(a+4026,i);
    a=mean(a,1);
    b=find(heatwaveDuration_95(4027:6039,i)>0);
    b=heatwaveDuration_95(b+4026,i);
    b=mean(b,1);
    HW_avg_duration_903(i,1)=a;
    HW_avg_duration_953(i,1)=b;
end


HW_avg_duration_904 =zeros(2256,1);
HW_avg_duration_954 =zeros(2256,1);

for i=1:2256
    a=find(heatwaveDuration_90(6040:8235,i)>0);
    a=heatwaveDuration_90(a+6039,i);
    a=mean(a,1);
    b=find(heatwaveDuration_95(6040:8235,i)>0);
    b=heatwaveDuration_95(b+6039,i);
    b=mean(b,1);
    HW_avg_duration_904(i,1)=a;
    HW_avg_duration_954(i,1)=b;
end








%Calculating Severity and Intensity of Heatwaves
HW_severity_90 = zeros(8235,2256);
HW_severity_95 = zeros(8235,2256);

HW_intensity_90 = zeros(8235,2256);
HW_intensity_95 = zeros(8235,2256);

for i=1:2256
    for j=1:8235
        if heatwaveDuration_90(j,i)>0
            HW_severity_90(j,i)   = max(hw_temp_90(j:(j+heatwaveDuration_90(j,i)-1),i));
            HW_intensity_90(j,i)  = HW_severity_90(j,i);
        end
        if heatwaveDuration_95(j,i)>0
            HW_severity_95(j,i)   = max(hw_temp_95(j:(j+heatwaveDuration_95(j,i)-1),i));
            HW_intensity_95(j,i)  = HW_severity_95(j,i);

          
        end
    end
end





%Analysis of Drought over Southwestern North America

%%Extracting CPC Precipitation data
lon = ncread("/Users/smonda15/Ongoing Work/Publication_Compound_Drought_Heatwave_Over_Phoneix/Datasets/CPC Global Precipitation/precip.1979.nc",'lon');
lat = ncread("/Users/smonda15/Ongoing Work/Publication_Compound_Drought_Heatwave_Over_Phoneix/Datasets/CPC Global Precipitation/precip.1979.nc",'lat');

folderPath = '/Users/smonda15/Ongoing Work/Publication_Compound_Drought_Heatwave_Over_Phoneix/Datasets/CPC Global Precipitation';
ncFiles = dir(fullfile(folderPath, '*.nc'));


b=zeros(259200,1);
d=zeros(259200,1);

c=1979;

for i = 1:45
    
    filePath = fullfile(folderPath, ncFiles(i).name);
    
    
    precipData = ncread(filePath, 'precip');
    a       = size(precipData);
    
  
    precipData = reshape(precipData,[a(1)*a(2),a(3)]);
    isLeapYear = (rem(c, 4) == 0 && rem(c, 100) ~= 0) || (rem(c, 400) == 0);
    if isLeapYear
        precipData1 = precipData(:,92:274);
        precipData2 = precipData;
    else 
        precipData1 = precipData(:,91:273);
        precipData2 = precipData;
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



SW_SPI12     =   reshape(SW_SPI12,[52,44,2256]);
SW_SPI3      =   reshape(SW_SPI3,[52,44,2256]);

%Extracting SPI3 for summer 

%SW_SPI3     =   SW_SPI3(14:39,:,:);
%SW_SPI3     =   reshape(SW_SPI3,[26*44,2256]);

%Using SPI=-1 as Drought Threshold moderate drought and -1.5 for extreme
%drought

SW_drought3 = zeros(2288,2256);

for i=1:2288
    for j=1:2256
        if SW_SPI3(i,j)<-1.0
            SW_drought3(i,j)=1;
        end
    end
end

SW_drought12 = zeros(2288,2256);

for i=1:2288
    for j=1:2256
        if SW_SPI12(i,j)<-1.0
            SW_drought12(i,j)=1;
        end
    end
end


%SW_drought = SW_drought(531:543,:);

% Assuming your matrix is called 'droughtData' and is of size 585x439
numWeeks = 2288;
numLocations = 2256;

% Initialize the matrices to store results
numDroughtEvents = zeros(1, numLocations);
DroughtSeverity = zeros(numWeeks, numLocations);
DroughtDuration = zeros(numWeeks, numLocations);
DroughtIntensity = zeros(numWeeks, numLocations);




for loc = 1:numLocations
    droughtPeriods = SW_drought12(:, loc);      %for 2020, 2081:2132 & %2237:2288
    droughtspi     = SW_SPI12(:,loc);
    % ... [any additional code for specific periods] ...

    % Find the start and end indices of each drought event
    startIndices = find(diff([0; droughtPeriods]) == 1);
    endIndices = find(diff([droughtPeriods; 0]) == -1);

    % Initialize the number of events for this location
    numEvents = 0;
    
    for i = 1:length(startIndices)
        eventDuration = endIndices(i) - startIndices(i) + 1;
        
        % Check if the event duration is at least two weeks
        if eventDuration >= 2
            numEvents = numEvents + 1; % Counting the valid event
            DroughtDuration(startIndices(i), loc) = eventDuration;
            
            eventSeverity = sum(droughtspi(startIndices(i):endIndices(i)));
            eventIntensity = eventSeverity / eventDuration;

            DroughtSeverity(startIndices(i), loc) = eventSeverity;
            DroughtIntensity(startIndices(i), loc) = eventIntensity;
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
    a  =  find(DroughtDuration(:,i)>0);
    b  =  DroughtDuration(a,i);
    b  =  mean(b,1);
    SW_dr_avg_duration(i,1)=b;
    c  =  DroughtSeverity(a,i);
    c  =  mean(c,1);
    SW_dr_avg_severity(i,1)=c;
    d  =  DroughtIntensity(a,i);
    d  =  mean(d,1);
    SW_dr_avg_intensity(i,1)=d;
end

    
a =find(isnan(SW_dr_avg_duration));

SW_dr_avg_duration(a,:) = 0;

SW_dr_avg_intensity(a,:) = 0;

SW_dr_avg_severity(a,:) = 0;


%Calculating Drought Characteristics for Each Year (For each year, 
% 
% we can consider two type of droughts: one which has started in that year or which
%has ended that year. Assuming we are focusing on droughts which started in
%each year or summer

SW_SPI3     =   reshape(SW_SPI3,[52,44,2256]);
SW_drought  =   reshape(SW_drought,[52,44,2256]);


% Assuming your matrix is called 'droughtData' and is of size 585x439
numWeeks = 2288;
numLocations = 2256;

% Initialize the matrices to store results
numDroughtEvents = zeros(1, numLocations);
DroughtSeverity = zeros(numWeeks, numLocations);
DroughtDuration = zeros(numWeeks, numLocations);
DroughtIntensity = zeros(numWeeks, numLocations);


for loc = 1:numLocations
    droughtPeriods = SW_drought(2237:2288, loc);      %for 2020, 2081:2132 & %2237:2288
    droughtspi     = SW_SPI3(2237:2288,loc);
    % ... [any additional code for specific periods] ...

    % Find the start and end indices of each drought event
    startIndices = find(diff([0; droughtPeriods]) == 1);
    endIndices = find(diff([droughtPeriods; 0]) == -1);

    % Initialize the number of events for this location
    numEvents = 0;
    
    for i = 1:length(startIndices)
        eventDuration = endIndices(i) - startIndices(i) + 1;
        
        % Check if the event duration is at least two weeks
        if eventDuration >= 2
            numEvents = numEvents + 1; % Counting the valid event
            DroughtDuration(startIndices(i), loc) = eventDuration;
            
            eventSeverity = sum(droughtspi(startIndices(i):endIndices(i)));
            eventIntensity = eventSeverity / eventDuration;

            DroughtSeverity(startIndices(i), loc) = eventSeverity;
            DroughtIntensity(startIndices(i), loc) = eventIntensity;
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
    a  =  find(DroughtDuration(:,i)>0);
    b  =  DroughtDuration(a,i);
    b  =  mean(b,1);
    SW_dr_avg_duration(i,1)=b;
    c  =  DroughtSeverity(a,i);
    c  =  mean(c,1);
    SW_dr_avg_severity(i,1)=c;
    d  =  DroughtIntensity(a,i);
    d  =  mean(d,1);
    SW_dr_avg_intensity(i,1)=d;
end

    
a =find(isnan(SW_dr_avg_duration));

SW_dr_avg_duration(a,:) = 0;

SW_dr_avg_intensity(a,:) = 0;

SW_dr_avg_severity(a,:) = 0;








%Checking Drought Characteristics over Different Decades

%SW_drought = SW_drought(531:543,:);

% Assuming your matrix is called 'droughtData' and is of size 585x439
numWeeks = 572;
numLocations = 2256;

% First 11 Years (1979-1990)
numDroughtEvents1 = zeros(1, numLocations);
DroughtSeverity = zeros(numWeeks, numLocations);
DroughtDuration = zeros(numWeeks, numLocations);
DroughtIntensity = zeros(numWeeks, numLocations);




for loc = 1:numLocations
    droughtPeriods = SW_drought3(1:572, loc);      %for 2020, 2081:2132 & %2237:2288
    droughtspi     = SW_SPI3(1:572,loc);
    % ... [any additional code for specific periods] ...

    % Find the start and end indices of each drought event
    startIndices = find(diff([0; droughtPeriods]) == 1);
    endIndices = find(diff([droughtPeriods; 0]) == -1);

    % Initialize the number of events for this location
    numEvents = 0;
    
    for i = 1:length(startIndices)
        eventDuration = endIndices(i) - startIndices(i) + 1;
        
        % Check if the event duration is at least two weeks
        if eventDuration >= 2
            numEvents = numEvents + 1; % Counting the valid event
            DroughtDuration(startIndices(i), loc) = eventDuration;
            
            eventSeverity = sum(droughtspi(startIndices(i):endIndices(i)));
            eventIntensity = eventSeverity / eventDuration;

            DroughtSeverity(startIndices(i), loc) = eventSeverity;
            DroughtIntensity(startIndices(i), loc) = eventIntensity;
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
    a  =  find(DroughtDuration(:,i)>0);
    b  =  DroughtDuration(a,i);
    b  =  mean(b,1);
    SW_dr_avg_duration1(i,1)=b;
    c  =  DroughtSeverity(a,i);
    c  =  mean(c,1);
    SW_dr_avg_severity1(i,1)=c;
    d  =  DroughtIntensity(a,i);
    d  =  mean(d,1);
    SW_dr_avg_intensity1(i,1)=d;
end

    
a =find(isnan(SW_dr_avg_duration1));

SW_dr_avg_duration1(a,:) = 0;

SW_dr_avg_intensity1(a,:) = 0;

SW_dr_avg_severity1(a,:) = 0;



% Second 11 Years (1991-2001)
numDroughtEvents2 = zeros(1, numLocations);
DroughtSeverity = zeros(numWeeks, numLocations);
DroughtDuration = zeros(numWeeks, numLocations);
DroughtIntensity = zeros(numWeeks, numLocations);




for loc = 1:numLocations
    droughtPeriods = SW_drought3(573:1144, loc);      %for 2020, 2081:2132 & %2237:2288
    droughtspi     = SW_SPI3(573:1144,loc);
    % ... [any additional code for specific periods] ...

    % Find the start and end indices of each drought event
    startIndices = find(diff([0; droughtPeriods]) == 1);
    endIndices = find(diff([droughtPeriods; 0]) == -1);

    % Initialize the number of events for this location
    numEvents = 0;
    
    for i = 1:length(startIndices)
        eventDuration = endIndices(i) - startIndices(i) + 1;
        
        % Check if the event duration is at least two weeks
        if eventDuration >= 2
            numEvents = numEvents + 1; % Counting the valid event
            DroughtDuration(startIndices(i), loc) = eventDuration;
            
            eventSeverity = sum(droughtspi(startIndices(i):endIndices(i)));
            eventIntensity = eventSeverity / eventDuration;

            DroughtSeverity(startIndices(i), loc) = eventSeverity;
            DroughtIntensity(startIndices(i), loc) = eventIntensity;
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
    a  =  find(DroughtDuration(:,i)>0);
    b  =  DroughtDuration(a,i);
    b  =  mean(b,1);
    SW_dr_avg_duration2(i,1)=b;
    c  =  DroughtSeverity(a,i);
    c  =  mean(c,1);
    SW_dr_avg_severity2(i,1)=c;
    d  =  DroughtIntensity(a,i);
    d  =  mean(d,1);
    SW_dr_avg_intensity2(i,1)=d;
end

    
a =find(isnan(SW_dr_avg_duration2));

SW_dr_avg_duration2(a,:) = 0;

SW_dr_avg_intensity2(a,:) = 0;

SW_dr_avg_severity2(a,:) = 0;


% Third 11 Years (2002-2012)
numDroughtEvents3 = zeros(1, numLocations);
DroughtSeverity = zeros(numWeeks, numLocations);
DroughtDuration = zeros(numWeeks, numLocations);
DroughtIntensity = zeros(numWeeks, numLocations);




for loc = 1:numLocations
    droughtPeriods = SW_drought3(1145:1716, loc);      %for 2020, 2081:2132 & %2237:2288
    droughtspi     = SW_SPI3(1145:1716,loc);
    % ... [any additional code for specific periods] ...

    % Find the start and end indices of each drought event
    startIndices = find(diff([0; droughtPeriods]) == 1);
    endIndices = find(diff([droughtPeriods; 0]) == -1);

    % Initialize the number of events for this location
    numEvents = 0;
    
    for i = 1:length(startIndices)
        eventDuration = endIndices(i) - startIndices(i) + 1;
        
        % Check if the event duration is at least two weeks
        if eventDuration >= 2
            numEvents = numEvents + 1; % Counting the valid event
            DroughtDuration(startIndices(i), loc) = eventDuration;
            
            eventSeverity = sum(droughtspi(startIndices(i):endIndices(i)));
            eventIntensity = eventSeverity / eventDuration;

            DroughtSeverity(startIndices(i), loc) = eventSeverity;
            DroughtIntensity(startIndices(i), loc) = eventIntensity;
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
    a  =  find(DroughtDuration(:,i)>0);
    b  =  DroughtDuration(a,i);
    b  =  mean(b,1);
    SW_dr_avg_duration3(i,1)=b;
    c  =  DroughtSeverity(a,i);
    c  =  mean(c,1);
    SW_dr_avg_severity3(i,1)=c;
    d  =  DroughtIntensity(a,i);
    d  =  mean(d,1);
    SW_dr_avg_intensity3(i,1)=d;
end

    
a =find(isnan(SW_dr_avg_duration3));

SW_dr_avg_duration3(a,:) = 0;

SW_dr_avg_intensity3(a,:) = 0;

SW_dr_avg_severity3(a,:) = 0;


% Last 11 Years (2013-2024)
numDroughtEvents4 = zeros(1, numLocations);
DroughtSeverity = zeros(numWeeks, numLocations);
DroughtDuration = zeros(numWeeks, numLocations);
DroughtIntensity = zeros(numWeeks, numLocations);




for loc = 1:numLocations
    droughtPeriods = SW_drought3(1717:2288, loc);      %for 2020, 2081:2132 & %2237:2288
    droughtspi     = SW_SPI3(1717:2288,loc);
    % ... [any additional code for specific periods] ...

    % Find the start and end indices of each drought event
    startIndices = find(diff([0; droughtPeriods]) == 1);
    endIndices = find(diff([droughtPeriods; 0]) == -1);

    % Initialize the number of events for this location
    numEvents = 0;
    
    for i = 1:length(startIndices)
        eventDuration = endIndices(i) - startIndices(i) + 1;
        
        % Check if the event duration is at least two weeks
        if eventDuration >= 2
            numEvents = numEvents + 1; % Counting the valid event
            DroughtDuration(startIndices(i), loc) = eventDuration;
            
            eventSeverity = sum(droughtspi(startIndices(i):endIndices(i)));
            eventIntensity = eventSeverity / eventDuration;

            DroughtSeverity(startIndices(i), loc) = eventSeverity;
            DroughtIntensity(startIndices(i), loc) = eventIntensity;
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
    a  =  find(DroughtDuration(:,i)>0);
    b  =  DroughtDuration(a,i);
    b  =  mean(b,1);
    SW_dr_avg_duration4(i,1)=b;
    c  =  DroughtSeverity(a,i);
    c  =  mean(c,1);
    SW_dr_avg_severity4(i,1)=c;
    d  =  DroughtIntensity(a,i);
    d  =  mean(d,1);
    SW_dr_avg_intensity4(i,1)=d;
end

    
a =find(isnan(SW_dr_avg_duration4));

SW_dr_avg_duration4(a,:) = 0;

SW_dr_avg_intensity4(a,:) = 0;

SW_dr_avg_severity4(a,:) = 0;



%%Analyzing Compound Drought and Heatwaves


%%Extracting CPC Precipitation data
lon = ncread("/Users/smonda15/Ongoing Work/Publication_Compound_Drought_Heatwave_Over_Phoneix/Datasets/CPC Global Precipitation/precip.1979.nc",'lon');
lat = ncread("/Users/smonda15/Ongoing Work/Publication_Compound_Drought_Heatwave_Over_Phoneix/Datasets/CPC Global Precipitation/precip.1979.nc",'lat');

folderPath = '/Users/smonda15/Ongoing Work/Publication_Compound_Drought_Heatwave_Over_Phoneix/Datasets/CPC Global Precipitation';
ncFiles = dir(fullfile(folderPath, '*.nc'));


b=zeros(259200,1);
d=zeros(259200,1);

c=1979;

for i = 1:length(ncFiles)
    
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



SW_SPI12     =   reshape(SW_SPI12,[52,44,2256]);
SW_SPI3      =   reshape(SW_SPI3,[52,44,2256]);

%Extracting SPI3 for summer 

%SW_SPI3     =   SW_SPI3(14:39,:,:);
%SW_SPI3     =   reshape(SW_SPI3,[26*44,2256]);

%Using SPI=-1 as Drought Threshold moderate drought and -1.5 for extreme
%drought

SW_drought3 = zeros(2288,2256);

for i=1:2288
    for j=1:2256
        if SW_SPI3(i,j)<-1.0
            SW_drought3(i,j)=1;
        end
    end
end

SW_drought12 = zeros(2288,2256);

for i=1:2288
    for j=1:2256
        if SW_SPI12(i,j)<-1.0
            SW_drought12(i,j)=1;
        end
    end
end

%converting heatwave event daily data to weekly data

heatwaveEvents_90_r = zeros(2256,183,45);

for i=1:45
    m = 183*(i-1)+1;
    n = 183*i;

    heatwaveEvents_90_r(:,:,i)= (heatwaveEvents_90(m:n,:))';
end

heatwaveEvents_90_r = heatwaveEvents_90_r(:,1:182,:);

heatwaveEvents_90_r = reshape(heatwaveEvents_90_r,[2256,182*45]);

heatwaveEvents_90_rw = zeros(2256,1170,7);

for i=1:1170
    k = 7*(i-1)+1;
    l = 7*i;
    heatwaveEvents_90_rw(:,i,:) = (heatwaveEvents_90_r(:,k:l));
end

heatwaveEvents_90_rwsum  =  sum(heatwaveEvents_90_rw,3);


%Extracting data from 1980 to match with drought

heatwaveEvents_90_rwsum      = heatwaveEvents_90_rwsum(:,27:1170);


%Changing Drought data to match the same dimension as the heatwaves
%SW_drought3_r       = SW_drought3';
%SW_drought3_r       = reshape(SW_drought3_r,[2256,52,44]);
%SW_drought3_r       = SW_drought3_r(:,14:39,:);
%SW_drought3_r       = reshape(SW_drought3_r,[2256,26*44]);


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

heatwaveEvents90_w = zeros(2256,1144);

for i=1:2256
    for j=1:1144
        if heatwaveEvents_90_rwsum(i,j)>=2
            heatwaveEvents90_w(i,j)=1;
        end
    end
end

heatwaveEvents90_wsum = (sum(heatwaveEvents90_w,1))';

SW_drought3_r_sum     = (sum(SW_drought3_r,1))';

Compound_dr_hw  = zeros(2256,1144);

for i=1:2256
    for j=1:1144
        if heatwaveEvents90_w(i,j)==1 && SW_drought3_r(i,j)==1
            Compound_dr_hw(i,j)=1;
        end
    end
end

Compound_dr_hw_sum    = (sum(Compound_dr_hw,1))';

heatwaveEvents90_wsum = reshape(heatwaveEvents90_wsum,[26,44]);
heatwaveEvents90_wsum = (mean(heatwaveEvents90_wsum,1))';

SW_drought3_r_sum     = reshape(SW_drought3_r_sum,[26,44]);
SW_drought3_r_sum     = (mean(SW_drought3_r_sum,1))';

Compound_dr_hw_sum    = reshape(Compound_dr_hw_sum,[26,44]);
Compound_dr_hw_sum    = (mean(Compound_dr_hw_sum,1))';


hw_drought_dependency = [heatwaveEvents90_wsum,SW_drought3_r_sum,Compound_dr_hw_sum];

% Assume matrix A contains your data, where A = [Vector1, Vector2, Vector3]
A = hw_drought_dependency; % Example data. Replace with your actual matrix.

% Creating the figure
figure;

% Setting up the plot for the first two vectors on the primary y-axis
yyaxis left;
plot(A(:,1), 'LineWidth', 2, 'Color', 'blue');
hold on;
plot(A(:,2), 'LineWidth', 2, 'Color', 'red');
ylabel('Primary Y-axis Label', 'FontSize', 12, 'FontWeight', 'bold');

% Setting up the plot for the third vector on the secondary y-axis
yyaxis right;
plot(A(:,3), 'LineWidth', 2, 'LineStyle', '--', 'Color', 'black');
ylabel('Secondary Y-axis Label', 'FontSize', 12, 'FontWeight', 'bold');

% Additional plot customizations
xlabel('X-axis Label', 'FontSize', 12, 'FontWeight', 'bold');
title('Title of Your Graph', 'FontSize', 14, 'FontWeight', 'bold');
legend({'Vector 1', 'Vector 2', 'Vector 3'}, 'Location', 'bestoutside', 'FontSize', 10);
grid on;
box on;
set(gca, 'FontSize', 10, 'LineWidth', 1.5);
set(gcf, 'Color', 'w'); % Sets background color to white

% Adjust these parameters as necessary to fit your specific data and visualization needs.



%Calculating Trend in frequency of compound drought and heatwave

%When dividing the years to reshape, keep the years 

%CDHW                  = zeros(2256,44);

%Compound_dr_hw1       = reshape(Compound_dr_hw,[2256,26,44]);

%Compound_dr_hwsum1    = sum(Compound_dr_hw1,2);
%Compound_dr_hwsum2    = reshape(Compound_dr_hwsum1,[2256,1*44]);

%Calculating the average spatiotemporal characteristics of compound drought
%and heatwave

%Duration & Frequency

% Assuming your data matrix is named 'compoundEvents'
numGrids = 2256; % Number of grid locations
numWeeks = 1144; % Number of weeks

% Initialize output matrices
averageFrequency = zeros(numGrids, 1); % Average frequency for each grid location
averageDuration = zeros(numGrids, 1); % Average duration for each grid location


for i = 1:numGrids
    % Extract the time series for the current grid location
    timeSeries = Compound_dr_hw(i, :);

    % Frequency calculation
    totalEvents = sum(timeSeries);
    averageFrequency(i) = totalEvents / 44;

    % Duration calculation
    durations = [];
    currentDuration = 0;
    for t = 1:numWeeks
        if timeSeries(t) == 1
            currentDuration = currentDuration + 1; % Continue counting the length of the current event
        elseif currentDuration > 0
            durations = [durations, currentDuration]; % End of an event, add its duration to the list
            currentDuration = 0; % Reset for the next event
        end
    end
    % Check if the last period is a continuous event
    if currentDuration > 0
        durations = [durations, currentDuration]; % Include the duration of the last event
    end

    % Calculate the average duration if there are any durations recorded
    if ~isempty(durations)
        averageDuration(i) = mean(durations);
    else
        averageDuration(i) = 0; % No events for this grid location
    end
end



%Calculating Precipitation Anomalies for April-May-June and
%July-August-September

%%Extracting CPC Precipitation data
lon = ncread("/Users/smonda15/Ongoing Work/Publication_Compound_Drought_Heatwave_Over_Phoneix/Datasets/CPC Global Precipitation/precip.1979.nc",'lon');
lat = ncread("/Users/smonda15/Ongoing Work/Publication_Compound_Drought_Heatwave_Over_Phoneix/Datasets/CPC Global Precipitation/precip.1979.nc",'lat');

folderPath = '/Users/smonda15/Ongoing Work/Publication_Compound_Drought_Heatwave_Over_Phoneix/Datasets/CPC Global Precipitation';
ncFiles = dir(fullfile(folderPath, '*.nc'));


b=zeros(259200,1);
d=zeros(259200,1);

c=1979;

for i = 1:length(ncFiles)
    
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

%Extracting Long term mean data

pcp_ltm         =  ncread('/Users/smonda15/Ongoing Work/Publication_Compound_Drought_Heatwave_Over_Phoneix/Datasets/CPC Global Long Term Mean/precip.day.1991-2020.ltm.nc','precip');

pcp_ltm         =  reshape(pcp_ltm,[720*360,365]);

pcp_ltm_summer  =  pcp_ltm(SWNA_index,91:273);


SW_pcp_anm  = zeros(8235,2256);

SW_summer_pcp1 = SW_summer_pcp';

for i=1:2256
    for j=1:45
        k =183*(j-1)+1;
        l = 183*j;

        SW_pcp_anm(k:l,i)=SW_summer_pcp1(k:l,i)-(pcp_ltm_summer(i,:))';
    end
end

SW_pcp_anm_2023_Apr_June       = (mean(SW_pcp_anm(8053:8143,:),1))';
SW_pcp_anm_2023_Jul_Aug       = (mean(SW_pcp_anm(8144:8235,:),1))';

SW_pcp_anm_2023_cum   = cumsum(SW_pcp_anm_2023);


%pcp_ltm  =  pcp_ltm(SWNA_index,1:364);

numLocations = 2256;
numDays = 364;
numYears = 1;

% Reshape the data to have one year per layer
numDaysPerYear = 364;
dataMatrixReshaped = reshape(pcp_ltm, numLocations, numDaysPerYear, numYears);

% Determine the number of weeks (assuming 12 weeks per summer, adjust as necessary)
numWeeksPerYear = 52;

% Initialize the weekly data matrix
weeklyDataMatrix_ltm = zeros(numLocations, numWeeksPerYear, numYears);

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
        weeklyDataMatrix_ltm(:, week, year) = sum(dataMatrixReshaped(:, startDay:endDay, year), 2);
    end
end

% Reshape the weekly data matrix back to two dimensions if needed
weeklyDataMatrix_rltm = reshape(weeklyDataMatrix_ltm, numLocations, numWeeksPerYear * numYears);
weeklyDataMatrix_rltm(find(weeklyDataMatrix_rltm<0))=0;

%Calculating weekly precipitation data anomaly

weeklyDataMatrix_anm = zeros(2256,52,45);

for i=1:45
    weeklyDataMatrix_anm(:,:,i)=weeklyDataMatrix(:,:,i)-weeklyDataMatrix_rltm;
end


%Calculating Anomalies for 2020 and 2023

%2020

weeklyDataMatrix_anm_2020_summer = weeklyDataMatrix_anm(:,14:39,42);
weeklyDataMatrix_anm_2023_summer = weeklyDataMatrix_anm(:,14:39,45);

%Getting anomalies corresponding to different decades
weeklyDataMatrix_anm_1981_1990 = weeklyDataMatrix_anm(:,14:39,3:12);
weeklyDataMatrix_anm_1991_2000 = weeklyDataMatrix_anm(:,14:39,13:22);
weeklyDataMatrix_anm_1991_2000(:,:,3) = [];
weeklyDataMatrix_anm_2001_2010 = weeklyDataMatrix_anm(:,14:39,23:32);
weeklyDataMatrix_anm_2011_2020 = weeklyDataMatrix_anm(:,14:39,33:42);

weeklyDataMatrix_anm_1981_1990 = mean(weeklyDataMatrix_anm_1981_1990,3);
weeklyDataMatrix_anm_1991_2000 = mean(weeklyDataMatrix_anm_1991_2000,3);
weeklyDataMatrix_anm_2001_2010 = mean(weeklyDataMatrix_anm_2001_2010,3);
weeklyDataMatrix_anm_2011_2020 = mean(weeklyDataMatrix_anm_2011_2020,3);


%Converting to Temporal Anomaly

weeklyDataMatrix_anm_1981_1990    = cumsum(mean(weeklyDataMatrix_anm_1981_1990,1))';
weeklyDataMatrix_anm_1991_2000    = cumsum(mean(weeklyDataMatrix_anm_1991_2000,1))';
weeklyDataMatrix_anm_2001_2010    = cumsum(mean(weeklyDataMatrix_anm_2001_2010,1))';
weeklyDataMatrix_anm_2011_2020    = cumsum(mean(weeklyDataMatrix_anm_2011_2020,1))';
weeklyDataMatrix_anm_2023         = (mean(weeklyDataMatrix_anm_2023_summer,1))';
weeklyDataMatrix_anm_2023_summer  = cumsum(weeklyDataMatrix_anm_2023(14:26));


t =1:26;
% Plotting the time series
figure; % Creates a new figure window
plot(t, weeklyDataMatrix_anm_2023 , 'LineWidth', 2); % Plot ts1 with a line width of 2
hold on; % Holds the current plot so that new plots are added to the same figure
%plot(t, cumMean_2020, '--', 'LineWidth', 2); % Plot ts2 with dashed lines
plot(t, weeklyDataMatrix_anm_1981_1990, '--', 'LineWidth', 2); % Plot ts2 with dashed lines
plot(t, weeklyDataMatrix_anm_1991_2000, ':', 'LineWidth', 2); % Plot ts3 with dotted lines
plot(t, weeklyDataMatrix_anm_2001_2010, '-.', 'LineWidth', 2); % Plot ts4 with dash-dot lines
plot(t, weeklyDataMatrix_anm_2011_2020, 'LineWidth', 2); % Plot ts5 with a different line style or color

% Enhancing the plot
title('Comparative Time Series Plot'); % Adds a title
xlabel('Time (s)'); % Adds an x-axis label
ylabel('Data Value'); % Adds a y-axis label
legend('2023','1981-1990', '1991-2000', '2001-2010', '2011-2020', 'Location', 'best'); % Adds a legend
grid on; % Adds a grid for better readability

% Optional: Set the axes properties for more control
ax = gca; % Gets the current axes
ax.FontSize = 12; % Sets the font size
ax.LineWidth = 1.5; % Sets the axes line width
ax.Box = 'on'; % Turns the box border on

% Adjust axis limits if necessary
% xlim([min_time max_time]);
% ylim([min_data_value max_data_value]);

hold off; % Releases the hold on the current figure




%Spatially averaged Cumuliative Precipitation Anomaly

weeklyDataMatrix_anm_2020_summer_savg_cum = (mean(weeklyDataMatrix_anm_2020_summer,1))'; 

weeklyDataMatrix_anm_2023_summer_savg_cum = (mean(weeklyDataMatrix_anm_2023_summer,1))';

runningMean2020 = movmean(weeklyDataMatrix_anm_2020_summer_savg_cum, [2 0]); % Includes the current and the previous 2 days
runningMean2023 = movmean(weeklyDataMatrix_anm_2023_summer_savg_cum , [2 0]);

dates =1:26;
figure;
% Plot original data in dashed lines
plot(dates, weeklyDataMatrix_anm_2020_summer_savg_cum, 'b--', 'LineWidth', 1); % 2020 data in blue dashed
hold on;
plot(dates, weeklyDataMatrix_anm_2023_summer_savg_cum, 'r--', 'LineWidth', 1); % 2023 data in red dashed

% Plot running mean in solid bold lines
plot(dates, runningMean2020, 'b-', 'LineWidth', 2); % 2020 running mean in blue solid
plot(dates, runningMean2023, 'r-', 'LineWidth', 2); % 2023 running mean in red solid
hold off;

% Enhancements for clarity and publication quality
xlabel('Date', 'FontSize', 12);
ylabel('Cumuliative Precipitation Anomaly (°C)', 'FontSize', 12);
title('Daily Cumuliative Precipitation Anomalies: 2020 vs 2023', 'FontSize', 14);
legend({'2020', '2023', '2020 3-day Mean', '2023 3-day Mean'}, 'Location', 'northeastoutside');
grid on;
datetick('x', 'mmm', 'keeplimits');

% Adjusting aesthetics
set(gca, 'FontSize', 10);
set(gcf, 'Color', 'w'); % Set background color to white


%Temporally averaged Spatial Anomaly of April-May-June &
%July-August-September

weeklyDataMatrix_anm_2020_summer_Apr_Jun = mean(weeklyDataMatrix_anm_2020_summer(:,1:13),2);
weeklyDataMatrix_anm_2020_summer_Jul_Aug = mean(weeklyDataMatrix_anm_2020_summer(:,14:26),2);


weeklyDataMatrix_anm_2023_summer_Apr_Jun = mean(weeklyDataMatrix_anm_2023_summer(:,1:13),2);
weeklyDataMatrix_anm_2023_summer_Jul_Aug = mean(weeklyDataMatrix_anm_2023_summer(:,14:26),2);

% Assuming x, y, and z are your 2256x1 vectors
x = SWNAheat_anm_2020_JUL_SEP_TMAX_mean;
y = SWNAheat_anm_2020_JUL_SEP_TMIN_mean;
z = weeklyDataMatrix_anm_2023_summer_Jul_Aug;

% Generate a figure
figure;

% Create a 3D scatter plot
scatter3(x, y, z, 36, z, 'filled'); % 36 is the size of the markers
hold on; % Keep the plot for further modifications

% Customize the colormap
colormap(jet); % You can choose from MATLAB's colormaps such as jet, hsv, hot, etc.

% Add colorbar
colorbar;
ylabel(colorbar, 'Value of Z'); % Label the colorbar

% Enhance the plot
title('3D Scatter Plot with Color Coding', 'Interpreter', 'latex', 'FontSize', 16);
xlabel('X-axis', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Y-axis', 'Interpreter', 'latex', 'FontSize', 14);
zlabel('Z-axis', 'Interpreter', 'latex', 'FontSize', 14);
grid on; % Add a grid to make data points easier to interpret
axis tight; % Adjust the axes limits to the data range
view(-30, 10); % Adjust the viewing angle for better visualization

% Set the figure's properties for publication quality
set(gcf, 'Color', 'w'); % Set background color to white
set(gca, 'FontSize', 12); % Adjust font size for readability
set(gca, 'LineWidth', 1.5); % Set the width of the axis lines
set(gca, 'FontName', 'Times New Roman'); % Use a font that is commonly used in publications

% Optionally, save the figure
% saveas(gcf, '3dscatterplot.png'); % Save the figure as a PNG file
% saveas(gcf, '3dscatterplot.eps', 'epsc'); % Save the figure as an EPS file for high-quality printing

hold off;


%Contour Plot with Projection
[X, Y] = meshgrid(linspace(min(x), max(x), 100), linspace(min(y), max(y), 100));
Zq = griddata(x, y, z, X, Y, 'natural');
contourf(X, Y, Zq, 20, 'LineStyle', 'none'); % 20 levels of contours
hold on;
scatter3(x, y, z, 36, z, 'filled');
colorbar;
xlabel('X-axis');
ylabel('Y-axis');
title('Contour Plot with 3D Scatter');
hold off;

%Creating 3-d plot of changing characteristics of daytime heatwaves,
%nighttime heatwaves and droughts

%Getting changing characteristics of Daytime Heatwaves
dayhw = readmatrix('/Users/smonda15/Ongoing Work/Publication_Compound_Drought_Heatwave_Over_Phoneix/Results/Results Southwest North America/SWNA_heatwave_avg_2020_2023_90percentile.xlsx');

dayhw_avg_freq            = dayhw(:,3);
dayhw_2023_freq           = dayhw(:,11);

dayhw_change_freq_2023    = zeros(2256,1);

for i=1:2256

dayhw_change_freq_2023(i,1)    = dayhw_2023_freq(i,1)/(dayhw_avg_freq(i,1)/45);
end


dayhw_avg_duration               = dayhw(:,4);
dayhw_2023_duration              = dayhw(:,12);

dayhw_change_dur_2023            = zeros(2256,1);

for i=1:2256
    dayhw_change_dur_2023(i,1)   = dayhw_2023_duration(i,1)/dayhw_avg_duration(i,1);
end



dayhw_avg_amplitude              = dayhw(:,5);
dayhw_2023_amplitude             = dayhw(:,13);

dayhw_change_amp_2023  = zeros(2256,1);

for i=1:2256
    dayhw_change_amp_2023(i,1)= dayhw_2023_amplitude(i,1)/dayhw_avg_amplitude(i,1);
end



%Getting the changing characteristics of nighttime heatwaves

nighthw = readmatrix('/Users/smonda15/Ongoing Work/Publication_Compound_Drought_Heatwave_Over_Phoneix/Results/Results Southwest North America/SWNA_nightime_heatwave_average_2020_2023_characteristics_90perc.csv');

nighthw_avg_freq        = nighthw(:,3);
nighthw_2023_freq       = nighthw(:,11);

nighthw_change_freq_2023    = zeros(2256,1);

for i=1:2256

nighthw_change_freq_2023(i,1)    = nighthw_2023_freq(i,1)/(nighthw_avg_freq(i,1)/45);
end


nighthw_avg_duration               = nighthw(:,4);
nighthw_2023_duration              = nighthw(:,12);

nighthw_change_dur_2023            = zeros(2256,1);

for i=1:2256
    nighthw_change_dur_2023(i,1)   = nighthw_2023_duration(i,1)/nighthw_avg_duration(i,1);
end

nighthw_avg_amplitude              = nighthw(:,5);
nighthw_2023_amplitude             = nighthw(:,13);

nighthw_change_amp_2023  = zeros(2256,1);

for i=1:2256
    nighthw_change_amp_2023(i,1)= nighthw_2023_amplitude(i,1)/nighthw_avg_amplitude(i,1);
end




%Getting the changing characteristics of Droughts

dr   =  readmatrix('/Users/smonda15/Ongoing Work/Publication_Compound_Drought_Heatwave_Over_Phoneix/Results/Results Southwest North America/SWNA_SPI12_Moderate_2020_2023_Drought_Characteristics_Anomaly.csv');

dr_avg_freq    = dr(:,3)/44;

dr_2023_freq   = dr(:,11);

dr_change_freq_2023 = zeros(2256,1);

for i=1:2256
    dr_change_freq_2023(i,1)=dr_2023_freq(i,1)/dr_avg_freq(i,1);
end


dr_avg_duration               = dr(:,4);
dr_2023_duration              = dr(:,12);

dr_change_dur_2023            = zeros(2256,1);

for i=1:2256
    dr_change_dur_2023(i,1)   = dr_2023_duration(i,1)/dr_avg_duration(i,1);
end


dr_avg_amplitude              = dr(:,6);
dr_2023_amplitude             = dr(:,14);

dr_change_amp_2023  = zeros(2256,1);

for i=1:2256
    dr_change_amp_2023(i,1)= dr_2023_amplitude(i,1)/dr_avg_amplitude(i,1);
end



x = dayhw_change_amp_2023;
y = nighthw_change_amp_2023;
z = dr_change_amp_2023;

% Generate a figure
figure;

% Create a 3D scatter plot
scatter3(x, y, z, 36, z, 'filled'); % 36 is the size of the markers
hold on; % Keep the plot for further modifications

% Customize the colormap
colormap(jet); % You can choose from MATLAB's colormaps such as jet, hsv, hot, etc.

% Add colorbar
colorbar;
ylabel(colorbar, 'Value of Z'); % Label the colorbar

% Enhance the plot
title('3D Scatter Plot with Color Coding', 'Interpreter', 'latex', 'FontSize', 16);
xlabel('X-axis', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Y-axis', 'Interpreter', 'latex', 'FontSize', 14);
zlabel('Z-axis', 'Interpreter', 'latex', 'FontSize', 14);
grid on; % Add a grid to make data points easier to interpret
axis tight; % Adjust the axes limits to the data range
view(-30, 10); % Adjust the viewing angle for better visualization

% Set the figure's properties for publication quality
set(gcf, 'Color', 'w'); % Set background color to white
set(gca, 'FontSize', 12); % Adjust font size for readability
set(gca, 'LineWidth', 1.5); % Set the width of the axis lines
set(gca, 'FontName', 'Times New Roman'); % Use a font that is commonly used in publications

% Optionally, save the figure
% saveas(gcf, '3dscatterplot.png'); % Save the figure as a PNG file
% saveas(gcf, '3dscatterplot.eps', 'epsc'); % Save the figure as an EPS file for high-quality printing

hold off;




%%Calculating Severity of Heatwaves, Droughts and Compound
%%Drought-Heatwaves

lon = ncread("/Users/smonda15/Ongoing Work/Publication_Compound_Drought_Heatwave_Over_Phoneix/Datasets/CPC Global Temperature/tmax.1979.nc",'lon');
lat = ncread("/Users/smonda15/Ongoing Work/Publication_Compound_Drought_Heatwave_Over_Phoneix/Datasets/CPC Global Temperature/tmax.1979.nc",'lat');

folderPath = '/Users/smonda15/Ongoing Work/Publication_Compound_Drought_Heatwave_Over_Phoneix/Datasets/CPC Global Temperature';
ncFiles = dir(fullfile(folderPath, '*.nc'));


b=zeros(259200,1);
d=zeros(259200,1);
c=1979;

for i = 1:length(ncFiles)
    
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
SWheat_75 = zeros(2256,183);
SWheat_25 = zeros(2256,183);


for i=1:2256
    for j=1:183
        SWheat_95(i,j)=prctile(SWNAheat_pool(i,j,:),95);
        SWheat_90(i,j)=prctile(SWNAheat_pool(i,j,:),90);
        SWheat_25(i,j)=prctile(SWNAheat_pool(i,j,:),25);
        SWheat_75(i,j)=prctile(SWNAheat_pool(i,j,:),75);
    end
end

%Calculating Spatiotemporal Data with heatwave occurrence as 1 and not
%occurrence as 0

SWheat_r    =  reshape(SWNAheat,[2256,183,45]);
SWheat_b    =  zeros(2256,183,45);
SWheat_d    =  zeros(2256,183,45);

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
        end
    end
end

SWheat_d1 = (reshape(SWheat_d,[2256,183*45]))';
SWheat_b1 = (reshape(SWheat_b,[2256,183*45]))';



%Calculating Daily Temperature - the corresponding calendar threshold
%(first step of heatwave severity calculation)

SWheat_90_r = SWheat_90';
SWheat_95_r = SWheat_95'; 

SWheat_severity_90 = zeros(8235,2256);
SWheat_severity_95 = zeros(8235,2256);

for i=1:45
    for j=1:2256
        k = 183*(i-1)+1;
        l = 183*i;
        SWheat_severity_90(k:l,j)=SWNAheat(j,k:l)'-(SWheat_90_r(:,j));
        SWheat_severity_95(k:l,j)=SWNAheat(j,k:l)'-(SWheat_95_r(:,j));
    end
end

%Calculating normalizing factor for each day (denomitor of the heatwave
%component in the compound drought and heatwave severity)

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


%Calculate the heatwave component of Compound Drought & Heatwave severity

CDHW_hw_sev = zeros(2256,8235);

for i=1:45
    k = 183*(i-1)+1;
    l = 183*i;
    for j=1:2256

    CDHW_hw_sev(j,k:l) =  SWheat_numerator(j,k:l)./SWheat_normf(j,:);
    end

end




CDHW_hw_sev_1980 = CDHW_hw_sev(:,184:8235);

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

SW_SPI3_summer = zeros(1144,2256);

for i=1:44
    m   = 52*(i-1)+1;
    n   = 52*i;

    p   = 26*(i-1)+1;
    q   = 26*i;
    l   = SW_SPI3(m:n,:);
    SW_SPI3_summer(p:q,:) = l(14:39,:);
end





CDHW_hw_dr3_sev  = zeros(2256,1144);

for i=1:2256
    for j=1:1144

        CDHW_hw_dr3_sev(i,j) = CDHW_hw_sev_rw(i,j)*SW_SPI3_summer(j,i)*(-1);
    end
end


%Calculating Spatiotemporal Characteristics of Compound Drought & Heatwave

Compound_dr_hw    = Compound_dr_hw';
CDHW_hw_dr3_sev   = CDHW_hw_dr3_sev';

CDHW_hw_dr3_sev((find(CDHW_hw_dr3_sev<0)))=0;


numLocations = 2256;
numWeeks     = 26;

numDroughtEvents4 = zeros(1, numLocations);
cDroughtSeverity = zeros(numWeeks, numLocations);
cDroughtDuration = zeros(numWeeks, numLocations);
cDroughtIntensity = zeros(numWeeks, numLocations);




for loc = 1:numLocations
    droughtPeriods = Compound_dr_hw(1119:1144, loc);      %for 2020, 2081:2132 & %2237:2288
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




 




%Calculating Occurrence of Heatwaves as Binary Time series        
heatwaveEvents_95 = findHeatwaves(SWheat_b);


heatwaveEvents_90 = findHeatwaves(SWheat_d);


heatwaveEvents_95 = (reshape(heatwaveEvents_95,[2256,183*45]))';
heatwaveEvents_90 = (reshape(heatwaveEvents_90,[2256,183*45]))';


hw_temp_95        = (heatwaveEvents_95.*SWNAheat'); 
hw_temp_90        = (heatwaveEvents_90.*SWNAheat');


%Calculating Duration of each Heatwaves

heatwaveDuration_95       = DurationHeatwaves(SWheat_b);

heatwaveDuration_90       = DurationHeatwaves(SWheat_d);

heatwaveDuration_95       = (reshape(heatwaveDuration_95,[2256,183*45]))';
heatwaveDuration_90       = (reshape(heatwaveDuration_90,[2256,183*45]))';


%Calculating Frequency of Heatwaves
heatwave_Frequency_95 = zeros(2256,1);
heatwave_Frequency_90 = zeros(2256,1);

for i=1:2256
    a=find(heatwaveDuration_90(:,i)>0);
    a=size(a);
    b=find(heatwaveDuration_95(:,i)>0);
    b=size(b);
    heatwave_Frequency_90(i,1)=a(1,1);
    heatwave_Frequency_95(i,1)=b(1,1);
end

isnan(heatwave_Frequency_90);


%Calculating Duration of Heatwaves

HW_avg_duration_90 =zeros(2256,1);
HW_avg_duration_95 =zeros(2256,1);

for i=1:2256
    a=find(heatwaveDuration_90(:,i)>0);
    a=heatwaveDuration_90(a,i);
    a=mean(a,1);
    b=find(heatwaveDuration_95(:,i)>0);
    b=heatwaveDuration_95(b,i);
    b=mean(b,1);
    HW_avg_duration_90(i,1)=a;
    HW_avg_duration_95(i,1)=b;
end



%Calculating Heatwave Severity


SWheat_hwseverity_90 = zeros(8235,2256);

SWheat_hwseverity_95 = zeros(8235,2256);


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

%Calculating Average Heatwave Severity for  whole summer

SWNAheat_avg_sev_90 = zeros(2256,1);
SWNAheat_avg_sev_95 = zeros(2256,1);

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


%Calculating average heatwave severity for different decades (1979-1989,
%1990-2000, 2001-2011, 2012-2023)


SWNAheat_avg_sev_90_1989 = zeros(2256,1);

SWNAheat_avg_sev_90_2000 = zeros(2256,1);

SWNAheat_avg_sev_90_2011 = zeros(2256,1);

SWNAheat_avg_sev_90_2023 = zeros(2256,1);


SWNAheat_avg_sev_95_1989 = zeros(2256,1);

SWNAheat_avg_sev_95_2000 = zeros(2256,1);

SWNAheat_avg_sev_95_2011 = zeros(2256,1);

SWNAheat_avg_sev_95_2023 = zeros(2256,1);

%Decadal Division would be 1:2013, 2014:4026, 4027:6039, 6040:8236


for i=1:2256
    a = find(heatwaveDuration_90(1:2013,i)>0);
    if a~=0
        b = SWheat_hwseverity_90(a,i);
        SWNAheat_avg_sev_90_1989(i,1) = mean(b);
    else 
        SWNAheat_avg_sev_90_1989(i,1) = 0;
    end
end

for i=1:2256
    a = find(heatwaveDuration_90(2014:4026,i)>0);
    if a~=0
        b = SWheat_hwseverity_90(a+2013,i);
        SWNAheat_avg_sev_90_2000(i,1) = mean(b);
    else 
        SWNAheat_avg_sev_90_2000(i,1) = 0;
    end
end

for i=1:2256
    a = find(heatwaveDuration_90(4027:6039,i)>0);
    if a~=0
        b = SWheat_hwseverity_90(a+4026,i);
        SWNAheat_avg_sev_90_2011(i,1) = mean(b);
    else 
        SWNAheat_avg_sev_90_2011(i,1) = 0;
    end
end

for i=1:2256
    a = find(heatwaveDuration_90(6040:8235,i)>0);
    if a~=0
        b = SWheat_hwseverity_90(a+6039,i);
        SWNAheat_avg_sev_90_2023(i,1) = mean(b);
    else 
        SWNAheat_avg_sev_90_2023(i,1) = 0;
    end
end


for i=1:2256
    a = find(heatwaveDuration_95(1:2013,i)>0);
    if a~=0
        b = SWheat_hwseverity_95(a,i);
        SWNAheat_avg_sev_95_1989(i,1) = mean(b);
    else 
        SWNAheat_avg_sev_95_1989(i,1) = 0;
    end
end

for i=1:2256
    a = find(heatwaveDuration_95(2014:4026,i)>0);
    if a~=0
        b = SWheat_hwseverity_95(a+2013,i);
        SWNAheat_avg_sev_95_2000(i,1) = mean(b);
    else 
        SWNAheat_avg_sev_95_2000(i,1) = 0;
    end
end

for i=1:2256
    a = find(heatwaveDuration_95(4027:6039,i)>0);
    if a~=0
        b = SWheat_hwseverity_95(a+4026,i);
        SWNAheat_avg_sev_95_2011(i,1) = mean(b);
    else 
        SWNAheat_avg_sev_95_2011(i,1) = 0;
    end
end

for i=1:2256
    a = find(heatwaveDuration_95(6040:8235,i)>0);
    if a~=0
        b = SWheat_hwseverity_95(a+6039,i);
        SWNAheat_avg_sev_95_2023(i,1) = mean(b);
    else 
        SWNAheat_avg_sev_95_2023(i,1) = 0;
    end
end

%Calculating Heatwave average severity for April-May-June,
%July-August-September respectively
 

%First, we have to seperate severity values for April-May-June &
%July-August-September

SWheat_hwseverity_95_Apr_June = zeros(4095,2256);
SWheat_hwseverity_95_Jul_Sep  = zeros(4140,2256);


SWheat_hwseverity_90_Apr_June = zeros(4095,2256);
SWheat_hwseverity_90_Jul_Sep  = zeros(4140,2256);


for i=1:45
    p = 183*(i-1)+1;
    q = 183*i;
    r = 91*(i-1)+1;
    s = 91*i;
    t = 92*(i-1)+1;
    v = 92*i;

    m = SWheat_hwseverity_95(p:q,:);
    n = SWheat_hwseverity_90(p:q,:);

    SWheat_hwseverity_95_Apr_June(r:s,:) = m(1:91,:);
    SWheat_hwseverity_95_Jul_Sep(t:v,:)  = m(92:183,:);

    SWheat_hwseverity_90_Apr_June(r:s,:) = n(1:91,:);
    SWheat_hwseverity_90_Jul_Sep(t:v,:)  = n(92:183,:);
end


SWNAheat_avg_sev_90_Apr_June = zeros(2256,1);
SWNAheat_avg_sev_90_Jul_Sep  = zeros(2256,1);
SWNAheat_avg_sev_95_Apr_June = zeros(2256,1);
SWNAheat_avg_sev_95_Jul_Sep = zeros(2256,1);

for i=1:2256
    a = find(SWheat_hwseverity_95_Apr_June(:,i)>0);
    if a~=0
        b = SWheat_hwseverity_95_Apr_June(a,i);
        SWNAheat_avg_sev_95_Apr_June(i,1) = mean(b);
    else 
        SWNAheat_avg_sev_95_Apr_June(i,1) = 0;
    end
end

for i=1:2256
    a = find(SWheat_hwseverity_95_Jul_Sep(:,i)>0);
    if a~=0
        b = SWheat_hwseverity_95_Jul_Sep(a,i);
        SWNAheat_avg_sev_95_Jul_Sep(i,1) = mean(b);
    else 
        SWNAheat_avg_sev_95_Jul_Sep(i,1) = 0;
    end
end

for i=1:2256
    a = find(SWheat_hwseverity_90_Apr_June(:,i)>0);
    if a~=0
        b = SWheat_hwseverity_90_Apr_June(a,i);
        SWNAheat_avg_sev_90_Apr_June(i,1) = mean(b);
    else 
        SWNAheat_avg_sev_90_Apr_June(i,1) = 0;
    end
end

for i=1:2256
    a = find(SWheat_hwseverity_90_Jul_Sep(:,i)>0);
    if a~=0
        b = SWheat_hwseverity_90_Jul_Sep(a,i);
        SWNAheat_avg_sev_90_Jul_Sep(i,1) = mean(b);
    else 
        SWNAheat_avg_sev_90_Jul_Sep(i,1) = 0;
    end
end




%Calculating average severity for 2023

SWNAheat_avg_sev_90_2023 = zeros(2256,1);
SWNAheat_avg_sev_95_2023 = zeros(2256,1);

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


%Calculating average severity for 2023 for April-May-June,
%July-August-September

SWNAheat_avg_sev_90_2023_Apr_June = zeros(2256,1);
SWNAheat_avg_sev_90_2023_Jul_Sep  = zeros(2256,1);
SWNAheat_avg_sev_95_2023_Apr_June = zeros(2256,1);
SWNAheat_avg_sev_95_2023_Jul_Sep  = zeros(2256,1);

for i=1:2256
    a = find(heatwaveDuration_90(8053:8143,i)>0);
    if a~=0
        b = SWheat_hwseverity_90(a+8052,i);
        SWNAheat_avg_sev_90_2023_Apr_June(i,1) = mean(b);
    else 
        SWNAheat_avg_sev_90_2023_Apr_June(i,1) = 0;
    end
end

for i=1:2256
    a = find(heatwaveDuration_90(8144:8235,i)>0);
    if a~=0
        b = SWheat_hwseverity_90(a+8143,i);
        SWNAheat_avg_sev_90_2023_Jul_Sep(i,1) = mean(b);
    else 
        SWNAheat_avg_sev_90_2023_Jul_Sep(i,1) = 0;
    end
end

for i=1:2256
    a = find(heatwaveDuration_95(8053:8143,i)>0);
    if a~=0
        b = SWheat_hwseverity_95(a+8052,i);
        SWNAheat_avg_sev_95_2023_Apr_June(i,1) = mean(b);
    else 
        SWNAheat_avg_sev_95_2023_Apr_June(i,1) = 0;
    end
end

for i=1:2256
    a = find(heatwaveDuration_95(8144:8235,i)>0);
    if a~=0
        b = SWheat_hwseverity_95(a+8143,i);
        SWNAheat_avg_sev_95_2023_Jul_Sep(i,1) = mean(b);
    else 
        SWNAheat_avg_sev_95_2023_Jul_Sep(i,1) = 0;
    end
end


%Considering Summer Droughts, Both initiating and terminating droughts


%Calculating Precipitation Anomalies for April-May-June and
%July-August-September

%%Extracting CPC Precipitation data
lon = ncread("/Users/smonda15/Ongoing Work/Publication_Compound_Drought_Heatwave_Over_Phoneix/Datasets/CPC Global Precipitation/precip.1979.nc",'lon');
lat = ncread("/Users/smonda15/Ongoing Work/Publication_Compound_Drought_Heatwave_Over_Phoneix/Datasets/CPC Global Precipitation/precip.1979.nc",'lat');

folderPath = '/Users/smonda15/Ongoing Work/Publication_Compound_Drought_Heatwave_Over_Phoneix/Datasets/CPC Global Precipitation';
ncFiles = dir(fullfile(folderPath, '*.nc'));


b=zeros(259200,1);
d=zeros(259200,1);

c=1979;

for i = 1:45
    
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



SW_SPI12     =   reshape(SW_SPI12,[52,44,2256]);
SW_SPI3      =   reshape(SW_SPI3,[52,44,2256]);

%Extracting SPI3 for summer 

%SW_SPI3     =   SW_SPI3(14:39,:,:);
%SW_SPI3     =   reshape(SW_SPI3,[26*44,2256]);

%Using SPI=-1 as Drought Threshold moderate drought and -1.5 for extreme
%drought

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


%SW_drought = SW_drought(531:543,:);

% Assuming your matrix is called 'droughtData' and is of size 585x439
numWeeks = 2288;
numLocations = 2256;
%Corresponding to Drought Initiation
% Initialize the matrices to store results
numDroughtEvents = zeros(1, numLocations);
DroughtSeverity = zeros(numWeeks, numLocations);
DroughtDuration = zeros(numWeeks, numLocations);
DroughtIntensity = zeros(numWeeks, numLocations);




for loc = 1:numLocations
    droughtPeriods = SW_drought12(:, loc);      %for 2020, 2081:2132 & %2237:2288
    droughtspi     = SW_SPI12(:,loc);
    % ... [any additional code for specific periods] ...

    % Find the start and end indices of each drought event
    startIndices = find(diff([0; droughtPeriods]) == 1);
    endIndices = find(diff([droughtPeriods; 0]) == -1);

    % Initialize the number of events for this location
    numEvents = 0;
    
    for i = 1:length(startIndices)
        eventDuration = endIndices(i) - startIndices(i) + 1;
        
        % Check if the event duration is at least two weeks
        if eventDuration >= 2
            numEvents = numEvents + 1; % Counting the valid event
            DroughtDuration(startIndices(i), loc) = eventDuration;
            
            eventSeverity = sum(droughtspi(startIndices(i):endIndices(i)));
            eventIntensity = eventSeverity / eventDuration;

            DroughtSeverity(startIndices(i), loc) = eventSeverity;
            DroughtIntensity(startIndices(i), loc) = eventIntensity;
        end
    end

    % Store the number of valid drought events for this location
    numDroughtEvents(loc) = numEvents;
end


%Extracting Droughts initiating only in summer

DroughtDuration_Summer_Init    = zeros(1144,2256);
DroughtIntensity_Summer_Init   = zeros(1144,2256);
DroughtSeverity_Summer_Init    = zeros(1144,2256); 

for i=1:44
    k   =  52*(i-1)+1;
    l   =  52*i;

    m   =  26*(i-1)+1;
    n   =  26*i;

    a   =  DroughtDuration(k:l,:);
    b   =  DroughtIntensity(k:l,:);
    c   =  DroughtSeverity(k:l,:);

    DroughtDuration_Summer_Init(m:n,:)   = a(14:39,:);
    DroughtIntensity_Summer_Init(m:n,:)  = b(14:39,:);
    DroughtSeverity_Summer_Init(m:n,:)   = c(14:39,:);
end

%Calculating Average Drought Duration, Severity and Intensity
numDroughtEvents         =  zeros(2256,1);
SW_dr_avg_duration       =  zeros(2256,1);
SW_dr_avg_severity       =  zeros(2256,1);
SW_dr_avg_intensity      =  zeros(2256,1);

for i=1:2256
    a  =  find(DroughtDuration_Summer_Init(:,i)>0);
    numDroughtEvents(i,1)= size(a,1);
    b  =  DroughtDuration_Summer_Init(a,i);
    b  =  mean(b,1);
    SW_dr_avg_duration(i,1)=b;
    c  =  DroughtSeverity_Summer_Init(a,i);
    c  =  mean(c,1);
    SW_dr_avg_severity(i,1)=c;
    d  =  DroughtIntensity_Summer_Init(a,i);
    d  =  mean(d,1);
    SW_dr_avg_intensity(i,1)=d;
end

    
a =find(isnan(SW_dr_avg_duration));

SW_dr_avg_duration(a,:) = 0;

SW_dr_avg_intensity(a,:) = 0;

SW_dr_avg_severity(a,:) = 0;

%Calculating Decadal change in spatiotemporal characteristics of initiating
%drought

%Calculating Average Drought Duration, Severity and Intensity

%1980-1990
numDroughtEvents1         =  zeros(2256,1);
SW_dr_avg_duration1       =  zeros(2256,1);
SW_dr_avg_severity1       =  zeros(2256,1);
SW_dr_avg_intensity1      =  zeros(2256,1);

for i=1:2256
    a  =  find(DroughtDuration_Summer_Init(1:286,i)>0);
    numDroughtEvents1(i,1)= size(a,1);
    b  =  DroughtDuration_Summer_Init(a,i);
    b  =  mean(b,1);
    SW_dr_avg_duration1(i,1)=b;
    c  =  DroughtSeverity_Summer_Init(a,i);
    c  =  mean(c,1);
    SW_dr_avg_severity1(i,1)=c;
    d  =  DroughtIntensity_Summer_Init(a,i);
    d  =  mean(d,1);
    SW_dr_avg_intensity1(i,1)=d;
end

    
a =find(isnan(SW_dr_avg_duration1));

SW_dr_avg_duration1(a,:) = 0;

SW_dr_avg_intensity1(a,:) = 0;

SW_dr_avg_severity1(a,:) = 0;


%1991-2001
numDroughtEvents2         =  zeros(2256,1);
SW_dr_avg_duration2       =  zeros(2256,1);
SW_dr_avg_severity2       =  zeros(2256,1);
SW_dr_avg_intensity2      =  zeros(2256,1);

for i=1:2256
    a  =  find(DroughtDuration_Summer_Init(287:572,i)>0);
    numDroughtEvents2(i,1)= size(a,1);
    b  =  DroughtDuration_Summer_Init(a+286,i);
    b  =  mean(b,1);
    SW_dr_avg_duration2(i,1)=b;
    c  =  DroughtSeverity_Summer_Init(a+286,i);
    c  =  mean(c,1);
    SW_dr_avg_severity2(i,1)=c;
    d  =  DroughtIntensity_Summer_Init(a+286,i);
    d  =  mean(d,1);
    SW_dr_avg_intensity2(i,1)=d;
end

    
a =find(isnan(SW_dr_avg_duration2));

SW_dr_avg_duration2(a,:) = 0;

SW_dr_avg_intensity2(a,:) = 0;

SW_dr_avg_severity2(a,:) = 0;


%2002-2012
numDroughtEvents3         =  zeros(2256,1);
SW_dr_avg_duration3       =  zeros(2256,1);
SW_dr_avg_severity3       =  zeros(2256,1);
SW_dr_avg_intensity3      =  zeros(2256,1);

for i=1:2256
    a  =  find(DroughtDuration_Summer_Init(573:858,i)>0);
    numDroughtEvents3(i,1)= size(a,1);
    b  =  DroughtDuration_Summer_Init(a+572,i);
    b  =  mean(b,1);
    SW_dr_avg_duration3(i,1)=b;
    c  =  DroughtSeverity_Summer_Init(a+572,i);
    c  =  mean(c,1);
    SW_dr_avg_severity3(i,1)=c;
    d  =  DroughtIntensity_Summer_Init(a+572,i);
    d  =  mean(d,1);
    SW_dr_avg_intensity3(i,1)=d;
end

    
a =find(isnan(SW_dr_avg_duration3));

SW_dr_avg_duration3(a,:) = 0;

SW_dr_avg_intensity3(a,:) = 0;

SW_dr_avg_severity3(a,:) = 0;

%2013-2023
numDroughtEvents4         =  zeros(2256,1);
SW_dr_avg_duration4       =  zeros(2256,1);
SW_dr_avg_severity4       =  zeros(2256,1);
SW_dr_avg_intensity4      =  zeros(2256,1);

for i=1:2256
    a  =  find(DroughtDuration_Summer_Init(859:1144,i)>0);
    numDroughtEvents4(i,1)= size(a,1);
    b  =  DroughtDuration_Summer_Init(a+858,i);
    b  =  mean(b,1);
    SW_dr_avg_duration4(i,1)=b;
    c  =  DroughtSeverity_Summer_Init(a+858,i);
    c  =  mean(c,1);
    SW_dr_avg_severity4(i,1)=c;
    d  =  DroughtIntensity_Summer_Init(a+858,i);
    d  =  mean(d,1);
    SW_dr_avg_intensity4(i,1)=d;
end

    
a =find(isnan(SW_dr_avg_duration4));

SW_dr_avg_duration4(a,:) = 0;

SW_dr_avg_intensity4(a,:) = 0;

SW_dr_avg_severity4(a,:) = 0;




%Extracting Droughts initiating only in summer,2023

%this should be done if only initiating droughts need to be considered

%DroughtDuration_Summer_Init_2023    = DroughtDuration_Summer_Init(1119:1144,:);
%DroughtIntensity_Summer_Init_2023   = DroughtIntensity_Summer_Init(1119:1144,:);
%DroughtSeverity_Summer_Init_2023    = DroughtSeverity_Summer_Init(1119:1144,:); 

numWeeks = 39;
numLocations = 2256;
%Corresponding to Drought Initiation
% Initialize the matrices to store results
numDroughtEvents = zeros(1, numLocations);
DroughtSeverity = zeros(numWeeks, numLocations);
DroughtDuration = zeros(numWeeks, numLocations);
DroughtIntensity = zeros(numWeeks, numLocations);




for loc = 1:numLocations
    droughtPeriods = SW_drought12(2250:2288, loc);      %for 2020, 2081:2132 & %2237:2288
    droughtspi     = SW_SPI12(2250:2288,loc);
    % ... [any additional code for specific periods] ...

    % Find the start and end indices of each drought event
    startIndices = find(diff([0; droughtPeriods]) == 1);
    endIndices = find(diff([droughtPeriods; 0]) == -1);

    % Initialize the number of events for this location
    numEvents = 0;
    
    for i = 1:length(startIndices)
        eventDuration = endIndices(i) - startIndices(i) + 1;
        
        % Check if the event duration is at least two weeks
        if eventDuration >= 2
            numEvents = numEvents + 1; % Counting the valid event
            DroughtDuration(startIndices(i), loc) = eventDuration;
            
            eventSeverity = sum(droughtspi(startIndices(i):endIndices(i)));
            eventIntensity = eventSeverity / eventDuration;

            DroughtSeverity(startIndices(i), loc) = eventSeverity;
            DroughtIntensity(startIndices(i), loc) = eventIntensity;
        end
    end

    % Store the number of valid drought events for this location
    numDroughtEvents(loc) = numEvents;
end

DroughtDuration_Summer_Init_2023    = DroughtDuration(1:26,:);
DroughtIntensity_Summer_Init_2023   = DroughtIntensity(1:26,:);
DroughtSeverity_Summer_Init_2023    = DroughtSeverity(1:26,:); 


%Calculating frequency, Duration, Severity and Intensity of initiating
%summer droughts in 2023
numDroughtEvents         =  zeros(2256,1);
SW_dr_avg_duration       =  zeros(2256,1);
SW_dr_avg_severity       =  zeros(2256,1);
SW_dr_avg_intensity      =  zeros(2256,1);

for i=1:2256
    a  =  find(DroughtDuration_Summer_Init_2023(:,i)>0);
    numDroughtEvents(i,1)= size(a,1);
    b  =  DroughtDuration_Summer_Init_2023(a,i);
    b  =  mean(b,1);
    SW_dr_avg_duration(i,1)=b;
    c  =  DroughtSeverity_Summer_Init_2023(a,i);
    c  =  mean(c,1);
    SW_dr_avg_severity(i,1)=c;
    d  =  DroughtIntensity_Summer_Init_2023 (a,i);
    d  =  mean(d,1);
    SW_dr_avg_intensity(i,1)=d;
end

    
a =find(isnan(SW_dr_avg_duration));

SW_dr_avg_duration(a,:) = 0;

SW_dr_avg_intensity(a,:) = 0;

SW_dr_avg_severity(a,:) = 0;



%Corresponding to Terminating Droughts

numDroughtEvents_T = zeros(1, numLocations);
DroughtSeverity_T = zeros(numWeeks, numLocations);
DroughtDuration_T = zeros(numWeeks, numLocations);
DroughtIntensity_T = zeros(numWeeks, numLocations);


for loc = 1:numLocations
    droughtPeriods = SW_drought12(:, loc);      %for 2020, 2081:2132 & %2237:2288
    droughtspi     = SW_SPI12(:,loc);
    % ... [any additional code for specific periods] ...

    % Find the start and end indices of each drought event
    startIndices = find(diff([0; droughtPeriods]) == 1);
    endIndices = find(diff([droughtPeriods; 0]) == -1);

    % Initialize the number of events for this location
    numEvents = 0;
    
    for i = 1:length(startIndices)
        eventDuration = endIndices(i) - startIndices(i) + 1;
        
        % Check if the event duration is at least two weeks
        if eventDuration >= 2
            numEvents = numEvents + 1; % Counting the valid event
            DroughtDuration_T(endIndices(i), loc) = eventDuration;
            
            eventSeverity = sum(droughtspi(startIndices(i):endIndices(i)));
            eventIntensity = eventSeverity / eventDuration;

            DroughtSeverity_T(endIndices(i), loc) = eventSeverity;
            DroughtIntensity_T(endIndices(i), loc) = eventIntensity;
        end
    end

    % Store the number of valid drought events for this location
    numDroughtEvents_T(loc) = numEvents;
end

%Extracting Droughts Terminating only in summer

DroughtDuration_Summer_Term    = zeros(1144,2256);
DroughtIntensity_Summer_Term   = zeros(1144,2256);
DroughtSeverity_Summer_Term    = zeros(1144,2256); 

for i=1:44
    k   =  52*(i-1)+1;
    l   =  52*i;

    m   =  26*(i-1)+1;
    n   =  26*i;

    a   =  DroughtDuration_T(k:l,:);
    b   =  DroughtIntensity_T(k:l,:);
    c   =  DroughtSeverity_T(k:l,:);

    DroughtDuration_Summer_Term(m:n,:)   = a(14:39,:);
    DroughtIntensity_Summer_Term(m:n,:)  = b(14:39,:);
    DroughtSeverity_Summer_Term(m:n,:)   = c(14:39,:);
end

%Calculating Average Drought Duration, Severity and Intensity
numDroughtEvents_T         =  zeros(2256,1);
SW_dr_avg_duration_T       =  zeros(2256,1);
SW_dr_avg_severity_T       =  zeros(2256,1);
SW_dr_avg_intensity_T      =  zeros(2256,1);

for i=1:2256
    a  =  find(DroughtDuration_Summer_Term(:,i)>0);
    numDroughtEvents_T(i,1)=size(a,1);
    b  =  DroughtDuration_Summer_Term(a,i);
    b  =  mean(b,1);
    SW_dr_avg_duration_T(i,1)=b;
    c  =  DroughtSeverity_Summer_Term(a,i);
    c  =  mean(c,1);
    SW_dr_avg_severity_T(i,1)=c;
    d  =  DroughtIntensity_Summer_Term(a,i);
    d  =  mean(d,1);
    SW_dr_avg_intensity_T(i,1)=d;
end

    
a =find(isnan(SW_dr_avg_duration_T));

SW_dr_avg_duration_T(a,:) = 0;

SW_dr_avg_intensity_T(a,:) = 0;

SW_dr_avg_severity_T(a,:) = 0;

%Calculating change in spatiotemporal characteristics of terminating
%droughts

%1980-1990
numDroughtEvents_T1         =  zeros(2256,1);
SW_dr_avg_duration_T1       =  zeros(2256,1);
SW_dr_avg_severity_T1       =  zeros(2256,1);
SW_dr_avg_intensity_T1      =  zeros(2256,1);

for i=1:2256
    a  =  find(DroughtDuration_Summer_Term(1:286,i)>0);
    numDroughtEvents_T1(i,1)=size(a,1);
    b  =  DroughtDuration_Summer_Term(a,i);
    b  =  mean(b,1);
    SW_dr_avg_duration_T1(i,1)=b;
    c  =  DroughtSeverity_Summer_Term(a,i);
    c  =  mean(c,1);
    SW_dr_avg_severity_T1(i,1)=c;
    d  =  DroughtIntensity_Summer_Term(a,i);
    d  =  mean(d,1);
    SW_dr_avg_intensity_T1(i,1)=d;
end

    
a =find(isnan(SW_dr_avg_duration_T1));

SW_dr_avg_duration_T1(a,:) = 0;

SW_dr_avg_intensity_T1(a,:) = 0;

SW_dr_avg_severity_T1(a,:) = 0;



%1991-2001
numDroughtEvents_T2         =  zeros(2256,1);
SW_dr_avg_duration_T2       =  zeros(2256,1);
SW_dr_avg_severity_T2       =  zeros(2256,1);
SW_dr_avg_intensity_T2      =  zeros(2256,1);

for i=1:2256
    a  =  find(DroughtDuration_Summer_Term(287:572,i)>0);
    numDroughtEvents_T2(i,1)=size(a,1);
    b  =  DroughtDuration_Summer_Term(a+286,i);
    b  =  mean(b,1);
    SW_dr_avg_duration_T2(i,1)=b;
    c  =  DroughtSeverity_Summer_Term(a+286,i);
    c  =  mean(c,1);
    SW_dr_avg_severity_T2(i,1)=c;
    d  =  DroughtIntensity_Summer_Term(a+286,i);
    d  =  mean(d,1);
    SW_dr_avg_intensity_T2(i,1)=d;
end

    
a =find(isnan(SW_dr_avg_duration_T2));

SW_dr_avg_duration_T2(a,:) = 0;

SW_dr_avg_intensity_T2(a,:) = 0;

SW_dr_avg_severity_T2(a,:) = 0;


%2002-2012
numDroughtEvents_T3         =  zeros(2256,1);
SW_dr_avg_duration_T3       =  zeros(2256,1);
SW_dr_avg_severity_T3       =  zeros(2256,1);
SW_dr_avg_intensity_T3      =  zeros(2256,1);

for i=1:2256
    a  =  find(DroughtDuration_Summer_Term(573:858,i)>0);
    numDroughtEvents_T3(i,1)=size(a,1);
    b  =  DroughtDuration_Summer_Term(a+572,i);
    b  =  mean(b,1);
    SW_dr_avg_duration_T3(i,1)=b;
    c  =  DroughtSeverity_Summer_Term(a+572,i);
    c  =  mean(c,1);
    SW_dr_avg_severity_T3(i,1)=c;
    d  =  DroughtIntensity_Summer_Term(a+572,i);
    d  =  mean(d,1);
    SW_dr_avg_intensity_T3(i,1)=d;
end

    
a =find(isnan(SW_dr_avg_duration_T3));

SW_dr_avg_duration_T3(a,:) = 0;

SW_dr_avg_intensity_T3(a,:) = 0;

SW_dr_avg_severity_T3(a,:) = 0;


%2013-2023
numDroughtEvents_T4         =  zeros(2256,1);
SW_dr_avg_duration_T4       =  zeros(2256,1);
SW_dr_avg_severity_T4       =  zeros(2256,1);
SW_dr_avg_intensity_T4      =  zeros(2256,1);

for i=1:2256
    a  =  find(DroughtDuration_Summer_Term(859:1144,i)>0);
    numDroughtEvents_T4(i,1)=size(a,1);
    b  =  DroughtDuration_Summer_Term(a+858,i);
    b  =  mean(b,1);
    SW_dr_avg_duration_T4(i,1)=b;
    c  =  DroughtSeverity_Summer_Term(a+858,i);
    c  =  mean(c,1);
    SW_dr_avg_severity_T4(i,1)=c;
    d  =  DroughtIntensity_Summer_Term(a+858,i);
    d  =  mean(d,1);
    SW_dr_avg_intensity_T4(i,1)=d;
end

    
a =find(isnan(SW_dr_avg_duration_T4));

SW_dr_avg_duration_T4(a,:) = 0;

SW_dr_avg_intensity_T4(a,:) = 0;

SW_dr_avg_severity_T4(a,:) = 0;





%Calculating frequency, Duration, Severity and Intensity of Terminating
%summer droughts in 2023, Here we are also considering droughts which are continuing
%through whole summer by considering the  

%DroughtDuration_Summer_Term_2023    = DroughtDuration_Summer_Term(1119:1144,:);
%DroughtIntensity_Summer_Term_2023   = DroughtIntensity_Summer_Term(1119:1144,:);
%DroughtSeverity_Summer_Term_2023    = DroughtSeverity_Summer_Term(1119:1144,:); 

% Assuming your matrix is called 'droughtData' and is of size 585x439
numWeeks = 2275;
numLocations = 2256;

numDroughtEvents_T = zeros(1, numLocations);
DroughtSeverity_T = zeros(numWeeks, numLocations);
DroughtDuration_T = zeros(numWeeks, numLocations);
DroughtIntensity_T = zeros(numWeeks, numLocations);


for loc = 1:numLocations
    droughtPeriods = SW_drought12(1:2275, loc);      %for 2020, 2081:2132 & %2237:2288
    droughtspi     = SW_SPI12(1:2275,loc);
    % ... [any additional code for specific periods] ...

    % Find the start and end indices of each drought event
    startIndices = find(diff([0; droughtPeriods]) == 1);
    endIndices = find(diff([droughtPeriods; 0]) == -1);

    % Initialize the number of events for this location
    numEvents = 0;
    
    for i = 1:length(startIndices)
        eventDuration = endIndices(i) - startIndices(i) + 1;
        
        % Check if the event duration is at least two weeks
        if eventDuration >= 2
            numEvents = numEvents + 1; % Counting the valid event
            DroughtDuration_T(endIndices(i), loc) = eventDuration;
            
            eventSeverity = sum(droughtspi(startIndices(i):endIndices(i)));
            eventIntensity = eventSeverity / eventDuration;

            DroughtSeverity_T(endIndices(i), loc) = eventSeverity;
            DroughtIntensity_T(endIndices(i), loc) = eventIntensity;
        end
    end

    % Store the number of valid drought events for this location
    numDroughtEvents_T(loc) = numEvents;
end


DroughtDuration_Summer_Term_2023     =   DroughtDuration_T(2250:2275,:);
DroughtIntensity_Summer_Term_2023    =   DroughtIntensity_T(2250:2275,:);
DroughtSeverity_Summer_Term_2023     =   DroughtSeverity_T(2250:2275,:);

numDroughtEvents         =  zeros(2256,1);
SW_dr_avg_duration       =  zeros(2256,1);
SW_dr_avg_severity       =  zeros(2256,1);
SW_dr_avg_intensity      =  zeros(2256,1);

for i=1:2256
    a  =  find(DroughtDuration_Summer_Term_2023(:,i)>0);
    numDroughtEvents(i,1)= size(a,1);
    b  =  DroughtDuration_Summer_Term_2023(a,i);
    b  =  mean(b,1);
    SW_dr_avg_duration(i,1)=b;
    c  =  DroughtSeverity_Summer_Term_2023(a,i);
    c  =  mean(c,1);
    SW_dr_avg_severity(i,1)=c;
    d  =  DroughtIntensity_Summer_Term_2023(a,i);
    d  =  mean(d,1);
    SW_dr_avg_intensity(i,1)=d;
end

    
a =find(isnan(SW_dr_avg_duration));

SW_dr_avg_duration(a,:) = 0;

SW_dr_avg_intensity(a,:) = 0;

SW_dr_avg_severity(a,:) = 0;





%For each year's summer (26 weeks from April 1 to September 30), for each year's grid location 
% I want to extract those drought events which has not started in summer but terminate in summer. Can you write the appropriate matlab code




