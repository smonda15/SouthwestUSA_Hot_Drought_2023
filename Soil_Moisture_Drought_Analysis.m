% Read data from the source file (file to be regridded)
ncdisp('/Volumes/Expansion/Gleam_Data/SMs_1980_GLEAM_v4.1a.nc');

lon = ncread('/Volumes/Expansion/Gleam_Data/SMs_1980_GLEAM_v4.1a.nc','lon');
lat = ncread('/Volumes/Expansion/Gleam_Data/SMs_1980_GLEAM_v4.1a.nc','lat');

%First subsetting for western USA
%500:800, 40.05N-10.05N
%445:895, -135.55 to -90.55

SMs_1980   = ncread('/Volumes/Expansion/Gleam_Data/SMs_1980_GLEAM_v4.1a.nc','SMs');
SMs_1980w  = SMs_1980(445:895,445:800,91:273);%considering 1980 is a leap year

lonw    = lon(445:895,1);
latw    = lat(445:800,1);

V= zeros(451*356,2);
for i=1:356
    V(((451*(i-1)+1):(451*i)),1)=lonw(1:451,1);
    V(((451*(i-1)+1):(451*i)),2)=latw(i,1);
end

SMs_1980w      = reshape(SMs_1980w,[451*356,183]);

SMs_1980w(isnan(SMs_1980w)) = 0;
SMs_1980w_sum = sum(SMs_1980w,2);


nullout  = find(SMs_1980w_sum~=0);

SMs_1980n   = SMs_1980w(nullout,:);

% Define your target years
years = 1980:2023;

% Preallocate a cell to store each year's data
SMs_all_years = cell(numel(years),1);

for i = 1:numel(years)
    yr = years(i);
    isLeap = ( (mod(yr,4) == 0 && mod(yr,100) ~= 0) || (mod(yr,400) == 0) );
    if isLeap
        dayrange = 92:274;   
    else
        dayrange = 91:273;  
    end
    filename = sprintf('/Volumes/Expansion/Gleam_Data/SMs_%d_GLEAM_v4.1a.nc', yr);
    SMs_data = ncread(filename, 'SMs');
    SMs_data = SMs_data(445:895, 445:800, dayrange);
    SMs_data(isnan(SMs_data)) = 0;
    SMs_data = reshape(SMs_data, [451*356, size(SMs_data,3)]);
    SMs_data     = SMs_data(nullout, :);
    
    
    SMs_all_years{i} = SMs_data;
    i
end

all_data_concat = cat(2, SMs_all_years{:});
 

all_data_yearly           =      reshape(all_data_concat,[63697,183,44]);
%all_data_yearly_az_2023   =      all_data_yearly(idxArizona,:,44);

%V_az                      =      V_nullout(idxArizona,:);
      

all_data_1991_2020 = all_data_yearly(:,:,12:41);

all_data_ltm       = mean(all_data_1991_2020,3);

all_data_anm       = zeros(63697,8052);

for i=1:44
    for j=1:63697
    m =183*(i-1)+1;
    n =183*i;
    all_data_anm(j,m:n)= all_data_concat(j,m:n)-all_data_ltm(j,:);

    end
    i
end

all_data_anmw = reshape(all_data_anm,[63697,183,44]);
all_data_anmw = all_data_anmw(:,1:182,:);
%Extracting soil moisture weekly scale for all years
all_data_anm1          = reshape(all_data_anmw,[63697,182*44]);
all_data_anm1          = reshape(all_data_anm1,[63697,1144,7]);
all_data_anm1          = mean(all_data_anm1,3);
all_data_anm1_mx       = all_data_anm1(idxMexico,:);
all_data_anm1_mx_mean  = (mean(all_data_anm1_mx,1))';
all_data_anm1_mx_mean_2023 = all_data_anm1_mx_mean(1119:1144,1);

all_data_anm1_mx_mean_y = mean(reshape(all_data_anm1_mx_mean,[44,26]),2);

all_data_anm1_mx_mean_y1 = mean(reshape(all_data_anm1_mx_mean_y,[11,4]),2);

all_data_anm1_az       = all_data_anm1(idxArizona,:);
all_data_anm1_az_mean  = (mean(all_data_anm1_az,1))';
all_data_anm1_az_mean_2023 = all_data_anm1_az_mean(1119:1144,1);


all_data_anm1_nm       = all_data_anm1(idxNM,:);
all_data_anm1_nm_mean  = (mean(all_data_anm1_nm,1))';
all_data_anm1_nm_mean_2023 = all_data_anm1_nm_mean(1119:1144,1);


%continuation before extracting soil moisture in weekly scale for all years
all_data_anmw = all_data_anmw(:,:,44);

all_data_anmw = reshape(all_data_anmw,[63697,26,7]);
all_data_anmw = mean(all_data_anmw,3);

CDHWstart     = mean(all_data_anmw(:,12:14),2);
CDHWmid       = mean(all_data_anmw(:,18:20),2);
CDHWend       = mean(all_data_anmw(:,24:26),2);

V_nullout     = V(nullout,:);

filename      = '/Volumes/Expansion/Publication_Compound_Drought_Heatwave_Over_Phoneix/Datasets/Shapefiles/Southwest_North_America.shp';
S = shaperead(filename, 'UseGeoCoords',true);

latPoints = V_nullout(:,2);
lonPoints = V_nullout(:,1);
nPoints = size(V_nullout,1);
pointRegionName = repmat({''},nPoints,1);

for i = 1:numel(S)
    % Extract the polygon’s boundary
    polyLat = S(i).Lat;
    polyLon = S(i).Lon;

    % For the polygon’s attribute (e.g., S(i).NAME or S(i).STATE_NAME),
    % use the correct field name from your shapefile structure
    thisName = S(i).NAME;  % or S(i).Name, S(i).STATE, etc.

    % Check which points fall inside this polygon
    in = inpolygon( ...
              lonPoints, ...  % Xs
              latPoints, ...  % Ys
              polyLon,  ...   % polygon X boundary
              polyLat);       % polygon Y boundary

    % Assign the polygon’s name to those points
    pointRegionName(in) = {thisName};
    i
end



idxArizona = find(strcmp(pointRegionName,'Arizona'));
idxMexico  = find(strcmp(pointRegionName,'Mexico'));
idxNM      = find(strcmp(pointRegionName,'New Mexico'));
idxCO      = find(strcmp(pointRegionName,'Colorado'));
idxUT      = find(strcmp(pointRegionName,'Utah'));




all_data_anmw_az = all_data_anmw(idxArizona,:);
all_data_anmw_mx = all_data_anmw(idxMexico,:);
all_data_anmw_nm = all_data_anmw(idxNM,:);
all_data_anmw_co = all_data_anmw(idxCO,:);
all_data_anmw_ut = all_data_anmw(idxUT,:);



all_data_anmw_azm = (mean(all_data_anmw_az,1))';
all_data_anmw_mxm = (mean(all_data_anmw_mx,1))';
all_data_anmw_nmm = (mean(all_data_anmw_nm,1))';
all_data_anmw_com = (mean(all_data_anmw_co,1))';
all_data_anmw_utm = (mean(all_data_anmw_ut,1))';


%Preparing the graphs for Soil Moisture, Compound Drought and Heatwave
%Severity
%Calculating the Difference in Soil Moisture

all_data_anmw_mxm_diff = diff(all_data_anmw_mxm);
all_data_anmw_mxm_diff = cumsum(all_data_anmw_mxm_diff);

all_data_anmw_azm_diff = cumsum(diff(all_data_anmw_azm));
all_data_anmw_nmm_diff = cumsum(diff(all_data_anmw_nmm));


% Now, 'all_data_concat' is [ (451*301) x total_days_across_years ]
%   where columns are consecutive days from 1980 to 2023 (with leap-year
%   day ranges included appropriately).

%Processing the precipitation data

lon = ncread('/Volumes/Expansion/Publication_Compound_Drought_Heatwave_Over_Phoneix/Datasets/CPC_Global_Precipitation/precip.1979.nc','lon');
lat = ncread('/Volumes/Expansion/Publication_Compound_Drought_Heatwave_Over_Phoneix/Datasets/CPC_Global_Precipitation/precip.1979.nc','lat');

folderPath = '/Volumes/Expansion/Publication_Compound_Drought_Heatwave_Over_Phoneix/Datasets/CPC_Global_Precipitation';
ncFiles = dir(fullfile(folderPath, '*.nc'));
ncFiles(1:46)= [];

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


SWNA             =  readmatrix('/Volumes/Expansion/Publication_Compound_Drought_Heatwave_Over_Phoneix/Datasets/Shapefiles/Southwest_north_America_Grid_Locations.csv');
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

SWNA             =  readtable('/Volumes/Expansion/Publication_Compound_Drought_Heatwave_Over_Phoneix/Datasets/Shapefiles/Southwest_north_America_Grid_Locations.csv');
%SWNA_latlon      =  SWNA(:,1:2);
%SWNA_index       =  SWNA(:,3);
States           =  SWNA(:,19);

index_AZ = find(strcmp(States.NAME, 'Arizona'));
index_UT = find(strcmp(States.NAME, 'Utah'));
index_NM = find(strcmp(States.NAME, 'New Mexico'));
index_CO = find(strcmp(States.NAME, 'Colorado'));
index_MX = find(strcmp(States.NAME, ''));

weeklyDataMatrix_r1                  = reshape(weeklyDataMatrix_r,[2256,52,45]); 

weeklyDataMatrix_ltm                 = mean(weeklyDataMatrix_r1(:,:,13:32),3);

weeklyDataMatrix_anm                 = zeros(2256,52,45);

for i=1:45
    weeklyDataMatrix_anm(:,:,i) = weeklyDataMatrix_r1(:,:,i) - weeklyDataMatrix_ltm;
end

weeklyDataMatrix_anm = weeklyDataMatrix_anm(:,14:39,:);
weeklyDataMatrix_anm_mx =mean(weeklyDataMatrix_anm(index_MX,:,:),2);
weeklyDataMatrix_anm_mx =reshape(weeklyDataMatrix_anm_mx,[700,1*45]);
weeklyDataMatrix_anm_mx = (mean(weeklyDataMatrix_anm_mx(:,2:45),1))';


weeklyDataMatrix_anm_az =mean(weeklyDataMatrix_anm(index_AZ,:,:),2);
weeklyDataMatrix_anm_az =reshape(weeklyDataMatrix_anm_az,[113,1*45]);
weeklyDataMatrix_anm_az = (mean(weeklyDataMatrix_anm_az(:,2:45),1))';


weeklyDataMatrix_anm_nm =mean(weeklyDataMatrix_anm(index_NM,:,:),2);
weeklyDataMatrix_anm_nm =reshape(weeklyDataMatrix_anm_nm,[122,1*45]);
weeklyDataMatrix_anm_nm = (mean(weeklyDataMatrix_anm_nm(:,2:45),1))';



weeklyDataMatrix_anm_mx = mean(reshape(weeklyDataMatrix_anm_mx,[11,4]),2);


weeklyDataMatrix_2023                = weeklyDataMatrix_r1(:,:,45);

weeklyDataMatrix_2023_anm            = weeklyDataMatrix_2023-weeklyDataMatrix_ltm;

weeklyDataMatrix_2023_anm_summer     = weeklyDataMatrix_2023_anm(:,14:39);
weeklyDataMatrix_2023_anm_summer_az  = (mean(weeklyDataMatrix_2023_anm_summer(index_AZ,:),1))';

weeklyDataMatrix_2023_anm_summer_mx  = (mean(weeklyDataMatrix_2023_anm_summer(index_MX,:),1))';
weeklyDataMatrix_2023_anm_summer_nm  = (mean(weeklyDataMatrix_2023_anm_summer(index_NM,:),1))';

sp_az  = [rescale(weeklyDataMatrix_2023_anm_summer_az),rescale(all_data_anmw_azm)];
sp_mx  = [rescale(weeklyDataMatrix_2023_anm_summer_mx),rescale(all_data_anmw_mxm)];
sp_nm  = [rescale(weeklyDataMatrix_2023_anm_summer_nm),rescale(all_data_anmw_nmm)];


%Extracting CDHW data statewise
load('Daytime_CDHW_severity_matrix_0702025.mat');

day_cdhw_severity           =       CDHW_hw_dr3_sev;
day_hw_severity             =       CDHW_hw_sev_rw;
day_cdhw_occurrence         =       Compound_dr_hw;
drought_index               =       SW_SPI3_summer;
drought_occurrence          =       SW_drought3_summer;

load('Nighttime_CDHW_severity_matrix_0702025.mat');

night_cdhw_severity         =       CDHW_hw_dr3_sev;
night_hw_severity           =       CDHW_hw_sev_rw;
night_cdhw_occurrence       =       Compound_dr_hw;




day_cdhw_severity_az_2023     =       (mean(day_cdhw_severity(index_AZ,1119:1144),1))';
day_cdhw_severity_mx_2023     =       (mean(day_cdhw_severity(index_MX,1119:1144),1))';
day_cdhw_severity_nm_2023     =       (mean(day_cdhw_severity(index_NM,1119:1144),1))';

night_cdhw_severity_az_2023   =       (mean(night_cdhw_severity(index_AZ,1119:1144),1))';
night_cdhw_severity_mx_2023   =       (mean(night_cdhw_severity(index_MX,1119:1144),1))';
night_cdhw_severity_nm_2023   =       (mean(night_cdhw_severity(index_NM,1119:1144),1))';

%extracting CDHW severity for all years for mexico, Arizona & New Mexico
day_cdhw_severity_mx          =       day_cdhw_severity(index_MX,:); 
day_cdhw_severity_mx_mean     =       (mean(day_cdhw_severity_mx,1))';
day_cdhw_severity_mx_meany    =       reshape(day_cdhw_severity_mx_mean,[44,26]);
day_cdhw_severity_mx_meany    =       mean(day_cdhw_severity_mx_meany,2);

night_cdhw_severity_mx        =       night_cdhw_severity(index_MX,:); 
night_cdhw_severity_mx_mean   =       (mean(night_cdhw_severity_mx,1))';
night_cdhw_severity_mx_meany    =       reshape(night_cdhw_severity_mx_mean,[44,26]);
night_cdhw_severity_mx_meany    =       mean(night_cdhw_severity_mx_meany,2);

day_cdhw_severity_az          =       day_cdhw_severity(index_AZ,:); 
day_cdhw_severity_az_mean     =       (mean(day_cdhw_severity_az,1))';
day_cdhw_severity_az_meany    =       reshape(day_cdhw_severity_az_mean,[44,26]);
day_cdhw_severity_az_meany    =       mean(day_cdhw_severity_az_meany,2);

night_cdhw_severity_az        =       night_cdhw_severity(index_AZ,:); 
night_cdhw_severity_az_mean   =       (mean(night_cdhw_severity_az,1))';
night_cdhw_severity_az_meany    =       reshape(night_cdhw_severity_az_mean,[44,26]);
night_cdhw_severity_az_meany    =       mean(night_cdhw_severity_az_meany,2);

day_cdhw_severity_nm          =       day_cdhw_severity(index_NM,:); 
day_cdhw_severity_nm_mean     =       (mean(day_cdhw_severity_nm,1))';
day_cdhw_severity_nm_meany    =       reshape(day_cdhw_severity_nm_mean,[44,26]);
day_cdhw_severity_nm_meany    =       mean(day_cdhw_severity_nm_meany,2);

night_cdhw_severity_nm        =       night_cdhw_severity(index_NM,:); 
night_cdhw_severity_nm_mean   =       (mean(night_cdhw_severity_nm,1))';
night_cdhw_severity_nm_meany    =       reshape(night_cdhw_severity_nm_mean,[44,26]);
night_cdhw_severity_nm_meany    =       mean(night_cdhw_severity_nm_meany,2);


%Counting the compound day and night
day_cdhw_occurrence_az        =       (sum(day_cdhw_occurrence(index_AZ,:),1)/1.13)';
day_cdhw_occurrence_mx        =       (sum(day_cdhw_occurrence(index_MX,:),1)/7)';
day_cdhw_occurrence_nm        =       (sum(day_cdhw_occurrence(index_NM,:),1)/1.22)';

day_cdhw_occurrence_az        =     reshape(day_cdhw_occurrence_az,[44,26]);
day_cdhw_occurrence_az        =     mean(day_cdhw_occurrence_az,2);
day_cdhw_occurrence_mx        =     reshape(day_cdhw_occurrence_mx,[44,26]);
day_cdhw_occurrence_mx        =     mean(day_cdhw_occurrence_mx,2);
day_cdhw_occurrence_nm        =     reshape(day_cdhw_occurrence_nm,[44,26]);
day_cdhw_occurrence_nm        =     mean(day_cdhw_occurrence_nm,2);


%Continuation before extracting CDHW...

sp_az  = [normalize(weeklyDataMatrix_2023_anm_summer_az),normalize(all_data_anm1_az_mean_2023), normalize(day_cdhw_severity_az_2023),normalize(night_cdhw_severity_az_2023)];

sp_nm  = [normalize(weeklyDataMatrix_2023_anm_summer_nm),normalize(all_data_anm1_nm_mean_2023),normalize(day_cdhw_severity_nm_2023),normalize(night_cdhw_severity_nm_2023)];


sp_mx  = [normalize(weeklyDataMatrix_2023_anm_summer_mx),normalize(all_data_anm1_mx_mean_2023),normalize(day_cdhw_severity_mx_2023),normalize(night_cdhw_severity_mx_2023)];


sp_az1  = [weeklyDataMatrix_2023_anm_summer_az,all_data_anm1_az_mean_2023, day_cdhw_severity_az_2023,night_cdhw_severity_az_2023];

sp_nm1  = [weeklyDataMatrix_2023_anm_summer_nm,all_data_anm1_nm_mean_2023,day_cdhw_severity_nm_2023,night_cdhw_severity_nm_2023];


sp_mx1  = [weeklyDataMatrix_2023_anm_summer_mx,all_data_anm1_mx_mean_2023,day_cdhw_severity_mx_2023,night_cdhw_severity_mx_2023];



%-----------------------------------------------------------------
%  Four‑series × Three‑region time‑series stack (Nature style)
%-----------------------------------------------------------------
%  INPUTS (26×4 each)
data_mx = sp_mx1;        % Mexico
data_az = sp_az1;        % Arizona
data_nm = sp_nm1;        % New Mexico

%-----------------------------------------------------------------
%  Figure & global aesthetics
%-----------------------------------------------------------------
f = figure('Units','centimeters','Position',[5 5 16 18], ...
           'Color','w','Renderer','painters');  % vector PDF‑friendly

set(f,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',8, ...
      'DefaultLineLineWidth',1.2,'DefaultAxesLineWidth',0.8, ...
      'DefaultAxesTickDir','out','DefaultAxesBox','off');

% ColourBrewer “Set2” for the four series (colour‑blind safe)
pal = [102 194 165;
       252 141  98;
       141 160 203;
       231 138 195] ./ 255;

regionNames = {'Mexico','Arizona','New Mexico'};

%-----------------------------------------------------------------
%  Tiled layout: 4 rows × 3 columns
%-----------------------------------------------------------------
tlo = tiledlayout(f,4,3,'TileSpacing','compact','Padding','compact');

for r = 1:4                     % row  = series index
    for c = 1:3                 % col  = region
        ax = nexttile;
        %--- select the appropriate column vector ----------------
        switch c
            case 1, y = data_mx(:,r);
            case 2, y = data_az(:,r);
            case 3, y = data_nm(:,r);
        end
        
        plot(1:26, y, ...
             'Color', pal(r,:), 'Marker','o', ...
             'MarkerFaceColor', pal(r,:), 'MarkerSize',6);

        %--- axes cosmetics --------------------------------------
        grid(ax,'on'); ax.GridAlpha = 0.15;
        
        % Row labels only on first column
        if c == 1
            ylabel(sprintf('Series %d',r), ...
                   'FontSize',8,'FontWeight','bold');
        end
        
        % Column titles only on top row
        if r == 1
            title(regionNames{c}, ...
                  'FontSize',16,'FontWeight','bold');
        end
        
        % Hide intermediate x‑tick labels
        if r < 4
            ax.XTickLabel = [];
        else
            xlabel('Time (weeks)','FontSize',16,'FontWeight','bold');
        end
    end
end

%-----------------------------------------------------------------
%  Figure‑level title (optional) & export
%-----------------------------------------------------------------
% Your data
% Your data

data = sp_az; % 26x4 matrix
varnames = {'Var1', 'Var2', 'Var3', 'Var4'}; % optional: your variable names
nvars = size(data,2);

% Create figure
figure('Units','inches','Position',[1 1 6 6],'Color','w');

% Get all pairs (i,j) where i<j
pairs = nchoosek(1:nvars,2); % all variable pairs
npairs = size(pairs,1);

% Setup a compact grid
nrows = ceil(sqrt(npairs));
ncols = ceil(npairs / nrows);

tiledlayout(nrows, ncols, 'Padding', 'compact', 'TileSpacing', 'compact');

for k = 1:npairs
    i = pairs(k,1);
    j = pairs(k,2);

    nexttile
    x = data(:,i);
    y = data(:,j);
    valid = ~isnan(x) & ~isnan(y);

    % Scatter plot
    scatter(x(valid), y(valid), 20, ...
            'filled', 'MarkerFaceAlpha', 0.7, ...
            'MarkerFaceColor', [0.5 0.5 0.5]);
    hold on

    % Regression line
    if sum(valid) >= 5
        p = polyfit(x(valid), y(valid), 1);
        xfit = linspace(min(x(valid)), max(x(valid)), 100);
        yfit = polyval(p, xfit);
        plot(xfit, yfit, '-', 'Color', [0 0.2 0.8], 'LineWidth', 3);

        % Correlation and p-value
        [R, Pval] = corr(x(valid), y(valid), 'rows', 'pairwise');
        text(0.5, 0.9, sprintf('r=%.2f\np=%.3f', R, Pval), ...
            'Units', 'normalized', ...
            'HorizontalAlignment', 'center', ...
            'FontSize', 16, ...
            'FontName', 'Arial', ...
            'Color', [0.1 0.1 0.1]);
    end
    hold off

    % Beautify
    set(gca, 'FontName', 'Arial', 'FontSize', 16, ...
             'XColor', [0.2 0.2 0.2], 'YColor', [0.2 0.2 0.2]);

    % Labels (optional, could remove if too crowded)
    xlabel(varnames{i}, 'FontSize', 20);
    ylabel(varnames{j}, 'FontSize', 20);
end

% Export to vector graphic
exportgraphics(gcf, 'scattermatrix_compact_final.pdf', 'ContentType', 'vector');














%Analysis of compounding drought over day

maxLagDays  = 3;                                 % search –3 … +3
dd_cdhw_coherence     = zeros(2256,2256);
dd_cdhw_lag           = zeros(2256,2256);

for i = 1:2256
    for j = 1:2256
        %–– choose the two time series ––%
        series1 = (day_cdhw_severity(i ,:))';       % daytime @ grid‑i
        series2 = (day_cdhw_severity(j,:))';      % ← use j if cross‑grid
        % series2 = night_hw_severity(i,:).';    % ← use i if same grid

        % cross‑correlation and associated lags
        [c,lags] = xcorr(series1,series2,maxLagDays,'coeff');

        % strongest positive / negative correlation and their lags
        [posCorr,posIdx] = max(c);               % ≥ 0
        [negCorr,negIdx] = min(c);               % ≤ 0
        posLag = lags(posIdx);
        negLag = lags(negIdx);

        % keep whichever has the larger |r|
        if abs(posCorr) >= abs(negCorr)
            dd_cdhw_coherence(i,j) =  posCorr;
            dd_cdhw_lag(i,j)       =  posLag;
        else
            dd_cdhw_coherence(i,j) =  negCorr;
            dd_cdhw_lag(i,j)       =  negLag;
        end
    end
    i
end

dd_cdhw_coherence_az = (mean(dd_cdhw_coherence(index_AZ,:),1))';
dd_cdhw_lag_az = (mean(dd_cdhw_lag(index_AZ,:),1))';

%Analysis of compounding drought over Night

maxLagDays  = 3;                                 % search –3 … +3
nn_cdhw_coherence     = zeros(2256,2256);
nn_cdhw_lag           = zeros(2256,2256);

for i = 1:2256
    for j = 1:2256
        %–– choose the two time series ––%
        series1 = (night_cdhw_severity(i ,:))';       % daytime @ grid‑i
        series2 = (night_cdhw_severity(j,:))';      % ← use j if cross‑grid
        % series2 = night_hw_severity(i,:).';    % ← use i if same grid

        % cross‑correlation and associated lags
        [c,lags] = xcorr(series1,series2,maxLagDays,'coeff');

        % strongest positive / negative correlation and their lags
        [posCorr,posIdx] = max(c);               % ≥ 0
        [negCorr,negIdx] = min(c);               % ≤ 0
        posLag = lags(posIdx);
        negLag = lags(negIdx);

        % keep whichever has the larger |r|
        if abs(posCorr) >= abs(negCorr)
            nn_cdhw_coherence(i,j) =  posCorr;
            nn_cdhw_lag(i,j)       =  posLag;
        else
            nn_cdhw_coherence(i,j) =  negCorr;
            nn_cdhw_lag(i,j)       =  negLag;
        end
    end
    i
end

nn_cdhw_coherence_az = (mean(`nn_cdhw_coherence(index_AZ,:),1))';
nn_cdhw_lag_az = (mean(nn_cdhw_lag(index_AZ,:),1))';





%Analysis of compounding of drought over day and night

maxLagDays  = 3;                                 % search –3 … +3
dn_cdhw_coherence     = zeros(2256,2256);
dn_cdhw_lag           = zeros(2256,2256);

for i = 1:2256
    for j = 1:2256
        %–– choose the two time series ––%
        series1 = (day_cdhw_severity(i ,:))';       % daytime @ grid‑i
        series2 = (night_cdhw_severity(j,:))';      % ← use j if cross‑grid
        % series2 = night_hw_severity(i,:).';    % ← use i if same grid

        % cross‑correlation and associated lags
        [c,lags] = xcorr(series1,series2,maxLagDays,'coeff');

        % strongest positive / negative correlation and their lags
        [posCorr,posIdx] = max(c);               % ≥ 0
        [negCorr,negIdx] = min(c);               % ≤ 0
        posLag = lags(posIdx);
        negLag = lags(negIdx);

        % keep whichever has the larger |r|
        if abs(posCorr) >= abs(negCorr)
            dn_cdhw_coherence(i,j) =  posCorr;
            dn_cdhw_lag(i,j)       =  posLag;
        else
            dn_cdhw_coherence(i,j) =  negCorr;
            dn_cdhw_lag(i,j)       =  negLag;
        end
    end
    i
end

dn_cdhw_coherence_az = (mean(dn_cdhw_coherence(index_AZ,:),1))';
dn_cdhw_lag_az = (mean(dn_cdhw_lag(index_AZ,:),1))';


%Calculating soil moisture temperature coupling

%day_cdhw_severity_mx_mean
%night_cdhw_severity_mx_mean

st_corr = zeros(44,1);

maxLagDays = 3;

for i=1:44
    j                =   26*(i-1)+1;
    k                =   26*i;


    series1          =  night_cdhw_severity_mx_mean(j:k,1);

    series2          =  all_data_anm1_mx_mean(j:k,1);

    [crossCorr, lags] = xcorr(series1, series2, maxLagDays, 'coeff');
        
        % Find the lag with the highest absolute correlation
     [maxCorr, maxIndex] = max(crossCorr);
     st_corr(i,1)    =  maxCorr;
end



st_corr    = zeros(44,1);
maxLagDays = 3;                       % ±3‑day window (total 7 coefficients)

for i = 1:44
    j       = 26*(i-1)+1;
    k       = 26*i;

    series1 = night_cdhw_severity_mx_mean(j:k,1);
    series2 = all_data_anm1_mx_mean   (j:k,1);

    [crossCorr, ~] = xcorr(series1, series2, maxLagDays, 'coeff');

    % strongest positive and strongest negative
    posCorr = max(crossCorr);         % ≥ 0
    negCorr = min(crossCorr);         % ≤ 0

    % keep the one with the larger absolute magnitude
    if abs(posCorr) >= abs(negCorr)
        st_corr(i) = posCorr;         % dominant positive correlation
    else
        st_corr(i) = negCorr;         % dominant negative correlation
    end
end

scorrelation = [normalize(st_corr),normalize(all_data_anm1_mx_mean_y)];


st_corr_az    = zeros(44,1);
maxLagDays = 3;                       % ±3‑day window (total 7 coefficients)

for i = 1:44
    j       = 26*(i-1)+1;
    k       = 26*i;

    series1 = night_cdhw_severity_az_mean(j:k,1);
    series2 = all_data_anm1_az_mean   (j:k,1);

    [crossCorr, ~] = xcorr(series1, series2, maxLagDays, 'coeff');

    % strongest positive and strongest negative
    posCorr = max(crossCorr);         % ≥ 0
    negCorr = min(crossCorr);         % ≤ 0

    % keep the one with the larger absolute magnitude
    if abs(posCorr) >= abs(negCorr)
        st_corr_az(i) = posCorr;         % dominant positive correlation
    else
        st_corr_az(i) = negCorr;         % dominant negative correlation
    end
end

scorrelation = [normalize(st_corr_az),normalize(all_data_anm1_mx_mean_y)];

st_corr_nm    = zeros(44,1);
maxLagDays = 3;                       % ±3‑day window (total 7 coefficients)

for i = 1:44
    j       = 26*(i-1)+1;
    k       = 26*i;

    series1 = night_cdhw_severity_nm_mean(j:k,1);
    series2 = all_data_anm1_nm_mean   (j:k,1);

    [crossCorr, ~] = xcorr(series1, series2, maxLagDays, 'coeff');

    % strongest positive and strongest negative
    posCorr = max(crossCorr);         % ≥ 0
    negCorr = min(crossCorr);         % ≤ 0

    % keep the one with the larger absolute magnitude
    if abs(posCorr) >= abs(negCorr)
        st_corr_nm(i) = posCorr;         % dominant positive correlation
    else
        st_corr_nm(i) = negCorr;         % dominant negative correlation
    end
end

a = [st_corr, st_corr_az, st_corr_nm];



st_corr_az_mx    = zeros(44,1);
maxLagDays = 3;                       % ±3‑day window (total 7 coefficients)

for i = 1:44
    j       = 26*(i-1)+1;
    k       = 26*i;

    series1 = night_cdhw_severity_mx_mean(j:k,1);
    series2 = all_data_anm1_az_mean   (j:k,1);

    [crossCorr, ~] = xcorr(series1, series2, maxLagDays, 'coeff');

    % strongest positive and strongest negative
    posCorr = max(crossCorr);         % ≥ 0
    negCorr = min(crossCorr);         % ≤ 0

    % keep the one with the larger absolute magnitude
    if abs(posCorr) >= abs(negCorr)
        st_corr_az_mx(i) = posCorr;         % dominant positive correlation
    else
        st_corr_az_mx(i) = negCorr;         % dominant negative correlation
    end
end

st_corr_mx_az    = zeros(44,1);
maxLagDays = 3;                       % ±3‑day window (total 7 coefficients)

for i = 1:44
    j       = 26*(i-1)+1;
    k       = 26*i;

    series1 = night_cdhw_severity_az_mean(j:k,1);
    series2 = all_data_anm1_mx_mean   (j:k,1);

    [crossCorr, ~] = xcorr(series1, series2, maxLagDays, 'coeff');

    % strongest positive and strongest negative
    posCorr = max(crossCorr);         % ≥ 0
    negCorr = min(crossCorr);         % ≤ 0

    % keep the one with the larger absolute magnitude
    if abs(posCorr) >= abs(negCorr)
        st_corr_mx_az(i) = posCorr;         % dominant positive correlation
    else
        st_corr_mx_az(i) = negCorr;         % dominant negative correlation
    end
end




st_corr_nm_mx    = zeros(44,1);
maxLagDays = 3;                       % ±3‑day window (total 7 coefficients)

for i = 1:44
    j       = 26*(i-1)+1;
    k       = 26*i;

    series1 = night_cdhw_severity_mx_mean(j:k,1);
    series2 = all_data_anm1_nm_mean   (j:k,1);

    [crossCorr, ~] = xcorr(series1, series2, maxLagDays, 'coeff');

    % strongest positive and strongest negative
    posCorr = max(crossCorr);         % ≥ 0
    negCorr = min(crossCorr);         % ≤ 0

    % keep the one with the larger absolute magnitude
    if abs(posCorr) >= abs(negCorr)
        st_corr_nm_mx(i) = posCorr;         % dominant positive correlation
    else
        st_corr_nm_mx(i) = negCorr;         % dominant negative correlation
    end
end


st_corr_mx_nm    = zeros(44,1);
maxLagDays = 3;                       % ±3‑day window (total 7 coefficients)

for i = 1:44
    j       = 26*(i-1)+1;
    k       = 26*i;

    series1 = night_cdhw_severity_nm_mean(j:k,1);
    series2 = all_data_anm1_mx_mean   (j:k,1);

    [crossCorr, ~] = xcorr(series1, series2, maxLagDays, 'coeff');

    % strongest positive and strongest negative
    posCorr = max(crossCorr);         % ≥ 0
    negCorr = min(crossCorr);         % ≤ 0

    % keep the one with the larger absolute magnitude
    if abs(posCorr) >= abs(negCorr)
        st_corr_mx_nm(i) = posCorr;         % dominant positive correlation
    else
        st_corr_mx_nm(i) = negCorr;         % dominant negative correlation
    end
end

b = [st_corr, st_corr_mx_az, st_corr_mx_nm];

c = [st_corr_az, st_corr_mx_az];


% Your data matrix a = [st_corr, st_corr_az, st_corr_nm];
pairs = [1 2; 1 3; 3 2];
titles = {'st\_corr vs st\_corr\_az', 'st\_corr vs st\_corr\_nm', 'st\_corr\_nm vs st\_corr\_az'};

figure('Position', [100, 100, 1200, 350]);

for i = 1:3
    idx = pairs(i,:);
    X = a(:, idx);  % Extract pair

    maxK = 6;
    sil_scores = NaN(maxK, 1);  % Preallocate
    all_labels = cell(maxK,1);

    % Loop over k from 2 to maxK
    for k = 2:maxK
        try
            [labels, ~] = kmeans(X, k, 'Replicates', 20, 'Start', 'plus', 'Display', 'off');
            s = silhouette(X, labels);
            sil_scores(k) = mean(s);
            all_labels{k} = labels;
        catch
            % If kmeans or silhouette fails, skip
            sil_scores(k) = NaN;
        end
    end

    % Get the best k (ignoring NaN)
    [max_sil, opt_k] = max(sil_scores(2:end));
    opt_k = opt_k + 1;  % since indexing starts at 2

    % Fallback if silhouette fails
    if isnan(max_sil) || max_sil < 0.1
        warning('Low silhouette or failure for plot %d — fallback to k=3.', i);
        opt_k = 3;
        [labels, C] = kmeans(X, opt_k, 'Replicates', 20, 'Start', 'plus', 'Display', 'off');
    else
        labels = all_labels{opt_k};
        [~, C] = kmeans(X, opt_k, 'Replicates', 20, 'Start', 'plus', 'Display', 'off');
    end

    % Plotting
    subplot(1, 3, i);
    gscatter(X(:,1), X(:,2), labels, lines(opt_k), 'o', 8, 'on');
    hold on;
    plot(C(:,1), C(:,2), 'kx', 'MarkerSize', 15, 'LineWidth', 2);
    hold off;
    xlabel(sprintf('a(:,%d)', idx(1)));
    ylabel(sprintf('a(:,%d)', idx(2)));
    title([titles{i} sprintf(' | Optimal k = %d', opt_k)]);
    legendStrings = arrayfun(@(c) sprintf('Cluster %d', c), 1:opt_k, 'UniformOutput', false);
    legend([legendStrings, {'Centroids'}], 'Location', 'best');
    grid on;
end

sgtitle('K-means Clustering with Robust Silhouette-based Optimal k');

title(sprintf('%s | Optimal k = %d | Silhouette = %.2f', titles{i}, opt_k, max_sil));




% Your matrix: a = [st_corr, st_corr_az, st_corr_nm];
pairs = [1 2; 1 3; 3 2];
titles = {'st\_corr vs st\_corr\_az', 'st\_corr vs st\_corr\_nm', 'st\_corr\_nm vs st\_corr\_az'};

figure('Position', [100, 100, 1200, 350]);

for i = 1:3
    idx = pairs(i,:);
    X = a(:, idx);  % Extract variable pair

    maxK = 6;
    sil_scores = NaN(maxK, 1);
    all_labels = cell(maxK, 1);

    % Try clustering from k=2 to k=6
    for k = 2:maxK
        try
            [labels, ~] = kmeans(X, k, 'Replicates', 20, 'Start', 'plus', 'Display', 'off');
            s = silhouette(X, labels);
            sil_scores(k) = mean(s);
            all_labels{k} = labels;
        catch
            sil_scores(k) = NaN;
        end
    end

    % Determine best k
    [max_sil, opt_k] = max(sil_scores(2:end));
    opt_k = opt_k + 1;

    if isnan(max_sil) || max_sil < 0.1
        warning('Low silhouette or failure for plot %d — fallback to k=3.', i);
        opt_k = 3;
        [labels, C] = kmeans(X, opt_k, 'Replicates', 20, 'Start', 'plus', 'Display', 'off');
        max_sil = NaN;
    else
        labels = all_labels{opt_k};
        [~, C] = kmeans(X, opt_k, 'Replicates', 20, 'Start', 'plus', 'Display', 'off');
    end

    % Plotting
    subplot(1, 3, i);
    gscatter(X(:,1), X(:,2), labels, lines(opt_k), 'o', 8, 'on');
    hold on;
    plot(C(:,1), C(:,2), 'kx', 'MarkerSize', 15, 'LineWidth', 2);
    hold off;
    xlabel(sprintf('a(:,%d)', idx(1)));
    ylabel(sprintf('a(:,%d)', idx(2)));
    if isnan(max_sil)
        sil_text = 'N/A';
    else
        sil_text = sprintf('%.2f', max_sil);
    end
    title(sprintf('%s | k = %d | Silhouette = %s', titles{i}, opt_k, sil_text));
    legendStrings = arrayfun(@(c) sprintf('Cluster %d', c), 1:opt_k, 'UniformOutput', false);
    legend([legendStrings, {'Centroids'}], 'Location', 'best');
    grid on;
end

sgtitle('K-means Clustering with Silhouette-based Optimal k and Scores');




% Your data
% st = ... (44x7 matrix)

% Your data
% st = ... (44x7 matrix)

% Your data
% st = ... (44x7 matrix)

%For Arizona & Mexico

st = [st_corr,st_corr_az,st_corr_az_mx,st_corr_mx_az];
p  = [normalize(weeklyDataMatrix_anm_mx),normalize(weeklyDataMatrix_anm_az)];

%all_data_anm1_mx_meany = mean(reshape(all_data_anm1_mx_mean,[44,26]),2);
%all_data_anm1_az_meany = mean(reshape(all_data_anm1_az_mean,[44,26]),2);

all_data_anm1_mx_meany  = reshape(all_data_anm1_mx_mean,[44,26]);
all_data_anm1_mx_meany  = mean(all_data_anm1_mx_meany(44,1:8),2);

all_data_anm1_az_meany  = reshape(all_data_anm1_az_mean,[44,26]);
all_data_anm1_az_meany  = mean(all_data_anm1_az_meany(44,1:8),2);
sm  = [normalize(all_data_anm1_mx_meany),normalize(all_data_anm1_az_meany)];

figure('Units','inches','Position',[1 1 6 5],'Color','w');
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

% Define colors
precip_color = [0.2 0.6 0.8]; % Blueish for Precipitation
soilmoist_color = [0.8 0.4 0.2]; % Orangish for Soil Moisture
ts_color1 = [0.1 0.4 0.8]; % Color for st line 1
ts_color2 = [0.8 0.2 0.2]; % Color for st line 2
ts_color3 = [0.2 0.6 0.3]; % Color for st line 3
ts_color4 = [0.5 0.2 0.7]; % Color for st line 4

bar_width = 0.4; % Width of bars

% X-axis values
x1 = 1:size(p,1);
x2 = 1:size(p,1); % assuming same size for simplicity

% --- First Panel ---
nexttile
hold on

yyaxis left
b1 = bar(x1, p(:,1), bar_width, 'FaceColor', precip_color, 'EdgeColor', 'none','FaceAlpha',0.5);
b2 = bar(x1, sm(:,1), bar_width, 'FaceColor', soilmoist_color, 'EdgeColor', 'none', 'FaceAlpha',0.5); % Set transparency for overlapping bars
ylabel('Precipitation / Soil Moisture', 'FontSize', 10, 'FontName', 'Arial');

yyaxis right
plot(x1, st(:,1), '-o', 'LineWidth', 1.8, 'MarkerSize', 4, 'Color', ts_color1);
plot(x1, st(:,3), '-s', 'LineWidth', 1.8, 'MarkerSize', 4, 'Color', ts_color2);
ylabel('Influence Metrics', 'FontSize', 10, 'FontName', 'Arial');

legend({'Precipitation (1)','Soil Moisture (1)','MX\_soil-MX\_night','MX\_night-AZ\_soil'}, ...
    'Location', 'best', 'FontSize', 8);

title('First Panel: Soil Moisture Influence & Hydroclimate Variables', 'FontWeight', 'normal', 'FontSize', 11);
set(gca, 'FontSize', 9, 'FontName', 'Arial');
grid on

% --- Second Panel ---
nexttile
hold on

yyaxis left
b3 = bar(x2, p(:,2), bar_width, 'FaceColor', precip_color, 'EdgeColor', 'none','FaceAlpha',0.5);
b4 = bar(x2, sm(:,2), bar_width, 'FaceColor', soilmoist_color, 'EdgeColor', 'none', 'FaceAlpha',0.5); % Transparent
ylabel('Precipitation / Soil Moisture', 'FontSize', 10, 'FontName', 'Arial');

yyaxis right
plot(x2, st(:,2), '-^', 'LineWidth', 1.8, 'MarkerSize', 4, 'Color', ts_color3);
plot(x2, st(:,4), '-d', 'LineWidth', 1.8, 'MarkerSize', 4, 'Color', ts_color4);
ylabel('Influence Metrics', 'FontSize', 10, 'FontName', 'Arial');

legend({'Precipitation (2)','Soil Moisture (2)','AZ\_soil-AZ\_night','MX\_soil-AZ\_night'}, ...
    'Location', 'best', 'FontSize', 8);

title('Second Panel: Arizona CDHW Influence & Hydroclimate Variables', 'FontWeight', 'normal', 'FontSize', 11);
set(gca, 'FontSize', 9, 'FontName', 'Arial');
grid on

% Common x-label
xlabel('Index', 'FontSize', 10, 'FontName', 'Arial');

% Export to vector graphic
exportgraphics(gcf, 'combined_barplot_timeseries_FIXED.pdf', 'ContentType', 'vector');

% For Mexico & New Mexico

st_nm = [st_corr,st_corr_nm,st_corr_nm_mx,st_corr_mx_nm];
p_nm  = [normalize(weeklyDataMatrix_anm_mx),normalize(weeklyDataMatrix_anm_nm)];

%all_data_anm1_mx_meany = mean(reshape(all_data_anm1_mx_mean,[44,26]),2);
%all_data_anm1_az_meany = mean(reshape(all_data_anm1_az_mean,[44,26]),2);

all_data_anm1_mx_meany  = reshape(all_data_anm1_mx_mean,[44,26]);
all_data_anm1_mx_meany  = mean(all_data_anm1_mx_meany,2);

all_data_anm1_nm_meany  = reshape(all_data_anm1_nm_mean,[44,26]);
all_data_anm1_nm_meany  = mean(all_data_anm1_nm_meany,2);
sm_nm  = [normalize(all_data_anm1_mx_meany),normalize(all_data_anm1_nm_meany)];

figure('Units','inches','Position',[1 1 6 5],'Color','w');
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

% Define colors
precip_color = [0.2 0.6 0.8]; % Blueish for Precipitation
soilmoist_color = [0.8 0.4 0.2]; % Orangish for Soil Moisture
ts_color1 = [0.1 0.4 0.8]; % Color for st line 1
ts_color2 = [0.8 0.2 0.2]; % Color for st line 2
ts_color3 = [0.2 0.6 0.3]; % Color for st line 3
ts_color4 = [0.5 0.2 0.7]; % Color for st line 4

bar_width = 0.4; % Width of bars

% X-axis values
x1 = 1:size(p_nm,1);
x2 = 1:size(p_nm,1); % assuming same size for simplicity

% --- First Panel ---
nexttile
hold on

yyaxis left
b1 = bar(x1, p_nm(:,1), bar_width, 'FaceColor', precip_color, 'EdgeColor', 'none','FaceAlpha',0.5);
b2 = bar(x1, sm_nm(:,1), bar_width, 'FaceColor', soilmoist_color, 'EdgeColor', 'none', 'FaceAlpha',0.5); % Set transparency for overlapping bars
ylabel('Precipitation / Soil Moisture', 'FontSize', 10, 'FontName', 'Arial');

yyaxis right
plot(x1, st_nm(:,1), '-o', 'LineWidth', 1.8, 'MarkerSize', 4, 'Color', ts_color1);
plot(x1, st_nm(:,3), '-s', 'LineWidth', 1.8, 'MarkerSize', 4, 'Color', ts_color2);
ylabel('Influence Metrics', 'FontSize', 10, 'FontName', 'Arial');

legend({'Precipitation (1)','Soil Moisture (1)','MX\_soil-MX\_night','MX\_night-AZ\_soil'}, ...
    'Location', 'best', 'FontSize', 8);

title('First Panel: Soil Moisture Influence & Hydroclimate Variables', 'FontWeight', 'normal', 'FontSize', 11);
set(gca, 'FontSize', 9, 'FontName', 'Arial');
grid on

% --- Second Panel ---
nexttile
hold on

yyaxis left
b3 = bar(x2, p_nm(:,2), bar_width, 'FaceColor', precip_color, 'EdgeColor', 'none','FaceAlpha',0.5);
b4 = bar(x2, sm_nm(:,2), bar_width, 'FaceColor', soilmoist_color, 'EdgeColor', 'none', 'FaceAlpha',0.5); % Transparent
ylabel('Precipitation / Soil Moisture', 'FontSize', 10, 'FontName', 'Arial');

yyaxis right
plot(x2, st_nm(:,2), '-^', 'LineWidth', 1.8, 'MarkerSize', 4, 'Color', ts_color3);
plot(x2, st_nm(:,4), '-d', 'LineWidth', 1.8, 'MarkerSize', 4, 'Color', ts_color4);
ylabel('Influence Metrics', 'FontSize', 10, 'FontName', 'Arial');

legend({'Precipitation (2)','Soil Moisture (2)','AZ\_soil-AZ\_night','MX\_soil-AZ\_night'}, ...
    'Location', 'best', 'FontSize', 8);

title('Second Panel: Arizona CDHW Influence & Hydroclimate Variables', 'FontWeight', 'normal', 'FontSize', 11);
set(gca, 'FontSize', 9, 'FontName', 'Arial');
grid on

% Common x-label
xlabel('Index', 'FontSize', 10, 'FontName', 'Arial');

% Export to vector graphic
exportgraphics(gcf, 'combined_barplot_timeseries_FIXED.pdf', 'ContentType', 'vector');


% figure('Units','inches','Position',[1 1 6 5],'Color','w');
% tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
% 
% % First panel: st(:,1) and st(:,4)
% nexttile
% hold on
% plot(st(:,1), '-o', 'LineWidth', 1.8, 'MarkerSize', 4, 'Color', [0.1 0.4 0.8]);
% plot(st(:,3), '-s', 'LineWidth', 1.8, 'MarkerSize', 4, 'Color', [0.8 0.2 0.2]);
% hold off
% ylabel('Values', 'FontSize', 10, 'FontName', 'Arial');
% legend({'MX_soil-MXnight','MX_soil-AZ_night'}, 'Location', 'best', 'FontSize', 8);
% title('Influence of Soil Moisture over Mexico', 'FontWeight', 'normal', 'FontSize', 11);
% set(gca, 'FontSize', 9, 'FontName', 'Arial');
% grid on
% 
% % Second panel: st(:,2) and st(:,3)
% nexttile
% hold on
% plot(st(:,2), '-^', 'LineWidth', 1.8, 'MarkerSize', 4, 'Color', [0.2 0.6 0.3]);
% plot(st(:,4), '-d', 'LineWidth', 1.8, 'MarkerSize', 4, 'Color', [0.5 0.2 0.7]);
% hold off
% ylabel('Values', 'FontSize', 10, 'FontName', 'Arial');
% legend({'AZ_soil-AZ_night','MX_soil-AZ_night'}, 'Location', 'best', 'FontSize', 8);
% title('Influence on Arizona CDHW_n_i_g_h_t', 'FontWeight', 'normal', 'FontSize', 11);
% set(gca, 'FontSize', 9, 'FontName', 'Arial');
% grid on
% 
% % Common x-label
% xlabel('Index', 'FontSize', 10, 'FontName', 'Arial');
% 
% % Export to vector graphic
% exportgraphics(gcf, 'stackedplot_2x1_customized.pdf', 'ContentType', 'vector');



nindmx = find(st_corr<0);
nindaz = find(st_corr_az>0);

nint   = intersect(nindmx,nindaz);

st_corr_d = st_corr_az-st_corr;
st_corr_d_int = st_corr_d(nint,:);
st_corr_d_indmx = st_corr_d(nindmx,:);
highestCorrelation_AZMX_int = highestCorrelation_AZMX(nint,:);
highestCorrelation_AZMX_indmx = highestCorrelation_AZMX(nindmx,:);
smpt  = [st_corr_d, highestCorrelation_AZMX];
smp   = [st_corr_d_int,highestCorrelation_AZMX_int];
smpdx = [st_corr_d_indmx,highestCorrelation_AZMX_indmx];


[r1,p1] = corrcoef(smp);
r1      = r1(2,1);
p1      = p1(2,1);

[r2,p2] = corrcoef(smpdx);
r2      = r2(2,1);
p2      = p2(2,1);

[r3,p3] = corrcoef(smpt);
r3      = r3(2,1);
p3      = p3(2,1);


% Input data
datasets = {smpdx, smp, smpt};
titles = {'(a) smpdx', '(b) smp', '(c) smpt'};

% Provided correlation and p-values
r_vals = [r2, r1, r3];
p_vals = [p2, p1, p3];

% Create figure
figure('Units', 'normalized', 'OuterPosition', [0 0 0.5 1]);
%tiledlayout(3,1, 'TileSpacing', 'compact', 'Padding', 'compact');
tiledlayout(3,1, 'TileSpacing', 'compact');
for i = 1:3
    nexttile;

    data = datasets{i};
    x = data(:,1);
    y = data(:,2);

    % Linear regression
    p = polyfit(x, y, 1);
    slope = p(1);
    intercept = p(2);
    y_fit = polyval(p, x);

    % Scatter plot
    scatter(x, y, 20, 'MarkerFaceColor', [0.2 0.4 0.7], ...
            'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.7);
    hold on;

    % Sorted line plot
    [sorted_x, sort_idx] = sort(x);
    sorted_y = polyval(p, sorted_x);
    plot(sorted_x, sorted_y, 'r-', 'LineWidth', 2);

    % Labels and title
    xlabel('x', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('y', 'FontSize', 16, 'FontWeight', 'bold');
    title(titles{i}, 'FontSize', 16, 'FontWeight', 'bold');

    % Equation and correlation info
    regression_eqn = sprintf('y = %.2fx + %.2f', slope, intercept);
    correlation_info = sprintf('r = %.2f, p = %.3g', r_vals(i), p_vals(i));
    full_text = {regression_eqn; correlation_info};

    % Position text
    x_pos = min(x) + 0.05 * range(x);
    y_pos = max(y) - 0.1 * range(y);
    text(x_pos, y_pos, full_text, ...
        'FontSize', 14, 'BackgroundColor', 'none', 'EdgeColor', 'none', ...
        'Margin', 6, 'FontWeight', 'bold');

    % Styling
    grid on;
    set(gca, 'FontSize', 14, 'LineWidth', 1.2, 'Box', 'off');
    axis tight;
end

% Save high-quality figure (optional)
% exportgraphics(gcf, 'fig_linear_regression_nature.pdf', 'ContentType', 'vector');


% Data: smpdx, smp, smpt should be (n x 2) matrices

datasets = {smp,smpdx, smpt};
titles = {'(a) Linear Regression: smpdx', '(b) Linear Regression: smp', '(c) Linear Regression: smpt'};

% Bootstrapping parameters
nBoot = 5000; % number of bootstrap resamples

% Create figure
figure('Units', 'normalized', 'OuterPosition', [0 0 0.5 1]);
tiledlayout(3,1, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:3
    nexttile;
    data = datasets{i};
    x = data(:,1);
    y = data(:,2);

    % Linear regression (simple)
    p = polyfit(x, y, 1);
    slope = p(1);
    intercept = p(2);
    y_fit = polyval(p, x);

    % Calculate Pearson correlation
    [r_obs, ~] = corr(x, y, 'Type', 'Pearson');

    % Bootstrap correlation
    r_boot = zeros(nBoot,1);
    n = length(x);
    for b = 1:nBoot
        idx = randi(n, n, 1); % random resampling with replacement
        xb = x(idx);
        yb = y(idx);
        r_boot(b) = corr(xb, yb, 'Type', 'Pearson');
    end

    % Bootstrapped p-value (two-sided)
    p_corr = mean(abs(r_boot) >= abs(r_obs)) * 2;

    % Scatter plot
    scatter(x, y, 20, 'MarkerFaceColor', [0.2 0.4 0.7], ...
            'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.7);
    hold on;
    [sorted_x, sort_idx] = sort(x);
    sorted_y = polyval(p, sorted_x);
    plot(sorted_x, sorted_y, 'r-', 'LineWidth', 2);

    % Labels and title
    xlabel('x', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('y', 'FontSize', 16, 'FontWeight', 'bold');
    title(titles{i}, 'FontSize', 16, 'FontWeight', 'bold', 'Interpreter', 'none');

    % Equation and correlation info
    regression_eqn = sprintf('y = %.2fx + %.2f', slope, intercept);
    correlation_info = sprintf('r = %.2f, p = %.3g', r_obs, p_corr);
    full_text = {regression_eqn; correlation_info};

    % Place text
    x_pos = min(x) + 0.05 * range(x);
    y_pos = max(y) - 0.1 * range(y);
    text(x_pos, y_pos, full_text, ...
        'FontSize', 14, 'BackgroundColor', 'none', 'Margin', 6);

    % Styling
    grid on;
    set(gca, 'FontSize', 14, 'LineWidth', 1.2, 'Box', 'off');
    axis tight;
end

% Save high-quality figure (optional)
% exportgraphics(gcf, 'linear_regression_bootstrap_corr.pdf', 'ContentType', 'vector');

nindmx1 = find(st_corr<0);
nindaz1 = find(st_corr_nm>0);

nint1   = intersect(nindmx1,nindaz1);

st_corr_d1 = st_corr_nm-st_corr;
st_corr_d_int1 = st_corr_d1(nint1,:);
st_corr_d_indmx1 = st_corr_d1(nindmx1,:);
highest_correlation_NM_MX_int1 = highest_correlation_NM_MX(nint1,:);
highest_correlation_NM_MX_indmx1 = highest_correlation_NM_MX(nindmx1,:);
smpt  = [st_corr_d1, highest_correlation_NM_MX];
smp   = [st_corr_d_int1,highest_correlation_NM_MX_int1];
smpdx = [normalize(st_corr_d_indmx1),normalize(highest_correlation_NM_MX_indmx1)];


[r1,p1] = corrcoef(smp);
r1      = r1(2,1);
p1      = p1(2,1);

[r2,p2] = corrcoef(smpdx);
r2      = r2(2,1);
p2      = p2(2,1);

[r3,p3] = corrcoef(smpt);
r3      = r3(2,1);
p3      = p3(2,1);


% Input data
datasets = {smpdx, smp, smpt};
titles = {'(a) smpdx', '(b) smp', '(c) smpt'};

% Provided correlation and p-values
r_vals = [r2, r1, r3];
p_vals = [p2, p1, p3];

% Create figure
figure('Units', 'normalized', 'OuterPosition', [0 0 0.5 1]);
tiledlayout(3,1, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:3
    nexttile;

    data = datasets{i};
    x = data(:,1);
    y = data(:,2);

    % Linear regression
    p = polyfit(x, y, 1);
    slope = p(1);
    intercept = p(2);
    y_fit = polyval(p, x);

    % Scatter plot
    scatter(x, y, 20, 'MarkerFaceColor', [0.2 0.4 0.7], ...
            'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.7);
    hold on;

    % Sorted line plot
    [sorted_x, sort_idx] = sort(x);
    sorted_y = polyval(p, sorted_x);
    plot(sorted_x, sorted_y, 'r-', 'LineWidth', 2);

    % Labels and title
    xlabel('x', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('y', 'FontSize', 16, 'FontWeight', 'bold');
    title(titles{i}, 'FontSize', 16, 'FontWeight', 'bold');

    % Equation and correlation info
    regression_eqn = sprintf('y = %.2fx + %.2f', slope, intercept);
    correlation_info = sprintf('r = %.2f, p = %.3g', r_vals(i), p_vals(i));
    full_text = {regression_eqn; correlation_info};

    % Position text
    x_pos = min(x) + 0.05 * range(x);
    y_pos = max(y) - 0.1 * range(y);
    text(x_pos, y_pos, full_text, ...
        'FontSize', 14, 'BackgroundColor', 'none', 'EdgeColor', 'none', ...
        'Margin', 6, 'FontWeight', 'bold');

    % Styling
    grid on;
    set(gca, 'FontSize', 14, 'LineWidth', 1.2, 'Box', 'off');
    axis tight;
end


% Your data
x = s(:,1);
y = s(:,2);

% Observed correlation
r_obs = corr(x, y, 'Type', 'Pearson');

% Bootstrap settings
nboot = 5000;  % number of bootstrap samples
r_boot = zeros(nboot,1);

% Bootstrap loop
n = length(x);
for i = 1:nboot
    idx = randi(n, n, 1);   % resample indices with replacement
    x_boot = x(idx);
    y_boot = y(idx);
    r_boot(i) = corr(x_boot, y_boot, 'Type', 'Pearson');
end

% Two-tailed p-value: how many bootstrap correlations are more extreme than observed
p_boot = mean(abs(r_boot) >= abs(r_obs));

% Display results
fprintf('Observed correlation (r) = %.3f\n', r_obs);
fprintf('Bootstrap p-value = %.4f\n', p_boot);



st_corr = -st_corr;

st_corr1 = zeros(11,1);

for i=1:11
    j                =   104*(i-1)+1;
    k                =   104*i;


    series1          =  night_cdhw_severity_mx_mean(j:k,1);

    series2          =  all_data_anm1_mx_mean(j:k,1);

    a = corrcoef(series1, series2);
        
        % Find the lag with the highest absolute correlation
    
     st_corr1(i,1)    =  a(2,1);
end

% Assuming dayMatrix and nightMatrix are your input matrices with dimensions 2256x1144
dayMatrix    = (mean(day_cdhw_severity(index_AZ,:),1))'; % Example data for day
nightMatrix  = (mean(night_cdhw_severity(index_MX,:),1))'; % Example data for night

% Initialize matrices to store the highest correlation and corresponding lag
highestCorrelation_AZMX    = zeros(44,1);
lagWithMaxCorrelation_AZMX = zeros(44,1);

% Time lag range from -3 weeks to +3 weeks, assuming 1 unit is 1 day
maxLagDays = 3; % 3 weeks

% Loop through each combination of columns (pairwise)

        for k=1:44
            m = 26*(k-1)+1;
            n = 26*k;

       
        % Extract the time series for the current pair
        series1 = dayMatrix(m:n,1);
        series2 = nightMatrix(m:n,1);
        
        % Compute cross-correlation for this pair with specified max lag
        [crossCorr, lags] = xcorr(series1, series2, maxLagDays, 'coeff');
        
        [posCorr,posIdx] = max(crossCorr);               % ≥ 0
        [negCorr,negIdx] = min(crossCorr);               % ≤ 0
        posLag = lags(posIdx);
        negLag = lags(negIdx);

        % keep whichever has the larger |r|
        if abs(posCorr) >= abs(negCorr)
            highestCorrelation_AZMX(k)      =  posCorr;
            lagWithMaxCorrelation_AZMX(k)   =  posLag;
        else
            highestCorrelation_AZMX(k)      =  negCorr;
            lagWithMaxCorrelation_AZMX(k)   =  negLag;
        end
        
        
        end
    
   
[a,b] = mann_kendall(highestCorrelation_AZMX);

% Compute mean across time (transpose to get column vectors)
dayMatrix   = (mean(day_cdhw_severity(index_NM,:), 1))';   % Daytime series
nightMatrix = (mean(night_cdhw_severity(index_MX,:), 1))'; % Nighttime series

% Compute mean across time (transpose to get column vectors)
dayMatrix   = (mean(day_cdhw_severity(index_NM,:), 1))';   % Daytime series
nightMatrix = (mean(night_cdhw_severity(index_NM,:), 1))'; % Nighttime series

% Initialize storage
numPairs = 44;                   % Number of AZ-MX pairs
daysPerPair = 26;                % Days per time series
maxLagDays = 3;                  % Max lag to test, in days
nBoot = 1000;                    % Number of bootstrap samples
alpha = 0.05;                    % Significance threshold

% Result containers
highestCorrelation_AZMX    = zeros(numPairs,1);
lagWithMaxCorrelation_AZMX = zeros(numPairs,1);
p_boot                     = zeros(numPairs,1);   % Bootstrapped p-values

% Loop through each AZ-MX pair
for k = 1:numPairs
    m = daysPerPair * (k-1) + 1;
    n = daysPerPair * k;

    % Extract time series
    series1 = dayMatrix(m:n);
    series2 = nightMatrix(m:n);

    % Cross-correlation
    [crossCorr, lags] = xcorr(series1, series2, maxLagDays, 'coeff');
    [posCorr, posIdx] = max(crossCorr);
    [negCorr, negIdx] = min(crossCorr);
    posLag = lags(posIdx);
    negLag = lags(negIdx);

    % Pick the stronger correlation
    if abs(posCorr) >= abs(negCorr)
        bestCorr = posCorr;
        bestLag  = posLag;
    else
        bestCorr = negCorr;
        bestLag  = negLag;
    end

    % Store actual result
    highestCorrelation_AZMX(k)      = bestCorr;
    lagWithMaxCorrelation_AZMX(k)   = bestLag;

    % ==== Bootstrapping for significance ====
    bootCorrs = zeros(nBoot,1);
    for b = 1:nBoot
        shuffled = series2(randperm(length(series2)));  % Shuffle series2
        [bootC, ~] = xcorr(series1, shuffled, maxLagDays, 'coeff');
        bootCorrs(b) = max(abs(bootC));  % Max absolute cross-corr under null
    end

    % Bootstrapped p-value (2-tailed)
    p_boot(k) = sum(abs(bootCorrs) >= abs(bestCorr)) / nBoot;
end


% === Mann-Kendall Trend Test ===
[tau, p_mk] = mann_kendall(highestCorrelation_AZMX);

% === Report ===
fprintf('\n--- Mann-Kendall Test on Peak Correlations ---\n');
fprintf('Kendall''s tau = %.4f\n', tau);
fprintf('p-value = %.4f\n', p_mk);
if p_mk < alpha
    fprintf('Trend is statistically significant (p < %.2f)\n', alpha);
else
    fprintf('No significant trend detected (p ≥ %.2f)\n', alpha);
end

% === Display significant cross-correlations ===
significantIdx = find(p_boot < alpha);
fprintf('\n--- Significant Cross-Correlations (p < %.2f) ---\n', alpha);
disp(significantIdx);

% === Optional Visualization ===
figure;
bar(p_boot);
yline(alpha, 'r--', 'Significance Threshold');
xlabel('Pair Index');
ylabel('Bootstrapped p-value');
title('Significance of Peak Cross-Correlations via Bootstrapping');





% Assuming dayMatrix and nightMatrix are your input matrices with dimensions 2256x1144
dayMatrix    = (mean(day_cdhw_severity(index_AZ,:),1))'; % Example data for day
nightMatrix  = (mean(night_cdhw_severity(index_AZ,:),1))'; % Example data for night

% Initialize matrices to store the highest correlation and corresponding lag
highestCorrelation_AZAZ    = zeros(44,1);
lagWithMaxCorrelation_AZAZ = zeros(44,1);

% Time lag range from -3 weeks to +3 weeks, assuming 1 unit is 1 day
maxLagDays = 3; % 3 weeks

% Loop through each combination of columns (pairwise)

        for k=1:44
            m = 26*(k-1)+1;
            n = 26*k;

       
        % Extract the time series for the current pair
        series1 = dayMatrix(m:n,1);
        series2 = nightMatrix(m:n,1);
        
        % Compute cross-correlation for this pair with specified max lag
        [crossCorr, lags] = xcorr(series1, series2, maxLagDays, 'coeff');
        
        [posCorr,posIdx] = max(crossCorr);               % ≥ 0
        [negCorr,negIdx] = min(crossCorr);               % ≤ 0
        posLag = lags(posIdx);
        negLag = lags(negIdx);

        % keep whichever has the larger |r|
        if abs(posCorr) >= abs(negCorr)
            highestCorrelation_AZAZ(k)      =  posCorr;
            lagWithMaxCorrelation_AZAZ(k)   =  posLag;
        else
            highestCorrelation_AZAZ(k)      =  negCorr;
            lagWithMaxCorrelation_AZAZ(k)   =  negLag;
        end
        
        
        end
    
   
[a,b] = mann_kendall(highestCorrelation_AZAZ);


% Assuming dayMatrix and nightMatrix are your input matrices with dimensions 2256x1144
dayMatrix    = (mean(day_cdhw_severity(index_MX,:),1))'; % Example data for day
nightMatrix  = (mean(night_cdhw_severity(index_MX,:),1))'; % Example data for night

% Initialize matrices to store the highest correlation and corresponding lag
highestCorrelation_MXMX    = zeros(44,1);
lagWithMaxCorrelation_MXMX = zeros(44,1);

% Time lag range from -3 weeks to +3 weeks, assuming 1 unit is 1 day
maxLagDays = 3; % 3 weeks

% Loop through each combination of columns (pairwise)

        for k=1:44
            m = 26*(k-1)+1;
            n = 26*k;

       
        % Extract the time series for the current pair
        series1 = dayMatrix(m:n,1);
        series2 = nightMatrix(m:n,1);
        
        % Compute cross-correlation for this pair with specified max lag
        [crossCorr, lags] = xcorr(series1, series2, maxLagDays, 'coeff');
        
        [posCorr,posIdx] = max(crossCorr);               % ≥ 0
        [negCorr,negIdx] = min(crossCorr);               % ≤ 0
        posLag = lags(posIdx);
        negLag = lags(negIdx);

        % keep whichever has the larger |r|
        if abs(posCorr) >= abs(negCorr)
            highestCorrelation_MXMX(k)      =  posCorr;
            lagWithMaxCorrelation_MXMX(k)   =  posLag;
        else
            highestCorrelation_MXMX(k)      =  negCorr;
            lagWithMaxCorrelation_MXMX(k)   =  negLag;
        end
        
        
        end
    
   
[a,b] = mann_kendall(highestCorrelation_MXMX);


% Assuming dayMatrix and nightMatrix are your input matrices with dimensions 2256x1144
dayMatrix    = (mean(day_cdhw_severity(index_NM,:),1))'; % Example data for day
nightMatrix  = (mean(night_cdhw_severity(index_MX,:),1))'; % Example data for night

% Initialize matrices to store the highest correlation and corresponding lag
highestCorrelation_NMMX    = zeros(44,1);
lagWithMaxCorrelation_NMMX = zeros(44,1);

% Time lag range from -3 weeks to +3 weeks, assuming 1 unit is 1 day
maxLagDays = 3; % 3 weeks

% Loop through each combination of columns (pairwise)

        for k=1:44
            m = 26*(k-1)+1;
            n = 26*k;

       
        % Extract the time series for the current pair
        series1 = dayMatrix(m:n,1);
        series2 = nightMatrix(m:n,1);
        
        % Compute cross-correlation for this pair with specified max lag
        [crossCorr, lags] = xcorr(series1, series2, maxLagDays, 'coeff');
        
        [posCorr,posIdx] = max(crossCorr);               % ≥ 0
        [negCorr,negIdx] = min(crossCorr);               % ≤ 0
        posLag = lags(posIdx);
        negLag = lags(negIdx);

        % keep whichever has the larger |r|
        if abs(posCorr) >= abs(negCorr)
            highestCorrelation_NMMX(k)      =  posCorr;
            lagWithMaxCorrelation_NMMX(k)   =  posLag;
        else
            highestCorrelation_NMMX(k)      =  negCorr;
            lagWithMaxCorrelation_NMMX(k)   =  negLag;
        end
        
        
        end
    
   
[a,b] = mann_kendall(highestCorrelation_NMMX);

% Assuming dayMatrix and nightMatrix are your input matrices with dimensions 2256x1144
dayMatrix    = (mean(day_cdhw_severity(index_NM,:),1))'; % Example data for day
nightMatrix  = (mean(night_cdhw_severity(index_NM,:),1))'; % Example data for night

% Initialize matrices to store the highest correlation and corresponding lag
highestCorrelation_NMNM    = zeros(44,1);
lagWithMaxCorrelation_NMNM = zeros(44,1);

% Time lag range from -3 weeks to +3 weeks, assuming 1 unit is 1 day
maxLagDays = 3; % 3 weeks

% Loop through each combination of columns (pairwise)

        for k=1:44
            m = 26*(k-1)+1;
            n = 26*k;

       
        % Extract the time series for the current pair
        series1 = dayMatrix(m:n,1);
        series2 = nightMatrix(m:n,1);
        
        % Compute cross-correlation for this pair with specified max lag
        [crossCorr, lags] = xcorr(series1, series2, maxLagDays, 'coeff');
        
        [posCorr,posIdx] = max(crossCorr);               % ≥ 0
        [negCorr,negIdx] = min(crossCorr);               % ≤ 0
        posLag = lags(posIdx);
        negLag = lags(negIdx);

        % keep whichever has the larger |r|
        if abs(posCorr) >= abs(negCorr)
            highestCorrelation_NMNM(k)      =  posCorr;
            lagWithMaxCorrelation_NMNM(k)   =  posLag;
        else
            highestCorrelation_NMNM(k)      =  negCorr;
            lagWithMaxCorrelation_NMNM(k)   =  negLag;
        end
        
        
        end
    
   
[a,b] = mann_kendall(highestCorrelation_NMNM);



% Assuming dayMatrix and nightMatrix are your input matrices with dimensions 2256x1144
dayMatrix    = (mean(day_cdhw_severity(index_AZ,:),1))'; % Example data for day
nightMatrix  = (mean(night_cdhw_severity(index_NM,:),1))'; % Example data for night

% Initialize matrices to store the highest correlation and corresponding lag
highestCorrelation_AZNM    = zeros(44,1);
lagWithMaxCorrelation_AZNM = zeros(44,1);

% Time lag range from -3 weeks to +3 weeks, assuming 1 unit is 1 day
maxLagDays = 3; % 3 weeks

% Loop through each combination of columns (pairwise)

        for k=1:44
            m = 26*(k-1)+1;
            n = 26*k;

       
        % Extract the time series for the current pair
        series1 = dayMatrix(m:n,1);
        series2 = nightMatrix(m:n,1);
        
        % Compute cross-correlation for this pair with specified max lag
        [crossCorr, lags] = xcorr(series1, series2, maxLagDays, 'coeff');
        
        [posCorr,posIdx] = max(crossCorr);               % ≥ 0
        [negCorr,negIdx] = min(crossCorr);               % ≤ 0
        posLag = lags(posIdx);
        negLag = lags(negIdx);

        % keep whichever has the larger |r|
        if abs(posCorr) >= abs(negCorr)
            highestCorrelation_AZNM(k)      =  posCorr;
            lagWithMaxCorrelation_AZNM(k)   =  posLag;
        else
            highestCorrelation_AZNM(k)      =  negCorr;
            lagWithMaxCorrelation_AZNM(k)   =  negLag;
        end
        
        
        end
    
   
[a,b] = mann_kendall(highestCorrelation_AZNM);

day_night_corr = {highestCorrelation_MXMX, highestCorrelation_AZAZ, highestCorrelation_NMNM, ...
                  highestCorrelation_AZMX, highestCorrelation_NMMX, highestCorrelation_AZNM};




titles = {'MXMX', 'AZAZ', 'NMNM', 'AZMX', 'NMMX', 'AZNM'};
subplot_positions = [1, 3, 5, 2, 4, 6];  % Column-wise layout
years = (1980:2023)';

figure('Color', 'w', 'Position', [100 100 1000 800]);

for i = 1:6
    subplot(3, 2, subplot_positions(i));

    % Clean NaNs
    x = day_night_corr{i};
    valid = ~isnan(x);
    x_clean = x(valid);
    years_clean = years(valid);

    % Linear trend
    coeffs = polyfit(years_clean, x_clean, 1);
    slope = coeffs(1);
    yfit = polyval(coeffs, years_clean);

    % Mann-Kendall trend
    [tau, p] = mann_kendall(x_clean);

    % Plot original data
    scatter(years_clean, x_clean, 40, 'k', 'filled'); hold on;

    % Plot trend line (color by slope sign)
    if slope > 0
        lineColor = [0.85 0.1 0.1]; % Deep red
    else
        lineColor = [0.1 0.1 0.85]; % Deep blue
    end
    plot(years_clean, yfit, '--', 'Color', lineColor, 'LineWidth', 2);

    % Annotations
    txt = sprintf('\\tau = %.2f\np = %.3f\n\\Delta = %.4f/yr', tau, p, slope);
    annotationX = years_clean(round(end*0.05));
    annotationY = max(x_clean) - 0.05; % near top
    text(annotationX, annotationY, txt, 'BackgroundColor', 'w', ...
        'EdgeColor', 'k', 'Margin', 5, 'FontSize', 20, 'FontName', 'Times New Roman');

    title(titles{i}, 'FontSize', 13, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
    xlim([1980 2023]);
    ylim([-1 1]);
    yticks(-1:0.5:1);
    xlabel('Year');
    ylabel('Correlation');
    grid on; box off;
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 20, 'LineWidth', 2, ...
        'TickDir', 'out', 'YGrid', 'on', 'XMinorTick', 'on', 'YMinorTick', 'on');
end

sgtitle('Trend of Day–Night Correlation (1980–2023)', ...
    'FontSize', 15, 'FontWeight', 'bold', 'FontName', 'Times New Roman');







% Assuming dayMatrix and nightMatrix are your input matrices with dimensions 2256x1144
dayMatrix    = day_cdhw_severity; % Example data for day
nightMatrix  = night_cdhw_severity; % Example data for night

% Initialize matrices to store the highest correlation and corresponding lag
highestCorrelation = zeros(2256, 2256,44);
lagWithMaxCorrelation = zeros(2256, 2256,44);

% Time lag range from -3 weeks to +3 weeks, assuming 1 unit is 1 day
maxLagDays = 3; % 3 weeks

% Loop through each combination of columns (pairwise)
parfor i = 1:2256
    for j = 1:2256
        for k=1:44
            m = 26*(k-1)+1;
            n = 26*k;

       
        % Extract the time series for the current pair
        series1 = (dayMatrix(i,m:n))';
        series2 = (nightMatrix(j,m:n))';
        
        % Compute cross-correlation for this pair with specified max lag
        [crossCorr, lags] = xcorr(series1, series2, maxLagDays, 'coeff');
        
        [posCorr,posIdx] = max(crossCorr);               % ≥ 0
        [negCorr,negIdx] = min(crossCorr);               % ≤ 0
        posLag = lags(posIdx);
        negLag = lags(negIdx);

        % keep whichever has the larger |r|
        if abs(posCorr) >= abs(negCorr)
            highestCorrelation(i,j,k)      =  posCorr;
            lagWithMaxCorrelation(i,j,k)   =  posLag;
        else
            highestCorrelation(i,j,k)      =  negCorr;
            lagWithMaxCorrelation(i,j,k)   =  negLag;
        end
        
        
        end
    
    end
    i
end

%save('yearly_HD_Dependence.mat','lagWithMaxCorrelation','highestCorrelation','-v7.3');

dayMatrix_AZ     = dayMatrix(index_AZ,:);
nightMatrix_MX   = nightMatrix(index_MX,:);
% Time lag range from -3 weeks to +3 weeks, assuming 1 unit is 1 day
maxLagDays = 3; % 3 weeks
highest_correlation_AZ_MX_rand = zeros(1000,44);
% Loop through each combination of columns (pairwise)
for p = 1:1000
rng(p);highestCorrelation_ran = zeros(113, 700,44);
lagWithMaxCorrelation_ran = zeros(113, 700,44);
for i = 1:113
    for j = 1:700
        for k=1:44
            m = 26*(k-1)+1;
            n = 26*k;

       
        % Extract the time series for the current pair
        series1   = (dayMatrix_AZ(i,m:n))';
        series2   = (nightMatrix_MX(j,m:n))';

        ser_size  = size(series1,1);
        T         = randperm(ser_size);

        [series1_ran, series2_ran] = deal(series1(T,:), series2(T,:));
        [crossCorr, lags] = xcorr(series1_ran, series2_ran, maxLagDays, 'coeff');

        
        [posCorr,posIdx] = max(crossCorr);               % ≥ 0
        [negCorr,negIdx] = min(crossCorr);               % ≤ 0
        posLag = lags(posIdx);
        negLag = lags(negIdx);

        % keep whichever has the larger |r|
        if abs(posCorr) >= abs(negCorr)
            highestCorrelation_ran(i,j,k)      =  posCorr;
            lagWithMaxCorrelation_ran(i,j,k)   =  posLag;
        else
            highestCorrelation_ran(i,j,k)      =  negCorr;
            lagWithMaxCorrelation_ran(i,j,k)   =  negLag;
        end
        
        
        end
    
    end
    
end


highest_correlation_AZ_ran    = nanmean(highestCorrelation_ran,1);
highest_correlation_AZ_MX_ran = nanmean(highest_correlation_AZ_ran,2);
highest_correlation_AZ_MX_rand(p,:) = (reshape(highest_correlation_AZ_MX_ran,[1*1,44]))';

p
end
save('highest_correlation_AZ_MX_Randomization','highest_correlation_AZ_MX_rand');

%Day-Night Hot Drought Coupling over Arizona and Mexico
highest_correlation_AZ    = nanmean(highestCorrelation(index_AZ,:,:),1);
highest_correlation_AZ_MX = nanmean(highest_correlation_AZ(:,index_MX,:),2);
highest_correlation_AZ_MX = (reshape(highest_correlation_AZ_MX,[1*1,44]))';



highest_correlation_AZ1    = nanmean(highestCorrelation(index_AZ,:,:),1);
highest_correlation_AZ_MX1 = nanmean(highest_correlation_AZ1(:,index_MX,:),2);
highest_correlation_AZ_MX1 = (reshape(highest_correlation_AZ_MX1,[1*1,44]))';

day_night = [highest_correlation_AZ_MX,highest_correlation_AZ_MX1];

% Example vectors
x = day_night(:,1); 
y = day_night(:,2);  % Correlated with noise

% Compute correlation and p-value
[R, P] = corr(x, y, 'Type', 'Pearson');  % You can also use 'Spearman' or 'Kendall'

% Scatter plot
figure;
scatter(x, y, 50, 'filled');
xlabel('x');
ylabel('y');
title('Scatter plot with correlation and p-value');
grid on;

% Annotate with correlation and significance
annotationText = sprintf('r = %.2f, p = %.3g', R, P);
text(min(x) + 0.05*range(x), max(y) - 0.1*range(y), annotationText, ...
    'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', 'w');




highest_correlation_AZ_AZ = nanmean(highest_correlation_AZ(:,index_AZ,:),2);
highest_correlation_AZ_AZ = (reshape(highest_correlation_AZ_AZ,[1*1,44]))';


highest_correlation_MX    =  nanmean(highestCorrelation(index_MX,:,:),1);
highest_correlation_MX_MX = nanmean(highest_correlation_MX(:,index_MX,:),2);
highest_correlation_MX_MX = (reshape(highest_correlation_MX_MX,[1*1,44]))';


%Day-Night Hot Drought Coupling over New Mexico and Mexico
highest_correlation_NM    = nanmean(highestCorrelation(index_NM,:,:),1);
highest_correlation_NM_MX = nanmean(highest_correlation_NM(:,index_MX,:),2);
highest_correlation_NM_MX = (reshape(highest_correlation_NM_MX,[1*1,44]))';

highest_correlation_NM_NM = nanmean(highest_correlation_NM(:,index_NM,:),2);
highest_correlation_NM_NM = (reshape(highest_correlation_NM_NM,[1*1,44]))';


highest_correlation_NM_AZ = nanmean(highest_correlation_NM(:,index_AZ,:),2);
highest_correlation_NM_AZ = (reshape(highest_correlation_NM_AZ,[1*1,44]))';




ONI_DJF = [0.6
-0.3
0
2.2
-0.6
-1
-0.5
1.2
0.8
-1.7
0.1
0.4
1.7
0.1
0.1
1
-0.9
-0.5
2.2
-1.5
-1.7
-0.7
-0.1
0.9
0.4
0.6
-0.9
0.7
-1.6
-0.8
1.5
-1.4
-0.9
-0.4
-0.4
0.5
2.5
-0.3
-0.9
0.7
0.5
-1
-1
-0.7];

ONI_JJA = [0.3
-0.3
0.8
0.3
-0.3
-0.5
0.2
1.5
-1.3
-0.3
0.3
0.7
0.4
0.3
0.4
-0.2
-0.3
1.6
-0.8
-1.1
-0.6
-0.1
0.8
0.1
0.5
-0.1
0.1
-0.6
-0.4
0.5
-1
-0.5
0.2
-0.4
0
1.5
-0.4
0.1
0.1
0.3
-0.4
-0.4
-0.8
1.1];

y  = highest_correlation_AZ_MX(:);            % 44-by-1 column
tt = datetime((1980:2023)',1,1);              % years → datetime
x  = year(tt);                                % numeric predictor


mdl   = fitlm(x, y, 'RobustOpts', 'on');      % bisquare weighting
b0    = mdl.Coefficients.Estimate(1);         % intercept
b1    = mdl.Coefficients.Estimate(2);         % slope (units yr-1)
CI    = coefCI(mdl);                          % 95 % CI
p_val = mdl.Coefficients.pValue(2);           % p-value for slope
yhat  = b0 + b1*x;                            % fitted trend


f = figure('Units','centimeters','Position',[5 5 16 8], ...
           'Color','w','Renderer','painters');
set(f,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',8, ...
      'DefaultLineLineWidth',1.2,'DefaultAxesLineWidth',0.8, ...
      'DefaultAxesTickDir','out','DefaultAxesBox','off');

plot(tt, y, '-o', ...
     'MarkerSize',3, 'MarkerFaceColor',[.35 .35 .35], ...
     'Color',[.35 .35 .35]);                      % observations
hold on
plot(tt, yhat, 'k--', 'LineWidth',1.6);           % robust trend
hold off

xlabel('Year');
ylabel('Correlation');
title('Annual Time Series with Robust Trend');
legend({'Observations','Robust linear trend'}, ...
       'Location','best','Box','off');

% Inline statistics
txt = sprintf(['Slope = %.3f ± %.3f yr^{-1}\\newline' ...
               'p = %.3g  (95%% CI)'], ...
               b1, 0.5*diff(CI(2,:)), p_val);
annotation('textbox',[.15 .72 .3 .15],'String',txt, ...
           'FitBoxToText','on','EdgeColor','none','FontSize',20);

% ------------------------------------------------------------
% Export
% ------------------------------------------------------------
print(f,'-dpdf','-painters','time_series_robust_trend.pdf');


% Prepare Data
x = datetime((1980:2023)',1,1);
years = year(x);

data = [highest_correlation_AZ_AZ(:), highest_correlation_MX_MX(:), highest_correlation_AZ_MX(:)];
labels = {"AZ Day–AZ Night", "MX Day–MX Night", "AZ Day–MX Night"};

% Colors
clr_linear = [0.3 0.3 0.3];   % Gray for linear trend
colors = [0.25,0.41,0.88; 0.2,0.6,0.2; 0.9,0.3,0.3];  % Series colors

% LOESS smoothing parameter
span = 0.3;

% Initialize storage
p_MK = zeros(3,1);
slope = zeros(3,1);
CI = zeros(3,2);
p_linear = zeros(3,1);
yhat_linear = zeros(44,3);
loess_trend = zeros(44,3);

for i = 1:3
    % Robust Linear Regression
    mdl = fitlm(years, data(:,i), 'RobustOpts','on');
    slope(i) = mdl.Coefficients.Estimate(2);
    CI_temp = coefCI(mdl);
    CI(i,:) = CI_temp(2,:);
    p_linear(i) = mdl.Coefficients.pValue(2);
    yhat_linear(:,i) = mdl.Coefficients.Estimate(1) + slope(i)*years;
    
    % LOESS Trend
    loess_trend(:,i) = smooth(years, data(:,i), span, 'loess');
    
    % Mann-Kendall Test (Assuming function exists)
    [~, p_MK(i)] = mann_kendall(data(:,i));   % Replace with actual function
end

% --- Plot ---
figure('Units','centimeters','Position',[5 5 15 18],'Color','w');
tiledlayout(3,1,'Padding','compact','TileSpacing','compact');

for i = 1:3
    ax = nexttile;
    hold on
    
    % Raw Data
    plot(x, data(:,i), 'o-', 'Color', colors(i,:), 'MarkerFaceColor', colors(i,:), 'MarkerSize', 3);
    
    % LOESS Trend
    plot(x, loess_trend(:,i), '-', 'Color', colors(i,:), 'LineWidth', 1.5);
    
    % Linear Trend
    plot(x, yhat_linear(:,i), '--', 'Color', clr_linear, 'LineWidth', 1.2);
    
    ylabel(labels{i});
    set(gca,'FontName','Helvetica','FontSize',8,'Box','off','TickDir','out');
    xlim([x(1) x(end)]);
    
    % Annotation
    txt = sprintf(['Linear: slope = %.3f ± %.3f (p = %.3g)\\newline' ...
                   'LOESS: span = %.2f, M-K p = %.3g'], ...
                   slope(i), 0.5*diff(CI(i,:)), p_linear(i), span, p_MK(i));
    
    text(0.02, 0.80, txt, 'Units','normalized', 'FontSize',7,'FontName','Helvetica');
    
    if i == 3
        xlabel('Year');
    else
        set(gca,'XTickLabel',[]);
    end
end

% --- Export ---
print('-dpdf','-painters','linear_loess_trends_significance.pdf');


% --- Prepare Data ---
x = datetime((1980:2023)',1,1);
years = year(x);

data = [highest_correlation_AZ_AZ(:),highest_correlation_AZ_MX(:), highest_correlation_MX_MX(:),highest_correlation_NM_MX(:), highest_correlation_NM_NM(:),highest_correlation_NM_AZ(:)];

labels = {"AZ Day–AZ Night","AZ Day–MX Night", "MX Day–MX Night","NM Day–MX Night", "NM Day–NM Night", "AZ Day–NM Night"};

% Colors
clr_linear = [0.3 0.3 0.3];   % Gray for linear trend
colors = [0.25,0.41,0.88;     % AZ-AZ
          0.2,0.6,0.2;        % MX-MX
          0.85,0.33,0.1;      % NM-NM
          0.9,0.3,0.3;        % AZ-MX
          0.5,0.2,0.8;        % NM-MX
          0.1,0.7,0.7];       % NM-AZ

% LOESS smoothing parameter
span = 0.3;

% Initialize storage
p_MK = zeros(6,1);
slope = zeros(6,1);
CI = zeros(6,2);
p_linear = zeros(6,1);
yhat_linear = zeros(44,6);
loess_trend = zeros(44,6);

for i = 1:6
    % Robust Linear Regression
    mdl = fitlm(years, data(:,i), 'RobustOpts','on');
    slope(i) = mdl.Coefficients.Estimate(2);
    CI_temp = coefCI(mdl);
    CI(i,:) = CI_temp(2,:);
    p_linear(i) = mdl.Coefficients.pValue(2);
    yhat_linear(:,i) = mdl.Coefficients.Estimate(1) + slope(i)*years;
    
    % LOESS Trend
    loess_trend(:,i) = smooth(years, data(:,i), span, 'loess');
    
    % Mann-Kendall Test (Assuming function exists)
    [~, p_MK(i)] = mann_kendall(data(:,i));   % Replace with actual function
end

% --- Plot ---
figure('Units','centimeters','Position',[5 5 18 18],'Color','w');
tiledlayout(3,2,'Padding','compact','TileSpacing','compact');

for i = 1:6
    ax = nexttile(i);
    hold on
    
    % Raw Data
    plot(x, data(:,i), 'o-', 'Color', colors(i,:), 'MarkerFaceColor', colors(i,:), 'MarkerSize', 3);
    
    % LOESS Trend
    plot(x, loess_trend(:,i), '-', 'Color', colors(i,:), 'LineWidth', 1.5);
    
    % Linear Trend
    plot(x, yhat_linear(:,i), '--', 'Color', clr_linear, 'LineWidth', 1.2);
    
    ylabel(labels{i});
    set(gca,'FontName','Helvetica','FontSize',8,'Box','off','TickDir','out');
    xlim([x(1) x(end)]);
    
    % Annotation
    txt = sprintf(['Linear: slope = %.3f ± %.3f (p = %.3g)\\newline' ...
                   'LOESS: span = %.2f, M-K p = %.3g'], ...
                   slope(i), 0.5*diff(CI(i,:)), p_linear(i), span, p_MK(i));
    
    text(0.02, 0.80, txt, 'Units','normalized', 'FontSize',7,'FontName','Helvetica');
    
    if i > 4
        xlabel('Year');
    else
        set(gca,'XTickLabel',[]);
    end
end

% --- Export ---
print('-dpdf','-painters','linear_loess_trends_significance_2col.pdf');






correlation = [st_corr, highest_correlation_AZ_MX];

f = figure('Color','w','Units','centimeters','Position',[5 5 14 7]);
set(f,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',8, ...
      'DefaultLineLineWidth',1.2,'DefaultAxesLineWidth',0.8, ...
      'DefaultAxesTickDir','out','DefaultAxesBox','off');
t=(1:11)';
yyaxis left                       % left y–axis
plot(t, correlation(:,1), '-o', 'MarkerSize',3);
ylabel('Series 1  (units)');
ylim padded                       % keeps a little margin

yyaxis right                      % right y–axis
plot(t, correlation(:,2), '-s', 'MarkerSize',3);
ylabel('Series 2  (units)');
ylim padded

xlabel('Time');
grid on
legend({'Series 1','Series 2'},'Location','northwest');




correlation = [st_corr , highest_correlation_AZ_MX];   % [11×2]
barSeries   = normalize(weeklyDataMatrix_anm_mx);                             % <-- new (length 11)

f = figure('Color','w','Units','centimeters','Position',[5 5 14 7]);
set(f,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',8, ...
      'DefaultLineLineWidth',1.2,'DefaultAxesLineWidth',0.8, ...
      'DefaultAxesTickDir','out','DefaultAxesBox','off');

t = (1:11)';

% -------- LEFT axis ------------------------------------------------------
yyaxis left
bar(t, barSeries, ...                       % <-- NEW: bar overlay
    'FaceColor',[0.75 0.75 0.75], ...       %   light‑gray bars
    'EdgeColor','none', ...
    'BarWidth',0.5);                        %   a little narrower
hold on                                     %   keep the axis active
plot(t, correlation(:,1), '-o', 'MarkerSize',3, 'LineWidth',1.2);
hold off
ylabel('Series 1 / Bar series  (units)');
ylim padded
set(gca,'Layer','top');                     % draw gridlines behind bars

% -------- RIGHT axis -----------------------------------------------------
yyaxis right
plot(t, correlation(:,2), '-s', 'MarkerSize',3, 'LineWidth',1.2);
ylabel('Series 2  (units)');
ylim padded

% -------- decorations ----------------------------------------------------
xlabel('Time');
grid on
legend({'Bar series','Series 1','Series 2'}, ...
       'Location','northwest');



% --- data ---------------------------------------------------------------
correlation  = [st_corr , highest_correlation_AZ_MX];   % [11×2]
barSeries1   = normalize(weeklyDataMatrix_anm_mx);      % first bar series
barSeries2   = normalize(all_data_anm1_mx_mean_y1);     % second bar series

t = (1:11)';

% --- figure & styling ---------------------------------------------------
f = figure('Color','w','Units','centimeters','Position',[5 5 14 7]);
set(f,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',8, ...
      'DefaultLineLineWidth',1.2,'DefaultAxesLineWidth',0.8, ...
      'DefaultAxesTickDir','out','DefaultAxesBox','off');

% -------- LEFT axis: two grouped bars + line ----------------------------
yyaxis left

% grouped bars: Y must be an N×2 matrix
b = bar(t, [barSeries1 , barSeries2], 'grouped');  % just the data + style
b(1).BarWidth = 0.7;
b(2).BarWidth = 0.7;
b(1).FaceColor = [0.75 0.75 0.75];   % light gray
b(1).EdgeColor = 'none';
b(2).FaceColor = [0.15 0.55 0.80];   % blue‑ish
b(2).EdgeColor = 'none';

hold on
plot(t, correlation(:,1), '-o', 'MarkerSize',3, 'LineWidth',1.2);
hold off

ylabel('Series 1  /  Bar series (normalized units)');
ylim padded
set(gca,'Layer','top');              % grid behind bars

% -------- RIGHT axis: second line ---------------------------------------
yyaxis right
plot(t, correlation(:,2), '-s', 'MarkerSize',3, 'LineWidth',1.2);
ylabel('Series 2  (units)');
ylim padded

% -------- decorations ----------------------------------------------------
xlabel('Time');
grid on
legend({'Bar A','Bar B','Series 1','Series 2'}, ...
       'Location','northwest');



%Caluclating Dependence using Mutual Information

% Initialize matrices to store the highest mutual information and corresponding lag
AZ_data = day_cdhw_severity(index_AZ,:);
MX_data = night_cdhw_severity(index_MX,:);

highestMI = zeros(113, 700, 11);
lagWithMaxMI = zeros(113, 700,11);

% Time lag range from -3 weeks to +3 weeks, assuming 1 unit is 1 day
maxLagDays = 3; % 3 weeks in days, assuming a weekly resolution (7 days per week)

% Define the range of lags
lagRange = -maxLagDays:maxLagDays;

% Loop through each combination of columns (pairwise)
parfor i = 1:113
    for j = 1:700
        for k = 1:11
            m = 104 * (k - 1) + 1;
            n = 104 * k;

            % Extract the time series for the current pair
            series1 = (AZ_data(i, m:n))';
            series2 = (MX_data(j, m:n))';

            % Initialize variables to track the maximum MI and corresponding lag
            maxMI = -inf;
            bestLag = 0;
            
            % Iterate over the lag range (-3 weeks to +3 weeks)
            for lag = lagRange
                % Shift series2 by the current lag
                if lag < 0
                    shiftedSeries2 = series2(1:end+lag); % lag < 0, shift left
                    truncatedSeries1 = series1(-lag+1:end);
                elseif lag > 0
                    shiftedSeries2 = series2(lag+1:end); % lag > 0, shift right
                    truncatedSeries1 = series1(1:end-lag);
                else
                    % No lag
                    shiftedSeries2 = series2;
                    truncatedSeries1 = series1;
                end
                
                % Calculate mutual information for the current lag
                mim = mi(truncatedSeries1, shiftedSeries2); % Use the appropriate MI function
                
                % Store the maximum MI and corresponding lag
                if mim > maxMI
                    maxMI = mim;
                    bestLag = lag;
                end
            end
            
            % Store the highest MI and corresponding lag for this pair
            highestMI(i, j, k) = maxMI;
            lagWithMaxMI(i, j, k) = bestLag;
        end
    end
    i
end

highestMI_AZ    = nanmean(highestMI,1);
highestMI_AZ_MX = nanmean(highestMI_AZ,2);
highestMI_AZ_MX = (reshape(highestMI_AZ_MX,[1*1,11]))';

ncorrelation = [st_corr,highestMI_AZ_MX];

%Analyzing NARR data

ncdisp('/Volumes/Expansion/NARR Data/Soil Moisture/soilm.1980.nc');

ncPath = '/Volumes/Expansion/NARR Data/Soil Moisture/soilm.1980.nc';
lats = ncread(ncPath,'lat');   % size: [349 x 277]
lons = ncread(ncPath,'lon');   % size: [349 x 277]
soilm = ncread(ncPath,'soilm');
time_s = ncread(ncPath,'time');
targetName = 'Mexico';
idx_mx = find(strcmp({S.NAME}, targetName));
if isempty(idx_mx)
    error('No polygon found for NAME="%s". Check shapefile attributes.', targetName);
end
polyX = S(idx_mx(1)).Lon; 
polyY = S(idx_mx(1)).Lat;  
[IN, ON] = inpolygon(lons, lats, polyX, polyY); 
mask = IN | ON;  
insideLinearIdx = find(mask);
[I_grid, J_grid] = find(mask);
soilm_mx = soilm;  
nTime = size(soilm, 3); 
for t = 1:nTime
    frame = soilm_mx(:,:,t);
    frame(~mask) = NaN; 
    soilm_mx(:,:,t) = frame;
end
soilm_mxr     = reshape(soilm_mx,[349*277,2928]);
soilm_mxr(isnan(soilm_mxr)) = 0;
soilm_mxr_sum = sum(soilm_mxr,2);

NARR_idmx = find(soilm_mxr_sum~=0);

ncPath = '/Volumes/Expansion/NARR Data/Soil Moisture/soilm.1980.nc';
lats   = ncread(ncPath,'lat');      % size: [349 x 277]
lons   = ncread(ncPath,'lon');      % size: [349 x 277]
soilm  = ncread(ncPath,'soilm');    % size: [349 x 277 x 2928] (example)
time_s = ncread(ncPath,'time');
nTime  = size(soilm, 3);

regionNames = {'Mexico','Arizona','New Mexico'};
shortNames  = {'MX','AZ','NM'};
NARR_id = struct(); 

for r = 1:numel(regionNames)
    targetName = regionNames{r};
    targetName1 = shortNames{r};
    idx_region = find(strcmp({S.NAME}, targetName));
    if isempty(idx_region)
        error('No polygon found for NAME="%s". Check shapefile attributes.', targetName);
    end
    polyX = S(idx_region(1)).Lon; 
    polyY = S(idx_region(1)).Lat;
    [IN, ON] = inpolygon(lons, lats, polyX, polyY);
    mask     = IN | ON;
    
    soilm_region = soilm;  
    for t = 1:nTime
        frame = soilm_region(:,:,t);
        frame(~mask) = NaN;
        soilm_region(:,:,t) = frame;
    end
    soilm_region_reshaped = reshape(soilm_region, [numel(lats), nTime]);
    soilm_region_reshaped(isnan(soilm_region_reshaped)) = 0;
    soilm_region_sum = sum(soilm_region_reshaped, 2);
    idx_nonzero = find(soilm_region_sum ~= 0);
    NARR_id.(targetName1) = idx_nonzero;
    
    fprintf('Found %d nonzero indices for %s.\n', numel(idx_nonzero), targetName);
end




%Soil Moisture Data Extraction for All Years
dataDir = '/Volumes/Expansion/NARR Data/Soil Moisture';
years = 1980:2023;
allSoilmDaily = [];   
allDates      = [];   

for yy = years
    ncPath = fullfile(dataDir, sprintf('soilm.%d.nc', yy));
    lats    = ncread(ncPath, 'lat');      
    lons    = ncread(ncPath, 'lon');      
    soilm3h = ncread(ncPath, 'soilm');    
    timeRaw = ncread(ncPath, 'time');     
    refDate  = datetime(1800, 1, 1, 0, 0, 0);
    timeVar  = refDate + hours(timeRaw);  
   
    startDate = datetime(yy, 4, 1);
    endDate   = datetime(yy, 9, 30);
    
    inRangeMask = (timeVar >= startDate) & (timeVar <= endDate);
    if ~any(inRangeMask)
        fprintf('No data for %d in Apr 1–Sep 30 range.\n', yy);
        continue
    end
    
    timeSub  = timeVar(inRangeMask);
    soilmSub = soilm3h(:,:,inRangeMask);  % [349 x 277 x #3-hourly slices]
    
  
    dateOnly = dateshift(timeSub, 'start', 'day');  
    [uniqueDays, ~, idxDay] = unique(dateOnly);
    nDays = numel(uniqueDays);
    
   
    dailySoilm = nan(size(soilmSub,1), size(soilmSub,2), nDays);
    
    for d = 1:nDays
        dayMask = (idxDay == d);
        dailySoilm(:,:,d) = mean(soilmSub(:,:,dayMask), 3, 'omitnan');
    end
    
    if isempty(allSoilmDaily)
        
        allSoilmDaily = dailySoilm;
        allDates      = uniqueDays;
    else
        allSoilmDaily = cat(3, allSoilmDaily, dailySoilm);
        allDates      = [allDates; uniqueDays];  % vertical concatenation
    end
    
    fprintf('Year %d: appended %d daily slices.\n', yy, nDays);
end
fprintf('Final dataset size: %s\n', mat2str(size(allSoilmDaily)));
fprintf('First date: %s, Last date: %s\n', datestr(allDates(1)), datestr(allDates(end)));

allSoilmDailyr = reshape(allSoilmDaily,[349*277,8052]);

allSoilmDailyr      = reshape(allSoilmDailyr,[96673,183,44]);
allSoilmDailyr_ltm  = mean(allSoilmDailyr(:,:,12:41),3);

allSoilmDailyr_anm = zeros(96673,183,44);

for i=1:44
    j   =   183*(i-1)+1;
    k   =   183*i;

    allSoilmDailyr_anm(:,:,i)=allSoilmDailyr(:,:,i)-allSoilmDailyr_ltm;
end


allSoilmDailyr_anm_2023    = allSoilmDailyr_anm(:,:,44);
allSoilmDailyr_anm_2023    = allSoilmDailyr_anm_2023(:,1:182);
allSoilmDailyr_anm_2023    = mean(reshape(allSoilmDailyr_anm_2023,[96673,26,7]),3);

allSoilmDailyr_anm_2023_az = (mean(allSoilmDailyr_anm_2023(NARR_id.AZ,:)))';

lonsv   = lons(:);
latsv   = lats(:);

latlon_NARR_az =[lonsv(NARR_id.AZ,:),latsv(NARR_id.AZ,:)];

lonPoint       = -110.5091;
latPoint       =  35.55;

% 1) Compute the Euclidean distance (in degrees) from the point to each row.
dists = sqrt( (latlon_NARR_az(:,1) - lonPoint).^2 + (latlon_NARR_az(:,2) - latPoint).^2 );

% 2) Find the index of the smallest distance (the nearest row).
[~, closestIndex] = min(dists);

latlon_NARR_az(218,:);

allSoilmDailyr_az    = allSoilmDailyr(NARR_id.AZ,:,44);
allSoilmDailyr_azs   = (allSoilmDailyr_az(218,:))';

%Getting the station soil moisture data
station_sm           = readmatrix('/Volumes/Expansion/Snowtel Station Data/AZ_soil_moisture_timeseries.csv');

station_index        = find(station_sm==3069);
station_smts         = station_sm(station_index,3);

station_model = [rescale(all_data_yearly_az_2023_s),rescale(allSoilmDailyr_azs),rescale(station_smts(91:273,1))];

% Define date range
startDate = datetime('01-Apr-2025','InputFormat','dd-MMM-yyyy');
endDate   = datetime('30-Sep-2025','InputFormat','dd-MMM-yyyy');
dates     = (startDate:endDate)';

% Generate sample data for demonstration
% Replace these with your actual data
numDays = length(dates);
Series_A = rescale(all_data_yearly_az_2023_s);
Series_B = rescale(allSoilmDailyr_azs);
Series_C = rescale(station_smts(91:273,1));

% Compute 7-day rolling means
Series_A_7d = movmean(Series_A, 7);
Series_B_7d = movmean(Series_B, 7);
Series_C_7d = movmean(Series_C, 7);

% Create a figure
figure;
hold on;  % Keep all plots on the same axes

% Plot daily data
plot(dates, Series_A, 'DisplayName', 'Series A (Daily)');
plot(dates, Series_B, 'DisplayName', 'Series B (Daily)');
plot(dates, Series_C, 'DisplayName', 'Series C (Daily)');

% Plot 7-day rolling means (dashed lines for distinction)
plot(dates, Series_A_7d, '--', 'DisplayName', 'Series A (7-day mean)');
plot(dates, Series_B_7d, '--', 'DisplayName', 'Series B (7-day mean)');
plot(dates, Series_C_7d, '--', 'DisplayName', 'Series C (7-day mean)');

% Enhance chart appearance
xlabel('Date');
ylabel('Value');
title('Three Time Series (April 1 – Sep 30) with 7-Day Running Mean');
grid on;
legend('Location','best');

% Optional: Format date labels nicely
datetick('x','mmm dd','keepticks','keeplimits');
xtickangle(45);  % rotate date labels for readability if desired

hold off;





ncPath = '/Volumes/Expansion/NARR Data/Sensible Heat Flux/shtfl.1980.nc';
lats   = ncread(ncPath,'lat');      % size: [349 x 277]
lons   = ncread(ncPath,'lon');      % size: [349 x 277]
shtfl  = ncread(ncPath,'shtfl');    % size: [349 x 277 x 2928] (example)
time_s = ncread(ncPath,'time');
nTime  = size(shtfl, 3);

regionNames = {'Mexico','Arizona','New Mexico'};
shortNames  = {'MX','AZ','NM'};
NARR_id = struct(); 

for r = 1:numel(regionNames)
    targetName = regionNames{r};
    targetName1 = shortNames{r};
    idx_region = find(strcmp({S.NAME}, targetName));
    if isempty(idx_region)
        error('No polygon found for NAME="%s". Check shapefile attributes.', targetName);
    end
    polyX = S(idx_region(1)).Lon; 
    polyY = S(idx_region(1)).Lat;
    [IN, ON] = inpolygon(lons, lats, polyX, polyY);
    mask     = IN | ON;
    
    shtfl_region = shtfl;  
    for t = 1:nTime
        frame = shtfl_region(:,:,t);
        frame(~mask) = NaN;
        shtfl_region(:,:,t) = frame;
    end
    shtfl_region_reshaped = reshape(shtfl_region, [numel(lats), nTime]);
    shtfl_region_reshaped(isnan(shtfl_region_reshaped)) = 0;
    shtfl_region_sum = sum(shtfl_region_reshaped, 2);
    idx_nonzero = find(shtfl_region_sum ~= 0);
    NARR_id.(targetName1) = idx_nonzero;
    
    fprintf('Found %d nonzero indices for %s.\n', numel(idx_nonzero), targetName);
end




%Sensible Heat Flux Data Extraction for All Years
dataDir = '/Volumes/Expansion/NARR Data/Sensible Heat Flux';
years = 1980:2023;
allshtflDaily = [];   
allDates      = [];   

for yy = years
    ncPath = fullfile(dataDir, sprintf('shtfl.%d.nc', yy));
    lats    = ncread(ncPath, 'lat');      
    lons    = ncread(ncPath, 'lon');      
    shtfl3h = ncread(ncPath, 'shtfl');    
    timeRaw = ncread(ncPath, 'time');     
    refDate  = datetime(1800, 1, 1, 0, 0, 0);
    timeVar  = refDate + hours(timeRaw);  
   
    startDate = datetime(yy, 4, 1);
    endDate   = datetime(yy, 9, 30);
    
    inRangeMask = (timeVar >= startDate) & (timeVar <= endDate);
    if ~any(inRangeMask)
        fprintf('No data for %d in Apr 1–Sep 30 range.\n', yy);
        continue
    end
    
    timeSub  = timeVar(inRangeMask);
    shtflSub = shtfl3h(:,:,inRangeMask);  % [349 x 277 x #3-hourly slices]
    
  
    dateOnly = dateshift(timeSub, 'start', 'day');  
    [uniqueDays, ~, idxDay] = unique(dateOnly);
    nDays = numel(uniqueDays);
    
   
    dailyshtfl = nan(size(shtflSub,1), size(shtflSub,2), nDays);
    
    for d = 1:nDays
        dayMask = (idxDay == d);
        dailyshtfl(:,:,d) = mean(shtflSub(:,:,dayMask), 3, 'omitnan');
    end
    
    if isempty(allshtflDaily)
        
        allshtflDaily = dailyshtfl;
        allDates      = uniqueDays;
    else
        allshtflDaily = cat(3, allshtflDaily, dailyshtfl);
        allDates      = [allDates; uniqueDays];  % vertical concatenation
    end
    
    fprintf('Year %d: appended %d daily slices.\n', yy, nDays);
end
fprintf('Final dataset size: %s\n', mat2str(size(allSoilmDaily)));
fprintf('First date: %s, Last date: %s\n', datestr(allDates(1)), datestr(allDates(end)));

allshtflDailyr = reshape(allshtflDaily,[349*277,8052]);

allshtflDailyr      = reshape(allshtflDailyr,[96673,183,44]);
allshtflDailyr_ltm  = mean(allshtflDailyr(:,:,12:41),3);

allshtflDailyr_anm = zeros(96673,183,44);

for i=1:44
    j   =   183*(i-1)+1;
    k   =   183*i;

    allshtflDailyr_anm(:,:,i)=allshtflDailyr(:,:,i)-allshtflDailyr_ltm;
end


allshtflDailyr_anm_2023    = allshtflDailyr_anm(:,:,44);
allshtflDailyr_anm_2023    = allshtflDailyr_anm_2023(:,1:182);
allshtflDailyr_anm_2023    = mean(reshape(allshtflDailyr_anm_2023,[96673,26,7]),3);

allshtflDailyr_anm_2023_az = (mean(allshtflDailyr_anm_2023(NARR_id.AZ,:)))';



