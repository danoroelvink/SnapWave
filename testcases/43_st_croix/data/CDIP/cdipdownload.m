% CDIP COMPENDIUM PLOT - MATLAB

% Plot Hs, Tp and Dp for a single station
% See http://cdip.ucsd.edu/themes/cdip?pb=1&bl=cdip?pb=1&d2=p1&d2=p70&u3=s:100:st:1:v:compendium for example Compendium plot.

%-----------------------------------------------------------------------
%%% PRE-CODE DOWNLOAD INSTRUCTIONS %%%
% You must follow these steps in order to access data directly from NetCDF files
%	- Download the nctoolbox from: https://github.com/nctoolbox
%		1. Download and extract into a folder
%		2. Open matlab and at the matlab prompt cd to the extracted folder where setup_nctoolbox lives. (Use pwd to see where you currently are).
%		3. Run setup_nctoolbox
%		4. Optionally edit startup.m and add lines in as described in the install instructions.
%-----------------------------------------------------------------------

% CHANGE MATLAB DIRECTORY TO 'NCTOOLBOX' DIRECTORY
cd nctoolbox

% OPEN NCTOOLBOX IN MATLAB
setup_nctoolbox

% USER ENTERS STATION NUMBER AND START/END DATES FOR PLOT
stn = '067';
startdate = '01/01/2011 00:00';
enddate = '01/31/2011 23:59';

% CONNECT TO THREDDS SERVER AND OPEN NETCDF FILE
urlbase = 'http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/archive/';  % Set 'base' THREDDS URL and pieces to concatenate with user-defined station number (from above)
p1 = 'p1/';
urlend = 'p1_historic.nc';
dsurl = strcat(urlbase,stn,p1,stn,urlend);
ds = ncdataset(dsurl);

% PRINT LIST OF VARIABLES IN NETCDF FILE
% ds.variables

% GET BUOY NAME AND TRANSPOSE TO HORIZONTAL STRING
buoyname = ds.data('metaStationName');
buoytitle = transpose(buoyname);

% EXTRACT MONTH/YEAR FROM 'STARTDATE' 
monthtitle = datestr(datenum(startdate,'mm/dd/yyyy'),'mmm yyyy'); % Create a string object of Month and Year from startdate


% CONVERT START/END DATES TO MATLAB SERIAL NUMBERS
startser = datenum(startdate, 'mm/dd/yyyy HH:MM');
endser = datenum(enddate, 'mm/dd/yyyy HH:MM');


% CALL 'TIME' VARIABLE 
timevar = ds.data('waveTime');
timeconvert = ds.time('waveTime',timevar); % Convert UNIX timestamps to Matlab serial units

% FIND INDEX NUMBERS OF CLOSEST TIME ARRAY VALUES TO START/END DATES
diffstart = abs(timeconvert-startser); % Compute the difference between the 'startdate' serial number and every serial number in 'timevar', to determine which difference is smallest (which number most closely matches 'startser')
[startidx startidx] = min(diffstart); % index of closest object to 'startser', based on the minimum value in 'diffstart' differences array
% closeststart = timeconvert(startidx);  % value of closest object to 'startser' (optional - not used in code)

diffend = abs(timeconvert-endser);
[endidx endidx] = min(diffend); % index of closest object to 'endser', based on the minimum value in 'diffstart' differences array
% closestend = timeconvert(endidx); % value of closest object to 'endser' (optional - not used in code)


% ASSIGN NAMES TO VARIABLES TO BE USED IN COMPENDIUM PLOT, AND SPECIFY TIMERANGE
timecut = ds.time('waveTime',ds.data('waveTime',startidx,endidx)); % Create a subset 'time' variable from startdate to enddate, and convert to Matlab serial units
timedt = datetime(timecut,'ConvertFrom','datenum'); % Convert Matlab serial units 'timecut' section to datetime objects

Hs = ds.data('waveHs',startidx,endidx);
Tp = ds.data('waveTp',startidx,endidx);
Dp = ds.data('waveDp',startidx,endidx);



%% SET UP FIGURE FOR PLOTTING
   
% MAKE SUBPLOTS AND ADD DATA

xvals = timedt; % define x-axis variable (time)
yvals1 = Hs; % define y-variable for Subplot 1
yvals2 = Tp; % define y-variable for Subplot 2
yvals3 = Dp; % define y-variable for Subplot 3

figure % Create new figure

% FIRST SUBPLOT
subplot(3,1,1)
plot(xvals,yvals1,'DatetimeTickFormat','dd')
ylim([0 8])
grid on
grid minor
title(monthtitle)
ylabel('Hs, m')

% SECOND SUBPLOT
subplot(3,1,2)
plot(xvals,yvals2,'DatetimeTickFormat','dd')
ylim([0 25])
grid on
grid minor
ylabel('Tp, s')

% THIRD SUBPLOT
subplot(3,1,3)
plot(xvals,yvals3,'o','DatetimeTickFormat','dd',... % Use 'plot' function (to be able to use datetime array) but set markers as scatter points
    'MarkerSize',2,...
    'MarkerEdgeColor','b',...
	'MarkerFaceColor','b')
ylim([0 360])
grid on
ylabel('Dp, deg')
xlabel('Day')

suptitle(buoytitle) % Set main plot title