%% Import Count Data
opts = spreadsheetImportOptions("NumVariables", 5);

% Specify sheet and range
opts.Sheet = "9495CO2T";
opts.DataRange = "A2:E1974";

% Specify column names and types
opts.VariableNames = ["TRAP_NUM", "MO", "WK", "YR", "CXT"];
opts.SelectedVariableNames = ["TRAP_NUM", "MO", "WK", "YR", "CXT"];
opts.VariableTypes = ["double", "double", "double", "double", "double"];

% Import the data
TrapCounts = readtable("H:\My Documents\9495co2t.xlsx", opts, "UseExcel", false);
TrapCounts=TrapCounts(~any(ismissing(TrapCounts),2),:);

clear opts

%% Import Trap Sites
opts = delimitedTextImportOptions("NumVariables", 4);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Lon", "Lat", "TRAP_NUM", "description"];
opts.VariableTypes = ["double", "double", "double", "string"];
opts = setvaropts(opts, 4, "WhitespaceRule", "preserve");
opts = setvaropts(opts, 4, "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
TrapSites = readtable("\\stuhome.psu.ds.pdx.edu\j\jocole\My Documents\TrapSites.csv", opts);
TrapSites = removevars(TrapSites,{'description'});

clear opts

%% Import Salton Sea Boundary
opts = delimitedTextImportOptions("NumVariables", 3);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Lon", "Lat", "Count"];
opts.VariableTypes = ["double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
bound = readtable("\\stuhome.psu.ds.pdx.edu\j\jocole\My Documents\MATLAB\saltonseaboundary.csv", opts);

clear opts

%% Setup Interpolation Area
% Coachella Valley Boundaries
minlon = -116.166666667;
minlat =   33.433333333;
maxlon = -115.833333333;
maxlat =   33.566666667;
[XQ,YQ] = meshgrid(linspace(minlon,maxlon,1000),linspace(minlat,maxlat,1000));
R = georefcells([minlat maxlat],[minlon maxlon],[1000 1000],'ColumnsStartFrom','north');

% Create an anonymous 2D interpolation function 'mygriddata'
mygriddata = @(x,y,c){griddata([x;bound.Lon],[y;bound.Lat],[c;bound.Count],XQ,YQ)};

%% Interpolate Groups by Date
% Setup the map coloring edges (uniformly spaced quantiles of the ecdf)
len_cmap = 64;
all_counts = TrapCounts.CXT;
sorted_counts = sort(all_counts);
ii = find(sorted_counts>0);
nonzero_sorted_counts = sorted_counts(ii(1):end);
min_color = nonzero_sorted_counts(1);
max_color = nonzero_sorted_counts(end);
[F,x] = ecdf(nonzero_sorted_counts);
quantiles = discretize(F,len_cmap);
quantile_derivitive = quantiles(2:end)-quantiles(1:end-1);
edges = unique([x(find(quantile_derivitive));max(nonzero_sorted_counts)]); %#ok<FNDSB>

% Setup the alpha gradient
alpha_max = 0.85;
alpha_min = 0.35;
alpha_step = (alpha_max-alpha_min)/(len_cmap-1);
alpha_vec = alpha_min:alpha_step:alpha_max;

% Find the groups of data for each frame
Traps = join(TrapCounts,TrapSites);
Times = findgroups(Traps.YR,Traps.WK);

% Feed lon, lat, and counts to the x,y,c placeholders in 'mygriddata'
AbundanceFrames = splitapply(mygriddata,Traps.Lon,Traps.Lat,Traps.CXT,Times);
AbundanceFrames = cat(3, AbundanceFrames{:});
AbundanceFrames = discretize(AbundanceFrames,edges);
null_counts = isnan(AbundanceFrames);
AbundanceFrames(null_counts) = 0;
AbundanceFrames = uint8(AbundanceFrames./max(AbundanceFrames(:)).*(len_cmap-1));

%% Create Video and Overlays
months = {'Jan.','Feb.','Mar.','Apr.','May ','Jun.','Jul.','Aug.','Sep.','Oct.','Nov.','Dec.'};
v = VideoWriter('abundance.avi');
v.FrameRate = 0.5;
open(v);
% Create a set of frames and write each frame to the file.
figure
cmap = colormap('parula');
TiffTags = struct('ExtraSamples', Tiff.ExtraSamples.AssociatedAlpha,...
                  'Photometric', Tiff.Photometric.RGB, ...
                  'Compression', Tiff.Compression.None);
for k = 1:size(AbundanceFrames,3)
   ind = Times==k;
   selection = Traps.WK(ind);
   week = selection(1);
   selection = Traps.MO(ind);
   month = selection(1);
   selection = Traps.YR(ind);
   year = selection(1);
   
   % Write the video frame
   image(AbundanceFrames(:,:,k));
   axis xy
   txt = [months{month},' ',num2str(year)];
   text(750,50,txt,'FontSize',14);
   frame = getframe;
   writeVideo(v,frame);
   
   % Save the overlay
   filename = ['CXTCounts_',num2str(year),num2str(week,'%02d'),'.tif'];
   Image_RGB = ind2rgb(flipud(AbundanceFrames(:,:,k)),cmap);
   Image_RGB = insertText(Image_RGB,[750,900],txt,'FontSize',36,'TextColor','white');
   alpha_map = alpha_vec(AbundanceFrames(:,:,k)+1);
   alpha_map(null_counts(:,:,k)) = 0;
   alpha_map = flipud(alpha_map);
   Image_RGB_4 = cat(3, Image_RGB, alpha_map);  % M-by-N-by-4 matrix with alpha data
   Image_RGB_4 = uint8(255.*Image_RGB_4);
   geotiffwrite2(filename, Image_RGB_4, R, 'TiffTags', TiffTags);
end
close(v);
