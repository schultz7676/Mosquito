basepath1 = 'H:\My Documents\MATLAB\frames\Coachella Valley Intro\footage';
basepath2 = 'H:\My Documents\MATLAB\frames';
basename1 = 'Coachella Valley Intro_';
basename2 = 'Coachella_Valley_Count_snapshot_';
mons = {'Jan.','Feb.','Mar.','Apr.','May ','Jun.','Jul.','Aug.','Sep.','Oct.','Nov.','Dec.'};
v = VideoWriter('CXTCounts.avi');
v.FrameRate = 30;
open(v);
fig = figure;

for k = 0:150
   tail1 = sprintf('%03d',k);
   filename = fullfile(basepath1,[basename1 tail1 '.jpeg']);
   framepic = imread(filename);
   
   % Write the video frame
   image(framepic);
   truesize(fig,[720 1280]);
   %axis square

   frame = getframe;
   writeVideo(v,frame);
end

for ii = 1:30
   writeVideo(v,frame);
end

year = 1994;
for week = 16:2:42
    tail2 = sprintf('%04d%02d',year,week);
    filename = fullfile(basepath2,[basename2 tail2 '.jpeg']);
    framepic = imread(filename);
    image(framepic);
    truesize(fig,[720 1280]);
    %axis square
    
    frame = getframe;
    for ii = 1:30
       writeVideo(v,frame);
    end
end

year = 1995;
for week = 9:2:45
    tail2 = sprintf('%04d%02d',year,week);
    filename = fullfile(basepath2,[basename2 tail2 '.jpeg']);
    framepic = imread(filename);
    image(framepic);
    truesize(fig,[720 1280]);
    %axis square
    
    frame = getframe;
    for ii = 1:30
       writeVideo(v,frame);
    end
end
close(v);

%txt = [mons{mon},' ',num2str(year)];
%text(750,50,txt,'FontSize',14);