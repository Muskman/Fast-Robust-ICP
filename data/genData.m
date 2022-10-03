%% Generating source and target clouds with known ground truth

addpath(genpath('/home/lavish/live-projects/scan-registration/'))

clc
clear

pcFileName = 'xyzrgb_statuette';

% create folder
mkdir(pcFileName); mkdir(strcat(pcFileName,'/res'));

% loading the point cloud
targetPtCloud = pcread(strcat(pcFileName,'.ply'));
pcwrite(targetPtCloud, strcat('./',pcFileName,'/target.ply'));

% ground truth transform
% tfGT = 0.5*[1 1 1 pi/2 pi/2 pi/2].*rand(1,6);
tGT = [0.1 0.1 0.2]; 
rGT = [pi/20 pi/6 pi/8];
tfGT = [tGT rGT];
solGT = invertTf(tfGT);

% transformed point cloud 
sourcePtCloud = transformPtCloud(targetPtCloud,tfGT);
pcwrite(sourcePtCloud, strcat('./',pcFileName,'/source.ply'));

% ground truth transfrom
GT = [eul2rotm(solGT(4:6)) solGT(1:3)'; 0 0 0 1]';

fileID = fopen(strcat('./',pcFileName,'/gt.txt'),'w');
nbytes = fprintf(fileID,'%1.5f \t %1.5f \t %1.5f \t %1.5f\n',GT);
fclose(fileID);

disp('GT = '); disp(GT)

% visualizing the point cloud
% figure
% pcshowpair(refPtCloud,currPtCloud); 
% title('Bunny');
% 
% transCurrPtCloud = transformPtCloud(currPtCloud,[solGT(1:3) solGT(4:6)]/0.99);
% figure
% pcshowpair(refPtCloud,transCurrPtCloud);
% title('Transformed Bunny');
