%%
% reorganize tracks_12k5_... format
% orig format: x,y,z,frame,trackID
% new format: trackID(from 0),frame(from 0),x,y,z

load tracks_12k5_coarse_250frames.mat
tracks_orig(:,[1,2,3,4,5]) = tracks(:,[5,4,1,2,3]);
tracks_orig = sortrows(tracks_orig,[1,2]);
tracks_orig(:,[1,2]) = tracks_orig(:,[1,2])-1;
% save("tracks_orig.mat","tracks_orig","-mat")
%%
clear;
load tracks_orig.mat

frame_start = 0;
frame_end = 49;

judge = tracks_orig(:,2) >= frame_start & tracks_orig(:,2) <= frame_end;
tracks_orig = tracks_orig(judge,:);

%%
tracks_code = table2array(...
    readtable('../../results/test_STB/Tracer_0/LongTrackActive_49.csv'));

% tracks_code_inactive = table2array(...
%     readtable('../../results/test_STB/Tracer_0/LongTrackInactive_49.csv'));
% tracks_code_inactive(:,1) = tracks_code_inactive(:,1) + ...
%                             size(unique(tracks_code(:,1)),1);
% tracks_code = [tracks_code; tracks_code_inactive];
%%
show_image = 0;
save_path = '../../results/test_STB/Tracer_0/Track_Quality.mat';

[fit_ratio, fit, correctness] = compareTracks(tracks_orig, tracks_code, show_image, save_path);

correct_ratio = correctRatio(tracks_code, correctness);
