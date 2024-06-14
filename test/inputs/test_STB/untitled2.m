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

tracks_code_inactive = table2array(...
    readtable('../../results/test_STB/Tracer_0/LongTrackInactive_49.csv'));
tracks_code_inactive(:,1) = tracks_code_inactive(:,1) + ...
                            size(unique(tracks_code(:,1)),1);
tracks_code = [tracks_code; tracks_code_inactive];

tracks_code_exit = table2array(...
    readtable('../../results/test_STB/Tracer_0/ExitTrack_49.csv'));
tracks_code_exit(:,1) = tracks_code_exit(:,1) + ...
                        size(unique(tracks_code(:,1)),1);
tracks_code = [tracks_code; tracks_code_exit];


judge = tracks_code(:,2) >= frame_start & tracks_code(:,2) <= frame_end;
tracks_code = tracks_code(judge, :);


tracks_test = tracks_code;

% tracks_test = removeShortTracks(tracks_code, 7);

%%
show_image = 0;
save_path = '../../results/test_STB/Tracer_0/Track_Quality.mat';

[fit_ratio, fit, correctness] = compareTracks(tracks_orig, tracks_test, show_image, save_path);
correct_ratio = correctRatio(tracks_test, correctness);

n_frag = sum(fit(:,2) > 1);
frag_ratio = n_frag / size(tracks_orig,1);

%%
tracksID_orig_notfind = find(fit(:,1) == 0)-1;

fig = figure(1);
for i = 1:size(tracksID_orig_notfind,1)
    judge = tracks_orig(:,1) == tracksID_orig_notfind(i);
    track = tracks_orig(judge,:);

    plotTracks(fig, track(:,3:5), 'b.-');
    hold on
end
