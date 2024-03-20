%%
load tracks_12k5_coarse_250frames.mat

frame_start = 1;
frame_end = 6;
judge = tracks(:,4) >= frame_start & tracks(:,4) <= frame_end;
tracks = tracks(judge,:);

%%
tracks_find = table2array(...
    readtable('../../results/test_STB/Tracer_0/LongTrackActive_4.csv'));

%%
tracks_id = unique(tracks(:,5));
tracks_find_id = unique(tracks_find(:,1));

%%
figure
for i = 1:size(tracks_id)
    judge = tracks(:,5) == i;
    track = tracks(judge, 1:3);
    plot3(track(:,1), track(:,2), track(:,3), 'k.-', MarkerSize=0.1)
    hold on
end

%%
for i = 1:size(tracks_find_id)
    judge = tracks(:,5) == i;
    track = tracks(judge, 1:3);
    plot3(track(:,1), track(:,2), track(:,3), 'k.-', MarkerSize=0.1)
    hold on
end