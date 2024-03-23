%%
load tracks_12k5_coarse_250frames.mat

frame_start = 1;
frame_end = 50;
judge = tracks(:,4) >= frame_start & tracks(:,4) <= frame_end;
tracks = tracks(judge,:);
tracks_id = unique(tracks(:,5));

%%
tracks_find = table2array(...
    readtable('../../results/test_STB/Tracer_0/LongTrackActive_49.csv'));
judge = tracks_find(:,2) >= frame_start-1 & tracks_find(:,2) <= frame_end-1;
tracks_find = tracks_find(judge,:);
tracks_find_id = unique(tracks_find(:,1));

% %%
% figure
% for i = 1:size(tracks_id)
%     judge = tracks(:,5) == i;
%     track = tracks(judge, 1:3);
%     plot3(track(:,1), track(:,2), track(:,3), 'k.-', MarkerSize=0.1)
%     hold on
% end
% 
% %%
% for i = 1:size(tracks_find_id)
% % for i = 0:10
%     judge = tracks_find(:,1) == tracks_find_id(i);
%     track = tracks_find(judge, 3:5);
%     plot3(track(:,1), track(:,2), track(:,3), 'r.-', MarkerSize=0.3)
%     hold on
% end

%%
tor = 0.3;
npts = size(unique(tracks(:,5)),1);
n_track_find = size(unique(tracks_find(:,1)),1);

is_correct = zeros(n_track_find,1);
is_long = zeros(n_track_find,1);

for track_find_id = 0:n_track_find-1
    judge = tracks_find(:,1) == track_find_id;
    track = tracks_find(judge,:);
    len = size(track,1);
    match_id = 0;
    
    % if len < 6
    %     continue;
    % end
    
    is_long(track_find_id+1) = 1;
    is_correct(track_find_id+1) = 1;
    err = 0;
    for i = 1:len
        j = track(i,2)+1; % frame_id: start from 0; j: start from 1
        pts = tracks((j-1)*npts+1:j*npts, 1:3);

        % find closest pt
        dist = sqrt(sum((pts-track(i,3:5)).^2, 2));
        [val, id] = min(dist);
        err = err + val;

        if match_id == 0
            match_id = id;
        elseif match_id ~= id
            is_correct(track_find_id+1) = 0;
            break;
        end
    end

    if is_correct(track_find_id+1) > 0
        if err / len > tor
            is_correct(track_find_id+1) = 0;
        end
    end
end


n_correct = sum(is_correct);
n_long = sum(is_long);