function tracks_new = removeShortTracks(tracks, n_threshold)
    n_track = size(unique(tracks(:,1)), 1);
    
    is_keep = zeros(size(tracks,1), 1);
    for i = 1:n_track
        judge = tracks(:,1) == i-1;
        len = sum(judge);

        if len >= n_threshold
            is_keep = is_keep | judge;
        end
    end

    tracks_temp = tracks(is_keep,:);
    tracks_new = tracks_temp;
    trackID_list = unique(tracks_new(:,1));
    n_track_new = size(trackID_list, 1);

    for i = 1:n_track_new
        judge = tracks_temp(:,1) == trackID_list(i);
        tracks_new(judge,1) = i-1;
    end
end