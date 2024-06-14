function correct_ratio = correctRatio(tracks, correctness)

num_tracks = size(unique(tracks(:,1)),1);
correct_ratio = zeros(num_tracks,1);

for i = 1:num_tracks
    judge = tracks(:,1) == i-1;
    len = sum(judge);
    correct_ratio(i) = correctness(i) / len;
end

end