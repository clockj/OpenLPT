function [fit_ratio, fit, correctness] = compareTracks(track0, track1, show_image, save_path)
% track0: the original track
% track1: the track to be compared
% track_id: start from 0

num_track0 = size(unique(track0(:,1)),1);
num_track1 = size(unique(track1(:,1)),1);

frameID_list = unique(track0(:,2));
num_frame = size(frameID_list,1);

tolerant = 3e-1;

fit = zeros(num_track0, 3); % the indicator for fitability for each track,
                            % the first value is the percentage of fitability, 
                            % the second value is the number of tracks to be fitted.
                            % Track length
correctness = zeros(num_track1, 1); % the number of particles which are in the right track


if show_image 
    fig = figure('Visible', 'on'); 
end
for i = 1 : num_track0      % go over every track
    track0_id = i-1;
    judge = track0(:,1) == track0_id;
    track_orig = track0(judge,:);

    fit_number = 0;
    fit_track_number = 0;
    if show_image
        plotTracks(fig, track_orig(:, 3:5), 'b.'); % plot the original track 
        hold on
    end
    j = 1;
    while (j <= num_frame)   % go over every frame
        frame_id = frameID_list(j);
        pt_orig = track_orig(j, 3:5);
        
        judge = track1(:,2) == frame_id;
        track1_search = track1(judge,:);
        index1_search = ...
            track1_search(:, 3) > pt_orig(1) - tolerant & ...
            track1_search(:, 3) < pt_orig(1) + tolerant & ...
            track1_search(:, 4) > pt_orig(2) - tolerant & ...
            track1_search(:, 4) < pt_orig(2) + tolerant & ...
            track1_search(:, 5) > pt_orig(3) - tolerant & ...
            track1_search(:, 5) < pt_orig(3) + tolerant;
        track1_search = track1_search(index1_search, :);

        dist = sqrt(sum((track1_search(:,3:5)-pt_orig).^2, 2));
        [min_dist, id] = min(dist);
        if (min_dist < tolerant)
            track1_id = track1_search(id,1);
            judge = track1(:,1) == track1_id;
            track1_candidate = track1(judge,:);
            
            track1_length = size(track1_candidate,1);
            k = find(track1_candidate(:,2)==frame_id);
            fit_track_number = fit_track_number + 1;

            if show_image 
                plotTracks(fig, track1_candidate(k, 3:5), 'r.'); 
                if k > 1
                    plotTracks(fig, track1_candidate(1:k-1, 3:5), 'k.');
                end
            end

            % search for later frames
            for m = j+1:num_frame
                k = k + 1;

                if k > track1_length
                    m = m-1;

                    break;
                end

                min_dist = sqrt(sum((track1_candidate(k,3:5)-track_orig(m,3:5)).^2,2));
                
                if (min_dist >= tolerant)
                    m = m-1;

                    if show_image 
                        plotTracks(fig, track1_candidate(k:end, 3:5), 'y^'); 
                        pause(0.1)
                    end

                    break;
                end

                fit_number = fit_number + 1;
                correctness(track1_id+1) = correctness(track1_id+1) + 1; % increase the number of correct particle in the track
                
                if show_image 
                    plotTracks(fig, track1_candidate(k, 3:5), 'g*'); 
                end
            end
            j = m;
        end

        j = j + 1;
    end

    
    fit(i, 1) = fit_number;
    fit(i, 2) = fit_track_number;
    fit(i, 3) = num_frame;

    if show_image 
        hold off; 
    end
    if ~(mod(i, 100)) 
        i
        save(save_path, 'fit', 'correctness','-mat');
    end
end
    fit_ratio = sum(fit(:, 1)) / sum(fit(:, 3));

end

