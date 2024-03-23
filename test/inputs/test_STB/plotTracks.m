function plotTracks(fig, track, style)
figure(fig);
plot3(track(:,1), track(:,2), track(:,3), style)
end