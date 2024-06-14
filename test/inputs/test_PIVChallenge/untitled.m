%%
cam_id = 0;
frame_start = 0;
frame_end = 30;


figure('Position',[50,50,800,500]);
imh = imshow([],[0,5000]);
for i = frame_start:frame_end
    img = imread(['D_img/D_cam',num2str(cam_id),'_',num2str(i,'%04d'),'.tif']);
    set(imh,'CData',img(200:500,600:800))
    drawnow()
    pause(0.5)
    % for j = 1:30
        % writeVideo(v, uint8(img));
    % end
end


%%
