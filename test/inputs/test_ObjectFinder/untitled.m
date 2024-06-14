%%
img = imread('test_function_1.tif');
pt_list = readtable('../../results/test_ObjectFinder/test_function_1.csv');
pt_list = table2array(pt_list);
%%
imshow(img)
hold on
plot(pt_list(:,1)+1, pt_list(:,2)+1,'r.')

%%
pt_list_matlab = Get2DPosOnImage(img);

file = fopen('../../solutions/test_ObjectFinder/test_function_1.csv', 'w');
fprintf(file, '%.8e,%.8e\n', pt_list_matlab');
fclose(file);
%%
function position2D = Get2DPosOnImage(img) 
[rows, cols] = size(img);
img = double(img);
threshold = 30;
position2D = [];
for i = 2 : rows - 1
    for j = 2 : cols - 1
        if (img(i, j) >= threshold && IsLocalMax(img, i, j))
            x1 = j - 1; x2 = j; x3 = j + 1;
            y1 = i - 1; y2 = i; y3 = i + 1;
            x1 = x1 - 1; x2 = x2 - 1; x3 = x3 - 1;
            y1 = y1 - 1; y2 = y2 - 1; y3 = y3 - 1;
            ln1 = NoInfLog(img(i, j - 1));
            ln2 = NoInfLog(img(i, j));
            ln3 = NoInfLog(img(i, j + 1));
        
            xc = -.5 * (ln1 * (x2 ^ 2 - x3 ^ 2) - ln2 * (x1 ^ 2 - x3 ^ 2) + ln3 * (x1 ^ 2 - x2 ^ 2)) / ...
            (ln1 * (x3 - x2) - ln3 * (x1 - x2) + ln2 * (x1 - x3));
            
            ln1 = NoInfLog(img(i - 1, j));
            ln2 = NoInfLog(img(i, j));
            ln3 = NoInfLog(img(i + 1, j));
            yc = -.5 * (ln1 * (y2 ^ 2 - y3 ^ 2) - ln2 * (y1 ^ 2 - y3 ^ 2) + ln3 * (y1 ^ 2 - y2 ^ 2)) / ...
            (ln1 * (y3 - y2) - ln3 * (y1 - y2) + ln2 * (y1 - y3));
            
            if (~isinf(xc) && ~isinf(yc))
                position2D = [position2D; xc, yc];
            end
        end
    end
end

end

function result = IsLocalMax(img, i, j)
   result = 1;
   if (img(i - 1, j) > img(i, j) || img(i + 1, j) > img(i, j) || ...
           img(i, j - 1) > img(i, j) || img(i, j + 1) > img(i, j))
       result = 0;
   end
end

function y = NoInfLog(x)
if x < 1e-6
    x = 1e-6;
end
y = log(double(x));
end
