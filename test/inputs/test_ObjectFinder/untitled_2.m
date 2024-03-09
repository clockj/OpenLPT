%%
pt_list_all = cell(4,1);
pt_list_all_res = cell(4,1);

for i = 1:4
    pt_list = readtable(['../../solutions/test_ObjectFinder/pt2d_list_cam',num2str(i),'.csv']);
    pt_list = table2array(pt_list);
    pt_list_all{i} = pt_list;
    
    pt_list_res = readtable(['../../results/test_ObjectFinder/pt2d_list_cam',num2str(i),'.csv']);
    pt_list_res = table2array(pt_list_res);
    pt_list_all_res{i} = pt_list_res;
end

%%
plot(pt_list_all{1}(:,1), pt_list_all{1}(:,2), 'k.')

hold on

plot(pt_list_all_res{1}(:,1), pt_list_all_res{1}(:,2), 'r.')