fid = fopen('P2093.txt');
data = textscan(fid, '%s%s%s%s%s%s%s %s %q %s');
fclose(fid);
vtype = data{1,9};
lv = zeros(length(vtype),1);
for i = 1:length(vtype)
    if sum(vtype{i} == 'large-vehicle')==13
        lv(i,1) = 1;
    end
end

pos = find(lv==0);

x1 = data{1,1};
x1 = str2double(x1);
y1 = data{1,2};
y1 = str2double(y1);
x2 = data{1,3};
x2 = str2double(x2);
y2 = data{1,4};
y2 = str2double(y2);
x3 = data{1,5};
x3 = str2double(x3);
y3 = data{1,6};
y3 = str2double(y3);
x4 = data{1,7};
x4 = str2double(x4);
y4 = data{1,8};
y4 = str2double(y4);
X = [x1,y1, x2,y2,x3,y3,x4,y4];
X(pos,:) = [];
clear x1 x2 x3 y1 y2 y3 y4
gt = [(X(:,1)+X(:,3)+X(:,5)+X(:,7))./4, (X(:,2)+X(:,4)+X(:,6)+X(:,8))./4];

figure;  imshow (colorcluster,'Border', 'loose'); 
title('Clustering result with ground truth')
hold on
[row,col] = size(colorcluster);
gt(:,2) = gt(:,2) - 586;
gt = gt.*DownSampleRate;
for j = 1: length(gt)
    hold on
    r1 = gt(j,1);
    c1 = gt(j,2);
    plot(r1,c1 ,'wo', 'MarkerFaceColor','red',...
        'MarkerSize', 5, 'LineWidth', 1);
end