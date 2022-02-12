
positive      = EdgeMask>0;
numbervector  = [1:n*m];
numbermatrix  = reshape(numbervector,[n,m]);
numbers       = numbermatrix(positive == 1);
plot_new      = double(EdgeMask);



for i = 1:length(IDX)
    plot_new(numbers(i)) = IDX(i);
end


colorcluster = zeros([size(plot_new),3]);
cn = max(max(plot_new))-min(min(plot_new))+2;
cncolor = hsv(cn);
randomindex = randperm(cn);

for i = 1:size(plot_new,1)
    for j = 1:size(plot_new,2)
        clusternm = plot_new(i,j);
        if clusternm == 0
            colorcluster(i,j,:) = [0,0,0];
        else
            colorcluster(i,j,:) = cncolor( randomindex(clusternm+1),:);
        end
    end
end
colorcluster = colorcluster(4:end-3,4:end-3,:);

figure;  imshow (colorcluster,'Border', 'tight');title('Clustering result without ground truth')

