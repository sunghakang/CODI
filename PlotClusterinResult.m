% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YPML110
% Project Title: Implementation of DBSCAN Clustering in MATLAB
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%

function PlotClusterinResult(X, IDX, d1, d2)

    k=max(IDX);

    Colors=hsv(k);

    Legends = {};

    for i=0:k
        Xi=X(IDX==i,:);
        if i~=0
            Style = 'x';
            MarkerSize = 8;
            Color = Colors(i,:);
%             Legends{end+1} = ['Cluster #' num2str(i)];
        else
            Style = 'o';
            MarkerSize = 6;
            Color = [0 0 0];
%             if ~isempty(Xi)
%                 Legends{end+1} = 'Noise';
%             end
        end
        if ~isempty(Xi)
            plot(Xi(:,d1),Xi(:,d2),Style,'MarkerSize',MarkerSize,'Color',Color);
        end
        hold on;
    end
    hold off;
    title('Projection of Clustering result in high dimensionable space')
%     axis equal;
%     grid on;
%     legend(Legends);
%     legend('Location', 'NorthEastOutside');
xlim([0  255])
ylim([0  255])
end