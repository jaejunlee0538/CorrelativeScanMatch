addpath('tools');
more off
% close all
clear all
% close all
%% intel dataset
load('../data/laser.mat')

% 100/101 ¼ö·Å ¾ÈÇÔ.
% 102/103 ¼ö·Å ¾ÈÇÔ
idx1 = 125;
idx2 = 128;
%%
xy1 = robotlaser_as_cartesian(laser(1, idx1));
xy2 = robotlaser_as_cartesian(laser(1, idx2));
dpos_true = t2v(inv(v2t(laser(1,idx1).pose))*v2t(laser(1,idx2).pose));
theta = dpos_true(3);

%%
cellSize = 0.03;
[lookUpTable, Xmin, Ymin, Xmax, Ymax] = computeModel(xy1(1:2,:), 0.3 ,cellSize);
figure(2)
imagesc(lookUpTable')
set(gca,'YDir','normal');
XTick = get(gca, 'XTick');
YTick = get(gca, 'YTick');
XTick = pixelsToPoints([XTick;XTick], cellSize, Xmin, Ymin, Xmax, Ymax);
YTick = pixelsToPoints([YTick;YTick], cellSize, Xmin, Ymin, Xmax, Ymax);

set(gca, 'XTickLabel', XTick(1,:));
set(gca, 'YTickLabel', YTick(2,:));
%%
Theta0 = theta;
X0 = 0;
Y0 = 0;

%%
figure(1); hold off;
plot(xy1(1,:),xy1(2,:),'r.');
hold on;
plot(xy2(1,:), xy2(2,:),'b.');
legend('scan1','scan2')
grid on; axis equal
%%
gridTheta = [-0.34:0.075:0.34] + Theta0;
gridX = [0,-2.0:cellSize:2.0] + X0;
gridY = [0,-2.0:cellSize:2.0] + Y0;

scores = [];
best_score = -10000;
best_align = [0 0];
for itheta = 1:length(gridTheta)
    probMap = zeros(length(gridX), length(gridY));
    xy2_r = v2t([0 0 gridTheta(itheta)]) * xy2;
    for ix = 1:length(gridX)
        for iy = 1:length(gridY)
            tmp = xy2_r(1:2,:) + repmat([gridX(ix);gridY(iy)], 1, size(xy2_r,2));
    %         figure(1);plot(tmp(1,:),tmp(2,:),'m.'); drawnow; return;
            lookup = pointsToPixels(tmp, cellSize, Xmin, Ymin, Xmax, Ymax);
            lookup(:,lookup(1,:)<1 | lookup(2,:)<1 | lookup(1,:)>size(lookUpTable,1)...
                | lookup(2,:)>size(lookUpTable,2)) = [];

            if ~isempty(lookup)
                scanInd = lookup(1,:) + (lookup(2,:) - 1).*size(lookUpTable,1);
                prob = sum(lookUpTable(scanInd));
                if best_score < prob
                    best_score = prob;
                    best_align = [gridX(ix) gridY(iy)];

                    figure(2);
                    imagesc(lookUpTable'); hold on; 
                    set(gca,'YDir','normal'); axis equal;
                    plot(lookup(1,:), lookup(2,:),'ro'); hold off; drawnow;
                end
            else
                prob = 0;
            end

            probMap(ix, iy) = prob;
        end
    end
    figure(3);
    imagesc(probMap');
    set(gca, 'YDir','normal');
end
figure(4);
figure(4); hold off;
plot(xy1(1,:),xy1(2,:),'r.');
hold on;
xy2_aligned = v2t([best_align(1) best_align(2) theta])* xy2;
plot(xy2_aligned(1,:), xy2_aligned(2,:),'m.');
legend('scan1','scan2 aligned')
grid on; axis equal


