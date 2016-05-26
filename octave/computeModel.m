function [lookUpTable, Xmin, Ymin, Xmax, Ymax]= computeModel(xy_pts, searchRadius, cellSize)
% xy_pts : xy points. 2 x N matrix. 
% searchRadius : 
% cellSize : 
Xmin = min(xy_pts(1,:)) - (searchRadius*2);
Xmax = max(xy_pts(1,:)) + (searchRadius*2);
Ymin = min(xy_pts(2,:)) - (searchRadius*2);
Ymax = max(xy_pts(2,:)) + (searchRadius*2);

windowSize = ceil(searchRadius / cellSize);

tableNumCellsX = ceil((Xmax - Xmin) / cellSize);
tableNumCellsY = ceil((Ymax - Ymin) / cellSize);

lookUpTable = zeros(tableNumCellsX, tableNumCellsY);

gaussianWindow = zeros(2 * windowSize + 1);
for ix=-windowSize:windowSize
    for iy=-windowSize:windowSize
        dist = cellSize * sqrt(ix^2 + iy^2);
        gaussianWindow(ix+windowSize+1, iy+windowSize+1) = normpdf(dist, 0, 0.1);
    end
end
max_prob =  normpdf(0, 0, 0.1);
% imagesc(gaussianWindow);

xy_indice = pointsToPixels(xy_pts, cellSize, Xmin, Ymin, Xmax, Ymax);
for idx=1:size(xy_pts,2)
    offset= xy_indice(:, idx);
    
    for ix=-windowSize:windowSize
        for iy=-windowSize:windowSize
            p = lookUpTable(ix+offset(1), iy+offset(2));
            q = gaussianWindow(ix+windowSize+1, iy+windowSize+1);
            
%             lookUpTable(ix+offset(1), iy+offset(2)) = p+q;
            lookUpTable(ix+offset(1), iy+offset(2)) = min(p+q, 5*max_prob);
        end
    end
end


end

