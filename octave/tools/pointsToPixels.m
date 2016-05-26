function pixeles = pointsToPixels( pts, cellSize, Xmin, Ymin, Xmax, Ymax)
% pts : 2xN matrix
    pts = pts - repmat([Xmin; Ymin], 1, size(pts,2));
    pixeles = floor(pts / cellSize)+ones(size(pts));
end

