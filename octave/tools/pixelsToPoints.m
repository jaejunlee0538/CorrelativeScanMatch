function pts= pixelsToPoints( pixels, cellSize, Xmin, Ymin, Xmax, Ymax )
% pixels : 2 x N matrix
%     pts = pts - repmat([Xmin; Ymin], 1, size(pts,2));
%     pixeles = floor(pts / cellSize)+ones(size(pts));
    
    pts = (pixels - ones(size(pixels))) * cellSize ...
        + repmat([Xmin; Ymin], 1, size(pixels,2));
end

