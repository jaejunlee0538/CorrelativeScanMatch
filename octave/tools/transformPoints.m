function [ pts_transformed ] = transformPoints( pts, se2 )
    % pts : 3 X N matrix
    % se2 : 3 X 3 transformation matrix
    pts = se2 * pts;
    pts_transformed = pts(1:2, :);

end

