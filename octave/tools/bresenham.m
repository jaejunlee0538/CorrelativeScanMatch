function [X, Y]=bresenham(p1,p2)
%Matlab optmized version of Bresenham line algorithm. No loops.
%Format:
%               [x y]=bham(x1,y1,x2,y2)
%
%Input:
%               (x1,y1): Start position
%               (x2,y2): End position
%
%Output:
%               x y: the line coordinates from (x1,y1) to (x2,y2)
%
%Usage example:
%               [x y]=bham(1,1, 10,-5);
%               plot(x,y,'or');
% http://www.mathworks.com/matlabcentral/fileexchange/28190-bresenham-optimized-for-matlab
x1 = p1(1); y1 = p1(2);
x2 = p2(1); y2 = p2(2);
x1=round(x1); x2=round(x2);
y1=round(y1); y2=round(y2);
dx=abs(x2-x1);
dy=abs(y2-y1);
steep=abs(dy)>abs(dx);
if steep %dy > dx
    t=dx;
    dx=dy;
    dy=t; 
end

%The main algorithm goes here.
if dy==0 
    q=zeros(dx+1,1);
else
    q=[0;diff(mod([floor(dx/2):-dy:-dy*dx+floor(dx/2)]',dx))>=0];
end

%and ends here.

if steep
    if y1<=y2 
        Y=[y1:y2]'; 
    else
        Y=[y1:-1:y2]'; 
    end
    if x1<=x2 
        X=x1+cumsum(q);
    else
        X=x1-cumsum(q);
    end
else
    if x1<=x2 
        X=[x1:x2]'; 
    else
        X=[x1:-1:x2]'; 
    end
    if y1<=y2 
        Y=y1+cumsum(q);
    else
        Y=y1-cumsum(q); 
    end
end
X = X';
Y = Y';