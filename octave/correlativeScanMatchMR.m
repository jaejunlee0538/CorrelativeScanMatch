function [ Tv, score ] = correlativeScanMatchMR(scan1, scan2, varargin)
% multi-resolution correlative scan match implementation
    global corrSMMR_verbose;
    global fig_handle;
    global struc_axes;
    
    args = inputParser;
    args.addParameter('verbose', []);
    args.addParameter('initGuess',[0 0 0]);
    args.parse(varargin{:});
    initGuess = args.Results.initGuess;
    if ~isempty(args.Results.verbose)
        fig_handle = args.Results.verbose;
        clf
        struc_axes.axis1 = subplot(2,2,1); %clf;
        struc_axes.axis2 = subplot(2,2,2); %clf;
        struc_axes.axis3 = subplot(2,2,3); %clf;
        struc_axes.axis4 = subplot(2,2,4); %clf;
        corrSMMR_verbose = true;
    else
        corrSMMR_verbose = false;
    end
    
    %%
    if corrSMMR_verbose
        axes(struc_axes.axis1);
        plot(scan1(1,:), scan1(2,:), 'k.'); hold on;
        scan2_init = transformScan(scan2, initGuess);
        plot(scan2_init(1,:), scan2_init(2,:), 'r.'); hold off;
        legend('scan1', 'scan2');
        grid on; axis equal;
    end
    
    %%
    lookupTableH = getLookupTableH(scan1, 0.1, 0.03);
%     figure(1);subplot(1,2,1);imagesc(lookupTableH.lookup');set(gca,'YDir','normal');
    lookupTableL = getLookupTableL(lookupTableH, 5);
%     figure(1);subplot(1,2,2);imagesc(lookupTableL.lookup');set(gca,'YDir','normal');

    theta_range = deg2rad([-10, 10]);
    theta_inc = deg2rad(1);
    [align_L, score_L] = searchBestAlign(lookupTableL, scan2, initGuess,...
     lookupTableL.Xrange, lookupTableL.Yrange,...
     theta_range, theta_inc)
 
    if corrSMMR_verbose
        axes(struc_axes.axis4);
        scan2_t = transformScan(scan2, align_L);
        plot(scan1(1,:), scan1(2,:), 'k.'); hold on;
        plot(scan2_t(1,:), scan2_t(2,:), 'r.'); hold off;
        legend('scan1', 'scan2 aligned LOW');
        grid on;  axis equal; drawnow;
    end
%     Tv = align_L;     score = score_L;
    searchRange = [-lookupTableL.cellSize, lookupTableL.cellSize];
    [Tv, score] = searchBestAlign(lookupTableH, scan2, align_L,...
        searchRange, searchRange,...
        [0 0], deg2rad(0.2));
    
    if corrSMMR_verbose
        axes(struc_axes.axis2);
        scan2_t = transformScan(scan2, Tv);
        plot(scan1(1,:), scan1(2,:), 'k.'); hold on;
        plot(scan2_t(1,:), scan2_t(2,:), 'r.'); hold off;
        legend('scan1', 'scan2 aligned');
        grid on;  axis equal; drawnow;
    end
end

function [best_align, best_score] = searchBestAlign(lookupTable, scan, initGuess, XRange, YRange, ThRange, dTh)
    % XRange : [Xstart, Xend]
    % YRange : [Ystart, Yend]
    % ThRange : [thetaStart, thetaEnd]
    % dTh
    global corrSMMR_verbose;
    global fig_handle;
    global struc_axes;
    
    cellSize = lookupTable.cellSize;
    thetas = [ThRange(1):dTh:ThRange(2)] + initGuess(3)
%     x_search = [XRange(1):cellSize:XRange(2)] + initGuess(1);
%     y_search = [YRange(1):cellSize:YRange(2)] + initGuess(2);
    x_search = floor(XRange(1)/cellSize):1:ceil(XRange(2)/cellSize)
    y_search = floor(YRange(1)/cellSize):1:ceil(YRange(2)/cellSize)
    scoreMap = zeros(size(x_search,2), size(y_search, 2));
    best_score = -10000;
    best_align = [0 0 0];
    for ith = 1:size(thetas, 2)
        scan_r = rotateScan(scan, thetas(ith)) + repmat([initGuess(1);initGuess(2)],1,size(scan,2));
        scan_r_idx = pointsToPixels(scan_r, lookupTable);

        for ix = 1:size(x_search, 2)
            idx2d = scan_r_idx + repmat([x_search(ix);y_search(1)-1],1,size(scan_r_idx,2));
            idx2d(:,idx2d(1,:)<1) = [];
            idx2d(:,idx2d(1,:)>size(lookupTable.lookup, 1)) = [];
            for iy = 1:size(y_search, 2)
                idx2d(2,:) = idx2d(2,:) + 1;
                valid_idx = idx2d;
                valid_idx(:,valid_idx(2,:)<1) = [];
                valid_idx(:,valid_idx(2,:)>size(lookupTable.lookup, 2)) = [];
                if ~isempty(valid_idx)
                    idx1d = valid_idx(1,:) + (valid_idx(2,:) -1).*size(lookupTable.lookup, 1);
                    score = sum(lookupTable.lookup(idx1d));
                    
                    if best_score < score
                        best_score = score;
                        best_align = [cellSize*x_search(ix)+initGuess(1)...
                            cellSize*y_search(iy)+initGuess(2)...
                            thetas(ith)];
                        
                        if corrSMMR_verbose
                            axes(struc_axes.axis3); hold off;
                            imagesc(lookupTable.lookup'); hold on;
                            set(gca, 'YDir', 'normal'); axis equal;
                            plot(idx2d(1,:), idx2d(2,:), 'r.','markersize',20); hold off; drawnow;
                        end
                    end
                else
                    score = 0;
                end
                scoreMap(ix, iy) = score;
            end
        end
    end
end

function [scan_r] = rotateScan(scan, theta)
    c = cos(theta);
    s = sin(theta);
    scan_r = zeros(size(scan));
    for i=1:size(scan,2)
        scan_r(1,i) = c * scan(1,i) - s * scan(2,i);
        scan_r(2,i) = s * scan(1,i) + c * scan(2,i);
    end
end

function [scan_t] = transformScan(scan, v)
    c = cos(v(3));     s = sin(v(3));
    dx = v(1); dy = v(2);
    scan_t = zeros(size(scan));
    for i=1:size(scan,2)
        scan_t(1,i) = c * scan(1,i) - s * scan(2,i) + dx;
        scan_t(2,i) = s * scan(1,i) + c * scan(2,i) + dy;
    end
end

function [lookupTable] = getLookupTableH(scan, kernelRadius, cellSize)
    % scan : xy scan data. 2 x N matrix. 
    % kernelRadius : 
    % cellSize : 
    Xmin = min(scan(1,:)) - (kernelRadius*2);
    Xmax = max(scan(1,:)) + (kernelRadius*2);
    Ymin = min(scan(2,:)) - (kernelRadius*2);
    Ymax = max(scan(2,:)) + (kernelRadius*2);

    windowSize = ceil(kernelRadius / cellSize);

    tableNumCellsX = ceil((Xmax - Xmin) / cellSize);
    tableNumCellsY = ceil((Ymax - Ymin) / cellSize);
    
    %% update Xmax, Ymax
    Xmax = Xmin + tableNumCellsX * cellSize; 
    Ymax = Ymin + tableNumCellsY * cellSize;
    
    %% init lookup table
    lookupTable.cellSize = cellSize;
    lookupTable.Xrange = [Xmin Xmax];
    lookupTable.Yrange = [Ymin Ymax];
    lookupTable.lookup = zeros(tableNumCellsX, tableNumCellsY);
    
    %%
    gaussianWindow = zeros(2 * windowSize + 1);
    for ix=-windowSize:windowSize
        for iy=-windowSize:windowSize
            dist = cellSize * sqrt(ix^2 + iy^2);
            gaussianWindow(ix+windowSize+1, iy+windowSize+1) = normpdf(dist, 0, 0.1);
        end
    end
    max_prob =  normpdf(0, 0, 0.1);
    % imagesc(gaussianWindow);
    
    %%
    xy_indice = pointsToPixels(scan, lookupTable);
    for idx=1:size(scan,2)
        offset= xy_indice(:, idx);

        for ix=-windowSize:windowSize
            for iy=-windowSize:windowSize
                p = lookupTable.lookup(ix+offset(1), iy+offset(2));
                q = gaussianWindow(ix+windowSize+1, iy+windowSize+1);

                %lookUpTable.lookup(ix+offset(1), iy+offset(2)) = p+q;
                lookupTable.lookup(ix+offset(1), iy+offset(2)) = min(p+q, 5*max_prob);
            end
        end
    end
end

function [lookupTableL] = getLookupTableL(lookupTableH, NCells)
    % lookupTableH : high-resolution lookup table to be low-resolutionized
    % NCells : each cell of low-resolution lookup table will be composed of
    %          NCells X NCells of high-resolution cells.
    
    szCellH = size(lookupTableH.lookup);
    cellSizeL = NCells * lookupTableH.cellSize;
    lookupTableL.cellSize = cellSizeL;
    szCellL = ceil(szCellH/NCells);
    lookupTableL.lookup = zeros(szCellL);
    lookupTableL.Xrange = [lookupTableH.Xrange(1) lookupTableH.Xrange(1)+cellSizeL*szCellL(1)];
    lookupTableL.Yrange = [lookupTableH.Yrange(1) lookupTableH.Yrange(1)+cellSizeL*szCellL(2)];
    for ixH = 1:size(lookupTableH.lookup, 1)
        ixL = ceil(ixH / NCells);
        for iyH = 1:size(lookupTableH.lookup, 2)
            iyL = ceil(iyH / NCells);
            if lookupTableL.lookup(ixL, iyL) < lookupTableH.lookup(ixH, iyH)
                lookupTableL.lookup(ixL, iyL) = lookupTableH.lookup(ixH, iyH);
            end
        end
    end
end

function pixeles = pointsToPixels( pts, lookupTable)
% pts : 2xN matrix
    pts = pts - repmat([lookupTable.Xrange(1); lookupTable.Yrange(1)], 1, size(pts,2));
    pixeles = floor(pts / lookupTable.cellSize)+ones(size(pts));
end

function pts= pixelsToPoints( pixels, lookupTable)
% pixels : 2 x N matrix
    
    pts = (pixels - ones(size(pixels))) * lookupTable.cellSize ...
        + repmat([lookupTable.Xrange(1); lookupTable.Yrange(1)], 1, size(pixels,2));
end


