function [EdgePoints] = imsegmentation(im_thin, Image,T1,T2)
%This function segment the fibers using DBSCAN Algorithm
%   im_thin is the thinned image, Image is the original SEM image
%   T1, and T2 are two tolerance limit of angle to merge clusters

    EdgePoints = {};
    tol1 = T1;
    tol2 = T2;
    Fibers = bwconncomp(im_thin);
    fib = cell(Fibers.NumObjects,1);
    for i = 1:Fibers.NumObjects
        temp = zeros(length(Fibers.PixelIdxList{i}),2);
        ImageTemp = zeros(size(Image));
        for k = 1:length(Fibers.PixelIdxList{i})
            [temp(k,1),temp(k,2)] = ind2sub(size(Image), Fibers.PixelIdxList{i}(k));
            ImageTemp(temp(k,1),temp(k,2)) = 1;
        end
        fib{i} = temp;

        [y, x] = find(ImageTemp);
        data = [x, y];

        [IDX, noise, D]=DBSCAN(data,1.5,4);

        K = find(IDX == 0);
        NOISE = data(K,:); % Junction points
        [IDXN, noiseN, DN]=DBSCAN(NOISE,1.5,3);

        SEDGES = {}; % Seperated Edges
        A = []; % Array for angels 
        RA = [];
        for j = 1:max(IDXN)
            L = find(IDXN == j);
            EP = NOISE(L,:); % EP for edges pixels (edges are seperated at juntions points)
            if size(EP,1) >= 10 % Skip the small edges
                X = FurthestPixel(EP); % Find the pixels at max distance of a cluster
                A = [A, atand((X(1,2)-X(2,2))/(X(1,1)-X(2,1)))];
                RA = [RA, atand((X(1,1)-X(2,1))/(X(1,2)-X(2,2)))];
                SEDGES = [SEDGES;EP];
            end
        end

        B = A;
        C = SEDGES;

        MEDGES = {}; % Merged Edges

        while isempty(B) == 0    
            [I,dir] = ClosestAngel(B,tol1);
            Merge = []; % Array of edges to merge
                Merge = [Merge;cell2mat(C(I(1)))];
                DI = 1; % Index of cluster to be deleted
                P1 = cell2mat(C(I(1)));
                for k = 2:length(I)
                    P2 = cell2mat(C(I(k)));
                    MP = [P1(round(size(P1,1)/2),:);P2(round(size(P2,1)/2),:)];
                    % MP for Middle Point of two edges
                    tempA = atand((MP(1,2)-MP(2,2))/(MP(1,1)-MP(2,1)));
                    if tempA >= dir-tol2 && tempA <= dir+tol2
                        Merge = [Merge;cell2mat(C(I(k)))];
                        DI = [DI;I(k)];
                    end
                end
            MEDGES = [MEDGES;Merge];
            B(DI) = [];
            C(DI) = [];
        end
        EdgePoints = [EdgePoints;MEDGES];
    end
end

function [I,dir] = ClosestAngel(B,tol)
    dir = B(1);
    I = find(B >= dir-tol & B <= dir+tol);

end

