function [detect_bound, centerxy, false_boundary, boundary, fiber_ori, fiber_orientation] = opening_method(processed_image,SE_L,DT,CT,TT)
    PI = processed_image;
    deg = -89:DT:90;
    LEN = SE_L;
    detect_bound = {};
    centerxy = [];
    fiber_ori = [];
    for i = 1:length(deg)
        SE = strel('line', LEN,deg(i));
        seDI = strel('disk',1);
        SE2 = imdilate(SE.getnhood,seDI);
        SE3 = strel(SE2);
        afteropening = imopen(PI,SE3);

        detect_bound = [detect_bound;bwboundaries(afteropening,'noholes')];
        stats = regionprops(afteropening,'Orientation','Centroid');
        centroid = cat(1,stats.Centroid);
        centerxy = [centerxy;centroid];
        ori = cat(1,stats.Orientation);
        fiber_ori = [fiber_ori;ori];
    end
    
    false_boundary = [];
    for i = 1:length(fiber_ori)
        ref_center = centerxy(i,:);
        for j = i+1:length(fiber_ori)
            angl_diff = abs(abs(fiber_ori(j))-abs(fiber_ori(i)));
            if angl_diff <= TT
                center_distance = sqrt((ref_center(1)-centerxy(j,1))^2 + (ref_center(2)-centerxy(j,2))^2);
                if center_distance <= CT
                    commonpoint = intersect(detect_bound{i,1},detect_bound{j,1},'rows');
                    if commonpoint ~= 0
                        false_boundary = [false_boundary; j fiber_ori(j)];
                    end
                end
            end
        end
    end
    
    fiber_ori_indx = (1:length(fiber_ori))';
%     false_boundary_u = unique(false_boundary,'rows','stable');
    fiber_indx = setdiff(fiber_ori_indx,false_boundary(:,1),'stable');
    
    fiber_orientation = [];
    boundary = {};
    for ii = 1:length(fiber_indx)
        fiber_orientation = [fiber_orientation;fiber_ori(fiber_indx(ii))];
        boundary = [boundary;detect_bound{fiber_indx(ii),1}];
    end

end