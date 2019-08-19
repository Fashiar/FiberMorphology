function nlines = apply_hough(Image, boundary,rho,theta,thresh,p,scale,gap,mlength)
    true_fiber = boundary;
    whole = zeros(size(Image));
    nlines = [];
    for i = 1 : length(true_fiber)
        fiber = zeros(size(Image));

        for j = 1: length(true_fiber{i})
           x(j)= true_fiber{i}(j,1);
           y(j) = true_fiber{i}(j,2);
           fiber(x(j),y(j)) = 1;
        end
        [H, Theta,Rho] = hough(fiber,'RhoResolution',rho,'ThetaResolution',theta);

        %set the NHoodSize parameter
        hood = floor(size(H)/scale);
        sizeN = hood - (1 - mod(hood,2));
        peaks = houghpeaks(H,p,'Threshold',thresh*max(H(:)),'NHoodSize',sizeN);
        lines = houghlines(fiber,Theta,Rho,peaks,'FillGap',gap,'MinLength',mlength);
        whole = whole + fiber;
        p1 = cat(1,lines.point1);
        p2 = cat(1,lines.point2);
        nlines = [nlines;p1 p2];
    end
    
end