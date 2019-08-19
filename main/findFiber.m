function [nlines,lines_index]=findFiber(bwImage,rho,theta,scale,p,thresh,gap,mlength)


[H,Theta,Rho] = hough(bwImage,'RhoResolution',rho,'Theta',-90:theta:89.5);
hood = floor(size(H)/scale);
sizeN = hood - (1 - mod(hood,2));
peaks = houghpeaks(H,p,'Threshold', thresh*max(H(:)),'NHoodSize',sizeN);
lines = houghlines(bwImage,Theta,Rho,peaks,'FillGap',gap,'MinLength',mlength);
nlines=length(lines);
lines_index=zeros(nlines,4);

for k = 1:length(lines)
   xy = [lines(k).point1; lines(k).point2];
   lines_index(k,:)=[xy(1,1),xy(1,2),xy(2,1),xy(2,2)];
end
end


  