function [num33] = findpb(num,deltatheta, Ori)
%This function find the probability distribution
%   num (90, 85, 80.....-80,-85,-90), deltatheta is the width of bar
%   orientation is the array of angle of each segmented fiber 

    num11 = zeros(180,2);
    num11(:,1) = -89:90;
    for ai = 1:length(Ori)
        temp = round(Ori(ai)*180/pi);
        if temp < -89
            temp = -89;
        end
        temp2 = find(num11(:,1) == temp);
        num11(temp2,2) = num11(temp2,2)+1;
    end

    num22 = zeros(180/deltatheta,1);
    for bi = 1:180
       temp = round(num11(bi,1)/deltatheta)*deltatheta;
       temp2 = find(num == temp);
       num22(temp2,1) = num22(temp2,1) + num11(bi,2);
    end

    numsum = sum(num22(:,1));
    num22(:,1) = num22(:,1)/numsum;

    num33 = zeros(180/deltatheta,1);
    for ci = 1:180/deltatheta
       num33(ci,1) = num22(180/deltatheta+1-ci,1);
    end
end

