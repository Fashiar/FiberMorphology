function angle=orientation(type)
% orientation of the simualtion nanofibers
switch type
    case 1
        angle = unifrnd(-pi/2,pi/2); %uniform
    case 2
        angle = normrnd(50,15)*pi/180;
    case 3
        angle = normrnd(-30,5)*pi/180;
    case 4
        p=unifrnd(0,1);
        if p<0.5
            angle=normrnd(30,10)*pi/180;
        else
            angle=normrnd(-60,10)*pi/180;
        end
    case 5
        p=unifrnd(0,1);
        if  p<=1/3
            angle=-70*pi/180;
        elseif p>1/3 && p<= 2/3
            angle=-30*pi/180;
        else
            angle=30*pi/180;
        end
%     if p<1/6
%         angle =-70*pi/180;
%     elseif p>=1/6 && p<1/3
%         angle =-30*pi/180;
%     elseif p>=1/3 && p<1/2
%         angle=  10*pi/180;
%     elseif p>=1/2 && p<2/3
%         angle= 45*pi/180;
%     elseif p>=2/3 && p<0.8
%         angle=60*pi/180;
%     else
%         angle=64*pi/180;
%     end
    case 6
        angle = pi/4;
end
end