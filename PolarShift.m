% this code will draw an arch around a point O' with angle ang_polar and
% radius r_polar that is not the origin of the coordinate system 
% O' is displaced from the origin of the coordinate system by angle
% angle_shift and distance r_shift in the O coordinate system
% ang_num is the num of data points to draw the arch as points
% addPi is to decide whether to add \pi to the phase

function bbbbb=PolarShift(ang_shift, r_shift, ang_polar, r_polar, ang_num, addPi)
aaa=zeros(ang_num,4);
ang_step=ang_polar/(ang_num-1);

for iii=1:1:ang_num
    aaa(iii,1)=(iii-1)*ang_step;
    aaa(iii,2)=r_polar;
    ttemp=aaa(iii,2)*exp(i*aaa(iii,1))*exp(i*ang_shift);  % add the shifting angle
    xtemp=real(ttemp)-r_shift*cos(ang_shift); % add the shifting distance in x
    ytemp=imag(ttemp)-r_shift*sin(ang_shift); % add the shifting distance in y
    if addPi==1
       aaa(iii,3)=atan(ytemp/xtemp)+pi(); %  % add \pi to the phase
    else
       aaa(iii,3)=atan(ytemp/xtemp);
    end
    aaa(iii,4)=hypot(ytemp,xtemp);   % calculate the magnitude or radius of the new vector
%    pp=polarplot(aaa(iii,3),aaa(iii,4),'k.','MarkerSize',2);
end
bbbbb=aaa;
