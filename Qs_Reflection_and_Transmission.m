%********
% Step1 import and display the original data, display the new data after
% electric delay correction
%********
close all;
load("R8As6p5726GHzSpan10MHz.mat");
%s12=a1st; % S21 is the input file, 1st column is real and 2nd column is imaginary.
%s12=s5p18_20
%fontsize_reso=10;
fontsize_reso=14;
%%%%%%%%%%%%%%%%%%
p_index=10;
s12=R8As6p5726GHzSpan10MHz(:,2*(p_index-1)+1:2*p_index);
w_center=6.5726E9;
w_span=10E6;
c_delay=50*10^-9;
fo_fitSpan=100;
%%%%%%%%%%%%%
rlim_upper1=0.02;
rlim_upper2=0.01;
% s12 has real and imaginary data as 1st and 2nd column
[m_row,n_column]=size(s12);
t12=zeros(m_row,8);
%c_delay=60*10^-9;      % in nano seconds and sets the electric delay due
%to the cable length. This value could easily be retrieved through the
%slope of phase of the s12 data. The right value will make the phase
%outside of the resonance flat instead of a downhill slope. 

pcolor=[0, 0.4470, 0.7410];
pltcolor_arr=[0 0.447 0.741; 0.85 0.325 0.098; 0.929 0.694 0.125];
c_xy=zeros(m_row,2);
W_all=zeros(m_row,1);
w_points=m_row;
w_step=w_span/(w_points-1);
W_all=w_center-w_span/2:w_step:w_center+w_span/2;
W_all_unit=W_all/w_center;
for indexer=1:1:m_row
    t12(indexer,1)=W_all(indexer);      % frequency array
    temp1=s12(indexer,1)+i*s12(indexer,2); % convert to complex number
    t12(indexer,3)=abs(temp1);          % get the amplitude
    t12(indexer,4)=angle(temp1);        % get the phase
    temp2=temp1*exp(2*pi*i*t12(indexer,1)*c_delay); % correct electic delay
    t12(indexer,6)=abs(temp2);          % amplitude after delay correction
    t12(indexer,7)=angle(temp2);        % phase after delay correction
end
[minval minindex]=min(t12(:,3));
%**********
% Polar plot of the original data
%**********
fig1=1;
if fig1==1
    figure(1);
%    set(gcf, 'PaperUnits', 'centimeters', 'PaperSize', [4.2 4.2], 'PaperPositionMode', 'manual', 'PaperPosition', [0 0 4.2 4.2]);
    ax = polaraxes;
    pp=polarplot(t12(:,4),t12(:,3),'+','MarkerSize',0.5,'MarkerEdgeColor',[0 0.447 0.741]);  % polar plot of the original data
    pp.Color=pcolor;
    ax.FontSize = fontsize_reso;
    rlim([0 rlim_upper1]);
end
%**********
% step 1 removing electric delay and display the corrected data
% coordinate system
%**********
fig2=1;
if fig2==1
    figure(2);
 %   set(gcf, 'PaperUnits', 'centimeters', 'PaperSize', [4.2 4.2], 'PaperPositionMode', 'manual', 'PaperPosition', [0 0 4.2 4.2]);
    ax = polaraxes;
    pp=polarplot(t12(:,7),t12(:,6),'+', 'MarkerSize', 0.5,'MarkerEdgeColor',[0 0.447 0.741]);  % polar plot of the data after delay correction
    pp.Color=pcolor;
    ax.FontSize = fontsize_reso;
    rlim([0 rlim_upper1]);
    print(figure(2),'polar_afterEdelay.png','-dpng','-r400');
end
%**********
% step 2 fit circle to the data and shift the data to the origin of the
% coordinate system
%**********
fig3=1;
if fig3==1
    figure(3);
    pp=polar(t12(:,7),t12(:,6),'c*');  % polar plot of the data after delay correction
    pp.Color=pltcolor_arr(1,:);
    hold on;
end
c_xy_after_delay=zeros(m_row,2); % array for circle fitting
for jjj=1:1:m_row
    c_xy_after_delay(jjj,1)=t12(jjj,6)*cos(t12(jjj,7));  % real part of the new data
    c_xy_after_delay(jjj,2)=t12(jjj,6)*sin(t12(jjj,7));  % imaginary part of the new data
end
c_num=CircleFitByTaubin(c_xy_after_delay(minindex-fo_fitSpan:minindex+fo_fitSpan,:));     %run circle fit here;
c_x=c_num(1);    % x coordinate of the fitted circle center
c_y=c_num(2);    % y coordinate of the fitted circle center
c_r=c_num(3);    % radius of the fitted circle
%%% generating a circle based on the parameters from the fitting circle 
ang=0:0.01:2*pi;     % the angle
xp=c_r*cos(ang);     % x
yp=c_r*sin(ang);     % y
xp_n=xp+c_x;         % shift the x origin of the x & y
yp_n=yp+c_y;         % shift the y origin of the x & y
[c_m,c_n]=size(ang); % get size of the cicle array
c_fit=zeros(c_n,3);
c_fit_shift=zeros(c_n,3);
for kk=1:1:c_n
    c_fit(kk,1)=xp_n(kk)+i*yp_n(kk);   % convert to complex number
    c_fit(kk,2)=abs(c_fit(kk,1));      % get amplitude
    c_fit(kk,3)=angle(c_fit(kk,1));    % get angle
    c_fit_shift(kk,1)=xp(kk)+i*yp(kk);   % convert to complex number
    c_fit_shift(kk,2)=abs(c_fit_shift(kk,1));      % get amplitude
    c_fit_shift(kk,3)=angle(c_fit_shift(kk,1));    % get angle  
end
if fig3==1
    ax(1)=polar(c_fit(:,3),c_fit(:,2),'-r');      % plot the fitted circle on the new data
end
Angle_center=angle(c_x+i*c_y);  % get the angle of the center of the circle
Angle_center_degree=Angle_center*180/pi();
ang_alpha=Angle_center;
rrr=sqrt(c_x^2+c_y^2);
theta_m=[0 ang_alpha];
    r_m=[0 rrr];
    [x y]=pol2cart(theta_m,r_m);
c_xy_after_delay_center=zeros(m_row,2);
c_xy_after_delay_center(:,1)=c_xy_after_delay(:,1)-c_x;  % shift the x to the orinin
c_xy_after_delay_center(:,2)=c_xy_after_delay(:,2)-c_y;  % shift the y to the origin
c_xy_after_delay_center_polar=zeros(2,m_row);
for lll=1:1:m_row
    c_xy_after_delay_center_polar(lll,1)=angle(c_xy_after_delay_center(lll,1)+i*c_xy_after_delay_center(lll,2)); % angle of shifted circle 
    c_xy_after_delay_center_polar(lll,2)=abs(c_xy_after_delay_center(lll,1)+i*c_xy_after_delay_center(lll,2));   % amplitude of the shifted circle
end
if fig3==1
    pp=polar(c_xy_after_delay_center_polar(:,1),c_xy_after_delay_center_polar(:,2),'+'); % plot the shifted circle. ...
    pp.Color=pltcolor_arr(2,:);
    set(gca,'FontSize',fontsize_reso);
    % The center of the circle should align with the origin of the coordinate system
    hold off;
end
set(findall(gcf,'-property','FontSize'),'FontSize',fontsize_reso);
c_xy_after_delay_center_rotate_polar=zeros(m_row,2);
c_xy_after_delay_center_rotate_polar_expanded=zeros(m_row,1);
for lll=1:1:m_row
    c_xy_after_delay_center_rotate_polar(lll,1)=angle((c_xy_after_delay_center(lll,1)+i*c_xy_after_delay_center(lll,2))*exp(-i*Angle_center));
    c_xy_after_delay_center_rotate_polar(lll,2)=abs((c_xy_after_delay_center(lll,1)+i*c_xy_after_delay_center(lll,2))*exp(-i*Angle_center));
end
%**********
% step 3 rotate the new data to align with the x axis of the
% coordinate system
%**********
fig4=1;
if fig4==1
    figure(4);  %  display the rotated circle 
    set(gcf, 'PaperUnits', 'centimeters', 'PaperSize', [4.2 4.2], 'PaperPositionMode', 'manual', 'PaperPosition', [0 0 4.2 4.2]);
    ax = polaraxes;
    pp=polarplot(c_xy_after_delay_center_rotate_polar(:,1),c_xy_after_delay_center_rotate_polar(:,2),'+','MarkerSize',0.3,'MarkerEdgeColor',[0 0.447 0.741]);
    pp.Color=pcolor;
    hold on
    pp=polarplot(c_fit_shift(:,3),c_fit_shift(:,2),'+','MarkerSize',0.2,'MarkerEdgeColor',[0.85 0.325 0.098]);
    hold off;  
    ax.FontSize = fontsize_reso;
    rlim([0 rlim_upper2]);
    print(figure(4),'polar_shiftToOrigin.png','-dpng','-r400');
end
%**********
% step 4 use linear regression to find W_o, Q_L and theta_o
%********** 
Angle_theta=zeros(m_row,1);
Angle_theta=c_xy_after_delay_center_rotate_polar(:,1);
Angle_theta_expanded=unwrap(Angle_theta);
diffmax=0;
diffmaxflag=0;  % the index of the fre_o
[mmm nnn] = QLinFit(W_all, W_all_unit, Angle_theta_expanded, 1);
fitval = coeffvalues(mmm);
Q_r_reg=fitval(2);
Theta_o=fitval(1);
%**********
% plot the phase and expanded or unwrapped phase and mark the resonance
% point
%**********
fre_o=0.1E9;
diffmaxflag=0;  % below is finding the index of the resonance point
for xxx=1:1:m_row-1
    if abs(W_all(xxx)-fitval(3))  < fre_o
       fre_o=abs(W_all(xxx)-fitval(3));
       diffmaxflag=xxx;
    end
end
% run linear regression again to tune up the fitting parameters
[mmm nnn] = QLinFit(W_all(1,diffmaxflag-fo_fitSpan:diffmaxflag+fo_fitSpan),...
   W_all_unit(1,diffmaxflag-fo_fitSpan:diffmaxflag+fo_fitSpan),...
   Angle_theta_expanded(diffmaxflag-fo_fitSpan:diffmaxflag+fo_fitSpan)',0);  
fitval = coeffvalues(mmm);
Q_r_reg=fitval(2);
Theta_o=fitval(1)-Angle_theta_expanded(diffmaxflag);
Theta_o_abs=fitval(1);
f_o_fit=fitval(3);
fig6=1;
    if fig6==1
    figure(6);
    hold on;
    plot(W_all,Angle_theta);   % display the phase as function of frequency
    scatter(W_all(diffmaxflag),Angle_theta(diffmaxflag),'r');  %
    hold off;
    end
fig7=1;
if fig7==1
    figure(7);
    hold on;
    plot(W_all,Angle_theta_expanded, 'Color',pcolor);   % display the phase as function of frequency
    scatter(W_all(diffmaxflag),Angle_theta_expanded(diffmaxflag));
    ylabel('Phase (unwrapped)');
    xlabel('Frequency (Hz)');
    hold off;
    set(gca,'FontSize',fontsize_reso);
end
%**********
% find the 90 degree points before and after the resonance points
%**********
mmm=Angle_theta_expanded;
angdiff=0.1;
angleft=0;
% find the left of the 45 degree point
for iii=1:1:diffmaxflag
    if abs(abs(mmm(iii)-mmm(diffmaxflag))-pi()/2)<angdiff
       angdiff=abs(abs(mmm(iii)-mmm(diffmaxflag))-pi()/2);
       angleft=iii;
    end
end
angdiff=0.1;
angright=0;
% find the right of the 45 degree point
for iii=diffmaxflag:1:m_row
    if abs(abs(mmm(iii)-mmm(diffmaxflag))-pi()/2)<angdiff
       angdiff=abs(abs(mmm(iii)-mmm(diffmaxflag))-pi()/2);
       angright=iii;  
    end
end
figure(8);  %% plot resonance and +-3dB points.
ax(1)=polar(c_xy_after_delay_center_rotate_polar(:,1),c_xy_after_delay_center_rotate_polar(:,2),'-');
ax(1).Color=pcolor;
hold on;
ax(2) = polar([0 0],[0 0.0025],'g');
ang_alpha=c_xy_after_delay_center_rotate_polar(diffmaxflag+1,1);
rrr=c_xy_after_delay_center_rotate_polar(diffmaxflag+1,2);
theta_m=[0 ang_alpha];
    r_m=[0 rrr];
    [x y]=pol2cart(theta_m,r_m);
set(ax(2),'xData',x,'yData',y);
ax(3) = polar([0 0],[0 0.0025],'r');
ang_alpha=c_xy_after_delay_center_rotate_polar(angleft,1);
rrr=c_xy_after_delay_center_rotate_polar(angleft,2);
theta_m=[0 ang_alpha];
    r_m=[0 rrr];
    [x y]=pol2cart(theta_m,r_m);
set(ax(3),'xData',x,'yData',y);
ax(4) = polar([0 0],[0 0.0025],'r');
ang_alpha=c_xy_after_delay_center_rotate_polar(angright,1);
rrr=c_xy_after_delay_center_rotate_polar(angright,2);
theta_m=[0 ang_alpha];
    r_m=[0 rrr];
    [x y]=pol2cart(theta_m,r_m);
set(ax(4),'xData',x,'yData',y);
hold off;
set(findall(gcf,'-property','FontSize'),'FontSize',fontsize_reso);
%*********************
% plot the phase with the fitted Q_L to check the quality of the regression
% a staight line is expected for best quality 
%*********************
fig9=1;
if fig9==1
   figure(9);
   plot(Angle_theta_expanded,2*atan(2*Q_r_reg*(1-W_all/fitval(3)))+Theta_o_abs,'LineWidth',1.5, 'Color',[0 0.447 0.741]);
   xlabel('Measured Phase (unwrapped)');
   ylabel('atan(2Q_{r}\delta_{L})');
   xlim([-6 0.15]);
   ylim([-6 0.15]);
   hold off;
   set(gca,'FontSize',fontsize_reso);
end
fig12=1;
if fig12==1
   figure(12);
   hold on;
   W_all_less=W_all(1:20:m_row);
   Angle_theta_expanded_less=2*atan(2*Q_r_reg*(1-W_all_less/fitval(3)))+Theta_o_abs;
   plot(W_all_less*1E-9, Angle_theta_expanded_less,'mo',...
    'MarkerEdgeColor',[0 0.447 0.741],...
    'MarkerFaceColor','none',...
    'MarkerSize',6,'LineWidth',1);
   plot(W_all*1E-9, 2*atan(2*Q_r_reg*(1-W_all/fitval(3)))+Theta_o_abs,'LineWidth',1.5,'Color','r');
   legend('measured','fitted');
   legend boxoff;
   ylabel('Phase (unwrapped)');
   xlabel('frequency (GHz)');
   hold off;
   set(gca,'FontSize',fontsize_reso);
end
%*********************
% plot the three angles to check the calculation
%*********************
ang_plot=1;   % set to 1 for plotting the three angles or not. 
if ang_plot==1
    figure(10);
    polarplot(t12(:,7),t12(:,6),'Color',pcolor,'LineWidth',2);
    hold on;
    Angle_center=angle(c_x+i*c_y);  % get the angle of the center of the circle
    Angle_center_degree=Angle_center*180/pi();
    ang_alpha=Angle_center;
    rrr=sqrt(c_x^2+c_y^2);
    pp=polarplot(ang_alpha,rrr,'r*','MarkerSize',5);
    polarplot(t12(diffmaxflag,7),t12(diffmaxflag,6),'r*','MarkerSize',10); % 
    polarplot([ang_alpha;t12(diffmaxflag,7)], [rrr; t12(diffmaxflag,6)],'r','LineWidth',2); % center of circle to resonance point  
    polarplot([ang_alpha;0], [rrr*1.4;0],'r','LineWidth',2); % origin to center of circle
    ang_num=101;
    ang_step=Angle_center/ang_num;
    rrr_ang=rrr/8;
    for iii=1:1:ang_num
        pp=polarplot((iii-1)*ang_step,rrr_ang,'k.','MarkerSize',2);
    end
    text(Angle_center/2,rrr/5,'\alpha','FontSize',18);
    ang_theta=acos((c_r^2+rrr^2-t12(diffmaxflag,6)^2)/(2*c_r*rrr));
    d_dcm=sqrt(c_r^2+rrr^2+2*c_r*rrr*cos(ang_theta));
    ang_dcm=acos((rrr^2+d_dcm^2-c_r^2)/(2*rrr*d_dcm));
    ang_phi=ang_theta-ang_dcm;
    pha_unwrap=unwrap(t12(:,7));
    pha_mean=mean(pha_unwrap);
    n_2pi=round(pha_mean/(2*pi()));
    if ang_alpha<pha_unwrap(diffmaxflag)-n_2pi*2*pi()  
        ang_sign=-1;
    else
        ang_sign=1;
    end
    polarplot([ang_alpha+ang_sign*ang_dcm;0], [d_dcm;0],'r','LineWidth',2); % off resonance point to origin
    polarplot([ang_alpha+ang_sign*ang_dcm;ang_alpha], [d_dcm;rrr],'r','LineWidth',2); % circle center to off resonance point
    amp_max=max(t12(:,6));
    text(ang_alpha+ang_sign*ang_dcm/1.2,d_dcm-amp_max/10,'\phi','FontSize',15);
    text(ang_alpha+ang_sign*ang_dcm/5,rrr+amp_max/10,'\theta_{o}','FontSize',15);
    aaaaa=PolarShift(Angle_center,rrr,-ang_theta,-c_r/5,100,0);
    for iii=1:1:ang_num
        pp=polarplot(aaaaa(:,3),aaaaa(:,4),'k.','MarkerSize',2);  %  arch for theta_o
    end
    aaaaa=PolarShift(ang_alpha+ang_sign*ang_dcm,d_dcm,-(ang_theta-ang_dcm),c_r/3,100,0);
    for iii=1:1:ang_num
        pp=polarplot(aaaaa(:,3),aaaaa(:,4),'k.','MarkerSize',2);  %  arch for phi
    end
end
%********************
% calculate the Qs
%********************
ang_phi=ang_theta-ang_dcm;
% for transmission type
is_tran=1;
if is_tran==1  % transmission type
    fprintf('transmission type');
    %%% 1. with +-90 points:
    Q_lc=W_all(diffmaxflag)/(W_all(angright)-W_all(angleft)); % total Q
    Q_c=((sqrt(c_x^2+c_y^2)+c_r)/(2*c_r))*Q_r_reg;    % coupling or external Q
    Q_i=1/(1/Q_lc-1/Q_c);       % internal Q
    %%% 2. with linear regression and diameter correction method
    c_r_dcm=abs(c_r*cos(ang_phi));
    Q_c_reg_dcm=((sqrt(c_x^2+c_y^2)+c_r_dcm)/(2*c_r_dcm))*Q_r_reg;  % coupling Q
    Q_i_reg_dcm=1/(1/Q_r_reg-1/Q_c_reg_dcm);   % internal Q
    Q_c_reg_dcm_img=Q_c_reg_dcm*abs(tan(ang_phi));   % imaginary part of coupling Q
else  % reflective type 
    fprintf('reflection type');
    %%% 1. with +-90 points:
    Q_lc=W_all(diffmaxflag)/(W_all(angright)-W_all(angleft)); % total Q
    Q_c=((sqrt(c_x^2+c_y^2)+c_r)/c_r)*Q_r_reg;    % coupling or external Q
    Q_i=1/(1/Q_lc-1/Q_c);       % internal Q
    %%% 2. with linear regression and diameter correction method
    c_r_dcm=abs(c_r*cos(ang_phi));
    Q_c_reg_dcm=((sqrt(c_x^2+c_y^2)+c_r_dcm)/c_r_dcm)*Q_r_reg;  % coupling Q
    Q_i_reg_dcm=1/(1/Q_r_reg-1/Q_c_reg_dcm);   % internal Q
    Q_c_reg_dcm_img=Q_c_reg_dcm*abs(tan(ang_phi));   % imaginary part of coupling Q

end