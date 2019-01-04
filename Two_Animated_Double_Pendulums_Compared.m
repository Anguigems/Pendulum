% Animation program for 2 double pendulums set off at different angles in same windows, for YOUTUBE.
domain = [0 15];      % Domain to solve ODE (in seconds), speeds up pendulum or slows it down.
N=1000;                % N = number of steps
Y0=[1.9 0 1.9 0 2 0 2 0];     % initial [theta1,theta1',theta2,theta2',theta3,theta3',theta4,theta4']
method =3 ;           % 1 = Euler, 2= Mod Euler, 3 = RK4 
h=(domain(2)-domain(1))/N; 
[m1,m2,l1,l2,g,u] = param();   % calling parameters
Y(1,:)=Y0; % tells it to set row 1 as y0 vector. 
t(1)=domain(1);



% Cosmetics for plot
set(gca,'xlim',[-1 1],'ylim',[-1 1], ...
  'Visible','on','NextPlot','add');
cla;
axis square              
bob1=line('color','r','Marker','.','markersize',40,...
    'xdata',[],'ydata',[]); 
rod1=line('color','b','LineStyle','-','LineWidth',3,...
    'xdata',[],'ydata',[]);
bob2=line('color','g','Marker','.','markersize',40,...
    'xdata',[],'ydata',[]);
rod2=line('color','k','LineStyle','-','LineWidth',3,...
    'xdata',[],'ydata',[]);  
bob3=line('color','c','Marker','.','markersize',40,...
    'xdata',[],'ydata',[]); 
rod3=line('color','b','LineStyle','-','LineWidth',3,...
    'xdata',[],'ydata',[]);
bob4=line('color','m','Marker','.','markersize',40,...
    'xdata',[],'ydata',[]);
rod4=line('color','k','LineStyle','-','LineWidth',3,...
    'xdata',[],'ydata',[]);
title('Two Double Pendulum Animations on Same Axis')

for k=1:N
  t(k+1)=t(k)+h;
  Y(k+1,:)=ODEstep(t(k),Y(k,:),h,method); %k+1th row is step up from kth row.
  xbob1=l1*sin(Y(k+1,1)); ybob1=-l1*cos(Y(k+1,1)); % bob1 position, i.e m1.
  xrod1=[0 xbob1]; yrod1=[0 ybob1]; % rod1 position. 
  xbob2=l1*sin(Y(k+1,1))+l2*sin(Y(k+1,3)); % bob2 position, i.e m2 1st Pend.
  ybob2=-l1*cos(Y(k+1,1))-l2*cos(Y(k+1,3)); % bob2 position, i.e m2 1st Pend.
  xrod2=[xbob2 xbob1]; yrod2=[ybob2 ybob1]; % rod2 position. 
  set(rod1,'xdata',xrod1,'ydata',yrod1)
  set(bob1,'xdata',xbob1,'ydata',ybob1)
  set(rod2,'xdata',xrod2,'ydata',yrod2)
  set(bob2,'xdata',xbob2,'ydata',ybob2)
  xbob3=l1*sin(Y(k+1,5)); ybob3=-l1*cos(Y(k+1,5)); % bob1 position, i.e m1 2nd Pend.
  xrod3=[0 xbob3]; yrod3=[0 ybob3]; % rod3 position. 
  xbob4=l1*sin(Y(k+1,5))+l2*sin(Y(k+1,7)); % bob4 position, i.e m2 2nd Pend.
  ybob4=-l1*cos(Y(k+1,5))-l2*cos(Y(k+1,7)); % bob4 position, i.e m2 2nd Pend.
  xrod4=[xbob4 xbob3]; yrod4=[ybob4 ybob3]; % rod4 position. 
  set(rod3,'xdata',xrod3,'ydata',yrod3)
  set(bob3,'xdata',xbob3,'ydata',ybob3)
  set(rod4,'xdata',xrod4,'ydata',yrod4)
  set(bob4,'xdata',xbob4,'ydata',ybob4)
  drawnow; pause(h)
end
% F=getframe(gfc);
myVideo = VideoWriter('Two_Double_Animation_Comparison_for_youtube.mp4','MPEG-4');
open(myVideo)
writeVideo(myVideo);  
%writeVideo(myVideo,F);
close(myVideo)

function Y=ODEstep(t,Y,h,method)
if method ==3
% RK4 stepper
%Input: current t, current Y (array),  step size h
%Output: approximate solution value at t+h
    hh=h/2;
    s1=dydx(t,Y);
    s2=dydx(t+hh,Y+hh*s1);
    s3=dydx(t+hh,Y+hh*s2);
    s4=dydx(t+h,Y+h*s3);
    Y=Y+h*(s1+2*(s2+s3)+s4)./6;
elseif method == 2
    yhat = Y+h*dydx(t,Y);
    Y=Y+h*(dydx(t,Y)+dydx(t+h,yhat))/2;
elseif method ==1
    Y=Y+h*dydx(t,Y);
end 

end

% Put the RHS of your ODE (system) here
function f=dydx(t,Y)
[m1,m2,l1,l2,g,u] = param();
% f and y are arrays
f(1) = Y(2);
f(2) = (-g*sin(Y(1))-u*l2*Y(4)*Y(4)*sin(Y(1)-Y(3))+u*g*sin(Y(3))*cos(Y(1)-Y(3))-u*l1*Y(2)*Y(2)*sin(Y(1)-Y(3))*cos(Y(1)-Y(3)))/(l1*(1-u*(cos(Y(1)-Y(3))*cos(Y(1)-Y(3)))));
f(3) = Y(4);
f(4) = (-g*sin(Y(3))+l1*Y(2)*Y(2)*sin(Y(1)-Y(3))+g*sin(Y(1))*cos(Y(1)-Y(3))+u*l2*Y(4)*Y(4)*sin(Y(1)-Y(3))*cos(Y(1)-Y(3)))/(l2*(1-u*(cos(Y(1)-Y(3))*cos(Y(1)-Y(3)))));
f(5) = Y(6);
f(6) = (-g*sin(Y(5))-u*l2*Y(8)*Y(8)*sin(Y(5)-Y(7))+u*g*sin(Y(7))*cos(Y(5)-Y(7))-u*l1*Y(6)*Y(6)*sin(Y(5)-Y(7))*cos(Y(5)-Y(7)))/(l1*(1-u*(cos(Y(5)-Y(7))*cos(Y(5)-Y(7)))));
f(7) = Y(8);
f(8) = (-g*sin(Y(7))+l1*Y(6)*Y(6)*sin(Y(5)-Y(7))+g*sin(Y(5))*cos(Y(5)-Y(7))+u*l2*Y(8)*Y(8)*sin(Y(5)-Y(7))*cos(Y(5)-Y(7)))/(l2*(1-u*(cos(Y(5)-Y(7))*cos(Y(5)-Y(7)))));
end

% Put the physical parameters here.
function [m1,m2,l1,l2,g,u] = param()
m1=1; m2=1; l1=0.45; l2=0.45; g= 9.8; u= m2/(m1+m2);
% l1 is the length of the rod between the pivot and bob 1.
% l2 is the length of the rod between bob 1 and bob 2.
% only need an l1,l2,m1,m2 since the masses and lengths of both pendulums
% will be the same.
% g is gravity.
% omega is the frequency of the pivots oscillation, it causes the pendulum 
% to move slightly around it's pivot point. 
% If A is non-zero then the peundulum pivot point will move slightly.
end