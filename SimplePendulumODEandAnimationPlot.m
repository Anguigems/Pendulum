% animation and ODE solver using RK4 method for Single Pendulum Combined.
% noe that when angle is minimum veloctity is maximum.

domain = [0,30];      % Domain to solve ODE (in seconds), speeds up pendulum or slows it down.
N=500;                % N = number of steps
Y0=[3 0];             % initial theta and theta' 
method =3 ;           % 1 = Euler, 2= Mod Euler, 3 = RK4 
h=(domain(2)-domain(1))/N; 
[l,g,A,omega] = param();   % calling parameters
Y(1,:)=Y0;  t(1)=domain(1);

set(gca,'xlim',[-1 1],'ylim',[-1 1], ...
  'Visible','on','NextPlot','add');
cla;
axis square              
bob=line('color','r','Marker','.','markersize',40,...
    'xdata',[],'ydata',[]);
rod=line('color','b','LineStyle','-','LineWidth',3,...
    'xdata',[],'ydata',[]);

for k=1:N
  t(k+1)=t(k)+h;
  Y(k+1,:)=ODEstep(t(k),Y(k,:),h,method); %k+1th row is step up from kth row.
  ypivot=A*cos(omega*t(k+1)); % the pivots position at time t.   
  xbob=l*sin(Y(k+1,1)); ybob= ypivot-l*cos(Y(k+1,1)); % bob position.
  xrod=[0 xbob]; yrod=[ypivot ybob]; % rod position. 
  set(rod,'xdata',xrod,'ydata',yrod) 
  set(bob,'xdata',xbob,'ydata',ybob)
  drawnow; pause(h)
  figure(2)
  plot(t,Y(:,1),'--b', 'Linewidth', 2.5, 'DisplayName','Angle theta') % plots angle theta
  plot(t,Y(:,2),'--r', 'Linewidth', 2.5, 'DisplayName','Velocity') % plots velocity theta'
    hold on
    legend('Angle theta','Velocity')
    % Cosmetics for graphs
set(gca,'FontSize',16)
%legend('Trapezium','Simpson', 'Location', 'northwest')
ylim([-10,10])
xlim([0,15])
xlabel('time t')
ylabel('angle \theta')
end

function Y=ODEstep(t,Y,h,method)
if method ==3
% RK4 stepper
%Input: current t, current Y (array),  step size h
%Output: approximate solution value at t+h
    hh=h/2;
    s1=dydt(t,Y);
    s2=dydt(t+hh,Y+hh*s1);
    s3=dydt(t+hh,Y+hh*s2);
    s4=dydt(t+h,Y+h*s3);
    Y=Y+h*(s1+2*(s2+s3)+s4)./6;
elseif method == 2
    yhat = Y+h*dydt(t,Y);
    Y=Y+h*(dydt(t,Y)+dydt(t+h,yhat))/2;
elseif method ==1
    Y=Y+h*dydt(t,Y);
end 

end

% Put the RHS of your ODE (system) here
function f=dydt(t,Y)
[l,g,A,omega] = param();
f(1)=Y(2);
f(2)= -(g/l - omega^2*A*cos(omega*t)/l)*sin(Y(1));
% MATLAB 2nd Year Lab 10 explains this ODE.
end

% Put the physical parameters here.
function [l,g,A,omega] = param()
l = 0.5; g= 9.8; A=0; omega=0.1;
% l is the length.
% g is gravity.
% A is the amplitude.
% omega is the frequency of the pivots oscillation, it causes the pendulum 
% to move slightly around it's pivot point. 
% If A is non-zero then the peundulum pivot point will move slightly.
end
