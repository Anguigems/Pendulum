% ODE solver using RK4 method for Single Pendulum
% In my project y=angle theta, x=time t, i.e, labels for axis

domain = [0,30];       % Domain to solve ODE of order k
Narray=[1000];         % All the number of steps to try. If this is an array, the code will plot them all
y0=[0.2 0]; % Initial Conditions:[theta velocity]
[l,g] = param(); % calling parameters
% the initial value of y up to its (k-1)th derivatives- possibly an array
% Specify the RHS of the ODE at the bottom of this code
% Output: t, solution  y

% In this code y will be treated as a matrix
% 1st column = y 
% 2nd column = y'  etc..
% Each row = values at different t.


for i =1:length(Narray)
    N=Narray(i);
    [t,y]= eulerODE(domain,y0,N);
    figure('Name','Simple Pendulum RK4 Solution','NumberTitle','off')
    plot(t,y(:,1),'--b', 'Linewidth', 2.5, 'DisplayName','RK4')  % plots angle theta w.r.t t
%     plot(t,y(:,2),'--r', 'Linewidth', 2.5, 'DisplayName','velocity') % plots velocity theta' 
    hold on
end


% using conditions y(0)=0, y'(0)=2
actual =  @(t) y0(:,1)*cos(t*sqrt(g/l))+y0(:,2)*sqrt(l/g)*sin(t*sqrt(g/l)); 
% actual = @(t) 2*sin(t);
% if you wish to compare with exact solution
T=linspace(0,30);
plot(T, actual(T), 'k', 'LineWidth', 1.5, 'DisplayName','Actual')
legend('RK4','Actual')
hold off


% Cosmetics for graphs
set(gca,'FontSize',16)
%legend('Trapezium','Simpson', 'Location', 'northwest')
ylim([-0.5,0.5])
xlim([0,30])
xlabel('time t (seconds)')
ylabel('angle \theta (radians)')
title('Simple Pendulum RK4 Solution')


function [t,y]= eulerODE(domain,y0,N)
% Euler's Method for Solving IVP dydt = f(t,y)
% y and f are possibly column vectors (i.e. arrays)
% Input: domain, initial value y0 (an array), number of steps N
% Output:  t and y for each step
% Example usage: eulerODE([0 1],1,10);
t(1)=domain(1); y(1,:)=y0;
h=(domain(2)-domain(1))/N;         % this is the step size
for i=1:N
  t(i+1)=t(i)+h;
  y(i+1,:)=ODEstep(t(i),y(i,:),h);
end
end


function y=ODEstep(t,y,h)
% RK4 method
% Input: t,  y (an array),  step size h
% Output: approximate solution value at x+h
s_1=dydt(t,y);
s_2=dydt(t+(h/2),y+(h/2)*s_1);
s_3=dydt(t+(h/2),y+(h/2)*s_2);
s_4=dydt(t+h,y+h*s_3);
y=y+(h/6)*(s_1+2*(s_2+s_3)+s_4);
% The current step is determined by the previous step added to the average of the previous step and the estimation of the current step.  
end


function f = dydt(t,y)
% Put the RHS of your ODE (system) here.
[l,g] = param();
% f and y are arrays
f(1) = y(2);
f(2) = -(g/l)*sin(y(1));
end


% Put the physical parameters here.
function [l,g] = param()
l = 10; g= 9.8;
% l is the length.
% g is gravity.
end
