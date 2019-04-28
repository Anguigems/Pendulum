% MIT License

% Copyright (c) 2018 Anguigems

% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:

% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.


% ODE solver using RK4 method to approximate the solution to two Double 
% Pendulums in order to plot the 2D and 3D Phase Space (Poncaire) diagrams 
% for each pendulum for direct visual comparison. Also calculates
% corresponding lyapunov exponent for mathematical comparison.

% t=time, y=theta.


domain = [0,11.98];       % Domain to solve ODE of order k
Narray=[1000];         % All the number of steps to try. If this is an array, the code will plot them all
y0=[pi/2, 0, 1.63, 0, pi/2, 0, 1.64, 0]; % Initial Conditions:[theta1 velocity1 theta2 velocity2] [y1 y1' y2 y2']
% manually update initial difference d on line 129 according to imput
% values used.
% both pendulums display different solutions for even 1 degree=0.017
% radians different.
[m1,m2,m3,m4,l1,l2,l3,l4,g,u1,u2] = param(); % calling parameters
% the initial value of y up to its (k-1)th derivatives- possibly an array
% Output: t, solution  y


for i =1:length(Narray)
    N=Narray(i);
    [t,y]= eulerODE(domain,y0,N);
% 2D Phase Space
    figure('Name','Comparing Double Pendulum 2D Phase Space','NumberTitle','off')
    plot(y(:,1),y(:,2),'-b', 'Linewidth',1, 'DisplayName','upper bob of 0004') % upper bob m1 of first (angle,velocity) plot
    hold on
    plot(y(:,3),y(:,4),'-r','Linewidth',1, 'DisplayName','lower bob of 0004') %lower bob m2 of first (angle,velocity) plot
    hold on
    plot(y(:,5),y(:,6),'-c', 'Linewidth',1, 'DisplayName','upper bob of 0006') % upper bob m1 (angle,velocity) plot
    hold on
    plot(y(:,7),y(:,8),'-m','Linewidth',1, 'DisplayName','lower bob of 0006') %lower bob m2 (angle,velocity) plot
    hold on
    legend('upper bob of first','lower bob of first','upper bob of second','lower bob of second','Location', 'northeast')
    % Cosmetics for 2D phase plot graph
    set(gca,'FontSize',14)
    ylim([-20,20])
    xlim([-5,5])
    xlabel('Angle \theta (radians)')
    ylabel('Velocity (radians per second)')
    title('Double Pendulum Phase Space Plot')
    % 3D Phase Space
    figure('Name','Comparing Double Pendulum 3D Phase Space','NumberTitle','off')
    plot3(t,y(:,1),y(:,2),'-b', 'Linewidth',1, 'DisplayName','upper bob of first') % upper bob m1 of first (angle,velocity) plot
    hold on
    plot3(t,y(:,3),y(:,4),'-r','Linewidth',1, 'DisplayName','lower bob of first') %lower bob m2 of first (angle,velocity) plot
    hold on
    plot3(t,y(:,5),y(:,6),'-c', 'Linewidth',1, 'DisplayName','upper bob of second') % upper bob m1 (angle,velocity) plot
    hold on
    plot3(t,y(:,7),y(:,8),'-m','Linewidth',1, 'DisplayName','lower bob of second') %lower bob m2 (angle,velocity) plot
    hold on
    legend('upper bob of first','lower bob of first','upper bob of second','lower bob of second','Location', 'northeast')
% Cosmetics for 3D phase plot graph
    set(gca,'FontSize',14)
    xlim([0,11.98])
    zlim([-20,20])
    ylim([-5,5])
    ylabel('Angle \theta (radians)')
    zlabel('Velocity  (radians per second)')
    xlabel('Time (seconds)')
    title('Double Pendulum 3D Phase Space Plot')
end


function [t,y]= eulerODE(domain,y0,N)
% Euler's Method for Solving IVP dydt = f(t,y)
% y and f are possibly column vectors (i.e. arrays)
% Input: domain, initial value y0 (an array), number of steps N
% Output:  t and y for each step
t(1)=domain(1); 
y(1,:)=y0; % tells it to set row 1 as y0 vector.
h=(domain(2)-domain(1))/N;         % this is the step size
for i=1:N
  t(i+1)=t(i)+h;
  y(i+1,:)=ODEstep(t(i),y(i,:),h);
% above line states i+1th row = ODEstep of ith row. 
end
end

function y=ODEstep(t,y,h)
% RK4 method
s_1=dydt(t,y);
s_2=dydt(t+(h/2),y+(h/2)*s_1);
s_3=dydt(t+(h/2),y+(h/2)*s_2);
s_4=dydt(t+h,y+h*s_3);
y=y+(h/6)*(s_1+2*(s_2+s_3)+s_4);
% The current step is determined by the previous step added to the average of the previous step and the estimation of the current step.  
end


function f = dydt(t,y)
% Put the RHS of your ODE (system) here.
[m1,m2,m3,m4,l1,l2,l3,l4,g,u1,u2] = param();
% f and y are arrays
f(1) = y(2);
f(2) = (-g*sin(y(1))-u1*l2*y(4)^2*sin(y(1)-y(3))+u1*g*sin(y(3))*cos(y(1)-y(3))-u1*l1*y(2)^2*sin(y(1)-y(3))*cos(y(1)-y(3)))/(l1*(1-u1*(cos(y(1)-y(3)))^2));
f(3) = y(4);
f(4) = (-g*sin(y(3))+l1*y(2)^2*sin(y(1)-y(3))+g*sin(y(1))*cos(y(1)-y(3))+u1*l2*y(4)^2*sin(y(1)-y(3))*cos(y(1)-y(3)))/(l2*(1-u1*(cos(y(1)-y(3)))^2));
% The first f is column 2 of y? Yes it is.
f(5) = y(6);
f(6) = (-g*sin(y(5))-u2*l4*y(8)^2*sin(y(5)-y(7))+u2*g*sin(y(7))*cos(y(5)-y(7))-u2*l3*y(6)^2*sin(y(5)-y(7))*cos(y(5)-y(7)))/(l3*(1-u2*(cos(y(5)-y(7)))^2));
f(7) = y(8);
f(8) = (-g*sin(y(7))+l3*y(6)^2*sin(y(5)-y(7))+g*sin(y(5))*cos(y(5)-y(7))+u2*l4*y(8)^2*sin(y(5)-y(7))*cos(y(5)-y(7)))/(l4*(1-u2*(cos(y(5)-y(7)))^2));
% 
% Difference:
D = y(:,1)-y(:,5) + y(:,3)-y(:,7) + y(:,2)-y(:,6) + y(:,4)-y(:,8)
d = 1.63-1.64 % Initial difference d. Manulally adjust this value according to 
% initial conditions on imput.
% 
% lambda = (1/T)*log|sqrt(D^2)/({initial angle 1 separation)+(intial angle 2 separation)
% +(initial velocity 1 separation)+(initial velocity 2 separation})|
lambda = (1/t)*log(abs(D/d))
end

% Put the physical parameters here.
function [m1,m2,m3,m4,l1,l2,l3,l4,g,u1,u2] = param()
m1 = 0.09702; 
m2 = 0.08119;
m3 = 0.09702; 
m4 = 0.08119; 
l1 = 0.254; 
l2 = 0.2286; 
l3 = 0.254; 
l4 = 0.2286; 
g = 9.8; 
u1=m2/(m1+m2); 
u2=m4/(m3+m4);
% l1 is the length of the rod between the pivot and bob 1.
% l2 is the length of the rod between bob 1 and bob 2.
% l3 is the length of the rod between the pivot and bob 3.
% l4 is the length of the rod between bob 3 and bob 4.
% m1 is the mass/bob between the two rods of pendulum 1.
% m2 is the mass/bob at the end of the second rod of pendulum 1.
% m3 is the mass/bob between the two rods of pendulum 2.
% m4 is the mass/bob at the end of the second rod of pendulum 2.
% g is gravity.
end
