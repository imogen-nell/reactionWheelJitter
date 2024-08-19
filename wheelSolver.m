dt = 0.01; 
tspan = 0:dt:5;
y0 = 0
[t,y] = ode45(@func,tspan,y0);


function dx = func(t,y)
    dx = 200  * 2*pi/60; %%units rad/s
end

figure;
hold on;
% 1 rpm = 2pi / 60 rad /s 
% First subplot
subplot(1,1,1); % 4 rows, 1 column, 1st plot
plot(t, y(:), '-m', 'DisplayName', 'RW 1 ');
title('RW 1: Wheel Speed / RPM ', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'g');

% Second subplot
subplot(2,1,2); % 4 rows, 1 column, 2nd plot
plot(t, rad2deg(mod(y(:),2*pi)), '-c', 'DisplayName', 'RW 2 ');
title('RW 2: Wheel Speed / RPM ', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'g');

