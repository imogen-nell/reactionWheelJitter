%parameters
m; %wheel mass

%torques
%rocking/radial
f_r;
z_r; % damping ??
w_r = 2 * pi * f_r;
k_r = I_rr * w_r ^ 2;
c_r = 4 * z * f * I;
theta_x;
theta_y;

%torques due to rocking
T = [k 0, k 0]*[theta_x theta_y];


%axial forces
f_ax;
w_ax = 2* pi * f_ax;
k_ax = m * w_ax^2;
c_ax = 4 * z_ax * f_ax * m;
