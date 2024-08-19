%%%%%%parameters2%%%%%%
%contants :  r_bc_b, masses , r_w_b, ihub_bc, irw_wc

%time dept :

%state vars: W_frame, r_b_nd , W, w_b_n


%masses/ kg
msc = 680 ; % space craft mass, estimate a
mrw = 12 ;% RW masses , used flywheel assembly mass in grams
mhub = 644 ; %estmate 

%hub
ihub_bc = [550 .1045 -0.084;.1045 650 .0001;-.084 .0001 650]; %random, hub inertia tensor abt com
r_bc_b = [.01 -.02 .10]'; %hub COM location wrt B, point on sc

%wheel properties
g   = [.7887 -.2113 -1/sqrt(3)]';
W_frame_init = [g, cross(g,[1 0 0]')/norm(cross(g,[1 0 0]')), cross(cross(g,[1 0 0]'),g )/norm(cross(cross(g,[1 0 0]'),g ))]; %wheel orientation matrix

U_s = 0.0000048 ; %wheel static imbal
U_d = 0 ;% .0000154 ; %dynamic imbal
d   =  U_s/mrw; % gc of rw to com  (not constant) (cm estimate)

%note irw_wc must take on this form due to assumptions made 
irw_wc = [1.5915 0 U_d ;0 0.8594 0; U_d 0 0.8594];%wheel inertial tensor 
r_w_b = [0.6309 -0.1691 0.4619]';%location vector of rw. B frame origin to rw origin / gc 
J11 = irw_wc(1,1);
J12 = 0;
J13 = U_d;
J22 = irw_wc(2,2);
J23 = irw_wc(2,3);
J33 = irw_wc(3,3);

%%FIX


%initial conditions

W_init = 200* 2*pi/60; % cap omega (random rad/s chosen) (initial wheel speed

us_0 = 0.005; 

%r_c_ndd = F / msc;

r_c_ndd = [0 0 0]'; % NO EXTERNAL FORCES ON SC

%%%%%%%%%%%%%%%%%%%%%%
save('parameters2.mat', 'g','r_c_ndd', 'msc', 'mrw', 'mhub', 'ihub_bc', 'r_bc_b', 'W_frame_init','W_init', 'd', 'U_s', 'U_d', 'irw_wc', 'r_w_b', 'J11', 'J12', 'J13', 'J22', 'J23', 'J33', 'us_0');
