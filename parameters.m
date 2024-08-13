%%%%%%parameters%%%%%%
%contants :  r_bc_b, masses , r_w_b, ihub_bc, irw_wc

%time dept :

%state vars: W_frame, r_b_nd , W, w_b_n


%masses/ kg
msc = 680 ; % space craft mass, estimate a
mrw = 12 ;% RW masses , used flywheel assembly mass in grams
mhub = 644 ; %estmate 

%hub
ihub_bc = [550 .1045 -0.084;.1045 650 .0001;-.084 .0001 650]; %random, hub inertia tensor abt com
r_bc_b = [.01 -.02 .10]; %hub COM location wrt B, point on sc

%wheel properties
W_frame_init = [.7887 -.2113 -.5774; -0.2113 .7887 -1/sqrt(3); 1/sqrt(3) 1/sqrt(3) 1/sqrt(3)]; %wheel orientation matrix
g = W_frame_init(:,1);
d = .0000004; % gc of rw to com  (not constant) (cm estimate)
U_s = 0.0000048 ; %wheel static imbal
U_d = 0.000154 ; %dynamic imbal
%note irw_wc must take on this form due to assumptions made 
irw_wc = [1.5915 0 1.54*10^(-6) ;0 0.8594 0; 1.54*10^(-6) 0 0.8594];%wheel inertial tensor 
r_w_b = [0.6309 -0.1691 0.4619];%location vector of rw. B frame origin to rw origin / gc 
J11 = irw_wc(1,1);
J12 = irw_wc(1,2);
J13 = irw_wc(1,3);
J22 = irw_wc(2,2);
J23 = irw_wc(2,3);
J33 = irw_wc(3,3);

%%FIX


%initial conditions

W_init = 0; %-558; % cap omega (random rad/s chosen) (initial wheel speed
%%FIX???
%Wd = 50 ; %angular acceleration set for RW 
%theta = 43 ; %(initial wheel angle deg )
us = 0.2 ; % motor torque SHOULD NOT BE A PARAM CHECK EQ 66 P2
%t_step = 0.1 ; % ms 
%jitter modeling

w_b_nd = [0 0 0]; %FIX THIS

%r_c_ndd = F / msc;

r_c_ndd = [0 0 0]'; % NO EXTERNAL FORCES ON SC

%%%%%%%%%%%%%%%%%%%%%%
save('parameters.mat', 'g','r_c_ndd','w_b_nd', 'msc', 'mrw', 'mhub', 'ihub_bc', 'r_bc_b', 'W_frame_init','W_init', 'd', 'U_s', 'U_d', 'irw_wc', 'r_w_b', 'J11', 'J12', 'J13', 'J22', 'J23', 'J33', 'us');
