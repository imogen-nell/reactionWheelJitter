%masses/ kg
msc = 680 ; % space craft mass, estimate a
mrw = [12 12 12 12] ;% RW masses , used flywheel assembly mass in grams
mhub = 644 ; %estmate 

%hub
ihub_bc = [550 .1045 -0.084;.1045 650 .0001;-.084 .0001 650]; %random, hub inertia tensor abt com
r_bc_b = [.01 -.02 .10]'; %hub COM location wrt B, point on sc

%wheel properties
W_frame_init1 = [.7887 -.2113 -1/sqrt(3); -0.2113 .7887 -1/sqrt(3); 1/sqrt(3) 1/sqrt(3) 1/sqrt(3)]; %wheel orientation matrix
W_frame_init2 = W_frame_init1; 
W_frame_init3 = W_frame_init1; 
W_frame_init4 = [1 0 0 ; 0 1 0 ; 0 0 1]; 

g = [W_frame_init1(:,1) W_frame_init2(:,1) W_frame_init3(:,1) W_frame_init4(:,1) ];
d = [.0000004 .0000004 .0000004 0.0000001]; 


U_s = 0.0000048 ; %wheel static imbal
U_d = 0.000154 ; %dynamic imbal

irw_wc = [1.5915 0 1.54*10^(-6) ;0 0.8594 0; 1.54*10^(-6) 0 0.8594] * 0.0001 ;%wheel inertial tensor 
r_w_b = [[0.6309 -0.1691 0.4619]' [0.6309 -0.1691 0.4619]' [0.6309 -0.1691 0.4619]' [0 0 0]'];%location vector of rw. B frame origin to rw origin / gc 
J11 = irw_wc(1,1);
J12 = irw_wc(1,2);
J13 = irw_wc(1,3);
J22 = irw_wc(2,2);
J23 = irw_wc(2,3);
J33 = irw_wc(3,3);

W_init = [0 50 -500 200]' ; %-558; % cap omega (random rad/s chosen) (initial wheel speed

us_0 = [0.06 ; 0; 0 ; 0] ; 

r_c_ndd = [0 0 0]'; % NO EXTERNAL FORCES ON SC

t1 = 0;
t2 = 0;
t3 = 0;
t4 = 0;


% L (u) = q* u q
%initial pointing direction
%let u0 = B1 of B frame (body fixed)
u0 = [1 0 0]';


%%%%%%%%%%%%%%%%%%%%%%
save('parameters3.mat', 'u0', 'W_init', 't1', 't2', 't3', 't4', 'g','r_c_ndd', 'msc', 'mrw', 'mhub', 'ihub_bc', 'r_bc_b', 'W_frame_init1','W_frame_init2','W_frame_init3','W_frame_init4', 'd', 'U_s', 'U_d', 'irw_wc', 'r_w_b', 'J11', 'J12', 'J13', 'J22', 'J23', 'J33', 'us_0');
