% ALL UNITS SI : KG, M, S, RAD, RAD/S, etc. unless otherwise stated 

%%%%%%%%%%%%%%%%%%%%%%%%%% INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%
%initial RPM of each rw INPUT RPM NOT RAD/s
rw1_w = 0 ; %RPM
rw2_w = 0 ;
rw3_w = 0 ;
rw4_w = 100 ;

W_init = [rw1_w rw2_w rw3_w rw4_w]'  * 2*pi/60 ; % initial wheel speed, converted to rad/s 

%initial comanded torques (Nm)
us_0 = [ 0 ; 0.001 ; 0 ; -0.005 ];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%masses/ kg
mrw = [12 12 12 12] ; % RW masses 
mhub = 644 ; %sc weight excluding rws
msc = mhub + mrw(1) + mrw(2) + mrw(3) + mrw(4) ; % space craft mass, 

%hub  
ihub_bc = [550 .1045 -0.084;.1045 650 .0001;-.084 .0001 650]; % hub inertia tensor about HUB COM
r_bc_b = [.01 -.02 .10]'; %hub COM location wrt Body-fixed (B) frame origin


g = [[.7887 -.2113 -1/sqrt(3)]', [-0.2113 .7887 -1/sqrt(3)]', [1/sqrt(3) 1/sqrt(3) 1/sqrt(3)]', [-1/sqrt(3) -1/sqrt(3) -1/sqrt(3)]'];  %columns of g are the axis of rotation of each rw UNITARY VECTORS 


%build frame components for each wheel . 
%W_frame(:, 1) = axis of rotation,
%W_frame(:,2) is perpendicular to axis of rotation, and frame(:,3) completes the right hand rule 
% make sure axis of rotations are not parallel to the vector used to make
% axis 2 & 3 ([1 0 0] chosen randomly) 
W_frame_init1  =  [g(:,1), cross(g(:,1),[1 0 0]')/norm(cross(g(:,1),[1 0 0]')), cross(cross(g(:,1),[1 0 0]'),g(:,1)) / norm(cross(cross(g(:,1),[1 0 0]'),g(:,1)) )]; %wheel orientation matrix
W_frame_init2  =  [g(:,2), cross(g(:,2),[1 0 0]')/norm(cross(g(:,2),[1 0 0]')), cross(cross(g(:,2),[1 0 0]'),g(:,2)) / norm(cross(cross(g(:,2),[1 0 0]'),g(:,2)) )]; %wheel orientation matrix
W_frame_init3  =  [g(:,3), cross(g(:,3),[1 0 0]')/norm(cross(g(:,3),[1 0 0]')), cross(cross(g(:,3),[1 0 0]'),g(:,3)) / norm(cross(cross(g(:,3),[1 0 0]'),g(:,3)) )]; %wheel orientation matrix
W_frame_init4  =  [g(:,4), cross(g(:,4),[0 1 0]')/norm(cross(g(:,4),[0 1 0]')), cross(cross(g(:,4),[0 1 0]'),g(:,4)) / norm(cross(cross(g(:,4),[0 1 0]'),g(:,4)) )]; %wheel orientation matrix

%wheel properties 
U_s = [0.0000048 0.0000048 0.0000048 0.0000048]; %wheel static imbal
U_d = [0.0000154 0.0000154 0.0000154 0.0000154]; %dynamic imbal
d = [U_s(1)/mrw(1) U_s(2)/mrw(2) U_s(3)/mrw(3) U_s(4)/mrw(4)]; %derived from wheel properites, d is the distance from the the rw com to its axis of rotation
irw_wc = [1.5915 0 U_d(1) ;0 0.8594 0; U_d(1) 0 0.8594] * 0.0001 ;%wheel inertial tensor, currently the same for all wheels 
r_w_b = [ [.01 -.02 .10]' [0.6309 -0.1691 0.4619]' [0.3309 -0.1691 0.2619]' [0.6309 -0.1691 0.4619]' ];%location vector of rw wrt B frame origin. eg r_w_b(:,3) = location of (geometric) center of rw #3 wrt to B frame, note this vector is contant

J11 = irw_wc(1,1);
J12 = 0;
J13 = irw_wc(1,3); % = U_d
J22 = irw_wc(2,2);
J23 = irw_wc(2,3);
J33 = irw_wc(3,3);


F = [0 0 0]' ; %EXTERNAL FORCE 
r_c_ndd = F / msc; % ASSUMED NO EXTERNAL FORCES ON SC, (see eq (1) from paper) 

%initial wheel angles of rws, doesnt matter since frame axis ae arbitrary,
% probably will never need to change these 
t1 = 0;
t2 = 0;
t3 = 0;
t4 = 0;

save('parameters3.mat', 'W_init', 't1', 't2', 't3', 't4', 'g','r_c_ndd', 'msc', 'mrw', 'mhub', 'ihub_bc', 'r_bc_b', 'W_frame_init1','W_frame_init2','W_frame_init3','W_frame_init4', 'd', 'U_s', 'U_d', 'irw_wc', 'r_w_b', 'J11', 'J12', 'J13', 'J22', 'J23', 'J33', 'us_0');
