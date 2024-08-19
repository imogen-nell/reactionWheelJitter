%mod 14 because mod13 is messy and wrong 
%now actually adding yaw pitch roll and MRP even though i do not know what
%mrp means physically 

load("parameters3.mat")
%%%%%%%%%%%%%%%%%%%%%%%%%% INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%
dt = 0.01; %doesnt matter, ODE45 will chose timestep 
tf = 6; %sim time, seconds
tspan = 0:dt:tf;
w_b_n_0 = [0 0 0]'; %inital angular velocity vector of space craft.
rd_0 = [0 0 0]'; %initial translational velocity vector of space craft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global us
us = [NaN, NaN, NaN, NaN];

%build ititial vector for ODE45 : 
theta_0  = [t1 t2 t3 t4]';  %rw angle displacement
pyr_0    = [0 0 0]'; %initial pitch yaw roll
rCOM_0   = rd_0 + (1/msc) * ( mhub * r_bc_b...
    + mrw(1) *(r_w_b(:,1) + d(1) * W_frame_init1(:,2)) ...
    + mrw(2) *(r_w_b(:,2) + d(2) * W_frame_init2(:,2)) ...
    + mrw(3) *(r_w_b(:,3) + d(3) * W_frame_init3(:,2))...
    + mrw(4) *(r_w_b(:,4) + d(4) * W_frame_init4(:,2))) ; %com of space craft + rw mass.
W_0     = W_init; %wheel speeds 
mrp_0   = [0 0 0]';

y0      = [ w_b_n_0 ; rd_0 ; W_0  ; theta_0 ; pyr_0 ; rCOM_0 ; mrp_0];
% 20x1  = [3x1      ; 3x1  ; 4x1  ; 4x1     ; 3x1   ; 3x1    ; 3x1  ] <- sizes of things 
% y0    = [1:3      ; 4:6  ; 7:10 ; 11:14   ; 15:17 ; 18:20  ; 21:23] <- indexing for reference   

%%%%%%%%%SOLVE%%%%%%%%%%
[t,y] = ode45(@func,tspan,y0);
%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%Graph%%%%%%%%%%
graphMRP(t, y(:,21:23))  % modified rodrigues parameters (quaternions) 
graphWBN(t, y)          % w_b_n -> graph componenets of spacecraft omega
graphPYR(t, y)           % graph pitch, yaw, roll
graphWHEELSPEED(t,y)     % graph wheel speeds 
graphRCOM(t, y(:,18:20)) % SANITY CHECK COM DOES NOT MOVE, unless you give sc initial translational speed
graph3D(y(1,18), y(1,19), y(1,20), y(:,1:3), y(1,7:10)) % see directions of omegas, black arrow: rws, colored arrows: spacecraft
%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%Functions%%%%%%%%%%
%ode function,
function dy = func(t,y)
%all equations have a number next to them (XX) corresponding to eqn number
%from the paper, if no number is present next to an eqn, I made it up and
%it might be wrong 
    load('parameters3.mat')
    global us
    if ismember(1, isnan(us))
        us = us_0;
        fprintf('US reset to US_0 %f\n', us);
        %should only print once per each reaction wheel 
    end

    %y0 = [     1:3 ;  4:6  ; 7:10  ; 11:14   ;        15:17   ; 18:20 ] 
    %y0 = [ w_b_n_0 ; rd_0  ; W_0   ; theta_0 ; pitch-yaw-roll ; mrp   ]
     
    %exctract data from y 
    w_b_n   = y(1:3);
    r_b_nd  = y(4:6);
    W       = y(7:10);
    thetas  = y(11:14);
    pyr     = y(15:17);
    mrp     = y(18:20);
   %find the orientation of the reaction wheels by rotating their initial
   %direction vectors by their angle of displacement. probably there is a
   %better way to do this 
   
   %columns of wframe_2 are the second axis of each rw (not the axis of
   %rotation)
   %cols of wframe_3 are the 3rd axis of each rw (not the axis of rotation)
   %eg wframe_2(:,2) = axis 2 of rw 2 
   %eg wframe_3(:,1) = axis 3 of rw 1 
   %size 3x4 
    wframe_2 = [rotateg(W_frame_init1(:,2), thetas(1),  g(:, 1)) , ...
        rotateg(W_frame_init2(:,2), thetas(2), g(:, 2)) , ...
        rotateg(W_frame_init3(:,2), thetas(3), g(:, 3)), ...
        rotateg(W_frame_init4(:,2), thetas(4), g(:, 4))];
    wframe_3 = [rotateg(W_frame_init1(:,3), thetas(1), g(:, 1)) , ...
        rotateg(W_frame_init2(:,3), thetas(2), g(:, 2)) ,...
        rotateg(W_frame_init3(:,3), thetas(3), g(:, 3)) , ...
        rotateg(W_frame_init4(:,3), thetas(4), g(:, 4)) ];
    
    %columns of r_wc_b are the poition vectors from the com of each rw to
    %the b frame origin'
    % eg r_wc_b(:,3) = position vector of com of rw 3 wrt B frame
    % size 3x4
    r_wc_b = [ r_w_b(:,1) + d(1) * wframe_2(:, 1)  , ...
        r_w_b(:,2) + d(2) * wframe_2(:, 2) , ...
        r_w_b(:,3) + d(3) * wframe_2(:, 3) , ...
        r_w_b(:,4) + d(4) * wframe_2(:, 4) ] ;  %(6)

    %r_wc_b_dash is the B frame (non inertial) derivative of r_wc_b
    % eg r_wc_b_dash(:,3) = B frame derivative of position vector of com of rw 3 wrt B frame
    % size 3x4
    r_wc_b_dash  = [d(1) * W(1) * wframe_3(:,1) ,...
        d(2) * W(2) * wframe_3(:,2) , ...
        d(3) * W(3) * wframe_3(:,3) , ...
        d(4) * W(4) * wframe_3(:,4) ] ;    %(7)

    %inertia tensor for hub about b frame origin (not about com)
    %size 3x3
    ihub_b = ihub_bc + mhub * skew(r_bc_b) * skew(r_bc_b)'; %(eqn 18 from paper1) but its just parallel axis thm 
    
    % irw_wc assumed all the same for each rw 
    % irw_b contains the inertia tensors for each rw about the B frame
    % origin
    % eg irw_b(1:3,1:3) = inertia tensor for rw 1 wrt B frame origin
    % eg irw_b(1:3,4:6) = inertia tensor for rw 2 ... etc.
    %size = 3x12
    irw_b = [irw_wc  + mrw(1)  * skew(r_wc_b(:,1)) * skew(r_wc_b(:,1))' , ...
             irw_wc  + mrw(2)  * skew(r_wc_b(:,2)) * skew(r_wc_b(:,2))' , ...
             irw_wc  + mrw(3)  * skew(r_wc_b(:,3)) * skew(r_wc_b(:,3))' , ...
             irw_wc  + mrw(4)  * skew(r_wc_b(:,4)) * skew(r_wc_b(:,4))' ]; %(26) p1 

    % space craft inertia tensor about b frame origin
    % size = 3x3 
    isc_b  = ihub_b  + irw_b(1:3,1:3) + irw_b(1:3,4:6) + irw_b(1:3,7:9) + irw_b(1:3,10:12);
    
    % Wheel frame time derivative of wheel inertia tensor about wheel com 
    % irw_wc_dash1 = inertia derivative for rw 1 
    % irw_wc_dash2 = ... for rw 2, etc. 
    % irw_wc_dash2 size = 3x3
    irw_wc_dash1 = (J12 *  wframe_3(:,1) - J13*wframe_2(:,1)) * W(1) * g(:,1)'...
        + (-J13 *  g(:,1) + J22 * wframe_3(:,1) -2 * J23 * wframe_2(:,1)...
        - J33 & wframe_3(:,1)) * W(1) * wframe_2(:,1)'...
        + (J12 * g(:,1) + J22 * wframe_2(:,1) + 2* J23 * wframe_3(:,1)...
        - J33 * wframe_2(:,1)) * W(1) * wframe_3(:,1)'; %29
    irw_wc_dash2 = (J12 *  wframe_3(:,2) - J13*wframe_2(:,2)) * W(2) * g(:,2)'...
        + (-J13 *  g(:,2) + J22 * wframe_3(:,2) -2 * J23 * wframe_2(:,2)...
        - J33 & wframe_3(:,2)) * W(2) * wframe_2(:,2)'...
        + (J12 * g(:,2) + J22 * wframe_2(:,2) + 2* J23 * wframe_3(:,2)...
        - J33 * wframe_2(:,2)) * W(2) * wframe_3(:,2)';
    irw_wc_dash3 = (J12 *  wframe_3(:,3) - J13*wframe_2(:,3)) * W(3) * g(:,3)'...
        + (-J13 *  g(:,3) + J22 * wframe_3(:,3) -2 * J23 * wframe_2(:,3)...
        - J33 & wframe_3(:,3)) * W(3) * wframe_2(:,3)'...
        + (J12 * g(:,3) + J22 * wframe_2(:,3) + 2* J23 * wframe_3(:,3)...
        - J33 * wframe_2(:,3)) * W(3) * wframe_3(:,3)';
    irw_wc_dash4 = (J12 *  wframe_3(:,4) - J13*wframe_2(:,4)) * W(4) * g(:,4)'...
        + (-J13 *  g(:,4) + J22 * wframe_3(:,4) -2 * J23 * wframe_2(:,4)...
        - J33 & wframe_3(:,4)) * W(4) * wframe_2(:,4)'...
        + (J12 * g(:,4) + J22 * wframe_2(:,4) + 2* J23 * wframe_3(:,4)...
        - J33 * wframe_2(:,4)) * W(4) * wframe_3(:,4)';

    % ihub_b does not change, therefore its derivative is 0. i just added
    % this for mysself 
    ihub_b_dash = zeros(3,3);

    % irw_b_dashX is the non inertial time derivative of rw inertia tensor about b frame 
    % irw_b_dash1 for rw, irw_b_dash2 for rw 3 etc.
    % size of each 3x3 
    irw_b_dash1  = irw_wc_dash1 + mrw(1) * skew(r_wc_b_dash(:,1))...
        *skew(r_wc_b(:,1))' ...
        + mrw(1) * skew(r_wc_b(:,1))*skew( r_wc_b_dash(:,1))';     %(35)
    irw_b_dash2  = irw_wc_dash2 + mrw(2) * skew(r_wc_b_dash(:,2))...
        *skew(r_wc_b(:,2))' ...
        + mrw(2) * skew(r_wc_b(:,2))*skew( r_wc_b_dash(:,2))';
    irw_b_dash3  = irw_wc_dash3 + mrw(3) * skew(r_wc_b_dash(:,3))...
        *skew(r_wc_b(:,3))' ...
        + mrw(3) * skew(r_wc_b(:,3))*skew( r_wc_b_dash(:,3))';
    irw_b_dash4  = irw_wc_dash4 + mrw(4) * skew(r_wc_b_dash(:,4))...
        *skew(r_wc_b(:,4))' ...
        + mrw(4) * skew(r_wc_b(:,4))*skew( r_wc_b_dash(:,4))';


    % isc_b_dash is the non inertial time derivative of space craft 
    % (incl. hub and  rw) inertia tensor about b frame
    % size 3x3
    isc_b_dash  = ihub_b_dash + irw_b_dash1+ irw_b_dash2 + irw_b_dash3+ irw_b_dash4;
    
    % c = vector from b frame origin to com of spacecraft (hub + rws)
    % cdash = b frame (non inertial) time derivative of c
    % size 3x1
    [c,cdash] = getcs(r_wc_b, r_wc_b_dash);

    % find A to find e, see paper for details 
    A = zeros(4,4) ;  
    for i = 1:4
        for j = 1:4
            if i == j
                A(i,j) = J11 + mrw(i)*(1 - mrw(i)/msc)*d(i)^2 ;
            else
                A(i,j) = - (mrw(i) * d(i) * wframe_3(:,i)'/msc) * ...
                    mrw(j) * d(j) * wframe_3(:,j);
            end
        end 
    end
    e = inv(A); % (71)

    % preset size then fill
    F = zeros(4,3); % (72)
    for i = 1:4
        F(i,:) = [- ((J11 + mrw(i) * d(i)^2) * g(:,i)' + J12*wframe_2(:,i)' + ...
        J13 * wframe_3(:,i)'+ (mrw(i)*d(i)* wframe_3(:,i)') * (skew(c) - ....
        skew(r_w_b(:,i))))]'; %(72)
    end
    
    % see paper
    i_lhs = isc_b + msc * skew(c)*skew(c); %(80)
    for i = 1:4

        i_lhs = i_lhs +( irw_wc*g(:,i) ...
        +mrw(i) * d(i) * (skew(r_wc_b(:,i)) ...
        - skew(c))*wframe_3(:,i))*e(i,:)*F;      %(80)
    end
    
    % get the wheel frame (magnitude) components of w_b_n, spacecraft angular velocity
    % w_1 = rw spin axis components of w_b_n for each rw, size = 1x4
    % w_2 = rw axis 2    components of w_b_n for each rw, size = 1x4
    % w_3 = rw axis 3    components of w_b_n for each rw, size = 1x4
    % eg w_2(3) = magnitude of component of spacecraft angular velocity in
    % axis 2 dirrection of rw 3.
    % eg w_1(4) = magnitude of component of spacecraft angular velocity in
    % axis 1 (spin axis) dirrection of rw 4.
    [w_1, w_2, w_3] = ws(w_b_n, g ,wframe_2, wframe_3);%(45)-(47)
    
    % see paper
    v = [0 0 0 0]'; %(73)
    
    for i = 1:4
        sum = 0;
        for j = 1:4
            if i ~=j
                sum = sum + mrw(j)*d(j)*W(j)^2*wframe_2(:,j);
            end
        end
        v(i) = -mrw(i)*d(i)*wframe_3(:,i)' * (r_c_ndd - 2 * skew(w_b_n) * cdash - ...
            skew(w_b_n) * skew(w_b_n)* c + (1/msc)*(sum))...
            +J23*(w_3(i)^2 - w_2(i)^2) + w_1(i)*(J12 * w_3(i) - J13 * w_2(i)) ...
            + w_2(i) * w_3(i) * ( J22 - J33 - mrw(i) * d(i)^2) - mrw(i) * d(i) * ...
            wframe_3(:,i)' * skew(w_b_n)*skew(w_b_n) * r_w_b(:,i) + us(i) ;      %(73)
    end

    %%FIX LB BECUASE CHOSING U,V CHANGES THINGS ACTUALLY AND I DONT THINK
    %%IT SHOULD 
    %(83) (84) ????????? I MADE THIS UP I DONT KNOW WHAT THIS SHOULD BE
    % ITS NOT DEFINED IN THE PAPER 
    %Lb is the torque on the spacecraf from rws 
    %size 3x1
    Lb = W(1)^2 * (U_d(1) * wframe_3(:,1) + U_s(1)*skew(r_w_b(:,1))*wframe_2(:,1))+ ...
         W(2)^2 * (U_d(2) * wframe_3(:,2) + U_s(2)*skew(r_w_b(:,2))*wframe_2(:,2))+...
         W(3)^2 * (U_d(3) * wframe_3(:,3) + U_s(3)*skew(r_w_b(:,3))*wframe_2(:,3))+...
         W(4)^2 * (U_d(4) * wframe_3(:,4) + U_s(4)*skew(r_w_b(:,4))*wframe_2(:,4));

    %see paper
    tau_rhs = - skew(w_b_n) * isc_b * w_b_n - isc_b_dash * w_b_n ...
        -msc * skew(c)*(r_c_ndd - 2 * skew(w_b_n) * cdash - ...
        skew(w_b_n)*skew(w_b_n)*c) + Lb; % (81) INCOMPLETE SEE BELOW 
    
    % collecting irw_wc_dashs for use in building tau_rhs 
    irw_wc_dash = [irw_wc_dash1 irw_wc_dash2 irw_wc_dash3 irw_wc_dash4];

    for i = 1:4
        tau_rhs = tau_rhs + mrw(i) * d(i) * W(i)^2 * skew(r_wc_b(:,i)) ...
            * wframe_2(:,i) -  irw_wc_dash(:,i*3 - 2 :i*3 ) * W(i) * g(:,i) - ...
            skew(w_b_n) * (irw_wc * W(i) * g(:,i) + mrw(i) * skew(r_wc_b(:,i))*...
            r_wc_b_dash(:,i)) - mrw(i) * d(i)* W(i)^2 * skew(c) * wframe_2(:,i) - ...
            (irw_wc*g(:,i) + mrw(i) * d(i) * (skew(r_wc_b(:,i)) - skew(c)) ...
            *wframe_3(:,i))* e(i,:) * v  ;                              %(81)
    end 

    %inertial time derivate of w_b_n
    % size 3x1
    w_b_nd = inv(i_lhs)*tau_rhs; %(79)

    %inertial time derivative of rw speeds (magnitude)
    % eg Wd(3) = rw 3 acceleration
    % size 4x1
    %Wd = e * F * w_b_nd + e * v %(75)
    Wd = A\ F * w_b_nd + A \ v; %(75)

    % non inertial b frame 2nd time derivate of r_wc_b
    % size 3x4
    %eg r_wc_b_ddash(:,3) =  acceleration of com of rw 3 wrt b frame origin
    r_wc_b_ddash = [d(1) * Wd(1) * wframe_3(:,1) - d(1) * W(1)^2 * wframe_2(:,1) ,...
        d(2) * Wd(2) * wframe_3(:,2) - d(2) * W(2)^2 * wframe_2(:,2) ,...
        d(3) * Wd(3) * wframe_3(:,3) - d(3) * W(3)^2 * wframe_2(:,3) ,...
        d(4) * Wd(4) * wframe_3(:,4) - d(4) * W(4)^2 * wframe_2(:,4) ]; %(8)
    
    %cddash = 2nd order B frame (non inertial) time derivative of c
    %size 3x1
    cddash  = (1/msc) *(mrw(1)* r_wc_b_ddash(:,1) +...
        mrw(2) * r_wc_b_ddash(:,2) +...
        mrw(3) * r_wc_b_ddash(:,3) +...
        mrw(4) * r_wc_b_ddash(:,4)) ; %(5)
    
    % cdd = 2nd order N frame (inertial) time derivative of c
    % size 3x1
    cdd     = cddash + 2 * cross(w_b_n, cdash) +...
        cross(w_b_nd, c) + cross(w_b_n, cross(w_b_n, c)) ; %(9)
    
    % r_b_ndd = inertial translational acceleration of spacecraft's B frame, will
    % change as spacecraft rotates (b frame not coincident with space craft
    % com) 
    % to check only rotation is occuring graph COM translation 
    r_b_ndd = r_c_ndd- cdd ; 
    
    %CHECK US IS CORRECT
    if norm(J11*(g' * w_b_nd + Wd) - us) > 0.00001
        fprintf("error")
        norm(us -  J11*(g' * w_b_nd + Wd))
       
    end
        
    %derivative of pitch, yaw, roll. this eqn is available online
    %size 3x1
    pyrd = ( 1/cos(pyr(2)) )  * [0 ,sin(pyr(3)),cos(pyr(3)) ; ...
        0, cos(pyr(3)) * cos(pyr(2)) ,sin(pyr(3)) * cos(pyr(2));...
        cos(pyr(2)),sin(pyr(3))*sin(pyr(2)) ,cos(pyr(3)) * sin(pyr(2))] * w_b_n;

    %i think i made this up, used transport thm, pretty sure its right???
    %inertial time derivative of com of spacecraft incl. hub and rw. 
    rcomd  =  r_b_nd + (1/msc) * ( d(1) * mrw(1) *  W(1) * wframe_3(:,1)  + cross(w_b_n,wframe_2(:,1)) + ...
        d(2) * mrw(2) *  W(2) * wframe_3(:,2)  + cross(w_b_n,wframe_2(:,2))+...
        d(3) * mrw(3) *  W(3) * wframe_3(:,3)  + cross(w_b_n,wframe_2(:,3))+...
        d(4) * mrw(4) *  W(4) * wframe_3(:,4)  + cross(w_b_n,wframe_2(:,4)));

    %mrp derivative, available online
    mrpd = ((1 + norm(mrp)^2)/4) * ( eye(3) + 2 * (skew(mrp) + skew(mrp) ...
        * skew(mrp)/(1 + norm(mrp)^2))) * w_b_n ; 
    %%%%%%%%%%%%%%%%%%%%%%%% COLLECT DERIVATIVES%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %y0 = [ w_b_n ; r_b_nd  ; W  ; thetas ; pyr  ; rcom  ; mrp  ] ;
    dy  = [w_b_nd ; r_b_ndd ; Wd ; W      ; pyrd ; rcomd ; mrpd ] ;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

%returns skew symetric matrix of vector v 
function S = skew(v)
    S = [  0   -v(3)  v(2); v(3)   0   -v(1); -v(2)  v(1)   0  ];
end

function [c,cdash] = getcs(r_wc_b, r_wc_b_dash)
    load('parameters3.mat')
    c      = (1/msc) * ( mhub * r_bc_b  ...
        + mrw(1) * r_wc_b(:,1)...
        + mrw(2) * r_wc_b(:,2) ...
        + mrw(3) * r_wc_b(:,3)...
        + mrw(4) * r_wc_b(:,4));  %(3)
    cdash  = (1/msc) * ( mrw(1) * r_wc_b_dash(:,1) ...
        +mrw(2) * r_wc_b_dash(:,2) ...
        +mrw(3) * r_wc_b_dash(:,3) ...
        +mrw(4) * r_wc_b_dash(:,4));            %(4) 
end

function [w_1, w_2, w_3] = ws(w_b_n, g ,wframe_2, wframe_3)
%w_1 = [ w_s1 rw1; w_s1 rw2; etc]
    w_1 = [g(:,1)' * w_b_n , g(:,2)' * w_b_n, g(:,3)' * w_b_n, g(:,4)' * w_b_n];
    w_2 = [wframe_2(:,1)' * w_b_n, wframe_2(:,2)' * w_b_n, wframe_2(:,3)' * w_b_n, wframe_2(:,4)' * w_b_n ];
    w_3 = [wframe_3(:,1)' * w_b_n, wframe_3(:,2)' * w_b_n, wframe_3(:,3)' * w_b_n, wframe_3(:,4)' * w_b_n ] ; 
    %(45)-(47)
end

% rodriguez formula i think
% rotates vector v an angle theta about axis of rotation g
function v_rot = rotateg(v, theta, g)
    g = g / norm(g); % just making sure, g shoudl aready be unitary 
    v_rot = dot(v, g) * g + cos(theta) * (v - dot(v, g) * g ) + sin(theta) * cross(g, v);
end


function graphMRP(t,mrp)
    figure;
    title('Modified Rodriguez Parameters', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'g');
    xlabel('Time - Seconds', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Angle - Degrees', 'FontSize', 12, 'FontWeight', 'bold');
    hold on;
    plot(t, mrp(:, 1), '-b', 'DisplayName', '1');
    plot(t, mrp(:, 2), '-g', 'DisplayName', '2');
    plot(t, mrp(:, 3), '-c', 'DisplayName', '3');

    hold off;

end


function graphWBN(t,y)
    figure;
    title('Space Craft Angular Velocity', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'g');
    xlabel('Time - Seconds', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Angular Velocity - Rad/s', 'FontSize', 12, 'FontWeight', 'bold');
    hold on;
    plot(t, y(:, 1), '-b', 'DisplayName', 'W_b_n (1)');
    plot(t, y(:, 2), '-g', 'DisplayName', 'W_b_n (2)');
    plot(t, y(:, 3), '-c', 'DisplayName', 'W_b_n (3)');
    hold off;
end


function graphPYR(t, y)
    figure;
    title('Pitch Yaw Roll / ArcSeconds', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'g');
    xlabel('Time - Seconds', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Angle - Arcsecond', 'FontSize', 12, 'FontWeight', 'bold');
    hold on;
    plot(t, rad2arc(y(:, 15)), '-m', 'DisplayName', 'Pitch');
    plot(t, rad2arc(y(:, 16)), '-g', 'DisplayName', 'Yaw');
    plot(t, rad2arc(y(:, 17)), '-b', 'DisplayName', 'Roll');
    hold off;
end 

%DONT BOTHER GRAPHING THIS IF TORQUE IS APPLIED WHEEL IS TOO FAST 
function graphWHEELANGLE(t,y)
    figure;
    
    % First subplot
    subplot(4,1,1); % 4 rows, 1 column, 1st plot
    plot(t, rad2deg(mod(y(:,11),2*pi)), '-m', 'DisplayName', 'RW 1 ');
    title('RW 1: Wheel angle / DEG ', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'g');
    
    
    % Second subplot
    subplot(4,1,2); % 4 rows, 1 column, 2nd plot
    plot(t, rad2deg(mod(y(:,12),2*pi)), '-c', 'DisplayName', 'RW 2');
    title('RW 2: Wheel angle / DEG ', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'g');
    
    % Third subplot
    subplot(4,1,3); % 4 rows, 1 column, 3rd plot
    plot(t, rad2deg(mod(y(:,13),2*pi)), '-r', 'DisplayName', 'RW 3');
    title('RW 3: Wheel angle / DEG ', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'g');
    
    % Fourth subplot
    subplot(4,1,4); % 4 rows, 1 column, 4th plot
    %plot(t, rad2deg(mod(y(:,14),2*pi)), '-g', 'DisplayName', 'RW 4');
    plot(t, rad2deg(mod(y(:,14),2*pi)), '-g', 'DisplayName', 'RW 4');
    title('RW 4: Wheel angle / DEG ', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'g');

end

function graphWHEELSPEED(t,y)
    figure;
    % First subplot
    subplot(4,1,1); % 4 rows, 1 column, 1st plot
    plot(t, y(:,7)*60/(2*pi), '-m', 'DisplayName', 'RW 1 '); % convert rads/s to rpm
    title('RW 1: Wheel Speed / RPM', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'm');

    
    % Second subplot
    subplot(4,1,2); % 4 rows, 1 column, 2nd plot
    plot(t, y(:,8)*60/(2*pi), '-c', 'DisplayName', 'RW 2 ');
    title('RW 2: Wheel Speed / RPM', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'c');

    % Third subplot
    subplot(4,1,3); % 4 rows, 1 column, 3rd plot
    plot(t, y(:,9)*60/(2*pi), '-r', 'DisplayName', 'RW 3 ');
    title('RW 3: Wheel Speed / RPM', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'r');

    % Fourth subplot
    subplot(4,1,4); % 4 rows, 1 column, 4th plot
    plot(t, y(:,10)*60/(2*pi), '-g', 'DisplayName', 'RW 4 ');
    title('RW 4: Wheel Speed / RPM ', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'g');


end


function graph3D(x0, y0, z0, omega, rwommegs)
    load('parameters3.mat')
    [dim, ~] = size(omega);

    cols = [ 'r' 'y'  'g' 'b' 'c' 'm' 'k'  ];
    % red yellow green  blue light-blue pink 
    [lim ~] = size(cols');
    figure
    title('Angular Velocity', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'g');

    scale = 5000;
    %%plot sc OMEGA 
    jump = round( dim / lim );
    for i =1: lim -1
        label = sprintf('Space craft, t = %d', i);
        quiver3(x0, y0, z0, omega(i*jump,1)*scale , omega(i*jump,2)*scale,  omega(i*jump,3)*scale, cols(i), 'LineWidth', 2,'DisplayName', label);
        title('Angular Velocity of Space Craft and RWs', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'g');

        hold on;
    end
    %plot rw omega 
    
    rwx1 = r_w_b(1,1);
    rwy1 = r_w_b(2,1);
    rwz1 = r_w_b(3,1);
    
    rwx2 = r_w_b(1,2);
    rwy2 = r_w_b(2,2);
    rwz2 = r_w_b(3,2);

    rwx3 = r_w_b(1,3);
    rwy3 = r_w_b(2,3);
    rwz3 = r_w_b(3,3);

    rwx4 = r_w_b(1,4);
    rwy4 = r_w_b(2,4);
    rwz4 = r_w_b(3,4);

    secondary = .0025;
    rw1scale = secondary * rwommegs(1); % rw scale
    rw2scale = secondary * rwommegs(2);
    rw3scale = secondary * rwommegs(3);
    rw4scale = secondary * rwommegs(4);
    

    extrascale  = 100; % for reaction wheels

    if rw1scale == 0 && us_0(1) ~= 0
        quiver3(rwx1, rwy1, rwz1, g(1, 1) *extrascale * us_0(1),g(2,1)*extrascale* us_0(1),  g(3,1)*extrascale* us_0(1),  'k', 'LineWidth', 1,'DisplayName', 'RW 1');
        hold on;
        fprintf('rw 1, initial torque only \n');
    elseif rw1scale ~=0 
        rw1scale = .1;
        quiver3(rwx1, rwy1, rwz1, g(1, 1) * rw1scale,g(2,1) * rw1scale,  g(3,1) * rw1scale,  'k', 'LineWidth', 1,'DisplayName', 'RW 1');
        fprintf('rw 1, initial speed only\n');
        hold on;
    end

    if rw2scale == 0 && us_0(2) ~= 0
        extrascale =100;
        quiver3(rwx2, rwy2, rwz2, g(1, 2)*extrascale* us_0(2) ,g(2,2)*extrascale* us_0(2),  g(3,2)*extrascale* us_0(2),  'k', 'LineWidth', 1,'DisplayName', 'RW 2');
        fprintf('rw 2, initial torque only \n');
        hold on;
    elseif rw2scale ~=0 
        quiver3(rwx2, rwy2, rwz2, g(1, 2) * rw2scale,g(2,2) * rw2scale,  g(3,2) * rw2scale,  'k', 'LineWidth', 1,'DisplayName', 'RW 2 ');
        fprintf('rw 2, initial speed only\n');
        hold on;
    end

    if rw3scale == 0 && us_0(3) ~= 0
        extrascale = 1000;
        fprintf('rw 3, initial torque only \n');
        quiver3(rwx3, rwy3, rwz3, g(1, 3)*extrascale* us_0(3) ,g(2,3)*extrascale* us_0(3),  g(3,3)*extrascale* us_0(3),  'k', 'LineWidth', 1,'DisplayName', 'RW 3');
        hold on;
    elseif rw3scale ~=0 
        fprintf('rw 3, initial speed only\n');
        quiver3(rwx3, rwy3, rwz3, g(1, 3) * rw3scale,g(2,3) * rw3scale,  g(3,3) * rw3scale,  'k', 'LineWidth', 1,'DisplayName', 'RW 3 ');
        hold on;
    end

    if rw4scale == 0 && us_0(4) ~= 0
        fprintf('rw 4, initial torque only \n');
        quiver3(rwx4, rwy4, rwz4, g(1, 4) *extrascale* us_0(4),g(2,4)*extrascale* us_0(4),  g(3,4)*extrascale* us_0(4),  'k', 'LineWidth', 1,'DisplayName', 'RW 4');
    elseif rw4scale ~=0 
        rw4scale = .1;
        quiver3(rwx4, rwy4, rwz4, g(1, 4) * rw4scale,g(2,4) * rw4scale,  g(3,4) * rw4scale,  'k', 'LineWidth', 1,'DisplayName', 'RW 4');
        fprintf('rw 4, initial speed only\n');
    end

end


function graphRCOM(t, rcom)
    figure;
    title('COM Translation', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'g');
    xlabel('Time - Seconds', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Center of Mass Poision (x,y,z)', 'FontSize', 12, 'FontWeight', 'bold');

    hold on;
    plot(t, rcom(:, 1), '-b', 'DisplayName', 'component (x)');
    plot(t, rcom(:, 2), '-m', 'DisplayName', 'component (y)'); % r 
    plot(t, rcom(:, 3), '-c', 'DisplayName', 'component (z)');
    hold off;
end 


function arcs = rad2arc(rad)
    arcs = rad * 60 ^2 * 360 / (2*pi) ; 
end