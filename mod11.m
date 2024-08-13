% mod 11 : added more rws, still need to add JS and IrwWc
load("parameters3.mat")
%set sizes
dt = 0.01; 
tspan = 0:dt:10;

global us
us = [NaN, NaN, NaN, NaN] ;
global angM
global Hrwlog
global Hhublog

Hhublog = [];
Hrwlog  = [];
angM    = [];



%%initial conditions
w_b_n_0 = [0 0 0]'; 
rd_0 = [0 0 0]';

W_0 = W_init;
theta_0 = [t1 t2 t3 t4]'; %reactionwheel angle of rotation in b frame 

%%%%%%%%%SOLVE%%%%%%%%%%
y0 = [ w_b_n_0 ; rd_0; W_0 ; theta_0 ];
[t,y] = ode45(@func,tspan,y0);
%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%POINTING%%%%%%%%%
[dim1,dim2] = size(t);
angles = zeros(1, dim1);

for i = 1:dim1-1
    time = t(i+1) - t(i);
    angles(i) = time * norm(y(i, 1:3));
end
%%%%%%%%%%%%%%%%%%%%%%%%

graph(t, y, angles');

%%%%%%%%GRAPHS%%%%%%%%
%figure;
%hold on;
%title('H rw', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'g');
%plot(Hrwlog)
%hold off;
%figure;
%hold on;
%title('H hub', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'g');
%plot(Hhublog)
%hold off;
%figure;
%hold on;
%title('Angular Momentum', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'g');
%plot(angM)
%hold off;
function graph(t,y, angles)
figure;
hold on;
title('Pointing Error', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'g');
plot(t, angles);
hold off;

figure;
title('Space Craft Angular Velocity', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'g');
hold on;
plot(t, y(:, 1), '-g', 'DisplayName', 'Component 1');
plot(t, y(:, 2), '-g', 'DisplayName', 'Component 1');
plot(t, y(:, 3 ), '-g', 'DisplayName', 'Component 1');
hold off;

figure;
title('Inertial Translational Velocity of sc', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'g');
hold on;
plot(t, y(:, 4), '-b', 'DisplayName', 'Component 1');
plot(t, y(:, 5), '-b', 'DisplayName', 'Component 2');
plot(t, y(:, 6), '-b', 'DisplayName', 'Component 3');
hold off;

figure;
hold on;
title('Reaction Wheel angular velocity', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'r');
plot(t, y(:, 7), '-r', 'DisplayName', 'RW 1');%%omega
plot(t, y(:,8), '-b', 'DisplayName', 'RW 2');%%theta
plot(t, y(:,9), '-g', 'DisplayName', 'RW 3');%%theta
plot(t, y(:,10), '-y', 'DisplayName', 'RW 4');%%theta

hold off;
end
%%%%%%%%%%%%%%%%%%%%%%


function dx = func(t,y)
    global us
    global Hrwlog
    global Hhublog
    global angM
    load('parameters3.mat')
    if ismember(1, isnan(us))
        us = us_0;
        fprintf('US reset to US_0 %f\n', us);
        %should only print once per each reaction wheel 
    end
    
    %         3         3   1x4     1x4
    %y0 = [ w_b_n_0 ; rd_0; W_0 ; theta_0 ];


    w_b_n   = y(1:3);
    r_b_nd  = y(4:6);
    W = y(7:10);
    thetas = y(11:14); %%angle of displacement for each reaction wheel

   %find the orientation of the reaction wheels by rotating their initial
   %direction vectors by their angle of displacement
   % wframe_2 = [w2-rw1; w2-rw2; w2-rw3; w2-rw4] etc. 
    wframe_2 = [rotateg(W_frame_init1(:,2), thetas(1),  g(:, 1)) , ...% rotates w2 about g by theta all for rw 1 
        rotateg(W_frame_init2(:,2), thetas(2), g(:, 2)) , ...
        rotateg(W_frame_init3(:,2), thetas(3), g(:, 3)), ...
        rotateg(W_frame_init4(:,2), thetas(4), g(:, 4))];
    wframe_3 = [rotateg(W_frame_init1(:,3), thetas(1), g(:, 1)) , ...
        rotateg(W_frame_init2(:,3), thetas(2), g(:, 2)) ,...
        rotateg(W_frame_init3(:,3), thetas(3), g(:, 3)) , ...
        rotateg(W_frame_init4(:,3), thetas(4), g(:, 4)) ];
    
    %r_wc_b = [r_wc_b for rw1; r_wc_b for rw2; etc]
    r_wc_b = [ r_w_b(:,1) + d(1) * wframe_2(:, 1)  , ...
        r_w_b(:,2) + d(2) * wframe_2(:, 2) , ...
        r_w_b(:,3) + d(3) * wframe_2(:, 3) , ...
        r_w_b(:,4) + d(4) * wframe_2(:, 4) ] ;  %(6)
    %non interial, b frame derivative of r_wc_b
    r_wc_b_dash  = [d(1) * W(1) * wframe_3(:,1) ,...
        d(2) * W(2) * wframe_3(:,2) , ...
        d(3) * W(3) * wframe_3(:,3) , ...
        d(4) * W(4) * wframe_3(:,4) ] ;    %(7)

    ihub_b = ihub_bc + mhub * skew(r_bc_b) * skew(r_bc_b)'; %(18) p1 ir parallel axis thm 
%%irw_wc all the same
    irw_b = [irw_wc  + mrw(1)  * skew(r_wc_b(:,1)) * skew(r_wc_b(:,1))' , ...
             irw_wc  + mrw(2)  * skew(r_wc_b(:,2)) * skew(r_wc_b(:,2))' , ...
             irw_wc  + mrw(3)  * skew(r_wc_b(:,3)) * skew(r_wc_b(:,3))' , ...
             irw_wc  + mrw(4)  * skew(r_wc_b(:,4)) * skew(r_wc_b(:,4))' ]; %(26) p1 


    isc_b  = ihub_b  + irw_b(1:3,1:3) + irw_b(1:3,4:6) + irw_b(1:3,7:9) + irw_b(1:3,10:12);

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


    ihub_b_dash = zeros(3,3);
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



    isc_b_dash  = ihub_b_dash + irw_b_dash1+ irw_b_dash2 + irw_b_dash3+ irw_b_dash4;
    
    [c,cdash] = getcs(r_wc_b, r_wc_b_dash);
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
    
    %%FIX????? WHY IS A DIAGONAL PLS 
    e = inv(A);% E = INV(A), (71)

    
    F = zeros(4,3);
    for i = 1:4
        F(i,:) = [- ((J11 + mrw(i) * d(i)^2) * g(:,i)' + J12*wframe_2(:,i)' + ...
        J13 * wframe_3(:,i)'+ (mrw(i)*d(i)* wframe_3(:,i)') * (skew(c) - ....
        skew(r_w_b(:,i))))]'; %(72)
    end
    
    i_lhs = isc_b + msc * skew(c)*skew(c); %(80)
    for i = 1:4
        i_lhs = i_lhs +( irw_wc*g(:,i) ...
        +mrw(i) * d(i) * (skew(r_wc_b(:,i)) ...
        - skew(c))*wframe_3(:,i))*e(i,:)*F;      %(80)
    end
    
    %wframe = [g ,wfame_2, wframe_3]
    [w_1, w_2, w_3] = ws(w_b_n, g ,wframe_2, wframe_3);
    
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

    
    Lb = W(1)^2 * (U_d * W_frame_init1(:,2) +U_s*skew(r_w_b(:,1))*wframe_2(:,1))+ ...
       W(2)^2 * (U_d * W_frame_init2(:,2) +U_s*skew(r_w_b(:,2))*wframe_2(:,2))+...
       W(3)^2 * (U_d * W_frame_init3(:,2) +U_s*skew(r_w_b(:,3))*wframe_2(:,3))+...
       W(4)^2 * (U_d * W_frame_init4(:,2) +U_s*skew(r_w_b(:,4))*wframe_2(:,4)); %(83) (84) ????????? 
 
    tau_rhs = - skew(w_b_n) * isc_b * w_b_n - isc_b_dash * w_b_n ...
        -msc * skew(c)*(r_c_ndd - 2 * skew(w_b_n) * cdash - ...
        skew(w_b_n)*skew(w_b_n)*c) + Lb;
    
    irw_wc_dash = [irw_wc_dash1 irw_wc_dash2 irw_wc_dash3 irw_wc_dash4];


    for i = 1:4
        tau_rhs = tau_rhs + mrw(i) * d(i) * W(i)^2 * skew(r_wc_b(:,i)) ...
            * wframe_2(:,i) -  irw_wc_dash(:,i*3 - 2 :i*3 ) * W(i) * g(:,i) - ...
            skew(w_b_n) * (irw_wc * W(i) * g(:,i) + mrw(i) * skew(r_wc_b(:,i))*...
            r_wc_b_dash(:,i)) - mrw(i) * d(i)* W(i)^2 * skew(c) * wframe_2(:,i) - ...
            (irw_wc*g(:,i) + mrw(i) * d(i) * (skew(r_wc_b(:,i)) - skew(c)) ...
            *wframe_3(:,i))* e(i,:) * v  ;                              %(81)
    end 

    w_b_nd = inv(i_lhs)*tau_rhs ;%(79) 
    %w_b_nd = i_lhs \ tau_rhs ;
    Wd = e * F * w_b_nd + e * v ;%(75)


    r_wc_b_ddash = [d(1) * Wd(1) * wframe_3(:,1) - d(1) * W(1)^2 * wframe_2(:,1) ,...
        d(2) * Wd(2) * wframe_3(:,2) - d(2) * W(2)^2 * wframe_2(:,2) ,...
        d(3) * Wd(3) * wframe_3(:,3) - d(3) * W(3)^2 * wframe_2(:,3) ,...
        d(4) * Wd(4) * wframe_3(:,4) - d(4) * W(4)^2 * wframe_2(:,4) ]; %(8)
    
    cddash  = (1/msc) *(mrw(1)* r_wc_b_ddash(:,1) +...
        mrw(2) * r_wc_b_ddash(:,2) +...
        mrw(3) * r_wc_b_ddash(:,3) +...
        mrw(4) * r_wc_b_ddash(:,4)) ; %(5)

    cdd     = cddash + 2 * cross(w_b_n, cdash) +...
        cross(w_b_nd, c) + cross(w_b_n, cross(w_b_n, c)) ;%(9)

    r_b_ndd = r_c_ndd- cdd ; 
    
    %here
    us = J11*(g' * w_b_nd + Wd);
    
    thetad = W ;
    
    %Pointing error - not sure about any of this at all
    

    x1 = w_b_nd   ;
    x2 = r_b_ndd  ;
    x3 = Wd       ;
    x4 = thetad   ;

    
    dx = [x1; x2; x3 ; x4] ;
    %         3         3   1x4     1x4
    %y0 = [ w_b_n_0 ; rd_0; W_0 ; theta_0 ];
end

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

function v_rot = rotateg(v, ang, g)
    g = g / norm(g);
    
    theta = deg2rad(ang);

    v_rot = dot(v, g) * g + cos(theta) * (v - dot(v, g) * g ) + sin(theta) * cross(g, v);

end
