% MOD 8 : add theta due to w frame converging 

load("parameters2.mat")
%%FOR GRAPHING PURPOSES%%
%set sizes
dt = 0.01; 
tspan = 0:dt:10 ;
%%initial conditions
w_b_n_0 = [0 0 0]'; 
rd_0 = [0 0 0]';
theta_0 = 0;
W_0 = [W_init theta_0 0 ]';
w2_0 = W_frame_init(:,2);
w3_0 = W_frame_init(:,3);
theta_0 = 0;
y0= [ w_b_n_0 , rd_0, W_0, w2_0, w3_0];

yin = y0 ;

[dim1,dim2] = size(tspan);
store_W       = zeros(3, dim2);
store_R       = zeros(3, dim2);

store_omega   = zeros(2, dim2);
store_Wframe2 = zeros(3, dim2);
store_Wframe3 = zeros(3, dim2);

store_forcedw2 = zeros(3, dim2);
store_forcedw3 = zeros(3, dim2);


%populate with initial conditions
store_W(:,1)        = y0(:,1);
store_R(:,1)        = y0(:,2);
store_omega(1)      = y0(1,3);
store_omega(2)      = y0(2,3);
%store_Wframe2(:,1)  = y0(:,4);
%store_Wframe3(:,1)  = y0(:,5);

store_forcedw2(:,1) = y0(:,4);
store_forcedw3(:,1) = y0(:,5);

%%%%%%%%%%%%%%%%%%%%%%%%%%





for i=1:tspan(end)/dt
    time = i*dt;
    yout = rk4singlestep(@(t,y)func(t,y),dt,time,yin);

    store_W(:,i+1)       = yout(:,1);
    store_R(:,i+1)       = yout(:,2);
    store_omega(1,i+1)   = yout(1,3); %not a vector
    theta                = mod(yout(2,3),360);
    store_omega(2,i+1)   = theta;
    %store_Wframe2(:,i+1) = yout(:,4);
    %store_Wframe3(:,i+1) = yout(:,5);

    store_forcedw2(:,i+1) = rotateg(w2_0, theta);
    store_forcedw3(:,i+1) = rotateg(w3_0, theta);


    yin = yout;
end

%%%%%%%%GRAPHS%%%%%%%%
%figure;
%plot(tspan, store_W(1, :), '-g', 'DisplayName', 'Component 1');
%hold on;
%title('Space Craft Angular Velocity', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'g');
%plot(tspan, store_W(2, :), '-g', 'DisplayName', 'Component 2');
%plot(tspan, store_W(3, :), '-g', 'DisplayName', 'Component 3');
%hold off;
%figure;
%plot(tspan, store_R(1, :), '-b', 'DisplayName', 'Component 1');
%hold on;
%title('Inertial Translational Velocity of sc', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'g');
%plot(tspan, store_R(2, :), '-b', 'DisplayName', 'Component 2');
%plot(tspan, store_R(3, :), '-b', 'DisplayName', 'Component 3');
%hold off;
%figure;
%hold on;
%title('Reaction Wheel angular velocity', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'r');
%plot(tspan, store_omega, '-r', 'DisplayName', 'Component 2');
%hold off;
figure;
title('(W2) Motor frame components', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'g');
hold on;
%plot(tspan, store_Wframe2(1, :), 'Color', [0.3, 0.75, 0.93], 'DisplayName', 'Component 2');
%plot(tspan, store_Wframe2(2, :), 'Color', [0.47, 0.67, 0.19], 'DisplayName', 'Component 2');
%plot(tspan, store_Wframe2(3, :), 'Color', [1, 0.71, 0.76], 'DisplayName', 'Component 2');

plot(tspan, store_forcedw2(1, :), 'Color', [0, 0, 0.55], 'DisplayName', 'Component 2');
plot(tspan, store_forcedw2(3, :), 'Color', [0, 0.5, 0], 'DisplayName', 'Component 2');
plot(tspan, store_forcedw2(2, :), 'Color', [0.55, 0, 0.55], 'DisplayName', 'Component 2');
hold off;

figure;
title('(W3) Motor frame components', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'g');
hold on;
%plot(tspan, store_Wframe3(1, :), 'Color', [0.3, 0.75, 0.93], 'DisplayName', 'Component 2');
%plot(tspan, store_Wframe3(2, :), 'Color', [0.47, 0.67, 0.19], 'DisplayName', 'Component 2');
%plot(tspan, store_Wframe3(3, :), 'Color', [1, 0.71, 0.76], 'DisplayName', 'Component 2');

plot(tspan, store_forcedw3(1, :), 'Color', [0, 0, 0.55], 'DisplayName', 'Component 2');
plot(tspan, store_forcedw3(3, :), 'Color', [0, 0.5, 0], 'DisplayName', 'Component 2');
plot(tspan, store_forcedw3(2, :), 'Color', [0.55, 0, 0.55], 'DisplayName', 'Component 2');
hold off;
%%%%%%%%%%%%%%%%%%%%%%


function xout = rk4singlestep(f, dt, tk, xk)
    load('parameters.mat')
    f1 = f(tk,xk) ;
    f2 = f(tk + dt/2, xk + (dt/2)*f1);
    f3 = f(tk + dt/2, xk + (dt/2)*f2);
    f4 = f(tk+dt, xk + dt*f3);

    xout = xk + (dt/6) * (f1 + 2*f2 + 2*f3 + f4);
    

end



function dx = func(t,y)
    
    load('parameters2.mat')
    w_b_n   = y(:,1);
    r_b_nd  = y(:,2);
    W       = y(1,3);
    theta   = mod( y(2,3), 360);
    wframe_2 = rotateg(W_frame_init(:,2), theta);
    wframe_3 = rotateg(W_frame_init(:,3), theta);
    W_frame = [g wframe_2 wframe_3];

   % W_frame = W_frame ./ vecnorm(W_frame); %FORCE NORMALISE
    %[g y(:,4) y(:,5)] - W_frame
    
    r_wc_b       = r_w_b + d * W_frame(:,2) ; %(6)
    r_wc_b_dash  = d * W * W_frame(:,3)      ; %(7)
    

    ihub_b = ihub_bc + mhub * skew(r_bc_b) * skew(r_bc_b)'; %(18) P1 - parallel axis thm 
    irw_b  = irw_wc  + mrw  * skew(r_wc_b) * skew(r_wc_b)'; %(26)
    isc_b  = ihub_b  + irw_b;  

    irw_wc_dash = (J12 *  W_frame(:,3) - J13*W_frame(:,2)) * W * W_frame(:,1)'...
        + (-J13 *  W_frame(:,1) + J22 * W_frame(:,3) -2 * J23 * W_frame(:,2)...
        - J33 & W_frame(:,3)) * W * W_frame(:,2)'...
        + (J12 * W_frame(:,1) + J22 * W_frame(:,2) + 2* J23 * W_frame(:,3)...
        - J33 * W_frame(:,2)) * W * W_frame(:,3)';

    ihub_b_dash = zeros(3,3);
    irw_b_dash  = irw_wc_dash + mrw * skew(r_wc_b_dash)*skew(r_wc_b)' ...
        + mrw * skew(r_wc_b)*skew( r_wc_b_dash)';                   %(35)
    isc_b_dash  = ihub_b_dash + irw_b_dash ;
    
    [c,cdash] = getcs(r_wc_b, W, r_wc_b_dash);

    e = inv([J11 + mrw * d^2 - mrw^2 * d^2/msc]);       % E = INV(A), (71)


    F = - ((J11 + mrw * d^2) * W_frame(:,1)' + J12*W_frame(:,2)' + ...
        J13 * W_frame(:,3)'+ (mrw*d*W_frame(:,3)') * (skew(c) - ....
        skew(r_w_b)));                                              %(72)

    i_lhs = isc_b + msc * skew(c)*skew(c) +( irw_wc*W_frame(:,1) ...
        +mrw * d * (skew(r_wc_b) - skew(c))*W_frame(:,3))*e*F;      %(80)

    [w_1, w_2, w_3] = ws(w_b_n, W_frame);
    
    v = -mrw*d*W_frame(:,3)' * (r_c_ndd - 2 * skew(w_b_n) * cdash - ...
        skew(w_b_n) * skew(w_b_n)* c + (mrw/msc)*d*W^2*W_frame(:,2))...
        +J23*(w_3^2 - w_2^2) + w_1*(J12 * w_3 - J13 * w_2) ...
        + w_2 * w_3 * ( J22 - J33 - mrw * d^2) - mrw * d * ...
        W_frame(:,3)' * skew(w_b_n)*skew(w_b_n) * r_w_b + us;       %(73)

    Lb = W^2 * (U_d+U_s);%WHERE IS THIS FROM PLS 

    tau_rhs = - skew(w_b_n) * isc_b * w_b_n - isc_b_dash * w_b_n ...
        -msc * skew(c)*(r_c_ndd - 2 * skew(w_b_n) * cdash - ...
        skew(w_b_n)*skew(w_b_n)*c) + mrw * d * W^2 * skew(r_wc_b) ...
        * W_frame(:,2) - irw_wc_dash * W * W_frame(:,1) - ...
        skew(w_b_n) * (irw_wc * W * W_frame(:,1) + mrw * skew(r_wc_b)*...
        r_wc_b_dash) - mrw * d* W^2 * skew(c) * W_frame(:,2) - ...
        (irw_wc*W_frame(:,1) + mrw * d * (skew(r_wc_b) - skew(c)) ...
        *W_frame(:,3))* e * v +Lb ;                                 %(81)

    w_b_nd = inv(i_lhs)*tau_rhs ;%(79) 
    %w_b_nd = i_lhs \ tau_rhs 
    Wd = e * F * w_b_nd + e * v ;%(75)
    
    r_wc_b_ddash = d * Wd * W_frame(:,3) - d * W^2 * W_frame(:,2) ; %(8)
    cddash       = (mrw/msc) * r_wc_b_ddash ; %(5)
    cdd          = cddash + 2 * cross(w_b_n, cdash) + cross(w_b_nd, c) + ....
        cross(w_b_n, cross(w_b_n, c)) ;

    r_b_ndd = - cdd ; 
    
    w_frame_2d = W * W_frame(:,3);
    w_frame_3d = -W * W_frame(:,2);

    thetad = W ;
    x1 = w_b_nd;
    x2 = r_b_ndd;
    x3 = [Wd thetad 0]';
    x4 = w_frame_2d;
    x5 = w_frame_3d;

    dx = [x1 x2 x3 x4 x5];

end

function S = skew(v)
    S = [  0   -v(3)  v(2); v(3)   0   -v(1); -v(2)  v(1)   0  ];
end

function [c,cdash] = getcs(r_wc_b, curr_W, r_wc_b_dash)
    load('parameters2.mat')
    c      = (1/msc) * ( mhub * r_bc_b + mrw * r_wc_b);  %(3)
    cdash  = (1/msc) * ( mrw * r_wc_b_dash);            %(4)
    %cddash = (1/msc) * ( mrw * r_wc_b_ddash);
    
end

function [w_1, w_2, w_3] = ws(w_b_n, W_frame)
    w_1 = W_frame(:,1)' * w_b_n ;
    w_2 = W_frame(:,2)' * w_b_n ;
    w_3 = W_frame(:,3)' * w_b_n ; 
    %(45)-(47)
end

function v_rot = rotateg(v, theta)
    load('parameters2.mat')
    g = g / norm(g);
    ang = deg2rad(theta);

    v_rot = dot(v, g) * g + cos(theta) * (v - dot(v, g) * g ) + sin(theta) * cross(g, v);
    
    %R = axang2rotm([g' theta]) ;
    %v_rotated = R * v;
    %v_rot - v_rotated
end