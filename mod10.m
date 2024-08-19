%mod10 - implement ODE 45 for dynamic time step
% add pointing error 
load("parameters2.mat")
%set sizes
dt = 0.01; 
tspan = 0:dt:5;
%%initial conditions
w_b_n_0 = [0 0 0]'; 
rd_0 = [0 0 0]';
theta_0 = 0;
W_0 = [W_init theta_0 0 ]';

w2_0 = W_frame_init(:,2);
w3_0 = W_frame_init(:,3);
global us
us = NaN ;

y0 = [ w_b_n_0 ; rd_0; [W_init theta_0 ]'];
[t,y] = ode45(@func,tspan,y0);


%%%%%%%%GRAPHS%%%%%%%%
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
plot(t, y(:, 7), '-r', 'DisplayName', 'RW Omega');%%omega
plot(t, y(:,8), '-b', 'DisplayName', 'Theta');%%theta
hold off;

%%%%%%%%%%%%%%%%%%%%%%


function xout = rk4singlestep(f, dt, tk, xk)
    load('parameters2.mat')
    f1 = f(tk,xk) ;
    f2 = f(tk + dt/2, xk + (dt/2)*f1);
    f3 = f(tk + dt/2, xk + (dt/2)*f2);
    f4 = f(tk+dt, xk + dt*f3);

    xout = xk + (dt/6) * (f1 + 2*f2 + 2*f3 + f4);
end


function dx = func(t,y)
    global us
    load('parameters2.mat')
    if isnan(us)
        us = us_0;
        fprintf('US reset to US_0 %f\n', us);
    end

    y = [y(1:3) y(4:6) [y(7:8) ;0]];

    w_b_n   = y(:,1);
    r_b_nd  = y(:,2);
    W       = y(1,3);
    theta   = mod(y(2,3), 360);
    wframe_2 = rotateg(W_frame_init(:,2), theta);
    wframe_3 = rotateg(W_frame_init(:,3), theta);
    
    r_wc_b       = r_w_b + d * wframe_2 ; %(6)
    r_wc_b_dash  = d * W * wframe_3      ; %(7)

    ihub_b = ihub_bc + mhub * skew(r_bc_b) * skew(r_bc_b)'; %(18) P1 - parallel axis thm 
    irw_b  = irw_wc  + mrw  * skew(r_wc_b) * skew(r_wc_b)'; %(26)
    isc_b  = ihub_b  + irw_b;  

    irw_wc_dash = (J12 *  wframe_3 - J13*wframe_2) * W * g'...
        + (-J13 *  g + J22 * wframe_3 -2 * J23 * wframe_2...
        - J33 & wframe_3) * W * wframe_2'...
        + (J12 * g + J22 * wframe_2 + 2* J23 * wframe_3...
        - J33 * wframe_2) * W * wframe_3';

    ihub_b_dash = zeros(3,3);
    irw_b_dash  = irw_wc_dash + mrw * skew(r_wc_b_dash)*skew(r_wc_b)' ...
        + mrw * skew(r_wc_b)*skew( r_wc_b_dash)';                   %(35)
    isc_b_dash  = ihub_b_dash + irw_b_dash ;
    
    [c,cdash] = getcs(r_wc_b, W, r_wc_b_dash);

    e = inv([J11 + mrw * d^2 - mrw^2 * d^2/msc]);       % E = INV(A), (71)


    F = - ((J11 + mrw * d^2) * g' + J12*wframe_2' + ...
        J13 * wframe_3'+ (mrw*d* wframe_3') * (skew(c) - ....
        skew(r_w_b)));                                              %(72)

    i_lhs = isc_b + msc * skew(c)*skew(c) +( irw_wc*g ...
        +mrw * d * (skew(r_wc_b) - skew(c))*wframe_3)*e*F;      %(80)

    %wframe = [g ,wfame_2, wframe_3]
    [w_1, w_2, w_3] = ws(w_b_n, [g ,wframe_2, wframe_3]);
    
    v = -mrw*d*wframe_3' * (r_c_ndd - 2 * skew(w_b_n) * cdash - ...
        skew(w_b_n) * skew(w_b_n)* c + (mrw/msc)*d*W^2*wframe_2)...
        +J23*(w_3^2 - w_2^2) + w_1*(J12 * w_3 - J13 * w_2) ...
        + w_2 * w_3 * ( J22 - J33 - mrw * d^2) - mrw * d * ...
        wframe_3' * skew(w_b_n)*skew(w_b_n) * r_w_b + us;       %(73)

    Lb = W^2 * (U_d+U_s);%WHERE IS THIS FROM PLS 

    tau_rhs = - skew(w_b_n) * isc_b * w_b_n - isc_b_dash * w_b_n ...
        -msc * skew(c)*(r_c_ndd - 2 * skew(w_b_n) * cdash - ...
        skew(w_b_n)*skew(w_b_n)*c) + mrw * d * W^2 * skew(r_wc_b) ...
        * wframe_2 - irw_wc_dash * W * g - ...
        skew(w_b_n) * (irw_wc * W * g + mrw * skew(r_wc_b)*...
        r_wc_b_dash) - mrw * d* W^2 * skew(c) * wframe_2 - ...
        (irw_wc*g + mrw * d * (skew(r_wc_b) - skew(c)) ...
        *wframe_3)* e * v +Lb ;                                 %(81)

    w_b_nd = inv(i_lhs)*tau_rhs ;%(79) 
    %w_b_nd = i_lhs \ tau_rhs 
    Wd = e * F * w_b_nd + e * v ;%(75)
    
    r_wc_b_ddash = d * Wd * wframe_3 - d * W^2 * wframe_2 ; %(8)
    cddash       = (mrw/msc) * r_wc_b_ddash ; %(5)
    cdd          = cddash + 2 * cross(w_b_n, cdash) + cross(w_b_nd, c) + ....
        cross(w_b_n, cross(w_b_n, c)) ;

    r_b_ndd = - cdd ; 
    

    us = J11*(g' * w_b_nd + Wd);
    

    thetad = W ;

    x1 = w_b_nd;
    x2 = r_b_ndd;
    x3 = [Wd thetad]';


    dx = [x1; x2; x3];

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

function v_rot = rotateg(v, ang)
    load('parameters2.mat')
    g = g / norm(g);
    
    theta = deg2rad(ang);

    v_rot = dot(v, g) * g + cos(theta) * (v - dot(v, g) * g ) + sin(theta) * cross(g, v);

end
