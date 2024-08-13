%%MODEL 5 : remove global variable, 

load('parameters.mat')
%%%%%%%%%%%%%%%%%%%%%%
%%only w2 and w3 chage  of W_Frame

%%W frame is the non inertial frame of the rw's spinning part, consisting
%%of unit vectors relative to body fixed frame of space craft 

W_frame = W_frame_init;

dt = 0.01;
tspan = 0:dt:5 ;
w_b_n_0 = [0 0 0]'; 
rd_0 = [ 0 0 0]';
W_0 = [W_init 0 0 ]';
w2_0 = W_frame(:,2);
w3_0 = W_frame(:,3);
y0= [ w_b_n_0 , rd_0, W_0, w2_0, w3_0];

%%FOR GRAPHING PURPOSES%%
%set sizes
[dim1,dim2] = size(tspan);
store_W       = zeros(3, dim2);
store_R       = zeros(3, dim2);
store_omega   = zeros(1, dim2);
store_Wframe2 = zeros(3, dim2);
store_Wframe3 = zeros(3, dim2);

%populate with initial conditions
store_W(:,1)        = y0(:,1);
store_R(:,1)        = y0(:,2);
store_omega(1)      = y0(1,3);
store_Wframe2(:,1)  = y0(:,4);
store_Wframe3(:,1)  = y0(:,5);
%%%%%%%%%%%%%%%%%%%%%%%%%%

yin = y0;

for i=1:tspan(end)/dt
    time = i*dt;
    yout = rk4singlestep(@(t,y)func(t,y,W_frame),dt,time,yin);

    %for graphing
    store_W(:,i+1)       = yout(:,1);
    store_R(:,i+1)       = yout(:,2);
    store_omega(i+1)     = yout(1,3); %not a vector
    store_Wframe2(:,i+1) = yout(:,4);
    store_Wframe3(:,i+1) = yout(:,5);

    yin = yout;
end

%%%%%%%%GRAPHS%%%%%%%%
figure;
plot(tspan, store_W(1, :), '-g', 'DisplayName', 'Component 1');
hold on;
title('Space Craft Angular Velocity', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'g');
plot(tspan, store_W(2, :), '-g', 'DisplayName', 'Component 2');
plot(tspan, store_W(3, :), '-g', 'DisplayName', 'Component 3');
hold off;
figure;
plot(tspan, store_R(1, :), '-b', 'DisplayName', 'Component 1');
hold on;
title('Inertial Translational Velocity of sc', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'g');
plot(tspan, store_R(2, :), '-b', 'DisplayName', 'Component 2');
plot(tspan, store_R(3, :), '-b', 'DisplayName', 'Component 3');
hold off;
figure;
hold on;
title('Reaction Wheel angular velocity', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'r');
plot(tspan, store_omega, '-r', 'DisplayName', 'Component 2');
hold off;
figure;
plot(tspan, store_Wframe3(1, :), '-b', 'DisplayName', 'Component 2');
hold on;
title('(W) Motor frame components', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'g');
plot(tspan, store_Wframe3(2, :), '-b', 'DisplayName', 'Component 2');
plot(tspan, store_Wframe3(3, :), '-b', 'DisplayName', 'Component 2');
plot(tspan, store_Wframe2(1, :), '-g', 'DisplayName', 'Component 2');
plot(tspan, store_Wframe2(2, :), '-g', 'DisplayName', 'Component 2');
plot(tspan, store_Wframe2(3, :), '-g', 'DisplayName', 'Component 2');
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

function dx = func(t,y, W_frame)
    load('parameters.mat')

    % y = [ w_b_n   , r_b_nd , W, wframe 2, wframe 3]
    %dy = [ w_b_n_d , rdd     , Wd , wframe2d , wframe3 d  ]  
    
    %unpack states : 
    w_b_n        = y(:,1);
    r_b_nd       = y(:,2);
    W            = y(1,3);
    W_frame(:,2) = y(:,4);
    W_frame(:,3) = y(:,5);
    
    F = W^2 * U_s;
    Lb = W^2 * (U_d+U_s);

    r_wc_b      =   r_w_b + d * W_frame(:,2)' ; %(6)
    r_wc_b_dash =   d * W * W_frame(:,3)      ; %(7)

    [c,cdash]       = getcs(r_wc_b, W, r_wc_b_dash);
    [w_1, w_2, w_3] = ws(w_b_n, W_frame)       ; %pass w_b_n current to get W frame omega components
    
    ihub_b = ihub_bc + mhub * skew(r_bc_b) * skew(r_bc_b)'; %(18) P1
    irw_b  = irw_wc  + mrw  * skew(r_wc_b) * skew(r_wc_b)'; %(26)
    isc_b  = ihub_b  + irw_b;                               % total I in B frame is sum of I of hub + RWs

    %build A first(>1 RW)
    irw_wc_dash = [0, -J13 ,J12 ; -J13, -2*J23 , J22 - J33; J12, J22 - J33, 2*J23] * W; %(23)

    e = inv([J11 + mrw * d^2 - mrw^2 * d^2/msc]);           % E = INV(A), (71)

    v = -mrw*d*W_frame(:,3)' * (r_c_ndd - 2 * skew(w_b_n) * cdash' - ...
        skew(w_b_n) * skew(w_b_n)* c' + (mrw/msc)*d*W^2*W_frame(:,2))...
        +J23*(w_3^2 - w_2^2) + w_1*(J12 * w_3 - J13 * w_2) ...
        + w_2 * w_3 * ( J22 - J33 - mrw * d^2) - ...
        mrw * d * W_frame(:,3)' * skew(w_b_n)*skew(w_b_n) * r_w_b' + us; %(73)

    tau_rhs = (- skew(w_b_n) * isc_b - isc_b') * w_b_n - msc * skew(c) * ...
        (r_c_ndd - 2 * skew(w_b_n) * cdash' - skew(w_b_n)*skew(w_b_n)* c') ...
        + mrw * d * W^2 * skew(r_wc_b)*W_frame(:,2) - irw_wc_dash * W * W_frame(:,1) ...
        -skew(w_b_n) * ( irw_wc * W * W_frame(:,1) + mrw * skew(r_wc_b)*r_wc_b_dash)...
        - mrw*d*W^2 * skew(c) * W_frame(:,2) ...
        - (irw_wc*W_frame(:,1) + mrw * d * (r_wc_b - skew(c)) * W_frame(:,3))* e' * v +Lb; %(81)

    f = - ((J11 + mrw * d^2) * W_frame(:,1)' + J12*W_frame(:,2)' + J13 * W_frame(:,3)'....
        + (mrw*d*W_frame(:,3)') * (skew(c) - skew(r_w_b))); %(72)

    ilhs = isc_b + msc * skew(c) * skew(c) + (irw_wc * W_frame(:,1) + mrw * d * (skew(r_wc_b)...
        - skew(c))*W_frame(:,3)) * e' * f; %(80)

    w_b_nd = tau_rhs \ ilhs ; %(79)

    Wd = e * f * w_b_nd' + e * v  ; %(74)

    x1 = (tau_rhs \ ilhs)';          % 79  (1 x 3)
    x2 = r_c_ndd - 2 * skew(w_b_n) * cdash' - ...
        skew(w_b_n)*skew(w_b_n)*c' + ...
        (mrw * d*W^2/msc) * W_frame(:,2)...
        +skew(c)*w_b_nd' - ...
        (mrw*d/msc)*W_frame(:,3)*Wd; %(10)
    x3 = [Wd 0 0]';              %(74)
    x4 =  W * W_frame(:,3) ;         %(27)
    x5 = -W * W_frame(:,2) ;         %(28)

    dx= [x1 ,x2 , x3 , x4, x5];
    % w_b_n d , r dd, W d

end

function [c,cdash] = getcs(r_wc_b, curr_W, r_wc_b_dash)
    load('parameters.mat')
    c      = (1/msc) * ( mhub * r_bc_b + mrw * r_wc_b); %(3)
    cdash  = (1/msc) * ( mrw * r_wc_b_dash');            %(4)
end

function S = skew(v)
    S = [  0   -v(3)  v(2); v(3)   0   -v(1); -v(2)  v(1)   0  ];
end


function [w_1, w_2, w_3] = ws(w_b_n, W_frame)
    w_1 = W_frame(:,1)' * w_b_n ;
    w_2 = W_frame(:,2)' * w_b_n;
    w_3 = W_frame(:,3)' * w_b_n; 
    %(45)-(47)
end