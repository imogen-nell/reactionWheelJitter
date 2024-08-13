global W_frame
load('parameters.mat')
%%%%%%%%%%%%%%%%%%%%%%
%%only w2 and w3 chage 

W_frame = W_frame_init;

dt = 0.01;
tspan = 0:dt:1 ;
w_b_n_0 = [0 0 0]';
rd_0 = [ 0 0 0]';
W_0 = [W 0 0 ]';
y0= [ w_b_n_0 , rd_0, W_0];

[dim1,dim2] = size(tspan);
store_W=zeros(3, dim2);
store_R=zeros(3, dim2);
store_omega=zeros(1, dim2);
store_W(:,1) = y0(:,1);
store_R(:,1) = y0(:,2);
store_omega(1)= y0(1,3);
yin = y0;

for i=1:tspan(end)/dt
    time = i*dt;
    yout = rk4singlestep(@(t,y)func(t,y),dt,time,yin);
    store_W(:,i+1) = yout(:,1);
    store_R(:,i+1) =  yout(:,2);
    store_omega(i+1)= yout(1,3);
    yin = yout;
end

figure;
plot(tspan, store_W(1, :), '-g', 'DisplayName', 'Component 1');
hold on;
plot(tspan, store_W(2, :), '-g', 'DisplayName', 'Component 2');
plot(tspan, store_W(3, :), '-g', 'DisplayName', 'Component 3');
hold off;
figure;
plot(tspan, store_R(1, :), '-b', 'DisplayName', 'Component 1');
hold on;
plot(tspan, store_R(2, :), '-b', 'DisplayName', 'Component 2');
plot(tspan, store_R(3, :), '-b', 'DisplayName', 'Component 3');
hold off;
figure;
hold on;
plot(tspan, store_omega, '-r', 'DisplayName', 'Component 2');
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
    global W_frame
    load('parameters.mat')
% y = [ w_b_n , r_b_nd , W]
%dy =  [ wd, rdd]    
    
    %%UPDATE HERE (OR AT END ???) : W FRAME 
    %(26) - (28)
    %%CHECK/ FIX THIS PLS 
   % W_frame(:,2) = y(1, 3) * W_frame(:,3)
   % W_frame(:,3) = -y(1,3) * (W_frame(:,2))

    r_wc_b      = r_w_b + d * W_frame(:,2)'; %(6)
    r_wc_b_dash = d *y(1,3)* W_frame(:,3) ; %(7)

    [c,cdash] = getcs(r_wc_b,y(1, 3));

    
    [w_1, w_2, w_3] = ws(y(:,1)); %pass w_b_n current to get W fram omega componenets
    
    ihub_b = ihub_bc + mhub * skew(r_bc_b) * skew(r_bc_b)'; %(18) P1
    irw_b = irw_wc + mrw * skew(r_wc_b)*skew(r_wc_b)';% (26)
    isc_b = ihub_b + irw_b; % total I in B frame is sum of I of hub + RWs

    %FIX US PLEASE


    %build A first(>1 RW)
    e = inv([J11 + mrw * d^2 - mrw^2 * d^2/msc]); % E = INV(A), (71)

    v = -mrw*d*W_frame(:,3)' * (r_c_ndd - 2 * skew(y(:,1)) * cdash' - ...
        skew(y(:,1)) * skew(y(:,1))* c' + (mrw/msc)*d*y(1,3)^2*W_frame(:,2))...
        +J23*(w_3^2 - w_2^2) + w_1*(J12 * w_3 - J13 * w_2) + ...
        w_2 * w_3 * ( J22 - J33 - mrw * d^2) - ...
        mrw * d * W_frame(:,3)' * skew(y(:,1))*skew(y(:,1)) * r_w_b' + us; %%FIXFIXFIX us %(73)

    tau_rhs = (- skew(y(:,1)) * isc_b - isc_b') * y(:,1) - msc * skew(c) * ...
        (r_c_ndd - 2 * skew(y(:,1)) * cdash' - skew(y(:,1))*skew(y(:,1))* c') ...
        + mrw * d * y(1,3)^2 * skew(r_wc_b)*W_frame(:,2) - irw_wc_dash * y(1,3) * W_frame(:,1) ...
        -skew(y(:,1)) * ( irw_wc * y(1,3) * W_frame(:,1) + mrw * skew(r_wc_b)*r_wc_b_dash)...
        - mrw*d*y(1,3)^2 * skew(c) * W_frame(:,2) ...
        - (irw_wc*W_frame(:,1) + mrw * d * (r_wc_b - skew(c)) * W_frame(:,3))* e' * v +Lb; %(81)

    f = - ((J11 + mrw * d^2) * W_frame(:,1)' + J12*W_frame(:,2)' + J13 * W_frame(:,3)'....
        + (mrw*d*W_frame(:,3)') * (skew(c) - skew(r_w_b))); %(72)

    ilhs = isc_b + msc * skew(c) * skew(c) + (irw_wc * W_frame(:,1) + mrw * d * (skew(r_wc_b)...
        - skew(c))*W_frame(:,3)) * e' * f; %(80)

    %%fix???
    sub_for_wbnd = tau_rhs \ ilhs ; %(79)
    %%also fix??? coupling ??????
    sub_for_Wd = e * f * sub_for_wbnd' + e * v  ; %(74)

    x1 = (tau_rhs \ ilhs)';  % 79  (1 x 3)
    x2 = r_c_ndd - 2 * skew(y(:,1)) * cdash' - skew(y(:,1))*skew(y(:,1))*c' + (mrw * d*W^2/msc) * W_frame(:,2)...
        +skew(c)*sub_for_wbnd' - (mrw*d/msc)*W_frame(:,3)*sub_for_Wd; %(10)
    x3 = [sub_for_Wd 0 0]'; %(74)
    dx= [x1 ,x2 , x3];
    % w_b_n d , r dd, W d

end

function [c,cdash] = getcs(r_wc_b, curr_W)
    load('parameters.mat')
    global W_frame
    r_wc_bdash  = d * curr_W  * W_frame(:,3)'; % (7)
    %r_wc_bddash = d * Wd * W_frame(:,3) - d * W^2 * W_frame(:,2);

    c      = (1/msc) * ( mhub * r_bc_b + mrw * r_wc_b); %(3)
    cdash  = (1/msc) * ( mrw * r_wc_bdash); %(4)
    %cddash = (1/msc) * ( mrw * r_wc_bddash);
end

function S = skew(v)
    S = [  0   -v(3)  v(2); v(3)   0   -v(1); -v(2)  v(1)   0  ];
end


function [w_1, w_2, w_3] = ws(w_b_n)
    load('parameters.mat');
    global W_frame
    w_1 = W_frame(:,1)' * w_b_n' ;
    w_2 = W_frame(:,2)' * w_b_n';
    w_3 = W_frame(:,3)' * w_b_n'; 
    %(45)-(47)
end