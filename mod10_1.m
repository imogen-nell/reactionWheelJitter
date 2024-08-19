%mod10 - same as mod 10 but with pointing error from mod 14
load("parameters2.mat")

%set sizes
dt = 0.001; 
tspan = 0:dt:5;

%%initial conditions
w_b_n_0 = [0 0 0]'; 
rd_0 = [0 0 0]'; %r_b_n
theta_0 = 0;
pyr_0 = [0 0 0]';
w2_0 = W_frame_init(:,2);
w3_0 = W_frame_init(:,3);
r_0 = rd_0 + (1/msc) * ( mhub * r_bc_b + mrw *(r_w_b + d * w2_0)) ;
global us
us = NaN ;

y0 = [ w_b_n_0 ; rd_0; W_init ; theta_0 ; pyr_0; r_0];

[t,y] = ode45(@func,tspan,y0);
%graph(t, y);
graphEUL(t, y)
graphMRP(t, y)
graphWBN(t, y)
graphPYR(t, y)
graphWHEELANGLE(t,y)
graphWHEELSPEED(t,y)
graph3D(y)
graphRCOM(t,y)



function dx = func(t,y)
    global us
    load('parameters2.mat')
    if isnan(us)
        us = us_0;
        fprintf('US reset to US_0 %f\n', us);
    end

    %y = [y(1:3) y(4:6) [y(7:8) ;0]];

    w_b_n    = y(1:3);
    r_b_nd   = y(4:6);
    W        = y(7);
    theta    = y(8);%rw displacement 
    pyr      = y(9:11);
    r        = y(12:14);
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
    
    [c,cdash] = getcs(r_wc_b, r_wc_b_dash);

    e = inv([J11 + mrw * d^2 - mrw^2 * d^2/msc]);       % E = INV(A), (71)


    F = - ((J11 + mrw * d^2) * g' + J12*wframe_2' + ...
        J13 * wframe_3'+ (mrw*d* wframe_3') * (skew(c) - ....
        skew(r_w_b)));                                             %(72)

    i_lhs = isc_b + msc * skew(c)*skew(c) +( irw_wc*g ...
        +mrw * d * (skew(r_wc_b) - skew(c))*wframe_3)*e*F;      %(80)

    %wframe = [g ,wfame_2, wframe_3]
    [w_1, w_2, w_3] = ws(w_b_n, [g ,wframe_2, wframe_3]); %
    

    v = -mrw*d*wframe_3' * (r_c_ndd - 2 * skew(w_b_n) * cdash - ...
        skew(w_b_n) * skew(w_b_n)* c + (mrw/msc)*d*W^2*wframe_2)...
        +J23*(w_3^2 - w_2^2) + w_1*(J12 * w_3 - J13 * w_2) ...
        + w_2 * w_3 * ( J22 - J33 - mrw * d^2) - mrw * d * ...
        wframe_3' * skew(w_b_n)*skew(w_b_n) * r_w_b + us;       %(73)

    Lb = W^2 * (U_d * wframe_3 + U_s * skew(r_w_b)*wframe_3);%not sure about tis thb 

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


    
    pyrd = ( 1/cos(pyr(2)) )  * [0 sin(pyr(3)) cos(pyr(3)) ; ...
        0 cos(pyr(3)) * cos(pyr(2)) sin(pyr(3)) * cos(pyr(2));...
        cos(pyr(2)) sin(pyr(3))*sin(pyr(2)) cos(pyr(3)) * sin(pyr(2))] * w_b_n;

%r  = r_b_n + (1/msc) * ( mhub * r_bc_b + mrw *(r_w_b + d * w2_0)) ;

    rd = r_b_nd + (mrw/msc) * d* (  W * wframe_3 + cross(w_b_n,wframe_2));


    thetad = W ;

    x1 = w_b_nd;
    x2 = r_b_ndd;
    x3 = Wd;
    x4 = thetad;
    x5 = pyrd ; 
    x6 = rd;
    dx = [x1; x2; x3; x4; x5 ; x6];

end

function S = skew(v)
    S = [  0   -v(3)  v(2); v(3)   0   -v(1); -v(2)  v(1)   0  ];
end

function [c,cdash] = getcs(r_wc_b, r_wc_b_dash)
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
    
    v_rot = dot(v, g) * g + cos(theta) * (v - dot(v, g) * g ) + sin(theta) * cross(g, v);

end
function graphMRP(t,y)

    [dim ~] = size(y);

    Ws = y(:,1:3);
    psiTrack = [0 0 0 0 ]';
    for i = 2:dim
        %axis of rotation
        W = Ws(i,:);
        %angle of rotation
        theta =( t(i)-t(i-1) )* norm(W);
        %theta =norm(W);
        psi = tan(theta/4)*W /theta;
        
        psiTrack = [psiTrack [psi(1) ,psi(2) ,psi(3) ,norm(psi)]'];
    end

    figure;
    title('MRP ', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'g');
    hold on;
    plot(t, psiTrack(1, :), '-b', 'DisplayName', 'psi 1');
    plot(t, psiTrack(2, :), '-b', 'DisplayName', 'psi 2');
    plot(t, psiTrack(3, : ), '-b', 'DisplayName', 'psi 3');
    plot(t, psiTrack(4, : ), '-g', 'DisplayName', 'y');
    hold off;

end

%y0 = [ w_b_n_0 ; rd_0; [W_init theta_0 ]'];


function graphEUL(t,y)
%%THIS IS DEFINITELY NOT RIGHT PLS 
    [dim ~] = size(y);

    Ws = y(:,1:3); % W_b_n
    zxy = [0 0 0 0]';
    for i = 2:dim
        %axis of rotation
        W = Ws(i,:) / norm(Ws(i,:));
        %angle of rotation
        theta = ( t(i)-t(i-1) ) * norm(W) ;%%%UNITS
        %theta = norm(W);
        K = skew(W);
        %rotation matrix of theta about axis 
        rot = eye(3) + sin(theta) * K + (1-cos(theta)) * K^2;
        %euler angles defined by rot matrix in order z x y 
        eul = rad2deg(rotm2eul(rot)) /60^2 * 1000;%convert to milli arcsec
        zxy = [zxy [eul'; rad2deg(theta)/60^2*1000 ]];
    end

    figure;
    title('Euler Angles / MilliArcsec ', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'g');
    hold on;
    plot(t, zxy(1, :), '-b', 'DisplayName', 'z');
    plot(t, zxy(2, :), '-b', 'DisplayName', 'x');
    plot(t, zxy(3, : ), '-b', 'DisplayName', 'y');
    plot(t, zxy(4, : ), '-g', 'DisplayName', 'y');
    hold off;

end

function graphWBN(t,y)
    figure;
    title('Space Craft Angular Velocity', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'g');
    hold on;
    plot(t, y(:, 1), '-g', 'DisplayName', 'W_b_n (1)');
    plot(t, y(:, 2), '-r', 'DisplayName', 'W_b_n (2)');
    plot(t, y(:, 3), '-b', 'DisplayName', 'W_b_n (3)');
    hold off;
end


function graphPYR(t, y)
    figure;
    title('Pitch Yaw Roll / Degrees', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'g');
    hold on;
    plot(t, rad2deg(y(:, 9)), '-g', 'DisplayName', 'W_b_n (1)');
    plot(t, rad2deg(y(:, 10)), '-g', 'DisplayName', 'W_b_n (2)');
    plot(t, rad2deg(y(:, 11)), '-g', 'DisplayName', 'W_b_n (3)');
    hold off;
end 


function graphWHEELANGLE(t,y)
    figure;
    plot(t, rad2deg(mod(y(:,8),2*pi)), '-m', 'DisplayName', 'RW 1 ');
    title('RW 1: Wheel angle / DEG ', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'g');

end

function graphWHEELSPEED(t,y)
    figure; % 7 
    
    plot(t, y(:,7)*60/(2*pi), '-m', 'DisplayName', 'RW 1 '); % convert TO RPM
    title('RW 1: Wheel Speed / RPM ', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'g');

end

function graph3D(y)
    load('parameters2.mat')

    x0 = y(1,12);
    y0 = y(1,13);
    z0 =  y(1,14);
    [dim, ~] = size(y);
    cols = [ 'r' 'g' 'b' 'c' 'm' 'y'];
    figure
    %plot sc omega 
    scale2 = 20000;
    for i =1:100:dim
        quiver3(x0, y0, z0, y(i,1)*scale2 , y(i,2)*scale2,  y(i,3)*scale2, 'm', 'LineWidth', 2);
        hold on;
    end
    %plot rw omega 

    rwx0 = r_w_b(1);
    rwy0 = r_w_b(2);
    rwz0 = r_w_b(3);
    scale = 1;
    
    %plot rw direction
    quiver3(rwx0, rwy0, rwz0, g(1) * scale,g(2) * scale,  g(3) * scale,  'c', 'LineWidth', 2);
  %  plot3(rwx0,rwy0,rwz0,'bo', 'MarkerSize', 4, 'MarkerFaceColor', 'b')
    %plot3(rwx0+g(1),rwy0+g(2), rwz0+g(3),'bo', 'MarkerSize', 4, 'MarkerFaceColor', 'k')




end


function graphRCOM(t, y)
    figure;
    title('COM Translation', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'g');
    hold on;
    plot(t, y(:, 12), '-b', 'DisplayName', 'component (1)');
    plot(t, y(:, 13), '-m', 'DisplayName', 'component (2)');
    plot(t, y(:, 14), '-c', 'DisplayName', 'component (3)');
    hold off;
end 