load('parameters.mat')
%%%%%%%%%%%%%%%%%%%%%%
dt = 0.0001;
tspan = [0:dt:3];
y0= [0 0 0];
Y(:,1) = y0;
yin = y0;
%for i=1:tspan(end)/dt
%    time = i*dt;
 %   yout = rk4singlestep(@(t,y)func(t,y),dt,time,yin);
%    Y =[Y yout];
%    yin = yout;
%end
%plot2(Y(1,:),Y(2,:))
%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%FUNCTIONS%%%%%%%


function xout = rk4singlestep(f, dt, tk, xk)
    f1 = f(tk,xk);
    f2 = f(tk + dt/2, xk + (dt/2)*f1);
    f3 = f(tk + dt/2, xk + (dt/2)*f2);
    f4 = f(tk+dt, xk + dt*f3);

    xout = xk + (dt/6) * (f1 + 2*f2 + 2*f3 + f4);
end

function dx = func(t,y)
    load('parameters.mat')

%statespace =  [ wd, rdd]
    [c,cdash, ~, ~] = getcs();
 
    e = inv([J11 + mrw * d^2 - mrw^2 * d^2/msc]);

    v = [-mrw*d*W_frame(:,3)' * (r_c_ndd - 2 * skew(w_b_n) * cdash -...
        skew(w_b_n) * skew(w_b_n)* c + (mrw/msc)*d*W^2*W_frame(:,2))...
        +J23*(w_3^2 - w_2^2) + w_1*(J12 * w_3 - J13 * w_2) + ...
        w_2 * w_3 * ( J22 - J33 - mrw * d^2) - ...
        mrw * d * W_frame(:,3)' * skew(w_b_n)*skew(w_b_n) * r_w_b + us];

    tau_rhs = (- skew(w_b_n) * isc_b - isc_b') * w_b_n -msc * skew(c) * ...
        (r_c_ndd - 2 * skew(w_b_n) * cdash - skew(w_b_n)*skew(w_b_n)* c) ...
        + mrw * d * W^2 * skew(r_wc_b)*W_frame(:,2) - irw_wc_dash * W * W_frame(:,1) ...
        -skew(w_b_n) * ( irw_wc * W * W_frame(:,1) + mrw * skew(r_wc_b)*r_wc_b_dash)...
        - mrw*d*W^2 * skew(c) * W_frame(:,2) ...
        - (irw_wc*W_frame(:,1) + mrw*d * (r_wc_b - skew(c)) * W_frame(:,3))* e * v +Lb;

    f = - [(J11 + mrw * d^2) * W_frame(:,1)' + J12*W_frame(:,2)' + J13 * W_frame(:,3)'....
        +(mrw*d*W_frame(:,3)') * (skew(c) - skew(r_w_b))];

    ilhs = isc_b + msc * skew(c) * skew(c) + (irw_wc * W_frame(:,1) + mrw * d * (skew(r_wc_b) - skew(c))*W_frame(:,3)) * e' * f;

    sub_for_wbnd = tau_rhs \ ilhs;

    sub_for_Wd = e' * f * sub_for_wbnd + e' * v ; 

    dx = [
    isc_b \ (rhs30  - lhs30); %w_b_nd%30
    F / msc; %rdd(2)
    
    tau_rhs \ ilhs; % P2 78 wd
    r_c_ndd - 2 * skew(w_b_n) * cdash - skew(w_b_n)*skew(w_b_n)*c + (mrw * d*W^2/msc) * W_frame(:,2)...
    +skew(c)*sub_for_wbnd - (mrw*d/msc)*W_frame(:,3)*sub_for_Wd;%(P2 10) rdd
        ];
end

function [c,cdash, cddash, cdd] = getcs()
    load('parameters.mat')
    r_wc_b = r_w_b' +  d * W_frame(:, 2) ; %(24)
    r_wc_b_dash = d *W* W_frame(:,3) ; % (6)
    r_wc_b_ddash = -d * W^2 * W_frame(:,2) + d * Wd * W_frame(:,3); %(7)

    c   = (1 / msc) * (mhub * r_bc_b + mrw * r_wc_b'); % (3)%wheel center of mass position vector relative to B 
    cdash = (1 / msc) * mrw * r_wc_b_dash; %(4)
    cddash = (1 / msc) * mrw * r_wc_b_ddash; %(4)
    cdd = cddash' + 2 * cross(w_b_n, cdash') + cross(w_b_nd, c) + cross(w_b_n, cross(w_b_n,c)); %(8) 

end

function [W, Wd] = wcomponents(W_frame)
%FIX THIS IT IS WRONG PLS W_B_N CHANGE 
    load('parameters.mat');
    % w_update  Returns the updates w_frame and its velocity components .
    % Input: w - 3x3 w_frame matrix 
    % Output: w_updated - 3x3 new matrix
    %W frame components (36) (35)
    w_1 = W_frame(:,1)' .* w_b_n ;
    w_2 = W_frame(:,2)' .* w_b_n;
    w_3 = W_frame(:,3)' .* w_b_n;
    
    w_1d = W_frame(:,1)' .* w_b_nd ;
    w_2d = W_frame(:,2)' .* w_b_nd + W * w_3;
    w_3d = W_frame(:,3)' .* w_b_nd - W * w_2;

    W = [w_1; w_2; w_3];
    Wd = [w_1d; w_2d; w_3d];
end

function irw_wc = inertupdate()
    load('parameters.mat');
    irw_wc = J11 * W_frame(:,1)* W_frame(:,1)' ...
        + J12 * W_frame(:,1)* W_frame(:,2)' ...
        + J13 * W_frame(:,1)* W_frame(:,3)' ...
        + J12 * W_frame(:,2)* W_frame(:,1)' ...
        + J22 * W_frame(:,2)* W_frame(:,2)' ...
        + J23 * W_frame(:,2)* W_frame(:,3)' ...
        + J13 * W_frame(:,3)* W_frame(:,1)' ...
        + J23 * W_frame(:,3)* W_frame(:,2)' ...
        + J33 * W_frame(:,3)* W_frame(:,3)'; %(21)
end


function [H_sc_bd, Hrw_bd,Hhub_bd, rhs30, isc_b] = hscbd(ihub_bc, r_bc_b, r_wc_b, r_wc_b_dash)
    load('parameters.mat');
    
    irw_wc_dash = [0, -J13 ,J12 ; -J13, -2*J23 , J22 - J33; J12, J22 - J33, 2*J23] * W; %(23)
    irw_b_dash = irw_wc_dash + mrw * skew(r_wc_b_dash) * skew(r_wc_b)' + mrw * skew(r_wc_b) * skew(r_wc_b_dash)' ; % (28)

    ihub_b_dash = 0 * eye(3);
    isc_b_dash = ihub_b_dash + irw_b_dash ;
    [ihub_b, irw_b, isc_b ] = getis(ihub_bc, r_bc_b, r_wc_b);

    H_sc_bd = isc_b * w_b_nd' + skew(w_b_n)* isc_b *w_b_n' + isc_b_dash * w_b_n'...
        + irw_wc_dash * W * W_frame(:,1) + irw_wc * Wd * W_frame(:,1) ...
        + skew(w_b_n) * (irw_wc * W * W_frame(:, 1) + mrw * skew(r_wc_b) * r_wc_b_dash)...
        + mrw * skew(r_wc_b) * ( d*Wd * W_frame(:,3) - d * W^2 * W_frame(:, 2)); %(29)

    Hrw_bd = irw_b_dash * w_b_n' + irw_b * w_b_nd' + skew(w_b_n) * irw_b * w_b_n' ...
        + irw_wc_dash * W * W_frame(:, 1) + irw_wc * Wd * W_frame(:,1) + skew(w_b_n) * irw_wc * W * W_frame(:,1)...
        + mrw *  cross (r_wc_b, (d * Wd * W_frame(:,3) -d*W^2 * W_frame(:, 2) )) + mrw * cross(w_b_n,cross(r_wc_b, r_wc_b_dash))'; %(27)

    Hhub_bd = ihub_b * w_b_nd' + skew(w_b_n) * ihub_b * w_b_n'; %(19)

    rhs30 = mrw  * skew(r_wc_b) * d * W^2 * W_frame(:,2) - skew(w_b_n) * (irw_wc * W * W_frame(:,1) + mrw * skew(r_wc_b) * r_wc_b_dash) - irw_wc_dash * W * W_frame(:,3) - skew(w_b_n) * isc_b * w_b_n' - isc_b_dash * w_b_n' + Lb;


end


function [ihub_b, irw_b, isc_b ]= getis(ihub_bc, r_bc_b, r_wc_b)
    load('parameters.mat');
    ihub_b = ihub_bc + mhub * skew(r_bc_b) * skew(r_bc_b)'; %(18)
    irw_b = irw_wc + mrw * skew(r_wc_b)*skew(r_wc_b)';% (26)
    isc_b = ihub_b + irw_b;

end


function [Hrw_b, Hhub_b]  = hbs(r_bc_b, r_bc_bd, r_wc_b, r_wc_bd)
    load('parameters.mat');
    Hrw_b = (irw_wc * (w_b_n' + W * W_frame(:,1)))' + mrw * cross(r_wc_b, r_wc_bd); %(12)
    Hhub_b = (ihub_bc * w_b_n')' + mhub * cross(r_bc_b, r_bc_bd); %(11)
end

function S = skew(v)

    S = [  0   -v(3)  v(2); v(3)   0   -v(1); -v(2)  v(1)   0  ];
end

