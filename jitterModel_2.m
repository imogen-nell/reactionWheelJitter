load('parameters.mat')

r_wc_b = r_w_b' +  d * W_frame(:, 2) ; %(24)
ihub_b = ihub_bc + mhub * skew(r_bc_b) * skew(r_bc_b)'; %(18)
irw_b = irw_wc + mrw * skew(r_wc_b)*skew(r_wc_b)';% (26)
isc_b = ihub_b + irw_b;
c   = (1 / msc) * (mhub * r_bc_b + mrw * r_wc_b'); % (3)%wheel center of mass position vector relative to B 
lhs30 = msc*skew(c) * r_b_ndd' + (irw_wc * W_frame(:,1) + mrw*d*skew(r_wc_b) * W_frame(:,3) ) * Wd;


w_b_nd =  isc_b \ (rhs30  - lhs30); %30
%r_b_ndd = -((J11 + mrw * d^2) * W_frame(:,1)' + J12 *  W_frame(:,2)' + J13 * W_frame(:,3)'...
 %   -mrw *d * W_frame(:,3)' * skew(r_w_b)) * w_b_nd - (J11+mrw*d^2)* Wd...
 %   + J23 * ( w_w(3)^2 - w_w(2)^2) + w_w(1)*(J12*w_w(3) - J13 * w_w(2)) ...
 %   + w_w(2)*w_w(3) * (J22 -J33 -mrw * d^2) - mrw * d * W_frame(:,3)' * skew(w_b_n) * skew(w_b_n) * r_w_b' + us; %(41)




function S = skew(v)
    S = [  0   -v(3)  v(2); v(3)   0   -v(1); -v(2)  v(1)   0  ];
end