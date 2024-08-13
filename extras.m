%r_wc_bdd = -d*W^2 * W_frame(:,2)+ d*Wd*W_frame(:,3) + cross(w_b_nd, r_wc_b) + 2 * cross(w_b_n, d*W*W_frame(:,3)) + cross( w_b_n, cross(w_b_n, r_wc_b)) ; %(25)
%r_bc_bdd = cross(w_b_nd, r_bc_b) + cross(w_b_n, cross(w_b_n, r_bc_b)); %(16)
%Hhub_bd = ihub_bc * w_b_nd' + cross( w_b_n, ihub_bc*w_b_n')' + mhub * cross(r_bc_b, r_bc_bdd)'; %(15)
%[Hrw_b, Hhub_b] = hbs(r_bc_b, r_bc_bd, r_wc_b, r_wc_bd);




%Equations of  Motion%
r_wc_b = r_w_b' +  d * W_frame(:, 2) ; %(24)
r_wc_b_dash = d *W* W_frame(:,3) ; % (6)
[c,cdash, cddash, cdd] = getcs();
r_b_ndd = F/msc - cdd ;  %(1) (2) --------------------------------------------
irw_wc = inertupdate();
%r_wc_bd = r_wc_b_dash' + cross(w_b_n,r_wc_b) ; 
%r_bc_bd = cross(w_b_n, r_bc_b);


[H_sc_bd, Hrw_bd, Hhub_bd, rhs30, isc_b] = hscbd(ihub_bc, r_bc_b, r_wc_b, r_wc_b_dash);


%(13) H_sc_bd = RHS
%rhs_13 = Lb + msc * cross(r_b_ndd,c);
%lhs_13 = Hhub_bd + Hrw_bd; 
c   = (1 / msc) * (mhub * r_bc_b + mrw * r_wc_b');
Lb  = Hhub_bd + Hrw_bd - msc * cross(r_b_ndd,c); 
us_check = J11 * (W_frame(:,1)' * w_b_nd' + Wd); %rw motor torque

%%%%%%%%%%%%%%%%%%%%%%
[w_w, w_wd] = wcomponents(W_frame);



%30
lhs30 = msc*skew(c) * r_b_ndd' + (irw_wc * W_frame(:,1) + mrw*d*skew(r_wc_b) * W_frame(:,3) ) * Wd;

%%%%%ss equations%%%%%
w_b_nd =  isc_b \ (rhs30  - lhs30); %30
r_b_ndd = -((J11 + mrw * d^2) * W_frame(:,1)' + J12 *  W_frame(:,2)' + J13 * W_frame(:,3)'...
    -mrw *d * W_frame(:,3)' * skew(r_w_b)) * w_b_nd - (J11+mrw*d^2)* Wd...
    + J23 * ( w_w(3)^2 - w_w(2)^2) + w_w(1)*(J12*w_w(3) - J13 * w_w(2)) ...
    + w_w(2)*w_w(3) * (J22 -J33 -mrw * d^2) - mrw * d * W_frame(:,3)' * skew(w_b_n) * skew(w_b_n) * r_w_b' + us; %(41)