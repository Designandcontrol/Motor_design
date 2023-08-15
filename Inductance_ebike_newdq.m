function  [Output,L, flux_phase, one_iter]=Inductance_ebike_newdq(coil,geo,seq,one_slot_area)
% ampere_per_density=coil.ampere_per_density;
D_stator_outer =  geo.D_os;
D_rotor_outer= geo.D_or;
airgap = geo.airgap;
w_teeth = geo.W_t;
th_core = geo.T_s;
shoe_1=geo.shoe1;
shoe_2=geo.shoe2;
slot_number = geo.slot_number;
pole_number= geo.pole_number;
Num_turns=coil.Num_turns;
Num_parallel=coil.Num_parallel;
R_phase = geo.R_phase;
rpm = geo.rpm;
wn = 2*pi/60 * rpm;

angle_offset =geo.angle_offset * pi/180; % 원하는 각도를 라디안으로 변환   코드추가

Num_turns=coil.Num_turns
Num_parallel=coil.Num_parallel
current_density=coil.current_density
coil_diameter=coil.coil_diameter
Num_stranded_wire=coil.Num_stranded_wire
one_coil_area=(coil_diameter^2*pi/4)*10^-6*Num_stranded_wire;
one_slot_current=one_coil_area*current_density*Num_turns;

phase_current=geo.A_coil_width*geo.A_coil_height*coil.current_density*Num_parallel/coil.Num_turns

slot_current_density=current_density;
J_max=slot_current_density*2;

iq=[];id=[];
num_iter=21;
for n=1:(num_iter)
    for m=1:1:(num_iter)
        Jd((n-1)*(num_iter)+m)=-J_max+J_max/(num_iter-1)*(m-1)*2;
        Jq((n-1)*(num_iter)+m)=-J_max+J_max/(num_iter-1)*(n-1)*2;
    end
end
%%%%%%%%%%%
% Jd=0;
% Jq=0;
% temp_J=[];J =[];elec_deg=[];
% for n=1:length(Jd)
%     J(n)=(Jq^2+Jd^2)^0.5;
%     temp_J(n)=J(n);
%     elec_deg(n)=atan2(Jq(n),Jd(n));
% end
%%%%%%%%55
temp_J=[];J =[];elec_deg=[]; theta=0;
for n=1:length(Jd)
    J(n)=(Jq(n)^2+Jd(n)^2)^0.5;
    temp_J(n)=J(n);
    elec_deg(n)=atan2(Jq(n),Jd(n));
end

%%
rot_angle=0;time_passed=0;
flux_phase_a=[];flux_phase_b=[];flux_phase_c=[];
for r=1:length(J)
    tic

    J_a(r)=temp_J(r)*cos(elec_deg(r)+(rot_angle));
    J_b(r)=temp_J(r)*cos(elec_deg(r)-2*pi/3+(rot_angle));
    J_c(r)=temp_J(r)*cos(elec_deg(r)+2*pi/3+(rot_angle));

    mi_modifymaterial('b+',4,-J_b(r))
    mi_modifymaterial('b-',4,+J_b(r))   
    mi_modifymaterial('a+',4,-J_a(r))
    mi_modifymaterial('a-',4,+J_a(r))
    mi_modifymaterial('c+',4,-J_c(r))
    mi_modifymaterial('c-',4,+J_c(r))

  mi_analyze(1)
  mi_loadsolution()
    
    flux_a=0;    flux_b=0;    flux_c=0;
    for v=1:length(seq)*0.5
        tic
        mo_clearcontour()
        x1=(D_stator_outer*0.5-(shoe_1+shoe_2)*2)*cos((w_teeth/2)/(D_stator_outer*0.5-(shoe_1+shoe_2)*2)*1.4+pi/2+(2*pi/slot_number)*(v-1));
        y1=(D_stator_outer*0.5-(shoe_1+shoe_2)*2)*sin((w_teeth/2)/(D_stator_outer*0.5-(shoe_1+shoe_2)*2)*1.4+pi/2+(2*pi/slot_number)*(v-1));
        
        x2=(D_stator_outer*0.5-(shoe_1+shoe_2)*2)*cos(-(w_teeth/2)/(D_stator_outer*0.5-(shoe_1+shoe_2)*2)*1.4+pi/2+(2*pi/slot_number)*(v-1));
        y2=(D_stator_outer*0.5-(shoe_1+shoe_2)*2)*sin(-(w_teeth/2)/(D_stator_outer*0.5-(shoe_1+shoe_2)*2)*1.4+pi/2+(2*pi/slot_number)*(v-1));
       
        rotated_x1 = x1 * cos(angle_offset) - y1 * sin(angle_offset);
        rotated_y1 = x1 * sin(angle_offset) + y1 * cos(angle_offset);
        rotated_x2 = x2 * cos(angle_offset) - y2 * sin(angle_offset);                        %%%%전기각 맞추기 코드%%%%
        rotated_y2 = x2 * sin(angle_offset) + y2 * cos(angle_offset);
         
            mo_addcontour(rotated_x2, rotated_y2)  
        mo_addcontour(rotated_x1, rotated_y1) 
     
        
    
        
        if seq(v*2)==-10
            flux_temp=mo_lineintegral(0);
            flux_b=flux_b-flux_temp(1);
            mo_clearcontour()
        elseif seq(v*2)==10
            flux_temp=mo_lineintegral(0);
            flux_b=flux_b+flux_temp(1);
            mo_clearcontour()
        elseif seq(v*2)==-11
            flux_temp=mo_lineintegral(0);
            flux_a=flux_a-flux_temp(1);
            mo_clearcontour()
        elseif seq(v*2)==11
            flux_temp=mo_lineintegral(0);
            flux_a=flux_a+flux_temp(1);
            mo_clearcontour()
        elseif seq(v*2)==-12
            flux_temp=mo_lineintegral(0);
            flux_c=flux_c-flux_temp(1);
            mo_clearcontour()
        elseif seq(v*2)==12
            flux_temp=mo_lineintegral(0);
            flux_c=flux_c+flux_temp(1);
            mo_clearcontour()
        end
    end
    
    flux_phase.flux_phase_a(r)=flux_a*Num_turns/Num_parallel;
    flux_phase.flux_phase_b(r)=flux_b*Num_turns/Num_parallel;
    flux_phase.flux_phase_c(r)=flux_c*Num_turns/Num_parallel;
    mo_clearblock()
    mo_groupselectblock(1)
    mo_groupselectblock(2)
    
    T(r) = mo_blockintegral(22)
    Fx(r) = mo_blockintegral(18);
    Fy(r) = mo_blockintegral(19);
    mo_close()
    mi_clearselected()
    
    length(T)
    one_iter(r)=toc
    time_require=one_iter(r)*length(J)
    time_passed=(one_iter(r)+time_passed)/3600
    time_left=(one_iter(r)*length(J)-time_passed)/3600

%      J_a(n)=temp_J(n)*sin(elec_deg(n)+(rot_angle));
%      J_b(n)=temp_J(n)*sin(elec_deg(n)-2*pi/3+(rot_angle));
%      J_c(n)=temp_J(n)*sin(elec_deg(n)+2*pi/3+(rot_angle));
     
     i1=J_a(r)*one_slot_area/Num_turns*Num_parallel*10^6;
     i2=J_b(r)*one_slot_area/Num_turns*Num_parallel*10^6;
     i3=J_c(r)*one_slot_area/Num_turns*Num_parallel*10^6;


    
    idq = abc2dq(i1,i2,i3,theta*pi/180);
    i_dd(r)=idq(1);
    i_qq(r)=idq(2);



    flux11=flux_phase.flux_phase_a(r);
    flux21=flux_phase.flux_phase_b(r);
    flux31=flux_phase.flux_phase_c(r);
    fluxdq1 = abc2dq(flux11,flux21,flux31,theta*pi/180);
    
    fluxd(r)=fluxdq1(1)
    fluxq(r)=fluxdq1(2)
end

num_iter=21
for p=1:1:num_iter
    for q=1:1:num_iter
        i_d(p,q)=i_dd((p-1)*num_iter+q);          %% 인덕턴스는 비선형이기 때문에 고려
        i_q(p,q)=i_qq((p-1)*num_iter+q);         
       flux_q(p,q)=fluxq((p-1)*num_iter+q);
       flux_d(p,q)=fluxd((p-1)*num_iter+q);

    
    end
end



for p=1:1:size(flux_d,1)
    L_dd(p,:)= gradient(flux_d(p,:))./gradient(i_d(p,:))*1000;
    L_qq(:,p)=gradient(flux_q(:,p))./gradient(i_q(:,p))*1000;
    L_dq(p,:)=gradient(flux_d(:,p))./gradient(i_q(:,p))*1000;
    L_qd(:,p)=gradient(flux_q(p,:))./gradient(i_d(p,:))*1000;
end

Output.T=T;
Output.Fx=Fx;
Output.Fy=Fy;

 figure(134)
 title('TORQUE');
 plot(T)
 xlabel('elec. angle[deg]')
 ylabel('Torque[Nm]')
 grid on


L.i_d=i_d;
L.i_q=i_q;
L.flux_q=flux_q;
L.flux_d=flux_d;
L.L_dd=L_dd;                         
L.L_qq=L_qq;
L.L_dq=L_dq;
L.L_qd=L_qd;
L.flux11=flux11;
L.flux21=flux21;
L.flux31=flux31;
%%%%% Vd Vq %%%%%%%%%

% Vd = R_phase * L.i_d - wn*L.L_qq * L.i_q
% Vq = R_phase * L.i_q + wn*flux_d
% Va = (Vd^2+Vq^2)^0.5;
% L.Vd = Vd;
% L.Vq = Vq;
% L.Va = Va;
% Ia = (i_d^2+i_q^2)^0.5;
% phasor;
%%%% plot %%%%%%%%%%






end

