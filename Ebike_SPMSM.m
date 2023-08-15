clear all;close all;clc;
%Measured
%fbm75n68 75A 68V
%Torque=1000W/(480rpm/60*2*pi)=19.98Nm
R_uv=0.275;R_vw=0.275;R_wu=0.275;
L_uv=0.416e-3;L_vw=0.479e-3;L_wu=0.485e-3;
R_s=(R_uv+R_vw+R_wu)/3*0.5 %0.1436 ohm
L_s=(L_uv+L_vw+L_wu)/3*0.5 %23 mH

K_e=48/3^0.5/(480/60*2*pi) %Epahse Vpk/(rad/sec)
K_t=K_e*3
Phase_current_est=19.98/K_t
geo.angle_offset = 6+0.6667;   %% 전기각 코드 변수 추가

%Resistancer4
Num_slot=27
Num_pole=30
Num_series=9          
Num_parallel=1
Num_turns=15
Num_stranded_wire=4
% with insulation
% AWG 23 0.574
coil_diameter_w_insu=0.500
%without insulation
coil_diameter=coil_diameter_w_insu/1.05


%coil
coil_depth=25
coil_width=5      
coil_f_radius=0.2
winding_temperature=25
resistivity_cu=(3.76*winding_temperature+873)*10^(-6)/55


%*1.05 is a factor to consider increasing winding length for crossing
%winding
for n=1:1:Num_turns
%     length_turn(n)=2*(coil_depth+coil_width)*(coil_f_radius+coil_diameter_w_insu/2+3^0.5/2*coil_diameter_w_insu*(n-1))
    length_turn(n)=2*coil_depth+2.4*coil_width+0.1;
end
end_winding=0;
length_tot=sum(length_turn)+end_winding;

coil_resistance=resistivity_cu*(length_tot)/(pi*coil_diameter^2/4*Num_stranded_wire)
%theoretical result
R_phase=coil_resistance*Num_series/Num_parallel
%experimental result
% R_phase_exp=0.115

current_density=0;
I_rated_est=coil_diameter^2*pi/4*Num_parallel*current_density*Num_stranded_wire
coil.current_density=current_density;
coil.Num_turns=Num_turns;
coil.Num_parallel=Num_parallel;
coil.Num_series=Num_series;
coil.coil_diameter=coil_diameter;
coil.Num_stranded_wire=Num_stranded_wire;
coil.length_tot=length_tot;
% stack length
geo.L_st=20;
% teeth height
%geo.T_h=13.4;   %%%%/입력 필요x
% teeth width
geo.W_t=4.2;
geo.pole_number=Num_pole;
geo.slot_number=Num_slot;
geo.airgap=1;
geo.T_m=2.5;
geo.pole_arc_ratio=0.99;
geo.slot_ratio = 0.82;
geo.D_os=105;
geo.D_is=60;
geo.D_or=111.2; 
geo.T_s=7.2;
geo.T_r=6;
geo.shoe1=1.0;
geo.shoe2=1.0;
geo.A_coil_width=1.5;
geo.A_coil_height=11.3;
rpm=480;
geo.R_phase=R_phase;
geo.rpm=rpm;

% % magnet with arc
% [seq,one_slot_area]=Draw_ebike(coil,geo);
% magnet with arc
% [seq,one_slot_area]=Draw_ebike(coil,geo);
% [Output,E,flux_e]=BackEMF_ebike(coil,geo,seq,one_slot_area,rpm)
% save('backEMF_FT0500')

%[seq,one_slot_area]=Draw_ebike(coil,geo)
%[T_i,L,flux_i,solve_time_i]=Inductance_ebike(coil,geo,seq,one_slot_area)
%save('Inductance_FT0500')

%%
current_density=1.5*2^0.5;
coil.current_density=current_density;
I_rated_est=geo.A_coil_width*geo.A_coil_height*coil.current_density*Num_parallel/coil.Num_turns
%%%%%               bemf        %%%%%%%%%%%%%
   [seq,one_slot_area]=Draw_ebike_newdq(coil,geo);
  [Output,E,flux_e]=BackEMF_ebike_newdq(coil,geo,seq,one_slot_area,rpm)
  save('backEMF_FT0500_newdq')


% % %%%%%     inductance          %%%%%%%%%%%%%55
%    [seq,one_slot_area]=Draw_ebike_newdq(coil,geo)
%   [T_i,L,flux_i,solve_time_i]=Inductance_ebike_newdq(coil,geo,seq,one_slot_area)
%    save('Inductance_FT0500_newdq')
% % % 
% %%
% phase_shift=pi+pi/12;
% [widning_factor]=winding_poleslot(pole_number,slot_number,phase_shift)