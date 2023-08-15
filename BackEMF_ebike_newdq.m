function [Output,E,flux_phase]=BackEMF_ebike_newdq(coil,geo,seq,one_slot_area,rpm)
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

angle_offset =geo.angle_offset * pi/180; % 원하는 각도를 라디안으로 변환   코드추가

current_density=coil.current_density
Num_turns=coil.Num_turns;
Num_parallel=coil.Num_parallel;
coil_diameter=coil.coil_diameter
Num_stranded_wire=coil.Num_stranded_wire  

one_coil_area=(coil_diameter^2*pi/4)*10^-6*Num_stranded_wire;    %%코일면적
one_slot_current=one_coil_area*current_density*Num_turns;       
slot_current_density=one_slot_current/one_slot_area
coil_fill_factor=slot_current_density/current_density

phase_current=one_slot_area*(current_density*10^6)*Num_parallel/coil.Num_turns  %%상전류

Jd=0;                         %% 자기장 최대하기 위해 d축전류 0?
Jq=current_density;
temp_J=[];I =[];elec_deg=[];

for n=1:length(Jq)
    J(n)=(Jq(n)^2+Jd(n)^2)^0.5;
    temp_J(n)=J(n);
    elec_deg(n)=atan2(Jq(n),Jd(n));
end

rot_angle=0;
time_passed=0;
flux_phase_a=[];
flux_phase_b=[];
flux_phase_c=[];

num_iter=36;          %% 회전각도 설정

%%%%%           자속 경로설정 대칭맞추기위한 코드 추가       %%%%%%%%%%%%%
   mi_clearselected()
    mi_selectcircle(0,0,(D_stator_outer*0.5+airgap*0.5),4)
    mi_selectgroup(2)
   mi_moverotate(0,0,-4*180/pole_number/num_iter)

%%%%%%%%%%%                                      %%%%%%%%%%%%%%%%%%55

for r=1:num_iter
%     J_a(r)=temp_J(1)*cos(elec_deg(1)+(rot_angle));
%     J_b(r)=temp_J(1)*cos(elec_deg(1)-2*pi/3+(rot_angle));
%     J_c(r)=temp_J(1)*cos(elec_deg(1)+2*pi/3+(rot_angle));
        J_a(r)=temp_J(1)*cos(elec_deg(1));
        J_b(r)=temp_J(1)*cos(elec_deg(1)-2*pi/3);
        J_c(r)=temp_J(1)*cos(elec_deg(1)+2*pi/3);
%     
    mi_modifymaterial('b+',4,-J_b(r))
    mi_modifymaterial('b-',4,+J_b(r))
    mi_modifymaterial('a+',4,-J_a(r))
    mi_modifymaterial('a-',4,+J_a(r))
    mi_modifymaterial('c+',4,-J_c(r))
    mi_modifymaterial('c-',4,+J_c(r))
    
    mi_analyze(1)
    mi_loadsolution()
    
    %     mo_showdensityplot(0,0,1.6,0,'bmag')
    %     path='C:/femm42/'
    %     name_fem=tostring(r)..".png"
    %     mo_savebitmap(path..name_fem)
    flux_a=0;flux_b=0;flux_c=0;


           %%%%%%%%%%% 전기각 맞추기위한 angle offset 추가

      
     for v=1:length(seq)*0.5
        tic
           mo_clearcontour()
             

        x1=(D_stator_outer*0.5-(shoe_1+shoe_2)*2)*cos(-r*4*pi/pole_number/num_iter+(w_teeth/2)/(D_stator_outer*0.5-(shoe_1+shoe_2)*2)*1.4+pi/2+(2*pi/slot_number)*(v-1));   
        y1=(D_stator_outer*0.5-(shoe_1+shoe_2)*2)*sin(-r*4*pi/pole_number/num_iter+(w_teeth/2)/(D_stator_outer*0.5-(shoe_1+shoe_2)*2)*1.4+pi/2+(2*pi/slot_number)*(v-1));
     
        x2=(D_stator_outer*0.5-(shoe_1+shoe_2)*2)*cos(-r*4*pi/pole_number/num_iter-(w_teeth/2)/(D_stator_outer*0.5-(shoe_1+shoe_2)*2)*1.4+pi/2+(2*pi/slot_number)*(v-1));
        y2=(D_stator_outer*0.5-(shoe_1+shoe_2)*2)*sin(-r*4*pi/pole_number/num_iter-(w_teeth/2)/(D_stator_outer*0.5-(shoe_1+shoe_2)*2)*1.4+pi/2+(2*pi/slot_number)*(v-1));
       
         rotated_x1 = x1 * cos(angle_offset) - y1 * sin(angle_offset);
        rotated_y1 = x1 * sin(angle_offset) + y1 * cos(angle_offset);
         rotated_x2 = x2 * cos(angle_offset) - y2 * sin(angle_offset);
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
    flux_phase_a(r)=flux_a*Num_turns/Num_parallel;
    flux_phase_b(r)=flux_b*Num_turns/Num_parallel;
    flux_phase_c(r)=flux_c*Num_turns/Num_parallel;
    mo_clearblock()
    mo_groupselectblock(1)
    mo_groupselectblock(2)
    
    T(r) = mo_blockintegral(22)
    Fx(r) = mo_blockintegral(18);
    Fy(r) = mo_blockintegral(19);
    
    mo_close()
    mi_clearselected()
    mi_selectcircle(0,0,(D_stator_outer*0.5+airgap*0.5),4)
    mi_selectgroup(2)
    
    mi_moverotate(0,0,-4*180/pole_number/num_iter)           %% if polenumber 30 > 2/3 씩 iter(36)[24도] 만큼 회전
    rot_angle=-2*pi/num_iter+rot_angle;
end
%%

Output.T=T;
Output.Fx=Fx;
Output.Fy=Fy;
Output.Ip=phase_current;
% Output.coil_fill_factor=coil_fill_factor;

flux_phase.flux_phase_a=flux_phase_a;
flux_phase.flux_phase_b=flux_phase_b;
flux_phase.flux_phase_c=flux_phase_c;
 figure(113)
 title('TORQUE')
 plot(T)
 xlabel('elec. angle[deg]')
 ylabel('Torque[Nm]')
 grid on
 ylim([-4 4])


%%
Npoints=num_iter;
Hz=pole_number/2*rpm/60;
elec_angle=0:360/num_iter:(360-360/num_iter);

for n=1:length(flux_phase_a)-1
    E.a(n)=-(flux_phase_a(n+1)-flux_phase_a(n))*Hz*Npoints;
    E.b(n)=-(flux_phase_b(n+1)-flux_phase_b(n))*Hz*Npoints;
    E.c(n)=-(flux_phase_c(n+1)-flux_phase_c(n))*Hz*Npoints;
    E.ab(n)=-((flux_phase_a(n+1)-flux_phase_b(n+1))-(flux_phase_a(n)-flux_phase_b(n)))*Hz*Npoints;
    E.bc(n)=-((flux_phase_b(n+1)-flux_phase_c(n+1))-(flux_phase_b(n)-flux_phase_c(n)))*Hz*Npoints;
    E.ca(n)=-((flux_phase_c(n+1)-flux_phase_a(n+1))-(flux_phase_c(n)-flux_phase_a(n)))*Hz*Npoints;
     figure(133)
     title('BEMF')
         plot(E.a)
         hold on
         plot(E.b)
         hold on 
         plot(E.c)
         hold on
         xlabel('elec. angle[deg]')
         ylabel('BEMF')
         grid on

end
average_T = mean(Output.T);
max_T = max(Output.T);
min_T = min(Output.T);
Torque_Riffle = (max_T - min_T)/average_T

E.elec_angle=elec_angle;
% E.k_e=max(E.a)/(rpm/60/2/pi)
E.k_e=max(E.a)/(rpm/1000)
% closefemm
end