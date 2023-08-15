function  [seq,one_slot_area]=Draw_ebikes_newdq(coil,geo)
%%
Num_turns=coil.Num_turns;
Num_series=coil.Num_series;
Num_parallel=coil.Num_parallel;

%%파일저장
openfemm;
newdocument(0);
path='C:/femm42/';
name_fem='ebike.fem';

%%%% 재료 불러오기
mi_saveas([path,name_fem]);
Hc=982000;    %% 보자력
idensity = 0;  %%??
u_r_mag=1.05;   %%자석 비투자율
u_r_epoxy=1;   
mi_addmaterial('NXX',u_r_mag,u_r_mag,Hc,0,0.56)        %%자석물성이게다? 포화등은?
mi_addmaterial('NoMag',u_r_epoxy,u_r_epoxy)
mi_getmaterial('Air')
mi_getmaterial('M-19 Steel')                         %% 라이브러리?
mi_getmaterial('Pure Iron')
mi_addmaterial('a+', 1, 1, 0, -idensity/2, 56)
mi_addmaterial('a-', 1, 1, 0, idensity/2, 56)
mi_addmaterial('b+', 1, 1, 0, idensity, 56)
mi_addmaterial('b-', 1, 1, 0, -idensity, 56)
mi_addmaterial('c+', 1, 1, 0, -idensity/2, 56)
mi_addmaterial('c-', 1, 1, 0, idensity/2, 56)
mi_addboundprop('A_0',0,0,0,0,0,0,0,0,0)          %%도화지의 태두리 기본 설정


stator_yoke_thickness = geo.T_s;     %% 고정자 요크두께
rotor_yoke_thickness=geo.T_r;        %% 회전자 요크 두께

% add parameter    %%form Ebike_SPMSM
depth = geo.L_st;                  %% stack length 
pole_number= geo.pole_number;      %%극수
slot_number = geo.slot_number;     %%슬롯수
airgap = geo.airgap;               %%공극 
m_thickness=geo.T_m;               %%자석 두께

% pole_arc_ratio = geo.pole_arc_ratio;  %%spmsm 에서 설정할땐 0.99? 입력이 어디
% >> pole arc 와 어떤 것의 비율?
pole_arc_ratio=geo.pole_arc_ratio;            
slot_ratio = geo.slot_ratio;   
D_stator_outer =  geo.D_os;                  
D_rotor_outer= geo.D_or;
D_rotor_inner=D_rotor_outer-m_thickness*2;
D_stator_inner=geo.D_is;
angle_offset = geo.angle_offset 
w_teeth = geo.W_t;                   %% teeth width 
shoe_1=geo.shoe1;                   
shoe_2=geo.shoe2;

A_coil_width=geo.A_coil_width;      
A_coil_height=geo.A_coil_height;     
one_slot_area=A_coil_width*A_coil_height*10^-6;    %%10^-6은 왜? 단위?

D_rotor_outer=D_stator_outer+(airgap+rotor_yoke_thickness+m_thickness)*2;    
D_rotor_inner=D_stator_outer+airgap*2;

teeth_height=D_rotor_outer/2-rotor_yoke_thickness-m_thickness-airgap-D_stator_inner/2-stator_yoke_thickness; %% stator_yoke_thickness

m=0.1;n=0.1;

aa = 360/slot_number;
bb = 360/pole_number;
cc = tan(pi/slot_number)*(D_stator_inner/2+stator_yoke_thickness);       

mi_probdef(0,'millimeters','planar',1e-008,depth,30);  %% (freq, units, planar: 2D평면, precision, depth, minangle) 

q_denominator=pole_number*3/gcd(slot_number,pole_number*3);   %% 3: 3상  LPM에서도 나왔는데 이해가 잘?
q_numerator=slot_number/gcd(slot_number,pole_number*3);
num_zeros=q_denominator-q_numerator;
seq_zero= num_zeros/q_numerator;
sequence=[];                                                 %% sequence? 
num_repeat=slot_number/q_numerator;
temp_num_zeros=ceil(seq_zero-0.5);

count_zeros=0;
count_ones=0;

for i = 1: num_repeat*q_numerator
    count_ones=count_ones+1;
    sequence(count_ones+count_zeros)=1;
    for x = 1:temp_num_zeros
        count_zeros=count_zeros+1;
        sequence(count_ones+count_zeros)=0;
    end
    if ((count_zeros/count_ones)>seq_zero)
        temp_num_zeros=floor(seq_zero);
    elseif((count_zeros/count_ones)<seq_zero)
        temp_num_zeros=ceil(seq_zero);
    end
end

seq=[];
Dist=[  -11, 10, -12,11,-10, 12];
seq_count=1;
seq_temp=0;

for j=1:length(sequence)
    if (sequence(j)==1)
        seq(seq_count)=Dist(mod(j,6)+1);
        seq_count=seq_count+1;
        seq(seq_count)=seq(seq_count-1)*-1;
        seq_count=seq_count+1;
    end
end

mi_probdef(0,'millimeters','planar',1e-008,depth,30)    %% 왜 또?

% square=D_rotor_outer*0.5*2.2;
% mi_drawrectangle(-square/2,-square/2,square/2,square/2)

% x1=-square/2;
% y1=-square/2;
% x2=square/2;
% y2=square/2;
% mi_selectsegment((x1+x2)/2,y1);
% mi_selectsegment((x1+x2)/2,y2);
% mi_selectsegment(x1,(y1+y2)/2);
% mi_selectsegment(x2,(y1+y2)/2);
% mi_setsegmentprop('A_0', 0, 1,0,0)
% mi_clearselected()
% mi_zoomnatural()

%mi_addnode(0,D_rotor_outer/2-rotor_yoke_thickness)
%mi_addnode(0,D_rotor_outer/2-rotor_yoke_thickness-m_thickness)
%mi_addnode((D_rotor_outer/2-rotor_yoke_thickness)*cos(pi/180*(-bb*pole_arc_ratio/2+90)),(D_rotor_outer/2-rotor_yoke_thickness)*sin(pi/180*((-bb*pole_arc_ratio/2+90))))
%%%%% mi_addnode((D_rotor_outer/2-rotor_yoke_thickness-m_thickness)*cos(pi/180*(-bb*pole_arc_ratio/2+90+m*180/pi/(D_rotor_outer/2-rotor_yoke_thickness-m_thickness))),(D_rotor_outer/2-rotor_yoke_thickness-m_thickness)*sin(pi/180*((-bb*pole_arc_ratio/2+90+m*180/pi/(D_rotor_outer/2-rotor_yoke_thickness-m_thickness)))))
%mi_addnode((D_rotor_outer/2-rotor_yoke_thickness-m_thickness)*cos(pi/180*(-bb*pole_arc_ratio/2+90)),(D_rotor_outer/2-rotor_yoke_thickness-m_thickness)*sin(pi/180*((-bb*pole_arc_ratio/2+90))))
%mi_addarc((D_rotor_outer/2-rotor_yoke_thickness-m_thickness)*cos(pi/180*(-bb*pole_arc_ratio/2+90+m*180/pi/(D_rotor_outer/2-rotor_yoke_thickness-m_thickness))),(D_rotor_outer/2-rotor_yoke_thickness-m_thickness)*sin(pi/180*((-bb*pole_arc_ratio/2+90+m*180/pi/(D_rotor_outer/2-rotor_yoke_thickness-m_thickness)))),0,D_rotor_outer/2-rotor_yoke_thickness-m_thickness,bb/2*pole_arc_ratio-m*180/pi/(D_rotor_outer/2-rotor_yoke_thickness-m_thickness),5)
%mi_addsegment((D_rotor_outer/2-rotor_yoke_thickness)*cos(pi/180*(-bb*pole_arc_ratio/2+90)),(D_rotor_outer/2-rotor_yoke_thickness)*sin(pi/180*((-bb*pole_arc_ratio/2+90))),(D_rotor_outer/2-rotor_yoke_thickness-m_thickness+m)*cos(pi/180*(-bb*pole_arc_ratio/2+90)),(D_rotor_outer/2-rotor_yoke_thickness-m_thickness+m)*sin(pi/180*((-bb*pole_arc_ratio/2+90))))

%%%%%%%%%         rotor 그리기       %%%%%%%%%%

r=D_rotor_outer/2-rotor_yoke_thickness;      %% 식이길어져서 밑에 있던 r= 위로 올림
mi_addnode(0,r)                             %%자석 그리기
mi_addnode(0,r-m_thickness)
mi_addnode((r)*cos(pi/180*(-bb*pole_arc_ratio/2+90)),(r)*sin(pi/180*((-bb*pole_arc_ratio/2+90))))
% mi_addnode((D_rotor_outer/2-rotor_yoke_thickness-m_thickness)*cos(pi/180*(-bb*pole_arc_ratio/2+90+m*180/pi/(D_rotor_outer/2-rotor_yoke_thickness-m_thickness))),(D_rotor_outer/2-rotor_yoke_thickness-m_thickness)*sin(pi/180*((-bb*pole_arc_ratio/2+90+m*180/pi/(D_rotor_outer/2-rotor_yoke_thickness-m_thickness)))))
mi_addnode((r-m_thickness)*cos(pi/180*(-bb*pole_arc_ratio/2+90)),(r-m_thickness)*sin(pi/180*((-bb*pole_arc_ratio/2+90))))
mi_addarc((r-m_thickness)*cos(pi/180*(-bb*pole_arc_ratio/2+90+m*180/pi/(r-m_thickness))),(r-m_thickness)*sin(pi/180*((-bb*pole_arc_ratio/2+90+m*180/pi/(r-m_thickness)))),0,r-m_thickness,bb/2*pole_arc_ratio-m*180/pi/(r-m_thickness),5)
mi_addsegment((r)*cos(pi/180*(-bb*pole_arc_ratio/2+90)),(r)*sin(pi/180*((-bb*pole_arc_ratio/2+90))),(r-m_thickness+m)*cos(pi/180*(-bb*pole_arc_ratio/2+90)),(r-m_thickness+m)*sin(pi/180*((-bb*pole_arc_ratio/2+90))))
mi_zoomnatural()
% mi_addarc((D_rotor_outer/2-rotor_yoke_thickness-m_thickness)*cos(pi/180*(-bb*pole_arc_ratio/2+90+m*180/pi/(D_rotor_outer/2-rotor_yoke_thickness-m_thickness))),(D_rotor_outer/2-rotor_yoke_thickness-m_thickness)*sin(pi/180*((-bb*pole_arc_ratio/2+90+m*180/pi/(D_rotor_outer/2-rotor_yoke_thickness-m_thickness)))),(D_rotor_outer/2-rotor_yoke_thickness-m_thickness+m)*cos(pi/180*(-bb*pole_arc_ratio/2+90)),(D_rotor_outer/2-rotor_yoke_thickness-m_thickness+m)*sin(pi/180*((-bb*pole_arc_ratio/2+90))),90,5)
mi_selectcircle(0,0,D_rotor_outer/2,4)
mi_mirror2(0,0,0,10,4)
mi_selectcircle(0,0,D_rotor_outer/2,4)
mi_copyrotate2(0,0,bb,pole_number-1,4)

x1=0;y1=0;                                      %% x1 y1 은 어떨때 ?
mi_addnode(x1+r,y1)
mi_addnode(x1-r,y1)
mi_addarc(x1+r,y1,x1-r,y1,180,5)
mi_addarc(x1-r,y1,x1+r,y1,180,5)

x1=0;y1=0;               %% 로터 YOKE 그리기
k=D_rotor_outer/2;
mi_addnode(x1+k,y1)
mi_addnode(x1-k,y1)
mi_addarc(x1+k,y1,x1-k,y1,180,5)
mi_addarc(x1-k,y1,x1+k,y1,180,5)

mi_clearselected();
mi_selectarcsegment(0,k);
mi_setarcsegmentprop(1,'A_0', 0,0);
mi_selectarcsegment(0,-k) ;
mi_setarcsegmentprop(1,'A_0', 0,0);

mi_zoomnatural()

mi_selectcircle(0,0,D_rotor_outer/2,4)
mi_moverotate(0,0,180/pole_number)

for z=1:pole_number/2             %% 자석 재질 입력
    mi_clearselected();
    mi_addblocklabel(cos(pi/180*(90+180/pole_number+(z-1)*720/pole_number))*(D_rotor_inner+m_thickness)*0.5,sin(pi/180*(90+180/pole_number+(z-1)*720/pole_number))*(D_rotor_inner+m_thickness)*0.5);
    mi_selectlabel(cos(pi/180*(90+180/pole_number+(z-1)*720/pole_number))*(D_rotor_inner+m_thickness)*0.5,sin(pi/180*(90+180/pole_number+(z-1)*720/pole_number))*(D_rotor_inner+m_thickness)*0.5);
    mi_setblockprop('NXX',1,0,'None',90+180/pole_number+(z-1)*720/pole_number,(10+z*2-1),1)
    mi_clearselected();
    mi_addblocklabel(cos(pi/180*(90-180/pole_number+(z-1)*720/pole_number))*(D_rotor_inner+m_thickness)*0.5,sin(pi/180*(90-180/pole_number+(z-1)*720/pole_number))*(D_rotor_inner+m_thickness)*0.5);
    mi_selectlabel(cos(pi/180*(90-180/pole_number+(z-1)*720/pole_number))*(D_rotor_inner+m_thickness)*0.5,sin(pi/180*(90-180/pole_number+(z-1)*720/pole_number))*(D_rotor_inner+m_thickness)*0.5);
    mi_setblockprop('NXX',1,0,'None',-90-180/pole_number+(z-1)*720/pole_number,(10+z*2),1);
end

% mi_selectcircle(0,0,(D_rotor_outer*0.5+airgap*0.5),4)
% for m=11:10+pole_number*2
%     mi_selectgroup(m);
% end
% mi_moverotate(0,0,-180/pole_number)

%%%%        stator 그리기          %%%%%%%%%%%
mi_addnode(cc,D_stator_inner/2+stator_yoke_thickness) 
mi_addnode(w_teeth/2,D_stator_inner/2+stator_yoke_thickness)
mi_addsegment(cc,D_stator_inner/2+stator_yoke_thickness,w_teeth/2,D_stator_inner/2+stator_yoke_thickness)
mi_addnode(0,D_stator_inner/2+stator_yoke_thickness+teeth_height)
mi_addnode((D_stator_inner/2+stator_yoke_thickness+teeth_height)*cos(pi/180*(90-180/slot_number*slot_ratio)),(D_stator_inner/2+stator_yoke_thickness+teeth_height)*sin(pi/180*(90-180/slot_number*slot_ratio)))
mi_addarc((D_stator_inner/2+stator_yoke_thickness+teeth_height)*cos(pi/180*(90-180/slot_number*slot_ratio)),(D_stator_inner/2+stator_yoke_thickness+teeth_height)*sin(pi/180*(90-180/slot_number*slot_ratio)),0,D_stator_inner/2+stator_yoke_thickness+teeth_height,(180/slot_number*slot_ratio),5)
mi_addnode((D_stator_inner/2+stator_yoke_thickness+teeth_height-shoe_1)*cos(pi/180*(90-180/slot_number*slot_ratio)),(D_stator_inner/2+stator_yoke_thickness+teeth_height-shoe_1)*sin(pi/180*(90-180/slot_number*slot_ratio)))
mi_addsegment((D_stator_inner/2+stator_yoke_thickness+teeth_height)*cos(pi/180*(90-180/slot_number*slot_ratio)),(D_stator_inner/2+stator_yoke_thickness+teeth_height)*sin(pi/180*(90-180/slot_number*slot_ratio)),(D_stator_inner/2+stator_yoke_thickness+teeth_height-shoe_1)*cos(pi/180*(90-180/slot_number*slot_ratio)),(D_stator_inner/2+stator_yoke_thickness+teeth_height-shoe_1)*sin(pi/180*(90-180/slot_number*slot_ratio)))

mi_addnode(w_teeth/2,D_stator_outer*0.5-shoe_1-shoe_2)
mi_addsegment(w_teeth/2,D_stator_outer*0.5-shoe_1-shoe_2,(D_stator_inner/2+stator_yoke_thickness+teeth_height-shoe_1)*cos(pi/180*(90-180/slot_number*slot_ratio)),(D_stator_inner/2+stator_yoke_thickness+teeth_height-shoe_1)*sin(pi/180*(90-180/slot_number*slot_ratio)))
mi_addsegment(w_teeth/2,D_stator_outer*0.5-shoe_1-shoe_2,w_teeth/2,D_stator_inner/2+stator_yoke_thickness)

mi_addnode(w_teeth/2+A_coil_width,D_stator_inner/2+stator_yoke_thickness)
mi_addnode(w_teeth/2+A_coil_width,D_stator_inner/2+stator_yoke_thickness+A_coil_height)
mi_addnode(w_teeth/2,D_stator_inner/2+stator_yoke_thickness+A_coil_height)
mi_addsegment(w_teeth/2+A_coil_width,D_stator_inner/2+stator_yoke_thickness+A_coil_height,w_teeth/2+A_coil_width,D_stator_inner/2+stator_yoke_thickness)
mi_addsegment(w_teeth/2+A_coil_width,D_stator_inner/2+stator_yoke_thickness+A_coil_height,w_teeth/2,D_stator_inner/2+stator_yoke_thickness+A_coil_height)

mi_selectcircle(0,0,D_stator_outer*0.5+airgap*0.5,4)
mi_mirror2(0,0,0,10,4)
mi_selectcircle(0,0,D_stator_outer*0.5+airgap*0.5,4)
mi_copyrotate2(0,0,aa,slot_number-1,4)



%%%%% shaft? 뭐라고 해야하나 그리기?

x1=0;y1=0;
r=D_stator_inner*0.5;                %% 위에 변수랑 같아서 헷갈릴수도? 
mi_addnode(x1+r,y1)
mi_addnode(x1-r,y1)
mi_addarc(x1+r,y1,x1-r,y1,180,5)
mi_addarc(x1-r,y1,x1+r,y1,180,5)


%%%%%%코일 a b c 상 입력  %%%
seq = circshift(seq,42);
zz=1;

for v=1:length(seq)
    xx=floor((v-1)/2);
    zz=zz*-1;
    
    mi_addblocklabel((D_stator_inner/2+stator_yoke_thickness+A_coil_height/2)*cos(pi/180*(zz*120/slot_number*0.9+90+360/slot_number*xx)),(D_stator_inner/2+stator_yoke_thickness+A_coil_height/2)*sin(pi/180*(zz*120/slot_number*0.9+90+360/slot_number*xx)));
    mi_selectlabel((D_stator_inner/2+stator_yoke_thickness+A_coil_height/2)*cos(pi/180*(zz*120/slot_number*0.9+90+360/slot_number*xx)),(D_stator_inner/2+stator_yoke_thickness+A_coil_height/2)*sin(pi/180*(zz*120/slot_number*0.9+90+360/slot_number*xx)));
    if seq(v)==-10
        mi_setblockprop('b-',1, 0, '<None>', 0, 2, 0);
        mi_clearselected();
    elseif seq(v)==10
        mi_setblockprop('b+',1, 0, '<None>', 0, 2, 0);
        mi_clearselected();
    elseif seq(v)==-11
        mi_setblockprop('a-',1, 0, '<None>', 0, 2, 0);
        mi_clearselected();
    elseif seq(v)==11
        mi_setblockprop('a+',1, 0, '<None>', 0, 2, 0);
        mi_clearselected();
    elseif seq(v)==-12
        mi_setblockprop('c-',1, 0, '<None>', 0, 2, 0);
        mi_clearselected();
    elseif seq(v)==12
        mi_setblockprop('c+',1, 0, '<None>', 0, 2, 0);
        mi_clearselected();
    end
end

% rot_slot_num=29;
% mi_selectcircle(0,0,(D_stator_outer*0.5+airgap*0.5),4)
% mi_moverotate(0,0,360/slot_number*(slot_number-rot_slot_num))
%%%%%% stator 철심 재질입력
mi_clearselected();
mi_addblocklabel(0,D_stator_outer*0.5-1);
mi_selectlabel(0,D_stator_outer*0.5-1);
mi_setblockprop('M-19 Steel',1, 0, '<None>', 0, 1, 0);
%%%%% rotor 철심 재질입력
mi_clearselected();
mi_addblocklabel(0,D_rotor_outer*0.5-1);
mi_selectlabel(0,D_rotor_outer*0.5-1);
mi_setblockprop('M-19 Steel',1, 0, '<None>', 0, 10, 0);
%%% airgap 공기 재질 입력
mi_clearselected();
mi_addblocklabel(0,D_stator_outer*0.5+airgap*0.5);
mi_selectlabel(0,D_stator_outer*0.5+airgap*0.5);
mi_setblockprop('Air',1, 0, '<None>', 0, 0, 0);

% mi_clearselected();
% mi_addblocklabel(0,D_rotor_outer*0.5+1);
% mi_selectlabel(0,D_rotor_outer*0.5+1);
% mi_setblockprop('Air',1, 0, '<None>', 0, 0, 0);



 %%%%    전기각 맞추기 코드 추가     %%%%%%%%5
     mi_clearselected()
      mi_selectcircle(0,0,(D_stator_outer*0.5+airgap*0.5),4)
      mi_selectgroup(2)
      mi_moverotate(0,0,angle_offset)    
% %%%%%%%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%55555

%%

%%%%%% shaft?? 공기 재질 입력
mi_clearselected();
mi_addblocklabel(0,0);
mi_selectlabel(0,0);
mi_setblockprop('Air',1, 0, '<None>', 0, 0, 0);


mi_zoomnatural();
mi_clearselected()
mi_refreshview()

end

