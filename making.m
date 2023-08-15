% -------------------------------------------------------------------------
%    srpm to fit into de optimization
% -------------------------------------------------------------------------

% ---> begin of the 1st part.
lerror = -Inf;
     
% ----------------------------------------------------------------------
%    stator
% ----------------------------------------------------------------------
Num_turns=coil.Num_turns;
Num_series=coil.Num_series;
Num_parallel=coil.Num_parallel;
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
angle_offset = geo.angle_offset;
w_teeth = geo.W_t;                   %% teeth width 
shoe_1=geo.shoe1;                   
shoe_2=geo.shoe2;

A_coil_width=geo.A_coil_width;      
A_coil_height=geo.A_coil_height;     
one_slot_area=A_coil_width*A_coil_height*10^-6;    %%10^-6은 왜? 단위?
coil.current_density=current_density;
coil.Num_turns=Num_turns;
coil.Num_parallel=Num_parallel;
coil.Num_series=Num_series;
coil.coil_diameter=coil_diameter;
coil.Num_stranded_wire=Num_stranded_wire;
coil.length_tot=length_tot;
coil_height=geo.A_coil_height;
th_core=geo.T_s;


R_phase=geo.R_phase;
length_tot=coil.length_tot;
g=geo.airgap;
l=geo.L_st;



ns = slot_number; 
S1 = ns;
r = geo.D_is;
us = geo.slot_ratio*2*pi*r/ns;

ds = shoe_1+shoe_2;
% lamt = 1 - lams;
hs = coil_height;
wt = w_teeth;   
wst = 2*pi*(r+th_core)/ns-wt;
wsb = 2*pi*(r+hs+th_core)/ns - wt; % <--- parallel tooth / tapered slot
wsa = 0.5*(wsb+wst);
rci = th_core+r+hs;
rco = rci + ds;

% -------------------------------------------------------------------------
%    geff
% -------------------------------------------------------------------------

sts = us;
tws = 2*pi*r/ns - sts;
alphas = atan(sts/(2*g));
ccs = g*(alphas*tan(alphas) - log (1/cos(alphas)));
geff = g*(sts+tws)/(sts+tws-ccs);

% -------------------------------------------------------------------------
%    Coil <--- check again
% 1st variables: nc N_Coil_Slot npar
% 2nd variables: na Lac Aac Ra Ib massa
% -------------------------------------------------------------------------

% na = (nc/npar)*N_Coil_Slot*ma*p; % series turns / phase
% nsfp = ns/(2*p); 
% nsct = nsfp - nst; 
% thetae = acos(lams);
% laz = pi*(r+ds+.5*hs)*nsct/ns;
% %le1 = laz/tan(thetae);
% le2 = laz/sin(thetae);
% Lac = 2*na*(l+2*le2); % length of series wire

% For concentrated winding
na = Num_turns*ns/3;


Lac = length_tot;
%Aac = hs*wsa*lama/((nc/npar)*N_Coil_Slot); % area of (parallel wires together)
Aac = pi*coil_diameter^2/4*Num_stranded_wire; 
Ib = I_rated_est;
%theoretical result
Ra=R_phase;
%%
% -------------------------------------------------------------------------
%    Rotor
% -------------------------------------------------------------------------

rag = r - g/2;
rro = r - g;
rri = lamrd*rro;
hm = ratm*(rro-rri);
hm1 = lamh*hm;
hm2 = (1-lamh)*hm;
rrm1 = rri + hm2 + 0.5*hm1 + (lammd1+lammd2)*(1-ratm)*(rro-rri);
rrm2 = rri + 0.5*hm2 + lammd2*(1-ratm)*(rro-rri);

alphm2 = ratm2*pi;
alpm2 = (1-ratm2)*pi/2;
ang02 = alpm2/p;
ang12 = alphm2/(4*p);
rrem2 = rrm2/cos(ang12);
len12 = locside(ang12,(rro-brw),rrem2);
ang22 = locangle(rrem2,(rro-brw),len12);
ang32 = pi - ang22 - ang12;
ang42 = 0.5*(pi/2 + ang22 + 2*ang12);
len22 = hm2/cos(ang22)/2;
len32 = hm2/sin(ang42)/2;
alpm2 = angle(pol2comp((rro-brw),ang02)-pol2comp(hm2/(2*cos(ang22)),ang02-ang22))*p;
brw2  = rro-hm2/2-abs(pol2comp((rro-brw),ang02)-pol2comp(hm2/(2*cos(ang22)),ang02-ang22));

ang11i = atan((rrm2*tan(ang12)-(rrm1-rrm2)/tan(ang42))/rrm1);
rrem1 = rrm1/cos(ang11i); 
ang31 = 2*pi-2*ang42-(pi/2-ang11i);
ang21 = losangle(rro-brw,ang31,rrem1);
ang11o = pi - ang31 - ang21;
alphm1 = 2*p*(ang11o+ang11i);
len11 = locside(ang11o,(rro-brw),rrem1); 
ratm1 = alphm1/pi;
ang01 = (1-ratm1)*pi/(2*p);
alpm1 = angle(pol2comp((rro-brw),ang01)-pol2comp(hm1/(2*cos(ang21)),ang01-ang21))*p;
brw1 = rro-hm1/2-abs(pol2comp((rro-brw),ang01)-pol2comp(hm1/(2*cos(ang21)),ang01-ang21));
ang41 = ang42;
len21 = hm1/cos(ang21)/2;
len31 = hm1/sin(ang41)/2;

%    changed for concentrated winding
%-------------------------------------------------------------
%    feasibility test by Original
% if ((p<1) || (ma<1) || (ns<24) || (ns>120) || (nst<0) || (nst>ma))
%      lerror=Inf; return;
% end;
%-------------------------------------------------------------
if (lamrd*(r-g)<rriMIN)
     lerror=Inf; return;
end;
if ((ratm2*(r-g)*pi/p <= (2*lammd1*(1-ratm)+(1+lamh)*ratm)*(r-g)*(1-lamrd)))
     lerror=Inf; return;
end;
%-------------------------------------------------------------
%    more feasibility test here
if (rrm1*tan(ang11i))<(crw+hm1/2+hm1/2/tan(ang41)), lerror=Inf; return; end
if (rrm2*tan(ang12) )<(crw+hm2/2+hm2/2/tan(ang42)), lerror=Inf; return; end
if len11<(len21+len31), lerror=Inf; return; end
if (hm1)<hmMIN, lerror=Inf; return; end
if (hm2)<hmMIN, lerror=Inf; return; end
if (lammd1*(1-ratm)*(rro-rri))<hmMIN, lerror=Inf; return; end; % dr2/g<1
if (lammd2*(1-ratm)*(rro-rri))<hmMIN, lerror=Inf; return; end; % dr3/g<1
%-------------------------------------------------------------

% -------------------------------------------------------------------------
%    Magenetic circuit (scanned)
% -------------------------------------------------------------------------

dalpm1 = pi/2 - alpm1;
dalpm2 = alpm1 - alpm2;
lm1 = 2*(rrm1*tan(ang11i) - crw ...
      + sqrt(rrem1^2 + (rro - brw)^2 - 2*rrem1*(rro - brw)*cos(ang11o)));
lm2 = 2*(rrm2*tan(ang12) - crw ...
      + sqrt(rrem2^2 + (rro - brw)^2 - 2*rrem2*(rro - brw)*cos(ang12)));
lry = 0.5*[0.5*(lm1+2*(brw+crw)),0.5*(lm1+lm2+4*(brw+crw)),2*rrm2*tan(alphm2/(4*p))];
lst = hs+ds+db/2;
lsy = (r+ds+hs+0.5*db)*[2*dalpm1,2*(dalpm2+dalpm1),pi-2*dalpm1]/(4*p);
Asy = db*l;
Am1 = (lm1+2*(crw+brw))*l;
Am2 = (lm2+2*(crw+brw))*l;
Aag = rro*l*pi/p;
Ag1 = rro*l*2*dalpm1/p;
Ag2 = rro*l*2*dalpm2/p;
Ag = 0.5*[Ag1,Ag2,rro*l*(pi-2*(dalpm1+dalpm2))/p];
Ast = lamt*Ag;
Ary = l*[(rro-rrm1-hm1/2),(rrm1-rrm2-hm1/2-hm2/2),(rrm2-rri-hm2/2)];
Ab1 = (brw+crw)*l;
Ab2 = (brw+crw)*l;
lb1 = len21;
lb2 = len22;
As = 2*pi*r*l/ns;
dalps = 2*pi*p/ns;

% -------------------------------------------------------------------------
%    Mass
% -------------------------------------------------------------------------
% 
% massa = 3*Lac*Aac*rhos;
% mcore = rhofe*l*pi*(rco^2-rci^2);			
% mtooth = rhofe*l*(ns*wt*hs+2*pi*r*ds-ns*ds*us);	
% mtooth = rhofe*l*ns*(wt*(hs+ds) + (ds+os)*(2*pi*r/ns-wt-us));
% marea = 2*p*(lm1*hm1+lm2*hm2);
% mmagnet = l*marea*rhom;
% mcorer = rhofe*l*(pi*(rro^2-rri^2)-marea);
% miron = mcore+mtooth+mcorer;
% Jm = (mcorer+mmagnet)*(rro^2+rri^2)/2;

% -------------------------------------------------------------------------
%    Winding Factor (scanned)
% -------------------------------------------------------------------------
% Changed to concentrated winding
%gama = 2*pi*p/ns;
%nssp = nsct/nsfp;
nssp = 2/3;
alfa = pi*nssp;

% kb = sin (0.5 * ma * gama) / (ma * sin ( 0.5 * gama));
% kb5 = sin (2.5 * ma * gama) / (ma * sin (2.5 * gama));
% kb7 = sin (3.5 * ma * gama) / (ma * sin (3.5 * gama));

kp = sin (alfa/2.0);
kp5 = sin (5.0 * alfa/2.0);
kp7 = sin (7.0 * alfa/2.0);
	
% ks = sin (gama/2.0) / (gama/2.0);
% ks5 = sin (5*gama/2.0) / (5*gama/2.0);
% ks7 = sin (7*gama/2.0) / (7*gama/2.0);
	
% ka = kb * kp * ks;
% ka5 = kb5 * kp5 * ks5;
% ka7 = kb7 * kp7 * ks7;
ka = kp;
ka5 = kp5;
ka7 = kp7;

np = ns/p+1;
nm = ns/p-1; 
% if ma>1,
%   kbp = sin (0.5 * ma * np * gama) / (ma * sin ( 0.5 * np * gama));
%   kbm = sin (0.5 * ma * nm * gama) / (ma * sin ( 0.5 * nm * gama));
%   kpp = sin (np * alfa/2.0);
%   kpm = sin (nm * alfa/2.0);
%   ksp = sin (gama * np/2.0) / (gama * np/2.0);
%   ksm = sin (gama * nm/2.0) / (gama * nm/2.0);	
%   kap = kbp * kpp * ksp;
%   kam = kbm * kpm * ksm;
% end

% -------------------------------------------------------------------------
%    Leakage component (scanned)
% -------------------------------------------------------------------------

Lag = (3/2)*(4/pi)*muzero*na^2*ka^2*l*r/(p^2*geff);

perm = muzero*l*(ds/us+(1/3)*hs/wst);
Lslot = 2*p*nc^2*perm*(4*ma-nst);
if (nssp > (2/3)),
	knssp = 1.5*nssp - 0.5;
elseif (nssp < (1/3)),
	knssp = 1.5*nssp - 1;
else,
	knssp = 3*nssp - 1.5;
end
permn = muzero*l*((2*os/us + 2*(ds-os)*log(wsa/us)/(wsa-us) + (5/6)*(hs/wsa)) ...
        + knssp*(os/us + (ds-os)*log(wsa/us)/(wsa-us) + hs/(4*wsa)));
Lslot = na^2*permn/(2*ma*p);

Lend = 0.5*(140/(4 * pi^2))*(3/2) * muzero * r * na^2 * (alfa/pi - 0.3)/p^2;

La5 = Lag * (ka5 / (ka * 5))^2;
La7 = Lag * (ka7 / (ka * 7))^2;

% if ma>1,
%     Lap = Lag * (kap / (ka * np))^2;
%     Lam = Lag * (kam / (ka * nm))^2;
% else,
% 	Lap = 0;
% 	Lam = 0;
% end
Lap = 0;
Lam = 0;
Ll = Lslot + Lend + La5 + La7 + Lap + Lam;

% -------------------------------------------------------------------------
%    d-axis by Vagati <--- understand this
% -------------------------------------------------------------------------

rm1 = hm1*As/(geff*0.5*Am1*Mur);
rm1 = par(hm1*As/(geff*(0.5*Am1-(brw+crw)*l)*Mur),0.5*(brw+crw)*As/(geff*Ab1));
rm2 = hm2*As/(geff*0.5*Am2*Mur);
rm2 = par(hm2*As/(geff*(0.5*Am2-(brw+crw)*l)*Mur),0.5*(brw+crw)*As/(geff*Ab2));
rg1 = dalps/dalpm1;
rg2 = dalps/dalpm2;

puf1 = (cos(alpm1) - 0)/dalpm1;
puf2 = (cos(alpm2) - cos(alpm1))/dalpm2;

pur1 = ((rm1+par(rm2,rg2))/(rm1+rg1+par(rm2,rg2)))*puf1 ...
       + (rm2*rg1/(rm2+rm1+rg1)/(rg2+par(rm2,rm1+rg1)))*puf2;
pur2 = (par(rm2,rg2)/(rm1+rg1+par(rm2,rg2)))*puf1 ...
       + (par(rm2,rm1+rg1)/(rg2+par(rm2,rm1+rg1)))*puf2;
ratdcqm = 1 - (4/pi) * (dalpm1*puf1^2 + dalpm2*puf2^2);
ratdtqm = (4/pi) * (puf1*(puf1-pur1)*dalpm1 + puf2*(puf2-pur2)*dalpm2);

Laqm = Lag;
Ladt = Laqm*ratdtqm;
Ladc = Laqm*ratdcqm;
Ladm = Ladt + Ladc;
Ld = (Ladm + Ll)*(FuzzyLD);

% -------------------------------------------------------------------------

Rm1 = hm1/(muzero*Mur*Am1);
Rm1 = par(hm1/(muzero*Mur*(Am1-2*(brw+crw)*l)),(brw+crw)/(muzero*2*Ab1));
Rm1 = hm1/(muzero*Mur*(Am1-2*(brw+crw)*l));
Rm2 = hm2/(muzero*Mur*Am2);
Rm2 = par(hm2/(muzero*Mur*(Am2-2*(brw+crw)*l)),(brw+crw)/(muzero*2*Ab2));
Rm2 = hm2/(muzero*Mur*(Am2-2*(brw+crw)*l));
Rg1 = geff/(muzero*Ag1);
Rg2 = geff/(muzero*Ag2);

Areq1 = (hm2*Rg2+hm1*(Rm2+Rg1))/(Rm2*Rg2+(Rm2+Rg2)*(Rm1+Rg1))/muzero/Mur;
Areq2 = (hm2*(Rm1+Rg1)-hm1*Rm2)/(Rm2*Rg2+(Rm2+Rg2)*(Rm1+Rg1))/muzero/Mur;
Aseq1 = (-Rm2*Rg2/(Rm2*(Rg2+Rg1+Rm1)+Rg2*(Rg1+Rm1)))*Ab2 ...
        + (-(Rg2+Rm2)*Rm1/(Rg2*Rm2+(Rm1+Rg1)*(Rm2+Rg2)))*Ab1;
Aseq2 = (-Rm2*(Rm1+Rg1)/(Rm2*(Rg2+Rg1+Rm1)+Rg2*(Rg1+Rm1)))*Ab2 ...
        + (Rm2*Rm1/(Rg2*Rm2+(Rm1+Rg1)*(Rm2+Rg2)))*Ab1;
     
Bg1 = (Br*Areq1 + Bs*2*Aseq1)/Ag1;
Bg2 = (Br*Areq2 + Bs*2*Aseq2)/Ag2;
B1 = (4/pi) * (Bg1*sin(dalpm1) + Bg2*(sin(dalpm1+dalpm2) - sin(dalpm1)));

lambda = 2*rro*l*B1*na*ka/p;
Lambda=(lambda/sqrt(2))*(FuzzyPM);
% Lambda1 = (na*ka)*(rro*l/p)*2*B1/sqrt(2); % sinusoidal
% Lambda2 = (na*ka)*(rro*l/p)*2*(Bg1*dalpm1+Bg2*dalpm2)/sqrt(2); % staircase

% -------------------------------------------------------------------------
%    solve Q-axis
% -------------------------------------------------------------------------

Laqm = Lag;
Lq = Laqm + Ll;
Lql = Lq;
Sr = Lq/Ld;

if (mflag==1)
	Iqx = linspace(.001, Ib);
	Lqy = Lq*linspace(1.001,.999);
	Lqy(1) = Lq;
elseif (mflag==3)
     qaxmodel3;
else
     qaxmodel4;
     if lerror>0; lerror=Inf; return; end; % <--- Improved?
end

if(ismono(Iqx)==0);
	lerror = Inf;
	return;
end

[Lqmax,Indmax] = max(Lqy);
iqmax = Iqx(Indmax)/Ib;
Indmin = length(Lqy);
while((Lqy(Indmin) < Ld) & (Indmin > Indmax))
	Indmin = Indmin - 1;
end
iqmin = min(1,Iqx(Indmin)/Ib);
if ((iqmax > 1) | (Indmin <= Indmax))
	lerror = Inf;
	return;
end

Sr = max(Lqy)/Ld;
      
% -------------------------------------------------------------------------
%    End of the 1st part
% -------------------------------------------------------------------------

Tb = 3*p*Lambda*Ib;
xd = Ld*Ib/Lambda;
lerror=max([lerror,Jm/JmMax-1]);
lerror=max([lerror,(max(speedm)*pi/30*p*Lambda)/BackEMFScreen-1]);
% ---> end of the 1st part.
return;


% -------------------------------------------------------------------------
%    The second phase of analysis
% -------------------------------------------------------------------------

Preq=Prat;

Tatspecifiedspeed=Trat;

f = p*rpm/60;
om = 2*pi*p*rpm/60;
omegam = 2*pi*rpm/60;
stip = rro*om/p;
Ea = om*Lambda;
Xd = om*Ld;
Xq = om*Lq;

if (gen == 0) % Motoring
	TAGM1;
else % Generating
     PAmax;
end % of mot/gen
lerror=max([lerror,id,ia-1,Va/mvolt-1.001]);
lerror=max([lerror,(Ia/Aac)/JaMaxm(spind)-1]);
lerror=max([lerror,Ia/IaScreen(spind)-1]);
lerror=max([lerror,gamidqMin/abs(gamidq)-1]);
IaMSave(spind)=Ia; IqMSave(spind)=Iq; IdMSave(spind)=Id;

% Coercive Force Estimate
psidelt = atan(Xq/Ra) + gen*pi;
Iashort = Ea/(Xd*sin(psidelt) + Ra*cos(psidelt));
Idshort = Iashort*sin(psidelt);
Cf = Xd*Idshort*Br/Ea;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if gen==0
       Bagv = abs(Va-((Ia*(cos(psi)+j*sin(psi)))*(Ra+j*om*Ll)))*p/(om*2*rro*l*na*ka); 
elseif gen==1
       Bagv = (Va+Ia*(Ra+om*Ll))*p/(om*2*rro*l*na*ka); % Original
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Bagva = 2*sqrt(2)*Bagv/pi;
Bt = Bagva/lamt;
Bb = Bagva*pi*rro/(2*p*db);
Brc = Bagva*pi*rro/(2*p*(rro-rri-hm1-hm2));
BavgmaxM(spind) = max([Bt Bb Brc]);
if (BavgmaxM(spind)>BmaxScreen), lerror=1; return; end;

% -------------------------------------------------------------------------
%    Current,Power,Torque,Efficiency,PowerFacter,ShearStress
% -------------------------------------------------------------------------

Iq = Ia*cos(psi+delt);
Id = Ia*sin(psi+delt);

% % Core Iron Loss  (account for the lamination factor)
% 
% cored = 0.5e-3;     % Thickness of lamination(m)
% cond = 1.46e6;          % Conductivity(S/m)
% Vt = (wt*(hs+ds+os)+ds*(ds+os))*l*ns;
% Vy = pi*((dout/2)^2-(rro+g+ds+hs)^2)*l;
% % Eddy current loss
% Pey = (pi^2*cored^2*cond)/6*(Bb^2*f^2)*Vy;
% Pet = (pi^2*cored^2*cond)/6*(Bt^2*f^2)*Vt;
% % Hysteresis loss
% kh = 40;
% beta = 1.93;
% % Teeth
% Pht = kh*om*Bt^beta*Vt;
% % Yoke
% Phy = kh*om*Bb^beta*Vy;
% % total Iron loss
% Pc = Pet+Pey+Pht+Phy;


% Core Iron Loss
Pcb = mcore*pb*abs(Bb/bb)^epb*abs(om/omb)^epf;
Pct = mtooth*pb*abs(Bt/bb)^epb*abs(om/omb)^epf;
Pc = Pcb+Pct;





% Harmonic Loss
sigc = 1.46e6;
theta0 = us*ns/r;
Bh = Bagv*(2/pi)*sin(theta0/2);
omh = om*ns/p;
Ez = rro.*omh.*Bh/p;
B0 = 0.75*Bs;
deltc = (3.79*B0./(sigc^2.*omh.*Ez)).^.3333;
Ph = 1.75.*rro.*l.*sigc.*deltc.*Ez.^2;

pwind = 0;
pfan = 0;

% Armature loss
Pa = 6*Ia^2*Ra;
if (gen == 0),
	Pin = Preq+Pc+Pa+Ph+pwind+pfan;
	if (rpm == 0),
		Tm = 3*p*Lambda*Ia*cos(psi);
	else,
		Tm = Preq/omegam;
	end
else,
	Pin = Preq-Pc-Pa-Ph-pwind-pfan;
	Tm = Pin/omegam;
end
eff = Preq/Pin;

%-------------------------------------------------------
pfact = cos(psi);
% if (eff < Effm(spind)), lerror=1; return; end;

shear = abs(Tm)/(2*pi*rro^2*l);
