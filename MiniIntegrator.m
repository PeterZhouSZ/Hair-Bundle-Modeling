function [T,Y]=MiniIntegrator(delta_t,t_end)

T(1:t_end/(delta_t))=0;                 %preallocate array
Y(1:t_end/(delta_t),1:5)=0;             %preallocate array

tspan = [0:delta_t:t_end];
start = [0 0 0.01 .2 .2];

[T,Y] = ode45(@MiniDerivatives,tspan,start,[]);

    function [dy] = MiniDerivatives(t,y)
    k_sf    =   150e-6;     % fiber stiffness
    drag_hb =   130e-9;     % drag coeff of hair bundle
    drag_sf =   100e-9;     % drag coeff of stimulus fiber
    N_gs    =   35;         % number of gating elements (channels)
    gamma   =   0.14;       % geometrical gain
    x_c     =   12e-9;      % resting extension of gating spring
    d       =   7e-9;       % distance gating spring shortens on MTC opening
    k_sp    =   200e-6;     % stereociliary pivot stiffness
    X_sp    =   251e-9;     % resting position of pivots
    boltzmann=  1.3806503*10^-23; %Boltzmann constant
    temp    =   300;        % Kelvin temperature
    k_es    =   140e-6;     % stiffness of extent spring
    x_es    =   0e-6;       % resting deflection of extent spring
    k_on_m  =   10e9;       % binding constant for Ca-myosin
    k_off_m =   80e3;      % unbinding constant for Ca-myosin
    k_on_re =   2e9;        % binding constant for Ca-release-element
    k_off_re=   75e3;       % unbinding constant for Ca-release-element
    m_hb    =   60e-18;     % mass of hair bundle
    m_sf    =   100e-18;    % mass of stimulus fiber
    
    dy = zeros(5,1);

    dy(1) = y(2);
    dy(2) = ( -(k_sf * y(1)) - (drag_hb + drag_sf)*y(2) - N_gs * gamma * K_gs(y(5)) * (gamma * y(1) - y(3) + x_c - Po(y(1),y(3),y(5),d,gamma,x_c,boltzmann,temp)*d) - k_sp*(y(1)-X_sp)) / (m_hb + m_sf);
    dy(3) = -C(y(4)) + S(y(4))*((K_gs(y(5))*(gamma*y(1) - y(3) + x_c - Po(y(1),y(3),y(5),d,gamma,x_c,boltzmann,temp)*d)) - k_es*(y(3) + x_es));
    dy(4) = k_on_m  * Ca(t) * (1 - y(4)) - k_off_m  * y(4);
    dy(5) = k_on_re * Ca(t) * (1 - y(5)) - k_off_re * y(5);



function [K_gs] = K_gs(y5)
K_reMax = 1600e-6;      % Max stiffness of reclosure element
K_reMin = 200e-6;       % Min stiffness of reclosure element
K_tl    = 4000e-6;      % stiffness of tip link
K_re = (1 - y5)*(K_reMax - K_reMin) + K_reMin;
K_gs = (K_tl * K_re)/(K_tl + K_re);

function [Po] = Po(y1,y3,y5,d,gamma,x_c,boltzmann,temp)
delta_e = 68e-21;       % Intrinsic E change on MTC opening
Po = 1 / (1 +  exp( (delta_e - (K_gs(y5)*d*(gamma * y1 - y3 + x_c - d/2))) / (boltzmann*temp) ) );

    
function [C] = C(y4)
Cmax    =   0.06e-6;    % Max myosin climbing rate
Cmin    =   0;          % Min myosin climbing rate
C = (1 - y4)*(Cmax - Cmin) + Cmin;

function [S] = S(y4)
Smin    =   0;          % Min myosin slip rate
Smax    =   110e3;      % Max myosin slip rate
S = y4*(Smax - Smin) + Smin;

function [Ca] = Ca(t)
    Ca = 1 + square(t);
%{
Pca     =   10e-18;     % Ca permeability through channel
zCa     =   2;          % Valence of Ca ion
electron=   1.6021766e-19;% charge of an electron
Faraday =   96485.3;    % Coulombs/mol; Faraday's constant
Vm      =   -55e-3;     % membrane potential
CaOut   =   250e-6;     % [Ca] outside cell; 250uM = 250e-6 M
Dca     =   8e-10;      % diffusion coeff for Ca
rm      =   20e-9;      % distance from channel to myosin
Po = Po(y1,y3,d,gamma,x_c,boltzmann,temp);
Ica = (Po * Pca * zCa^2 * electron * Faraday * Vm * CaOut) / (boltzmann*temp*(1 - exp((zCa*electron*Vm)/boltzmann*temp)));
CaIN = -Ica / (2*pi*zCa*Faraday*Dca*rm);
    %}
    


