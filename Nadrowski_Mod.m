function [T,Y]=Nadrowski_Mod(Amp_Input,Freq_Input,delta_t,t_end)
%final version, Chagit Braiman
%delta_t is the change in step for integration, t_end is the end time for
%the simulation, so for example in below, delta_t is .0025 seconds, and
%t_end is 10.24 seconds
%[T,Y] = ode45(@Hair_Bundle_ODE45,[0:.0025:10.24],[0 0 0 0],options,Param_Input);

T(1:t_end./(delta_t))=0;                 %preallocate array
Y(1:t_end./(delta_t),1:4)=0;             %preallocate array

options = odeset('RelTol',1e-12,'AbsTol',[1e-12 1e-12 1e-12 1e-12]);
[T,Y] = ode45(@Hair_Bundle_ODE45,[0:delta_t:t_end],[0 0 0 0 1],options,Amp_Input,Freq_Input);

    
    function dy = Hair_Bundle_ODE45(t,y,Amp_Input,Freq_Input)

    lambda=28*10^-7;                               %friction coefficient of hair bundle
	lambda_a=10*10^-6 ;                            %friction coefficient of adaptation motors
	K_gs= 750*10^-6;                               %combined gating-spring stiffness
	K_sp=600*10^-6;                                %combined stiffness of stereociliary pivots and loads
	gamma=(1./7);                                  %gamma, geometrical gain of stereociliary shear motion
	tau=1*10^-4;                                   %time constant of calcium feedbac5
	C_0=0;                                         %Intracellular calcium concentration with channels closed
	N=50;                                           %Number of stereocilia
    k_on_re =   2e9;                                % binding constant for Ca-release-element
    k_off_re=   75e3;                               % unbinding constant for Ca-release-element
    
	Boltzman_Constant=1.3806503*10^-23;             %Boltzman's constant k_b
	Temp=300;                                       %temperature
	delta_G=10*Boltzman_Constant*Temp;              %intrinisic energy change on channel opening        
	C_M=.1*10^-3;                                   %C_M for equation of calcium, value from Daibhid's supplemental info 


%/**********************************Now Define Model Parameters Initialized in Other Pages*********************************************/

	 D=61*10^-9;                                   %displacement of gating spring   (page 12195)
	 A=exp( (delta_G+((K_gs*D*D)./(2*N)))./(Boltzman_Constant*Temp));     %intrinsic energy difference between open and closed states of a transduction channel (page 12195)
     
     delta_Po_term = (N*Boltzman_Constant*Temp)/(K_gs*D);                         %delta for the P_0 term (equation 1)  page 12195
   
     
%/******************************************************Here using Nadrowski's Thesis***************************************************/

	f_max=550*(10^-12);                              %//term needed in equation 2.31,  value given in Appendix B of thesis, column G
	S=.88;                                 %// term used in equation 2.31, value given in Appendix B of thesis, column B 
	K_E=0;                            %//  Term used in equation 2.31, for now leaving out second term by saying K_E = 0 (I think this account for incomplete adaptation
	X_E=0;                            %//  Term used in equation 2.31, initialized in equation 2.34, but I don't know what their value of small xe is
    
	X_Noise=0;
	Xa_Noise=0;
	Calcium_Noise=0;

dy = zeros(5,1);   


Amp=Amp_Input * 10^-12; %.2*10^-12;                              %this is to cycle through freq with fixed
freq=Freq_Input;

dy(1)=(1/lambda)*(-K_gs(y(5))*(y(1)-y(2)-D*PO(y(1),y(2),A,delta_Po_term))-K_sp*y(1)+Amp*sin(2*pi*freq*y(4))+X_Noise); 
dy(2) =(1/lambda_a)*(K_gs(y(5))*(y(1)-y(2)-D*PO(y(1),y(2),A,delta_Po_term))+(gamma*f_max*  ((S*(y(3)/C_M))-1))+Xa_Noise);
dy(3) = (1/tau)*(C_0-y(3)+C_M*PO(y(1),y(2),A,delta_Po_term)+Calcium_Noise);
dy(4)=1;
dy(5) = k_on_re * Ca(t) * (1 - y(5)) - k_off_re * y(5);



function [K_gs] = K_gs(y5)
K_reMax = 1600e-6;      % Max stiffness of reclosure element
K_reMin = 200e-6;       % Min stiffness of reclosure element
K_tl    = 4000e-6;      % stiffness of tip link
K_re = (1 - y5)*(K_reMax - K_reMin) + K_reMin;
K_gs = (K_tl * K_re)/(K_tl + K_re);

function[POValue] =PO(X,Xa,A,delta_Po_term)
            POValue=1./(1+(A*exp(((-X+Xa)/delta_Po_term))));
            
function [Ca] = Ca(t)
    Ca = 1 + square(t);

