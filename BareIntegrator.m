function [T,Y]=BareIntegrator(delta_t,t_end)

T(1:t_end/(delta_t))=0;                 %preallocate array
Y(1:t_end/(delta_t),1)=0;             %preallocate array

tspan = [0:delta_t:t_end];
start = [0.1 0];

[T,Y] = ode45(@BareDerivatives,tspan,start,[]);



function [dy] = BareDerivatives(t,y)
    k_on_re =   2e9;        % binding constant for Ca-release-element
    k_off_re=   75e3;       % unbinding constant for Ca-release-element
    dy = zeros(2,1);
    dy(1) = k_on_re * CaINmod4(y(2),1,2,30e-9,2e-6,1e-6) * (1 - y(1)) - k_off_re * y(1); 
    dy(2) = 1;

    
function [CaINmod4] = CaINmod4(t,onset,offset,I,r,bathCa)
region1 = t < onset;
CaINmod4(region1) = bathCa;

region2 = (t > onset) & (t < offset);
CaINmod4(region2) = bathCa +  3;

region3 =  t > offset;
CaINmod4(region3) = bathCa;



% function [Ca] = Ca(t)
%     bathCa = 1e-6;
%     I = 30e-9;
%     r = 2e-6;
%     Ca =  bathCa +  (((0.12 * I) / (4 * pi * 8e-10 * 1.602e-19 * r)) *   erfc(r./(2*(sqrt(8e-10*t)))) / 6.022e23 / 1000);
%         
