function [ dydt ] = kCaDeriv(t,y,k_on_re,k_off_re,k_on_m,k_off_m)
Ca = y(1);
RE = y(2);
M = y(3);
CaRE = y(4);
CaM = y(5);

dydt = [ -k_on_re * Ca * RE - k_on_m * Ca * M + k_off_re * CaRE + k_off_m * CaM;
        -k_on_re * RE + k_off_re * CaRE;
        -k_on_m * M + k_off_m * CaM;
        k_on_re * Ca * RE - CaRE * k_off_re;
        k_on_m * Ca * M - CaM * k_off_m;
        ];
    
        