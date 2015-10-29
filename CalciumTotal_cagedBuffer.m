function CalciumTotal_cagedBuffer(FreeCalcium)
%Enter the desired Free Calcium concentration in nM
%Enter the buffer concentration below
Bt=20;                  % NP-EGTA concentration in mM (millimoles)
Ca = FreeCalcium*10^-6; % in nM (micromoles)       
Kd=80e-6;               % in mM (micromoles); Kd = 80nM, which is 80e-6 mM

CaTotal = (Kd*Ca+Bt*Ca+Ca^2)/(Ca+Kd);    % in mM (millimoles)
CaNeeded = ['Add ' num2str(CaTotal) ' mM of Calcium']