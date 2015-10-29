function CalciumTotal(FreeCalcium)
%Enter the desired Free Calcium concentration in uM
%Enter the buffer concentration
Bt=1000;            % BAPTA concentration in uM (micromoles)
Ca = FreeCalcium;    % in uM (micromoles)       
Kd=.14;             % in uM (micromoles); 0.14 is the Kd for BAPTA ; 0.067 for EGTA

CaTotal = (Kd*Ca+Bt*Ca+Ca^2)/(Ca+Kd);    % in uM (micromoles)
CaVolume = CaTotal/10;
CaNeeded = ['Add ' num2str(CaTotal) ' uM of Calcium']
VolumeNeeded = ['This is equivalent to ' num2str(CaVolume) ' uL of 1M CaCl2 for 100mL of solution']
