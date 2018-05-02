# FreeFluidsC
Physical properties calculation using EOS and correlations.
Cubic EOS:
Implements the Peng Robinson (PR) and Soave Redlich Kwong (SRK) types, 
with different theta function calculations, and volume traslation. An special 
Peng Robinson type (PRFIT4) uses fitted parameters, not related to real Tc and Pc, 
for a,b,w and Tc, obtaining good results in density,Vp. and Cp calculations. 
At the moment VdW, and Panagiotopoulos Reid, mixing rules are working. 
The Mathias Klotz Prausnitz one needs some revision, and gE rules need an update. 
PCSAFT eos:
Includes association term. Polar PCSAFT is not working properly. 
CR-1 (Wolbach and Sandler modified) is used as mixing rule for mixtures.
Multiparameter eos as per Span and Wagner:
Includes IAWPS95.The ideal gas contribution is implemented,except for IAPWS95, 
in a non standard way, using ideal gas Cp instead of ideal Helmholtz energy, 
and using always as reference state that of an ideal gas. No work has been done in mixtures.
Correlations:
All the normal ones, mainly DIPPR.
Optimization tools:
For parameter finding, both for physical properties correlations and cubic 
and PCSAFT eos.
