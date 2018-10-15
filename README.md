# FreeFluidsC
Physical properties calculation using EOS, activity coefficients, and correlations.

Cubic EOS:
Implements the Peng Robinson (PR) and Soave Redlich Kwong (SRK) types, 
with different theta function calculations, and volume traslation. An special 
Peng Robinson type (PRFIT4) uses fitted parameters, not related to real Tc and Pc, 
for a,b,w and Tc, obtaining good results in density,Vp. and Cp calculations.
For mixtures, VdW, Panagiotopoulos and Reid, and MKP mixing rules are implemented. 
Also gE mixing rules (HV,MHV1,MHV2,LCVM and UMR) based on activity models.

PCSAFT eos:
Includes PCSAFT with association term. Also Polar PCSAFT with dipolar term, both
Gross and Vrabeck and Jog and Chapman models. but only the Gross and Vrabeck is implemented for mixtures. 
CR-1 (Wolbach and Sandler modified) is used as mixing rule for mixtures.
SAFT VR with Mie potential is also implemented, but only for pure substances.

Multiparameter eos as per Schmidt and Wagner:
Includes IAWPS95.The ideal gas contribution is implemented,except for IAPWS95, 
in a non standard way, using ideal gas Cp instead of ideal Helmholtz energy, 
and using always as reference state that of an ideal gas. No work has been done in mixtures.

Activity coefficients:
Includes Wilson, NRTL, UNIQUAC and UNIFAC models. The last one comprising: standard, PSRK, 
Dortmund and NIST modifications. For polymers: Unifac ZM and Entropic FV models.
The files included in the data folder must be accesible to the program if UNIFAC models 
are to be used.

Correlations:
All the normal ones, mainly DIPPR.

Optimization tools:
For parameter finding, both for physical properties correlations and cubic 
and PCSAFT eos.
