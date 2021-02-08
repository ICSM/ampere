pro calculate_crosssection

readcol, '/home/zeegers/git_ampere/ampere/small_silicates/olivines/broadened_freq_in_cm-1/10/broadened_freq_cm_1', freq, int_abs_coeff

readcol, '/home/zeegers/git_ampere/ampere/small_silicates/olivines/broadened_freq_in_microns/10/broadened_freq_microns_1', wav, int_abs_coeff_wav


Navogadro=6.02214076e23 ; units per mol

conversion_value_eps_sig=alog(10.)/Navogadro

; putting int_abs_coeff in the right units
int_abs_coeff_cm_mol=int_abs_coeff*1./(42.2561)*100000.

cross_section=conversion_value_eps_sig*int_abs_coeff
cross_section=conversion_value_eps_sig*int_abs_coeff*10000.

; just a test to calculate Qext
Qext=cross_section/(1e-7)^2.*!pi




stop

end
