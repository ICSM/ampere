This directory contains the files associated with the fit presented in Kemper et al. (2002, A&A, 394, 679).

+ crystallinity.ia:
   - This is the IDL script used to make/plot the fit, according to the chi-by-eye principle. Note that this is not a standard IDL script, but it assumes that you run the IA3 data reduction suite, which was written for ISO-SWS data.

+ sh_calcaar.pro
+ ck_modbb.pro
+ pl.pro
+ lezen.pro
   - these are IDL/IA3 routines called by crystallinity.ia

+ NGC6302_nolines.tab
   - data file (wavelength (um), flux (Jy)) with the spectral lines removed.
+ NGC6302_100.tab
   - data file (wavelength (um), flux (Jy)) rebinned to a spectral resolution of 100.

+ calcite_improved.dat
+ dolomite_wh.mac
   - carbonate opacity files, in kappa (mass absorption coefficient), in units of g cm-2. I must say that I don't know how calcite_improved.dat is improved though.

+ am_oliv_CDE.dat
+ am_oliv.dat
+ cr_ice.dat
+ cr_oliv_CDE.dat
+ fe.dat
+ forst_m.q
    - various dust opacities calculated by modust from n,k values. these are in dimensionless Q_abs.

+ Koike1999_c_enst.q
+ Koike1999_o_enst.q
+ Koike2000_diopside.q_a
   - Q/a opacities measured by Koike et al. in the lab. Units cm-1.



