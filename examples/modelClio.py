# This routine is written to reproduce the model used by Clio Gielen.
# The model is described by Gielen et al. (2008, A&A 490, 725) for post-AGB
# stars. It has also been used by Gielen et al. (2011, A&A 533, A99) to
# model and compare Galactic and Magellanic post-AGB stars.
# An erratum was issued (Gielen et al. 2010, A&A 515, C2), but this seems to
# be only for the code which contained a bug.
#
# The model emission is given by:
#
# F_lambda ~ ( sum_i alpha_i * kappa_i ) x ( sum_j beta_j B_lambda(T_j) )
# Note: this formula as given in the paper is actually not exactly right, see description below!
#
# with
# kappa_i = mass absorption coefficient of dust component i
# alpha_i = the (mass?) fraction of dust component i
# B_lambda (T_j) = the Planck function at temperature T_j
# beta_j = the (mass?) fraction of dust at temperature T_j
#
# The temperature is assumed to be the independent of grain size and grain shape
#
# Constraints on the parameters used by Gielen et al. 2008
# Only two temperature components: j=1 and j=2, between 100 and 1000 K, with
# their relative fractions
# Four silicate species (amorphous and crystalline olivine and pyroxene),
# each with two dust sizes, and their relative fractions (7 free parameters)
# Thus we get a total of 7 * 2 (cold and warm) + 1 (relative fraction
# between cold and warm) = 15 free parameters. 


