from Dusty import DustySpectrum
import numpy as np
foo = DustySpectrum(np.linspace(0.1,100,500),
                    spectralshape="blackbody",
                    spectralshape_blackbody_temperatures=5000,
                    spectralscale_inner_temperature=1500.,
                    density_distribution="powerlaw",
                    density_transition_radii=[1,100],
                    density_powers=-2,
                    grain_fractional_abundances=[1,0,0,0,0,0],
                    tau_min=10,
                    tau_max=10
                    
)

print(foo.geometry)

foo.dusty_write_input_file()
