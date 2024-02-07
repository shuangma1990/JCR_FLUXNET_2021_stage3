
max_infil=4
liquid_in=seq(0,40,0.1)
infil = max_infil*(1 - exp(-liquid_in/max_infil))

plot(liquid_in,infil)
