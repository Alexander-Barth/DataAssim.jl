using Random
using Test
using DataAssim: talagrand_diagram

# Example from
# ABC of ensemble forecasting
# by Peter Houtekamer, ARMA
# https://collaboration.cmc.ec.gc.ca/cmc/cmoi/product_guide/docs/lib/ens_en.pdf
# https://web.archive.org/web/20231118202650/https://collaboration.cmc.ec.gc.ca/cmc/cmoi/product_guide/docs/lib/ens_en.pdf

x = [1.5;; 2.3;; 0.8;; 4.1;; 0.3]
y = [2.95]


freq = talagrand_diagram(x,y)
@test freq == [0,0,0,0,1,0]
