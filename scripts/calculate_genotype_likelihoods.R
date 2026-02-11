#!/usr/bin/env R

genotype_likelihood <- function(m,g,e,ref,alt){
  (((m-g)*e+g*(1-e))^alt * ((m-g)*(1-e)+g*e)^ref)/(m^(ref+alt))
}

print(-10*log10(genotype_likelihood(m = 2, g= 0, e = 0.01, ref = 22, alt = 4)))
print(-10*log10(genotype_likelihood(m = 2, g= 1, e = 0.01, ref = 22, alt = 4)))
print(-10*log10(genotype_likelihood(m = 2, g= 2, e = 0.01, ref = 22, alt = 4)))
    
