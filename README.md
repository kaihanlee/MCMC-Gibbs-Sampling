# MCMC Gibbs Sampling

main code is in BI_code.R

Daily sales of a store are known to be affected by whether the store has announced super discounts on that day or not. From past experience, it is known that the fluctuation of the amount of daily sales around its mean value is normally distributed; if Y denotes the amount of daily sales, then Y ~ N(μ, Ψ) where μ and Ψ denote the unknown mean and precision parameter respectively. Super discount days may change the mean level, but leave the precision unaffected. To investigate the extent to which super discount days influence μ, a random sample of n = 100 daily sales amounts is obtained. The sample is represented as y = (y1, y2, ..., yn) where yi ~ N(μ_i, Ψ), independently for i = 1, 2, ..., n. An indicator variable di, for i = 1, 2, ..., n, denotes whether the i-th observation was taken on a day when super discount was announced (corresponding to di = 1). The following equation gives a model for μ_i that captures the effect of di: 

μ_i = alf + γ * d_i

for i = 1, 2, ..., n. The store provides you with data on y but is unable to provide the data on d = (d1, d2, ..., dn). You gather information from them that super discounts days are held on 100p% (0 < p < 1) of days in a year, independently of each other. The departmetal store believes that p is around 0.4.
