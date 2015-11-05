"The binary Kullback-Leibler divergence."
function KL(q, p)
    (0.0 ≤ q ≤ 1.0 && 0.0 ≤ p ≤ 1.0) || error("q and p must in [0,1]")
    kl = 0.0
    if (q > 0.0 && p > 0.0)
        kl += q*log(q/p)
    end
    if q < 1.0 && p < 1.0
        kl += (1.0-q)*log((1.0-q)/(1.0-p))
    end
    kl
end
