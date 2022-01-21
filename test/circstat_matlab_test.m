% test data
av = [0:0.2:pi,0.4*pi:0.05:0.6*pi,pi:0.4:2*pi]'
am = [av av]
aa = cat(3,am,am)
wv = [0:0.1:0.5*pi,0.4*pi:0.05:0.6*pi,0:0.2:0.5*pi]'
wm = [wv wv]
wa = cat(3,wm,wm)

circ_r(av)
circ_r(av,wv)
circ_r(av,wv,0.1)

[mu, ul, ll] = circ_mean(av)
[mu, ul, ll] = circ_mean(av,wv)

[mu, ul, ll] = circ_mean(am)
[mu, ul, ll] = circ_mean(am,wm)
[mu, ul, ll] = circ_mean(am,wm,2)

circ_median(av)

[s, s0] = circ_std(av)

[S, s] = circ_var(av)


circ_dist(av,pi)
circ_dist(av,flipud(av))

a = [1,2,3]'
b = [2,3]'
circ_dist2(a)
circ_dist2(a,b)

[mp, rp, mup] = circ_moment(av)
[mp, rp, mup] = circ_moment(av,[],[],true)

[b, b0] = circ_skewness(av)
[b, b0] = circ_skewness(am)

[k, k0] = circ_kurtosis(av)
[k, k0] = circ_kurtosis(am)

[p,z] = circ_rtest(av)
[p,m] = circ_otest(av)

idx = [ones(18,1);2*ones(19,1)]
[p,F] = circ_wwtest(av,idx)

idp = [repmat([1;2;3;4],9,1);1]
idq = [repmat([1;2],18,1);1]
[p,s] = circ_hktest(av,idp,idq)
[p,f] = circ_ktest(av,flipud(av))
p = circ_symtest(av)
[p,k,K]=circ_kuipertest(av,flipud(av))
[p, u, UC] = circ_raotest(av)
[p, v] = circ_vtest(av,0)
p = circ_medtest(av,0)
[h, mu, ul, ll] = circ_mtest(av,0)


[r,p] = circ_corrcc(av,flipud(av))
[r,p] = circ_corrcl(av,flipud(av))

circ_kappa(av)

[phi,c] = circ_samplecdf(av,8)


