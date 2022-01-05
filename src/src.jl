"""
Computes mean resultant vector length for circular data.

	1. α: sample of angles in radians
    - w: number of incidences in case of binned angle data
    - d: spacing of bin centers for binned data, used to correct for bias in estimation of r, in radians
    - dims: compute along this dimension(default=1)

	return: mean resultant length
"""
function circ_r(α; w = ones(size(α)), d = 0, dims = 1)
  # compute weighted sum of cos and sin of angles
  r = sum(w .* exp.(im * α); dims)

  # obtain length
  r = abs.(r) ./ sum(w; dims)

  # for data with known spacing, apply correction factor to correct for bias in the estimation of r (see Zar, p. 601, equ. 26.16)
  if d != 0
    c = d / 2 / sin(d / 2)
    r *= c
  end
  length(r) == 1 ? r[1] : r
end

"""
Computes the mean direction for circular data.

	1. α: sample of angles in radians
    - w: weightings in case of binned angle data
    - dims: compute along this dimension(default=1)

	return:
	- μ: mean direction
	- ul: upper 95% confidence limit
	- ll: lower 95% confidence limit
"""
function circ_mean(α; w = ones(size(α)), dims = 1)
  # compute weighted sum of cos and sin of angles
  r = sum(w .* exp.(im * α); dims)

  # obtain mean by
  μ = angle.(r)
  μ = length(μ) == 1 ? μ[1] : μ

  # confidence limits
  t = circ_confmean(α; xi = 0.05, w, d = 0, dims)
  ul = μ .+ t
  ll = μ .- t

  return (; μ, ul, ll)
end

"""
Transforms p-axial data to a common scale.

  1. α: sample of angles in radians
  2. p: number of modes(default=1)

  return: transformed data
"""
circ_axial(α, p = 1) = mod.(α * p, 2π)

"""
Computes the median direction for circular data.

  1. α: sample of angles in radians

  return: median direction
"""
function circ_median(α::AbstractVector)
  α = mod.(α, 2π)
  n = length(α)

  dd = circ_dist2(α)
  m1 = dropdims(count(>=(0), dd; dims = 1), dims = 1)
  m2 = dropdims(count(<=(0), dd; dims = 1), dims = 1)

  dm = abs.(m1 .- m2)
  m = minimum(dm)
  md = circ_mean(α[dm.==m]).μ
  αμ = circ_mean(α).μ
  if abs(circ_dist(αμ, md)) > abs(circ_dist(αμ, md + π))
    md = mod(md + π, 2π)
  end
  return md
end

"""
Computes circular standard deviation for circular data (equ. 26.20, Zar).

	1. α: sample of angles in radians
   - w: weightings in case of binned angle data
   - d: spacing of bin centers for binned data, used to correct for bias in estimation of r, in radians
   - dims: compute along this dimension(default=1)

	return:
  - s: angular deviation
  - s0: circular standard deviation
"""
function circ_std(α; w = ones(size(α)), d = 0, dims = 1)
  # compute mean resultant vector length
  r = circ_r(α; w, d, dims)

  s = sqrt.(2 .* (1 .- r))      # 26.20
  s0 = sqrt.(-2 .* log.(r))   # 26.21
  return (; s, s0)
end

"""
Computes circular variance for circular data (equ. 26.17/18, Zar).

	1. α: sample of angles in radians
  - w: number of incidences in case of binned angle data
  - d: spacing of bin centers for binned data, used to correct for bias in estimation of r, in radians
	- dims: compute along this dimension(default=1)

	return:
	- S: circular variance 1-r
	- s: angular variance 2(1-r)
"""
function circ_var(α; w = ones(size(α)), d = 0, dims = 1)
  # compute mean resultant vector length
  r = circ_r(α; w, d, dims)

  # apply transformation to var
  S = 1 .- r
  s = 2 * S
  return (; S, s)
end


"""
Computes the confidence limits on the mean for circular data.

	1. α: sample of angles in radians
   - xi: (1-xi) confidence limits are computed, default=0.05
   - w: number of incidences in case of binned angle data
   - d: spacing of bin centers for binned data, used to correct for bias in estimation of r, in radians
   - dims: compute along this dimension(default=1)

	return: mean ± d yields upper/lower (1-xi)% confidence limit
"""
function circ_confmean(α; xi = 0.05, w = ones(size(α)), d = 0, dims = 1)
  # compute ingredients for conf. lim.
  r = circ_r(α; w, d, dims)
  n = sum(w; dims)
  R = n .* r
  c2 = quantile(Chisq(1), 1 - xi)

  # check for resultant vector length and select appropriate formula
  t = zeros(size(r))

  for i = 1:length(r)
    if r[i] < 0.9 && r[i] > sqrt(c2 / 2 / n[i])
      t[i] = sqrt((2 * n[i] * (2 * R[i]^2 - n[i] * c2)) / (4 * n[i] - c2))  # equ. 26.24
    elseif r[i] >= 0.9
      t[i] = sqrt(n[i]^2 - (n[i]^2 - R[i]^2) * exp(c2 / n[i]))      # equ. 26.25
    else
      t[i] = NaN
      @warn "Requirements for confidence levels not met."
    end
  end

  # apply final transform
  t = acos.(t ./ R)
  t = length(t) == 1 ? t[1] : t
end


"""
Calculates a measure of angular skewness.

	1. α: sample of angles in radians
    - w: weightings in case of binned angle data
    - dims: compute along this dimension(default=1)

	return:
	- b: skewness (from Pewsey)
	- b0: alternative skewness measure (from Fisher)
"""
function circ_skewness(α; w = ones(size(α)), dims = 1)
  # compute neccessary values
  R = circ_r(α; w, d = 0, dims)
  θ = circ_mean(α; w, dims).μ
  _, ρ₂, μ₂ = circ_moment(α; w, p = 2, cent = false, dims)

  # compute skewness
  if ndims(θ) > 0
    θ₂ = repeat(θ, outer = Int.(size(α) ./ size(θ)))
  else
    θ₂ = θ
  end
  b = sum(w .* sin.(2 * circ_dist(α, θ₂)); dims) ./ sum(w; dims)
  b0 = ρ₂ .* sin.(circ_dist(μ₂, 2θ)) ./ (1 .- R) .^ (3 / 2)   # (formula 2.29)
  b = length(b) == 1 ? b[1] : b
  b0 = length(b0) == 1 ? b0[1] : b0
  return (; b, b0)
end


"""
Calculates a measure of angular kurtosis.

  1. α: sample of angles in radians
    - w: weightings in case of binned angle data
    - dims: compute along this dimension(default=1)

  return:
  - k: kurtosis (from Pewsey)
  - k0: kurtosis (from Fisher)
"""
function circ_kurtosis(α; w = ones(size(α)), dims = 1)
  # compute mean direction
  R = circ_r(α; w, d = 0, dims)
  θ = circ_mean(α; w, dims).μ
  _, ρ₂, _ = circ_moment(α; w, p = 2, cent = true, dims)
  _, _, μ₂ = circ_moment(α; w, p = 2, cent = false, dims)

  # compute skewness
  if ndims(θ) > 0
    θ₂ = repeat(θ, outer = Int.(size(α) ./ size(θ)))
  else
    θ₂ = θ
  end
  k = sum(w .* cos.(2 * circ_dist(α, θ₂)); dims) ./ sum(w; dims)
  k0 = (ρ₂ .* cos.(circ_dist(μ₂, 2θ)) .- R .^ 4) ./ (1 .- R) .^ 2    # (formula 2.30)
  k = length(k) == 1 ? k[1] : k
  k0 = length(k0) == 1 ? k0[1] : k0
  return (; k, k0)
end


"""
Calculates the complex p-th centred or non-centred moment of the angular data in angle.

  1. α: sample of angles in radians
   - w: number of incidences in case of binned angle data
   - p: p-th moment to be computed, default=1
   - cent: if true, central moments are computed, default = false
   - dims: compute along this dimension(default=1)

  return:
  - mp: complex p-th moment
  - ρp: magnitude of the p-th moment
  - μp: angle of th p-th moment
"""
function circ_moment(α; w = ones(size(α)), p = 1, cent = false, dims = 1)
  if cent
    θ = circ_mean(α; w, dims).μ
    if ndims(θ) > 0
      v = Int.(size(α) ./ size(θ))
      θ = repeat(θ, outer = v)
    end
    α = circ_dist(α, θ)
  end

  n = size(α, dims)
  c̄ = sum(cos.(p * α) .* w; dims) / n
  s̄ = sum(sin.(p * α) .* w; dims) / n
  mp = c̄ .+ im * s̄
  ρp = abs.(mp)
  μp = angle.(mp)
  mp = length(mp) == 1 ? mp[1] : mp
  ρp = length(ρp) == 1 ? ρp[1] : ρp
  μp = length(μp) == 1 ? μp[1] : μp
  return (; mp, ρp, μp)
end



"""
Pairwise difference α-β around the circle computed efficiently.

  1. α: sample of linear random variable
  2. β: sample of linear random variable or one single angle

  return: matrix with differences
"""
circ_dist(α, β) = angle.(exp.(im * α) ./ exp.(im * β))


"""
All pairwise difference α-β around the circle computed efficiently.

  1. α: sample of linear random variable
  2. β: sample of linear random variable

  return: matrix with pairwise differences
"""
circ_dist2(α) = circ_dist2(α, α)
circ_dist2(α, β) = angle.(exp.(im * α) ./ exp.(im * β'))


"""
Computes descriptive statistics for circular data.

  1. α: sample of angles in radians
    - w: weightings in case of binned angle data
    - d: spacing of bin centers for binned data, used to correct for bias in estimation of r, in radians

  return: descriptive statistics
"""
function circ_stats(α::AbstractVector; w = ones(size(α)), d = 0)
  # mean
  mean = circ_mean(α; w).μ
  # median
  median = circ_median(α)
  # variance
  var = circ_var(α; w, d)
  # standard deviation
  std, std0 = circ_std(α; w, d)
  # skewness
  skewness, skewness0 = circ_skewness(α; w)
  # kurtosis
  kurtosis, kurtosis0 = circ_kurtosis(α; w)

  return (;
    mean,
    median,
    var,
    std,
    std0,
    skewness,
    skewness0,
    kurtosis,
    kurtosis0,
  )
end



"""
Computes Rayleigh test for non-uniformity of circular data.

	H0: the population is uniformly distributed around the circle
  HA: the populatoin is not distributed uniformly around the circle
	Assumption: the distribution has maximally one mode and the data is
							sampled from a von Mises distribution!

	1. α: sample of angles in radians
    - w: number of incidences in case of binned angle data
    - d: spacing of bin centers for binned data, used to correct for bias in estimation of r, in radians

	return:
	- p: p-value of Rayleigh's test
	- z: value of the z-statistic
"""
function circ_rtest(α; w = ones(size(α)), d = 0)
  r = circ_r(α; w, d)
  n = sum(w)

  # compute Rayleigh's R (equ. 27.1)
  R = n * r

  # compute Rayleigh's z (equ. 27.2)
  z = R^2 / n

  # compute p value using approxation in Zar, p. 617
  p = exp(sqrt(1 + 4n + 4(n^2 - R^2)) - (1 + 2n))
  return (; p, z)
end


"""
Computes Omnibus or Hodges-Ajne test for non-uniformity of circular data.

	H0: the population is uniformly distributed around the circle
	HA: the population is not distributed uniformly around the circle
	Alternative to the Rayleigh and Rao's test. Works well for unimodal,
	bimodal or multimodal data. If requirements of the Rayleigh test are met, the latter is more powerful.

	1. α: sample of angles in radians
    - sz: step size for evaluating distribution, default 1 degree
    - w: number of incidences in case of binned angle data

	return:
	- p: p-value
	- m: minimum number of samples falling in one half of the circle
"""
function circ_otest(α; sz = deg2rad(1), w = ones(size(α)))
  α = mod.(α, 2π)
  n = sum(w)
  dg = 0:sz:π

  m1 = zeros(size(dg))
  m2 = zeros(size(dg))
  for i = 1:length(dg)
    m1[i] = sum((α .> dg[i]) .& (α .< π + dg[i]) .* w)
    m2[i] = n - m1[i]
  end
  m = min(minimum(m1), minimum(m2))

  if n > 50
    # approximation by Ajne (1968)
    A = π * sqrt(n) / 2 / (n - 2m)
    p = sqrt(2π) / A * exp(-π^2 / 8 / A^2)
  else
    # exact formula by Hodges (1955)
    # p = 2^(1-n) * (n-2m) * nchoosek(n,m)  # revised below for numerical stability
    p = exp(
      (1 - n) * log(2) + log(n - 2m) + loggamma(n + 1) - loggamma(m + 1) -
      loggamma(n - m + 1),
    )
  end
  return (; p, m)
end


# circ_raotest Rao""s spacing test for nonuniformity
# circ_vtest V-Test for nonuniformity with known mean direction
# circ_medtest Test for median angle
# circ_mtest One-sample test for specified mean direction


"""
Parametric Watson-Williams multi-sample test for equal means. Can be used as a one-way ANOVA test for circular data.

  H0: the s populations have equal means
  HA: the s populations have unequal means

  !!! note
  Use with binned data is only advisable if binning is finer than 10 deg.
  In this case, α is assumed to correspond to bin centers.

  The Watson-Williams two-sample test assumes underlying von-Mises distributrions.
  All groups are assumed to have a common concentration parameter k.

  1. α: angles in radians
  2. idx: indicates which population the respective angle in α comes from, 1:s
    - w: number of incidences in case of binned angle data

  return:
  - p: p-value of the Watson-Williams multi-sample test. Discard H0 if p is small.
  - F: F statistics
"""
function circ_wwtest(α, idx; w = ones(size(α)))
  # number of groups
  u = unique(idx)
  s = length(u)

  # number of samples
  n = sum(w)

  # compute relevant quantitites
  pn = zeros(s)
  pr = zeros(s)
  for t = 1:s
    pidx = idx .== u[t]
    pn[t] = sum(pidx .* w)
    pr[t] = circ_r(α[pidx]; w = w[pidx])
  end

  r = circ_r(α; w)
  rw = sum(pn .* pr) / n

  # make sure assumptions are satisfied
  checkAssumption =
    (rw, n) -> begin
      if n >= 11 && rw < 0.45
        @warn "Test not applicable. Average resultant vector length < 0.45."
      elseif n < 11 && n >= 7 && rw < 0.5
        @warn "Test not applicable. Average number of samples per population 6 < x < 11 and average resultant vector length < 0.5."
      elseif n >= 5 && n < 7 && rw < 0.55
        @warn "Test not applicable. Average number of samples per population 4 < x < 7 and average resultant vector length < 0.55."
      elseif n < 5
        @warn "Test not applicable. Average number of samples per population < 5."
      end
    end
  checkAssumption(rw, mean(pn))

  # test statistic
  kk = circ_kappa(rw)
  β = 1 + 3 / (8 * kk)    # correction factor
  A = sum(pr .* pn) - r * n
  B = n - sum(pr .* pn)

  F = β * (n - s) * A / (s - 1) / B
  # p = 1 - fcdf(F, s - 1, n - s)
  p = 1 - cdf(FDist(s - 1, n - s), F)

  # fprintf('\nANALYSIS OF VARIANCE TABLE (WATSON-WILLIAMS TEST)\n\n');
  # fprintf('%s\t\t\t\t%s\t%s\t\t%s\t\t%s\t\t\t%s\n', ' ' ,'d.f.', 'SS', 'MS', 'F', 'P-Value');
  # fprintf('--------------------------------------------------------------------\n');
  # fprintf('%s\t\t\t%u\t\t%.2f\t%.2f\t%.2f\t\t%.4f\n', 'Columns', s-1 , A, A/(s-1), F, pval);
  # fprintf('%s\t\t%u\t\t%.2f\t%.2f\n', 'Residual ', n-s, B, B/(n-s));
  # fprintf('--------------------------------------------------------------------\n');
  # fprintf('%s\t\t%u\t\t%.2f', 'Total   ',n-1,A+B);
  # fprintf('\n\n')

  # table = ["Source",'d.f.','SS','MS','F','P-Value';
  #          'Columns', s-1 , A, A/(s-1), F, pval;
  #          'Residual ', n-s, B, B/(n-s), [], [];
  #          'Total',n-1,A+B,[],[],[]]

  return (; p, F)
end


"""
Parametric two-way ANOVA for circular data with interations.

  !!!note
  The test assumes underlying von-Mises distributrions.
  All groups are assumed to have a common concentration parameter k, between 0 and 2.

  1. α: angles in radians
  2. idp: indicates the level of factor 1 (1:p)
  3. idq: indicates the level of factor 2 (1:q)
    - inter: whether to include effect of interaction
    - fn: string array containing names of the factors

  return:
  - p: pvalues for factors and interaction
  - s: statistic of each p value
"""
function circ_hktest(α, idp, idq; inter = true, fn = ['A', 'B'])
  # number of groups for every factor
  pu = unique(idp)
  p = length(pu)
  qu = unique(idq)
  q = length(qu)

  # number of samples
  n = length(α)

  # compute important sums for the test statistics
  cn = zeros(p, q)
  cr = zeros(p, q)
  pm = zeros(p)
  pr = zeros(p)
  pn = zeros(p)
  qm = zeros(q)
  qr = zeros(q)
  qn = zeros(q)
  for pp = 1:p
    p_id = idp .== pu[pp] # indices of factor1 = pp
    for qq = 1:q
      q_id = idq .== qu[qq] # indices of factor2 = qq
      idx = p_id .& q_id
      cn[pp, qq] = sum(idx)     # number of items in cell
      cr[pp, qq] = cn[pp, qq] * circ_r(α[idx]) # R of cell
    end
    # R and mean angle for factor 1
    pr[pp] = sum(p_id) * circ_r(α[p_id])
    pm[pp] = circ_mean(α[p_id]).μ
    pn[pp] = sum(p_id)
  end

  # R and mean angle for factor 2
  for qq = 1:q
    q_id = idq .== qu[qq]
    qr[qq] = sum(q_id) * circ_r(α[q_id])
    qm[qq] = circ_mean(α[q_id]).μ
    qn[qq] = sum(q_id)
  end

  # R and mean angle for whole sample (total)
  tr = n * circ_r(α)

  # estimate kappa
  kk = circ_kappa(tr / n)

  # different formulas for different width of the distribution
  if kk > 2
    # large kappa

    # effect of factor 1
    eff_1 = sum(pr .^ 2 ./ dropdims(sum(cn, dims = 2), dims = 2)) - tr .^ 2 / n
    df_1 = p - 1
    ms_1 = eff_1 / df_1

    # effect of factor 2
    eff_2 = sum(qr .^ 2 ./ dropdims(sum(cn, dims = 1), dims = 1)) - tr .^ 2 / n
    df_2 = q - 1
    ms_2 = eff_2 / df_2

    # total effect
    eff_t = n - tr .^ 2 / n
    df_t = n - 1

    m = mean(cn)

    if inter
      # correction factor for improved F statistic
      β = 1 / (1 - 1 / (5kk) - 1 / (10kk^2))

      # residual effects
      eff_r = n - sum(cr .^ 2 ./ cn)
      df_r = p * q * (m - 1)
      ms_r = eff_r / df_r

      # interaction effects
      eff_i =
        sum(cr .^ 2 ./ cn) - sum(qr .^ 2 ./ qn) - sum(pr .^ 2 ./ pn) +
        tr .^ 2 / n
      df_i = (p - 1) * (q - 1)
      ms_i = eff_i / df_i

      # interaction test statistic
      FI = ms_i / ms_r
      pI = 1 - cdf(FDist(df_i, df_r), FI)

    else
      # residual effect
      eff_r = n - sum(qr .^ 2 ./ qn) - sum(pr .^ 2 ./ pn) + tr .^ 2 / n
      df_r = (p - 1) * (q - 1)
      ms_r = eff_r / df_r

      # interaction effects
      eff_i = []
      df_i = []
      ms_i = []

      # interaction test statistic
      FI = []
      pI = NaN
      β = 1
    end

    # compute all test statistics as F = β * MS(A) / MS(R)
    F1 = β * ms_1 / ms_r
    p1 = 1 - cdf(FDist(df_1, df_r), F1)

    F2 = β * ms_2 / ms_r
    p2 = 1 - cdf(FDist(df_2, df_r), F2)
    s = [F1, F2, FI]
  else
    # small kappa

    # correction factor
    rr = besseli(1, kk) / besseli(0, kk)
    f = 2 / (1 - rr^2)

    chi1 = f * (sum(pr .^ 2 ./ pn) - tr .^ 2 / n)
    df_1 = 2 * (p - 1)
    p1 = 1 - cdf(Chisq(df_1), chi1)

    chi2 = f * (sum(qr .^ 2 ./ qn) - tr .^ 2 / n)
    df_2 = 2 * (q - 1)
    p2 = 1 - cdf(Chisq(df_2), chi2)

    chiI =
      f * (
        sum(cr .^ 2 ./ cn) - sum(pr .^ 2 ./ pn) - sum(qr .^ 2 ./ qn) +
        tr .^ 2 / n
      )
    df_i = (p - 1) * (q - 1)
    pI = 1 - cdf(Chisq(df_i), chiI)
    s = [chi1, chi2, chiI]
  end
  p = [p1, p2, pI]
  #
  # if kk>2
  #
  #   fprintf('\nANALYSIS OF VARIANCE TABLE (HIGH KAPPA MODE)\n\n');
  #
  #   fprintf('%s\t\t\t\t%s\t%s\t\t%s\t\t%s\t\t\t%s\n', ' ' ,'d.f.', 'SS', 'MS', 'F', 'P-Value');
  #   fprintf('--------------------------------------------------------------------\n');
  #   fprintf('%s\t\t\t\t%u\t\t%.2f\t%.2f\t%.2f\t\t%.4f\n', fn{1}, df_1 , eff_1, ms_1, F1, p1);
  #   fprintf('%s\t\t\t\t%u\t\t%.2f\t%.2f\t%.2f\t\t%.4f\n', fn{2}, df_2 , eff_2, ms_2, F2, p2);
  #   if (inter)
  #       fprintf('%s\t\t%u\t\t%.2f\t%.2f\t%.2f\t\t%.4f\n', 'Interaction', df_i , eff_i, ms_i, FI, pI);
  #   end
  #   fprintf('%s\t\t%u\t\t%.2f\t%.2f\n', 'Residual ', df_r, eff_r, ms_r);
  #   fprintf('--------------------------------------------------------------------\n');
  #   fprintf('%s\t\t%u\t\t%.2f', 'Total   ',df_t,eff_t);
  #   fprintf('\n\n')
  # else
  #   fprintf('\nANALYSIS OF VARIANCE TABLE (LOW KAPPA MODE)\n\n');
  #
  #   fprintf('%s\t\t\t\t%s\t%s\t\t\t%s\n', ' ' ,'d.f.', 'CHI2', 'P-Value');
  #   fprintf('--------------------------------------------------------------------\n');
  #   fprintf('%s\t\t\t\t%u\t\t%.2f\t\t\t%.4f\n', fn{1}, df_1 , chi1, p1);
  #   fprintf('%s\t\t\t\t%u\t\t%.2f\t\t\t%.4f\n', fn{2}, df_2 , chi2, p2);
  #   if (inter)
  #       fprintf('%s\t\t%u\t\t%.2f\t\t\t%.4f\n', 'Interaction', df_i , chiI, pI);
  #   end
  #   fprintf('--------------------------------------------------------------------\n');
  #   fprintf('\n\n')
  #
  # end
  return (; p, s)
end



# circ_ktest Test for equal concentration parameter
# circ_symtest Test for symmetry around median angle
# circ_kuipertest Test whether two distributions are identical (like KS test)


"""
Circular correlation coefficient for two circular random variables.

  1. α₁: sample of angles in radians
  2. α₂: sample of angles in radians

  return:
  - ρ: correlation coefficient
  - p: p-value
"""
function circ_corrcc(α₁, α₂)
  # compute mean directions
  n = length(α₁)
  ᾱ₁ = circ_mean(α₁).μ
  ᾱ₂ = circ_mean(α₂).μ

  # compute correlation coeffcient from p. 176
  num = sum(sin.(α₁ .- ᾱ₁) .* sin.(α₂ .- ᾱ₂))
  den = sqrt(sum(sin.(α₁ .- ᾱ₁) .^ 2) .* sum(sin.(α₂ .- ᾱ₂) .^ 2))
  ρ = num / den

  # compute pvalue
  l20 = mean(sin.(α₁ .- ᾱ₁) .^ 2)
  l02 = mean(sin.(α₂ .- ᾱ₂) .^ 2)
  l22 = mean((sin.(α₁ .- ᾱ₁) .^ 2) .* (sin.(α₂ .- ᾱ₂) .^ 2))

  ts = sqrt((n * l20 * l02) / l22) * ρ
  p = 2 * (1 - cdf(Normal(), abs(ts)))
  return (; ρ, p)
end


"""
Correlation coefficient between one circular and one linear random variable.

  1. α: sample of angles in radians
  2. x: sample of linear random variable

  return:
  - ρ: correlation coefficient
  - p: p-value
"""
function circ_corrcl(α, x)
  n = length(α)

  # compute correlation coefficent for sin and cos independently
  rxs = cor(x, sin.(α))
  rxc = cor(x, cos.(α))
  rcs = cor(sin.(α), cos.(α))

  # compute angular-linear correlation (equ. 27.47)
  ρ = sqrt((rxc^2 + rxs^2 - 2 * rxc * rxs * rcs) / (1 - rcs^2))

  # compute pvalue
  p = 1 - cdf(Chisq(2), n * ρ^2)
  return (; ρ, p)
end


"""
Computes an approximation to the ML estimate of the concentration parameter kappa of the von Mises distribution.

  1. α: angles in radians OR α is length resultant
    - w: number of incidences in case of binned angle data

  return: estimated value of kappa
"""
function circ_kappa(α; w = ones(size(α)))
  N = length(α)

  if N > 1
    R = circ_r(α; w)
  else
    R = α
  end

  if R < 0.53
    κ = 2R + R^3 + 5R^5 / 6
  elseif R >= 0.53 && R < 0.85
    κ = -0.4 + 1.39R + 0.43 / (1 - R)
  else
    κ = 1 / (R^3 - 4R^2 + 3R)
  end

  if N < 15 && N > 1
    if κ < 2
      κ = max(κ - 2 * (N * κ)^-1, 0)
    else
      κ = (N - 1)^3 * κ / (N^3 + N)
    end
  end
  κ
end




# circ_plot Visualization for circular data
# circ_clust Simple clustering for circular data
# circ_samplecdf Evaluate CDF of a sample of angles
