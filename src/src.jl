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

  return (; mean, median, var,std,std0, skewness,skewness0,kurtosis,kurtosis0)
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

"""
Calculates Rao's spacing test by comparing distances between points on a circle to those expected from a uniform distribution.
  H0: Data is distributed uniformly around the circle.
  H1: Data is not uniformly distributed around the circle.

  Alternative to the Rayleigh test and the Omnibus test. Less powerful
  than the Rayleigh test when the distribution is unimodal on a global scale but uniform locally.

  Due to the complexity of the distributioin of the test statistic, we resort to the tables published by
  Russell, Gerald S. and Levitin, Daniel J.(1995) An expanded table of probability values for rao's spacing test, Communications in Statistics - Simulation and Computation
  Therefore the reported p-value is the smallest α level at which the test would still be significant.
  If the test is not significant at the α=0.1 level, we return the critical value for α = 0.05 and p = 0.5.

  1. α: sample of angles in radians

  return:
  - p: smallest p-value at which test would be significant
  - u: computed value of the test-statistic u
  - UC: critical value of the test statistic at sig-level
"""
function circ_raotest(α)
  # for the purpose of the test, convert to angles
  α = sort(rad2deg.(α))
  n = length(α)

  # compute test statistic
  u = 0
  λ = 360 / n
  for j = 1:n-1
    ti = α[j+1] - α[j]
    u += abs(ti - λ)
  end

  tn = 360 - α[n] + α[1]
  u = 0.5(u + abs(tn - λ))

  getVal = (N, u) -> begin
      # Table II from Russel and Levitin, 1995
      α = [0.001, 0.01, 0.05, 0.10]
      table = [
        4 247.32 221.14 186.45 168.02
        5 245.19 211.93 183.44 168.66
        6 236.81 206.79 180.65 166.30
        7 229.46 202.55 177.83 165.05
        8 224.41 198.46 175.68 163.56
        9 219.52 195.27 173.68 162.36
        10 215.44 192.37 171.98 161.23
        11 211.87 189.88 170.45 160.24
        12 208.69 187.66 169.09 159.33
        13 205.87 185.68 167.87 158.50
        14 203.33 183.90 166.76 157.75
        15 201.04 182.28 165.75 157.06
        16 198.96 180.81 164.83 156.43
        17 197.05 179.46 163.98 155.84
        18 195.29 178.22 163.20 155.29
        19 193.67 177.08 162.47 154.78
        20 192.17 176.01 161.79 154.31
        21 190.78 175.02 161.16 153.86
        22 189.47 174.10 160.56 153.44
        23 188.25 173.23 160.01 153.05
        24 187.11 172.41 159.48 152.68
        25 186.03 171.64 158.99 152.32
        26 185.01 170.92 158.52 151.99
        27 184.05 170.23 158.07 151.67
        28 183.14 169.58 157.65 151.37
        29 182.28 168.96 157.25 151.08
        30 181.45 168.38 156.87 150.80
        35 177.88 165.81 155.19 149.59
        40 174.99 163.73 153.82 148.60
        45 172.58 162.00 152.68 147.76
        50 170.54 160.53 151.70 147.05
        75 163.60 155.49 148.34 144.56
        100 159.45 152.46 146.29 143.03
        150 154.51 148.84 143.83 141.18
        200 151.56 146.67 142.35 140.06
        300 148.06 144.09 140.57 138.71
        400 145.96 142.54 139.50 137.89
        500 144.54 141.48 138.77 137.33
        600 143.48 140.70 138.23 136.91
        700 142.66 140.09 137.80 136.59
        800 142.00 139.60 137.46 136.33
        900 141.45 139.19 137.18 136.11
        1000 140.99 138.84 136.94 135.92]

      ridx = findfirst(table[:, 1] .>= N)
      cidx = findfirst(table[ridx, 2:end] .< u)

      if isnothing(cidx)
        UC = table[ridx, end-1]
        p = 0.5
      else
        UC = table[ridx, cidx+1]
        p = α[cidx]
      end
      return p, UC
    end

  # get critical value from table
  p, UC = getVal(n, u)

  return (; p, u, UC)
end

"""
Computes V test for non-uniformity of circular data with a specified mean direction.

  H0: the population is uniformly distributed around the circle
  HA: the population is not distributed uniformly around the circle but has a mean of `m`.

  !!!note
    Not rejecting H0 may mean that the population is uniformly distributed around the circle OR
    that it has a mode but that this mode is not centered at `m`.

  The V test has more power than the Rayleigh test and is preferred if there is reason to believe in a specific mean direction.

  1. α: sample of angles in radians
  - m: suspected mean direction, default=0
  - w: number of incidences in case of binned angle data
  - d: spacing of bin centers for binned data, used to correct for bias in estimation of r, in radians

  return:
  - p: p-value of V test
  - v: value of the V statistic
"""
function circ_vtest(α; m = 0, w = ones(size(α)), d = 0)
  # compute some ingredients
  r = circ_r(α; w, d)
  μ, = circ_mean(α; w)
  n = sum(w)

  # compute Rayleigh's R (equ. 27.1)
  R = n * r

  # compute the V statistic (equ. 27.5)
  v = R * cos(μ - m)

  # compute u (equ. 27.6)
  u = v * sqrt(2 / n)

  # compute p-value from one tailed normal approximation
  p = 1 - cdf(Normal(), u)

  return (; p, v)
end


"""
Tests for significance of the median.

  H0: the population has median angle `md`
  HA: the population has not median angle `md`

  1. α: sample of angles in radians
  - md: median to test, default=0

  return: p-value
"""
function circ_medtest(α; md = 0)
  n = length(α)

  # compute deviations from median
  d = circ_dist(α, md)

  n1 = sum(d .< 0)
  n2 = sum(d .> 0)

  # compute p-value with binomial test
  sum(pdf.(Binomial(n, 0.5), [0:min(n1, n2); max(n1, n2):n]))
end


"""
One-Sample test for the mean angle.
  H0: the population has mean `m`.
  HA: the population has not mean `m`.

  !!!note: This is the equvivalent to a one-sample t-test with specified mean direction.

  1. α: sample of angles in radians
  - m: assumed mean direction, default=0
  - w: number of incidences in case of binned angle data
  - d: spacing of bin centers for binned data, used to correct for bias in estimation of r, in radians
  - xi: alpha level of the test

  return:
  - h: false if H0 can not be rejected, true otherwise
  - μ: mean
  - ul: upper (1-xi) confidence level
  - ll: lower (1-xi) confidence level
"""
function circ_mtest(α; xi = 0.05, m = 0, w = ones(size(α)), d = 0)
  # compute ingredients
  μ, = circ_mean(α; w)
  t = circ_confmean(α; xi, w, d)
  ul = μ + t
  ll = μ - t

  # compute test via confidence limits (example 27.3)
  h = abs(circ_dist2(m, μ)) > t

  return (; h, μ, ul, ll)
end



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


"""
A parametric two-sample test to determine whether two concentration parameters are different.

  H0: The two concentration parameters are equal.
  HA: The two concentration parameters are different.

  Assumptions: both samples are drawn from von Mises type distributions and their joint resultant vector length should be > .7

  1. α₁: sample of angles in radians
  2. α₂: sample of angles in radians

  return:
  - p: p-value that samples have different concentrations
  - f: f-statistic calculated
"""
function circ_ktest(α₁, α₂)
  n1 = length(α₁)
  n2 = length(α₂)

  R1 = n1 * circ_r(α₁)
  R2 = n2 * circ_r(α₂)

  # make sure that r̄ > .7
  r̄ = (R1 + R2) / (n1 + n2)

  if r̄ < 0.7
    @warn "Resultant vector length should be > 0.7"
  end

  # calculate test statistic
  f = ((n2 - 1) * (n1 - R1)) / ((n1 - 1) * (n2 - R2))
  if f > 1
    p = 2 * (1 - cdf(FDist(n1, n2), f))
  else
    f = 1 / f
    p = 2 * (1 - cdf(FDist(n2, n1), f))
  end
  return (; p, f)
end


"""
Tests for symmetry about the median.
  H0: the population is symmetrical around the median
  HA: the population is not symmetrical around the median

  1. α: sample of angles in radians

  return: p-value
"""
function circ_symtest(α)
  # compute median
  md = circ_median(α)

  # compute deviations from median
  d = circ_dist(α, md)

  # compute wilcoxon sign rank test
  pvalue(SignedRankTest(d))
end

const kuipertable= [5	1.45800000000000	1.56500000000000	1.68200000000000	1.76300000000000	1.83800000000000	1.92000000000000	1.97000000000000;
                    6	1.47100000000000	1.58200000000000	1.71100000000000	1.79300000000000	1.86700000000000	1.95700000000000	2.02000000000000;
                    7	1.48300000000000	1.59800000000000	1.72700000000000	1.81400000000000	1.89400000000000	1.98700000000000	2.05100000000000;
                    8	1.49300000000000	1.60800000000000	1.74100000000000	1.83000000000000	1.91100000000000	2.00900000000000	2.07700000000000;
                    9	1.50000000000000	1.61800000000000	1.75200000000000	1.84300000000000	1.92600000000000	2.02700000000000	2.09700000000000;
                    10	1.50700000000000	1.62500000000000	1.76100000000000	1.85400000000000	1.93800000000000	2.04100000000000	2.11300000000000;
                    11	1.51300000000000	1.63100000000000	1.76900000000000	1.86200000000000	1.94800000000000	2.05300000000000	2.12600000000000;
                    12	1.51700000000000	1.63700000000000	1.77600000000000	1.87000000000000	1.95700000000000	2.06200000000000	2.13700000000000;
                    13	1.52200000000000	1.64200000000000	1.78200000000000	1.87600000000000	1.96400000000000	2.07100000000000	2.15400000000000;
                    14	1.52500000000000	1.64600000000000	1.78700000000000	1.88200000000000	1.97000000000000	2.07800000000000	2.15400000000000;
                    15	1.52900000000000	1.65000000000000	1.79100000000000	1.88700000000000	1.97600000000000	2.08500000000000	2.16100000000000;
                    16	1.53200000000000	1.65300000000000	1.79500000000000	1.89200000000000	1.98100000000000	2.09000000000000	2.16800000000000;
                    17	1.53400000000000	1.65700000000000	1.79900000000000	1.89600000000000	1.98600000000000	2.09600000000000	2.17300000000000;
                    18	1.53700000000000	1.65900000000000	1.80200000000000	1.89900000000000	1.99000000000000	2.10000000000000	2.17800000000000;
                    19	1.53900000000000	1.66200000000000	1.80500000000000	1.90300000000000	1.99300000000000	2.10400000000000	2.18300000000000;
                    20	1.54100000000000	1.66400000000000	1.80800000000000	1.90600000000000	1.99700000000000	2.10800000000000	2.18700000000000;
                    21	1.54300000000000	1.66700000000000	1.81000000000000	1.90800000000000	2	                2.11200000000000	2.19100000000000;
                    22	1.54500000000000	1.66900000000000	1.81300000000000	1.91100000000000	2.00300000000000	2.11500000000000	2.19400000000000;
                    23	1.54700000000000	1.67000000000000	1.81500000000000	1.91300000000000	2.00500000000000	2.11800000000000	2.19800000000000;
                    24	1.54900000000000	1.67200000000000	1.81700000000000	1.91600000000000	2.00800000000000	2.12100000000000	2.20100000000000;
                    25	1.55000000000000	1.67400000000000	1.81900000000000	1.91800000000000	2.01000000000000	2.12300000000000	2.20300000000000;
                    30	1.55600000000000	1.68100000000000	1.82600000000000	1.92600000000000	2.01900000000000	2.13400000000000	2.21500000000000;
                    35	1.56100000000000	1.68600000000000	1.83200000000000	1.93300000000000	2.02600000000000	2.14100000000000	2.22300000000000;
                    40	1.56500000000000	1.69000000000000	1.83700000000000	1.93800000000000	2.03200000000000	2.14800000000000	2.23000000000000;
                    45	1.56800000000000	1.69400000000000	1.84100000000000	1.94200000000000	2.03600000000000	2.15200000000000	2.23500000000000;
                    50	1.57100000000000	1.69700000000000	1.84400000000000	1.94600000000000	2.04000000000000	2.15700000000000	2.23900000000000;
                    100	1.58800000000000	1.71400000000000	1.86200000000000	1.96500000000000	2.06000000000000	2.17800000000000	2.26200000000000;
                    200	1.60000000000000	1.72600000000000	1.87600000000000	1.97900000000000	2.07500000000000	2.19400000000000	2.27900000000000;
                    500	1.61000000000000	1.73700000000000	1.88700000000000	1.99000000000000	2.08700000000000	2.20700000000000	2.29200000000000;
                    501	1.62000000000000	1.74700000000000	1.89800000000000	2.00100000000000	2.09800000000000	2.21800000000000	2.30300000000000]


"""
The Kuiper two-sample test tests whether the two samples differ significantly.
The difference can be in any property, such as mean location and dispersion.
It is a circular analogue of the Kolmogorov-Smirnov test.

  H0: The two distributions are identical.
  HA: The two distributions are different.

  1. α₁: sample of angles in radians
  2. α₂: sample of angles in radians
    - n: resolution at which the cdf are evaluated

  return:
  - p: p-value; the smallest of .10, .05, .02, .01, .005, .002, .001,
        for which the test statistic is still higher than the respective critical value.
        this is due to the use of tabulated values. if p>.1, pval is set to 1.
  - k: test statistic
  - K: critical value
"""
function circ_kuipertest(α₁, α₂; n = 100)
  n = length(α₁)
  m = length(α₂)

  # create cdfs of both samples
  ϕ₁, cdf1 = circ_samplecdf(α₁; n)
  _, cdf2 = circ_samplecdf(α₂; n)

  # maximal difference between sample cdfs
  dplus, gdpi = findmax([0; cdf1 .- cdf2])
  dminus, gdmi = findmax([0; cdf2 .- cdf1])

  # calculate k-statistic
  k = n * m * (dplus + dminus)

  # find p-value
  p, K = kuiperlookup(min(n, m), k / sqrt(n * m * (n + m)))
  K *= sqrt(n * m * (n + m))

  return (; p, k, K)
end

function kuiperlookup(n, k)
  α = [0.10, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001]
  nn = kuipertable[:, 1]

  # find correct row of the table
  row = findfirst(n.==nn)
  if isnothing(row)
    # find closest value if no entry is present
    row = length(nn) - sum(n .< nn)
    if row == 0
      error("n too small.")
    else
      @warn "n=$n not found in table, using closest n=$(nn[row]) present."
    end
  end

  # find minimal p-value and test-statistic
  idx = findlast(kuipertable[row, 2:end] .< k)
  if isnothing(idx)
    p = 1
    K = NaN
  else
    p = α[idx]
    K = kuipertable[row, idx+1]
  end

  return (; p, K)
end


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

"""
Performs a simple agglomerative clustering of angular data.

1. α: sample of angles in radians
  - k: number of clusters desired, default=2

return:
- cid: cluster id for each entry of α
- α: sorted angles, matched with cid
- μ: mean direction of angles in each cluster
"""
function circ_clust(α; k = 2)

  n = length(α)
  n < k && error("Not enough data for clusters.")

  # prepare data
  cid = 1:n

  # main clustering loop
  num_unique = length(unique(cid))

  while (num_unique > k)

    # find centroid means...

    # calculate the means for each putative cluster
    μ = fill(NaN, n)
    for j = 1:n
      if sum(cid .== j) > 0
        μ[j], = circ_mean(α(cid .== j))
      end
    end

    # find distance between centroids...
    μdist = abs.(circ_dist2(μ))

    # find closest pair of clusters/datapoints
    mindist = minimum(μdist[tril(ones(size(μdist), size(μdist)), -1)==1])
    row, col = findall(μdist .== mindist)

    # update cluster id's
    cid[cid.==maximum(row)] = minimum(col)

    # update stop criteria
    num_unique = length(unique(cid))

  end

  # renumber cluster ids (so cids [1 3 7 10] => [1 2 3 4])
  cid2 = deepcopy(cid)
  uniquecids = unique(cid)
  for j = 1:length(uniquecids)
    cid[cid2.==uniquecids[j]] = j
  end

  # compute final cluster means
  μ = fill(NaN, num_unique)
  for j = 1:num_unique
    if sum(cid .== j) > 0
      μ[j], = circ_mean(α[cid.==j]')
    end
  end

  return (; cid, α, μ)
end

"""
Helper function for circ_kuipertest. Evaluates CDF of sample.

  1. α: sample of angles in radians
    - n: resolution at which the cdf are evaluated

  return:
  - ϕ: angles at which CDF are evaluated
  - c: CDF values at ϕ
"""
function circ_samplecdf(α; n = 100)
  ϕ = range(0, 2π, length = n + 1)

  # ensure all points in α are on interval [0, 2pi)
  α = sort(mod.(α,2π))

  dp = 1 / length(α) # incremental change in probability
  c = cumsum(map((l, r) -> dp * sum(l .<= α .< r), ϕ[1:end-1], ϕ[2:end]))
  ϕ = ϕ[1:end-1]

  return (; ϕ, c)
end
