using Test, CircStats

@testset "CircStats" begin
    av = [0:0.2:π;0.4π:0.05:0.6π;π:0.4:2π]
    am = [av av]
    aa = [am;;;am]
    wv = [0:0.1:0.5π;0.4π:0.05:0.6π;0:0.2:0.5π]
    wm = [wv wv]
    wa = [wm;;;wm]
    d = 0.1

    @test circ_r(av) ≈ 0.481125184946105
    @test circ_r(av;w=wv)  ≈ 0.634569454038054
    @test circ_r(av;w=wv,d)  ≈ 0.634833935115386

    circ_r(am)
    circ_r(am;w=wm)
    circ_r(am;w=wm,d)
    circ_r(am;w=wm,d,dims=2)

    circ_r(aa)
    circ_r(aa;w=wa)
    circ_r(aa;w=wa,d)
    circ_r(aa;w=wa,d,dims=3)

    μ, ul, ll = circ_mean(av)
    @test μ ≈ 1.568889360630885
    @test ul ≈ 2.037005335388447
    @test ll ≈ 1.100773385873324
    μ, ul, ll = circ_mean(av;w=wv)
    @test μ ≈ 1.691236041279729
    @test ul ≈ 2.018457532767024
    @test ll ≈ 1.364014549792434

    circ_mean(am)
    circ_mean(am;w=wm)

    circ_mean(aa)
    circ_mean(aa;w=wa)

    @test circ_median(av) ≈ 1.556637061435917

    s, s0 = circ_std(av)
    @test s ≈ 1.018699970603607
    @test s0 ≈ 1.209651009981690

    circ_std(am)
    circ_std(am;dims=2)
    circ_std(aa)
    circ_std(am;dims=3)

    S, s = circ_var(av)
    @test S ≈ 0.518874815053895
    @test s ≈ 1.037749630107790

    circ_var(am)
    circ_var(am;dims=2)
    circ_var(aa)
    circ_var(am;dims=3)

    α = [1,2,3];β=[2,3,4]
    circ_dist(α,π)
    circ_dist(α,β)

    circ_dist2(α)
    circ_dist2(α,β)

    mp, ρp, μp = circ_moment(av)
    @test mp ≈ 0.000917488892268 + 0.481124310135703im
    @test ρp ≈ 0.481125184946105
    @test μp ≈ 1.568889360630885
    mp, ρp, μp = circ_moment(av;cent=true)
    @test mp ≈ 0.481125184946105 - 0.000000000000000im
    @test ρp ≈ 0.481125184946105
    @test μp ≈ 3.741981750042832e-17 atol=1e-16

    mp, ρp, μp = circ_moment(am)
    mp, ρp, μp = circ_moment(am;cent=true)

    b, b0 = circ_skewness(av)
    @test b ≈ -0.005585405160200
    @test b0 ≈ -0.014943791328267
    b, b0 = circ_skewness(am)

    k, k0 = circ_kurtosis(av)
    @test k ≈ 0.315477455279792
    @test k0 ≈ 0.972747287142949
    k, k0 = circ_kurtosis(am)

    circ_stats(av)

    p,z = circ_rtest(av)
    @test p ≈ 1.247327496731614e-04
    @test z ≈ 8.564813412808679

    p,m = circ_otest(av)
    @test p ≈ 0.003445833222941
    @test m ≈ 7

    idx = [fill(1,18);fill(2,19)]
    p,F = circ_wwtest(av,idx)
    @test p ≈ 0.489205605248093
    @test F ≈ 0.488523443700675

    idp = [repeat([1,2,3,4],outer=9);1]
    idq = [repeat([1,2],outer=18);1]
    p,s = circ_hktest(av,idp,idq)
    @test p[1] ≈ 0.988563843496855
    @test p[2] ≈ 0.865090769456459
    @test isnan(p[3])
    @test s[1] ≈ 0.917000658130361
    @test s[2] ≈ 0.289841683535636
    @test isnan(s[3])

    p,f = circ_ktest(av,reverse(av))
    @test p ≈ 0.999999999999993
    @test f ≈ 1

    p = circ_symtest(av)
    # @test p ≈ 0.850454511827570 # diff in signrank in matlab?

    p,k,K=circ_kuipertest(av,reverse(av))
    @test p ≈ 1
    @test k ≈ 0

    p, u, UC = circ_raotest(av)
    @test p ≈ 0.5
    @test u ≈ 121.1858500639162
    @test UC ≈ 153.82

    p, v = circ_vtest(av)
    @test p ≈ 0.496851365629160
    @test v ≈ 0.033947089013911

    p = circ_medtest(av)
    @test p ≈ .0004719867429230360

    h, μ, ul, ll = circ_mtest(av)
    @test h
    @test μ ≈ 1.568889360630885
    @test ul ≈ 2.037005335388447
    @test ll ≈ 1.100773385873324

    ρ,p = circ_corrcc(av,reverse(av))
    @test ρ ≈ 0.204639214081096
    @test p ≈ 0.187332538620452

    ρ,p = circ_corrcl(av,reverse(av))
    @test ρ ≈ 0.591325482630456
    @test p ≈ 0.001551058328371

    @test circ_kappa(av) ≈ 1.095105628679724

    ϕ,c = circ_samplecdf(av,n=8)
    @test ϕ[4] ≈ 2.356194490192345
    @test c[4] ≈ 0.783783783783784

end
