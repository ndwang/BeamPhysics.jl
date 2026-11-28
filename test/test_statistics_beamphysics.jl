using AtomicAndPhysicalConstants: massof, C_LIGHT
using LinearAlgebra: det
import Statistics
import Statistics: mean as _mean
import StatsBase
import StatsBase: weights
import BeamPhysics:
    norm_emit_calc, twiss_calc, twiss_dispersion,
    A_mat_calc, A_inverse_mat_calc, amplitude_calc,
    particle_amplitude, normalized_particle_coordinate,
    delta, ptp

const _STAT_BEAM = ParticleGroup(joinpath(@__DIR__, "initial_particles.h5"))

@testset "Statistics and Beam Physics" begin

    @testset "Basic statistics" begin
        @testset "mean" begin
            pg = _STAT_BEAM
            @test abs(Statistics.mean(pg, :x)) < 1e-4

            pg1 = single_particle(x=42.0)
            @test Statistics.mean(pg1, :x) == 42.0
        end

        @testset "std" begin
            pg = _STAT_BEAM
            @test Statistics.std(pg, :x) > 0
            @test Statistics.std(pg, :pz) > 0

            pg1 = single_particle(x=1.0)
            @test Statistics.std(pg1, :x) == 0.0
        end

        @testset "minimum / maximum" begin
            pg = _STAT_BEAM
            @test Base.minimum(pg, :x) == minimum(pg.x)
            @test Base.maximum(pg, :x) == maximum(pg.x)
            @test Base.minimum(pg, :x) < Base.maximum(pg, :x)
        end

        @testset "weighted moments with stored and computed columns" begin
            pg = single_particle(x=1.0, px=1e6, pz=2e6, weight=1.0) +
                 single_particle(x=5.0, px=3e6, pz=2e6, weight=3.0, status=0)
            beta = beta_x(pg)
            expected = 0.75 * (beta[2] - beta[1])
            @test Statistics.mean(pg, :x) == 4.0
            @test Statistics.std(pg, p -> p.x) ≈ sqrt(3.0)
            @test Statistics.cov(pg, :x, beta_x) ≈ expected
            @test covariance_matrix(pg, :x, beta_x)[1, 2] ≈ expected
        end

        @testset "masked weighted statistics" begin
            pg = single_particle(x=1.0, weight=1.0) +
                 single_particle(x=20.0, weight=5.0, status=0) +
                 single_particle(x=5.0, weight=3.0)
            alive = pg.status .== 1
            double_x = p -> 2 .* p.x
            @test mean(pg, :x; mask=alive) == 4.0
            @test std(pg, double_x; mask=alive) ≈ sqrt(12.0)
            @test cov(pg, :x, double_x; mask=alive) ≈ 6.0
            @test covariance_matrix(pg, :x, double_x; mask=alive) ≈ [3 6; 6 12]
            @test mean(pg, p -> delta(p, :x); mask=alive) ≈ 0.0
            subset = view(pg, [3, 2])
            @test mean(subset, :x; mask=subset.status .== 1) == 5.0
            @test_throws DimensionMismatch mean(pg, :x; mask=[true])
            @test_throws ArgumentError mean(pg, :x; mask=falses(length(pg)))
        end

        @testset "delta" begin
            pg = _STAT_BEAM
            d = delta(pg, :x)
            @test length(d) == length(pg)
            @test abs(_mean(d, weights(pg.weight))) < 1e-12
        end

        @testset "ptp" begin
            @test ptp([1.0, 5.0, 3.0]) == 4.0

            pg = _STAT_BEAM
            @test ptp(pg, :x) ≈ maximum(pg.x) - minimum(pg.x)
        end
    end

    @testset "Emittance" begin
        @testset "norm_emit_calc single plane" begin
            pg = _STAT_BEAM
            ex = norm_emit_calc(pg, (:x,))
            ey = norm_emit_calc(pg, (:y,))
            @test ex > 0
            @test ey > 0
        end

        @testset "norm_emit_calc 4D" begin
            pg = _STAT_BEAM
            e4 = norm_emit_calc(pg, (:x, :y))
            @test e4 > 0
        end

        @testset "convenience wrappers consistency" begin
            pg = _STAT_BEAM
            @test norm_emit_x(pg) ≈ norm_emit_calc(pg, (:x,))
            @test norm_emit_y(pg) ≈ norm_emit_calc(pg, (:y,))
            @test norm_emit_4d(pg) ≈ norm_emit_calc(pg, (:x, :y))
        end
    end

    @testset "Twiss" begin
        @testset "twiss_calc known sigma" begin
            beta0  = 10.0
            alpha0 = -1.5
            emit0  = 2e-6
            gamma0 = (1 + alpha0^2) / beta0

            sigma = emit0 * [beta0  -alpha0;
                             -alpha0 gamma0]

            tw = twiss_calc(sigma)

            @test tw.beta  ≈ beta0   rtol=1e-12
            @test tw.alpha ≈ alpha0  rtol=1e-12
            @test tw.gamma ≈ gamma0  rtol=1e-12
            @test tw.emit  ≈ emit0   rtol=1e-12
        end

        @testset "twiss_calc identity sigma" begin
            tw = twiss_calc(Float64[1 0; 0 1])
            @test tw.emit  ≈ 1.0
            @test tw.beta  ≈ 1.0
            @test tw.gamma ≈ 1.0
            @test tw.alpha ≈ 0.0
        end

        @testset "twiss_calc shape assertion" begin
            @test_throws AssertionError twiss_calc(ones(3, 3))
        end

        @testset "twiss_dispersion matrix — zero dispersion" begin
            beta0  = 5.0
            alpha0 = -0.5
            emit0  = 1e-6
            gamma0 = (1 + alpha0^2) / beta0
            delta2 = 1e-4

            sigma3 = zeros(3, 3)
            sigma3[1:2, 1:2] .= emit0 * [beta0 -alpha0; -alpha0 gamma0]
            sigma3[3, 3] = delta2

            tw = twiss_dispersion(sigma3)
            @test tw.beta  ≈ beta0   rtol=1e-10
            @test tw.alpha ≈ alpha0  rtol=1e-10
            @test tw.emit  ≈ emit0   rtol=1e-10
            @test tw.eta   ≈ 0.0     atol=1e-14
            @test tw.etap  ≈ 0.0     atol=1e-14
        end

        @testset "twiss_dispersion matrix — nonzero dispersion" begin
            beta0  = 5.0
            alpha0 = 0.0
            emit0  = 1e-6
            gamma0 = (1 + alpha0^2) / beta0
            delta2 = 1e-4
            eta_val = 0.5
            etap_val = 0.1

            sigma3 = zeros(3, 3)
            sigma3[1, 1] = emit0 * beta0 + eta_val^2 * delta2
            sigma3[2, 2] = emit0 * gamma0 + etap_val^2 * delta2
            sigma3[1, 2] = -emit0 * alpha0 + eta_val * etap_val * delta2
            sigma3[2, 1] = sigma3[1, 2]
            sigma3[1, 3] = eta_val * delta2
            sigma3[3, 1] = sigma3[1, 3]
            sigma3[2, 3] = etap_val * delta2
            sigma3[3, 2] = sigma3[2, 3]
            sigma3[3, 3] = delta2

            tw = twiss_dispersion(sigma3)
            @test tw.eta   ≈ eta_val   rtol=1e-10
            @test tw.etap  ≈ etap_val  rtol=1e-10
            @test tw.beta  ≈ beta0     rtol=1e-10
            @test tw.emit  ≈ emit0     rtol=1e-10
        end

        @testset "Twiss on a ParticleGroup" begin
            pg = _STAT_BEAM
            tw = twiss(pg; plane=:x)
            @test tw.beta > 0
            @test tw.emit > 0
            @test hasproperty(tw, :norm_emit)
            @test tw.norm_emit > 0

            shifted_ref = twiss(pg; plane=:x, p0c=2 * Statistics.mean(pg, momentum))
            @test shifted_ref.norm_emit ≈ tw.norm_emit
            @test shifted_ref.eta ≈ 2 * tw.eta
            @test shifted_ref.etap ≈ tw.etap
        end

        @testset "twiss lowest-action fraction" begin
            pg = _STAT_BEAM[1:101]
            core = pg[sortperm(Jx(pg))[1:90]]
            @test all(isapprox.(values(twiss(pg; fraction=0.9)), values(twiss(core)); rtol=1e-12))
        end

        @testset "Twiss plane selection" begin
            pg = _STAT_BEAM
            p0c = Statistics.mean(pg, momentum)
            sigma_y = covariance_matrix(pg, :y, p -> p.py ./ p0c, p -> momentum(p) ./ p0c)
            @test twiss(pg; plane=:y).beta ≈ twiss_dispersion(sigma_y).beta
        end
    end

    @testset "Normal form matrices" begin
        @testset "A_mat_calc" begin
            beta  = 10.0
            alpha = -1.5

            A = A_mat_calc(beta, alpha)
            @test size(A) == (2, 2)
            @test A[1, 1] ≈ sqrt(beta)
            @test A[1, 2] ≈ 0.0
            @test A[2, 1] ≈ -alpha / sqrt(beta)
            @test A[2, 2] ≈ 1.0 / sqrt(beta)
            @test det(A) ≈ 1.0
        end

        @testset "A_inverse_mat_calc" begin
            beta  = 10.0
            alpha = -1.5

            A    = A_mat_calc(beta, alpha)
            Ainv = A_inverse_mat_calc(beta, alpha)

            @test A * Ainv ≈ Float64[1 0; 0 1]  atol=1e-12
            @test Ainv * A ≈ Float64[1 0; 0 1]  atol=1e-12
        end

        @testset "A_mat identity Twiss" begin
            A = A_mat_calc(1.0, 0.0)
            @test A ≈ Float64[1 0; 0 1]
        end
    end

    @testset "Amplitude" begin
        @testset "amplitude_calc identity Twiss" begin
            x = [1.0, 0.0, 1.0]
            p = [0.0, 1.0, 1.0]
            J = amplitude_calc(x, p; beta=1.0, alpha=0.0)
            @test J ≈ [0.5, 0.5, 1.0]
        end

        @testset "amplitude_calc with Twiss" begin
            beta  = 4.0
            alpha = -2.0
            gamma_tw = (1 + alpha^2) / beta

            x = [1.0]
            p = [0.0]
            J = amplitude_calc(x, p; beta=beta, alpha=alpha)
            @test J[1] ≈ gamma_tw / 2 * x[1]^2
        end

        @testset "particle_amplitude" begin
            pg = _STAT_BEAM
            J = particle_amplitude(pg; plane=:x)
            @test length(J) == length(pg)
            @test all(J .>= 0)
        end

        @testset "particle_amplitude with provided Twiss" begin
            pg = _STAT_BEAM
            sigma2 = StatsBase.cov(
                hcat(pg.x, pg.px ./ massof(Species("electron"))),
                weights(pg.weight), 1; corrected=false
            )
            tw = twiss_calc(sigma2)
            J1 = particle_amplitude(pg; plane=:x)
            J2 = particle_amplitude(pg; plane=:x, twiss=tw)
            @test J1 ≈ J2
        end
    end

    @testset "Normalized particle coordinate" begin
        @testset "position" begin
            pg = _STAT_BEAM
            xn = normalized_particle_coordinate(pg, :x)
            @test length(xn) == length(pg)
        end

        @testset "momentum" begin
            pg = _STAT_BEAM
            pxn = normalized_particle_coordinate(pg, :px)
            @test length(pxn) == length(pg)
        end

        @testset "Transverse actions" begin
            pg = _STAT_BEAM
            @test isapprox(2 .* Jx(pg), x_bar(pg).^2 .+ px_bar(pg).^2)
            @test isapprox(2 .* Jy(pg), y_bar(pg).^2 .+ py_bar(pg).^2)
        end
    end

    @testset "Equal-count slice views" begin
        n = 10
        z = collect(Float64, n:-1:1)
        pg = ParticleGroup(; x=zeros(n), px=zeros(n), y=zeros(n), py=zeros(n),
                           z, pz=ones(n), t=zeros(n), weight=ones(n),
                           status=ones(Int, n), species=Species("electron"))
        groups = slices(pg; n=4)
        @test map(length, groups) == [3, 2, 2, 3]
        @test vcat((s.z for s in groups)...) == reverse(z)
        groups[1].x[1] = 2.0
        @test pg.x[end] == 2.0
        @test map(s -> Statistics.mean(s, :t), slices(_STAT_BEAM; n=4)) ==
              map(s -> Statistics.mean(s, :t), slices(_STAT_BEAM; by=:t, n=4))
    end

end  # Statistics and Beam Physics
