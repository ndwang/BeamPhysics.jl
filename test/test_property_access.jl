using AtomicAndPhysicalConstants: Species, massof, E_CHARGE

const _PROP_BEAM = ParticleGroup(joinpath(@__DIR__, "initial_particles.h5"))

@testset "Property Access" begin

    @testset "Standalone functions" begin
        @testset "nalive / ndead" begin
            pg = single_particle() + single_particle() + single_particle()
            @test nalive(pg) == 3
            @test ndead(pg) == 0

            pg.status[1] = 0
            pg.status[2] = 0
            @test nalive(pg) == 1
            @test ndead(pg) == 2
        end

        @testset "mass / species_charge / charge" begin
            pg = _PROP_BEAM
            @test mass(pg) == massof(Species("electron"))
            @test species_charge(pg) == -E_CHARGE
            @test charge(pg) ≈ sum(pg.weight)
        end

        @testset "energy / momentum / kinetic_energy" begin
            pz = 1e9
            pg = single_particle(pz=pz)
            m = massof(Species("electron"))

            @test momentum(pg)[1] ≈ abs(pz)
            @test energy(pg)[1]   ≈ sqrt(pz^2 + m^2)
            @test kinetic_energy(pg)[1] ≈ energy(pg)[1] - m
        end

        @testset "gamma / beta / beta_x / beta_y / beta_z" begin
            pz = 1e9
            pg = single_particle(px=1e5, py=-2e5, pz=pz)
            m = massof(Species("electron"))
            E = energy(pg)[1]
            p = momentum(pg)[1]

            @test gamma(pg)[1]  ≈ E / m
            @test beta(pg)[1]   ≈ p / E
            @test beta_x(pg)[1] ≈ 1e5 / E
            @test beta_y(pg)[1] ≈ -2e5 / E
            @test beta_z(pg)[1] ≈ pz / E
        end

        @testset "xp / yp" begin
            pg = single_particle(px=1e6, py=-2e6, pz=1e9)
            @test xp(pg)[1] ≈ 1e6 / 1e9
            @test yp(pg)[1] ≈ -2e6 / 1e9
        end

        @testset "r / theta" begin
            pg = single_particle(x=3.0, y=4.0)
            @test r(pg)[1]     ≈ 5.0
            @test theta(pg)[1] ≈ atan(4.0, 3.0)

            pg0 = single_particle(x=1.0, y=0.0)
            @test theta(pg0)[1] ≈ 0.0
        end

        @testset "pr / ptheta" begin
            # Particle on x-axis: pr = px, ptheta = py
            pg = single_particle(x=1.0, y=0.0, px=3e6, py=4e6)
            @test pr(pg)[1]     ≈ 3e6  atol=1e-6
            @test ptheta(pg)[1] ≈ 4e6  atol=1e-6

            # Particle on y-axis: pr = py, ptheta = -px
            pg2 = single_particle(x=0.0, y=1.0, px=3e6, py=4e6)
            @test pr(pg2)[1]     ≈ 4e6  atol=1e-6
            @test ptheta(pg2)[1] ≈ -3e6 atol=1e-6
        end

        @testset "Lz" begin
            pg = single_particle(x=1.0, py=2.0, y=3.0, px=4.0)
            @test Lz(pg)[1] ≈ 1.0 * 2.0 - 3.0 * 4.0
        end

        @testset "set_charge!" begin
            pg = single_particle(weight=1e-15) + single_particle(weight=2e-15)
            new_q = 5e-15
            set_charge!(pg, new_q)
            @test charge(pg) ≈ new_q

            @test_throws ArgumentError set_charge!(pg, -1.0)
        end

        @testset "average_current" begin
            pg = _PROP_BEAM
            I = average_current(pg)
            @test I ≈ charge(pg) / ptp(pg.t)
        end
    end

    @testset "Emittance and normalized coordinates" begin
        @testset "norm_emit_x / norm_emit_y / norm_emit_4d" begin
            pg = _PROP_BEAM
            ex = norm_emit_x(pg)
            ey = norm_emit_y(pg)
            e4 = norm_emit_4d(pg)
            @test ex > 0
            @test ey > 0
            @test e4 > 0
        end

        @testset "x_bar, px_bar, Jx, Jy" begin
            pg = _PROP_BEAM

            xb  = x_bar(pg)
            pxb = px_bar(pg)
            @test length(xb)  == length(pg)
            @test length(pxb) == length(pg)

            jx = Jx(pg)
            jy = Jy(pg)
            @test all(jx .>= 0)
            @test all(jy .>= 0)
        end
    end

    @testset "Coordinate checks" begin
        @test in_z_coordinates(_PROP_BEAM)
        @test !in_t_coordinates(_PROP_BEAM)

        pg_t = single_particle(t=0.0, z=0.1) + single_particle(t=0.0, z=0.2)
        @test in_t_coordinates(pg_t)
        @test !in_z_coordinates(pg_t)
    end

end  # Property Access
