using AtomicAndPhysicalConstants: C_LIGHT

@testset "ParticleGroupView basic behavior" begin
    # Construct a small ParticleGroup with three particles
    pg1 = single_particle(x = 0.0)
    pg2 = single_particle(x = 1.0, px = 1e5, py = 2e5, pz = 1e6)
    pg3 = single_particle(x = 2.0, px = -2e5, py = 1e5, pz = 2e6)
    pg = pg1 + pg2 + pg3

    @test length(pg) == 3

    # Create a view onto a subset of particles
    view_inds = 2:3
    pgv = view(pg, view_inds)

    @test length(pgv) == length(view_inds)
    @test all(pgv.x .== pg.x[view_inds])
    @test all(pgv.y .== pg.y[view_inds])
    @test all(pgv.z .== pg.z[view_inds])

    # Mutations through the view should reflect in the parent
    pgv.x .= 42.0
    @test all(pg.x[view_inds] .== 42.0)

    # Derived quantities should be consistent between parent and view
    @test all(momentum(pgv) .== momentum(pg)[view_inds])
    @test all(energy(pgv) .== energy(pg)[view_inds])

    # Drift only the viewed subset and ensure other particles are unchanged
    for dt in (1e-12, [1e-12, 2e-12])
        before = copy(pgv)
        drift!(pgv, dt)
        @test pgv.x .- before.x ≈ beta_x(before) .* C_LIGHT .* dt
        @test pgv.y .- before.y ≈ beta_y(before) .* C_LIGHT .* dt
        @test pgv.z .- before.z ≈ beta_z(before) .* C_LIGHT .* dt
        @test pgv.t ≈ before.t .+ dt
        @test pg[1] == pg1
    end

    before = copy(pg)
    @test_throws ArgumentError drift!(pgv, [1e-12])
    @test pg == before
end

