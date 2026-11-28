using AtomicAndPhysicalConstants: Species, massof

@testset "ParticleGroup Construction" begin

    @testset "Constructors" begin
        @testset "named columns with explicit IDs" begin
            n = 5
            x  = randn(n)
            px = randn(n) .* 1e6
            y  = randn(n)
            py = randn(n) .* 1e6
            z  = randn(n)
            pz = ones(n) .* 1e9
            t  = zeros(n)
            status = ones(Int, n)
            weight = fill(1e-15, n)
            sp = Species("electron")
            id = collect(1:n)

            pg = ParticleGroup(; x, px, y, py, z, pz, t, status, weight, species=sp, id)

            @test length(pg) == n
            @test pg.x === x
            @test pg.px === px
            @test pg.y === y
            @test pg.py === py
            @test pg.z === z
            @test pg.pz === pz
            @test pg.t === t
            @test pg.status === status
            @test pg.weight === weight
            @test pg.species === sp
            @test pg.id === id
        end

        @testset "named columns with generated IDs" begin
            n = 3
            x  = [0.0, 1.0, 2.0]
            px = [1e6, 2e6, 3e6]
            y  = zeros(n)
            py = zeros(n)
            z  = zeros(n)
            pz = fill(1e9, n)
            t  = zeros(n)
            status = ones(Int, n)
            weight = fill(1e-15, n)
            sp = Species("electron")

            pg = ParticleGroup(; x, px, y, py, z, pz, t, status, weight, species=sp)

            @test length(pg) == n
            @test pg.x == x
            @test pg.id == [1, 2, 3]
        end

        @testset "column lengths match x" begin
            v = zeros(2)
            columns = (; x=v, px=v, y=v, py=v, z=v, pz=v, t=v,
                       status=[1, 1], weight=ones(2), species=Species("electron"))
            @test_throws DimensionMismatch ParticleGroup(; columns..., px=[3.0])
            @test_throws DimensionMismatch ParticleGroup(; columns..., status=[1])
            @test_throws DimensionMismatch ParticleGroup(; columns..., id=[1])
            @test_throws DimensionMismatch ParticleGroup(; columns..., extra=Dict(:spin_x => zeros(3)))
            spin = [0.5, -0.5]
            @test ParticleGroup(; columns..., extra=Dict(:spin_x => spin)).extra[:spin_x] === spin
        end

        @testset "single_particle with defaults" begin
            pg = single_particle()

            @test length(pg) == 1
            @test pg.x  == [0.0]
            @test pg.px == [0.0]
            @test pg.y  == [0.0]
            @test pg.py == [0.0]
            @test pg.z  == [0.0]
            @test pg.pz == [0.0]
            @test pg.t  == [0.0]
            @test pg.status == [1]
            @test pg.weight == [1.0]
            @test pg.id == [1]
            @test nameof(pg.species) == "electron"
        end

        @testset "single_particle with custom values" begin
            pg = single_particle(
                x=1, px=5e5, y=2e-3, py=-3e5,
                z=0.01, pz=1e9, t=1e-12,
                weight=2e-15, status=0, species=Species("proton")
            )

            @test pg.x  == [1.0]
            @test pg.px == [5e5]
            @test pg.y  == [2e-3]
            @test pg.py == [-3e5]
            @test pg.z  == [0.01]
            @test pg.pz == [1e9]
            @test pg.t  == [1e-12]
            @test pg.status == [0]
            @test pg.weight == [2e-15]
            @test nameof(pg.species) == "proton"
        end

        @testset "type parameters" begin
            pg64 = single_particle(x=1.0, pz=1e9)
            @test pg64 isa ParticleGroup{Float64,Int64}

            pg32 = ParticleGroup(;
                x=Float32[0, 1], px=Float32[1e6, 2e6], y=zeros(Float32, 2), py=zeros(Float32, 2),
                z=zeros(Float32, 2), pz=fill(1f9, 2), t=zeros(Float32, 2), status=Int32[1, 1],
                weight=fill(1f-15, 2), species=Species("electron")
            )
            @test pg32 isa ParticleGroup{Float32,Int32}
            @test length(pg32) == 2
        end
    end

    @testset "Species handling" begin
        for (name, charge_C) in (("electron", -1.602176634e-19), ("proton", 1.602176634e-19))
            pg = single_particle(species=Species(name))
            @test nameof(pg.species) == name
            @test mass(pg) == massof(Species(name))
            @test species_charge(pg) == charge_C
        end
    end

    @testset "Length and indexing" begin
        @testset "length" begin
            @test length(single_particle()) == 1

            pg3 = single_particle() + single_particle() + single_particle()
            @test length(pg3) == 3
        end

        @testset "integer-range indexing" begin
            pg = single_particle(x=1.0) + single_particle(x=2.0) + single_particle(x=3.0)

            sub = pg[1:2]
            @test sub isa ParticleGroup
            @test length(sub) == 2
            @test sub.x == [1.0, 2.0]
            @test sub.pz == pg.pz[1:2]

            pg.id[2] = 42
            particle = pg[2]
            @test particle isa ParticleGroup
            @test particle.id == [42]
            @test particle.x == [2.0]
            particle.x[1] = 5.0
            @test pg.x[2] == 2.0

            particle_view = view(pg, 1:3)[2]
            @test particle_view isa ParticleGroupView
            particle_view.x[1] = 5.0
            @test pg.x[2] == 5.0
        end

        @testset "vector indexing" begin
            pg = single_particle(x=10.0) + single_particle(x=20.0) + single_particle(x=30.0)

            sub = pg[[1, 3]]
            @test sub isa ParticleGroup
            @test length(sub) == 2
            @test sub.x == [10.0, 30.0]
        end

        @testset "BitVector indexing" begin
            pg = single_particle(x=1.0) + single_particle(x=2.0) + single_particle(x=3.0)

            sub = pg[BitVector([true, false, true])]
            @test length(sub) == 2
            @test sub.x == [1.0, 3.0]
        end
    end

    @testset "Joining (+)" begin
        @testset "two-way addition" begin
            pg1 = single_particle(x=1.0, pz=1e9)
            pg2 = single_particle(x=2.0, pz=2e9)
            pg = pg1 + pg2

            @test length(pg) == 2
            @test pg.x  == [1.0, 2.0]
            @test pg.pz == [1e9, 2e9]
            @test pg.species === pg1.species
        end

        @testset "three-way addition" begin
            pg = single_particle(x=1.0) + single_particle(x=2.0) + single_particle(x=3.0)
            @test length(pg) == 3
            @test pg.x == [1.0, 2.0, 3.0]
        end
    end

    @testset "Equality and comparison" begin
        @testset "==" begin
            pg1 = single_particle(x=1.0, pz=1e9)
            pg2 = single_particle(x=1.0, pz=1e9)
            @test pg1 == pg2

            pg3 = single_particle(x=2.0, pz=1e9)
            @test pg1 != pg3

            pg4 = single_particle(x=1.0, pz=1e9, species=Species("proton"))
            @test pg1 != pg4
        end

        @testset "isapprox" begin
            pg1 = single_particle(x=1.0, pz=1e9)
            pg2 = single_particle(x=1.0 + 1e-14, pz=1e9)
            @test isapprox(pg1, pg2)

            pg3 = single_particle(x=2.0, pz=1e9)
            @test !isapprox(pg1, pg3)
        end
    end

    @testset "show / display" begin
        pg = single_particle(weight=1e-15)
        s = sprint(show, pg)
        @test occursin("ParticleGroup", s)
        @test occursin("1 particles", s)
    end

end  # ParticleGroup Construction
