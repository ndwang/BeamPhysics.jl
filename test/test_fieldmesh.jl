using BeamPhysics
using Test
import Plots

const FIELDMESH_FILE = joinpath(@__DIR__, "lcls2_solenoid_fieldmesh.h5")

@testset "FieldMesh" begin

    @testset "Linear interpolation" begin
        fm = FieldMesh(FIELDMESH_FILE)
        itp = interpolator(fm, :Bz)
        r, theta, z = coord_vecs(fm)
        bz = fm.Bz
        @test itp(0, theta[1], z[end]) ≈ bz[1, 1, end]
        @test itp((r[2] + r[3]) / 2, theta[1], (z[200] + z[201]) / 2) ≈
            sum(bz[2:3, 1, 200:201]) / 4
        @test_throws BoundsError itp(r[1], theta[1] + 0.1, z[1])
        @test_throws BoundsError itp(r[1], theta[1], z[1] - 0.001)

        ex = [x + 2y + 3z + im * (x - z)
              for x in (0, 1), y in (0, 2, 4), z in (-2, 0, 2, 4)]
        rf = FieldMesh(; geometry="rectangular", grid_origin=(0.0, 0.0, -2.0),
                       grid_spacing=(1.0, 2.0, 2.0), grid_size=size(ex), ele_anchor_pt="center",
                       harmonic=1, fundamental_frequency=1.3e9, components=Dict(:Ex => ex),
                       field_scale=2.0, rf_phase=0.25)
        @test interpolator(rf, rf.Ex)(0.25, 3, 1.5) == interpolator(rf, :Ex)(0.25, 3, 1.5)
        @test_throws DimensionMismatch interpolator(rf, zeros(2, 2, 2))
        @test_throws ArgumentError interpolator(rf, :grid_size)
        @test interpolate(rf, :Ex, (0.25, 3, 1.5)) ≈ -2.5 - 21.5im
        @test interpolate(rf, abs.(rf.Ex), (0.25, 0, -2)) ≈
            0.75 * abs(rf.Ex[1, 1, 1]) + 0.25 * abs(rf.Ex[2, 1, 1])
    end

    @testset "Construction from HDF5 file" begin
        fm = FieldMesh(FIELDMESH_FILE)
        @test fm isa FieldMesh
        @test fm.geometry == "cylindrical"
        @test fm.axis_labels == (:r, :theta, :z)
        @test fm.grid_size == (31, 1, 481)
        @test fm.ele_anchor_pt == "center"
        @test fm.harmonic == 0
        @test fm.Bz[1, 1, 241] ≈ 0.06
    end

    @testset "Rectangular slice plots" begin
        x, y, z = 0:2, -1:2, 0:1
        bz = [xi + 10zi for xi in x, yi in y, zi in z]
        fm = FieldMesh(; geometry="rectangular", grid_origin=(0.0, -1.0, 0.0),
                       grid_spacing=(1.0, 1.0, 1.0), grid_size=size(bz), ele_anchor_pt="center",
                       harmonic=0, fundamental_frequency=0.0, components=Dict(:Bz => bz))

        xy = fieldmesh_plot(fm; component=:Bz, z=0.4)
        @test xy.series_list[1][:z].surf ≈ [xi + 4 for yi in y, xi in x]

        yz = fieldmesh_plot(fm; component=:Bz, x=0.25)
        @test yz.series_list[1][:z].surf ≈ [0.25 + 10zi for yi in y, zi in z]
    end

    @testset "Grid accessors" begin
        rf = FieldMesh(joinpath(@__DIR__, "fixtures", "openpmd_python_rf.h5"))
        @test mins(rf) == map(first, coord_vecs(rf))
        @test maxs(rf) == map(last, coord_vecs(rf))
        @test interpolate(rf, :Ex, maxs(rf)) ≈ rf.Ex[end, end, end]

        fm = FieldMesh(FIELDMESH_FILE)

        # mins, maxs, deltas
        m = mins(fm)
        @test m == (0.0, 0.0, -0.24)
        d = deltas(fm)
        @test d == (0.001, 0.0, 0.001)
        mx = maxs(fm)
        @test mx[1] ≈ 0.03
        @test mx[3] ≈ 0.24

        # axis_index
        @test axis_index(fm, :r) == 1
        @test axis_index(fm, :theta) == 2
        @test axis_index(fm, :z) == 3
        @test_throws ArgumentError axis_index(fm, :x)

        # coord_vec
        r = coord_vec(fm, :r)
        @test length(r) == 31
        @test r[1] ≈ 0.0
        @test r[end] ≈ 0.03

        z = coord_vec(fm, :z)
        @test length(z) == 481
        @test z[1] ≈ -0.24
        @test z[end] ≈ 0.24

        # coord_vecs
        vecs = coord_vecs(fm)
        @test length(vecs) == 3
        @test length(vecs[1]) == 31
        @test length(vecs[3]) == 481

        # Per-axis accessors
        @test dr(fm) == 0.001
        @test dz(fm) == 0.001
        @test rmin(fm) == 0.0
        @test rmax(fm) ≈ 0.03
        @test zmin(fm) == -0.24
        @test zmax(fm) ≈ 0.24
    end

    @testset "Boolean queries" begin
        fm = FieldMesh(FIELDMESH_FILE)

        @test is_static(fm)
        @test frequency(fm) == 0.0
        @test is_pure_magnetic(fm)
        @test !is_pure_electric(fm)
    end

    @testset "Native physical components" begin
        fm = FieldMesh(FIELDMESH_FILE)
        @test hasproperty(fm, :Br) && hasproperty(fm, :Ez)
        @test !hasproperty(fm, :Bx)
        @test_throws ArgumentError fm.Bx
        @test_throws ArgumentError fm.weight
        @test size(fm.Bz) == fm.grid_size
        @test fm.components[:Bz][1, 1, 241] ≈ 1.0
        @test fm.Bz[1, 1, 241] ≈ 0.06
        @test field_magnitude(fm, :B) ≈ sqrt.(abs2.(fm.Br) .+ abs2.(fm.Btheta) .+ abs2.(fm.Bz))
        @test all(iszero, field_magnitude(fm, :E))
        @test !haskey(fm.components, :Ez)
        untouched = copy(fm)
        @test all(iszero, fm.Ez)
        fm.Ez[1] = 0
        @test !haskey(fm.components, :Ez)
        @test fm == untouched
        fm.Ez[1] = 9.0
        @test fm.Ez[1] == 9.0
        @test field_magnitude(fm, :E)[1] == 9.0
        @test !is_pure_magnetic(fm)
        @test_throws ArgumentError field_magnitude(fm, :Q)

        grid = (; geometry="rectangular", grid_origin=(0, 0, 0), grid_spacing=(1, 1, 1),
                  grid_size=(2, 2, 2), ele_anchor_pt="center", harmonic=1, fundamental_frequency=10.0)
        mesh = FieldMesh(; grid..., components=Dict(:Ex => fill(4im, grid.grid_size)))
        @test is_pure_electric(mesh)
        @test frequency(mesh) == 10.0
        @test all(==(4), field_magnitude(mesh, :E))
        mesh.By[1] = im
        @test mesh.By[1] == im
        @test_throws DimensionMismatch FieldMesh(; grid..., components=Dict(:Ex => zeros(1, 1, 1)))
        @test_throws ArgumentError FieldMesh(; grid..., components=Dict(:Er => zeros(2, 2, 2)))
    end

    @testset "Changing the solenoid scale" begin
        fm = FieldMesh(FIELDMESH_FILE)
        bz = fm.Bz
        axis_profile = @view bz[1, 1, :]
        snapshot = copy(bz)
        center = CartesianIndex(1, 1, 241)
        @test fm.field_scale ≈ 0.06
        @test snapshot[center] ≈ 0.06

        fm.field_scale = 0.12
        @test bz[center] ≈ 0.12
        @test axis_profile[241] ≈ 0.12
        @test snapshot[center] ≈ 0.06
        @test interpolate(fm, :Bz, (0, 0, 0)) ≈ 0.12
        @test field_magnitude(fm, :B)[center] ≈ 0.12

        bz[center] = 0.24
        @test fm.components[:Bz][center] ≈ 2.0
        fm.field_scale = 0
        @test size(bz) == fm.grid_size
        @test all(iszero, bz)
        @test all(iszero, field_magnitude(fm, :B))
        bz[center] = 0
        @test_throws ArgumentError (bz[center] = 1.0)
        fm.field_scale = 0.06
        @test bz[center] ≈ 0.12

        fm.rf_phase = 0.25
        @test eltype(bz) === Float64
        @test iszero(axis_profile[241])
        @test eltype(fm.Bz) === ComplexF64
        @test collect(fm.Bz)[center] ≈ -0.12im
        @test only(fieldmesh_plot_onaxis(fm).series_list)[:y][241] ≈ -0.12
        @test interpolate(fm, :Bz, (0, 0, 0)) ≈ -0.12im
    end

    @testset "Setters" begin
        fm = FieldMesh(FIELDMESH_FILE)
        orig_zmin = zmin(fm)

        # Shift z origin
        set_min!(fm, :z, -0.5)
        @test zmin(fm) == -0.5
        @test zmax(fm) ≈ -0.5 + 0.001 * 480

        # Shift z max
        fm2 = FieldMesh(FIELDMESH_FILE)
        set_max!(fm2, :z, 1.0)
        @test zmax(fm2) ≈ 1.0

    end

    @testset "show and copy" begin
        fm = FieldMesh(FIELDMESH_FILE)
        s = sprint(show, fm)
        @test occursin("cylindrical", s)
        @test occursin("(31, 1, 481)", s)

        fm2 = copy(fm)
        @test fm == fm2
        fm2.Bz[1] = 99.0
        @test fm.Bz[1] != fm2.Bz[1]
    end

    @testset "HDF5 round-trip" begin
        fm = FieldMesh(FIELDMESH_FILE)

        tmpfile = tempname() * ".h5"
        try
            write_fieldmesh(tmpfile, fm)
            fm2 = FieldMesh(tmpfile)
            @test fm2.geometry == fm.geometry
            @test fm2.axis_labels == fm.axis_labels
            @test fm2.grid_origin == fm.grid_origin
            @test fm2.grid_spacing == fm.grid_spacing
            @test fm2.grid_size == fm.grid_size
            @test fm2.ele_anchor_pt == fm.ele_anchor_pt
            @test fm2.harmonic == fm.harmonic
            @test fm2.field_scale == fm.field_scale
            @test fm2.rf_phase == fm.rf_phase
            @test fm2.metadata == fm.metadata
            @test keys(fm2.components) == keys(fm.components)
            for k in keys(fm.components)
                @test fm2.components[k] ≈ fm.components[k]
            end
            # Detect an invalid replacement before opening the existing file for writing.
            fm.components[:Bz] = zeros(1, 1, 1)
            @test_throws DimensionMismatch write_fieldmesh(tmpfile, fm)
            @test FieldMesh(tmpfile) == fm2
        finally
            isfile(tmpfile) && rm(tmpfile)
        end
    end
end

@testset "On-axis profiles and RF quadratures" begin
    x, y, z = (-1, 1), (-1, 1), 0:2
    field = [10 + 3xi + 2yi + zi for xi in x, yi in y, zi in z]
    fm = FieldMesh(; geometry="rectangular", grid_origin=(-1., -1., 0.), grid_spacing=(2., 2., 1.),
                   grid_size=size(field), ele_anchor_pt="center", harmonic=1, fundamental_frequency=1e9,
                   components=Dict(:Ex => im .* field, :Bz => (1+2im) .* field))
    profile = fieldmesh_plot_onaxis(fm)
    @test length(profile.series_list) == 3
    ex_im, bz_re, bz_im = profile.series_list
    @test ex_im[:y] ≈ 10 .+ z
    @test bz_re[:y] ≈ 10 .+ z
    @test bz_im[:y] ≈ 2 .* (10 .+ z)
    @test startswith(ex_im[:label], "Im ") && occursin("V/m", ex_im[:label])
    @test startswith(bz_re[:label], "Re ") && startswith(bz_im[:label], "Im ")

    @test fm.field_scale == 1 && fm.rf_phase == 0
    for op in (real, imag, abs)
        cut = fieldmesh_plot(fm; component=:Bz, transform=op, z=1)
        @test cut.series_list[1][:z].surf ≈ [op((1+2im)*(11 + 3xi + 2yi)) for yi in y, xi in x]
    end

    bz = fm.Bz
    fm.field_scale = 2
    fm.rf_phase = 0.25      # one quarter turn, multiplier -2im
    @test bz ≈ (4-2im) .* field
    @test field_magnitude(fm, :B) ≈ abs.(bz)
    @test interpolate(fm, :Bz, (0, 0, 1)) ≈ 11 * (4-2im)
    rf_cut = fieldmesh_plot(fm; component=:Bz, transform=imag, z=1)
    @test rf_cut.series_list[1][:z].surf ≈ [-2*(11 + 3xi + 2yi) for yi in y, xi in x]
    ex_re, bz_re, bz_im = fieldmesh_plot_onaxis(fm).series_list
    @test ex_re[:y] ≈ 2 .* (10 .+ z)
    @test bz_re[:y] ≈ 4 .* (10 .+ z)
    @test bz_im[:y] ≈ -2 .* (10 .+ z)
    bz[1] = 6 - 8im
    @test fm.components[:Bz][1] ≈ 4 + 3im

    # Opposite phases at neighboring samples have the same magnitude.
    fm.Ex[1, :, :] .= -im
    fm.Ex[2, :, :] .= im
    magnitude = fieldmesh_plot(fm; component=:Ex, transform=abs, x=0)
    @test magnitude.series_list[1][:z].surf == ones(length(y), length(z))

    cylindrical = FieldMesh(; geometry="cylindrical", grid_origin=(0., 0., 0.), grid_spacing=(1., 0., 1.),
                            grid_size=(2, 1, 3), ele_anchor_pt="center", harmonic=1, fundamental_frequency=1e9,
                            components=Dict(:Ez => [im*(1+r+zi) for r in 0:1, t in 1:1, zi in z]))
    cylindrical.Bz  # inspecting an omitted component should not add a trace
    cyl_series = only(fieldmesh_plot_onaxis(cylindrical).series_list)
    @test cyl_series[:y] ≈ 1 .+ z
    @test startswith(cyl_series[:label], "Im ")
    @test fieldmesh_plot(cylindrical; component=:Ez, transform=imag).series_list[1][:z].surf ≈
          [1+r+zi for r in 0:1, zi in z]
    @test fieldmesh_plot(cylindrical).series_list[1][:z].surf ≈
          [1+r+zi for r in 0:1, zi in z]

    cylindrical.Br[2, 1, :] .= im .* (1 .+ z)
    radial = fieldmesh_plot(cylindrical; component=:Br, transform=imag, mirror=true)
    @test radial.series_list[1][:z].surf[1, :] ≈ -(1 .+ z)
    radial_abs = fieldmesh_plot(cylindrical; component=:Br, transform=abs, mirror=true)
    @test radial_abs.series_list[1][:z].surf[1, :] ≈ 1 .+ z
end
