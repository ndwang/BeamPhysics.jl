using BeamPhysics
using Test
import Plots

const FIELDMESH_FILE = joinpath(@__DIR__, "lcls2_solenoid_fieldmesh.h5")

@testset "FieldMesh" begin

    @testset "Linear interpolation" begin
        fm = FieldMesh(FIELDMESH_FILE)
        itp = interpolator(fm, "Bz")
        r, theta, z = coord_vecs(fm)
        bz = fm["Bz"]
        @test itp(0, theta[1], z[end]) ≈ bz[1, 1, end]
        @test itp((r[2] + r[3]) / 2, theta[1], (z[200] + z[201]) / 2) ≈
            sum(bz[2:3, 1, 200:201]) / 4
        @test_throws BoundsError itp(r[1], theta[1] + 0.1, z[1])
        @test_throws BoundsError itp(r[1], theta[1], z[1] - 0.001)

        ex = [x + 2y + 3z + im * (x - z)
              for x in (0, 1), y in (0, 2, 4), z in (-2, 0, 2, 4)]
        rf = FieldMesh("rectangular", ("x", "y", "z"), (0.0, 0.0, -2.0),
                       (1.0, 2.0, 2.0), size(ex), "center", 1, 1.3e9, 0.25, 2.0,
                       Dict("electricField/x" => ex))
        @test interpolate(rf, "Ex", (0.25, 3, 1.5)) ≈ -2.5 - 21.5im
        @test interpolate(rf, "abs_Ex", (0.25, 0, -2)) ≈
            0.75 * abs(rf["Ex"][1, 1, 1]) + 0.25 * abs(rf["Ex"][2, 1, 1])
    end

    @testset "Construction from HDF5 file" begin
        fm = FieldMesh(FIELDMESH_FILE)
        @test fm isa FieldMesh
        @test fm.geometry == "cylindrical"
        @test fm.axis_labels == ("r", "theta", "z")
        @test fm.grid_size == (31, 1, 481)
        @test fm.ele_anchor_pt == "center"
        @test fm.harmonic == 0
        @test fm.field_scale == 0.06
    end

    @testset "Rectangular slice plots" begin
        x, y, z = 0:2, -1:2, 0:1
        bz = [xi + 10zi for xi in x, yi in y, zi in z]
        fm = FieldMesh("rectangular", ("x", "y", "z"), (0.0, -1.0, 0.0),
                       (1.0, 1.0, 1.0), size(bz), "center", 0, 0.0, 0.0, 1.0,
                       Dict("magneticField/z" => bz))

        xy = fieldmesh_plot(fm; component="Bz", z=0.4)
        @test xy.series_list[1][:z].surf ≈ [xi + 4 for yi in y, xi in x]

        yz = fieldmesh_plot(fm; component="Bz", x=0.25)
        @test yz.series_list[1][:z].surf ≈ [0.25 + 10zi for yi in y, zi in z]
    end

    @testset "Grid accessors" begin
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
        @test axis_index(fm, "r") == 1
        @test axis_index(fm, "theta") == 2
        @test axis_index(fm, "z") == 3
        @test_throws ArgumentError axis_index(fm, "x")

        # coord_vec
        r = coord_vec(fm, "r")
        @test length(r) == 31
        @test r[1] ≈ 0.0
        @test r[end] ≈ 0.03

        z = coord_vec(fm, "z")
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

    @testset "Field scaling" begin
        fm = FieldMesh(FIELDMESH_FILE)

        @test scale(fm) == 0.06
        @test phase(fm) == 0.0
        @test factor(fm) == 0.06
        @test is_static(fm)
        @test frequency(fm) == 0.0
    end

    @testset "Boolean queries" begin
        fm = FieldMesh(FIELDMESH_FILE)

        @test is_static(fm)
        @test is_pure_magnetic(fm)
        @test !is_pure_electric(fm)
    end

    @testset "Component access" begin
        fm = FieldMesh(FIELDMESH_FILE)

        # Raw components
        @test haskey(fm.components, "magneticField/z")
        @test haskey(fm.components, "magneticField/r")
        @test size(fm.components["magneticField/z"]) == (31, 1, 481)

        # Scaled component
        bz = scaled_component(fm, "Bz")
        @test size(bz) == (31, 1, 481)
        @test bz ≈ fm.field_scale .* fm.components["magneticField/z"]

        # Missing component returns zeros
        ez = scaled_component(fm, "Ez")
        @test all(iszero, ez)
    end

    @testset "getindex" begin
        fm = FieldMesh(FIELDMESH_FILE)

        # Coordinate access
        z = fm["z"]
        @test z[1] ≈ -0.24
        @test length(z) == 481

        # Scaled component
        bz = fm["Bz"]
        @test size(bz) == (31, 1, 481)
        @test bz ≈ fm.field_scale .* fm.components["magneticField/z"]

        raw = fm.components["magneticField/z"]
        @test_deprecated r"unscaled field" begin
            stored = fm["magneticField/z"]
            @test stored === raw
        end
        @test raw[1, 1, 241] ≈ 1.0       # r=0, z=0
        @test bz[1, 1, 241] ≈ 0.06
        @test scaled_component(fm, "magneticField/z") ≈ bz

        # Operator prefix
        re_bz = fm["re_Bz"]
        @test re_bz ≈ real.(bz)

        # Field magnitude
        b = fm["B"]
        @test size(b) == (31, 1, 481)
    end

    @testset "Setters" begin
        fm = FieldMesh(FIELDMESH_FILE)
        orig_zmin = zmin(fm)

        # Shift z origin
        set_min!(fm, "z", -0.5)
        @test zmin(fm) == -0.5
        @test zmax(fm) ≈ -0.5 + 0.001 * 480

        # Shift z max
        fm2 = FieldMesh(FIELDMESH_FILE)
        set_max!(fm2, "z", 1.0)
        @test zmax(fm2) ≈ 1.0

        # Direct field mutation
        fm3 = FieldMesh(FIELDMESH_FILE)
        fm3.field_scale = 2.0
        @test scale(fm3) == 2.0
        @test factor(fm3) == 2.0
    end

    @testset "show and copy" begin
        fm = FieldMesh(FIELDMESH_FILE)
        s = sprint(show, fm)
        @test occursin("cylindrical", s)
        @test occursin("(31, 1, 481)", s)

        fm2 = copy(fm)
        @test fm == fm2
        fm2.field_scale = 99.0
        @test fm.field_scale != fm2.field_scale
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
            @test keys(fm2.components) == keys(fm.components)
            for k in keys(fm.components)
                @test fm2.components[k] ≈ fm.components[k]
            end
        finally
            isfile(tmpfile) && rm(tmpfile)
        end
    end
end
