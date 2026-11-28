import Plots

@testset "Slice statistics in particle plots" begin
    pg = reduce(+, [single_particle(z=i, x=i * 1e-6) for i in 1:4])
    my_func(g) = mean(g, :x)
    raw = [1.5e-6, 3.5e-6]

    for (plotter, coords, axis) in ((slice_plot, (), 1), (density_and_slice_plot, (:z, :x), 2))
        for (statistic, label) in ((my_func, "my_func"), ("my_func (1/J)" => my_func, "my_func (1/J)"))
            p = plotter(pg, coords..., statistic; n_slice=2)
            curve = p[axis].series_list[1]
            @test curve[:x] ≈ [1.5, 3.5]
            @test curve[:y] ≈ raw
            @test curve[:label] == label
            @test p[axis][:yaxis][:guide] == label
        end

        known = "mean_x" => my_func
        p = plotter(pg, coords..., known; n_slice=2, tex=false)
        @test p[axis].series_list[1][:y] ≈ raw ./ 1e-6
        @test p[axis].series_list[1][:label] == "mean_x (µm)"
        @test p[axis][:yaxis][:guide] == "mean_x (µm)"

        p = plotter(pg, coords..., known, my_func; n_slice=2)
        @test p[axis].series_list[1][:y] ≈ raw
        @test p[axis].series_list[2][:y] ≈ raw
    end
    @test_throws MethodError density_and_slice_plot(pg, :z, :x)
end
