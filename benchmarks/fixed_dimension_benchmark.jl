# activates the script environment and instantiates it
import Pkg
Base.active_project() != joinpath(@__DIR__, "Project.toml") && Pkg.activate(@__DIR__)
haskey(Pkg.project().dependencies, "QRWithoutQ") ||
    Pkg.develop(path = joinpath(@__DIR__, "../"))
#isfile(joinpath(@__DIR__, "Manifest.toml")) && Pkg.resolve()
#Pkg.instantiate()

using LinearAlgebra, Revise, QRWithoutQ

using BenchmarkTools, CairoMakie

includet("utils.jl")



function fixed_dimension_benchmark(T, ns, m)

    fixed_dim_bm = [lyapunov_benchmark(T, n, m) for n in ns]
    qr_results = [bm[1] for bm in fixed_dim_bm]
    qrwoq_results = [bm[2] for bm in fixed_dim_bm]

    nan2mili = 1e-6
    byte2KiB = 0.000976562

    # execution times
    qr_min_times = nan2mili * [minimum(time(trial)) for trial in qr_results]
    qr_med_times = nan2mili * [median(time(trial)) for trial in qr_results]

    qrwoq_min_times = nan2mili * [minimum(time(trial)) for trial in qrwoq_results]
    qrwoq_med_times = nan2mili * [median(time(trial)) for trial in qrwoq_results]


    fig_times = with_theme(Theme(Lines = (cycle = Cycle([:linestyle]),))) do

        fig = Figure(size = (1200, 600 / 1.61))
        axmin = Axis(
            fig[1, 1],
            title = "Minimum",
            ylabel = "ms",
            xlabel = "number of factorizations",
        )
        lines!(axmin, ns, qr_min_times, label = "qr")
        lines!(axmin, ns, qrwoq_min_times, label = "qrwoq")
        axislegend(axmin; position = :lt)

        axmed = Axis(
            fig[1, 2],
            title = "Median",
            ylabel = "ms",
            xlabel = "number of factorizations",
        )
        lines!(axmed, ns, qr_med_times)
        lines!(axmed, ns, qrwoq_med_times)

        Label(fig[0, :], "Execution times - fixed dimension")
        fig
    end
    display(fig_times)

    # allocations
    qr_nallocs = [allocs(trial) for trial in qr_results]
    qr_mem = byte2KiB * [memory(trial) for trial in qr_results]

    qrwoq_nallocs = [allocs(trial) for trial in qrwoq_results]
    qrwoq_mem = byte2KiB * [memory(trial) for trial in qrwoq_results]

    fig_allocs = with_theme(Theme(Lines = (cycle = Cycle([:linestyle]),))) do

        fig = Figure(size = (1200, 600 / 1.61))
        axnallocs = Axis(
            fig[1, 1],
            title = "number of allocations",
            ylabel = "#",
            xlabel = "number of factorizations",
        )
        lines!(axnallocs, ns, qr_nallocs, label = "qr")
        lines!(axnallocs, ns, qrwoq_nallocs, label = "qrwoq")
        axislegend(axnallocs; position = :lt)

        axmem = Axis(
            fig[1, 2],
            title = "memory",
            ylabel = "KiB",
            xlabel = "number of factorizations",
        )
        lines!(axmem, ns, qr_mem)
        lines!(axmem, ns, qrwoq_mem)

        Label(fig[0, :], "Allocations - fixed dimension")
        fig
    end
    display(fig_allocs)

    save(joinpath(@__DIR__, "fixed_dim_times.pdf"), fig_times, pt_per_unit=1.0)
    save(joinpath(@__DIR__,"fixed_dim_allocs.pdf"), fig_allocs, pt_per_unit=1.0)


    return qr_results, qrwoq_results, fig_times, fig_allocs
end

T = Float64
ns = 50:50:1000
m = 10
qr_results, qrwoq_results, fig_times, fig_allocs = fixed_dimension_benchmark(T, ns, m)



