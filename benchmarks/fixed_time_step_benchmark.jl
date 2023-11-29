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



function fixed_time_step_benchmark(T, n, ms)

    fixed_dim_bm = [lyapunov_benchmark(T, n, m) for m in ms]
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
            xlabel = "dimension of problem",
        )
        lines!(axmin, ms, qr_min_times, label = "qr")
        lines!(axmin, ms, qrwoq_min_times, label = "qrwoq")
        axislegend(axmin; position = :lt)

        axmed = Axis(
            fig[1, 2],
            title = "Median",
            ylabel = "ms",
            xlabel = "dimension of problem",
        )
        lines!(axmed, ms, qr_med_times)
        lines!(axmed, ms, qrwoq_med_times)

        Label(fig[0, :], "Execution times - fixed step")
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
            xlabel = "dimension of problem",
        )
        lines!(axnallocs, ms, qr_nallocs, label = "qr")
        lines!(axnallocs, ms, qrwoq_nallocs, label = "qrwoq")
        axislegend(axnallocs; position = :lt)

        axmem = Axis(
            fig[1, 2],
            title = "memory",
            ylabel = "KiB",
            xlabel = "dimension of problem",
        )
        lines!(axmem, ms, qr_mem)
        lines!(axmem, ms, qrwoq_mem)

        Label(fig[0, :], "Allocations - fixed step")
        fig
    end
    display(fig_allocs)

    save(joinpath(@__DIR__, "fixed_time_step_times.pdf"), fig_times, pt_per_unit=1.0)
    save(joinpath(@__DIR__, "fixed_time_step_allocs.pdf"), fig_times, pt_per_unit=1.0)


    return qr_results, qrwoq_results, fig_times, fig_allocs
end

T = Float64
n = 500
ms = 1:10
qr_results, qrwoq_results, fig_times, fig_allocs = fixed_time_step_benchmark(T, n, ms)



