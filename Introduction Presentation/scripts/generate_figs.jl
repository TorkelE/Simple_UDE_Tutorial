### Preparations ###

# Fetch packages.
begin
    using Catalyst
    using OrdinaryDiffEqDefault
    using Plots
    import Plots: mm
end

# Plotting defaults.
default(xlimit = (0.0, 30.0), ylimit = (0.0, 100.0), framestyle = :box, legend = :none, xguide = "",
    formatter = :none, bg = :white, bg_inside = :transparent)

# Declare the model.
sir = @reaction_network begin
    β, S + I --> 2I
    γ, I --> R
end

# Create a basic ODEProblem
begin
    u0 = [:S => 99, :I => 1, :R => 0]
    ps = [:β => 0.01, :γ => 0.1]
    tspan = (0.0, 30.0)
    prob = ODEProblem(sir, u0, tspan, ps)
    plot(solve(prob); lw = 5, la = 0.8)
end

# For saving figures.
function save_figure(plt, filename)
    savefig(plt, joinpath(@__DIR__, "figure_base", filename * ".svg"))
    plt
end

### Create Figures ###

# Figure 1 (initial data).
begin
    sol = solve(remake(prob, p = [:β => 0.01, :γ => 0.1]))
    sample_ts = 0:1:6
    sample_I = [I * (0.75 + 0.5rand()) + 5*(rand()-0.5) for I in sol(sample_ts; idxs = :I).u]
    fig_1 = plot(sample_ts, sample_I; seriestype = :scatter, ms = 8, color = 2)
    save_figure(fig_1, "figure_1")
end

# Figure 1 (potential future data).
begin
    sample_ts_cont = 7:1:30
    plt_1 = begin
        sample_I_cont1 = [1.05 * (I * (0.75 + 0.5rand()) + 8*(rand()-0.5)) for I in sol(sample_ts_cont; idxs = :I).u]
        plot(sample_ts, sample_I; seriestype = :scatter, ms = 8, xlimit = (0.0, 30.0), ylimit = (0.0, 100.0), color = 2)
        plt1 = plot!(sample_ts_cont, sample_I_cont1; seriestype = :scatter, alpha = 0.7, marker = :diamond, ms = 8, xlimit = (0.0, 30.0), ylimit = (0.0, 100.0), color = 2)
    end
    plt_2 = begin
        sample_I_cont2 = [I * (0.75 + 0.5rand()) + 5*(rand()-0.5) for I in sol(sample_ts_cont; idxs = :I).u]
        plot(sample_ts, sample_I; seriestype = :scatter, ms = 8, xlimit = (0.0, 30.0), ylimit = (0.0, 100.0), color = 2)
        plt2 = plot!(sample_ts_cont, sample_I_cont2; seriestype = :scatter, alpha = 0.7, marker = :rect, ms = 8, xlimit = (0.0, 30.0), ylimit = (0.0, 100.0), color = 2)
    end
    plt_3 = begin
        sample_I_cont3 = [(I * (0.75 + 0.5rand()) + 5*(rand()-0.5))/1.2 for I in sol(sample_ts_cont; idxs = :I).u]
        plot(sample_ts, sample_I; seriestype = :scatter, ms = 8, xlimit = (0.0, 30.0), ylimit = (0.0, 100.0), color = 2)
        plt3 = plot!(sample_ts_cont, sample_I_cont3; seriestype = :scatter, alpha = 0.7, marker = :star5, ms = 8, xlimit = (0.0, 30.0), ylimit = (0.0, 100.0), color = 2)
    end
    fig_1_alts = plot(plt1, plt2, plt3; layout = (1, 3), size = (1400, 400))
    save_figure(fig_1_alts, "figure_1_alts")
end

# Figure 2 (basic simulation).
begin
    plot(solve(prob); lw = 8, la = 0.8)
    scatter!((0.2, 98.8), ms = 12, color = 1, alpha = 0.8)
    scatter!((0.2, 1.2), ms = 12, color = 2, alpha = 0.9)
    scatter!((0.2, 0.2), ms = 12, color = 3, alpha = 0.6)
    fig_2 = plot!(xguide = "")
    save_figure(fig_2, "figure_2")
end

# Figure 3 (varying parameters => different simulation outcomes).
begin
    prob1 = remake(prob, p = [:β => 0.01, :γ => 0.16])
    prob2 = remake(prob, p = [:β => 0.0125, :γ => 0.12])
    prob3 = remake(prob, p = [:β => 0.015, :γ => 0.08])
    plt1 = plot(solve(prob1); lw = 8, la = 0.8)
    plt2 = plot(solve(prob2); lw = 8, la = 0.8)
    plt3 = plot(solve(prob3); lw = 8, la = 0.8)
    fig_3 = plot(plt1, plt2, plt3; layout = (1, 3), size = (1400, 400), xguide = "")
    save_figure(fig_3, "figure_3")
end

# Figure 4 (lots of simulations for different parameter sets).
begin
    oprobs = [remake(prob, p = [:β => β, :γ => γ]) for γ in 0.05:0.05:0.25, β in 0.005:0.005:0.025]
    sols = solve.(oprobs)
    function make_plot(sol)
        plot(sample_ts, sample_I; seriestype = :scatter, ms = 8, xlimit = (0.0, 30.0), ylimit = (0.0, 100.0), color = 2)
        plot!(sol; lw = 5, la = 0.9, xguide = "", bottom_margin = -4.5mm, left_margin = -5mm,
            color = [1 2 3])
    end
    fig_4 = plot(make_plot.(sols)..., layout = (5,5), size = (2000, 1200))
    save_figure(fig_4, "figure_4")
end

# Figure 5 (incorrect development).
begin
    prob_best = remake(prob, p = [:β => 0.01, :γ => 0.15])
    plot(solve(prob_best); lw = 8, la = 0.8, xguide = "")
    plot!(sample_ts, sample_I; seriestype = :scatter, ms = 8, xlimit = (0.0, 30.0), ylimit = (0.0, 100.0), color = 2)
    fig_5 = plot!(sample_ts_cont, sample_I_cont2; seriestype = :scatter, alpha = 0.7, marker = :rect, ms = 8, xlimit = (0.0, 30.0), ylimit = (0.0, 100.0), color = 2)
    save_figure(fig_5, "figure_5")
end

# Figure 6 (similar parameters gives different outcome).
begin
    prob_trueish = remake(prob, p = [:β => 0.0105, :γ => 0.10])
    plot(solve(prob_trueish); lw = 8, la = 0.8, xguide = "")
    plot!(sample_ts, sample_I; seriestype = :scatter, ms = 8, xlimit = (0.0, 30.0), ylimit = (0.0, 100.0), color = 2)
    fig_6_true = plot!(sample_ts_cont, sample_I_cont2; seriestype = :scatter, alpha = 0.7, marker = :rect, ms = 8, xlimit = (0.0, 30.0), ylimit = (0.0, 100.0), color = 2)
    fig_6 = plot(fig_5, fig_6_true; layout = (1, 2), size = (1200, 400))
    save_figure(fig_6, "figure_6")
end

# Figure 7 (profile likelihood intro).
begin
    get_loss(sol) = sum(abs2, sol(sample_ts, idxs = :I).u .- sample_I)
    losses = get_loss.(sols)
    loss_mins = [minimum(losses[i,:]) for i in 1:5]
    loss_args = [argmin(losses[i,:]) for i in 1:5]
    oprobs_dense = [remake(prob, p = [:β => β, :γ => γ]) for γ in 0.05:0.005:0.25, β in 0.005:0.005:0.025]
    sols_dense = solve.(oprobs_dense)
    losses_dense = get_loss.(sols_dense)
    loss_mins_dense = [minimum(losses_dense[i,:]) for i in 1:size(losses_dense, 1)]
    loss_args_dense = [argmin(losses_dense[i,:]) for i in 1:size(losses_dense, 1)]
    plot(0.05:0.005:0.25, -loss_mins_dense; ylimit = :auto, lw = 12, color = 4, la = 0.8, size = (1600,500), xlimit = (0.05, 0.25))
    plot_7_upper = plot!(0.05:0.05:0.25, -loss_mins; seriestype = :scatter, ms = 15, color = 4, marker = :oct, xlimit = (0.049, 0.251), xticks = [0.1,0.15,0.2], yticks = [])
    save_figure(plot_7_upper, "plot_7_upper")
end
begin
    plot_7_lower = plot(make_plot.(sols)..., layout = (5,5), size = (2000, 900))
    save_figure(plot_7_lower, "plot_7_lower")
end
