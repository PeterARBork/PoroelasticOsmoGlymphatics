using Plots

function diskshape(x, y, r, disccount)
    θ = LinRange(0, 2*π, disccount)
    x .+ r * sin.(θ), y .+ r * cos.(θ)
end

function plotannularparamdeformation(parameters, baselineparams; varargs...)
    parameters = merge(parameters, varargs)
    tissuecolor = RGB(254/255,234/255,204/255)

    undeformed_rs = range(parameters.a, parameters.b, length=5)
    baseline_rs = [r + calcur(r, baselineparams) for r in undeformed_rs]
    osmotic_rs = [calcur(r, parameters) for r in undeformed_rs]

    mOsm = round(convertPsiPatomOsm(parameters.Ψ), digits=1)
    comp = round(calcOsmcompression(params), digits=1)

    n = 100
    p = plot!(
        diskshape(0, 0, 2*baseline_rs[end], n),
        seriestype=[:shape, ],
        color=:lightblue,
        grid=false,
        linecolor=:lightblue,
        label="",#CSF",
        title="Osmotic brain $comp% compression due to $mOsm mOsm",
        titlefontsize=12,
        aspect_ratio=1;
        varargs...
    )
    plot!(
        diskshape(0, 0, baseline_rs[end], n),
        seriestype=[:shape, ],
        color=tissuecolor,
        alpha=1.0,
        linecolor=tissuecolor,
        label="";
        varargs...
    )
    plot!(
        diskshape(0, 0, baseline_rs[1], n),
        label="",
        color=:lightblue,
        seriestype=[:shape, ],
        linecolor=:lightblue;
        varargs...
    )
    for r in baseline_rs
        label = r == baseline_rs[1] ? "baseline" : ""
        circlearc = Plots.partialcircle(0.0, π, n, r)
        plot!(
            circlearc,
            label=label,
            color=:green,
            linewidth=2;
            varargs...
        )
    end
    for (r, u_r) in zip(undeformed_rs, osmotic_rs)
        label = r == undeformed_rs[1] ? "osmotic deformation" : ""
        circlearc = Plots.partialcircle(0.0, π, n, r + u_r)
        plot!(
            circlearc,
            label=label,
            color=2,
            linewidth=2;
            varargs...
        )
    end
    plot!(
        xlims=(-0.2e-3, 1e-3),
        ylims=(0.9e-3, 2.1e-3),
        xlabel="radius [mm]",
        ylabel="radius [mm]",
        #xticks=[0, 0.5, 1.0],
        #yticks=[1, 1.5, 2.0],
        guidefontsize=10,
        xformatter= x -> 1e3*x,
        yformatter= y -> 1e3*y;
        varargs...
    )
    annotate!(0.1e-3, 0.93e-3, text("ventricle", 9, :right); varargs...)
    annotate!(0., 2.05e-3, text("SAS", 9, :right); varargs...)

    return p
end

function addurplot!(plotparams, calcfun, units, latexsymbol; filllabel="", minparams=nothing, maxparams=nothing, yfac=1, plotvarargs...)
    p = current()
    colorindex = length(p.series_list) > 1 ? p.series_list[end][:series_index] : 1
    ylabel = L"$%$latexsymbol$ %$units"
    ylabel = length(ylabel) > 3 ? ylabel : ""

    plot!(
        x -> yfac * calcfun(x, plotparams),
        ylabel=ylabel, label="",
        xlims=(plotparams.a, plotparams.b),
        color=colorindex;
        plotvarargs...
    )
    if !isnothing(minparams) && !isnothing(maxparams)
        xs = range(plotparams.a, plotparams.b, length=40)
        minys = [yfac * calcfun(x, plotparams; minparams...) for x in xs]
        maxys = [yfac * calcfun(x, plotparams; maxparams...) for x in xs]
        plot!(
            xs,
            minys,
            fillbetween=maxys,
            alpha=0.5, linewidth=0,
            color=colorindex,
            xlims=(plotparams.a, plotparams.b),
            label=filllabel;
            plotvarargs...,
        )
    end

    return p
end

function addgammaplot!(parameters, func; filllabel="", yfac=1, varargs...)
    params = merge(parameters, varargs)
    gammas = range(0.1, 10, length=500)
    comps = yfac * [func(params; γ=γ) for γ in gammas]
    uppercomps = yfac * [func(params; γ=γ, E=params.E / 2) for γ in gammas]
    lowercomps = yfac * [func(params; γ=γ, E=params.E * 2) for γ in gammas]
    uppercomps = uppercomps[end] > lowercomps[end] ? uppercomps : lowercomps

    mOsm = round(convertPsimmHgtomOsm(params.Ψ), digits=1)
    plot!(
        1 ./ gammas,
        comps,
        label="",
        xlabel=L"$\gamma$ [mm$^{-1}$]";
        #ribbon=(lowercomps, uppercomps),;
        varargs...
    )
    plot!(
        1 ./ gammas,
        lowercomps, fillbetween=uppercomps,
        label=filllabel,
        alpha=0.5;
        varargs...
    )
end

function plotvbversusc!(params; plotPlogMestre=true, minparams=nothing, maxparams=nothing, varargs...)
    
    plot!(
        c -> -60^2 * 1e3 * calcSASv(params; Ψ=-params.RT_Pa_per_Osm * 1e-3 * c),
        xlims=(0, 50),
        xlabel=L"$\Delta c$ [mOsm]",
        ylabel=L"$|v(b)|$ [mm/h]",
        label="";
        varargs...
    )
    if !isnothing(minparams)
        plot!(
            c -> -60^2 * 1e3 * calcSASv(params; Ψ=-params.RT_Pa_per_Osm * 1e-3 * c, lowparams...),
            fillbetween=c -> -60^2 * 1e3 * calcSASv(params; Ψ=-params.RT_Pa_per_Osm * 1e-3 * c, highparams...),
            xlims=(0, 50),
            label="", color=1, alpha=0.5;
            varargs...
        )
    end
    if plotPlogMestre
        PlogMestrev_mmperh = params.PlogMestrev_mumpermin * 60 / 1000
        uncertainty = 10 * 60 / 1000
        plot!(x -> PlogMestrev_mmperh, color=:grey, label="", xlims=(0, 50); varargs...)
        plot!(
            x -> PlogMestrev_mmperh - uncertainty,
            fillbetween=x -> PlogMestrev_mmperh + uncertainty,
            xlims=(0, 50),
            color=:grey, alpha=0.5, label="";
            varargs...,
        )
    end
end

function plotcompressionversusc!(params; plotLilius=false, minparams=nothing, maxparams=nothing, varargs...)
    plot!(
        c -> calcOsmcompression(params; Ψ=-params.RT_Pa_per_Osm * 1e-3 * c),
        xlims=(0, 50),
        xlabel=L"$\Delta c$ [mOsm]",
        ylabel="compression [%]",
        label="";
        varargs...
    )
    if !isnothing(minparams) && !isnothing(maxparams)
        plot!(
            c -> calcOsmcompression(params; Ψ=-params.RT_Pa_per_Osm * 1e-3 * c, minparams...),
            fillbetween=c -> calcOsmcompression(params; Ψ=-params.RT_Pa_per_Osm * 1e-3 * c, maxparams...),
            xlims=(0, 50),
            label="",
            color=1, alpha=0.5;
            varargs...
        )
    end
    if plotLilius
        plot!(
            x->params.Liliuscompression,
            xlims=(0, 50),
            color=:grey, label="";
            varargs...,
        )
        plot!(
            x->params.Liliuscompression - 1.0,
            fillbetween=x->params.Liliuscompression + 1.0,
            xlims=(0, 50),
            label="", color=:grey, alpha=0.5;
            varargs...,
        )
    end
end

function makepolargrid(params)
    @unpack a, b = params

    θs = 0:π/16:π/2
    rs = 1.1*a:(b-a)/5:(b - 0.1*a)
    angles = θs * ones(length(rs))'
    Rs = ones(length(θs)) * rs'
    X = Rs[:] .* cos.(angles[:])
    Y = Rs[:] .* sin.(angles[:])

    return X, Y
end

function df(x, y; params, scale=1000.0)
    r = √(x^2 + y^2)
    v = calcv(r, params) * scale
    return (v / r) * x, (v / r) * y
end

function normalize(vecs; maxlength=1)
    rawlengths = [√(x^2 + y^2) for (x, y) in vecs]
    rawmaxlength = maximum(abs.(rawlengths))
    return [maxlength .* v ./ rawmaxlength for v in vecs]
end

function addflowquivers!(params; plotvarargs...)
    X, Y = makepolargrid(params)
    vs = df.(X[:], Y[:]; params)
    vs = normalize(vs; maxlength=0.15e-3)
    quiver!(X[:], Y[:], quiver=vs, color=:black, arrow=:closed; plotvarargs...)
end