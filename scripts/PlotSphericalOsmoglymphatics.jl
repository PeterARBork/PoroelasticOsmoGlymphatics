
using Plots, LaTeXStrings
Plots.default(fontfamily = ("computer modern"), guidefontsize=9, legend_foreground_color=:transparent)
include("../src/Osmoglymphatics.jl")
include("../src/SphericalOsmoglymphatics.jl")
include("../src/plotfunctions.jl")

params = loadparameters("data/baselineparameters50mOsm_meters_Pascal.csv")

mOsm = convertPsiPatomOsm(params.Ψ)

calcur(params.a, params)
calcur(params.b, params)
critur0 = calcur(params.b, params; type="fixedSAS")
round(critur0, digits=10) == 0.0 || error("brain not fixed at SAS")

calcur(params.a, params; Ψ=-0.05*params.RT_Pa_per_Osm, type="bothfree")
calcur(params.b, params; Ψ=-0.05*params.RT_Pa_per_Osm, type="bothfree")

physiological_Δc_mM = 3
v_b_at_1mOsm_mm_per_h = 1e3 * 60^2 * calcSASv(params; Ψ = params.RT_Pa_per_Osm * 1e-3 * physiological_Δc_mM)
v_b_at_50mOsm_mm_per_h = 1e3 * 60^2 * calcSASv(params;)
compression_in_percent = calcOsmcompression(params)
calcQnumeric_nlpermin(; params...)

calcOsmcompression(params)

calcsigmarr(params.b, params)

p = plot(); plotannularparamdeformation(params, merge(params, (Ψ=-0,)))
savefig("img/shell deformation.pdf")

begin
    lowparams = (E=params.E * 2, γ=params.γ * 0.5)
    highparams = (E=params.E / 2, γ=params.γ * 1.5)
    l = @layout [
        grid(3, 1){0.33w} [a{0.68h}
                        grid(1, 2)]
    ]
    p = plot(layout=l, size=(600, 500))
    plotannularparamdeformation(params, merge(params, (Ψ=-0,)); subplot=4, title="")
    addflowquivers!(params; subplot=4)

    plotcompressionversusc!(params; plotLilius=false, minparams=lowparams, maxparams=highparams, xticks=[0, 20, 40], yticks=[0, 4, 8], subplot=5)
    plotvbversusc!(params; plotPlogMestre=false, minparams=lowparams, maxparams=highparams, xticks=[0, 20, 40], yticks=[0, 1, 2], subplot=6)
    addurplot!(params, calcpr, "[mmHg]", "p(r)"; yfac=params.mmHg_per_Pa, minparams=lowparams, maxparams=highparams, subplot=1, color=1, ylims=(-4, 0), xformatter=_->"", yticks=[0, -2, -4])
    addurplot!(params, calcv, "[mm/h]", "v(r)"; minparams=lowparams, maxparams=highparams, subplot=2, color=1, ylims=(-4, 4), yticks=[-2, 0, 2], yfac=1e3*60^2, xformatter=_->"")
    addurplot!(params, calcur, "[mm]", "u_r(r)"; minparams=lowparams, maxparams=highparams, yfac=1e3, color=1, subplot=3, xlabel=L"$r$ [mm]", yticks=[0, 0.1, 0.2, 0.3], ylims=(-0.01, 0.2), xformatter=x->1e3*x)

    cs_mM = range(0, 5, length=20)
    vs = abs.(60^2 * 1e3 * [calcSASv(params; Ψ=c * 1e-3 * params.RT_Pa_per_Osm) for c in cs_mM])
    plot!(
        cs_mM, vs,
        xlabel="", ylabel="",
        xticks=[2, 4, ],
        yticks=[0, 0.1],
        label="",
        subplot=7,
        inset =(6, bbox(0.25, -0.1, 0.5, 0.35, :top, :left));
    )

    annotate!((-0.3, 1.05), "a", subplot=1)
    annotate!((-0.3, 1.05), "b", subplot=2)
    annotate!((-0.3, 1.05), "c", subplot=3)
    annotate!((-0.2, 1.), "d", subplot=4)
    annotate!((-0.25, 1.05), "e", subplot=5)
    annotate!((-0.2, 1.05), "f", subplot=6)
end
savefig("img/Figure 2.pdf")

begin
    filtrationparams = merge(params, (Ψ = 0.004 * params.RT_Pa_per_Osm, ))
    lowparams = (E=params.E * 2, γ=params.γ * 0.5)
    highparams = (E=params.E / 2, γ=params.γ * 1.5)
    l = @layout [
        grid(3, 1){0.33w} [a{0.68h}
                        grid(1, 2)]
    ]
    p = plot(layout=l, size=(600, 500))
    plotannularparamdeformation(filtrationparams, merge(filtrationparams, (Ψ=-0,)); subplot=4, title="")
    addflowquivers!(filtrationparams; subplot=4)
    plotcompressionversusc!(filtrationparams; minparams=lowparams, maxparams=highparams, plotLilius=false, xlims=(0, 10), xticks=[0, 4, 8], yticks=[0, 1, 2], subplot=5)
    plotvbversusc!(filtrationparams; minparams=lowparams, maxparams=highparams, plotPlogMestre=false, xlims=(0, 10), xticks=[0, 4, 8], yticks=[0, 0.2, 0.4], subplot=6)
    addurplot!(filtrationparams, calcpr, "[mmHg]", "p(r)"; minparams=lowparams, maxparams=highparams,yfac=params.mmHg_per_Pa, subplot=1, color=1, ylims=(0, 0.5), yticks=[0, 0.2, 0.4], xformatter=_->"",)
    addurplot!(filtrationparams, calcv, "[mm/h]", "v(r)"; minparams=lowparams, maxparams=highparams,subplot=2, color=1, yfac=1e3*60^2, ylims=(-0.3, 0.3), yticks=[-0.2, 0, 0.2], xformatter=_->"")
    addurplot!(filtrationparams, calcur, "[μm]", "u_r(r)"; minparams=lowparams, maxparams=highparams, yfac=1e6, color=1, subplot=3, xlabel=L"$r$ [mm]", ylims=(-10, 10), yticks=[-5, 0, 5], xformatter=x->1e3*x)

    annotate!((-0.3, 1.05), "a", subplot=1)
    annotate!((-0.3, 1.05), "b", subplot=2)
    annotate!((-0.3, 1.05), "c", subplot=3)
    annotate!((-0.2, 1.), "d", subplot=4)
    annotate!((-0.2, 1.05), "e", subplot=5)
    annotate!((-0.2, 1.05), "f", subplot=6)
end
savefig("img/filtrationfig.pdf")
