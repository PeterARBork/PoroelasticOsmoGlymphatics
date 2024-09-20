using Parameters, CSV
using LinearAlgebra, QuadGK

function calcpr(r, parameters; varargs...)
    params = merge(parameters, varargs)
    @unpack a, b, Ps, Pv, Ψ, γ = params
    Pa = Pv - Ψ
    Pb = Ps - Ψ
    aterm = a * Pa * (exp.(γ * (a .- 2 * b .+ r)) - exp.(γ * (a .- r)))
    bterm = b * Pb * (exp.(γ * (2 * a .- b .- r)) - exp.(-γ * (b .- r)))
    denominator = (exp(2 * γ * (a - b)) - 1) * r
    frac = (aterm + bterm) ./ denominator

    return frac .+ Ψ
end

function calcQnumeric_cubicmeterpersec(; a, b, kappa, Ψ, η, γ, varargs...)
    #integrand(r) = 4 * π * r^2 * kappa * (calcpr(r, params) - Ψ) / η
    integrand(r) = 4 * π * r^2 * γ^2 * kappa * (calcpr(r, params) - Ψ) / η
    integral, error = quadgk(integrand, a, b)
    println("integral= $integral with error = $error")

    return integral
end

function calcQnumeric_nlpermin(; params...)
    nlpercubicmeter = 1e12
    secondspermin = 60
    return calcQnumeric_cubicmeterpersec(; params...) * nlpercubicmeter * secondspermin
end

function calcdPx(x, parameters; varargs...)
    parameters = merge(parameters, varargs)
    @unpack A, B, N = makePparams(parameters)
    @unpack γ = parameters

    term1 = (B * γ * x .+ A) .* sinh.(γ * x)
    term2 = -(A * γ * x .+ B) .* cosh.(γ * x)

    return (term1 .+ term2) ./ (N * x.^2)
end

function calcdpr(r, parameters; varargs...)
    params = merge(parameters, varargs)
    @unpack a, b, Ps, Pv, Ψ, γ = params
    Pa = Pv - Ψ
    Pb = Ps - Ψ
    aPa_terms = (γ*r .- 1) .* exp.(γ*(a .- 2*b .+ r)) .+ (γ*r .+ 1) .* exp.(γ * (a .- r))
    bPb_terms = (γ*r .+ 1) .* exp.(γ*(2*a .- b .- r)) .+ (γ*r .- 1) .* exp.(-γ * (b .- r))
    denominator = (exp(2 * γ * (a - b)) - 1) * r.^2

    return (a * Pa * aPa_terms .- b * Pb * bPb_terms) ./ denominator
end

function calcd2Px(x, parameters; varargs...)
    parameters = merge(parameters, varargs)
    @unpack A, B, N = makePparams(parameters)
    @unpack γ = parameters

    term1 = (B * γ^2 * x.^2 .+ 2 * A * γ * x .+ 2 * B) .* cosh.(γ * x)
    term2 = -(A * γ^2 * x.^2 .+ 2 * B * γ * x .+ 2 * A) .* sinh.(γ * x)

    return (term1 + term2) ./ (N * x.^3)
end

function calcd2pr(r, parameters; varargs...)
    params = merge(parameters, varargs)
    @unpack a, b, Ps, Pv, Ψ, γ = params
    Pa = Pv - Ψ
    Pb = Ps - Ψ
    gammafacminus = γ^2 * r.^2 .- 2 * γ * r .+ 2
    gammafacplus = γ^2 * r.^2 .+ 2 * γ * r .+ 2
    aPa_terms = gammafacminus .* exp.(γ*(a .- 2*b .+ r)) .- gammafacplus .* exp.(γ*(a .- r))
    bPb_terms = gammafacplus .* exp.(γ*(2*a .- b .- r)) .- gammafacminus .* exp.(-γ*(b .- r))
    denominator = (exp(2*γ*(a-b)) - 1) * r.^3

    return (a*Pa * aPa_terms .+ b*Pb * bPb_terms) ./ denominator
end

function calcv(r, parameters; varargs...)
    parameters = merge(parameters, varargs)
    @unpack η, kappa = parameters

    v = -kappa * calcdpr(r, parameters) / η
    return v
end

function calcSASv(parameters; varargs...)
    return calcv(parameters.b, parameters; varargs...)
end

function calcC1_bothfree(parameters; varargs...)
    parameters = merge(parameters, varargs)
    @unpack a, b, γ, E, ν = parameters

    M = calcM(;E, ν)
    lambda = calcLameλ(;E, ν)

    pa = calcpr(a, parameters)
    pb = calcpr(b, parameters)
    dpa = calcdpr(a, parameters)
    ddpa = calcd2pr(a, parameters)
    dpb = calcdpr(b, parameters)
    ddpb = calcd2pr(b, parameters)

    # Copied from Maple
    C1 = -(-pa * a^3 * γ^2 * M + pb * b^3 * γ^2 * M + M * ddpa * a^3 - M * ddpb * b^3 + 2 * dpa * a^2 * lambda - 2 * dpb * b^2 * lambda) / (a^3 * M - M * b^3 + 2 * a^3 * lambda - 2 * b^3 * lambda) / γ^2 / M
    return C1
end

function calcC2_bothfree(parameters; varargs...)
    parameters = merge(parameters, varargs)
    @unpack a, b, γ, E, ν = parameters

    M = calcM(;E, ν)
    lambda = calcLameλ(;E, ν)

    pa = calcpr(a, parameters)
    pb = calcpr(b, parameters)
    dpa = calcdpr(a, parameters)
    ddpa = calcd2pr(a, parameters)
    dpb = calcdpr(b, parameters)
    ddpb = calcd2pr(b, parameters)

    # copied from Maple
    C2 = -(-M * pa * a * b * γ^2 + M * pb * a * b * γ^2 + M * ddpa * a * b - M * ddpb * a * b - 2 * dpb * a * lambda + 2 * dpa * b * lambda) * b^2 * a^2 / (a^3 - b^3) / γ^2 / M / (M - lambda) / 2
    return C2
end

function calcC1_fixedSAS(parameters; varargs...)
    parameters = merge(parameters, varargs)
    @unpack a, b, γ, E, ν = parameters

    M = calcM(;E, ν)
    lambda = calcLameλ(;E, ν)
    m = lambda / M

    Pa = calcpr(a, parameters) / M
    dPa = calcdpr(a, parameters) / M
    ddPa = calcd2pr(a, parameters) / M
    dPb = calcdpr(b, parameters) / M

    cg0 = (Pa * a^3 * γ^2 - ddPa * a^3 - 2 * dPa * a^2 * m + 2 * dPb * b^2 * m - 2 * dPb * b^2) / γ^2 / (2 * a^3 * m - 2 * b^3 * m + a^3 + 2 * b^3)

    return cg0
end

function calcC2_fixedSAS(parameters; varargs...)
    parameters = merge(parameters, varargs)
    @unpack a, b, γ, E, ν = parameters

    M = calcM(;E, ν)
    lambda = calcLameλ(;E, ν)
    m = lambda / M

    Pa = calcpr(a, parameters) / M
    dPa = calcdpr(a, parameters) / M
    ddPa = calcd2pr(a, parameters) / M
    dPb = calcdpr(b, parameters) / M

    cg1 = -b^2 * a^2 * (Pa * a * b * γ^2 - ddPa * a * b + 2 * dPb * a * m - 2 * dPa * b * m + dPb * a) / γ^2 / (2 * a^3 * m - 2 * b^3 * m + a^3 + 2 * b^3)

    return cg1
end

function calcur(r, parameters; type="fixedSAS", varargs...)
    parameters = merge(parameters, varargs)
    @unpack E, ν, γ = parameters

    M = calcM(;E, ν)
    dPr = calcdpr(r, parameters) / M

    if type == "bothfree"
        C1 = calcC1_bothfree(parameters)
        C2 = calcC2_bothfree(parameters)
    elseif type == "fixedSAS"
        C1 = calcC1_fixedSAS(parameters)
        C2 = calcC2_fixedSAS(parameters)
    else
        error("type = $type must be either \"bothfree\" or \"fixedSAS\".")
    end

    return C1 * r + C2 / r^2 + γ^(-2) * dPr
end

function calcdurdr(r, parameters; type="fixedSAS", varargs...)
    parameters = merge(parameters, varargs)
    @unpack E, ν, γ = parameters

    M = calcM(;E, ν)
    d2Pr = calcd2pr(r, parameters) / M

    if type == "bothfree"
        C1 = calcC1_bothfree(parameters)
        C2 = calcC2_bothfree(parameters)
    elseif type == "fixedSAS"
        C1 = calcC1_fixedSAS(parameters)
        C2 = calcC2_fixedSAS(parameters)
    else
        error("type = $type must be either \"bothfree\" or \"fixedSAS\".")
    end

    return C1 - 2 * C2 / r^3 + γ^(-2) * d2Pr
end

function calcsigmarr(r, parameters; varargs...)
    parameters = merge(parameters, varargs)
    @unpack E, ν = parameters

    M = calcM(;E, ν)
    λ = calcLameλ(;E, ν)

    durdr = calcdurdr(r, parameters)
    ur = calcur(r, parameters)
    pr = calcpr(r, parameters)

    return M * durdr + 2 * λ * ur / r - pr
end

function calcV(parameters; varargs...)
    params = merge(parameters, varargs)
    @unpack a, b = params

    aa = a + calcur(a, params)
    bb = b + calcur(b, params)

    return (4/3) * π * (bb^3 - aa^3)
end

function calcOsmcompression(parameters; varargs...)
    params = merge(parameters, varargs)
    V0 = calcV(params; Ψ=0)
    V = calcV(params;)

    return -100 * (V - V0) / V0
end

function calcrange(outputfunc, baseparameters, paramtovary, paramdomain; varargs...)
    baseparameters = merge(baseparameters, varargs)
    paramtovary in keys(baseparameters) || throw(ArgumentError("paramtovary ($paramtovary) is not in given baseparameters."))
    results = []
    for pvalue in paramdomain
        tmpparams = merge(
            baseparameters,
            NamedTuple(zip((paramtovary, ), (pvalue, ))),
        )
        tmpresult = outputfunc(tmpparams)
        tmpresult = isa(tmpresult, NamedTuple) ? [tmpresult, ] : tmpresult
        append!(results, tmpresult)
    end

    return results
end
