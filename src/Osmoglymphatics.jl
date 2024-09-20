using Parameters, CSV

function loadparameters(filepath)
    return NamedTuple((Symbol(row.name), row.value) for row in CSV.File(filepath))
end

function convertPsimmHgtomOsm(ΨmmHg; RT_mmHg_per_mOsm=18.7)
    return abs(ΨmmHg / RT_mmHg_per_mOsm)
end

function convertPsiPatomOsm(ΨPa; RT_Pa_per_Osm=2.5e6)
    # The gas constant is 8314 Pa L / K / mol = 8314 Pa / K / Osm
    # With temperature T = 300, RT = 2.5e6.
    # The corresponding mOsm will be factor 1e3 greater than the Osm.
    return 1e3 * abs(ΨPa / RT_Pa_per_Osm)
end

function convertfromLametoYoungs(;λ, μ)
    E = μ * (3 * λ + 2 * μ) / (λ + μ)
    ν = λ / (2 * (λ + μ))

    return (E=E, ν=ν)
end

function calcM(; E, ν)
    return E * (1 - ν) / ((1+ν) * (1-2*ν))
end

function calcLameλ(; E, ν)
    return E * ν / ((1 + ν) * (1 - 2 * ν))
end

function calcLameμ(; E, ν)
    return E / (2 * (1 + ν))
end
