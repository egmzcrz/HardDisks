using Random
include("box.jl")

function ρ(r::Float64, dr::Float64, N::Int, config::Vector{Particle})
    rmin = r - dr
    rmax = r + dr
    n = 0
    for i = 1:N-1
        p = config[i]
        for j = i+1:N
            q = config[j]
            dx = q.rx - p.rx
            dy = q.ry - p.ry
            if rmin < sqrt(dx^2 + dy^2) < rmax
                n += 2
            end
        end
    end
    return n
end

function gdr(r::Vector{Float64}, dr::Float64,
        ensembleSize::Int, numPerType::Vector{Int},
        growthPerType::Vector{Float64}, ϕ::Float64)
    g::Vector{Float64} = zeros(length(r))
    N::Int = sum(numPerType)
    nbins::Int = length(r)
    pairDistance::Vector{Float64} = []
    for i in 1:ensembleSize
        box = Box(numPerType, growthPerType, 0.0)
        updateBox!(box, ϕ)
        # Measure g(r)
        for k in 1:nbins
            g[k] += ρ(r[k], dr, N, box.config)
        end
        print("\t", i, " ")
    end
    println()
    # Averaging.
    for k in eachindex(r); g[k] /= ensembleSize*N*(2*π*r[k]); end

    return g
end

#function plotgdr(ϕ::Float64, n::Vector{Int}, g::Vector{Float64}, gapLength::Float64)
function plotgdr()
    ensembleSize = 1
    N = 1000
    growthPerType = [1.0,2.0] # 1:1.4 size ratio
    ϕ = 0.85

    for n = 0:10:1000
        numPerType = [n,N-n] # 50:50 mixture ratio

        # define range r, and shell width dr
        box = Box(numPerType, growthPerType, 0.0)
        updateBox!(box, ϕ)
        println(n, "\t", box.packingRate*box.clock^2)
        σ = 2*minimum(growthPerType)*box.clock # 1000 particles: σ = 0.031917212440780114
        r = [r for r = σ/2:σ/100:σ*5]
        dr = σ/25 # shell width

        g = gdr(r, dr, ensembleSize, numPerType, growthPerType, ϕ)

        gp = open(`gnuplot -p`, "w")
        println(gp, """set terminal epslatex standalone '12pt' background 'white' color background 'white' header '\\usepackage[utf8]{inputenc} \\usepackage{amsmath}'
                set output 'tmp-$n.tex'
                set xlabel '\$r/\\sigma\$'
                set ylabel '\$g(r)\\mathrm{d}r\$'
                set xtics 1
                set yrange [0:25]
                set arrow from 0,1 to 10,1 nohead dashtype 2 lc 'gray90'
                plot '-' w l noti
                """)
        for k in eachindex(r)
            println(gp, r[k]/σ, ",", g[k])
        end
        println(gp, "e")
        close(gp)
    end
end

function plotTimeIncrement()
    box = Box([1000], [0.5], 0.0)
    gp = open(`gnuplot`, "w")
    println(gp, "set term qt persist noraise size 700,700 pos 700,50 font 'Helvetica,8'")
    println(gp, "set logscale xy")
    println(gp, "plot '-' w l")
    i = 0
    while i < 10000000
        if mod(i, 100) == 0
            println(gp, i, ",", box.clock)
        end
        if updateBox!(box); i+=1; end
    end
    println(gp, "e")
    close(gp)
    println(lpad("ϕ = ", 10), round(box.packingRate*box.clock^2, digits=5))
end

seed = rand(1:10^10)
Random.seed!(seed)
println(lpad("seed = ", 10), seed)

@time plotgdr()
#plotTimeIncrement()
