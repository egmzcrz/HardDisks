using Random
include("box.jl")

function gpPrint(gp::Base.Process, X...)
    for x in X[1:end-1]
        print(gp, x, " ")
    end
    println(gp, X[end])
end

function gpPrintParticles(gp::Base.Process, box::Box)
    for p in box.config
        r = box.clock * p.growthRate
        dt = box.clock - p.clock
        m = p.mass
        rx = p.rx + p.vx * dt
        ry = p.ry + p.vy * dt
        if ry + r > 1.0 + 1e-10
            ry -= 1.0
        elseif ry - r < -1e-10
            ry += 1.0
        end
        gpPrint(gp, rx, ry, r, m)
    end
    println(gp, "e")
end

function gpPrintVelocity(gp::Base.Process, box::Box)
    for p in box.config
        dt = box.clock - p.clock
        rx = p.rx + p.vx * dt
        ry = p.ry + p.vy * dt
        r = box.clock * p.growthRate
        v = sqrt(p.vx^2 + p.vy^2)
        gpPrint(gp, rx, ry, r * p.vx / v, r * p.vy / v)
    end
    println(gp, "e")
end

function gpPrintLabels(gp::Base.Process, box::Box)
    for i in eachindex(box.config)
        p = box.config[i]
        dt = box.clock - p.clock
        rx = p.rx + p.vx * dt
        ry = p.ry + p.vy * dt
        gpPrint(gp, rx, ry, i)
    end
    println(gp, "e")
end

function boxStyle(cellWidth::Float64, gapLength::Float64)
    w = cellWidth
    L = w / 2.0
    x = 0.5 - gapLength / 2.0
    X = 0.5 + gapLength / 2.0
    return """
    set colorsequence podo
    set datafile separator ' '
    unset colorbox
    set cbrange [0:1]
    set palette rgb 23,28,-3
    set style fill solid 0.4 border 1
    set size ratio -1
    set tics 0,$w,1 scale 0 format ''
    set grid lw 1
    set xrange [-$L:1+$L]
    set yrange [-$L:1+$L]
    set obj rect from -$L,1   to 1+$L,1+$L  fc 'gray90' front fs noborder
    set obj rect from -$L,0   to 1+$L, -$L  fc 'gray90' front fs noborder
    set obj rect from 0,-$L   to -$L, 1+$L  fc 'gray90' front fs noborder
    set obj rect from 1,-$L   to 1+$L,1+$L  fc 'gray90' front fs noborder
    set obj rect from $x,1 to $X,1+$L fc '#8848A14D' front fs noborder
    set obj rect from $x,0 to $X,-$L  fc '#88B33F40' front fs noborder
    set obj rect from 0,0     to 1,1        fc 'white' front fs transparent solid 0.1
    set obj rect from -$L,-$L to 1+$L,1+$L  fc 'white' front fs transparent solid 0.1
    """
end

function vizBox(ϕ₀::Float64, n::Vector{Int}, g::Vector{Float64}, gapLength::Float64,
    frames::Int, frameSkip::Int; pause::Bool=false)
    # create box
    box = Box(n, g, gapLength)
    updateBox!(box, ϕ₀)
    println(lpad("ϕ₀ = ", 10), round(box.packingRate * box.clock^2, digits=5))
    # start visualization
    gp = open(`gnuplot`, "w")
    println(gp, boxStyle(box.cellWidth, box.gapLength))
    println(gp, "set term qt persist noraise size 700,700 pos 700,50 font 'Helvetica,8'")
    p = box.config[40]
    for k = 0:frameSkip:frames
        println(minKey(box.events), " ", p.rx, " ", p.ry)
        println(
            gp,
            """unset label
        set label 'i = $k' at 0.3,1.025 front center
        set label 't = $(box.clock)' at 0.5,1.025 front center
        set label 'ϕ = $(box.packingRate*box.clock^2)' at 0.7,1.025 front center
        """
        )

        println(
            gp,
            """plot '-' w circles lc pal lw 0.5 noti,\\
        '-' w vectors filled head lw 3 noti,\\
        '-' w labels noti"
        """
        )
        gpPrintParticles(gp, box)
        gpPrintVelocity(gp, box)
        gpPrintLabels(gp, box)

        if pause
            println(gp, "pause mouse any")
        end

        for x = 1:frameSkip
            updateBox!(box)
        end
    end
    close(gp)
end

function drawBox(ϕ::Float64, n::Vector{Int}, g::Vector{Float64}, gapLength::Float64)
    # create and update box
    box = Box(n, g, gapLength)
    updateBox!(box, ϕ)
    println(lpad("ϕ = ", 10), round(box.packingRate * box.clock^2, digits=20))
    # draw box
    gp = open(`gnuplot`, "w")
    println(gp, boxStyle(box.cellWidth, box.gapLength))
    println(gp, "set term qt persist noraise size 700,700 pos 700,50")
    println(gp, "plot '-' w circles lc pal lw 0.5 noti")
    gpPrintParticles(gp, box)
    close(gp)
end

# define seed
seed = rand(1:10^10)
Random.seed!(seed)
println(lpad("seed = ", 10), seed)
#vizBox(0.816348045, [100], [0.5], 0.0, 5000, 10, pause=false)
@time drawBox(0.85, [500, 500], [1.0, 2.0], 0.0)
