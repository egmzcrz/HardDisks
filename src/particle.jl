export Particle

mutable struct Particle
    rx::Float64
    ry::Float64
    vx::Float64
    vy::Float64
    mass::Float64
    growthRate::Float64
    clock::Float64

    function Particle(mass::Float64, growthRate::Float64)
        rx = rand()
        ry = rand()
        vx = -1.0 + 2.0 * rand()
        vy = -1.0 + 2.0 * rand()
        norm = sqrt(vx^2 + vy^2)

        p = new()
        p.rx = rx
        p.ry = ry
        p.vx = vx / norm
        p.vy = vy / norm
        p.mass = mass
        p.growthRate = growthRate # radius growth rate
        p.clock = 0
        return p
    end
end

function move!(p::Particle, t::Float64)
    # Move in a straight line.
    dt = t - p.clock
    p.rx += p.vx * dt
    p.ry += p.vy * dt

    # Update particle's internal clock
    p.clock = t
end

function collide!(p::Particle, q::Particle)
    dx::Float64 = q.rx - p.rx
    dy::Float64 = q.ry - p.ry
    dvx::Float64 = q.vx - p.vx
    dvy::Float64 = q.vy - p.vy
    G::Float64 = p.growthRate + q.growthRate

    drdr::Float64 = dx^2 + dy^2
    dr::Float64 = sqrt(drdr)
    drdv::Float64 = dx * dvx + dy * dvy

    M::Float64 = p.mass + q.mass

    # Velocities update rules:
    # Δpᵢ = mᵢmⱼ/M × 2(dv⋅n - G)n
    # vᵢ' = vᵢ + Δpᵢ/mᵢ
    # vⱼ' = vⱼ - Δpᵢ/mⱼ
    c::Float64 = 2 * (drdv / drdr - G / dr)
    cp::Float64 = q.mass / M * c
    cq::Float64 = p.mass / M * c

    p.vx += cp * dx
    p.vy += cp * dy
    q.vx -= cq * dx
    q.vy -= cq * dy
end

# TODO: weird behaviour when the collision is a back-collision
function reflectLeft!(p::Particle)
    p.vx = -(abs(p.vx) + p.growthRate)
end
function reflectRight!(p::Particle)
    p.vx = abs(p.vx) + p.growthRate
end
function reflectDown!(p::Particle)
    p.vy = -(abs(p.vy) + p.growthRate)
end
function reflectUp!(p::Particle)
    p.vy = abs(p.vy) + p.growthRate
end

function predictCollisionTime(p::Particle, q::Particle)
    if p === q
        return Inf
    end

    pClock = p.clock
    qClock = q.clock
    t = max(pClock, qClock)
    dt_p = t - pClock
    dt_q = t - qClock
    rx_p = p.rx
    ry_p = p.ry
    rx_q = q.rx
    ry_q = q.ry
    vx_p = p.vx
    vy_p = p.vy
    vx_q = q.vx
    vy_q = q.vy

    dx::Float64 = (rx_q + vx_q * dt_q) - (rx_p + vx_p * dt_p)
    dy::Float64 = (ry_q + vy_q * dt_q) - (ry_p + vy_p * dt_p)
    dvx::Float64 = vx_q - vx_p
    dvy::Float64 = vy_q - vy_p
    G2::Float64 = (p.growthRate + q.growthRate)^2

    drdr::Float64 = dx^2 + dy^2
    dvdv::Float64 = dvx^2 + dvy^2
    drdv::Float64 = dx * dvx + dy * dvy

    # en la ecuación 4.12, G es el crecimiento del diámetro!
    # aquí G es el crecimiento del radio
    a::Float64 = dvdv - G2
    b::Float64 = drdv - G2 * t
    c::Float64 = drdr - G2 * t^2

    Δ::Float64 = b^2 - a * c
    if c < 1e-12
        c = 0
    end # particles touch or overlap
    if (b <= 0 || a < 0) && Δ >= 0
        return t + max(-(b + sqrt(Δ)) / a, 0)
    else
        return Inf
    end
end

function predictTransferTime(p::Particle, w::Float64)
    ϵ::Float64 = 1e-13

    rx = p.rx
    ry = p.ry
    vx = p.vx
    vy = p.vy
    x = Int(cld(rx, w))
    y = Int(cld(ry, w))
    xmin = (x - 1) * w - ϵ
    xmax = x * w + ϵ
    ymin = (y - 1) * w - ϵ
    ymax = y * w + ϵ
    dtx = vx > 0.0 ? (xmax - rx) / vx : (xmin - rx) / vx
    dty = vy > 0.0 ? (ymax - ry) / vy : (ymin - ry) / vy
    dt = min(dtx, dty)
    x´ = Int(cld(rx + vx * dt, w))
    y´ = Int(cld(ry + vy * dt, w))
    # TODO: check if it happens that x´==x and y´==y
    return p.clock + dt, y, x, y´, x´
end
