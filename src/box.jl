export Box, predictCollisionInCell, predictCollisionInNeighborhood, updateBox!, handleCollisionEvent!, handleTransferEvent!

include("event.jl")
include("particle.jl")

@inline real2grid(x::Float64, N::Int) = Int(ceil(x * N))
@inline real2grid(x::Float64, w::Float64) = Int(cld(x, w))
@inline image(x::Float64, n::Int, N::Int) = n < 1 ? x + 1 : (n > N ? x - 1 : x)
@inline pbcCell(n::Int, N::Int) = mod(n - 1, N) + 1
@inline pbcPosition(x::Float64) = x - floor(x)

mutable struct Box
    gapLength::Float64        # particles cross the gap with pbc's
    events::IndexMinPQ{Event} # priority queue of collision and transfer events
    config::Vector{Particle}  # vector of particles
    grid::Matrix{Vector{Int}} # space partitioning data structure
    gridSize::Int             # number of rows (cols) in the grid
    cellWidth::Float64        # width of each cell
    flux::Int                 # counter for the number of particles that jump the gap
    clock::Float64            # global time
    packingRate::Float64      # packing fraction: ϕ(t) = ∑nₖπrₖ² = (π∑nₖgₖ²)t² = ϕ´t²

    function Box(n::Vector{Int}, g::Vector{Float64}, gapLength::Float64)
        #
        # Define grid
        #
        N = 0    # = ∑nₖ        number of particles
        G = 0.0  # = maximum(g) max growth rate
        N´ = 0.0 # = ∑nₖgₖ²/G²  number of area-normalized particles
        ϕ´ = 0.0 # = (π∑nₖgₖ²)  packing rate
        for k in eachindex(n)
            N += n[k]
            ϕ´ += n[k] * g[k]^2
            G = max(g[k], G) # apparently more performant than an if-statement
        end
        N´ = ϕ´ / G^2
        ϕ´ *= π
        # gridSize = Int(max(1, floor(sqrt(N´)) - 1))
        gridSize = Int(floor(sqrt(N´)))
        cellWidth = 1.0 / gridSize
        grid = [Int[] for y in 1:gridSize, x in 1:gridSize] # in this notation y runs faster than x

        #
        # Populate grid
        #
        config = Vector{Particle}(undef, N)
        offset = 0 # keeps track of already created particles
        for k in eachindex(n)
            nₖ = n[k]
            gₖ = g[k]
            for i in 1+offset:nₖ+offset
                # create and insert particle into config vector
                p = Particle(gₖ / G, gₖ)
                config[i] = p
                # insert particle reference into corresponding grid cell
                y = real2grid(p.ry, gridSize)
                x = real2grid(p.rx, gridSize)
                push!(grid[y, x], i)
            end
            offset += nₖ
        end

        #
        # Predict first events
        #
        events = IndexMinPQ{Event}(N)
        box = new()
        box.config = config
        box.grid = grid
        box.gridSize = gridSize
        for x in 1:gridSize, y in 1:gridSize # in this notation y runs faster than x
            currCell = grid[y, x]    # current cell
            nCell = length(currCell) # number of particles inside current cell

            for k in 1:nCell
                i = currCell[k] # current particle index
                p = config[i] # current particle
                tCollision = Inf
                partner = 0
                border = NONE
                # check for particle collisions inside current cell
                for l in k+1:nCell
                    j = currCell[l] # partner particle index
                    q = config[j] # partner particle
                    t = predictCollisionTime(p, q)
                    if tCollision > t
                        tCollision = t
                        partner = j
                    end
                end
                # check for particle collisions in an L-shaped neighborhood
                t, j = predictCollisionInNeighborhood(box, p)
                if tCollision > t
                    tCollision = t
                    partner = j
                end
                # check for border collisions
                t, j = predictCollisionTime(p, gridSize)
                if tCollision > t
                    tCollision = t
                    border = j
                    partner = 0
                end
                tTransfer, y, x, y´, x´ = predictTransferTime(p, cellWidth)
                type = tCollision < tTransfer ? COLLISION : TRANSFER
                event = Event(tCollision, tTransfer, i, partner, border, y, x, y´, x´, type)
                insert!(events, i, event)
            end
        end

        box.gapLength = 0.0
        box.events = events
        box.cellWidth = cellWidth
        box.flux = 0
        box.clock = 0.0
        box.packingRate = ϕ´
        return box
    end
end


"""
    predictCollisionInCell(box, p, row, col, currTime, currPartner)

Computes the next collision time and partner of particle `p` against the particles inside the **cell at `(row,col)`**.

# Arguments

- `box::Box`: the box containing the particle configuration and grid
- `p::Particle`: the particle for which the prediction is computed
- `row::Int`: the row position of the cell to check against
- `col::Int`: the column position of the cell to check against
- `currTime::Float64`: the current time of an existing prediction for particle `p`
- `currPartner::Int`: the current partner of an existing prediction for particle `p`

# Returns

- `currTime::Float64`: the time of the next collision
- `currPartner::Int`: the index of the collision partner
"""
@inline function predictCollisionInCell(box::Box,
    p::Particle, row::Int, col::Int,
    currTime::Float64, currPartner::Int)

    gridSize = box.gridSize
    grid = box.grid
    config = box.config

    # (row,col)-cell might be out-of-bounds, handle PBCs
    row´ = pbcCell(row, gridSize)
    col´ = pbcCell(col, gridSize)

    # particle p might need to be translated to its image cell
    ry = p.ry
    rx = p.rx
    if row´ != row || col´ != col
        p.ry = image(ry, row, gridSize)
        p.rx = image(rx, col, gridSize)
    end

    # check against particles in (row´,col´)-cell
    for partner in grid[row´, col´]
        time = predictCollisionTime(p, config[partner])
        if currTime > time # TODO 1: AND if currTime > partner's collision time
            currTime = time
            currPartner = partner
        end
    end

    # put particle p back to its original position
    p.ry = ry
    p.rx = rx

    return currTime, currPartner
end


"""
    predictCollisionInNeighborhood(box, p)

Computes the next collision time and partner of particle `p` against the particles inside an **L-shaped neighborhood**:

    |   | X | X |
    |   | p | X |
    |   |   | X |

# Arguments

- `box::Box`: the box containing the particle configuration and grid
- `p::Particle`: the particle for which the prediction is computed

# Returns

- `time::Float64`: the time of the next collision
- `partner::Int`: the index of the collision partner
"""
function predictCollisionInNeighborhood(box::Box, p::Particle)
    col = real2grid(p.rx, box.gridSize)
    row = real2grid(p.ry, box.gridSize)

    time = Inf
    partner = 0

    col´ = col + 1
    for row´ in row-1:row+1
        time, partner = predictCollisionInCell(box, p, row´, col´, time, partner)
    end
    time, partner = predictCollisionInCell(box, p, row + 1, col, time, partner)

    return time, partner
end

"""
    predictCollisionInNeighborhood(box, p, dy, dx)

Computes the next collision time and partner of particle `p` against the particles inside a **proper neighborhood**:

    |   |   | X |
    |   | p | X |
    |   |   | X |
    
    | X | X | X |
    |   | p |   |
    |   |   |   |
    
    | X | X | X |
    |   | p | X |
    |   |   | X |
    
    | X | X | X |
    | X | X | X |
    | X | X | X |

# Arguments

- `box::Box`: the box containing the particle configuration and grid
- `p::Particle`: the particle for which the prediction is computed
- `dy::Int`: the direction (sign) of the step the particle took while changing rows
- `dx::Int`: the direction (sign) of the step the particle took while changing columns

# Returns

- `time::Float64`: the time of the next collision
- `partner::Int`: the index of the collision partner
"""
function predictCollisionInNeighborhood(box::Box, p::Particle, dy::Int, dx::Int)
    # REMEMBER that return statements inside IFs make the function allocate a lot more memory!

    col = real2grid(p.rx, box.gridSize)
    row = real2grid(p.ry, box.gridSize)

    time = Inf
    partner = 0

    if dx == 0 && dx == 0
        for row´ in row-1:row+1, col´ in col-1:col+1
            time, partner = predictCollisionInCell(box, p, row´, col´, time, partner)
        end
    elseif dx != 0 && dy != 0
        # TODO: check if this ever happens
        # (y + dy, x + dx), (y, x + dx), (y - dy, x + dx),
        # (y + dy, x), (y + dy, x - dx)
        col´ = col + dx
        for row´ in row-dy:row+dy
            time, partner = predictCollisionInCell(box, p, row´, col´, time, partner)
        end
        row´ = row + dy
        for col´ in col:col-dx
            time, partner = predictCollisionInCell(box, p, row´, col´, time, partner)
        end
    elseif dx != 0
        # (y + 1, x + dx), (y, x + dx), (y - 1, x + dx)
        col´ = col + dx
        for row´ in row-1:row+1
            time, partner = predictCollisionInCell(box, p, row´, col´, time, partner)
        end
    else
        # (y + dy, x + 1), (y + dy, x), (y + dy, x - 1)
        row´ = row + dy
        for col´ in col-1:col+1
            time, partner = predictCollisionInCell(box, p, row´, col´, time, partner)
        end
    end
    return time, partner
end


"""
    handleCollisionEvent!(box, event)

Handles a wall collision or a particle collision. The involved particles get their state updated.

# Arguments

- `box::Box`: the box containing the particle configuration, grid and event queue
- `event::Event`: the collision event to handle; it can be a wall or particle collision
"""
@inline function handleCollisionEvent!(box::Box, event::Event)
    i = event.particle
    j = event.partner # if wall collision then partner doesn't exist (j==0)
    p = box.config[i]
    t = event.timestamp
    if j > 0 # particle collision
        collide!(p, box.config[j], t)
        # update partner event to a CHECK event
        changeKey!(box.events, j, Event(t, Inf, j, i, NONE, 0, 0, 0, 0, CHECK))
    else # border collision
        border = event.border
        if border == RIGHT
            collideWithRightWall!(p, t)
        elseif border == TOP
            collideWithTopWall!(p, t)
        elseif border == LEFT
            collideWithLeftWall!(p, t)
        else
            collideWithBottomWall!(p, t)
        end
    end
    box.clock = t
end

"""
    handleTransferEvent!(box, event)

Handles a transfer event.

# Arguments

- `box::Box`: the box containing the particle configuration, grid and event queue
- `event::Event`: the collision event to handle; it can be a wall or particle collision

"""
@inline function handleTransferEvent!(box::Box, event::Event)
    i = event.particle
    p = box.config[i]
    t = event.timestamp
    move!(p, t)
    y = real2grid(p.ry, box.gridSize)
    x = real2grid(p.rx, box.gridSize)
    if y == event.y && x == event.x
        #println("particle moved but did not transfer")
        t = t + 1e-17
        move!(p, t)
        #println(event.timestamp)
        #println(p.ry, " ", p.rx)
        #println(p.ry, " ", p.rx)
    end
    # delete from old cell
    cellPrev = box.grid[event.y, event.x]
    deleteat!(cellPrev, findfirst(isequal(i), cellPrev))
    # add to new cell
    y´ = pbcCell(event.y´, box.gridSize)
    x´ = pbcCell(event.x´, box.gridSize)
    cellNext = box.grid[y´, x´]
    push!(cellNext, i)

    box.clock = t
end


"""
    updateBox!(box)

Handles the next event in line and updates the state of the event queue and the particles involved.

The algorithm follows this steps:

1. Gets the next event from the _Indexed Minimum Priority Queue_
2. Handles the event according to its type
  1. In a particle collision both particles get their velocities updated. The next event of particle `p` gets recomputed against all neighboring cells. The event of its partner gets cancelled by changing its type to a CHECK
  2. In a wall collision the particle involved gets its velocity updated. Its next event gets recomputed against all neighboring cells
  3. In a transfer event the particle involved gets its grid position updated and its event recomputed according to the next visible neighboring cells
  4. In a check event the particle gets its event recomputed against all neighboring cells

# Arguments

- `box::Box`: the box containing the particle configuration, the grid and the event queue

"""
function updateBox!(box::Box)
    # find next event in event queue
    event = minKey(box.events)
    isBoxUpdated = event.type == CHECK ? false : true

    i = event.particle
    p = box.config[i]
    dy = 0 # used to compute the new visible neighborhood
    dx = 0 # used to compute the new visible neighborhood

    # handle different types of events
    if event.type == COLLISION
        handleCollisionEvent!(box, event)
    elseif event.type == TRANSFER
        handleTransferEvent!(box, event)

        # get the new visible part of the neighborhood
        dy = sign(event.y´ - event.y)
        dx = sign(event.x´ - event.x)

        # check if transfer contributes to flux
        if event.y´ > box.gridSize
            box.flux += 1
            p.rx = pbcPosition(p.rx)
            p.ry = pbcPosition(p.ry)
        end
    else # event.type == CHECK
        # Nothing to be handled
    end

    # compute particle's next transfer time
    tTransfer, y, x, y´, x´ = predictTransferTime(p, box.cellWidth)

    # compute particle's next collision time within appropriate neighborhood
    tParticleCollision, partner = predictCollisionInNeighborhood(box, p, dy, dx)
    tBorderCollision, border = predictCollisionTime(p, box.gridSize)

    # update particle 
    tCollision = tParticleCollision
    if tParticleCollision < tBorderCollision
        border = NONE
        # TODO-2 (see TODO-1): Update partner's event and change the event type of partner's old parter to a CHECK
        #event´ = keyOf(box.events, partner)
        #partner´ = event´.partner # partner's old partner may not exist
        ## adjust event for newly found partner
        #type´ = tCollision < event´.tTransfer ? COLLISION : TRANSFER
        #event´ = Event(tCollision, event´.tTransfer,
        #    partner, i,
        #    event´.border,
        #    event´.y, event´.x, event´.y´, event´.x´,
        #    type´)
        #changeKey!(box.events, partner, event´)

        #update event type to a CHECK for old partner of newly found partner
        #if partner´ > 0
        #	event´´ = keyOf(box.events, partner´)
        #	event´´ = Event(event´´.tCollision, event´´.tTransfer,
        #									partner´, event´´.partner, event´´.border,
        #                                   event´´.y, event´´.x, event´´.y´, event´´.x´, CHECK)
        #	changeKey!(box.events, partner´, event´´)
        #end
    else
        partner = 0
        tCollision = tBorderCollision
    end


    # TODO: this check may be implemented at the moment of recomputing against new visible neighborhood
    # adjust particle event
    if event.type == TRANSFER
        # check if new partner is better than current
        if event.tCollision < tCollision
            tCollision = event.tCollision
            partner = event.partner
            border = event.border
        end
    end

    event = Event(tCollision, tTransfer,
        i, partner,
        border,
        y, x, y´, x´,
        tTransfer < tCollision ? TRANSFER : COLLISION)
    changeKey!(box.events, i, event)
    return isBoxUpdated
end

"""
    updateBox!(box, ϕMax, steps)

This is the main program loop. It calls [`updateBox!(box)`](@ref) until the packing fraction reaches an upper limit or its relative increment gets lower than a given threshold

# Arguments

- `box::Box`: the box containing the particle configuration, grid and event queue
- `ϕMax::Float64`: the target packing fraction
- `n::Int`: the number of measurements to average to get the average relative increment of the packing fraction
"""
function updateBox!(box::Box, ϕMax::Float64; n::Int=1000000)
    ϕCurr = 0.0
    δϕ = 0.0    # ∑(tCurr/tPrev)²
    Δϕ = Inf    # <Δϕ> = ∑[(ϕCurr-ϕPrev)/ϕPrev]/n = ∑[(tCurr/tPrev)²-1]/n = ∑(δϕ)/n - 1
    k = 0
    ϕRate = box.packingRate
    while ϕCurr < ϕMax && Δϕ > 1.0e-6
        tPrev = box.clock
        updateBox!(box)
        tCurr = box.clock

        ϕCurr = ϕRate * tCurr^2

        k += 1
        δϕ += (tCurr / tPrev)^2
        if mod(k, n) == 0
            Δϕ = δϕ / n - 1
            δϕ = 0.0
            k = 0
        end
    end
end
