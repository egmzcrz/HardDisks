export Event

#/*---------------------------------------------------------------------------+
# |                              Event Structure                              |
# +---------------------------------------------------------------------------*/

@enum Border::Int NONE RIGHT TOP LEFT BOTTOM
@enum EventType::Int COLLISION TRANSFER CHECK

struct Event
    timestamp::Float64  # min(tCollision,tTransfer)
    tCollision::Float64 # collision timestamp
    tTransfer::Float64  # transfer timestamp
    particle::Int       # particle index
    partner::Int        # partner index
    border::Border      # border index
    y::Int              # current y position in grid (row)
    x::Int              # current x position in grid (col)
    y´::Int             # next y position in grid (row)
    x´::Int             # next x position in grid (col)
    type::EventType     # event type
    function Event(tCollision::Float64, tTransfer::Float64, particle::Int,
            partner::Int, border::Border, y::Int, x::Int, y´::Int, x´::Int, type::EventType)
        return new(min(tCollision,tTransfer), tCollision, tTransfer,
                   particle, partner,
                   border,
                   y, x, y´, x´,
                   type)
    end
end

import Base: <
function <(a::Event, b::Event)
    return a.timestamp < b.timestamp
end

#/*---------------------------------------------------------------------------+
# |                              Priority Queue                               |
# +---------------------------------------------------------------------------*/
#/*---------------------------------------------------------------------------+
# |         Implementation based on:                                          |
# |         https://algs4.cs.princeton.edu/24pq/                              |
# |         https://algs4.cs.princeton.edu/24pq/IndexMinPQ.java.html          |
# +---------------------------------------------------------------------------*/

mutable struct IndexMinPQ{T}
    maxN::Int       # maximum number of elements in PQ
    n::Int          # number of elements in PQ
    pq::Vector{Int} # binary heap using 1-based indexing
    qp::Vector{Int} # inverse of pq: qp[pq[i]] = pq[qp[i]] = i
    keys::Vector{T} # keys[i] = priority of i
    function IndexMinPQ{T}(size::Int) where {T}
        maxN = size
        n = 0
        keys = Vector{T}(undef, maxN)
        pq = Vector{Int}(undef, maxN)
        qp = Vector{Int}(undef, maxN)
        fill!(qp, -1)
        return new(maxN, n, pq, qp, keys)
    end
end

function insert!(pq::IndexMinPQ{T}, i::Int, key::T) where {T}
    pq.n += 1
    pq.qp[i] = pq.n
    pq.pq[pq.n] = i
    pq.keys[i] = key
    swim!(pq, pq.n)
end

function minIndex(pq::IndexMinPQ{T}) where {T}
    return pq.pq[1]
end

function minKey(pq::IndexMinPQ{T}) where {T}
    return pq.keys[pq.pq[1]]
end

function keyOf(pq::IndexMinPQ{T}, i::Int) where {T}
    return pq.keys[i]
end

function changeKey!(pq::IndexMinPQ{T}, i::Int, key::T) where {T}
    pq.keys[i] = key
    swim!(pq, pq.qp[i])
    sink!(pq, pq.qp[i])
end

function increaseKey!(pq::IndexMinPQ{T}, i::Int, key::T) where {T}
    pq.keys[i] = key
    sink!(pq, pq.qp[i])
end

function greater(pq::IndexMinPQ{T}, i::Int, j::Int) where {T}
    keys = pq.keys
    return keys[pq.pq[j]] < keys[pq.pq[i]]
end

function exch!(pq::IndexMinPQ{T}, i::Int, j::Int) where {T}
    pq.pq[i], pq.pq[j] = pq.pq[j], pq.pq[i]
    pq.qp[pq.pq[i]] = i
    pq.qp[pq.pq[j]] = j
end

half(i::Int) = i >> 1
dble(i::Int) = i << 1

function swim!(pq::IndexMinPQ{T}, k::Int) where {T}
    while k > 1 && greater(pq, half(k), k)
        exch!(pq, k, half(k))
        k = half(k)
    end
end

function sink!(pq::IndexMinPQ{T}, k::Int) where {T}
    n = pq.n
    while dble(k) <= n
        j = dble(k)
        if j < n && greater(pq, j, j+1); j += 1 end
        if !greater(pq, k, j); break; end
        exch!(pq, k, j)
        k = j
    end
end
