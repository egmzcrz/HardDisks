# Hard Disks Project Documentation

```@docs
updateBox!(box::Box, ϕMax::Float64; steps::Int=1000000)
updateBox!(box::Box)
computeNextCollisionInsideCell(box::Box, p::Particle, row::Int, col::Int, currTime::Float64, currPartner::Int)
computeNextCollisionInsideLNeighborhood(box::Box, p::Particle)
computeNextCollisionInsideVisibleNeighborhood(box::Box, p::Particle, dy::Int, dx::Int)
handleCollisionEvent!(box::Box, event::Event)
handleTransferEvent!(box::Box, event::Event)
```
