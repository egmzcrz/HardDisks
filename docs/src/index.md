# Hard Disks Project Documentation

```@docs
predictCollisionInCell(box::Box, p::Particle, row::Int, col::Int, currTime::Float64, currPartner::Int)
predictCollisionInNeighborhood(box::Box, p::Particle)
predictCollisionInNeighborhood(box::Box, p::Particle, dy::Int, dx::Int)
handleCollision!(box::Box, event::Event)
updateBox!(box::Box)
updateBox!(box::Box, ϕMax::Float64; steps::Int=1000000)
```
