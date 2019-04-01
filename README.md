# NAFF.jl

A [Julia](https://julialang.org) interpretation of the NAFF algorithm for the numerical analysis of the fundamental frequencies by J. Laskar.


## Description
Function that takes an Array (Complex or Float) and runs the NAFF [J. Laskar]. Implementation of the PyNAFF code [https://github.com/nkarast/PyNAFF](https://github.com/nkarast/PyNAFF)
algorithm.

## Usage & Arguments
- `data::Array{Float64|ComplexF64}` : The data on which to perform NAFF
- `turns::Int64` : Number of turns (i.e. indices) to use from data
- `nterms::Int64`: Number of frequencies to return. The code can return up to the maximum number of frequencies found in the signal, which can be less than nterms provided.
- `skipTurns::Int64` : Number of turns to skip. Data array starts from skipTurns+1
- `window::Int64` : Power of Hann window to use.

Returns a matrix. Each colum includes 5 rows of coded elements as
- order of harmonic
- frequency
- Amplitude
- Re{Amplitude}
- Im{Amplitude}

## Example:
```julia
julia> data = sin.(collect(0:999)*2*pi*0.123) .+ sin.(collect(0:999)*2*pi*0.234);

julia> NAFF.naff(data, 300, 4, 0, 1)
5×4 Array{Float64,2}:
  1.0           2.0           3.0         4.0
  0.123001     -0.123001      0.234      -0.234
  0.5           0.5           0.5         0.5
 -0.000431821  -0.000434902  -0.0003409  -1.45175e-7
 -0.5           0.5          -0.5         0.5
```

## Contact
N. Karastathis, nkarast .at. cern .dot. ch

## References:
J. Laskar, "The chaotic motion of the solar system: A numerical estimate of the size of the chaotic zones", Icarus Vol. 88, 2, Dec.1990, p.266-291
