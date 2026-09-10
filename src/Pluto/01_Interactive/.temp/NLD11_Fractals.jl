### A Pluto.jl notebook ###
# v0.19.5

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 489110e0-a3fa-47a5-a2c5-603af5563b5b
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", "..", ".."))
    Pkg.instantiate()
    using Plots, PlutoUI, LinearAlgebra
end

# ╔═╡ 4aec301f-3507-4627-b1f6-a885bfb7b593
gr()

# ╔═╡ 710c9fe1-6d3b-5987-b559-ed4a116038fa
md"""
This discrete map generates the Julia-set escape-time image; the code counts iterations until |z| > 2.

```math
\begin{aligned}
z_{n+1} &= z_n^2+c, \qquad z_n,c\in\mathbb{C}.
\end{aligned}
```
"""

# ╔═╡ 5df1d017-97da-4085-96a5-ea0d46f1b436
function juliaSetPixel(z0, c)
    z = z0
    niter = 255
    for i in 1:niter
        abs2(z)> 4.0 && return (i - 1)%UInt8
        z = z*z + c
    end
    return niter%UInt8
end

# ╔═╡ 9623f878-e816-49d9-a545-718ceccd8eab
function calcColumn!(pic, c, n, j)
    x = -2.0 + (j-1)*4.0/(n-1)
    for i in 1:n
        y = -2.0 + (i-1)*4.0/(n-1)
        @inbounds pic[i,j] = juliaSetPixel(x+im*y, c)
    end
    nothing
end

# ╔═╡ 7efd2dc2-ac71-4b69-aca1-22d4fefe9589
function juliaSetCalcRecSpawn!(pic, c, n, lo=1, hi=n, ntasks=16)
    if hi - lo > n/ntasks-1
        mid = (lo+hi)>>>1
        finish = Threads.@spawn juliaSetCalcRecSpawn!(pic, c, n, lo, mid, ntasks)
        juliaSetCalcRecSpawn!(pic, c, n, mid+1, hi, ntasks)
        wait(finish)
        return
    end
    for j in lo:hi
        calcColumn!(pic, c, n, j)
    end
    nothing
end

# ╔═╡ cbb00399-b821-4ab2-8c4b-a153d19b101c
function juliaSet(x, y, n=1000, method = juliaSetCalcRecSpawn!, extra...)
    c = x + y*im
    pic = Array{UInt8,2}(undef,n,n)
    method(pic, c, n, extra...)
    return pic
end

# ╔═╡ 98aab381-c40c-412a-af14-0c422be0b508
md"""
real(z) $(@bind rz Slider(-1:0.02:1.0,default=-0.74;show_value=true)) \
imag(z) $(@bind iz Slider(0.0:0.05:1.0,default=0.15;show_value=true)) \
N $(@bind N Slider(200:200:2000,default=400;show_value=true))
"""

# ╔═╡ ce001864-869d-4019-9947-8d1381585162
begin
	frac = juliaSet(rz,iz,N,juliaSetCalcRecSpawn!)
	plot(heatmap(frac, color=:Spectral),colorbar=false,size=(600,600))
end	

# ╔═╡ Cell order:
# ╠═489110e0-a3fa-47a5-a2c5-603af5563b5b
# ╠═4aec301f-3507-4627-b1f6-a885bfb7b593
# ╟─710c9fe1-6d3b-5987-b559-ed4a116038fa
# ╠═5df1d017-97da-4085-96a5-ea0d46f1b436
# ╠═9623f878-e816-49d9-a545-718ceccd8eab
# ╠═7efd2dc2-ac71-4b69-aca1-22d4fefe9589
# ╠═cbb00399-b821-4ab2-8c4b-a153d19b101c
# ╟─ce001864-869d-4019-9947-8d1381585162
# ╠═98aab381-c40c-412a-af14-0c422be0b508
