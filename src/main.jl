### A Pluto.jl notebook ###
# v0.12.4

using Markdown
using InteractiveUtils

# ╔═╡ 7b279ba0-6ddf-11eb-241a-a38f6c0f64d1
begin
	using Pkg
	Pkg.add("LinearAlgebra")
	
end

# ╔═╡ 63718250-6ddf-11eb-265f-a9d6289e0bc1
begin
	using Gadfly
	using Pipe: @pipe
	using LinearAlgebra
end

# ╔═╡ 3a798ba0-6de3-11eb-1545-692843f2573e
begin
	sol = (pos=Float64[0, 0, 0], mass=1.0, speed=Float64[0, 0, 0])
	earth = (pos=Float64[1, 0, 0], mass=0.01, speed=Float64[0, 1, 0])
	
	G = 1.
end

# ╔═╡ f8f3cb32-6de4-11eb-0a3d-f37fa4b3aed6
begin
	doc"""
	Computes gravity applied to each object
	
	$\vec{F_{g12}} = G \cdot \frac{m_2}{\Vert\vec{x_2}-\vec{x_1}\Vert^2} \cdot \frac{\vec{x_2}-\vec{x_1}}{\Vert\vec{x_2}-\vec{x_1}\Vert}$ 
	"""
	function gravity(obj1, obj2)
		F = G / sum((obj2.pos - obj1.pos) .^ 2)
		u⃗12 = (obj2.pos - obj1.pos) ./ norm(obj2.pos - obj1.pos)
		(u⃗12 .* F .* obj2.mass, .-u⃗12 .* F .* obj1.mass)
	end
	
	gravity(sol, earth)

end

# ╔═╡ 806743a0-6ded-11eb-0e5c-65664ca23ccc
begin
	doc"""
	Computes acceleration of each body
	
	$\vec{F} = m\vec{a}$
	
	$\vec{a} = \frac{\vec{F}}{m}$
	"""
	function acceleration(obj1, obj2)
		g⃗ = gravity(obj1, obj2)
		
		a⃗ = g⃗ ./ [o.mass for o in (obj1, obj2)]
	end
	
	acceleration(sol, earth)
end

# ╔═╡ 72ae6d50-6dee-11eb-324f-f7057c68da16
begin
	doc"""
	Apply acceleration step during dt time and moves object
	Uses Euler method
	"""
	function step!(obj, a⃗, dt)
		obj.pos .+= dt * (obj.speed .+ a⃗*dt/2)
		obj.speed .+= a⃗ .* dt
		
		obj
	end
end

# ╔═╡ 28897240-6dfa-11eb-0972-fdac4b2e863c
begin
	a⃗ = acceleration(sol, earth)
	step!.([sol, earth], a⃗, 0.01)
end

# ╔═╡ Cell order:
# ╠═7b279ba0-6ddf-11eb-241a-a38f6c0f64d1
# ╠═63718250-6ddf-11eb-265f-a9d6289e0bc1
# ╠═3a798ba0-6de3-11eb-1545-692843f2573e
# ╠═f8f3cb32-6de4-11eb-0a3d-f37fa4b3aed6
# ╠═806743a0-6ded-11eb-0e5c-65664ca23ccc
# ╠═72ae6d50-6dee-11eb-324f-f7057c68da16
# ╠═28897240-6dfa-11eb-0972-fdac4b2e863c
