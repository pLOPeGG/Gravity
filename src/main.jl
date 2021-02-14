### A Pluto.jl notebook ###
# v0.12.20

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 7b279ba0-6ddf-11eb-241a-a38f6c0f64d1
begin
	using Pkg
	Pkg.add("PlutoUI")
	
end

# ╔═╡ 63718250-6ddf-11eb-265f-a9d6289e0bc1
begin
	using Gadfly
	using Pipe: @pipe
	using LinearAlgebra
	
	using PlutoUI
end

# ╔═╡ 3a798ba0-6de3-11eb-1545-692843f2573e
begin
	sol = (pos=Float64[0, 0, 0], mass=1.0, speed=Float64[0, 0, 0])
	earth = (pos=Float64[1, 0, 0], mass=0.01, speed=Float64[0, 9, 0])
	
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
		u12 = (obj2.pos - obj1.pos) ./ norm(obj2.pos - obj1.pos)
		[u12 .* F .* obj2.mass, .-u12 .* F .* obj1.mass]
	end
	
	gravity(sol, earth)

end

# ╔═╡ 6c756960-6e5d-11eb-1217-b5b6a5f7b90d
begin
	function kinetic_energy(obj)
		0.5 * obj.mass * norm(obj.speed) ^ 2
	end
	
	function potential_energy(obj1, obj2)
		norm.(gravity(obj1, obj2)) .* norm(obj2.pos - obj1.pos) |> sum
	end
	
end

# ╔═╡ 806743a0-6ded-11eb-0e5c-65664ca23ccc
begin
	doc"""
	Computes acceleration of each body
	
	$\vec{F} = m\vec{a}$
	
	$\vec{a} = \frac{\vec{F}}{m}$
	"""
	function acceleration(obj1, obj2)
		g = gravity(obj1, obj2)
		
		a = g ./ [o.mass for o in (obj1, obj2)]
	end
	
	acceleration(sol, earth)
end

# ╔═╡ 72ae6d50-6dee-11eb-324f-f7057c68da16
begin
	doc"""
	Apply acceleration step during dt time and moves object
	Uses Euler method
	"""
	function step!(obj, a, dt)
		obj.pos .+= dt * (obj.speed .+ a*dt/2)
		obj.speed .+= a .* dt
		
		obj
	end
end

# ╔═╡ 0b609aee-6e58-11eb-0541-078a3980f960
begin
	n = 10000
	@bind step_number Slider(1:n, show_value=true)
	
end

# ╔═╡ 28897240-6dfa-11eb-0972-fdac4b2e863c
begin
	positions = zeros(Float64, 3, 2, n)
	energy = zeros(Float64, n)
	for i in 1:n
		a = acceleration(sol, earth)
		positions[:, :, i] = vcat(sol.pos, earth.pos)
		step!.([sol, earth], a, 0.01)
		
		energy[i] = sum(kinetic_energy.((sol, earth))) + potential_energy(sol, earth)
	end
end

# ╔═╡ b76fd152-6e5b-11eb-052c-9177071b4bfe
plot(x=positions[1, :, step_number], y=positions[2, :, step_number], Geom.point, Coord.Cartesian(xmin=-5, xmax=5, ymin=-5, ymax=5))

# ╔═╡ 3d9eca90-6e5e-11eb-07f5-9925802f60c2
plot(x=1:n, y=energy, xintercept=[step_number], Geom.line, Geom.vline)

# ╔═╡ Cell order:
# ╠═7b279ba0-6ddf-11eb-241a-a38f6c0f64d1
# ╠═63718250-6ddf-11eb-265f-a9d6289e0bc1
# ╠═3a798ba0-6de3-11eb-1545-692843f2573e
# ╠═6c756960-6e5d-11eb-1217-b5b6a5f7b90d
# ╠═f8f3cb32-6de4-11eb-0a3d-f37fa4b3aed6
# ╠═806743a0-6ded-11eb-0e5c-65664ca23ccc
# ╠═72ae6d50-6dee-11eb-324f-f7057c68da16
# ╠═0b609aee-6e58-11eb-0541-078a3980f960
# ╠═28897240-6dfa-11eb-0972-fdac4b2e863c
# ╠═b76fd152-6e5b-11eb-052c-9177071b4bfe
# ╠═3d9eca90-6e5e-11eb-07f5-9925802f60c2
