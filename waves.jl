using Plots, IJulia

# first default source function:
# 	replicates an oscillator in the middle point of the space
Φ₁(u,Δt,t) = begin
	o = u
	s = size(o)[1]
	middle = Int(ceil(s/2))
	o[middle,middle] = Δt^2*sin(2*pi*t)
	return o
end

# second default source function that introduces no source
Φ₂(u,Δt,t) = u

# default initial state:
# 	sets all values to 0
Ψ₁(n) = zeros(n,n)

# first default boundary function:
# 	implements reflecting boundaries	
η₁(u,u₁,u₂,CFL) = begin
	o = u₁
	o[1,:] .= 0
	o[:,1] .= 0
	o[end,:] .= 0
	o[:,end] .= 0
	return u, o
end

# second default boundary function:
# 	implements absorbing boundaries	using Mur's absorbing boundary conditions
η₂(u,u₁,u₂,CFL) = begin
	o = u
	o[1,:] = u₁[2,:] + ((CFL-1)/(CFL+1))*(u[2,:] - u₁[1,:])
	o[:,1] = u₁[:,2] + ((CFL-1)/(CFL+1))*(u[:,2] - u₁[:,1])
	o[end,:] = u₁[end-1,:] + ((CFL-1)/(CFL+1))*(u[end-1,:] - u₁[end,:])
	o[:,end] = u₁[:,end-1] + ((CFL-1)/(CFL+1))*(u[:,end-1] - u₁[:,end])
	return o,u₁
end


# Sampler struct used to record the value at a given point in the space over the simulation
# 	s denotes the recorded values
# 	x,y denote position
struct Sampler
	x::Int64
	y::Int64
	s::Array{Float64,1}

	Sampler(x,y,s) = new(x,y,s)
	Sampler(x,y) = Sampler(x,y,[])
end

# FieldObj struct used to place objects in the environment
# NOTE: objects placed on a boundary, of size (1,x)/(x,1) or overlap other objects have undefined behaviour
# 	t enumerates objects tendicies as: { 1:reflect, 2:absorb, 3:noise }
# 	x₁,x₂,y₁,y₂ denote position and size of object
struct FieldObj
	t::Int64
	x₁::Int64
	x₂::Int64
	y₁::Int64
	y₂::Int64

	FieldObj(t,p,s) = new(t,p[1],p[1]+s[1],p[2],p[2]+s[2])
end


function _obj_helper(u,u₁,Δt,CFL,objs)
	for o in objs
		r = @view u[o.y₁-1:o.y₂+1,o.x₁-1:o.x₂+1]
		r₁ = @view u₁[o.y₁-1:o.y₂+1,o.x₁-1:o.x₂+1]
		if o.t == 1
			r[2:end-1,2:end-1] .= 0
			r₁[2:end-1,2:end-1] .= 0
		elseif o.t == 2
			r[2,2:end-1] = r₁[1,2:end-1] + ((CFL-1)/(CFL+1))*(r[1,2:end-1] - r₁[2,2:end-1])
			r[2:end-1,2] = r₁[2:end-1,1] + ((CFL-1)/(CFL+1))*(r[2:end-1,1] - r₁[2:end-1,2])
			r[end-1,2:end-1] = r₁[end,2:end-1] + ((CFL-1)/(CFL+1))*(r[end,2:end-1] - r₁[end-1,2:end-1])
			r[2:end-1,end] = r₁[2:end-1,end] + ((CFL-1)/(CFL+1))*(r[2:end-1,end] - r₁[2:end-1,end-1])
			r[3:end-2,3:end-2] .= 0
		elseif o.t == 3
			Δnoise = Δt^2*(rand(abs(o.y₂-o.y₁)+1,abs(o.x₂-o.x₁)+1) .- 0.5)
			r[2:end-1,2:end-1] += Δnoise
		end
	end
	return u, u₁
end


"""
	two_d_sim()

Runs a numerical simulation of the wave equation in 2-d by using a time discretization of the wave equation,

(∂^2 u/ ∂ t^2) = c^2 (∂^2 u/ ∂ x^2) + c^2(∂^2 u/ ∂ y^2) + f,

with the assumption that Δx = Δy.

# Examples:
julia-repl \\
julia> include("waves.jl") \\
julia> two _ d _ sim() \\
Output:  \\
	"2d-wave-3d.gif": this is an animation of the simulation with default parameters using a 3-d surface plot \\
	"2d-wave-hm.gif": this is an animation of the simulation with default parameters using a heatmap plot \\
	"2d-wave-[i].gif": this is a plot of the values of the sampler i over the simulation \\
"""
function two_d_sim(;fn::String="2d-wave", # filename to prepend output
		B::Function=η₁, # boundary function
		I::Function=Ψ₁, # initial state function
		f::Function=Φ₁, # source funtion
		Ts::Real=1.0, # number of T increments in a second 
		ns::Real=10.0, # number of seconds
		L::Real=10, # length of x and y-axis, the unit is arbitrary
		smprs::Array{Sampler,1}=Sampler[], # records input at different points of the 2-d space;
		objs::Array{FieldObj,1}=FieldObj[], # array of objects in the environment
		# the following parameters I do not recommend changing
		CFL::Real=0.5, # CFL (Courant-Friedrichs-Lewy) constant
		c::Real=1.0, # wave propagation rate 
		Δl::Real=0.1 # distance between points of mesh, unit corresponds to L
		)
	
	T = Ts*ns # determine the total number of time-steps in the simulation
	
	x = collect(0:Δl:L) # create mesh points
	n = size(x,1) # number of mesh points in the x and y-axis

	u = I(n) # u denotes the wave at time point t+Δt; I(n) creates initial state of wave
	u₁ = zeros(n,n) # u₁ denotes the wave at time point t
	u₂ = zeros(n,n) # u₂ denotes the wave at time point t-Δt

	change = zeros(n,n) # denotes the total change at a location throughout the simulation

	Δt = CFL*Δl/c # compute the timestep that will gaurentee stability of the simulation
	
	fps = ceil((T/Δt)/(60*ns)) # determine every *fps* timesteps will be taken for the animation
	l_frame = 0 # record how many timesteps since last frame
	
	# initialize animations
	anim1 = Plots.Animation()
	anim2 = Plots.Animation()
	anim3 = Plots.Animation()
	
	for t=Δt:Δt:T
		u,u₁ = B(u,u₁,u₂,CFL) # enforce the boundary conditions
		u,u₁ = _obj_helper(u,u₁,Δt,CFL,objs) # enforce conditions on objects

		step_change = abs.(u₁-u₂)
		change += step_change
		u₂ = u₁ # update data at timestep t-Δt by setting it to the previous data at t
		u₁ = u # update data at timestep t-Δt by setting it to the previous data at t+Δt
		
		u₁ = f(u,Δt,t) # add any source to the space


		# compute the data at timestep t+Δt
		u = 2*u₁ - u₂
		c = -4*u₁[2:end-1,2:end-1]
		c += u₁[3:end,2:end-1]
		c += u₁[1:end-2,2:end-1]
		c += u₁[2:end-1,1:end-2]
		c += u₁[2:end-1,3:end]
		c = CFL^2*c
		u[2:n-1,2:n-1] += c

		u,u₁ = _obj_helper(u,u₁,Δt,CFL,objs) # enforce conditions on objects
	
		# record the data at the sampling points
		old_smprs = smprs
		smprs = Sampler[]
		for s in old_smprs
			s = Sampler(s.x,s.y,[s.s; u[s.y,s.x]])
			push!(smprs,s)
		end
		
		# plot the data at timestep t
		if l_frame == 0
			@show t
			plot(x,x,u₁,st=:surface,c=cgrad(:grays),zlims=(-Δt^2,Δt^2), clims=(-Δt^2,Δt^2),title=fn*" surface plot")
			Plots.frame(anim1)
			heatmap(u₁,c=cgrad(:grays),clims=(-Δt^2,Δt^2),title=fn*" heatmap plot")
			Plots.frame(anim2)	
			heatmap(step_change,c=cgrad(:grays),clims=(0,Δt^2),title=fn*" magnitude of derivative plot")
			Plots.frame(anim3)	
		end
		l_frame = (l_frame + 1) % fps
	end
	
	# save the animation created
	gif(anim1,fn*"-3d.gif",fps=60)
	gif(anim2,fn*"-hm.gif",fps=60)
	gif(anim3,fn*"-chm.gif",fps=60)
	
	# save the data the samplers received
	smp_num = 1
	for s in smprs
		plot(s.s,ylims=(-Δt^2,Δt^2),legend=false,title=fn*" sampler at x: "*string(s.x)*", y: "*string(s.y))
		png(fn*"-s"*string(smp_num))
		smp_num += 1
	end
	heatmap(change,c=cgrad(:grays))
	png(fn*"-tc")
end


# same as above function but only computes the changed heatmap
function two_d_sim_ch(;fn::String="2d-wave", # filename to prepend output
		B::Function=η₁, # boundary function
		I::Function=Ψ₁, # initial state function
		f::Function=Φ₁, # source funtion
		Ts::Real=1.0, # number of T increments in a second 
		ns::Real=10.0, # number of seconds
		L::Real=-1, # length of x and y-axis, the unit is arbitrary
		smprs::Array{Sampler,1}=Sampler[], # records input at different points of the 2-d space;
		objs::Array{FieldObj,1}=FieldObj[], # array of objects in the environment
		# the following parameters I do not recommend changing
		CFL::Real=0.5, # CFL (Courant-Friedrichs-Lewy) constant
		c::Real=1.0, # wave propagation rate 
		Δl::Real=0.1 # distance between points of mesh, unit corresponds to L
		)
	
	if L == -1
		L = 100*Δl
	end

	T = Ts*ns # determine the total number of time-steps in the simulation
	
	x = collect(0:Δl:L) # create mesh points
	n = size(x,1) # number of mesh points in the x and y-axis

	u = I(n) # u denotes the wave at time point t+Δt; I(n) creates initial state of wave
	u₁ = zeros(n,n) # u₁ denotes the wave at time point t
	u₂ = zeros(n,n) # u₂ denotes the wave at time point t-Δt

	total_change = zeros(n,n) # denotes the total change at a location throughout the simulation
	change = zeros(n,n) # denotes the change at a location throughout the simulation

	Δt = CFL*Δl/c # compute the timestep that will gaurentee stability of the simulation
	
	for t=Δt:Δt:T

		u,u₁ = B(u,u₁,u₂,CFL) # enforce the boundary conditions
		u,u₁ = _obj_helper(u,u₁,Δt,CFL,objs) # enforce conditions on objects

		change += abs.(u₁-u₂)
		u₂ = u₁ # update data at timestep t-Δt by setting it to the previous data at t
		u₁ = u # update data at timestep t-Δt by setting it to the previous data at t+Δt
		
		u₁ = f(u,Δt,t) # add any source to the space

		# compute the data at timestep t+Δt
		u = 2*u₁ - u₂ 
		c = -4*u₁[2:n-1,2:n-1]
		c += u₁[3:n,2:n-1]
		c += u₁[1:n-2,2:n-1]
		c += u₁[2:n-1,1:n-2]
		c += u₁[2:n-1,3:n]
		c = CFL^2*c
		u[2:n-1,2:n-1] += c

		u,u₁ = _obj_helper(u,u₁,Δt,CFL,objs) # enforce conditions on objects
		
		diff_t = 16 # denotes the difference in t for when saving the change
		if t % diff_t == 0
			heatmap(change)#,c=cgrad(:grays))
			png(fn*"-tc-"*string(t))
			total_change += change
			change = zeros(n,n)
		end

		# record the data at the sampling points
		for s in smprs
			s = Sampler(s.x,s.y,[s.s; u[s.y,s.x]])
		end
	end
	
	# save the data the samplers received
	smp_num = 1
	for s in smprs
		plot(s.s,legend=false,title=fn*" sampler at x: "*string(s.x)*", y: "*string(s.y))
		png(fn*"-s"*string(smp_num))
		smp_num += 1
	end
	
	heatmap(total_change)#,c=cgrad(:grays))
	png(fn*"-total")
	
	return total_change
end


