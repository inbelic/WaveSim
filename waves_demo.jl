### total script takes about 15 mins to run

include("waves.jl")

### initial demo ###
two_d_sim()

ns = 20
objs = FieldObj[]
o1 = FieldObj(1,(65,65),(10,10))
o2 = FieldObj(2,(25,25),(10,10))
push!(objs,o1)
push!(objs,o2)

two_d_sim(fn="objects",objs=objs,ns=ns)



# microphone visualization
Φᵪ(u,Δt,t) = begin
	o = u
	s = size(o)[1]
	middle = Int(ceil(s/2))
	o[end-1,middle] = Δt^2*sin(2*pi*t)
	return o
end

smprs = Sampler[]
s = Sampler(31,31)
push!(smprs,s1)

two_d_sim(fn="rfl_mic",f=Φᵪ,smprs=smprs,ns=ns,L=5)
two_d_sim(fn="abs_mic",f=Φᵪ,smprs=smprs,B=η₂,ns=ns,L=5)



# initialize parameters for timestep to be 1/60 to allow for real-time animations
ns = 15
c = 5
Δl = 1/6
L = 10
conv_num = 1



# base case of 1hz from center with reflecting boundaries
two_d_sim(fn="bc",ns=ns,Δl=Δl,L=L,c=c)
println("bc done")



# reflecting boundaries with one absorbing edge
objs = FieldObj[]

o1 = FieldObj(2,(2,2),(56,0))
push!(objs,o1)

two_d_sim(fn="abs_edge",ns=ns,Δl=Δl,L=L,objs=objs,c=c)
println("abs_edge done")



# absorbing boundaries with one reflecting edge
objs = FieldObj[]

o1 = FieldObj(1,(2,2),(56,0))
push!(objs,o1)

two_d_sim(fn="rfl_edge",ns=ns,Δl=Δl,L=L,objs=objs,B=η₂,c=c)
println("rfl_edge done")



# clamped center and source from edge

objs = FieldObj[]

o1 = FieldObj(1,(28,28),(6,6))
push!(objs,o1)

two_d_sim(fn="clmp_ctr",ns=ns,Δl=Δl,L=L,f=Φᵪ,objs=objs,c=c)
println("clmp_ctr done")



# source placed at a non-center/edge location
Φᵪ(u,Δt,t) = begin
	o = u
	s = size(o)[1]
	middle = Int(ceil(s/2))
	o[middle+3,middle-7] = Δt^2*sin(2*pi*t)
	return o
end

two_d_sim(fn="odd_src",ns=ns,Δl=Δl,L=L,f=Φᵪ,c=c)
println("odd_src done")



# source with a limited source
Φᵪ(u,Δt,t) = begin
	o = u
	s = size(o)[1]
	middle = Int(ceil(s/2))
	if t <=3
		o[middle,middle] = Δt^2*sin(2*pi*t)
	end
	return o
end

two_d_sim(fn="lim_src",ns=ns,Δl=Δl,L=L,f=Φᵪ,c=c)
println("lim_src done")



# source with a 4hz source
Φᵪ(u,Δt,t) = begin
	o = u
	s = size(o)[1]
	middle = Int(ceil(s/2))
	o[middle,middle] = Δt^2*sin(4*2*pi*t)
	return o
end

two_d_sim(fn="4hz_src",ns=ns,Δl=Δl,L=L,f=Φᵪ,c=c)
println("4hz_src done")



# pseudo 1-d wave equation
objs = FieldObj[]

o1 = FieldObj(1,(24,2),(5,58))
o2 = FieldObj(1,(33,2),(5,58))
push!(objs,o1)
push!(objs,o2)

two_d_sim(fn="1d",ns=ns,Δl=Δl,L=L,objs=objs,B=η₂,c=c)
println("1d done")
