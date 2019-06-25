module AutoVectors
using LinearAlgebra
using DSP		# for conv = convolution

import Base.deepcopy, Base.getindex, Base.setindex!

mutable struct AutoVector{T}  			#   <: AbstractVector{T}   doesn't work
	mini::Int64
	maxi::Int64
	miniloc::Int64		# the location/index of mini-1 in dat
	def::T
	#dat::AbstractVector{T}
	dat::Vector{T}
    function AutoVector{T}(def::T,mini::Integer,maxi::Integer,miniloc::Integer,
			   				dat::Vector{T}) where T
	length(dat) < maxi-mini+1 && error("dat too small in AutoVector constructor")
	new(mini,maxi,miniloc,deepcopy(def),dat)
    end
end

AutoVector(def::T,mini::Integer,maxi::Integer,miniloc::Integer,dat::Vector{T}) where T = 
	    AutoVector{T}(def,mini,maxi,miniloc,dat)
AutoVector(def::T,mini::Integer,maxi::Integer,miniloc::Integer) where T = 
    AutoVector(def,mini,maxi,miniloc,T[deepcopy(def) for j=1:maxi-mini+1])
AutoVector(def::T) where T = AutoVector(deepcopy(def),1,0,0)

function AutoVector(f::Function,mini::Integer,maxi::Integer)
    dat = [f(i) for i=mini:maxi]
    AutoVector(zero(typeof(dat[1])),mini,maxi,0,dat)
end
AutoVector(v::Vector{T},mini::Integer,maxi::Integer) where T =
    AutoVector(zero(typeof(v[1])),mini,maxi,0,v)
AutoVector(v::Vector{T}) where T = AutoVector(zero(typeof(v[1])),1,length(v),0,v)

# Usage:  v = AutoVector(0.0)

mini(v::AutoVector) = v.mini
maxi(v::AutoVector) = v.maxi
import Base.length
length(v::AutoVector) = v.maxi-v.mini+1
arange(v::AutoVector) = mini(v):maxi(v)
olaprange(v::AutoVector,w::AutoVector) = max(mini(v),mini(w)):min(maxi(v),maxi(w))

avlocation(v::AutoVector,i) = i-v.mini+v.miniloc+1
avlocmin(v::AutoVector) = v.miniloc+1
avlocmax(v::AutoVector) = v.maxi-v.mini+v.miniloc+1
avvec(v::AutoVector) = v.dat[avlocmin(v):avlocmax(v)]

function getindex(v::AutoVector,i::Integer)
	if length(v.dat) == 0 || i < v.mini || i > v.maxi 
		return v.def
	end
	v.dat[avlocation(v,i)]
#       v.dat[i-v.mini+v.miniloc+1]
end

function fast(v::AutoVector,i)
	v.dat[i-v.mini+v.miniloc+1]
end

function clear!(v::AutoVector{T}) where {T}
	v.mini = 1
	v.maxi = 0
	v.miniloc = 0
	v.dat = T[deepcopy(v.def) for j=1:v.maxi-v.mini+1]
end

function setindex!(v::AutoVector{T},x,i::Integer) where {T}
	if(length(v.dat) == 0)
		v.dat = T[deepcopy(v.def) for j=1:19]
		v.miniloc = 10
		v.mini = v.maxi = i
	elseif(i < v.mini || i > v.maxi)
		newmini = min(i,v.mini)
		newmaxi = max(i,v.maxi)
		j = i - v.mini + v.miniloc
		if j >= 0 && j < length(v.dat)
			v.miniloc += newmini - v.mini;
		else
			newlen = (newmaxi-newmini+20)*2
			newminiloc=div(newlen-(newmaxi-newmini+1),2)
			newdat = T[deepcopy(v.def) for j=1:newlen]
			oldbegin = v.miniloc
			oldend = v.miniloc+v.maxi-v.mini
			newbegin = newminiloc - newmini + v.mini
			for k = oldbegin:oldend
				newdat[newbegin+k-oldbegin+1] = v.dat[k+1]
			end
			v.dat = newdat
			v.miniloc = newminiloc
		end
		v.mini = newmini
		v.maxi = newmaxi
	end
	v.dat[i-v.mini+v.miniloc+1] = x
end

import Base.copy
function copy(x::AutoVector)
    deepcopy(x)
end

import Base.+,Base.-
function +(v::AutoVector,w::AutoVector)
    AutoVector(i->v[i]+w[i],min(v.mini,w.mini),max(v.maxi,w.maxi))
end
function -(v::AutoVector,w::AutoVector)
    AutoVector(i->v[i]-w[i],min(v.mini,w.mini),max(v.maxi,w.maxi))
end
function +(v::AutoVector,w::Float64)
    AutoVector(i->v[i] .+ w,v.mini,v.maxi)
end
function -(v::AutoVector,w::Float64)
    AutoVector(i->v[i] .- w,v.mini,v.maxi)
end

import Base.*
function *(f::Real,v::AutoVector)
    newv = deepcopy(v)
    newv.dat *= f
    newv
end
function *(v::AutoVector,f::Real) f*v end

function *(coef::Vector{Float64},v::Vector{AutoVector{Float64}})
    n = length(v)
    mi = minimum([mini(v[i]) for i=1:n])
    ma = maximum([maxi(v[i]) for i=1:n])
    len = ma-mi+1
    res = zeros(len)
    for i=1:n
	ra = arange(v[i]) .+ (1-mi)
	res[ra] += coef[i] * avvec(v[i])
    end
    AutoVector(res,mi,ma)
end

import Base.broadcast
#Base.broadcast(::typeof(*), ...)"
#function broadcast(*,x::AutoVector,y::AutoVector)
#import Base.(.*)
#function .*(x::AutoVector,y::AutoVector)
function broadcast(::typeof(*),x::AutoVector,y::AutoVector)
    a = max(mini(x),mini(y))
    b = min(maxi(x),maxi(y))
    @inbounds xy = Float64[fast(x,i)*fast(y,i) for i=a:b]
    newv = AutoVector(xy,a,b)
    #=
    oldv = AutoVector(0.0)
    for i=a:b
	oldv[i] = x[i] * y[i]
    end
    @show vecnorm((oldv-newv).dat)
    =#
    newv
end

function avdot(x::AutoVector,y::AutoVector)
    a = max(mini(x),mini(y))
    b = min(maxi(x),maxi(y))
    a > b && return 0.0

    xoff = -x.mini + x.miniloc + 1
    yoff = -y.mini + y.miniloc + 1
    xa = a +xoff
    #@views rf = LinearAlgebra.dot(x.dat[a+xoff:b+xoff],y.dat[a+yoff:b+yoff])
    @views rf = dot(x.dat[a+xoff:b+xoff],y.dat[a+yoff:b+yoff])
#=
    res = 0.0
    for j = max(mini(x),mini(y)):min(maxi(x),maxi(y))
@inbounds	res += fast(x,j) * fast(y,j)
    end
    if abs(rf-res) > 1.0e-10
	@show (rf,res)
    end
=#
    rf
end

function avtriple(x::AutoVector,y::AutoVector,z::AutoVector)
    res = 0.0
    for j = max(mini(x),mini(y),mini(z)):min(maxi(x),maxi(y),maxi(z))
@inbounds	res += fast(x,j) * fast(y,j) * fast(z,j)
    end
    res
end

function doprint(v::AutoVector; spacing = 1)
    for i=mini(v):maxi(v)
        println(i*spacing," ",v[i])
    end
end
function doprint(FI,v::AutoVector; spacing = 1)
    for i=mini(v):maxi(v)
        println(FI,i*spacing," ",v[i])
    end
end

import LinearAlgebra.axpy!
function axpy!(y::AutoVector,a::Float64,x::AutoVector)	# y += a * x
    for k = mini(x):maxi(x)
	y[k] += a * x[k]
    end
end

function axpy!(y::AutoVector,a::Float64,x::AutoVector, cutoff::Float64)	# y += a * x
    for k = mini(x):maxi(x)
	r = a * x[k]
	abs(r) > cutoff && (y[k] += a * x[k])
    end
end

function dotrip(ud,gd,vd,ua,ub,ga,va,gb,vb,uoff,goff,voff)
    res = 0.0
		    #j+k >= va  so k >= va - j
		    #j+k <= vb  so k <= vb - j
    for j = ua:ub, k = max(ga,va-j):min(gb,vb-j)
	@inbounds res += ud[j-uoff] * gd[k-goff] * vd[j+k-voff]
    end
    res
end

#  Same as avdot(convolve(u,g),v)
function avtripconv(u::AutoVector,g::AutoVector,v::AutoVector)
    ua,ub = mini(u),maxi(u)
    ga,gb = mini(g),maxi(g)
    va,vb = mini(v),maxi(v)
    ub+gb < va && return 0.0
    ua+ga > vb && return 0.0
    dotrip(u.dat,g.dat,v.dat,ua,ub,ga,va,gb,vb,u.mini-u.miniloc-1,g.mini-g.miniloc-1,v.mini-v.miniloc-1)
end

function tripconv(u::Vector{Float64},g::Vector{Float64},v::Vector{Float64})
    k,m,n = length(u),length(g),length(v)
    res = 0.0
    for i=1:min(k,n), j=1:min(m,n-i+1)
	@inbounds res += u[i] * g[j] * v[i+j-1]
    end
    res
end

function convolvecheck(x::Vector{Float64},y::Vector{Float64})
    mx = length(x)
    my = length(y)
    res = zeros(Float64,mx+my-1)
    for j = 1:mx
@inbounds	for k = 1:my
	    res[j+k-1] += x[j] * y[k]
	end
    end
    res
end

function myconv(u::Vector{Float64},v::Vector{Float64})
    m,n = length(u),length(v)
    if m < n
        return myconv(v,u)
    end
    if n > 40
        return conv(u,v)
    end
    s = m+n-1
    res = zeros(s)
    for i=1:m
        @inbounds res[i:i+n-1] += u[i] * v
    end

    res
end

function convolve(x::AutoVector,y::AutoVector,cut=1.0e-14)		# use absolute cutoff
    res = AutoVector(x.def)
#println(maxi(x)-mini(x),"  ",maxi(y)-mini(y))
    if typeof(x.def) == Float64
	mix = mini(x)
	mx = maxi(x)
	mx-mix < 0 && return res
	#xx = Float64[fast(x,i+mix-1) for i=1:mx-mix+1]
	xx = avvec(x)
	miy = mini(y)
	my = maxi(y)
	my-miy < 0 && return res
	#yy = Float64[fast(y,i+miy-1) for i=1:my-miy+1]
	yy = avvec(y)
#	@time vres2 = convolvecheck(xx,yy)
	vres = myconv(xx,yy)
#	@show vecnorm(vres-vres2)
	mir = mix+miy
	minj = 1
	for j = 1:length(vres)
	    abs(vres[j]) > cut && (minj = j; break)
	end
	maxj = length(vres)
	for j = length(vres):-1:1
	    abs(vres[j]) > cut && (maxj = j; break)
	end
	res.dat = vres[minj:maxj]
	res.mini=minj+mir-1
	res.maxi=maxj+mir-1
	res.miniloc=0
    else
	for j = mini(x):maxi(x)
	    for k = mini(y):maxi(y)
		res[j+k] += x[j] * y[k]
	    end
	end
    end
    res
end

function makeautotake(v::Vector{Float64},offset::Integer)		# take v as the data; very fast
    AutoVector(0.0,1-offset,length(v)-offset,0,v)
end

function makeauto(v::Vector{Float64},offset::Integer)
    res = AutoVector(0.0)
    for i=1:length(v)
	res[i-offset] = v[i]
    end
    res
end
function makeauto(v::Vector{Float64},offset::Integer,cutoff::Float64)
    res = AutoVector(0.0)
    for i=1:length(v)
	vi = v[i]
	abs(vi) > cutoff && ( res[i-offset] = v[i] )
    end
    res
end

avnorm(x::AutoVector) = vecnorm(x.dat)

# This does not make a new dat array
function applyshift(x::AutoVector,offset::Integer)
    AutoVector(x.def,x.mini+offset,x.maxi+offset,x.miniloc,x.dat)
end

function symmetrize!(m::Array{Float64,2})
    n = size(m,1)
    for i=1:n
	for j=i:n
	    m[i,j] = m[j,i] = 0.5 * (m[i,j]+m[j,i])
	end
    end
end

function domin(xd,st,ma,cut)
    @inbounds for i=st:ma
	if(abs(xd[i]) > cut)
	    return i
	end
    end
    return length(xd)+1
end
function domax(xd,st,ma,cut)
    @inbounds for i=ma:-1:st
	if(abs(xd[i]) > cut)
	    return i
	end
    end
    return 0
end

function shrink!(x::AutoVector,cut)
    ami,ama = avlocmin(x),avlocmax(x)
    mi = domin(x.dat,ami,ama,cut)
    newmin = mini(x)+mi-ami
    if newmin > maxi(x) 
	clear!(x)
	return
    end
    ma = domax(x.dat,ami,ama,cut)
    x.maxi = maxi(x) + ma - ama 
    x.miniloc += newmin - x.mini
    x.mini = newmin
end
function shrinkold!(x::AutoVector,cut)
    newmin = maxi(x)+1
    for i=arange(x)
	if(abs(fast(x,i)) > cut)
	    newmin = i
	    break
	end
    end
    if newmin == maxi(x)+1 
	clear!(x)
	return
    end
    newmax = newmin-1
    for i=maxi(x):-1:newmin
	if(abs(fast(x,i)) > cut)
	    newmax = i
	    break
	end
    end
    x.maxi = newmax
    x.miniloc += newmin - x.mini
    x.mini = newmin
end

function eigsym(Marg)
    M = copy(Marg)
    symmetrize!(M)
    F = eigen(M)
    F.values,F.vectors
end

const autovector = AutoVector

export AutoVector, autovector, mini, maxi, clear!, copy, avdot, doprint, axpy!, convolve, 
		makeauto,makeautotake,applyshift,avtriple, fast, arange, symmetrize!, 
		avlocation, avlocmin,avlocmax,avvec, shrink!, avnorm, avtripconv,eigsym

end
