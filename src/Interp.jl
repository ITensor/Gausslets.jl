module Interp

using AutoVectors

# 8th order polynomial fit to evenly spaced grid points, with constraining deriv continuity
function interpgetcoefs(f) 
    Float64[(-5*f[1]+49*f[2]-245*f[3]+1225*f[4]+1225*f[5]-245*f[6]+49*f[7]-5*f[8])/2048.,
	(75*f[1]-1029*f[2]+8575*f[3]-128625*f[4]+128625*f[5]-8575*f[6]+1029*f[7]-75*f[8])/107520.,
	(259*f[1]-2495*f[2]+11691*f[3]-9455*f[4]-9455*f[5]+11691*f[6]-2495*f[7]+259*f[8])/23040.,
	(-37*f[1]+499*f[2]-3897*f[3]+9455*f[4]-9455*f[5]+3897*f[6]-499*f[7]+37*f[8])/11520.,
	(-7*f[1]+59*f[2]-135*f[3]+83*f[4]+83*f[5]-135*f[6]+59*f[7]-7*f[8])/1152.,
	(5*f[1]-59*f[2]+225*f[3]-415*f[4]+415*f[5]-225*f[6]+59*f[7]-5*f[8])/2880.,
	(f[1]-5*f[2]+9*f[3]-5*f[4]-5*f[5]+9*f[6]-5*f[7]+f[8])/1440.,
	(-f[1]+7*f[2]-21*f[3]+35*f[4]-35*f[5]+21*f[6]-7*f[7]+f[8])/5040.]
end

mutable struct Interpolator
    coefset::AutoVector{Vector{Float64}}
    ispacing::Float64
    function Interpolator(f::AutoVector{Float64},spacing)
	cs = AutoVector(zeros(8),0,0,0)
	for i=arange(f)
	    v = [f[i+j] for j=-3:4]
	    cs[i] = interpgetcoefs(v)
	end
	new(cs,1/spacing)
    end
end

function Interpolator() 
    f = AutoVector(0.0)
    f[0] = 0.0
    f[1] = 0.0
    Interpolator(f,1.0)
end

function poly(xi::Float64,co::Vector{Float64})
    @inbounds co[1]+xi*(co[2]+xi*(co[3]+xi*(co[4]+xi*(co[5]+xi*(co[6]+xi*(co[7]+xi*co[8]))))))
end
function (p::Interpolator)(x::Real)
    xx = x * p.ispacing
    if xx < mini(p.coefset) || xx > maxi(p.coefset)
	return 0.0
    end
    i = floor(Int64,xx)
    xi = xx-i-0.5
    ind = avlocation(p.coefset,i)
    #co = fast(p.coefset,i)
    #poly(xi,co)
    poly(xi,p.coefset.dat[ind])
end

export Interpolator

end

