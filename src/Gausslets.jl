module Gausslets

using Reexport

include("Interp.jl")
@reexport using .Interp
@reexport using AutoVectors

include("gausslet.jl")

export gausslet, fg, transc,
	c220, c212, c322, c432,  c652, c872, c1092, 
	c212left, c220left, c432left, c652left, c872left, c1092left, cf432, cf652,
	cf872, cf1092,  c652right, c1092right,  c872right, c1092mid,c1092high,
	c652mid,c652high,c432mid,c432high,c432left,c432right,
	gaussletFamilyLR, gaussletFamilyMH

end

