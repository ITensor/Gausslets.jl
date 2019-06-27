using Gausslets,
      Test

g = gausslet(cf652)
x = 0.234828785

g10 = gausslet(cf1092)

g1010set = (g1010,g1010mid,g1010high) = gaussletFamilyMH(10,10)

g101010set = (g1,g2l,g2r,g3l,g3r) = gaussletFamilyMH(10,10,10)

@testset "Gausslets" begin
  olap(f1,f2) = 0.001 * sum(j->f1(j*0.001)*f2(j*0.001),-20000:20000)

  n = length(g1010set)
  for i1 in 1:n, i2 in 1:n
    f1 = g1010set[i1]
    f2 = g1010set[i2]
    if i1 == i2
      @test olap(f1,f2) ≈ 1
    else
      @test olap(f1,f2) ≈ 0 atol = 1e-15
    end
  end

  n = length(g101010set)
  for i1 in 1:n, i2 in 1:n
    f1 = g101010set[i1]
    f2 = g101010set[i2]
    if i1 == i2
      @test olap(f1,f2) ≈ 1
    else
      @test olap(f1,f2) ≈ 0 atol = 1e-15
    end
  end
end
