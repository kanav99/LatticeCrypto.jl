
function dot(v::Array{CType}, u::Array{CType}) where {CType<:AbstractFloat}
  s = zero(CType)
  for (a, b) in zip(v, u)
    s += a * b
  end
  return s
end

function normsquare(v::Array{CType}) where {CType<:AbstractFloat}
  s = zero(CType)
  for i in v
    s += i^2
  end
  return s
end

function norm(v::Array{CType}) where {CType<:AbstractFloat}
  return sqrt(normsquare(v))
end

function μ(i::Array{CType}, j::Array{CType}) where {CType<:AbstractFloat}
  return dot(i, j) / normsquare(i)
end
