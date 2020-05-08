
# computes dot product of two vectors, works only for real valued vectors
function dot(v::Array{CType}, u::Array{CType}) where {CType<:AbstractFloat}
  s = zero(CType)
  for (a, b) in zip(v, u)
    s += a * b
  end
  return s
end

# calculates ||v||^2
function normsquare(v::Array{CType}) where {CType<:AbstractFloat}
  s = zero(CType)
  for i in v
    s += i^2
  end
  return s
end

# calculates ||v||
function norm(v::Array{CType}) where {CType<:AbstractFloat}
  return sqrt(normsquare(v))
end

# computes μ(v, u) = (v.u)/||v||^2
# used in gram-schmidt orthogonalization and LLL Algorithm
# used in calculating projection of `j` on `i`
function μ(i::Array{CType}, j::Array{CType}) where {CType<:AbstractFloat}
  return dot(i, j) / normsquare(i)
end
