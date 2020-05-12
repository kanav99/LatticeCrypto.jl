
# calculates ||v||^2
normsquare(v) = LinearAlgebra.dot(v, v)

# calculates ||v||
norm(v) = sqrt(normsquare(v))

# computes μ(v, u) = (v.u)/||v||^2
# used in gram-schmidt orthogonalization and LLL Algorithm
# used in calculating projection of `j` on `i`
@inline μ(i, j) = LinearAlgebra.dot(i, j) / normsquare(i)

function allocate_btilde(basis::Array{<:Array{CType}}) where {CType}
  n = length(basis)
  return [ similar(basis[1], float(CType)) for i in 1:n ]
end

function determinant(basis::Array{<:Array{CType}}) where {CType}
  det(hcat(basis...))
end

# BigInt

function add!(a::BigInt, b::BigInt, c::BigInt)
   ccall((:__gmpz_add, :libgmp), Nothing, (Ref{BigInt}, Ref{BigInt}, Ref{BigInt}), a, b, c)
   return a
end

function mul!(a::BigInt, b::BigInt, c::BigInt)
   ccall((:__gmpz_mul, :libgmp), Nothing, (Ref{BigInt}, Ref{BigInt}, Ref{BigInt}), a, b, c)
   return a
end

function zero!(a::BigInt)
   ccall((:__gmpz_set_si, :libgmp), Nothing, (Ref{BigInt}, Int), a, 0)
   return a
end

function addeq!(a::BigInt, b::BigInt)
   ccall((:__gmpz_add, :libgmp), Nothing, (Ref{BigInt}, Ref{BigInt}, Ref{BigInt}), a, a, b)
   return a
end

function addmul!(a::BigInt, b::BigInt, c::BigInt) # special case, no temporary required
   ccall((:__gmpz_addmul, :libgmp), Nothing, (Ref{BigInt}, Ref{BigInt}, Ref{BigInt}), a, b, c)
   return a
end

# BigFloat

function zero!(a::BigFloat)
   ccall((:mpfr_set_si, :libmpfr), Nothing,
         (Ref{BigFloat}, Int, Int32), a, 0, Base.MPFR.ROUNDING_MODE[])
   return a
end

function mul!(a::BigFloat, b::BigFloat, c::BigFloat)
   ccall((:mpfr_mul, :libmpfr), Nothing,
         (Ref{BigFloat}, Ref{BigFloat}, Ref{BigFloat}, Int32),
                 a, b, c, Base.MPFR.ROUNDING_MODE[])
   return a
end

function add!(a::BigFloat, b::BigFloat, c::BigFloat)
   ccall((:mpfr_add, :libmpfr), Nothing,
         (Ref{BigFloat}, Ref{BigFloat}, Ref{BigFloat}, Int32),
                 a, b, c, Base.MPFR.ROUNDING_MODE[])
   return a
end

function addeq!(a::BigFloat, b::BigFloat)
   ccall((:mpfr_add, :libmpfr), Nothing,
         (Ref{BigFloat}, Ref{BigFloat}, Ref{BigFloat}, Int32),
                 a, a, b, Base.MPFR.ROUNDING_MODE[])
   return a
end

function addmul!(a::BigFloat, b::BigFloat, c::BigFloat) # special case, no temporary required
   ccall((:mpfr_fma, :libmpfr), Nothing,
         (Ref{BigFloat}, Ref{BigFloat}, Ref{BigFloat}, Ref{BigFloat}, Int32),
                 a, b, c, a, Base.MPFR.ROUNDING_MODE[])
   return a
end

