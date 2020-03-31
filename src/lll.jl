mutable struct LLLCache{btType,tmpType}
  btilde::btType
  tmp::tmpType
end

function generate_btilde!(btilde::Array{<:Array{CType}}, basis::Array{<:Array{CType}}) where {CType<:AbstractFloat}
  len = length(basis)
  for j in 2:len
    btilde[j] .= basis[j]
    for i in 1:j-1
      btilde[j] .= btilde[j] .- btilde[i] .* μ(btilde[i], basis[j])
    end
  end
end

function size_reduce!(basis::Array{<:Array{CType}}, btilde::Array{<:Array{CType}}) where {CType<:AbstractFloat}
  len = length(basis)
  for j in 2:len
    for i in j-1:-1:1
      basis[j] .= basis[j] .- round(μ(btilde[i], basis[j])) .* basis[i]
    end
  end
end

function lovasz(basis::Array{<:Array{CType}}, btilde::Array{<:Array{CType}}, tmp::Array{CType}) where {CType<:AbstractFloat}
  len = length(basis)
  for i in 1:len-1
    tmp .= μ(btilde[i], basis[i+1]) .* btilde[i] .+ btilde[i+1]
    if 0.75normsquare(btilde[i]) > normsquare(tmp)
      return i
    end
  end
  return 0
end

function lll!(basis::Array{<:Array{CType}}, cache) where {CType<:AbstractFloat}
  btilde = cache.btilde
  while true
    generate_btilde!(btilde, basis)
    size_reduce!(basis, btilde)
    i = lovasz(basis, btilde, cache.tmp)
    if i == 0
      break
    else
      tmp = basis[i]
      basis[i] = basis[i+1]
      basis[i+1] = tmp 
    end
  end
  return basis
end

function run_lll!(basis)
  cache = LLLCache(deepcopy(basis), similar(basis[1]))
  lll!(basis, cache)
end
