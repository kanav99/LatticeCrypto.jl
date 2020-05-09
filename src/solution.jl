mutable struct Solution{vecType, normType, tmpType, bType}
	minvec::vecType
	minnorm::normType
	basis::bType
	tmp::tmpType
end

function update_solution!(solution, v)
	tmp = solution.tmp
	basis = solution.basis
	tmp .= false
	for i in 1:length(v)
		@. tmp += basis[i] * v[i]
	end
	nrm = norm(tmp)
	if !iszero(nrm) && nrm < solution.minnorm
		solution.minnorm = nrm
		solution.minvec .= v
	end
end

(s::Solution)(v) = update_solution!(s, v)

function Base.show(io::IO, s::Solution)
  println(io, "coefficients: ")
  println(io, s.minvec)
  println(io, "norm: ")
  println(io, s.minnorm)
end

