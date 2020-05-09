"""
	enumerate_lattice(basis, R)

Generate neccesary data and start enumeration DFS for basis `basis` which visits all the nodes
with norm less than `R`
"""
function enumerate_lattice(basis, R; lll=true)
	# number of basis vectors
	n = length(basis)
	lll && run_lll!(basis)
	# calculate bstar
	bstar = allocate_btilde(basis)
	generate_btilde!(bstar, basis)
	# calculate norms and normsquare of bstar
	bstarnormsquare = [ normsquare(v) for v in bstar ]
	bstarnorm = [ norm(v) for v in bstar ]
	# we don't want to calculate μij all the time, lets precalculate them
	μarray = zeros(eltype(bstar[1]), n, n)
	for i in 2:n
		for j in 1:(i-1)
			μarray[i, j] = μ(bstar[j], basis[i])
		end
	end
	# start leaf vector with all zeros
	v = [ 0 for _ in 1:n ]
	# callback, what to do when a leaf node is encountered
	solution = Solution(zeros(eltype(basis[1]), n), Inf, basis, similar(bstar[1]))
	# start DFS with depth 1
	@time enumeration_tree_bfs_step(1, v, n, μarray, bstarnormsquare, bstarnorm, R, solution)
	solution
end

"""
	enumeration_tree_bfs_step(depth, v, n, μ, bstarnormsquare, bstarnorm, R, callback)

Run enumeration DFS at depth `depth`. Calculates interval of possible values of `v[n + 1 - depth]`
using `get_interval` and performs DFS starting from `center`, increasing its distance.

# Arguments:
- `depth`: depth where interval is calculated
- `v`: vector with correct values from index `n + 2 - depth` till `n`
- `n`: number of  of lattice
- `μ`: `μ` matrix for the basis
- `bstarnormsquare`: square of norms of all the vectors of B* matrix
- `bstarnorm`: norms of all the vectors of B* matrix
- `R`: upper limit of size of vector
- `callback`: if a leaf is encountered, run the `callback` function
"""
function enumeration_tree_bfs_step(depth, v, n, μ, bstarnormsquare, bstarnorm, R, callback)
	if depth == n + 1
		callback(v)
		return
	end

	# get the interval of values of `v[n + 1 - depth]`
	center, radius = get_interval(depth, v, n, μ, bstarnormsquare, bstarnorm, R)

	# enumerate at center first
	real_center = eltype(v)(round(center))
	v[n - depth + 1] = real_center
	enumeration_tree_bfs_step(depth + 1, v, n, μ, bstarnormsquare, bstarnorm, R, callback)
	# now enumerate in order `center + 1`, `center - 1`, `center + 2`, `center - 2`....
	i = 1
	while i <= radius - 1
		v[n - depth + 1] = real_center + i
		enumeration_tree_bfs_step(depth + 1, v, n, μ, bstarnormsquare, bstarnorm, R, callback)
		v[n - depth + 1] = real_center - i
		enumeration_tree_bfs_step(depth + 1, v, n, μ, bstarnormsquare, bstarnorm, R, callback)
		i += 1
	end

	# ends of radius need to be cautiously handled 
	if real_center + i <= center + radius
		v[n - depth + 1] = real_center + i
		enumeration_tree_bfs_step(depth + 1, v, n, μ, bstarnormsquare, bstarnorm, R, callback)
	end

	if real_center - i >= center - radius
		v[n - depth + 1] = real_center - i
		enumeration_tree_bfs_step(depth + 1, v, n, μ, bstarnormsquare, bstarnorm, R, callback)
	end

	i += 1
	if real_center + i <= center + radius
		v[n - depth + 1] = real_center + i
		enumeration_tree_bfs_step(depth + 1, v, n, μ, bstarnormsquare, bstarnorm, R, callback)
	end

	if real_center - i >= center - radius
		v[n - depth + 1] = real_center - i
		enumeration_tree_bfs_step(depth + 1, v, n, μ, bstarnormsquare, bstarnorm, R, callback)
	end
	nothing
end

"""
	get_interval(k, v, μ, bstarnormsquare, bstarnorm, R)

Calculates interval of possible values of `v[n + 1 - k]`. Returns a tuple `(center, radius)`
which means that `| v[n + 1 - k] - center | <= radius `

# Arguments:
- `k`: depth where interval is calculated
- `v`: vector with correct values from index `n + 2 - k` till `n`
- `n`: dimensions of lattice
- `μ`: `μ` matrix for the basis
- `bstarnormsquare`: square of norms of all the vectors of B* matrix
- `bstarnorm`: norms of all the vectors of B* matrix
- `R`: upper limit of size of vector
"""
function get_interval(k, v, n, μ, bstarnormsquare, bstarnorm, R)
	# calculating radius of interval
	if k == 1
		center = R / (2 * bstarnorm[n])
		return (center, center)
	end

	rhs = R * R
	for j in (n + 2 - k):n
		s = v[j]
		for i in (j + 1):n
			s += μ[i, j] * v[i]
		end
		s = s * s * bstarnormsquare[j]
		rhs = rhs - s
	end
	j = n + 1 - k
	radius = sqrt(rhs) / bstarnorm[j]
	# calculating center of interval
	center = zero(eltype(v))
	for i in (j + 1):n
		center = center - μ[i, j] * v[i]
	end
	return (center, radius)
end
