using LinearAlgebra
using Arpack
using SparseArrays
using Combinatorics

# A julia file that computes the upper bound on the connective constant at level k on the square lattice.

"""
rotate_point(x, y, n)

Rotate a point (x, y) by n * 90 degrees (pi/2 radians) counterclockwise.

Arguments
x: The x-coordinate of the point.
y: The y-coordinate of the point.
n: The number of 90-degree rotations to apply.

Returns
A tuple (x_rot, y_rot) representing the rotated point as integers.
"""
function rotate_point(x, y, n)
    return (Int(round(x*cos(n*pi/2) - y*sin(n*pi/2))),
            Int(round(x*sin(n*pi/2) + y*cos(n*pi/2))))
end


vertex_1 = Set([Set([(-1, 0),(0, 1)]), Set([(1, 0),(0, -1)])]) 

vertex_2 = Set([Set([(-1, 0),(1, 0)])]) 


allowed_vertices = Set([vertex_1,vertex_2])
for vertex_set in copy(allowed_vertices)
    for a in 1:3
        rotated_coords = Set([Set([rotate_point(x, y, a) for (x, y) in s]) for s in vertex_set])
        push!(allowed_vertices, rotated_coords)
    end

    for a in 0:3
        rotated_coords = Set([Set([rotate_point(x, -y, a) for (x, y) in s]) for s in vertex_set])
        push!(allowed_vertices, rotated_coords)
    end
end


"""
has_valid_matching(l, S)

Check if there exists a permutation of the collection S such that each element of l is a subset of the corresponding element in S.

Arguments
l: A collection of sets to match.
S: A collection of sets to match against.

Returns
true if there exists a permutation of S where each element of l is a subset of the corresponding element; false otherwise.
"""

function has_valid_matching(l, S)

    l_list = collect(l)
    S_list = collect(S)
    if length(l_list)>length(S_list)
        return false
    end
    for perm in permutations(S_list)
        if all(issubset(li, si) for (li, si) in zip(l_list, perm))
            return true
        end
    end
    return false
end


"""
    shifted_allowed_vertices(allowed_vertices, x, y)

Shift all vertex sets in allowed_vertices from origin (0,0) to (x,y).

# Arguments
- allowed_vertices: A nested collection of sets of vertex coordinates.
  Each innermost set contains tuples (i, j) representing points.
- x: Horizontal shift amount.
- y: Vertical shift amount.

# Returns
- A new Set containing all shifted vertices, preserving the original nested structure.
"""

function shifted_allowed_vertices(allowed_vertices, x, y)

    # shift the allowed_vertices location from (0,0) to (x,y)

    shifted = Set([])
    for vertex_set in allowed_vertices
        shifted_vertex_set = Set([])
        for set in vertex_set
            push!(shifted_vertex_set, Set([(loc[1]+x,loc[2]+y) for loc in set]))
        end
        push!(shifted, shifted_vertex_set)
    end
    return shifted
end


"""
    is_allowed_vertex(path, new_x, new_y)

Determine if a path together with a candidate new step `(new_x, new_y)` is allowed.

# Arguments
- `path`: A collection of vertices or coordinates representing the current path.
- `new_x`: The x-coordinate of the candidate new vertex.
- `new_y`: The y-coordinate of the candidate new vertex.

# Returns
- `true` if the new vertex can be added to the path according to the rules.
- `false` otherwise.
"""
function is_allowed_vertex(path,new_x,new_y)

    x,y = path[end]

    indices = findall(z -> z == (x,y), path)[1:end-1]
    indices_new = findall(z -> z == (new_x,new_y), path)

    allowed_vertices_xy = shifted_allowed_vertices(allowed_vertices,x,y)

    allowed_vertices_new_x_new_y = shifted_allowed_vertices(allowed_vertices,new_x,new_y)

    if (new_x,new_y) == path[end-1] #No going back
        return false
    end

    if count(p -> p == (new_x, new_y), path) == 2   # A vertex cannot be visited >2 times
        return false
    end

    #Check allowed paths wrt (x,y)

    vertex_set = Set{Any}([Set([path[end-1],(new_x,new_y)])]) #vertex_set related to (x,y)
    for idx in indices
        if idx>1
            push!(vertex_set, Set([path[idx-1],path[idx+1]]))
        else
            push!(vertex_set, Set([path[idx+1]])) #origin
        end
    end

    list = collect(allowed_vertices_xy)

    is_valid = any(has_valid_matching(vertex_set, allowed_vertex_set) for allowed_vertex_set in allowed_vertices_xy)

    #Check allowed paths wrt (new_x,new_y)

    vertex_set_new = Set([Set([(x,y)])]) #vertex_set related to (new_x,new_y)

    for idx in indices_new
        if idx>1
            push!(vertex_set_new, Set([path[idx-1],path[idx+1]]))
        else
            push!(vertex_set_new, Set([path[idx+1]])) #origin
        end
    end


    is_valid *= any(has_valid_matching(vertex_set_new, allowed_vertex_set) for allowed_vertex_set in allowed_vertices_new_x_new_y)


    return is_valid
end


"""
    find_loops(k)

Find all loops of size k (up to reflection/rotation)

# Arguments
- `k`: Size of the loop

# Returns
- List of loops of size k
"""
function find_loops(k)

    if k == 2
        return [[(0, 0), (1, 0), (0, 0)]]
    end

    moves =  [(1, 0), (-1,0),(0,1),(0,-1)]

    solutions = []

    function backtrack(path)
        if length(path) == k + 1
            push!(solutions, copy(path))
            return
        end
        x, y = path[end]  # Get current position
        subpath_less_k = path[max(1, length(path) - (k - 2)):end] #partial path

        for (dx, dy) in moves
            new_x, new_y = x + dx, y + dy

            if (x,y) == (length(path)-1,0) && dy<0 # Skip paths that can be reflected along x to speed up code later
                continue
            end

            if !is_allowed_vertex(subpath_less_k, new_x, new_y) #Skip all invalid paths wrt to the partial path
                continue
            end

            
            # Keep loops of size k
            if length(path) >= k && is_allowed_vertex(path, new_x, new_y)
                continue
            end

            
            push!(path, (new_x, new_y))
            backtrack(path)
            pop!(path)  # Backtrack
        end
    end

    # Start backtracking from different initial directions
    backtrack([(0.0, 0.0), (1.0, 0.0)])

    return solutions
end


"""
evolve(path)

Generate all possible one-step evolutions of a self-osculating walk (SOW) from the current path.

Arguments
path: A collection of coordinates representing the current walk.

Returns
A collection of new paths, each representing a valid one-step extension of the current path.
"""

function evolve(path)
    moves =  [(1, 0), (-1,0),(0,1),(0,-1)]
    solutions = []  # List of valid paths
    x, y = path[end]  # Get current position
    
    path_length = length(path)

    # Precompute the index of the current position in path for later use
    idx_current_pos = findfirst(p -> p == (x, y), path)

    for (dx, dy) in moves
        new_x, new_y = x + dx, y + dy

        if !is_allowed_vertex(path, new_x, new_y) #Skip all invalid paths
                continue
        end

        # Create new path and add to solutions
        new_path = copy(path)  # Create a copy of the current path
        push!(new_path, (new_x, new_y))  # Add the new position to the path
        push!(solutions, new_path)  # Store the new valid path
    end

    return solutions
end


"""
chop_loop(loop)

Generate all proper prefixes (subpaths keeping the starting point) of a loop.

Arguments
loop: A collection representing a loop (e.g., an array of coordinates or elements).

Returns
A vector of sub-loops, where each sub-loop is a prefix of the original loop.
"""
function chop_loop(loop)
    chopped_loop = Vector{Any}()
    for n in 2:length(loop)-1 
        push!(chopped_loop, loop[1:n]) 
    end
    return chopped_loop
end

"""
rotations_reflection(path)

Return all unique rotations and reflections of a given path.

Arguments
path: A collection of coordinates representing the original path (each coordinate is a tuple (x, y)).

Returns
A collection of paths, each being a unique rotation or reflection of the original path.
"""
function rotations_reflection(path)
    rotation_path = [path]

    # Add rotations by 90, 180, 270 degrees
    for a in 1:3
        push!(rotation_path, [rotate_point(x, y, a) for (x, y) in path])
    end
    
    # Add reflections across the x-axis and their rotations
    for a in 0:3
        push!(rotation_path, [rotate_point(x, -y, a) for (x, y) in path])
    end

    return unique(rotation_path)
end

"""
equivalence_k(path, subloops)

Determine which sub-loop a given path is equivalent to.

Arguments
path: A collection of coordinates representing the current path.
subloops: A collection of sub-loops against which to check equivalence.

Returns
The index of the sub-loop in `subloops` to which `path` is equivalent. Returns `nothing` if no equivalent sub-loop is found.
"""
function equivalence_k(path, subloops)
    for n in 1:length(path)-1
        offset = path[n]
        path_translated = [(x - offset[1], y - offset[2]) for (x, y) in path[n:end]]
        
        for rotated in rotations_reflection(path_translated)
            if haskey(subloops, rotated)
                return rotated 
            end
        end
    end
    return nothing
end

"""
automata(k)

Build the transfer matrix to count paths that avoid loops of size less than or equal to `k`.

Arguments
k: An integer specifying the maximum loop size to avoid.

Returns
A sparse matrix used to count all valid paths without loops of size ≤ k efficiently.
"""
function automata(k)
    loop_list = [loop for n in 5:2:k for loop in find_loops(n)]


    chopped_loop_list = unique(vcat([chop_loop(loop) for loop in loop_list]...))

    subloops = Dict(loop => idx for (idx, loop) in enumerate(chopped_loop_list))
    
    matrix = spzeros(UInt8, length(subloops), length(subloops))

    for index in 1:length(chopped_loop_list)
        chopped_loop = chopped_loop_list[index]
        for path in evolve(chopped_loop)
            equivalent_state = equivalence_k(path, subloops)
            if equivalent_state !== nothing
                index_equivalent_state = subloops[equivalent_state]
                matrix[index_equivalent_state, index] +=1
            end
        end
    end

    return matrix
end


"""
approximate_mu(k)

Estimate the connectivity (growth constant) for walks that avoid loops of size less than or equal to `k`.

Arguments
k: An integer specifying the maximum loop size to avoid.

Returns
A numerical value representing the approximate connectivity of the walk.
"""
function approximate_mu(k)
    #return maximum(abs.(eigvals(automata(k))))
    vals, _ = eigs(automata(k), nev=1, which=:LM, maxiter=10000, tol=1e-9)

    return vals[1]
    
end