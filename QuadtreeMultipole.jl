"""

Multipole math

Upward: P2M, M2M
Downward: M2L, L2L, ??, neighbor contribution

"""

include("ArrayQuadtree.jl")

# number of additional coefficients of multipole expansion to keep
# total number of coefficients will be P + 1
const P = 50


"""

Upward pass functions

"""


""" P2M : Particle to Multipole """

function P2M(quadtree::Array{Box, 1}, points::Array{ComplexF64, 1}, masses::Array{Float64, 1}, tree_depth::Int)
  leaf_offset::Int = getOffsetOfDepth(tree_depth)
  # determine multipole expansion at all leaf boxes
  for global_idx in leaf_offset+1:length(quadtree)
    box::Box = quadtree[global_idx]
    # builtin sum uses cumsum so using the builtin for better error
    if (boxHasPoints(box))
      relevant_points = @view points[box.start_idx:box.final_idx]
      relevant_masses = @view masses[box.start_idx:box.final_idx]
      box.a[1] = sum(relevant_masses)
      for i in 2:P+1
        # subtract box center to form multipole expansion about box center
        k = i - 1
        box.a[i] = -1/k*sum(((relevant_points .- box.center).^k).*relevant_masses)
      end
    else 
      box.a = zeros(P+1)
    end
  end
end


""" M2M : Multipole to multipole"""

function M2M
  # propogate up multipole expansions from the leaves
  # to the highest boxes at depth 2
  for depth in max_depth-1:-1:1
    for idx in 1:4^depth
      global_idx::Int = depth_offsets[depth] + idx
      parent_box::Box = quadtree[global_idx]
      # zero out as will be adding on contributions from children
      # do not want data from previous timestep
      parent_box.a .= zero(parent_box.a(1))
      for child_idx in parent_box.children_idxs
        child_global_idx::Int = depth_offsets[depth+1] + child_idx
        child_box::Box = quadtree[child_global_idx]
        parent_box.a(1) += child_box.a(1)
        for l in 1:P
          i = l + i
          parent_box.a(i) -= 1/l*child_box.a(1)*(child_box.center - parent_box.center).^l
          for k in 1:l
            j = k + 1
            parent_box.a(i) += binomial(l, k)*child_box.a(j)*(child_box.center - parent_box.center).^(l-k)
          end
        end
      end
    end
  end
end




