
function simulate(curr_pos, prev_pos, tmp_pos, mass, acc, timesteps=TIMESTEPS)
  """
  # using Plots + GR to create a gif
  gr() # GR backend
  anim = Animation()
  for i ∈ 1:TIMESTEPS
    updateAcc!(acc, curr_pos, mass)
    # CHECK/VERIFY: verlet integration
    tmp_pos = curr_pos
    curr_pos .= 2*curr_pos - prev_pos + acc * Δt^2
    prev_pos = tmp_pos
    assert(curr_pos != prev_pos)
    x = @view curr_pos[:, 1]
    y = @view curr_pos[:, 2]
    scatter(x, y, xlim = (-10, 10), ylim = (-10, 10))
    frame(anim)
  end
  gif(anim, "fmm/n_body.gif", fps = 5)
  """

  """
  # using Plots + GR realtime display plotting all previous positions
  #gr(show=true)
  gr(reuse=true, size = (1000, 1000))
  x = @view curr_pos[:, 1]
  y = @view curr_pos[:, 2]
  plot = scatter(x, y, xlim = (-2, 10), ylim = (-2, 2), legend=false)
  for i ∈ 1:timesteps
    updateAcc!(acc, curr_pos, mass)
    # CHECK/VERIFY: verlet integration
    tmp_pos = copy(curr_pos)
    curr_pos .= 2*curr_pos - prev_pos + acc * Δt^2
    prev_pos = tmp_pos
    @assert curr_pos != prev_pos
    x = @view curr_pos[:, 1]
    y = @view curr_pos[:, 2]
    scatter!(plot, x, y, xlim = (-2, 2), ylim = (-2, 2), legend=false)
    gui() 
  end
  """

  # using Plots + GR realtime display plotting current position
  #gr(show=true)
  gr(reuse=true, size = (1000, 1000))
  for i ∈ 1:timesteps
    updateAcc!(acc, curr_pos, mass)
    # CHECK/VERIFY: verlet integration
    tmp_pos = copy(curr_pos)
    curr_pos .= 2*curr_pos - prev_pos + acc * Δt^2
    prev_pos = tmp_pos
    @assert curr_pos != prev_pos
    x = @view curr_pos[:, 1]
    y = @view curr_pos[:, 2]
    scatter((x, y), xlim = (-2, 2), ylim = (-2, 2), color=1, label="", legend=false)
    gui() 
  end


end
