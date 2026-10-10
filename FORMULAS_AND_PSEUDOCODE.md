# Formula definitions and executable-code inventory

F01--F27 identify computational definitions, assumptions, units and pseudocode. Source positions are generated from the executable Python syntax tree.

## F01: Domains and volume

$$\Omega=\prod_{k=1}^{3}[a_k,b_k],\qquad |\Omega|=\prod_{k=1}^{3}(b_k-a_k)$$

Units: Bounds: m; volume: m cubed; normalized weights: dimensionless.

Source: code/formula_reference.py:13 (domain_volume); code/run_experiments.py:28 (quad)

Implemented and checked against mathematical contracts.

1. Check that each upper bound exceeds its lower bound.
2. Multiply the three positive side lengths.
3. For target averaging, normalize each one-dimensional quadrature weight by its side length.

## F02: Gaussian density and physical truncation

$$\phi(v)=\frac{e^{-v^2/2}}{\sqrt{2\pi}},\quad Z=\Phi(\beta)-\Phi(\alpha),\quad f(x)=\frac{\phi((x-\mu)/\sigma)}{\sigma Z}\mathbf1_{[a,b]}(x)$$

Units: Density: inverse m; sigma: m; Z: dimensionless.

Source: code/formula_reference.py:19 (normal_pdf); code/formula_reference.py:23 (truncation_normalizer); code/formula_reference.py:30 (truncated_density)

Implemented and checked against mathematical contracts.

1. Set alpha = (a - mu)/sigma and beta = (b - mu)/sigma.
2. Compute the Gaussian probability mass between the bounds.
3. When alpha is positive, use Phi(-alpha) - Phi(-beta) to avoid cancellation.
4. Evaluate the parent density and divide by sigma times the normalizer; return zero outside the bounds.

## F03: Actual-depth mean and additive error

$$\mathbb E[S_z]=s_z+\sigma_z\frac{\phi(\alpha_z)-\phi(\beta_z)}{Z_z},\qquad\mathbb E[\varepsilon_z]=\mathbb E[S_z]-s_z$$

Units: All depths and errors: m.

Source: code/formula_reference.py:35 (truncated_mean)

Implemented and checked against mathematical contracts.

1. Evaluate the density at standardized lower and upper bounds.
2. Divide their difference by the truncation normalizer and multiply by sigma.
3. Add the parent centre to obtain the actual mean.
4. Subtract the nominal centre to obtain the mean additive error.

## F04: Joint coordinate density

$$f_i(\mathbf x)=\prod_{k\in\{x,y,z\}}f_{ik}(x_k)$$

Units: Joint density: inverse m cubed.

Source: code/formula_reference.py:40 (joint_density)

Implemented and checked against mathematical contracts.

1. Evaluate the three bounded marginal densities.
2. Multiply them only under the declared coordinate-independence assumption.

## F05: Spherical sensing probability

$$p_i(\mathbf q)=\int_{\mathcal F}\mathbf1_{\{\|\mathbf x-\mathbf q\|\le R\}}f_i(\mathbf x)\,d\mathbf x$$

Units: Probability: dimensionless; distance and radius: m.

Source: code/run_experiments.py:41 (pointwise); code/formula_reference.py:44 (spherical_hit)

Implemented and checked against mathematical contracts.

1. Define the Euclidean sphere event, not a bounding cuboid.
2. Intersect the sphere with physical depth bounds.
3. Evaluate the equivalent disk-CDF/depth integral described in F09 and F10.

## F06: Pointwise network union

$$P(\mathbf q)=1-\prod_i(1-p_i(\mathbf q))$$

Units: Dimensionless; independent sensors, or conditional independence given a common current.

Source: code/formula_reference.py:48 (pointwise_union); code/run_experiments.py:57 (objective)

Implemented and checked against mathematical contracts.

1. At the same target point, collect every sensor probability.
2. Multiply their complements across sensors.
3. Subtract the product from one.
4. Repeat for every target point before any spatial averaging.

## F07: Expected covered volume and numerical target averaging

$$C=\frac1{|\Omega|}\int_\Omega P(\mathbf q)\,d\mathbf q=\frac{\mathbb E[|\Omega\cap\bigcup_i B(\mathbf S_i,R)|]}{|\Omega|}\simeq\sum_jw_jP(\mathbf q_j)$$

Units: Coverage fraction: dimensionless; multiply by 100 for percent.

Source: code/run_experiments.py:28 (quad); code/formula_reference.py:52 (spatial_average); code/run_experiments.py:57 (objective)

Implemented and checked against mathematical contracts.

1. Build tensor-product Gauss--Legendre target points and weights.
2. Check that the weights sum to one.
3. Compute a network union separately at each point.
4. Return the weighted sum of these unions.

## F08: Aggregation-order covariance audit

$$C=\bar p_1+\bar p_2-\mathbb E_Q[p_1p_2],\quad 1-(1-\bar p_1)(1-\bar p_2)-C=\operatorname{Cov}_Q(p_1,p_2)$$

Units: Dimensionless.

Source: code/formula_reference.py:58 (covariance_audit)

Implemented and checked against mathematical contracts.

1. Compute each spatial mean and the spatial mean of their product.
2. Evaluate correct pointwise-union averaging.
3. Evaluate the product-of-spatial-means shortcut.
4. Check that shortcut minus correct equals the spatial covariance.

## F09: Horizontal noncentrality and disk probability

$$\eta_i(\mathbf q)=\frac{(q_x-s_{ix})^2+(q_y-s_{iy})^2}{\sigma_h^2},\quad F_{\chi^2_2(\eta_i)}(r^2/\sigma_h^2)$$

Units: Dimensionless.

Source: code/formula_reference.py:67 (horizontal_noncentrality); code/run_experiments.py:41 (pointwise)

Implemented and checked against mathematical contracts.

1. Divide squared horizontal separation by parent horizontal variance.
2. Evaluate a noncentral chi-square CDF with two degrees of freedom.
3. Normalize by the horizontal physical-conditioning mass.

## F10: Depth integration and zero-depth limit

$$\ell=\max(a_z,q_z-R),\quad h=\max(\ell,\min(b_z,q_z+R)),\quad p_i=\frac1{Z_{ix}Z_{iy}}\int_\ell^h F_{\chi^2_2(\eta_i)}\!\left(\frac{R^2-(q_z-z)^2}{\sigma_h^2}\right)f_{iz}(z)\,dz$$

Units: Integration depth: m; final probability: dimensionless.

Source: code/formula_reference.py:71 (depth_limits); code/run_experiments.py:41 (pointwise)

Implemented and checked against mathematical contracts.

1. Clip the depth interval to the sphere and physical bounds.
2. Map Gauss--Legendre nodes into this interval.
3. At each depth compute the horizontal disk CDF times the bounded depth density.
4. Multiply by interval-scaled quadrature weights and horizontal normalization.
5. If sigma_z is zero, use the disk at nominal depth directly; retain volume averaging.

## F11: Structured grid and nominal separation

$$\mathbf s_{ij}=(c_x+id,c_y+jd,z_0),\quad i,j\in\{-1,0,1\},\quad d\ge40\ \mathrm m$$

Units: Coordinates and pitch: m; actual positions are not collision-constrained.

Source: code/run_experiments.py:39 (layout)

Implemented and checked against mathematical contracts.

1. Enumerate the nine index pairs (-1, 0, 1) by (-1, 0, 1).
2. Construct coordinates at the common nominal depth.
3. The smallest distinct grid separation equals the pitch.

## F12: Bounded parameterization and objective

$$d=40+460\theta_3,\quad z_0=120+150\theta_4,\quad c_{x,y}=(2\theta_{1,2}-1)(500-d),\quad\max_{\theta\in[0,1]^4}C(\mathcal S(\theta))$$

Units: Optimization variables: dimensionless; decoded coordinates: m.

Source: code/run_experiments.py:35 (decode); code/run_experiments.py:61 (optimize)

Implemented and checked against mathematical contracts.

1. Map unit-box variables to pitch and depth.
2. Map the first two variables to feasible grid-centre coordinates.
3. Construct all nine nominal coordinates.
4. Evaluate the shared deterministic objective and increment the call counter.
5. Retain the best evaluated candidate; stop at the declared limit or Powell tolerance.

## F13: Monte Carlo union observation and batch mean

$$I_{rm}=\mathbf1_{\{\min_i\|\mathbf S_{irm}-\mathbf Q_{rm}\|\le R\}},\qquad\widehat C_r=M^{-1}\sum_{m=1}^{M}I_{rm}$$

Units: Bernoulli observations and batch fraction: dimensionless.

Source: code/formula_reference.py:44 (spherical_hit); code/run_experiments.py:118 (main); code/run_ocean_experiments.py:82 (run)

Implemented and checked against mathematical contracts.

1. Draw an independent target uniformly in the monitoring volume.
2. Draw the entire array according to the same bounded position model.
3. For current replay, draw one held-out current day shared across sensors.
4. Return one if at least one sensor is within the spherical radius.
5. Average observations in each independent batch; archive packed bits and counts.

## F14: Inverse-CDF actual-position sampling

$$S=\mu+\sigma\Phi^{-1}\{\Phi(\alpha)+U[\Phi(\beta)-\Phi(\alpha)]\},\quad U\sim\mathcal U(0,1)$$

Units: Position: m; uniform variates: dimensionless.

Source: code/formula_reference.py:76 (sample_actual); code/run_experiments.py:118 (main)

Implemented and checked against mathematical contracts.

1. Generate independent uniform variates for each coordinate, sensor and trial.
2. Map each variate into the parent-CDF interval between physical bounds.
3. Apply the inverse Gaussian CDF and rescale.
4. Use the equivalent survival-tail identity for a far-left parent centre.

## F15: Independent-batch confidence intervals

$$\bar C\ \pm\ t_{0.975,n-1}\frac{s_C}{\sqrt n},\quad n=30$$

Units: Fraction; display conversion to percent is explicit.

Source: code/formula_reference.py:83 (batch_interval); code/run_experiments.py:118 (main)

Implemented and checked against mathematical contracts.

1. Treat the 30 independent batch estimates as replicates.
2. Calculate their sample mean and sample standard deviation.
3. Multiply the standard error by the Student-t critical value.
4. For paired gains, apply the same calculation to within-batch differences.

## F16: Coverage gain and expected uncovered fraction

$$G_{\rm pp}=100(C_o-C_f),\quad G_{\rm rel}=100(C_o/C_f-1),\quad G_{\rm unc}=100(C_o-C_f)/(1-C_f)$$

Units: First output: percentage points; remaining outputs: relative percent.

Source: code/formula_reference.py:90 (coverage_gains)

Implemented and checked against mathematical contracts.

1. Subtract fixed coverage from optimized coverage for absolute gain.
2. Divide by fixed coverage for relative coverage gain.
3. Divide the decrease in uncovered fraction by initial uncovered fraction.
4. Multiply by 100 and label percentage points versus relative percent correctly.

## F17: PSO update and bounds

$$\mathbf v\leftarrow0.7298\mathbf v+1.49618\mathbf r_1\odot(\mathbf p-\mathbf u)+1.49618\mathbf r_2\odot(\mathbf g-\mathbf u)$$

Units: Unit-box coordinates and velocities: dimensionless.

Source: code/formula_reference.py:96 (pso_update); code/run_experiments.py:61 (optimize)

Implemented and checked against mathematical contracts.

1. Draw independent componentwise uniform r1 and r2.
2. Update and clip velocity to [-0.2, 0.2].
3. Add velocity to position and project to the unit box.
4. Zero velocity components whose proposed position crossed a bound.
5. Evaluate candidates and update personal/global best records; count initialization.

## F18: AFSA schedule and specified variant

$$\rho=n_{\rm eval}/B,\quad v=0.5(1-\rho)+0.05,\quad h=0.1(1-\rho)+0.01$$

Units: Unit-box distances: dimensionless.

Source: code/formula_reference.py:104 (afsa_schedule); code/run_experiments.py:61 (optimize)

Implemented and checked against mathematical contracts.

1. Set visual and movement scales from the fraction of consumed calls.
2. Find neighbors within the visual radius.
3. With 1--12 neighbors, consider an improving best neighbor and evaluate the centroid.
4. Try at most three projected uniform-ball prey candidates, stopping at the first improvement.
5. Move toward the best eligible target by a capped random step, or use a projected random movement.
6. Evaluate the moved point; every exploratory call updates the best-observed bulletin.

## F19: Prescribed synthetic bias and ideal compensation

$$\mathbf b_k=(80\sin(k\pi/6),40\sin(k\pi/3)),\quad\mathbf s_{i,k}=\mathbf s_i^*-(b_{kx},b_{ky},0)$$

Units: Bias and position: m; steps have no assigned physical duration.

Source: code/formula_reference.py:109 (prescribed_bias); code/run_experiments.py:118 (main)

Implemented and checked against mathematical contracts.

1. Generate the declared sinusoidal horizontal bias for steps 0--12.
2. Shift actual parent centres for the static layout.
3. For ideal compensation, subtract the known bias from nominal coordinates.
4. Reevaluate the identical physically truncated coverage law.

## F20: Movement-distance proxy

$$D=\sum_{k,i}\|\mathbf s_{i,k}-\mathbf s_{i,k-1}\|_2$$

Units: Metres; this is not energy.

Source: code/formula_reference.py:113 (movement_distance); code/run_experiments.py:118 (main)

Implemented and checked against mathematical contracts.

1. Take consecutive position differences for each node.
2. Compute Euclidean distances and sum across nodes and steps.

## F21: Energy-aware extension interface

$$J=C-\lambda E$$

Units: E: joules; lambda: inverse joules; J: dimensionless.

Source: code/formula_reference.py:117 (energy_penalized_objective)

Prospective interface only. No empirical energy advantage is reported.

1. Require a separately calibrated energy input measured in joules.
2. If it is missing, raise an error instead of generating an energy result.
3. Subtract weighted energy from coverage only when that external input is available.

## F22: Frozen-current displacement and parent centre

$$\mathbf b_t=\Delta t(u_t,v_t),\qquad\boldsymbol\mu_{it}=\mathbf s_i+(b_{tx},b_{ty},0)$$

Units: Velocity: m/s; horizon: s; derived displacement: m.

Source: code/formula_reference.py:122 (current_displacement); code/run_ocean_experiments.py:24 (conditional_coverage)

Current forcing is real reanalysis; the displacement rule and residual law are assumptions.

1. Read real eastward and northward daily reanalysis velocities.
2. Multiply by the declared horizon in seconds.
3. Add the same horizontal displacement to every sensor parent centre.
4. Condition residual coordinates on the same physical bounds.

## F23: Coverage averaged over common environmental forcing

$$C_{\rm env}=T^{-1}\sum_{t=1}^{T}\frac1{|\Omega|}\int_\Omega\left[1-\prod_i(1-p_i(\mathbf q\mid\mathbf b_t))\right]d\mathbf q$$

Units: Dimensionless.

Source: code/formula_reference.py:126 (empirical_current_mean); code/run_ocean_experiments.py:24 (conditional_coverage); code/run_ocean_experiments.py:82 (run)

Implemented and checked against mathematical contracts.

1. Condition every sensor probability on the same current state.
2. Form the conditional union at each target point.
3. Integrate over the target volume.
4. Average these complete conditional coverage values across current states.
5. Do not multiply sensor complements after unconditional current averaging.

## F24: Temporal block-bootstrap interval

$$\widehat G^*_b=T^{-1}\sum_{j=1}^{T}g_{I_{bj}},\qquad\mathrm{CI}=[Q_{.025}(\widehat G^*),Q_{.975}(\widehat G^*)]$$

Units: Fraction; temporal-record variability is distinct from independent Monte Carlo error.

Source: code/formula_reference.py:133 (block_bootstrap_interval)

Implemented and checked against mathematical contracts.

1. Compute the paired daily coverage difference on held-out dates.
2. Sample circular block starts with replacement.
3. Concatenate 14-day blocks and retain the first 365 indices.
4. Repeat 2000 times with a fixed seed and take the percentile endpoints.

## F25: Conditional spherical integral with drift

$$\eta_{it}=\|\mathbf q_h-\boldsymbol\mu_{it,h}\|^2/\sigma_h^2,\qquad Z_{it,k}=\Phi((b_k-\mu_{it,k})/\sigma_k)-\Phi((a_k-\mu_{it,k})/\sigma_k)$$

Units: Dimensionless probability and normalizers.

Source: code/run_ocean_experiments.py:24 (conditional_coverage)

Implemented and checked against mathematical contracts.

1. Construct displaced horizontal parent centres for a bounded chunk of current days.
2. Recompute both horizontal physical normalizers for each centre.
3. Use F09/F10 with the displaced noncentrality and unchanged depth density.
4. Aggregate conditionally using F06/F07; return one coverage value per day.

## F26: Training-state reduction and weighted objective

$$w_k=n_k/T_{\rm train},\quad\bar{\mathbf b}_k=n_k^{-1}\sum_{t\in\mathcal C_k}\mathbf b_t,\quad C_K=\sum_kw_kC(\mathcal S\mid\bar{\mathbf b}_k)$$

Units: Centres: m; weights and coverage: dimensionless.

Source: code/run_ocean_experiments.py:53 (training_states); code/run_ocean_experiments.py:82 (run)

Implemented and checked against mathematical contracts.

1. Use only training-year displacements to fit 16 k-means++ representatives.
2. Run 60 centroid iterations with seed 6200; assign final nearest-centre labels.
3. Divide cluster counts by the training-date count to obtain weights.
4. Optimize the weighted conditional coverage with a matched budget.
5. Select by refined training coverage; evaluate all held-out dates without tuning.
6. Compare 16 and 32 representatives against full-year training evaluation.

## F27: Random design sample and coordinate plots

$$c_x,c_y\sim\mathcal U(-160,160),\quad z_0\sim\mathcal U(120,270),\quad d=40\ \mathrm m$$

Units: Coordinates: m; coverage: fraction.

Source: code/run_experiments.py:118 (main)

Implemented and checked against mathematical contracts.

1. With seed 43210, draw 600 independent feasible design vectors.
2. Construct the fixed-pitch grid and evaluate the deterministic objective.
3. Group by coordinate intervals for descriptive plots; do not infer significance from meshes.
