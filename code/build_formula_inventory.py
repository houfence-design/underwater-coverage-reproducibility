"""Generate the complete equation/pseudocode/executable-source concordance."""
from pathlib import Path
import ast,json,re

ROOT=Path(__file__).resolve().parents[1];CODE=ROOT/'code'
entries=[]
def entry(number,title,location,formula,functions,pseudocode,unit,note='Implemented and checked against mathematical contracts.'):
    entries.append({'id':f'F{number:02d}','title':title,'manuscript_location':location,
        'formula_latex':formula,'functions':functions,'pseudocode':pseudocode,'units':unit,'status':note})

entry(1,'Domains and volume','Main Eq. (1), domain dimensions',
 r'\Omega=\prod_{k=1}^{3}[a_k,b_k],\qquad |\Omega|=\prod_{k=1}^{3}(b_k-a_k)',
 ['formula_reference.domain_volume','run_experiments.quad'],
 ['Check that each upper bound exceeds its lower bound.','Multiply the three positive side lengths.','For target averaging, normalize each one-dimensional quadrature weight by its side length.'],
 'Bounds: m; volume: m cubed; normalized weights: dimensionless.')
entry(2,'Gaussian density and physical truncation','Main Eq. (2)',
 r'\phi(v)=\frac{e^{-v^2/2}}{\sqrt{2\pi}},\quad Z=\Phi(\beta)-\Phi(\alpha),\quad f(x)=\frac{\phi((x-\mu)/\sigma)}{\sigma Z}\mathbf1_{[a,b]}(x)',
 ['formula_reference.normal_pdf','formula_reference.truncation_normalizer','formula_reference.truncated_density'],
 ['Set alpha = (a - mu)/sigma and beta = (b - mu)/sigma.','Compute the Gaussian probability mass between the bounds.','When alpha is positive, use Phi(-alpha) - Phi(-beta) to avoid cancellation.','Evaluate the parent density and divide by sigma times the normalizer; return zero outside the bounds.'],
 'Density: inverse m; sigma: m; Z: dimensionless.')
entry(3,'Actual-depth mean and additive error','Main Eq. (3), actual-minus-nominal definition',
 r'\mathbb E[S_z]=s_z+\sigma_z\frac{\phi(\alpha_z)-\phi(\beta_z)}{Z_z},\qquad\mathbb E[\varepsilon_z]=\mathbb E[S_z]-s_z',
 ['formula_reference.truncated_mean'],
 ['Evaluate the density at standardized lower and upper bounds.','Divide their difference by the truncation normalizer and multiply by sigma.','Add the parent centre to obtain the actual mean.','Subtract the nominal centre to obtain the mean additive error.'], 'All depths and errors: m.')
entry(4,'Joint coordinate density','Main text after Eq. (3)',
 r'f_i(\mathbf x)=\prod_{k\in\{x,y,z\}}f_{ik}(x_k)',
 ['formula_reference.joint_density'],
 ['Evaluate the three bounded marginal densities.','Multiply them only under the declared coordinate-independence assumption.'],
 'Joint density: inverse m cubed.')
entry(5,'Spherical sensing probability','Main Eq. (4)',
 r'p_i(\mathbf q)=\int_{\mathcal F}\mathbf1_{\{\|\mathbf x-\mathbf q\|\le R\}}f_i(\mathbf x)\,d\mathbf x',
 ['run_experiments.pointwise','formula_reference.spherical_hit'],
 ['Define the Euclidean sphere event, not a bounding cuboid.','Intersect the sphere with physical depth bounds.','Evaluate the equivalent disk-CDF/depth integral described in F09 and F10.'],
 'Probability: dimensionless; distance and radius: m.')
entry(6,'Pointwise network union','Main Eq. (5)',
 r'P(\mathbf q)=1-\prod_i(1-p_i(\mathbf q))',
 ['formula_reference.pointwise_union','run_experiments.objective'],
 ['At the same target point, collect every sensor probability.','Multiply their complements across sensors.','Subtract the product from one.','Repeat for every target point before any spatial averaging.'],
 'Dimensionless; independent sensors, or conditional independence given a common current.')
entry(7,'Expected covered volume and numerical target averaging','Main Eq. (6) and Fubini interpretation',
 r'C=\frac1{|\Omega|}\int_\Omega P(\mathbf q)\,d\mathbf q=\frac{\mathbb E[|\Omega\cap\bigcup_i B(\mathbf S_i,R)|]}{|\Omega|}\simeq\sum_jw_jP(\mathbf q_j)',
 ['run_experiments.quad','formula_reference.spatial_average','run_experiments.objective'],
 ['Build tensor-product Gauss--Legendre target points and weights.','Check that the weights sum to one.','Compute a network union separately at each point.','Return the weighted sum of these unions.'],
 'Coverage fraction: dimensionless; multiply by 100 for percent.')
entry(8,'Aggregation-order covariance audit','Main Eq. (7)',
 r'C=\bar p_1+\bar p_2-\mathbb E_Q[p_1p_2],\quad 1-(1-\bar p_1)(1-\bar p_2)-C=\operatorname{Cov}_Q(p_1,p_2)',
 ['formula_reference.covariance_audit'],
 ['Compute each spatial mean and the spatial mean of their product.','Evaluate correct pointwise-union averaging.','Evaluate the product-of-spatial-means shortcut.','Check that shortcut minus correct equals the spatial covariance.'], 'Dimensionless.')
entry(9,'Horizontal noncentrality and disk probability','Main Eq. (8)',
 r'\eta_i(\mathbf q)=\frac{(q_x-s_{ix})^2+(q_y-s_{iy})^2}{\sigma_h^2},\quad F_{\chi^2_2(\eta_i)}(r^2/\sigma_h^2)',
 ['formula_reference.horizontal_noncentrality','run_experiments.pointwise'],
 ['Divide squared horizontal separation by parent horizontal variance.','Evaluate a noncentral chi-square CDF with two degrees of freedom.','Normalize by the horizontal physical-conditioning mass.'], 'Dimensionless.')
entry(10,'Depth integration and zero-depth limit','Main Eq. (9), following zero-scale expression',
 r'\ell=\max(a_z,q_z-R),\quad h=\max(\ell,\min(b_z,q_z+R)),\quad p_i=\frac1{Z_{ix}Z_{iy}}\int_\ell^h F_{\chi^2_2(\eta_i)}\!\left(\frac{R^2-(q_z-z)^2}{\sigma_h^2}\right)f_{iz}(z)\,dz',
 ['formula_reference.depth_limits','run_experiments.pointwise'],
 ['Clip the depth interval to the sphere and physical bounds.','Map Gauss--Legendre nodes into this interval.','At each depth compute the horizontal disk CDF times the bounded depth density.','Multiply by interval-scaled quadrature weights and horizontal normalization.','If sigma_z is zero, use the disk at nominal depth directly; retain volume averaging.'],
 'Integration depth: m; final probability: dimensionless.')
entry(11,'Structured grid and nominal separation','Main Eq. (10)',
 r'\mathbf s_{ij}=(c_x+id,c_y+jd,z_0),\quad i,j\in\{-1,0,1\},\quad d\ge40\ \mathrm m',
 ['run_experiments.layout'],
 ['Enumerate the nine index pairs (-1, 0, 1) by (-1, 0, 1).','Construct coordinates at the common nominal depth.','The smallest distinct grid separation equals the pitch.'], 'Coordinates and pitch: m; actual positions are not collision-constrained.')
entry(12,'Bounded parameterization and objective','Main Eq. (11) and maximization statement',
 r'd=40+460\theta_3,\quad z_0=120+150\theta_4,\quad c_{x,y}=(2\theta_{1,2}-1)(500-d),\quad\max_{\theta\in[0,1]^4}C(\mathcal S(\theta))',
 ['run_experiments.decode','run_experiments.optimize'],
 ['Map unit-box variables to pitch and depth.','Map the first two variables to feasible grid-centre coordinates.','Construct all nine nominal coordinates.','Evaluate the shared deterministic objective and increment the call counter.','Retain the best evaluated candidate; stop at the declared limit or Powell tolerance.'], 'Optimization variables: dimensionless; decoded coordinates: m.')
entry(13,'Monte Carlo union observation and batch mean','Main Eq. (12)',
 r'I_{rm}=\mathbf1_{\{\min_i\|\mathbf S_{irm}-\mathbf Q_{rm}\|\le R\}},\qquad\widehat C_r=M^{-1}\sum_{m=1}^{M}I_{rm}',
 ['formula_reference.spherical_hit','run_experiments.main','run_ocean_experiments.run'],
 ['Draw an independent target uniformly in the monitoring volume.','Draw the entire array according to the same bounded position model.','For current replay, draw one held-out current day shared across sensors.','Return one if at least one sensor is within the spherical radius.','Average observations in each independent batch; archive packed bits and counts.'],
 'Bernoulli observations and batch fraction: dimensionless.')
entry(14,'Inverse-CDF actual-position sampling','Main Monte Carlo sampling paragraph',
 r'S=\mu+\sigma\Phi^{-1}\{\Phi(\alpha)+U[\Phi(\beta)-\Phi(\alpha)]\},\quad U\sim\mathcal U(0,1)',
 ['formula_reference.sample_actual','run_experiments.main'],
 ['Generate independent uniform variates for each coordinate, sensor and trial.','Map each variate into the parent-CDF interval between physical bounds.','Apply the inverse Gaussian CDF and rescale.','Use the equivalent survival-tail identity for a far-left parent centre.'], 'Position: m; uniform variates: dimensionless.')
entry(15,'Independent-batch confidence intervals','Main Monte Carlo confidence interval paragraph',
 r'\bar C\ \pm\ t_{0.975,n-1}\frac{s_C}{\sqrt n},\quad n=30',
 ['formula_reference.batch_interval','run_experiments.main'],
 ['Treat the 30 independent batch estimates as replicates.','Calculate their sample mean and sample standard deviation.','Multiply the standard error by the Student-t critical value.','For paired gains, apply the same calculation to within-batch differences.'], 'Fraction; display conversion to percent is explicit.')
entry(16,'Coverage gain and expected uncovered fraction','Main Results RQ3 and uncovered-volume interpretation',
 r'G_{\rm pp}=100(C_o-C_f),\quad G_{\rm rel}=100(C_o/C_f-1),\quad G_{\rm unc}=100(C_o-C_f)/(1-C_f)',
 ['formula_reference.coverage_gains'],
 ['Subtract fixed coverage from optimized coverage for absolute gain.','Divide by fixed coverage for relative coverage gain.','Divide the decrease in uncovered fraction by initial uncovered fraction.','Multiply by 100 and label percentage points versus relative percent correctly.'], 'First output: percentage points; remaining outputs: relative percent.')
entry(17,'PSO update and bounds','Supplementary Methods S1 PSO expression',
 r'\mathbf v\leftarrow0.7298\mathbf v+1.49618\mathbf r_1\odot(\mathbf p-\mathbf u)+1.49618\mathbf r_2\odot(\mathbf g-\mathbf u)',
 ['formula_reference.pso_update','run_experiments.optimize'],
 ['Draw independent componentwise uniform r1 and r2.','Update and clip velocity to [-0.2, 0.2].','Add velocity to position and project to the unit box.','Zero velocity components whose proposed position crossed a bound.','Evaluate candidates and update personal/global best records; count initialization.'], 'Unit-box coordinates and velocities: dimensionless.')
entry(18,'AFSA schedule and specified variant','Supplementary Methods S1 AFSA expressions',
 r'\rho=n_{\rm eval}/B,\quad v=0.5(1-\rho)+0.05,\quad h=0.1(1-\rho)+0.01',
 ['formula_reference.afsa_schedule','run_experiments.optimize'],
 ['Set visual and movement scales from the fraction of consumed calls.','Find neighbors within the visual radius.','With 1--12 neighbors, consider an improving best neighbor and evaluate the centroid.','Try at most three projected uniform-ball prey candidates, stopping at the first improvement.','Move toward the best eligible target by a capped random step, or use a projected random movement.','Evaluate the moved point; every exploratory call updates the best-observed bulletin.'], 'Unit-box distances: dimensionless.')
entry(19,'Prescribed synthetic bias and ideal compensation','Main Discussion; Supplementary Analysis S4',
 r'\mathbf b_k=(80\sin(k\pi/6),40\sin(k\pi/3)),\quad\mathbf s_{i,k}=\mathbf s_i^*-(b_{kx},b_{ky},0)',
 ['formula_reference.prescribed_bias','run_experiments.main'],
 ['Generate the declared sinusoidal horizontal bias for steps 0--12.','Shift actual parent centres for the static layout.','For ideal compensation, subtract the known bias from nominal coordinates.','Reevaluate the identical physically truncated coverage law.'], 'Bias and position: m; steps have no assigned physical duration.')
entry(20,'Movement-distance proxy','Main Discussion; Supplementary Analysis S4',
 r'D=\sum_{k,i}\|\mathbf s_{i,k}-\mathbf s_{i,k-1}\|_2',
 ['formula_reference.movement_distance','run_experiments.main'],
 ['Take consecutive position differences for each node.','Compute Euclidean distances and sum across nodes and steps.'], 'Metres; this is not energy.')
entry(21,'Energy-aware extension interface','Main Discussion energy objective',
 r'J=C-\lambda E',
 ['formula_reference.energy_penalized_objective'],
 ['Require a separately calibrated energy input measured in joules.','If it is missing, raise an error instead of generating an energy result.','Subtract weighted energy from coverage only when that external input is available.'],
 'E: joules; lambda: inverse joules; J: dimensionless.',
 'Prospective interface only. No empirical energy advantage is reported.')
entry(22,'Frozen-current displacement and parent centre','Main current-replay Methods',
 r'\mathbf b_t=\Delta t(u_t,v_t),\qquad\boldsymbol\mu_{it}=\mathbf s_i+(b_{tx},b_{ty},0)',
 ['formula_reference.current_displacement','run_ocean_experiments.conditional_coverage'],
 ['Read real eastward and northward daily reanalysis velocities.','Multiply by the declared horizon in seconds.','Add the same horizontal displacement to every sensor parent centre.','Condition residual coordinates on the same physical bounds.'],
 'Velocity: m/s; horizon: s; derived displacement: m.',
 'Current forcing is real reanalysis; the displacement rule and residual law are assumptions.')
entry(23,'Coverage averaged over common environmental forcing','Main current-replay objective',
 r'C_{\rm env}=T^{-1}\sum_{t=1}^{T}\frac1{|\Omega|}\int_\Omega\left[1-\prod_i(1-p_i(\mathbf q\mid\mathbf b_t))\right]d\mathbf q',
 ['formula_reference.empirical_current_mean','run_ocean_experiments.conditional_coverage','run_ocean_experiments.run'],
 ['Condition every sensor probability on the same current state.','Form the conditional union at each target point.','Integrate over the target volume.','Average these complete conditional coverage values across current states.','Do not multiply sensor complements after unconditional current averaging.'], 'Dimensionless.')
entry(24,'Temporal block-bootstrap interval','Main current-replay statistical Methods',
 r'\widehat G^*_b=T^{-1}\sum_{j=1}^{T}g_{I_{bj}},\qquad\mathrm{CI}=[Q_{.025}(\widehat G^*),Q_{.975}(\widehat G^*)]',
 ['formula_reference.block_bootstrap_interval'],
 ['Compute the paired daily coverage difference on held-out dates.','Sample circular block starts with replacement.','Concatenate 14-day blocks and retain the first 365 indices.','Repeat 2000 times with a fixed seed and take the percentile endpoints.'],
 'Fraction; temporal-record variability is distinct from independent Monte Carlo error.')
entry(25,'Conditional spherical integral with drift','Main current-replay Methods, F10 with shifted horizontal centres',
 r'\eta_{it}=\|\mathbf q_h-\boldsymbol\mu_{it,h}\|^2/\sigma_h^2,\qquad Z_{it,k}=\Phi((b_k-\mu_{it,k})/\sigma_k)-\Phi((a_k-\mu_{it,k})/\sigma_k)',
 ['run_ocean_experiments.conditional_coverage'],
 ['Construct displaced horizontal parent centres for a bounded chunk of current days.','Recompute both horizontal physical normalizers for each centre.','Use F09/F10 with the displaced noncentrality and unchanged depth density.','Aggregate conditionally using F06/F07; return one coverage value per day.'], 'Dimensionless probability and normalizers.')
entry(26,'Training-state reduction and weighted objective','Main current-replay optimization settings',
 r'w_k=n_k/T_{\rm train},\quad\bar{\mathbf b}_k=n_k^{-1}\sum_{t\in\mathcal C_k}\mathbf b_t,\quad C_K=\sum_kw_kC(\mathcal S\mid\bar{\mathbf b}_k)',
 ['run_ocean_experiments.training_states','run_ocean_experiments.run'],
 ['Use only training-year displacements to fit 16 k-means++ representatives.','Run 60 centroid iterations with seed 6200; assign final nearest-centre labels.','Divide cluster counts by the training-date count to obtain weights.','Optimize the weighted conditional coverage with a matched budget.','Select by refined training coverage; evaluate all held-out dates without tuning.','Compare 16 and 32 representatives against full-year training evaluation.'], 'Centres: m; weights and coverage: dimensionless.')
entry(27,'Random design sample and coordinate plots','Main design sampling; Supplementary Analysis S3',
 r'c_x,c_y\sim\mathcal U(-160,160),\quad z_0\sim\mathcal U(120,270),\quad d=40\ \mathrm m',
 ['run_experiments.main'],
 ['With seed 43210, draw 600 independent feasible design vectors.','Construct the fixed-pitch grid and evaluate the deterministic objective.','Group by coordinate intervals for descriptive plots; do not infer significance from meshes.'], 'Coordinates: m; coverage: fraction.')

def source_index(identifier):
    module,function=identifier.split('.');path=CODE/(module+'.py')
    tree=ast.parse(path.read_text(encoding='utf-8'))
    node=next(n for n in tree.body if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef)) and n.name==function)
    return {'path':'code/'+path.name,'function':function,'start_line':node.lineno,'end_line':node.end_lineno}

def main():
    for item in entries:
        item['source'] = list(map(source_index, item['functions']))
    (ROOT / 'formula_registry.json').write_text(json.dumps(entries, indent=2), encoding='utf-8')
    lines = ['# Formula definitions and executable-code inventory', '',
             'F01--F27 identify computational definitions, assumptions, units and pseudocode. Source positions are generated from the executable Python syntax tree.', '']
    for item in entries:
        sources = '; '.join(f"{x['path']}:{x['start_line']} ({x['function']})" for x in item['source'])
        lines.extend([f"## {item['id']}: {item['title']}", '', '$$' + item['formula_latex'] + '$$', '',
                      f"Units: {item['units']}", '', f"Source: {sources}", '', item['status'], '',
                      *[f'{i + 1}. {step}' for i, step in enumerate(item['pseudocode'])], ''])
    (ROOT / 'FORMULAS_AND_PSEUDOCODE.md').write_text('\n'.join(lines), encoding='utf-8')
    print('Generated', len(entries), 'formula definitions, pseudocode and source mappings')

if __name__ == '__main__':
    main()
