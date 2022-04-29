# CAPRICA

Continuous Analysis and Probabilistic Inference for Cancer

Caprica is a Global Cumulative Treatment Analysis (GCTA) simulation toolkit.
It is used to model perpetual trials and evaluate active learning strategies.

## Tasks

1. Initialize
   * define the scope (e.g., treatment options, $A$, biomarkers, $X$)
   * fix true parameter values, $\theta_\mathrm{true}$ (effect sizes, noise, etc.)
   * define priors, $P(\theta)$, i.e., fix hyperparameters
   * define the patient selection distribution, $P(X)$
   * set decision policy, $\pi: \{\hat{y}_k\} \rightarrow a_k$
2. Predict 
   * inputs: 
	 * patient biomarkers, $x_i \in X$
	 * posterior distribution, $P(\theta \ | \ D)$ where $D = (\{y_j\}, \{x_j\}, \{a_j\})$
	 * treatment options, $A = \{a_1, a_2, ..., a_K\}$
   * output:
	 * $\{\hat{y}_{ik}\} \sim P(\hat{y}_{ik} \ | \ \theta, x_i, a_k)$
3. Decide
   * inputs:
	 * $K$ predictive distributions, $P(\hat{y}_{ik} \ | \ \theta, x_i, a_k)$
   * output:
	 * $a_i$
4. Observe
   * inputs:
	 * biomarkers, $x_i$
	 * chosen treatment, $a_i$
   * output:
	 * outcome, $y_i \sim P(y_i \ | \ a_i, x_i, \theta_\mathrm{true})$
5. Infer
   * inputs:
	 * data, $D$
   * outputs:
	 * $\theta \sim P(\theta \ | \ D)$
