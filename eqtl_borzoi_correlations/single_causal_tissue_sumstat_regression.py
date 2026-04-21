import numpy as np


class SINGLE_CAUSAL_TISSUE_SUMSTAT_REGRESSION(object):
	def __init__(self):
		self.best_tissue_index = None
		self.beta_mu = None
		self.beta_cov = None
		self.component_variances = None
		self.objective = None
		self.objective_per_tissue = None

	def fit(self, tau, tau_cov):
		tau = np.asarray(tau)
		tau_cov = np.asarray(tau_cov)

		if tau.ndim != 1:
			raise ValueError('tau must be a 1D vector')
		if tau_cov.ndim != 2:
			raise ValueError('tau_cov must be a 2D matrix')
		if tau_cov.shape[0] != tau_cov.shape[1]:
			raise ValueError('tau_cov must be square')
		if tau_cov.shape[0] != len(tau):
			raise ValueError('tau and tau_cov have incompatible shapes')

		n_tissues = len(tau)
		objective_per_tissue = np.full(n_tissues, np.nan)
		posterior_means = np.zeros(n_tissues)
		posterior_vars = np.zeros(n_tissues)

		for tissue_iter in range(n_tissues):
			tau_var = tau_cov[tissue_iter, tissue_iter]
			if tau_var <= 0.0:
				raise ValueError('tau_cov diagonal entries must be positive')

			# Single-tissue Gaussian regression using only the marginal mean/variance
			posterior_means[tissue_iter] = tau[tissue_iter]
			posterior_vars[tissue_iter] = tau_var
			objective_per_tissue[tissue_iter] = (tau[tissue_iter] ** 2.0) / tau_var

		self.best_tissue_index = int(np.argmax(objective_per_tissue))
		self.objective_per_tissue = objective_per_tissue
		self.objective = objective_per_tissue[self.best_tissue_index]

		self.beta_mu = np.zeros(n_tissues)
		self.beta_mu[self.best_tissue_index] = posterior_means[self.best_tissue_index]

		self.beta_cov = np.zeros((n_tissues, n_tissues))
		self.beta_cov[self.best_tissue_index, self.best_tissue_index] = posterior_vars[self.best_tissue_index]

		self.component_variances = np.zeros(n_tissues)
		self.component_variances[self.best_tissue_index] = (
			posterior_means[self.best_tissue_index] ** 2.0 +
			posterior_vars[self.best_tissue_index]
		)

		return self
