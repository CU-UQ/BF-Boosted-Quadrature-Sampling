# This script creates the plots in Figure 1 of our paper

from numpy import arange, corrcoef, sqrt, sum, transpose, zeros
from numpy.linalg import lstsq, svd, norm
from numpy.random import choice, normal, seed
import numpy as np
from scipy.linalg import null_space
import matplotlib.pyplot as plt

# Settings
N, d = 1000, 50  # Size of A
kappa_list = np.linspace(0.1,0.9,9)
phi_list = np.linspace(0.1,0.9,9)
nu_list = []
mu_dif_list = []
no_sketches = 10
m = 100  # Embedding dimension
seed(2)  # Set seed for reproducability

# Generate A, Q
A = normal(size=(N, d))
Q = svd(A, full_matrices=False)[0]

# Generate the b components
z1 = normal(size=(d, 1))
b1 = Q @ z1
b1 = b1 / norm(b1)
z2 = normal(size=(N, 1))
b2 = z2 - Q @ (transpose(Q) @ z2)
b2 = b2 / norm(b2)

# Generate the random component for bt (called \tilde{b} in the paper)
z3 = normal(size=(N, 1))

# Draw leverage score sampling sketches
# Storing these as dense matrices since they're small
p = sum(Q**2, axis=1) / d  # Leverage score sampling distribution

# Look through all pairs of kappa, phi and generate plots
for kappa in kappa_list:
    for phi in phi_list:
        # Draw all Gaussian sketches
        G = [normal(scale=1/sqrt(m), size=(m, N)) for k in range(no_sketches)]
        S = []
        rows = arange(m)
        for k in range(no_sketches):
            sketch = zeros(shape=(m, N))
            cols = choice(N, size=m, replace=True, p=p)
            vals = 1/sqrt(m * p[cols])
            sketch[rows, cols] = vals
            S.append(sketch)
        # Finish constructing b
        b = kappa * b1 + sqrt(1 - kappa**2) * b2

        # Finish constructing bt
        bt2 = z3 - b @ (transpose(b) @ z3)  # Component orthogonal to b
        bt2 = bt2 / norm(bt2)
        bt = phi * b + sqrt(1 - phi**2) * bt2

        # Check generated data
        print("Checking generated data for kappa = {:.4f}, phi = {:.4f}...".
              format(kappa, phi))
        kappa_emp = norm(Q @ transpose(Q) @ b) / norm(b)
        print("\tEmprirical kappa: {:.4f}".format(kappa_emp))
        phi_emp = float(transpose(b) @ bt)
        print("\tEmprirical phi: {:.4f}".format(phi_emp))
        Qper = null_space(Q.T)
        PQperb = Qper @ transpose(Qper) @ b 
        PQperbt = Qper @ transpose(Qper) @ bt
        nu_emp = float((transpose(PQperb) @ PQperbt)/(norm(PQperb)*norm(PQperbt)))
        print("\tEmprirical nu: {:.4f}".format(nu_emp))
        nu_list.append(nu_emp)

        # Compute all optimality coefficients for Gaussian sketch
        mu_Gaussian = zeros(shape=no_sketches)  # \mu^2(b, S)
        mu_Gaussian_t = zeros(shape=no_sketches)  # \mu^2(\tilde{b}, S)
        x_opt_b = lstsq(A, b, rcond=None)[0]
        x_opt_bt = lstsq(A, bt, rcond=None)[0]
        r = norm(A @ x_opt_b - b)
        rt = norm(A @ x_opt_bt - bt)
        new_mu_diff = np.empty(shape=2)
        min_mu_Gaussian_t = np.inf
        for sketch in range(no_sketches):
            x_hat_b = lstsq(G[sketch] @ A, G[sketch] @ b, rcond=None)[0]
            x_hat_bt = lstsq(G[sketch] @ A, G[sketch] @ bt, rcond=None)[0]
            rS = norm(A @ x_hat_b - b)
            rSt = norm(A @ x_hat_bt - bt)
            mu_Gaussian[sketch] = rS**2 / r**2 - 1
            mu_Gaussian_t[sketch] = rSt**2 / rt**2 - 1
            if mu_Gaussian_t[sketch] < min_mu_Gaussian_t:
                best_sketch = sketch
                min_mu_Gaussian_t = mu_Gaussian_t[best_sketch]
        new_mu_diff[0] = np.sqrt(mu_Gaussian[best_sketch]) - np.sqrt(np.min(mu_Gaussian))

        # Compute all optimality coefficients for leverage score sketch
        mu_levscore = zeros(shape=no_sketches)
        mu_levscore_t = zeros(shape=no_sketches)
        min_mu_levscore_t = np.inf
        for sketch in range(no_sketches):
            x_hat_b = lstsq(S[sketch] @ A, S[sketch] @ b, rcond=None)[0]
            x_hat_bt = lstsq(S[sketch] @ A, S[sketch] @ bt, rcond=None)[0]
            rS = norm(A @ x_hat_b - b)
            rSt = norm(A @ x_hat_bt - bt)
            mu_levscore[sketch] = rS**2 / r**2 - 1
            mu_levscore_t[sketch] = rSt**2 / rt**2 - 1
            if mu_levscore_t[sketch] < min_mu_levscore_t:
                best_sketch = sketch
                min_mu_levscore_t = mu_levscore_t[best_sketch]
        
        new_mu_diff[1] = np.sqrt(mu_levscore[best_sketch]) - np.sqrt(np.min(mu_levscore))
        mu_dif_list.append(new_mu_diff)

        print("\tOptimality coefficient correlation...")
        print("\t\tGaussian: {:.2f}".
              format(corrcoef(mu_Gaussian_t, mu_Gaussian)[0, 1]))
        print("\t\tLeverage score: {:.2f}".
              format(corrcoef(mu_levscore_t, mu_levscore)[0, 1]))

nu_list, mu_dif_list = np.array(nu_list),np.array(mu_dif_list)
plt.scatter(nu_list,mu_dif_list[:,0],label=r'Gaussian',alpha=0.7)
plt.scatter(nu_list,mu_dif_list[:,1],label=r'Leverage score',alpha=0.7)
plt.plot(np.sort(nu_list),0.6*np.sqrt(6)*np.sqrt(1.0-np.sort(nu_list)),c='g',label=r'$2\sqrt{6(1-\nu)\epsilon},\epsilon=0.09$')
plt.xlabel(r'$\nu$', fontsize=10)
plt.ylabel(r'$\mu(b,S_{\ell^*})-\mu(b,S_{\ell^{**}})$',fontsize=10)
plt.legend()
plt.savefig("synthetic_experiment_thm33.eps", bbox_inches="tight",format='eps')
plt.show()

print("Done")
