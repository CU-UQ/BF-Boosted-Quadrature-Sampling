# This script creates the plots in Figure 2 of our paper

from numpy import arange, corrcoef, sqrt, sum, transpose, zeros
from numpy.linalg import lstsq, svd, norm
from numpy.random import choice, normal, seed
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

# Settings
N, d = 1000, 50  # Size of A
kappa_list = [0.2, 0.95]  # Low and high kappa values
phi_list = [0.3, 0.95]  # Low and high phi values
no_sketches = 100
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

# Draw all Gaussian sketches
G = [normal(scale=1/sqrt(m), size=(m, N)) for k in range(no_sketches)]

# Draw leverage score sampling sketches
# Storing these as dense matrices since they're small
p = sum(Q**2, axis=1) / d  # Leverage score sampling distribution
S = []
rows = arange(m)
for k in range(no_sketches):
    sketch = zeros(shape=(m, N))
    cols = choice(N, size=m, replace=True, p=p)
    vals = 1/sqrt(m * p[cols])
    sketch[rows, cols] = vals
    S.append(sketch)

# Create plot window
fig, axs = plt.subplots(2, 4, figsize=(12, 5))
plt_col = 0

# Look through all pairs of kappa, phi and generate plots
for kappa in kappa_list:
    for phi in phi_list:
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

        # Compute all optimality coefficients for Gaussian sketch
        mu_Gaussian = zeros(shape=no_sketches)  # \mu^2(b, S)
        mu_Gaussian_t = zeros(shape=no_sketches)  # \mu^2(\tilde{b}, S)
        x_opt_b = lstsq(A, b, rcond=None)[0]
        x_opt_bt = lstsq(A, bt, rcond=None)[0]
        r = norm(A @ x_opt_b - b)
        rt = norm(A @ x_opt_bt - bt)
        for sketch in range(no_sketches):
            x_hat_b = lstsq(G[sketch] @ A, G[sketch] @ b, rcond=None)[0]
            x_hat_bt = lstsq(G[sketch] @ A, G[sketch] @ bt, rcond=None)[0]
            rS = norm(A @ x_hat_b - b)
            rSt = norm(A @ x_hat_bt - bt)
            mu_Gaussian[sketch] = rS**2 / r**2 - 1
            mu_Gaussian_t[sketch] = rSt**2 / rt**2 - 1

        # Compute all optimality coefficients for leverage score sketch
        mu_levscore = zeros(shape=no_sketches)
        mu_levscore_t = zeros(shape=no_sketches)
        for sketch in range(no_sketches):
            x_hat_b = lstsq(S[sketch] @ A, S[sketch] @ b, rcond=None)[0]
            x_hat_bt = lstsq(S[sketch] @ A, S[sketch] @ bt, rcond=None)[0]
            rS = norm(A @ x_hat_b - b)
            rSt = norm(A @ x_hat_bt - bt)
            mu_levscore[sketch] = rS**2 / r**2 - 1
            mu_levscore_t[sketch] = rSt**2 / rt**2 - 1

        axs[0, plt_col].scatter(mu_Gaussian_t, mu_Gaussian, color="tab:red")
        axs[0, plt_col].title.set_text(r"$\kappa = {:.2f}$, $\varphi = {:.2f}$"
                                       .format(kappa, phi))
        axs[1, plt_col].scatter(mu_levscore_t, mu_levscore, color="tab:blue")
        plt_col += 1

        print("\tOptimality coefficient correlation...")
        print("\t\tGaussian: {:.2f}".
              format(corrcoef(mu_Gaussian_t, mu_Gaussian)[0, 1]))
        print("\t\tLeverage score: {:.2f}".
              format(corrcoef(mu_levscore_t, mu_levscore)[0, 1]))

plt.setp(axs[-1, :], xlabel=r"Low-fidelity $\mu^2$")
plt.setp(axs[:, 0], ylabel=r"High-fidelity $\mu^2$")

for ax in axs.flat:
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))

plt.savefig("synthetic_experiment.pdf", bbox_inches="tight")
plt.show()

print("Done")