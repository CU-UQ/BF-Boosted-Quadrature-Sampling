# Bi-fidelity Boosted Quadrature Sampling

This repo contains code for our paper 
> N. Cheng, O. A. Malik, Y. Xu, S. Becker, A. Doostan, A. Narayan.
> *Quadrature Sampling of Parametric Models with Bi-fidelity Boosting*.
> **arXiv preprint arXiv:2209.05705**,
> 2022.

The paper is available at [arXiv](https://arxiv.org/abs/2209.05705).

## Referencing this code

If you use this code in any of your own work, please reference our paper:
```
@misc{cheng2022quadrature,
      title={Quadrature Sampling of Parametric Models with Bi-fidelity Boosting}, 
      author={Nuojin Cheng and Osman Asif Malik and Yiming Xu and Stephen Becker and Alireza Doostan and Akil Narayan},
      year={2022},
      eprint={2209.05705},
      archivePrefix={arXiv},
      primaryClass={math.NA}
}
```

## Description of code

### Main files
- **cavity\_flow\_bifidelity\_HC.m:** Codes generate boosted quadrature sampling results for cavity flow data;

- **beam\_bifidelity\_HC:** Codes generate boosted quadrature sampling results for beam data;

### Sampling functions for <ins>hyperbolic cross</ins> polynomial spaces 

- **hyperbolic\_cross\_sampling.m:** This function samples and solve the least square problem with many different sampling options. The function utilizes the tensor structure and solves sampled least square problems with less time. The sampling methods provided are uniform sampling and leverage score sampling. 

- **det\_rejection\_sampling.m:** This function allows sampling according to leveraged volume sampling. This is done via the determinantal rejection sampling algorith in *[Derezinski et al., arXiv:1802.06749, 2018]*. Moreover, it is implemented to allow both boosted and unboosted sampling for both unstructured and Kronecker structured matrices.

### Help functions

- **RegVol.m:** This algorithm does volume sampling. It was presented in *[Derezinski and Warmuth, JMLR 19, 2018]*.

- **FastRegVol.m:** This algorithm is a faster version of RegVol.m, and uses RegVol.m as a substep. It was also presented in *[Derezinski and Warmuth, JMLR 19, 2018]*.

- **data_loader.m:** This function is used to load the various datasets.

- **sub\_tp\_idx\_set.m:** This function computes the sub-index sets corresponding to tensor product, total degree, and hyperbolic cross polynomial spaces. It can be expanded to handle other polynomial spaces if necessary.

- **my\_legendre\_1d.m:** This function generates 1d Gauss-Legendre node between -1 and 1;

- **mu.m:** This function computes optimality coefficient of the given sketching and QoI vector b.
