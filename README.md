# Code for Simulation-based validation of Bayes factor computation

Preprint: https://arxiv.org/abs/2508.11814

In R + Quarto, using renv to capture dependencies.

Once all libraries are installed, you can run `process.sh` to render all outputs.

Individual files describing the simulations are:

Toy examples:
- `binary_example.qmd`
- `poisson_nb_example.qmd`

Realistic examples:
- `turtles.qmd` - testing ranef presence in the turtles examples + discovering bad normalization constant with `bridgesampling`
- `lmbf_ranef_presence_post_sbc.qmd` - testing for ranef presence and using posterior SBC with the `BayesFactor` package
- `ttestBF.qmd` - JZS t-test with posterior SBC via the `BayesFactor` package
- `lmm_multiple.qmd` - comparing multiple linear mixed models with `brms` and `bridgesampling` as well as `brms::hypothesis`

Appendices:
-  `good_check_examples.qmd` examples showing problematic convergence for the Good check while showing that SBC and binary calibration do well
- `dap_metrics.qmd` comparing various tests for the data-averaged posterior criterion
- `large-linear-rjags.qmd` sampling a large model space directly with a Gibbs sampler


- `BF_SBC_paper.qmd` produces all of the final figures, assuming all of the above files were compiled and produce the necessary intermediate outputs
