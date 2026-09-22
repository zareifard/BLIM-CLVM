# Bayesian Modeling Framework for ChIA-PET Data Analysis

This repository contains the implementation of the statistical framework proposed in the paper: **"A Gibbs Sampler for Detecting Significant chromatin Interactions from ChIA-PET data."**

## Abstract
Chromatin Interaction Analysis with Paired-End Tag (ChIA-PET) sequencing is a genome-wide, high-throughput technology for detecting chromatin interactions... [اینجا خلاصه مقاله را کپی کنید].

## Repository Structure
This repository provides implementations for two distinct Bayesian latent variable models:

*   **[BLIM (Binary Latent Indicator Model)](./BLIM/README.md):** A discrete classification approach to distinguish signal from background noise.
*   **[CLVM (Continuous Latent Variable Model)](./CLVM/README.md):** A flexible approach that captures heterogeneity in signal strength across interaction pairs.

Each model is implemented with efficient MCMC sampling schemes. Detailed instructions, model parameters, and usage examples are available within each respective folder.

---
*For questions or collaborations, please contact [Your Name/Email].*
