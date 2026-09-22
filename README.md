# Bayesian Modeling Framework for ChIA-PET Data Analysis

This repository contains the implementation of the statistical framework proposed in the paper: **"A Gibbs Sampler for Detecting Significant chromatin Interactions from ChIA-PET data."**

## Abstract
Chromatin Interaction Analysis with Paired-End Tag (ChIA-PET) sequencing is a genome-wide‎, ‎high-throughput technology for detecting chromatin interactions associated with specific proteins of interest‎. ‎Despite its effectiveness‎, ‎ChIA-PET data often contain substantial background noise resulting from random ligation events‎, ‎which can obscure true chromatin interaction signals and lead to unreliable inferences‎. ‎To address this challenge‎, ‎we propose a Bayesian modeling framework for selecting significant chromatin interactions from candidate interactions provided by other ChIA-PET loop callers‎. ‎Specifically‎, ‎we introduce a latent variable modeling approach that distinguishes signal from noise through two complementary perspectives‎. ‎The proposed methods differ in how the transition between noise and signal interaction intensities is characterized‎. ‎The first approach‎, ‎referred to as the Binary Latent Indicator Model (BLIM)‎, ‎adopts a discrete classification of interactions‎, ‎whereas the second approach employs a more flexible Continuous Latent Variable Model (CLVM) to represent heterogeneity in signal strength across interaction pairs‎. ‎Parameter estimation is performed using an efficient Markov Chain Monte Carlo (MCMC) sampling scheme‎, ‎yielding posterior probabilities that quantify the strength and credibility of detected interactions. ‎We evaluate the performance of the proposed framework using both simulated and real ChIA-PET datasets‎, ‎GM12878 and HCT116‎, ‎and compare it with widely used tools‎: ‎ChIA-PIPE and ChIA-PET2‎. ‎The results demonstrate that our Bayesian approach consistently achieves high accuracy and sensitivity‎. ‎The high-probability inteactions identified by our method among the candidates selected from ChIA-PIPE and ChIA-PET2 outputs show improved quality in GM12878 and comparable quality in HCT116‎, ‎as indicated by aggregated peak analysis‎.
‎By integrating statistical rigor with biological interpretability‎, ‎the proposed framework provides a robust and extensible foundation for the analysis of high-throughput chromatin conformation data.

## Repository Structure
This repository provides implementations for two distinct Bayesian latent variable models:

*   **[BLIM (Binary Latent Indicator Model)](./BLIM/README.md):** A discrete classification approach to distinguish signal from background noise.
*   **[CLVM (Continuous Latent Variable Model)](./CLVM/README.md):** A flexible approach that captures heterogeneity in signal strength across interaction pairs.

Each model is implemented with efficient MCMC sampling schemes. Detailed instructions, model parameters, and usage examples are available within each respective folder.
