# SelfConsistentPerception
Matlab code demonstrates the self-consistent/conditioned Bayesian observer model for the experiments in my paper. 

Long Luu and Alan A Stocker (2018). Post-decision biases reveal a
self-consistency principle in perceptual inference. eLife.
https://elifesciences.org/articles/33334

## My experiments
The code is selfConsistentModel.m

For experiment 1 and 2, set flagDecisionGiven = 0.

For experiment 3, set flagDecisionGiven = 1.

## Zamboni et.al. experiments
The codes are condEstimMem_ShiftBoundary.m and condEstimMem_DisappearBoundary.m.

This is to deal with the situation when the decision boundary is shifted or disappear and I show that the conditioned Bayesian model can account well for both situations.

