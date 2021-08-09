# Overall description
This repository contains the code of machine learning model (Bayesian inference) that I use to gain insight from the data I collected on how humans make a sequence of visual judgments. This is basically my whole PhD thesis.

# Model fits the data well with a few interpretable and meaningful parameters
In this project, I collected, cleaned, analyzed and then used Bayesian inference to gain a deep understanding on how humans make sequential judgments of simple visual information. The full description and detail is in my published paper https://elifesciences.org/articles/33334. The code is in directory "BasicModel_eLife". Here I just sketch out some important points with figures taken from the paper.

## Experiment 1
Below is a simple layout of the experiment (panel b) in which human subjects had to make a sequence of 2 judgments: whether the mean orientation of the stimulus is clockwise or counterclockwise of a decision boundary (black lines); then they had to adjust a probe line to estimate the stimulus orientation. In panels c and d, I show the data for both tasks. Note that the distribution of estimation data in panel d shows an interesting bimodal pattern, which indicates repulsive biases away from the decision boundary.

![](/figures/fig1.png)

To gain insight on the decision strategy people used in this task, I employed a conditioned Bayesian observer model. The key idea is that after making the first judgment, the observer considers that as correct and use it to constrain the subsequent judgment. This was implemented by a conditional prior that is consistent with the first judgment. The model can actually jointly account for data in both tasks with the same set of a few parameters. Notably, the parameters can be easily interepreted and intuitively meaningful (check out the paper if you are interested!).

![](/figures/fig2.png)

In fact, the model not only accounts for the data pooled from all subjects but can also account for the invidual differences across subjects. Below I show how the model can fit the data of individual subjects by meaningful variations of the parameters. For example, we can see that subjects S2 and S3 have big biases in their judgments and that is due to the high memory noise (i.e. they tend to quickly forget what they saw).

![](/figures/fig2_2.png)

## Experiment 2, 3
In fact, the model can account for the data across a wide range of situations. Here's one more example when the same model can explain the data in 2 experiments jointly with the same parameters.

![](/figures/fig3.png)

In my paper, you can find several other situations that the model can explain the data well (both my data and other people's data).

# Model predicts data in different conditions
Ok now we see that the model can fit the data across many different situation with a few meaningful parameters. However, can the model generalize (i.e. predict) the data in a totally new conditions? I just show 2 examples here but you can read my full published paper here https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008968.

## Fit correct trials, predict incorrect trials
First, in this situation, we fit the model to data in correct trials and predict the data in incorrect trials. The  results show that the model can generalize well to another condition.

![](/figures/fig4.png)

## Fit one experiment, predict another one
Now, let's bring it up another level. I fit the model to the data in one experiment (condition 1) and try to predict the data in a totally different experiment (condition 2). Again, the results demonstrate that the model can generalize very well to a totally different situation. More details can be found here https://github.com/longluu/High-to-low-visual-decoding.

![](/figures/fig5.png)

So, the power of the model is that it can generalize well to a different situation such as different stimuli, tasks, etc. Just to demonstrate that, I have the codes in the directories "MultiAttributeDecision" and "SunkCostFallacy" to illustrate how the model will predict human behaviors in totally different situations such as making judgments on different attributes/features or incorporating additional evidence. Note that those are purely prediction and there is no data yet but it can illustrate an important point that generative models like Bayesian model inherently possess the power to generalize very well to novel situations.

