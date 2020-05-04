# Overall description
This repository contains the code of machine learning model (Bayesian inference) that I use to gain insight from the data I collected on how humans make a sequence of visual judgments. This is basically my whole PhD thesis.

# Model fits the data well with a few interpretable and meaningful parameters
In this project, I collected, cleaned, analyzed and then used Bayesian inference to gain a deep understanding on how humans make sequential judgments of simple visual information. The full description and detail is in my published paper https://elifesciences.org/articles/33334. The code is in directory "BasicModel_eLife". Here I just sketch out some important points with figures taken from the paper.

## Experiment 1
Below is a simple layout of the experiment (panel b) in which human subjects had to make a sequence of 2 judgments: whether the mean orientation of the stimulus is clockwise or counterclockwise of a decision boundary (black lines); then they had to adjust a probe line to estimate the stimulus orientation. In panels c and d, I show the data for both tasks. Note that the distribution of estimation data in panel d shows an interesting bimodal pattern, which indicates repulsive biases away from the decision boundary.

![](/figures/fig1.png)

To gain insight on the decision strategy people used in this task, I employed a conditioned Bayesian observer model. The model can actually jointly account for data in both tasks with the same set of a few parameters. Notably, the parameters can be easily interepreted and intuitively meaningful (check out the paper if you are interested!).
![](/figures/fig2.png)

## Experiment 2, 3
In fact, the model can account for the data across a wide range of situations. Here's one more example when the same model can explain the data in 2 experiments jointly with the same parameters.
![](/figures/fig3.png)

In my paper, you can find several other situations that the model can explain the data well (both my data and other people's data).

# Model predicts data in different conditions
Ok now we see that the model can fit the data across many different situation with a few meaningful parameters. However, can the model generalize (i.e. predict) the data in a totally new conditions? Because I'm still working on the papers for those situations, I just show 2 examples here.

First, in this situation, we fit the model to data in correct trials and predict the data in incorrect trials. The  results show that the model can generalize well to another condition.
![](/figures/fig4.png)

Now, let's bring it up another level. I fit the model to the data in one experiment and try to predict the data in a totally different experiment. Again, the results demonstrate that the model can generalize very well to a totally different situation.
![](/figures/fig5.png)



