# Regression Analysis
This repository depends on https://github.com/PittSMLlab/labTools

# Abstract 

Humans can adapt their gait to compensate for changes in environmental demands and generalize learned movements from one situation to another.  One way to study the processes  that underlie locomotor adaptation is by exposing participants to split-belt treadmill walking; this paradigm induces robust changes in gait kinematics and muscle activation patterns. Additionally, we can study generalization by contrasting the adaptation effects (i.e., aftereffects) that participants exhibit in the same (treadmill) or different (overground) contexts from the adaptation. Here, we propose a novel approach for characterizing locomotor adaptation by analyzing individual muscle activation patterns. We ask what processes underlie the adaptation of muscle activation patterns and what aspects are observed during generalization (overground). 

We hypothesize that at least two processes with distinct dynamics underlie the changes in muscle activity during locomotor adaptation and generalization. Specifically, we posit that a fast reactive process will recruit a neuromuscular pattern to maintain balance at each transition between walking environments. A slow adaptive process will forge a contextual pattern that meets the demands of the novel split environment, but this pattern will be slowly disengaged during post-adaptation (tied walking) on the treadmill. It will not be used during overground post-adaptation due to environmental differences.  We recorded the activity of 28 leg muscles of twenty-four young adults (<40 yrs. old) who experienced split-belt walking during their adaptation, and their de-adaptation was measured on either the treadmill (n=12) or overground (n=12) walking. We used a data-driven approach to measure individual muscles’ reactive and contextual patterns to reproduce the evolution of muscle activity during the split-belt walking paradigm.

Our analysis showed that the reactive and contextual processes contribute to the adaptation and post-adaptation of muscle activity on the treadmill. However, during overground post-adaptation, 2 out of the 28 muscles generalized the contextual pattern, and all other muscles exhibited the reactive pattern, suggesting that the kinematic effects previously reported overground are mostly induced by reactive processes in response to a small number of muscles generalizing the split pattern. These findings provide insights into locomotor adaptation features beyond those drawn from traditional kinetic or kinematic analyses, which can be leveraged to study the effect of aging and brain lesions on the carryover of muscle activity.


# **How to process the data (step-by-step)**
This repository depends on labtools. Thus, the functions and scripts assume they work with an AdaptationData object. Ensure that your data is in the correct format. 

To view the quality of the data, you can check the EMGnorm of all the muscles or individual muscles. For this: 
* To plot the EMGnorm, please refer to AddingNormToAdaptData.m

Once you confirm that the quality of your data is good and select muscles that need to be removed (if necessary), you can start to estimate the dynamics of the regressors. 

1. To compute the EMG dynamics, you need:
2. Use preProcessingLinearModel.m to extract the data from the param file.
* Then you can estimate the dynamics using RegressionByMuscle.m
* To Plot the time courses for the reactive and contextual, please refer to ReactiveandContextualTimeCourses.m 
