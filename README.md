# Regression Analysis
This repository depends on https://github.com/PittSMLlab/labTools

# Abstract 

Humans can adapt their gait to compensate for changes in environmental demands and generalize learned movements from one situation to another.  One way to study the processes  that underlie locomotor adaptation is by exposing participants to split-belt treadmill walking; this paradigm induces robust changes in gait kinematics and muscle activation patterns. Additionally, we can study generalization by contrasting the adaptation effects (i.e., aftereffects) that participants exhibit in the same (treadmill) or different (overground) contexts from the adaptation. Here, we propose a novel approach for characterizing locomotor adaptation by analyzing individual muscle activation patterns. We ask what processes underly the adaptation of muscle activation patterns and what aspects are observed during generalization (overground). 

We hypothesize that at least two processes with distinct dynamics underlie the changes in muscle activity during locomotor adaptation and generalization. Specifically, we posit that a fast reactive process will recruit a neuromuscular pattern to maintain balance at each transition between walking environments. Whereas a slow adaptive process will forge a contextual pattern meeting the demands of the novel split environment, this pattern will be slowly disengaged during post-adaptation (tied walking) on the treadmill. It will not be used during overground post-adaptation due to environmental differences.  We recorded the activity of 28 leg muscles of twenty-four young adults (<40 yrs. old) who experienced split-belt walking during their adaptation, and their de-adaptation was measured on either the treadmill (n=12) or overground (n=12) walking. We used a data-driven approach to measure individual muscles’ reactive and contextual patterns to reproduce the evolution of muscle activity during the split-belt walking paradigm.

Our analysis showed that the reactive and contextual processes contribute to the adaptation and post-adaptation of muscle activity on the treadmill. However, during overground post-adaptation, 2 out of the 28 muscles generalized the contextual pattern, and all other muscles exhibited the reactive pattern, suggesting that the kinematic effects previously reported overground are mostly induced by reactive processes in response to a small number of muscles generalizing the split pattern. These findings provide insights into locomotor adaptation features beyond those drawn from traditional kinetic or kinematic analyses, which can be leveraged to study the effect of aging and brain lesions on the carryover of muscle activity.


**How to process the data (step-by-step)**
1. To plot the EMGnorm please refer to AddingNormToAdaptData.m
2. To compute the EMG dynamics, you need:
3. Use preProcessingLinearModel.m to extract the data from the param files 
  3.a Then you can estimate the dynamics using model_fit_individual_muscles_PerSubject.m
  3.b To Plot the time courses for the reactive and contextual, please refer to ReactiveandContextualTimeCourses.m 
