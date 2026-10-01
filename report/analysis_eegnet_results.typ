=== Exemplary-task classification

#include "tables/eegnet_default_source.typ"

EEGNet discriminates the conditions in the exemplary tasks for subjects whose data were not used to train the classifier, with higher mean accuracy for oddball than gambling.

=== Gameplay transfer

The table reports the classifier results on gameplay events. For example, for "Rare" and "Missile hit", we have 27.1% of missile-hit trials classified as rare, whereas 17% of the random control trials were classified as rare, which gives a difference of 9.8%.

#include "tables/eegnet_default_whole_epoch_transfer.typ"

All three game events receive rare labels more frequently than random gameplay, although their mean rare-label fractions remain below 50%. This means the classifier seems to have a shift towards the rare category.

Missile hits also receive more win labels than random gameplay. The same direction is observed for enemy and wall crashes, which means we cannot interpret the gambling classifier as a predictor of reward. The original paper similarly could not make this inference.

The mean event-minus-control differences are positive for all six comparisons. These are descriptive results from the fitted models.

The original paper specifically identified specific shorter intervals of the epochs to be responsible for the decisions, rather than the complete 800 ms response. In the following we will again train classifiers, using only the time intervals isolated in the original study.