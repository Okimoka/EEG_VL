=== Exemplary-task classification

#include "tables/eegnet_default_source.typ"

EEGNet discriminates the conditions in the exemplary tasks across held-out participants, with higher mean accuracy for oddball than gambling. The 95% bootstrap interval summarizes uncertainty around the mean score across participants. It is calculated by repeatedly resampling participants’ scores and taking the central 95% of the resulting averages. The three-run average for each participant is resampled; the networks are not refitted during the bootstrap.

=== Gameplay transfer

The table reports whole-epoch EEGNet predictions. Control means use the same participants as the corresponding event comparison; differences are expressed in percentage points (pp).

#include "tables/eegnet_default_whole_epoch_transfer.typ"

All three game events receive rare labels more frequently than random gameplay, although their mean rare-label fractions remain below 50%. This supports a relative shift towards the rare category, but does not reproduce the paper’s time-specific predominance of target classifications.

Missile hits also receive more win labels than random gameplay. The same direction is observed for enemy and wall crashes, limiting the interpretation of the gambling classifier as a selective measure of reward. The original paper similarly did not establish the predicted loss-like response to crashes.

The event-minus-control differences are positive for all six comparisons, and their descriptive 95% intervals exclude zero. These intervals summarize variation across participants conditional on the fitted models. Because the models share training participants, they should not be interpreted as independent confirmatory tests.

Using the complete 800 ms response does not localize the information responsible for these shifts. The shorter intervals below assess whether the pattern is also present in the response periods relevant to the original study.
