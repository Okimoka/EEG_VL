== Comparison in the P3 and feedback intervals

We also train EEGNet using only 300–600 ms for oddball and 200–350 ms for gambling. These intervals focus on the P3 and feedback-related responses. The classifiers use the same default model settings, participant assignments and training procedure as the whole-epoch analysis. This comparison is exploratory.

The original paper used short moving windows and participant-specific classifiers. Our broader, fixed intervals and cross-participant training address the general transfer pattern without reproducing the original time courses. We examine both the proportion of rare or win labels and its difference from random gameplay.

=== Exemplary-task classification and gameplay transfer

Within these intervals, EEGNet reaches 74.9% balanced accuracy for oddball and 60.3% for gambling. The descriptive 95% intervals are 70.1–79.3% and 56.7–64.1%, respectively.

The table reports transfer from the direct rare/standard and win/loss classifiers. Event and control means use matched participant sets.

#include "tables/eegnet_default_fixed_window_transfer.typ"

For oddball, all three events shift towards rare relative to random gameplay, but their mean rare-label fractions remain below 50%. The predominantly target-like result is therefore not recovered by the cross-participant models in the selected 300–600 ms interval. For gambling, hits receive predominantly win labels at 200–350 ms, consistent with the direction of the paper’s direct win-versus-loss result. Both crash conditions also receive predominantly win labels.

Missile hits receive 61.5% win labels, with a descriptive 95% interval of 55.4–67.8%. The difference from controls is +16.5 percentage points, with an interval of +10.2 to +22.5 points. As with the whole-epoch analysis, all six event-minus-control intervals exclude zero. These are descriptive summaries of the fitted models.

#figure(image("../figures/eegnet_default_transfer.png", width: 100%),
  caption: [EEGNet transfer using the whole response (top) and the shorter response intervals (bottom). Points show the mean within-participant difference between event and control label fractions after averaging three training runs; bars are descriptive 95% participant bootstrap intervals. Positive values indicate more rare or win labels than for random gameplay.],
)

=== Exemplary tasks versus gameplay

The paper also compared each exemplary category with random gameplay. We therefore train separate EEGNet classifiers for rare versus gameplay, standard versus gameplay, win versus gameplay and loss versus gameplay in the same shorter intervals. Training uses exemplary trials and random gameplay epochs from the training participants only; the test participant remains excluded.

The table shows held-out classification performance for these contrasts and the fraction of game-event trials assigned the exemplary label. Each row represents a separate binary classifier, so the label fractions across rows are not complementary probabilities.

#include "tables/eegnet_default_background.typ"

Although these classifiers discriminate exemplary trials from gameplay, all mean exemplary-label fractions for the game events remain below 50%. Missile hits receive 26.6% win labels in the win-versus-gameplay model, compared with 14.4% for random gameplay. The loss-versus-gameplay model gives hits 19.2% loss labels, compared with 13.0% for controls. Both exemplary-label fractions therefore increase relative to controls, while gameplay remains the majority classification.

Agreement with the original hit result is consequently limited to the direct win-versus-loss contrast in this interval. The additional models do not recover a predominance of exemplary classifications over gameplay. They may also exploit differences between the exemplary-task and gameplay contexts rather than a specific event response.
