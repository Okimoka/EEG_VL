== Comparison in the P3 and feedback intervals

We now train EEGNet using only 300–600 ms for oddball and 200–350 ms for gambling. These intervals focus on the P3 and feedback-related responses. The classifiers use the same default model settings, participant assignments and training procedure as the whole-epoch analysis.

=== Gameplay transfer

Within these intervals, EEGNet reaches 74.9% balanced accuracy for oddball and 60.3% for gambling, which is ~5% worse for both classifiers.

The table reports transfer results for these new classifiers

#include "tables/eegnet_default_fixed_window_transfer.typ"

For oddball, all three events shift towards rare relative to random gameplay, but their mean rare-label fractions remain below 50%. We can therefore again not interpret an oddball rare as a "target-like" result.
For gambling, all of the events now receive predominantly win labels, which matches the papers ambiguous result here as well.


