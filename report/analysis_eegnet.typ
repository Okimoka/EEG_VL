== Replication objective

The original study examined whether EEG patterns associated with exemplary events could be identified during continuous and complex gameplay. To do this, the authors trained participant-specific classifiers on oddball and gambling trials, then transferred them to game events. Missile hits, enemy crashes and wall crashes showed a bias towards rare-oddball classifications at later post-event times (P3). Missile hits also resembled gambling wins. However, the expected loss-like response to crashes could not be established, and crashes were partly also classified as win-like. The authors interpreted these findings as evidence that signatures from the exemplary tasks associated with salience and reward generalize to gameplay. @cavanagh2016

Our objective is to see whether these relationships remain observable after our different preprocessing pipeline and a different analysis approach (using EEGNet). We therefore examine whether gameplay events are assigned the same exemplary categories as in the paper, and whether these classifications differ from those of ordinary gameplay.

== Analysis design

For each test participant, we train a classifier using exemplary trials from all other participants. We then apply this classifier to the test participant’s exemplary trials and gameplay events.
In the original paper, classifiers were trained and transfered within the same participant instead. This means in addition to testing the generalizability across tasks, we are now testing how well the training from other subjects generalizes to other subjects as well.

The training inputs consist of the complete 0–800 ms response and two shorter intervals: 300–600 ms for oddball and 200–350 ms for gambling. The shorter intervals focus on the P3 and feedback-related responses, which were the deciding intervals for the original paper.

=== Classifier choice

We use Braindecode 1.5.1 and its implementation of the EEGNet model with default parameters.
EEGNET is a compact convolutional network designed to learn spatial and temporal EEG features from voltage samples. We chose it due to it's ease of use and small size (which suits our relatively small amount of training data). We are only supplying the input dimensions, sampling frequency and number of output classes. Training settings are specified separately in the methods appendix. @eegnet @braindecode

=== Input data and participant inclusion

The cleaned dataset contains 2,786 oddball trials from 17 participants: 2,247 standard and 539 rare. Fourteen participants contribute 1,009 gambling trials: 494 losses and 515 wins. These labelled exemplary trials train the classifiers that are supposed to differentiate between oddball standards and oddball rares, or between gambling wins and gambling losses.

"Cleaned dataset" refers to the data that has been resampled to 250 Hz, filtered 0.1–40 Hz with −200 to 0 ms baseline. Bad channels have been interpolated to obtain common 63 EEG inputs, EOG is excluded.

== Whole-epoch classification and transfer: 0–800 ms

We use a training/test split of *16/1 participants for oddball* and *13/1 participants for gambling*. All trials from a participant stay in that person’s set. Each model is trained for *20 epochs with a batch size of 64*. The test participant remains excluded throughout. This procedure is repeated until each participant has served as the test participant once.

We run each training/test split three times with fixed random seeds. Performance and transfer fractions are calculated separately for each run and then averaged within each test participant. This reduces dependence on a single random training run.

EEG values are standardized using means and standard deviations calculated from the training data only. Class weights are inversely proportional to the number of training trials in each class. Detailed settings are provided in the methods appendix.

For the test participant’s exemplary trials, we compare the predicted labels with the known event labels. We calculate the percentage of correctly classified trials separately for each class, then average the two percentages. This is balanced accuracy: standard and rare trials, or wins and losses, contribute equally even when their trial counts differ. Random guessing has an expected balanced accuracy of 50%.

For gameplay events, there is no known exemplary-task label against which to measure accuracy. Instead, we report the percentage of each event type assigned to rare by the oddball classifier or win by the gambling classifier. A model score of at least 0.5 gives the rare or win label; a lower score gives standard or loss, respectively. We compare these percentages with the same classifier’s predictions on epochs sampled at random times during gameplay. These control epochs represent ordinary gameplay and may contain recorded events.

We calculate the percentages separately for each participant and then average them, giving every participant equal weight. Each event comparison requires at least 20 retained trials of that event type and 20 control epochs per participant. Only subject 007 falls below this minimum, with ten retained wall-crash trials, and is therefore excluded from the wall-crash comparison. Their other gameplay comparisons are retained.

#include "analysis_eegnet_results.typ"
