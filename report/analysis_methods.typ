#pagebreak(weak: true)
= Appendix: classifier settings and checks

== Inputs and gameplay controls

Analysis copies of the native repaired epochs retain the 0.1–40 Hz band, −200 to 0 ms baseline, reference and local trial repairs. Globally bad channels are spherically interpolated to give 63 common scalp inputs; EOG is omitted. Polyphase resampling gives 250 Hz, with 200 samples from 0 to 796 ms for the whole-epoch analysis. Neither the +40 ms event correction nor the 500 ms repeat-event rule is applied again.

Up to 400 control anchors per participant are sampled on randomly phased 2 s grids within event-supported gameplay spans; gaps above 30 s split spans. Controls may overlap recorded events and represent ongoing gameplay. AutoReject estimators recovered from the ICA-only short epochs reproduce the original rejection masks and checked repaired voltages. The saved estimators repair or reject controls before the same baseline, interpolation and resampling. This retains 6,049 of 6,800 candidates. Estimators, repair logs and trial IDs are retained.

Every source model includes all eligible participants with both training labels. Gameplay contrasts require at least 20 retained event and control trials per participant. Subject 007 has only ten wall crashes, giving 16 oddball and 13 gambling participants for wall-crash comparisons.

== EEGNet architecture and training

We instantiate Braindecode 1.5.1 EEGNet as `EEGNet(n_chans=63, n_outputs=2, n_times=n_times, sfreq=250)`. All optional model arguments use their defaults: F1 = 8 temporal filters, depth multiplier D = 2, F2 = 16 separable filters, a 64-sample temporal kernel and dropout 0.25. The whole-epoch network has 2,306 trainable parameters. These are implementation defaults, not parameters shown to be optimal for this dataset. @eegnet @braindecode

Training uses Braindecode’s `EEGClassifier` with Adam, learning rate 0.001, batch size 64 and 20 training epochs. Cross-entropy loss uses inverse-frequency class weights calculated from the training labels. Channel means and standard deviations come only from the training data. A model tested on a participant receives none of that person’s trials for fitting, scaling or batch-normalization updates.

Each participant serves once as the test participant while all other eligible participants form the training set. We repeat each split with three fixed seeds, which control model initialization and training randomness. Scores and label fractions are calculated for each fitted model and then averaged within participants. We do not combine model scores before assigning labels. The seed for every model is saved with its weights.

For the shorter intervals, inputs contain 75 samples from 300–600 ms for oddball and 38 samples from 200–350 ms for gambling. Endpoints are half-open on the 250 Hz grid. Default architecture settings are retained, giving 2,178 and 2,146 parameters. The 64-sample kernel exceeds the gambling input duration; the implementation handles this with zero padding, without receiving EEG from outside the cropped interval. Earlier filtering and baseline correction can still affect temporal interpretation.

The shorter-interval analysis includes three contrasts per task: the direct distinction between exemplary categories, its positive category versus ordinary gameplay, and its negative category versus ordinary gameplay. Positive means rare or win; negative means standard or loss. Classifiers comparing exemplary tasks with gameplay use random controls only from their training participants. The held-out random controls also appear in their source evaluation and transfer summaries, so these are not independent pieces of evidence.

== Summaries and reproducibility

Exemplary-task summaries first average the three runs within each participant, then average participant balanced accuracies. Gameplay summaries use the same order of averaging. Gameplay labels use a fixed 0.5 threshold. Each person’s event fraction is paired with their ordinary-gameplay fraction before averaging; control means use the same participants as the event comparison. Intervals use 10,000 participant bootstrap resamples and describe uncertainty conditional on fitted models. They are not simultaneous intervals, do not refit networks and do not account for all dependence arising from shared training participants. No formal confirmatory significance claim is made.

The analysis comprises 31 whole-epoch and 93 shorter-interval train/test splits, each fitted three times, giving 372 networks. Saved-model checks independently reproduce every exemplary-task and gameplay prediction and verify trial identities, participant separation, default architecture, training-only scaling, class weights and training duration. These checks establish implementation consistency, not psychological specificity.

The analysis uses three scripts: `prepare.py` constructs EEG arrays and random gameplay controls from the post-ICA files; `train.py` fits the networks; and `results.py` creates the summaries, report tables and transfer figure. Settings are specified directly in these scripts. Generated EEG arrays, cleaning records, model weights and numeric results are stored outside the submission in `artifacts/analysis_final/`. The submission contains the scripts and the tables and figure used in the report.
