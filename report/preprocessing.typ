// Preprocessing snapshot; analysis chapter will be added to the final report separately.
#set document(
  title: "Project Report: Identification of canonical neural events during continuous gameplay of an 8-bit style video game",
  author: ("Sana Hafeez", "Hadi Ismail", "Okan Mazlum"),
)

#title()

Sana Hafeez, Hadi Ismail, Okan Mazlum #h(1fr) SS 2026

= Introduction

This report documents an attempt to reproduce the main findings of _Identification of canonical neural events during continuous gameplay of an 8-bit style video game_ by Cavanagh and Castellanos (2016), for the module _Signal processing and Analysis of human brain potentials (EEG)_. @cavanagh2016 We use the study’s publicly shared EEG dataset, OpenNeuro ds003517. @openneuro The accompanying project repository provides the code and reproduction instructions. @repository

Our main question is whether the relationship between laboratory-task responses and gameplay events remains observable with a different preprocessing pipeline. Following the workbench guidance, we compare the original and implemented methods, justify differences, inspect intermediate signals and document both successful and unsuccessful steps. This tests the robustness of the findings to reasonable analysis choices while making the workflow reproducible. @wbs-replication The analysis also extends the original approach with pooled classifiers evaluated on held-out subjects; differences from the paper therefore cannot be attributed to preprocessing alone.

== Original paper

The paper addresses a problem of generalizability: EEG findings established with isolated stimuli in carefully controlled tasks can be difficult to interpret in a complex, continuous environment. The authors developed the space-shooter game _Escape from Asteroid Axon_ to study goal-oriented events while retaining precise event timing. Before playing, subjects completed two established exemplar tasks: oddball target detection and gambling wins/losses. The authors compared the EEG responses statistically and trained subject-specific classifiers on the exemplar tasks, then transferred their learned patterns to gameplay. @cavanagh2016

The main finding was that wall crashes, enemy crashes and missile hits showed a transfer bias toward oddball-target patterns during the P3 period. Missile hits also resembled gambling wins, whereas crashes did not reliably resemble gambling losses. The authors interpreted these results as evidence that some laboratory signatures associated with salience and reward generalize to gameplay; the expected punishment-related result was not established. @cavanagh2016

== Dataset description

The dataset contains 17 subjects (11 men and 6 women; ages 18–39, mean 20.94 years), reported as right handed. Each subject has two runs under the BIDS task name `ContinuousVideoGamePlay`. Run 01 contains oddball and gambling; run 02 contains gameplay, with the exception described below. Full recording durations average 777.78 s and 1,976.49 s respectively and include periods outside the active tasks. EEG was sampled at 500 Hz with CPz as acquisition reference. There are 63 recorded EEG channels plus bipolar VEOG and HEOG; CPz is not independently recorded. The supplied electrode positions match the BESA `standard-10-5-cap385.elp` lookup montage and are shared template positions, rather than individual head measurements. @openneuro

=== Oddball

Participants mentally counted rare enemy-spacecraft images without a button response. The intended design contained 160 standard and 40 rare images. Each image was shown for 500 ms, followed by a random 1–2 s interval. Rare events were intended to occur after at least three but before ten standards. @cavanagh2016 In the supplied EEG event tables, the interval from the first to the last oddball-task event averages 6.69 min across the 17 subjects. This measures the recorded event span, excluding instructions and time outside those events, and closely agrees with the paper’s 6.68 min.

=== Gambling

Participants chose one of two doors using the left/right gamepad controls. Each trial resulted in a win of 75 credits or a loss of 60 credits; the intended design contained 40 wins and 40 losses. The doors remained on screen for up to 4 s or until a choice. Feedback was then presented immediately for 500 ms, followed by a 1 s interval. The paper gives a mean task duration of 2.86 min. @cavanagh2016 In the available EEG files, the first-to-last outcome-event span averages 3.18 min across the 14 subjects with gambling events, including the unusual trial counts described below. The full task duration and the span between recorded outcomes are different measures.

=== Escape from Asteroid Axon

Participants steered a ship and fired missiles, avoided wall/enemy crashes, collected stars to restore health and collected ammunition. An empty health bar ended the round; another round could be started by a button press. Game events include shots, star collections, ammunition collection, wall crashes, enemy crashes and missile hits. Gameplay was intended to last 30 min. @cavanagh2016

=== Dataset exceptions and recording boundaries

- All 17 subjects have 200 oddball events, but some early recordings contain 37–44 rare events rather than exactly 40.
- Subjects 002, 003 and 005 have no gambling-outcome events, consistent with the three missing gambling recordings reported in the paper.
- Subject 001 has only 18 win and 17 loss events; subject 006 has 120 wins and 120 losses. The remaining subjects with gambling events have the intended 40 wins and 40 losses.
- Subject 001 has the extended gameplay session: approximately 45.9 min through the last logged game event. The paper describes one 45 min session and reports an overall mean of 31 min. @cavanagh2016 The behavioural logs can continue after gameplay ends, so their final timestamp is not automatically the task duration.
- For the other 15 subjects, the first game-round marker occurs a median of 90.44 s after the beginning of run 02. Alignment with the behavioural logs shows that approximately 40.9 s of initial gameplay is absent for subject 003. For subject 008, approximately 11.2 s of initial gameplay is recorded at the end of run 01, followed by a gap of approximately 34.8 s before run 02. Thus run 02 begins about 46.0 s into that game session; it is not simply missing its first 35 s.

The subject-008 run-01 tail contains eight shots, two enemy crashes, one star and one game-over-labelled event. No oddball or gambling events occur in run 02. The current analyses retain eligible gameplay events from both runs; the tail has not been cropped or silently excluded. The recording/log alignment and duration calculations are saved in `reference/report_notes_revision/`.

= Preprocessing

We use MNE-BIDS-Pipeline because it is the workflow used in the course; the original authors used EEGLAB, a MATLAB toolbox. Differences between these software ecosystems affect many implementation details and defaults. This chapter focuses on the main methodological choices, their reasons and the checks used to assess them. It does not attempt an exhaustive comparison of every solver or interpolation implementation.

== Initial loading attempt

=== Channel metadata

With a minimal configuration, loading the original dataset fails during the pipeline’s initial data-quality report: `picks (None, treated as "data") yielded no channels, consider passing picks explicitly`. All entries in the original `*_channels.tsv` type and units columns are `n/a`, so MNE-BIDS imports the channels as miscellaneous and finds no EEG data channels for the spectrum calculation.

This motivated `scripts/prepare_bids.py`, which creates a separate prepared dataset with corrected metadata. The 63 scalp channels receive type `EEG`, VEOG and HEOG receive type `EOG`, and all receive units `uV`, following the EEGLAB storage convention. The corresponding `*_eeg.json` counts change from EEG/EOG/miscellaneous = 64/0/65 to 63/2/0. Metadata files are copied independently; the original signal files are linked without changing their samples. Assigning units does not resolve the auxiliary channels’ uncertain physiological calibration.

=== Electrode coordinates

With channel types corrected, another warning appears: “Other is not an MNE-Python coordinate frame for EEG data and so will be set to 'unknown'”. The original `*_coordsystem.json` declares `EEGCoordinateSystem = "Other"`, which BIDS permits when adequately defined, but its description does not establish the correct frame here. @bids-eeg

The supplied MATLAB scripts, `Convert2BIDS_ExAAx.m` and `ExAAx_Preprocess.m`, identify `standard_BESA/standard-10-5-cap385.elp`. We therefore extend preparation to correct `*_electrodes.tsv` and `*_coordsystem.json`, retaining the supplied positions and transforming them with the matching template landmarks into MNE head coordinates. The corrected coordinates use metres and the `CapTrak` coordinate convention; this describes the coordinate frame, not subject-specific CapTrak digitization. Details remain in `docs/montage_comparison.md`.

=== First inspection and visual delay

The prepared dataset allows an initial pipeline execution and inspection of minimally filtered EEG. The rare-oddball overview below shows an early posterior positivity near 0.138 s, a broader positivity near 0.410 s and a much larger frontal deflection near 0.870 s, consistent with visual P1, P3/P300 and ocular activity respectively.

#figure(
  image("figures/initial_oddball_rare.png", width: 100%),
  caption: [Initial rare-oddball butterfly plot with minimal filtering, before the visual-delay correction. The marked times are 0.138, 0.410 and 0.870 s.],
)

The early peak is consistent with a visual P1 near 100 ms after allowing for the paper’s independently measured 40 ms display delay. @cavanagh2016 We shift the ten analysis-event types *40 ms forward* in the prepared `*_events.tsv` files. This moves the 138 ms EEG peak to *98 ms relative to the corrected event*, without moving the EEG samples. Acquisition/status markers are unchanged. Peak timing alone does not establish the delay, and P3 or blink timing need not match one fixed latency. The installed pipeline has no dedicated visual-delay setting, so the correction is recorded once in preparation rather than repeated in downstream analyses.

== Pipeline configuration

The table summarizes the main methodological differences. Event-delay and repeated-event rules follow the paper and are described in the second-pass section rather than listed as differences.

#table(columns: (22%, 37%, 41%),
  table.header([*Aspect*], [*Original paper*], [*This project*]),
  [Bad-channel selection], [Four ventral electrodes excluded, then FASTER], [PyPREP including RANSAC; all threshold-based flags accepted after inspection; no default ventral exclusions],
  [Frequency filtering], [0.01–100 Hz acquisition; 0.1–20 Hz ERP/classification], [0.1–100 Hz continuous EEG with a 60 Hz notch; ICA fit at 1–100 Hz; ERP/classification at 0.1–40 Hz],
  [Rereferencing], [Initial average reference and CPz reconstruction; subsequent bad-channel interpolation], [Average reference over the common non-bad EEG set; no CPz reconstruction],
  [Combining runs], [Supplied scripts suggest separate exemplar/gameplay decompositions; pooling is not explicit in the paper], [One ICA per subject across both runs; common union of bad-channel flags],
  [ICA training], [Methods describe −2…+2 s event epochs before ICA], [Evenly spaced 4 s windows, downsampled to 250 Hz, without baseline or channel interpolation before fitting],
  [ICA component selection], [VEOG/HEOG correlations with manual verification], [ICLabel proposals plus visual review; automatic EOG-based selection disabled],
  [Bad-epoch handling], [FASTER, channel interpolation and epoch rejection], [AutoReject selects ICA fitting windows; a fresh post-ICA AutoReject fit repairs channels and rejects analysis trials],
)

=== Bad-channel selection

We chose PyPREP because it works with MNE data and provides several complementary channel-quality checks. We use the installed `NoisyChannels` detector suite with its standard thresholds, including NaN/flat, deviation, high-frequency noise, correlation/dropout, SNR, PSD and RANSAC checks. @pyprep In exploratory runs, the flagged channels varied substantially. We therefore fixed the random seed at *2026* and increased RANSAC’s `n_samples` to *1,000 random channel subsets* to improve reproducibility and sampling stability. All non-RANSAC checks run first; RANSAC is then called explicitly with that subset setting.

To inspect the flags, we built a small tool, `./review-channels` in `01_review_bad_channels/`, using MNE’s data browser. It shows the flagged channel alongside nearby electrodes in original or frequency-filtered EEG, together with detector scores. Correlation/dropout and RANSAC provide one-second and five-second window evidence respectively, allowing affected intervals to be inspected. Whole-record amplitude or spectral flags need not have one unique offending interval.

Our visual inspection led to three observations:

- Many flagged channels were visibly noisier than nearby electrodes.
- Some flags were difficult to verify visually, especially RANSAC flags: prediction from multiple electrodes is not necessarily apparent in a simple trace comparison.
- Some unflagged channels looked problematic but had detector scores just below the rejection threshold.

// Replace these three SVG placeholders with the selected annotated EEG screenshots.
// Keep the three image calls in this grid to preserve the side-by-side arrangement.
#figure(
  grid(
    columns: (1fr, 1fr, 1fr),
    gutter: 1em,
    image("figures/channel-example-a.svg", width: 100%),
    image("figures/channel-example-b.svg", width: 100%),
    image("figures/channel-example-c.svg", width: 100%),
  ),
  caption: [Placeholders for the three review examples: (A) an obviously noisy flagged channel compared with its neighbours; (B) a flag that is difficult to explain visually, particularly from RANSAC; (C) a suspicious channel just below its detector threshold. Annotated EEG screenshots will replace these panels.],
)

Given the project’s scope and our limited experience with manual channel rejection, we adopted all PyPREP threshold-based decisions without manual overrides. There are *68 flagged channel–recording cases*, not 68 distinct electrodes. The median is two flags per recording (range 0–6); six recordings have none. CP2 is flagged most often, in seven recordings. The paper also reports a median of two interpolated electrodes, although that is not an identical selection measure. @cavanagh2016

#figure(
  image("figures/bad_channel_topomap.svg", width: 100%),
  caption: [Frequency of PyPREP flags over all 34 recordings and separately by run. All panels use the same count scale. Labels identify electrodes flagged in at least five recordings overall.],
)

The paper excluded FT9, TP9, TP10 and FT10 by default, corresponding to 100% exclusion at those sites. PyPREP flags 21/136 ventral channel–recording cases (15.4%), compared with 47/2,006 (2.3%) at the other 59 electrodes. We apply the same automatic criteria to both groups.

PyPREP runs in a *separate script*, `scripts/run_pyprep.py`, after the initial preparation and preliminary frequency-filtering pass, before the main ICA pipeline execution. It is not an integrated pipeline detector. The accepted flags are then written into the analysis-ready prepared sidecars. Because joint ICA fitting requires identical bad-channel lists in both runs, preparation uses their union for each subject and preserves the original run-specific flags in an additional column. Thus the actual sequence is metadata preparation → preliminary filtering → standalone PyPREP → preparation of the shared bad-channel policy → main ICA pipeline.

=== Frequency filtering

Continuous EEG is filtered to 0.1–100 Hz with a 60 Hz notch at the original 500 Hz sampling rate. The ERP/classification branch uses a continuously filtered 0.1–40 Hz copy. The 40 Hz low-pass follows the filtering lecture’s general recommendation to retain at least 40 Hz; it is the bandwidth of these saved analysis epochs, not only a plotting setting. @filter-lecture Wider-band continuous files remain available.

ICA training uses a separate 1–100 Hz copy to follow ICLabel’s documented input convention, together with average reference and extended Infomax. @icalabel This is a recommended input specification rather than a software-enforced prohibition on other bands. The paper specifies 0.01–100 Hz acquisition and 0.1–20 Hz ERP/classification filtering, but no distinct ICA fitting filter. @cavanagh2016 Filtering is performed before epoching.

=== Rereferencing

The paper describes an initial average reference and CPz reconstruction before ventral-channel removal and FASTER cleaning; the supplied scripts also contain later rereferencing after interpolation. @cavanagh2016 To avoid known bad channels contributing to the average, this project references over the common accepted EEG set after bad-channel selection, excluding EOG and globally bad EEG channels. The same contributing channels and reference are used for ICA fitting and application. We do not reconstruct CPz. The paper’s description does not establish that untreated bad channels remained in its final reference.

=== Combining runs within subjects

We retain the pipeline’s normal behavior of fitting one ICA per subject/task across all included runs. Our configuration explicitly includes runs 01 and 02, equivalent here to selecting all runs. The supplied MATLAB scripts suggest separate exemplar and gameplay processing, but the dataset warns that these scripts may not match the final published analysis, and the paper does not resolve the pooling question.

Joint fitting is supported by the same-session acquisition and common channel layout. Subject 008’s gameplay tail illustrates why run labels alone do not perfectly separate the tasks. A common channel/reference policy and ICA transform also make source-task and gameplay data more comparable for transfer analysis. Pooling provides more training EEG when one task loses many windows; it does not guarantee better separation, and assumes sufficiently stable spatial mixing across the session. It cannot restore EEG missing between recordings.

This choice requires the union of each subject’s run-specific bad channels and produces one saved ICA model per subject. The installed pipeline fails if the runs’ bad-channel lists differ. Each run is filtered and windowed separately before the fitting epochs are concatenated, so windows never cross recording boundaries. For the run-01 reference this additionally excludes PO7 in 009 and CP1 in 012, which were flagged in run 02.

=== ICA training

An initial pilot used −2…+2 s epochs around selected events, following the paper’s epoching description. Frequent gameplay events produced extensive overlap, repeatedly presenting the same samples to ICA and increasing processing time. Local AutoReject also removed uneven amounts of exemplar data. We explored fixed windows with a 500 µV peak-to-peak screen, then returned to fresh local AutoReject for the selected native workflow. These tests changed more than one setting, so their differences cannot be assigned to a single cause.

The final configuration uses `task_is_rest = True` to create evenly spaced, non-overlapping four-second fitting windows and downsamples to 250 Hz with antialiasing. This option selects a time grid; the experiment is not resting-state EEG. The paper also used 250 Hz for classification, but does not clearly specify ICA downsampling. The final approach was much faster and retained recognizable brain and artifact components, although some pilot comparisons showed weaker alpha structure or imperfect component correspondence. It therefore does not establish unchanged decomposition quality in every subject.

The grid is independent of task events: it can include pauses or task transitions within a run. We chose this simpler native behavior instead of custom inference of task blocks, gaps and padding; it does not guarantee balanced task coverage. The earlier workflow comparisons are historical; this submission packages the selected final preprocessing workflow. The supplied scripts contain a commented whole-epoch mean-removal step before ICA. Our fitting windows receive no baseline correction and no channel interpolation; the −200…0 ms ERP baseline is applied later.

A controlled subject-001 comparison used identical retained data for Picard extended Infomax and regular extended Infomax. Picard took 16 min 47 s, compared with 8 min 49 s for regular Infomax, with closely corresponding components and identical matched automatic exclusions. We therefore retained regular extended Infomax. This single-subject timing comparison is not a general solver benchmark.

Across all 17 subjects, *9,714/11,655 fitting windows (83.35%)* were retained, providing 24.3–54.5 min of training EEG per subject. There are 1,003 components, 56–62 per subject according to available EEG rank after average reference. Post-hoc comparison with event-inferred task spans gives 85.7% oddball, 87.2% gambling and 87.3% gameplay retention. These approximate spans did not control the native grid and do not cover every fitting window.

#figure(
  image("figures/ica_training_rejection.svg", width: 100%),
  caption: [AutoReject rejection of ICA fitting windows by subject. The dotted line is the median across the 17 subject-specific percentages (12.1%); the separately calculated pooled rejection rate is 16.65%.],
)

Rejection varied strongly: subjects 005 and 014 lost 48.2% and 39.8% of fitting windows respectively. Subject 014 retained only 4/33 windows attributed to gambling, but still had more than 24 min of training data overall. Replaying AutoReject reproduced their masks exactly. A counterfactual 500 µV cutoff would retain 94.7% and 96.7%, including 97/99 oddball windows in 005 and all 33 gambling windows in 014. This check evaluated retention only: no alternative ICA fit was performed for those subjects, so it cannot establish better or worse decomposition. Visual review found their existing decompositions acceptable, and we retained them. Fitting-window losses are distinct from later ERP trial losses.

=== ICA component selection

#figure(
  image("figures/ica_component_proposals.svg", width: 100%),
  caption: [Automatic ICLabel exclusion proposals by subject and category, before manual additions.],
)

ICLabel proposed *178 exclusions*: 118 muscle, 47 eye, 12 channel-noise and one heart component. Among the first ten components (ICA000–ICA009), there are 37 automatic eye proposals across the cohort: *2.18 per subject*. Across all component indices, the mean is 47/17 = 2.76. The paper reports one blink component per subject, plus two horizontal-eye components in 13 subjects and one in four: also 2.76 ocular components per subject overall, predominantly among its earliest components. @cavanagh2016 ICLabel’s combined eye category does not distinguish blink from horizontal-eye sources, so these are descriptive counts, not a one-to-one validation.

Visual review found recognizable ocular, muscle and brain components, particularly among early components, but also 15 additional artifacts that did not pass the automatic exclusion threshold. We inspected maps, spectra and time courses and manually excluded these additional sources. The following panel illustrates one such decision.

#figure(
  image("figures/sub-008_ICA005_native_properties.webp", width: 85%),
  caption: [Manual exclusion example: subject 008, ICA005. A focal posterior map and disproportionately large early single-window contribution accompany a 1,641 µV P6 transient. Its channel-noise probability of 0.622 was below the automatic 0.8 threshold; visual review supported exclusion. The fixed-window component average is not a task ERP.],
)

With the manual additions, *193 components* were approved for removal, with a *median of eight per subject*.

#figure(
  image("figures/ica_blink_attenuation.svg", width: 100%),
  caption: [Effect of applying the reviewed ICA models on subjects 005 and 014. Peak-to-peak amplitude in blink-like windows decreases from approximately 273 to 42 µV and 287 to 63 µV respectively. Both curves use 0.1–40 Hz EEG and the same average reference, with individual medians removed for display],
)

The first pipeline pass ends here with saved ICA models, component decisions and reports. A second pass applies those models to EEG prepared for event-related analysis.

==== Note about EOG channels

[merge this section with Low-variation VEOG channels and significantly shorten it. all the reader needs to know is:
1. EOG is excluded from PyPREP, the EEG ICA fit and the average reference.
2. the original study does use eog for ica classification
3. we could have done it too using ica_use_eog_detection = True, but we decided against it
4. the reasons are (do like 3 main points. e.g. veog showed no correlation with blinks for some subjects. amplitudes/offsets/excursions were too large, ...)

]

[summarize all this]
EOG is excluded from scalp-channel PyPREP, the EEG ICA fit and the scalp average reference. Automatic EOG-based component selection remains disabled (`ica_use_eog_detection = False`): the detector first chooses an EOG channel using filtered energy, so heterogeneous channel gains can affect that selection even though a correlation coefficient itself is invariant to fixed positive scaling. There is no arbitrary auxiliary amplitude cutoff or invented rescaling.
Enabling `ica_use_eog_detection = True` would be the pipeline counterpart to selecting components by their relationship to EOG, followed by manual verification as in the paper. We instead use ICLabel and visual EEG-component review. In subject 014, visual inspection additionally suggested possibly exchanged EOG roles, with blink-like spikes apparent in HEOG. This is an inspection-based suspicion, not a verified wiring diagnosis; the channel labels and samples were left unchanged. EOG was used only as auxiliary diagnostic evidence, not as an automatic cleaning input.
The supplied auxiliary values show large offsets and excursions in many recordings, with substantial variation between files. The following summary uses the conventional EEGLAB-to-SI scale applied by the loader. It does not establish the original physiological gain. For each recording, the 99th percentile is calculated from the absolute deviation from that recording's median; this measures excursions without allowing a DC offset to dominate.

#figure(
  table(columns: (14%, 28%, 29%, 29%),
    table.header([*Channel*], [*Recording-median offset range (V)*], [*Median 99th-percentile excursion (mV)*], [*99th-percentile excursion range (mV)*]),
    [VEOG], [−1.603 to 2.992], [135.3], [0.571 to 514.4],
    [HEOG], [−1.515 to 1.126], [101.3], [0.246 to 3,390.3],
  ),
  caption: [Auxiliary values across 34 recordings per channel under the supplied scaling convention. The excursion columns summarize one 99th-percentile value per recording; they are not pooled-sample percentiles.],
)

These values are not uniformly large: some auxiliary channels have much smaller variation. Recording peak-to-peak values reach 4.22 V for VEOG and 5.48 V for HEOG, while the lowest values are approximately 1.95 and 1.04 mV. Thus neither a single representative trace nor one correction factor describes the entire cohort.
A direct sample audit found that MNE's values agree with the native `.fdt` values after the expected factor of $10^(-6)$. This does not support a second, accidental scaling error in the current loader. Acquisition headers and a documented auxiliary gain are absent. The manufacturer's guidance describes auxiliary EOG arrangements whose output scaling can differ from scalp EEG; a gain mismatch is plausible, but the available files do not establish the precise hardware configuration or a universal correction. The values remain unchanged. @brainproducts-eog

[summarize this as well for above]
===== Low-variation VEOG channels

In the audit's selected 120 s intervals, VEOG for sub-003, sub-004 and sub-006 showed both near-zero association with the frontal comparator Fp1−Pz and unusually small filtered variation relative to other VEOG channels. Sub-001 instead showed a strong frontal association. The audit used a 0.5–15 Hz bandpass on those intervals, independently of the main pipeline. This can reflect poor auxiliary signal, an uninformative interval, or another acquisition difference. It does not by itself prove that a channel is unusable throughout both recordings. A matched visualization against sub-001 and the frontal EEG comparator is retained in the working project’s TODO list.
ICLabel provides EEG component proposals for inspection of component maps, spectra and time courses. Auxiliary traces remain available for normalized, per-channel inspection; they are useful timing evidence where a physiological association is visible, without treating their current amplitudes as calibrated EOG. The full audit, segment definitions and per-recording results are in `docs/eog_investigation.md`.

== Second pipeline pass: from ICA models to analysis epochs

After the computation of the ICA models, we have continuous 250 Hz data, 1–100 Hz filtered and epoched to fixed windows. The analysis instead needs 500 Hz, 0.1–40 Hz EEG aligned to specific events, followed by trial cleaning and a prestimulus baseline. We therefore use a separate pipeline configuration, `config_postica.py`, reusing the reviewed models from the previous pipeline execution.

=== 1. Select eligible events from the prepared dataset

The original paper notes that some of the game events occur in quick succession, not giving each of the game events a meaningful epoch. A look at the data confirms this issue:
[picture of a segment where many events happen within the window, and the eeg reaction seems to be most tied to the first of those events. also it should be evident that a 500ms window is enough to circumvent this issue]

The issue affects wall crashes (91) and star collection (93).
To circumvent this issue each of these events is compared to the last retained event of the same type: if the time between these events is below 500 ms, this repeat event gets marked as excluded. This is done by adding a `EXCLUDED_REPEAT__…` prefix to the these affected wall/star events, allowing e.g. ERPs to ignore these events. 
The same effective step was done in the original paper as well and affects 3,010 / 4,443 wall events and 6,556 of 13,776 star events. 

=== 2. Filter, reference and epoch the EEG

The original 500 Hz runs have to be filtered continuously to 0.1–40 Hz. Then, with `task_is_rest = False`, the new pipeline creates −0.2…+0.8 s epochs around the ten eligible event labels to supply the intervals used by our classifiers.

=== 3. Apply the reviewed ICA exclusions

The pipeline loads the saved models and native component-status tables and applies the 193 approved exclusions to compatible continuous EEG and event epochs. EOG samples remain unchanged.

=== 4. Repair or reject analysis epochs

AutoReject is now run again. Previously, it was used before ICA to select training windows. This time, using `reject = "autoreject_local"`, it can repair individual channels in retained epochs as well as reject complete trials. The following counts describe this second stage.

Of 36,258 events that are eligible for the classifiers, 13 windows extend beyond a recording boundary. AutoReject rejects 4,049/36,244 (11.2%) epochs. Among the retained epochs, *28,316 (88.0%)* receive at least one channel repair. The pooled median is four repaired channels per retained epoch; subject medians range from one to nine. Selected interpolation limits are 4, 8 or 16 (default).

[put subject names on the ticks for the top half as well]
#table(columns: (1.8fr, 1fr, 1fr),
  table.header([*Condition*], [*Eligible Events*], [*Final trials*]),
  [Oddball Standard], [2,719], [2,247],
  [Oddball Rare], [681], [539],
  [Gambling Win], [618], [515],
  [Gambling Loss], [617], [494],
  [Player Crash Wall], [1,433], [1,235],
  [Collect Star], [7,220], [6,457],
  [Missile Hit Enemy], [6,495], [5,923],
  [Player Crash Enemy], [1,810], [1,556],
  [Collect Ammo], [1,103], [1,007],
  [Shoot Button], [13,562], [12,222],
)

#figure(
  image("figures/postica_rejection.png", width: 100%),
  caption: [Post-ICA trial rejection and the percentage of retained trials with at least one channel repair.],
)

The rejection is uneven across subjects: 63.2% for 004 and 57.6% for 007. Subject 007 retains only 5/44 rare oddball trials and 34/156 standards; 004 retains 11/40 gambling losses and 17/40 wins. Subject 001 retains only 87/200 oddball trials despite 4.6% rejection pooled across all conditions (this is due to gameplay being a much larger contributor of epochs than oddball or gambling).

We again experimented with lowering rejection to retain more trials, but did not see an improvement in data quality.

=== 5. Baseline correction and saved analysis outputs

After repair/rejection, the pipeline applies the −0.2…0 s baseline (same as the original paper) and finishes the final `*_proc-clean_epo.fif` files reports. Global bad channels remain marked in these epochs. 

== ERP checks before and after cleaning

The following figure compares the ERPs for exemplary trials before and after our main cleaning. For the Oddball task, it can be seen that the Standard event stays much closer to the baseline after cleaning, and both events no stay much more stable after 0.6 s.
For the Gambling task, cleaning mostly corrects the dip seen at \~0.5 s from -4 uV to 0 uV.

[remove the baked in text at the bottom of the image]
#figure(
  image("figures/erp_cleaning_stages.svg", width: 100%),
  caption: [Grand-average condition ERPs before ICA (top) and after ICA plus final AutoReject repairs (bottom), at Pz for oddball (left, N = 17) and Cz for gambling (right, N = 14). Both rows use identical final retained trial IDs, 0.1–40 Hz filtering and −0.2…0 s baseline. Shading is the pointwise 95% t confidence interval],
)


#bibliography("references.bib", style: "ieee", title: [References])
