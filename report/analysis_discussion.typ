#pagebreak(weak: true)
= Discussion

== Comparison with the original findings

EEGNet discriminates the conditions in the exemplary tasks across held-out participants and shows event-related shifts in gameplay classifications. The comparison in the P3 and feedback intervals gives different outcomes for oddball and gambling.

For oddball, gameplay events shift towards rare labels relative to controls, but the mean fractions remain below 50% at 300–600 ms. We therefore do not recover the original predominance of target classifications in this analysis. This result concerns the selected interval and cross-participant models; it does not establish an absence of attention-related responses during gameplay.

For gambling, hits receive predominantly win labels at 200–350 ms and more win labels than random gameplay. This agrees with the direction of the original direct win-versus-loss result. However, the separate win-versus-gameplay classifier assigns most hit trials to gameplay, so agreement does not extend to the full set of original contrasts. Crashes also shift towards win, limiting a selective reward interpretation. The original paper similarly failed to establish its predicted loss-like crash response and reported some win-like classifications for crashes. @cavanagh2016

These findings constitute partial agreement in a reanalysis of the same dataset. The differences in preprocessing, classifier type, training population and temporal intervals prevent attribution of the discrepancies to any one methodological choice.

== Limitations

Transferred labels indicate similarity according to the learned classifier, not direct evidence of attention, reward or punishment. EEGNet may use shared amplitude features, sensory or motor activity, differences between tasks, or residual artifacts. The analysis does not establish which features drive its decisions or whether they identify the same psychological process in both settings.

The participant sample is limited to 17 for oddball and 14 for gambling. Large gameplay trial counts improve within-participant estimates but do not increase the number of independent participants or labelled exemplary training examples. Several participants have few retained trials in one exemplary category, and rejection and interpolation vary across recordings.

Participant exclusion applies to classifier fitting and scaling. Preprocessing was completed retrospectively using each participant’s whole recording, so the results are conditional on that cleaning procedure. They do not represent a fully independent evaluation from unseen raw data to prediction. Furthermore, the held-out models share training participants. Results average three fixed training runs per split. The participant bootstrap intervals remain conditional on those fitted networks and do not capture all training variability or dependence between model results.

The analysis is exploratory and uses the original dataset rather than an independent replication sample. The fixed intervals do not provide a complete time course, and neighbouring gameplay events may contribute to the measured responses. Positive transfer shifts therefore support further investigation rather than a definitive identification of canonical neural events.

== Further work

A matched within-participant EEGNet analysis would bring the training design closer to the paper and help assess the effect of pooling participants. The paper’s shorter moving windows and derivative representation could then be examined as additional methodological comparisons. Stronger claims about psychological specificity would require controls that distinguish sensory and motor responses from salience or reward, followed by evaluation in an independent sample.

= Conclusion

Using EEGNet with default model parameters, we identify EEG distinctions between exemplary categories that generalize across participants and observe shifts in gameplay classifications. The feedback-interval analysis agrees with the direction of the paper’s win-versus-loss result for missile hits. It does not recover predominantly target-like classifications in the selected oddball interval or predominantly exemplary classifications against gameplay controls.

The results therefore show partial agreement with the original transfer findings. Win-like crash responses and uncertainty about the features used by EEGNet limit their interpretation as specific measures of attention, reward or punishment.
