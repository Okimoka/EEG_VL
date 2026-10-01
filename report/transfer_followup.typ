== Exploratory follow-up

Following the whole-epoch results, we examined temporal information, event overlap and direct differences between game-event types. These analyses were motivated by the initial outcomes and are reported as exploratory.

=== Temporal intervals

Separate linear classifiers were fitted to four fixed intervals, retaining the original participant splits. These models assess where discriminative information is available; they do not identify the features used by the whole-epoch EEGNet.

#table(columns: (1.4fr, 1fr, 1fr),
  table.header([*Time after event*], [*Oddball balanced accuracy*], [*Gambling balanced accuracy*]),
  [0–200 ms], [69.6%], [55.8%],
  [200–350 ms], [69.7%], [60.7%],
  [300–600 ms], [73.8%], [61.9%],
  [600–800 ms], [67.1%], [52.6%],
)

Oddball accuracy is numerically highest at 300–600 ms, the interval used for the P3 comparison. Discrimination is also present at 0–200 ms, so the classification results cannot be attributed exclusively to a late P3 response. Gambling performance is strongest in the two middle intervals and approaches 50% at 600–800 ms. The intervals overlap and differ in duration and feature count, limiting direct comparisons of their performance.

For the 200–350 ms gambling model, win-label fractions increase over random gameplay by 15.5 percentage points for hits, 25.0 for enemy crashes and 16.2 for wall crashes. The transfer pattern is therefore also present in a classifier restricted to the feedback-response interval.

Fixed ERP contrasts provide a complementary check of the laboratory responses. Rare minus standard at Pz over 300–600 ms is +4.60 µV, and win minus loss at Cz over 200–350 ms is +2.08 µV. Both pass Holm correction across these two exploratory tests. These differences support the presence of recognizable laboratory responses after cleaning, without establishing the specificity of their transfer to gameplay.

=== Sensitivity to neighbouring events

To assess overlapping responses, we repeated the existing prediction summaries for events with no other recorded hit, crash, star, ammunition outcome or game-over in the −200 to +800 ms interval. Control epochs were required to contain none of these outcomes. Shots were permitted in both cases.

The restriction retains 4,793 of 5,923 hits, 1,093 of 1,556 enemy crashes and 118 of 1,235 wall crashes. Eligibility requires at least 20 event trials and 20 controls per participant. No participant meets this criterion for isolated wall crashes, so that contrast cannot be estimated under the specified rule.

#table(columns: (1.2fr, 1.1fr, .5fr, 1fr, 1fr),
  table.header([*EEGNet label*], [*Event*], [*People*], [*All trials, same people (pp)*], [*Isolated events (pp)*]),
  [Rare], [Missile hit], [17], [+10.9], [+14.1],
  [Rare], [Enemy crash], [15], [+25.2], [+29.0],
  [Win], [Missile hit], [14], [+15.8], [+15.7],
  [Win], [Enemy crash], [12], [+22.7], [+21.1],
)

The rare- and win-label shifts persist for hits and enemy crashes. Neighbouring recorded outcomes therefore do not fully account for these effects. The restriction does not remove contributions from shots, earlier response tails or unrecorded events, and it changes the set of gameplay situations represented in the analysis.

=== Direct event contrasts

For whole-epoch EEGNet, the hit-minus-enemy-crash difference in win-label fractions is −5.5 percentage points, while the hit-minus-wall-crash difference is +6.3 points. After Holm correction across four contrasts covering both model families, the first does not reach the exploratory 0.05 threshold (p = 0.088), whereas the second does (p = 0.022). The numerical ordering of the event means therefore does not establish reliable EEGNet discrimination between hits and enemy crashes.

Mean continuous classifier scores preserve the directions of the original event-minus-control differences, indicating that the observations are not restricted to the 0.5 label threshold. These scores remain uncalibrated with respect to psychological states. Complete results and supplementary figures are saved in `artifacts/transfer_followup/`.
