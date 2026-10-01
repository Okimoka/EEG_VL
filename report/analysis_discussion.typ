#pagebreak(weak: true)
= Discussion

This project tested whether the original findings remain observable with a different preprocessing and analysis pipeline. Preparing the shared dataset required corrections to channel metadata, electrode coordinates and event timing. PyPREP, reviewed ICA removal and AutoReject reduced visible artifacts while retaining the main ERP distinctions. However, trial loss was uneven, particularly for subjects 004 and 007, and 88% of retained analysis epochs required channel repair. Successful pipeline execution therefore does not establish that cleaning was equally effective across subjects or tasks. Despite the large differences in preprocessing, many of the results, such as average number of eye components, number of electrodes that were marked as bad, and the resulting ERP were very similar.

The analysis shows partial agreement with the original findings. EEGNet distinguishes the exemplary categories for subjects whose data were not used to train the classifier. Gameplay events receive more rare labels than random gameplay, but their mean rare-label fractions remain below 50% in the P3 interval. Missile hits receive predominantly win labels in the feedback interval, matching the direction reported in the original paper, but crashes also shift towards win, limiting the extent of interpretation @cavanagh2016

Future work should vary only one choice at a time when varying the analysis setup. In our analysis, results varied vastly even for small changes such as different random seeds. One focus should be to identify where these differences come from and whether they can be better controlled and accounted for. It should also be a goal to be able to assess that the classifier actually bases its decision off of real psychological properties, and not any sort of artifact or other irrelevant data.

= Conclusion

The alternative workflow was able to train classifiers that distinguished exemplary task conditions, and was able to reproduce some of the same gameplay transfer observations. However, it could not reproduce all of the original findings, showing that changes in analysis methodology can be very sensible to changes, and should only be done in a controlled way.
