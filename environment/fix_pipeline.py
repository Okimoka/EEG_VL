"""
Apply upstream pipeline fix #1164 to the pinned course version.
The branch still has the bug with ICA only fitting on one run
"""
from importlib.metadata import distribution

OLD = '''            if cfg.ica_l_freq is not None or h_freq is not None:
                logger.info(**gen_log_kwargs(message=msg))
                raw.filter(l_freq=cfg.ica_l_freq, h_freq=h_freq, n_jobs=1)
            del nyq, h_freq
'''

NEW = '''            if msg:
                logger.info(**gen_log_kwargs(message=msg))
            del nyq

        if cfg.ica_l_freq is not None or h_freq is not None:
            raw.filter(l_freq=cfg.ica_l_freq, h_freq=h_freq, n_jobs=1)
'''

path = distribution("mne-bids-pipeline").locate_file("mne_bids_pipeline/steps/preprocessing/_06a1_fit_ica.py")
path.write_text(path.read_text().replace(OLD, NEW))
