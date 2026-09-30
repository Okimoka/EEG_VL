"""Run only the installed pipeline's continuous ICA application function.

The installed stage main() unconditionally processes epochs first. This small
adapter selects its raw function without patching the pipeline or inventing
epochs. Filtering, reference, TSV exclusion loading, application, native caching,
raw saving and standard reports retain the installed pipeline implementation.
"""
from pathlib import Path
from mne_bids_pipeline._config_import import _import_config
from mne_bids_pipeline._config_utils import _get_ssrt
from mne_bids_pipeline._parallel import get_parallel_backend,parallel_func
from mne_bids_pipeline._run import save_logs
from mne_bids_pipeline.steps.preprocessing import _08a_apply_ica as step
ROOT=Path(__file__).resolve().parents[1]

def main():
    config=_import_config(config_path=ROOT/'config_broadband.py')
    assert config.l_freq==.1 and config.h_freq==100 and config.raw_resample_sfreq is None
    assert config.reject is None and config.ica_use_icalabel
    jobs=_get_ssrt(config=config,which=('runs',));assert len(jobs)==34
    with get_parallel_backend(config.exec_params):
        parallel,run=parallel_func(step.apply_ica_raw,exec_params=config.exec_params,n_iter=len(jobs))
        logs=parallel(run(cfg=step.get_config(config=config,subject=subject),exec_params=config.exec_params,
                         subject=subject,session=session,run=run_id,task=task)
                      for subject,session,run_id,task in jobs)
    # Supply the native step identity expected by its execution-log writer.
    __mne_bids_pipeline_step__=Path(step.__file__)
    save_logs(config=config,logs=logs)
    print('Completed native continuous ICA application for 34 recordings.',flush=True)

if __name__=='__main__':main()
