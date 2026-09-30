# Floating exploration — not proof evidence

This separate, orthonormal-tensor construction guided the choice of the grouped-pair inequality. All 35 exploratory sector calculations were completed; their maximum estimated overlap was approximately 0.217711. No estimate here is used as a rigorous lower gap bound or a rigorous upper overlap bound.

The development source was copied here with only portable path changes and deterministic tie-breaking for equivalent choices of branch order. `kvbs.py` is copied from the original submission, under its MIT notice. The copied portable exploration has not been rerun after these path/tie-breaking changes; the stored JSON records the earlier development execution. It may return a different ordering of equivalent groupings when rerun. The authoritative deterministic code is `../code/certified/`, with a fresh full run in `../evidence/new/`.

To explore in a working copy, run `OPENBLAS_NUM_THREADS=2 python exploratory/run_angles.py`. This is optional and not part of the proof replay.
