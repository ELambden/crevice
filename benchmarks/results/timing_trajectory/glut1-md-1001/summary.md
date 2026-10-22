<!-- System: glut1-md-1001. Mutant GLUT1 membrane MD (prod.gro + prod2_short.xtc, 1001 frames, CHARMM naming) with a prepared candidate-site region (reference frame 0, 0.2 A grid, probe 1.4 A, margin 2 A). GLUT1 has no through-channel, so no profile cases. Machine: Intel(R) Core(TM) Ultra 7 268V, 8 logical CPUs, Linux-7.0.0-31-generic-x86_64-with-glibc2.39, Python 3.13.13; load average at start [1.11083984375, 1.453125, 1.68115234375]. CPUs allowed per run: 1 (max MHz 5000); frames: 1000; repeats: 3. Timings depend on machine load (see the load columns). -->

| Case | Tool | grid (Å) | CPUs / workers | frames (with result) | median wall (s) | median CPU (s) | frames/s | peak RSS, largest process / tree (MB) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| crevice_region | CREVICE trajectory --region | 0.2 (volume grid) | 1 / 1 | 1000 (1000) | 1665 | 1665 | 0.601 | 7570 / 7553 |
| crevice_region_0.35 | CREVICE trajectory --cavities (reference region at 0.35 A) | 0.35 (volume grid) | 1 / 1 | 1000 (1000) | 607 | 606 | 1.65 | 7450 / 7449 |
| crevice_region_0.5 | CREVICE trajectory --cavities (reference region at 0.5 A) | 0.5 (volume grid) | 1 / 1 | 1000 (1000) | 432 | 432 | 2.32 | 7403 / 7403 |
| crevice_region_0.7 | CREVICE trajectory --cavities (reference region at 0.7 A) | 0.7 (volume grid) | 1 / 1 | 1000 (1000) | 370 | 370 | 2.7 | 7310 / 7310 |
