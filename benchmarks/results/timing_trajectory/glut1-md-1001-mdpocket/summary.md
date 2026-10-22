<!-- System: glut1-md-1001-mdpocket. Mutant GLUT1 membrane MD (prod.gro + prod2_short.xtc, 1001 frames, CHARMM naming) with a prepared candidate-site region (reference frame 0, 0.2 A grid, probe 1.4 A, margin 2 A). GLUT1 has no through-channel, so no profile cases. Machine: Intel(R) Core(TM) Ultra 7 268V, 8 logical CPUs, Linux-7.0.0-31-generic-x86_64-with-glibc2.39, Python 3.13.13; load average at start [0.86181640625, 1.52490234375, 1.72314453125]. CPUs allowed per run: 2 (max MHz 4900); frames: 1000; repeats: 1. Timings depend on machine load (see the load columns). -->

| Case | Tool | grid (Å) | CPUs / workers | frames (with result) | median wall (s) | median CPU (s) | frames/s | peak RSS, largest process / tree (MB) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| mdpocket_input | MDpocket input (MDAnalysis fit and write) | – | 1 / 1 | 1000 (1000) | 9.26 | 9.24 | 108 | 186.5 / 185.4 |
| mdpocket | MDpocket 4.1.5 (selected pocket) | 1 (MDpocket grid) | 1 / 1 | 1000 (1000) | 9326 | 9323 | 0.107 | 114.9 / 115 |
