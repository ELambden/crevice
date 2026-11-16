<!-- System: ompf-trimer-pore1. OmpF trimer membrane MD (1 M KCl, 0.05 V/nm), 1001 frames over 100 ns (100 ps per frame); pore 1, point at the widest point (2.84 A clearance) of the constriction plane in frame 0, 5.5 A from the Arg42/Arg82/Arg132 CZ centroid; private data, not distributed Machine: Intel(R) Core(TM) Ultra 7 268V, 8 logical CPUs, Linux-7.0.0-42-generic-x86_64-with-glibc2.39, Python 3.13.13; load average at start [0.484375, 0.93994140625, 1.29150390625]. CPUs allowed per run: 1 (max MHz 5000); frames: all; repeats: 1. Timings depend on machine load (see the load columns). -->

CREVICE row from the second run (10 Nov 2026, run with `--max-frames 0`, which lifts the frame limit of the profile command); HOLE row from the first run, minutes earlier, same machine and settings.

| Case | Tool | grid (Å) | CPUs / workers | frames (with result) | median wall (s) | median CPU (s) | frames/s | peak RSS, largest process / tree (MB) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| crevice_profile | CREVICE trajectory (channel profile) | 0.5 (section grid) | 1 / 1 | 1001 (933) | 986 | 986 | 1.01 | 658.4 / 606.7 |
| hole | HOLE 2.3.1 (mdahole2, per frame) | 0.2 (plane spacing) | 1 / 1 | 1001 (1001) | 1187 | 1187 | 0.843 | 697.8 / 696.9 |
