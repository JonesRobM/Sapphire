# Tutorial data

Four bimetallic 14 ns MD data sets (two independent runs each, ASE `.traj`):

| System | Morphology | Protocol |
|---|---|---|
| AgNi | random alloy, 80/20 | freezing 1000 → 300 K, −50 K / ns |
| AuPd | Janus, 80/20 | freezing 1000 → 300 K, −50 K / ns |
| AuPt | random alloy, 80/20 | melting 300 → 1000 K, +50 K / ns |
| CuNi | Janus, 80/20 | melting 300 → 1000 K, +50 K / ns |

The full tarballs (~58 MB each) are **not in git**. They are published as GitHub
release assets and fetched on demand:

```python
from Sapphire.Tutorials import data
path = data.fetch("AuPt")           # ~/.cache/sapphire/AuPt_Traj/
xyz  = data.sample("AuPt", ".")     # bundled every-20th-frame sample, offline
traj = data.synthetic(n_frames=50)  # seeded EMT toy trajectory, no download
```

`samples/` holds one down-sampled run per system (every 20th frame, gzip extxyz).
