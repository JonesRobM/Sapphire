# Shrinking `.git` from 721 MB (runbook for Rob)

Everything here is run by **Rob**, not Claude. Do it on a *fresh clone* so the
working repository is never at risk; the original stays as a backup until you
are satisfied.

## 0. Preconditions
- The four tarballs are safely stored somewhere outside git (e.g. `~/SapphireData/`):
  `AgNi.tar.gz AuPd.tar.gz AuPt.tar.gz CuNi.tar.gz` — SHA-256s are in
  `main/Sapphire/Tutorials/data.py`.
- `brew install git-filter-repo`
- All open branches/PRs merged or noted — the rewrite invalidates every existing clone.

## 1. Publish the data as a release asset (once)
```bash
gh release create data-v1 --title "Tutorial MD data v1" \
  --notes "Four bimetallic 14 ns trajectories used by the tutorials. See Tutorials/Data/ReadMe.md." \
  main/Sapphire/Tutorials/Data/AgNi.tar.gz main/Sapphire/Tutorials/Data/AuPd.tar.gz \
  main/Sapphire/Tutorials/Data/AuPt.tar.gz main/Sapphire/Tutorials/Data/CuNi.tar.gz
```
Then verify one download end-to-end: `python -c "from Sapphire.Tutorials import data; print(data.fetch('AuPt'))"`.

## 2. Rewrite on a fresh clone
```bash
cd ~/tmp && git clone --mirror https://github.com/JonesRobM/Sapphire.git Sapphire-rewrite.git
cd Sapphire-rewrite.git
# strip: big tarballs, every build artefact ever committed, eggs, copied .git internals
git filter-repo --strip-blobs-bigger-than 2M \
  --path build --path main/build --path main/dist --path main/.eggs \
  --path-glob '*.egg-info' --path-glob '*.egg' --path-glob '**/__pycache__' \
  --path-glob '**/.ipynb_checkpoints' \
  --path HEAD --path config --path description --path main/HEAD --path main/config \
  --path main/description --path main/hooks --path main/info \
  --invert-paths
git reflog expire --expire=now --all && git gc --prune=now --aggressive
du -sh objects        # ~10 MB after the 2026-08-28 rewrite
```

## 3. Push and re-clone
```bash
git remote add origin https://github.com/JonesRobM/Sapphire.git   # filter-repo removes it
git push --force --all origin
git push --force --tags origin
# then, from a clean directory:
git clone https://github.com/JonesRobM/Sapphire.git
```
Everyone else with a clone must re-clone (or `git fetch && git reset --hard origin/main`).

## 4. Guard against recurrence
- `.gitignore` already blocks `*.traj` and `Tutorials/Data/*.tar.gz`.
- `.pre-commit-config.yaml` refuses files over 2 MB (`pre-commit install` once per clone).
