"""Command line front end to Sapphire.

The analysis itself is :func:`Sapphire.api.run_config`; this module only turns
argv into a :class:`Sapphire.api.Config` and reports what happened. Two rules
shape it, both of which come from running on a remote batch node where nobody
is watching the terminal:

* a run that could not compute what was asked exits non-zero;
* nothing optional (a quote lookup, a missing plotting extra) may fail a run.

Usage::

    sapphire run movie.xyz -o out/ -q pdf,adj,agcn,cna_sigs --frames 0:1000:10
    sapphire run --config analysis.toml
    sapphire run movie.xyz --dry-run          # resolve and print the plan, compute nothing
    sapphire quantities                       # what can be asked for
"""
from __future__ import annotations

import argparse
import os
import sys

from Sapphire.api import NOT_PERSISTED, READER_ALIASES, REQUIRES


def parse_frames(spec):
    """``START:END:STEP`` -> ``(start, end, step)``; empty fields take defaults."""
    if spec is None:
        return (0, None, 1)
    parts = spec.split(':')
    if len(parts) > 3:
        raise argparse.ArgumentTypeError(f"--frames takes START:END:STEP, got {spec!r}")
    parts += [''] * (3 - len(parts))
    try:
        start = int(parts[0]) if parts[0] else 0
        end = int(parts[1]) if parts[1] else None
        step = int(parts[2]) if parts[2] else 1
    except ValueError:
        raise argparse.ArgumentTypeError(f"--frames fields must be integers, got {spec!r}")
    if step < 1:
        raise argparse.ArgumentTypeError(f"--frames step must be >= 1, got {step}")
    if end is not None and end <= start:
        raise argparse.ArgumentTypeError(f"--frames end ({end}) must exceed start ({start})")
    return (start, end, step)


def _csv(value):
    """``a,b , c`` -> ``['a', 'b', 'c']``; an empty string means an empty list."""
    return [item.strip() for item in value.split(',') if item.strip()]


def build_parser():
    parser = argparse.ArgumentParser(
        prog='sapphire',
        description='Post-process metallic nanoalloy trajectories.',
    )
    from Sapphire import __version__
    parser.add_argument('--version', action='version', version=f'sapphire {__version__}')
    sub = parser.add_subparsers(dest='command', required=True)

    run = sub.add_parser('run', help='analyse a trajectory')
    run.add_argument('trajectory', nargs='?', help='trajectory file (xyz, extxyz, traj, ...)')
    run.add_argument('--config', metavar='TOML',
                     help='load settings from a TOML file written by a previous run')
    run.add_argument('-o', '--out-dir', default='sapphire_run/', metavar='DIR')
    run.add_argument('-q', '--quantities', type=_csv, metavar='A,B,C',
                     help='comma separated; see `sapphire quantities`')
    run.add_argument('--homo', type=_csv, metavar='A,B', help='per-species quantities')
    run.add_argument('--hetero', type=_csv, metavar='A,B', help='between-species quantities')
    run.add_argument('--species', type=_csv, metavar='Au,Pt',
                     help='default: inferred from the first frame')
    run.add_argument('--frames', type=parse_frames, metavar='START:END:STEP',
                     help='frame slice, e.g. 0:1000:10 (default: every frame)')
    run.add_argument('--band', type=float, default=0.05, metavar='A',
                     help='KDE bandwidth in angstrom for the cutoff (default: 0.05)')
    run.add_argument('-j', '--jobs', type=int, default=1, metavar='N',
                     help='analyse frames across N processes (default: 1)')
    run.add_argument('--adj-format', choices=('npz', 'text'), default='npz',
                     help="adjacency on disk: sparse 'npz' (default) or dense 'text'. "
                          "Dense costs O(N^2) bytes per frame -- 8 MB at N=2000, 160 GB "
                          "over 20k frames. Use `sapphire expand` to get text on demand.")
    run.add_argument('--no-quote', action='store_true',
                     help='skip the quote lookup in the log header (it needs network access)')
    run.add_argument('--strict', action='store_true',
                     help='stop at the first failing quantity instead of continuing')
    run.add_argument('--no-overwrite', action='store_true',
                     help='keep results from a previous run in the output directory')
    run.add_argument('--dry-run', action='store_true',
                     help='resolve and print the plan, then exit without computing')

    expand = sub.add_parser(
        'expand', help='write sparse adjacency matrices out as dense text')
    expand.add_argument('run_dir', help='a Sapphire output directory')
    expand.add_argument('-o', '--out-dir', metavar='DIR',
                        help='where to write (default: alongside the npz files)')
    expand.add_argument('--frames', type=parse_frames, metavar='START:END:STEP',
                        help='only these frame numbers (default: all)')
    expand.add_argument('--stdout', action='store_true',
                        help='write to standard output instead of files, for piping')

    sub.add_parser('quantities', help='list the supported quantities')
    return parser


def cmd_expand(args):
    import pathlib

    import numpy as np
    import scipy.sparse as spa

    from Sapphire.Post_Process.Adjacent import _write_int_matrix
    from Sapphire.parallel import _PER_FRAME

    run_dir = pathlib.Path(args.run_dir)
    if not run_dir.is_dir():
        print(f"sapphire: no such directory: {run_dir}", file=sys.stderr)
        return 2

    # Full and homo matrices live under Adjacency/, hetero ones under Time_Dependent/.
    # Search by name so every flavour is covered.
    files = sorted(p for p in run_dir.rglob('*File*.npz') if _PER_FRAME.fullmatch(p.name))
    if not files:
        plain = [p for p in run_dir.rglob('*File*') if _PER_FRAME.fullmatch(p.name)]
        if plain:
            print(f"sapphire: {run_dir} already holds dense text ({len(plain)} matrices); "
                  f"nothing to expand", file=sys.stderr)
            return 0
        print(f"sapphire: no adjacency matrices under {run_dir}", file=sys.stderr)
        return 2

    def frame_of(path):
        return int(path.name[:-4].split('File')[-1])

    if args.frames is not None:
        start, end, step = args.frames
        wanted = range(start, end if end is not None else 1 << 30, step)
        files = [f for f in files if frame_of(f) in wanted]
        if not files:
            print(f"sapphire: no adjacency matrices in frame range {start}:{end}:{step}",
                  file=sys.stderr)
            return 2

    out_dir = pathlib.Path(args.out_dir) if args.out_dir else None
    written = 0
    for f in sorted(files, key=lambda p: (p.parent.name, p.name)):
        dense = spa.load_npz(f).toarray().astype(np.int32)
        if args.stdout:
            _write_int_matrix(sys.stdout, dense)
        else:
            # Alongside the npz by default; under -o keep the run's own layout so
            # Adjacency/ and Time_Dependent/ matrices stay distinguishable.
            if out_dir is None:
                target = f.with_suffix('')
            else:
                target = out_dir / f.relative_to(run_dir).with_suffix('')
            target.parent.mkdir(parents=True, exist_ok=True)
            with open(target, 'w') as handle:
                _write_int_matrix(handle, dense)
        written += 1
    if not args.stdout:
        where = out_dir if out_dir is not None else run_dir
        print(f"sapphire: expanded {written} matrices into {where}")
    return 0


def cmd_quantities(_args):
    from Sapphire.Utilities.Supported import Supported
    sup = Supported()
    for title, items in (('Full', sup.Full()), ('Homo', sup.Homo()), ('Hetero', sup.Hetero())):
        print(f"\n{title}:")
        for name in sorted(items):
            deps = REQUIRES.get(name)
            note = f"   (requires {', '.join(deps)})" if deps else ''
            print(f"  {name}{note}")
    print()
    return 0


def cmd_run(args):
    from Sapphire.api import Config, run_config

    if args.no_quote:
        os.environ['SAPPHIRE_NO_QUOTE'] = '1'

    if args.config:
        cfg = Config.from_toml(args.config)
        if args.trajectory:                       # an explicit path overrides the stored one
            cfg.trajectory = args.trajectory
        if args.out_dir != 'sapphire_run/':
            cfg.out_dir = args.out_dir
    else:
        if not args.trajectory:
            print("sapphire: give a trajectory, or --config with a TOML file", file=sys.stderr)
            return 2
        cfg = Config(trajectory=args.trajectory, out_dir=args.out_dir)
        if args.quantities is not None:
            cfg.quantities = args.quantities
        if args.frames is not None:
            cfg.frames = args.frames
        cfg.band = args.band

    # These apply whether or not a config file was loaded; an explicit flag beats a
    # stored setting, but a default must not silently override one.
    if args.homo is not None:
        cfg.homo = args.homo
    if args.hetero is not None:
        cfg.hetero = args.hetero
    if args.species is not None:
        cfg.species = args.species
    cfg.strict = args.strict
    cfg.overwrite = not args.no_overwrite
    if args.adj_format != 'npz' or not args.config:
        cfg.adj_format = args.adj_format
    if args.jobs != 1 or not args.config:
        cfg.jobs = args.jobs

    try:
        cfg.validate()                       # resolves prerequisites into cfg.quantities
    except FileNotFoundError as exc:
        print(f"sapphire: no such trajectory: {exc}", file=sys.stderr)
        return 2
    except ValueError as exc:
        print(f"sapphire: {exc}", file=sys.stderr)
        return 2

    if cfg.resolved_dependencies:
        print(f"sapphire: also computing {', '.join(sorted(set(cfg.resolved_dependencies)))} "
              f"(required by what you asked for)", file=sys.stderr)

    start, end, step = cfg.frames
    if args.dry_run:
        print(f"trajectory : {cfg.trajectory}")
        print(f"out_dir    : {cfg.out_dir}")
        print(f"frames     : {start}:{end if end is not None else 'end'}:{step}")
        print(f"quantities : {', '.join(cfg.quantities)}")
        print(f"homo       : {', '.join(cfg.homo) if cfg.homo else '-'}")
        print(f"hetero     : {', '.join(cfg.hetero) if cfg.hetero else '-'}")
        print(f"band       : {cfg.band}")
        print(f"adj_format : {cfg.adj_format}")
        print(f"jobs       : {cfg.jobs}")
        return 0

    reader = run_config(cfg)

    # Process logs failures per quantity and carries on, so a zero exit status on
    # its own means very little. Ask the Reader what actually landed on disk.
    produced = set(reader.available())

    def landed(q):
        """Did quantity ``q`` actually reach disk?

        Two wrinkles: the Reader files some quantities under another name, and homo/hetero
        results carry a species suffix (hocomdist -> hocomdistAu).
        """
        key = READER_ALIASES.get(q, q)
        return any(k == key or k.startswith(key) for k in produced)

    checkable = [q for q in cfg.quantities if q not in NOT_PERSISTED]
    missing = [q for q in checkable if not landed(q)]
    if missing and produced:
        print(f"sapphire: no output for {', '.join(missing)} - see "
              f"{os.path.join(cfg.out_dir, 'Sapphire_Errors.log')}", file=sys.stderr)
        return 1
    if not produced:
        print(f"sapphire: the run produced no output - see "
              f"{os.path.join(cfg.out_dir, 'Sapphire_Errors.log')}", file=sys.stderr)
        return 1
    print(f"sapphire: wrote {len(produced)} quantities to {cfg.out_dir}")
    return 0


def main(argv=None):
    args = build_parser().parse_args(argv)
    handler = {'run': cmd_run, 'quantities': cmd_quantities, 'expand': cmd_expand}[args.command]
    try:
        return handler(args)
    except KeyboardInterrupt:
        print("\nsapphire: interrupted", file=sys.stderr)
        return 130
    except BrokenPipeError:
        # A downstream reader (head, less) closed the pipe. That is ordinary shell usage,
        # not an error: point stdout at devnull so the interpreter does not complain on
        # shutdown, and report the status a signal-killed process would.
        os.dup2(os.open(os.devnull, os.O_WRONLY), sys.stdout.fileno())
        return 141


if __name__ == '__main__':
    sys.exit(main())
