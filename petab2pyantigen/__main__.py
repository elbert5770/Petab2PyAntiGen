"""python -m petab2pyantigen <command>

  antimony  SBML -> clean Antimony (Tellurium, no `is` clutter)
  emit      write parameters.csv, InitialConditions.csv and the wide data csv
  verify    compare a PyAntiGen project with the PEtab reference (see verify.py)
  stats     size of the PEtab problem vs the PyAntiGen project
"""
import argparse
import sys
from pathlib import Path


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="petab2pyantigen", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    a = sub.add_parser("antimony", help="SBML -> clean Antimony")
    a.add_argument("sbml")
    a.add_argument("-o", "--out", help="write here instead of stdout")

    e = sub.add_parser("emit", help="write the data files of a PyAntiGen project")
    e.add_argument("petab", help="PEtab v1 directory")
    e.add_argument("--root", default=".", help="repository root (holds antimony_models/ and data/)")
    e.add_argument("--model", required=True, help="PyAntiGen model name, e.g. Zheng_PNAS2012")

    v = sub.add_parser("verify", help="compare with the PEtab reference")
    v.add_argument("rest", nargs=argparse.REMAINDER)

    s = sub.add_parser("stats", help="PEtab vs PyAntiGen size")
    s.add_argument("petab")
    s.add_argument("--root", default=".")
    s.add_argument("--model", required=True)
    s.add_argument("--modules", nargs="+", required=True, help="files under antimony_modules/")

    args = ap.parse_args(argv)
    if args.cmd == "antimony":
        from .sbml_to_antimony import sbml_to_antimony
        r = sbml_to_antimony(args.sbml)
        if args.out:
            Path(args.out).write_text(r.text, encoding="utf-8")
        else:
            sys.stdout.write(r.text)
        print(f"removed {r.removed}", file=sys.stderr)
        if r.leftover_is:
            print(f"WARNING: identifier `is` statements remain: {r.leftover_is}", file=sys.stderr)
    elif args.cmd == "emit":
        from . import emit
        from .petab_reader import load_petab
        p, root = load_petab(args.petab), Path(args.root)
        base = root / "antimony_models" / args.model
        base.mkdir(parents=True, exist_ok=True)
        (root / "data").mkdir(exist_ok=True)
        emit.write_parameters_csv(p, base / f"{args.model}_parameters.csv")
        emit.write_initial_conditions_csv(p, base / f"{args.model}_InitialConditions.csv")
        emit.write_wide_data_csv(p, root / "data" / f"{args.model}.csv")
        print(f"wrote data files for {args.model} under {root}")
    elif args.cmd == "verify":
        from .verify import main as verify_main
        return verify_main(args.rest)
    elif args.cmd == "stats":
        from .petab_reader import load_petab
        from .sbml_to_antimony import sbml_to_antimony
        from .stats import report
        p = load_petab(args.petab)
        print(report(p, args.root, args.model, args.modules, sbml_to_antimony(p.sbml_path).text))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
