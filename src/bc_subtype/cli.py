import argparse
import json

from bc_subtype.data import prepare
from bc_subtype.download import download
from bc_subtype.evaluate import evaluate
from bc_subtype.selection import rank_folds


def main():
    parser = argparse.ArgumentParser(description="TCGA-BRCA data and baseline setup")
    commands = parser.add_subparsers(dest="command", required=True)
    fetch = commands.add_parser("download", help="Fetch and verify pinned public source files")
    fetch.add_argument("--raw-dir", default="data/raw")
    prep = commands.add_parser("prepare", help="Match, log-transform and audit the cohort")
    prep.add_argument("--raw-dir", default="data/raw")
    prep.add_argument("--out-dir", default="data/processed")
    run = commands.add_parser("baseline", help="Run reference baseline validation")
    run.add_argument("--data-dir", default="data/processed")
    run.add_argument("--config", default="configs/baselines.json")
    run.add_argument("--output-dir", default="results/baselines")
    rank = commands.add_parser("rank", help="Rank genes (ANOVA, L1) on every training fold")
    rank.add_argument("--data-dir", default="data/processed")
    rank.add_argument("--config", default="configs/baselines.json")
    rank.add_argument("--output-dir", default="results/selection")
    rank.add_argument("--methods", nargs="+", default=["anova", "l1"])
    rank.add_argument("--n-rank", type=int, default=500)
    rank.add_argument("--jobs", type=int, default=4)
    args = parser.parse_args()
    if args.command == "download":
        download(args.raw_dir)
    elif args.command == "prepare":
        print(json.dumps(prepare(args.raw_dir, args.out_dir), indent=2))
    elif args.command == "rank":
        rankings = rank_folds(
            args.data_dir, args.output_dir, args.config, args.methods, args.n_rank, args.jobs
        )
        print(f"Saved {len(rankings)} ranked genes to {args.output_dir}/rankings.tsv")
    else:
        print(evaluate(args.data_dir, args.output_dir, args.config).to_string(index=False))


if __name__ == "__main__":
    main()
