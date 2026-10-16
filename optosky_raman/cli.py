"""Command-line contract: explicit modes and nonzero exit on partial failure."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys
from .pipeline import run_pipeline
from .report import generate_report


def main(argv: list[str]|None=None) -> int:
    parser=argparse.ArgumentParser(description='Auditable Optosky RamanSPy analysis')
    parser.add_argument('--mode',required=True,choices=['audit','analyze','compare','report'])
    parser.add_argument('--input',nargs='+',type=Path,help='CSV batch directories')
    parser.add_argument('--config',type=Path,help='YAML experiment/analysis settings')
    parser.add_argument('--output',required=True,type=Path,help='Output directory outside source directories')
    parser.add_argument('--results',type=Path,help='Saved analysis directory for report mode')
    args=parser.parse_args(argv)
    try:
        if args.mode=='report':
            if args.input or args.config: parser.error('report mode reads --results only; input/config would imply recomputation')
            outputs=generate_report(args.results or args.output,args.output)
            print(json.dumps(outputs,ensure_ascii=False)); return 0
        if not args.input: parser.error('--input is required for audit/analyze/compare')
        result=run_pipeline(args.input,args.output,args.config,args.mode)
        print(json.dumps({'status':result['status'],'counts':result['counts'],'output':str(args.output.resolve())},ensure_ascii=False))
        return 0 if result['status']=='complete' else 2
    except (ValueError,FileNotFoundError,OSError) as exc:
        print(f'{type(exc).__name__}: {exc}',file=sys.stderr); return 1


if __name__=='__main__': raise SystemExit(main())
