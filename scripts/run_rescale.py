#!/usr/bin/env python
"""
run_rescale.py — rescale アプローチによる end-to-end mock 生成スクリプト。

処理フロー:
    1. compute_all_params()  … PFS ng に合うよう Ac/As を rescale
    2. save_params()         … hod_params_rescale.yaml を書き出し
    3. run_all()             … 全スナップショットの mock を生成・保存
       （または --params-only で 1-2 のみ実行）

使い方:
    # パラメータ計算 + YAML 保存のみ（mock は生成しない）
    python -m hodmock.scripts.run_rescale --params-only

    # 特定の z だけ mock を生成
    python -m hodmock.scripts.run_rescale --z 0.94 1.03

    # 全スナップショット mock を生成
    python -m hodmock.scripts.run_rescale

    # z_list / nthreads など Config を上書き
    python -m hodmock.scripts.run_rescale --z-list 0.94 1.03 1.12 --nthreads 8

    # 出力先を変える
    python -m hodmock.scripts.run_rescale --mock-outdir /tmp/mock_test
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import replace
from pathlib import Path


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="rescale アプローチで HOD mock を生成する",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--params-only", action="store_true",
        help="Ac/As の rescaling + YAML 保存だけ行い、mock は生成しない",
    )
    p.add_argument(
        "--z", nargs="+", type=float, metavar="Z", default=None,
        help="mock を生成するスナップショット z（省略時は cfg.z_list 全体）",
    )
    p.add_argument(
        "--z-list", nargs="+", type=float, metavar="Z", default=None,
        help="rescaling の対象 z を上書きする（--z と区別するため）",
    )
    p.add_argument(
        "--nthreads", type=int, default=None,
        help="ハローカタログ並列読み込みのスレッド数",
    )
    p.add_argument(
        "--mock-outdir", type=Path, default=None,
        help="mock カタログの出力先ディレクトリ",
    )
    p.add_argument(
        "--halobase", type=Path, default=None,
        help="ハローカタログのベースディレクトリ",
    )
    p.add_argument(
        "--yaml-dir", type=Path, default=None,
        help="per-z HODDIES YAML の保存先（省略時は hod_params_yaml の隣の yaml_per_z/）",
    )
    p.add_argument(
        "--seed", type=int, default=None,
        help="乱数シード",
    )
    p.add_argument(
        "--no-verbose", action="store_true",
        help="進捗表示を抑制する",
    )
    return p.parse_args()


def main() -> None:
    args = _parse_args()

    # ── Config 構築 ──────────────────────────────────────────────────────────
    from hodmock.config import HodMockConfig

    cfg = HodMockConfig()

    overrides: dict = {}
    if args.z_list is not None:
        overrides["z_list"] = args.z_list
    if args.nthreads is not None:
        overrides["nthreads"] = args.nthreads
    if args.mock_outdir is not None:
        overrides["mock_outdir"] = Path(args.mock_outdir)
    if args.halobase is not None:
        overrides["halobase"] = Path(args.halobase)
    if args.seed is not None:
        overrides["seed"] = args.seed
    if overrides:
        cfg = replace(cfg, **overrides)

    verbose = not args.no_verbose

    # ── Step 1-2: rescaling + YAML 保存 ──────────────────────────────────────
    from hodmock.params.rescale import compute_all_params, save_params

    print("=" * 60, flush=True)
    print("Step 1/2: Ac/As rescaling", flush=True)
    print("=" * 60, flush=True)

    all_params = compute_all_params(cfg, verbose=verbose)
    yaml_path  = save_params(all_params, cfg)

    print(f"\nhod_params saved → {yaml_path}", flush=True)

    if args.params_only:
        print("\n--params-only が指定されたため mock 生成をスキップ。", flush=True)
        return

    # ── Step 3: mock 生成 ─────────────────────────────────────────────────────
    from hodmock.mock.maker import make_mock, run_all

    print("\n" + "=" * 60, flush=True)
    print("Step 2/2: mock 生成", flush=True)
    print("=" * 60, flush=True)

    if args.z is not None:
        # 指定した z だけ生成
        missing = [z for z in args.z if z not in all_params]
        if missing:
            print(f"[警告] 以下の z は rescaling の対象外（z_list に含まれていない）: {missing}", flush=True)
            print("  --z-list も合わせて指定してください。", flush=True)
            sys.exit(1)

        for z in args.z:
            make_mock(z, all_params[z], cfg, yaml_dir=args.yaml_dir)
    else:
        # cfg.z_list 全体
        run_all(all_params, cfg, yaml_dir=args.yaml_dir)

    print("\n完了。", flush=True)


if __name__ == "__main__":
    main()
