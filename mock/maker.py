"""
maker.py — HODDIES を使って mock カタログを生成する。

使い方:
    from hodmock.params.rescale import compute_all_params
    from hodmock.mock.maker import make_mock
    from hodmock.config import HodMockConfig

    cfg    = HodMockConfig()
    params = compute_all_params(cfg)          # {z: {Ac, As, density, ...}}
    output = make_mock(0.94, params[0.94], cfg)  # → Path to .fits
"""

from __future__ import annotations

import glob
import multiprocessing
import time
from functools import partial
from pathlib import Path

import h5py
import numpy as np

from hodmock.config import DEFAULT_CONFIG, HodMockConfig

# HODDIES / mpytools は hoddies conda 環境でのみ利用可能なため遅延インポートする。
# load_hcat / make_mock 内で必要時にインポートする。

# ─────────────────────────────────────────────────────────────────────────────
# YAML 生成（テンプレート → per-z HODDIES YAML）
# ─────────────────────────────────────────────────────────────────────────────

def generate_yaml(
    z: float,
    hod_params: dict,
    cfg: HodMockConfig = DEFAULT_CONFIG,
    yaml_dir: Path | None = None,
) -> Path:
    """HODDIES YAML テンプレートに Ac/As/density/z を埋め込んで per-z YAML を生成する。

    テンプレート内のプレースホルダ:
        __AC__      → hod_params["Ac"]
        __AS__      → hod_params["As"]
        __DENSITY__ → hod_params["density"]
        __Z__       → z

    Parameters
    ----------
    z : float
    hod_params : dict
        rescale_params() または同等の関数が返す dict。
        必須キー: "Ac", "As", "density"
    cfg : HodMockConfig
    yaml_dir : Path, optional
        生成 YAML の保存ディレクトリ。省略時は cfg.hod_params_yaml.parent / "yaml_per_z"。

    Returns
    -------
    Path
        生成された YAML ファイルのパス。
    """
    if yaml_dir is None:
        yaml_dir = cfg.hod_params_yaml.parent / "yaml_per_z"
    yaml_dir = Path(yaml_dir)
    yaml_dir.mkdir(parents=True, exist_ok=True)

    template_path = cfg.template_yaml
    if not template_path.exists():
        raise FileNotFoundError(f"テンプレート YAML が見つからない: {template_path}")

    tmpl = template_path.read_text()
    out = (
        tmpl
        .replace("__AC__",      str(hod_params["Ac"]))
        .replace("__AS__",      str(hod_params["As"]))
        .replace("__DENSITY__", str(hod_params["density"]))
        .replace("__Z__",       str(z))
    )

    out_path = yaml_dir / f"Uchuu_HOD_z{z:.2f}.yaml"
    out_path.write_text(out)
    return out_path


# ─────────────────────────────────────────────────────────────────────────────
# ハローカタログ読み込み（Uchuu / Rockstar h5 形式）
# ─────────────────────────────────────────────────────────────────────────────

def _read_one_h5(fname: str, mass_cut: float) -> dict:
    """1 つの h5 ファイルからハローデータを読み込む（並列ワーカー用）。"""
    cols = ["x", "y", "z", "vx", "vy", "vz",
            "M200c", "Rvir", "rs", "vrms", "id", "upid"]
    with h5py.File(fname, "r") as f:
        d = {c: f[c][:] for c in cols}

    mask = (d["upid"] == -1) & (np.log10(d["M200c"]) > mass_cut)
    Mh = d["M200c"][mask]
    Rh = d["Rvir"][mask]
    Rs = d["rs"][mask]
    return {
        "x":        d["x"][mask],
        "y":        d["y"][mask],
        "z":        d["z"][mask],
        "vx":       d["vx"][mask],
        "vy":       d["vy"][mask],
        "vz":       d["vz"][mask],
        "Mh":       Mh,
        "Rh":       Rh,
        "Rs":       Rs,
        "c":        Rh / Rs,
        "Vrms":     d["vrms"][mask],
        "halo_id":  d["id"][mask],
        "log10_Mh": np.log10(Mh),
    }


def load_hcat(
    halodir: Path | str,
    mass_cut: float,
    nthreads: int,
):
    """Uchuu ハローカタログ（halolist_*.h5）を並列で読み込み mpytools.Catalog を返す。

    Parameters
    ----------
    halodir : Path or str
        ハローカタログ h5 ファイルが入ったディレクトリ。
    mass_cut : float
        log10(M200c / [Msun/h]) の下限カット。
    nthreads : int
        並列読み込みのスレッド数。

    Returns
    -------
    mpytools.Catalog
    """
    from mpytools import Catalog  # hoddies 環境でのみ利用可能

    halodir = Path(halodir)
    fnames = sorted(glob.glob(str(halodir / "halolist_*.h5")))
    if not fnames:
        raise FileNotFoundError(f"h5 ファイルが見つからない: {halodir}/halolist_*.h5")

    print(f"[hcat] {len(fnames)} files in {halodir}", flush=True)
    t0 = time.time()

    fn = partial(_read_one_h5, mass_cut=mass_cut)
    with multiprocessing.Pool(min(nthreads, len(fnames))) as pool:
        results = pool.map(fn, fnames)

    merged = {k: np.concatenate([r[k] for r in results]) for k in results[0]}
    print(
        f"[hcat] done in {time.time() - t0:.1f}s  |  "
        f"halos: {len(merged['Mh']):,}  |  "
        f"log10_Mh: {merged['log10_Mh'].min():.2f}–{merged['log10_Mh'].max():.2f}",
        flush=True,
    )
    return Catalog.from_dict(merged)


# ─────────────────────────────────────────────────────────────────────────────
# mock 生成
# ─────────────────────────────────────────────────────────────────────────────

def make_mock(
    z: float,
    hod_params: dict,
    cfg: HodMockConfig = DEFAULT_CONFIG,
    hcat=None,
    yaml_dir: Path | None = None,
) -> Path:
    """指定した z の mock カタログを生成して保存する。

    ① HODDIES YAML を生成（generate_yaml）
    ② ハローカタログを読み込む（load_hcat）— hcat を渡した場合はスキップ
    ③ HODDIES HOD を初期化して make_mock_cat を実行
    ④ カタログを cfg.mock_output(z) に保存

    Parameters
    ----------
    z : float
    hod_params : dict
        rescale_params() の戻り値。必須キー: "Ac", "As", "density"
    cfg : HodMockConfig
    hcat : mpytools.Catalog, optional
        事前に読み込んだハローカタログ。None の場合は自動で読み込む。
    yaml_dir : Path, optional
        生成 YAML の保存先（generate_yaml に渡す）。

    Returns
    -------
    Path
        保存した mock カタログのパス。
    """
    from HODDIES import HOD  # hoddies 環境でのみ利用可能

    print(f"=== make_mock  z={z:.2f} ===", flush=True)

    # ① HODDIES YAML 生成
    yaml_path = generate_yaml(z, hod_params, cfg, yaml_dir=yaml_dir)
    print(f"[yaml] {yaml_path}", flush=True)

    # ② ハローカタログ読み込み
    if hcat is None:
        hcat = load_hcat(cfg.halodir(z), cfg.mass_cut, cfg.nthreads)

    # ③ HODDIES HOD 初期化 + mock 生成
    print("[HOD] initializing ...", flush=True)
    hod = HOD(param_file=str(yaml_path), hcat_file=hcat)

    print(f"[HOD] making mock (tracer={cfg.tracer}) ...", flush=True)
    t0 = time.time()
    cat = hod.make_mock_cat(tracers=cfg.tracer, fix_seed=cfg.seed, verbose=True)
    print(f"[HOD] done in {time.time() - t0:.1f}s", flush=True)

    # サマリ
    n_cen = int(cat["Central"].sum())
    n_sat = int((1 - cat["Central"]).sum())
    print(
        f"[result] N={len(cat):,}  cen={n_cen:,}  sat={n_sat:,}  "
        f"sat_frac={n_sat / len(cat):.3f}",
        flush=True,
    )

    # ④ 保存
    output = cfg.mock_output(z)
    output.parent.mkdir(parents=True, exist_ok=True)
    print(f"[save] {output}", flush=True)
    cat.write(str(output))

    return output


def run_all(
    all_params: dict[float, dict],
    cfg: HodMockConfig = DEFAULT_CONFIG,
    yaml_dir: Path | None = None,
) -> dict[float, Path]:
    """全スナップショットの mock を順次生成する。

    Parameters
    ----------
    all_params : dict[float, dict]
        compute_all_params() の戻り値 {z: hod_params}。
    cfg : HodMockConfig
    yaml_dir : Path, optional

    Returns
    -------
    dict[float, Path]
        {z: 保存した mock カタログのパス}
    """
    outputs: dict[float, Path] = {}
    skipped: list[float] = []

    for z, params in all_params.items():
        # Z_TO_HALODIR が未設定（"halodir_???"）の z はスキップ
        halodir_name = cfg.z_to_halodir.get(z, "")
        if "???" in halodir_name:
            print(f"[run_all] skip z={z:.2f}: halodir 未設定", flush=True)
            skipped.append(z)
            continue

        halodir = cfg.halobase / halodir_name
        if not halodir.exists():
            print(f"[run_all] skip z={z:.2f}: {halodir} が存在しない", flush=True)
            skipped.append(z)
            continue

        outputs[z] = make_mock(z, params, cfg, yaml_dir=yaml_dir)

    if skipped:
        print(f"[run_all] skipped z = {skipped}", flush=True)

    return outputs
