"""
config.py — hodmock パイプライン全体の設定値を一元管理する。

ここにある値を変更すれば、パイプライン全体に反映される。
"""

from pathlib import Path

import numpy as np

# ─────────────────────────────────────────────────────────────────────────────
# Uchuu スナップショット赤方偏移リスト
# ─────────────────────────────────────────────────────────────────────────────
Z_LIST: list[float] = [
    0.63, 0.70, 0.78, 0.86, 0.94, 1.03,
    1.12, 1.22, 1.32, 1.43, 1.54, 1.65,
    1.77, 1.90, 2.03, 2.17, 2.31,
]

# z → ハローカタログディレクトリ名の対応表
# （実際のスナップショット番号を確認して埋めること）
Z_TO_HALODIR: dict[float, str] = {
    0.63: "halodir_???",  # TODO
    0.70: "halodir_???",  # TODO
    0.78: "halodir_???",  # TODO
    0.86: "halodir_???",  # TODO
    0.94: "halodir_034",
    1.03: "halodir_???",  # TODO
    1.12: "halodir_???",  # TODO
    1.22: "halodir_???",  # TODO
    1.32: "halodir_???",  # TODO
    1.43: "halodir_???",  # TODO
    1.54: "halodir_???",  # TODO
    1.65: "halodir_???",  # TODO
    1.77: "halodir_???",  # TODO（Uchuu にこの z があるか要確認）
    1.90: "halodir_???",  # TODO（同上）
    2.03: "halodir_???",  # TODO（同上）
    2.17: "halodir_???",  # TODO（同上）
    2.31: "halodir_???",  # TODO（同上）
}

# ─────────────────────────────────────────────────────────────────────────────
# PFS ELG ターゲット銀河数密度
# ─────────────────────────────────────────────────────────────────────────────
# 単位: (Mpc/h)^-3
PFS_BINS: list[tuple[float, float]] = [
    (0.6, 0.8),
    (0.8, 1.0),
    (1.0, 1.2),
    (1.2, 1.4),
    (1.4, 1.6),
    (1.6, 2.0),
    (2.0, 2.4),
]

PFS_NG: np.ndarray = np.array([
    5.4e-4,
    11.2e-4,
    10.6e-4,
    14.4e-4,
    10.8e-4,
    4.3e-4,
    4.8e-4,
])

# ─────────────────────────────────────────────────────────────────────────────
# HOD パラメータ bin（rescaling 用）
# ─────────────────────────────────────────────────────────────────────────────
# この 2 bin に対して別々のベースライン HOD を定義する。
# 範囲外の z は端の bin に外挿する（rescale.py の extrapolate オプション参照）。
HOD_BINS: list[tuple[float, float]] = [
    (0.8, 1.1),   # z1 bin
    (1.1, 1.6),   # z2 bin
]

# ─────────────────────────────────────────────────────────────────────────────
# HOD ベースラインパラメータ
# ─────────────────────────────────────────────────────────────────────────────
# halomodel_module での ng 計算に使うベースライン HOD。
# Ac=0.1 のときの ng が基準値 ng_hod になる。
# HOD_BINS の各 bin に対応する dict のリスト。
HOD_BASELINE: list[dict] = [
    {   # z1 bin (0.8 < z < 1.1)
        "model":      "mhmq",
        "logMc":      11.62,
        "sigma_M":    0.21,
        "Ac":         0.1,
        "gamma":      6.49,
        "logM0":      11.10,
        "logM1":      13.0,
        "alpha":      0.74,
        "As":         0.38,
        "poff":       0.0,
        "Roff":       0.0,
        "conformity": True,
        "sat_model":  "elg",
    },
    {   # z2 bin (1.1 < z < 1.6)
        "model":      "mhmq",
        "logMc":      11.62,
        "sigma_M":    0.21,
        "Ac":         0.1,
        "gamma":      6.49,
        "logM0":      11.10,
        "logM1":      13.0,
        "alpha":      0.74,
        "As":         0.38,
        "poff":       0.0,
        "Roff":       0.0,
        "conformity": True,
        "sat_model":  "elg",
    },
]

# HODDIES の Ac/As ベースライン値（HOD_BINS に対応）。
# HODDIES では Ac と As が独立した絶対振幅。ng に合わせて両者を同率でスケールする。
#   Ac_scaled = HOD_AC_BASE * scale
#   As_scaled = HOD_AS_BASE * scale  （HODDIES YAML にそのまま書く値）
# halomodel_module 検証用の「絶対 As」は HOD_AS_BASE * HOD_AC_BASE で求まる。
HOD_AC_BASE: np.ndarray = np.array([0.1,  0.1 ])
HOD_AS_BASE: np.ndarray = np.array([0.38, 0.47])  # HODDIES baseline As

# ─────────────────────────────────────────────────────────────────────────────
# Halo model 計算グリッド
# ─────────────────────────────────────────────────────────────────────────────
K_ARR: np.ndarray = np.logspace(-4, 3, 1000)  # [h/Mpc]
M_ARR: np.ndarray = np.logspace(10, 15, 200)  # [Msun/h]

# ─────────────────────────────────────────────────────────────────────────────
# シミュレーション / ジョブ設定
# ─────────────────────────────────────────────────────────────────────────────
MASS_CUT: float = 11.0  # log10(M200c / [Msun/h]) の下限カット
NTHREADS: int   = 27    # ハローカタログ読み込みの並列数
SEED: int       = 42
TRACER: str     = "ELG"

# ─────────────────────────────────────────────────────────────────────────────
# パス設定
# ─────────────────────────────────────────────────────────────────────────────
# Uchuu ハローカタログのベースディレクトリ（サブディレクトリは Z_TO_HALODIR で決まる）
HALOBASE: Path = Path("/data/PFS/Uchuu/RockstarExtendedM200c1e11")

# mock カタログの出力先
MOCK_OUTDIR: Path = Path("/home/honke/data/HOD_mock")

# HODDIES に渡す YAML テンプレート
TEMPLATE_YAML: Path = Path(__file__).parent / "mock" / "template" / "Uchuu_ELG_template.yaml"

# rescaling で生成するコンパクト HOD パラメータ YAML の保存先
HOD_PARAMS_YAML: Path = Path(__file__).parent / "hod_params_rescale.yaml"
