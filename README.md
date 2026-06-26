# hodmock

HOD mock カタログ生成パイプライン。[HODDIES](https://github.com/HodorFM/HODDIES) を用いて Uchuu シミュレーションから ELG mock を作成する。

## 概要

PFS などのサーベイターゲット銀河数密度 `ng` に合わせた HOD パラメータ（Ac, As）を求め、HODDIES で mock カタログを生成する。

現在実装済みのアプローチ:

| アプローチ | 説明 | 状態 |
|---|---|---|
| **rescale** | ベースライン HOD の ng を計算し、PFS ターゲット ng に合うよう Ac/As を一括スケール | 実装済み |
| **mcmc** | DESI 等の xi からパラメータを MCMC フィット | 未実装 |

## インストール / 環境

HODDIES と halomodel_module はそれぞれ特定の conda 環境を必要とする。各モジュールは必要な関数内で遅延インポートされているため、環境を分けて使用できる。

```
hoddies 環境  — HODDIES, mpytools, CAMB
               → mock 生成 / ng 計算に必要
```

```bash
# パッケージを編集可能モードでインストール（hodmock/ 直下で）
pip install -e .
```

または `sys.path` に `/home/honke/code` を追加する方法でも動作する。

## ディレクトリ構成

```
hodmock/
├── config.py               # HodMockConfig — 全設定値を一元管理
├── params/
│   └── rescale.py          # Ac/As rescaling ロジック
├── mock/
│   ├── maker.py            # HODDIES ラッパー（YAML 生成 + mock 生成）
│   └── template/
│       └── Uchuu_ELG_template.yaml   # HODDIES YAML テンプレート
└── scripts/
    └── run_rescale.py      # rescale アプローチの CLI ランナー
```

## 設定 — `HodMockConfig`

全パラメータは `config.py` の `HodMockConfig` dataclass で管理している。変えたいフィールドだけ指定すればよい。

```python
from hodmock.config import HodMockConfig
from dataclasses import replace

# デフォルト（Uchuu + PFS ELG 設定）
cfg = HodMockConfig()

# z_list だけ変える（別のハローカタログを使うとき）
cfg = HodMockConfig(z_list=[0.5, 0.8, 1.0, 1.5])

# 複数フィールドを上書き
cfg = replace(cfg, nthreads=8, mock_outdir=Path("/scratch/mock"))
```

### 主要フィールド

| フィールド | デフォルト | 説明 |
|---|---|---|
| `z_list` | Uchuu スナップショット 17 点 | 対象赤方偏移リスト |
| `z_to_halodir` | `{z: "halodir_???"}` | z → ハローカタログサブディレクトリ名 |
| `pfs_bins` | 7 z ビン (0.6–2.4) | PFS ELG z ビン定義 |
| `pfs_ng` | PFS ターゲット ng | ターゲット銀河数密度 [(Mpc/h)^-3] |
| `hod_bins` | `[(0.8,1.1), (1.1,1.6)]` | rescaling 用 HOD パラメータ bin |
| `hod_ac_base` | `[0.1, 0.1]` | ベースライン Ac（ng 計算基準） |
| `hod_as_base` | `[0.38, 0.47]` | ベースライン As（HODDIES 絶対振幅） |
| `halobase` | `/data/PFS/Uchuu/RockstarExtendedM200c1e11` | ハローカタログのベースディレクトリ |
| `mock_outdir` | `/home/honke/data/HOD_mock` | mock 出力先 |
| `nthreads` | `27` | 並列読み込みスレッド数 |
| `seed` | `42` | 乱数シード |
| `tracer` | `"ELG"` | トレーサー名（HODDIES に渡す） |

> **注意**: `z_to_halodir` の `"halodir_???"` は実際の Uchuu スナップショット番号に要更新。現在 z=0.94 (`halodir_034`) のみ設定済み。

## HODDIES の As 規約について

HODDIES では `Ac` と `As` は独立した**絶対振幅**。halomodel_module の `As` は比率（`Ns = Ac * As_ratio * f_sat`）とは異なる。

rescaling での計算式:

```
scale = ng_pfs / ng_hod
Ac    = hod_ac_base * scale
As    = hod_as_base * scale    ← HODDIES YAML に書く値
```

halomodel_module で検証する場合:

```
As_ratio = As / Ac    ← halomodel_module に渡す値
```

## 使い方

### CLI（推奨）

```bash
# Ac/As パラメータ計算 + YAML 保存のみ（CAMB が必要、mock は生成しない）
python -m hodmock.scripts.run_rescale --params-only

# 特定の z だけ mock を生成
python -m hodmock.scripts.run_rescale --z 0.94 1.03

# 全スナップショットの mock を生成
python -m hodmock.scripts.run_rescale

# z_list やスレッド数を上書き
python -m hodmock.scripts.run_rescale --z-list 0.94 1.03 1.12 --nthreads 8 --mock-outdir /tmp/mock_test
```

### Python API

```python
from hodmock.config import HodMockConfig
from hodmock.params.rescale import compute_all_params, save_params
from hodmock.mock.maker import make_mock, run_all

cfg = HodMockConfig()

# Step 1: Ac/As を rescale
params = compute_all_params(cfg)   # {z: {"Ac", "As", "density", ...}}

# Step 2: パラメータを YAML に保存
save_params(params, cfg)           # → hodmock/hod_params_rescale.yaml

# Step 3: 特定の z の mock を生成
output = make_mock(0.94, params[0.94], cfg)   # → Path(".../HODmock_ELG_z0.94.fits")

# または全スナップショット一括
outputs = run_all(params, cfg)     # {z: Path}
```

### 個別関数

```python
from hodmock.params.rescale import compute_baseline_ng, rescale_params
from hodmock.mock.maker import generate_yaml, load_hcat

# ベースライン ng だけ計算
ng_dict = compute_baseline_ng(cfg)   # {z: ng}

# 1 点だけ rescale
r = rescale_params(0.94, ng_dict, cfg)
# → {"Ac": ..., "As": ..., "density": ..., "scale": ..., ...}

# HODDIES YAML だけ生成（mock は作らない）
yaml_path = generate_yaml(0.94, r, cfg)

# ハローカタログだけ読み込む
hcat = load_hcat(cfg.halodir(0.94), cfg.mass_cut, cfg.nthreads)
```

## テンプレート YAML の編集

`mock/template/Uchuu_ELG_template.yaml` が HODDIES に渡す設定ファイルのテンプレート。
以下 4 つのプレースホルダが `generate_yaml()` で置き換えられる。

| プレースホルダ | 置き換え後 |
|---|---|
| `__AC__` | Ac（rescale 後） |
| `__AS__` | As（rescale 後、HODDIES 絶対振幅） |
| `__DENSITY__` | ng_pfs |
| `__Z__` | スナップショット赤方偏移 |

## 今後の予定

- `params/mcmc.py` — DESI xi からの MCMC フィット
- `z_to_halodir` の全スナップショット番号を埋める
- PBS / qsub ジョブ投入スクリプト（`scripts/submit_pbs.py`）
