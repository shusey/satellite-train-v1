# 保存データと再現手順

## 計算を再利用する方針

基準1倍の30/60日物理run、既存Kd解析の残留状態、同じ外乱の全LOW軌道を再利用する。新しい物理計算は0.3/3倍だけ。出力先が存在する場合は上書きを拒否するので、再現するときは新しい親ディレクトリを指定する。

今回の保存先は`outputs/kp_residual_20261010`。巨大な時系列・終了状態NPZと全局所JSONはローカル保存し、Gitには定義・設定・集計・イベント・図を保存する。

既存のKd解析スクリプトへ`--gain-axis kp`を追加して再利用した。省略時は従来のKd解析として動く。以下はリポジトリルートで実行する。`python`は今回`.analysis-venv/Scripts/python.exe`を使用。

## 30日までの新規計算

```powershell
python scripts/run_kp_case.py --kp-factor 0.3 --days 30 --output outputs/kp_new/runs/kp_0.3_30d
python scripts/run_kp_case.py --kp-factor 3 --days 30 --output outputs/kp_new/runs/kp_3_30d
```

CLIは新規計算の0.3/3倍だけを許す。凍結24時間設定からcontroller.kp_s2_invと観測時間だけを変更し、従来のsimulate_caseを実行する。Kp=1倍を重複計算しない。

## 保存runからの後処理

```powershell
python scripts/analyze_kd_residual.py --gain-axis kp --run outputs/kp_new/runs/kp_0.3_30d --all-low outputs/long_horizon_q_20261008/runs/low_30d --output outputs/kp_new/analysis/kp_0.3_30d
python scripts/analyze_kd_residual.py --gain-axis kp --run outputs/long_horizon_q_20261008/runs/central_30d --all-low outputs/long_horizon_q_20261008/runs/low_30d --cached-residual outputs/kd_residual_20261010/analysis_release/kd_1_30d --output outputs/kp_new/analysis/kp_1_30d
python scripts/analyze_kd_residual.py --gain-axis kp --run outputs/kp_new/runs/kp_3_30d --all-low outputs/long_horizon_q_20261008/runs/low_30d --output outputs/kp_new/analysis/kp_3_30d
python scripts/summarize_kd_residual.py --gain-axis kp --analysis outputs/kp_new/analysis/kp_0.3_30d outputs/kp_new/analysis/kp_1_30d outputs/kp_new/analysis/kp_3_30d --output outputs/kp_new/comparison_30d
```

baselineのcached-residualはKd=1倍解析を指定して、その保存済みe/v/P/D終了状態を再利用する。解析CLIは元run・設定・期間・D一定判定許容差の一致を検証する。新しい指標である起動時の位置誤差だけは保存物理状態から計算する。

全LOW比較の一致検査はKp以外の全設定と時間格子・外乱フラグを確認する。無効な制御器のKpは軌道へ作用しないので、この基準軌道を共有できる。

## 結果を確認後の継続

30日比較を先に確認してから、必要な条件だけ次を実行する。30日までの状態を再積分しない。

```powershell
python scripts/continue_observation_case.py --source outputs/kp_new/runs/kp_0.3_30d --days 60 --output outputs/kp_new/runs/kp_0.3_60d
python scripts/continue_observation_case.py --source outputs/kp_new/runs/kp_3_30d --days 60 --output outputs/kp_new/runs/kp_3_60d
```

継続関数は観測時間以外の変更を拒否する。長半径・位相・姿勢状態・保留指令・dwell時計を復元する。60日後処理は上の30dを60dへ変更する。

## 図・小さい成果物の保存

plot_manifest.jsonは各条件のrun・解析・全LOW参照先を指定する。既存の日本語フォントを使い、依存追加は不要だった。図の描画だけ間引き、数値指標は全10秒サンプルを使用する。

```powershell
python scripts/plot_kd_residual.py --gain-axis kp --manifest docs/kp_residual_analysis/results/plot_manifest.json --output outputs/kp_new/figures
python scripts/export_kd_results.py --gain-axis kp --comparison outputs/kp_new/comparison --output outputs/kp_new/curated_results
```

Git用local_cycle_samples.csvは各衛星の最初/最後3回だけであり、全分布ではない。全局所CSV・P/D起動CSV・帰還判定CSV・符号付き終了状態NPZ・物理時系列NPZは新規outputsに保存し、local_files.jsonから参照できる。全配列図の再生成には元NPZが必要。

## 回帰・入力・非破壊確認

```powershell
python -m pytest -q
python scripts/check_baseline_regression.py --existing-suite outputs/long_horizon_q_20261008/regression_24h/suite --output outputs/kp_check_new/regression
python scripts/check_kp_baseline_reuse.py --old outputs/kd_residual_20261010/analysis_release --new outputs/kp_residual_20261010/analysis --output outputs/kp_check_new/baseline_reuse
python scripts/check_kd_inputs.py --gain-axis kp --runs outputs/kp_residual_20261010/runs --output outputs/kp_check_new/inputs
python scripts/verify_long_horizon_preservation.py --manifest outputs/kp_residual_20261010/preserved_inputs_sha256.json --runs outputs/kp_residual_20261010/runs --output outputs/kp_check_new/preservation
```

24時間は保存済み再実行suiteと凍結referenceの全イベント・主要配列・集計値を照合する。別途短時間hybridテストで省略Kpと明示同値Kpの物理結果一致を確認する。新しい解析が既存30/60日baselineの残留各行・要約・全終了状態配列を変えないことも照合する。

作業前の保護manifestは既存879ファイルのSHA-256。完了時に欠落・変更がないこと、新規30→60日の全11配列とイベントprefix一致を確認する。mainへmergeしない。
