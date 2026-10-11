# 再現・データ再利用

物理core、既存設定・runは変更していない。Pythonは既存の `.analysis-venv`、pymsis 0.12.0 / NRLMSIS **2.0** を使用する。新しい依存の追加はない。作業branchは `astra/real-environment-applicability`。元の物理coreは `e759c99a66cd15b5d3dfd6b0d676a064d246f4f7` と同じ。

## 保存場所

- 公開する小さい成果：このディレクトリの `data/`、`figures/`、報告書。実際のF10.7/Ap入力は `data/quiet_drivers.csv`、`active_drivers.csv`、`strong_drivers.csv`。
- ローカル大配列：`outputs/real_environment_20261011/`。新しい60日runは `runs/active_60d` と `runs/strong_60d`。
- 再利用した静穏期：`outputs/nrlmsis2_20261010/runs/all_low_60d`。再計算していない。
- 再利用するGRACE観測：`outputs/observed_density_20261010/raw/`。新規ダウンロードしていない。
- その実際の2003年F10.7/Ap入力も `data/grace_2003_drivers.csv` へ保存した。背景NRLMSIS比較はこの固定系列を使う。
- 宇宙天気の生データ：既存 `outputs/nrlmsis2_20261010/SW-All.csv` と、新しく一度取得したNASA `outputs/real_environment_20261011/raw/omni2_2024.dat`。

物理run・解析の出力は新しい場所だけに作成する。既に完成した出力を再実行で置き換えない。大きいNPZ、生観測アーカイブ、年度OMNIファイルはGitへ含めない。別マシンで完全な時系列を再現する場合は、元runを移すか、保存したdriverを用いて計算する。既存runがあるマシンでは再計算しない。

## 今回実行した主要手順

以下はリポジトリrootから実行する。出力が既に存在する場合は再実行しない。

```powershell
$env:PYTHONUTF8='1'
.analysis-venv/Scripts/python.exe scripts/prepare_environment_comparison.py
.analysis-venv/Scripts/python.exe scripts/run_nrlmsis_case.py --environment configs/real_environment_applicability/active.yaml --weather outputs/real_environment_20261011/weather/active/drivers.csv --case all_low --duration-s 5184000 --cache-seconds 120 --output outputs/real_environment_20261011/runs/active_60d
.analysis-venv/Scripts/python.exe scripts/run_nrlmsis_case.py --environment configs/real_environment_applicability/strong.yaml --weather outputs/real_environment_20261011/weather/strong/drivers.csv --case all_low --duration-s 5184000 --cache-seconds 120 --output outputs/real_environment_20261011/runs/strong_60d
.analysis-venv/Scripts/python.exe scripts/analyze_environment_comparison.py
.analysis-venv/Scripts/python.exe scripts/analyze_grace_integrated_difference.py
.analysis-venv/Scripts/python.exe scripts/plot_environment_comparison.py
.analysis-venv/Scripts/python.exe scripts/write_environment_report.py
```

`prepare_environment_comparison.py` は間隔誤差を一切読まず、宇宙天気だけで選定する。元CSV snapshotのSHA256と、全候補・除外理由・採用期間は `data/weather_selection.json` に保存した。後から別の期間を選び直すことはしない。

driver CSVは前日observed F10.7、当日のcentered81 F10.7、7要素のApを3時間ごとに記録する。平均への置換や予測値での補完は行わない。supportに含まれる2024-07-30のF10.7=400.7は今回のMay1–June30の直接入力には使われない。実際に使うdaily/centered fluxだけを適正範囲チェックする。July/September候補は必要supportにINT値もあり不採用となる。

## 密度補間の検証（物理runを再計算しない）

```powershell
.analysis-venv/Scripts/python.exe scripts/validate_environment_density_integrals.py --environment configs/real_environment_applicability/strong.yaml --weather outputs/real_environment_20261011/weather/strong/drivers.csv --output outputs/real_environment_20261011/strong_density_convergence.json
.analysis-venv/Scripts/python.exe scripts/audit_environment_integrated_cache.py --run outputs/nrlmsis2_20261010/runs/all_low_60d --environment configs/real_environment_applicability/quiet.yaml --weather outputs/real_environment_20261011/weather/quiet/drivers.csv --output outputs/real_environment_20261011/quiet_integrated_cache_audit.json
```

前者は最大3時間apから選んだ約3軌道の**指定位置上の密度だけ**を直接2秒評価と比較する。後者は保存済み10秒位置で直接密度を再評価し、密度キャッシュ誤差の積算差を確認する。別の60日物理計算ではない。軌道へのfeedbackを再伝播していないため、厳密な60日物理解の誤差保証とは区別する。

同じ監査をactive/strongの新規保存runにも実施した。`--run`、`--environment`、`--weather`、`--output` のquiet部分を各環境の保存先へ合わせる。監査の追加基準が不合格なら、JSONを残して非ゼロ終了する。これは集計値や回帰を都合よく書き換える機能ではない。

新しいマシンでは、公開driverを `--weather docs/real_environment_applicability/data/active_drivers.csv` 等として物理runnerへ直接渡せる。後処理scriptは上記の保存run位置を期待するため、元のrunを移すか、そのマシンで不足分だけ作成する。元の宇宙天気snapshotがない場合、現在の配布ファイルに差し替えて凍結期間選定を再現したと主張しない。観測rawの取得方法・URL・ハッシュは既存 `docs/observed_density_design/reproduction.md` と本研究の `data/grace_integrated_summary.json` を参照する。

## テスト・回帰・旧成果の検査

```powershell
.analysis-venv/Scripts/python.exe -m pytest -q
.analysis-venv/Scripts/python.exe scripts/check_baseline_regression.py --existing-suite outputs/relative_motion_20261011/regression_24h/suite --output outputs/real_environment_20261011/regression_24h
.analysis-venv/Scripts/python.exe scripts/check_saved_data_preservation.py --manifest outputs/real_environment_20261011/preserved_sha256.json --output outputs/real_environment_20261011/preservation_check
```

回帰は保存済みの再計算suiteを再利用する。新規処理は物理coreを呼んで改変するものではなく、元coreのハッシュも保持検査に含めている。回帰では主要時系列、全イベント、集計値を凍結referenceと照合する。初期の保存データ一覧8404ファイルのSHA256を比較する。mainへのmergeは行わない。
