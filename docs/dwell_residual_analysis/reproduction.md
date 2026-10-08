# 再現と保存データの扱い

## 環境

Python 3.11、既存のNumPy/Matplotlib/PyYAML/pytestを使用。依存追加なし。
図は環境のYu Gothicを使用し日本語ラベルと単位を付ける。
以下はリポジトリ直下で実行。ローカルではpythonを
`.analysis-venv/Scripts/python.exe`に読み替える。
**出力先が既に存在する場合は処理を停止するので、再実行時は新しい名前を指定する。**
旧runや中間解析を削除・置換して実行しない。

## 再利用したデータ

- 300 s 30/60日：outputs/long_horizon_q_20261008/runs/central_30d、central_60d
- 対応する全LOW：同じ親のlow_30d、low_60d
- 既存300秒集計・群区分：outputs/long_horizon_q_20261008/analysis/30d、60d
- 共通初期LOW：outputs/long_horizon_q_20261008/regression_24h/suite/controlled_disturbed
- 群区分の既定値：outputs/evaluation_20261007/runs/central_14d のイベント

300秒・全LOWの物理再計算はしていない。24時間の前回実再実行済みsuiteを再照合。
150/600秒は初回指令20430 sの1制御周期前、20420 sまで共通の保存配列を利用した。
この時点まで全LOW、指令なし、初期next_allowed=0で、dwellはまだ一度も作用していない。
それ以降は新しい設定による既存continue_caseで計算。
60日は各30日状態からdurationだけを延長した。

## 新規30日、必要性を確認した後の60日

```powershell
python scripts/run_dwell_case.py --dwell 150 --days 30 --output outputs/dwell_residual_20261008/runs/dwell_150_30d
python scripts/run_dwell_case.py --dwell 600 --days 30 --output outputs/dwell_residual_20261008/runs/dwell_600_30d
python scripts/continue_observation_case.py --source outputs/dwell_residual_20261008/runs/dwell_150_30d --days 60 --output outputs/dwell_residual_20261008/runs/dwell_150_60d
python scripts/continue_observation_case.py --source outputs/dwell_residual_20261008/runs/dwell_600_30d --days 60 --output outputs/dwell_residual_20261008/runs/dwell_600_60d
```

両延長の判断は30日最大誤差・活動数の悪化と、観測終了後の横断予測による。
150秒の最後の全LOW境界から33.406日、600秒では30.200日が次回指令の予測。
最初から60日を無条件に実行していない。

## 保存データの後処理

```powershell
python scripts/analyze_dwell_residual.py --run outputs/dwell_residual_20261008/runs/dwell_150_30d --all-low outputs/long_horizon_q_20261008/runs/low_30d --output outputs/dwell_residual_20261008/analysis_final/dwell_150_30d
python scripts/analyze_dwell_residual.py --run outputs/long_horizon_q_20261008/runs/central_30d --all-low outputs/long_horizon_q_20261008/runs/low_30d --cached-analysis outputs/long_horizon_q_20261008/analysis/30d --output outputs/dwell_residual_20261008/analysis_final/dwell_300_30d
python scripts/analyze_dwell_residual.py --run outputs/dwell_residual_20261008/runs/dwell_600_30d --all-low outputs/long_horizon_q_20261008/runs/low_30d --output outputs/dwell_residual_20261008/analysis_final/dwell_600_30d
```

60日も上記30dを60dへ置換。300秒だけ既存metrics/wavesを読んで集計を再利用する。
`--relative-D-tolerance`は後処理の品質設定。変更してもシミュレーション再実行は不要。
metrics.jsonの旧指標・コアは変更しておらず、新しいresidual_metrics.json/CSVに追加指標を分離した。

```powershell
python scripts/summarize_dwell_residual.py --analysis outputs/dwell_residual_20261008/analysis_final/dwell_150_30d outputs/dwell_residual_20261008/analysis_final/dwell_300_30d outputs/dwell_residual_20261008/analysis_final/dwell_600_30d outputs/dwell_residual_20261008/analysis_final/dwell_150_60d outputs/dwell_residual_20261008/analysis_final/dwell_300_60d outputs/dwell_residual_20261008/analysis_final/dwell_600_60d --output outputs/dwell_residual_20261008/comparison_final
python scripts/plot_dwell_residual.py --manifest docs/dwell_residual_analysis/results/plot_manifest.json --output outputs/dwell_residual_20261008/figures_release
python scripts/check_baseline_regression.py --existing-suite outputs/long_horizon_q_20261008/regression_24h/suite --output outputs/dwell_residual_20261008/regression_24h
python scripts/check_dwell_preservation.py --output outputs/dwell_residual_20261008/checks_final
python -m pytest -q
```

## Git成果物とローカル成果物

Gitへ保存するのは評価コード、テスト、文書、比較CSV/JSON、局所・群・全LOW境界CSV、
集計・設定・イベント・検証結果・図。
results.npz、cycle_end_states.npz、旧runと試作中間解析はローカルに保持。
最終解析はanalysis_final、最終比較はcomparison_final、最終図はfigures_release。
長時間曲線は描画時だけ約340 sへ間引く。集計値・指令イベント・パルス拡大図は元の10 s標本を使用。
CSV end_state_rowは同じ解析ディレクトリのcycle_end_states.npzを参照する。

NPZは元々Git対象外なので、クローンだけでは既存長時間データを得られない。
保存NPZを入手できない環境のみ、凍結24時間suiteを
`check_baseline_regression.py --output <新しいパス>`で生成し、
300秒も`continue_observation_case.py`で30→60日へ再構成する必要がある。
新規dwell実行の`--source`には、その24時間controlled_disturbedを指定できる。
群区分参照14日データもなければ、先に同じ300秒基準で生成する。
これは新しい環境での再現手順であり、本実験では既存300秒を再計算していない。

## 変更ファイル

- src/sat_string/residual.py：終了状態、3種境界、P/D予測、有効性、実パルス幅。
- scripts/run_dwell_case.py：物理コアを変えず、初期共通LOWからdwell介入。
- scripts/analyze_dwell_residual.py：既存集計を利用する後処理、CSV/NPZ出力。
- scripts/summarize_dwell_residual.py：比較表と分布・予測統計。
- scripts/plot_dwell_residual.py：7枚の日本語比較図。
- scripts/check_dwell_preservation.py：旧ファイルhash、条件、初期/延長prefix照合。
- tests/test_residual.py：14件の指標・予測・境界・設定・再利用検査。
- docs/dwell_residual_analysis/：定義、結果、手順、図、保存可能な小規模成果物。
