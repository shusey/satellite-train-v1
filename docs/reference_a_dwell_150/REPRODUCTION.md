# A-ref D300/D150：7日介入実験の再現

## 必要な既存入力

- `outputs/reference_dispersion_20261012/batch_30d/`：D300=member0、対応all LOW=member3。先頭7日を読み取る。
- `docs/reference_dispersion/data/profiles.json`：固定A-ref profile。乱数を再生成しない。
- `configs/reference_dispersion.yaml`：active環境・weather・cache methodへの参照。
- `configs/real_environment_applicability/active.yaml`、`docs/real_environment_applicability/data/active_drivers.csv`。
- `outputs/reference_dispersion_20261012/zero_regression_24h/suite/`：以前に生成した24時間回帰suite。

初回物理実行はD150のみ。既に `completion.json` があるrunを重複計算しない。再現目的で再計算する際も、存在しない新しい保存先を明示する。

```powershell
$env:PYTHONUTF8='1'
.analysis-venv/Scripts/python.exe -X utf8 -m pytest tests/test_reference_a_dwell.py tests/test_pulse_diagnostics.py -q
.analysis-venv/Scripts/python.exe -X utf8 scripts/run_reference_a_dwell.py --output outputs/reference_a_dwell_reproduce/D150_7d
.analysis-venv/Scripts/python.exe -X utf8 scripts/analyze_reference_a_dwell.py --root outputs/reference_a_dwell_reproduce
.analysis-venv/Scripts/python.exe -X utf8 scripts/verify_reference_a_dwell.py --root outputs/reference_a_dwell_reproduce --output outputs/reference_a_dwell_reproduce/verified
.analysis-venv/Scripts/python.exe -X utf8 scripts/plot_reference_a_dwell.py --input outputs/reference_a_dwell_reproduce/analysis --output outputs/reference_a_dwell_reproduce/figures
.analysis-venv/Scripts/python.exe -X utf8 -m pytest -q
```

`run`は7日・150秒だけを許可するrunner。HIGH/LOW両側のdwellが変わる。既存density cacheは軌道追従の計算法を再利用するので、分岐後のD150密度値そのものは既存軌道からコピーしない。

## 保存物

- `D150_7d/`：元10秒a、lambda、density、area、mode、metadata、normalized config、environment、input difference、completion、再開可能なcheckpoint。
- `analysis/*_metrics.json`：全metrics。粒度の異なる残留速度統計は別キー。
- `analysis/*_cycles.csv`：局所cycle、パルス量、LOW→次HIGH、D寄与、打切り。
- `analysis/*_incident_gaps.csv`：符号付きe/v、増減・反転、e×v前後、寄与収支。
- `analysis/*_cycle_stages.csv`：HIGH開始/HIGH完了/LOW指令/LOW完了の接続gap状態。
- `analysis/*_four_stage_states.npz`：上記4時点の全20gapのe/v、全21機のqP/qD。配列はJSONへ重複しない。
- `analysis/*_endpoint_states.npz`：既存residual解析のLOW復帰・全LOW境界状態。
- `analysis/*_quiet_intervals.csv`：全LOW期間・既存D一定予測。
- 公開集計の `data/*_quiet_diagnostics.csv`：実際の次起動衛星について開始時/起動時のqP/qDを明示。
- `analysis/paired_initial_pulses.csv`：0秒に同時起動した7機だけの対応比較。
- `analysis/*_satellites.csv`：衛星別cycles、遷移、HIGH時間、総高度損失、同機all LOW比損失、最終高度。
- `analysis/*_history_60s.csv` と `*_plot_series.npz`：描画用だけ60秒へ間引く。指標集計には使用しない。
- `data/comparison_summary.csv`、`comparison.json`、`verification.json`：主要表、全集計、検証。

CSV内の辞書・配列は既存の出力規約に従う。Pythonの辞書表記で保存されている列は、必要なら `ast.literal_eval` で読み、`eval` を使わない。通常の集計はJSON側を参照する。

## 変更ファイル

- `scripts/run_reference_a_dwell.py`：同一入力を確認しD150の1本を実行。
- `scripts/analyze_reference_a_dwell.py`：既存D300/all LOWとD150の全件比較。既存residual・metrics・寄与分解を再利用。
- `scripts/verify_reference_a_dwell.py`：物理収支、4段階状態、回帰、入力非破壊確認と小規模成果の保存。
- `scripts/plot_reference_a_dwell.py`：日本語図5枚。
- `src/sat_string/dwell_intervention_metrics.py`：分布・観測打切りの後処理のみ。
- `tests/test_reference_a_dwell.py`：新規3件。
- 本レポート、再現手順、`data/`、`figures/`。

既存の力学・制御・積分器・姿勢state machineは変更しない。291件テスト合格、保存24時間suiteの再照合も合格。今回の新規物理計算は7日D150のみ（372.60秒）。既存run、旧Reference profile、過去実験は変更しない。
