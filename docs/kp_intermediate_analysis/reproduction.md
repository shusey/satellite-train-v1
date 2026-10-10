# Kp中間2条件の保存先と再現方法

## 再利用するもの

基点は前回Kp解析branchの`510c19234d321b90ff579efb3551edbf29438b0e`。
今回の出力先は`outputs/kp_intermediate_20261010`。既存の1倍・3倍の物理runと後処理結果、全LOW軌道を再利用する。

- 1倍物理run：`outputs/long_horizon_q_20261008/runs/central_{30,60}d`
- 3倍物理run：`outputs/kp_residual_20261010/runs/kp_3_{30,60}d`
- 1倍・3倍後処理：`outputs/kp_residual_20261010/analysis/kp_{1,3}_{30,60}d`
- 旧比較：`outputs/kp_residual_20261010/comparison/comparison.json`
- 全LOW：`outputs/long_horizon_q_20261008/runs/low_{30,60}d`

巨大なNPZはローカル保存。出力ディレクトリは既存なら拒否するので、再現時は以下の`kp_mid_new`を新しい名前にする。
以下はリポジトリルートで実行。今回のPythonは`.analysis-venv/Scripts/python.exe`。

## 30日を先に実行

```powershell
python scripts/run_kp_case.py --kp-factor 2 --days 30 --output outputs/kp_mid_new/runs/kp_2_30d
python scripts/run_kp_case.py --kp-factor 2.5 --days 30 --output outputs/kp_mid_new/runs/kp_2.5_30d
python scripts/check_kd_inputs.py --gain-axis kp --factors 2 2.5 --runs outputs/kp_mid_new/runs --output outputs/kp_mid_new/inputs_30d
```

凍結24時間設定からKpと観測時間だけを変更。物理core、制御則、積分器、保存間隔は従来のまま。

## 新規2条件だけを後処理

```powershell
python scripts/analyze_kd_residual.py --gain-axis kp --run outputs/kp_mid_new/runs/kp_2_30d --all-low outputs/long_horizon_q_20261008/runs/low_30d --output outputs/kp_mid_new/analysis/kp_2_30d
python scripts/analyze_kd_residual.py --gain-axis kp --run outputs/kp_mid_new/runs/kp_2.5_30d --all-low outputs/long_horizon_q_20261008/runs/low_30d --output outputs/kp_mid_new/analysis/kp_2.5_30d
python scripts/summarize_kd_residual.py --gain-axis kp --analysis outputs/kp_mid_new/analysis/kp_2_30d outputs/kp_mid_new/analysis/kp_2.5_30d --output outputs/kp_mid_new/new_30d
python scripts/compare_kp_tradeoff.py --comparison outputs/kp_residual_20261010/comparison/comparison.json outputs/kp_mid_new/new_30d/comparison.json --days 30 --output outputs/kp_mid_new/comparison_30d
```

1倍・3倍の既存比較行を読み取り、新規2条件の行と結合する。旧解析を再計算しない。全4条件の設定がKp・観測時間以外は凍結入力と一致することも結合時に検査する。

## 結果確認後だけ60日へ継続

```powershell
python scripts/continue_observation_case.py --source outputs/kp_mid_new/runs/kp_2_30d --days 60 --output outputs/kp_mid_new/runs/kp_2_60d
python scripts/continue_observation_case.py --source outputs/kp_mid_new/runs/kp_2.5_30d --days 60 --output outputs/kp_mid_new/runs/kp_2.5_60d
```

既存30日状態を復元し、その部分を再積分しない。60日後処理は上の30dを60dへ変更。
`compare_kp_tradeoff.py`には旧比較・新規30日比較・新規60日比較を渡し、`--days 30 60`で全8行を結合する。

## 5図と小さい成果物

```powershell
python scripts/plot_kp_tradeoff.py --comparison outputs/kp_mid_new/comparison/comparison.json --output outputs/kp_mid_new/figures
python scripts/export_kd_results.py --gain-axis kp --factors 2 2.5 --comparison outputs/kp_mid_new/comparison --output outputs/kp_mid_new/curated_results
```

図は新規5枚だけ。日本語フォントは既存環境を使用し、依存追加なし。
1倍・3倍の既存成果ファイルは前回ディレクトリを参照し、今回の設定・イベント・残留CSVの書き出しは新規2条件だけにする。比較表には全4条件を残す。
図の時系列は描画用に間引くが、集計は全10秒保存サンプルを使う。
Git用の局所往復CSV見本は衛星ごとに最初/最後3回のみ。全データの所在は各条件の`local_files.json`に残す。

## 回帰と非破壊確認

```powershell
python -m pytest -q
python scripts/check_baseline_regression.py --existing-suite outputs/long_horizon_q_20261008/regression_24h/suite --output outputs/kp_mid_check/regression
python scripts/verify_long_horizon_preservation.py --manifest outputs/kp_intermediate_20261010/preserved_inputs_sha256.json --runs outputs/kp_intermediate_20261010/runs --output outputs/kp_mid_check/preservation
```

保存済み24時間suiteを凍結referenceと照合し、短時間hybridテストで省略Kpと同値明示Kpを比較する。新たにbaselineの長時間軌道を計算しない。
作業開始時に既存`outputs/`・`runs/`・`docs/`の1,555ファイルのSHA-256を保存した。完了時に変更・欠落と30→60日の全配列・イベントprefixを検査する。
