# データ再利用・実行・検証

## 保存場所

今回の新規runは `outputs/kd_residual_20261010/runs/kd_{0.3,3}_{30,60}d`。baselineは既存 `outputs/long_horizon_q_20261008/runs/central_{30,60}d`、全LOW基準は同じ親ディレクトリの `low_{30,60}d`。Kd=1倍の物理計算は行っていない。

正式な解析は `outputs/kd_residual_20261010/analysis_release/kd_{0.3,1,3}_{30,60}d`。初期検証用のanalysis/analysis_finalは正式成果ではない。正式解析のsigned endpoint配列は `cycle_end_states.npz`。元の全時系列は各runの `results.npz`。

巨大NPZはローカルに保存し、Gitには条件・イベント・集計・定義・図・再現手順を保存する。解析provenanceには元runと全LOW、設定、再利用したbaseline解析、群の分類間隔を記録する。Git上の小さい成果だけから全時系列図を再生成する場合は、元NPZを用意するか該当する新規runを再実行する必要がある。

すべての実行先は新しいディレクトリとする。スクリプトは既存出力の上書きを拒否する。

## 新規Kd run

リポジトリルートで実行。以下の `python` は依存がインストールされたPython 3.11以上（今回 `.analysis-venv/Scripts/python.exe`）。

```powershell
python scripts/run_kd_case.py --kd-factor 0.3 --days 30 --output outputs/kd_new/runs/kd_0.3_30d
python scripts/run_kd_case.py --kd-factor 3 --days 30 --output outputs/kd_new/runs/kd_3_30d
```

baseline=1を誤って再計算しないよう、run_kd_caseのCLIは0.3/3倍のみ受け付ける。最初の外乱応答からKdに依存するため、制御済みbaselineのprefixは借用しない。積分・制御・姿勢state machineは既存simulate_caseをそのまま使用。

30日結果を解析・確認してから、必要な条件だけ保存状態を継続する。

```powershell
python scripts/continue_observation_case.py --source outputs/kd_new/runs/kd_0.3_30d --output outputs/kd_new/runs/kd_0.3_60d --days 60
python scripts/continue_observation_case.py --source outputs/kd_new/runs/kd_3_30d --output outputs/kd_new/runs/kd_3_60d --days 60
```

今回は0.3倍で次回切替予測が30.06日を示し、3倍で30日末にも多数の切替が続いていたため、双方を延長した。既存30日prefixを再積分せず、長半径・位相・姿勢状態・dwell時計・pendingイベントを復元する。duration以外の変更は継続関数が拒否する。

## 後処理

```powershell
python scripts/analyze_kd_residual.py --run outputs/kd_new/runs/kd_0.3_30d --all-low outputs/long_horizon_q_20261008/runs/low_30d --output outputs/kd_new/analysis/kd_0.3_30d
python scripts/analyze_kd_residual.py --run outputs/long_horizon_q_20261008/runs/central_30d --all-low outputs/long_horizon_q_20261008/runs/low_30d --cached-residual outputs/dwell_residual_20261008/analysis_final/dwell_300_30d --output outputs/kd_new/analysis/kd_1_30d
python scripts/analyze_kd_residual.py --run outputs/kd_new/runs/kd_3_30d --all-low outputs/long_horizon_q_20261008/runs/low_30d --output outputs/kd_new/analysis/kd_3_30d
python scripts/summarize_kd_residual.py --analysis outputs/kd_new/analysis/kd_0.3_30d outputs/kd_new/analysis/kd_1_30d outputs/kd_new/analysis/kd_3_30d --output outputs/kd_new/comparison_30d
```

60日の解析も同じ引数の30dを60dへ変更する。Kdだけが異なるcontroller無効の全LOW物理結果は共用できる。解析CLIは他の全条件と外乱フラグ・時間格子の一致を検証する。baselineでは既存dwell=300秒の残留解析・endpoint NPZも再利用する。既存パラメータの科学的な妥当性を仮定するものではない。

予測品質設定は `--relative-D-tolerance 0.1`、P/D分類設定は `--dominance-ratio 2`。物理入力ではなく後処理設定。予測許容差を変更する場合、異なる設定のcacheは拒否するため `--cached-residual` を外す（物理再計算は不要）。

## 図

results/plot_manifest.jsonは3条件の60日run・解析・全LOWを指す。日本語フォントは環境に既存のYu Gothicを使用し、新しい依存は追加していない。

```powershell
python scripts/plot_kd_residual.py --manifest docs/kd_residual_analysis/results/plot_manifest.json --output outputs/kd_new/figures
```

図だけ保存点を約2万点へ間引く。数値集計は全10秒サンプルを使用。既存のstate metrics、HIGH equivalent time、波の外側フィット、constant-D予測を再利用する。

## 回帰・テスト・非破壊確認

```powershell
python -m pytest -q
python scripts/check_baseline_regression.py --existing-suite outputs/long_horizon_q_20261008/regression_24h/suite --output outputs/kd_check_new/regression
python scripts/check_kd_metrics_compatibility.py --run outputs/long_horizon_q_20261008/runs/central_60d --all-low outputs/long_horizon_q_20261008/runs/low_60d --cached-residual outputs/dwell_residual_20261008/analysis_final/dwell_300_60d --output outputs/kd_check_new/metrics
python scripts/check_kd_inputs.py --runs outputs/kd_residual_20261010/runs --output outputs/kd_check_new/inputs
python scripts/verify_long_horizon_preservation.py --manifest outputs/kd_residual_20261010/preserved_inputs_sha256.json --runs outputs/kd_residual_20261010/runs --output outputs/kd_check_new/preservation
```

24時間については保存済みの実再実行suiteを凍結referenceと照合し、物理計算を重複しない。追加の単体テストで暗黙Kdと明示的な同値Kdの短いhybrid simulationが全配列・イベント一致することを確認する。

頻繁な切替で旧残留後処理の全時系列/全イベント走査が高コストになったため、検索をソート済み時刻とsearchsorted、局所区間slice、衛星cycle番号の辞書へ置き換えた。計算式・サンプル選択・境界の包含条件は変更していない。既存60日baselineの旧残留JSON各行・要約・全endpoint配列との一致を別途確認する。

保護manifestは作業開始時点の既存runs/outputs全725ファイルのSHA-256。最終確認では欠落・変更なしと、各30→60日の全11配列/イベントprefix一致を確認する。入力検査ではKdと観測時間以外が凍結条件と同じことを確認する。mainへmergeしない。
