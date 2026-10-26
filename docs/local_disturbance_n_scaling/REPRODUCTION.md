# 再現手順

リポジトリルートでPython 3.11以上とプロジェクト依存パッケージを使用する。
今回の背景密度はexponentialなのでpymsisの再取得・宇宙天気データ取得は不要。

```powershell
$env:PYTHONUTF8='1'
.analysis-venv/Scripts/python.exe -u scripts/run_local_disturbance_n_scaling.py
.analysis-venv/Scripts/python.exe scripts/analyze_local_disturbance_n_scaling.py --publish docs/local_disturbance_n_scaling
.analysis-venv/Scripts/python.exe -m pytest -q > outputs/local_disturbance_n_scaling_20261026/tests_verified.txt
.analysis-venv/Scripts/python.exe scripts/check_baseline_regression.py --output outputs/local_disturbance_n_scaling_20261026/regression_24h
.analysis-venv/Scripts/python.exe scripts/verify_local_disturbance_n_scaling.py --destination outputs/local_disturbance_n_scaling_20261026/verification
```

runnerは既存の完了runを再利用し、不完全runは保持して停止する。
上書き・削除・暗黙の再実行は行わない。
解析は完了済みのmetricsを再利用する。
`--publish`は出力先の`data/`が存在すると停止する。
既に公開済みの成果を作り直す場合は、別の新規ディレクトリを指定する。
回帰確認先も既存の場合は新しい名前を指定する。

入力を変更して別実験を行う場合は設定ファイルを複製し、runnerのSPEC/OUTPUTを
その新規実験へ明示的に変更する。固定結果をそのまま再利用して別条件と称さない。

## ローカル保存

`outputs/local_disturbance_n_scaling_20261026/Nxxx/`の各Nについて：

- time_s.npy：保存時刻、10 s。
- a_m.npy、lambda_rad.npy、density_kg_m3.npy、area_m2.npy、mode.npy：
  (時刻, 2ケース, N)。ケース0=nominal、1=disturbed。
- metadata.json：全入力、node入力仕様、個体差ゼロのprofile。
- checkpoint.json：終了状態、全イベント、node入力仕様。既存prefixは変更しない。
- normalized_config.yaml、input_checks.json、completion.json。
- analysis/raw_gap_error_m.npy：(時刻, 2ケース, N-1)。
- analysis/delta_gap_error_m.npy：(時刻, N-1)、nominal差分。
- analysis/gaps.csv：左右、graph distance、最大追加誤差と時刻、初回1 m到達。
- analysis/spatial.csv：左右別peak、平均、最大値。
- analysis/metrics.json：全体・終盤指標、切替/HIGH量、時間窓別評価。
- analysis/history_60s.csv：図用の代表時系列。指標計算は元の10 s全時刻を使用。

raw NPZ/NPYはGitに追加しない。Gitには設定、集計、イベント、図、文書、テストを保存。
公開用時系列CSVはさらに600 s間隔へ間引き、Git容量を抑える。
この間引きはpeak計算・図作成へ適用していない。

## 検証

- baselineの`formation.satellite_count`、`simulation.duration_s`、
  `disturbance.enabled`以外が不変であることを毎runでassert。
- nominal/disturbedのinitial stateとB倍率は同一、どちらもcontrol ON。
- node入力なし・ゼロ入力の結果一致、単位・符号、波形積分値を単体検証。
- 既存Gaussianを用いた切替を含む区間で、batchと元simulatorの
  全保存軌道・密度・面積・mode・イベントが一致することを検証。
- 既存24 h software referenceを新しく3ケース再実行して照合。
- 開始前に既存runs/outputsのサイズ・更新時刻と凍結fixtureのSHA256を保存し、
  終了後に不変を検証。今回ディレクトリは照合対象から除く。

ブランチ：`astra/local-disturbance-n-scaling`。mainへmergeしない。
正確なcommitは`git rev-parse HEAD`で確認でき、完了時の会話にも記載する。
