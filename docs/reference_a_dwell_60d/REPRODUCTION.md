# A-ref D150の60日観測

既存7日runを先頭区間として保存したまま、同じstate・姿勢・イベント履歴から再開する。D300と同一profileの全LOWは再計算しない。

## 実行手順

リポジトリ直下の既存解析環境を使用する。すべて新規保存先へ書く。既存の完了済みsuffixは再計算しない。

```powershell
$env:PYTHONUTF8='1'
.analysis-venv/Scripts/python.exe -X utf8 -m pytest tests/test_reference_a_dwell.py -q
.analysis-venv/Scripts/python.exe -X utf8 scripts/run_reference_a_dwell.py --days 30 --prefix outputs/reference_a_dwell_20261017/D150_7d --output outputs/reference_a_dwell_60d_20261017/D150_30d
.analysis-venv/Scripts/python.exe -X utf8 scripts/analyze_reference_a_dwell_horizon.py --days 30
```

30日完了後、completion、有限な状態、正のsigned gap、順序維持を確認する。以下のrunnerもcheckpointの正常性・入力一致を検査する。

```powershell
.analysis-venv/Scripts/python.exe -X utf8 scripts/run_reference_a_dwell.py --days 60 --prefix outputs/reference_a_dwell_60d_20261017/D150_30d --output outputs/reference_a_dwell_60d_20261017/D150_60d
.analysis-venv/Scripts/python.exe -X utf8 scripts/analyze_reference_a_dwell_horizon.py --days 60
.analysis-venv/Scripts/python.exe -X utf8 scripts/plot_reference_a_dwell_horizon.py
.analysis-venv/Scripts/python.exe -X utf8 -m pytest -q
```

## 必要な既存データ

- `outputs/reference_a_dwell_20261017/D150_7d/`：7日状態とcheckpoint。
- `outputs/reference_dispersion_20261012/batch_30d/` と `batch_60d/`：D300、対応all LOW。
- `outputs/reference_dispersion_20261012/A_plot_30d.npz` と `A_plot_60d.npz`：既存D300の描画用保存履歴。
- `outputs/reference_pulse_mechanism_20261017_v3/A_pulse_cycles.csv` と `A_incident_gaps.csv`：既存D300のパルス・残留速度診断。
- `docs/reference_dispersion/data/30d/`、`60d/`：既存D300指標と、各観測長で打ち切った全LOW区間。
- `docs/reference_dispersion/data/profiles.json`、既存active weather/environment、保存済み24時間回帰suite。

巨大runはGit対象外。入力がない場合、コードは勝手に別profileを生成したり、既存runを再計算したりしない。

## 接続と集計

7日・30日境界は各1時点だけを残す。軌道半径、位相、面積、姿勢modeは完全一致を検査する。densityだけ既存trajectory cache再アンカーに伴うrtol=1e-8を許容し、差を記録する。イベント履歴もprefixの完全一致を検査する。

30日・60日の解析出力はそれぞれ `analysis_30d/`、`analysis_60d/`。D150の全cycle・接続gap状態・全LOW区間・再切替・区間別コストを、7日実験と同じ指標関数で評価する。D300は既存指標と既存パルス診断を読み取る。

すべての指標は元の10秒保存値を使う。描画用CSVだけ60秒へ間引く。全RMSは全保存時点・全20gapを均等に二乗平均したもの。追加高度損失は同一profileの全LOW平均半径−制御平均半径。現在半径の物理弧長、固定基準半径の誤差、全ペアchordを分ける。

7日と30日について、後から観測した次回切替を混入させない。各観測長の終了を超える次回イベントは打切りとする。D300の全LOW予測も対応する既存30日／60日集計を使う。

## 今回の変更

- `run_reference_a_dwell.py`：同一D150の30/60日checkpoint延長を許可。
- `dwell_horizon.py`：保存区間の非破壊接続。
- `analyze_reference_a_dwell_horizon.py`：D150延長解析と既存D300再利用。
- `plot_reference_a_dwell_horizon.py`：日本語比較図5枚。
- `test_reference_a_dwell.py`：観測長だけの変更、境界重複除去・状態不一致検出を確認。

既存7日レポートの結論・保存データは更新せず、長期観測の結論は本ディレクトリのレポートに記す。
