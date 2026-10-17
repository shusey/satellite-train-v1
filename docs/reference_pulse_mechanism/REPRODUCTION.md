# Referenceパルス診断の再現

## 入力

物理モデル・設定・乱数profileを変更せず、以下を読み取る。

- `outputs/reference_dispersion_20261012/batch_30d/` と `batch_60d/`：元10秒配列・checkpointイベント。
- `docs/reference_dispersion/data/profiles.json`：固定21機profile。
- `docs/reference_dispersion/data/60d/`：既存の局所cycle・全LOW予測・集計。
- `outputs/reference_dispersion_20261012/zero_regression_24h/suite/`：以前に生成したゼロ分散回帰suite。

大配列はGit対象外。存在しない環境では今回のコードが勝手に再シミュレーションや再取得をすることはない。

## 実行

リポジトリ直下、既存解析用Python環境を使う。出力先は存在しない新しいディレクトリを指定する。既存の成果へ上書きするためのオプションはない。

```powershell
$env:PYTHONUTF8='1'
.analysis-venv/Scripts/python.exe -X utf8 -m pytest tests/test_pulse_diagnostics.py -q
.analysis-venv/Scripts/python.exe -X utf8 scripts/analyze_reference_pulses.py --output outputs/reference_pulse_reproduce
.analysis-venv/Scripts/python.exe -X utf8 scripts/verify_reference_pulses.py --input outputs/reference_pulse_reproduce --output outputs/reference_pulse_reproduce_verified
.analysis-venv/Scripts/python.exe -X utf8 scripts/plot_reference_pulses.py --input outputs/reference_pulse_reproduce --output outputs/reference_pulse_reproduce_figures
.analysis-venv/Scripts/python.exe -X utf8 -m pytest -q
```

`verify`は既存24時間データを再照合する。新規物理計算ではない。新しい解析モジュールは既存シミュレーションからimportされない。

## 出力の定義

- `*_pulse_cycles.csv`：全16,983サイクルの局所前後速度、HIGH換算時間、qP/qD、自己/近傍/自然抗力の寄与、次の同機HIGHまでの寄与。未観測の次回イベントは空欄。
- `*_incident_gaps.csv`：接続間隔ごとの符号付きe、v、半径差、積算速度寄与。独立標本とは見なさない。
- `*_quiet_intervals.csv`：全機LOW期間のP蓄積、D変化、B差寄与、既存予測との比較。
- `*_first_pulse_trace.csv` 等：最初／30日以降最初の局所cycle・全LOW区間、前後600秒を含む10秒時系列。
- `summary.json`、`interval_summary.csv`：全件集計。
- `examples.json`：事前定義した代表区間。`own_excess_g_impulse`は固定基準半径上の速度寄与。
- `input_manifest.json`：入力のハッシュまたはサイズ・時刻。旧データは変更しない。

`next_qD_neighbor_control`は対象機がLOWの間の近傍操作の寄与。`next_delta_qP`と区別する。全LOW区間ではこの操作寄与はゼロであることを検査する。

## 追加したファイル

- `src/sat_string/pulse_diagnostics.py`：観測軌道に沿う加法分解と積分。
- `scripts/analyze_reference_pulses.py`：保存runの全件診断。
- `scripts/verify_reference_pulses.py`：収支、求積、24時間保存suite、非破壊の確認。
- `scripts/plot_reference_pulses.py`：日本語図5枚。
- `tests/test_pulse_diagnostics.py`：新規5テスト。
- `docs/reference_pulse_mechanism/REPORT.md`、本書、`data_final/`、`figures_final/`：成果。

288件合格時のテスト実行時間68.22秒。検証・図生成は上記コマンドで別途確認済み。mainへのmerge、ゲインやdwell変更、追加物理runは行っていない。
