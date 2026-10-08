# 再現・継続・保存

すべてリポジトリ直下から実行する。依存は既存pyprojectと同じで、追加なし。
各`--output`には未使用のディレクトリを指定する。スクリプトは既存の保存先を拒否する。
今回の新規runは`outputs/long_horizon_q_20261008/`に保存し、旧runは読み取りのみ。

## 最小の後処理

```powershell
python scripts/analyze_long_horizon.py --run outputs/evaluation_20261007/runs/central_14d --all-low outputs/evaluation_20261007/runs/low_14d --output outputs/q_analysis_new/14d
```

metrics、波一覧、P/D分解、トリガー表、局所cycle間隔、7日区間統計、settling感度、
7–10日の10 s間隔の詳細データを保存する。物理計算は行わない。

## 再計算を避ける継続

```powershell
python scripts/continue_observation_case.py --source outputs/evaluation_20261007/runs/central_14d --output outputs/q_runs_new/central_21d --days 21
python scripts/continue_observation_case.py --source outputs/evaluation_20261007/runs/low_14d --output outputs/q_runs_new/low_21d --days 21
python scripts/analyze_long_horizon.py --run outputs/q_runs_new/central_21d --all-low outputs/q_runs_new/low_21d --output outputs/q_analysis_new/21d
```

次の段階はsourceを21日、outputを新しい30日ディレクトリ、daysを30へ変える。
必要性を判断してから同じ手順で45日、60日へ進む。制御と全LOWを対で継続する。
各runのNPZは元の0時刻からの全時系列を含み、continuation.jsonに継続元を残す。
条件が異なる保存runをつなげることはできない。
今回の14/21/30/45/60日データが残っていれば、同じ条件を再計算せず読み出す。

Gitのクローンだけには巨大なNPZは含まれない。保存14日runがない別環境の場合だけ、
下記で継続元を新しく用意する。今回の作業ではこの14日再計算を行っていない。

```powershell
python scripts/run_observation_case.py --hours 336 --output outputs/bootstrap_new/central_14d
python scripts/run_observation_case.py --hours 336 --case all_low_drag_baseline --output outputs/bootstrap_new/low_14d
```

その場合は継続のsourceをbootstrap_new内へ変え、解析の
`--reference-14d outputs/bootstrap_new/central_14d`も明示する。
24時間の新規確認は、保存suiteがない場合に限りcheck_baseline_regression.pyの
`--existing-suite`を省略して実行できる。

## 図

```powershell
python scripts/plot_long_horizon.py --run outputs/q_runs_new/central_21d --all-low outputs/q_runs_new/low_21d --analysis outputs/q_analysis_new/21d --output outputs/q_figures_new
```

日本語フォントYu Gothic/Meiryo/Noto Sans CJK JP/IPAexGothicのいずれかを使う。
基本4枚、45日を超えるrunでは後半q図を加えた5枚の合成図に、
長時間誤差・速度・追加損失、イベントと空間パターン、波の傾向、qの原因をまとめる。
7–10日図は既存第2波を詳細表示するもので、長時間データの間引きでcommandを見落とさない。
元時系列と詳細NPZはいずれも10 s出力。command判定も10 s間隔。

## 検証

```powershell
python -m pytest tests/test_continuation.py tests/test_long_horizon.py -q
python -m pytest -q
python scripts/check_baseline_regression.py --existing-suite outputs/long_horizon_q_20261008/regression_24h/suite --output outputs/q_regression_check_new
python scripts/verify_long_horizon_preservation.py --manifest outputs/long_horizon_q_20261008/preserved_inputs_sha256.json --runs outputs/long_horizon_q_20261008/runs --output outputs/q_preservation_check_new
```

今回の24時間確認は、旧3ケースの0–20500 sを再利用し、制御ケースの衛星12が
LOW→HIGHのslew中である20500 sから24時間まで実際に継続した。
固定fixtureと集計値・全イベント・主要時系列を比較した。短時間ユニットテストでは
制御あり外乱、制御あり無外乱、全LOW外乱について、全保存配列と全イベントを
途中再開なしの実行と比較する。遅延・slew・dwellの復元も別途確認する。

`results/regression_24h.json`に今回の判定、`results/preservation.json`に旧442ファイルの
ハッシュ検証と継続prefix検証を保存する。固定24時間fixtureは変更しない。
NPZの巨大な全履歴は新規outputs内のローカル保存、Gitには小さな集計・イベント・図を保存する。
