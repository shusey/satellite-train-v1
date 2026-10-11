# Git・検証・変更一覧

研究判断用サマリーは[README.md](README.md)を参照してください。ここは検証と再現用の付録です。

## 作業系

- 基点：`astra/initial-dispersion-screening`、`0e32b68d3b5eb07024873616caef91e2d844e12c`
- 新規branch：`astra/dispersion-control-transfer`
- mainへのmerge、既存branchへの強制pushは行いません。
- 最終commitは取得したbranch上で `git rev-parse HEAD`、commit一覧は `git log 0e32b68d..HEAD --oneline` で確認できます。文書に自分自身のcommit hashを埋め込むための追加commitは作りません。

## 変更一覧

すべて新規ファイルです。既存の物理モデル・制御則・積分器・設定・過去成果は編集していません。

|ファイル|目的|
|---|---|
|`src/sat_string/dispersion_control.py`|衛星固有B倍率をLOW/HIGH/slewで維持し、既存の運動方程式・制御・姿勢機械・RK4を利用して独立ケースを計算。全状態のcheckpointを保存。|
|`src/sat_string/dispersion_transfer_metrics.py`|同一配置の全LOW比較、仮想許容誤差、運用区間、残留状態、順序維持、高度コストを評価。|
|`configs/dispersion_control_transfer.yaml`|変更禁止の環境・前回配置・分散・観測時間を指定。|
|`scripts/run_dispersion_control.py`|前回代表seedを選び、新規ON3ケースと必要なAB全LOWを計算。30日レビュー後のcheckpoint延長。|
|`scripts/extend_dispersion_control_OFF.py`|既存A/B全LOWの30日までを再計算せず、30–60日のみ延長。|
|`scripts/analyze_dispersion_control.py`|保存データから既存残留速度・再横断予測を再利用し、イベント・gap・HIGH・高度・P/Dを集計。|
|`scripts/audit_dispersion_control_density.py`|新しい軌道上の補間密度をNRLMSIS2.0直接評価とスポット比較。|
|`scripts/plot_dispersion_control.py`|保存データから日本語6図を作成。|
|`scripts/report_dispersion_control.py`|保存した集計から技術付録を生成。判断サマリーは実際の結果を読んで記述。|
|`scripts/reproduce_dispersion_control.py`|巨大な旧runを持たないGit取得直後の再現専用。凍結した3配置ON/OFFのみを新規保存先へ計算。今回の追加実験としては未実行。|
|`scripts/verify_dispersion_control.py`|回帰・境界・入力一致・既存ファイル非破壊・新規追加だけのGit差分を確認。|
|`tests/test_dispersion_control.py`|B倍率、独立性、旧core一致、slew/dwellを含む再開、配置固定の11検証。|
|`tests/test_dispersion_transfer_metrics.py`|精度・高度・順序・静穏候補・閾値超過の6検証。|
|`docs/dispersion_control_transfer/`|小規模集計、凍結設定、入力配置、図、判断と再現手順。|

## 検証

- 既存256件＋新規17件＝**273件合格**。Windows子プロセスも `PYTHONUTF8=1` とし、警告なし。
- 標準24時間3ケースについて、集計・全切替イベント・主要時系列の回帰が合格。
- 新しい分散対応計算でも分散ゼロ24時間を照合。a・平均経度・高度・間隔誤差・速度・密度・面積・姿勢・イベントは旧coreと完全一致。後処理qの演算順序差だけは約1.06×10⁻²² m/s²。
- 再開境界の軌道・姿勢・面積と30日までのイベントは完全一致を確認。密度キャッシュは再初期化のため境界密度が約6.12×10⁻¹⁰の相対差となり、境界診断は相対1×10⁻⁸、絶対0を許容。既存の密度補間基準1×10⁻³より十分小さい差です。
- 作業前に8570ファイル（約9.13 GB）のSHA256を記録し、全ファイルが変更・欠損なし。機械判定は[data/verification.json](data/verification.json)、[data/preservation.json](data/preservation.json)。
- 大きな時系列NPY、終了状態NPZ、checkpointはローカル `outputs/dispersion_control_20261011/` にのみ保存し、Gitには加えません。

## 主要な実行手順

既存保存データを再利用する今回の手順と、Git取得直後の再現方法はREADMEの「再現方法」に分けて記載しています。同じ物理ケースを重複実行しないでください。

```powershell
$env:PYTHONUTF8="1"
.analysis-venv/Scripts/python.exe -X utf8 -m pytest -q
.analysis-venv/Scripts/python.exe -X utf8 scripts/check_saved_data_preservation.py --manifest outputs/dispersion_control_20261011/preserved_inputs_sha256.json --output outputs/dispersion_control_20261011/preservation_final
.analysis-venv/Scripts/python.exe -X utf8 scripts/verify_dispersion_control.py
```

最後の検証は既存runがある研究環境向けで、新規ファイルをGitへstageした状態で実行します。保存済み結果を削除・上書きして再実行する手順ではありません。
