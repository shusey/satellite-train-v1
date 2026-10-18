# A-ref HIGH/LOW dwell 2×2実験の再現

既存runを再利用する手順を優先する。同じ出力先が存在する場合、runnerは停止し、削除・上書きしない。
新規の30日物理計算は交差2条件だけ。60日へ延長する処理はない。

## 入力

- 固定A profile: `docs/reference_dispersion/data/profiles.json` の `profiles[0]`。21個の実値を使用し、再生成しない。
- 半径分散のsample std: ddof=1、0.5 m。B倍率は全衛星1。
- Reference conditionは実機母集団の観測標準偏差ではなく、工学的な感度評価条件。
- 初期半径と初期位置を区別し、初期mean longitudeは100 km間隔のまま。
- 環境: `configs/real_environment_applicability/active.yaml`、歴史的指標: `docs/real_environment_applicability/data/active_drivers.csv`。
- 2024-03-01 00:00 UTC、NRLMSIS 2.0 / pymsis 0.12.0、inclination 51.6度。
- 120秒密度キャッシュの方法・指標は共通。軌道が変われば評価位置も変わるので、他caseの密度時系列をコピーしない。

## 後方互換

`controller.dwell_high_time_s` と `controller.dwell_low_time_s` を任意指定する。
未指定の側は従来の `controller.dwell_time_s` を使用する。従来の設定ファイルと正規化結果を保つため、未指定キーは保存しない。

## 実行

リポジトリの既存 `.analysis-venv` を使用する。テスト記録用の出力親フォルダは事前に新規作成する。以下の順序は計算を初めて再現する場合だけ実施する。
既に同じrunがある場合は下の後処理へ進む。新規出力先へ変更するときは各scriptのOUTPUT定数も整合させる。

```powershell
.analysis-venv/Scripts/python.exe -m pytest -q | Tee-Object outputs/reference_dwell_factorial_20261018/tests.txt
.analysis-venv/Scripts/python.exe scripts/check_baseline_regression.py --output outputs/reference_dwell_factorial_20261018/regression_24h
.analysis-venv/Scripts/python.exe scripts/check_dwell_diagonal.py
.analysis-venv/Scripts/python.exe scripts/run_reference_dwell_factorial.py --case H150-L300
.analysis-venv/Scripts/python.exe scripts/run_reference_dwell_factorial.py --case H300-L150
```

## 保存runのみからの後処理

```powershell
.analysis-venv/Scripts/python.exe scripts/analyze_reference_dwell_factorial.py --case H300-L300
.analysis-venv/Scripts/python.exe scripts/analyze_reference_dwell_factorial.py --case H150-L300
.analysis-venv/Scripts/python.exe scripts/analyze_reference_dwell_factorial.py --case H300-L150
.analysis-venv/Scripts/python.exe scripts/analyze_reference_dwell_factorial.py --case H150-L150
.analysis-venv/Scripts/python.exe scripts/analyze_reference_dwell_factorial.py --summarize
.analysis-venv/Scripts/python.exe scripts/plot_reference_dwell_factorial.py
.analysis-venv/Scripts/python.exe scripts/publish_reference_dwell_factorial.py --output outputs/reference_dwell_factorial_20261018/published_reproduction
```

H300-L300は `outputs/reference_dispersion_20261012/batch_30d` のmember 0を使用する。
同一A-refの全LOWは同じbatchのmember 3。
H150-L150は `outputs/reference_a_dwell_20261017/D150_7d` と
`outputs/reference_a_dwell_60d_20261017/D150_30d` を結合する。
7日境界は状態・面積・modeの完全一致を確認し、重複sampleを1個除く。
既存の密度キャッシュ再始動差はrtol=1e-8以内で検査する。

全10秒sampleを評価へ使用し、図だけ60秒に間引く。
cycle終端・四段階の全状態配列は `outputs/reference_dwell_factorial_20261018/analysis/<case>/` のNPZへ保存する。
再出版先は既存のdocsを上書きしない新規フォルダにする。既存REPORTの数値は新しいcomparison.csvと照合できる。
Gitには生軌道・状態NPZを入れず、小さい表・設定・診断・図・文書を保存する。

## 検証と原資料保全

24時間の人工外乱3ケースを新規再計算し、凍結fixtureの全イベント・主要時系列・集計値と照合。
さらにHIGH/LOWを明示的に300/300または150/150とした6時間runを、保存済みscalar dwell runと全配列・イベントで完全一致照合する。
原資料はread-onlyで使用し、大ファイルはsize/mtime、小さいファイルはSHA256も実行前後で照合する。
`data/source_manifest.json` と `data/verification.json` を参照。
