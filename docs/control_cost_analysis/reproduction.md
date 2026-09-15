# 再現手順とGit情報

## Git

- 対象: shusey/satellite-train-v1
- default branch: main
- 解析開始コミット: 9449a78 (num1)
- 作業branch: astra/control-cost-analysis
- mainへのmergeは行っていない。

主要コミット：

1. d9dfdfb — Document implemented satellite model and align NumPy requirement
2. 5b7b110 — Add event-exact drag cost metrics, bounded experiments and PD analysis
3. この成果物・再現手順を含むコミット — Record validated control cost results and horizon dependence

正確な一覧は次で取得できる。

~~~powershell
git log --oneline 9449a78..astra/control-cost-analysis
git diff --stat 9449a78..astra/control-cost-analysis
~~~

## セットアップ

リポジトリのルートで実行する。追跡済み.venvを上書きせず解析用環境を作る。

~~~powershell
git clone --branch astra/control-cost-analysis https://github.com/shusey/satellite-train-v1.git
cd satellite-train-v1
python -m venv .analysis-venv
.\.analysis-venv\Scripts\Activate.ps1
python -m pip install -r docs/control_cost_analysis/requirements-reproduction.txt
python -m pip install -e . --no-deps
python -m pytest -q
~~~

環境記録: Python 3.11.3 / Windows、実行依存の正確な版はrequirements-reproduction.txt。
coreの数値計算は元から不変。異なるOS・数値ライブラリでのbitwise一致までは保証しない。
リポジトリに元からあるnp.trapezoidのため、NumPyの最低版を2.0に修正した。

## 主実験

出力先は存在しない新しいフォルダーを指定する。既存出力は上書きしない。

~~~powershell
python scripts/run_control_cost_analysis.py --config configs/baseline.yaml --output outputs/control_cost_analysis --stage baseline
python scripts/run_control_cost_analysis.py --config configs/baseline.yaml --output outputs/control_cost_analysis --stage sweep
~~~

または --stage all で標準3ケースの成功後に7条件を実行できる。
現在設定1.5e-7は標準ケース結果を再利用し、追加計算は6条件。
物理条件はbaselineから継承し、変更するのは閾値だけ。
h_offは元の幅h_on-h_offを維持する。今回の幅は1e-8。

各ケースにはresults.npz（全時系列）、switch_events.csv、normalized_config.yaml、
既存metrics、control_cost.json、satellite_metrics.csv、gap_metrics.csvを保存する。
同一キー名で符号付き高度/長半径/エネルギー変化を読める。
JSONの未定義値はnull、CSVでは空欄。
図と集計だけ再生成する場合はシミュレーションを実行しない。

~~~powershell
python scripts/run_control_cost_analysis.py --output outputs/control_cost_analysis --stage postprocess
~~~

## 72時間の追加確認

~~~powershell
python scripts/check_control_cost_horizon.py --output outputs/control_cost_analysis
~~~

期間は既存runs/standard_suite/normalized_config.yamlのduration_s=259200から読む。
同ファイルの別のゲイン・閾値・積分刻みは採用しない。
無制御、現在設定、1.8e-7、2.0e-7の4条件だけを延長する。
出力: horizon_check/comparison.csv、各ケース全時系列、plots/11_horizon_cost.png。
72時間の追加試験も新しい出力先を必要とする。

## 変更前結果との照合・数値収束

本作業では解析追加前に以下を実行し、元コードの時系列を保存した。

~~~powershell
python scripts/run_standard_suite.py --config configs/baseline.yaml --output outputs/control_cost_before --skip-animation
~~~

現在もcoreは変更していないので同コマンドを利用できる。
次は3条件のdt=2→1 sチェックをまとめて行う。検証用サブフォルダーも既存なら停止する。

~~~powershell
python scripts/validate_control_cost_analysis.py --output outputs/control_cost_analysis --before outputs/control_cost_before --refine baseline/controlled_disturbed --refine sweep/h_on_1.380e-07 --refine sweep/h_on_2.000e-07
~~~

検証スクリプトは次を保存する：

- 元コミット9449a78からのcore差分の有無
- 変更前後3ケースの全配列の一致・最大差
- 保存済みexperiment_13_1.5の指標JSON・イベントログとの一致
- 積分刻み半減による配列差、ピーク差、コスト差、切替記録の一致

本作業では既存64テスト通過後に変更し、追加後の全72テスト、
関連する8テストの再確認、代表図の目視確認を行った。

## 成果物と再現の範囲

docs/control_cost_analysis/results と figures にレビュー用の表・設定・イベント・図・
検証証拠を格納した。元のruns内のファイルは変更していない。
大きい生時系列はローカルoutputs/control_cost_analysis（約85 MB）に保持し、
Gitには重複登録していない。上記コマンドで再生成できる。

provenance.jsonのsource_commitは最初の計算開始時のHEAD=d9dfdfbを記録している。
当時は解析コード追加中だった。最終解析コードは作業branchにコミット済みであり、
計算に用いた元のcoreはprovenanceのSHA256とvalidationの比較で確認できる。
日本語の研究上の解釈はreport.md、数式の詳細はpd_analysis.md。
