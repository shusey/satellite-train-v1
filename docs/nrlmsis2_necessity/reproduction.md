# NRLMSIS 2.0：再現手順

## 環境

リポジトリ直下から実行する。Python3.11、NumPy2.4.6、Matplotlib3.11.2、PyYAML6.0.3、pytest9.1.1、pymsis0.12.0で検証。図にはWindowsのYu Gothicを使用し、新しい日本語フォント依存は追加しなかった。

```powershell
python -m pip install -e ".[test,nrlmsis]"
python -m pytest -q
```

毎回、未使用の出力先を指定する。既存ディレクトリの上書きを許可するオプションはない。

## 歴史的入力

Gitに保存した `docs/nrlmsis2_necessity/data/drivers.csv` を使えば再取得不要。`observed_support.csv` とprovenanceも保存する。元データを再取得する場合、公開CSVは更新されるためSHAが同じとは限らない。保存した支持データから同じ481入力を再生成できる。

```powershell
python scripts/prepare_nrlmsis_weather.py --raw docs/nrlmsis2_necessity/data/observed_support.csv --output outputs/nrl_weather_new
```

この場合provenanceの元ファイルSHAは支持データのSHAになるが、計算入力driversの内容は同じである。

## 本格計算前の短時間検証

```powershell
python scripts/run_nrlmsis_case.py --weather docs/nrlmsis2_necessity/data/drivers.csv --case all_low --duration-s 16670 --cache-seconds 0 --output outputs/nrl_direct_low_new
python scripts/run_nrlmsis_case.py --weather docs/nrlmsis2_necessity/data/drivers.csv --case baseline --duration-s 16670 --cache-seconds 0 --output outputs/nrl_direct_control_new
python scripts/validate_nrlmsis.py --weather docs/nrlmsis2_necessity/data/drivers.csv --case all_low --direct-run outputs/nrl_direct_low_new --output outputs/nrl_convergence_low_new
python scripts/validate_nrlmsis.py --weather docs/nrlmsis2_necessity/data/drivers.csv --case baseline --direct-run outputs/nrl_direct_control_new --output outputs/nrl_convergence_control_new
```

120 sを選んだ根拠は技術付録と保存した収束結果。キャッシュを0にすると直接NRLMSIS計算になる。物理のRK4刻みはどの設定でも2 sのまま。

## 60日全LOW

```powershell
python scripts/run_nrlmsis_case.py --weather docs/nrlmsis2_necessity/data/drivers.csv --case all_low --duration-s 5184000 --cache-seconds 120 --prefix outputs/nrl_convergence_low_new/all_low_120s --output outputs/nrl_low_60d_new
```

`--prefix` では観測終了時刻だけ変更でき、キャッシュ・幾何・宇宙天気・全物理／制御設定を確認する。前の1日を再計算しない。

## 60日基準制御

非作動区間を省かず独立に計算する方法：

```powershell
python scripts/run_nrlmsis_case.py --weather docs/nrlmsis2_necessity/data/drivers.csv --case baseline --duration-s 5184000 --cache-seconds 120 --prefix outputs/nrl_convergence_control_new/baseline_120s --output outputs/nrl_control_60d_new
```

計算を節約する方法：

```powershell
python scripts/reuse_nrlmsis_inactive_prefix.py --all-low outputs/nrl_low_60d_new --output outputs/nrl_control_inactive_new
```

`inactivity_proof.json` が60日全区間の非作動を確認した場合、その出力が完了した基準制御ケースである。途中でHIGH条件に達した場合は、出力されたprefixから実際の基準制御を60日まで継続する。

```powershell
python scripts/run_nrlmsis_case.py --weather docs/nrlmsis2_necessity/data/drivers.csv --case baseline --duration-s 5184000 --cache-seconds 120 --prefix outputs/nrl_control_inactive_new --output outputs/nrl_control_60d_new
```

## 後処理・密度監査

```powershell
python scripts/audit_nrlmsis_cache.py --run outputs/nrl_low_60d_new --weather docs/nrlmsis2_necessity/data/drivers.csv --cache-seconds 120 --output outputs/nrl_low_60d_new/cache_audit.json
python scripts/audit_nrlmsis_cache.py --run outputs/nrl_control_60d_new --weather docs/nrlmsis2_necessity/data/drivers.csv --cache-seconds 120 --output outputs/nrl_control_60d_new/cache_audit.json
python scripts/analyze_nrlmsis_comparison.py --all-low outputs/nrl_low_60d_new --controlled outputs/nrl_control_60d_new --output outputs/nrl_comparison_new
python scripts/plot_nrlmsis_validation.py --run outputs/nrl_direct_low_new --output outputs/nrl_validation_plot_new
```

全区間再利用時は、上記のcontrolledの入力パスを `outputs/nrl_control_inactive_new` にする。大きなNPZを持たない環境でも、Gitのhourly_summary、gap_metrics、tolerance_exceedanceと図で今回の主要判断を読める。

## 人工外乱の回帰

```powershell
python scripts/check_baseline_regression.py --output outputs/nrl_legacy_regression_new
```

3ケースを再実行し、旧基準の集計・イベント・主要時系列を比較する。同じ再実行済みsuiteを再検査する場合は `--existing-suite` で重複計算を省ける。
