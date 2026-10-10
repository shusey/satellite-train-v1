# 観測に基づく外乱解析の再現

## 既存データを読む場合

新しい出力先を使う。既存フォルダが存在すると処理を拒否する。

```powershell
python -m pip install -e ".[test,observations]"
python scripts/compare_grace_density.py --raw-dir outputs/observed_density_20261010/raw --drivers docs/observed_density_design/data/drivers.csv --output outputs/observed_reproduce/grace_ab
python scripts/summarize_density_design.py --analysis outputs/observed_density_20261010/analysis_v2 --grace-ab outputs/observed_reproduce/grace_ab --output outputs/observed_reproduce/assessment
```

この方法はCHAMP/GRACEの保存NPZを再利用し、主データのNRLMSIS密度を再計算しない。GRACE A/B比較のみ小さい新規背景計算を行う。A/B保存比較もあるなら、後者の `--grace-ab` に既存 `outputs/observed_density_20261010/grace_ab` を指定する。

## データ取得から始める場合

PANGAEAはZIP全体を取得しない。選択月だけをHTTP rangeで読む。

```powershell
python scripts/download_density_zip_members.py --url https://download.pangaea.de/dataset/931347/files/champ_TND-IGG_RL01.zip --cache outputs/observed_reproduce/raw --member champ_TND-IGG_RL01_2003-10.nc --member champ_TND-IGG_RL01_2003-11.nc --receipt outputs/observed_reproduce/champ_download.json
python scripts/download_density_zip_members.py --url https://download.pangaea.de/dataset/931347/files/grace_TND-IGG_RL01.zip --cache outputs/observed_reproduce/raw --member grace_TND-IGG_RL01_2003-10.nc --member grace_TND-IGG_RL01_2003-11.nc --receipt outputs/observed_reproduce/grace_download.json
python scripts/download_observed_density.py --manifest docs/observed_density_design/data/aux_download_manifest.json --cache outputs/observed_reproduce/raw --receipt outputs/observed_reproduce/aux_download.json
python scripts/analyze_observed_density.py --raw-dir outputs/observed_reproduce/raw --weather-raw docs/observed_density_design/data/observed_support.csv --output outputs/observed_reproduce/analysis
python scripts/compare_grace_density.py --raw-dir outputs/observed_reproduce/raw --drivers outputs/observed_reproduce/analysis/drivers.csv --output outputs/observed_reproduce/grace_ab
python scripts/summarize_density_design.py --analysis outputs/observed_reproduce/analysis --grace-ab outputs/observed_reproduce/grace_ab --output outputs/observed_reproduce/assessment
```

`data/download_receipts.json` に凍結したSHA256を、PANGAEAの新取得receiptのselected項目と照合する。異なれば同一製品版ではないので、そのまま同じ基準解析と扱わない。補助download manifestはSHA256検証を自動実施する。全期間ZIP・大型NC・年tar・NPZはGit外に置く。

F10.7/apは `observed_support.csv` を入力にするため、宇宙天気を再取得する必要はない。元のSW-All.csvとの出典・checksum関係は `data/provenance.json` に保存。GPS→UTCは1999–2005年の13秒差のみを明示対応しており、年代変更時は実装を見直す必要がある。今回のスクリプトは選択した2003年期間用で、汎用の全球解析器ではない。

日本語図にはWindowsのMeiryoを使った。別環境では日本語フォントを設定する。japanized-matplotlibは追加していない。フォントcacheへの書込みが許されない環境では、新しい書込み可能なローカルフォルダを `MPLCONFIGDIR` に指定する。

## テスト・回帰

```powershell
python -m pytest tests/test_observed_density.py -q
python -m pytest -q
python scripts/check_baseline_regression.py --existing-suite outputs/nrlmsis2_20261010/regression_24h/suite --output outputs/observed_reproduce/regression
```

今回は物理coreを変更していないので、保存済みの24時間再実行suiteを照合に再利用した。suiteが別PCに存在しない場合のみ、最後のコマンドから `--existing-suite` を除き、新しい3ケースを実行する。

## 保存物

- `data/window_statistics.csv`、`scale_statistics.csv`：全6組・全5距離スケールの結果。
- `data/rho400_pair_relative_difference.csv`：規格化物理密度の2点差。
- `data/local_events_with_weather.csv`：正負全候補窓、fit品質、地方時／緯度／高度、Kp/ap/Dst/AE。
- `data/representative_segments.csv`：3候補の小さい実データ区間。
- `data/design_candidates.json`：A・幅・背景と、寿命／運動未確定の明示。
- `data/grace_ab_statistics.csv`、`grace_ab_sample.csv`：粗いA/B比較。
- `data/quality_control.csv`、`drivers.csv`、`space_weather_hourly.csv`、`observed_support.csv`：品質と実入力。
- `figures/`：日本語4図。
- `data/regression.json`、`preservation_check.json`、`completion.json`：検証情報。

生データと全処理NPZの保存先は `outputs/observed_density_20261010`。最初の入力照合失敗フォルダ `analysis` は上書きせず残し、正規解析は `analysis_v2`。Gaussian偽陽性対策後の最終集計は `assessment`。ここでのfitは密度再計算を伴わない。

## Gitの確認

```powershell
git branch --show-current
git log -2 --oneline
git status --short
```

作業branchは `astra/observed-density-design`。mainへmergeせず、通常のpushだけを行う。
