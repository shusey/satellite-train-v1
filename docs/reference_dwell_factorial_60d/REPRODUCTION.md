# 2×2 dwell実験：30日checkpointから60日への延長

## 再利用するもの

30日研究: `docs/reference_dwell_factorial_30d/`、基点commit `e735788be82fff9a26454c342e9d8845c0c78b4a`。
初期半径分散sample std(ddof=1)=0.5 m、B差なしの固定A profileを使用する。再抽選しない。
Reference条件は実機母集団の観測標準偏差ではなく、工学的な感度条件である。

- H300-L300: `outputs/reference_dispersion_20261012/batch_30d` と `batch_60d` のmember 0。
- 全LOW高度基準: 同じ2 batchのmember 3。
- H150-L150: `outputs/reference_a_dwell_20261017/D150_7d`、`outputs/reference_a_dwell_60d_20261017/D150_30d`、`D150_60d`。
- 交差2条件のprefix: `outputs/reference_dwell_factorial_20261018/H150-L300`、`H300-L150`。

大きな生データはローカル保存で、Gitへは追加しない。新規物理計算は交差2条件の30–60日suffixのみ。

## 延長

以下はsuffixが未計算の場合だけ実行する。同じ出力フォルダがあればrunnerは停止し、上書きしない。

```powershell
.analysis-venv/Scripts/python.exe scripts/extend_reference_dwell_factorial.py --case H150-L300
.analysis-venv/Scripts/python.exe scripts/extend_reference_dwell_factorial.py --case H300-L150
```

30日終了のcheckpoint、正規化設定、最終保存軌道、姿勢、固定profile、環境、weather SHA256、pymsis versionを照合する。
設定差分は観測終了30→60日だけ。姿勢の完了後dwell時刻を復元し、30日境界の指令を重複実行しない。
積分器・制御・密度評価coreは変更していない。

## 後処理

```powershell
.analysis-venv/Scripts/python.exe scripts/analyze_reference_dwell_factorial.py --days 60 --case H300-L300
.analysis-venv/Scripts/python.exe scripts/analyze_reference_dwell_factorial.py --days 60 --case H150-L300
.analysis-venv/Scripts/python.exe scripts/analyze_reference_dwell_factorial.py --days 60 --case H300-L150
.analysis-venv/Scripts/python.exe scripts/analyze_reference_dwell_factorial.py --days 60 --case H150-L150
.analysis-venv/Scripts/python.exe scripts/analyze_reference_dwell_factorial.py --days 60 --summarize
.analysis-venv/Scripts/python.exe scripts/plot_reference_dwell_factorial.py --days 60
```

接続境界の半径・位相・面積・modeは完全一致、密度は既存キャッシュの再始動による微差をrtol=1e-8以内で照合する。
境界の重複sampleを1個除き、0–60日を10秒格子で評価。図だけ60秒に間引く。
全累積イベントのprefix完全一致と、保存済み30日最大誤差・RMSの再一致も確認する。
全LOW高度基準も同じ保存runを結合し、違う大気・違うprofileを基準にしない。

区間別性能・負荷は0–7日、7–30日、30–60日。cycleはLOW完了時刻に割り当てる。
`late_30_60d` はLOW完了が30–60日に入るcycleの残留速度・パルス・再起動・近傍寄与。
再起動未観測は右打切り。全LOW時間率は30日境界をまたぐ区間も切り詰めて時間積算する。

## 検証

```powershell
.analysis-venv/Scripts/python.exe -m pytest -q
```

物理coreの変更がないため、24時間baselineは前回新規再生成したsuiteを現在の照合コードで再検証する。
`outputs/reference_dwell_factorial_20261018/regression_24h/suite` と `tests/fixtures/regression_24h` を
`sat_string.regression.compare_baseline` で照合する。suiteの再計算は不要。
HIGH/LOW個別指定の300/300・150/150冒頭6時間完全一致検証も前回の結果を再利用する。

公開用小規模成果物は、検証済みの `regression_24h/regression.json`、`diagonal_regression.json`、原資料manifestを
新規suffixと同じ出力rootへ保存した後、次で新しいフォルダへ書き出せる。

```powershell
.analysis-venv/Scripts/python.exe scripts/publish_reference_dwell_factorial.py --days 60 --output outputs/reference_dwell_factorial_60d_20261019/published_reproduction
```

保存済みのdocsや30日結果を上書きしない。四段階の全状態NPZはローカルanalysisフォルダに保存し、
Gitには設定・小さい集計・cycle/event表・図・文書を保存する。
