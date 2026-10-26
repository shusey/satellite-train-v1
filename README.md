# sat-string

差動大気抵抗で有限開放列の 3U CubeSat を維持し、局所的な単発密度パルスに対する沿軌道間隔応答を評価する Python 3.11 以上向けシミュレータです。連続軌道状態には固定刻み RK4、姿勢・制御にはイベント駆動の離散状態機械を使います。

## セットアップ

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -e ".[test]"
```

実行時依存は NumPy、Matplotlib、PyYAML、imageio-ffmpeg です。MP4 は imageio-ffmpeg 付属の FFmpeg と H.264 (`libx264`) を使います。

## 実行

標準3ケース、全指標、300 dpi PNG、MP4をまとめて生成します。

```powershell
python scripts/run_standard_suite.py --config configs/baseline.yaml --output runs/standard_suite
```

単一ケースを実行する場合:

```powershell
python scripts/run_case.py --config configs/baseline.yaml --output runs/baseline
```

保存済み時系列だけから指標・図・アニメーションを再生成する場合（シミュレーションは再実行しません）:

```powershell
python scripts/make_report.py --run runs/standard_suite
```

既存出力先は保護されます。置換を明示する場合だけ `--overwrite`（または仕様互換の `-overwrite`）を付けてください。数値結果だけを先に確認したい場合、標準スイートに `--skip-animation` を指定できます。

## 標準スイート

- `controlled_disturbed`: 制御あり、密度パルスあり
- `controlled_undisturbed`: 制御あり、密度パルスなし
- `all_low_drag_baseline`: 全機低抗力固定、密度パルスあり

各ケースには `results.npz`、`metadata.json`、`normalized_config.yaml`、`switch_events.csv`、`metrics.json`、`metrics.csv`、`run.log` が保存されます。スイート直下には比較配列、集計指標、8枚の必須PNG、`formation_animation.mp4` が保存されます。回復判定の時間分解能は出力保存周期で、metadata と metrics に記録されます。

## 単位・符号・モデル仮定

- 設定と内部計算は SI 単位（m、s、kg、rad）です。図では時間を hour、距離を m または km に変換して軸へ明記します。
- 衛星番号は初期順序に固定し、平均経度は unwrap 済み連続角として保持します。並べ替えは行いません。
- 正の間隔誤差は前方が空きすぎ、負は前方へ詰まりすぎを表します。
- 軌道は平均化された近円軌道です。高度依存密度には局所指数近似を使います。
- 外乱は実在TID/TADの再現ではなく、基準軌道回転座標に置くテスト入力です。`propagation_speed_m_s` は正負とも設定できます。
- 制御は前後最近傍の相対情報だけを使う二値差動抗力です。高抗力化は片側入力であり、姿勢変更中の面積は cubic smoothstep で連続補間されます。
- 同じ制御時刻の全指令は、更新前の状態から一括決定されます。制御・保存・遅延開始・遷移完了のイベント時刻を積分が跨ぐことはありません。

## 評価上の注意と既知の限界

本実装は J2、離心率、面外運動、太陽・月摂動、大気共回転、風、通信遅延/途絶、観測誤差、実行中のグラフ変更、PID、実時刻NRLMSIS呼出しを含みません。衝突半径や恣意的な分断閾値も導入せず、最小隣接間隔と順序逆転を報告します。

ピーク比は数学的な string stability の証明ではなく「実験的伝播増幅率」です。伝播速度は片側3点以上かつ線形回帰の $R^2\ge0.8$ の場合だけ採用します。外乱なしでも高度依存密度による軌道減衰があるため、外乱影響は無外乱ケースとの差、高度コストは同じ外乱を受ける全機低抗力ケースとの差として解釈してください。

## テスト

```powershell
python -m pytest
```

単体試験に加え、非整数イベント時刻、3ケース独立実行、保存データ再生成、1.0 s/0.5 s刻み収束、移動外乱MP4生成を検証します。

## 切替波伝播の閾値走査

`h_on` の粗い走査（`1.0e-7`～`2.0e-7`, `1.0e-8`刻み）と遷移領域の追加走査
（`1.3e-7`～`1.5e-7`, `2.0e-9`刻み）を、重複を除いた19ケースとして実行します。
各ケースの `h_off` は `h_on - 1.0e-8` です。

```powershell
python scripts/run_switching_wave_sweep.py `
  --config configs/baseline.yaml `
  --output runs/switching_wave_sweep
```

出力先には、各閾値を1行にした `switching_wave_sweep.csv`、5枚の集約グラフ、
代表閾値（既定値は `1.0e-7`, `1.4e-7`, `2.0e-7`）の距離–切替開始時刻フィット図、
`sweep_metadata.json`、`run.log` が保存されます。CSVの `delta_t_left/right` は各隣接伝播時間、
`delta_t_pairs_left/right` は対応する内側・外側衛星IDを含むJSON配列です。

左右は衛星添字について `left: i < i_center`、`right: i > i_center` と定義します。
外乱中心衛星は外乱開始時の中心座標に最も近い衛星から自動決定し、必要なら
`--center-satellite`（1始まり）で明示できます。フィットは中心から4衛星以上離れた点のみを使い、
3点未満なら速度、R²、RMSEを `NaN` とします。標準偏差は母標準偏差（`ddof=0`）です。

## 制御価値・抗力コスト解析

現在の標準条件から、制御なし／制御あり／無外乱を比較し、7条件の閾値実験を行います。

```powershell
python scripts/run_control_cost_analysis.py --config configs/baseline.yaml --output outputs/control_cost_analysis
```

[解析結果・図・研究上の判定](docs/control_cost_analysis/report.md)、
[実装モデル](docs/control_cost_analysis/model.md)、
[PD制御の導出と保証範囲](docs/control_cost_analysis/pd_analysis.md)、
[検証と再現手順](docs/control_cost_analysis/reproduction.md)を参照してください。
24時間の見かけのコスト削減は、72時間では到達遅れと分かった候補があります。

## 評価基盤・24時間回帰・観測時間

[評価指標の定義](docs/evaluation_observation_horizon/definitions.md)、
[再現方法](docs/evaluation_observation_horizon/reproduction.md)、
[観測時間の結果と推奨](docs/evaluation_observation_horizon/report.md)を参照してください。
固定された24時間基準はsoftware regression referenceであり、暫定ゲインの科学的妥当性を意味しません。

```powershell
python scripts/check_baseline_regression.py --output outputs/regression_new
python scripts/evaluate_saved_case.py --run outputs/regression_new/suite/controlled_disturbed --all-low outputs/regression_new/suite/all_low_drag_baseline --output outputs/evaluation_new --epsilon-e 50 1 5 10
python scripts/run_observation_case.py --hours 168 --output outputs/observation_7d_new
```

新しいスクリプトは既存の出力先への書込みを拒否します。保存データだけで判定閾値・時間窓・観測終了時刻を変えて比較できます。

## 長時間再切替・qのP/D解析

[長時間結果と次段階への判断](docs/long_horizon_q_analysis/report.md)、
[波・q・高度損失の定義](docs/long_horizon_q_analysis/definitions.md)、
[保存runの継続と再現手順](docs/long_horizon_q_analysis/reproduction.md)を参照してください。
物理・制御・積分coreを変更せず、既存14日の保存状態から段階的に継続します。
波の時間的群分けと、衛星ごとのL-H-L cycleを区別します。

```powershell
python scripts/continue_observation_case.py --source outputs/evaluation_20261007/runs/central_14d --output outputs/long_horizon_new/central_21d --days 21
python scripts/continue_observation_case.py --source outputs/evaluation_20261007/runs/low_14d --output outputs/long_horizon_new/low_21d --days 21
python scripts/analyze_long_horizon.py --run outputs/long_horizon_new/central_21d --all-low outputs/long_horizon_new/low_21d --output outputs/long_horizon_analysis_new/21d
```

## dwell・cycle終了時の残留状態

[150/300/600秒の比較と因果仮説の評価](docs/dwell_residual_analysis/report.md)、
[残留速度・P/D予測の定義](docs/dwell_residual_analysis/definitions.md)、
[データ再利用・実行・検証手順](docs/dwell_residual_analysis/reproduction.md)を参照してください。
各衛星のL-H-L終了と、全衛星LOWになる境界を区別して評価します。
物理・制御・積分coreと凍結24時間基準は維持しています。

## Kd単独変更・残留速度・長時間再切替

[0.3/1/3倍の30・60日比較](docs/kd_residual_analysis/report.md)、
[速度分布・P/D分類・予測の定義](docs/kd_residual_analysis/definitions.md)、
[再利用と再現手順](docs/kd_residual_analysis/reproduction.md)を参照してください。
Kpとdwell等を固定してKdだけを指定できます。1倍baselineは保存runを再利用します。
残留速度最大値だけでなく、RMS・符号・空間分布・再切替コストを比較します。

## Kp単独変更・位置誤差の蓄積と再切替

[0.3/1/3倍の比較と平易な解説](docs/kp_residual_analysis/report.md)、
[位置誤差・起動・残留状態・予測の定義](docs/kp_residual_analysis/definitions.md)、
[データ再利用と再現手順](docs/kp_residual_analysis/reproduction.md)を参照してください。
Kd・dwell・閾値等を固定してKpだけを指定できます。1倍baselineは保存runを再利用します。
判定量q_Pと物理的な位置誤差を区別し、HIGH起動時にどれだけ誤差が蓄積していたかも評価します。

## Kp中間倍率の精度・コスト比較

Kp=2倍・2.5倍を追加し、既存1倍・3倍を再利用した30/60日比較。
普段読む[研究判断用サマリー](docs/kp_intermediate_analysis/report.md)と、
必要時に読む[詳細技術付録](docs/kp_intermediate_analysis/technical_appendix.md)を分離した。
数値条件・保存先・再現方法は[再現手順](docs/kp_intermediate_analysis/reproduction.md)を参照。

## NRLMSIS 2.0：経験大気モデル下での制御必要性

[研究判断用サマリー](docs/nrlmsis2_necessity/report.md)、
[詳細技術付録](docs/nrlmsis2_necessity/technical_appendix.md)、
[再現手順](docs/nrlmsis2_necessity/reproduction.md)を参照。
`density_model=nrlmsis2` は任意依存 `pymsis==0.12.0` と明示的なモデル2.0指定を使用する。
2024年1月1日〜3月1日UTCの歴史的F10.7/Apを固定し、全LOWと既存基準制御を比較する。
元の指数密度・Gaussian／half-cosine人工外乱モデルと24時間回帰基準は維持する。

```powershell
python -m pip install -e ".[test,nrlmsis]"
python scripts/run_nrlmsis_case.py --weather docs/nrlmsis2_necessity/data/drivers.csv --case all_low --cache-seconds 120 --output outputs/nrl_low_new
```

lambdaを地理経度として直接使わず、固定軌道面・地球自転・WGS84の測地座標へ変換する。
密度補間の検証と非作動区間の再利用方法は付録に記載する。NRLMSISは実測密度そのものではない。

## CHAMP・GRACE観測に基づく設計外乱の検討

[研究判断用サマリー](docs/observed_density_design/report.md)、
[データ仕様・定量結果・留保を含む詳細技術付録](docs/observed_density_design/technical_appendix.md)、
[選択データ取得と再現手順](docs/observed_density_design/reproduction.md)を参照。
NRLMSIS 2.0比の観測偏差、400 km規格化、100 kmの短時間差推定、GRACE A/B比較を解析する。
既存simulation coreと人工外乱は維持。候補の寿命・相対移動が未確定なので、標準外乱は置き換えていない。

## A-ref：HIGH側・LOW側dwellの2×2分離実験

[30日結果と研究判断](docs/reference_dwell_factorial_30d/REPORT.md)、[再現手順](docs/reference_dwell_factorial_30d/REPRODUCTION.md)。
HIGH側短縮は主に操作量・高度コストを減らし、LOW側短縮は精度と負荷の交換条件へ作用する。
`dwell_high_time_s` / `dwell_low_time_s` は独立指定でき、未指定の側は既存 `dwell_time_s` を使う。
対角2条件は保存runを再利用し、新規の交差2条件だけ30日計算した。60日へは延長していない。

## A-ref：HIGH/LOW dwell分離実験の60日検証

[60日結果と30日判断の再評価](docs/reference_dwell_factorial_60d/REPORT.md)、[checkpoint延長・再現手順](docs/reference_dwell_factorial_60d/REPRODUCTION.md)。

## 中央局所入力に対する衛星数スケーリング

既存300秒dwellの制御則を固定し、N=11/21/51/101/201/301でnominalと中央nodeパルスを各14日比較しました。
21機以上の最大追加gap errorは約88.92 mで横ばいでしたが、再切替と時間方向のpeak更新は残っています。
有限時間の探索結果であり、string stabilityの証明ではありません。
[研究判断と6図](docs/local_disturbance_n_scaling/REPORT.md)、[入力定義と実行前固定条件](docs/local_disturbance_n_scaling/DESIGN.md)、[再現手順](docs/local_disturbance_n_scaling/REPRODUCTION.md)。
交差2条件の30–60日だけを新規計算し、対角条件・全LOW・最初の30日は再利用した。
HIGH側短縮の低負荷効果は持続したが、LOW側短縮のRMS改善は後半には持続しなかった。
