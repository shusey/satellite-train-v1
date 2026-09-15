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
