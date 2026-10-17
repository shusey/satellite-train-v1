# A-ref D150：60日観測への延長

ユーザーの60日観測指示に基づく。7日結果の再計算や追加のdwell条件は実行しない。

- 初期平均軌道半径sample std=0.5 m、B分散=0。固定済み21機A-ref profileを使用。
- dwell=150秒（HIGH/LOW両側）。Kp/Kd、閾値、slew120秒、制御周期10秒、環境・積分器・保存間隔は7日runから変更しない。
- NRLMSIS 2.0、epoch2024-03-01、historical active weather、120秒trajectory density cache、RK4最大2秒、保存10秒。
- 既存7日checkpoint→新規30日checkpoint→新規60日checkpoint。初期状態へ戻して同じ7日を再計算しない。
- 30日時点で完了フラグ、状態有限、軌道半径>地球半径、正のsigned gap・順序維持を確認してから60日へ進む。
- D300と同一profileの全LOWは既存30/60日runを再利用。D300の既存metrics・パルス診断も再利用する。
- 主比較は7/30/60日の最大誤差、RMS、終了時誤差、残留速度、cycles、HIGH換算量、同じ全LOW比の追加平均高度損失、再切替と閾値超過。
- 区間別0–7/7–30/30–60日のRMS、切替/日、cycles/日、HIGH量/日、高度損失/日を確認する。
- 7日・30日境界は同じ時刻の1点だけを残す。a/lambda/area/modeは完全一致、densityは既存cache再アンカーのrtol=1e-8を許容し、実差を保存する。
- 完全な有限時間観測であり、漸近安定・最適dwell・実機適合の証明には用いない。HIGH保持短縮とLOW待機短縮の総効果である。

既存7日成果は `docs/reference_a_dwell_150/`、既存D300 60日成果は `docs/reference_dispersion/`。今回の新規runは `outputs/reference_a_dwell_60d_20261017/`、公開用小規模成果は本ディレクトリへ保存する。
