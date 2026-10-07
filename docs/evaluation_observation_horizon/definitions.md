# 評価指標と24時間ソフトウェア回帰基準

## 対象と変更範囲

前回の解析成果 `a81426a` を引き継いだ。default branchはmain、物理・制御コードの元の研究版は `9449a78`。今回、運動方程式、制御則、姿勢状態機械、積分器、既存の設定・保存結果を変更していない。

代表設定は21機、100 km間隔、400 km高度、4 kg、CD=2.2、LOW/HIGH面積0.01/0.03 m²、密度3.584968e-12 kg/m³、スケール高度58,310 m。軌道長半径aと非折返し経度λをRK4で積分する円軌道近似。初期状態は全LOW、全機同じa。制御10 s、slew120 s、完了後dwell300 s、command delay0 s、積分最大2 s、保存10 s。外乱は10,800–14,400 s、密度増加50%、FWHM100 kmの回転座標系Gaussian。

間隔誤差は `e_i=a0(λ_(i+1)-λ_i)-d`、速度は `de_i/dt=a0(n_(i+1)-n_i)`。`p_i=Kp e_i+Kd de_i/dt`として、端部 `q_1=p_1, q_N=-p_(N-1)`、内部 `q_i=p_i-p_(i-1)`。LOWかつq>=h_onでHIGH、HIGHかつq<=h_offでLOWを指令する。slew中・dwell中・pending中は指令不可。

Kp=1/21600²、Kd=2/21600、h_on=1.5e-7、h_off=1.4e-7 m/s²は現行コードの暫定値。妥当性・最適性を保証しない。この関係を最終的な制御設計として正当化しない。

## 編隊誤差

`formation.max_abs_error_m = max_(k,i)|e_i(t_k)|`。

正式RMSは **全保存時刻・全隣接間隔に等しい重み**を与える：`sqrt(sum_(k,i) e_i(t_k)^2/(M(N-1)))`。各間隔RMSは `sqrt(sum_k e_i(t_k)^2/M)`。開始・終了サンプルも含む。非等間隔保存の場合にもこれはサンプルRMSであり、時間平均とは異なる。

既存control_costの時間積分RMSは変更せず、`formation.time_weighted_rms_error_m`にも併記する。24時間制御ケースではサンプルRMS=5.991741618 m、従来の台形積分RMS=5.991770632 m。差は定義によるもので物理結果の変更ではない。

`final_max_abs_error_m`は最後の保存時刻の最大絶対誤差であり、収束指標ではない。

## 切替・cycle

標準切替回数は、既存互換の **実際に開始された姿勢遷移数**。command、start、completionを別集計し、同じ順序の衛星別配列を保存する。衛星番号は1始まり、配列indexは番号-1。

イベントCSVはcommand時点で未来のstart/completionを記録する。従って時刻<=実際の観測終了時刻のイベントだけをそれぞれ数える。観測末尾の未開始指令・未完了slewは `pending_or_active_transitions_at_tmax` に保存する。

`last_switch_time_s`は全機の最後の **実際の遷移完了**。イベントなしはnull。未完了がある場合、この時刻だけで切替停止と解釈しない。

`min_switch_interval_s`は同一衛星の連続する **start間隔**の最小値。command/complete間隔も別保存。比較対象は定数slew+dwell、標準ケースの下限420 s。各衛星のイベントが1以下なら候補なし、そのケース全体でも候補がなければnull。

全LOW開始で、**完了した**L→Hとその後の完了したH→Lを1 L-H-L cycleとする。slew中に観測が終わったcycleは数えない。非交互・重複イベントはエラーとして検出する。

- secondary_switching_detected：いずれかの衛星で完了cycle>=2。
- multiple_switching_detected：いずれかの衛星で完了cycle>=3。
- secondary_response_started：未完了も含め、いずれかで2回目のL→Hが実際に開始。完了cycle数だけでは見落とす観測末尾の二次反応を補う。

HIGH exposureは既存のevent-exact積分を利用する。stable HIGH時間はslewを除外、equivalent HIGH時間は `∫(A-AL)/(AH-AL)dt` としてslewも含む。

## RecoveryとSettling

`EvaluationSettings(epsilon_e_m,epsilon_v_m_s,hold_window_s)`は後処理設定。省略時は既存config.recovery（50 m、0.001 m/s、1800 s）から取得し、物理条件の妥当性は主張しない。

許容条件は全間隔の `|e|<=epsilon_e` かつ `|de/dt|<=epsilon_v`。等号は既存recoveryとの互換性のため含める。標準評価開始は外乱終了14,400 s。無外乱ケースも比較のため同じ開始時刻を使う。

**Recovery**：許容条件をW以上維持した最初の区間について、`recovery_entry_time_s`は入口、`recovery_time_s`は確認完了時刻（絶対simulation時刻）。`recovery_time_after_disturbance_s`は確認時刻-外乱終了。後の逸脱・切替はこの判定を取り消さない。既存の `recovery` JSON構造・定義は保持。

**Settling**：外乱終了後かつ最後の実際のslew完了以後で、終了まで全保存時刻で許容条件を満たす最初の時刻。この最後の区間がW以上続き、観測末尾に未完了・pendingの姿勢遷移がない場合だけ成立する。定義上、HIGHのまま静止して許容条件を満たす場合も成立しうる。

`settled_by_tmax`はこの有限観測での成否。`verification_horizon_s=Tmax-settling_time_s`、未成立ならnull。10 sサンプルの間やTmax以後について保証しない。後に観測を延長すると成立から不成立へ変わりうる。入口時刻を使うためSettlingがRecoveryの確認時刻より早くても矛盾しない。

1/5/10/50 mなどの評価は `evaluate_saved_case.py --epsilon-e 50 1 5 10` で再積分なく取得できる。WはRecoveryのholdとSettlingの最低verificationの両方に適用する。

## 切替波

各衛星の初回 **L→Hのstart** を到達時刻とする。左右別に中心からの番号差×目標間隔を距離として `distance=v*t+b` の最小二乗を行う。既存互換で中心から3 hop以内をfitから除外。3点以上、外向きの時刻が単調増加、到達の抜けがなく、R²>=0.8、正の速度を満たすときだけ正式な速度を返す。

max_hopsは到達衛星の最大番号差。last_arrival_time_sは初回到達群の最遅時刻であり、後の二次切替の最終時刻とは異なる。左右の最遠連続到達衛星・最初の未到達衛星、fit状態も残す。未到達は有限観測の打切りであり、永久に伝播が停止した証明ではない。

外乱位置は番号をGaussianの回転座標位置に変換する：`s_center=(番号-(N+1)/2)*d`。21機なら11番0 m、6番-500,000 m、1番-1,000,000 m。これは外乱開始時の衛星への追従ではなく、既存の静止Gaussian中心を平行移動する。

## 高度損失

円軌道近似のため `h=a-R_E`、高度変化とa変化は同じ。

- mean_altitude_loss_m：`mean_i[a_i(0)-a_i(T)]`。正が損失。
- extra_loss_vs_all_low_m：**同じ初期状態・外乱・観測終了時刻**の全LOW基準に対し `mean_i[a_LOW,i(T)-a_CASE,i(T)]`。正が制御による追加損失。全時間平均ではなく、Tでの衛星算術平均差。

24時間の追加損失0.821149887 mは既存 `final_mean_altitude_loss_vs_all_low_m` / control_cost `extra_mean_a_loss_vs_low_m` と同じ計算。無外乱制御を外乱ありLOWと比較して制御コストと呼ぶことは避け、その値はnullにする。基準historyがない場合も追加損失はnull。

## 固定された回帰基準

`tests/fixtures/regression_24h/` は元の変更前24時間historyから作成。設定・全保存時刻の誤差/速度/a/高度/面積/姿勢状態・イベント・集計値・SHA256を保持する。今後configs/baseline.yamlのゲイン等が変更されてもfixtureは自動更新しない。

fixtureのテキスト改行は.gitattributesでLFに固定し、NPZはbinaryとして扱う。Windows/Unix間のcheckoutでも参照SHA256が変わらないようにした。

比較は3段階：集計、イベントの番号/向き/command/start/complete、主要時系列。許容誤差はmanifestに明示し、状態modeは厳密一致、時間1e-9 s、イベント1e-6 s、位置・高度1e-4 m、面積1e-12 m²、速度1e-9 m/s、rtol0。bit一致は要求しない。

これは **software regression reference**。新しい科学的に良い制御条件を拒否する基準ではなく、固定入力に対する実装の変化を検出するもの。通常pytestは保存fixtureの整合性とcheckerが改変を検出することを試験し、フル24時間再実行は専用コマンドで明示的に実施する。
