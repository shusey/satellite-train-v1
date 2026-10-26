# 14日から30日への固定条件延長

ユーザーの指示によりN=21と301だけを延長する。無外乱nominalと中央nodeパルスdisturbedの両方で、同じ現行制御を用いる。

- 新規の物理積分はDay14〜30のみ。Day0〜14は保存済みrunを再利用。
- 元run：`outputs/local_disturbance_n_scaling_20261026/N021`、`N301`。
- 新規suffix：`outputs/local_disturbance_n_scaling_30d_20261026/N021`、`N301`。
- 制御：Kp=1/21600²、Kd=2/21600、h_on=1.5e-7、h_off=1.4e-7。
  HIGH/LOW dwell=300 s、slew=120 s、control=10 s、delay=0。
- 元のexponential高度依存密度、初期高度400 km、spacing100 km、mass4 kg、Cd2.2、LOW/HIGH面積.01/.03 m²を維持。
- Gaussian、NRLMSIS空間変動、B分散、初期半径分散、位相分散、ノイズ、新しい力は追加しない。
- 入力は元の一般化nodeパルスを維持。d_j=(d/dt)[a_ref*n(a_j)]_external。
  中央のみ、最大1e-6 m/s²、開始3時間、継続1時間、raised cosine。
  時刻はepochからの絶対経過秒。したがって延長部分ではd_j=0で、パルスを再生しない。
- RK4最大刻み2 s、保存10 sを維持。
- 主要指標はdelta e=e_disturbed-e_nominal。raw eも保持。
- 接続境界のa、lambda、density、area、modeを厳密一致で検証し、重複する14日目の保存点を1つ除く。
- 元14日のdelta eとイベント列prefixが完全一致することを後処理でも検証。
- 最後の3日（Day27〜30）を終盤窓とする。旧E_ultはDay11〜14で、窓幅は同一。
- 回復候補は全gapの|delta e|<=5 mが観測終了まで続き、確認期間が1日以上ある場合のみ。
- 14/21/30日の累積peakとDay14〜21、21〜27、27〜30の区間peak/RMS/負荷を保存する。
- Nは2点だけなので、べき乗則や一般的なN-scaling lawはfitしない。
- 今回の成功条件は時間方向のpeak更新・終盤残留・切替負荷を明らかにすること。
  有界性、漸近収束、string stabilityを有限30日だけで証明したとは表現しない。

変更は継続実行の検証・保存区間を結合する後処理・図・テストに限定する。物理・制御・積分coreは変更しない。
