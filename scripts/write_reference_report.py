"""Render the two-layer Japanese report from verified saved metrics only."""
import json
from run_reference_dispersion import DOC,DATA


def read(name):return json.loads((DATA/name).read_text())
def f(x,d=2):return '—' if x is None else f'{x:,.{d}f}'
def day(x):return f(None if x is None else x/86400,4)
def table(headers,rows):
    escape=lambda x:str(x).replace('|','\\|')
    return '\n'.join(['| '+' | '.join(map(escape,headers))+' |','| '+' | '.join(['---']*len(headers))+' |']+
        ['| '+' | '.join(map(escape,row))+' |' for row in rows])


def main():
    verify=read('verification.json')
    if not verify['passed']:raise ValueError('Verified results required')
    cases=['A','B','AB'];on={d:{c:read(f'{d}d/{c}_metrics.json') for c in cases} for d in (30,60)}
    off={d:{c:read(f'{d}d/{c}_OFF_metrics.json') for c in cases} for d in (30,60)}
    diag={c:read(f'60d/{c}_diagnostics.json') for c in cases};s=read('60d/reference_summary.json')
    profiles=read('profiles.json');meta=read('profile_metadata.json');a,b,ab=[on[60][c] for c in cases]
    passed500=all(x['max_error_m']<=500 and x['first_label_order_reversal_s'] is None for x in (a,b,ab))
    summary_table=table(['Case','全LOW・60日','制御最大誤差 [m]','完了L-H-L [回]','追加平均高度損失 [m]'],[
        [c, f(off[60][c]['max_error_m']/1000)+' km' if off[60][c]['first_label_order_reversal_s'] is None
         else day(off[60][c]['first_label_order_reversal_s'])+'日で順序逆転',
         f(on[60][c]['max_error_m']),on[60][c]['completed_cycles'],f(on[60][c]['extra_mean_altitude_loss_m'])] for c in cases])
    summary=f'''# 研究判断用サマリー

## 1. 結論

Reference 0.5 m / 0.5%条件で、現行controllerは{'A/B/ABすべての衛星順序を保ち、60日間500 m以内を維持できた' if passed500 else '500 m以内の維持に限界が確認された'}。ただし、精度を保つための姿勢切替は60日末まで続き、制御が停止して落ち着いたとは言えない。無制御ではAもkm級へdriftし、B・ABは約1か月で番号順序が逆転した。したがって、この固定条件は今後の制御機構解析に使える一方、現行controllerの効率が十分とはまだ判断できない。

> {profiles['statement']}。

## 2. 主要表

{summary_table}

誤差は従来の固定基準半径での100 km間隔からのずれ。B・ABの逆転後を正常な精度として比較していない。500 mはReference評価閾値であり、実ミッション要求ではない。3ケースとも100 m閾値は超過したが、500/1000 m閾値の超過は{'なかった' if passed500 else '付録で区別した'}。

## 3. A-ref：初期残差のみ

無制御の最大誤差は30日{f(off[30]['A']['max_error_m']/1000)} km、60日{f(off[60]['A']['max_error_m']/1000)} km。制御は{f(a['max_error_m'])} mへ抑えたが、初期補正後も再切替した。「最初に直せば終わる」挙動は確認できず、残る速度差と二値パルスの影響を調べる必要がある。

## 4. B-ref：持続個体差のみ

無制御は{day(off[60]['B']['first_label_order_reversal_s'])}日で順序逆転した。制御は{f(b['max_error_m'])} mまで抑え、持続する抗力差へ継続的に反応した。原因自体が消えないため、切替が続くことを初期条件補正と同一視できない。

## 5. AB-ref：両方

同じA配置と同じB配置を組み合わせた結果、制御最大誤差は{f(ab['max_error_m'])} m。A/Bの最大誤差を足した値にはならない。各gapでの打消し・強め合いと切替時刻の変化を含む、二値制御の応答として読むべきである。

## 6. Low conditionとの比較

旧Aのsample標準偏差は0.1134 m、旧Bは0.2304%。今回の実配列倍率は約4.41倍・2.17倍で、名目値比とは異なる。旧ABは配置も異なるため、純粋な分散倍率比較には使えない。制御下の誤差・cycle数・高度コストは単純比例せず、旧結果より誤差が小さくなる場合もあった。

## 7. 次の判断

1. このprofile・大気・controllerを固定し、最小HIGHパルスが残す半長軸差・速度差を先に調べる。
2. その結果を踏まえ、dwellの役割を別実験で評価する。今回のrunでは調整していない。
3. ゲイン探索やcontroller再設計は、パルス由来の負荷を理解した後に判断する。
'''
    parts=[summary,'# 詳細技術付録\n']
    parts.append('## A. Reference profile definition\n\n保存済みA seed 4、B生成seed 10002を再使用した。新しい乱数探索はない。B旧ケース名の末尾「2」は生成seed 10002を表す。ABの半長軸offsetはAと、固定抗力倍率はBと完全一致する。21機の数値は[profiles.json](data/profiles.json)および[reference_profiles.csv](data/reference_profiles.csv)へ固定した。\n')
    parts.append('## B. realized standard deviation\n\n今回の「0.5」は最大値ではなく21要素のsample標準偏差である。平均を引いてから規格化し、ddof=1を使用した。\n\n'+
        table(['Case','保存offsetのsample std [m]','初期aから再算出 [m]','εのsample std','初期最大|e| [m]'],[
            [x['case'],f(x['specified_a_sample_std_m'],9),f(x['initialized_a_sample_std_m'],12),f(x['B_sample_std'],9),f(x['initial_max_error_m'],12)] for x in meta['reference']])+
        '\n\n半長軸へ加算した後の約5×10^-11 mの標準偏差差は浮動小数点の丸め。位置は変えていない。A/AB初期最大gap速度は約3.118 mm/sで、円軌道一次近似との差は約2.0×10^-10 m/s。\n\n`δa = 0.5 (z-mean(z))/std(z,ddof=1)`、`ε = 0.005 (z-mean(z))/std(z,ddof=1)`。これは元の配置を保ちながら、平均0と指定sample標準偏差をそろえる操作である。\n')
    parts.append('## C. all LOW results\n\nAは60日でも順序を保つがdriftが蓄積する。B/ABは順序逆転するため、その後の固定番号誤差を正常な編隊精度と解釈しない。\n\n'+
        table(['Case','日数','逆転前max [m]','全時間の固定番号max（診断）[m]','初回逆転 [日]','最小signed gap [m]'],[
            [c,d,f(off[d][c]['max_error_before_order_reversal_m']),f(off[d][c]['max_error_m']),day(off[d][c]['first_label_order_reversal_s']),f(off[d][c]['min_signed_reference_gap_m'])] for d in (30,60) for c in cases])+'\n\n逆転の判定分解能は10 s。衝突・衛星間相互作用を含まない平均軌道モデルなので、逆転後OFFは反実仮想の診断軌道としてのみ保存する。\n')
    parts.append('## D. control results\n\n制御ONの精度は全観測時間で評価できる。最大誤差が小さくても、終了時速度差や後半の切替が残ることは収束とは異なる。\n\n'+
        table(['Case','日数','最大 [m]','全体RMS [m]','終了時max [m]','最大gap速度 [mm/s]','最終完了 [日]'],[
            [c,d,f(on[d][c]['max_error_m']),f(on[d][c]['rms_error_m']),f(on[d][c]['final_max_error_m']),f(on[d][c]['max_gap_rate_m_s']*1000,3),day(on[d][c]['last_switch_time_s'])] for d in (30,60) for c in cases])+
        '\n\nRMS = sqrt(mean(e_i(t)^2))：全保存時刻・全20gapを同じ重みで平均する。gap別RMS、最小全衛星対距離、現在半径の弧長誤差は各metrics.jsonに保存した。主要誤差はa0で換算し、現在半径上の物理弧長とは混同しない。\n')
    parts.append('### 同一の順序逆転前時間窓での改善\n\nB/ABの全60日OFF誤差による削減率を主指標にしない。OFFが逆転する直前の保存サンプルまで、ONも同じ時間窓へ切り詰めた。\n\n'+
        table(['Case','共通窓の終端 [日]','OFF max [m]','ON max [m]','max削減率 [%]','OFF RMS [m]','ON RMS [m]'],[
            [c,day(on[60][c]['matched_pre_allLOW_reversal_improvement']['last_valid_sample_s']),
             f(on[60][c]['matched_pre_allLOW_reversal_improvement']['all_low_m']),f(on[60][c]['matched_pre_allLOW_reversal_improvement']['controlled_m']),
             f(on[60][c]['matched_pre_allLOW_reversal_improvement']['reduction_percent']),
             f(on[60][c]['matched_pre_allLOW_reversal_RMS_improvement']['all_low_m']),f(on[60][c]['matched_pre_allLOW_reversal_RMS_improvement']['controlled_m'])] for c in cases])+'\n')
    threshold_rows=[]
    for c in cases:
        for label,m in [('ON',on[60][c]),('OFF',off[60][c])]:
            for x in m['exceedances']:
                threshold_rows.append([c,label,int(x['tolerance_m']),day(x['first_exceedance_time_s']),f(x['exceedance_fraction']*100),f(x['longest_continuous_exceedance_s']/86400,4)])
    parts.append('## E. threshold statistics\n\n100/500/1000 mはstrict/reference/loose評価閾値。実ミッション要求ではない。ONの500 m維持は保存時刻上の結果であり、保存間の連続時間保証ではない。\n\n'+
        table(['Case','状態','閾値 [m]','初回超過 [日]','超過時間率 [%]','最長連続超過 [日]'],threshold_rows)+
        '\n\nOFFの全60日時間率には逆転後の診断区間を含む。通常の精度として使える逆転前のON/OFF閾値統計は各metricsの`matched_pre_allLOW_reversal_thresholds`へ別保存した。期間長に依存する時間率を混同しない。\n')
    parts.append('## F. switching statistics\n\n精度維持の裏側で、繰返し姿勢変更を要している。回数はslew開始、完了L-H-LはLOW帰還完了を数える。\n\n'+
        table(['Case','L→H','H→L','全開始','全完了','完了cycle','stable HIGH [衛星h]','HIGH相当 [衛星h]','同時HIGH最大 [機]'],[
            [c,on[60][c]['l_to_h'],on[60][c]['h_to_l'],on[60][c]['transitions'],on[60][c]['completions'],on[60][c]['completed_cycles'],
             f(on[60][c]['high_stable_total_s']/3600),f(on[60][c]['high_equivalent_total_s']/3600),on[60][c]['max_simultaneous_stable_HIGH']] for c in cases])+
        '\n\nHIGH相当時間はintegral((A-A_LOW)/(A_HIGH-A_LOW)dt)を全衛星で合計し、slew中も含む。stable HIGHは姿勢変更完了後のHIGHだけ。衛星別回数・cycle数・最終姿勢・高度は各`satellites.csv`に保存した。command/開始/完了の全時刻は`events.csv`にある。\n')
    parts.append('### 0–7 / 7–30 / 30–60日の運用負荷\n\n後半の活動が停止するか、精度と負荷がどう変わるかを区間別に確認する。率の大小だけで定常性や増大を証明するものではない。\n\n'+
        table(['Case','区間 [日]','RMS [m]','遷移/日','cycle/日','HIGH相当 [衛星h/日]','追加高度損失 [m/日]'],[
            [x['case'],f"{x['start_day']}-{int(x['end_day'])}",f(x['error_rms_m']),f(x['transitions_per_day']),f(x['cycles_per_day']),
             f(x['high_equivalent_satellite_hours_per_day']),f(x['extra_mean_altitude_loss_per_day_m'],3)] for x in s['segments']])+ '\n')
    parts.append('## G. altitude cost\n\n追加損失は通常の大気抵抗で全衛星が下がった量とは別である。同一profileの全LOWとONを同じ時刻で比較する。\n\n'+
        table(['Case','追加平均損失 [m]','全LOW自然平均損失 [m]','ON最終平均高度 [km]','ON最終高度幅 [m]','OFF最終高度幅 [m]'],[
            [c,f(on[60][c]['extra_mean_altitude_loss_m']),f(off[60][c]['mean_altitude_loss_m']),f(on[60][c]['final_mean_altitude_m']/1000,3),
             f(on[60][c]['final_altitude_span_m']),f(on[60][c]['final_allLOW_altitude_span_m'])] for c in cases])+
        '\n\nΔh_extra=mean(a_OFF-a_ON)。正なら制御が余計に軌道半径を下げた。OFF逆転後も高度は番号ごとの反実仮想軌道として計算できるが、安全な運用基準と見なさない。高度コストをそのまま寿命短縮日数へ換算していない。\n')
    parts.append('## H. A/B/AB comparison・q診断\n\nA/ABは初期速度差があるため時刻0でHIGH条件を満たす。BはLOW中に抗力差から半長軸差・速度差が育つ。この違いを、最初のtriggerとLOW帰還後の状態で区別する。\n\n'+
        table(['Case','初回HIGH [s]','代表衛星','初回qP [m/s²]','初回qD [m/s²]','初回局所LOW完了 [s]','その衛星の次HIGH [s]'],[
            [c,f(diag[c]['first_HIGH_trigger']['command_time_s'],0),diag[c]['first_HIGH_trigger']['satellite_id'],
             f"{diag[c]['first_HIGH_trigger']['command_q_P_m_s2']:.4e}",f"{diag[c]['first_HIGH_trigger']['command_q_D_m_s2']:.4e}",
             f(s['correction_diagnostics'][j]['first_local_LOW_return_s'],0),f(s['correction_diagnostics'][j]['first_local_next_HIGH_s'],0)] for j,c in enumerate(cases)])+
        '\n\n上表の「次HIGH」は最初にLOWへ戻った衛星自身の次回で、列全体の次回ではない。初回に同時起動した他衛星もイベントCSVに保持した。Aでも速度差が完全には消えず再切替する。Bには固定倍率差による再成長があり、一定D予測をそのまま適用できない区間がある。\n\n'+
        table(['Case','P優勢起動','D優勢起動','混合','dwell300 s帰還','300 s後の帰還','840 s再起動間隔'],[
            [c,*[diag[c]['command_contribution_counts'][k] for k in ('P_dominant','D_dominant','mixed')],
             diag[c]['returns_dwell_limited'],diag[c]['returns_after_dwell'],diag[c]['cycles_at_840s']] for c in cases])+
        '\n\nP/D優勢は絶対値2:1以上の寄与分類で、符号と打消しも保存する。q=Kp·C e+Kd·C vで、Cは開いた列の隣接gap信号を前側−後側へ写す行列。端では片側だけを使う。ここでCはballistic coefficient Bとは別物。\n\n120 sの往路slew＋300 sのHIGH dwell＋120 sの復路slew＋300 sのLOW dwellから、同一衛星の再起動下限840 sが生じる。最小パルスの詳細解析・dwell調整は今回実施しなかった。初回および観測後半の代表q、面積、姿勢、dwell解除時刻は`representative_triggers.csv`と小さいq時系列CSVに保存した。\n')
    parts.append('### 静穏・残留状態・予測\n\n一時的に全LOWへ戻ることと、制御が停止することは異なる。\n\n'+
        table(['Case','最長観測全LOW区間 [h]','局所cycle後・全gap最大速度の中央値 [mm/s]','関連gap最大速度中央値 [mm/s]','最後の完了 [日]','最終静穏candidate'],[
            [c,f(s['correction_diagnostics'][j]['longest_observed_allLOW_interval_s']/3600 if s['correction_diagnostics'][j]['longest_observed_allLOW_interval_s'] is not None else None),
             f(diag[c]['completed_local_Vres_median_m_s']*1000,3),f(diag[c]['completed_local_incident_Vres_median_m_s']*1000,3),
             day(on[60][c]['last_switch_time_s']),str(on[60][c]['quiet_candidate']['candidate_start_s'])] for j,c in enumerate(cases)])+
        '\n\n局所LOW帰還時に他衛星が制御中の場合も含むため、全gap最大速度の中央値だけで補正性能を決めない。符号付きe/v/P/D配列はローカルのend_states NPZへ保存した。最終静穏candidateは既存の50 m・1 mm/s・1800 sを満たす全LOWの最終連続区間であり、有限観測の定義である。\n\n一定D予測はq(t)≈qP(t0)+qD(t0)+(Kp/Kd)qD(t0)(t−t0)で、dwellと10 s指令格子を加える。残留速度でP項が増えるという意味だが、B差や大気変動がDを変えると外挿が外れる。許容D変化10%などの有効性フラグは予測精度の保証ではない。各diagnosticsに実測との差・無効理由を保持した。\n')
    parts.append('## I. low/reference comparison\n\nA/Bは同じ空間配置を規格化した比較、旧ABは別配置を含む参考比較である。制御下で分散を大きくしても誤差や負荷が同じ倍率になる保証はない。\n\n'+
        table(['Case','旧条件（名目）','Reference','旧max [m]','Ref max [m]','旧cycle','Ref cycle','旧追加損失 [m]','Ref追加損失 [m]'],[
            [x['case'],{'A':'0.1 m','B':'0.3%','AB':'0.1 m + 0.3%'}[x['case']],{'A':'0.5 m','B':'0.5%','AB':'0.5 m + 0.5%'}[x['case']],
             f(x['old_control_max_error_m']),f(x['reference_control_max_error_m']),x['old_cycles'],x['reference_cycles'],
             f(x['old_extra_altitude_loss_m']),f(x['reference_extra_altitude_loss_m'])] for x in s['low_reference']])+
        '\n\n旧A/Bの実現sample標準偏差は0.11342855 m / 0.23042064%。旧ABは0.08509592 m / 0.28987992%で、配置も違う。新しいOFFの1–7日は、NRLMSIS-onlyの誤差を差し引くと実配列倍率にほぼ比例した。7日RMS比例差はA約0.0083%、B約0.0020%。[比例性確認](data/early_OFF_proportionality.json)は保存済みrunだけで計算した。順序逆転後や長期までの完全比例を主張しない。\n')
    parts.append('## J. verification / regression\n\n解析の追加によって既存baselineを変えていないことを確認した。\n\n'+
        f"- 全テスト **{verify['tests_passed']}件合格**。新規10件は規格化、共有profile、不正入力、初期間隔、共通時間窓、cycle区間計数を検証。\n"+
        '- 分散ゼロの新しい実行経路で24時間3ケースを新規照合。最大誤差・RMS・全イベント・高度損失・左右切替波速度・主要時系列がsoftware referenceと一致。\n'+
        '- 30→60日境界のa/lambda/area/modeは完全一致。30日までのイベント列も完全一致。密度キャッシュ再起動の差は既存の厳しい照合範囲に収まった。\n'+
        '- ON/OFFの初期状態・profile・大気・宇宙天気・積分条件を機械比較し、物理上の違いが制御ON/OFFだけであることを確認。\n'+
        '- 新軌道上で直接NRLMSISと密度cacheをspot照合。30/60日とも既存の相対差0.1%基準に合格。\n'+
        '- 旧追跡ファイルのSHA256、旧ローカルrunのサイズ・mtimeが不変。大きな旧runの全内容を再ハッシュした検証ではない。\n\n[検証JSON](data/verification.json)、[回帰JSON](data/regression_24h.json)、[非破壊JSON](data/preservation.json)。24時間referenceはsoftware regression用で、制御の科学的妥当性を証明するものではない。\n')
    parts.append('## K. reproduction・図\n\n[再現手順](REPRODUCTION.md)に実行順序、評価定義、checkpoint継続、後処理のみの再利用方法を記載した。新規物理計算は6ケースの0–30日と同じcheckpointからの30–60日だけであり、旧run・旧profileを再計算していない。\n\n'+
        '\n'.join(f'![図{j}](figures/{name})' for j,name in enumerate(['01_error_comparison.png','02_A_gaps.png','03_B_gaps.png','04_AB_gaps.png','05_switching_activity.png','06_extra_altitude_loss.png'],1))+
        '\n\n図1はOFFを逆転前まで表示。図2–4は全20gapと代表gap速度。図5は累積遷移と同時HIGH数。図6は同一profileの全LOW比追加平均高度損失。描画用は60 s間引きにイベント時刻・最大誤差時刻を追加し、metricsは10 s全保存データから計算した。\n')
    parts.append('## L. Git・変更内容\n\n'+f"作業branch：`{verify['branch']}`。起点commit：`{verify['base_commit']}`。mainへmergeしない。今回のcommit hashは最終応答と`git log -1`で確認できる。\n\n"+
        '- `src/sat_string/reference_condition.py`：保存profileの規格化と凍結後の検査。\n'+
        '- `src/sat_string/reference_metrics.py`：既存評価への共通窓RMS・閾値統計・区間cycle率追加。\n'+
        '- `configs/reference_dispersion.yaml`と`run_reference_dispersion.py`：固定入力とcheckpoint実行。\n'+
        '- `analyze_reference_dispersion.py`、`audit_reference_density.py`：既存後処理・密度監査を再利用する薄い接続。\n'+
        '- `supplement_reference_results.py`、`reference_profile_metadata.py`、`check_reference_proportionality.py`：共有profile、旧条件比較、代表q、初期比例性。\n'+
        '- `check_reference_regression.py`、`preserve_reference_inputs.py`、`verify_reference_results.py`：24h分散ゼロ回帰、旧成果非破壊、境界・入力・監査確認。\n'+
        '- `plot_reference_dispersion.py`、`write_reference_report.py`、`tests/test_reference_condition.py`：6図、報告、10テスト。\n'+
        '- docsにはprofile、config、集計JSON/CSV、小さいcycle sample、イベント、6図、再現手順を保存。巨大な配列・全幅cycle表・全HIGH trigger表はローカル保存。\n\n'+
        '### 観測・解釈・未確定\n\n**観測**：表の精度・順序・切替・高度損失は固定Referenceの数値結果。**解釈**：現行制御はこの入力に対し精度を保つが、繰返し切替と追加抗力を要する。**未確定**：実機分散の統計的代表性、500 mの実ミッション要求、永久安定性、許容できる姿勢負荷・寿命損失、他のprofileや宇宙天気への一般性。今回の結果だけで解決したとはしない。\n')
    path=DOC/'REPORT.md'
    if path.exists():raise FileExistsError(path)
    path.write_text('\n'.join(parts),encoding='utf-8')
    print('Saved',path,'summary characters',len(summary))


if __name__=='__main__':main()
