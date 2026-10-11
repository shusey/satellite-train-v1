"""Build a plain-Japanese decision report from saved metrics, not new physics."""
import json
from run_dispersion_control import DOC,DATA,OUT,configuration,frozen_profiles


def main():
    if (DOC/'README.md').exists():
        raise FileExistsError('A reviewed report already exists. Reuse it or select a new SAT_CONTROL_DOC; do not overwrite it.')
    days=60 if (DATA/'60d/comparison.json').exists() else 30
    rows={d:json.loads((DATA/f'{d}d/comparison.json').read_text())['cases'] for d in (30,60) if (DATA/f'{d}d/comparison.json').exists()}
    spec,cfg=configuration(days);profiles=frozen_profiles(spec)[:3]
    get=lambda d,c,on:next(r for r in rows[d] if r['case']==c and r['control_enabled']==on)
    text=['# 研究判断用サマリー','','<!-- Summary is finalized after scientific review; technical appendix below is generated from actual outputs. -->','','# 詳細技術付録','','## 1. 比較条件と乱数配置','',
        '今回変えたのは初期半長軸残差、衛星固有の抗力倍率、制御ON/OFFだけです。ゲイン・閾値・姿勢変更時間・大気場・積分器は固定しました。分散量は暫定スクリーニング値であり、実衛星の分布と確認された値ではありません。','',
        '|Case|指定した半長軸標準偏差|指定したB相対標準偏差|a乱数番号|B乱数番号|実現したa標準偏差|実現したB標準偏差|',
        '|---|---:|---:|---:|---:|---:|---:|']
    for p in profiles:
        ase=p['seed'] if p['case']!='B' else '未使用'; bse=p['B_seed'] if p['case']!='A' else '未使用'
        text.append(f"|{p['case']}|{.1 if p['case']!='B' else 0:g} m|{.3 if p['case']!='A' else 0:g}%|{ase}|{bse}|{p['realized_sigma_a_m']:.6f} m|{p['realized_sigma_B_fraction']*100:.6f}%|")
    text+=['','平均を引いたGaussian標本を再規格化していないため、実現した標準偏差は指定値と一致しません。Aは前回7日最大誤差の中央値に対応するseed4、Bはseed2を再利用しました。ABは前回の組合せseed0/10000を縮小した配置です。ABはA/Bと同じ個別配置ではなく、差の単純な加算や相乗効果を論じる比較ではありません。','',
        '各衛星の値は[data/satellite_profiles.csv](data/satellite_profiles.csv)、元配置への参照は[data/profiles.json](data/profiles.json)に保存しました。初期の平均経度は従来通り100 km間隔で、半長軸だけをずらしています。','',
        '```text\nKp = 1/21600² = 2.143347050754458e-9 s^-2\nKd = 2/21600 = 9.259259259259259e-5 s^-1\nh_on = 1.5e-7 m/s²; h_off = 1.4e-7 m/s²\n制御周期10 s / slew120 s / 完了後dwell300 s / command delay0 s\nN=21 / 目標間隔100 km / 初期高度400 km\n質量4 kg / Cd2.2 / LOW0.01 m² / HIGH0.03 m²\nNRLMSIS2.0 / pymsis0.12.0 / 質量密度kg/m³\n2024-03-01 UTC開始 / 傾斜角51.6° / RAAN0° / 初期位相0°\n既存120 s密度補間 / RK4最大2 s / 保存10 s\n```','',
        'B=C_D A/m（逆数の定義ではない）です。各衛星の倍率1+εをLOW/HIGH/姿勢変更中の面積すべてに掛けています。元のspacecraftのスカラー設定や既存coreは変更していません。NRLMSISは各ケースの実際のa・平均経度から既存幾何変換で評価し、旧密度時系列の単純再生はしていません。宇宙天気CSV、epoch、軌道面、キャッシュ方式は前回と同じです。','',
        '## 2. 精度と制御効果','',
        '無制御と制御は同じ配置・同じ観測時間で比較しています。RMSは全保存時刻と全20間隔の二乗誤差を等しい重みで平均し、平方根を取ったものです。終了時最大誤差は収束判定とは別です。']
    for d in rows:
        text+=['',f'### {d}日','', '|Case・制御|最大誤差 m|RMS m|終了時最大誤差 m|最大間隔速度 mm/s|最大誤差削減率|','|---|---:|---:|---:|---:|---:|']
        for c in ['A','B','AB']:
            for on in (False,True):
                r=get(d,c,on);red=(f"{r['error_improvement']['reduction_percent']:.3f}%"+('†' if not r['full_horizon_precision_reduction_valid'] else '')) if on else '—'
                text.append(f"|{c}・{'制御' if on else '全LOW'}|{r['max_error_m']:.3f}|{r['rms_error_m']:.3f}|{r['final_max_error_m']:.3f}|{r['max_gap_rate_m_s']*1000:.6f}|{red}|")
        text+=['','制御による削減率だけでなく、左側の実際の誤差を合わせて読む必要があります。†は全LOWに順序逆転があり、全期間の削減率が固定番号の参考診断であることを示します。逆転後の最大誤差は隣接編隊の精度として解釈しません。逆転前の同一時間窓による改善量も各metrics.jsonのmatched_pre_allLOW_reversal_improvementへ別保存しています。','',
            '|Case|姿勢遷移開始/完了|L→H / H→L|完了L-H-L|安定HIGH 衛星時間 h|HIGH換算 衛星時間 h|追加平均高度損失 m|','|---|---:|---:|---:|---:|---:|---:|']
        for c in ['A','B','AB']:
            r=get(d,c,True)
            text.append(f"|{c}|{r['transitions']} / {r['completions']}|{r['l_to_h']} / {r['h_to_l']}|{r['completed_cycles']}|{r['high_stable_total_s']/3600:.3f}|{r['high_equivalent_total_s']/3600:.3f}|{r['extra_mean_altitude_loss_m']:.3f}|")
        text+=['','安定HIGH時間は姿勢変更中を除きます。HIGH換算時間は∫(A−A_LOW)/(A_HIGH−A_LOW)dtで、姿勢変更中の部分面積も数えます。ここでの時間は21機の合計であり、経過時間とは異なります。姿勢遷移は実際に始まった変更を数え、L-H-LはHIGHとLOWへの変更が両方完了した1機単位の往復です。']
    text+=['','## 3. 仮想許容誤差と順序','',
        '100/500/1000 mは実ミッション要求ではありません。保存サンプルで絶対誤差が厳密に許容値を超えたかを判定し、時間割合は各10秒区間の左端値から求めました。未観測のサンプル間挙動を保証しません。','',
        '|日数|Case・制御|許容 m|初回超過 日|超過時間 %|最大連続超過 日|','|---:|---|---:|---:|---:|---:|']
    for d in rows:
        for c in ['A','B','AB']:
            for on in (False,True):
                for x in get(d,c,on)['exceedances']:
                    first='なし' if x['first_exceedance_time_s'] is None else f"{x['first_exceedance_time_s']/86400:.5f}"
                    text.append(f"|{d}|{c}・{'制御' if on else '全LOW'}|{x['tolerance_m']:g}|{first}|{x['exceedance_fraction']*100:.3f}|{x['longest_continuous_exceedance_s']/86400:.5f}|")
    text+=['','主要誤差は従来の固定基準半径a_ref上の誤差 e=a_ref(λ[j+1]−λ[j]−目標角間隔)です。現在半径上の物理弧長は隣接2機の半径平均×角差として別保存し、フィードバックには使いません。','',
        '|日数|Case・制御|最小固定基準の符号付きgap km|最小現在半径の符号付き弧長 km|現在半径上の最大誤差 m|全ペア最小直線距離 km|初回順序逆転 日|','|---:|---|---:|---:|---:|---:|---:|']
    for d in rows:
        for c in ['A','B','AB']:
            for on in (False,True):
                r=get(d,c,on);bad=r['first_label_order_reversal_s'];when='なし' if bad is None else f'{bad/86400:.5f}'
                text.append(f"|{d}|{c}・{'制御' if on else '全LOW'}|{r['min_signed_reference_gap_m']/1000:.6f}|{r['min_signed_current_radius_arc_m']/1000:.6f}|{r['max_current_radius_arc_error_m']:.3f}|{r['min_physical_pair_chord_m']/1000:.6f}|{when}|")
    text+=['','全ペア直線距離は同一軌道面の円軌道幾何から算出し、保存時刻に限定した距離です。衝突確率や運用安全性を保証する指標ではありません。順序逆転時刻以降は固定番号誤差を通常のformation precisionとして扱いません。','',
        '## 4. 初期過渡と長期運用','',
        '切替が初期だけで終わったかを見るため、0–7日、7–30日、30–60日を分けました。追加高度損失は同じ個体差の全LOWに対する差なので、自然な高度低下やB差自体による高度差とは区別できます。','',
        '|Case|区間 日|遷移/日|HIGH換算 衛星時間 h/日|区間RMS m|追加平均高度損失 m/日|','|---|---|---:|---:|---:|---:|']
    for c in ['A','B','AB']:
        for s in get(days,c,True)['segments']:
            text.append(f"|{c}|{s['start_day']}–{s['end_day']:g}|{s['transitions_per_day']:.3f}|{s['high_equivalent_satellite_s_per_day']/3600:.3f}|{s['error_rms_m']:.3f}|{s['extra_mean_altitude_loss_per_day_m']:.3f}|")
    text+=['','期間末の指令は最後の区間だけに含めます。RMSは各区間の両端保存サンプルを含みます。追加平均高度損失=mean(a_OFF(T)−a_ON(T))で、正の値ほど制御によって平均高度が余分に下がったことを表します。全LOW側の順序が崩れても平均半径差としては計算できますが、その時点の番号gap精度と同じ意味ではありません。','',
        '## 5. 切替原因・残留状態・予測の限界','',
        'qPは位置誤差、qDは速度差を切替判定へ写した量です。端部 q[1]=Kp e[1]+Kd v[1]、内部 q[j]=Kp(e[j]−e[j−1])+Kd(v[j]−v[j−1])、反対端 q[N]=−Kp e[N−1]−Kd v[N−1]です。q≥h_onでLOW→HIGH、q≤h_offでHIGH→LOWですが、姿勢変更・dwell中の指令は禁止されます。','',
        '|Case|初回HIGH 日・衛星|P優勢起動|D優勢起動|混合起動|300秒で帰還した往復|帰還が遅れた往復|最短同一機HIGH開始間隔 s|840秒に一致した間隔|','|---|---|---:|---:|---:|---:|---:|---:|---:|']
    for c in ['A','B','AB']:
        d=json.loads((DATA/f'{days}d/{c}_diagnostics.json').read_text());f=d['first_HIGH_trigger'];cc=d['command_contribution_counts']
        first='なし' if f is None else f"{f['command_time_s']/86400:.6f}・{f['satellite_id']}"
        text.append(f"|{c}|{first}|{cc['P_dominant']}|{cc['D_dominant']}|{cc['mixed']}|{d['returns_dwell_limited']}|{d['returns_after_dwell']}|{d['min_same_satellite_HIGH_start_interval_s']}|{d['cycles_at_840s']}|")
    text+=['','優勢は絶対値が相手の2倍以上という分類であり、P/Dが逆符号で打ち消す場合もあります。初回HIGHのqP/qD、関連gapの符号付き誤差・速度は*_HIGH_triggers.csvに残しています。HIGH完了から帰還指令まで300秒ならdwell支配、300秒より長ければ閾値条件による待ち時間も発生しています。最短840秒は120+300+120+300の2回の姿勢変更と両側dwellを足した値です。','',
        '残留速度V_res=max_gap|v|はLOW復帰完了直後の全gap速度です。他衛星の活動がある局所完了と、全LOWになる完了境界を分け、各状態の符号付きe/v/qP/qDをローカルの*_end_states_*.npzへ保存しています。要約CSVでは最大・RMS・絶対値中央値とe×vの悪化方向割合を残しました。','',
        '静穏中の簡易予測は q(t)≈qP(t0)+qD(t0)+(Kp/Kd)qD(t0)(t−t0)です。残った速度が一定なら、その速度に応じてP項が増え、閾値へ達する時刻を計算できるという意味です。実際にはNRLMSISの変動とB差が継続するため、qD一定が破れることがあります。後続切替、他機のHIGH、D変化10%以上などを無効理由として記録し、無効な予測を成功例に含めません。','',
        '|Case|局所完了時V_res中央値 mm/s|全LOW完了境界数|有効な全LOW次回予測数（実測・起点一致）|静穏候補開始 日|その後の監視 日|','|---|---:|---:|---:|---:|---:|']
    for c in ['A','B','AB']:
        r=get(days,c,True);d=json.loads((DATA/f'{days}d/{c}_diagnostics.json').read_text());q=r['quiet_candidate']
        start='なし' if q['candidate_start_s'] is None else f"{q['candidate_start_s']/86400:.6f}"
        horizon='—' if q['verification_horizon_s'] is None else f"{q['verification_horizon_s']/86400:.6f}"
        vr='—' if d['completed_local_Vres_median_m_s'] is None else f"{d['completed_local_Vres_median_m_s']*1000:.6f}"
        text.append(f"|{c}|{vr}|{d['quiet_packets']}|{d['quiet_valid_origin_matched_predictions']}|{start}|{horizon}|")
    text+=['','予測の一致を、閾値が将来横断される予測と、既に閾値を超えていてLOW側dwell解除を待つ予測に分けました。後者の時刻一致は待ち時間の再現であり、残留速度によるP項蓄積の証拠には数えません。以下はD変化10%以内・実測あり・予測した起点が実測と一致した区間に限定します。','',
        '|Case|将来横断の有効件数|時刻誤差の絶対値中央値 s|待ち時間に対する誤差中央値 %|既に閾値以上で待機した件数|','|---|---:|---:|---:|---:|']
    for c in ['A','B','AB']:
        d=json.loads((DATA/f'{days}d/{c}_diagnostics.json').read_text())['quiet_prediction_status_summary']
        f=d.get('future_upward_crossing',{});gate=d.get('already_at_or_above_on',{})
        text.append(f"|{c}|{f.get('count',0)}|{f.get('median_absolute_error_s',0):.1f}|{100*f.get('median_absolute_relative_error',0):.3f}|{gate.get('count',0)}|")
    text+=['','D項の大きさがP項の2倍以上でも、次の切替へ近づく変化をP項が担う場合があります。大小の分類は起動時の内訳であり、閾値へ近づいた原因の証明ではありません。D一定を10%以内とした有効判定も、予測時刻の小さな誤差を保証しません。BとABでは持続する個体差でDが変わり、単発外乱後と同じ予測精度を前提にできません。','',
        '静穏候補は既存評価条件50 m・1 mm/s・1800 sの全LOW最終区間です。数学的収束とは呼びません。この条件で候補がないことと、誤差を500 m以内に保てないことは同じではありません。*_mechanism_hourly.csvにはB差→半長軸変化率差→半長軸差→速度差→誤差の符号付き履歴と、同じ瞬間状態でB倍率だけを1へ戻した補助計算を保存しました。','',
        '## 6. 図','',
        '履歴図は1分抽出に切替時刻・極値を加えています。HIGH衛星数は全10秒保存の状態変化を省略せず表示し、初回q診断も10秒です。集計は全10秒保存データ・正確なイベントを使い、描画用間引きは集計へ適用しません。']
    for f,caption in [('01_precision_comparison.png','3ケースの全LOWと制御の最大誤差'),('02_A_gaps.png','Aの符号付き全gap誤差・速度'),('03_B_gaps.png','Bの符号付き全gap誤差・速度'),('04_control_activity.png','BのHIGH衛星数と各ケースの累積切替'),('05_extra_altitude_loss.png','同じ個体差の全LOWに対する追加平均高度損失'),('06_first_trigger_PD.png','A/B初回HIGHの判定量と関連gap速度')]:
        text+=['',f'### {caption}','',f'![{caption}](figures/{f})']
    text+=['','## 7. 観測・解釈・未確定','',
        '### 観測されたこと','', '上の表は指定した3つの配置とNRLMSISの1期間に対する数値結果です。誤差・イベント・高度コストは別々に保存し、順序逆転と予測無効も隠さず記録しました。','',
        '### そこから言えること','', 'Aでは初期半長軸差だけでなく、制御パルスが新しい半長軸差を作ることと、通常の大気空間差が継続することも考慮します。BではB倍率が制御後も残るため、LOWへ戻っても半長軸変化率の差が再発し得ます。再切替をすべて単発外乱の残留速度だけで説明することはできません。','',
        '### まだ言えないこと','', '0.1 m・0.3%が実衛星の典型的分散か、永久にboundedか、切替と高度コストが運用上許容できるかは未確定です。1配置・1大気期間の結果から一般的な性能保証や最適ゲインを主張しません。実機の姿勢変更制限・寿命予算がないため、負荷の「過大」判定も数値を示した上で留保します。','',
        '## 8. 再現方法','',
        '本研究で実際に実行した手順は下記です。A/B全LOWの前半30日は既存保存データを再利用しました。巨大NPYはGitに含めません。既存結果がある場合は再利用し、必要な再現計算は新しい保存先を指定してください。既存runを消して再実行する手順ではありません。','',
        '```powershell\n.analysis-venv/Scripts/python.exe -X utf8 scripts/run_dispersion_control.py --days 30\n.analysis-venv/Scripts/python.exe -X utf8 scripts/analyze_dispersion_control.py --days 30\n# data/extension_decision.jsonの30日結果をレビューしてからgoを記録\n.analysis-venv/Scripts/python.exe -X utf8 scripts/run_dispersion_control.py --days 60 --prefix outputs/dispersion_control_20261011/batch_30d\n.analysis-venv/Scripts/python.exe -X utf8 scripts/extend_dispersion_control_OFF.py\n.analysis-venv/Scripts/python.exe -X utf8 scripts/analyze_dispersion_control.py --days 60\n.analysis-venv/Scripts/python.exe -X utf8 scripts/audit_dispersion_control_density.py --days 60\n.analysis-venv/Scripts/python.exe -X utf8 scripts/plot_dispersion_control.py --days 60\n.analysis-venv/Scripts/python.exe -X utf8 scripts/report_dispersion_control.py\n.analysis-venv/Scripts/python.exe -X utf8 -m pytest -q\n```','',
        '30日後は姿勢状態・dwellの時刻・pending command・イベント履歴を含むcheckpointから再開します。境界サンプルは解析時に1つだけ残し、0–30日を再計算しません。','',
        '### 巨大な旧runがないGit取得直後の再現','',
        'scripts/reproduce_dispersion_control.pyは凍結した3配置の制御ON/OFFだけを計算する再現専用です。この研究では追加実験として実行していません。旧スクリーニング全37条件を計算し直す必要はありません。30日結果の順序・精度・負荷を確認してから、必要なら60日へ進めます。','',
        '```powershell\n$env:PYTHONUTF8="1"\n.analysis-venv/Scripts/python.exe -X utf8 scripts/reproduce_dispersion_control.py --days 30 --output-root outputs/control_reproduction\n$env:SAT_CONTROL_OUTPUT="outputs/control_reproduction"\n$env:SAT_CONTROL_DOC="docs/control_reproduction"\n.analysis-venv/Scripts/python.exe -X utf8 scripts/analyze_dispersion_control.py --days 30\n# 30日の結果を確認した後のみ実行\n.analysis-venv/Scripts/python.exe -X utf8 scripts/reproduce_dispersion_control.py --days 60 --output-root outputs/control_reproduction --prefix outputs/control_reproduction/batch_30d\n.analysis-venv/Scripts/python.exe -X utf8 scripts/analyze_dispersion_control.py --days 60\n.analysis-venv/Scripts/python.exe -X utf8 scripts/audit_dispersion_control_density.py --days 60\n.analysis-venv/Scripts/python.exe -X utf8 scripts/plot_dispersion_control.py --days 60\n```','',
        'Python環境は既存pyproject.tomlのtest/nrlmsis追加依存（pip install -e ".[test,nrlmsis]"）に従ってください。この解析で使用した版はnumpy2.4.6、matplotlib3.11.2、pytest9.1.1、pymsis0.12.0です。特にNRLMSIS version2.0を維持してください。宇宙天気CSVはconfigs/dispersion_control_transfer.yamlが指す既存ファイルを使います。SAT_CONTROL_OUTPUT/SAT_CONTROL_DOCは別保存先を選ぶための環境変数です。再現後は解除すれば元の研究成果を参照できます。','',
        '## 9. 検証・Git','', '[検証結果](data/verification.json)と[Git情報](GIT.md)を参照してください。標準3ケースの24時間回帰に加え、分散ゼロの新しい実験用コードを24時間計算し、物理時系列・イベントが旧coreの保存結果と一致することを確認しました。']
    (DOC/'README.md').write_text('\n'.join(text)+'\n',encoding='utf-8')
    print('Technical appendix written; summary requires final scientific review')


if __name__=='__main__':main()
