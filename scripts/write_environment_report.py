"""Publish a short research-decision summary and numerical technical appendix."""
import json
from pathlib import Path
import sys
import numpy as np
import yaml
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from sat_string.nrlmsis import OrbitGeometry
from analyze_environment_comparison import ROOT,OUT,DOC


def main():
    results=json.loads((DOC/'data'/'comparison.json').read_text(encoding='utf-8'))['results']
    protocol=json.loads((OUT/'protocol.json').read_text());obs=json.loads((DOC/'data'/'grace_integrated_summary.json').read_text())
    maximum=max(r['max_error_m'] for r in results)
    if maximum<100:
        classification='D';decision='制御最適化を主題から外す'
        interpretation='今回の60日・固定基準半径の定義では、3環境とも100 m未満でした。「100 km列の間隔維持に制御が必須」という包括的な主張は撤回すべきです。'
    elif maximum<300:
        classification='C';decision='制御最適化を主題から外す（要求境界・適用限界の確認を先行）'
        interpretation='100 mを超える環境はありますが、km級ではありません。「100 km列には制御が必須」という一般的主張は支持されません。精度要求によって必要性が変わる領域なので、制御性能の最適化より適用限界を先に調べるべきです。'
    elif maximum<1000:
        classification='B';decision='条件付きで戻る'
        interpretation='数百m級の環境があり、要求精度との比較が必要です。制御必須とは断定せず、実要求と観測側の積算差を確認してから戻ります。'
    else:
        classification='A';decision='戻る（経験モデルでの支持と実観測の検証を区別）'
        interpretation='経験モデルでkm級が出ました。間隔維持の問題はモデル内で存在しますが、実観測による検証は別途必要です。'
    names=dict(quiet='静穏側',active='活動期',strong='強嵐期')
    table='| 環境・開始日 | 最大間隔誤差 | max|ΔJρ| [kg s/m³] |\n|---|---:|---:|\n'
    # Escape the inner pipe so Markdown keeps a three-column table.
    table=table.replace('max|ΔJρ|','max abs(ΔJρ)')
    for m in results:table+=f"| {names[m['group']]}・{m['epoch'][:10]} | {m['max_error_m']:.1f} m | {m['max_cumulative_adjacent_Jrho']:.2e} |\n"
    audit=json.loads((OUT/'quiet_integrated_cache_audit.json').read_text())
    audits={g:json.loads((OUT/f'{g}_integrated_cache_audit.json').read_text()) for g in names}
    worst_cache=max(r['max_estimated_spacing_error_m'] for r in audits.values())
    summary=f'''# 研究判断用サマリー

## 1. 一番重要な結論

NRLMSIS 2.0・全LOWの60日最大誤差は、静穏側{results[0]['max_error_m']:.1f} m、活動期{results[1]['max_error_m']:.1f} m、強嵐期{results[2]['max_error_m']:.1f} mでした。{interpretation}実観測から100 km間隔の60日誤差まで確定したわけではありません。

## 2. PHASE 1：TADの時間変化

波の周期や移動速度は確認できましたが、個別構造の振幅が4〜5分でどれだけ変わるかを観測で拘束できませんでした。「5分で密度増加」という重要な例は数値モデルの引用でした。不変と証明したのではなく、任意の寿命を置く根拠がないため、新しい時間変化Gaussian計算は見送りました。

## 3. PHASE 2：宇宙天気だけで事前選定

{table}
ΔJρは隣接衛星の累積密度差で、表は全時刻の最大です。終了値だけで位置ずれは決まらず、その差が続いた時間も効きます。太陽活動・季節・初期地方時も異なるため、地磁気活動だけの効果とは言いません。補間の追加監査では最大約{worst_cache:.2f} mの履歴差が見積もられ、全履歴1%基準を満たさない期間がありました。細かい桁を物理的精度として扱いません。

## 4. PHASE 3：観測から言えた範囲

キャッシュ済みGRACE A/Bから、地磁気指標で選んだ各1日の同時刻差・近傍同地点通過差・有効区間の積算差を取得しました。ただし各値は付属密度誤差の積算包絡内で、欠測もあります。持続する実密度差としては識別できません。約220 km対・3度平均のデータを100 kmの直接証拠には使いません。

## 5. 制御研究へ戻るべきか

**{decision}。** 事前基準では{classification}に近い判定です。研究の中心は当面「熱圏密度変動の影響と差動抗力制御の適用限界」とします。実環境で永久に制御不要と証明した、という意味ではありません。

## 6. 次の作業

1. 固定基準間隔か現在半径の弧長かを決め、実際の許容誤差を確認する。
2. 60日補間の積算精度を改善し、必要なら代表ケースのみ180日へ延長する。
3. GRACE A/Bの相互校正と欠測の扱いを整え、実積算差を識別する。
'''
    (DOC/'report.md').write_text(summary,encoding='utf-8')
    appendix='# 詳細技術付録\n\n'
    appendix+='## 1. 観測・解釈・未確定\n\n**計算で得られたこと：** 以下は既存平均軌道モデルとNRLMSIS 2.0の保存軌道の後処理値である。観測密度は加速度計観測から推定された量である。\n\n**そこから言えること：** 選定した環境・軌道・期間での差動抗力と間隔誤差の規模を判断できる。積算差と速度の関係は元ODEの恒等関係でもあるため、独立な観測による因果証明と区別する。\n\n**まだ言えないこと：** 全軌道面・実大気の上限、永久安定性、実ミッション要求、実測100 km対の60日誤差は未確定。NRLMSISは平均的な経験モデルで、個々のTADを再現しない。[原論文](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2020ea001321)。\n\n'
    appendix+='## 2. PHASE 1の停止判断\n\n[原典別の確認結果](PHASE1.md)を参照。新しい成長・減衰モデルは0ケース。波周期と包絡成長を混同しない。既存1回通過の結果を覆すために任意の時間窓は置いていない。隣接100 kmの通過時差は旧1D比較で約13〜14秒、列全体では約4〜5分であり、両者も区別する。\n\n'
    appendix+='## 3. 宇宙天気選定・固定入力\n\n選定は2024年Jan/Mar/May/Jul/Sep/Novの1日から各60日を候補とし、各クラスで観測入力が完全な最初の期間を採用。quietは最大Kp<5かつ最低Dst>−50、activeはKp>=5の日が2日以上かつ−250<最低Dst<=−100、strongは最低Dst<=−250 nT。これは研究上の区分で、公式嵐分類やSYM-HをDstへ置き換えた区分ではない。[全候補と凍結driver](data/weather_selection.json)。\n\n| 環境 | UTC開始〜終了（終了を含まない） | 最大Kp | 最低Dst [nT] | Kp>=5の日数 | 平均観測F10.7 [sfu] | 初期中央地方時 [h] |\n|---|---|---:|---:|---:|---:|---:|\n'
    for m in results:
        g=m['group'];r=protocol['selected'][g];spec=yaml.safe_load((ROOT/f'configs/real_environment_applicability/{g}.yaml').read_text())
        geo=OrbitGeometry(**{k:spec[k] for k in ('epoch_utc','inclination_deg','raan_deg','initial_phase_deg')});lon,_,_=geo.coordinates(0,6778137.,0.)
        appendix+=f"| {names[g]} | {r['start']}〜{r['end_exclusive']} | {r['max_kp']:.1f} | {r['min_dst_nT']:.0f} | {r['days_kp_ge5']} | {r['mean_f107_observed']:.1f} | {(float(lon)/15)%24:.2f} |\n"
    appendix+='\nJanuaryは既存runを再利用した。normalized物理・数値設定、初期状態、保存時刻は機械的に一致確認した。変わるのはepochと宇宙天気系列。N21、400 km、100 km、4 kg、Cd2.2、LOW0.01 m²、RK4最大2秒・保存10秒、傾斜51.6度・RAAN0を固定し、制御・人工外乱はOFF。F10.7/centered81と7成分Apは[CelesTrak公式仕様](https://celestrak.org/SpaceData/SpaceWx-format.php)、Dst/Kpは[NASA OMNI仕様](https://omniweb.gsfc.nasa.gov/html/ow_data.html)を確認。入力を公開CSVへ凍結し、実際の系列に欠損補完はない。\n\n'
    appendix+='## 4. 60日結果と誤差定義\n\n主要指標は従来の固定基準半径誤差。現在半径の弧長は共通高度低下も含む別の補助量である。[定義・式・限界](metrics_and_limits.md)。\n\n| 環境 | 最大 [m] | RMS [m] | 終了最大 [m] | 100km比 [%] | 最大速度 [µm/s] | 終了max abs(ΔJρ) [kg s/m³] | 弧長最大 [m] | 平均高度損失 [m] |\n|---|---:|---:|---:|---:|---:|---:|---:|---:|\n'
    for m in results:appendix+=f"| {names[m['group']]} | {m['max_error_m']:.3f} | {m['rms_error_m']:.3f} | {m['final_max_error_m']:.3f} | {m['max_error_percent_of_100km']:.4f} | {m['max_gap_speed_m_s']*1e6:.3f} | {m['final_max_adjacent_Jrho']:.3e} | {m['current_radius_arc']['max_error_m']:.3f} | {m['mean_altitude_loss_m']:.1f} |\n"
    appendix+='\nこの表の細かい桁は保存軌道の再現用で、実環境や直接モデル解の精度を意味しない。高度損失は初期半径−終了半径を21機平均した通常の全LOW損失であり、「制御の追加高度損失」ではない。\n\n![間隔比較](figures/02_spacing.png)\n\n'
    appendix+='| 環境 | 終了max abs(ΔJdrag) [m/s] | 終了最大gap速度 [m/s] | 最後14日・最大誤差の傾き [m/day] |\n|---|---:|---:|---:|\n'
    for m in results:appendix+=f"| {names[m['group']]} | {m['final_max_adjacent_Jdrag_m_s']:.3e} | {m['final_max_gap_speed_m_s']:.3e} | {m['tail']['end14day_error_growth_m_s']*86400:.3f} |\n"
    appendix+='\nNRLMSISは継続して作用しており、終了時速度を「外乱後一定の残留速度」として60日以降へ直線外挿しない。初期地方時・季節・太陽活動が違うため、強い嵐ほど必ず最大誤差が大きいとの法則も主張しない。最後14日の正の傾きがある期間は180日を調べる価値があるが、今回180日・1年runは実施せず、まず数値精度・要求・観測の識別性を優先する。\n\n![密度差と状態](figures/03_mechanism.png)\n\n'
    appendix+='## 5. 仮想許容誤差の超過\n\nこれは実ミッション要求ではない。10秒標本のstrict >、時間割合は左端標本保持、最後の標本は長さ0で計算した。\n\n| 環境 | 許容 [m] | 初回超過 [day] | 超過 [%] | 最大連続超過 [day] |\n|---|---:|---:|---:|---:|\n'
    for m in results:
        for r in m['virtual_tolerances']:
            first='なし' if r['first_exceedance_time_s'] is None else f"{r['first_exceedance_time_s']/86400:.3f}"
            appendix+=f"| {names[m['group']]} | {r['tolerance_m']:.0f} | {first} | {r['exceedance_fraction']*100:.3f} | {r['longest_continuous_exceedance_s']/86400:.3f} |\n"
    appendix+='\n現在半径弧長での超過統計は各 `data/*_metrics.json` に別保存。共通の軌道収縮を相対密度差による誤差と呼ばない。\n\n![許容値](figures/04_tolerances.png)\n\n'
    appendix+='## 6. PHASE 3のgo判断と観測支持区間\n\nPHASE 2だけでは個別構造を含む実環境を判断できないため、キャッシュ済み観測の小規模後処理へ進んだ。新しい21機観測入力シミュレーションは0件。[Sutton v2.3製品仕様](https://zenodo.org/records/13287614)の3度平均GRACE A/Bを利用。加速度計観測から推定された密度であり、100 kmを直接解像したデータではない。気候学的分布や上位percentileは推定しない。\n\n選定は既存driver supportの2003-10-25〜11-23内で、地磁気指標だけから各クラスの最初の1日を採用。[選定一覧](data/grace_window_selection.json)。quietは地磁気のみを意味し太陽活動は高い。disturbedはHalloween嵐の回復期であり、3つの独立嵐ではない。\n\n| 観測日・区分 | 被覆 [h] | 最長連続 [min] | 同時刻J_A−J_B [kg s/m³] | 同地点近似J差 | 誤差包絡 | 一定校正ずれ相当 [%] |\n|---|---:|---:|---:|---:|---:|---:|\n'
    for m in obs['results']:
        r=m['simultaneous_obs'];appendix+=f"| {m['window']} | {r['coverage_s']/3600:.2f} | {r['max_contiguous_coverage_s']/60:.2f} | {r['integral_observed_support']:.3e} | {m['matched_obs']['integral_observed_support']:.3e} | {r['fully_correlated_error_envelope']:.3e} | {m['constant_relative_calibration_offset_equivalent']*100:.3f} |\n"
    appendix+='\nこれらは有効支持区間の和で、欠測を含む1日全体や60日の持続積算差ではない。σの包絡、位置合わせ後の残る数kmの差、約220 kmの対間隔、300〜350 km平均、400 km規格化の仮定がある。観測由来値が同じ位置のモデルより大きくても、持続する100 km差動密度を識別したとは言えない。各区間の符号を残し、平均を引いて積算差を強制的に0にしていない。詳細は `grace_integrated_summary.json` と `grace_segments.csv`。\n\n![観測支持と不確かさ](figures/05_grace_support.png)\n\n'
    appendix+='## 7. 数値検証・回帰・テスト\n\n密度120秒キャッシュは各実runの保存位置と内部時刻で直接NRLMSIS2.0に照合した。すべて0.001相対誤差基準内。最大3時間apの約3軌道について10/30/60/120秒と直接2秒密度の収束も確認し、120秒の間隔への影響診断は約0.000508 mだった。ただし、これは60日の物理軌道の収束証明ではない。\n\n'
    appendix+=f"静穏期の全{audit['all_saved_points']:,}密度点の直接再評価では、最大相対密度誤差{audit['max_relative_density_error']:.3e}、最大累積隣接誤差{audit['max_cumulative_adjacent_density_error_kg_s_m3']:.3e} kg s/m³、凍結軌道による間隔差診断{audit['max_estimated_spacing_error_m']:.3f} m。追加の全履歴1%基準（{audit['acceptance_estimated_spacing_error_m']:.3f} m）は**不合格**。これは単体テストや24時間回帰の失敗ではない。既存runを勝手に置き換えず、この精度限界を保持した。軌道・密度feedbackを再伝播した誤差上界ではなく、100 m近傍の厳密境界を議論する前に改善が必要。\n\n"
    appendix+='新しい積算恒等式・符号・欠測非補間・誤差包絡・宇宙天気分類の14テストを含め、**245テスト合格**。24時間referenceは保存済み再計算suiteから集計・全イベント・主要時系列を照合し、3ケースとも合格。core変更はない。旧成果8404ファイルのSHA256保持検査を別保存する。\n\n'
    appendix+='| 全保存位置の補間監査 | 間隔差診断 [m] | 追加基準 [m] | 判定 |\n|---|---:|---:|---|\n'
    for g,r in audits.items():
        appendix+=f"| {names[g]} | {r['max_estimated_spacing_error_m']:.3f} | {r['acceptance_estimated_spacing_error_m']:.3f} | {'合格' if r['passed'] else '不合格'} |\n"
    appendix+='\n各監査JSONとソフトウェア検証記録は `data/` に保存した。補間を含む長期精度の厳密な合格とは扱わない。特に境界近くの要求判定には、より細かい密度キャッシュまたは直接評価軌道との収束確認が必要。\n\n'
    appendix+='## 8. 最終10問\n\n'
    answers=[
      '4〜5分で全構造が不変とは言えない。局所の包絡変化率を観測で拘束できず、任意の成長・減衰modelは作らなかった。',
      '単発・時間不変TADの1回通過を大きな長期誤差の主要因とする根拠は今回残らない。継続forcingや急速生成を一般に否定するものではない。',
      f"静穏側は{results[0]['max_error_m']:.3f} m。",f"活動期は{results[1]['max_error_m']:.3f} m。",f"強嵐期は{results[2]['max_error_m']:.3f} m。",
      '元ODEでは積算密度から半径・速度への関係が導ける。位置誤差にはその速度履歴の積算が必要で、最終ΔJρや3期間の最大値だけでは説明しない。',
      '観測由来支持区間の差は推定できたが、付属密度誤差の積算包絡内で欠測もあり、NRLMSISを上回る持続した実差動積算密度は識別できていない。',
      f'この平均モデルで得た最大は{maximum:.1f} m。100 km対・400 kmの実観測から数百m〜km級を確定したわけではない。',
      interpretation,
      decision+'。追加の制御パラメータ計算は0件。']
    for i,a in enumerate(answers,1):appendix+=f'{i}. {a}\n'
    appendix+='\n## 9. Git・再現・変更内容\n\nbranch：`astra/real-environment-applicability`。mainへのmergeなし。新しい物理runは活動期・強嵐期の2件のみ、既存静穏期を再利用。新しい追加は `environment_integrals.py`、選定・集計・密度精度監査・GRACE積算・描画・文書出力scripts、3条件の環境YAML、14 tests、この報告書・小さいCSV/JSON/図。既存core/config/runは変更していない。\n\n[再現手順と主要コマンド](reproduction.md)。生観測・大NPZはローカルのみ。最終commit hashは `git log -2 --oneline` と完了回答で確認できる。\n'
    (DOC/'technical_appendix.md').write_text(appendix,encoding='utf-8')
    print('Summary characters',len(summary),'max_error_m',maximum,'decision',decision,flush=True)


if __name__=='__main__':main()
