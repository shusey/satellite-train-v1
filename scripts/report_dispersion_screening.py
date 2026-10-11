"""Write a two-layer Japanese research report from compact saved diagnostics."""
import json
import numpy as np
from run_dispersion_screening import ROOT,OUT,DOC,DATA,configuration


def main():
    spec,cfg=configuration(7)
    groups=json.loads((DATA/'group_summary.json').read_text())['groups']
    rows=json.loads((DATA/'metrics.json').read_text())['results']
    mech=json.loads((DATA/'mechanism_summary.json').read_text())
    audit=json.loads((DATA/'density_spot_audit.json').read_text())
    regression=json.loads((OUT/'regression_24h'/'regression.json').read_text())
    preservation=json.loads((OUT/'preservation_check'/'preservation.json').read_text())
    def group(f,s):return next(g for g in groups if g['factor']==f and g['sigma']==s)
    def val(f,s,d):return group(f,s)['horizons'].get(str(d))
    def median(f,s,d):return val(f,s,d)['median_max_error_m']
    def horizon(h):
        if h is None:return '省略：7日で十分大きい'
        suffix=f"（{h['n']}通り）"
        return f"{h['median_max_error_m']:,.1f} [{h['min_max_error_m']:,.1f}–{h['max_max_error_m']:,.1f}] {suffix}" if h['n']>1 else f"{h['median_max_error_m']:,.1f}（代表1通り）"
    def day(t):return '未超過' if t is None else f'{t/86400:.3f}日'
    base={int(r['days']):r for r in rows if r['factor']=='zero'}
    n=cfg.derived.reference_mean_motion_rad_s
    text=f'''# 研究判断用サマリー

## 1. 結論

**今回の暫定範囲では、初期半長軸のずれと衛星ごとの抗力の違いは、両方とも本格調査する価値があります。** 初期半長軸差は開始直後から間隔をずらし、抗力の個体差は日数とともに影響を増しました。1日中央値では半長軸0.1 mが約{median('a',.1,1):.0f} m、抗力比0.3%が約{median('B',.003,1):.0f} mと近いものの、7日では後者がより大きくなりました。実衛星での分散は未確認です。

## 2. 主要比較表

全ケースは21機・100 km間隔・全LOW・制御なしです。単位はm、7日は乱数5通りの中央値、括弧内は最小–最大です。30日は必要と判断した条件だけ延長しました。

| 要因 | 指定分散 | 7日最大誤差 | 30日最大誤差 | 暫定判定 |
|---|---:|---:|---:|---|
'''
    for g in groups:
        label='初期半長軸' if g['factor']=='a' else 'C_D A/m'
        sigma=f"{g['sigma']:g} m" if g['factor']=='a' else f"{g['sigma']*100:g}%"
        seven=horizon(g['horizons']['7']).replace(' （5通り）','')
        text+=f"| {label} | {sigma} | {seven} | {horizon(g['horizons'].get('30'))} | A※ |\n"
    text+='\n※暫定レンジ内でkm級に達したという判定で、入力値の現実性は未確認です。\n'
    text+=f'''
## 3. 初期半長軸分散

0.1 mでも7日に約{median('a',.1,7):.0f} m、1 mでは約{median('a',1,7)/1000:.2f} km、10 mでは約{median('a',10,7)/1000:.1f} kmになりました。最初の100 km間隔を乱した結果ではありません。軌道の大きさが異なることで周回速度が少し違い、その差が時間とともに位置ずれになります。0.1 m条件も30日にkm級へ達しました。

## 4. 抗力の個体差

0.3%でも7日に約{median('B',.003,7)/1000:.2f} km、1%で約{median('B',.01,7)/1000:.2f} kmでした。5–10%はさらに大きく、延長を省きました。同じ密度でも抗力が違うと高度の低下量が異なり、やがて周回速度と間隔がずれます。

## 5. 大気変動だけとの比較

分散なしの同じ7日は{base[7]['max_error_m']:.2f} m、30日は{base[30]['max_error_m']:.2f} m、60日は{base[60]['max_error_m']:.2f} mでした。今回の分散による誤差は、期間を揃えてもこれを大きく上回ります。「NRLMSISの密度差だけでは小さい」という結果から、個体差を含む編隊全体にも制御が不要とは言えません。一方、この実験だけで制御必要性が証明されたとも言えません。

## 6. 次の判断

**両方を本格調査へ進めます。** 長期の個体差としては C_D A/m を、編隊形成後に残すべき精度としては相対半長軸を調べます。次は、①形成後の半長軸残差・推定精度、②同型衛星の質量・投影面積・抗力係数のばらつき、③実際の許容間隔誤差、の根拠を確認します。それまで今回の乱数分布を卒論の正式な現実モデルとはしません。

# 詳細技術付録

## A. 観測事実・解釈・未確定事項

### 観測されたこと

ここでいう観測事実は、この数値モデルで得た計算結果です。衛星の実観測値ではありません。両要因について、暫定入力の小さい側でも数百m–km級の間隔誤差が生じました。全て全LOWで、姿勢遷移数は0です。

### そこから言えること

全衛星の半長軸と C_D A/m が完全に同じという理想化は、今回のモデルで長期間隔誤差を小さくする重要な条件でした。数百m–km級のずれを作る要因候補として両方を残す合理性があります。

### まだ言えないこと

0.1 mや1 m、0.3%や1%が実衛星の典型的な標準偏差であること、許容誤差に反すること、制御で抑制できることは確認していません。有限30日での増加を数学的な「発散」とは呼びません。5乱数は母集団の確率や上位百分位を推定するには不足しています。

今回は衛星ごとのB倍率を時間一定としました。全衛星に共通する抗力モデルの誤差や、平均0へ戻る短時間の姿勢変動を、そのまま独立で永続するB差とは扱えません。また軌道推定の不確かさと、実際に残っている半長軸差は別の量です。次の文献調査ではこの区別を確認します。

## B. 入力固定と計算方法

元の active 60日runの正規化設定を読み、観測終了時刻だけを1/7/30日へ変えました。設定辞書を機械比較し、他の物理・数値値が一致することを確認しています。個体差は別の配列として与え、各時点の実際の位置・高度でNRLMSIS 2.0を評価します。元runの密度履歴を全ケースへ固定して流用してはいません。

- epoch: 2024-03-01 00:00 UTC、既存active期間の先頭部分。
- N=21、基準半径 {cfg.derived.reference_radius_m:.0f} m、高度400 km、目標間隔100 km。
- m=4 kg、C_D=2.2、LOW面積0.01 m²、B0=0.0055 m²/kg。
- 元の平均軌道ODE `orbital_rates` とRK4、最大刻み2 s、保存10 s。
- 元の軌道幾何：傾斜角51.6°、RAAN0°、初期位相0°。J2等の追加なし。
- NRLMSIS 2.0 / pymsis 0.12.0、凍結済みF10.7/Ap CSV、既存120 s密度キャッシュ方式。
- 人工Gaussian外乱はOFF。制御器・姿勢状態機械を呼ばずLOW面積を常時使用。

独立ケース37本（乱数35本＋交互配置2本）をまとめて計算しました。他ケースと力学的に結合しません。各ケース21機の密度は変更された軌道位置から求めます。Bの分散は元の抗力式に B_j/B0 を掛けることで実装しており、「物理的に C_D だけが違う」という仮定にはしていません。

## C. 乱数・定義・全評価値

乱数番号は0–4、半長軸とBの標準化パターンは独立です。B側の実際の乱数番号は10000–10004です。同じ要因の大小条件では同じパターンを倍率変更しました。生成後に平均を引くだけで、実現した標準偏差を指定値へ再調整していません。N=21なので実現値にはばらつきがあります。

- `delta_a_m` は基準半長軸からの差、`epsilon` はB0からの相対差。
- 初期lambdaは元の100 km配列を維持。全衛星平均の半長軸/Bは元の値を維持。
- 間隔誤差は従来の固定基準半径定義：e_i=a_ref[(lambda_(i+1)-lambda_i)-100km/a_ref]。
- 全時間・全20間隔RMSは sqrt(mean(e_i(t)^2))。10 s保存時刻を等しく重み付け。
- 最大間隔速度は max_(i,t)|de_i/dt|。差分微分ではなく元の角速度から算出。
- 初回超過は max_i|e_i| が100/500/1000 mを**厳密に上回る最初の保存時刻**。10 sの判定分解能。
- 百分率は誤差/100000 m×100。終了時最大誤差は収束判定ではない。

[入力と全衛星ベクトル](data/profiles.json)、[衛星別CSV](data/satellite_profiles.csv)、[全metrics JSON](data/metrics.json)、[全metrics CSV](data/metrics.csv)、[集約値](data/group_summary.json) に全結果があります。

## D. 初期半長軸差と理論の照合

微小な半長軸差に対して、n=sqrt(mu/a³)を一次近似すると、

de_i/dt ≈ −(3/2)n0[delta_a_(i+1)−delta_a_i]。

これは、後ろ側番号の衛星の半長軸が大きいと、その衛星の角速度が小さくなり、符号付き間隔が縮むことを意味します。符号は元のlambda順序に合わせています。n0={n:.10g} rad/s、隣接半長軸差1 mあたり約{1.5*n*1000:.4f} mm/s、抗力を無視した1日あたり約{1.5*n*86400:.2f} mです。sigma_aと隣接差は同じ量ではありません。

保存初期状態の厳密なKepler角速度との差は、最大絶対誤差 {mech['theory_max_absolute_error_m_s']:.3g} m/s、最大相対誤差 {mech['theory_max_relative_error']:.3g} でした。符号・単位が整合しています。長期値は抗力と密度変化も含むので、この一次式の直線と完全一致するとはしません。[全gap理論照合](data/initial_drift_theory.csv)。

## E. C_D A/m の作用

元の平均軌道式は da_j/dt=−B_j rho_j sqrt(mu a_j) です。共通初期半長軸では、B分散だけによる初期の変化量は −B0 epsilon_j rho_j sqrt(mu a0) になります。

抗力比が大きい衛星ほど半長軸が早く小さくなり、その後角速度が大きくなります。この速度差がさらに間隔誤差へ積算されます。一定密度・小さい差という近似なら、初期半長軸差による位置ずれは時間にほぼ比例し、B差が半長軸差を育てる場合は初期には時間の二乗で育つ構造です。本計算の密度は一定ではなく、二乗則は全期間の厳密な法則ではありません。

[代表Bケースの毎時・全gap履歴](data/representative_B_hourly.csv) に ΔB、分散なしを差し引いた隣接半長軸差、間隔速度、間隔誤差があります。[機構要約](data/mechanism_summary.json) は初期B由来の半長軸変化率と7日状態を保存しています。

## F. 超過時刻と残留速度

以下は各指定水準で7日最大誤差が中央値に位置した乱数の、**7日までに限った**結果です。7日内の最大速度を記載し、瞬間の終端速度とは区別しています。

| 要因 | 分散 | 代表乱数 | 100m初回 | 500m初回 | 1km初回 | 7日最大間隔速度 [mm/s] |
|---|---:|---:|---:|---:|---:|---:|
'''
    for g in groups:
        h=g['horizons']['7'];r=next(r for r in rows if r['case_id']==h['representative_case_id'] and r['days']==7)
        s=f"{g['sigma']:g} m" if g['factor']=='a' else f"{g['sigma']*100:g}%"
        crosses=[day(r['first_exceedance_s'][str(k)]) for k in (100,500,1000)]
        text+=f"| {g['factor']} | {s} | {r['seed']} | {' | '.join(crosses)} | {r['max_gap_rate_m_s']*1000:.3f} |\n"
    text+='\n30日へ延長した代表乱数では、次の結果でした。Aは乱数4、Bは乱数2です。全乱数のRMS・終了値・超過時刻はmetrics CSVに保存されています。\n\n| 要因・分散 | 30日最大/終了誤差 [m] | 全時間RMS [m] | 100m初回 | 500m初回 | 1km初回 |\n|---|---:|---:|---:|---:|---:|\n'
    for f,s in [('a',.1),('B',.003)]:
        h=val(f,s,7);r=next(r for r in rows if r['case_id']==h['representative_case_id'] and r['days']==30)
        crosses=[day(r['first_exceedance_s'][str(k)]) for k in (100,500,1000)]
        text+=f"| {f} {s:g} | {r['max_error_m']:,.2f} / {r['final_max_error_m']:,.2f} | {r['rms_error_m']:,.2f} | {' | '.join(crosses)} |\n"
    text+='''
## G. 延長判断・交互配置・組合せ

7日で100–500 mの乱数は30日へ継続し、7日で十分大きい水準の全乱数延長を省きました。ただし最小B分散は、長期の成長を確認するため中央値の乱数を1本だけ30日へ継続しました。保存された7日の物理状態から続きを計算し、最初の7日を再計算していません。50 m条件・追加乱数・60日計算・制御ありケースは実施していません。判定内容は [extension_decision.json](data/extension_decision.json) にあります。

交互配置は平均を引いた±パターンで、半長軸1 m、B1%を各1本実行しました。確率分布でも最悪ケースの保証でもありません。ランダム配置では隣接差が指定sigmaの2倍を超えることもあるため、交互配置が最大値になるとは限りません。

| 負荷確認 | 7日最大誤差 [m] | RMS [m] |
|---|---:|---:|
'''
    for r in rows:
        if r['kind']=='stress' and r['days']==7:text+=f"| {r['factor']}：{r['nominal_sigma']:g} | {r['max_error_m']:,.2f} | {r['rms_error_m']:,.2f} |\n"
    combo=[r for r in rows if r['factor']=='combined' and r['days']==7]
    if combo:
        r=combo[0];deviation=json.loads((DATA/'combination_comparison.json').read_text())['superposition_deviation_max_m']
        text+=f'''
両要因が有意だったため、半長軸1 m＋B1%を独立パターンの1組だけ追加しました。7日最大誤差 {r['max_error_m']:,.2f} m、RMS {r['rms_error_m']:,.2f} mでした。各単独ケースを線形重ね合わせした近似からの最大差は {deviation:.3f} mです。各gapの符号次第で打消し・強め合いがあり、最大誤差同士を単純加算できません。1組のみなので組合せ分布は評価していません。[入力](data/combination_profile.json)。
'''
    reversed_rows=[r for r in rows if r['days'] in (7,30) and r.get('first_label_order_reversal_s') is not None]
    if reversed_rows:
        text+='\n大きい側の次のケースは、固定番号間の符号付き間隔が0以下となり、元の衛星順序を失いました。以後の値は元の衛星番号間の誤差であり、並べ替え後の実際の近傍間隔ではありません。上限スクリーニングの結果として扱い、健全な編隊維持状態の指標とはしません。\n\n'
        for r in reversed_rows:text+=f"- {r['case_id']}（{r['days']:g}日）：初回 {day(r['first_label_order_reversal_s'])}、最小符号付き間隔 {r['min_labelled_spacing_m']:,.1f} m。\n"
    else:
        text+='\n全解析期間で元の衛星番号順序は保たれました。固定した衛星番号間の間隔を自動的に並べ替えてはいません。\n'
    text+=f'''
## H. 同じ観測期間での比較

| 分散なし・NRLMSISのみ | 最大誤差 [m] | RMS [m] |
|---|---:|---:|
'''
    for d,r in base.items():text+=f"| {d}日 | {r['max_error_m']:.3f} | {r['rms_error_m']:.3f} |\n"
    text+='\n| 要因・指定分散 | 7日中央値 / 分散なし7日 | 7日中央値 / 113.13m（60日） | 同じ7日での規模分類 |\n|---|---:|---:|---|\n'
    for g in groups:
        value=g['horizons']['7']['median_max_error_m'];ratio=value/base[7]['max_error_m']
        category='2桁以上大きい' if ratio>=100 else ('1桁大きい' if ratio>=10 else ('同程度' if ratio>=.1 else '小さい'))
        text+=f"| {g['factor']} {g['sigma']:g} | {ratio:.1f}倍 | {value/base[60]['max_error_m']:.1f}倍 | {category} |\n"
    text+='''
60日の113 mは背景として併記します。分散ケースの7日値を113 mだけで比較せず、11.75 mの同じ7日とも比べてください。条件ごとのRankは入力の現実性を認定するものではなく、この暫定範囲で生じた大きさを選別するものです。

## I. 図の読み方

1. [半長軸分散と最大誤差](figures/01_initial_a.png)：5乱数の中央値と範囲。30日の試行数は各水準で異なります。
2. [抗力比分散と最大誤差](figures/02_B_dispersion.png)：同じ期間・同じ乱数水準で比較。30日のBは選んだ1乱数の成長例です。
3. [初期半長軸0.1 mの履歴](figures/03_initial_drift.png)：速度差が最初からあり、位置ずれが積み上がる様子。破線は初期速度を一定とした一次予測です。
4. [B0.3%の機構](figures/04_B_mechanism.png)：半長軸差→速度差→位置差の順で影響が育ちます。上段は同時刻の分散なしを差し引いています。
5. [小さい入力と分散なしの比較](figures/05_baseline_comparison.png)：観測期間を揃え、密度差だけより個体差の影響が大きいことを示します。

## J. 10の問いへの回答

'''
    for num,f,s in [(1,'a',.1),(2,'a',1),(3,'a',10),(4,'B',.003),(5,'B',.01)]:
        text+=f"- Q{num}：7日中央値 {median(f,s,7):,.1f} m。1日は {median(f,s,1):,.1f} m。\n"
    text+=f'''- Q6：B5%は7日中央値 {median('B',.05,7):,.1f} m、10%は {median('B',.1,7):,.1f} m。
- Q7：両要因とも小さい指定水準でもNRLMSIS-onlyを上回ります。同期間での比較は表H、全乱数の比は元metricsから算出できます。
- Q8：この暫定レンジなら500 m–1 km級は生じます。ただしレンジの実現頻度は未確認で、「現実の個体差で必ず起きる」とは言えません。
- Q9：両方を本格文献調査へ。個体差としてはB、形成後の精度要求としては相対半長軸。
- Q10：編隊維持の必要性を再検討する根拠候補です。実際の残差・個体差・許容値を確認してから、制御研究へ戻るか決めます。

## K. 精度・テスト・回帰・既存成果の保全

NRLMSISを実際の変更軌道で直接評価した密度と、保存キャッシュ密度を約6時間ごと、さらにその10/50/90 s後の計 {audit['points']} 点で照合しました。キャッシュ区切り時刻だけでなく、補間区間の途中を含みます。最大相対差 {audit['max_relative_error']:.3g}、基準0.001未満で合格です。これは疎な点検であり、全軌道伝播誤差の上限保証ではありません。元研究の120 sキャッシュ方式を使い、今回のkm級効果を解釈する際にこの限界を残します。[監査](data/density_spot_audit.json)。

全テストの件数とログは [verification.json](data/verification.json) にあります。新規テストでは平均0の生成、正のB倍率、初期lambda固定、Kepler理論、一定密度の解析解、独立ケースのまとめ計算と単独計算の一致、保存状態からの継続を検証しました。

既存24時間software referenceは3ケースとも {'合格' if regression['passed'] else '不合格'}。既存の再実行済み時系列を再利用し、最大誤差/RMS/高度差/左右波速度・全切替イベント・主要時系列を凍結referenceと照合しています。新しい解析は既存の制御・物理・積分器ファイルを変更していません。[回帰照合](data/regression_24h.json)。

作業前にSHA256を保存した旧ファイル {preservation['old_files_checked']} 個について、変更・欠損 {'0' if preservation['passed'] else 'あり'} を確認しました。新runは独立ディレクトリだけへ保存し、旧runへの追加や上書きはありません。[保全](data/preservation.json)。

## L. 変更ファイル・再現・Git

- `src/sat_string/passive_dispersion.py`：元ODE/RK4を再利用する全LOWの個体差オプションと評価。
- `tests/test_passive_dispersion.py`：独立性・初期条件・理論・継続等の検証。
- `configs/initial_dispersion_screening.yaml`：暫定範囲、乱数、保存・延長ルール。
- `scripts/run_dispersion_screening.py`：既存baseline再利用、7日・継続・組合せ・集計。
- `scripts/analyze_dispersion_details.py`：理論、因果履歴、直接密度監査。
- `scripts/plot_dispersion_screening.py`：日本語5図。
- `scripts/report_dispersion_screening.py`：保存指標からこの文書を生成。
- `scripts/verify_dispersion_screening.py`：設定差分・テスト・回帰・旧ファイル保全の証拠。
- `docs/initial_dispersion_screening/`：小さい設定・集計・図・検証・文書。巨大NPYはローカル `outputs/initial_dispersion_20261011/` に保存しGitへ追加しません。

このローカル環境で保存runから後処理する主要コマンド：

```powershell
.analysis-venv/Scripts/python.exe scripts/run_dispersion_screening.py --stage analyze
.analysis-venv/Scripts/python.exe scripts/analyze_dispersion_details.py
.analysis-venv/Scripts/python.exe scripts/plot_dispersion_screening.py
.analysis-venv/Scripts/python.exe -m pytest -q
```

物理計算を再現する必要がある場合のみ、既存active60日run・凍結weather CSV・同じ依存環境を用意し、新しい保存先を設定します。既存runを重複計算するためには実行しないでください。

```powershell
$env:SAT_DISPERSION_OUTPUT='C:/absolute/path/to/new_outputs'
$env:SAT_DISPERSION_DOC='C:/absolute/path/to/new_report'
.analysis-venv/Scripts/python.exe scripts/run_dispersion_screening.py --stage initial
.analysis-venv/Scripts/python.exe scripts/run_dispersion_screening.py --stage analyze
.analysis-venv/Scripts/python.exe scripts/run_dispersion_screening.py --stage extend
.analysis-venv/Scripts/python.exe scripts/run_dispersion_screening.py --stage combination
.analysis-venv/Scripts/python.exe scripts/run_dispersion_screening.py --stage analyze
.analysis-venv/Scripts/python.exe scripts/analyze_dispersion_details.py
.analysis-venv/Scripts/python.exe scripts/plot_dispersion_screening.py
```

新しい保存先で再計算する場合は7日結果を確認してからextend/combinationへ進みます。既存published図は上書き防止のため再作図時に新しい `SAT_DISPERSION_DOC` を指定してください。元baselineデータがないcheckoutでは、先に既存研究の [再現手順](../real_environment_applicability/reproduction.md) で必要なローカルデータを準備してください。既存データの大半はGit対象外であり、この手順にはそれらの準備が必要です。

branch：`astra/initial-dispersion-screening`。mainへのmergeは行いません。commitとpush情報は [GIT.md](GIT.md) に記載します。
'''
    (DOC/'README.md').write_text(text,encoding='utf-8')
    print('Saved',DOC/'README.md',flush=True)


if __name__=='__main__':main()
