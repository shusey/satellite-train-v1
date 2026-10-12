# 研究判断用サマリー

## 1. 一番重要な結論

**形成後の半長軸標準偏差0.1 mは、非現実的と判明した値ではありませんが、かなり厳しく形成した状態を想定する値です。同型3U衛星の抗力倍率標準偏差0.3%も物理的にはあり得ますが、実機集団の代表値と確認できませんでした。** 両方とも、今後は「観測された標準偏差」ではなく「文献で確認した機構・桁を踏まえた暫定条件」と記載してください。100/500/1000 mについても、100 km列に共通する実要求は見つかりません。主要比較は暫定的に500 m、厳しい側100 m、緩い側1 kmを併記し、ミッション要求とは区別することを推奨します。

## 2. 半長軸残差

CanX-4/5では位相調整・ドリフト停止後、約2.95 km離れた状態の**平均半長軸差0.5 m**が報告されています。これは形成後の具体的な推定状態で、標準偏差ではありません。TanDEM-Xでは半長軸差の推定精度約0.2 mが報告されていますが、こちらは状態を知る精度であり、実際に0.2 m以内へ制御したという意味ではありません。[CanX原論文](https://www.utias-sfl.net/wp-content/uploads/SSC15-I-4-CanX_4andCanX_5PrecisionFormationFlight_MissionAccomplished.pdf)、[TanDEM-X原論文](https://skoge.folk.ntnu.no/prost/proceedings/ifac11-proceedings/data/html/papers/2374.pdf)

つまり、0.1〜1 mを次の感度確認の中心に置く理由はあります。一方、0.01 mを通常の形成後残差として採用する直接根拠は弱く、極めて楽観的な下限確認に限ります。精密GNSS・推進機を使った実績を、推進機なし3Uの達成能力へそのまま移すこともできません。

## 3. 抗力個体差

3Uに近い直接例として、ANSERではLOW姿勢でも二機の自然降下率に違いがあり、運用側は姿勢制御性能の違いを原因候補としています。ただし0.3%という恒常差を測定した論文ではありません。大型同型衛星のTanDEM-Xでは、同じ面積・抗力係数を仮定した個体差が初期1%未満、燃料消費で質量が離れた後は約3%でした。**質量・恒常形状の差が持続外乱になることは支持されますが、3Uの0.3%という分布までは支持されません。** [ANSER](https://doi.org/10.1016/j.actaastro.2026.04.057)、[TanDEM-X運用報告](https://elib.dlr.de/225026/1/Kahle_et_al_ISSFD2026.pdf)

4 kg衛星で0.3%の質量差は約12 gに相当します。これは換算例であり、製造公差の実測値ではありません。質量差なら60日固定は合理的ですが、姿勢・大気組成・推定密度の誤差を全部まとめた倍率を60日一定とする根拠は弱くなります。

## 4. 要求精度

100 km級の実例は、むしろ広い維持範囲を許容します。ANSERの運用範囲は100〜200 km、GRACE-FOは220 km±50 kmです。精密SARのTanDEM-Xには沿軌道200 m RMSの制御要求がありますが、観測目的も距離も異なり、今回の最大誤差へ直接適用できません。[GRACE-FO公式資料](https://www.jpl.nasa.gov/news/press_kits/grace-fo/mission/launch-events/)、[TanDEM-X制御報告](https://issfd.org/ISSFD_2011/S4-Flight.Dynamics.Operations.1-FDOP1/S4_P1_ISSFD22_PF_037.pdf)

## 5. 今後使う条件

|項目|現在値|文献で確認された範囲・意味|次の推奨3水準|入力値への根拠の強さ|
|---|---:|---|---|---|
|衛星別半長軸標準偏差 σa|0.1 m|CanXの形成後二機差0.5 m。GRACEの維持運用には意図的なm級差もある。集団σは未確認|0.01 / 0.1 / 1 m|Engineering assumption。桁の選定はModerate evidence|
|衛星別抗力倍率標準偏差 σB/B|0.3%|3Uの集団σは未確認。大型同型機では1%未満〜約3%の二機差|0.1 / 0.3 / 1%|Engineering assumption。持続差の存在はModerate evidence|
|100 km列の許容最大誤差|100/500/1000 m|ミッションによりm級精密編隊から数十kmの維持範囲まで|厳しい100 / 比較基準500 / 緩い1000 m|Engineering assumption。普遍要求にはしない|

現在の318 m・236 mは、この500 m比較基準なら良好な誤差抑制ですが、100 mを常時満たす性能ではありません。**現行制御研究は条件付きで継続できます。ただし「実際の3U列に必須」とする根拠はまだ不足します。** 次は①想定観測ミッションから精度・運用予算を決める、②質量・LOW/HIGH姿勢別の実効抗力差を確認する、③形成完了の条件を平均半長軸差・相対速度で明示する、の順を推奨します。今回は追加シミュレーションを行っていません。

# 詳細技術付録

調査日：2026-10-12。対象モデル：3U級、400 km、21機、固定基準半径上の100 km間隔。数値の出典は本文・表・図説明まで確認し、抄録だけの数値は推奨値の直接根拠にしていません。ページはPDFの先頭を1ページとして数えます。

## A. Relative semi-major-axis literature

### A.1 まず何が分かったか

形成を終えた後でも、半長軸差を常にゼロにする運用ばかりではありません。小さな差を意図的に残して、許された間隔範囲内をゆっくり往復させる運用もあります。**実際の差・意図した目標・差を推定する精度を分けないと、0.1 mの妥当性を判断できません。**

|資料|機体・高度・形状|段階・方式|確認した平均半長軸関連量|数値の種類と注意|
|---|---|---|---|---|
|A1 CanX-4/5, Bonin et al. 2015|各約6 kg、20 cm立方体。高度は採用箇所に数値なし|位相調整後、冷ガス、差分GPS|2.95 kmの駐機編隊で差0.5 m|p9–10/Table 3。形成後の推定状態1点。50 m精密編隊の残差σではない|
|A2 TanDEM-X, Kahle et al. 2011|大型SAR二機、約505 km、ヘリックス|日毎の推進機による編隊維持|第2マヌーバで差をゼロ近傍へ。平均沿軌道誤差を調整する小オフセットを使用|p8–11/§3.2。Δaの一律上下限なし。図だけから残差のσを作らない|
|A3 TanDEM-X TAFF, Ardaens & Fischer 2011|同上|実機航法評価、差分GPS|Δa推定の精度約0.20 m|p5/§4。航法誤差であり物理残差ではない。 commissioning時の自律制御権限にも制限|
|A4 PRISMA, D’Amico & Larsson 2011|Mango150 kg/Tango40 kg、750 km、近接編隊|形成維持、GNSS、ヒドラジン|Table 1の制御目標 aδa=0 m。航法のみ−0.5 m、誘導5 m|p4/Table 1の値は設定目標。実測σではない。aδaは無次元δaにaを掛けた長さ|
|A5 PRISMA, D’Amico et al. 2012|同上|SAFE実機実験、形成維持・再配置|ゼロ平均半長軸差を目標に維持し、再配置には意図的な差を利用|§III–IV。推定・制御・再配置を区別。A4と同じ実験系で独立標本に数えない|
|A6 GRACE, Spiridonova et al. 2014|精密大型二機、近円軌道、沿軌道列|2012年の維持マヌーバ前後、GNSS軌道解から平均要素計算|マヌーバ前7.28±0.20 m、後−9.66 m|pp10–12/§4/Table 3。推定・予測値。符号反転は意図した次の維持周期。形成失敗の残差ではない|
|A7 GRACE運用資料, CSR|同軌道二機、170–220 km範囲の例|運用中の形成維持|維持周期全体でΔaが平均的にゼロとなる運用|公式Configuration。平均ゼロ≠全時刻ゼロ。残差σ数値は提示なし|

A1の0.5 mは、今回の問いに最も近い具体的な形成後の値です。ただし二機の差と21機の標準偏差は異なり、冷ガス・GNSSを備えた機体の成果でもあります。A3は0.1 mと近い航法精度の桁を示しますが、制御精度の証明にはなりません。A6は、間隔に広い許容幅を持つミッションならm級の差を意図して運用できることを示します。

### A.2 推定精度と実残差

A6は、平均要素算出方法を評価する研究でもあります。Table 2の追加平均化後のΔa変動±0.01〜0.03 mは、伝播した軌道から計算した値の短周期変動を抑えた結果です。**真値に対する独立なGNSS誤差検証でも、実衛星が1 cm以内に制御された結果でもありません。** 本文のGRACE間隔を170〜270「meters」とする記述は、公式資料のkm単位と整合しないため、要求間隔の根拠には採用しません。

A4/Table 4の相対位置3D RMS約8.15 cm、相対速度0.24 mm/sも航法精度です。これを半長軸残差σへ換算していません。半長軸は位置だけでなく速度にも依存し、位置精度をそのままΔa精度とすることはできません。

### A.3 現在のσa=0.1 mへの判断

- **現実性：** 精密に形成した後の小残差という条件としては排除されない。ただし3U・推進機なし・今回の二値制御で取得できるという直接実績は未確認。
- **厳しさ：** CanXの二機差0.5 mより小さな桁を狙う、楽観的で厳しい形成条件に近い。
- **保守性：** 分散を大きくして制御必要性を作った値と断定する根拠はない。逆に、より大きな残差を見落とす可能性もある。
- **代表性：** 21機の形成後残差の分布、独立性、Gaussian性は未検証。
- **次の水準：** 0.01 mは楽観的下限、0.1 mは現研究との接続点、1 mは控えめな上限感度確認。名称nominalは「解析の中心」であり「実機の平均的性能」ではない。

ここで採用するのは**平均半長軸**です。瞬時の軌道要素や局所高度差には軌道周期の変動が含まれます。文献の平均化・フィルタと今回の平均軌道モデルを対応させ、単なる高度差を無断でΔaへ置き換えません。

## B. Ballistic coefficient literature

### B.1 まず何が分かったか

持続的な質量・形状・姿勢差が相対ドリフトを生むことは、実機運用から支持されます。しかし公開論文が推定した「実効抗力係数」には、大気モデル誤差や運用姿勢も含まれます。**その時系列のばらつきを、そのまま衛星固有の固定倍率の標準偏差にすることはできません。**

今回の定義は B*=Cd A/m [m²/kg]、基準B0*=0.0055 m²/kgです。文献がβ=m/(Cd A) [kg/m²]を使うならB*=1/β。二機差の厳密な換算は B*i/B*j−1=βj/βi−1 であり、大きな差では符号反転だけでは不十分です。TLEのBstarはSGP4用のパラメータで、この物理量と同一ではありません。

|資料|対象・方法|本文で確認した数値・現象|今回の固定個体差へ使える範囲|
|---|---|---|---|
|B1 Foster et al. 2015|Planet 3U、軌道決定・抗力推定|BCを数日以下の軌道弧から推定。姿勢により面積が変化|衛星別の校正が必要という実機の方法。固定σ=0.3%は報告なし|
|B2 Foster et al. 2017 / arXiv2018|12機Planet、510 km、飛行実績|幾何面積比約10に対し実効抗力比約3|pp9–10。LOW/HIGH運用・姿勢の違い。個体間σではない|
|B3 Brandt et al. 2020|Planet Flock 3Pの10機、TLEとMSISE-00|実効面積を2〜3日区間で調整。Table 1の中央値例697.77と938.62 cm²|§2.2/Table 1。大気・姿勢・質量仮定を吸収した推定値。約35%の差を恒常製造差と読まない|
|B4 Mason et al. 2017|MinXSS-1 3U、XACT姿勢制御|軸別指向精度0.0042〜0.0117°、3σ|Table 1/§2。特定の太陽指向機体の実測。抗力姿勢の恒常二機差へは直接換算不可|
|B5 Castro-Fernández et al. 2026|ANSER 3U二機、約575 km、実機差動抗力|LOWでも自然降下率が異なり、姿勢制御差が原因候補。LOW側姿勢を−6°に変更した運用|§5.2。持続する非対称抗力の実例。大きさのσは未確定|
|B6 Kahle et al. 2026|TSX/TDX、約505 km、質量・面積から評価|同じCd・面積の仮定で二機差初期<1%、2023–25年約3%|p3/§I.C。質量差に由来する長期差。大型衛星からの一般化|
|B7 Bowman 2002|HASDM用校正衛星、同型機体の比較|同型と判断したペアに<1%、<0.1%、0.2%以内などの一致例|p3/Table 2・本文。モデルを用いた校正。3Uの実測分布ではない|
|B8 Dey et al. 2025, ESS|GRACE形状・空力モデル|姿勢・表面相互作用・組成によりCdと面積の両方が変化|§5。物理機構。モデル不確かさを個体差へ変換しない|
|B9 Dey et al. 2025, Space Weather|嵐時の空力計算|2024年5月例でCdの時間変化約8〜12%|§3/Fig.4。嵐時のモデル値。同型衛星が8〜12%異なる証拠ではない|
|B10 van der Zwaard et al. 2025|3P PocketQube、差動抗力による空力同定提案|異なる板姿勢の二機を利用してGSIを推定する方法|§2–3。3Pは3Uではない。シミュレーションであり固定個体差の飛行測定ではない|
|B11 Cal Poly CDS Rev14.1, 2022|CubeSat規格|3Uの典型的最大質量6 kg|§2.2.10/Table 1。許容最大質量であり同型機の質量マッチング公差ではない|

直接の3U資料はB1〜B5です。それらは姿勢と運用を無視した単一の精密BC値が不十分なことを示しますが、0.3%や1%の母集団分布を与えません。B6/B7は小さな恒常差が物理的に不自然ではないという補助根拠で、CubeSatへの転用の強さは下がります。

### B.2 0.3%と1%の物理的意味

小さな変化では、

    δB*/B* ≈ δCd/Cd + δA/A − δm/m

となります。これは、質量だけ違えば重い衛星ほど同じ面積に対する抗力加速度が小さく、面積やCdが大きければ大きくなる、という意味です。各項の誤差が共通か衛星固有かを分けて扱います。

4 kgで質量差だけを考えると、0.1/0.3/1%は約4/12/40 gに相当します。**これは実機の公差データではなく一次換算**です。最終搭載品の構成・配線・重心調整・質量マッチングによって、この範囲より小さくも大きくもなり得ます。CDSにはその集団公差の保証はありません。

今回の付加物なし直方体を仮定すると、LOW面を流れへ向けた状態から一軸にθだけ傾いた幾何面積は

    A(θ)=0.01 cosθ + 0.03 |sinθ| [m²]
    (A(θ)−A(0))/A(0) ≈ 3 |θ|     （θはrad）

です。したがって面積のみなら0.3%は約0.057°、1%は約0.191°の傾きで生じます。小さな角度でも長い側面が見え始めるためです。ただしCd一定・付加物なしの幾何換算であり、実機の姿勢誤差分布を測った結果ではありません。HIGH側に同じ倍率が掛かるとも限りません。

### B.3 60日一定とできるもの

|要因|60日間の扱い|現在の固定倍率モデルとの整合|
|---|---|---|
|推進剤消費のない機体の質量差|ほぼ恒常。相対質量変化が小さいことを確認|LOW/HIGHの双方へ同じ倍率を掛ける物理的根拠が最も明確|
|恒常形状・付加物の個体差|固定形状なら持続|姿勢によって面積差が変わるためLOW/HIGH共通倍率は別途検証|
|系統的な姿勢ずれ・校正差|持続し得るが運用モード・日照で変化|姿勢モード別の平均面積差として扱う方が適切|
|表面・組成依存Cd|表面差は持続し得るが入射流・温度・組成が変わる|60日完全固定のCd差は未検証|
|ゼロ平均の姿勢角ジッタ|短周期でも平均面積増加が残る場合あり|角度平均ゼロから抗力差平均ゼロを結論しない。非線形な投影面積のため|
|大気モデルの共通バイアス|共通尺度を変える|衛星固有固定倍率のσとしては使わない|

現モデルの0.3%固定差は、**質量・恒常ハードウェア差を想定した小さな工学的感度条件**に近いものです。representativeとは認定できず、現実性の裏付けが弱いという意味ではunsupportedな統計仮定を含みます。過大・aggressiveと断定する証拠も、上限を覆うconservativeという証拠もありません。1%も物理的にあり得る感度条件ですが、これを実衛星の標準偏差と呼びません。3〜5%は今回の直接証拠から通常3U群のnominalとして採用しません。

### B.4 支配要因の判断

**一般に質量・面積・Cd・姿勢のどれが最大かは、今回の公開資料だけでは決められません。** 展開パネルを持つANSER/Planetでは運用姿勢と展開状態が大きく効きます。精密姿勢制御・同じ形状なら、測定できる質量差が単純な固定倍率として最も扱いやすい要因です。同じ大気・材料に対するCdの大きな共通不確かさは、二機間のCd差が大きいことを意味しません。

## C. 3U CubeSat-specific evidence

### C.1 3Uの直接資料と限界

ANSERは今回に最も近い、3U・二値に近い差動抗力・100 km級の運用例です。ただし機体は大きな展開翼を持ち、今回の0.01/0.03 m²の単純な機体と同一ではありません。Table 1の質量は両機3.357 kg、翼面積は1800/2100 cm²。FLW1の一部展開不良が原因で、通常の製造ばらつきの標準偏差ではありません。比は約14.3%違いますが、LOW時に同じ差が掛かるとは限りません。本文§5.3の係数も固定参照面積0.21 m²とNRLMSISE-00を使う推定量なので、係数だけを今回のCd=2.2と比較しません。[B5原論文](https://doi.org/10.1016/j.actaastro.2026.04.057)

Planetの実効抗力推定、MinXSSの指向実績、CDSの機体規格は、それぞれ異なる種類の根拠です。指向精度が良い一機の実績から全21機の恒常BC差を作ることも、最大質量の規格から実機質量の標準偏差を作ることもできません。

### C.2 次に実機情報で埋める項目

1. 同型機の最終実測質量一覧、質量マッチングの有無、燃料・展開機器。
2. LOW/HIGH別の実測平均姿勢ずれ・ジッタ分布と実効投影面積。
3. 同時期・近い環境で校正した衛星別の抗力倍率。大気モデル誤差との共分散も記録。

これらがないうちは、今回の3水準は**分散の可能性を検討する入力**です。実機データの代替として扱わないでください。

## D. Spacing / formation requirement literature

### D.1 まず何が分かったか

実ミッションの要求には、近接して精密な観測基線を維持するものと、数十〜数百kmの観測間隔を大きな運用範囲で維持するものがあります。**同じ「formation」という名前でも、必要精度は桁違いです。** 100 km列の最大誤差に、SARの航法精度やレーザー測距精度をそのまま適用してはいけません。

|ミッション・原資料|目標配置|確認した要求・実績|理由・種類|100 km列への対応|
|---|---|---|---|---|
|GRACE-FO, JPL 2018|沿軌道220 km|170〜270 km維持、約5か月毎のstation keeping計画|重力観測、運用維持範囲。測距精度とは別|中心220 kmに対する±50 km。100 m級の維持要求の根拠にならない|
|ANSER, B5/§5.2|沿軌道100 km級|100≤S≤200 kmのstation keeping範囲|3U実機の運用閾値・編隊保持|今回に近い距離だが許容範囲は広い。100 km±100 mではない|
|CloudSat–CALIPSO, Richon & Case 2005|通過時間差12.5秒、約705 km|±2.5秒の計画維持範囲|観測協調・安全、打上げ前の運用設計|軌道速度換算で約94 km±19 km。地上投影距離とは別|
|TanDEM-X, A2/Table 4|近接SARヘリックス、数百m〜km|沿軌道制御要求200 m RMS、実績27.3 m RMS|SAR基線・編隊制御|RMS要求。今回の最大誤差500 mとは異なる基準|
|TanDEM-X, B6/§I.B|同上|36時間先まで沿軌道位置のknowledge <200 m|同期アンテナ選択・レーダ照射回避、予測精度|物理的な間隔誤差200 mの一律上限と混同しない|
|CanX, A1/Table 5|50〜1000 mのATO/PCO|実績3D位置誤差RMSの列に約0.345〜0.602 m|近接編隊技術実証、達成性能|要求でなく達成値。100 km列の必要精度には転用不可|
|PRISMA, A4/Table 2|数百m〜数km|形成D/E/Fの沿軌道最大誤差5.2/3.9/3.9 m|精密近接制御の実績|要求でなく達成値。GNSS航法精度とも別|
|Planet, B2/運用結果|510 km、12機の沿軌道配置|例の角間隔3.6°。均等配置・相対速度停止を目標|地上観測の配置規則性。公開数値の硬い維持deadbandは未確認|目標間隔と許容誤差を分ける。3.6°は約433 kmの換算値|

GRACE-FO/ANSERは大きな維持範囲でも運用可能な実例です。TanDEM-X/CanX/PRISMAは精密制御の実例ですが、今回と要求理由が異なります。今回の21機100 km列の用途が決まらない限り、100 m・500 m・1 kmのどれが「問題」かは文献だけでは一意に決まりません。

### D.2 時間要求からの換算

円軌道では近似的に e≈v Δt です。つまり同じ地点を通る観測時刻差の誤差を距離に直すには、衛星の軌道速度を掛けます。705 kmの円軌道で約7.50 km/sなら、CloudSatの12.5±2.5秒は約93.8±18.8 kmです。元論文p2の画像を確認し、OCRの「12.5 a.5」を±0.5秒と読み替えていません。

今回の400 kmで、100/500/1000 mは概ね13/65/130 msの通過時刻差に相当します。これだけの同期が観測に必要かは、計測方式・時間同期・地上処理と併せて決める必要があります。要求根拠がないまま数値だけ厳しくしません。

|最大誤差|100 km比|今回の役割|
|---:|---:|---|
|10 m|0.01%|非常に精密な側の参考値。今回の主基準にはしない|
|100 m|0.1%|厳しい仮想要求|
|500 m|0.5%|暫定的な主要比較基準|
|1 km|1%|緩い仮想要求|
|5 km|5%|広い配置要求の追加参考。自動的に失敗とは呼ばない|

500 mを中央に置くのは、既存結果を同じ尺度で比較し、100 mと1 kmの間の判断を見せるためです。**文献が500 mを推奨したという意味ではありません。** GRACE-FOやANSERに近い広い維持範囲を想定するなら、1 kmより緩い基準も必要です。逆に干渉計の用途なら100 mより厳しくなる場合があります。

## E. Evidence quality

### E.1 根拠の強さは何に対するものか

- **Strong evidence：** 本文で直接確認した当該ミッションの飛行値・公式要求。例えばCanXの0.5 m、GRACE-FOの維持範囲。
- **Moderate evidence：** 近い機構・近いサイズの実例から、影響の存在や入力の桁を支持する。ただし今回の機体・条件への転用に限界がある。
- **Weak evidence：** 大型機・別形状からの一般化、単一事例、推定モデルとの分離が不十分。
- **Engineering assumption：** 公開された分布がないため、感度確認として決めた入力値や仮想要求。

原資料の0.5 mがStrongでも、「21機のσaは0.1 m」はStrongになりません。同様に同型大型機の3%という根拠がStrongでも、「3U群のσBは0.3%」はEngineering assumptionです。

### E.2 採用しなかった読み替え

|読み替え|採用しない理由|
|---|---|
|展開速度・投入半長軸誤差→形成後の残差|形成操作前後が異なる|
|GNSS相対位置精度→平均半長軸残差|位置と速度・平均化が必要|
|平均要素計算のcm級変動→cm級制御成功|計算方法の安定性と実制御は異なる|
|Cdモデルの数%〜数十%不確かさ→衛星間に同じσ|共通誤差と差動誤差を混同|
|軌道推定した実効面積の分散→製造時の恒常面積σ|密度・姿勢・運用誤差を吸収|
|最大質量6 kg→4 kg機の質量公差|制約と製造分布は異なる|
|3P PocketQube→3U CubeSatの実証|規格も実験種別も異なる|
|SARの200 m RMS→全時刻最大誤差200 m|統計量が異なる|
|測距のnm/µm級精度→編隊を同じ精度で制御|測定精度と維持要求は異なる|

### E.3 文献表

全て本文または公式資料の該当節を確認した資料です。参考文献の件数は独立した飛行標本数を意味しません。著者が多いものは第一著者et al.を表記します。

|ID|題名・著者・年・掲載先|DOI / 原文URL|対象・採用箇所|数値の種類|
|---|---|---|---|---|
|A1|CanX–4 and CanX–5 Precision Formation Flight: Mission Accomplished! — Grant Bonin, Niels Roth, Scott Armitage, Josh Newman, Ben Risi, Robert E. Zee; 2015, SmallSat SSC15-I-4|[原文](https://www.utias-sfl.net/wp-content/uploads/SSC15-I-4-CanX_4andCanX_5PrecisionFormationFlight_MissionAccomplished.pdf)|CanX、pp9–10/Table3、p13/Table5|形成後推定差0.5 m、実測制御性能|
|A2|TerraSAR-X / TanDEM-X Formation Control – First Results from Commissioning and Routine Operations — Ralph Kahle, Benjamin Schlepp, Michael Kirschner; 2011, 22nd ISSFD|[原文](https://issfd.org/ISSFD_2011/S4-Flight.Dynamics.Operations.1-FDOP1/S4_P1_ISSFD22_PF_037.pdf)|TSX/TDX、§3.2/Table4|制御要求200 m RMSと実績27.3 m RMS|
|A3|TanDEM-X Autonomous Formation Flying System: Flight Results — Jean-Sébastien Ardaens, Daniel Fischer; 2011, IFAC World Congress|[原文](https://skoge.folk.ntnu.no/prost/proceedings/ifac11-proceedings/data/html/papers/2374.pdf)|TDX、p5/§4|Δa航法精度0.2 m|
|A4|Navigation and Control of the PRISMA Formation: In-Orbit Experience — Simone D’Amico, Robin Larsson; 2011, IFAC World Congress|[原文](https://skoge.folk.ntnu.no/prost/proceedings/ifac11-proceedings/data/html/papers/2494.pdf)|PRISMA、Tables1/2/4|目標、達成誤差、航法誤差を別々に採用|
|A5|Spaceborne Autonomous Formation-Flying Experiment on the PRISMA Mission — Simone D’Amico, Jean-Sébastien Ardaens, Robin Larsson; 2012, JGCD35(3),834–850|[DOI](https://doi.org/10.2514/1.55638)、[著者公開本文](https://www.researchgate.net/publication/225024880_Spaceborne_Autonomous_Formation_Flying_Experiment_on_the_PRISMA_Mission)|PRISMA、§III–IV|実機実験、再配置と維持を区別|
|A6|Precise Mean Orbital Elements Determination for LEO Monitoring and Maintenance — Sofya Spiridonova, Michael Kirschner, Urs Hugentobler; 2014, 24th ISSFD|[DLR原文](https://elib.dlr.de/103814/1/Spiridonova_ISSFD_2014_upd.pdf)、[本文ミラー](https://www.researchgate.net/publication/281638840_Precise_mean_orbital_elements_determination_for_LEO_monitoring_and_maintenance)|GRACE、§4/Tables2–3|運用軌道からの推定・伝播・平均化による変動|
|A7|GRACE Operations: Orbital Configuration — University of Texas CSR, 公開年記載なし|[公式](https://www2.csr.utexas.edu/grace/operations/configuration.html)|GRACE、Configuration本文・図|運用範囲、半長軸差の平均ゼロ|
|B1|Orbit Determination and Differential-Drag Control of Planet Labs Cubesat Constellations — Cyrus Foster, Henry Hallam, James Mason; 2015, AAS15-524|[原文](https://arxiv.org/pdf/1509.03270)|Planet3U、軌道決定・BC推定節|実機運用・推定方法。固定σ数値なし|
|B2|Differential Drag Control Scheme for Large Constellation of Planet Satellites and On-Orbit Results — Cyrus Foster, James Mason, Vivek Vittaldev, Lawrence Leung, Vincent Beukelaers, Leon Stepan, Rob Zimmerman; 2017, IWSCFF17-33; 公開稿2018|[原文](https://arxiv.org/pdf/1806.01218)|Planet3U、pp9–12|実効HIGH/LOW比約3、配置角3.6°。二機恒常差ではない|
|B3|A Simple Method for Correcting Empirical Model Densities During Geomagnetic Storms Using Satellite Orbit Data — Daniel A. Brandt, Charles D. Bussy-Virat, Aaron J. Ridley; 2020, Space Weather18,e2020SW002565|[本文](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2020SW002565)|Planet3U、§2.2/Table1|校正された実効面積。純粋な面積測定ではない|
|B4|MinXSS-1 CubeSat On-Orbit Pointing and Power Performance: The First Flight of the Blue Canyon Technologies XACT 3-axis Attitude Determination and Control System — James Paul Mason et al.; 2017, Journal of Small Satellites6(3),651–662|[著者稿](https://arxiv.org/pdf/1706.06967)|3U、Table1/§2|3σ指向精度の飛行実測|
|B5|Formation Flying of the ANSER Cluster Using Differential-Drag Control — Iván Castro-Fernández et al.; 2026, Acta Astronautica247,77–91|[DOI](https://doi.org/10.1016/j.actaastro.2026.04.057)、[著者公開本文](https://www.researchgate.net/publication/404260778_Formation_flying_of_the_ANSER_cluster_using_differential-drag_control)|3U、Table1/§5.2–5.3|質量・形状、運用範囲、姿勢差の解釈、実効Cd推定|
|B6|Operational Experience of Close Formation Flight at 500 km Altitude during Severe Geomagnetic Storms — Ralph Kahle, Benjamin Schlepp, Simon Gruber, Lars Sprengelmeyer, Martin Wermuth, Stefan Hackel, Florian Meissner, Frederic Schoutetens; 2026, 30th ISSFD|[DLR原文](https://elib.dlr.de/225026/1/Kahle_et_al_ISSFD2026.pdf)|TSX/TDX、Table1/§I.B–C|質量に由来するBC差、位置予測要求|
|B7|True Satellite Ballistic Coefficient Determination for HASDM — Bruce R. Bowman; 2002, AIAA2002-4887|[原文](https://sol.spacenvironment.net/JB2008/pubs/JB2006_AIAA_2002_4887.pdf)|同型校正物体、p3/Table2|モデル校正・形状同等ペアの比較|
|B8|Investigation of Gas-Surface Interactions and Neutral Atmospheric Properties and Their Impact on Satellite Drag Coefficients — Soumyajit Dey, Phillip C. Anderson, Aaron L. Bukowski, Marcin D. Pilinski; 2025, Earth and Space Science12,e2024EA003998|[本文](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2024EA003998)|GRACE形状等、§5|物理空力モデル、個体差測定ではない|
|B9|An Extreme Events Study of the Storm-Time Variations in Satellite Drag Coefficients — Soumyajit Dey et al.; 2025, Space Weather|[本文](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2025SW004634)|嵐時空力、§3/Fig4|時間変化の計算値|
|B10|Using Differential Satellite Drag as an Experiment to Characterise Gas–Surface Interactions — Rens van der Zwaard, Stefano Speretta, Christian Siemes, Eberhard Gill; 2025, CEAS Space Journal|[本文](https://link.springer.com/article/10.1007/s12567-025-00657-0)|3P PocketQube、§2–3|提案とシミュレーション|
|B11|CubeSat Design Specification Rev14.1 — Cal Poly; 2022-02-09, 公式規格|[原文](https://www.nasa.gov/wp-content/uploads/2018/01/cubesatdesignspecificationrev14_12022-02-09.pdf)|3U、§2.2.10/Table1|最大質量仕様。製造標準偏差なし|
|C1|GRACE-FO Launch Press Kit: Launch Events and Mission Phases — NASA/JPL; 2018|[公式](https://www.jpl.nasa.gov/news/press_kits/grace-fo/mission/launch-events/)|GRACE-FO、Mission Phases/Science Phase|計画維持範囲220±50 km|
|C2|An Overview of the Earth Science Afternoon Constellation Contingency Procedures — Karen V. Richon, Warren Case; 2005, AAS05-369|[NASA原文](https://ntrs.nasa.gov/api/citations/20050210076/downloads/20050210076.pdf)|A-Train、p2/Concept of the Control Box|打上げ前運用設計、時間間隔12.5±2.5秒|

## F. Recommended simulation matrix

### F.1 数値は何を指定するか

今回は新規runを作らず、次の入力候補を `recommended_conditions.json` に保存しました。既存設定を上書きしません。

|水準|σa [m]|σB/B|意味|
|---|---:|---:|---|
|optimistic|0.01|0.1%|極めて良い形成・個体差校正を仮定する下限|
|nominal sensitivity|0.1|0.3%|既存研究の中心。実機集団の典型値とは呼ばない|
|conservative sensitivity|1|1%|形成残差・恒常差が一桁程度大きい場合の感度確認。確率上限ではない|

**3水準の同時組合せ9条件を推奨しているわけではありません。** 当面は片方だけ変える既存の方法を維持し、観測ミッションと実機根拠を固めた後に必要な条件だけ実行します。0.01/0.1/1 mは平均半長軸の衛星別Gaussian生成時のσ、0.1/0.3/1%は衛星固有倍率の生成時のσです。

### F.2 二機差と衛星別σの換算

独立な衛星別Gaussianを生成し全体平均を引く現在の方式では、隣接差の標準偏差は

    std(δa[j+1]−δa[j])=√2 σa
    std(ε[j+1]−ε[j])=√2 σB

です。平均を引いても差は変わりません。つまり衛星別0.1 mは典型的な二機差約0.141 m、衛星別0.3%は二機差約0.424%に相当します。**文献の二機差0.3%と今回のσB=0.3%は同じ入力ではありません。** 二機差の標準偏差0.3%を仮に目標とするなら、独立生成σは約0.212%です。

平均補正だけでは21機の実現標準偏差は指定σと完全一致しません。現在の代表Aは指定0.1 mに対し実現約0.111 m、Bは指定0.3%に対し約0.225%です。次の報告でも指定σ・実現σ・最大隣接差を別々に保存してください。同一衛星群に共通する校正誤差や空間相関があるなら、独立Gaussianの仮定も再検討します。

### F.3 理論的な桁の確認

400 kmの基準円軌道では

    Δn≈−(3/2)n Δa/a
    de/dt≈−(3/2)n Δa

となります。半長軸が少し大きい衛星は平均的に遅く進むため、相対配置は徐々に変わります。二機差0.1 mなら、抗力差を無視したドリフトの大きさは約14.7 m/日です。衛星別σa=0.1 mなら隣接差の典型値√2σに対して約20.7 m/日です。20gapの最大値はさらに大きくなり得るため、単一二機の理論値を既存の全列最大誤差419 m/7日と直接同一視しません。

これは新たな物理計算ではなく、既存の基準半径・重力定数を使った一次近似の換算です。`audit.py` が計算し、JSONへ保存します。

### F.4 現在の制御結果をどう評価するか

既存保存結果では、Aは60日無制御3.774 km→制御317.8 m、Bは無制御で40.94日頃に順序逆転→制御235.7 mでした。この改善は**入力条件を仮定したモデル上の成果**として維持できます。入力値が実機の代表値であるという新たな証明にはなりません。

|評価条件|A：317.8 m|B：235.7 m|解釈|
|---|---|---|---|
|仮想最大誤差100 m|不適合|不適合|全時刻の厳しい要求は満たしていない|
|仮想最大誤差500 m|60日内で適合|60日内で適合|暫定比較尺度では有効|
|仮想最大誤差1 km|60日内で適合|60日内で適合|さらに緩い要求なら精度面では十分|
|実ミッション要求|未判定|未判定|想定用途・寿命・姿勢変更許容・高度予算が未決定|

TanDEM-Xの200 m RMSに、今回の最大値318 mを比較して不適合とは判定しません。A/Bの追加高度損失約276/367 mと数千回の姿勢変更も評価が必要です。100 kmに対する誤差割合が小さくても観測が要求すれば意味がありますが、要求が広ければこの負担を支払う理由は弱くなります。

現時点の判断は、ユーザーのCase 1〜3では**「Case 1を条件付きで扱う」**です。物理的にあり得る入力として研究は続けられますが、実機分布として確認済みとはしません。Case 2（大きすぎるので下げる）やCase 3（非現実的なので除外）を裏付ける証拠は得られませんでした。特にBの固定σについて、文献から実機代表性を確定することはできません。

## G. Unresolved points / 必須10問への回答

|問い|今回の回答|
|---|---|
|Q1：形成後0.1 mは現実的か|厳しく形成した小残差の条件としてはあり得る。3Uの通常残差という直接根拠はない|
|Q2：0.01/0.1/1 mの桁|0.1〜1 mを中心とする感度範囲が合理的。0.01 mは楽観的な下限。集団σは工学的仮定|
|Q3：3UのBC差は何%か|0.1〜1%の恒常差は質量・姿勢の換算から物理的にあり得る。直接確認した3Uの固定差分布はなく、通常何%とは断定できない|
|Q4：0.3%を60日固定|質量由来なら合理的。LOW/HIGH共通倍率の統計分布や姿勢由来の60日一定性は未検証|
|Q5：1%は現実的か|物理的にはあり得る。3Uの代表σとする証拠はない。上限感度側として採用|
|Q6：mass/area/Cd/attitudeの支配性|機体依存。3U運用例は姿勢・展開状態が重要。共通Cd不確かさを個体差と混同しない|
|Q7：100 km列の許容精度|共通要求はない。今回確認した実例は精密近接編隊のm級性能から100 km級列の数十km維持幅まで|
|Q8：100/500/1000 mのどれ|500 mを暫定比較の中心、100 m厳しい側、1 km緩い側。全て仮想要求。文献認定の主要求とはしない|
|Q9：318/236 mは良いか|500 m/1 km尺度なら60日で有効、100 m尺度なら不足。実要求・高度・切替予算の合格は未判定|
|Q10：次の3水準|σa=0.01/0.1/1 m、σB/B=0.1/0.3/1%、許容最大誤差100/500/1000 m。根拠の強さは各表のとおり|

未解決なのは「制御で誤差を抑えられるか」よりも、**どの精度を何のために維持し、形成後にどの状態を達成でき、衛星ごとの平均抗力をどれだけ校正できるか**です。追加ゲイン探索の前にこれを決める必要があります。今回は初期位置誤差・投入誤差・新しい外乱・制御パラメータ変更を導入しませんでした。

## H. 再現・非破壊確認・Git

- 新規成果：この報告、根拠台帳 `evidence.json`、候補条件 `recommended_conditions.json`、換算確認 `audit.py`、結果 `audit_result.json`。
- PDFと抽出本文はGit除外の `outputs/formation_input_literature/` へキャッシュ。公開論文の全文をGitへ転載しない。
- 既存の物理・制御・積分・設定・run・回帰referenceは変更していない。新規シミュレーションなし。
- 換算と台帳整合性の確認：`.analysis-venv/Scripts/python.exe -X utf8 docs/formation_input_evidence/audit.py`。
- 既存テスト・Git情報は `VERIFICATION.md` に保存する。
- 本報告の3水準は、実機母集団の確率分布を確定したものではない。生成seedを改めて選んで有利な結果を作ることも推奨していない。
