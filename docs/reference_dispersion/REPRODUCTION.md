# Reference conditionの再現方法

0.5 m・0.5%は実機母集団から推定された統計的nominal値ではなく、文献で確認された現実的な桁を踏まえた工学的Reference sensitivity conditionである。

## 固定する入力

- `data/profiles.json`の21要素を以後再利用する。乱数を再生成しない。
- 平均を引いた旧A seed 4、旧B生成seed 10002を、それぞれsample標準偏差 **ddof=1** で0.5 m・0.005へ規格化した。
- A/ABの初期半長軸offsetは完全一致。B/ABの固定倍率は完全一致。ON/OFFでも一致。
- 固定倍率はLOW/HIGH/slew中すべての`Cd*A/m`へ乗ずる。質量差等を代表する簡略化であり、姿勢誤差の精密モデルではない。
- active環境のepoch・NRLMSIS 2.0・歴史的宇宙天気・軌道幾何・120 s密度キャッシュ方式を継承する。
- 新しい軌道では衛星位置が異なるため、旧runの密度時系列をそのまま複製しない。既存の宇宙天気CSVとキャッシュ方式を再利用し、必要な位置で密度を計算する。
- 制御・積分coreは変更しない。RK4最大2 s、保存10 s。人工外乱はOFF。

## 実行順序

リポジトリrootで実行する。既存runがあれば後処理だけを行い、物理計算を繰り返さない。

```powershell
$env:PYTHONUTF8='1'
.analysis-venv/Scripts/python.exe scripts/preserve_reference_inputs.py
.analysis-venv/Scripts/python.exe scripts/run_reference_dispersion.py --days 30
.analysis-venv/Scripts/python.exe scripts/reference_profile_metadata.py
.analysis-venv/Scripts/python.exe scripts/analyze_reference_dispersion.py --days 30
.analysis-venv/Scripts/python.exe scripts/supplement_reference_results.py --days 30
.analysis-venv/Scripts/python.exe scripts/check_reference_proportionality.py
.analysis-venv/Scripts/python.exe scripts/audit_reference_density.py --days 30
```

30日結果の数値正常性・各ONケースの正のgap・制御の継続性・精度・コストを確認し、`data/extension_decision.json`へ判断理由を記入して`reviewed=true`とする。安全・有用性を確認した場合だけ`extend_to_60d=true`にする。

```powershell
.analysis-venv/Scripts/python.exe scripts/run_reference_dispersion.py --days 60 --prefix outputs/reference_dispersion_20261012/batch_30d
.analysis-venv/Scripts/python.exe scripts/analyze_reference_dispersion.py --days 60
.analysis-venv/Scripts/python.exe scripts/supplement_reference_results.py --days 60
.analysis-venv/Scripts/python.exe scripts/audit_reference_density.py --days 60
.analysis-venv/Scripts/python.exe scripts/plot_reference_dispersion.py
.analysis-venv/Scripts/python.exe scripts/check_reference_regression.py
.analysis-venv/Scripts/python.exe -m pytest -q
.analysis-venv/Scripts/python.exe scripts/preserve_reference_inputs.py --check
```

出力先は設定ファイル`configs/reference_dispersion.yaml`に固定する。実行器は既存run directoryを上書きしない。再現計算が必要なら新しいoutput先を設定し、完成した旧出力を保持する。後処理の主要出力directoryも既存なら拒否する。

## 評価量の読み方

主要誤差は従来と同じ固定基準半径の定義：`e_i=a0*(lambda_(i+1)-lambda_i)-100000`。現在半径の物理弧長誤差とは区別する。初期半長軸offsetを与えても初期mean longitudeは変更せず、e(0)=0である。

RMSは全保存時刻・全20gapの二乗誤差を平均して平方根を取る。区間RMSも同じ定義。100/500/1000 mはstrict/reference/looseの**仮想評価閾値**であり、実ミッション要求ではない。超過は厳密に`>閾値`。時間率は左端の保存サンプルを次のサンプルまで保持した近似で計算する。

遷移数はslew開始を数える。command、開始、完了をイベントCSVで別々に追跡する。command delay=0なのでcommandと開始は同時刻。完了L-H-Lは各衛星のLOW帰還完了まで観測された1往復。区間cycle数はLOW帰還完了時刻で割り当てる。

HIGH-equivalent timeは`integral((A(t)-A_LOW)/(A_HIGH-A_LOW) dt)`を全衛星で合計した量。姿勢変更中の面積増分も含む。stable HIGH timeは姿勢変更中を含まない。

追加平均高度損失は、**同じprofile・同じ大気場の全LOW**を基準として`mean(a_OFF-a_ON)`。正ならONが余計に高度を失ったことを意味する。初期高度からの通常の全体低下とは異なる。

順序逆転は固定番号の隣接signed gapが初めて0以下になった保存時刻。分解能は10 s。OFFが逆転した場合、ON/OFF精度改善の主比較は、その直前の共通時間窓で行う。逆転後OFFの誤差は診断値だけで、通常のformation accuracyではない。逆転後OFF半長軸を用いる高度コストは、衝突や衛星間相互作用を含まないモデル上の反実仮想基準である。

最小物理距離は同一軌道面の全衛星対の直線距離を保存時刻で計算する。これだけで衝突安全性を保証するものではない。

残留相対速度はLOW帰還後のgap速度。P/D分解・符号付きgap状態は既存残留評価を使用し、大きな配列はローカルNPZへ保存する。一定Dの再横断予測は仮定を伴う外挿であり、持続B差や継続大気変動でDが変わる場合は有効性フラグに従う。

## 保存先

- `data/profiles.json` / `reference_profiles.csv`：今後使う実数値21個。
- `data/30d/`, `data/60d/`：各ON/OFF metrics、共通時間窓比較、区間評価、イベント、cycle、代表q診断。
- 幅広い全cycle表・全HIGH trigger表はローカル解析用。Gitには明示的に抽出した小さなcycle sample表と集計・イベントを保存する。sample表だけから総数を再計算しない。
- `outputs/reference_dispersion_20261012/`：大きな状態履歴、checkpoint、残留状態配列、描画用cache。Gitへ追加しない。
- `figures/`：日本語の研究図6枚。
- `REPORT.md`：研究判断用サマリーと詳細技術付録。

Python環境は既存の依存関係と`pymsis==0.12.0`を使用する。新しい依存追加はない。
