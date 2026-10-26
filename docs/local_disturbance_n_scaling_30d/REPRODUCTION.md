# 再現手順

元の14日runがローカルに保存されている状態で、リポジトリルートから実行する。
なければ先に元の14日実験を再現する。既存runがある場合は再計算しない。

```powershell
$env:PYTHONUTF8='1'
.analysis-venv/Scripts/python.exe -u scripts/extend_local_disturbance_n_scaling.py --count 21
.analysis-venv/Scripts/python.exe -u scripts/extend_local_disturbance_n_scaling.py --count 301
.analysis-venv/Scripts/python.exe -m pytest -q > outputs/local_disturbance_n_scaling_30d_20261026/tests.txt
.analysis-venv/Scripts/python.exe scripts/analyze_local_disturbance_n_scaling_30d.py
```

延長runnerは、完全一致する完了suffixがあれば再利用する。入力不一致・不完全な既存suffixは保持して停止する。上書きしない。
公開先のdata/、figures/が既に存在する場合も停止する。別の再解析には新しい出力・公開ディレクトリを明示する。
既存のテストログを保持したい場合は新しいログ名を使い、解析スクリプトの読込先も合わせる。

## 保存内容

- 新規Nxxx直下のtime/a/lambda/rho/area/mode.npyは**Day14〜30だけ**。元のprefixを複製しない。
- checkpoint.jsonのイベントは累積、時刻は全て開始epochからの秒。
- input_checks.jsonに元checkpoint SHA256と唯一の設定変更（duration）を保存。
- analysis/time_s.npy、raw_gap_error_m.npy、delta_gap_error_m.npyは**Day0〜30**。
  元の状態全体を複製せず、prefix/suffixをchunk単位で読み、誤差だけを結合する。
- join_checks.jsonは接続境界の完全一致と重複点除外の記録。
- metrics.jsonは共通評価、extended_metrics.jsonは区間別・14日対比・イベント数等も含む。
- horizon_comparison.csvは14/21/30日比較、windows.csvは時刻区間比較。
- gaps.csv、spatial.csv、disturbed_events.csv、history_60s.csvも保存。
- Gitへは小さいCSV/JSON、設定、図、文書、コード、テストのみ。NPYは追加しない。
- 図は60 s間隔、公開用時系列CSVは600 s間隔。全指標は10 sの元データで計算。

物理・制御・積分coreが不変のため、24 h software referenceは前タスクで新規実行した3ケースを再利用し、現行評価コードで再照合する。
新しい継続テストでは、slew途中で分割した結果が一括積分と全保存状態・全イベントで厳密一致することを確認する。

branch：`astra/local-disturbance-n-scaling-30d`。mainへmergeしない。
