# Git・検証・実行記録

## Branchとcommit

- Branch：`astra/initial-dispersion-screening`。
- 研究の引継ぎ元：`b727c440ebf10d5f3bcb77a55060750d2d263154`。
- コード・設定・11新規テスト：`ea54359ffe6ac1b6197aeb212ad6863197ce7a1c` — `feat: add passive initial-a and ballistic-dispersion screening`。
- 成果commit：`docs: report dispersion screening results and regression evidence`。このファイルを含むcommitのhashは `git log -1 --format=%H` で取得できます。
- mainへのmergeやforce pushは行いません。

## 実行した主要操作

1. 既存active60日設定、保存結果、歴史的F10.7/Ap CSVを確認。旧ファイルのSHA256を作業前に保存。
2. `pytest -q tests/test_passive_dispersion.py`：11件合格。
3. `pytest -q`：既存245件＋新規11件＝256件合格、43.13 s。
4. `check_baseline_regression.py --existing-suite outputs/relative_motion_20261011/regression_24h/suite --output outputs/initial_dispersion_20261011/regression_24h`：既存24時間3ケースを凍結referenceと照合、全て合格。物理baselineの再計算なし。
5. `run_dispersion_screening.py --stage initial`：独立37ケースを0–7日。保存時系列の先頭1日と7日で評価。
6. 7日結果を確認し、`--stage extend`：半長軸0.1 mの5乱数＋B0.3%の代表乱数2を、保存7日状態から30日へ継続。
7. `--stage combination`：両要因の影響が大きいため、半長軸1 m＋B1%の独立パターン1組だけ0–7日。
8. `--stage analyze`、`analyze_dispersion_details.py`、`plot_dispersion_screening.py`、`verify_dispersion_screening.py`、`report_dispersion_screening.py`：集計・理論・密度照合・日本語5図・報告書。
9. 解析後に旧ファイルのSHA256を再照合：8,496ファイル、変更・欠損0。
10. 小さい集計・図・文書だけGitへ追加。旧ファイルの変更なし。巨大NPYと実行ログはGit対象外。

独立した延長と組合せは別プロセスで実行しました。モデル・外乱・初期配置を変更する並行作業は行っていません。全新規ケースは全LOW・制御OFFです。

## 計算量

- 最初の37ケース：7日をまとめて計算、約1,537 s。
- 6ケースの継続：7–30日、約1,379 s。0–7日の再計算なし。
- 組合せ1ケース：7日、約211 s。
- 既存のNRLMSIS-only0分散、60日データ、24時間回帰データは再利用。
- 50 m・追加乱数・60日・制御ケース・パラメータ再探索は省略。

上記は今回の環境での経過時間であり、他環境の性能を保証しません。後処理の密度監査はキャッシュ内部も含む136,668点、最大相対差4.71×10^-5で合格です。疎な点検と全伝播精度保証は区別しています。

## 保存先

- 公開する小さい成果：このディレクトリの `data/`、`figures/`、`README.md`。
- ローカル全時系列：`outputs/initial_dispersion_20261011/{initial_7d,extended_30d,combination_7d}/`。
- ローカル実行ログ・保全manifest：`outputs/initial_dispersion_20261011/`。

詳しい再現方法と入力前提は [READMEの再現節](README.md#l-変更ファイル再現git) を参照してください。保存済みデータがある環境では物理計算を重複実行しません。
