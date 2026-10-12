# 文献調査の確認記録

確認日：2026-10-12。branch：`astra/formation-input-literature`。
作業開始点：`8bd03b050cb6c4de2a426c9acce650a28d96a1be`（既存分散制御実験）。

## 実施した確認

1. 原文20資料の本文・表・図説明・公式運用資料を確認。該当箇所と数値の種類を `evidence.json` に保存。
2. A-Trainの12.5±2.5秒をPDFページ画像で照合。OCRのみでは符号・小数が不明瞭だった。
3. `audit.py` の検算が合格。基準B、逆数定義、質量・幾何面積換算、隣接差の意味、時間/距離換算を確認。
4. `.analysis-venv/Scripts/python.exe -X utf8 -m pytest -q`：**273 passed, 1 warning, 35.64 s**。
5. `git diff --exit-code HEAD -- src configs tests scripts docs/dispersion_control_transfer docs/initial_dispersion_screening`：差分なし。

警告は `test_analysis_cli_roundtrip_and_output_protection` 中の子プロセス出力読み取りスレッドの `UnicodeDecodeError`。UTF-8として読んだ出力に別エンコードのバイトが含まれる警告で、テストは合格した。今回の文献調査と無関係な既存実装は修正していない。警告がなかったとは報告しない。

## 既存baseline・runの扱い

この作業は文献調査と換算だけで、物理runの再実行はゼロ。24時間software reference、既存長時間run、制御・物理core、設定を編集していない。今回独立に24時間全ケースの物理runを再計算したという主張はしない。既存の回帰テストを含む全テスト合格と、基準commitからの既存ファイル差分ゼロを確認した。

新規成果は `docs/formation_input_evidence/` に限定。論文キャッシュはGit除外の `outputs/formation_input_literature/` に保存した。原論文の全文・大容量PDF・既存NPZはcommitしない。

## 再現

リポジトリルートで以下を実行する。

```text
.analysis-venv/Scripts/python.exe -X utf8 docs/formation_input_evidence/audit.py
.analysis-venv/Scripts/python.exe -X utf8 -m pytest -q
```

`audit.py --output NEW_PATH.json` は新規ファイルだけを書き、既存ファイルの置換を拒否する。新たなシミュレーションは起動しない。利用する定数は既存 `configs/baseline.yaml` から読む。

根拠再確認は `evidence.json` のURLと該当節・表番号を使う。A5/A6/B5は出版社・機関URLに加え、報告書文献表に示した著者公開全文を実際の本文確認に利用した。PDFの全再ダウンロードを必須としない。

## 成果物

- `README.md`：約2,355文字（Markdown・リンクを含む）の研究判断用サマリーと技術付録。
- `evidence.json`：20資料の書誌・数値・種類・限界・確認箇所。
- `recommended_conditions.json`：暫定3水準。実行指示・実機代表分布ではない。
- `audit.py`：定義と換算の再現。
- `audit_result.json`：確認時の数値。
- この確認記録。

mainへのmergeは実施しない。commit hashはGit履歴と作業完了メッセージで確認できる。
