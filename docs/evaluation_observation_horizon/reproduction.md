# 再現手順

## 環境とテスト

Python>=3.11、NumPy>=2、matplotlib、PyYAML、pytestを使用。新規依存なし。今回の環境は前回の `.analysis-venv` を再利用した。

```powershell
python -m pytest tests/test_evaluation.py tests/test_regression.py -q
python -m pytest -q
```

## 24時間の3ケース再実行・回帰

```powershell
python scripts/check_baseline_regression.py --output outputs/regression_new
```

固定fixtureのconfigを使用し、制御＋外乱、制御＋無外乱、全LOW＋外乱を1回ずつ計算。新しいsuiteとregression.jsonを保存する。通常のpytestは毎回この長い計算を行わない。すでに再実行したsuiteを二重計算せず比較する場合：

```powershell
python scripts/check_baseline_regression.py --existing-suite outputs/regression_new/suite --output outputs/regression_recheck_new
```

失敗時は非ゼロ終了。参照値を自動更新する機能は設けていない。

## 保存データの後処理

```powershell
python scripts/evaluate_saved_case.py --run outputs/regression_new/suite/controlled_disturbed --all-low outputs/regression_new/suite/all_low_drag_baseline --output outputs/evaluation_new --epsilon-e 50 1 5 10 --epsilon-v 0.001 --hold-window 1800
```

新しいディレクトリに後方互換metrics.json/metrics.csv、閾値ごとのsettling_sensitivity.json、元runと設定を指すprovenance.jsonを保存。元データは読み取り専用。省略した判定値は元config.recoveryから取得。

長い保存historyから短い観測期間を評価する場合は `--horizon-hours 24` 等を指定。最後に含めた保存サンプルの時刻を実際のTmaxとして保存する。範囲外のhorizon、物理条件や外乱有無が異なるLOW基準は拒否する。

## 段階的な観測時間

24 hの回帰と72 h保存データを確認後、必要なら次へ進む。

```powershell
python scripts/run_observation_case.py --hours 168 --output outputs/central_7d_new
python scripts/run_observation_case.py --hours 168 --case all_low_drag_baseline --output outputs/low_7d_new
python scripts/evaluate_saved_case.py --run outputs/central_7d_new --all-low outputs/low_7d_new --output outputs/evaluation_7d_new --epsilon-e 50 1 5 10
```

7日結果で再逸脱・未完了・新規切替があれば336 hへ延長する。現行coreはcheckpoint再開を持たないので延長はt=0から積分する。coreを変更して再開機能を作る代わりに、既存24/72 hを再利用し、必要な7/14 dを各条件で1回だけ計算する。保存済み長時間runがあれば後処理だけで再評価する。

## 外乱位置

```powershell
python scripts/run_observation_case.py --hours 168 --disturbance-satellite 6 --output outputs/mid_7d_new
python scripts/run_observation_case.py --hours 168 --disturbance-satellite 1 --output outputs/end_7d_new
```

他の物理・制御条件は固定し、既存Gaussian中心を沿軌道座標へ変換して位置だけ変更する。中央対称ケースの結論を非対称ケースに一般化しない。LOWコスト比較が必要なら同じ位置指定で `--case all_low_drag_baseline` を計算する。

## 出力とデータ保護

今回の新規数値historyは `outputs/evaluation_20261007/` に保存。既存runs/、前回outputs/を上書きしていない。既存299ファイルのSHA256を確認し、検証結果を成果物へ保存した。

Gitには固定回帰fixtureと小さい集計・図・文書を保存し、長時間の大きいNPZは既存のoutputs除外規則を維持する。集計成果と元runの対応は結果manifest/provenanceから確認できる。再現用コマンドの出力先は、実行ごとに未使用の名前にする。
