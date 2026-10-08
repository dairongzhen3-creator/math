# 数学成果审读 / Independent reading ledger

开始日期：2026-10-07。最新更新：2026-10-08。本目录是仓库所有者的独立阅读层，不是 OpenAI 官方认证，也不是原作者的勘误公告。

## 当前可核对的完成范围

- 2026-10-07的历史台账登记722篇目录摘要已读；本次恢复未重新逐篇重放该摘要阅读。
- **F325的两篇正文、内附证明、附录及文献表已按源码通读；两篇各一幅提供的页面图已查看。** 完整记录见 [F325全文阅读与建议](notes/F325_full_reading_20261008.md)。
- 当前有可交付全文源码读后记录的候选家族：**1/25（F325，共2篇）**。其余24个家族没有在此登记为全部完成。
- PEPS、无噪高斯回归和熵光子数三篇暂保留最后存档的部分阅读状态。中断运行中曾报告更进一步的阅读，但没有完整可复核记录的，不借该报告升级。
- 全局证明独立验收0；Lean独立构建0。源码通读与完整PDF逐页核对分开，当前两篇没有做全PDF/源码一致性核查。

**阅读状态的唯一更新入口是 [review_state.json](review_state.json)。** 自动生成的总表、CSV、JSONL和index_audit.json目前仍是2026-10-07快照，尚未按此次状态更新重新生成；其中旧数量不代表本页的最新进度。

## 直接入口

- [领域地图与全文阅读顺序](03_domain_map_and_reading_order.md)。
- [全文阅读后的建议模板](04_post_reading_record_template.md)。
- [F325两篇完整源码读后记录](notes/F325_full_reading_20261008.md)。
- [25个优先家族](02_high_priority_results.md)：这是摘要筛选后的候选队列，不是可直接使用的定理清单。
- [722篇逐篇台账](00_722_master_index.md)、[CSV台账](722_master_index.csv)、[JSONL记录](catalogue_rows.jsonl)：目录快照和导航。
- [37个摘要包](722_abstract_batches.md)与[家族包](372_family_batches.md)。
- 较早的部分阅读札记：[PEPS](notes/peps_representation_scope.md)、[有限记忆回归](notes/gaussian_regression_memory_scope.md)、[熵光子数](notes/entropy_photon_number_interface.md)、[Crouzeix初读](notes/crouzeix_nonnormal_error_interface.md)。Crouzeix最新阅读范围以2026-10-08完整记录为准，初读笔记保留为历史。

## 固定来源

原始论文基线commit：`adc7f1241b42e322a6451854ab7e4b4c146bf78a`。

原始目录[CONTENTS.md](../CONTENTS.md)的Git blob为`8148d2c9077e611d5665a2ded93da68e5b608f0f`，640371 bytes；overview.tex的blob为`fc40ebdcd0fb9db026a2e44bde66f2e469939fd4`。既有目录检查记录为722个不同PDF路径、372个家族和722份非空摘要。

原始家族编号最高到377，空号045、061、070、123、163。保留原始编号，不将数量372误写成连续编号001–372。单篇材料在[preprints](../preprints/)，形式化材料在[lean](../lean/)。存在形式化链接不等于本轮已编译验证。

## 阅读状态不相互替代

`UNREAD`：没有逐篇阅读记录。

`ABSTRACT_READ`：已读该篇目录摘要，不代表正文已读。

`PARTIAL_TEXT_READ`：只读状态表注明的正文范围。

`FULL_TEXT_READ`：正文及内附证明全部读完，须标明阅读介质、图表覆盖与未解决的外部依赖；不自动表示独立证明验收。

`PROOF_CHECKED_PARTIAL`：只复核明确列出的证明步骤。全文状态与proof_verified字段分别记录。

取得目录、下载文件、自动分类都不计作正文阅读。主稿引用的companion证明未读时仍是依赖缺口；同一家族的多篇不自动等于多次独立验证。

## 原始材料与隐私保护

新增阅读材料只放在`our_review/`。已有[目录生成工作流](../.github/workflows/build-math-review-index.yml)只执行固定版本的目录解析器，不运行原论文的构建或证明代码。原始论文、源码、LICENSE和官方目录保持不变。

本仓库公开。这里只保存公开文献的阅读札记与通用数学标签。未公开论文编号映射、内部研究索引、投稿状态与内部评审不上传到这里。

## 复现台账

`python our_review/build_review_index.py --root .` 从完整目录和人工状态生成台账。自动标签仅用于导航，不冒充语义审读。

GitHub工作流不是后台持续阅读任务。目录解析器或局部校准通过，不代表数学证明通过。所有结论须回到固定版本原文。
