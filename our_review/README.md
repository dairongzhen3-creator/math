# 数学成果审读 / Independent reading ledger

开始与本轮更新日期：2026-10-07。本目录是仓库所有者的独立阅读层，不是 OpenAI 官方认证，也不是原作者的勘误公告。

## 本轮完成情况

- **722篇目录摘要全部逐篇审读**：按37个摘要包顺序读取；被工具截断的尾段另外取回。覆盖全部372个家族。
- **4篇正文部分已读**：PEPS、无噪高斯回归记忆界、熵光子数不等式、complete Crouzeix。每篇实际阅读范围见状态表及札记。
- 完整论文已读0；整篇证明验收0；Lean独立构建0。摘要已读与正文部分已读是重叠范围，不把4篇另外加到722之上。
- 选出[25个优先家族](02_high_priority_results.md)，其中10个先核；不因此判其他家族无关。

## 直接入口

- [722篇逐篇台账](00_722_master_index.md)：精确原文路径、固定版本、学科、自动候选标签、实际阅读状态。
- [722篇CSV台账](722_master_index.csv)及[完整JSONL记录](catalogue_rows.jsonl)。
- [37个摘要阅读包](722_abstract_batches.md)与[家族阅读包](372_family_batches.md)。
- [本轮优先队列](02_high_priority_results.md)。
- 正文札记：[PEPS](notes/peps_representation_scope.md)、[有限记忆回归](notes/gaussian_regression_memory_scope.md)、[熵光子数](notes/entropy_photon_number_interface.md)、[非正规函数误差](notes/crouzeix_nonnormal_error_interface.md)。
- [人工阅读状态](review_state.json)和[目录一致性检查](index_audit.json)。

## 固定来源

基线commit：`adc7f1241b42e322a6451854ab7e4b4c146bf78a`。

原始目录[CONTENTS.md](../CONTENTS.md)的Git blob为`8148d2c9077e611d5665a2ded93da68e5b608f0f`，640371 bytes；overview.tex的blob为`fc40ebdcd0fb9db026a2e44bde66f2e469939fd4`。脚本核到722条记录、722个不同PDF路径、372个家族和722份非空摘要。

原始家族编号最高到377，空号045、061、070、123、163。保留这些原始编号，不将数量372误写成连续编号001–372。单篇材料在[preprints](../preprints/)，形式化材料在[lean](../lean/)。目录中有形式化链接不等于本轮已编译验证。

## 阅读状态不相互替代

`UNREAD`：没有逐篇阅读记录。

`ABSTRACT_READ`：已读该篇目录摘要，不代表正文已读。

`PARTIAL_TEXT_READ`：只读状态表注明的正文范围。

`FULL_TEXT_READ`：正文及内附证明全部读完，但不自动表示独立证明验收。

`PROOF_CHECKED_PARTIAL`：只复核明确列出的证明步骤。

形式化文件存在、成功编译、声明与论文同义，分别记录。取得目录、下载文件、自动关键词分类都不计作正文阅读。主稿引用的companion证明未读时仍是依赖缺口；一个家族中的多篇不自动等于多次独立证据。

## 原始材料与隐私保护

新增阅读材料只放在`our_review/`。唯一额外文件是[目录生成工作流](../.github/workflows/build-math-review-index.yml)：它只执行本次自写并固定版本的目录解析器，不运行任何原论文构建或证明代码。原有论文、源码、LICENSE和官方目录保持不变。

本仓库公开。这里只保存已公开文献的阅读札记与通用数学标签；未公开论文编号映射、内部研究索引、投稿状态与内部评审不上传到这里。

## 复现台账

`python our_review/build_review_index.py --root .` 从完整目录和人工状态生成台账。自动候选标签只用于导航，不冒充语义审读。

GitHub工作流只在手动触发或工作流文件被明确修改时运行；不是后台持续阅读任务。更新人工阅读状态后需要重新生成索引。目录一致性检查与解析器测试通过，也不代表数学证明通过。

所有结论须回到固定版本原文。每次新增实际阅读范围，不以旧对话印象替代本轮证据。
