# 数学成果审读 / Independent reading ledger

开始日期：2026-10-07。此目录是仓库所有者的独立阅读层，不是 OpenAI 官方认证，也不是原作者的勘误公告。

## 固定来源

- 阅读基线 commit：`adc7f1241b42e322a6451854ab7e4b4c146bf78a`。
- 原始总目录：[CONTENTS.md](../CONTENTS.md)。该版本声明 722 篇手稿、372 个成果家族。
- 总目录 Git blob：`8148d2c9077e611d5665a2ded93da68e5b608f0f`（640371 bytes）。
- 学科与家族概览：[overview.tex](../overview.tex)，Git blob `fc40ebdcd0fb9db026a2e44bde66f2e469939fd4`。
- 单篇正文与代码：[preprints](../preprints/)。形式化材料：[lean](../lean/)。

## 阅读状态不相互替代

`UNREAD` 未逐篇阅读；`ABSTRACT_READ` 已读该篇目录摘要；`PARTIAL_TEXT_READ` 已读明确列出的正文范围；`FULL_TEXT_READ` 已读所有正文及其内附证明；`PROOF_CHECKED_PARTIAL` 仅复核列出的步骤。形式化文件存在、成功编译、检查声明与论文同义，是三个独立字段。

**取得目录、下载文件、自动关键词分类，都不计作正文已读。** 主稿引用的 companion 证明未读时，必须单列依赖缺口。没有阅读记录的条目默认 `UNREAD`，而非默认无关。

## 原始材料保护

本次所有新增文件只放在 `our_review/`；不改原有论文、源码、LICENSE 或官方目录。不运行来源不明的构建脚本，不宣称已经验收任何未复核的证明。

本仓库公开。本目录仅保存已公开文献的阅读记录与通用数学问题标签；未公开论文编号映射、内部研究索引、投稿状态与内部评审不上传到这里。

## 文件入口

- `00_722_master_index.md`：完整目录的阅读入口与当前完成范围。
- `review_state.json`：按精确论文路径记录实际阅读进度。
- `build_review_index.py`：从本仓库的完整 CONTENTS.md 生成逐篇台账；读取默认未读状态，不将解析计为人工阅读。
- `02_high_priority_results.md`：已读材料支持的初步优先队列，不是已确认可迁移定理清单。
- `notes/`：注明原文位置、前提、结论、依赖和核查状态的阅读札记。

所有结论都应可回到固定版本的原文。每次更新追加具体阅读范围，不以对话里的旧评价代替本轮证据。
