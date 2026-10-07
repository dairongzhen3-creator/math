# Polynomial PEPS approximation：正文条件初审

日期：2026-10-07。状态：PARTIAL_TEXT_READ。没有验收整篇证明。

## 固定原文与实际阅读范围

基线 `adc7f1241b42e322a6451854ab7e4b4c146bf78a`。

论文：[Polynomial PEPS approximation of gapped square-grid ground states](../../preprints/Polynomial-PEPS-approximation-of-gapped-square-grid-ground-states-September-24-2026/paper.pdf)。

已读源码：
- [00-introduction.tex](../../preprints/Polynomial-PEPS-approximation-of-gapped-square-grid-ground-states-September-24-2026/build/sections/00-introduction.tex)，全文件，blob `0f6f4dcdd46d9e764b4f1e52718bf330cdc4313a`。
- [08-companion-interface.tex](../../preprints/Polynomial-PEPS-approximation-of-gapped-square-grid-ground-states-September-24-2026/build/sections/08-companion-interface.tex)，全文件，blob `f9ee6dd9ea8a4b4bde5618010f474f9a170390fe`。

未完成：其余主证明、area-law companion 正文、PDF与源码一致性、任何形式化构建。

## 原文定理声明

固定局域维数 q≥2、相互作用范数上界 J、全系统谱隙下界 Δ>0。模型是开放 L×L 正方格上的 Hermitian 单点项与最近邻项之和，各项范数≤J。要求唯一归一化基态 Ω，并满足全系统不等式 H−E0I≥Δ(I−|Ω><Ω|)。

原文声称存在同一格上的非零 PEPS Φ，最大键维数≤C L^c，归一化后到 Ω（允许整体相位）的向量范数误差≤L^−1。C、c 可依赖 q、J、Δ，但不依赖单个格长和模型实例。

原文不要求子区域限制 Hamiltonian 都有隙，也不要求 frustration free、平移不变、项彼此对易或从积态出发的有隙路径。

**关键：原文明确不声称张量能被高效找到或高效收缩。** 证明允许一方起初持有完整 Ω，因为目标是表示存在性，不是寻找基态的算法。

## 配套论证的依赖

附录显式使用 *A two-dimensional area law from a global spectral gap* 的面积律、正算子替换、collar bound、amplification 和 two-family entropy cancellation 接口。引言随后说明本文另做 tapering、patch encoding、density compression 与 routing。

因此，读完本篇的定理声明，不等于验证了面积律前提，也不等于独立核实整条表示构造。配套论文是承重依赖，不是另一篇独立重复证据。

附录最后特别区分：不能仅对普通 mutual-information 不等式增加条件变量，就推得所需 conditional-information 界；相应局部变体需要另证。

## 迁移前必须检查（审读者判断）

1. 目标模型能否在不改变问题的情况下满足开放方格、有限局域维数、局域作用与唯一基态条件？
2. 使用的 q、J、Δ 是否在所讨论的尺寸或截断极限中固定？每个有限截断有某个正隙，不等于共同的正下界；截断使局域维数或项范数增长时，不能忽略 C、c 对它们的依赖。
3. 目标是基态表示、局域观测量、激发谱还是非正规生成元的弛豫率？定理只直接声明第一项，跨对象误差传递需另证。
4. 不能以“多项式键维数”代替可实施算法和实际计算预算。

## 一个条件推论（本札记的初等推导，不是论文的新颖性认领）

若归一化态 ψ、φ 在选定整体相位后满足 ||ψ−φ||≤δ，则对任意有界 O，分两项并用 Cauchy–Schwarz 得

|<ψ,Oψ>−<φ,Oφ>|≤2||O||δ。

因此，**条件于论文存在性定理为真**，存在性结论可给有界局域观测的相应误差界；它仍不告诉我们如何找出该态，也不自动给出激发能或谱隙估计。

结论：这是重要的表示论候选工具，当前不是已验收的计算加速方案。
