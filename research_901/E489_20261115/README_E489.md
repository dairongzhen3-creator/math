# 901 E489 authoritative cloud handoff

- [E489 project continuation](https://drive.google.com/file/d/1W979Ktwfxq_Pdia1gnyFaQ6Ht4EMO58F/view)
- [E489 report](https://drive.google.com/file/d/1Rft5AJg9_yaimvVwgK1BNfnklPRUE6xV/view)
- [E489 complete 25-member original source and QA zip](https://drive.google.com/file/d/1ud5OpEpmxesAVR54Cy-31NbdJxQni28z/view)
- [E489 data archive](https://drive.google.com/drive/folders/1UragVzjU5bf5vZ1MRugSqKAj42c8myxv)
- [GitHub Actions source acquisition receipt](https://github.com/dairongzhen3-creator/math/actions/runs/43704328062)

### Full stage report

# 901号 E489 正式科研阶段报告

## 公开FLCS单分子检测与后续诱导响应：真实数值、独立日期验证及物理联接边界

**研究日期：2026-11-15。** 当前项目为生命内部／外部信息数量和物理成本的实证量化；本期接续E488，新增可复核的**同细胞「切换前分子检测状态→切换后诱导快慢」**统计信息，不替换原有母女荧光信息账。**主研究论文仍为 `PUBLICATION_NO_GO`。**

## 一、原作者论文与公开原件的真实接收

原始论文：Julou T, Zweifel L, Blank D, Fiori A, van Nimwegen E (2020), *PLoS Biology* 18:e3000952, DOI https://doi.org/10.1371/journal.pbio.3000952。原作者代码仓库 https://github.com/julou/MoM_lacInduction/tree/plos ；对应标签 `plos` 的树SHA `976f62f31433d315ae2337a72d99f6deba4165ee`。FLCS分析脚本为 `src/MoM_lacInduction_FLIM.Rmd` Git blob `bfc059418ebcc88f326c9ba722abb1a6d5ef52fe`。GitHub Actions真实源回执 https://github.com/dairongzhen3-creator/math/actions/runs/43704328062 。

本阶段**真正取得7份已发表的FLCS分析原件**（2份XLSX、2份CSV和3份DAT），逐件从作者Git LFS实体下载，按Git LFS正式公布的原字节SHA-256和尺寸核验7/7 PASS；XLSX ZIP CRC通过。原件及逐件校验收据保留于本期复现ZIP。未取得的原始 `.ptu` 单光子到达记录不得冒充已下载；作者给的FLCS估计表属于仪器数据的**分析后原件**。所有源哈希及链接见 `source_audit.json` 和 `E489_artifact_manifest.json`。

方法语义来自作者Methods及公开Rmd：采用荧光寿命相关光谱（FLCS），用GFP与细胞自发荧光的寿命差分离信号；用拟合扩散时间、相对单分子光子计数以及检测体积内的有效荧光粒子数作检测判据。原论文的阈值包括CPM信号5σ、3–30ms的扩散窗口及对超过20个粒子的异常拟合值限制。这里的粒子数是**探测体积内拟合值**，通常可能以LacZ-GFP四聚体为单位，不等于全细胞绝对蛋白拷贝数，亦非跨分裂实际分子转运。

本次原件比历史E487源的区别是实质性的：E487的日期来自**2015-12-04至2019-06-14**，在Nikon Ti-E常规宽场显微和MoMA处理的母机实验中追踪母女。FLCS工作簿记录**2020-01-22、2020-01-23、2020-01-28**三批点式FLCS微流控实验；原论文的FLCS由Olympus IX73、光子计数与荧光寿命成像流程完成。按物理实验日期，`E487 ∩ E489 = 0`，且E489的FCS文件名／位置未包含E487的母体ID、女儿ID和分裂帧联合键。因此**两套原件不存在可认证的同细胞母女物理联表**，不能把FLCS单分子阈值直接归因于E488已经计算的42,658次母女分裂事件。

## 二、逐源键匹配，明确保留测量歧义

`FLCS_lacInduction.xlsx` → `MotherMachine`表有74条真实记录；`FLCS_lags.xlsx` → `Sheet1`表亦有74条。以完整`(date, file)`严格连接：**73条精确匹配**。2020-01-23有1条源文件名冲突：FLCS `FCS_Ch9_bottom_2.ptu`，延迟源 `FCS_Ch9_bottom_1.ptu`。同日近似名称不是同一原始测量的有效身份证明，本轮不给其人工别名、不入主计算；两条原源名完整保留。

精确连接后，59条具有作者S（短）或L（长）响应类别及0、1或-1的FLCS检测标记；作者`gfp=-1`表示**dubious / 质量有疑的荧光判据**，此类共3条（S一条、L两条），不擅自转成0或1。故主估计使用56条确定0/1检测标记、S/L结局均明确的细胞，覆盖**3个实验日期**，各22、15、19条。其余源行、异常日期、`NA`/非生长值和疑似重复文件名保留在逐行联接表，不从「66 cells」的发表叙述反推未在源表出现的记录。

## 三、真正新增的逐细胞信息数量

定义 `X=1` 为母机FLCS在乳糖切换前检测到LacZ-GFP，`X=0` 为未检到；`Y=S/L` 为同一细胞切换后的短／长诱导延迟；分类直接采用**作者原始标记**，不利用后验最优阈值重新调节。

| FLCS切换前状态 | 短延迟 S | 长延迟 L | 合计 |
| --- | ---: | ---: | ---: |
| 检测到 LacZ-GFP（1） | **13** | **1** | **14** |
| 未检测到 LacZ-GFP（0） | **8** | **34** | **42** |
| 总计 | **21** | **35** | **56** |

以真实56条构造经验联合分布，计算

\[
\widehat I(X;Y)=\sum_{x,y}\widehat p(x,y)\log_2\frac{\widehat p(x,y)}{\widehat p(x)\widehat p(y)}
=\boxed{0.3347760077969534\ \mathrm{bit/valid\ cell}}.
\]

该子样本内 `P(short|detected)=13/14=92.8571%`，`P(short|undetected)=8/42=19.0476%`，差**73.8095个百分点**；结局本身的经验熵为**0.9544340029 bit/细胞**，观察到FLCS状态使该样本二分类结局的经验不确定性减少约**35.08%**。它是**状态—后续响应的统计信息数量**，并非细胞间消息、外界刺激随机编码通道容量、真实蛋白拷贝数传递，亦非bit/s、bit/J。

### 跨实验日期留出（预先固定Laplace-1解码器）

每次整实验日期留出，将该日期全部细胞保留作为测试；仅以其余两个日期训练类别无关基线`P(S)`，以及基于FLCS检测0/1的`P(S|X)`，均对二元比例施加固定Laplace(1,1)平滑。对每个测试细胞计算 `log2[P_true(S/L|X)/P_base_true(S/L)]`，没有测试集调参。

| 留出真实实验日期 | 细胞n | 折外对数得分新增（bit/细胞） |
| --- | ---: | ---: |
| 2020-01-22 | 22 | +0.281251106 |
| 2020-01-23 | 15 | +0.505451984 |
| 2020-01-28 | 19 | +0.216645351 |
| **全部56条加权** | **56** | **+0.319385103** |

3/3日期均获得正增益，说明在同实验室的这三批中该关联有折外预测价值。**只有三批**，且可能同培养通道记录不独立，不能将个体56条视作56个独立生物学重复；也不能将三个折正向直接声称跨实验室泛化。

技术诊断：在「各日期检测阳性总数及短延迟总数均固定、细胞标签在日期内可交换」这一**额外技术零假设**下，三个日期超几何分布卷积的观测阳性短延迟数量13的单侧条件参考尾为 `16/5,496,491 ≈ 2.91×10⁻⁶`。该参考并非真实随机分配FLCS状态的因果检验，也不纠正同一微流控通道内生物依赖；因此不作为主研究的生物学P值使用。

## 四、关键反事实与质量门

**三条可疑检测的穷举**：若限定本次59条S/L有效响应个体，允许3个原`-1` FLCS状态在真实二元0/1间任意赋值，共8种可能；其经验 `I(X;Y)`位于**[0.2246262843,0.3528913319] bit/59条样本的每一观察细胞**。这只是固定59条已观察反应并将3个不确定检测标签穷尽分配的**条件识别集合**，不是抽样置信区间，更不是分子传输的上下界。不能把59条区间和56条主估计误作同一分母的误差带。

**独立快照对照**：同一FLCS分析工作簿的`Snapshot`表中，原作者记录150个可确认0/1的单细胞：检测到43、未检测到107，检测比例43/150=**28.6667%**；另有3个`-1`未纳入。它与作者论文报告的约28%检测率相容，但为**不同液体培养快照实验**，不能与3批母机FLCS响应个体合并为一个联合分布。

**作者已发表贡献与本次贡献分开**：原作者首先指出少量预存LacZ/LacY与短延迟之间的关系，并以FLCS论证阈值接近单分子量级；本阶段新增的是针对公开逐细胞原件的精确源键匹配、三日期保留验证、0.334776-bit经验信息量、0.319385-bit折外对数得分以及可疑标签的完全枚举界，不宣称原来的生物机理是本次发现。

## 五、独立QA、文件保存与下一期精确任务

- `extract_flcs_originals.py`经`artifact_tool`逐表读取作者原XLSX并输出全部结构保留CSV；`analyze_e489_flcs.py`重建73条精确联表、观察MI、日期留出、三条可疑标记穷举及日期内固定边缘技术参考；`independent_qa_e489.py`**从两份正式XLSX重新读取**，不调用主MI和留出代码，用60位Decimal高精度复算MI与另一实现的全部56条留出得分，PASS。两算法使用相同的XLSX解码器，不冒称独立测序或另一个实验室。
- E487冻结的36日期源身份表原样保存在包中，仅用于日期交集和实验来源门验证；不重训E488的42,658次母女研究、不改变任何旧阴性结果。E489没有付费资源、没有湿实验或对外邮件。
- 最小下一研究门（E490）：获取与**同一2020 FLCS细胞**精确配对的`.ptu`光子到达记录、同态重复FLCS及有效背景/体积标定、真实膜片/微流控培养通道和`cell_id→parent_id→daughter_ids`逐分裂联接，且能判断条件施加时钟。若这些证据未出现，当前真正可报告的仍是**切换前分子检测状态对切换后诱导结果的预测信息量**，不是细胞内部全部分子信息或自然跨细胞通信总量。
- 投稿前仍需逐篇完整原文下载、逐条引文比对及PaperQA正式记录；本轮是可复现实证阶段，不触发投稿。

**科学状态：** `E489_7_LFS_REAL_ORIGINALS_SHA_PASS / 74_AND_74_SHEETS_73_EXACT_KEYS / ONE_RAW_ID_MISMATCH_EXCLUDED / 56_VALID_CELLS_3_EXPERIMENT_DATES / OBSERVED_FLCS_TO_INDUCTION_MI_0p334776 / DATE_HELDOUT_LOG_GAIN_0p319385 / BINARY_DUBIOUS_SET_REPORTED / CROSS_TO_E488_SAME_CELL_JOIN_ZERO / MOLECULAR_INHERITANCE_BIT_NOT_IDENTIFIED / MAIN_PUBLICATION_NO_GO`。
