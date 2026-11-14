# 00_901_E487_最新接续卡｜2026-11-15

**E487为901目前最新完成阶段。** E486和E485旧结果原样冻结，不重复。V5.3与项目README已读；旧README首行E347过时，只看当前E487正式报告并参照原始研究合同。

## 真正新增的原件与科学事实
- 作者Julou et al 2020，PLOS Biology DOI `10.1371/journal.pbio.3000952`，Zenodo正式数据DOI `10.5281/zenodo.3894719`，源云端GitHub独立验收 https://github.com/dairongzhen3-creator/math/actions/runs/43616297733 。
- 官方111,442,020 B `GL_Preproc.tar.gz` MD5 `02c8ce66d2603b8dbd6211794696aedd`、SHA256 `4b63dfee10ad9839060186671f888178e10d3962a80481c752175fa37d23ccd6`。36个实验日期、1057个逐lane框架（729双通道12列、328单通道9列）、117,565细胞、3,014,419细胞帧、111,267真实可解析父指针、55,454严格认证母体双女儿分裂事件；另1个假阳性日期20190605 `parent_type=exit`且出生时间不连续，保留原件并排除。
- 遵作者R脚本 `discarded_dates` 排9日期200lane，剩27候选日期857lane、42,694严格分裂。`20180515 mislabelled 20180514`来源别名保留双字段。27个日期不等于同处理独立重复。
- 新观察对象为按日期分箱 `I(mother_last_fluorescence; daughters_first_sum_fluorescence | date)` 经验 **1.3362566424874205 bit/事件**；同lane同阶段错配20/40/80帧null中位0.930141/0.758736/0.510639。纯技术错配不能用于扣除因果信息。第二套不依赖主代码的标准库重算完全一致。
- 8片原字节已经全部上传到Google Drive**独立901大数据仓**并独立回读SHA通过，附恢复脚本/manifest/55,454逐事件记录。原111MB单件被Drive Library尺寸门禁阻止，故分8个15MB以内片段，完整官方MD5复合再次PASS。

## 继续E488的准确断点
先回原作者R`dataLoad.R`冻结不同日期/格式的真实荧光-尺寸-仪器通道、真实采样时钟，优先在**生物学独立date、同培养条件与同lane时间基线**比较观测信息量，防止分割及观测噪声、时间漂移及家系共享混杂。若未得到同步双曝光和分子绝对校准，`I(真实母分子;真实女子分子)`与物理因果bit、bit/s、bit/J、全生命周期去重总量一律NOT_IDENTIFIED。

正式报告：`901_E487_正式科研阶段报告_20261115.md`。
完整复现：`901_E487_母女身份链与真实荧光信息_完整研究过程包_20261115.zip`。
大数据永久源目录：https://drive.google.com/drive/folders/1TmXJZJzxgLHTrSv2uLPqSoapCp-O2jyH。
云端源码：https://github.com/dairongzhen3-creator/math/tree/cloud-901-e487-julou-source-20261115。

**论文状态**：主论文PUBLICATION_NO_GO；正式投稿级引文全文下载/PaperQA未完成，不投稿、不邮件作者、不重跑旧深度检索；E1–E477历史文件未逐一新全文读透，不冒充已全面审核。
