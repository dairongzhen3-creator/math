# 901 E493 cloud continuation | 2026-11-15

Source is the immutable 901 project report and independent cloud data archives.

- [Formal Chinese report](https://drive.google.com/file/d/1nXi-Antyt2e_pa9BeFSRQPzz-7nhgp1M/view)
- [Frozen current continuation](https://drive.google.com/file/d/1pcTXY0uL5gH3UknjdLu29U46JISsUwP9/view)
- [Full 15-member original-data, code, QA archive](https://drive.google.com/file/d/1XhgmO1uT2VnQuIETE8j4QIdxRaVwddtN/view)
- [E493 independent big-data folder including original 842092-byte Table S2](https://drive.google.com/drive/folders/1xy0MZkebu4DuwtEK4FW5sY-cgwNQy6Bc)
- [Actual source acquisition, verified publisher MD5](https://github.com/dairongzhen3-creator/math/actions/runs/43712125191)

## Frozen official continuation

# 901号 E493 最新接续卡（2026-11-15）

**总目标：** 定量记录生命内部、外部真实信息和时间/能耗；不将概率预测分数、细胞荧光均值、蛋白稳定性、真正分子转运混为一体。E1–E492既有结果、原件和阴性冻结。研究状态`PUBLICATION_NO_GO`。

## 已取得并完成（不要重复计算）

1. 主论文 Auboiron et al 2021 *iScience* DOI `10.1016/j.isci.2021.102075`，正式补充 `mmc2.xlsx` Table S2原件（842092 bytes；MD5 `4d2104fd2a40a51bf21b8dddfc0f3428`；SHA256 `2ea57bedff6ccdd71ea982e695eea638cb6a815a1d2c8347b7ae91d8d8c21893`）。原件共**1338蛋白ID**；含CBF1与HHO1两种菌株背景，0–180min两背景SiLAC H/L及作者拟合`kd/half life`。云端Github Actions成功: https://github.com/dairongzhen3-creator/math/actions/runs/43712125191。
2. 正正配对`n=1335`。CBF1/HHO1中位半衰期`86.506734 / 60.940944 min`，Spearman `0.31801614`；83.5206%的蛋白在HHO1更短；93.1086%两背景差值≤2倍；HHO1/CBF1中位比0.693658。3个负拟合速率的蛋白：CBF1 DBP8、HHO1 LRP1/GYP5，原源记录保留而不并入正半衰期。
3. 全部1338行`kd*T_half≈ln2`PASS；作者源表`fitted_slope-kd`每背景恒定：CBF1 `0.002539000661392/min`，HHO1 `0.002423591540419/min`。数值证明修正口径一致，不将其自动解释为单细胞生长率。
4. 真正同一母体及两女儿的2021 STAR Protocols H4-Dendra2 Table1共145条，在第二女儿出现后有**34帧母+两女同步观测**。新同帧原始红绿光谱对照 `C=log2[(R_m/G_m)/(R_d/G_d)]`：母–女儿1中位+0.289519，32/34为正；母–女儿2中位+0.306796，31/34为正。第一女儿出生后至第二女儿前16帧中位+0.050540。总84母女同帧配对是**1个母系多次采样**，不是84独立实验。光学对照不等于分子数比。
5. 数学识别门：即使把背景、增益都设为有利的零/一致，只有原始每细胞`FL mean`而无`A_m,A_d`，任何母体总荧光份额`A_m I_m/(A_m I_m+A_d I_d)`均可由未知面积比对应；无法对真实分子份额给非平凡边界。实际还缺ROI面积、背景、光转换效率和重复。
6. 两条独立程序：主`artifact_tool`原XLSX读取；QA使用纯Python OOXML原生ZIP+XML重读1338行，并用Decimal48重读STAR Protocols XML。结果**PASS**。完成原件归档、官方MD5/SHA和代码清单。

## E494优先动作

找真正已开放且有`physical_mother_id`,`daughter_ids`,`division_time`,`old/new protein molecules before/after`,`same-state technical repeats`,`optical calibration`,`culture batch`的实验原件。TrIPP主论文Data Availability称显微图像与Perl/R程序须向通讯作者请求，尚未取得18蛋白完整BudJ多母系CSV；原作者已发表群体图表不作个体真实原件。若无新来源，停止同数据调阈值与换分箱。**没有邮件授权，不发作者邮件；不投稿、不开付费云**。

详见本期正式报告、完整复现ZIP、SiLAC原始Excel、source qualification CSV和独立GitHub分支。`PUBLICATION_NO_GO`。
