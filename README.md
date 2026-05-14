# medicalAsset
**Clinical Special Disease Data Asset Value Assessment R Package**

## English Introduction

medicalAsset is an open-source R package for value assessment of clinical disease-specific data assets. Built on a compliance-driven and machine learning-enhanced multi-dimensional framework, it provides standardized, reproducible, one-click valuation for clinical data assets.
Features
✅ Multi-disease support (diabetes, cardiovascular, oncology)
✅ Full-lifecycle compliance review
✅ Multi-algorithm intelligent feature selection
✅ Five-dimensional scoring
✅ Comprehensive score & grade (A/B/C/D)
✅ Four professional visualizations


    ## 中文介绍
    `medicalAsset` 是一款面向**临床专科疾病数据资产**价值评估的开源 R 工具包。本工具基于**合规驱动 + 机器学习增强**的多维评估框架，集成数据治理、全生命周期合规审查、智能指标筛选与多维价值量化全流程，实现临床专科数据资产的**标准化、可重复、一键式价值评估**。

### 功能
- ✅ 糖尿病、心血管病、肿瘤等专科数据资产评估
- ✅ 全生命周期合规审查
- ✅ GBDT/RFE/VIF/互信息多算法特征筛选
- ✅ 五维评分：数据质量、成本、应用价值、安全合规、指标敏感性
- ✅ 综合评分 + 资产等级（A/B/C/D）
- ✅ 4种专业可视化：柱状图、雷达图、特征图、交互式仪表盘

### 安装
```r
# 安装
devtools::install_github("supersicz/medicalAsset")

# 加载
library(medicalAsset)
